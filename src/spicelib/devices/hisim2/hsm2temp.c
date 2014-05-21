/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2temp.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "hsm2def.h"
#include "hsm2evalenv.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

#define  Nsubmin     (1e15 / C_cm2m_p3) 
#define  Nsubmin_dlt (0.01 / C_cm2m_p3)
#define  lpext_dlt   (1e-8 / C_m2cm) 

#define RANGECHECK(param, min, max, pname)          \
  if ( model->HSM2_coerrrep && ((param) < (min) || (param) > (max)) ) { \
    printf("warning: (%s = %g) range [%g , %g].\n", \
           (pname), (param), (min), (max) );  \
  }
 
/*---------------------------------------------------*
* smoothZero: flooring to zero.
*      y = 0.5 ( x + sqrt( x^2 + 4 delta^2 ) )
*-----------------*/

#define Fn_SZtemp( y , x , delta ) { \
    T1 = sqrt ( ( x ) *  ( x ) + 4.0 * ( delta ) * ( delta) ) ; \
    y = 0.5 * ( ( x ) + T1 ) ; \
  }

#define Fn_SUtemp( y , x , xmax , delta ) { \
    T1 = ( xmax ) - ( x ) - ( delta ) ; \
    T2 = 4.0 * ( xmax ) * ( delta) ; \
    T2 = T2 > 0.0 ?  T2 : - ( T2 ) ; \
    T2 = sqrt ( T1 * T1 + T2 ) ; \
    y = ( xmax ) - 0.5 * ( T1 + T2 ) ; \
  }

#define Fn_SLtemp( y , x , xmin , delta ) { \
    T1 = ( x ) - ( xmin ) - ( delta ) ; \
    T2 = 4.0 * ( xmin ) * ( delta ) ; \
    T2 = T2 > 0.0 ?  T2 : - ( T2 ) ; \
    T2 = sqrt ( T1 * T1 + T2 ) ; \
    y = ( xmin ) + 0.5 * ( T1 + T2 ) ; \
  }

int HSM2temp(
     GENmodel *inModel,
     CKTcircuit *ckt)
{
  HSM2model *model = (HSM2model *)inModel ;
  HSM2instance *here ;
  HSM2binningParam *pParam ;
  HSM2modelMKSParam *modelMKS ;
  HSM2hereMKSParam  *hereMKS ;
  double mueph ;
  double Leff, dL , LG, Weff, dW , WG , WL , Lgate , Wgate;
  double Nsubpp, Nsubps, Nsub, q_Nsub, Nsubb, Npext ;
  double Lod_half, Lod_half_ref ;
  double MUEPWD = 0.0 ;
  double MUEPLD = 0.0 ;
  double GDLD = 0.0 ;
  double T1, T2, T3 ;
  /* temperature-dependent variables */
  double Eg ,TTEMP, beta, Nin;
  double js, jssw, js2, jssw2 ;
  int i;

  /* declarations for the sc3 clamping part */
  double A, beta_inv, c_eox, cnst0, cnst1, Cox, Cox_inv;
  double Denom, dPpg, dVth, dVthLP, dVthLP_dVb, dVthSC, dVthW;
  double dVth0, dVth0_dVb, fac1, limVgp_dVbs, Pb20, Ps0, Ps0_dVbs;
  double Ps0_min, Qb0, sc3lim, sc3Vbs, sc3Vgs, term1, term2, term3, term4;
  double Tox, T0, T3_dVb, T4, T5, T6, T6_dVb, T8, T8_dVb;
  double T9, T9_dVb, Vgp, Vgs_min, Vfb, Vthp, Vth0;


  for ( ;model ;model = model->HSM2nextModel ) {

    modelMKS = &model->modelMKS ;

    for ( here = model->HSM2instances; here; here = here->HSM2nextInstance ) {
      pParam = &here->pParam ;
      hereMKS = &here->hereMKS ;

      here->HSM2_lgate = Lgate = here->HSM2_l + model->HSM2_xl ;
      Wgate = here->HSM2_w / here->HSM2_nf  + model->HSM2_xw ;

      LG = Lgate * C_m2um ;
      here->HSM2_wg = WG = Wgate * C_m2um ; 
      WL = WG * LG ;
      MUEPWD = model->HSM2_muepwd * C_m2um ;
      MUEPLD = model->HSM2_muepld * C_m2um ;

      /* Band gap */
      here->HSM2_egtnom = pParam->HSM2_eg0 - model->HSM2_ktnom    
        * ( 90.25e-6 + model->HSM2_ktnom * 1.0e-7 ) ;
  
      /* C_EOX */
      here->HSM2_cecox = C_VAC * model->HSM2_kappa ;
  
      /* Vth reduction for small Vds */
      here->HSM2_msc = model->HSM2_scp22  ;

      /* Poly-Si Gate Depletion */
      if ( pParam->HSM2_pgd1 == 0.0 ) {
        here->HSM2_flg_pgd = 0 ;
      } else {
        here->HSM2_flg_pgd = 1 ;
      }


      /* CLM5 & CLM6 */
      here->HSM2_clmmod = 1e0 + pow( LG , model->HSM2_clm5 ) * model->HSM2_clm6 ;

      /* Half length of diffusion */
      T1 = 1.0 / (model->HSM2_saref + 0.5 * here->HSM2_l)
         + 1.0 / (model->HSM2_sbref + 0.5 * here->HSM2_l);
      Lod_half_ref = 2.0 / T1 ;

      if (here->HSM2_sa > 0.0 && here->HSM2_sb > 0.0 &&
	  (here->HSM2_nf == 1.0 ||
           (here->HSM2_nf > 1.0 && here->HSM2_sd > 0.0))) {
        T1 = 0.0;
        for (i = 0; i < here->HSM2_nf; i++) {
          T1 += 1.0 / (here->HSM2_sa + 0.5 * here->HSM2_l
                       + i * (here->HSM2_sd + here->HSM2_l))
              + 1.0 / (here->HSM2_sb + 0.5 * here->HSM2_l
                       + i * (here->HSM2_sd + here->HSM2_l));
        }
        Lod_half = 2.0 * here->HSM2_nf / T1;
      } else {
        Lod_half = 0.0;
      }

      Npext = modelMKS->HSM2_npext * ( 1.0 + model->HSM2_npextw / pow( WG, model->HSM2_npextwp ) ); /* new */
      here->HSM2_mueph1 = pParam->HSM2_mueph1 ;
      here->HSM2_nsubp  = pParam->HSM2_nsubp ;
      here->HSM2_nsubc  = pParam->HSM2_nsubc ;

      /* DFM */
      if ( model->HSM2_codfm == 1 && here->HSM2_nsubcdfm_Given ) {
	RANGECHECK(here->HSM2_nsubcdfm,   1.0e16,   1.0e19, "NSUBCDFM") ;
 	here->HSM2_mueph1 *= here->HSM2_mphdfm
	  * ( log(hereMKS->HSM2_nsubcdfm) - log(here->HSM2_nsubc) ) + 1.0 ;
	here->HSM2_nsubp += hereMKS->HSM2_nsubcdfm - here->HSM2_nsubc ;
 	Npext += hereMKS->HSM2_nsubcdfm - here->HSM2_nsubc ;
 	here->HSM2_nsubc = hereMKS->HSM2_nsubcdfm ;
      }

	/* WPE */
        T0 = modelMKS->HSM2_nsubcwpe *
              ( here->HSM2_sca
                + model->HSM2_web * here->HSM2_scb
                + model->HSM2_wec * here->HSM2_scc ) ;
        here->HSM2_nsubc +=  T0 ;
        Fn_SLtemp( here->HSM2_nsubc , here->HSM2_nsubc , Nsubmin , Nsubmin_dlt ) ;
        T0 = modelMKS->HSM2_nsubpwpe *
              ( here->HSM2_sca
                + model->HSM2_web * here->HSM2_scb
                + model->HSM2_wec * here->HSM2_scc ) ;
        here->HSM2_nsubp +=  T0 ;
        Fn_SLtemp( here->HSM2_nsubp , here->HSM2_nsubp , Nsubmin , Nsubmin_dlt ) ;
        T0 = modelMKS->HSM2_npextwpe *
              ( here->HSM2_sca
                + model->HSM2_web * here->HSM2_scb
                + model->HSM2_wec * here->HSM2_scc ) ;
        Npext +=  T0 ;
        Fn_SLtemp( Npext , Npext , Nsubmin , Nsubmin_dlt ) ;
	/* WPE end */

      /* Coulomb Scattering */
      here->HSM2_muecb0 = pParam->HSM2_muecb0 * pow( LG, model->HSM2_muecb0lp );
      here->HSM2_muecb1 = pParam->HSM2_muecb1 * pow( LG, model->HSM2_muecb1lp );

      /* Phonon Scattering (temperature-independent part) */
      mueph = here->HSM2_mueph1  
        * (1.0e0 + (model->HSM2_muephw / pow( WG + MUEPWD , model->HSM2_muepwp))) 
        * (1.0e0 + (model->HSM2_muephl / pow( LG + MUEPLD , model->HSM2_mueplp))) 
        * (1.0e0 + (model->HSM2_muephw2 / pow( WG, model->HSM2_muepwp2))) 
        * (1.0e0 + (model->HSM2_muephl2 / pow( LG, model->HSM2_mueplp2))) 
        * (1.0e0 + (model->HSM2_muephs / pow( WL, model->HSM2_muepsp)));  
      if (Lod_half > 0.0) {
        T1 = 1.0e0 / (1.0e0 + pParam->HSM2_muesti2) ;
        T2 = pow (pParam->HSM2_muesti1 / Lod_half, pParam->HSM2_muesti3) ;
        T3 = pow (pParam->HSM2_muesti1 / Lod_half_ref, pParam->HSM2_muesti3) ;
        here->HSM2_mueph = mueph * (1.0e0 + T1 * T2) / (1.0e0 + T1 * T3); 
      } else {
        here->HSM2_mueph = mueph;
      }
      
      /* Surface Roughness Scattering */
      here->HSM2_muesr = model->HSM2_muesr0 
        * (1.0e0 + (model->HSM2_muesrl / pow (LG, model->HSM2_mueslp))) 
        * (1.0e0 + (model->HSM2_muesrw / pow (WG, model->HSM2_mueswp))) ;

      /* Coefficients of Qbm for Eeff */
      T1 = pow( LG, model->HSM2_ndeplp ) ;
      T2 = pow( WG, model->HSM2_ndepwp ) ; /* new */
      T3 = T1 + model->HSM2_ndepl ;
      T4 = T2 + model->HSM2_ndepw ;
      if( T3 < 1e-8 ) { T3 = 1e-8; } 
      if( T4 < 1e-8 ) { T4 = 1e-8; } 
      here->HSM2_ndep_o_esi = ( pParam->HSM2_ndep * T1 ) / T3  * T2 / T4
	/ C_ESI ;
      here->HSM2_ninv_o_esi = pParam->HSM2_ninv / C_ESI ;

      /* Metallurgical channel geometry */
      dL = model->HSM2_xld 
        + (modelMKS->HSM2_ll / pow (Lgate + model->HSM2_lld, model->HSM2_lln)) ;
      dW = model->HSM2_xwd 
        + (modelMKS->HSM2_wl / pow (Wgate + model->HSM2_wld, model->HSM2_wln)) ;  
    
      Leff = Lgate - 2.0e0 * dL ;
      if ( Leff <= 1.0e-9 ) {   
        SPfrontEnd->IFerrorf
          ( 
           ERR_FATAL, 
           "HiSIM2: MOSFET(%s) MODEL(%s): effective channel length is smaller than 1nm", 
           model->HSM2modName, here->HSM2name);
        return (E_BADPARM);
      }
      here->HSM2_leff = Leff ;

      /* Wg dependence for short channel devices */
      here->HSM2_lgatesm = Lgate + model->HSM2_wl1 / pow( WL , model->HSM2_wl1p ) ;
      here->HSM2_dVthsm = pParam->HSM2_wl2 / pow( WL , model->HSM2_wl2p ) ;

      /* Lg dependence of wsti */
      T1 = 1.0e0 + model->HSM2_wstil / pow( here->HSM2_lgatesm * C_m2um , model->HSM2_wstilp ) ;
      T2 = 1.0e0 + model->HSM2_wstiw / pow( WG , model->HSM2_wstiwp ) ;
      here->HSM2_wsti = pParam->HSM2_wsti * T1 * T2 ;

      here->HSM2_weff = Weff = Wgate - 2.0e0 * dW ;
      if ( Weff <= 0.0 ) {   
        SPfrontEnd->IFerrorf
          ( 
           ERR_FATAL, 
           "HiSIM2: MOSFET(%s) MODEL(%s): effective channel width is negative or 0", 
           model->HSM2modName, here->HSM2name);
        return (E_BADPARM);
      }
      here->HSM2_weff_nf = Weff * here->HSM2_nf ;

      /* Surface impurity profile */
      /* Nsubp */
      if(model->HSM2_nsubpfac < 1.0) {
      T1 = 2.0 * ( 1.0 - model->HSM2_nsubpfac ) / model->HSM2_nsubpl * LG + 2.0 * model->HSM2_nsubpfac - 1.0 ;
      Fn_SUtemp( T1 , T1 , 1 , model->HSM2_nsubpdlt ) ;
      Fn_SLtemp( T1 , T1 , model->HSM2_nsubpfac  , model->HSM2_nsubpdlt ) ;
      here->HSM2_nsubp *= T1 ;
      }

      /* Note: Sign Changed --> */
      Nsubpp = here->HSM2_nsubp  
        * (1.0e0 + (model->HSM2_nsubpw / pow (WG, model->HSM2_nsubpwp))) ;
      /* <-- Note: Sign Changed */

      if (Lod_half > 0.0) {
        T1 = 1.0e0 / (1.0e0 + pParam->HSM2_nsubpsti2) ;
        T2 = pow (pParam->HSM2_nsubpsti1 / Lod_half, pParam->HSM2_nsubpsti3) ;
        T3 = pow (pParam->HSM2_nsubpsti1 / Lod_half_ref, pParam->HSM2_nsubpsti3) ;
        Nsubps = Nsubpp * (1.0e0 + T1 * T2) / (1.0e0 + T1 * T3) ;
      } else {
        Nsubps = Nsubpp ;
      }

      T2 = 1.0e0 + ( model->HSM2_nsubcw / pow ( WG, model->HSM2_nsubcwp )) ;
      T2 *= 1.0e0 + ( model->HSM2_nsubcw2 / pow ( WG, model->HSM2_nsubcwp2 )) ;
      T3 = modelMKS->HSM2_nsubcmax / here->HSM2_nsubc ;

      Fn_SUtemp( T1 , T2 , T3 , 0.01 ) ;
      here->HSM2_nsubc *= T1 ;

      if (Lod_half > 0.0) {
        T1 = 1.0e0 / (1.0e0 + pParam->HSM2_nsubcsti2) ;
        T2 = pow (pParam->HSM2_nsubcsti1 / Lod_half, pParam->HSM2_nsubcsti3) ;
        T3 = pow (pParam->HSM2_nsubcsti1 / Lod_half_ref, pParam->HSM2_nsubcsti3) ;
        here->HSM2_nsubc = here->HSM2_nsubc * (1.0e0 + T1 * T2) / (1.0e0 + T1 * T3) ;
      }

      if(model->HSM2_coerrrep && (here->HSM2_nsubc <= 0.0)) {
        fprintf ( stderr , "*** warning(HiSIM): actual NSUBC value is negative -> reset to 1E+15.\n" ) ;
        fprintf ( stderr , "    The model parameter  NSUBCW/NSUBCWP and/or NSUBCW2/NSUBCW2P might be wrong.\n" ) ;
        here->HSM2_nsubc = 1e15 / C_cm2m_p3 ;
      }
      if(model->HSM2_coerrrep && (Npext < here->HSM2_nsubc || Npext > here->HSM2_nsubp)) {
        fprintf ( stderr , "*** warning(HiSIM): actual NPEXT value is smaller than NSUBC and/or greater than NSUBP.\n" ) ;
        fprintf ( stderr , "    ( Npext = %e , NSUBC = %e , NSUBP = %e ) \n",Npext,here->HSM2_nsubc,here->HSM2_nsubp);
        fprintf ( stderr , "    The model parameter  NPEXTW and/or NPEXTWP might be wrong.\n" ) ;
      }

      if( Lgate > model->HSM2_lp ){
        Nsub = (here->HSM2_nsubc * (Lgate - model->HSM2_lp) 
                +  Nsubps  * model->HSM2_lp) / Lgate ;
      } else {
        Nsub = Nsubps
          + (Nsubps - here->HSM2_nsubc) * (model->HSM2_lp - Lgate) 
          / model->HSM2_lp ;
      }
      T3 = 0.5e0 * Lgate - model->HSM2_lp ;
      Fn_SZtemp( T3 , T3 , lpext_dlt ) ;
      T1 = Fn_Max(0.0e0, model->HSM2_lpext ) ;
      T2 = T3 * T1 / ( T3 + T1 ) ;

      here->HSM2_nsub = 
        Nsub = Nsub + T2 * (Npext - here->HSM2_nsubc) / Lgate ;
      here->HSM2_qnsub = q_Nsub  = C_QE * Nsub ;
      here->HSM2_qnsub_esi = q_Nsub * C_ESI ;
      here->HSM2_2qnsub_esi = 2.0 * here->HSM2_qnsub_esi ;

      /* Pocket Overlap (temperature-independent part) */
      if ( Lgate <= 2.0e0 * model->HSM2_lp ) {
        Nsubb = 2.0e0 * Nsubps 
          - (Nsubps - here->HSM2_nsubc) * Lgate 
          / model->HSM2_lp - here->HSM2_nsubc ;
        here->HSM2_ptovr0 = log (Nsubb / here->HSM2_nsubc) ;
      } else {
        here->HSM2_ptovr0 = 0.0e0 ;
      }

      /* costi0 and costi1 for STI transistor model (temperature-independent part) */
      here->HSM2_costi00 = sqrt (2.0 * C_QE * pParam->HSM2_nsti * C_ESI ) ;
      here->HSM2_nsti_p2 = 1.0 / ( pParam->HSM2_nsti * pParam->HSM2_nsti ) ;

      /* Velocity Temperature Dependence (Temperature-dependent part will be multiplied later.) */
      here->HSM2_vmax0 = (1.0e0 + (model->HSM2_vover / pow (LG, model->HSM2_voverp)))
        * (1.0e0 + (model->HSM2_vovers / pow (WL, model->HSM2_voversp))) ;

      /* 2 phi_B (temperature-independent) */
      /* @300K, with pocket */
      here->HSM2_pb20 = 2.0e0 / C_b300 * log (Nsub / C_Nin0) ;
      /* @300K, w/o pocket */
      here->HSM2_pb2c = 2.0e0 / C_b300 * log (here->HSM2_nsubc / C_Nin0) ;


      /* constant for Poly depletion */
      here->HSM2_cnstpgd = pow ( 1e0 + 1e0 / LG , model->HSM2_pgd4 ) 
        * pParam->HSM2_pgd1 ;




      /* Gate resistance */
      if ( here->HSM2_corg == 1 ) {
        T1 = here->HSM2_xgw + Weff / (3.0e0 * here->HSM2_ngcon);
        T2 = Lgate - here->HSM2_xgl;
        here->HSM2_grg = model->HSM2_rshg * T1 / (here->HSM2_ngcon * T2 * here->HSM2_nf);
        if (here->HSM2_grg > 1.0e-3) here->HSM2_grg = here->HSM2_m / here->HSM2_grg;
        else {
          here->HSM2_grg = here->HSM2_m * 1.0e3;
          if(model->HSM2_coerrrep) 
          printf("warning(HiSIM2): The gate conductance reset to 1.0e3 mho.\n");
        }
      }

      /* Process source/drain series resistamce */
      here->HSM2_rd = 0.0;
      if ( model->HSM2_rsh > 0.0 ) {
        here->HSM2_rd += model->HSM2_rsh * here->HSM2_nrd ;
      } 
      if ( model->HSM2_rd > 0.0 ) {
       here->HSM2_rd += model->HSM2_rd / here->HSM2_weff_nf ;
     }

      here->HSM2_rs = 0.0;
      if ( model->HSM2_rsh > 0.0 ) {
        here->HSM2_rs += model->HSM2_rsh * here->HSM2_nrs ;
      }
      if ( model->HSM2_rs > 0.0 ) {
        here->HSM2_rs += model->HSM2_rs / here->HSM2_weff_nf ; 
      }

      if (model->HSM2_corsrd < 0) {
        if ( here->HSM2_rd > 0.0 ) {
          here->HSM2drainConductance = here->HSM2_m / here->HSM2_rd ;
        } else {
          here->HSM2drainConductance = 0.0;
        }
        if ( here->HSM2_rs > 0.0 ) {
          here->HSM2sourceConductance = here->HSM2_m / here->HSM2_rs ;
        } else {
          here->HSM2sourceConductance = 0.0;
        }
      } else if (model->HSM2_corsrd > 0) {
        here->HSM2drainConductance = 0.0 ;
        here->HSM2sourceConductance = 0.0 ;
        if ( here->HSM2_rd > 0.0 && model->HSM2_cothrml != 0 ) {
          here->HSM2internalGd = here->HSM2_m / here->HSM2_rd ;
        } else {
          here->HSM2internalGd = 0.0;
        }
        if ( here->HSM2_rs > 0.0 && model->HSM2_cothrml != 0 ) {
          here->HSM2internalGs = here->HSM2_m / here->HSM2_rs ;
        } else {
          here->HSM2internalGs = 0.0;
        }
      } else {
        here->HSM2drainConductance = 0.0 ;
        here->HSM2sourceConductance = 0.0 ;
      }


      /* Body resistance */
      if ( here->HSM2_corbnet == 1 ) {
        if (here->HSM2_rbdb < 1.0e-3) here->HSM2_grbdb = here->HSM2_m * 1.0e3 ; /* in mho */
        else here->HSM2_grbdb = here->HSM2_m * ( model->HSM2_gbmin + 1.0 / here->HSM2_rbdb ) ;

        if (here->HSM2_rbpb < 1.0e-3) here->HSM2_grbpb = here->HSM2_m * 1.0e3 ;
        else here->HSM2_grbpb = here->HSM2_m * ( model->HSM2_gbmin + 1.0 / here->HSM2_rbpb ) ;

        if (here->HSM2_rbps < 1.0e-3) here->HSM2_grbps = here->HSM2_m * 1.0e3 ;
        else here->HSM2_grbps = here->HSM2_m * ( model->HSM2_gbmin + 1.0 / here->HSM2_rbps ) ;

        if (here->HSM2_rbsb < 1.0e-3) here->HSM2_grbsb = here->HSM2_m * 1.0e3 ;
        else here->HSM2_grbsb = here->HSM2_m * ( model->HSM2_gbmin + 1.0 / here->HSM2_rbsb ) ;

        if (here->HSM2_rbpd < 1.0e-3) here->HSM2_grbpd = here->HSM2_m * 1.0e3 ;
        else here->HSM2_grbpd = here->HSM2_m * ( model->HSM2_gbmin + 1.0 / here->HSM2_rbpd ) ;
      }

      /* Vdseff */
      T1 = model->HSM2_ddltslp * LG + model->HSM2_ddltict ;
      here->HSM2_ddlt = T1 * model->HSM2_ddltmax / ( T1 + model->HSM2_ddltmax ) + 1.0 ;

      /* Isub */
      T2 = pow( Weff , model->HSM2_svgswp ) ;
      here->HSM2_vg2const = pParam->HSM2_svgs
         * ( 1.0e0
           + modelMKS->HSM2_svgsl / pow( here->HSM2_lgate , model->HSM2_svgslp ) )
         * ( T2 / ( T2 + modelMKS->HSM2_svgsw ) ) ; 

      here->HSM2_xvbs = pParam->HSM2_svbs 
         * ( 1.0e0
           + modelMKS->HSM2_svbsl / pow( here->HSM2_lgate , model->HSM2_svbslp ) ) ;
      here->HSM2_xgate = modelMKS->HSM2_slg  
         * ( 1.0
         + modelMKS->HSM2_slgl / pow( here->HSM2_lgate , model->HSM2_slglp ) ) ;

      here->HSM2_xsub1 = pParam->HSM2_sub1 
         * ( 1.0 
         + modelMKS->HSM2_sub1l / pow( here->HSM2_lgate , model->HSM2_sub1lp ) ) ;

      here->HSM2_xsub2 = pParam->HSM2_sub2
         * ( 1.0 + modelMKS->HSM2_sub2l / here->HSM2_lgate ) ;

      /* Fringing capacitance */
      here->HSM2_cfrng = C_EOX / ( C_Pi / 2.0e0 ) * here->HSM2_weff_nf 
         * log( 1.0e0 + model->HSM2_tpoly / model->HSM2_tox ) ; 

      /* Additional term of lateral-field-induced capacitance */
      here->HSM2_cqyb0 = C_m2um * here->HSM2_weff_nf
	* model->HSM2_xqy1 / pow( LG , model->HSM2_xqy2 ) ;

      /* Parasitic component of the channel current */
      GDLD = model->HSM2_gdld * C_m2um ;
      here->HSM2_ptl0 = model->HSM2_ptl * pow( LG        , - model->HSM2_ptlp ) ;
      here->HSM2_pt40 = model->HSM2_pt4 * pow( LG        , - model->HSM2_pt4p ) ;
      here->HSM2_gdl0 = model->HSM2_gdl * pow( LG + GDLD , - model->HSM2_gdlp ) ;


      /*-----------------------------------------------------------*
       * Temperature dependent constants. 
       *-----------------*/
        TTEMP = ckt->CKTtemp ;
        if ( here->HSM2_temp_Given ) TTEMP = here->HSM2_ktemp ;
        if ( here->HSM2_dtemp_Given ) {
          TTEMP = TTEMP + here->HSM2_dtemp ;
          here->HSM2_ktemp = TTEMP ;
        }

        /* Band gap */
        T1 = TTEMP - model->HSM2_ktnom ;
        T2 = TTEMP * TTEMP - model->HSM2_ktnom * model->HSM2_ktnom ;
        here->HSM2_eg = Eg = here->HSM2_egtnom - pParam->HSM2_bgtmp1 * T1
          - pParam->HSM2_bgtmp2 * T2 ;
        here->HSM2_sqrt_eg = sqrt( Eg ) ;


	T1 = 1.0 / TTEMP ;
	T2 = 1.0 / model->HSM2_ktnom ;
	T3 = here->HSM2_egtnom + model->HSM2_egig
	  + model->HSM2_igtemp2 * ( T1 - T2 )
	  + model->HSM2_igtemp3 * ( T1 * T1 - T2 * T2 ) ;
 	here->HSM2_egp12 = sqrt ( T3 ) ;
 	here->HSM2_egp32 = T3 * here->HSM2_egp12 ;

        
        /* Inverse of the thermal voltage */
        here->HSM2_beta = beta = C_QE / (C_KB * TTEMP) ;
        here->HSM2_beta_inv = 1.0 / beta ;
        here->HSM2_beta2 = beta * beta ;
        here->HSM2_betatnom = C_QE / (C_KB * model->HSM2_ktnom) ;

        /* Intrinsic carrier concentration */
        here->HSM2_nin = Nin = C_Nin0 * pow (TTEMP / model->HSM2_ktnom, 1.5e0) 
          * exp (- Eg / 2.0e0 * beta + here->HSM2_egtnom / 2.0e0 * here->HSM2_betatnom) ;


        /* Phonon Scattering (temperature-dependent part) */
        T1 =  pow (TTEMP / model->HSM2_ktnom, pParam->HSM2_muetmp) ;
        here->HSM2_mphn0 = T1 / here->HSM2_mueph ;
        here->HSM2_mphn1 = here->HSM2_mphn0 * model->HSM2_mueph0 ;


        /* Pocket Overlap (temperature-dependent part) */
        here->HSM2_ptovr = here->HSM2_ptovr0 / beta ;


        /* Velocity Temperature Dependence */
        T1 =  TTEMP / model->HSM2_ktnom ;
        here->HSM2_vmax = here->HSM2_vmax0 * pParam->HSM2_vmax 
          / (1.8 + 0.4 * T1 + 0.1 * T1 * T1 - pParam->HSM2_vtmp * (1.0e0 - T1)) ;


        /* Coefficient of the F function for bulk charge */
        here->HSM2_cnst0 = sqrt ( 2.0 * C_ESI * C_QE * here->HSM2_nsub / beta ) ;
	here->HSM2_cnst0over = here->HSM2_cnst0 * sqrt( pParam->HSM2_nover / here->HSM2_nsub ) ;     

        /* 2 phi_B (temperature-dependent) */
        /* @temp, with pocket */
        here->HSM2_pb2 =  2.0e0 / beta * log (here->HSM2_nsub / Nin) ;
        if ( pParam->HSM2_nover != 0.0) {
	   here->HSM2_pb2over = 2.0 / beta * log( pParam->HSM2_nover / Nin ) ;

           /* (1 / cnst1 / cnstCoxi) for Ps0LD_iniB */
           T1 = here->HSM2_cnst0over * model->HSM2_tox / here->HSM2_cecox  ;
           T2 = pParam->HSM2_nover / Nin ;
           T1 = T2 * T2 / ( T1 * T1 ) ;
           here->HSM2_ps0ldinib  = T1 ; /* (1 / cnst1 / cnstCoxi) */

        }else {
           here->HSM2_pb2over = 0.0 ;
           here->HSM2_ps0ldinib  = 0.0 ;
        }


        /* Depletion Width */
        T1 = 2.0e0 * C_ESI / C_QE ;
        here->HSM2_wdpl = sqrt ( T1 / here->HSM2_nsub ) ;
        here->HSM2_wdplp = sqrt( T1 / ( here->HSM2_nsubp ) ) ; 

        /* cnst1: n_{p0} / p_{p0} */
        T1 = Nin / here->HSM2_nsub ;
        here->HSM2_cnst1 = T1 * T1 ;


        /* for substrate-source/drain junction diode. */
        js   = pParam->HSM2_js0
          * exp ((here->HSM2_egtnom * here->HSM2_betatnom - Eg * beta
                  + model->HSM2_xti * log (TTEMP / model->HSM2_ktnom)) / pParam->HSM2_nj) ;
        jssw = pParam->HSM2_js0sw
          * exp ((here->HSM2_egtnom * here->HSM2_betatnom - Eg * beta 
                  + model->HSM2_xti * log (TTEMP / model->HSM2_ktnom)) / model->HSM2_njsw) ;

        js2  = pParam->HSM2_js0
          * exp ((here->HSM2_egtnom * here->HSM2_betatnom - Eg * beta
                  + model->HSM2_xti2 * log (TTEMP / model->HSM2_ktnom)) / pParam->HSM2_nj) ;  
        jssw2 = pParam->HSM2_js0sw
          * exp ((here->HSM2_egtnom * here->HSM2_betatnom - Eg * beta
                  + model->HSM2_xti2 * log (TTEMP / model->HSM2_ktnom)) / model->HSM2_njsw) ; 
      
        here->HSM2_isbd = here->HSM2_ad * js + here->HSM2_pd * jssw ;
        here->HSM2_isbd2 = here->HSM2_ad * js2 + here->HSM2_pd * jssw2 ;
        here->HSM2_isbs = here->HSM2_as * js + here->HSM2_ps * jssw ;
        here->HSM2_isbs2 = here->HSM2_as * js2 + here->HSM2_ps * jssw2 ;

        here->HSM2_vbdt = pParam->HSM2_nj / beta 
          * log (pParam->HSM2_vdiffj * (TTEMP / model->HSM2_ktnom) * (TTEMP / model->HSM2_ktnom) 
                 / (here->HSM2_isbd + 1.0e-50) + 1) ;
        here->HSM2_vbst = pParam->HSM2_nj / beta 
          * log (pParam->HSM2_vdiffj * (TTEMP / model->HSM2_ktnom) * (TTEMP / model->HSM2_ktnom) 
                 / (here->HSM2_isbs + 1.0e-50) + 1) ;

        here->HSM2_exptemp = exp (((TTEMP / model->HSM2_ktnom) - 1) * model->HSM2_ctemp) ;
        here->HSM2_jd_nvtm_inv = 1.0 / ( pParam->HSM2_nj / beta ) ;
        here->HSM2_jd_expcd = exp (here->HSM2_vbdt * here->HSM2_jd_nvtm_inv ) ;
        here->HSM2_jd_expcs = exp (here->HSM2_vbst * here->HSM2_jd_nvtm_inv ) ;

	/* costi0 and costi1 for STI transistor model (temperature-dependent part) */
	here->HSM2_costi0 = here->HSM2_costi00 * sqrt(here->HSM2_beta_inv) ;
	here->HSM2_costi0_p2 = here->HSM2_costi0 * here->HSM2_costi0 ; 
	here->HSM2_costi1 = here->HSM2_nin * here->HSM2_nin * here->HSM2_nsti_p2 ;

        /* check if SC3 is too large */
        if (pParam->HSM2_sc3 && model->HSM2_sc3Vbs < 0.0) {

          beta     = here->HSM2_beta ;
          beta_inv = here->HSM2_beta_inv ;
          Weff     = here->HSM2_weff ;
          Vfb      = pParam->HSM2_vfbc ;
          Pb20     = here->HSM2_pb20 ; 
          cnst0    = here->HSM2_cnst0 ;
          cnst1    = here->HSM2_cnst1 ;
          c_eox    = here->HSM2_cecox ;
          Tox      = model->HSM2_tox ;
          Cox      = c_eox / Tox ;
          Cox_inv  = 1.0 / Cox ;
          fac1     = cnst0 * Cox_inv ;
          Vgs_min  = model->HSM2_type * model->HSM2_Vgsmin ;
          sc3Vbs   = model->HSM2_sc3Vbs ;
          sc3Vgs   = 2.0 ;
          Ps0_min  = 2.0 * beta_inv * log(-Vgs_min/fac1) ;

          /* approximate solution of Poisson equation for large
             Vgs and negative Vbs (3 iterations!)*/
          Vgp = sc3Vgs - Vfb;
          Denom = fac1*sqrt(cnst1);
          Ps0 = 2.0 * beta_inv * log(Vgp/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);
          Ps0_dVbs = 0.0;

          T1   = here->HSM2_2qnsub_esi ;
          Qb0  = sqrt ( T1 ) ;
          Vthp = Ps0 + Vfb + Qb0 * Cox_inv + here->HSM2_ptovr ;

          T1     = 2.0 * C_QE * here->HSM2_nsubc * C_ESI ;
          T2     = sqrt( T1 ) ;
          Vth0   = Ps0 + Vfb + T2 * Cox_inv ;
          T1     = C_ESI * Cox_inv ;
          T2     = here->HSM2_wdplp ;
          T4     = 1.0e0 / ( model->HSM2_lp * model->HSM2_lp ) ;
          T3     = 2.0 * ( model->HSM2_vbi - Pb20 ) * T2 * T4 ;
          T5     = T1 * T3 ;
          T6     = Ps0 - sc3Vbs ;
          T6_dVb = Ps0_dVbs - 1.0 ;
          dVth0  = T5 * sqrt( T6 ) ;
          dVth0_dVb = T5 * 0.5 / sqrt( T6 ) * T6_dVb;
          T1     = Vthp - Vth0 ;
          T9     = Ps0 - sc3Vbs ;
          T9_dVb = Ps0_dVbs - 1.0 ;
          T3     = pParam->HSM2_scp1 + pParam->HSM2_scp3 * T9 / model->HSM2_lp;
          T3_dVb = pParam->HSM2_scp3 * T9_dVb / model->HSM2_lp ;
          dVthLP = T1 * dVth0 * T3 ;
          dVthLP_dVb = T1 * dVth0_dVb * T3 + T1 * dVth0 * T3_dVb;

          T3 = here->HSM2_lgate - model->HSM2_parl2 ;
          T4 = 1.0e0 / ( T3 * T3 ) ;
          T0 = C_ESI * here->HSM2_wdpl * 2.0e0 * ( model->HSM2_vbi - Pb20 ) * T4 ;
          T2 = T0 * Cox_inv ;
          T5 = pParam->HSM2_sc3 / here->HSM2_lgate ;
          T6 = pParam->HSM2_sc1 + T5 * ( Ps0 - sc3Vbs ) ;
          T1 = T6 ;
          A  = T2 * T1 ;
          T9 = Ps0 - sc3Vbs + Ps0_min ;
          T9_dVb = Ps0_dVbs - 1.0 ;
          T8     = sqrt( T9 ) ;
          T8_dVb = 0.5 * T9_dVb / T8 ;
          dVthSC = A * T8 ;

          T1 = 1.0 / Cox ;
          T3 = 1.0 / ( Cox + pParam->HSM2_wfc / Weff ) ;
          T5 = T1 - T3 ;
          dVthW = Qb0 * T5 + pParam->HSM2_wvth0 / here->HSM2_wg ;

          dVth = dVthSC + dVthLP + dVthW + here->HSM2_dVthsm ;
          dPpg = 0.0 ;
          Vgp = sc3Vgs - Vfb + dVth - dPpg ;

          /* Recalculation of Ps0, using more accurate Vgp */
          Ps0 = 2.0 * beta_inv * log(Vgp/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);
          Ps0 = 2.0 * beta_inv * log((Vgp-Ps0)/Denom);

          term1 = Vgp - Ps0;
          term2 = sqrt(beta*(Ps0-sc3Vbs)-1.0);
          term3 = term1 + fac1 * term2;
          term4 = cnst1 * exp(beta*Ps0);
          limVgp_dVbs = - beta * (term3 + 0.5*fac1 * term4/term2)
                  / (2.0*term1/fac1/fac1*term3 - term4);
          T2     = T0 * Cox_inv ;
          sc3lim = here->HSM2_lgate / T2 
               * (limVgp_dVbs - dVthLP_dVb - T2*pParam->HSM2_sc1*T8_dVb)
               / ((Ps0-sc3Vbs)*T8_dVb +(Ps0_dVbs-1.0)*T8);
          if (sc3lim < 1.0e-20)
              sc3lim = 1e-20 ;
          if (sc3lim < pParam->HSM2_sc3 * 0.999) {
            pParam->HSM2_sc3 = sc3lim;
          }
        }
    } /* End of instance loop */
  }
  return(OK);
}
