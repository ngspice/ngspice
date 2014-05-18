/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 1  SUBVERSION : 2  REVISION : 4 )
 Model Parameter VERSION : 1.23
 FILE : hsmhvld.c

 DATE : 2013.04.30

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "hisimhv.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

#define SHOW_EPS_QUANT 1.0e-15

static void ShowPhysVals
(
 HSMHVinstance *here,
 HSMHVmodel *model,
 int isFirst,
 double vds,
 double vgs,
 double vbs,
 double vgd,
 double vbd,
 double vgb
 )
{

  NG_IGNORE(vgd);
  NG_IGNORE(vbd);

  /*                                                            */
  /*                                                            */
  /* Note: This function is not yet adapted to the flat version */
  /*                                                            */
  /*                                                            */

  /* regard the epsilon-quantity as 0.0 */
  vds = (fabs(vds) < SHOW_EPS_QUANT) ? 0.0 : vds;
  vgs = (fabs(vgs) < SHOW_EPS_QUANT) ? 0.0 : vgs;
  vbs = (fabs(vbs) < SHOW_EPS_QUANT) ? 0.0 : vbs;
  vgb = (fabs(vgb) < SHOW_EPS_QUANT) ? 0.0 : vgb;
  switch (model->HSMHV_show) {
  case 1:
    if (isFirst) printf("Vds        Ids\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_mode*here->HSMHV_ids);
    break;
  case 2:
    if (isFirst) printf("Vgs        Ids\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_mode*here->HSMHV_ids);
    break;
  case 3:
    if (isFirst) printf("Vgs        log10(|Ids|)\n");
    printf("%e %e\n", model->HSMHV_type*vgs, log10(here->HSMHV_ids));
    break;
  case 4:
    if (isFirst) printf("log10(|Ids|)    gm/|Ids|\n");
    if (here->HSMHV_ids == 0.0)
      printf("I can't show gm/Ids - log10(Ids), because Ids = 0.\n");
    else
      printf("%e %e\n",  log10(here->HSMHV_ids), here->HSMHV_gm/here->HSMHV_ids);
    break;
  case 5:
    if (isFirst) printf("Vds        gds\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_gds);
    break;
  case 6:
    if (isFirst) printf("Vgs        gm\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_gm);
    break;
  case 7:
    if (isFirst) printf("Vbs        gbs\n");
    printf("%e %e\n", model->HSMHV_type*vbs, here->HSMHV_gmbs);
    break;
  case 8:
    if (isFirst) printf("Vgs        Cgg\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_cggb);
    break;
  case 9:
    if (isFirst) printf("Vgs        Cgs\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_cgsb);
    break;
  case 10:
    if (isFirst) printf("Vgs        Cgd\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_cgdb);
    break;
  case 11:
    if (isFirst) printf("Vgs        Cgb\n");
    printf("%e %e\n", model->HSMHV_type*vgs, -(here->HSMHV_cggb+here->HSMHV_cgsb+here->HSMHV_cgdb));
    break;
  case 12:
    if (isFirst) printf("Vds        Csg\n");
    printf("%e %e\n", model->HSMHV_type*vds, -(here->HSMHV_cggb+here->HSMHV_cbgb+here->HSMHV_cdgb));
    break;
  case 13:
    if (isFirst) printf("Vds        Cdg\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_cdgb);
    break;
  case 14:
    if (isFirst) printf("Vds        Cbg\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_cbgb);
    break;
  case 15:
    if (isFirst) printf("Vds        Cgg\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_cggb);
    break;
  case 16:
    if (isFirst) printf("Vds        Cgs\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_cgsb);
    break;
  case 17:
    if (isFirst) printf("Vds        Cgd\n");
    printf("%e %e\n", model->HSMHV_type*vds, here->HSMHV_cgdb);
    break;
  case 18:
    if (isFirst) printf("Vds        Cgb\n");
    printf("%e %e\n", model->HSMHV_type*vds, -(here->HSMHV_cggb+here->HSMHV_cgsb+here->HSMHV_cgdb));
    break;
  case 19:
    if (isFirst) printf("Vgs        Csg\n");
    printf("%e %e\n", model->HSMHV_type*vgs, -(here->HSMHV_cggb+here->HSMHV_cbgb+here->HSMHV_cdgb));
    break;
  case 20:
    if (isFirst) printf("Vgs        Cdg\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_cdgb);
    break;
  case 21:
    if (isFirst) printf("Vgs        Cbg\n");
    printf("%e %e\n", model->HSMHV_type*vgs, here->HSMHV_cbgb);
    break;
  case 22:
    if (isFirst) printf("Vgb        Cgb\n");
    printf("%e %e\n", model->HSMHV_type*vgb, -(here->HSMHV_cggb+here->HSMHV_cgsb+here->HSMHV_cgdb));
    break;
  case 50:
    if (isFirst) printf("Vgs  Vds  Vbs  Vgb  Ids  log10(|Ids|)  gm/|Ids|  gm  gds  gbs  Cgg  Cgs  Cgb  Cgd  Csg  Cbg  Cdg\n");
    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
            model->HSMHV_type*vgs, model->HSMHV_type*vds, model->HSMHV_type*vbs, model->HSMHV_type*vgb, here->HSMHV_mode*here->HSMHV_ids,
            log10(here->HSMHV_ids), here->HSMHV_gm/here->HSMHV_ids, here->HSMHV_gm, here->HSMHV_gds, here->HSMHV_gmbs, here->HSMHV_cggb,
            here->HSMHV_cgsb, -(here->HSMHV_cggb+here->HSMHV_cgsb+here->HSMHV_cgdb), here->HSMHV_cgdb,
            -(here->HSMHV_cggb+here->HSMHV_cbgb+here->HSMHV_cdgb), here->HSMHV_cbgb, here->HSMHV_cdgb);
    break;
  default:
    /*
      printf("There is no physical value corrsponding to %d\n", flag);
    */
    break;
  }
}

int HSMHVload(
     GENmodel *inModel,
     register CKTcircuit *ckt)
     /* actually load the current value into the 
      * sparse matrix previously provided 
      */
{
  register HSMHVmodel *model = (HSMHVmodel*)inModel;
  register HSMHVinstance *here;
  HSMHVbinningParam *pParam;
/*  HSMHVmodelMKSParam *modelMKS ;*/
/*  HSMHVhereMKSParam  *hereMKS ;*/



  /* node voltages */
#define XDIM 14
  double x[XDIM] ;

  /* branch voltages */
  double vbd=0.0,  vbs=0.0,  vds=0.0,  vgb=0.0, vgd=0.0, vgdo=0.0, vgs=0.0 ;
  double vdbs=0.0, vsbs=0.0, vdbd=0.0 ;
  double vges=0.0, vged=0.0, vgedo=0.0 ;
  double vbs_jct=0.0, vbd_jct=0.0;
  double deltemp = 0.0 , deltemp_old = 0.0 ;
  double vggp=0, vddp=0, vssp=0, vbpb=0, vbpsb=0, vbpdb=0 ;
  double vdse=0.0, vgse=0.0, vbse=0.0 ;
  double vsubs=0.0 ; /* substrate bias */
 
  double ivds=0.0,  ivgs=0.0,  ivbs=0.0 ;
  double ivgse=0.0, ivdse=0.0, ivbse=0.0 ;

  /* delta of branch voltages */
  double delvbs=0.0, delvds=0.0, delvgs=0.0, delvsbs=0.0, delvdbd=0.0, deldeltemp = 0.0 ;
  double delvdse=0.0, delvgse=0.0, delvbse=0.0, delvddp=0.0, delvssp=0.0, delvggp=0.0, delvbpb=0.0, delvbpdb=0.0, delvbpsb =0.0 ;
  double delvsubs=0.0; /* substrate bias change */

  /* branch currents */
  double Ids=0.0, gds=0.0,     gm=0.0,     gmbs=0.0,    gmT=0.0,
                  gds_ext=0.0, gm_ext=0.0, gmbs_ext=0.0 ;
  double Igd=0.0, dIgd_dVd=0.0, dIgd_dVg=0.0, dIgd_dVb=0.0, dIgd_dVs=0.0, dIgd_dT=0.0 ;
  double Igs=0.0, dIgs_dVd=0.0, dIgs_dVg=0.0, dIgs_dVb=0.0, dIgs_dVs=0.0, dIgs_dT=0.0 ;
  double Igb=0.0, dIgb_dVd=0.0, dIgb_dVg=0.0, dIgb_dVb=0.0, dIgb_dVs=0.0, dIgb_dT=0.0 ;
  double Isub=0.0, dIsub_dVds=0.0, dIsub_dVgs=0.0, dIsub_dVbs=0.0, dIsub_dT=0.0 ;
  double Isubs=0.0, dIsubs_dVds=0.0, dIsubs_dVgs=0.0, dIsubs_dVbs=0.0, dIsubs_dT=0.0 ;
  double dIsub_dVdse=0.0, dIsubs_dVdse=0.0 ;
  double Igidl=0.0, dIgidl_dVds=0.0, dIgidl_dVgs=0.0, dIgidl_dVbs=0.0, dIgidl_dT=0.0 ;
  double Igisl=0.0, dIgisl_dVds=0.0, dIgisl_dVgs=0.0, dIgisl_dVbs=0.0, dIgisl_dT=0.0 ;
  double Ibd=0.0, Gbd=0.0, Gbdt=0.0 ;
  double Ibs=0.0, Gbs=0.0, Gbst=0.0 ;
  double Iddp=0.0, dIddp_dVddp=0.0, dIddp_dVdse=0.0, dIddp_dVgse=0.0, dIddp_dVbse=0.0, dIddp_dVsubs=0.0, dIddp_dT =0.0 ;
  double Issp=0.0, dIssp_dVssp=0.0, dIssp_dVdse=0.0, dIssp_dVgse=0.0, dIssp_dVbse=0.0, dIssp_dVsubs=0.0, dIssp_dT =0.0 ;
  double Iggp=0.0, dIggp_dVggp =0.0 ;
  double Ibpb=0.0, dIbpb_dVbpb =0.0 ;
  double Ibpdb=0.0, dIbpdb_dVbpdb =0.0 ;
  double Ibpsb=0.0, dIbpsb_dVbpsb =0.0 ;
  double Ith=0.0,   dIth_dT =0.0 ;

  /* displacement currents */
  double cq_d=0.0,  cq_dP=0.0,  cq_g=0.0,  cq_gP=0.0,  cq_s=0.0,  cq_sP=0.0,  cq_bP=0.0,  cq_b=0.0,  cq_db=0.0,    cq_sb=0.0,  cq_t=0.0  ; 
  double cq_dE=0.0, cq_gE=0.0,  cq_sE=0.0,  cq_bE=0.0 ;

  /* node currents */
  double cur_d=0.0, cur_dP=0.0, cur_g=0.0, cur_gP=0.0, cur_s=0.0, cur_sP=0.0, cur_bP=0.0, cur_b=0.0, cur_db=0.0,   cur_sb=0.0, cur_t=0.0 ;
  double i_d=0.0,   i_dP=0.0,   i_g=0.0,   i_gP=0.0,   i_s=0.0,   i_sP=0.0,   i_bP=0.0,   i_b=0.0,   i_db=0.0,     i_sb=0.0,   i_t=0.0   ;

  /* resistances and conductances */
  double Rd=0.0, dRd_dVdse=0.0, dRd_dVgse=0.0, dRd_dVbse=0.0, dRd_dVsubs=0.0, dRd_dT=0.0 ;
  double Rs=0.0, dRs_dVdse=0.0, dRs_dVgse=0.0, dRs_dVbse=0.0, dRs_dVsubs=0.0, dRs_dT=0.0 ;

  double GD=0.0, GD_dVds=0.0, GD_dVgs=0.0, GD_dVbs=0.0, GD_dVsubs=0.0, GD_dT=0.0 ;
  double GS=0.0, GS_dVds=0.0, GS_dVgs=0.0, GS_dVbs=0.0, GS_dVsubs=0.0, GS_dT=0.0 ;
  double Gth=0.0 ;
  double GG=0.0,  GRBPD=0.0, GRBPS=0.0, GRBPB=0.0; 
 
  /* charges */
  double Qd=0.0, dQd_dVds=0.0, dQd_dVgs=0.0, dQd_dVbs=0.0, dQd_dT=0.0 ;
  double Qg=0.0, dQg_dVds=0.0, dQg_dVgs=0.0, dQg_dVbs=0.0, dQg_dT=0.0 ;
  double Qs=0.0, dQs_dVds=0.0, dQs_dVgs=0.0, dQs_dVbs=0.0, dQs_dT=0.0 ;
  double Qb=0.0, dQb_dVds=0.0, dQb_dVgs=0.0, dQb_dVbs=0.0, dQb_dT=0.0 ;
  double Qbd=0.0, Cbd=0.0, Cbdt=0.0,
         Qbs=0.0, Cbs=0.0, Cbst=0.0 ;
  double Qth=0.0, Cth =0.0 ;
  double Qfd=0.0, dQfd_dVdse=0.0, dQfd_dVgse=0.0, dQfd_dVbse=0.0, dQfd_dT=0.0 ;
  double Qfs=0.0, dQfs_dVdse=0.0, dQfs_dVgse=0.0, dQfs_dVbse=0.0, dQfs_dT=0.0 ;

  double Qdext=0.0, dQdext_dVdse=0.0, dQdext_dVgse=0.0, dQdext_dVbse=0.0, dQdext_dT=0.0 ;
  double Qgext=0.0, dQgext_dVdse=0.0, dQgext_dVgse=0.0, dQgext_dVbse=0.0, dQgext_dT=0.0 ;
  double /*Qsext=0.0,*/ dQsext_dVdse=0.0, dQsext_dVgse=0.0, dQsext_dVbse=0.0, dQsext_dT=0.0 ;
  double Qbext=0.0, dQbext_dVdse=0.0, dQbext_dVgse=0.0, dQbext_dVbse=0.0, dQbext_dT=0.0 ;
  /* 5th substrate node */
  int flg_subNode = 0 ;
  
  /* self heating */
  double Veffpower=0.0, dVeffpower_dVds=0.0, dVeffpower_dVdse =0.0 ;
  double P=0.0, dP_dVds=0.0,  dP_dVgs=0.0,  dP_dVbs=0.0, dP_dT =0.0,
                dP_dVdse=0.0, dP_dVgse=0.0, dP_dVbse =0.0 ;
  int flg_tempNode = 0 ;
  double T0 , T1 , T2 ;

#define SHE_MAX_dlt 0.1 

#define C_RTH_MIN 1.0e-4

  double TMF1 , TMF2  ;
/*---------------------------------------------------*
* smoothUpper: ceiling.
*      y = xmax - 0.5 ( arg + sqrt( arg^2 + 4 xmax delta ) )
*    arg = xmax - x - delta
*-----------------*/

#define Fn_SU( y , x , xmax , delta , dx ) { \
    TMF1 = ( xmax ) - ( x ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmax ) * ( delta) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 :  -( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 *  TMF1 + TMF2 ) ; \
    dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    y = ( xmax ) - 0.5 * ( TMF1 + TMF2 ) ; \
  }

  /* NQS related variables */
  int flg_nqs =0 ;
  double Qi_nqs=0.0, Qb_nqs=0.0, delQi_nqs=0.0, delQb_nqs=0.0, i_qi=0.0, i_qb=0.0, cq_qi=0.0, cq_qb=0.0, cur_qi=0.0, cur_qb =0.0 ;
  double Iqi_nqs=0.0, dIqi_nqs_dVds=0.0, dIqi_nqs_dVgs=0.0, dIqi_nqs_dVbs=0.0, dIqi_nqs_dT=0.0, dIqi_nqs_dQi_nqs =0.0 ;
  double Iqb_nqs=0.0, dIqb_nqs_dVds=0.0, dIqb_nqs_dVgs=0.0, dIqb_nqs_dVbs=0.0, dIqb_nqs_dT=0.0, dIqb_nqs_dQb_nqs =0.0 ;
  double Qd_nqs=0.0,  dQd_nqs_dVds=0.0,  dQd_nqs_dVgs=0.0,  dQd_nqs_dVbs=0.0,  dQd_nqs_dT=0.0,  dQd_nqs_dQi_nqs  =0.0 ;
  double Qs_nqs=0.0,  dQs_nqs_dVds=0.0,  dQs_nqs_dVgs=0.0,  dQs_nqs_dVbs=0.0,  dQs_nqs_dT=0.0,  dQs_nqs_dQi_nqs  =0.0 ;
  double Qg_nqs=0.0,  dQg_nqs_dQi_nqs=0.0, dQg_nqs_dQb_nqs =0.0 ;
  double tau=0.0,     dtau_dVds=0.0,     dtau_dVgs=0.0,     dtau_dVbs=0.0,     dtau_dT   =0.0 ;
  double taub=0.0,    dtaub_dVds=0.0,    dtaub_dVgs=0.0,    dtaub_dVbs=0.0,    dtaub_dT  =0.0 ;
  double Qdrat=0.0,   dQdrat_dVds=0.0,   dQdrat_dVgs=0.0,   dQdrat_dVbs=0.0,   dQdrat_dT =0.0 ;
  double Qi=0.0,      dQi_dVds=0.0,      dQi_dVgs=0.0,      dQi_dVbs=0.0,      dQi_dT    =0.0 ;
  double Qbulk=0.0,   dQbulk_dVds=0.0,   dQbulk_dVgs=0.0,   dQbulk_dVbs=0.0,   dQbulk_dT =0.0 ;

  /* output related variables */
  double dQi_nqs_dVds=0.0, dQi_nqs_dVgs=0.0, dQi_nqs_dVbs=0.0,
         dQb_nqs_dVds=0.0, dQb_nqs_dVgs=0.0, dQb_nqs_dVbs=0.0 ;
  double cgdb_nqs=0.0,     cggb_nqs=0.0,     cgsb_nqs=0.0,
         cbdb_nqs=0.0,     cbgb_nqs=0.0,     cbsb_nqs=0.0,
         cddb_nqs=0.0,     cdgb_nqs=0.0,     cdsb_nqs=0.0 ;
  double cgdb=0.0,         cggb=0.0,         cgsb=0.0,
         cbdb=0.0,         cbgb=0.0,         cbsb=0.0,
         cddb=0.0,         cdgb=0.0,         cdsb=0.0 ;

  /* rows of conductance and capacitance matrix stamp */
  double ydc_d[XDIM],   ydc_dP[XDIM],  ydc_g[XDIM],   ydc_gP[XDIM],  ydc_s[XDIM],  ydc_sP[XDIM] ; 
  double ydc_bP[XDIM],  ydc_b[XDIM],   ydc_db[XDIM],  ydc_sb[XDIM],  ydc_t[XDIM],  ydc_qi[XDIM],  ydc_qb[XDIM]  ;
  double ydyn_d[XDIM],  ydyn_dP[XDIM], ydyn_g[XDIM],  ydyn_gP[XDIM], ydyn_s[XDIM], ydyn_sP[XDIM] ; 
  double ydyn_bP[XDIM], ydyn_b[XDIM],  ydyn_db[XDIM], ydyn_sb[XDIM], ydyn_t[XDIM], ydyn_qi[XDIM], ydyn_qb[XDIM] ;

  /* limiter, bypass, and convergence */
  int ByPass=0, Check=0, Check1=0, Check2=0, Check3=0 ;
  double von=0.0, limval =0.0 ;
  double i_dP_hat=0.0, i_gP_hat=0.0, i_sP_hat=0.0, i_db_hat=0.0, i_sb_hat =0.0         ;

#define LIM_TOL 1.0e0
#define LIM_TOL2 1.0e0

  /* predictor and numerical integration stuff */
  double ag0=0.0, xfact=0.0 ;
  double ceq=0.0, geq=0.0 ;
  int ChargeComputationNeeded =  
    ((ckt->CKTmode & (MODEAC | MODETRAN | MODEINITSMSIG)) ||
     ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)))
    ? 1 : 0;
  int showPhysVal=0 ;
  int isConv=0 ;
  double vds_pre=0.0 ;
  int i=0, noncon_old=0 ;

  

#define dNode        0
#define dNodePrime   1
#define gNode        2
#define gNodePrime   3
#define sNode        4
#define sNodePrime   5
#define bNodePrime   6
#define bNode        7 
#define dbNode       8
#define sbNode       9
#define subNode     10
#define tempNode    11
#define qiNode      12
#define qbNode      13
#define lastNode    13        /* must be the last node! */


#define SPICE_rhs 1           /* set to 0 if rhs to be loaded for standard Newton */
                              /* set to 1 if rhs_eq to be loaded, e.g. for SPICE  */

  if (SPICE_rhs) {
    for (i=0; i<XDIM; i++) {
      x[i] = 0.0 ;
    }
  }

  /*  loop through all the HSMHV device models */
  for ( ; model != NULL; model = model->HSMHVnextModel ) {
    /* loop through all the instances of the model */

/*    modelMKS = &model->modelMKS ;*/

    for (here = model->HSMHVinstances; here != NULL ;
         here = here->HSMHVnextInstance) {
      
/*      hereMKS = &here->hereMKS ;*/
      pParam = &here->pParam ;
      showPhysVal = 0;
      Check=1;
      ByPass = 0;
      vsubs = 0.0 ; /* substrate bias */
      deltemp = 0.0 ;
      noncon_old = ckt->CKTnoncon;
      flg_nqs = model->HSMHV_conqs ;
      flg_subNode = here->HSMHVsubNode ; /* if flg_subNode > 0, external(/internal) substrate node exists */
      flg_tempNode = here->HSMHVtempNode ; /* if flg_tempNode > 0, external/internal temperature node exists */

#ifdef DEBUG_HISIMHVLD_VX
      printf("mode = %x\n", ckt->CKTmode);
      printf("Vd Vg Vs Vb %e %e %e %e\n", *(ckt->CKTrhsOld+here->HSMHVdNodePrime),
             *(ckt->CKTrhsOld+here->HSMHVgNodePrime),
             *(ckt->CKTrhsOld+here->HSMHVsNodePrime),
             *(ckt->CKTrhsOld+here->HSMHVbNodePrime));
#endif

      if ( ckt->CKTmode & MODEINITSMSIG ) {
        vbs = *(ckt->CKTstate0 + here->HSMHVvbs);
        vgs = *(ckt->CKTstate0 + here->HSMHVvgs);
        vds = *(ckt->CKTstate0 + here->HSMHVvds);

        vges = *(ckt->CKTstate0 + here->HSMHVvges);
        vdbd = *(ckt->CKTstate0 + here->HSMHVvdbd);
        vsbs = *(ckt->CKTstate0 + here->HSMHVvsbs);
	if (flg_subNode > 0) vsubs = *(ckt->CKTstate0 + here->HSMHVvsubs);
 	if( flg_tempNode > 0 ){
	  deltemp = *(ckt->CKTstate0 + here->HSMHVdeltemp);
	}
 	vdse = *(ckt->CKTstate0 + here->HSMHVvdse) ;
	vgse = *(ckt->CKTstate0 + here->HSMHVvgse) ;
	vbse = *(ckt->CKTstate0 + here->HSMHVvbse) ;
        if ( flg_nqs ) {
          Qi_nqs = *(ckt->CKTstate0 + here->HSMHVqi_nqs) ;
          Qb_nqs = *(ckt->CKTstate0 + here->HSMHVqb_nqs) ;
        } else {
          Qi_nqs = 0.0 ;
          Qb_nqs = 0.0 ;
        }
      /* printf("HSMHV_load: (from state0) vds.. = %e %e %e %e %e %e\n",
                                              vds,vgs,vbs,vdse,vgse,vbse); */
      } 
      else if ( ckt->CKTmode & MODEINITTRAN ) {
/* #include "printf_ld_converged.inc" */
        vbs = *(ckt->CKTstate1 + here->HSMHVvbs);
        vgs = *(ckt->CKTstate1 + here->HSMHVvgs);
        vds = *(ckt->CKTstate1 + here->HSMHVvds);

        vges = *(ckt->CKTstate1 + here->HSMHVvges);
        vdbd = *(ckt->CKTstate1 + here->HSMHVvdbd);
        vsbs = *(ckt->CKTstate1 + here->HSMHVvsbs);
	if (flg_subNode > 0) vsubs = *(ckt->CKTstate1 + here->HSMHVvsubs);
 	if( flg_tempNode > 0 ){
	  deltemp = *(ckt->CKTstate1 + here->HSMHVdeltemp);
	}
 	vdse = *(ckt->CKTstate1 + here->HSMHVvdse) ;
        vgse = *(ckt->CKTstate1 + here->HSMHVvgse) ;
        vbse = *(ckt->CKTstate1 + here->HSMHVvbse) ;
        if ( flg_nqs ) {
          Qi_nqs = *(ckt->CKTstate1 + here->HSMHVqi_nqs) ;
          Qb_nqs = *(ckt->CKTstate1 + here->HSMHVqb_nqs) ;
        } else {
          Qi_nqs = 0.0 ;
          Qb_nqs = 0.0 ;
        }
      } 
      else if ( (ckt->CKTmode & MODEINITJCT) && !here->HSMHV_off ) {
        vds = model->HSMHV_type * here->HSMHV_icVDS;
        vgs = vges = model->HSMHV_type * here->HSMHV_icVGS;
        vbs = vsbs = model->HSMHV_type * here->HSMHV_icVBS;
        vdbd = 0.0 ;
        if ( (vds == 0.0) && (vgs == 0.0) && (vbs == 0.0) && 
             ( (ckt->CKTmode & (MODETRAN|MODEAC|MODEDCOP|MODEDCTRANCURVE)) ||
               !(ckt->CKTmode & MODEUIC) ) ) { 
          /* set biases for starting analysis */
          vbs = vdbd = vsbs = 0.0;
          vgs = vges = 0.1;
          vds = 0.1;
        }
	if (flg_subNode > 0) vsubs = 0.0;
        if( flg_tempNode > 0 ) deltemp=0.0;
        vdse = vds ;
        vgse = vgs ;
        Qi_nqs = Qb_nqs = 0.0 ;
      } 
      else if ( ( ckt->CKTmode & (MODEINITJCT | MODEINITFIX) ) && 
                here->HSMHV_off ) {
        vbs = vgs = vds = 0.0; vges = 0.0; vdbd = vsbs = 0.0;
	if (flg_subNode > 0) vsubs = 0.0;
        if( flg_tempNode > 0 ) deltemp=0.0;
        vdse = vds ;
        vgse = vgs ;
        Qi_nqs = Qb_nqs = 0.0 ;
      } 
      else {
#ifndef PREDICTOR /* BSIM3 style */
        if (ckt->CKTmode & MODEINITPRED) {
/* #include "printf_ld_converged.inc" */
          /* if (here->HSMHV_mode > 0) {
             gds_ext = here->HSMHV_dIds_dVdse ;
          } else {
             gds_ext = + (here->HSMHV_dIds_dVdse + here->HSMHV_dIds_dVgse + here->HSMHV_dIds_dVbse) ;
          }
          printf("zzz %e %e\n",ckt->CKTtime,gds_ext) ; */
          xfact = ckt->CKTdelta / ckt->CKTdeltaOld[1];
          *(ckt->CKTstate0 + here->HSMHVvbs) = 
            *(ckt->CKTstate1 + here->HSMHVvbs);
          vbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvbs))
            -(xfact * (*(ckt->CKTstate2 + here->HSMHVvbs)));
          *(ckt->CKTstate0 + here->HSMHVvgs) = 
            *(ckt->CKTstate1 + here->HSMHVvgs);
          vgs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvgs))
            -(xfact * (*(ckt->CKTstate2 + here->HSMHVvgs)));
          *(ckt->CKTstate0 + here->HSMHVvds) = 
            *(ckt->CKTstate1 + here->HSMHVvds);
          vds = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvds))
            -(xfact * (*(ckt->CKTstate2 + here->HSMHVvds)));
          *(ckt->CKTstate0 + here->HSMHVvbd) = 
            *(ckt->CKTstate0 + here->HSMHVvbs)-
            *(ckt->CKTstate0 + here->HSMHVvds);

          *(ckt->CKTstate0 + here->HSMHVvges) = 
            *(ckt->CKTstate1 + here->HSMHVvges);
          vges = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvges))
            -(xfact * (*(ckt->CKTstate2 + here->HSMHVvges)));
          *(ckt->CKTstate0 + here->HSMHVvdbd) =
            *(ckt->CKTstate1 + here->HSMHVvdbd);
          vdbd = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvdbd))
            - (xfact * (*(ckt->CKTstate2 + here->HSMHVvdbd)));
          *(ckt->CKTstate0 + here->HSMHVvdbs) =
            *(ckt->CKTstate0 + here->HSMHVvdbd)
            + *(ckt->CKTstate0 + here->HSMHVvds);
          *(ckt->CKTstate0 + here->HSMHVvsbs) =
            *(ckt->CKTstate1 + here->HSMHVvsbs);
          vsbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvsbs))
            - (xfact * (*(ckt->CKTstate2 + here->HSMHVvsbs)));
	  if (flg_subNode > 0){
            vsubs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvsubs))
              - ( xfact * (*(ckt->CKTstate2 + here->HSMHVvsubs)));
	  }
          if( flg_tempNode > 0 ){
            deltemp = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVdeltemp))
              - ( xfact * (*(ckt->CKTstate2 + here->HSMHVdeltemp)));
            
            *(ckt->CKTstate0 + here->HSMHVdeltemp) = 
              *(ckt->CKTstate1 + here->HSMHVdeltemp);
          }
          *(ckt->CKTstate0 + here->HSMHVvdse) = 
            *(ckt->CKTstate1 + here->HSMHVvdse);
	  vdse = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvdse))
	    -(xfact * (*(ckt->CKTstate2 + here->HSMHVvdse)));
	  *(ckt->CKTstate0 + here->HSMHVvgse) = 
	    *(ckt->CKTstate1 + here->HSMHVvgse);
	  vgse = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvgse))
	    -(xfact * (*(ckt->CKTstate2 + here->HSMHVvgse)));
	  *(ckt->CKTstate0 + here->HSMHVvbse) = 
	    *(ckt->CKTstate1 + here->HSMHVvbse);
	  vbse = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVvbse))
	    -(xfact * (*(ckt->CKTstate2 + here->HSMHVvbse)));
          if (flg_nqs) {
	    *(ckt->CKTstate0 + here->HSMHVqi_nqs) = 
	      *(ckt->CKTstate1 + here->HSMHVqi_nqs);
	    Qi_nqs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVqi_nqs))
	    -(xfact * (*(ckt->CKTstate2 + here->HSMHVqi_nqs)));
	    *(ckt->CKTstate0 + here->HSMHVqb_nqs) = 
	      *(ckt->CKTstate1 + here->HSMHVqb_nqs);
	    Qb_nqs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSMHVqb_nqs))
	      -(xfact * (*(ckt->CKTstate2 + here->HSMHVqb_nqs)));
          } else {
            Qi_nqs = Qb_nqs = 0.0 ;
          }
        } 
        else {
#endif /* PREDICTOR */
/* 	  here->HSMHV_time = ckt->CKTtime; /\* for debug print *\/ */
          /* get biases from CKT */
          vbs = model->HSMHV_type * 
            (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) -
             *(ckt->CKTrhsOld+here->HSMHVsNodePrime));
          vgs = model->HSMHV_type * 
            (*(ckt->CKTrhsOld+here->HSMHVgNodePrime) -
             *(ckt->CKTrhsOld+here->HSMHVsNodePrime));
          vds = model->HSMHV_type * 
            (*(ckt->CKTrhsOld+here->HSMHVdNodePrime) -
             *(ckt->CKTrhsOld+here->HSMHVsNodePrime));

          vges = model->HSMHV_type * 
            (*(ckt->CKTrhsOld+here->HSMHVgNode) -
             *(ckt->CKTrhsOld+here->HSMHVsNodePrime));
          vdbd = model->HSMHV_type
            * (*(ckt->CKTrhsOld + here->HSMHVdbNode)
               - *(ckt->CKTrhsOld + here->HSMHVdNode));
          vsbs = model->HSMHV_type
            * (*(ckt->CKTrhsOld + here->HSMHVsbNode)
               - *(ckt->CKTrhsOld + here->HSMHVsNode));
	  if (flg_subNode > 0){
	    vsubs = model->HSMHV_type
	      * (*(ckt->CKTrhsOld + here->HSMHVsubNode)
		 - *(ckt->CKTrhsOld + here->HSMHVsNode));
	  }
	  if( flg_tempNode > 0 ){
	      deltemp = *(ckt->CKTrhsOld + here->HSMHVtempNode);
	  }
	  vbse = model->HSMHV_type * 
	      (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) -
	       *(ckt->CKTrhsOld+here->HSMHVsNode));
	  vgse = model->HSMHV_type * 
	      (*(ckt->CKTrhsOld+here->HSMHVgNodePrime) -
	       *(ckt->CKTrhsOld+here->HSMHVsNode));
	  vdse = model->HSMHV_type * 
	      (*(ckt->CKTrhsOld+here->HSMHVdNode) -
	       *(ckt->CKTrhsOld+here->HSMHVsNode));
          if ( flg_nqs ) {
            Qi_nqs = *(ckt->CKTrhsOld + here->HSMHVqiNode);
            Qb_nqs = *(ckt->CKTrhsOld + here->HSMHVqbNode);
          } else {
            Qi_nqs = Qb_nqs = 0.0 ;
          }	    
#ifndef PREDICTOR
        }
#endif /* PREDICTOR */

        /* printf("HSMHV_load: (from rhs   ) vds.. = %e %e %e %e %e %e\n",
                                                 vds,vgs,vbs,vdse,vgse,vbse); */

        vbd = vbs - vds;
        vgd = vgs - vds;
        vged = vges - vds;
        vdbs = vdbd + vdse;
        vgdo = *(ckt->CKTstate0 + here->HSMHVvgs) - *(ckt->CKTstate0 + here->HSMHVvds);
        vgedo = *(ckt->CKTstate0 + here->HSMHVvges) - *(ckt->CKTstate0 + here->HSMHVvds);

        vds_pre = vds;

#ifndef NOBYPASS
        /* start of bypass section
           ... no bypass in case of selfheating  */
        if ( !(ckt->CKTmode & MODEINITPRED) && ckt->CKTbypass && !model->HSMHV_coselfheat) {
          delvds  = vds  - *(ckt->CKTstate0 + here->HSMHVvds) ;
          delvgs  = vgs  - *(ckt->CKTstate0 + here->HSMHVvgs) ;
          delvbs  = vbs  - *(ckt->CKTstate0 + here->HSMHVvbs) ;
          delvdse = vdse - *(ckt->CKTstate0 + here->HSMHVvdse) ;
          delvgse = vgse - *(ckt->CKTstate0 + here->HSMHVvgse) ;
          delvbse = vbse - *(ckt->CKTstate0 + here->HSMHVvbse) ;
          delvdbd = vdbd - *(ckt->CKTstate0 + here->HSMHVvdbd) ;
          delvsbs = vsbs - *(ckt->CKTstate0 + here->HSMHVvsbs) ;
          if (flg_subNode > 0) delvsubs = vsubs - *(ckt->CKTstate0 + here->HSMHVvsubs) ; /* substrate bias change */
          deldeltemp = deltemp - *(ckt->CKTstate0 + here->HSMHVdeltemp) ;
          if (flg_nqs) {
            delQi_nqs = Qi_nqs - *(ckt->CKTstate0 + here->HSMHVqi_nqs) ;
            delQb_nqs = Qb_nqs - *(ckt->CKTstate0 + here->HSMHVqb_nqs) ;
          } else {
            delQi_nqs = delQb_nqs = 0.0 ;
          }

          /* now let's see if we can bypass                     */
          /* ... first perform the easy cheap bypass checks ... */
 /*          1 2     3       3                       3    4    4     4 5                               543                   2    1 */
          if ( ( fabs(delvds ) < ckt->CKTreltol * MAX(fabs(vds ),fabs(*(ckt->CKTstate0 + here->HSMHVvds ))) + ckt->CKTvoltTol ) &&
               ( fabs(delvgs ) < ckt->CKTreltol * MAX(fabs(vgs ),fabs(*(ckt->CKTstate0 + here->HSMHVvgs ))) + ckt->CKTvoltTol ) &&
               ( fabs(delvbs ) < ckt->CKTreltol * MAX(fabs(vbs ),fabs(*(ckt->CKTstate0 + here->HSMHVvbs ))) + ckt->CKTvoltTol ) &&
               ( fabs(delvdse) < ckt->CKTreltol * MAX(fabs(vdse),fabs(*(ckt->CKTstate0 + here->HSMHVvdse))) + ckt->CKTvoltTol ) &&
               ( fabs(delvgse) < ckt->CKTreltol * MAX(fabs(vgse),fabs(*(ckt->CKTstate0 + here->HSMHVvgse))) + ckt->CKTvoltTol ) &&
               ( fabs(delvbse) < ckt->CKTreltol * MAX(fabs(vbse),fabs(*(ckt->CKTstate0 + here->HSMHVvbse))) + ckt->CKTvoltTol ) &&
               ( fabs(delvdbd) < ckt->CKTreltol * MAX(fabs(vdbd),fabs(*(ckt->CKTstate0 + here->HSMHVvdbd))) + ckt->CKTvoltTol ) &&
               ( fabs(delvsbs) < ckt->CKTreltol * MAX(fabs(vsbs),fabs(*(ckt->CKTstate0 + here->HSMHVvsbs))) + ckt->CKTvoltTol ) &&
               ( fabs(delvsubs) < ckt->CKTreltol * MAX(fabs(vsubs),fabs(*(ckt->CKTstate0 + here->HSMHVvsubs))) + ckt->CKTvoltTol ) &&
               ( fabs(delQi_nqs) < ckt->CKTreltol *   fabs(Qi_nqs) + ckt->CKTchgtol*ckt->CKTabstol + 1.0e-20  ) &&
               ( fabs(delQb_nqs) < ckt->CKTreltol *   fabs(Qb_nqs) + ckt->CKTchgtol*ckt->CKTabstol + 1.0e-20  )    )
                                                                                  /* 1.0e-20: heuristic value, must be small enough     */
                                                                                  /* to ensure that bypass does not destroy convergence */
          { /* ... the first bypass checks are passed -> now do the more expensive checks ...*/
            if ( here->HSMHV_mode > 0 ) { /* forward mode */
              Ids        = here->HSMHV_ids ;
              gds        = here->HSMHV_dIds_dVdsi ;
              gm         = here->HSMHV_dIds_dVgsi ;
              gmbs       = here->HSMHV_dIds_dVbsi ;
              gmT        = (flg_tempNode > 0) ? here->HSMHV_dIds_dTi : 0.0  ;
	      gmbs_ext   = here->HSMHV_dIds_dVbse;
	      gds_ext    = here->HSMHV_dIds_dVdse ;
	      gm_ext     = here->HSMHV_dIds_dVgse;
	      Isub         = here->HSMHV_isub ;
	      dIsub_dVds   = here->HSMHV_dIsub_dVdsi ; 
	      dIsub_dVgs   = here->HSMHV_dIsub_dVgsi ; 
	      dIsub_dVbs   = here->HSMHV_dIsub_dVbsi ;
              dIsub_dT     = (flg_tempNode > 0) ? here->HSMHV_dIsub_dTi : 0.0  ; 
	      dIsub_dVdse  = here->HSMHV_dIsub_dVdse ;
	      Isubs        = 0.0 ;
	      dIsubs_dVds  = 0.0 ; 
	      dIsubs_dVgs  = 0.0 ; 
	      dIsubs_dVbs  = 0.0 ;
              dIsubs_dT    = 0.0 ;
	      dIsubs_dVdse = 0.0 ;
              Igidl        = here->HSMHV_igidl ;
              dIgidl_dVds  = here->HSMHV_dIgidl_dVdsi ;
              dIgidl_dVgs  = here->HSMHV_dIgidl_dVgsi ;
              dIgidl_dVbs  = here->HSMHV_dIgidl_dVbsi ;
              dIgidl_dT    = (flg_tempNode > 0) ? here->HSMHV_dIgidl_dTi : 0.0  ;
              Igisl        = here->HSMHV_igisl ;
              dIgisl_dVds  = here->HSMHV_dIgisl_dVdsi ;
              dIgisl_dVgs  = here->HSMHV_dIgisl_dVgsi ;
              dIgisl_dVbs  = here->HSMHV_dIgisl_dVbsi ;
              dIgisl_dT    = (flg_tempNode > 0) ? here->HSMHV_dIgisl_dTi : 0.0  ;
              Igd          = here->HSMHV_igd ;
              dIgd_dVd   = here->HSMHV_dIgd_dVdsi ;
              dIgd_dVg   = here->HSMHV_dIgd_dVgsi ;
              dIgd_dVb   = here->HSMHV_dIgd_dVbsi ;
              dIgd_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgd_dTi : 0.0  ;
              Igs          = here->HSMHV_igs ;
              dIgs_dVd   = here->HSMHV_dIgs_dVdsi ;
              dIgs_dVg   = here->HSMHV_dIgs_dVgsi ;
              dIgs_dVb   = here->HSMHV_dIgs_dVbsi ;
              dIgs_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgs_dTi : 0.0  ;
              Igb          = here->HSMHV_igb ;
              dIgb_dVd   = here->HSMHV_dIgb_dVdsi ;
              dIgb_dVg   = here->HSMHV_dIgb_dVgsi ;
              dIgb_dVb   = here->HSMHV_dIgb_dVbsi ;
              dIgb_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgb_dTi : 0.0  ;
	      Ibd = here->HSMHV_ibd ;
              Gbd = here->HSMHV_gbd ;
	      Gbdt = (flg_tempNode > 0) ? here->HSMHV_gbdT : 0.0 ;
	      Ibs = here->HSMHV_ibs ;
	      Gbs = here->HSMHV_gbs ;
	      Gbst = (flg_tempNode > 0) ? here->HSMHV_gbsT : 0.0 ;
            } else { /* reverse mode */
              Ids       = - here->HSMHV_ids ;
              gds       = + (here->HSMHV_dIds_dVdsi + here->HSMHV_dIds_dVgsi + here->HSMHV_dIds_dVbsi) ;
              gm        = - here->HSMHV_dIds_dVgsi ;
              gmbs      = - here->HSMHV_dIds_dVbsi ;
              gmT       = (flg_tempNode > 0) ? - here->HSMHV_dIds_dTi : 0.0  ;
              gds_ext   = + (here->HSMHV_dIds_dVdse + here->HSMHV_dIds_dVgse + here->HSMHV_dIds_dVbse) ;
  	      gm_ext    = - here->HSMHV_dIds_dVgse;
	      gmbs_ext  = - here->HSMHV_dIds_dVbse;
	      Isub         = 0.0 ;
	      dIsub_dVds   = 0.0 ; 
	      dIsub_dVgs   = 0.0 ; 
	      dIsub_dVbs   = 0.0 ;
              dIsub_dT     = 0.0 ; 
	      dIsub_dVdse  = 0.0 ;
	      Isubs        =   here->HSMHV_isub ;
	      dIsubs_dVds  = - (here->HSMHV_dIsub_dVdsi + here->HSMHV_dIsub_dVgsi + here->HSMHV_dIsub_dVbsi) ; 
	      dIsubs_dVgs  =   here->HSMHV_dIsub_dVgsi ; 
	      dIsubs_dVbs  =   here->HSMHV_dIsub_dVbsi ;
              dIsubs_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIsub_dTi : 0.0 ;
	      dIsubs_dVdse = - here->HSMHV_dIsub_dVdse ; /* = - (dIsub_dVdse + dIsub_dVbse + dIsub_dVgse) */
              Igidl        =   here->HSMHV_igisl ;
              dIgidl_dVds  = - (here->HSMHV_dIgisl_dVdsi + here->HSMHV_dIgisl_dVgsi + here->HSMHV_dIgisl_dVbsi) ;
              dIgidl_dVgs  =   here->HSMHV_dIgisl_dVgsi ;
              dIgidl_dVbs  =   here->HSMHV_dIgisl_dVbsi ;
              dIgidl_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIgisl_dTi : 0.0  ;
              Igisl        =   here->HSMHV_igidl ;
              dIgisl_dVds  = - (here->HSMHV_dIgidl_dVdsi + here->HSMHV_dIgidl_dVgsi + here->HSMHV_dIgidl_dVbsi) ;
              dIgisl_dVgs  =   here->HSMHV_dIgidl_dVgsi ;
              dIgisl_dVbs  =   here->HSMHV_dIgidl_dVbsi ;
              dIgisl_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIgidl_dTi : 0.0  ;
              Igd          =   here->HSMHV_igd ;
              dIgd_dVd   = - (here->HSMHV_dIgs_dVdsi + here->HSMHV_dIgs_dVgsi + here->HSMHV_dIgs_dVbsi) ;
              dIgd_dVg   =   here->HSMHV_dIgs_dVgsi ;
              dIgd_dVb   =   here->HSMHV_dIgs_dVbsi ;
              dIgd_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgs_dTi : 0.0  ;
              Igs          =   here->HSMHV_igs ;
              dIgs_dVd   = - (here->HSMHV_dIgd_dVdsi + here->HSMHV_dIgd_dVgsi + here->HSMHV_dIgd_dVbsi) ;
              dIgs_dVg   =   here->HSMHV_dIgd_dVgsi ;
              dIgs_dVb   =   here->HSMHV_dIgd_dVbsi ;
              dIgs_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgd_dTi : 0.0  ;
              Igb          =   here->HSMHV_igb ;
              dIgb_dVd   = - (here->HSMHV_dIgb_dVdsi + here->HSMHV_dIgb_dVgsi + here->HSMHV_dIgb_dVbsi) ;
              dIgb_dVg   =   here->HSMHV_dIgb_dVgsi ;
              dIgb_dVb   =   here->HSMHV_dIgb_dVbsi ;
              dIgb_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgb_dTi : 0.0  ;
	      Ibd = here->HSMHV_ibd ;
	      Gbd = here->HSMHV_gbd ;
	      Gbdt = (flg_tempNode > 0) ? here->HSMHV_gbdT : 0.0 ;
	      Ibs = here->HSMHV_ibs ;
	      Gbs = here->HSMHV_gbs ;
	      Gbst = (flg_tempNode > 0) ? here->HSMHV_gbsT : 0.0 ;
            } /* end of reverse mode */

            /* for bypass control, only nonlinear static currents are considered: */
            i_dP     = Ids  + Isub  + Igidl - Igd ;
            i_dP_hat = i_dP + gm         *delvgs + gds        *delvds + gmbs       *delvbs + gmT      *deldeltemp
                            + dIsub_dVgs *delvgs + dIsub_dVds *delvds + dIsub_dVbs *delvbs + dIsub_dT *deldeltemp
                            + dIsub_dVdse*delvdse 
                            + dIgidl_dVgs*delvgs + dIgidl_dVds*delvds + dIgidl_dVbs*delvbs + dIgidl_dT*deldeltemp
                            -(dIgd_dVg   *delvgs + dIgd_dVd   *delvds + dIgd_dVb   *delvbs + dIgd_dT  *deldeltemp)
                            + gm_ext     *delvgse + gds_ext    *delvdse + gmbs_ext   *delvbse ;

            i_gP     = Igd  + Igs   + Igb ;
            i_gP_hat = i_gP + dIgd_dVg   *delvgs + dIgd_dVd   *delvds + dIgd_dVb   *delvbs + dIgd_dT  *deldeltemp
                            + dIgs_dVg   *delvgs + dIgs_dVd   *delvds + dIgs_dVb   *delvbs + dIgs_dT  *deldeltemp
                            + dIgb_dVg   *delvgs + dIgb_dVd   *delvds + dIgb_dVb   *delvbs + dIgb_dT  *deldeltemp ;

            i_sP     =-Ids  + Isubs + Igisl - Igs ;
            i_sP_hat = i_sP -(gm         *delvgs + gds        *delvds + gmbs       *delvbs + gmT      *deldeltemp)
                            + dIsubs_dVgs*delvgs + dIsubs_dVds*delvds + dIsubs_dVbs*delvbs + dIsubs_dT*deldeltemp
                            + dIsubs_dVdse*delvdse 
                            + dIgisl_dVgs*delvgs + dIgisl_dVds*delvds + dIgisl_dVbs*delvbs + dIgisl_dT*deldeltemp
                            -(dIgs_dVg   *delvgs + dIgs_dVd   *delvds + dIgs_dVb   *delvbs + dIgs_dT  *deldeltemp)
                            -(gm_ext     *delvgse + gds_ext    *delvdse + gmbs_ext   *delvbse) ;

            i_db     = Ibd ;
            i_db_hat = i_db + Gbd*delvdbd + Gbdt*deldeltemp ;

            i_sb     = Ibs ;
            i_sb_hat = i_sb + Gbs*delvsbs + Gbst*deldeltemp ;

            /* ... second part of bypass checks: */
 /*            1 2     3               3                       3    4        4     4    43                  2    1 */
            if ( ( fabs(i_dP_hat - i_dP) < ckt->CKTreltol * MAX(fabs(i_dP_hat),fabs(i_dP)) + ckt->CKTabstol ) && 
                 ( fabs(i_gP_hat - i_gP) < ckt->CKTreltol * MAX(fabs(i_gP_hat),fabs(i_gP)) + ckt->CKTabstol ) &&
                 ( fabs(i_sP_hat - i_sP) < ckt->CKTreltol * MAX(fabs(i_sP_hat),fabs(i_sP)) + ckt->CKTabstol ) &&
                 ( fabs(i_db_hat - i_db) < ckt->CKTreltol * MAX(fabs(i_db_hat),fabs(i_db)) + ckt->CKTabstol ) &&
                 ( fabs(i_sb_hat - i_sb) < ckt->CKTreltol * MAX(fabs(i_sb_hat),fabs(i_sb)) + ckt->CKTabstol )    )
            {
              /* bypass code */
              vds  = *(ckt->CKTstate0 + here->HSMHVvds );
              vgs  = *(ckt->CKTstate0 + here->HSMHVvgs );
              vbs  = *(ckt->CKTstate0 + here->HSMHVvbs );
              vdse = *(ckt->CKTstate0 + here->HSMHVvdse);
              vgse = *(ckt->CKTstate0 + here->HSMHVvgse);
              vbse = *(ckt->CKTstate0 + here->HSMHVvbse);
              vdbd = *(ckt->CKTstate0 + here->HSMHVvdbd);
              vsbs = *(ckt->CKTstate0 + here->HSMHVvsbs);
              vsubs = *(ckt->CKTstate0 + here->HSMHVvsubs);
              deltemp = *(ckt->CKTstate0 + here->HSMHVdeltemp);
              if ( flg_nqs ) {
                Qi_nqs = *(ckt->CKTstate0 + here->HSMHVqi_nqs);
                Qb_nqs = *(ckt->CKTstate0 + here->HSMHVqb_nqs);
              }           

              vges = *(ckt->CKTstate0 + here->HSMHVvges);

              vbd  = vbs - vds;
              vgd  = vgs - vds;
              vgb  = vgs - vbs;
              vged = vges - vds;

              vbs_jct = vsbs;
              vbd_jct = vdbd;

              /* linear branch currents */
              vddp = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVdNode) - *(ckt->CKTrhsOld+here->HSMHVdNodePrime));
              vggp = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVgNode) - *(ckt->CKTrhsOld+here->HSMHVgNodePrime));
              vssp = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVsNode) - *(ckt->CKTrhsOld+here->HSMHVsNodePrime));
              vbpb  = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) - *(ckt->CKTrhsOld+here->HSMHVbNode));
              vbpdb = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) - *(ckt->CKTrhsOld+here->HSMHVdbNode));
              vbpsb = model->HSMHV_type * (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) - *(ckt->CKTrhsOld+here->HSMHVsbNode));

              ByPass = 1;
              goto line755;
            }
          }
        }  /* end of Bypass section */
#endif /*NOBYPASS*/

#ifdef DEBUG_HISIMHVLD_VX        
        printf( "vbd_p    = %12.5e\n" , vbd );
        printf( "vbs_p    = %12.5e\n" , vbs );
        printf( "vgs_p    = %12.5e\n" , vgs );
        printf( "vds_p    = %12.5e\n" , vds );
#endif

        /* start limiting of nonlinear branch voltages */

        von = here->HSMHV_von; 
        Check3 = 0 ;
        if(*(ckt->CKTstate0 + here->HSMHVvds) >= 0.0) { /* case vds>=0 for limiting */
          limval = DEVfetlim(vgs, *(ckt->CKTstate0 + here->HSMHVvgs), von);
          if (vgs != limval) {
            vgs = limval ;
            Check3 = 1 ;
          }
          if (Check3) vds = vgs - vgd;
          limval = DEVlimvds(vds, *(ckt->CKTstate0 + here->HSMHVvds));
          if (vds != limval) {
            vds = limval ;
            Check3 = 2 ;
          }
          vgd = vgs - vds;

          if (here->HSMHV_corg == 1) {
            limval = DEVfetlim(vges, *(ckt->CKTstate0 + here->HSMHVvges), von);
            if (vges != limval) {
              vges = limval ;
              Check3 = 3 ;
            }
            vged = vges - vds;
          }

        } else { /* case vds < 0 for limiting */
          limval = DEVfetlim(vgd, vgdo, von);
          if (vgd != limval) {
            vgd = limval ;
            Check3 = 4 ;
          }
          if (Check3) vds = vgs - vgd;
          limval = -DEVlimvds(-vds, -(*(ckt->CKTstate0 + here->HSMHVvds)));
          if (vds != limval) {
            vds = limval ;
            Check3 = 5 ;
          }
          vgs = vgd + vds;

          if (here->HSMHV_corg == 1) {
            limval = DEVfetlim(vged, vgedo, von);
            if (vged != limval) {
              vged = limval ;
              Check3 = 6 ;
            }
            vges = vged + vds;
          }
        } /* end of case vds< 0 for limiting */

        if (vds >= 0.0) { /* case vds >=0 for limiting of junctions */ 
          vbs = DEVpnjlim(vbs, *(ckt->CKTstate0 + here->HSMHVvbs),
                          CONSTvt0, model->HSMHV_vcrit, &Check1);
          if (Check1) Check3 = 7 ;
          vbd = vbs - vds;
          if (here->HSMHV_corbnet) {
            vsbs = DEVpnjlim(vsbs, *(ckt->CKTstate0 + here->HSMHVvsbs),
                             CONSTvt0, model->HSMHV_vcrit, &Check2);
            if (Check2) Check3 = 8 ;
          }
        } 
        else { /* case vds < 0 for limiting of junctions */
          vbd = DEVpnjlim(vbd, *(ckt->CKTstate0 + here->HSMHVvbd),
                          CONSTvt0, model->HSMHV_vcrit, &Check1);
          if (Check1) Check3 = 9 ;
          vbs = vbd + vds;
          if (here->HSMHV_corbnet) {
            vdbd = DEVpnjlim(vdbd, *(ckt->CKTstate0 + here->HSMHVvdbd),
                             CONSTvt0, model->HSMHV_vcrit, &Check2);
            if (Check2) {
               Check3 = 10 ;
               vdbs = vdbd + vdse;
            }
          }
        }

        if( flg_tempNode > 0 ){
           /* Logarithmic damping of deltemp beyond LIM_TOL */
           deltemp_old = *(ckt->CKTstate0 + here->HSMHVdeltemp);
           if (deltemp > deltemp_old + LIM_TOL)
             {deltemp = deltemp_old + LIM_TOL + log10((deltemp-deltemp_old)/LIM_TOL);
               Check3 = 11;}
           else if (deltemp < deltemp_old - LIM_TOL)
             {deltemp = deltemp_old - LIM_TOL - log10((deltemp_old-deltemp)/LIM_TOL);
               Check3 = 12;}
        }

        /* if (Check3) printf("HSMHV_load: Check3=%d\n",Check3) ; */

        /* limiting completed */
        if (Check3 == 0 ) Check = 0 ;
      } /* loading and limiting of nonlinear branch voltages is completed */

      
      vbd = vbs - vds;
      vgd = vgs - vds;
      vgb = vgs - vbs;
      vged = vges - vds;

      vbs_jct = vsbs;
      vbd_jct = vdbd;

      /* linear branch voltages */
      vddp = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVdNode) -
              *(ckt->CKTrhsOld+here->HSMHVdNodePrime));

      vggp = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVgNode) -
              *(ckt->CKTrhsOld+here->HSMHVgNodePrime));

      vssp = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVsNode) -
              *(ckt->CKTrhsOld+here->HSMHVsNodePrime));

      vbpdb = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) -
              *(ckt->CKTrhsOld+here->HSMHVdbNode));

      vbpb = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) -
              *(ckt->CKTrhsOld+here->HSMHVbNode));


      vbpsb = model->HSMHV_type * 
           (*(ckt->CKTrhsOld+here->HSMHVbNodePrime) -
              *(ckt->CKTrhsOld+here->HSMHVsbNode));


#ifdef DEBUG_HISIMHVLD_VX
      printf( "vbd    = %12.5e\n" , vbd );
      printf( "vbs    = %12.5e\n" , vbs );
      printf( "vgs    = %12.5e\n" , vgs );
      printf( "vds    = %12.5e\n" , vds );
#endif

      /* After loading (and limiting of branch voltages: Start model evaluation */

      /* printf("HSMHV_load: vds=%e vgs=%e vbs=%e vsd=%e vgd=%e vbd=%e\n",
                          vds,vgs,vbs,-vds,vgs-vds,vbs-vds); */

      if (vds >= 0) { /* normal mode */
        here->HSMHV_mode = 1;
        ivds = vds;
        ivgs = vgs;
        ivbs = vbs;
	ivdse = vdse;
	ivgse = vgse;
	ivbse = vbse;
      } else { /* reverse mode */
        here->HSMHV_mode = -1;
        ivds = -vds;
        ivgs = vgd;
        ivbs = vbd;
	ivdse = -vdse;
	ivgse = vgse - vdse;
	ivbse = vbse - vdse;
      }

      if ( model->HSMHV_info >= 5 ) { /* mode, bias conditions ... */
        printf( "--- variables given to HSMHVevaluate() ----\n" );
        printf( "type   = %s\n" , (model->HSMHV_type>0) ? "NMOS" : "PMOS" );
        printf( "mode   = %s\n" , (here->HSMHV_mode>0) ? "NORMAL" : "REVERSE" );
        
        printf( "vbse vbs    = %12.5e %12.5e\n" , vbse, ivbs );
        printf( "vdse vds    = %12.5e %12.5e\n" , vdse, ivds );
        printf( "vgse vgs    = %12.5e %12.5e\n" , vgse, ivgs );
      }
      if ( model->HSMHV_info >= 6 ) { /* input flags */
        printf( "corsrd = %s\n" , (model->HSMHV_corsrd)  ? "true" : "false" ) ;
        printf( "coadov = %s\n" , (model->HSMHV_coadov)  ? "true" : "false" ) ;
        printf( "coisub = %s\n" , (model->HSMHV_coisub)  ? "true" : "false" ) ;
        printf( "coiigs = %s\n" , (model->HSMHV_coiigs)  ? "true" : "false" ) ;
        printf( "cogidl = %s\n" , (model->HSMHV_cogidl)  ? "true" : "false" ) ;
        printf( "coovlp = %s\n" , (model->HSMHV_coovlp)  ? "true" : "false" ) ;
        printf( "coovlps = %s\n" , (model->HSMHV_coovlps) ? "true" : "false" ) ;
        printf( "coflick = %s\n", (model->HSMHV_coflick) ? "true" : "false" ) ;
        printf( "coisti = %s\n" , (model->HSMHV_coisti)  ? "true" : "false" ) ;
        printf( "conqs  = %s\n" , (model->HSMHV_conqs)   ? "true" : "false" ) ;
        printf( "cothrml = %s\n", (model->HSMHV_cothrml) ? "true" : "false" ) ;
        printf( "coign = %s\n"  , (model->HSMHV_coign)   ? "true" : "false" ) ;
        printf( "cosym   = %s\n" , (model->HSMHV_cosym) ? "true" : "false" ) ;
        printf( "coselfheat = %s\n" , (model->HSMHV_coselfheat) ? "true" : "false" ) ;
      }
      /* print inputs ------------AA */

 #ifdef DEBUG_HISIMHVCGG
      /* Print convergence flag */
      printf("isConv %d ", isConv );
      printf("CKTtime %e ", ckt->CKTtime );
      printf("Vb %1.3e ", (model->HSMHV_type>0) ? vbs:-vbs );
      printf("Vd %1.3e ", (model->HSMHV_type>0) ? vds:-vds );
      printf("Vg %1.3e ", (model->HSMHV_type>0) ? vgs:-vgs );
 #endif

      /* call model evaluation */
      if ( HSMHVevaluate(ivdse,ivgse,ivbse,ivds, ivgs, ivbs, vbs_jct, vbd_jct, vsubs, deltemp, here, model, ckt) == HiSIM_ERROR ) 
        return (HiSIM_ERROR);


#ifdef DEBUG_HISIMHVCGG
      printf("HSMHV_ids %e ", here->HSMHV_ids ) ;
      printf("HSMHV_cggb %e ", here->HSMHV_cggb ) ;
      printf("\n") ;
#endif

      here->HSMHV_called += 1;

line755: /* standard entry if HSMHVevaluate is bypassed */
         /* (could be shifted a bit forward ...)       */
      if ( here->HSMHV_mode > 0 ) { /* forward mode */
        Rd         = here->HSMHV_Rd ;
        dRd_dVdse  = here->HSMHV_dRd_dVdse  ;
        dRd_dVgse  = here->HSMHV_dRd_dVgse  ;
        dRd_dVbse  = here->HSMHV_dRd_dVbse  ;
        dRd_dVsubs = (flg_subNode > 0) ? here->HSMHV_dRd_dVsubs : 0.0 ; /* derivative w.r.t. Vsubs */
        dRd_dT     = (flg_tempNode > 0) ? here->HSMHV_dRd_dTi : 0.0  ;
        Rs         = here->HSMHV_Rs ;
        dRs_dVdse  = here->HSMHV_dRs_dVdse  ;
        dRs_dVgse  = here->HSMHV_dRs_dVgse  ;
        dRs_dVbse  = here->HSMHV_dRs_dVbse  ;
        dRs_dVsubs = (flg_subNode > 0) ? here->HSMHV_dRs_dVsubs : 0.0 ; /* derivative w.r.t. Vsubs */
        dRs_dT     = (flg_tempNode > 0) ? here->HSMHV_dRs_dTi : 0.0  ;
        Ids        = here->HSMHV_ids ;
        gds        = here->HSMHV_dIds_dVdsi ;
        gm         = here->HSMHV_dIds_dVgsi ;
        gmbs       = here->HSMHV_dIds_dVbsi ;
        gmT        = (flg_tempNode > 0) ? here->HSMHV_dIds_dTi : 0.0  ;
	gmbs_ext   = here->HSMHV_dIds_dVbse ;
	gds_ext    = here->HSMHV_dIds_dVdse ;
	gm_ext     = here->HSMHV_dIds_dVgse ;

        Qd        = here->HSMHV_qd ;
        dQd_dVds  = here->HSMHV_dQdi_dVdsi ;
        dQd_dVgs  = here->HSMHV_dQdi_dVgsi ;
        dQd_dVbs  = here->HSMHV_dQdi_dVbsi ;
        dQd_dT    = (flg_tempNode > 0) ? here->HSMHV_dQdi_dTi : 0.0  ;
        Qg         = here->HSMHV_qg ;
        dQg_dVds   = here->HSMHV_dQg_dVdsi ;
        dQg_dVgs   = here->HSMHV_dQg_dVgsi ;
        dQg_dVbs   = here->HSMHV_dQg_dVbsi ;
        dQg_dT     = (flg_tempNode > 0) ? here->HSMHV_dQg_dTi : 0.0  ;
        Qs        = here->HSMHV_qs ;
        dQs_dVds  = here->HSMHV_dQsi_dVdsi ;
        dQs_dVgs  = here->HSMHV_dQsi_dVgsi ;
        dQs_dVbs  = here->HSMHV_dQsi_dVbsi ;
        dQs_dT    = (flg_tempNode > 0) ? here->HSMHV_dQsi_dTi : 0.0  ;
        Qb         = - (here->HSMHV_qg + here->HSMHV_qd + here->HSMHV_qs) ;
        dQb_dVds   = here->HSMHV_dQb_dVdsi ;
        dQb_dVgs   = here->HSMHV_dQb_dVgsi ;
        dQb_dVbs   = here->HSMHV_dQb_dVbsi ;
        dQb_dT     = (flg_tempNode > 0) ? here->HSMHV_dQb_dTi : 0.0  ;
        Qfd        = here->HSMHV_qdp ;
        dQfd_dVdse = here->HSMHV_dqdp_dVdse ;
        dQfd_dVgse = here->HSMHV_dqdp_dVgse ;
        dQfd_dVbse = here->HSMHV_dqdp_dVbse ;
        dQfd_dT    = (flg_tempNode > 0) ? here->HSMHV_dqdp_dTi : 0.0  ;
        Qfs        = here->HSMHV_qsp ;
        dQfs_dVdse = here->HSMHV_dqsp_dVdse ;
        dQfs_dVgse = here->HSMHV_dqsp_dVgse ;
        dQfs_dVbse = here->HSMHV_dqsp_dVbse ;
        dQfs_dT    = (flg_tempNode > 0) ? here->HSMHV_dqsp_dTi : 0.0  ;

        Qdext        = here->HSMHV_qdext ;
        dQdext_dVdse = here->HSMHV_dQdext_dVdse ;
        dQdext_dVgse = here->HSMHV_dQdext_dVgse ;
        dQdext_dVbse = here->HSMHV_dQdext_dVbse ;
        dQdext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQdext_dTi : 0.0  ;
        Qgext        = here->HSMHV_qgext ;
        dQgext_dVdse = here->HSMHV_dQgext_dVdse ;
        dQgext_dVgse = here->HSMHV_dQgext_dVgse ;
        dQgext_dVbse = here->HSMHV_dQgext_dVbse ;
        dQgext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQgext_dTi : 0.0  ;
/*        Qsext        = here->HSMHV_qsext ;*/
        dQsext_dVdse = here->HSMHV_dQsext_dVdse ;
        dQsext_dVgse = here->HSMHV_dQsext_dVgse ;
        dQsext_dVbse = here->HSMHV_dQsext_dVbse ;
        dQsext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQsext_dTi : 0.0  ;
        Qbext        = - (here->HSMHV_qgext + here->HSMHV_qdext + here->HSMHV_qsext) ;
        dQbext_dVdse = here->HSMHV_dQbext_dVdse ;
        dQbext_dVgse = here->HSMHV_dQbext_dVgse ;
        dQbext_dVbse = here->HSMHV_dQbext_dVbse ;
        dQbext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQbext_dTi : 0.0  ;
	Isub         = here->HSMHV_isub ;
	dIsub_dVds   = here->HSMHV_dIsub_dVdsi ; 
	dIsub_dVgs   = here->HSMHV_dIsub_dVgsi ; 
	dIsub_dVbs   = here->HSMHV_dIsub_dVbsi ;
        dIsub_dT     = (flg_tempNode > 0) ? here->HSMHV_dIsub_dTi : 0.0  ;
        dIsub_dVdse  = here->HSMHV_dIsub_dVdse ;
	Isubs        = 0.0 ;
	dIsubs_dVds  = 0.0 ; 
	dIsubs_dVgs  = 0.0 ; 
	dIsubs_dVbs  = 0.0 ;
        dIsubs_dT    = 0.0 ;
	dIsubs_dVdse = 0.0 ;
        Igidl        = here->HSMHV_igidl ;
        dIgidl_dVds  = here->HSMHV_dIgidl_dVdsi ;
        dIgidl_dVgs  = here->HSMHV_dIgidl_dVgsi ;
        dIgidl_dVbs  = here->HSMHV_dIgidl_dVbsi ;
        dIgidl_dT    = (flg_tempNode > 0) ? here->HSMHV_dIgidl_dTi : 0.0  ;
        Igisl        = here->HSMHV_igisl ;
        dIgisl_dVds  = here->HSMHV_dIgisl_dVdsi ;
        dIgisl_dVgs  = here->HSMHV_dIgisl_dVgsi ;
        dIgisl_dVbs  = here->HSMHV_dIgisl_dVbsi ;
        dIgisl_dT    = (flg_tempNode > 0) ? here->HSMHV_dIgisl_dTi : 0.0  ;
        Igd          = here->HSMHV_igd ;
        dIgd_dVd   = here->HSMHV_dIgd_dVdsi ;
        dIgd_dVg   = here->HSMHV_dIgd_dVgsi ;
        dIgd_dVb   = here->HSMHV_dIgd_dVbsi ;
        dIgd_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgd_dTi : 0.0  ;
        Igs          = here->HSMHV_igs ;
        dIgs_dVd   = here->HSMHV_dIgs_dVdsi ;
        dIgs_dVg   = here->HSMHV_dIgs_dVgsi ;
        dIgs_dVb   = here->HSMHV_dIgs_dVbsi ;
        dIgs_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgs_dTi : 0.0  ;
        Igb          = here->HSMHV_igb ;
        dIgb_dVd   = here->HSMHV_dIgb_dVdsi ;
        dIgb_dVg   = here->HSMHV_dIgb_dVgsi ;
        dIgb_dVb   = here->HSMHV_dIgb_dVbsi ;
        dIgb_dT      = (flg_tempNode > 0) ? here->HSMHV_dIgb_dTi : 0.0  ;

        /*---------------------------------------------------* 
         * Junction diode.
         *-----------------*/ 
	Ibd = here->HSMHV_ibd ;
	Gbd = here->HSMHV_gbd ;
	Gbdt = (flg_tempNode > 0) ? here->HSMHV_gbdT : 0.0 ;
	
	/* Qbd = here->HSMHV_qbd ; */
	Qbd = *(ckt->CKTstate0 + here->HSMHVqbd) ;
	Cbd = here->HSMHV_capbd ;
	Cbdt = (flg_tempNode > 0) ? here->HSMHV_gcbdT : 0.0 ;

	Ibs = here->HSMHV_ibs ;
	Gbs = here->HSMHV_gbs ;
	Gbst = (flg_tempNode > 0) ? here->HSMHV_gbsT : 0.0 ;

	/* Qbs = here->HSMHV_qbs ; */
	Qbs = *(ckt->CKTstate0 + here->HSMHVqbs) ;
	Cbs = here->HSMHV_capbs ;
	Cbst = (flg_tempNode > 0) ? here->HSMHV_gcbsT : 0.0 ;

        if (flg_nqs) {
          tau         = here->HSMHV_tau       ;
          dtau_dVds   = here->HSMHV_tau_dVdsi ;
          dtau_dVgs   = here->HSMHV_tau_dVgsi ;
          dtau_dVbs   = here->HSMHV_tau_dVbsi ;
          dtau_dT     = here->HSMHV_tau_dTi   ;
          taub        = here->HSMHV_taub      ;
          dtaub_dVds  = here->HSMHV_taub_dVdsi ;
          dtaub_dVgs  = here->HSMHV_taub_dVgsi ;
          dtaub_dVbs  = here->HSMHV_taub_dVbsi ;
          dtaub_dT    = here->HSMHV_taub_dTi   ;
          Qdrat       = here->HSMHV_Xd         ;
          dQdrat_dVds = here->HSMHV_Xd_dVdsi   ;
          dQdrat_dVgs = here->HSMHV_Xd_dVgsi   ;
          dQdrat_dVbs = here->HSMHV_Xd_dVbsi   ;
          dQdrat_dT   = here->HSMHV_Xd_dTi     ;
          Qi          = here->HSMHV_Qi         ;
          dQi_dVds    = here->HSMHV_Qi_dVdsi   ;
          dQi_dVgs    = here->HSMHV_Qi_dVgsi   ;
          dQi_dVbs    = here->HSMHV_Qi_dVbsi   ;
          dQi_dT      = here->HSMHV_Qi_dTi     ;
          Qbulk       = here->HSMHV_Qbulk       ;
          dQbulk_dVds = here->HSMHV_Qbulk_dVdsi ;
          dQbulk_dVgs = here->HSMHV_Qbulk_dVgsi ;
          dQbulk_dVbs = here->HSMHV_Qbulk_dVbsi ;
          dQbulk_dT   = here->HSMHV_Qbulk_dTi   ;
        }

      } else { /* reverse mode */
        /* note: here->HSMHV_Rd and here->HSMHV_Rs are already subjected to mode handling,
           while the following derivatives here->HSMHV_Rd_dVdse, ... are not! */
        Rd        = here->HSMHV_Rd ;
        dRd_dVdse = here->HSMHV_dRd_dVdse ;
        dRd_dVgse = here->HSMHV_dRd_dVgse ;
        dRd_dVbse = here->HSMHV_dRd_dVbse ;
        dRd_dVsubs= (flg_subNode > 0) ? here->HSMHV_dRd_dVsubs : 0.0 ; /* derivative w.r.t. Vsubs */
        dRd_dT    = (flg_tempNode > 0) ? here->HSMHV_dRd_dTi : 0.0  ; 
        Rs        = here->HSMHV_Rs ;
        dRs_dVdse = here->HSMHV_dRs_dVdse ;
        dRs_dVgse = here->HSMHV_dRs_dVgse ;
        dRs_dVbse = here->HSMHV_dRs_dVbse ;
        dRs_dVsubs= (flg_subNode > 0) ? here->HSMHV_dRs_dVsubs : 0.0 ; /* derivative w.r.t. Vsubs */
        dRs_dT    = (flg_tempNode > 0) ? here->HSMHV_dRs_dTi : 0.0  ;
        Ids       = - here->HSMHV_ids ;
        gds       = + (here->HSMHV_dIds_dVdsi + here->HSMHV_dIds_dVgsi + here->HSMHV_dIds_dVbsi) ;
        gm        = - here->HSMHV_dIds_dVgsi ;
        gmbs      = - here->HSMHV_dIds_dVbsi ;
        gmT       = (flg_tempNode > 0) ? - here->HSMHV_dIds_dTi : 0.0  ;
	gds_ext   = + (here->HSMHV_dIds_dVdse + here->HSMHV_dIds_dVgse + here->HSMHV_dIds_dVbse) ;
	gm_ext    = - here->HSMHV_dIds_dVgse;
	gmbs_ext  = - here->HSMHV_dIds_dVbse;

        Qd        =   here->HSMHV_qs ;
        dQd_dVds  = - (here->HSMHV_dQsi_dVdsi + here->HSMHV_dQsi_dVgsi + here->HSMHV_dQsi_dVbsi) ;
        dQd_dVgs  =   here->HSMHV_dQsi_dVgsi ;
        dQd_dVbs  =   here->HSMHV_dQsi_dVbsi ;
        dQd_dT    =   (flg_tempNode > 0) ? here->HSMHV_dQsi_dTi : 0.0  ;
        Qg         =   here->HSMHV_qg ;
        dQg_dVds   = - (here->HSMHV_dQg_dVdsi + here->HSMHV_dQg_dVgsi + here->HSMHV_dQg_dVbsi) ;
        dQg_dVgs   =   here->HSMHV_dQg_dVgsi ;
        dQg_dVbs   =   here->HSMHV_dQg_dVbsi ;
        dQg_dT     =   (flg_tempNode > 0) ? here->HSMHV_dQg_dTi : 0.0  ;
        Qs        =   here->HSMHV_qd ;
        dQs_dVds  = - (here->HSMHV_dQdi_dVdsi + here->HSMHV_dQdi_dVgsi + here->HSMHV_dQdi_dVbsi) ;
        dQs_dVgs  =   here->HSMHV_dQdi_dVgsi ;
        dQs_dVbs  =   here->HSMHV_dQdi_dVbsi ;
        dQs_dT    =   (flg_tempNode > 0) ? here->HSMHV_dQdi_dTi : 0.0  ;
        Qb         = - (here->HSMHV_qg + here->HSMHV_qd + here->HSMHV_qs) ;
        dQb_dVds   = - (here->HSMHV_dQb_dVdsi + here->HSMHV_dQb_dVgsi + here->HSMHV_dQb_dVbsi) ;
        dQb_dVgs   =   here->HSMHV_dQb_dVgsi ;
        dQb_dVbs   =   here->HSMHV_dQb_dVbsi ;
        dQb_dT     =   (flg_tempNode > 0) ? here->HSMHV_dQb_dTi : 0.0  ;
        Qfd        =   here->HSMHV_qsp ;
        dQfd_dVdse = - (here->HSMHV_dqsp_dVdse + here->HSMHV_dqsp_dVgse + here->HSMHV_dqsp_dVbse) ;
        dQfd_dVgse =   here->HSMHV_dqsp_dVgse ;
        dQfd_dVbse =   here->HSMHV_dqsp_dVbse ;
        dQfd_dT    =   (flg_tempNode > 0) ? here->HSMHV_dqsp_dTi : 0.0  ;
        Qfs        =   here->HSMHV_qdp ;
        dQfs_dVdse = - (here->HSMHV_dqdp_dVdse + here->HSMHV_dqdp_dVgse + here->HSMHV_dqdp_dVbse) ;
        dQfs_dVgse =   here->HSMHV_dqdp_dVgse ;
        dQfs_dVbse =   here->HSMHV_dqdp_dVbse ;
        dQfs_dT    =   (flg_tempNode > 0) ? here->HSMHV_dqdp_dTi : 0.0  ;

        Qdext        = here->HSMHV_qsext ;
        dQdext_dVdse = - (here->HSMHV_dQsext_dVdse + here->HSMHV_dQsext_dVgse + here->HSMHV_dQsext_dVbse);
        dQdext_dVgse = here->HSMHV_dQsext_dVgse ;
        dQdext_dVbse = here->HSMHV_dQsext_dVbse ;
        dQdext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQsext_dTi : 0.0  ;
        Qgext        = here->HSMHV_qgext ;
        dQgext_dVdse = - (here->HSMHV_dQgext_dVdse + here->HSMHV_dQgext_dVgse + here->HSMHV_dQgext_dVbse);
        dQgext_dVgse = here->HSMHV_dQgext_dVgse ;
        dQgext_dVbse = here->HSMHV_dQgext_dVbse ;
        dQgext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQgext_dTi : 0.0  ;
/*        Qsext        = here->HSMHV_qdext ;*/
        dQsext_dVdse = - (here->HSMHV_dQdext_dVdse + here->HSMHV_dQdext_dVgse + here->HSMHV_dQdext_dVbse);
        dQsext_dVgse = here->HSMHV_dQdext_dVgse ;
        dQsext_dVbse = here->HSMHV_dQdext_dVbse ;
        dQsext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQdext_dTi : 0.0  ;
        Qbext        = - (here->HSMHV_qgext + here->HSMHV_qdext + here->HSMHV_qsext) ;
        dQbext_dVdse = - (here->HSMHV_dQbext_dVdse + here->HSMHV_dQbext_dVgse + here->HSMHV_dQbext_dVbse);
        dQbext_dVgse = here->HSMHV_dQbext_dVgse ;
        dQbext_dVbse = here->HSMHV_dQbext_dVbse ;
        dQbext_dT    = (here->HSMHV_coselfheat > 0) ? here->HSMHV_dQbext_dTi : 0.0  ;
	Isub         = 0.0 ;
	dIsub_dVds   = 0.0 ; 
	dIsub_dVgs   = 0.0 ; 
	dIsub_dVbs   = 0.0 ;
        dIsub_dT     = 0.0 ;
        dIsub_dVdse  = 0.0 ;
	Isubs        =   here->HSMHV_isub ;
	dIsubs_dVds  = - (here->HSMHV_dIsub_dVdsi + here->HSMHV_dIsub_dVgsi + here->HSMHV_dIsub_dVbsi) ; 
	dIsubs_dVgs  =   here->HSMHV_dIsub_dVgsi ; 
	dIsubs_dVbs  =   here->HSMHV_dIsub_dVbsi ;
        dIsubs_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIsub_dTi : 0.0 ;
	dIsubs_dVdse = - here->HSMHV_dIsub_dVdse ; /* = - (dIsub_dVdse + dIsub_dVbse + dIsub_dVgse) */
        Igidl        =   here->HSMHV_igisl ;
        dIgidl_dVds  = - (here->HSMHV_dIgisl_dVdsi + here->HSMHV_dIgisl_dVgsi + here->HSMHV_dIgisl_dVbsi) ;
        dIgidl_dVgs  =   here->HSMHV_dIgisl_dVgsi ;
        dIgidl_dVbs  =   here->HSMHV_dIgisl_dVbsi ;
        dIgidl_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIgisl_dTi : 0.0  ;
        Igisl        =   here->HSMHV_igidl ;
        dIgisl_dVds  = - (here->HSMHV_dIgidl_dVdsi + here->HSMHV_dIgidl_dVgsi + here->HSMHV_dIgidl_dVbsi) ;
        dIgisl_dVgs  =   here->HSMHV_dIgidl_dVgsi ;
        dIgisl_dVbs  =   here->HSMHV_dIgidl_dVbsi ;
        dIgisl_dT    =   (flg_tempNode > 0) ? here->HSMHV_dIgidl_dTi : 0.0  ;
        /* note: here->HSMHV_igd and here->HSMHV_igs are already subjected to mode handling,
           while the following derivatives here->HSMHV_dIgd_dVdsi, ... are not! */
        Igd          =   here->HSMHV_igd ;
        dIgd_dVd   = - (here->HSMHV_dIgs_dVdsi + here->HSMHV_dIgs_dVgsi + here->HSMHV_dIgs_dVbsi) ;
        dIgd_dVg   =   here->HSMHV_dIgs_dVgsi ;
        dIgd_dVb   =   here->HSMHV_dIgs_dVbsi ;
        dIgd_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgs_dTi : 0.0  ;
        Igs          =   here->HSMHV_igs ;
        dIgs_dVd   = - (here->HSMHV_dIgd_dVdsi + here->HSMHV_dIgd_dVgsi + here->HSMHV_dIgd_dVbsi) ;
        dIgs_dVg   =   here->HSMHV_dIgd_dVgsi ;
        dIgs_dVb   =   here->HSMHV_dIgd_dVbsi ;
        dIgs_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgd_dTi : 0.0  ;
        Igb          =   here->HSMHV_igb ;
        dIgb_dVd   = - (here->HSMHV_dIgb_dVdsi + here->HSMHV_dIgb_dVgsi + here->HSMHV_dIgb_dVbsi) ;
        dIgb_dVg   =   here->HSMHV_dIgb_dVgsi ;
        dIgb_dVb   =   here->HSMHV_dIgb_dVbsi ;
        dIgb_dT      =   (flg_tempNode > 0) ? here->HSMHV_dIgb_dTi : 0.0  ;

	/*---------------------------------------------------* 
	 * Junction diode.
	 *-----------------*/ 
	Ibd = here->HSMHV_ibd ;
	Gbd = here->HSMHV_gbd ;
	Gbdt = (flg_tempNode > 0) ? here->HSMHV_gbdT : 0.0 ;
	
	/* Qbd = here->HSMHV_qbd ; */
	Qbd = *(ckt->CKTstate0 + here->HSMHVqbd) ;
	Cbd = here->HSMHV_capbd ;
	Cbdt = (flg_tempNode > 0) ? here->HSMHV_gcbdT : 0.0 ;

	Ibs = here->HSMHV_ibs ;
	Gbs = here->HSMHV_gbs ;
	Gbst = (flg_tempNode > 0) ? here->HSMHV_gbsT : 0.0 ;

	/* Qbs = here->HSMHV_qbs ; */
	Qbs = *(ckt->CKTstate0 + here->HSMHVqbs) ;
	Cbs = here->HSMHV_capbs ;
	Cbst = (flg_tempNode > 0) ? here->HSMHV_gcbsT : 0.0 ;

        if (flg_nqs) {
          tau         =   here->HSMHV_tau       ;
          dtau_dVds   = -(here->HSMHV_tau_dVdsi + here->HSMHV_tau_dVgsi + here->HSMHV_tau_dVbsi) ;
          dtau_dVgs   =   here->HSMHV_tau_dVgsi ;
          dtau_dVbs   =   here->HSMHV_tau_dVbsi ;
          dtau_dT     =   here->HSMHV_tau_dTi   ;
          taub        =   here->HSMHV_taub      ;
          dtaub_dVds  = -(here->HSMHV_taub_dVdsi + here->HSMHV_taub_dVgsi + here->HSMHV_taub_dVbsi);
          dtaub_dVgs  =   here->HSMHV_taub_dVgsi ;
          dtaub_dVbs  =   here->HSMHV_taub_dVbsi ;
          dtaub_dT    =   here->HSMHV_taub_dTi   ;
          Qdrat       =   1.0 - here->HSMHV_Xd         ;
          dQdrat_dVds = +(here->HSMHV_Xd_dVdsi + here->HSMHV_Xd_dVgsi + here->HSMHV_Xd_dVbsi) ;
          dQdrat_dVgs = - here->HSMHV_Xd_dVgsi   ;
          dQdrat_dVbs = - here->HSMHV_Xd_dVbsi   ;
          dQdrat_dT   = - here->HSMHV_Xd_dTi     ;
          Qi          =   here->HSMHV_Qi         ;
          dQi_dVds    = -(here->HSMHV_Qi_dVdsi + here->HSMHV_Qi_dVgsi + here->HSMHV_Qi_dVbsi) ;
          dQi_dVgs    =   here->HSMHV_Qi_dVgsi   ;
          dQi_dVbs    =   here->HSMHV_Qi_dVbsi   ;
          dQi_dT      =   here->HSMHV_Qi_dTi     ;
          Qbulk       =   here->HSMHV_Qbulk       ;
          dQbulk_dVds = -(here->HSMHV_Qbulk_dVdsi + here->HSMHV_Qbulk_dVgsi + here->HSMHV_Qbulk_dVbsi) ;
          dQbulk_dVgs =   here->HSMHV_Qbulk_dVgsi ;
          dQbulk_dVbs =   here->HSMHV_Qbulk_dVbsi ;
          dQbulk_dT   =   here->HSMHV_Qbulk_dTi   ;
        }
      } /* end of reverse mode */

      if (flg_tempNode > 0) {
        if (pParam->HSMHV_rth > C_RTH_MIN) {
	  Gth = 1.0/pParam->HSMHV_rth ;
        } else {
	  Gth = 1.0/C_RTH_MIN ;
        }
        Ith         = Gth * deltemp ;
        dIth_dT     = Gth ;
        Cth         = pParam->HSMHV_cth ;
        Qth         = Cth * deltemp ;
        /*     P = Ids * (Vdsi + param * ( Vdse - Vdsi)) */
        /*       = Ids * Veffpower                       */
        if ( vds * (vdse - vds) >= 0.0) {
	  if ( pParam->HSMHV_powrat == 1.0 ) {
	    Veffpower        = vdse ;
	    dVeffpower_dVds  = 0.0 ;
	    dVeffpower_dVdse = 1.0 ;
	  } else {
	    Veffpower        = vds + here->HSMHV_powratio * (vdse - vds) ;
	    dVeffpower_dVds  = (1.0 - here->HSMHV_powratio) ;
	    dVeffpower_dVdse = here->HSMHV_powratio ;
	  }




        } else {
           Veffpower        = vds ;
           dVeffpower_dVds  = 1.0 ;
           dVeffpower_dVdse = 0.0 ;
        }
        P           = Ids  * Veffpower ;
        dP_dVds     = gds  * Veffpower + Ids * dVeffpower_dVds;
        dP_dVgs     = gm   * Veffpower ;
        dP_dVbs     = gmbs * Veffpower ;
        dP_dT       = gmT  * Veffpower ;
        dP_dVdse    = gds_ext  * Veffpower + Ids * dVeffpower_dVdse ;
        dP_dVgse    = gm_ext   * Veffpower ;
        dP_dVbse    = gmbs_ext * Veffpower ;

        /* Clamping the maximum rise tempaerarure (SHEMAX) */
        T1 = model->HSMHV_shemax * Gth ;
        Fn_SU( T2 , P , T1 , SHE_MAX_dlt * Gth , T0 ) ;
        P = T2 ;
        dP_dVds     = T0 * dP_dVds ;
        dP_dVgs     = T0 * dP_dVgs ;
        dP_dVbs     = T0 * dP_dVbs ;
        dP_dT       = T0 * dP_dT ;
        dP_dVdse    = T0 * dP_dVdse ; 
        dP_dVgse    = T0 * dP_dVgse ;
        dP_dVbse    = T0 * dP_dVbse ;

      } else {
        Gth         = 0.0 ;
        Ith         = 0.0 ;
        dIth_dT     = 0.0 ;        
        Cth         = 0.0 ;
        Qth         = 0.0 ;
        P           = 0.0 ;
        dP_dVds     = 0.0 ;
        dP_dVgs     = 0.0 ;
        dP_dVbs     = 0.0 ;
        dP_dT       = 0.0 ;
        dP_dVdse    = 0.0 ;
        dP_dVgse    = 0.0 ;
        dP_dVbse    = 0.0 ;
      }

      /* in case of nqs: construct static contributions to the nqs equations (Iqi_nqs, Iqb_nqs)       */
      /*   and nqs charge contributions to inner drain, gate and source node (Qd_nqs, Qg_nqs, Qs_nqs) */
                         
      if (flg_nqs) {
        /* .. tau, taub must be > 0 */
        if (tau < 1.0e-18) {
          tau = 1.e-18 ;
          dtau_dVds = dtau_dVgs = dtau_dVbs = dtau_dT = 0.0 ;
        }
        if (taub < 1.0e-18) {
          taub = 1.0e-18 ;
          dtaub_dVds = dtaub_dVgs = dtaub_dVbs = dtaub_dT = 0.0 ;
        }

        Iqi_nqs          = (Qi_nqs - Qi) / tau ;
        dIqi_nqs_dVds    = - (dQi_dVds + Iqi_nqs * dtau_dVds) / tau ;
        dIqi_nqs_dVgs    = - (dQi_dVgs + Iqi_nqs * dtau_dVgs) / tau ;
        dIqi_nqs_dVbs    = - (dQi_dVbs + Iqi_nqs * dtau_dVbs) / tau ;
        dIqi_nqs_dT      = - (dQi_dT   + Iqi_nqs * dtau_dT  ) / tau ;
        dIqi_nqs_dQi_nqs = 1.0 / tau ;

        Iqb_nqs          = (Qb_nqs - Qbulk) / taub ;
        dIqb_nqs_dVds    = - (dQbulk_dVds + Iqb_nqs * dtaub_dVds) / taub ;
        dIqb_nqs_dVgs    = - (dQbulk_dVgs + Iqb_nqs * dtaub_dVgs) / taub ;
        dIqb_nqs_dVbs    = - (dQbulk_dVbs + Iqb_nqs * dtaub_dVbs) / taub ;
        dIqb_nqs_dT      = - (dQbulk_dT   + Iqb_nqs * dtaub_dT  ) / taub ;
        dIqb_nqs_dQb_nqs = 1.0 / taub ;

        Qd_nqs           = Qi_nqs * Qdrat ;
        dQd_nqs_dVds     = Qi_nqs * dQdrat_dVds ;
        dQd_nqs_dVgs     = Qi_nqs * dQdrat_dVgs ;
        dQd_nqs_dVbs     = Qi_nqs * dQdrat_dVbs ;
        dQd_nqs_dT       = Qi_nqs * dQdrat_dT   ;
        dQd_nqs_dQi_nqs  = Qdrat ;

        Qg_nqs           = - Qi_nqs - Qb_nqs ;
        dQg_nqs_dQi_nqs  = - 1.0 ;
        dQg_nqs_dQb_nqs  = - 1.0 ;

        Qs_nqs           =   Qi_nqs * (1.0 - Qdrat) ;
        dQs_nqs_dVds     = - Qi_nqs * dQdrat_dVds ;
        dQs_nqs_dVgs     = - Qi_nqs * dQdrat_dVgs ;
        dQs_nqs_dVbs     = - Qi_nqs * dQdrat_dVbs ; 
        dQs_nqs_dT       = - Qi_nqs * dQdrat_dT   ;
        dQs_nqs_dQi_nqs  =   1.0 - Qdrat ;
      } else {
        Iqi_nqs = Iqb_nqs = Qd_nqs = Qg_nqs = Qs_nqs = 0.0 ;
      }

      dIgd_dVs = - (dIgd_dVd + dIgd_dVg + dIgd_dVb) ;
      dIgs_dVs = - (dIgs_dVd + dIgs_dVg + dIgs_dVb) ;
      dIgb_dVs = - (dIgb_dVd + dIgb_dVg + dIgb_dVb) ;

     /*---------------------------------------------------* 
      * External Resistances
      *-----------------*/ 
      if(model->HSMHV_corsrd == 1 || model->HSMHV_corsrd == 3 ) {
	if(Rd > 0){
	  GD = 1.0/Rd;
	  GD_dVgs = - dRd_dVgse /Rd/Rd;
	  GD_dVds = - dRd_dVdse /Rd/Rd;
	  GD_dVbs = - dRd_dVbse /Rd/Rd;
	  GD_dVsubs = - dRd_dVsubs /Rd/Rd;
	  GD_dT = - dRd_dT /Rd/Rd;
	}else{
	  GD=0.0;
	  GD_dVgs=0.0;
	  GD_dVds=0.0;
	  GD_dVbs=0.0;
	  GD_dVsubs=0.0;
          GD_dT  =0.0;
	}
	if(Rs > 0){
	  GS = 1.0/Rs;
	  GS_dVgs = - dRs_dVgse /Rs/Rs;
	  GS_dVds = - dRs_dVdse /Rs/Rs;
	  GS_dVbs = - dRs_dVbse /Rs/Rs;
	  GS_dVsubs = - dRs_dVsubs /Rs/Rs;
	  GS_dT = - dRs_dT /Rs/Rs;
	}else{
	  GS=0.0;
	  GS_dVgs=0.0;
	  GS_dVds=0.0;
	  GS_dVbs=0.0;
	  GS_dVsubs=0.0;
          GS_dT  =0.0;
	}
      }
      Iddp        = GD * vddp;
      dIddp_dVddp = GD;
      dIddp_dVdse = GD_dVds * vddp;
      dIddp_dVgse = GD_dVgs * vddp;
      dIddp_dVbse = GD_dVbs * vddp;
      dIddp_dVsubs= GD_dVsubs * vddp;
      dIddp_dT    = GD_dT * vddp;

      Issp        = GS * vssp;
      dIssp_dVssp = GS;
      dIssp_dVdse = GS_dVds * vssp ;
      dIssp_dVgse = GS_dVgs * vssp;
      dIssp_dVbse = GS_dVbs * vssp;
      dIssp_dVsubs= GS_dVsubs * vssp;
      dIssp_dT    = GS_dT * vssp;

      if( model->HSMHV_corg > 0.0 ){
        GG = here->HSMHV_grg ;
      }else{
        GG = 0.0 ;
      }
      Iggp        = GG * vggp;
      dIggp_dVggp = GG;

      if(model->HSMHV_corbnet == 1 && here->HSMHV_rbpb > 0.0 ){
        GRBPB = here->HSMHV_grbpb ;
      }else{
        GRBPB = 0.0 ;
      }
      Ibpb        = GRBPB * vbpb;
      dIbpb_dVbpb = GRBPB;
     
      if(model->HSMHV_corbnet == 1 && here->HSMHV_rbpd > 0.0 ){
        GRBPD = here->HSMHV_grbpd ;
      }else{
        GRBPD = 0.0 ;
      }
      Ibpdb         = GRBPD * vbpdb;
      dIbpdb_dVbpdb = GRBPD;
     
      if(model->HSMHV_corbnet == 1 && here->HSMHV_rbps > 0.0 ){
        GRBPS = here->HSMHV_grbps ;
      }else{
        GRBPS = 0.0 ;
      }
      Ibpsb         = GRBPS * vbpsb;
      dIbpsb_dVbpsb = GRBPS;

      /* printf("HSMHV_load: ByPass=%d\n",ByPass) ; */

      if (!ByPass) { /* no convergence check in case of Bypass */
        /*
         *  check convergence
         */
        isConv = 1;
        if ( (here->HSMHV_off == 0) || !(ckt->CKTmode & MODEINITFIX) ) {
          if (Check == 1) {
            ckt->CKTnoncon++;
            isConv = 0;
#ifndef NEWCONV
          } else {
          /* convergence check for branch currents is done in function HSMHVconvTest */
#endif /* NEWCONV */
          }
        }

        *(ckt->CKTstate0 + here->HSMHVvbs) = vbs;
        *(ckt->CKTstate0 + here->HSMHVvbd) = vbd;
        *(ckt->CKTstate0 + here->HSMHVvgs) = vgs;
        *(ckt->CKTstate0 + here->HSMHVvds) = vds;
        *(ckt->CKTstate0 + here->HSMHVvsbs) = vsbs; 
        *(ckt->CKTstate0 + here->HSMHVvdbs) = vdbs;
        *(ckt->CKTstate0 + here->HSMHVvdbd) = vdbd;
        *(ckt->CKTstate0 + here->HSMHVvges) = vges;
        *(ckt->CKTstate0 + here->HSMHVvsubs) = vsubs;
        *(ckt->CKTstate0 + here->HSMHVdeltemp) = deltemp;
        *(ckt->CKTstate0 + here->HSMHVvdse) = vdse;
        *(ckt->CKTstate0 + here->HSMHVvgse) = vgse;
        *(ckt->CKTstate0 + here->HSMHVvbse) = vbse;
        if ( flg_nqs ) {
          *(ckt->CKTstate0 + here->HSMHVqi_nqs) = Qi_nqs;
          *(ckt->CKTstate0 + here->HSMHVqb_nqs) = Qb_nqs;
        }           
        /* printf("HSMHV_load: (into state0) vds.. = %e %e %e %e %e %e\n",
                                               vds,vgs,vbs,vdse,vgse,vbse); */

        if ((ckt->CKTmode & MODEDC) && 
          !(ckt->CKTmode & MODEINITFIX) && !(ckt->CKTmode & MODEINITJCT)) 
          showPhysVal = 1;
        if (model->HSMHV_show_Given && showPhysVal && isConv) {
          static int isFirst = 1;
          if (vds != vds_pre) 
            ShowPhysVals(here, model, isFirst, vds_pre, vgs, vbs, vgd, vbd, vgb);
          else 
            ShowPhysVals(here, model, isFirst, vds, vgs, vbs, vgd, vbd, vgb);
          if (isFirst) isFirst = 0;
        }
      }

#include "hsmhvld_info_eval.h"

      /* For standard Newton method (SPICE_rhs == 0):                          */
      /*   if currents (and charges) are limited -> extrapolate onto x-values  */
      /* in SPICE mode (SPICE_rhs == 1):                                       */
      /*   extrapolate onto x = 0 (-> rhs_eq)                                  */
      /*                                                                       */
      /* note that                                                             */
      /* the charge extrapolation is replaced by extrapolation of displacement */
      /* currents, see below                                                   */
 
      /* ...... just for easier handling: collect node voltages in vector x:   */
      if (!SPICE_rhs) {
        x[dNode]      = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVdNode));
        x[dNodePrime] = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVdNodePrime));
        x[gNode]      = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVgNode));
        x[gNodePrime] = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVgNodePrime));
        x[sNode]      = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVsNode));
        x[sNodePrime] = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVsNodePrime));
        x[bNodePrime] = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVbNodePrime));
        x[bNode]      = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVbNode));
        x[dbNode]     = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVdbNode));
        x[sbNode]     = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVsbNode));
	if (flg_subNode > 0)
	  x[subNode]  = model->HSMHV_type *( *(ckt->CKTrhsOld+here->HSMHVsubNode)); /* previous vsub */
	else
	  x[subNode]  = 0.0;
        if (flg_tempNode > 0) 
          x[tempNode] =  *(ckt->CKTrhsOld+here->HSMHVtempNode);
        else
          x[tempNode] = 0.0;
        if ( flg_nqs ) {
          x[qiNode] =  *(ckt->CKTrhsOld+here->HSMHVqiNode);
          x[qbNode] =  *(ckt->CKTrhsOld+here->HSMHVqbNode);
        } else {
          x[qiNode] = 0.0;
          x[qbNode] = 0.0;
        }
      }

      delvgs = (x[gNodePrime] - x[sNodePrime]) - vgs;
      delvds = (x[dNodePrime] - x[sNodePrime]) - vds;
      delvbs = (x[bNodePrime] - x[sNodePrime]) - vbs;
      deldeltemp = x[tempNode] - deltemp;

      if (delvgs || delvds || delvbs ||deldeltemp) {
        Ids  += gm         *delvgs + gds        *delvds + gmbs       *delvbs + gmT      *deldeltemp ;
        Isub += dIsub_dVgs *delvgs + dIsub_dVds *delvds + dIsub_dVbs *delvbs + dIsub_dT *deldeltemp ;
        Isubs+= dIsubs_dVgs*delvgs + dIsubs_dVds*delvds + dIsubs_dVbs*delvbs + dIsubs_dT*deldeltemp ;
        Igd  += dIgd_dVg   *delvgs + dIgd_dVd   *delvds + dIgd_dVb   *delvbs + dIgd_dT  *deldeltemp ;
        Igs  += dIgs_dVg   *delvgs + dIgs_dVd   *delvds + dIgs_dVb   *delvbs + dIgs_dT  *deldeltemp ;
        Igb  += dIgb_dVg   *delvgs + dIgb_dVd   *delvds + dIgb_dVb   *delvbs + dIgb_dT  *deldeltemp ;
        Igidl+= dIgidl_dVgs*delvgs + dIgidl_dVds*delvds + dIgidl_dVbs*delvbs + dIgidl_dT*deldeltemp ;
        Igisl+= dIgisl_dVgs*delvgs + dIgisl_dVds*delvds + dIgisl_dVbs*delvbs + dIgisl_dT*deldeltemp ;
        P    += dP_dVgs    *delvgs + dP_dVds    *delvds + dP_dVbs    *delvbs + dP_dT    *deldeltemp ;
        if (flg_nqs) {
          Iqi_nqs += dIqi_nqs_dVgs*delvgs + dIqi_nqs_dVds*delvds + dIqi_nqs_dVbs*delvbs + dIqi_nqs_dT*deldeltemp ;
          Iqb_nqs += dIqb_nqs_dVgs*delvgs + dIqb_nqs_dVds*delvds + dIqb_nqs_dVbs*delvbs + dIqb_nqs_dT*deldeltemp ;
        }
      }

      delvgse = (x[gNodePrime] - x[sNode]) - vgse;
      delvdse = (x[dNode]      - x[sNode]) - vdse;
      delvbse = (x[bNodePrime] - x[sNode]) - vbse;
      if (flg_subNode > 0) delvsubs = (x[subNode] - x[sNode]) - vsubs; /* substrate bias change */

      if (delvgse || delvdse || delvbse ) {
        Ids  += gm_ext     *delvgse + gds_ext    *delvdse + gmbs_ext   *delvbse ;
        Isub += dIsub_dVdse*delvdse ;
        Isubs+= dIsubs_dVdse*delvdse ;
        P    += dP_dVgse   *delvgse + dP_dVdse   *delvdse + dP_dVbse   *delvbse ;
        Iddp += dIddp_dVgse*delvgse + dIddp_dVdse*delvdse + dIddp_dVbse*delvbse ;
        Issp += dIssp_dVgse*delvgse + dIssp_dVdse*delvdse + dIssp_dVbse*delvbse ;
      }

      if (delvsubs) {
        Iddp += dIddp_dVsubs*delvsubs ;
        Issp += dIssp_dVsubs*delvsubs ;
      }

      if (deldeltemp) {
        Iddp += dIddp_dT*deldeltemp ;
        Issp += dIssp_dT*deldeltemp ;
        Ith  += dIth_dT *deldeltemp ;
      }

      delvdbd = (x[dbNode] - x[dNode]) - vbd_jct ;
      if (delvdbd || deldeltemp) {
        Ibd += Gbd*delvdbd + Gbdt*deldeltemp ;
      }

      delvsbs = (x[sbNode] - x[sNode]) - vbs_jct ;
      if (delvsbs || deldeltemp) {
        Ibs += Gbs*delvsbs + Gbst*deldeltemp ;
      }

      delvddp = (x[dNode] - x[dNodePrime]) - vddp ;
      if (delvddp) {
        Iddp += dIddp_dVddp * delvddp ;
      }

      delvssp = (x[sNode] - x[sNodePrime]) - vssp ;
      if (delvssp) {
        Issp += dIssp_dVssp * delvssp ;
      }

      delvggp = (x[gNode] - x[gNodePrime]) - vggp ;
      if (delvggp) {
        Iggp += dIggp_dVggp * delvggp ;
      }

      delvbpb = (x[bNodePrime] - x[bNode]) - vbpb ;
      if (delvbpb) {
        Ibpb += dIbpb_dVbpb * delvbpb ;
      }

      delvbpdb = (x[bNodePrime] - x[dbNode]) - vbpdb ;
      if (delvbpdb) {
        Ibpdb += dIbpdb_dVbpdb * delvbpdb ;
      }

      delvbpsb = (x[bNodePrime] - x[sbNode]) - vbpsb ;
      if (delvbpsb) {
        Ibpsb += dIbpsb_dVbpsb * delvbpsb ;
      }

      if (flg_nqs) {
        delQi_nqs = x[qiNode] - Qi_nqs ;
        if (delQi_nqs) {
          Iqi_nqs += dIqi_nqs_dQi_nqs * delQi_nqs ;
        }
        delQb_nqs = x[qbNode] - Qb_nqs ;
        if (delQb_nqs) {
          Iqb_nqs += dIqb_nqs_dQb_nqs * delQb_nqs ;
        }
      }


      /* Assemble currents into nodes */
      /* ... static part              */
     

      /*  drain node  */
      i_d = Iddp - Ibd ;
      /*  intrinsic drain node */
      i_dP = -Iddp + Ids + Isub + Igidl - Igd ;
      /*  gate node */
      i_g = Iggp ;
      /*  intrinsic gate node */
      i_gP = - Iggp + Igd + Igs + Igb ;
      /*  source node  */
      i_s = Issp - Ibs ;
      /*  intrinsic source node  */
      i_sP = - Issp - Ids + Isubs + Igisl - Igs ;
      /*  intrinsic bulk node */
      i_bP = - Isub - Isubs - Igidl -Igb - Igisl  + Ibpdb + Ibpb + Ibpsb ;
      /*  base node */
      i_b = - Ibpb ;
      /*  drain bulk node  */
      i_db = Ibd - Ibpdb ;
      /*  source bulk node  */
      i_sb = Ibs - Ibpsb ;
      /*  temp node  */
      if (flg_tempNode > 0){
        i_t = Ith - P ;
      } else {
        i_t = 0.0;
      }
      /* nqs equations */
      i_qi = Iqi_nqs ;
      i_qb = Iqb_nqs ;

      for (i = 0; i < XDIM; i++) {
        ydc_d[i] = ydc_dP[i] = ydc_g[i] = ydc_gP[i] = ydc_s[i] = ydc_sP[i] = ydc_bP[i] = ydc_b[i]
	 = ydc_db[i] = ydc_sb[i] = ydc_t[i] = 0.0;
      }
      if (flg_nqs) {
        for (i = 0; i < XDIM; i++) {
           ydc_qi[i] = ydc_qb[i] = 0.0;
        }
      }
     
      /*  drain node  */
      ydc_d[dNode] = dIddp_dVddp + dIddp_dVdse + Gbd ;
      ydc_d[dNodePrime] = -dIddp_dVddp ;
      /* ydc_d[gNode] = 0.0 ; */
      ydc_d[gNodePrime] = dIddp_dVgse ;
      ydc_d[sNode] = - ( dIddp_dVdse + dIddp_dVgse + dIddp_dVbse ) - dIddp_dVsubs ;
      /* ydc_d[sNodePrime] = 0.0 ; */
      ydc_d[bNodePrime] =  dIddp_dVbse ;
      /* ydc_d[bNode] = 0.0 ; */
      ydc_d[dbNode] = - Gbd ;
      /* ydc_d[sbNode] = 0.0 ; */
      ydc_d[subNode] = dIddp_dVsubs ;
      ydc_d[tempNode] = dIddp_dT - Gbdt ;

      /*  intrinsic drain node  */
      ydc_dP[dNode] = - (dIddp_dVddp + dIddp_dVdse) + gds_ext + dIsub_dVdse ;
      ydc_dP[dNodePrime] = dIddp_dVddp + gds + dIsub_dVds + dIgidl_dVds - dIgd_dVd ;
      /* ydc_dP[gNode] = 0.0; */
      ydc_dP[gNodePrime] = -dIddp_dVgse + gm_ext
	 + gm + dIsub_dVgs + dIgidl_dVgs - dIgd_dVg ;
      ydc_dP[sNode] =  dIddp_dVdse + dIddp_dVgse + dIddp_dVbse + dIddp_dVsubs + (-gds_ext -gm_ext -gmbs_ext) - dIsub_dVdse;
      ydc_dP[sNodePrime] = -( gds + dIsub_dVds + dIgidl_dVds ) 
	 - ( gm + dIsub_dVgs + dIgidl_dVgs )
	 - ( gmbs + dIsub_dVbs + dIgidl_dVbs ) - dIgd_dVs ;
      ydc_dP[bNodePrime] = - dIddp_dVbse + gmbs_ext
	 + gmbs + dIsub_dVbs + dIgidl_dVbs - dIgd_dVb;
      /* ydc_dP[bNode] = 0.0; */
      /* ydc_dP[dbNode] = 0.0 ; */
      /* ydc_dP[sbNode] = 0.0 ; */
      ydc_dP[subNode] = - dIddp_dVsubs ;
      ydc_dP[tempNode] = - dIddp_dT 
	 + gmT + dIsub_dT + dIgidl_dT - dIgd_dT ;

      /*  gate node  */
      /* ydc_g[dNode] = 0.0 ; */
      /* ydc_g[dNodePrime] = 0.0 ; */
      ydc_g[gNode] = dIggp_dVggp ;
      ydc_g[gNodePrime] = - dIggp_dVggp ;
      /* ydc_g[sNode] = 0.0 ; */
      /* ydc_g[sNodePrime] = 0.0 ; */
      /* ydc_g[bNodePrime] = 0.0 ; */
      /* ydc_g[bNode] = 0.0 ; */
      /* ydc_g[dbNode] = 0.0 ; */
      /* ydc_g[sbNode] = 0.0 ; */
      /* ydc_g[tempNode] = 0.0 ; */

      /*  intrinsic gate node  */
      /* ydc_gP[dNode] = 0.0 ; */
      ydc_gP[dNodePrime] = dIgd_dVd + dIgs_dVd + dIgb_dVd ;
      ydc_gP[gNode] = - dIggp_dVggp ; 
      ydc_gP[gNodePrime] = dIggp_dVggp + dIgd_dVg + dIgs_dVg + dIgb_dVg ;
      /* ydc_gP[sNode] = 0.0 ; */
      ydc_gP[sNodePrime] = dIgd_dVs + dIgs_dVs + dIgb_dVs ;
      ydc_gP[bNodePrime] = dIgd_dVb + dIgs_dVb + dIgb_dVb ;
      /* ydc_gP[bNode] = 0.0 ; */
      /* ydc_gP[dbNode] = 0.0 ; */ 
      /* ydc_gP[sbNode] = 0.0 ; */
      ydc_gP[tempNode] = dIgd_dT + dIgs_dT + dIgb_dT ;

      /*  source node */
      ydc_s[dNode] = dIssp_dVdse;
      /* ydc_s[dNodePrime] = 0.0 */
      /* ydc_s[gNode] = 0.0 */
      ydc_s[gNodePrime] = dIssp_dVgse;
      ydc_s[sNode] = dIssp_dVssp - ( dIssp_dVgse + dIssp_dVdse + dIssp_dVbse ) - dIssp_dVsubs + Gbs;
      ydc_s[sNodePrime] = - dIssp_dVssp;
      ydc_s[bNodePrime] = dIssp_dVbse ;
      /* ydc_s[bNode] = 0.0 */
      /* ydc_s[dbNode] = 0.0 */
      ydc_s[sbNode]     = - Gbs ;
      ydc_s[subNode] = dIssp_dVsubs;
      ydc_s[tempNode] = dIssp_dT - Gbst;

      /*  intrinsic source node */
      ydc_sP[dNode] = - dIssp_dVdse -gds_ext + dIsubs_dVdse ;
      ydc_sP[dNodePrime] = - gds + dIsubs_dVds + dIgisl_dVds - dIgs_dVd ;
      /* ydc_sP[gNode] = 0.0 ; */
      ydc_sP[gNodePrime] = -dIssp_dVgse -gm_ext
	 - gm + dIsubs_dVgs + dIgisl_dVgs - dIgs_dVg ;
      ydc_sP[sNode] = - dIssp_dVssp - ( - dIssp_dVdse - dIssp_dVgse - dIssp_dVbse ) + dIssp_dVsubs +(gds_ext + gm_ext + gmbs_ext) - dIsubs_dVdse;
      ydc_sP[sNodePrime] = dIssp_dVssp - ( - gds + dIsubs_dVds + dIgisl_dVds )
	 - ( - gm + dIsubs_dVgs + dIgisl_dVgs )
	 - ( - gmbs + dIsubs_dVbs + dIgisl_dVbs ) - dIgs_dVs ;
      ydc_sP[bNodePrime] = -dIssp_dVbse -gmbs_ext
	 - gmbs + dIsubs_dVbs + dIgisl_dVbs - dIgs_dVb ; 
      /* ydc_sP[bNode] = 0.0 ; */
      /* ydc_sP[dbNode] = 0.0 ; */
      /* ydc_sP[sbNode] = 0.0 ; */
      ydc_sP[subNode] = - dIssp_dVsubs;
      ydc_sP[tempNode] = -dIssp_dT 
			  - gmT + dIsubs_dT + dIgisl_dT - dIgs_dT;
     
      /*  intrinsic bulk node */
      ydc_bP[dNode] = - dIsub_dVdse - dIsubs_dVdse ; 
      ydc_bP[dNodePrime] = - dIsub_dVds - dIsubs_dVds - dIgidl_dVds - dIgb_dVd - dIgisl_dVds ;
      /* ydc_bP[gNode] = 0.0 ; */
      ydc_bP[gNodePrime] = - dIsub_dVgs - dIsubs_dVgs - dIgidl_dVgs - dIgb_dVg - dIgisl_dVgs ;
      ydc_bP[sNode] = dIsub_dVdse + dIsubs_dVdse;
      ydc_bP[sNodePrime] = - ( - dIsub_dVds - dIsubs_dVds - dIgidl_dVds - dIgisl_dVds )
       - ( - dIsub_dVgs - dIsubs_dVgs - dIgidl_dVgs - dIgisl_dVgs )
       - ( - dIsub_dVbs - dIsubs_dVbs - dIgidl_dVbs - dIgisl_dVbs ) - dIgb_dVs ; 
      ydc_bP[bNodePrime] = - dIsub_dVbs - dIsubs_dVbs - dIgidl_dVbs - dIgb_dVb - dIgisl_dVbs + dIbpdb_dVbpdb + dIbpb_dVbpb + dIbpsb_dVbpsb ; 
      ydc_bP[bNode] = - dIbpb_dVbpb ; 
      ydc_bP[dbNode] = - dIbpdb_dVbpdb ;
      ydc_bP[sbNode] =  - dIbpsb_dVbpsb ;
      ydc_bP[tempNode] = - dIsub_dT - dIsubs_dT - dIgidl_dT - dIgb_dT - dIgisl_dT ;
     
      /*  bulk node */
      /* ydc_b[dNode] = 0.0 ; */
      /* ydc_b[dNodePrime] = 0.0 ; */
      /* ydc_b[gNode] = 0.0 ; */
      /* ydc_b[gNodePrime] = 0.0 ; */
      /* ydc_b[sNode] = 0.0 ; */
      /* ydc_b[sNodePrime] = 0.0 ; */
      ydc_b[bNodePrime] = - dIbpb_dVbpb ;
      ydc_b[bNode] = dIbpb_dVbpb ;
      /* ydc_b[dbNode] = 0.0 ; */
      /* ydc_b[sbNode] = 0.0 ; */
      /* ydc_b[tempNode] = 0.0 ; */

      /*  drain bulk node */
      ydc_db[dNode] = - Gbd ;
      /* ydc_db[dNodePrime] = 0.0 ; */
      /* ydc_db[gNode] = 0.0 ; */
      /* ydc_db[gNodePrime] = 0.0 ; */
      /* ydc_db[sNode] = 0.0 ; */
      /* ydc_db[sNodePrime] = 0.0 ; */
      ydc_db[bNodePrime] = - dIbpdb_dVbpdb ;
      /* ydc_db[bNode] = 0.0 ; */
      ydc_db[dbNode] = Gbd + dIbpdb_dVbpdb ;
      /* ydc_db[sbNode] = 0.0 ; */
      ydc_db[tempNode] = Gbdt ;
     
      /* source bulk node */
      /* ydc_sb[dNode] = 0.0 ; */
      /* ydc_sb[dNodePrime] = 0.0 ; */
      /* ydc_sb[gNode] = 0.0 ; */
      /* ydc_sb[gNodePrime] = 0.0 ; */
      ydc_sb[sNode] = - Gbs ;
      /* ydc_sb[sNodePrime] = 0.0 ; */
      ydc_sb[bNodePrime] = - dIbpsb_dVbpsb ;
      /* ydc_sb[bNode] = 0.0 ; */
      /* ydc_sb[dbNode] = 0.0 ; */
      ydc_sb[sbNode] = Gbs + dIbpsb_dVbpsb ;
      ydc_sb[tempNode] = Gbst ;
     
      /*  temp node */
      ydc_t[dNode] = - dP_dVdse ;
      ydc_t[dNodePrime] = - dP_dVds ;
      /* ydc_t[gNode] = 0.0 ; */
      ydc_t[gNodePrime] = - dP_dVgs - dP_dVgse ;
      ydc_t[sNode] = - ( - dP_dVdse - dP_dVgse - dP_dVbse ) ;
      ydc_t[sNodePrime] = - ( - dP_dVds - dP_dVgs - dP_dVbs ) ;
      ydc_t[bNodePrime] = - dP_dVbs - dP_dVbse ;
      /* ydc_t[bNode] = 0.0 ; */
      /* ydc_t[dbNode] = 0.0 ; */
      /* ydc_t[sbNode] = 0.0 ; */
      ydc_t[tempNode] = dIth_dT - dP_dT ;

      /* additional entries for flat nqs handling */
      if ( flg_nqs ) {
        ydc_qi[dNodePrime]   =   dIqi_nqs_dVds ;
        ydc_qi[gNodePrime]   =   dIqi_nqs_dVgs ;
        ydc_qi[sNodePrime]   = -(dIqi_nqs_dVds + dIqi_nqs_dVgs + dIqi_nqs_dVbs) ;
        ydc_qi[bNodePrime]   =   dIqi_nqs_dVbs ;  
        ydc_qi[tempNode] =   dIqi_nqs_dT   ;
        ydc_qi[qiNode]   =   dIqi_nqs_dQi_nqs ;    
        /* ydc_qi[qbNode]=   0.0 ;  */
        ydc_qb[dNodePrime]   =   dIqb_nqs_dVds ;
        ydc_qb[gNodePrime]   =   dIqb_nqs_dVgs ;
        ydc_qb[sNodePrime]   = -(dIqb_nqs_dVds + dIqb_nqs_dVgs + dIqb_nqs_dVbs) ;
        ydc_qb[bNodePrime]   =   dIqb_nqs_dVbs ;
        ydc_qb[tempNode] =   dIqb_nqs_dT   ;
        /* ydc_qb[qiNode]=   0.0 ; */
        ydc_qb[qbNode]   =   dIqb_nqs_dQb_nqs ;
      }
     

      /* Preset vectors and matrix for dynamic part */

      cq_d = cq_dP = cq_g = cq_gP = cq_s = cq_sP = cq_bP = cq_b = cq_db = cq_sb = cq_t = cq_qi = cq_qb = 0.0 ;
      for (i = 0; i < XDIM ; i++) {
        ydyn_d[i] = ydyn_dP[i] = ydyn_g[i] = ydyn_gP[i] = ydyn_s[i] = ydyn_sP[i] = ydyn_bP[i] = ydyn_b[i]
	 = ydyn_db[i] = ydyn_sb[i] = ydyn_t[i] = 0.0;
      }
      if (flg_nqs) {
         for (i = 0; i < XDIM ; i++) {
          ydyn_qi[i] = ydyn_qb[i] = 0.0;
        }
      }

      ag0 = ckt->CKTag[0];
     
      if (ChargeComputationNeeded) { /* start handling of dynamic part */

        if (!ByPass) { /* loading of state vector not necessary in case of Bypass */

          /*  intrinsic gate node (without fringing charges) */
          *(ckt->CKTstate0 + here->HSMHVqg) = Qg + Qg_nqs + Qgext;

          /*  intrinsic drain node */
          *(ckt->CKTstate0 + here->HSMHVqd) = Qd + Qd_nqs ; 

          /* intrinsic bulk node */
          *(ckt->CKTstate0 + here->HSMHVqb) = Qb + Qb_nqs + Qbext;

          /*  drain bulk node */
          *(ckt->CKTstate0 + here->HSMHVqbd) = Qbd ;

          /* source bulk node */
          *(ckt->CKTstate0 + here->HSMHVqbs) = Qbs ;

          /* temp node */
          *(ckt->CKTstate0 + here->HSMHVqth) = Qth ;

          /* fringing charges */
          *(ckt->CKTstate0 + here->HSMHVqfd) = Qfd ;
          *(ckt->CKTstate0 + here->HSMHVqfs) = Qfs ;

          /* external drain node */
          *(ckt->CKTstate0 + here->HSMHVqdE) = Qdext;

          /* nqs charges Qi_nqs, Qb_nqs: already loaded above */
          /* if ( flg_nqs ) {                                 */
          /*   *(ckt->CKTstate0 + here->HSMHVqi_nqs) = Qi_nqs; */
          /*   *(ckt->CKTstate0 + here->HSMHVqb_nqs) = Qb_nqs; */
          /* }                                                */
        }

        /* ... assemble capacitance matrix */

        /* ...... drain node */
        ydyn_d[dNode] = dQfd_dVdse + Cbd + dQdext_dVdse ;
        /* ydyn_d[dNodePrime] = 0.0 ; */
        /* ydyn_d[gNode] = 0.0 ; */
        ydyn_d[gNodePrime] = dQfd_dVgse + dQdext_dVgse ;
        ydyn_d[sNode] = - (dQfd_dVdse + dQfd_dVgse+ dQfd_dVbse) - ( dQdext_dVdse + dQdext_dVgse + dQdext_dVbse ) ;
        /* ydyn_d[sNodePrime ] = 0.0 ; */
        ydyn_d[bNodePrime] = dQfd_dVbse + dQdext_dVbse;
        /* ydyn_d[bNode ] = 0.0 ; */
        ydyn_d[dbNode] = - Cbd ;
        /* ydyn_d[sbNode ] = 0.0 ; */
        ydyn_d[tempNode] = dQfd_dT - Cbdt + dQdext_dT ;

        /* ...... intrinsic drain node */
        /* ydyn_dP[dNode] = 0.0 ; */
        ydyn_dP[dNodePrime] = dQd_dVds ;
        /* ydyn_dP[gNode] = 0.0 ; */
        ydyn_dP[gNodePrime] = dQd_dVgs ;
        /* ydyn_dP[sNode] = 0.0 ; */
        ydyn_dP[sNodePrime] = - ( dQd_dVds + dQd_dVgs + dQd_dVbs ) ;
        ydyn_dP[bNodePrime] = dQd_dVbs ;
        /* ydyn_dP[bNode] = 0.0 ; */
        /* ydyn_dP[dbNode] = 0.0 ; */
        /* ydyn_dP[sbNode] = 0.0 ; */
        ydyn_dP[tempNode] = dQd_dT ;

        /* ...... gate node  */
        /*        (no entry) */

        /* ...... intrinsic gate node */
        ydyn_gP[dNode] = -dQfd_dVdse - dQfs_dVdse + dQgext_dVdse ;
        ydyn_gP[dNodePrime] = dQg_dVds ; 
        /* ydyn_gP[gNode] = 0.0 ; */
        ydyn_gP[gNodePrime] = dQg_dVgs -dQfd_dVgse - dQfs_dVgse + dQgext_dVgse ; 
        ydyn_gP[sNode] = dQfd_dVdse + dQfs_dVdse + dQfd_dVgse + dQfs_dVgse + dQfd_dVbse + dQfs_dVbse
                       - ( dQgext_dVdse + dQgext_dVgse + dQgext_dVbse ) ;
        ydyn_gP[sNodePrime] = -( dQg_dVds + dQg_dVgs + dQg_dVbs ) ;
        ydyn_gP[bNodePrime] = dQg_dVbs -dQfd_dVbse - dQfs_dVbse + dQgext_dVbse ;
        /* ydyn_gP[bNode] = 0.0 ; */
        /* ydyn_gP[dbNode] = 0.0 ; */
        /* ydyn_gP[sbNode] = 0.0 ; */
        ydyn_gP[tempNode] = dQg_dT - dQfd_dT - dQfs_dT + dQgext_dT ;

        /* ...... source node */
        ydyn_s[dNode] = dQfs_dVdse + dQsext_dVdse ;
        /* ydyn_d[dNodePrime ] = 0.0 ; */
        /* ydyn_d[gNode ] = 0.0 ; */
        ydyn_s[gNodePrime] = dQfs_dVgse + dQsext_dVgse ;
        ydyn_s[sNode] = Cbs - (dQfs_dVdse + dQfs_dVgse+ dQfs_dVbse) - ( dQsext_dVdse + dQsext_dVgse + dQsext_dVbse ) ;
        /* ydyn_d[sNodePrime ] = 0.0 ; */
        ydyn_s[bNodePrime] = dQfs_dVbse + dQsext_dVbse ;
        /* ydyn_d[bNode ] = 0.0 ; */
        /* ydyn_d[dbNode ] = 0.0 ; */
        ydyn_s[sbNode] = - Cbs ;
        ydyn_s[tempNode] = dQfs_dT - Cbst + dQsext_dT ;

        /* ...... intrinsic source node */
        /* ydyn_sP[dNode] = 0.0 ; */
        ydyn_sP[dNodePrime] = dQs_dVds ;
        /* ydyn_sP[gNode] = 0.0 ; */
        ydyn_sP[gNodePrime] =  dQs_dVgs ; 
        /* ydyn_sP[sNode] = 0.0 ; */
        ydyn_sP[sNodePrime] = - ( dQs_dVds + dQs_dVgs + dQs_dVbs );
        ydyn_sP[bNodePrime] = dQs_dVbs ;
        /* ydyn_sP[bNode] = 0.0 ; */
        /* ydyn_sP[dbNode] = 0.0 ; */
        /* ydyn_sP[sbNode] = 0.0 ; */
        ydyn_sP[tempNode] = dQs_dT ;

        /* ...... intrinsic bulk node */
        ydyn_bP[dNode] = dQbext_dVdse ;
        ydyn_bP[dNodePrime] =  dQb_dVds ;
        /* ydyn_bP[gNode] = 0.0 ; */
        ydyn_bP[gNodePrime] = dQb_dVgs + dQbext_dVgse ; 
        ydyn_bP[sNode] = - ( dQbext_dVdse + dQbext_dVgse + dQbext_dVbse ) ;  
        ydyn_bP[sNodePrime] = - ( dQb_dVds + dQb_dVgs + dQb_dVbs );
        ydyn_bP[bNodePrime] = dQb_dVbs + dQbext_dVbse ;
        /* ydyn_bP[bNode] = 0.0 ; */
        /* ydyn_bP[dbNode] = 0.0 ; */
        /* ydyn_bP[sbNode] = 0.0 ; */
        ydyn_bP[tempNode] = dQb_dT + dQbext_dT ;

        /* ...... bulk node  */
        /*        (no entry) */

        /* ...... drain bulk node */
        ydyn_db[dNode] = - Cbd ;
        /* ydyn_db[dNodePrime] = 0.0 ; */
        /* ydyn_db[gNode] = 0.0 ; */
        /* ydyn_db[gNodePrime] = 0.0 ; */
        /* ydyn_db[sNode] = 0.0 ; */
        /* ydyn_db[sNodePrime] = 0.0 ; */
        /* ydyn_db[bNodePrime] = 0.0 ; */
        /* ydyn_db[bNode] = 0.0 ; */
        ydyn_db[dbNode] = Cbd ;
        /* ydyn_db[sbNode] = 0.0 ; */
        ydyn_db[tempNode] = Cbdt ;

        /* ...... source bulk  node */
        /* ydyn_sb[dNode] = 0.0 ; */
        /* ydyn_sb[dNodePrime] = 0.0 ; */
        /* ydyn_sb[gNode] = 0.0 ; */
        /* ydyn_sb[gNodePrime] = 0.0 ; */
        ydyn_sb[sNode] = - Cbs ;
        /* ydyn_sb[sNodePrime] = 0.0 ; */
        /* ydyn_sb[bNodePrime] = 0.0 ; */
        /* ydyn_sb[bNode] = 0.0 ; */
        /* ydyn_sb[dbNode] = 0.0 ; */
        ydyn_sb[sbNode] = Cbs ;
        ydyn_sb[tempNode] = Cbst ;
     
        /* ...... temp  node */
        /* ydyn_t[dNode] = 0.0 ; */
        /* ydyn_t[dNodePrime] = 0.0 ; */
        /* ydyn_t[gNode] = 0.0 ; */
        /* ydyn_t[gNodePrime] = 0.0 ; */
        /* ydyn_t[sNode] = 0.0 ; */
        /* ydyn_t[sNodePrime] = 0.0 ; */
        /* ydyn_t[bNodePrime] = 0.0 ; */
        /* ydyn_t[bNode] = 0.0 ; */
        /* ydyn_t[dbNode] = 0.0 ; */
        /* ydyn_t[sbNode] = 0.0 ; */
        ydyn_t[tempNode] = Cth ;

        /* additional entries for flat nqs handling */
        if (flg_nqs) {
          /* ...... intrinsic drain node */
          /* ydyn_dP[dNode] += 0.0 ; */
          ydyn_dP[dNodePrime] += dQd_nqs_dVds ;
          /* ydyn_dP[gNode] += 0.0 ; */
          ydyn_dP[gNodePrime] += dQd_nqs_dVgs ;
          /* ydyn_dP[sNode] += 0.0 ; */
          ydyn_dP[sNodePrime] += - ( dQd_nqs_dVds + dQd_nqs_dVgs + dQd_nqs_dVbs ) ;
          ydyn_dP[bNodePrime] += dQd_nqs_dVbs ;
          /* ydyn_dP[bNode] += 0.0 ; */
          /* ydyn_dP[dbNode] += 0.0 ; */
          /* ydyn_dP[sbNode] += 0.0 ; */
          ydyn_dP[tempNode] += dQd_nqs_dT ;
          ydyn_dP[qiNode]    = dQd_nqs_dQi_nqs ;

          /* ...... intrinsic gate node */
          /* ydyn_gP[dNode] += 0.0 ; */
          /* ydyn_gP[dNodePrime] += 0.0 ; */
          /* ydyn_gP[gNode] += 0.0 ; */
          /* ydyn_gP[gNodePrime] += 0.0 ; */
          /* ydyn_gP[sNode] += 0.0 ; */
          /* ydyn_gP[sNodePrime] += 0.0 ; */
          /* ydyn_gP[bNodePrime] += 0.0 ; */
          /* ydyn_gP[bNode] += 0.0 ; */
          /* ydyn_gP[dbNode] += 0.0 ; */
          /* ydyn_gP[sbNode] += 0.0 ; */
          /* ydyn_gP[tempNode] += 0.0 ; */
          ydyn_gP[qiNode] = dQg_nqs_dQi_nqs ;
          ydyn_gP[qbNode] = dQg_nqs_dQb_nqs ;

          /* ...... intrinsic source node */
          /* ydyn_sP[dNode] += 0.0 ; */
          ydyn_sP[dNodePrime] += dQs_nqs_dVds ;
          /* ydyn_sP[gNode] += 0.0 ; */
          ydyn_sP[gNodePrime] +=  dQs_nqs_dVgs ; 
          /* ydyn_sP[sNode] += 0.0 ; */
          ydyn_sP[sNodePrime] += - ( dQs_nqs_dVds + dQs_nqs_dVgs + dQs_nqs_dVbs );
          ydyn_sP[bNodePrime] += dQs_nqs_dVbs ;
          /* ydyn_sP[bNode] += 0.0 ; */
          /* ydyn_sP[dbNode] += 0.0 ; */
          /* ydyn_sP[sbNode] += 0.0 ; */
          ydyn_sP[tempNode] += dQs_nqs_dT ;
          ydyn_sP[qiNode]    = dQs_nqs_dQi_nqs ;

          /* ...... intrinsic bulk node */
          /* ydyn_bP[dNode] += 0.0 ; */
          /* ydyn_bP[dNodePrime] += 0.0 ; */
          /* ydyn_bP[gNode] += 0.0 ; */
          /* ydyn_bP[gNodePrime] += 0.0 ; */
          /* ydyn_bP[sNode] += 0.0 ; */
          /* ydyn_bP[sNodePrime] += 0.0 ; */
          /* ydyn_bP[bNodePrime] += 0.0 ; */
          /* ydyn_bP[bNode] += 0.0 ; */
          /* ydyn_bP[dbNode] += 0.0 ; */
          /* ydyn_bP[sbNode] += 0.0 ; */
          /* ydyn_bP[tempNode] += 0.0 ; */
          ydyn_bP[qbNode] = 1.0 ;

          /* ...... qi node */
          /* ydyn_qi[dNodePrime]   = 0.0 ; */
          /* ydyn_qi[gNodePrime]   = 0.0 ; */
          /* ydyn_qi[sNodePrime]   = 0.0 ; */
          /* ydyn_qi[bNodePrime]   = 0.0 ; */  
          /* ydyn_qi[tempNode] = 0.0 ; */
          ydyn_qi[qiNode]      = 1.0 ;    
          /* ydyn_qi[qbNode]   = 0.0 ; */

          /* ...... qb node */
          /* ydyn_qb[dNodePrime]   = 0.0 ; */
          /* ydyn_qb[gNodePrime]   = 0.0 ; */
          /* ydyn_qb[sNodePrime]   = 0.0 ; */
          /* ydyn_qb[bNodePrime]   = 0.0 ; */
          /* ydyn_qb[tempNode] = 0.0 ; */
          /* ydyn_qb[qiNode]   = 0.0 ; */
          ydyn_qb[qbNode]      = 1.0 ;
        }
     
        if (!ByPass) { /* integrate etc. only necessary if not in Bypass mode! */

          /* store small signal parameters */
          if (ckt->CKTmode & MODEINITSMSIG) {
            /* printf("HSMHV_load: (small signal) ByPass=%d\n",ByPass); */
/*          printf("HSMHV_load: ydc_dP=%e %e %e %e %e %e %e %e\n",
                    ydc_dP[0],ydc_dP[1],ydc_dP[2],ydc_dP[3],ydc_dP[4],ydc_dP[5],ydc_dP[6],ydc_dP[7]);
            printf("HSMHV_load: ych_dP=%e %e %e %e %e %e %e %e\n",
                    ydyn_dP[0],ydyn_dP[1],ydyn_dP[2],ydyn_dP[3],ydyn_dP[4],ydyn_dP[5],ydyn_dP[6],ydyn_dP[7]);
*/
            /* dc matrix into structure 0724*/
            for (i = 0; i < XDIM; i++) {
              here->HSMHV_ydc_d[i] = ydc_d[i]; 
              here->HSMHV_ydc_dP[i] = ydc_dP[i]; 
              here->HSMHV_ydc_g[i] = ydc_g[i]; 
              here->HSMHV_ydc_gP[i] = ydc_gP[i]; 
              here->HSMHV_ydc_s[i] = ydc_s[i]; 
              here->HSMHV_ydc_sP[i] = ydc_sP[i]; 
              here->HSMHV_ydc_bP[i] = ydc_bP[i]; 
              here->HSMHV_ydc_b[i] = ydc_b[i]; 
              here->HSMHV_ydc_db[i] = ydc_db[i]; 
              here->HSMHV_ydc_sb[i] = ydc_sb[i]; 
              here->HSMHV_ydc_t[i] = ydc_t[i]; 
            }
            /* capacitance matrix into structure 0724*/
            for (i = 0; i < XDIM; i++) {
              here->HSMHV_ydyn_d[i] = ydyn_d[i]; 
              here->HSMHV_ydyn_dP[i] = ydyn_dP[i]; 
              here->HSMHV_ydyn_g[i] = ydyn_g[i]; 
              here->HSMHV_ydyn_gP[i] = ydyn_gP[i]; 
              here->HSMHV_ydyn_s[i] = ydyn_s[i]; 
              here->HSMHV_ydyn_sP[i] = ydyn_sP[i]; 
              here->HSMHV_ydyn_bP[i] = ydyn_bP[i]; 
              here->HSMHV_ydyn_b[i] = ydyn_b[i]; 
              here->HSMHV_ydyn_db[i] = ydyn_db[i]; 
              here->HSMHV_ydyn_sb[i] = ydyn_sb[i]; 
              here->HSMHV_ydyn_t[i] = ydyn_t[i]; 
            }
            if (flg_nqs) {
              for (i = 0; i < XDIM; i++) {
                here->HSMHV_ydc_qi[i] = ydc_qi[i]; 
                here->HSMHV_ydc_qb[i] = ydc_qb[i]; 
                here->HSMHV_ydyn_qi[i] = ydyn_qi[i]; 
                here->HSMHV_ydyn_qb[i] = ydyn_qb[i];
              }
            }
            goto line1000; /* that's all for small signal analyses */
          }

          /* Continue handling of dynamic part: */
          /* ... calculate time derivatives of node charges */

          if (ckt->CKTmode & MODEINITTRAN) {
            /* at the very first iteration of the first timepoint:
	       copy charges into previous state -> the integrator may use them ... */
            *(ckt->CKTstate1 + here->HSMHVqb) = *(ckt->CKTstate0 + here->HSMHVqb);
            *(ckt->CKTstate1 + here->HSMHVqg) = *(ckt->CKTstate0 + here->HSMHVqg);
            *(ckt->CKTstate1 + here->HSMHVqd) = *(ckt->CKTstate0 + here->HSMHVqd);
            *(ckt->CKTstate1 + here->HSMHVqth) = *(ckt->CKTstate0 + here->HSMHVqth);
            *(ckt->CKTstate1 + here->HSMHVqbs) = *(ckt->CKTstate0 + here->HSMHVqbs);
            *(ckt->CKTstate1 + here->HSMHVqbd) = *(ckt->CKTstate0 + here->HSMHVqbd);

            *(ckt->CKTstate1 + here->HSMHVqfd) = *(ckt->CKTstate0 + here->HSMHVqfd);
            *(ckt->CKTstate1 + here->HSMHVqfs) = *(ckt->CKTstate0 + here->HSMHVqfs);

            *(ckt->CKTstate1 + here->HSMHVqdE) = *(ckt->CKTstate0 + here->HSMHVqdE);

            if (flg_nqs) {
              *(ckt->CKTstate1 + here->HSMHVqi_nqs) = *(ckt->CKTstate0 + here->HSMHVqi_nqs);
              *(ckt->CKTstate1 + here->HSMHVqb_nqs) = *(ckt->CKTstate0 + here->HSMHVqb_nqs);
            }
          }
     
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqb));
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqg));
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqd));
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqbs));
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqbd));

          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqth));

          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqfd));
          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqfs));

          return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqdE));

          if (flg_nqs) {
            return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqi_nqs));
            return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSMHVqb_nqs));
          }

          if (ckt->CKTmode & MODEINITTRAN) {
            /* at the very first iteration of the first timepoint:
               copy currents into previous state -> the integrator may use them ... */
            *(ckt->CKTstate1 + here->HSMHVcqb) = *(ckt->CKTstate0 + here->HSMHVcqb);
            *(ckt->CKTstate1 + here->HSMHVcqg) = *(ckt->CKTstate0 + here->HSMHVcqg);
            *(ckt->CKTstate1 + here->HSMHVcqd) = *(ckt->CKTstate0 + here->HSMHVcqd);

            *(ckt->CKTstate1 + here->HSMHVcqth) = *(ckt->CKTstate0 + here->HSMHVcqth);

            *(ckt->CKTstate1 + here->HSMHVcqbs) = *(ckt->CKTstate0 + here->HSMHVcqbs);
            *(ckt->CKTstate1 + here->HSMHVcqbd) = *(ckt->CKTstate0 + here->HSMHVcqbd);
            *(ckt->CKTstate1 + here->HSMHVcqfd) = *(ckt->CKTstate0 + here->HSMHVcqfd);
            *(ckt->CKTstate1 + here->HSMHVcqfs) = *(ckt->CKTstate0 + here->HSMHVcqfs);

            *(ckt->CKTstate1 + here->HSMHVcqdE) = *(ckt->CKTstate0 + here->HSMHVcqdE);
            if (flg_nqs) {
              *(ckt->CKTstate1 + here->HSMHVdotqi_nqs) = *(ckt->CKTstate0 + here->HSMHVdotqi_nqs);
              *(ckt->CKTstate1 + here->HSMHVdotqb_nqs) = *(ckt->CKTstate0 + here->HSMHVdotqb_nqs);
            }
          }
        }


        /* ... finally gather displacement currents from data structures */

        cq_dP = *(ckt->CKTstate0 + here->HSMHVcqd);
        cq_gP = *(ckt->CKTstate0 + here->HSMHVcqg);
        cq_bP = *(ckt->CKTstate0 + here->HSMHVcqb);
        cq_sP = - *(ckt->CKTstate0 + here->HSMHVcqg)
                - *(ckt->CKTstate0 + here->HSMHVcqb)
                - *(ckt->CKTstate0 + here->HSMHVcqd);
        cq_dE = *(ckt->CKTstate0 + here->HSMHVcqdE);
        cq_db = *(ckt->CKTstate0 + here->HSMHVcqbd);
        cq_sb = *(ckt->CKTstate0 + here->HSMHVcqbs);
        cq_g = 0.0 ;
        cq_b = 0.0 ;

        /* displacement currents at outer drain/source node (fringing part only!) */
        cq_d = *(ckt->CKTstate0 + here->HSMHVcqfd);
        cq_s = *(ckt->CKTstate0 + here->HSMHVcqfs);

        cq_t = *(ckt->CKTstate0 + here->HSMHVcqth); 

        /* displacement currents due to nqs */
        if (flg_nqs) {
          cq_qi = *(ckt->CKTstate0 + here->HSMHVdotqi_nqs);
          cq_qb = *(ckt->CKTstate0 + here->HSMHVdotqb_nqs);
        } else {
          cq_qi = cq_qb = 0.0 ;
        }

        /* ... and, if necessary: Extrapolate displacement currents onto x-values */

        if (delvdbd || deldeltemp) {
          cq_db += ag0*(Cbd*delvdbd + Cbdt*deldeltemp) ;
        }
        if (delvsbs || deldeltemp) {
          cq_sb += ag0*(Cbs*delvsbs + Cbst*deldeltemp) ;
        }
        if (delvgs || delvds || delvbs || deldeltemp) {
          cq_gP += ag0*(dQg_dVgs*delvgs + dQg_dVds*delvds + dQg_dVbs*delvbs + dQg_dT*deldeltemp) ;
          cq_dP += ag0*(dQd_dVgs*delvgs + dQd_dVds*delvds + dQd_dVbs*delvbs + dQd_dT*deldeltemp) ;
          cq_sP += ag0*(dQs_dVgs*delvgs + dQs_dVds*delvds + dQs_dVbs*delvbs + dQs_dT*deldeltemp) ;
          cq_bP  = - ( cq_gP + cq_dP + cq_sP );
        }
        if (deldeltemp) {
          cq_t += ag0*Cth    *deldeltemp ;
          cq_d += ag0*dQfd_dT*deldeltemp ;
          cq_s += ag0*dQfs_dT*deldeltemp ;

          cq_gE += ag0*dQgext_dT*deldeltemp ;
          cq_dE += ag0*dQdext_dT*deldeltemp ;
          cq_bE += ag0*dQbext_dT*deldeltemp ;
          cq_sE  = - ( cq_gE + cq_dE + cq_bE );
        }
        if (delvgse || delvdse || delvbse) {
          cq_d += ag0*(dQfd_dVgse*delvgse + dQfd_dVdse*delvdse + dQfd_dVbse*delvbse) ;
          cq_s += ag0*(dQfs_dVgse*delvgse + dQfs_dVdse*delvdse + dQfs_dVbse*delvbse) ;

          cq_gE += ag0*(dQgext_dVgse*delvgse + dQgext_dVdse*delvdse + dQgext_dVbse*delvbse) ;
          cq_dE += ag0*(dQdext_dVgse*delvgse + dQdext_dVdse*delvdse + dQdext_dVbse*delvbse) ;
          cq_bE += ag0*(dQbext_dVgse*delvgse + dQbext_dVdse*delvdse + dQbext_dVbse*delvbse) ;
          cq_sE  = - ( cq_gE + cq_dE + cq_bE );
        }

        if (flg_nqs) {
          if (delvgs || delvds || delvbs || deldeltemp) {
            cq_dP += ag0*(dQd_nqs_dVgs*delvgs + dQd_nqs_dVds*delvds + dQd_nqs_dVbs*delvbs + dQd_nqs_dT*deldeltemp) ;
            cq_sP += ag0*(dQs_nqs_dVgs*delvgs + dQs_nqs_dVds*delvds + dQs_nqs_dVbs*delvbs + dQs_nqs_dT*deldeltemp) ;
            cq_bP  = - ( cq_gP + cq_dP + cq_sP ); /* should be superfluous ? */
          }
          if (delQi_nqs) {
            cq_dP += ag0*dQd_nqs_dQi_nqs*delQi_nqs ;
            cq_gP += ag0*dQg_nqs_dQi_nqs*delQi_nqs ;
            cq_sP += ag0*dQs_nqs_dQi_nqs*delQi_nqs ;
            cq_qi += ag0*     1.0       *delQi_nqs ;
          }
          if (delQb_nqs) {
            cq_gP += ag0*dQg_nqs_dQb_nqs*delQb_nqs ;
            cq_bP += ag0*     1.0       *delQb_nqs ;
            cq_qb += ag0*     1.0       *delQb_nqs ;
          }
        }
      } /* End of handling dynamic part */


      /* Assemble total node currents (type-handling shifted to stamping) */

      cur_d  = i_d  + cq_d - cq_db + cq_dE ;
      cur_dP = i_dP + cq_dP ;
      cur_g  = i_g  + cq_g  ;
      cur_gP = i_gP + cq_gP - cq_d - cq_s + cq_gE ;
      cur_s  = i_s  + cq_s - cq_sb + cq_sE ;
      cur_sP = i_sP + cq_sP;
      cur_bP = i_bP + cq_bP + cq_bE ;
      cur_b  = i_b  + cq_b;
      cur_db = i_db + cq_db;
      cur_sb = i_sb + cq_sb;
      cur_t  = i_t  + cq_t;
      cur_qi = i_qi + cq_qi;
      cur_qb = i_qb + cq_qb;


      /* Now we can start stamping ...                                    */
      /* ... right hand side: subtract total node currents                */

      *(ckt->CKTrhs + here->HSMHVdNode)      -= model->HSMHV_type * cur_d;
      *(ckt->CKTrhs + here->HSMHVdNodePrime) -= model->HSMHV_type * cur_dP;
      *(ckt->CKTrhs + here->HSMHVgNode)      -= model->HSMHV_type * cur_g;
      *(ckt->CKTrhs + here->HSMHVgNodePrime) -= model->HSMHV_type * cur_gP;
      *(ckt->CKTrhs + here->HSMHVsNode)      -= model->HSMHV_type * cur_s;
      *(ckt->CKTrhs + here->HSMHVsNodePrime) -= model->HSMHV_type * cur_sP;
      *(ckt->CKTrhs + here->HSMHVbNodePrime) -= model->HSMHV_type * cur_bP;
      *(ckt->CKTrhs + here->HSMHVbNode)      -= model->HSMHV_type * cur_b;
      *(ckt->CKTrhs + here->HSMHVdbNode)     -= model->HSMHV_type * cur_db;
      *(ckt->CKTrhs + here->HSMHVsbNode)     -= model->HSMHV_type * cur_sb;
      if( flg_tempNode > 0) { 
        *(ckt->CKTrhs + here->HSMHVtempNode) -= cur_t;  /* temp node independent of model type! */
      }
      if (flg_nqs) {
        *(ckt->CKTrhs + here->HSMHVqiNode) -= cur_qi;
        *(ckt->CKTrhs + here->HSMHVqbNode) -= cur_qb;
      }


      /* ... finally stamp matrix */

      /*drain*/
      *(here->HSMHVDdPtr)  += ydc_d[dNode]         + ag0*ydyn_d[dNode];
      *(here->HSMHVDdpPtr) += ydc_d[dNodePrime]    + ag0*ydyn_d[dNodePrime];
      *(here->HSMHVDgpPtr) += ydc_d[gNodePrime]    + ag0*ydyn_d[gNodePrime];
      *(here->HSMHVDsPtr)  += ydc_d[sNode]         + ag0*ydyn_d[sNode];
      *(here->HSMHVDspPtr) += ydc_d[sNodePrime]    + ag0*ydyn_d[sNodePrime];
      *(here->HSMHVDbpPtr) += ydc_d[bNodePrime]    + ag0*ydyn_d[bNodePrime];
      *(here->HSMHVDdbPtr) += ydc_d[dbNode]        + ag0*ydyn_d[dbNode];
      if (flg_subNode > 0) {
	*(here->HSMHVDsubPtr) += ydc_d[subNode];
      }
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVDtempPtr) += model->HSMHV_type * (ydc_d[tempNode] + ag0*ydyn_d[tempNode]);
      }

      /*drain prime*/
      *(here->HSMHVDPdPtr)  +=  ydc_dP[dNode]      + ag0*ydyn_dP[dNode];
      *(here->HSMHVDPdpPtr) +=  ydc_dP[dNodePrime] + ag0*ydyn_dP[dNodePrime];
      *(here->HSMHVDPgpPtr) +=  ydc_dP[gNodePrime] + ag0*ydyn_dP[gNodePrime];
      *(here->HSMHVDPsPtr)  +=  ydc_dP[sNode]      + ag0*ydyn_dP[sNode];
      *(here->HSMHVDPspPtr) +=  ydc_dP[sNodePrime] + ag0*ydyn_dP[sNodePrime];
      *(here->HSMHVDPbpPtr) +=  ydc_dP[bNodePrime] + ag0*ydyn_dP[bNodePrime];
      if (flg_subNode > 0) {
	*(here->HSMHVDPsubPtr) += ydc_dP[subNode];
      }
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVDPtempPtr) +=  model->HSMHV_type * (ydc_dP[tempNode] + ag0*ydyn_dP[tempNode]);
      }
      if (flg_nqs) {
        *(here->HSMHVDPqiPtr) +=  model->HSMHV_type * (ydc_dP[qiNode] + ag0*ydyn_dP[qiNode]);
      }

      /*gate*/     
      *(here->HSMHVGgPtr)   +=  ydc_g[gNode]       + ag0*ydyn_g[gNode];
      *(here->HSMHVGgpPtr)  +=  ydc_g[gNodePrime]  + ag0*ydyn_g[gNodePrime];

      /*gate prime*/
      *(here->HSMHVGPdPtr)  +=  ydc_gP[dNode]      + ag0*ydyn_gP[dNode];
      *(here->HSMHVGPdpPtr) +=  ydc_gP[dNodePrime] + ag0*ydyn_gP[dNodePrime];
      *(here->HSMHVGPgPtr)  +=  ydc_gP[gNode]      + ag0*ydyn_gP[gNode];
      *(here->HSMHVGPgpPtr) +=  ydc_gP[gNodePrime] + ag0*ydyn_gP[gNodePrime];
      *(here->HSMHVGPsPtr)  +=  ydc_gP[sNode]      + ag0*ydyn_gP[sNode];
      *(here->HSMHVGPspPtr) +=  ydc_gP[sNodePrime] + ag0*ydyn_gP[sNodePrime];
      *(here->HSMHVGPbpPtr) +=  ydc_gP[bNodePrime] + ag0*ydyn_gP[bNodePrime];
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVGPtempPtr) +=  model->HSMHV_type * (ydc_gP[tempNode] + ag0*ydyn_gP[tempNode]);
      }
      if (flg_nqs) {
	*(here->HSMHVGPqiPtr) +=  model->HSMHV_type * (ydc_gP[qiNode] + ag0*ydyn_gP[qiNode]);
	*(here->HSMHVGPqbPtr) +=  model->HSMHV_type * (ydc_gP[qbNode] + ag0*ydyn_gP[qbNode]);
      }

      /*source*/
      *(here->HSMHVSdPtr)  += ydc_s[dNode]         + ag0*ydyn_s[dNode];
      *(here->HSMHVSsPtr)  += ydc_s[sNode]         + ag0*ydyn_s[sNode];
      *(here->HSMHVSdpPtr) += ydc_s[dNodePrime]    + ag0*ydyn_s[dNodePrime];
      *(here->HSMHVSgpPtr) += ydc_s[gNodePrime]    + ag0*ydyn_s[gNodePrime];
      *(here->HSMHVSspPtr) += ydc_s[sNodePrime]    + ag0*ydyn_s[sNodePrime];
      *(here->HSMHVSbpPtr) += ydc_s[bNodePrime]    + ag0*ydyn_s[bNodePrime];
      *(here->HSMHVSsbPtr) += ydc_s[sbNode]        + ag0*ydyn_s[sbNode];
      if (flg_subNode > 0) {
	*(here->HSMHVSsubPtr) += ydc_s[subNode];
      }
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVStempPtr) += model->HSMHV_type * (ydc_s[tempNode]+ ag0*ydyn_s[tempNode]);
      }

      /*source prime*/
      *(here->HSMHVSPdPtr)  +=  ydc_sP[dNode]      + ag0*ydyn_sP[dNode];
      *(here->HSMHVSPdpPtr) +=  ydc_sP[dNodePrime] + ag0*ydyn_sP[dNodePrime];
      *(here->HSMHVSPgpPtr) +=  ydc_sP[gNodePrime] + ag0*ydyn_sP[gNodePrime];
      *(here->HSMHVSPsPtr)  +=  ydc_sP[sNode]      + ag0*ydyn_sP[sNode];
      *(here->HSMHVSPspPtr) +=  ydc_sP[sNodePrime] + ag0*ydyn_sP[sNodePrime];
      *(here->HSMHVSPbpPtr) +=  ydc_sP[bNodePrime] + ag0*ydyn_sP[bNodePrime];
      if (flg_subNode > 0) {
	*(here->HSMHVSPsubPtr) += ydc_sP[subNode];
      }
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVSPtempPtr) +=  model->HSMHV_type * (ydc_sP[tempNode] + ag0*ydyn_sP[tempNode]);
      }
      if (flg_nqs) {
	*(here->HSMHVSPqiPtr) +=  model->HSMHV_type * (ydc_sP[qiNode] + ag0*ydyn_sP[qiNode]);
      }

      /*bulk prime*/
      *(here->HSMHVBPdPtr)  +=  ydc_bP[dNode]      + ag0*ydyn_bP[dNode]; 
      *(here->HSMHVBPdpPtr) +=  ydc_bP[dNodePrime] + ag0*ydyn_bP[dNodePrime];
      *(here->HSMHVBPgpPtr) +=  ydc_bP[gNodePrime] + ag0*ydyn_bP[gNodePrime];
      *(here->HSMHVBPspPtr) +=  ydc_bP[sNodePrime] + ag0*ydyn_bP[sNodePrime];
      *(here->HSMHVBPsPtr)  +=  ydc_bP[sNode]      + ag0*ydyn_bP[sNode];
      *(here->HSMHVBPbpPtr) +=  ydc_bP[bNodePrime] + ag0*ydyn_bP[bNodePrime];
      *(here->HSMHVBPbPtr)  +=  ydc_bP[bNode]      + ag0*ydyn_bP[bNode];
      *(here->HSMHVBPdbPtr) +=  ydc_bP[dbNode]     + ag0*ydyn_bP[dbNode];
      *(here->HSMHVBPsbPtr) +=  ydc_bP[sbNode]     + ag0*ydyn_bP[sbNode];
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVBPtempPtr) +=  model->HSMHV_type * (ydc_bP[tempNode] + ag0*ydyn_bP[tempNode]);
      }
      if (flg_nqs) {
	*(here->HSMHVBPqbPtr) +=  model->HSMHV_type * (ydc_bP[qbNode] + ag0*ydyn_bP[qbNode]);
      }

      /*bulk*/
      *(here->HSMHVBbpPtr) +=  ydc_b[bNodePrime]   + ag0*ydyn_b[bNodePrime];
      *(here->HSMHVBbPtr)  +=  ydc_b[bNode]        + ag0*ydyn_b[bNode];

      /*drain bulk*/
      *(here->HSMHVDBdPtr)  +=  ydc_db[dNode]      + ag0*ydyn_db[dNode];
      *(here->HSMHVDBbpPtr) +=  ydc_db[bNodePrime] + ag0*ydyn_db[bNodePrime];
      *(here->HSMHVDBdbPtr) +=  ydc_db[dbNode]     + ag0*ydyn_db[dbNode];
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVDBtempPtr) +=  model->HSMHV_type * (ydc_db[tempNode] + ag0*ydyn_db[tempNode]);
      }

      /*source bulk*/
      *(here->HSMHVSBsPtr)  +=  ydc_sb[sNode]      + ag0*ydyn_sb[sNode];
      *(here->HSMHVSBbpPtr) +=  ydc_sb[bNodePrime] + ag0*ydyn_sb[bNodePrime];
      *(here->HSMHVSBsbPtr) +=  ydc_sb[sbNode]     + ag0*ydyn_sb[sbNode];
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVSBtempPtr) +=  model->HSMHV_type * (ydc_sb[tempNode] + ag0*ydyn_sb[tempNode]);
      }

      /*temp*/
      if( flg_tempNode > 0) { 
        /* temp entries in matrix dependent on model type */
        *(here->HSMHVTempdPtr)  +=  model->HSMHV_type * (ydc_t[dNode]      + ag0*ydyn_t[dNode]     );
        *(here->HSMHVTempdpPtr) +=  model->HSMHV_type * (ydc_t[dNodePrime] + ag0*ydyn_t[dNodePrime]);
        *(here->HSMHVTempgpPtr) +=  model->HSMHV_type * (ydc_t[gNodePrime] + ag0*ydyn_t[gNodePrime]);
        *(here->HSMHVTempsPtr)  +=  model->HSMHV_type * (ydc_t[sNode]      + ag0*ydyn_t[sNode]     );
        *(here->HSMHVTempspPtr) +=  model->HSMHV_type * (ydc_t[sNodePrime] + ag0*ydyn_t[sNodePrime]);
        *(here->HSMHVTempbpPtr) +=  model->HSMHV_type * (ydc_t[bNodePrime] + ag0*ydyn_t[bNodePrime]);
        /* no type factor at main diagonal temp entry! */
        *(here->HSMHVTemptempPtr) +=  ydc_t[tempNode] + ag0*ydyn_t[tempNode];
      }
     
      /* additional entries for flat nqs handling */
      if ( flg_nqs ) {
        /*qi*/
	*(here->HSMHVQIdpPtr) +=  model->HSMHV_type * (ydc_qi[dNodePrime] + ag0*ydyn_qi[dNodePrime]);
	*(here->HSMHVQIgpPtr) +=  model->HSMHV_type * (ydc_qi[gNodePrime] + ag0*ydyn_qi[gNodePrime]);
	*(here->HSMHVQIspPtr) +=  model->HSMHV_type * (ydc_qi[sNodePrime] + ag0*ydyn_qi[sNodePrime]);
	*(here->HSMHVQIbpPtr) +=  model->HSMHV_type * (ydc_qi[bNodePrime] + ag0*ydyn_qi[bNodePrime]);
	*(here->HSMHVQIqiPtr) +=                     (ydc_qi[qiNode] + ag0*ydyn_qi[qiNode]);
        if ( flg_tempNode > 0 ) { /* self heating */
	  *(here->HSMHVQItempPtr) +=                     (ydc_qi[tempNode] + ag0*ydyn_qi[tempNode]);
        }

        /*qb*/
	*(here->HSMHVQBdpPtr) +=  model->HSMHV_type * (ydc_qb[dNodePrime] + ag0*ydyn_qb[dNodePrime]);
	*(here->HSMHVQBgpPtr) +=  model->HSMHV_type * (ydc_qb[gNodePrime] + ag0*ydyn_qb[gNodePrime]);
	*(here->HSMHVQBspPtr) +=  model->HSMHV_type * (ydc_qb[sNodePrime] + ag0*ydyn_qb[sNodePrime]);
	*(here->HSMHVQBbpPtr) +=  model->HSMHV_type * (ydc_qb[bNodePrime] + ag0*ydyn_qb[bNodePrime]);
	*(here->HSMHVQBqbPtr) +=                     (ydc_qb[qbNode] + ag0*ydyn_qb[qbNode]);
        if ( flg_tempNode > 0 ) { /* self heating */
	  *(here->HSMHVQBtempPtr) +=                     (ydc_qb[tempNode] + ag0*ydyn_qb[tempNode]);
        }
      }


line1000:

      if (ckt->CKTnoncon != noncon_old) {
        ckt->CKTtroubleElt = (GENinstance *) here;
      }
     
      
    } /* End of MOSFET Instance */
  } /* End of Model Instance */

  return(OK);
}
