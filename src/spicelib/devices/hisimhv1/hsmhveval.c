/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 1  SUBVERSION : 2  REVISION : 4 )
 Model Parameter VERSION : 1.23
 FILE : hsmhveval.c

 DATE : 2013.04.30

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

/**********************************************************************

The following source code, and all copyrights, trade secrets or other
intellectual property rights in and to the source code in its entirety,
is owned by the Hiroshima University and the STARC organization.

All users need to follow the "HISIM_HV Distribution Statement and
Copyright Notice" attached to HiSIM_HV model.

-----HISIM_HV Distribution Statement and Copyright Notice--------------
 Hiroshima University and/or Semiconductor Technology Academic Research 
Center ("STARC") grants to Licensee a worldwide, royalty-free, 
sub-licensable, retroactive, perpetual, irrevocable license to make, 
have made, offer to sell, sell, import, export, use, copy, redistribute, 
perform, display and incorporate HiSIM_HV intellectual property (the 
"License") subject to the conditions set forth below. 
 This License includes rights to use and modify copyrighted material 
for any purpose if the copyright is acknowledged as well as the use of 
related patents which STARC owns, has applied for or will apply for 
in connection with HiSIM_HV intellectual property for the purpose of 
implementing and using he HiSIM_HV intellectual property in connection
with the standard. This license applies to all past and future versions 
of HiSIM_HV.

1. HiSIM_HV intellectual property is offered "as is" without any warranty, 
explicit or implied, or service support. Hiroshima University, STARC, 
its University staff and employees assume no liability for the quality 
and performance of HiSIM_HV intellectual property.

2. As the owner of the HiSIM_HV intellectual property, and all other 
related rights, Hiroshima University and/or STARC grant the License 
as set forth above.

3. A Licensee may not charge an end user a fee for the HiSIM_HV source 
code, which Hiroshima University and STARC own, by itself, however, 
a Licensee may charge an end user a fee for alterations or additions 
to the HiSIM_HV source code or for maintenance service.

4. A Licensee of HiSIM_HV intellectual property agrees that Hiroshima 
University and STARC are the developers of HiSIM_HV in all products 
containing HiSIM_HV and the alteration thereof (subject to Licensee's 
ownership of the alterations). 
If future versions of HiSIM_HV incorporate elements of other CMC models 
the copyrights of those elements remains with the original developers. 
For this purpose the copyright notice as shown below shall be used.

"The HiSIM_HV source code, and all copyrights, trade secrets or other 
intellectual property rights in and to the source code, is owned 
by Hiroshima University and/or STARC."

5. A Licensee of HiSIM_HV intellectual property will comply with the 
export obligations pertaining to the export of the HiSIM_HV intellectual 
property.

6.  By using HiSIM_HV intellectual property owned by Hiroshima University
and/or STARC, Licensee agrees not to prosecute any patents or patent  
held by Licensee that are required for implementation of HiSIM_HV against
any party who is infringing those patents solely by implementing and/or
using the HiSIM_HV standard.


Toshimasa Asahara, President, Hiroshima University
Mitiko Miura-Mattausch, Professor, Hiroshima University
Katsuhiro Shimohigashi, President&CEO, STARC
June. 2008  (revised in June 2011)
*************************************************************************/

/*********************************************************************
* Memorandum on programming
* 
* (1) Bias (x: b|d|g)
*     . vxs : Input arguments: Outer branch voltages.
*     . vxsi: Input arguments: Inner branch voltages.
*     . deltemp: Input argument: delta temperature.
*     . Vxse: Internal name for outer branch voltages.
*     . Vxs:  Internal name for inner branch voltages.
*     . Vbscl:Inner bulk source voltage is clamped within a specified region. 
*     . Y_dVxs denotes the partial derivative of Y w.r.t. Vxs.
*     . Y_dVxse denotes the partial derivative of Y w.r.t. Vxse.
*     . Y_dT denotes derivatives with respect to deltemp.
* 
* (2) Device Mode
*     . Normal mode (Vds>=0 for nMOS) is assumed.
*     . The sign of Vdse is assumed to be changed simultaneously, if the sign of Vds is changed;
*       hence Vdse may become negative even thogh Vds >=0.
*     . In case of reverse mode, parent routines have to properly 
*       transform or interchange inputs and outputs except ones 
*       related to junction diodes, which are regarded as being 
*       fixed to the nodal S/D.
*
* (3) Modification for symmetry at Vds=0
*     . Vxsz: Modified bias. (x: b|d|g)
*     . Ps0z: Modified Ps0.
*     . The following variables are calculated as a function of 
*       modified biases or potential.
*         Tox, Cox, (-- with quantum effect)
*         Vth*, dVth*, dPpg, Igate, Igidl, Igisl. 
*     . The following variables are calculated using a transform
*       function.
*         Lred
*     
* (4) Zones and Cases (terminology)
* 
*       Chi:=beta*(Ps0-Vbs)=       0    3    5
*
*                      Zone:    A  | D1 | D2 | D3
*                                  |
*                    (accumulation)|(depletion)
*                                  |
*                      Vgs =     Vgs_fb                Vth
*                                              /       /
*                      Case:    Nonconductive / Conductive
*                                            /
*             VgVt:=Qn0/Cox=             VgVt_small
*
*     . Ids is regarded as zero in zone-A.
*     . Procedure to calculate Psl and dependent variables is 
*       omitted in the nonconductive case. Ids and Qi are regarded
*       as zero in this case.
*
*********************************************************************/

/*===========================================================*
* Preamble.
*=================*/
/*---------------------------------------------------*
* Header files.
*-----------------*/
#include "ngspice/ngspice.h"
#ifdef __STDC__
/* #include <ieeefp.h> */
#endif


/*-----------------------------------*
* HiSIM macros
*-----------------*/
#include "hisimhv.h"
#include "hsmhvevalenv.h"
#define C_IDD_MIN    1.0e-15
#define C_sub_delta 0.1 /* CHECK! */
#define C_sub_delta2 1.0e-9 /* CHECK! */
#define C_gidl_delta 0.5

/* local variables used in macro functions */
static double TMF0 , TMF1 , TMF2 , TMF3 , TMF4 , TMF5 , TMF6 ; 
/*===========================================================*
* pow
*=================*/
#ifdef POW_TO_EXP_AND_LOG
#define Fn_Pow( x , y )  exp( y * log( x )  ) 
#else
#define Fn_Pow( x , y )  pow( x , y )
#endif

/*===========================================================*
* Exp() for PGD.
* - ExpLim(-3)=0
*=================*/

#define Fn_ExpLim( y , x , dx ) { \
    if ( (x) < -3.0 ) { \
      dx = 0.0 ; \
      y = 0.0 ; \
    } else if ( (x) < 0.0 ) { \
      dx =  1.0 + (x) * ( 2 * (1.0/3.0) + (x) * 3 * (1.0/27.0) )  ; \
      y = 1.0 + (x) * ( 1.0 + (x) * ( (1.0/3.0) + (x) * (1.0/27.0) ) ) ; \
    } else { \
      dx =  1.0 + (x) * ( 2 * (1.0/3.0) + (x) * ( 3 * 0.0402052934513951 \
                + (x) * 4 * 0.148148111111111 ) ) ; \
      y = 1.0 + (x) * ( 1.0 + (x) * ( (1.0/3.0) + (x) * ( 0.0402052934513951 \
              + (x) * 0.148148111111111 ) ) ) ; \
    } \
}

/*===========================================================*
* Macro Functions for ceiling/flooring/symmetrization.
*=================*/

/*---------------------------------------------------*
* smoothUpper: ceiling.
*      y = xmax - 0.5 ( arg + sqrt( arg^2 + 4 xmax delta ) )
*    arg = xmax - x - delta
*-----------------*/

#define Fn_SU( y , x , xmax , delta , dx ) { \
    TMF1 = ( xmax ) - ( x ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmax ) * ( delta) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    y = ( xmax ) - 0.5 * ( TMF1 + TMF2 ) ; \
  }

#define Fn_SU2( y , x , xmax , delta , dy_dx , dy_dxmax ) { \
    TMF1 = ( xmax ) - ( x ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmax ) * ( delta) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dy_dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    dy_dxmax = 0.5 * ( 1.0 - ( TMF1 + 2.0 * delta ) / TMF2 ) ; \
    y = ( xmax ) - 0.5 * ( TMF1 + TMF2 ) ; \
  }


/*---------------------------------------------------*
* smoothLower: flooring.
*      y = xmin + 0.5 ( arg + sqrt( arg^2 + 4 xmin delta ) )
*    arg = x - xmin - delta
*-----------------*/

#define Fn_SL( y , x , xmin , delta , dx ) { \
    TMF1 = ( x ) - ( xmin ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmin ) * ( delta ) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    y = ( xmin ) + 0.5 * ( TMF1 + TMF2 ) ; \
  }

#define Fn_SL2( y , x , xmin , delta , dy_dx, dy_dxmin ) { \
    TMF1 = ( x ) - ( xmin ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmin ) * ( delta ) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dy_dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    dy_dxmin = 0.5 * ( 1.0 - ( TMF1 - 2.0 * delta ) / TMF2 ) ; \
    y = ( xmin ) + 0.5 * ( TMF1 + TMF2 ) ; \
  }

/*---------------------------------------------------*
* smoothZero: flooring to zero.
*      y = 0.5 ( x + sqrt( x^2 + 4 delta^2 ) )
*-----------------*/

#define Fn_SZ( y , x , delta , dx ) { \
    TMF2 = sqrt ( ( x ) *  ( x ) + 4.0 * ( delta ) * ( delta) ) ; \
    dx = 0.5 * ( 1.0 + ( x ) / TMF2 ) ; \
    y = 0.5 * ( ( x ) + TMF2 ) ; \
    if( y < 0.0 ) { y=0.0; dx=0.0; } \
  }


/*---------------------------------------------------*
* CeilingPow: ceiling for positive x, flooring for negative x.
*      y = x * xmax / ( x^{2m} + xmax^{2m} )^{1/(2m)}
* note: 
*   - xmax has to be positive.
*   - -xmax < y < xmax.
*   - dy/dx|_{x=0} = 1.
*-----------------*/

#define Fn_CP( y , x , xmax , pw , dx ) { \
  double x2 = (x) * (x) ; \
  double xmax2 = (xmax) * (xmax) ; \
  double xp = 1.0 , xmp = 1.0 ; \
  int   m =0, mm =0; \
  double arg =0.0, dnm =0.0; \
  for ( m = 0 ; m < pw ; m ++ ) { xp *= x2 ; xmp *= xmax2 ; } \
  arg = xp + xmp ; \
  dnm = arg ; \
  if ( pw == 1 || pw == 2 || pw == 4 || pw == 8 ) { \
    if ( pw == 1 ) { mm = 1 ; \
    } else if ( pw == 2 ) { mm = 2 ; \
    } else if ( pw == 4 ) { mm = 3 ; \
    } else if ( pw == 8 ) { mm = 4 ; } \
    for ( m = 0 ; m < mm ; m ++ ) { dnm = sqrt( dnm ) ; } \
  } else { dnm = Fn_Pow( dnm , 1.0 / ( 2.0 * pw ) ) ; } \
  dnm = 1.0 / dnm ; \
  y = (x) * (xmax) * dnm ; \
  dx = (xmax) * xmp * dnm / arg ; \
}

/*---------------------------------------------------*
* Declining function using a polynomial.
*-----------------*/

#define Fn_DclPoly4( y , x , dx ) { \
  TMF2 = (x) * (x) ; \
  TMF3 = TMF2 * (x) ; \
  TMF4 = TMF2 * TMF2 ; \
  y = 1.0 / ( 1.0 + (x) + TMF2 + TMF3 + TMF4 ) ; \
  dx = - ( 1.0 + 2.0 * (x) + 3.0 * TMF2 + 4.0 * TMF3 )  * y * y  ; \
} 

#define Fn_DclPoly6( y , x , dx ) { \
  TMF2 = (x) * (x) ; \
  TMF3 = TMF2 * (x) ; \
  TMF4 = TMF2 * TMF2 ; \
  TMF5 = TMF2 * TMF3 ; \
  TMF6 = TMF3 * TMF3 ; \
  y = 1.0 / ( 1.0 + (x) + TMF2 + TMF3 + TMF4 + TMF5 + TMF6 ) ; \
  dx = - ( 1.0 + 2.0 * (x) + 3.0 * TMF2 + 4.0 * TMF3 \
         + 5.0 * TMF4 + 6.0 * TMF5 )  * y * y  ; \
}

/*---------------------------------------------------*
* "smoothUpper" using a polynomial
*-----------------*/

#define Fn_SUPoly4( y , x , xmax , dx ) { \
 TMF1 = (x) / xmax ; \
 Fn_DclPoly4( y , TMF1 , dx ) ; \
 y = xmax * ( 1.0 - y ) ; \
 dx = - dx ; \
}
 

#define Fn_SUPoly4m( y , x , xmax , dx , dxmax ) { \
 TMF1 = (x) / xmax ; \
 Fn_DclPoly4( TMF0 , TMF1 , dx ) ; \
 y = xmax * ( 1.0 - TMF0 ) ; \
 dxmax = 1.0 - TMF0 + TMF1 * dx ; \
 dx = - dx ; \
}


#define Fn_SUPoly6m( y , x , xmax , dx , dxmax ) { \
 TMF1 = (x) / xmax ; \
 Fn_DclPoly6( TMF0 , TMF1 , dx ) ; \
 y = xmax * ( 1.0 - TMF0 ) ; \
 dxmax = 1.0 - TMF0 + TMF1 * dx ; \
 dx = - dx ; \
}

/*---------------------------------------------------*
* SymAdd: evaluate additional term for symmetry.
*-----------------*/

#define Fn_SymAdd( y , x , add0 , dx ) \
{ \
    TMF1 = 2.0 * ( x ) / ( add0 ) ; \
    TMF2 = 1.0 + TMF1 * ( (1.0/2) + TMF1 * ( (1.0/6) \
               + TMF1 * ( (1.0/24) + TMF1 * ( (1.0/120) \
               + TMF1 * ( (1.0/720) + TMF1 * (1.0/5040) ) ) ) ) ) ; \
    TMF3 = (1.0/2) + TMF1 * ( (1.0/3) \
               + TMF1 * ( (1.0/8) + TMF1 * ( (1.0/30) \
               + TMF1 * ( (1.0/144) + TMF1 * (1.0/840) ) ) ) ) ; \
    y = add0 / TMF2 ; \
    dx = - 2.0 * TMF3 / ( TMF2 * TMF2 ) ; \
}

/*===========================================================*
* Function hsmhvevaluate.
*=================*/
int HSMHVevaluate
(
 double        Vdse, /* external branch voltage ( Vds >= 0 are assumed -> Vdse might be negative.) */
 double        Vgse, /* external branch voltage */
 double        Vbse, /* external branch voltage */
 double        Vds, /* inner branch voltage */
 double        Vgs, /* inner branch voltage */
 double        Vbs, /* inner branch voltage */
 double        vbs_jct,
 double        vbd_jct,
 double        Vsubs, /* substrate-source voltage */
 double        deltemp,
 HSMHVinstance *here,
 HSMHVmodel    *model,
 CKTcircuit   *ckt
 ) 
{
  HSMHVbinningParam *pParam = &here->pParam ;
  HSMHVmodelMKSParam *modelMKS = &model->modelMKS ;
  HSMHVhereMKSParam  *hereMKS = &here->hereMKS ;
  /*-----------------------------------*
   * Constants for Smoothing functions
   *---------------*/
  const double vth_dlt = 1.0e-3 ;
  /*  const double cclmmdf = 1.0e-2 ;*/
  const double cclmmdf = 1.0e-1 ;
  const double C_cm2m  = 1.0e-2 ;
  const double qme_dlt = 1.0e-9 * C_cm2m ;
  const double rdsl2_dlt = 10.0e-3 * C_cm2m ; 
  const double rdsu2_dlt = 50.0e-6 * C_cm2m ;
  const double rdsz_dlt  = 0.1e-3 * C_cm2m ;
  const double qme2_dlt = 5.0e-2 ;
  const double eef_dlt = 1.0e-2 / C_cm2m ;
  const double sti2_dlt = 2.0e-3 ;
  const double pol_dlt = 5.0e-2 ; 
  const double psia_dlt = 1.0e-3 ;
  const double psia2_dlt = 5.0e-3 ;
  const double psisti_dlt = 5.0e-3 ;

  /*---------------------------------------------------*
   * Local variables. 
   *-----------------*/
  /* Constants ----------------------- */
  const int lp_s0_max = 20 ;
  const int lp_sl_max = 40 ;
  const double dP_max  = 0.1e0 ;
  const double ps_conv = 1.0e-12 ;
  /* double  ps_conv = 1.0e-13 ;*/
  const double gs_conv = 1.0e-8 ;
  /** depletion **/
  const double znbd3 = 3.0e0 ;
  const double znbd5 = 5.0e0 ;
  const double cn_nc3 = C_SQRT_2 / 108e0 ;
  /* 5-degree, contact:Chi=5 */
  const double cn_nc51 =  0.707106781186548 ; /* sqrt(2)/2 */
  const double cn_nc52 = -0.117851130197758 ; /* -sqrt(2)/12 */
  const double cn_nc53 =  0.0178800506338833 ; /* (187 - 112*sqrt(2))/1600 */
  const double cn_nc54 = -0.00163730162779191 ; /* (-131 + 88*sqrt(2))/4000 */
  const double cn_nc55 =  6.36964918866352e-5 ; /* (1509-1040*sqrt(2))/600000 */
  /** inversion **/
  /* 3-degree polynomial approx for ( exp[Chi]-1 )^{1/2} */
  const double cn_im53 =  2.9693154855770998e-1 ;
  const double cn_im54 = -7.0536542840097616e-2 ;
  const double cn_im55 =  6.1152888951331797e-3 ;
  /** initial guess **/
  const double c_ps0ini_2 = 8.0e-4 ;
  const double c_pslini_1 = 0.3e0 ;
  const double c_pslini_2 = 3.0e-2 ;

  const double VgVt_small = 1.0e-12 ;
  /* const double Vbs_min = model->HSMHV_vbsmin ; */
  /* const double Vds_max = 10.5e0 ; */
  /* const double Vgs_max = 10.5e0 ; */
  const double epsm10 = 10.0e0 * C_EPS_M ;
  const double small = 1.0e-50 ;
  const double small2= 1e-12 ;      /* for Ra(Vdse) dependence */
  const double c_exp_2 = 7.38905609893065 ;
  const double large_arg = 80 ; // log(1.0e100) ;

  double Vbs_max = 0.8e0, Vbs_max_dT =0.0 ; 
  double Vbs_bnd = 0.4e0, Vbs_bnd_dT =0.0 ; /* start point of positive Vbs bending */
  double Gdsmin = 0.0 ;
  double Gjmin = ckt->CKTgmin ;

  /* Internal flags  --------------------*/
  int flg_err = 0 ;   /* error level */
  int flg_rsrd = 0 ;  /* Flag for handling Rs and Rd */
  /* int flg_iprv = 0 ;  */ /* Flag for initial guess of Ids  -> not necessary any more */
  int flg_pprv = 0 ;  /* Flag for initial guesses of Ps0 and Pds */
  int flg_noqi =0;    /* Flag for the cases regarding Qi=Qd=0 */
  /* int flg_vbsc = 0 ; */ /* Flag for Vbs confining  -> not necessary any more */
  int flg_info = 0 ; 
  int flg_conv = 0 ;  /* Flag for Poisson loop convergence */
  int flg_qme = 0 ;   /* Flag for QME */
  int flg_nqs=0 ;     /* Flag for NQS calculation */
  
  /* Important Variables in HiSIM -------*/
  /* confined bias */
  double Vbscl=0.0, Vbscl_dVbs=0.0, Vbscl_dT=0.0, Vbscl_dVbs_dVbs = 0.0, Vbscl_dVbs_dT = 0.0 ;
  double Vgp =0.0, Vgp_dVbs =0.0, Vgp_dVds =0.0, Vgp_dVgs =0.0, Vgp_dT =0.0 ;
  double Vgs_fb =0.0 ;
  /* Ps0 : surface potential at the source side */
  double Ps0 =0.0, Ps0_dVbs =0.0, Ps0_dVds =0.0, Ps0_dVgs =0.0, Ps0_dT =0.0 ;
  double Ps0_ini =0.0 ;
  double Ps0_iniA =0.0, Ps0_iniA_dVxb =0.0, Ps0_iniA_dVgb =0.0, Ps0_iniA_dT =0.0 ;
  double Ps0_iniB =0.0;/*, Ps0_iniB_dVxb =0.0, Ps0_iniB_dVgb =0.0, Ps0_iniB_dT =0.0 ;*/
  /* Psl : surface potential at the drain side */
  double Psl =0.0, Psl_dVbs =0.0, Psl_dVds =0.0, Psl_dVgs =0.0, Psl_dT =0.0 ;
  double Psl_lim =0.0, dPlim =0.0 ;
  /* Pds := Psl - Ps0 */
  double Pds = 0.0, Pds_dVbs = 0.0, Pds_dVds = 0.0, Pds_dVgs  =0.0, Pds_dT =0.0 ;
  double Pds_ini =0.0 ;
  double Pds_max =0.0 ;
  /* iteration numbers of Ps0 and Psl equations. */
  int lp_s0 = 0 , lp_sl = 0 ;
  /* Xi0 := beta * ( Ps0 - Vbs ) - 1. */
  double Xi0 =0.0, Xi0_dVbs =0.0, Xi0_dVds =0.0, Xi0_dVgs =0.0, Xi0_dT =0.0 ;
  double Xi0p12 =0.0, Xi0p12_dVbs =0.0, Xi0p12_dVds =0.0, Xi0p12_dVgs =0.0, Xi0p12_dT =0.0 ;
  double Xi0p32 =0.0, Xi0p32_dVbs =0.0, Xi0p32_dVds =0.0, Xi0p32_dVgs =0.0, Xi0p32_dT =0.0 ;
  /* Xil := beta * ( Psl - Vbs ) - 1. */
  double Xil =0.0,    Xil_dVbs =0.0,    Xil_dVds =0.0,    Xil_dVgs =0.0,    Xil_dT =0.0 ;
  double Xilp12 =0.0, Xilp12_dVbs =0.0, Xilp12_dVds =0.0, Xilp12_dVgs =0.0, Xilp12_dT =0.0 ;
  double Xilp32 =0.0, Xilp32_dVbs =0.0, Xilp32_dVds =0.0, Xilp32_dVgs =0.0, Xilp32_dT =0.0 ;
  /* modified bias and potential for sym.*/
  double Vbsz =0.0,  Vbsz_dVbs =0.0,  Vbsz_dVds =0.0,  Vbsz_dT =0.0 ;
  double Vdsz =0.0,  Vdsz_dVbs =0.0,  Vdsz_dVds =0.0,  Vdsz_dT =0.0 ;
  double Vgsz =0.0,  Vgsz_dVbs =0.0,  Vgsz_dVds =0.0,  Vgsz_dVgs =0.0, Vgsz_dT =0.0  ;
  double Vzadd =0.0, Vzadd_dVbs =0.0, Vzadd_dVds = 0.0, Vzadd_dT = 0.0 ;
  double Ps0z =0.0,  Ps0z_dVbs =0.0,  Ps0z_dVds =0.0,  Ps0z_dVgs =0.0,  Ps0z_dT =0.0 ;
  double Pzadd =0.0, Pzadd_dVbs =0.0, Pzadd_dVds =0.0, Pzadd_dVgs =0.0, Pzadd_dT =0.0 ;

  /* IBPC */
  double dVbsIBPC =0.0, dVbsIBPC_dVbs =0.0, dVbsIBPC_dVds =0.0, dVbsIBPC_dVgs =0.0, dVbsIBPC_dT =0.0 ;
  double dG3 =0.0,      dG3_dVbs =0.0,      dG3_dVds =0.0,      dG3_dVgs =0.0,      dG3_dT =0.0 ;
  double dG4 =0.0,      dG4_dVbs =0.0,      dG4_dVds =0.0,      dG4_dVgs =0.0,      dG4_dT =0.0 ;
  double dIdd =0.0,     dIdd_dVbs =0.0,     dIdd_dVds =0.0,     dIdd_dVgs =0.0,     dIdd_dT =0.0 ;

  double betaWL =0.0, betaWL_dVbs =0.0, betaWL_dVds =0.0, betaWL_dVgs =0.0, betaWL_dT =0.0 ;

  /* Chi := beta * ( Ps{0/l} - Vbs ) */
  double Chi =0.0, Chi_dVbs =0.0, Chi_dVds =0.0, Chi_dVgs =0.0, Chi_dT =0.0 ;
  /* Rho := beta * ( Psl - Vds ) */
  double Rho =0.0, Rho_dT =0.0 ;
  /* threshold voltage */
  double Vth =0.0 ;
  double Vth0 =0.0, Vth0_dVb =0.0, Vth0_dVd =0.0, Vth0_dVg =0.0, Vth0_dT =0.0 ;
  /* variation of threshold voltage */
  double dVth =0.0, dVth_dVb =0.0, dVth_dVd =0.0, dVth_dVg =0.0, dVth_dT = 0.0 ;
  double dVth0 =0.0 ;
  double dVth0_dVb =0.0, dVth0_dVd =0.0, dVth0_dVg =0.0, dVth0_dT =0.0 ;
  double dVthSC =0.0 ;
  double dVthSC_dVb =0.0, dVthSC_dVd =0.0, dVthSC_dVg =0.0, dVthSC_dT =0.0 ;
  double delta0 = 5.0e-3 ;
  double Psi_a =0.0, Psi_a_dVg =0.0, Psi_a_dVb =0.0, Psi_a_dVd =0.0, Psi_a_dT =0.0 ;
  double Pb20a =0.0, Pb20a_dVg =0.0, Pb20a_dVb =0.0, Pb20a_dVd =0.0, Pb20a_dT =0.0 ;
  double Pb20b =0.0, Pb20b_dVg =0.0, Pb20b_dVb =0.0, Pb20b_dVd =0.0, Pb20b_dT =0.0 ;
  double dVthW =0.0, dVthW_dVb =0.0, dVthW_dVd =0.0, dVthW_dVg =0.0, dVthW_dT =0.0 ;
  /* Alpha and related parameters */
  double Alpha =0.0, Alpha_dVbs =0.0, Alpha_dVds =0.0, Alpha_dVgs =0.0, Alpha_dT =0.0 ;
  double Achi =0.0,  Achi_dVbs =0.0,  Achi_dVds =0.0,  Achi_dVgs =0.0,  Achi_dT =0.0 ;
  double VgVt =0.0,  VgVt_dVbs =0.0,  VgVt_dVds =0.0,  VgVt_dVgs =0.0,  VgVt_dT =0.0 ;
  double Pslsat =0.0 ;
  double Vdsat =0.0 ;
  double VdsatS =0.0, VdsatS_dVbs =0.0, VdsatS_dVds =0.0, VdsatS_dVgs =0.0, VdsatS_dT =0.0 ;
  double Delta =0.0 ;
  /* Q_B and capacitances */
  double Qb =0.0,  Qb_dVbs =0.0,  Qb_dVds =0.0,  Qb_dVgs =0.0,  Qb_dT=0.0 ;
  double Qbu =0.0, Qbu_dVbs =0.0, Qbu_dVds =0.0, Qbu_dVgs =0.0, Qbu_dT =0.0 ;
  /* Q_I and capacitances */
  double Qi =0.0,  Qi_dVbs =0.0,  Qi_dVds =0.0,  Qi_dVgs =0.0,  Qi_dT=0.0 ;
  double Qiu =0.0, Qiu_dVbs =0.0, Qiu_dVds =0.0, Qiu_dVgs =0.0, Qiu_dT =0.0 ;
  /* Q_D and capacitances */
  double Qd =0.0,  Qd_dVbs =0.0,  Qd_dVds =0.0,  Qd_dVgs =0.0,  Qd_dT =0.0 ;
  /* channel current */
  double Ids =0.0,  Ids_dVbs =0.0,  Ids_dVds =0.0,  Ids_dVgs =0.0,  Ids_dT =0.0, Ids_dRa =0.0 ;
  double Ids0 =0.0, Ids0_dVbs =0.0, Ids0_dVds =0.0, Ids0_dVgs =0.0, Ids0_dT =0.0 ;
  /* STI */
  double dVthSCSTI =0.0, dVthSCSTI_dVg =0.0, dVthSCSTI_dVd =0.0, dVthSCSTI_dVb =0.0, dVthSCSTI_dT =0.0 ;
  double Vgssti =0.0,    Vgssti_dVbs =0.0,   Vgssti_dVds =0.0,   Vgssti_dVgs =0.0,   Vgssti_dT =0.0 ;
  double costi0 =0.0 ;
  double costi1 =0.0, costi1_dT =0.0 ;
  double costi3 =0.0, costi3_dVb =0.0,    costi3_dVd =0.0,    costi3_dVg =0.0,    costi3_dT =0.0 ;
  double              costi3_dVb_c3 =0.0, costi3_dVd_c3 =0.0, costi3_dVg_c3 =0.0;/*, costi3_dT_c3 =0.0 ;*/
  double costi4 =0.0, costi4_dT =0.0 ;
  double costi5 =0.0, costi5_dT =0.0 ;
  double costi6 =0.0, costi6_dT =0.0 ;
  double costi7 =0.0, costi7_dT =0.0 ;
  double Psasti =0.0, Psasti_dVbs =0.0, Psasti_dVds =0.0, Psasti_dVgs =0.0, Psasti_dT =0.0 ;
  double Psbsti =0.0, Psbsti_dVbs =0.0, Psbsti_dVds =0.0, Psbsti_dVgs =0.0, Psbsti_dT =0.0 ;
  double Psab =0.0,   Psab_dVbs =0.0,   Psab_dVds =0.0,   Psab_dVgs =0.0,   Psab_dT =0.0 ;
  double Psti =0.0,   Psti_dVbs =0.0,   Psti_dVds =0.0,   Psti_dVgs =0.0,   Psti_dT =0.0 ;
  double sq1sti =0.0, sq1sti_dVbs =0.0, sq1sti_dVds =0.0, sq1sti_dVgs =0.0, sq1sti_dT =0.0 ;
  double sq2sti =0.0, sq2sti_dVbs =0.0, sq2sti_dVds =0.0, sq2sti_dVgs =0.0, sq2sti_dT =0.0 ;
  double Qn0sti =0.0, Qn0sti_dVbs =0.0, Qn0sti_dVds =0.0, Qn0sti_dVgs =0.0, Qn0sti_dT =0.0 ;
  double Idssti =0.0, Idssti_dVbs =0.0, Idssti_dVds =0.0, Idssti_dVgs =0.0, Idssti_dT=0.0 ;
  /* constants */
  double beta =0.0,     beta_dT =0.0 ;
  double beta_inv =0.0, beta_inv_dT =0.0 ;
  double beta2 =0.0 ;
  double Pb2 =0.0,      Pb2_dT =0.0 ;
  double Pb20 =0.0 ;
  double Pb2c =0.0 ;
  double Vfb =0.0 ;
  double c_eox =0.0 ;
  double Leff =0.0, Weff =0.0, WeffLD_nf =0.0, Ldrift =0.0 ;
  double Ldrift0 =0.0 ;
  double q_Nsub =0.0 ;
  /* PART-1 */
  /* Accumulation zone */
  double Psa =0.0 ;
  double Psa_dVbs =0.0, Psa_dVds =0.0, Psa_dVgs =0.0, Psa_dT =0.0 ;
  /* CLM*/
  double Psdl =0.0, Psdl_dVbs =0.0, Psdl_dVds =0.0, Psdl_dVgs =0.0, Psdl_dT =0.0 ;
  double Lred =0.0, Lred_dVbs =0.0, Lred_dVds =0.0, Lred_dVgs =0.0, Lred_dT =0.0 ;
  double Lch =0.0,  Lch_dVbs =0.0,  Lch_dVds =0.0,  Lch_dVgs =0.0,  Lch_dT =0.0 ;
  double Wd =0.0,   Wd_dVbs =0.0,   Wd_dVds =0.0,   Wd_dVgs =0.0,   Wd_dT =0.0 ;
  double Aclm =0.0 ;
  /* Pocket Implant */
  double Vthp=0.0,   Vthp_dVb=0.0,   Vthp_dVd=0.0,   Vthp_dVg =0.0,   Vthp_dT =0.0 ;
  double dVthLP=0.0, dVthLP_dVb=0.0, dVthLP_dVd=0.0, dVthLP_dVg =0.0, dVthLP_dT =0.0 ;
  double bs12=0.0,   bs12_dVb=0.0,   bs12_dVd =0.0,  bs12_dVg =0.0, bs12_dT =0.0 ;
  double Qbmm=0.0,   Qbmm_dVb=0.0,   Qbmm_dVd =0.0,  Qbmm_dVg =0.0, Qbmm_dT =0.0 ;
  double dqb=0.0,    dqb_dVb=0.0,    dqb_dVg=0.0,    dqb_dVd =0.0,  dqb_dT =0.0 ;
  double Vdx=0.0,    Vdx_dVbs=0.0;/*,   Vdx_dT=0.0 ;*/
  double Vdx2=0.0,   Vdx2_dVbs=0.0;/*,  Vdx2_dT=0.0 ;*/    
  double Pbsum=0.0,  Pbsum_dVb=0.0,  Pbsum_dVd=0.0,  Pbsum_dVg =0.0,  Pbsum_dT =0.0 ;
  double sqrt_Pbsum =0.0 ;
  /* Poly-Depletion Effect */
  const double pol_b = 1.0 ;
  double dPpg =0.0,  dPpg_dVb =0.0,  dPpg_dVd =0.0,  dPpg_dVg =0.0,   dPpg_dT = 0.0 ; 
  /* Quantum Effect */
  double Tox =0.0,     Tox_dVb =0.0,     Tox_dVd =0.0,     Tox_dVg =0.0,  Tox_dT =0.0 ;
  double dTox =0.0,    dTox_dVb =0.0,    dTox_dVd =0.0,    dTox_dVg =0.0, dTox_dT =0.0  ;
  double Cox =0.0,     Cox_dVb =0.0,     Cox_dVd =0.0,     Cox_dVg =0.0,  Cox_dT =0.0 ;
  double Cox_inv =0.0, Cox_inv_dVb =0.0, Cox_inv_dVd =0.0, Cox_inv_dVg =0.0, Cox_inv_dT =0.0 ;
  double Tox0 =0.0,    Cox0 =0.0,        Cox0_inv =0.0 ;
  double Vthq=0.0,     Vthq_dVb =0.0,    Vthq_dVd =0.0 ;
  /* Igate , Igidl , Igisl */
  const double igate_dlt = 1.0e-2 / C_cm2m ;
  const double gidlvds_dlt = 1.0e-5 ;
  const double gidla = 100.0 ;
  double Psdlz =0.0, Psdlz_dVbs =0.0, Psdlz_dVds =0.0, Psdlz_dVgs =0.0, Psdlz_dT =0.0 ;
  double Egp12 =0.0, Egp12_dT =0.0 ; 
  double Egp32 =0.0, Egp32_dT =0.0 ;
  double E1 =0.0,    E1_dVb =0.0,     E1_dVd =0.0,     E1_dVg =0.0,     E1_dT =0.0 ;
  double Etun =0.0,  Etun_dVbs =0.0,  Etun_dVds =0.0,  Etun_dVgs =0.0,  Etun_dT =0.0 ;
  double Vdsp=0.0,   Vdsp_dVd =0.0 ;
  double Egidl =0.0, Egidl_dVb =0.0,  Egidl_dVd =0.0,  Egidl_dVg =0.0,  Egidl_dT =0.0 ;
  double Egisl =0.0, Egisl_dVb =0.0,  Egisl_dVd =0.0,  Egisl_dVg =0.0,  Egisl_dT =0.0 ;
  double Igate =0.0, Igate_dVbs =0.0, Igate_dVds =0.0, Igate_dVgs =0.0, Igate_dT =0.0 ;
  double Igs =0.0,   Igs_dVbs =0.0,   Igs_dVds =0.0,   Igs_dVgs =0.0,   Igs_dT =0.0 ;
  double Igd =0.0,   Igd_dVbs =0.0,   Igd_dVds =0.0,   Igd_dVgs =0.0,   Igd_dT =0.0 ;
  double Igb =0.0,   Igb_dVbs =0.0,   Igb_dVds =0.0,   Igb_dVgs =0.0,   Igb_dT =0.0 ;
  double Igidl =0.0, Igidl_dVbs =0.0, Igidl_dVds =0.0, Igidl_dVgs =0.0, Igidl_dT =0.0 ;
  double Igisl =0.0, Igisl_dVbs =0.0, Igisl_dVds =0.0, Igisl_dVgs =0.0, Igisl_dT =0.0 ;
  double Vdb =0.0, Vsb =0.0 ;
  /* connecting function */
  double FD2 =0.0,    FD2_dVbs =0.0,    FD2_dVds =0.0,    FD2_dVgs =0.0,    FD2_dT =0.0 ;
  double FMDVDS =0.0, FMDVDS_dVbs =0.0, FMDVDS_dVds =0.0, FMDVDS_dVgs =0.0, FMDVDS_dT =0.0 ;
  double FMDVGS =0.0, FMDVGS_dVgs =0.0 ;
  double FMDPG =0.0,  FMDPG_dVbs =0.0,  FMDPG_dVds =0.0,  FMDPG_dVgs =0.0,  FMDPG_dT =0.0 ;

  double cnst0 =0.0,  cnst0_dT =0.0;
  double cnst1 =0.0,  cnst1_dT =0.0;
  double cnstCoxi =0.0, cnstCoxi_dVb =0.0, cnstCoxi_dVd =0.0, cnstCoxi_dVg =0.0, cnstCoxi_dT =0.0 ;
  double fac1 =0.0,     fac1_dVbs =0.0,    fac1_dVds =0.0,    fac1_dVgs =0.0,    fac1_dT =0.0 ;
  double fac1p2 =0.0,   fac1p2_dT =0.0 ;
  double fs01 =0.0, fs01_dVbs =0.0, fs01_dVds =0.0, fs01_dVgs =0.0, fs01_dT =0.0, fs01_dPs0 =0.0 ;
  double fs02 =0.0, fs02_dVbs =0.0, fs02_dVds =0.0, fs02_dVgs =0.0, fs02_dT =0.0, fs02_dPs0 =0.0 ;
  double fsl1 =0.0, fsl1_dVbs =0.0, fsl1_dVds =0.0, fsl1_dVgs =0.0, fsl1_dT =0.0, fsl1_dPsl =0.0 ;
  double fsl2 =0.0, fsl2_dVbs =0.0, fsl2_dVds =0.0, fsl2_dVgs =0.0, fsl2_dT =0.0, fsl2_dPsl =0.0 ;
  double cfs1 =0.0, cfs1_dT =0.0 ;
  double fb =0.0,   fb_dChi =0.0 ;
  double fi =0.0,   fi_dChi =0.0 ;
  double exp_Chi =0.0,     exp_Chi_dT =0.0 ;
  double exp_Rho =0.0,     exp_Rho_dT =0.0 ;
  double exp_bVbs =0.0,    exp_bVbs_dT =0.0 ; 
  double exp_bVbsVds =0.0, exp_bVbsVds_dT =0.0 ;
  double exp_bPs0 =0.0,    exp_bPs0_dT =0.0 ;
  double Fs0 =0.0,         Fs0_dPs0 =0.0 ; 
  double Fsl =0.0,         Fsl_dPsl =0.0 ;
  double dPs0 =0.0,        dPsl =0.0 ;
  double Qn0 =0.0,   Qn0_dVbs =0.0,   Qn0_dVds =0.0,   Qn0_dVgs =0.0,   Qn0_dT =0.0 ;
  double Qb0 =0.0,   Qb0_dVb =0.0,    Qb0_dVd =0.0,    Qb0_dVg =0.0,    Qb0_dT =0.0 ;
  double Qbnm =0.0,  Qbnm_dVbs =0.0,  Qbnm_dVds =0.0,  Qbnm_dVgs =0.0,  Qbnm_dT =0.0 ;
  double DtPds =0.0, DtPds_dVbs =0.0, DtPds_dVds =0.0, DtPds_dVgs =0.0, DtPds_dT =0.0 ;
  double Qinm =0.0,  Qinm_dVbs =0.0,  Qinm_dVds =0.0,  Qinm_dVgs =0.0,  Qinm_dT =0.0 ;
  double Qidn =0.0,  Qidn_dVbs =0.0,  Qidn_dVds =0.0,  Qidn_dVgs =0.0,  Qidn_dT =0.0 ;
  double Qdnm =0.0,  Qdnm_dVbs =0.0,  Qdnm_dVds =0.0,  Qdnm_dVgs =0.0,  Qdnm_dT =0.0 ;
  double Qddn =0.0,  Qddn_dVbs =0.0,  Qddn_dVds =0.0,  Qddn_dVgs =0.0,  Qddn_dT =0.0 ;
  double Quot =0.0 ;
  double Qdrat =0.5, Qdrat_dVbs =0.0, Qdrat_dVds =0.0, Qdrat_dVgs =0.0, Qdrat_dT =0.0 ;
  double Idd =0.0,   Idd_dVbs =0.0,   Idd_dVds =0.0,   Idd_dVgs =0.0,   Idd_dT =0.0 ;
  double Fdd =0.0,   Fdd_dVbs =0.0,   Fdd_dVds =0.0,   Fdd_dVgs =0.0,   Fdd_dT =0.0 ;
  double Eeff =0.0,  Eeff_dVbs =0.0,  Eeff_dVds =0.0,  Eeff_dVgs =0.0,  Eeff_dT =0.0 ;
  double Rns =0.0,   Rns_dT =0.0 ;
  double Mu = 0.0,   Mu_dVbs =0.0,    Mu_dVds =0.0,    Mu_dVgs =0.0,    Mu_dT =0.0 ;
  double Muun =0.0,  Muun_dVbs =0.0,  Muun_dVds =0.0,  Muun_dVgs =0.0,  Muun_dT =0.0 ;
  double Ey =0.0,    Ey_dVbs =0.0,    Ey_dVds =0.0,    Ey_dVgs =0.0,    Ey_dT =0.0 ;
  double Em =0.0,    Em_dVbs =0.0,    Em_dVds =0.0,    Em_dVgs =0.0,    Em_dT =0.0 ;
  double Vmax =0.0,  Vmax_dT =0.0 ;
  double Eta =0.0,   Eta_dVbs =0.0,   Eta_dVds =0.0,   Eta_dVgs =0.0,   Eta_dT =0.0 ;
  double Eta1 =0.0,    Eta1_dT =0.0 ;
  double Eta1p12 =0.0, Eta1p12_dT =0.0 ;
  double Eta1p32 =0.0, Eta1p32_dT =0.0 ; 
  double Eta1p52 =0.0, Eta1p52_dT =0.0 ;
  double Zeta12 =0.0,  Zeta12_dT =0.0 ;
  double Zeta32 =0.0,  Zeta32_dT =0.0 ;
  double Zeta52 =0.0,  Zeta52_dT =0.0 ;
  double F00 =0.0,   F00_dVbs =0.0,   F00_dVds =0.0,   F00_dVgs =0.0,   F00_dT =0.0 ;
  double F10 =0.0,   F10_dVbs =0.0,   F10_dVds =0.0,   F10_dVgs =0.0,   F10_dT =0.0 ;
  double F30 =0.0,   F30_dVbs =0.0,   F30_dVds =0.0,   F30_dVgs =0.0,   F30_dT =0.0 ;
  double F11 =0.0,   F11_dVbs =0.0,   F11_dVds =0.0,   F11_dVgs =0.0,   F11_dT =0.0 ;
  double Ps0_min =0.0, Ps0_min_dT =0.0 ;
  double Acn =0.0,   Acn_dVbs =0.0,   Acn_dVds =0.0,   Acn_dVgs =0.0,   Acn_dT =0.0 ;
  double Acd =0.0,   Acd_dVbs =0.0,   Acd_dVds =0.0,   Acd_dVgs =0.0,   Acd_dT =0.0 ;
  double Ac1 =0.0,   Ac1_dVbs =0.0,   Ac1_dVds =0.0,   Ac1_dVgs =0.0,   Ac1_dT =0.0 ;
  double Ac2 =0.0,   Ac2_dVbs =0.0,   Ac2_dVds =0.0,   Ac2_dVgs =0.0,   Ac2_dT =0.0 ;
  double Ac3 =0.0,   Ac3_dVbs =0.0,   Ac3_dVds =0.0,   Ac3_dVgs =0.0,   Ac3_dT =0.0 ;
  double Ac4 =0.0,   Ac4_dVbs =0.0,   Ac4_dVds =0.0,   Ac4_dVgs =0.0,   Ac4_dT =0.0 ;
  double Ac31 =0.0,  Ac31_dVbs =0.0,  Ac31_dVds =0.0,  Ac31_dVgs =0.0,  Ac31_dT =0.0 ;
  double Ac41 =0.0,                                                     Ac41_dT =0.0 ;
  double ninvd_dT =0.0 ;
  /* PART-2 (Isub) */
  double Isub =0.0,      Isub_dVbs =0.0,     Isub_dVds =0.0,     Isub_dVgs =0.0,     Isub_dT=0.0 ;
  double Isub_dVdse = 0.0 ;
  double Psislsat =0.0,  Psislsat_dVb =0.0,  Psislsat_dVd =0.0,  Psislsat_dVg =0.0,  Psislsat_dT =0.0 ;
  double Psisubsat =0.0, Psisubsat_dVb =0.0, Psisubsat_dVd =0.0, Psisubsat_dVg =0.0, Psisubsat_dT =0.0 ;
  double Ifn =0.0,       Ifn_dVb =0.0,       Ifn_dVd=0.0,        Ifn_dVg=0.0,        Ifn_dT = 0.0 ;
  double Eg12=0.0,       Eg32 =0.0 ;
  /* PART-3 (overlap) */
  /* const double cov_dlt = 1.0e-1 ; */
  /* const double covvgmax = 5.0 ;   */
  double cov_slp =0.0, cov_mag =0.0 ;
  double Qgos =0.0,    Qgos_dVbs =0.0,  Qgos_dVds =0.0,  Qgos_dVgs =0.0,  Qgos_dT =0.0 ;
  double Qgod =0.0,    Qgod_dVbs =0.0,  Qgod_dVds =0.0,  Qgod_dVgs =0.0,  Qgod_dT =0.0 ;
  double Qgbo =0.0,    Qgbo_dVbs =0.0,  Qgbo_dVds =0.0,  Qgbo_dVgs =0.0,  Qgbo_dT = 0.0 ;
  double Cgdo =0.0,    Cgso =0.0,       Cgbo_loc =0.0 ;
  double Qgso =0.0,    Qgso_dVbse =0.0, Qgso_dVdse =0.0, Qgso_dVgse =0.0 ;
  double Qgdo =0.0,    Qgdo_dVbse =0.0, Qgdo_dVdse =0.0, Qgdo_dVgse =0.0 ;
  /* fringe capacitance */
  double Qfd =0.0,     Cfd =0.0 ;
  double Qfs =0.0,     Cfs =0.0 ;
  /* Cqy */
  double Ec =0.0,      Ec_dVbs =0.0,    Ec_dVds =0.0,    Ec_dVgs =0.0,    Ec_dT =0.0 ;
  double Pslk =0.0,    Pslk_dVbs =0.0,  Pslk_dVds =0.0,  Pslk_dVgs =0.0,  Pslk_dT =0.0 ;
  double Qy =0.0,      Qy_dVbs =0.0,    Qy_dVds =0.0,    Qy_dVgs =0.0,    Qy_dT =0.0  ;
  /* PART-4 (junction diode) */
  double Ibs =0.0, Gbs =0.0, Gbse =0.0, Ibs_dT =0.0 ;
  double Ibd =0.0, Gbd =0.0, Gbde =0.0, Ibd_dT =0.0 ;
/*  double Nvtm =0.0 ;*/
  /* junction capacitance */
  double Qbs =0.0, Capbs =0.0, Capbse =0.0, Qbs_dT =0.0 ;
  double Qbd =0.0, Capbd =0.0, Capbde =0.0, Qbd_dT =0.0 ;
  double czbd =0.0,    czbd_dT=0.0 ;
  double czbdsw =0.0,  czbdsw_dT=0.0 ;
  double czbdswg =0.0, czbdswg_dT=0.0 ;
  double czbs =0.0,    czbs_dT=0.0 ;
  double czbssw =0.0,  czbssw_dT=0.0 ;
  double czbsswg =0.0, czbsswg_dT=0.0 ;
  double arg =0.0,     sarg =0.0 ;
  /* PART-5 (NQS) */
  double tau =0.0,  tau_dVbs=0.0,  tau_dVds=0.0,  tau_dVgs =0.0,  tau_dT=0.0  ;
  double taub =0.0, taub_dVbs=0.0, taub_dVds=0.0, taub_dVgs =0.0, taub_dT =0.0 ;
  /* PART-6 (noise) */
  /* 1/f */
  double NFalp =0.0, NFtrp =0.0, Cit =0.0, Nflic =0.0 ;
  /* thermal */
  double Eyd =0.0, Mu_Ave= 0.0, Nthrml =0.0, Mud_hoso =0.0 ;
  /* induced gate noise ( Part 0/3 ) */
  double kusai00 =0.0, kusaidd =0.0, kusaiL =0.0, kusai00L =0.0 ;
  int flg_ign = 0 ;
  double sqrtkusaiL =0.0, kusai_ig =0.0, gds0_ign =0.0, gds0_h2 =0.0, GAMMA =0.0, crl_f =0.0 ;
  const double c_sqrt_15 =3.872983346207417e0 ; /* sqrt(15) */
  const double Cox_small =1.0e-6 ;
  const double c_16o135 =1.185185185185185e-1 ; /* 16/135 */
  double Nign0 =0.0, MuModA =0.0, MuModB =0.0, correct_w1 =0.0 ;

  /* usage of previously calculated values */
  double vtol_pprv =1.01e-1 ;
  double Vbsc_dif  =0.0, Vdsc_dif =0.0,  Vgsc_dif =0.0,  sum_vdif =0.0 ;
  double Vbsc_dif2 =0.0, Vdsc_dif2 =0.0, Vgsc_dif2 =0.0, sum_vdif2 =0.0 ;
  double dVbs =0.0,      dVds =0.0,      dVgs =0.0 ;

  /* temporary vars. & derivatives*/
  double TX =0.0, TX_dVbs =0.0, TX_dVds =0.0, TX_dVgs =0.0, TX_dT =0.0 ;
  double TY =0.0, TY_dVbs =0.0, TY_dVds =0.0, TY_dVgs =0.0, TY_dT =0.0 ;
  double T0 =0.0, T0_dVb =0.0, T0_dVd =0.0, T0_dVg =0.0, T0_dT =0.0 ;
  double T1 =0.0, T1_dVb =0.0, T1_dVd =0.0, T1_dVg =0.0, T1_dT =0.0, T1_dVdse_eff =0.0 ;
  double T2 =0.0, T2_dVb =0.0, T2_dVd =0.0, T2_dVg =0.0, T2_dT =0.0 ;
  double T3 =0.0, T3_dVb =0.0, T3_dVd =0.0, T3_dVg =0.0, T3_dT =0.0 ;
  double T4 =0.0, T4_dVb =0.0, T4_dVd =0.0, T4_dVg =0.0, T4_dT =0.0 ;
  double T5 =0.0, T5_dVb =0.0, T5_dVd =0.0, T5_dVg =0.0, T5_dT =0.0 ;
  double T6 =0.0, T6_dVb =0.0, T6_dVd =0.0, T6_dVg =0.0, T6_dT =0.0 ;
  double T7 =0.0, T7_dVb =0.0, T7_dVd =0.0, T7_dVg =0.0, T7_dT =0.0 ;
  double T8 =0.0, T8_dVb =0.0, T8_dVd =0.0, T8_dVg =0.0, T8_dT =0.0 ;
  double T9 =0.0, T9_dVb =0.0, T9_dVd =0.0, T9_dVg =0.0, T9_dT =0.0, T9_dVdse_eff =0.0 ;
  double T10 =0.0, T10_dVb =0.0, T10_dVd =0.0, T10_dVg =0.0, T10_dT =0.0 ;
  double T11 =0.0,                                           T11_dT =0.0 ;
  double T12 =0.0,                                           T12_dT =0.0 ;
  double T15 =0.0, T16 =0.0, T17 =0.0 ;
  double T2_dVdse = 0.0, T5_dVdse = 0.0 ;
  double T4_dVb_dT, T5_dVb_dT, T6_dVb_dT, T7_dVb_dT ;

  int    flg_zone =0 ;
  double Vfbsft =0.0, Vfbsft_dVbs =0.0, Vfbsft_dVds =0.0, Vfbsft_dVgs =0.0, Vfbsft_dT =0.0 ;

  /* Vdseff */
  double Vdseff =0.0, Vdseff_dVbs =0.0, Vdseff_dVds =0.0, Vdseff_dVgs =0.0, Vdseff_dT =0.0 ;
  double Vdsorg =0.0 ;

  /* D/S Overlap Charges: Qovd/Qovs */
  double CVDSOVER =0.0 ;
  double Qovdext =0.0, Qovdext_dVbse =0.0, Qovdext_dVdse =0.0, Qovdext_dVgse =0.0, Qovdext_dT =0.0 ;
  double Qovsext =0.0, Qovsext_dVbse =0.0, Qovsext_dVdse =0.0, Qovsext_dVgse =0.0, Qovsext_dT =0.0 ;
  double Qovd =0.0,    Qovd_dVbs =0.0,     Qovd_dVds =0.0,     Qovd_dVgs =0.0,     Qovd_dT =0.0 ;
  double Qovs =0.0,    Qovs_dVbs =0.0,     Qovs_dVds =0.0,     Qovs_dVgs =0.0,     Qovs_dT =0.0 ;
  double QbuLD =0.0,   QbuLD_dVbs =0.0,    QbuLD_dVds =0.0,    QbuLD_dVgs =0.0,    QbuLD_dT =0.0 ;
  double QbdLD =0.0,   QbdLD_dVbs =0.0,    QbdLD_dVds =0.0,    QbdLD_dVgs =0.0,    QbdLD_dT =0.0 ;
  double QbsLD =0.0;/*,   QbsLD_dVbs =0.0,    QbsLD_dVds =0.0,    QbsLD_dVgs =0.0,    QbsLD_dT =0.0 ;*/
  double QbdLDext =0.0, QbdLDext_dVbse =0.0, QbdLDext_dVdse =0.0, QbdLDext_dVgse =0.0, QbdLDext_dT =0.0 ;
  double QbsLDext =0.0;/*, QbsLDext_dVbse =0.0, QbsLDext_dVdse =0.0, QbsLDext_dVgse =0.0, QbsLDext_dT =0.0 ;*/

  /* Vgsz for SCE and PGD */
  double dmpacc =0.0,  dmpacc_dVbs =0.0,   dmpacc_dVds =0.0,   dmpacc_dVgs =0.0 ;
  double Vbsz2 =0.0,   Vbsz2_dVbs =0.0,    Vbsz2_dVds =0.0,    Vbsz2_dVgs =0.0 , Vbsz2_dT =0.0;

  /* Multiplication factor * number of gate fingers */
  double Mfactor = here->HSMHV_m ;


  /*-----------------------------------------------------------*
   * HiSIM-HV
   *-----------------*/
  /* bias-dependent Rd, Rs */
  double Rdrift =0.0,  Rdrift_dVbse =0.0,  Rdrift_dVdse =0.0,  Rdrift_dVgse =0.0,  Rdrift_dT =0.0 ;
  double Rsdrift =0.0, Rsdrift_dVbse =0.0, Rsdrift_dVdse =0.0, Rsdrift_dVgse =0.0, Rsdrift_dT =0.0 ;
  double Rd =0.0,      Rd_dVbse =0.0,      Rd_dVdse =0.0,      Rd_dVgse =0.0,      Rd_dT =0.0 ;
  double Rs =0.0,      Rs_dVbse =0.0,      Rs_dVdse =0.0,      Rs_dVgse =0.0,      Rs_dT =0.0 ;
  double Ra =0.0,      Ra_dVbse =0.0,      Ra_dVdse =0.0,      Ra_dVgse =0.0 ; 
  double               Ra_dVbs =0.0,       Ra_dVds =0.0,       Ra_dVgs =0.0 ;
  double               Ra_dVdse_eff =0.0 ;
  const double delta_rd = 10e-3 * C_cm2m ;
  const double Ra_N = 20.0;	/* smoothing parameter for Ra */
  const double Res_min = 1.0e-4 ;
  double Rd0_dT =0.0,  Rs0_dT =0.0,        Rdvd_dT =0.0,       Rsvd_dT =0.0 ;
  double Vdse_eff =0.0, Vdse_eff_dVbse =0.0, Vdse_eff_dVdse =0.0, Vdse_eff_dVgse =0.0,
                        Vdse_eff_dVbs  =0.0, Vdse_eff_dVds  =0.0, Vdse_eff_dVgs  =0.0 ;
  double VdseModeNML =0.0, VdseModeRVS =0.0 ;
  double Vbsegmt =0.0, Vdsegmt =0.0, Vgsegmt =0.0 ;
  double Vbserev =0.0, Vdserev =0.0, Vgserev =0.0 ;
  double Ra_alpha,     Ra_beta ;

  /* modified external biases for symmetry */
  double /*Vzadd_ext = 0.0,*/ Vzadd_ext_dVd = 0.0 ;
  double Vdserevz = 0.0, Vdserevz_dVd = 0.0 ;
  double Vgserevz = 0.0, Vgserevz_dVd = 0.0 ;
  double Vbserevz = 0.0, Vbserevz_dVd = 0.0 ;

  /* Substrate Effect */
  const double RDVSUB =  model->HSMHV_rdvsub ; 
  const double RDVDSUB = model->HSMHV_rdvdsub ;
  const double DDRIFT =  model->HSMHV_ddrift ; 
  const double VBISUB =  model->HSMHV_vbisub ;
  const double NSUBSUB = modelMKS->HSMHV_nsubsub ;
  double Vsubsrev = 0.0 ;
  double Wdep = 0.0, Wdep_dVdserev = 0.0, Wdep_dVsubsrev = 0.0 ;
  double T1_dVdserev = 0.0, T1_dVsubsrev = 0.0, T6_dVdserev = 0.0, T6_dVsubsrev = 0.0 ;
  double Rs_dVsubs = 0.0, Rd_dVsubs = 0.0, Rdrift_dVsubs = 0.0, Rsdrift_dVsubs = 0.0 ;


  /* temperature-dependent variables for SHE model */
  double TTEMP =0.0, TTEMP0 =0.0 ;
  double/* Tdiff0 = 0.0, Tdiff0_2 = 0.0,*/ Tdiff = 0.0, Tdiff_2 = 0.0 ;
  double Eg =0.0,    Eg_dT =0.0 ;
  double Nin =0.0,   Nin_dT =0.0 ;
  double js =0.0,    js_dT =0.0 ;
  double jssw =0.0,  jssw_dT =0.0 ;
  double js2 =0.0,   js2_dT =0.0 ;
  double jssw2 =0.0, jssw2_dT =0.0 ;
  
  /* Qover 5/1 ckt-bias use */
  double Vgbgmt =0.0, Vgbgmt_dVbs =0.0, Vgbgmt_dVds =0.0, Vgbgmt_dVgs =0.0 ;
  double Vxbgmt =0.0, Vxbgmt_dVbs =0.0, Vxbgmt_dVds =0.0, Vxbgmt_dVgs =0.0 ;
  double Vxbgmtcl =0.0, Vxbgmtcl_dVxbgmt =0.0, Vxbgmtcl_dT =0.0 ;

  double ModeNML =0.0, ModeRVS =0.0 ;

  double QsuLD =0.0, QsuLD_dVbs =0.0, QsuLD_dVds =0.0, QsuLD_dVgs =0.0, QsuLD_dT =0.0 ;
  double QiuLD =0.0, QiuLD_dVbs =0.0, QiuLD_dVds =0.0, QiuLD_dVgs =0.0, QiuLD_dT =0.0 ;
  double /*QidLD =0.0,*/ QidLD_dVbs =0.0, QidLD_dVds =0.0, QidLD_dVgs =0.0, QidLD_dT =0.0 ;
  double /*QisLD =0.0,*/ QisLD_dVbs =0.0, QisLD_dVds =0.0, QisLD_dVgs =0.0, QisLD_dT =0.0 ;
  double /*QidLDext =0.0,*/ QidLDext_dVbse =0.0, QidLDext_dVdse =0.0, QidLDext_dVgse =0.0, QidLDext_dT =0.0 ;
  double /*QisLDext =0.0,*/ QisLDext_dVbse =0.0, QisLDext_dVdse =0.0, QisLDext_dVgse =0.0, QisLDext_dT =0.0 ;

  /* Self heating */
  double mphn0_dT =0.0 ;
  double ps0ldinib_dT =0.0, cnst0over_dT =0.0 ;
  double ps0ldinibs_dT =0.0, cnst0overs_dT =0.0 ;
  double Temp_dif =0.0 ;
  /* for SCE */
  double ptovr_dT =0.0 ;
 
  /* IBPC */
  double IdsIBPC =0.0, IdsIBPC_dVbs =0.0, IdsIBPC_dVds =0.0, IdsIBPC_dVgs =0.0, IdsIBPC_dT =0.0 ;

  /* Qover */
  int flg_ovzone = 0 ;
  double VgpLD =0.0,     VgpLD_dVgb =0.0 ;
  double /*VthLD =0.0,*/ Vgb_fb_LD =0.0 ;
  double Ac31_dVgb =0.0, Ac31_dVxb =0.0 ;
  double Ac1_dVgb =0.0, Ac1_dVxb =0.0 ;
  double Ac2_dVgb =0.0, Ac2_dVxb =0.0 ;
  double Ac3_dVgb =0.0, Ac3_dVxb =0.0 ;
  double Acn_dVgb =0.0, Acn_dVxb =0.0 ;
  double Acd_dVgb =0.0, Acd_dVxb =0.0 ;
  double Chi_dVgb =0.0, Chi_dVxb =0.0 ;
  double Psa_dVgb =0.0, Psa_dVxb =0.0 ;
  double QsuLD_dVgb =0.0, QsuLD_dVxb =0.0 ;
  double QbuLD_dVgb =0.0, QbuLD_dVxb =0.0 ;
  double fs02_dVgb =0.0 ;/*, fs02_dVxb =0.0 ;*/
  double TX_dVgb =0.0, TX_dVxb =0.0 ;
  double TY_dVgb =0.0, TY_dVxb =0.0 ;
  double Ps0LD =0.0,   Ps0LD_dVgb =0.0, Ps0LD_dVxb =0.0, Ps0LD_dT =0.0 ;
  double /*Ps0LD_dVbs =0.0,*/ Ps0LD_dVds =0.0; /*Ps0LD_dVgs =0.0 ;*/
  double Pb2over =0.0,                                   Pb2over_dT =0.0 ;

  int flg_overgiven =0 ;
  int Coovlps =0,   Coovlpd =0 ;
  double Lovers =0.0, Loverd =0.0 ;
  double Novers =0.0, Noverd =0.0 ;
  double Nover_func =0.0 ;
/*  double ps0ldinib_func =0.0, ps0ldinib_func_dT =0.0 ;*/
  double cnst0over_func =0.0, cnst0over_func_dT =0.0 ;
  double cnst1over =0.0, cnst1over_dT =0.0;
  /* Qover Analytical Model */
  int lp_ld;
  double Ta = 9.3868e-3, Tb = -0.1047839 ;
  double Tc,                   Tc_dT ; 
  double Td, Td_dVxb, Td_dVgb, Td_dT ;
  double Tv, Tv_dVxb, Tv_dVgb, Tv_dT ;
  double Tu, Tu_dVxb, Tu_dVgb, Tu_dT ;
  double Tp,                   Tp_dT ; 
  double Tq, Tq_dVxb, Tq_dVgb, Tq_dT ;
  double     T1_dVxb, T1_dVgb ;
/*  double     T2_dVxb, T2_dVgb ;*/
/*  double     T3_dVxb, T3_dVgb ;*/
  double     T5_dVxb, T5_dVgb ;
  double VgpLD_shift, VgpLD_shift_dT ;
  double VgpLD_shift_dVgb, VgpLD_shift_dVxb, exp_bVbs_dVxb ;
  double gamma, gamma_dVxb, gamma_dT ;
  double psi  , psi_dVgb  , psi_dVxb  , psi_dT ;
/*  double psi_B, arg_B ;*/
  double Chi_1, Chi_1_dVgb, Chi_1_dVxb ,Chi_1_dT ;
  double Chi_A, Chi_A_dVgb, Chi_A_dVxb, Chi_A_dT ;
  double Chi_B, Chi_B_dVgb, Chi_B_dVxb, Chi_B_dT;/*, Chi_B_dpsi , Chi_B_dgamma ;*/

  /* X_dT for leakage currents & junction diodes */
  double isbd_dT =0.0,      isbs_dT =0.0 ;
  double isbd2_dT =0.0,     isbs2_dT =0.0 ;
  double vbdt_dT =0.0,      vbst_dT = 0.0 ;
  double jd_expcd_dT =0.0 , jd_expcs_dT =0.0 ;
  double jd_nvtm_inv_dT =0.0 ;
  double exptemp_dT = 0.0 ;
  double tcjbd =0.0,    tcjbs =0.0,
         tcjbdsw =0.0,  tcjbssw =0.0,
         tcjbdswg =0.0, tcjbsswg =0.0 ;

  /*================ Start of executable code.=================*/


  if (here->HSMHV_mode == HiSIM_NORMAL_MODE) {
    ModeNML = 1.0 ;
    ModeRVS = 0.0 ;
  } else {
    ModeNML = 0.0 ;
    ModeRVS = 1.0 ;
  }

  T1 = Vdse + Vgse + Vbse + Vds + Vgs + Vbs + vbd_jct + vbs_jct ;
  if ( ! finite (T1) ) {
    fprintf (stderr ,
       "*** warning(HiSIM_HV): Unacceptable Bias(es).\n" ) ;
    fprintf (stderr , "----- bias information (HiSIM_HV)\n" ) ;
    fprintf (stderr , "name: %s\n" , here->HSMHVname ) ;
    fprintf (stderr , "states: %d\n" , here->HSMHVstates ) ;
    fprintf (stderr , "Vdse= %.3e Vgse=%.3e Vbse=%.3e\n"
            , Vdse , Vgse , Vbse ) ;
    fprintf (stderr , "Vdsi= %.3e Vgsi=%.3e Vbsi=%.3e\n"
            , Vds , Vgs , Vbs ) ;
    fprintf (stderr , "vbs_jct= %12.5e vbd_jct= %12.5e\n"
            , vbs_jct , vbd_jct ) ;
    fprintf (stderr , "vd= %.3e vs= %.3e vdp= %.3e vgp= %.3e vbp= %.3e vsp= %.3e\n" 
            , *( ckt->CKTrhsOld + here->HSMHVdNode ) 
            , *( ckt->CKTrhsOld + here->HSMHVsNode ) 
            , *( ckt->CKTrhsOld + here->HSMHVdNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVgNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVbNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVsNodePrime )  ) ;
    fprintf (stderr , "----- bias information (end)\n" ) ;
    return ( HiSIM_ERROR ) ;
  }

  flg_info = model->HSMHV_info ;
  flg_nqs = model->HSMHV_conqs ;
  
  /*-----------------------------------------------------------*
   * Start of the routine. (label)
   *-----------------*/
/* start_of_routine: */

  /*-----------------------------------------------------------*
   * Temperature dependent constants. 
   *-----------------*/
  if ( here->HSMHVtempNode > 0 && pParam->HSMHV_rth0 != 0.0 ) {

#define HSMHVEVAL
#include "hsmhvtemp_eval.h"

  } else {
    beta = here->HSMHV_beta ;
    TTEMP       = ckt->CKTtemp ;
    if ( here->HSMHV_dtemp_Given ) { TTEMP = TTEMP + here->HSMHV_dtemp ; }
    Eg_dT       = 0.0 ;
    beta_dT     = 0.0 ;
    beta_inv_dT = 0.0 ;
    mphn0_dT    = 0.0 ;
    ptovr_dT    = 0.0 ;
    Vmax_dT     = 0.0 ;
    Pb2_dT      = 0.0 ;
    cnst0_dT    = 0.0 ;
    cnst1_dT    = 0.0 ;
  }

  /* Inverse of the thermal voltage */
  beta_inv = here->HSMHV_beta_inv ;
  beta2 = here->HSMHV_beta2 ;

  /* Bandgap */
  Egp12 = here->HSMHV_egp12 ;
  Egp32 = here->HSMHV_egp32 ;

  /* Metallurgical channel geometry */
  Leff = here->HSMHV_leff ;
  Weff = here->HSMHV_weff ;
  WeffLD_nf = here->HSMHV_weff_ld * here->HSMHV_nf ;

  Ldrift0 = here->HSMHV_ldrift1 +  here->HSMHV_ldrift2 ;
  Ldrift = (model->HSMHV_coldrift) ? Ldrift0 
                                   : Ldrift0 + here->HSMHV_loverld ;

  /* Flat band voltage */
  Vfb = pParam->HSMHV_vfbc ;

  /* Surface impurity profile */
  q_Nsub = here->HSMHV_qnsub ;
  
  /* Velocity Temperature Dependence */
  Vmax = here->HSMHV_vmax ;
   
  /* 2 phi_B */
  Pb2 = here->HSMHV_pb2 ;
  Pb20 = here->HSMHV_pb20 ; 
  Pb2c = here->HSMHV_pb2c ;

  /* Coefficient of the F function for bulk charge */
  cnst0 = here->HSMHV_cnst0 ;

  /* cnst1: n_{p0} / p_{p0} */
  cnst1 = here->HSMHV_cnst1 ;

  /* c_eox: Permitivity in ox  */
  c_eox = here->HSMHV_cecox ;

  /* Tox and Cox without QME */
  Tox0 = model->HSMHV_tox ;
  Cox0 = c_eox / Tox0 ;
  Cox0_inv = 1.0 / Cox0 ;

  /*---------------------------------------------------*
   * Determine clamping limits for too large Vbs (internal).
   *-----------------*/

  Fn_SU( T1 , Pb2  - model->HSMHV_vzadd0 , Vbs_max , 0.1 , T0 ) ;
  Vbs_max = T1 ;
  Vbs_max_dT = Pb2_dT * T0 ;
  if ( Pb20 - model->HSMHV_vzadd0 < Vbs_max ) {
    Vbs_max = Pb20 - model->HSMHV_vzadd0 ;
    Vbs_max_dT = 0.0 ;
  }
  if ( Pb2c - model->HSMHV_vzadd0 < Vbs_max ) {
    Vbs_max = Pb2c - model->HSMHV_vzadd0 ;
    Vbs_max_dT = 0.0 ;
  }

  if ( Vbs_bnd > Vbs_max * 0.5 ) {
    Vbs_bnd = 0.5 * Vbs_max ;
    Vbs_bnd_dT = 0.5 * Vbs_max_dT ;
  }


  if (here->HSMHV_rs > 0.0 || here->HSMHV_rd > 0.0) {
    if ( model->HSMHV_corsrd == 1 ) flg_rsrd  = 1 ;
    if ( model->HSMHV_corsrd == 2 ) flg_rsrd  = 2 ;
    if ( model->HSMHV_corsrd == 3 ) flg_rsrd  = 3 ;
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * PART-1: Basic device characteristics. 
   *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
   * Prepare for potential initial guesses using previous values
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  flg_pprv = 0 ;

  if ( here->HSMHV_called >= 1 ) {

    Vbsc_dif = Vbs - here->HSMHV_vbsc_prv ;
    Vdsc_dif = Vds - here->HSMHV_vdsc_prv ;
    Vgsc_dif = Vgs- here->HSMHV_vgsc_prv ;

    sum_vdif  = fabs( Vbsc_dif ) + fabs( Vdsc_dif ) 
              + fabs( Vgsc_dif ) ;

    if ( model->HSMHV_copprv >= 1 && sum_vdif <= vtol_pprv  &&
         here->HSMHV_mode * here->HSMHV_mode_prv > 0 ) { flg_pprv = 1 ;}

    if ( here->HSMHV_called >= 2 && flg_pprv == 1 ) {
      Vbsc_dif2 = here->HSMHV_vbsc_prv - here->HSMHV_vbsc_prv2 ;
      Vdsc_dif2 = here->HSMHV_vdsc_prv - here->HSMHV_vdsc_prv2 ;
      Vgsc_dif2 = here->HSMHV_vgsc_prv - here->HSMHV_vgsc_prv2 ;
      sum_vdif2 = fabs( Vbsc_dif2 ) + fabs( Vdsc_dif2 ) 
                  + fabs( Vgsc_dif2 ) ;
      if ( epsm10 < sum_vdif2 && sum_vdif2 <= vtol_pprv &&
           here->HSMHV_mode_prv * here->HSMHV_mode_prv2 > 0 ) { flg_pprv = 2 ; }
    }
    Temp_dif = TTEMP - here->HSMHV_temp_prv ;

  } else {

    Vbsc_dif = 0.0 ;
    Vdsc_dif = 0.0 ;
    Vgsc_dif = 0.0 ;
    sum_vdif = 0.0 ;
    Vbsc_dif2 = 0.0 ;
    Vdsc_dif2 = 0.0 ;
    Vgsc_dif2 = 0.0 ;
    sum_vdif2 = 0.0 ;
    flg_pprv = 0 ;
    Temp_dif = 0.0 ;
  }

  dVbs = Vbsc_dif ;
  dVds = Vdsc_dif ;
  dVgs = Vgsc_dif ;

  if ( flg_pprv >= 1 ) {
    Ps0 = here->HSMHV_ps0_prv ;
    Ps0_dVbs = here->HSMHV_ps0_dvbs_prv ;
    Ps0_dVds = here->HSMHV_ps0_dvds_prv ;
    Ps0_dVgs = here->HSMHV_ps0_dvgs_prv ;
  
    Pds = here->HSMHV_pds_prv ;
    Pds_dVbs = here->HSMHV_pds_dvbs_prv ;
    Pds_dVds = here->HSMHV_pds_dvds_prv ;
    Pds_dVgs = here->HSMHV_pds_dvgs_prv ;
  }


  if ( flg_rsrd == 1 || flg_rsrd == 3 ) {

    /*----------------------------------------------------------*
     * Considering these special cases:
     * ( here->HSMHV_mode == HiSIM_NORMAL_MODE  && Vdse < 0.0 )
     * ( here->HSMHV_mode == HiSIM_REVERSE_MODE && Vdse < 0.0 )
     *----------------------------------------------------------*/
    Vdsegmt = here->HSMHV_mode * Vdse ; /* geometrical outer bias */
    Vgsegmt = Vgse - ModeRVS * Vdse ;   /* geometrical outer bias */
    Vbsegmt = Vbse - ModeRVS * Vdse ;   /* geometrical outer bias */
    if ( Vdsegmt >= 0.0 ) { /* vdse normal mode */
      VdseModeNML = 1 ;
      VdseModeRVS = 0 ;
      Vdserev = Vdsegmt ;
      Vgserev = Vgsegmt ;
      Vbserev = Vbsegmt ;
      Vsubsrev = Vsubs ;
    } else { /* vdse reverse mode */
      VdseModeNML = 0 ;
      VdseModeRVS = 1 ;
      Vdserev = - Vdsegmt ;
      Vgserev = Vgsegmt - Vdsegmt ;
      Vbserev = Vbsegmt - Vdsegmt ;
      Vsubsrev = Vsubs - Vdsegmt ;
    }


    if ( here->HSMHV_rdvd > 0.0 || here->HSMHV_rsvd > 0.0 || pParam->HSMHV_rdvg11 > 0.0 || pParam->HSMHV_rdvb > 0.0 || here->HSMHVsubNode >= 0 ) {
      /*-----------------------------------------------------------*
       * Vxserevz: Modified bias introduced to realize symmetry at Vds=0.
       *-----------------*/
      Fn_SymAdd( Vzadd , Vdserev / 2 , model->HSMHV_vzadd0 , T2 ) ;
      Vzadd_ext_dVd = T2 / 2 ;
      if ( Vzadd < ps_conv ) {
	Vzadd = ps_conv ;
	Vzadd_ext_dVd = 0.0 ;
      }
      Vdserevz = Vdserev + 2.0 * Vzadd ;
      Vdserevz_dVd = 1.0 + 2.0 * Vzadd_ext_dVd ;
      Vgserevz = Vgserev + Vzadd ;
      Vgserevz_dVd = Vzadd_ext_dVd ;
      Vbserevz = Vbserev + Vzadd ;
      Vbserevz_dVd = Vzadd_ext_dVd ;


      /* bias-dependent Rdrift for HVMOS/LDMOS */

      if ( model->HSMHV_cosym == 1 || VdseModeNML == 1 ) { /* HVMOS or normal mode LDMOS: */
      /* ... Vdse dependence             */
      T1     = VdseModeNML * here->HSMHV_rd + VdseModeRVS * here->HSMHV_rs ;
      T1_dT  = VdseModeNML * Rd0_dT         + VdseModeRVS * Rs0_dT ;
      T0     = VdseModeNML * here->HSMHV_rdvd + VdseModeRVS * here->HSMHV_rsvd ;
      T0_dT  = VdseModeNML * Rdvd_dT          + VdseModeRVS * Rsvd_dT ;
      T4     = T1 + T0 * Vdserevz ;
      T4_dVd = T0 * Vdserevz_dVd ;
      T4_dT  = T1_dT + T0_dT * Vdserevz ;

      /* ... Vgse dependence             */
      T10 = model->HSMHV_rdvg12 + small ;
      T1     =   T4     * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) ) ;
      T1_dVd =   T4_dVd * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) )
	       + T4     *         pParam->HSMHV_rdvg11 * (     - Vgserevz_dVd / T10 ) ;
      T1_dVg =   T4     *         pParam->HSMHV_rdvg11 * (     - 1.0     ) / T10 ;
      T1_dT  =   T4_dT  * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) ) ;
      Fn_SL2( T2 , T1 , T4 , rdsl2_dlt , T0 , T5 ) ;
      T2_dVd = T0 * T1_dVd + T5 * T4_dVd ;
      T2_dVg = T0 * T1_dVg ;
      T2_dT  = T0 * T1_dT  + T5 * T4_dT ;

      T3     = T4     * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
      T3_dVd = T4_dVd * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
      T3_dT  = T4_dT  * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
      Fn_SU2( Rdrift , T2 , T3 , rdsu2_dlt , T0, T5 ) ;
      Rdrift_dVdse = T0 * T2_dVd + T5 * T3_dVd ;
      Rdrift_dVgse = T0 * T2_dVg ;
      Rdrift_dT  = T0 * T2_dT  + T5 * T3_dT ;

      /* ... Vbse dependence             */
      T1 = 1.0 - pParam->HSMHV_rdvb * Vbserevz ;
      T1_dVb = - pParam->HSMHV_rdvb ;
      T1_dVd = - pParam->HSMHV_rdvb * Vbserevz_dVd ;
      Fn_SZ( T3 , T1 , rdsz_dlt , T4 ) ;
      T3_dVb = T4 * T1_dVb ;
      T3_dVd = T4 * T1_dVd ;
      T0           = Rdrift ;
      Rdrift       = Rdrift       * T3 ;
      Rdrift_dVdse = Rdrift_dVdse * T3 + T0 * T3_dVd ;
      Rdrift_dVgse = Rdrift_dVgse * T3 ;
      Rdrift_dVbse =                   + T0 * T3_dVb ;
      Rdrift_dT    = Rdrift_dT    * T3 ;

      } else { /* reverse mode LDMOS: */
	Rdrift =       here->HSMHV_rs ;
	Rdrift_dVdse = 0.0 ;
	Rdrift_dVgse = 0.0 ;
	Rdrift_dVbse = 0.0 ;
	Rdrift_dT =    Rs0_dT ;
      }
     
      
      /* Rsdrift */
      T4    = ( VdseModeNML * here->HSMHV_rs + VdseModeRVS * here->HSMHV_rd ) ;
      T4_dT = VdseModeNML * Rs0_dT         + VdseModeRVS * Rd0_dT ;
      T4_dVd = 0.0 ;

      if ( model->HSMHV_cosym == 1 || VdseModeRVS == 1 ) { /* HVMOS or reverse mode LDMOS: */
	/* ... Vdse dependence             */
	T0     = VdseModeNML * here->HSMHV_rsvd + VdseModeRVS * here->HSMHV_rdvd ;
	T0_dT  = VdseModeNML * Rsvd_dT          + VdseModeRVS * Rdvd_dT ;
/* 	if ( model->HSMHV_cosym == 2 ) { /\* latest case with bugfix: *\/ */
/* 	T4     = T4 + T0 * Vdserevz ; */
/* 	T4_dVd =      T0 * Vdserevz_dVd ; */
/* 	T4_dT  = T4_dT + T0_dT * Vdserevz ; */
/* 	} else { /\* HiSIM_HV 1.1.1 compatible case *\/ */
	  T4 = T4 + T0 * ( 2.0 * model->HSMHV_vzadd0 ) ; /* 2.0 * Fn_SymAdd( x=0, add0=model->HSMHV_vzadd0 ) */
	  T4_dT = T4_dT + T0_dT * ( 2.0 * model->HSMHV_vzadd0 ) ;
/* 	} */

	/* ... Vgse dependence             */
        T10 = model->HSMHV_rdvg12 + small ;
	T1     =   T4     * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) ) ;
	T1_dVd =   T4_dVd * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) )
	         + T4     *         pParam->HSMHV_rdvg11 * (     - Vgserevz_dVd / T10 ) ;
	T1_dVg =   T4     *         pParam->HSMHV_rdvg11 * (     - 1.0     ) / T10 ;
	T1_dT  =   T4_dT  * ( 1.0 + pParam->HSMHV_rdvg11 * ( 1.0 - Vgserevz   / T10 ) ) ;
	Fn_SL2( T2 , T1 , T4 , rdsl2_dlt , T0 , T5 ) ;
	T2_dVd = T0 * T1_dVd + T5 * T4_dVd ;
	T2_dVg = T0 * T1_dVg ;
	T2_dT  = T0 * T1_dT  + T5 * T4_dT ;

	T3     = T4     * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
	T3_dVd = T4_dVd * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
	T3_dT  = T4_dT  * ( 1.0 + pParam->HSMHV_rdvg11 ) ;
       Fn_SU2( Rsdrift , T2 , T3 , rdsu2_dlt , T0, T5 ) ;
	Rsdrift_dVdse = T0 * T2_dVd + T5 * T3_dVd ;
	Rsdrift_dVgse = T0 * T2_dVg ;
	Rsdrift_dT  = T0 * T2_dT  + T5 * T3_dT ;

	/* ... Vbse dependence             */
	T1 = 1.0 - pParam->HSMHV_rdvb * Vbserevz ;
	T1_dVb = - pParam->HSMHV_rdvb ;
	T1_dVd = - pParam->HSMHV_rdvb * Vbserevz_dVd ;
	Fn_SZ( T3 , T1 , rdsz_dlt , T4 ) ;
	T3_dVb = T4 * T1_dVb ;
	T3_dVd = T4 * T1_dVd ;
	T0            = Rsdrift ;
	Rsdrift       = Rsdrift       * T3 ;
	Rsdrift_dVdse = Rsdrift_dVdse * T3 + T0 * T3_dVd ;
	Rsdrift_dVgse = Rsdrift_dVgse * T3 ;
	Rsdrift_dVbse =                   + T0 * T3_dVb ;
	Rsdrift_dT    = Rsdrift_dT    * T3 ;
      } else { /* LDMOS normal mode: */
	  Rsdrift =       here->HSMHV_rs ;
	  Rsdrift_dVdse = 0.0 ;
	  Rsdrift_dVgse = 0.0 ;
	  Rsdrift_dVbse = 0.0 ;
	  Rsdrift_dT =    Rs0_dT ;
      }
      
      
      if ( here->HSMHVsubNode >= 0 && model->HSMHV_cosym == 0 && 
         ( pParam->HSMHV_nover * ( NSUBSUB + pParam->HSMHV_nover ) ) > 0 ) {
       /* external substrate node exists && LDMOS case: */
	/* Substrate Effect */
	T0 = VBISUB - RDVDSUB * Vdserevz - RDVSUB * Vsubsrev ;
	Fn_SZ( T1, T0, 10.0, T2 ) ;
	T1_dVdserev  = - RDVDSUB * Vdserevz_dVd * T2 ;  
	T1_dVsubsrev = - RDVSUB  * T2 ;

	T0 = NSUBSUB / ( pParam->HSMHV_nover * ( NSUBSUB + pParam->HSMHV_nover ) ) ;
	T4 = 2 * C_ESI / C_QE * T0 ;
	Wdep = sqrt ( T4 * T1 ) + small ;
	Wdep_dVdserev =  0.5 * T4 * T1_dVdserev  / Wdep ;
	Wdep_dVsubsrev = 0.5 * T4 * T1_dVsubsrev / Wdep ;

	Fn_SU( Wdep, Wdep, DDRIFT, C_sub_delta * DDRIFT, T0 ) ;
	Wdep_dVdserev *=  T0 ;
	Wdep_dVsubsrev *= T0 ;

	T0 = DDRIFT - Wdep ;
	Fn_SZ( T0, T0, C_sub_delta2, T2 ) ;
	T6 = Ldrift0 / T0 ;
	T6_dVdserev =  T2 * Wdep_dVdserev  * T6 / T0 ; 
	T6_dVsubsrev = T2 * Wdep_dVsubsrev * T6 / T0 ;

	if ( VdseModeNML == 1 ) { /* Vdse normal mode: */
	  T0 = Rdrift ;
	  Rdrift       = T0           * T6 ;
	  Rdrift_dVdse = Rdrift_dVdse * T6  +  T0 * T6_dVdserev ;
	  Rdrift_dVgse = Rdrift_dVgse * T6 ;
	  Rdrift_dVbse = Rdrift_dVbse * T6 ;
	  Rdrift_dVsubs=                       T0 * T6_dVsubsrev ;
	  Rdrift_dT    = Rdrift_dT    * T6 ;
	} else { /* Vdse reverse mode: */
	  T0 = Rsdrift ;
	  Rsdrift =       T0            * T6 ;
	  Rsdrift_dVdse = Rsdrift_dVdse * T6  +  T0 * T6_dVdserev ;
	  Rsdrift_dVgse = Rsdrift_dVgse * T6 ;
	  Rsdrift_dVbse = Rsdrift_dVbse * T6 ;
	  Rsdrift_dVsubs =                       T0 * T6_dVsubsrev ;
	  Rsdrift_dT =    Rsdrift_dT    * T6 ;
	}
      }

      Rd = Rdrift ;
      Rd_dVgse = Rdrift_dVgse ;
      Rd_dVdse = Rdrift_dVdse ;
      Rd_dVbse = Rdrift_dVbse ;
      Rd_dVsubs = Rdrift_dVsubs ;
      Rd_dT    = Rdrift_dT ; 
      Rs = Rsdrift ;
      Rs_dVgse = Rsdrift_dVgse ;
      Rs_dVdse = Rsdrift_dVdse ;
      Rs_dVbse = Rsdrift_dVbse ;
      Rs_dVsubs = Rsdrift_dVsubs ;
      Rs_dT    = Rsdrift_dT ; 

    } else { /* bias-independent Rs/Rd */
      Rd = VdseModeNML * here->HSMHV_rd + VdseModeRVS * here->HSMHV_rs ;
      Rd_dT = VdseModeNML * Rd0_dT      + VdseModeRVS * Rs0_dT ;
      Rs = VdseModeNML * here->HSMHV_rs + VdseModeRVS * here->HSMHV_rd ;
      Rs_dT = VdseModeNML * Rs0_dT      + VdseModeRVS * Rd0_dT ;
    }

    /* Weff dependence of the resistances */
    Rd = Rd  /  WeffLD_nf ;
    Rd_dVgse /= WeffLD_nf ;
    Rd_dVdse /= WeffLD_nf ;
    Rd_dVbse /= WeffLD_nf ;
    Rd_dVsubs /= WeffLD_nf ;
    Rd_dT    /= WeffLD_nf ;
    Rs =  Rs /  WeffLD_nf ;
    Rs_dVgse /= WeffLD_nf ;
    Rs_dVdse /= WeffLD_nf ;
    Rs_dVbse /= WeffLD_nf ;
    Rs_dVsubs /= WeffLD_nf ;
    Rs_dT    /= WeffLD_nf ;

    /* Sheet resistances are added. */
    Rd += VdseModeNML * here->HSMHV_rd0 + VdseModeRVS * here->HSMHV_rs0 ;
    Rs += VdseModeNML * here->HSMHV_rs0 + VdseModeRVS * here->HSMHV_rd0 ;

    /* Re-stamps for hsmhvnoi.c */
    /* Please see hsmhvnoi.c */
    T0 = VdseModeNML * Rd + VdseModeRVS * Rs ; /* mode-dependent --> geometrical */
    if ( T0 > 0.0 && model->HSMHV_cothrml != 0 ) here->HSMHVdrainConductance = Mfactor / T0 ;
    else here->HSMHVdrainConductance = 0.0 ;
    T0 = VdseModeNML * Rs + VdseModeRVS * Rd ; /* mode-dependent --> geometrical */
    if ( T0 > 0.0 && model->HSMHV_cothrml != 0 ) here->HSMHVsourceConductance = Mfactor / T0 ;
    else here->HSMHVsourceConductance = 0.0 ;

  } /* end of case flg_rsrd=1 or flg_rsrd=3 */




    /* Clamping for Vbs > Vbs_bnd */
    if ( Vbs > Vbs_bnd ) {
      T1 = Vbs - Vbs_bnd ;
      T2 = Vbs_max - Vbs_bnd ;
      T1_dT = - Vbs_bnd_dT ;
      T2_dT = Vbs_max_dT - Vbs_bnd_dT ;

      Fn_SUPoly4m( TY , T1 , T2 , Vbscl_dVbs , T0 ) ;
      TY_dT = T1_dT * Vbscl_dVbs + T2_dT * T0 ;

      Vbscl    = Vbs_bnd    + TY ;
      Vbscl_dT = Vbs_bnd_dT + TY_dT ;

      T3 = 1 / T2 ;

      /* x/xmax */
      T4 = T1 * T3 ;
      T4_dVb = T3 ;
      T4_dT = T1_dT * T3 - T1*T3*T3*T2_dT;
      T4_dVb_dT = -T3*T3*T2_dT ;
 
      T5 = T4 * T4;
      T5_dVb = 2 * T4_dVb * T4 ;
      T5_dT = 2.0*T4*T4_dT;
      T5_dVb_dT = 2 * T4_dVb_dT * T4 + 2 * T4_dVb * T4_dT ;
      T15 = 2 * T4_dVb * T4_dVb ; /* T15 = T5_dVb_dVb */
 
      T6 = T4 * T5 ;
      T6_dVb = T4_dVb * T5 + T4 * T5_dVb ;
      T6_dT = T4_dT * T5 + T4 * T5_dT ;
      T6_dVb_dT = T4_dVb_dT * T5 + T4_dVb * T5_dT + T4_dT * T5_dVb + T4*T5_dVb_dT ;
      T16 = T4_dVb * T5_dVb + T4_dVb * T5_dVb + T4 * T15 ; /* T16 = T6_dVb_dVb */
 
      /* T7 = Z  T7_dVb = dZ_dVb  T17 = dZ_dVb_dVb */
      T7 = 1 + T4 + T5 + T6 + T5 * T5 ;
      T7_dVb = T4_dVb + T5_dVb + T6_dVb + 2 * T5_dVb * T5 ;
      T7_dT = T4_dT + T5_dT + T6_dT + 2 * T5_dT * T5 ;
      T7_dVb_dT = T4_dVb_dT + T5_dVb_dT + T6_dVb_dT + 2 * T5_dVb_dT * T5 + 2 * T5_dVb * T5_dT ;
      T17 = T15 + T16 + 2 * T15 * T5 + 2 * T5_dVb * T5_dVb ;
 
      T8 = T7 * T7 ;
      T8_dVb = 2 * T7_dVb * T7 ;
      T8_dT = 2 * T7_dT * T7 ;
 
      T9 = 1 / T8 ;
      T9_dVb = - T8_dVb * T9 * T9 ;
      T9_dT = - T8_dT * T9 * T9 ;
 
      Vbscl_dVbs = T2 * T7_dVb * T9 ;
      Vbscl_dVbs_dT = T2_dT * T7_dVb * T9 + T2*(T7_dVb_dT * T9+ T7_dVb * T9_dT);
      Vbscl_dVbs_dVbs = T2 * ( T17 * T9 + T7_dVb * T9_dVb ) ;
    }  else {
      Vbscl      = Vbs ;
      Vbscl_dVbs = 1.0 ;
      Vbscl_dT   = 0.0 ;
      Vbscl_dVbs_dVbs = 0.0 ;
    }

    /*-----------------------------------------------------------*
     * Vxsz: Modified bias introduced to realize symmetry at Vds=0.
     *-----------------*/


    T1 = Vbscl_dVbs * Vds / 2 ;
    Fn_SymAdd(  Vzadd , T1 , model->HSMHV_vzadd0 , T2 ) ;
    Vzadd_dVbs = T2 * Vbscl_dVbs_dVbs * Vds / 2 ;
    Vzadd_dT = T2 * Vbscl_dVbs_dT * Vds / 2 ;
    T2 *= Vbscl_dVbs / 2 ;
    Vzadd_dVds = T2 ;

    if ( Vzadd < ps_conv ) {
      Vzadd = ps_conv ;
      Vzadd_dVds = 0.0 ;
      Vzadd_dVbs = 0.0 ;
      Vzadd_dT = 0.0 ;
    }
 
    Vbsz      = Vbscl + Vzadd ;
    Vbsz_dVbs = Vbscl_dVbs + Vzadd_dVbs ;
    Vbsz_dVds = Vzadd_dVds ;
    Vbsz_dT   = Vbscl_dT + Vzadd_dT;
 
    Vdsz = Vds + 2.0 * Vzadd ;
    Vdsz_dVbs = 2.0 * Vzadd_dVbs ;
    Vdsz_dVds = 1.0 + 2.0 * Vzadd_dVds ;
    Vdsz_dT = 2.0 * Vzadd_dT ;
 
    Vgsz = Vgs + Vzadd ;
    Vgsz_dVbs = Vzadd_dVbs ;
    Vgsz_dVgs = 1.0 ;
    Vgsz_dVds = Vzadd_dVds ;
    Vgsz_dT = Vzadd_dT ;

    /*---------------------------------------------------*
     * Factor of modification for symmetry.
     *-----------------*/

    T1 = here->HSMHV_qnsub_esi * Cox0_inv * Cox0_inv ;
    T2 = Vgs - Vfb ;
    T3 = 1 + 2.0 / T1 * ( T2 - 1.0 / here->HSMHV_betatnom - Vbscl ) ;

    Fn_SZ( T4 , T3 , 1e-3 , T5 ) ;
    TX = sqrt( T4 ) ;
    Pslsat = T2 + T1 * ( 1.0 - TX ) ; 
    VdsatS = Pslsat - Pb2c ;
    Fn_SL( VdsatS , VdsatS , 0.1 , 5e-2 , T6 ) ;

    VdsatS_dVbs = ( TX ? (T6 * T5 / TX * Vbscl_dVbs) : 0.0 ) ;
    VdsatS_dVds = 0.0 ;
    VdsatS_dVgs = ( TX ? (T6 * ( 1.0 - T5 / TX )) : 0.0 ) ;
    VdsatS_dT = (TX ? (T6* T5/TX * Vbscl_dT) : 0) ;

    T1 = Vds / VdsatS ;
    Fn_SUPoly4( TX , T1 , 1.0 , T0 ) ; 
    FMDVDS = TX * TX ;
    T2 = 2 * TX * T0 ;
    T3 = T2 / ( VdsatS * VdsatS ) ;
    FMDVDS_dVbs = T3 * ( - Vds * VdsatS_dVbs ) ;
    FMDVDS_dVds = T3 * ( 1.0 * VdsatS - Vds * VdsatS_dVds ) ;
    FMDVDS_dVgs = T3 * ( - Vds * VdsatS_dVgs ) ;
    FMDVDS_dT = T3 * ( - Vds * VdsatS_dT ) ;

    /*-----------------------------------------------------------*
     * Quantum effect
     *-----------------*/
    if ( model->HSMHV_flg_qme == 0 ) {
      flg_qme = 0 ;
    } else {
      flg_qme = 1 ;
    }

    T1 = here->HSMHV_2qnsub_esi ;
    T2 = sqrt( T1 * Pb20 ) ;
    Vthq = Pb20 + Vfb + T2 * Cox0_inv ;
    Vthq_dVb = 0.0 ;
    Vthq_dVd = 0.0 ;

    if ( flg_qme == 0 ) {
      Tox = Tox0 ;
      Tox_dVb = 0.0 ;
      Tox_dVd = 0.0 ;
      Tox_dVg = 0.0 ;
  
      Cox = Cox0 ;
      Cox_dVb = 0.0 ;
      Cox_dVd = 0.0 ;
      Cox_dVg = 0.0 ;
  
      Cox_inv  = Cox0_inv ;
      Cox_inv_dVb = 0.0 ;
      Cox_inv_dVd = 0.0 ;
      Cox_inv_dVg = 0.0 ;
  
      T0 = cnst0 * cnst0 * Cox_inv ;
      cnstCoxi = T0 * Cox_inv ;
      cnstCoxi_dVb = 0.0 ;
      cnstCoxi_dVd = 0.0 ;
      cnstCoxi_dVg = 0.0 ;
      cnstCoxi_dT = 2.0 * cnst0 * cnst0_dT * Cox_inv * Cox_inv ;

    } else {

      T1 = - model->HSMHV_qme2 ;
      T5 = Vgsz - Vthq - model->HSMHV_qme2  ;
      T5_dVb = Vgsz_dVbs - Vthq_dVb ;
      T5_dVd = Vgsz_dVds - Vthq_dVd ;
      T5_dVg = Vgsz_dVgs ;
      T5_dT = Vgsz_dT ;
      Fn_SZ( T2 , - T5 , qme2_dlt, T3) ;
      T2 = T2 + small ;
      T2_dVb = - T3 * T5_dVb ;
      T2_dVd = - T3 * T5_dVd ;
      T2_dVg = - T3 * T5_dVg ;
      T2_dT = - T3 * T5_dT ;
      T3 = model->HSMHV_qme12 * T1 * T1 ;
      T4 = model->HSMHV_qme12 * T2 * T2 + model->HSMHV_qme3 ;
      Fn_SU( dTox , T4 , T3 , qme_dlt , T6 ) ;
      T7 = 2 * model->HSMHV_qme12 * T2 * T6 ;
      dTox_dVb = T7 * T2_dVb ;
      dTox_dVd = T7 * T2_dVd ;
      dTox_dVg = T7 * T2_dVg ;
      dTox_dT = T7 * T2_dT ;

  
      if ( dTox * 1.0e12 < Tox0 ) {
        dTox = 0.0 ;
        dTox_dVb = 0.0 ;
        dTox_dVd = 0.0 ;
        dTox_dVg = 0.0 ;
        dTox_dT = 0.0 ;
        flg_qme = 0 ;
      }
 
      Tox = Tox0 + dTox ;
      Tox_dVb = dTox_dVb ;
      Tox_dVd = dTox_dVd ;
      Tox_dVg = dTox_dVg ;
      Tox_dT = dTox_dT ;
 
      Cox = c_eox / Tox ;
      T1  = - c_eox / ( Tox * Tox ) ;
      Cox_dVb = T1 * Tox_dVb ;
      Cox_dVd = T1 * Tox_dVd ;
      Cox_dVg = T1 * Tox_dVg ;
      Cox_dT = T1 * Tox_dT ;
 
      Cox_inv  = Tox / c_eox ;
      T1  = 1.0 / c_eox ;
      Cox_inv_dVb = T1 * Tox_dVb ;
      Cox_inv_dVd = T1 * Tox_dVd ;
      Cox_inv_dVg = T1 * Tox_dVg ;
      Cox_inv_dT = T1 * Tox_dT ;
 
      T0 = cnst0 * cnst0 * Cox_inv ;
      cnstCoxi = T0 * Cox_inv ;
      T1 = 2.0 * T0 ;
      cnstCoxi_dVb = T1 * Cox_inv_dVb ;
      cnstCoxi_dVd = T1 * Cox_inv_dVd ;
      cnstCoxi_dVg = T1 * Cox_inv_dVg ;
      cnstCoxi_dT = 2.0 * cnst0 * cnst0_dT * Cox_inv * Cox_inv + T1 * Cox_inv_dT;
    }

    /*---------------------------------------------------*
     * Vbsz2 : Vbs for dVth
     *-----------------*/
    Vbsz2 = Vbsz ;
    Vbsz2_dVbs =  Vbsz_dVbs ;
    Vbsz2_dVds = Vbsz_dVds ;
    Vbsz2_dVgs = 0.0  ;
    Vbsz2_dT = Vbsz_dT ;
 
    /*---------------------------------------------------*
     * Vthp : Vth with pocket.
     *-----------------*/
    T1 = here->HSMHV_2qnsub_esi ;
    Qb0 = sqrt (T1 * (Pb20 - Vbsz2)) ;
    T2 = 0.5 * T1 / Qb0 ;
    Qb0_dVb = T2 * (- Vbsz2_dVbs) ;
    Qb0_dVd = T2 * (- Vbsz2_dVds) ;
    Qb0_dVg = T2 * (- Vbsz2_dVgs) ;
    Qb0_dT = T2 * (- Vbsz2_dT) ;
 
    Vthp = Pb20 + Vfb + Qb0 * Cox_inv + here->HSMHV_ptovr;
    Vthp_dVb = Qb0_dVb * Cox_inv + Qb0 * Cox_inv_dVb ;
    Vthp_dVd = Qb0_dVd * Cox_inv + Qb0 * Cox_inv_dVd ;
    Vthp_dVg = Qb0_dVg * Cox_inv + Qb0 * Cox_inv_dVg ;
    Vthp_dT = Qb0_dT * Cox_inv + Qb0 * Cox_inv_dT + ptovr_dT ;
    
    if ( pParam->HSMHV_pthrou != 0.0 ) {
      /* Modify Pb20 to Pb20b */
      T11 = beta * 0.25 ;
      T10 = beta_inv - cnstCoxi * T11 + small ;
      T10_dVg = - T11 * cnstCoxi_dVg ;
      T10_dVd = - T11 * cnstCoxi_dVd ;
      T10_dVb = - T11 * cnstCoxi_dVb ;
      T10_dT = beta_inv_dT - ( T11 * cnstCoxi_dT + beta_dT * 0.25 * cnstCoxi ) ;

      T1 = Vgsz - T10 - psia2_dlt ;
      T1_dVg = Vgsz_dVgs - T10_dVg ;
      T1_dVd = Vgsz_dVds - T10_dVd ;
      T1_dVb = Vgsz_dVbs - T10_dVb ;
      T1_dT = Vgsz_dT - T10_dT ;
      T0 = Fn_Sgn (T10) ;
      T2 = sqrt (T1 * T1 + T0 * 4.0 * T10 * psia2_dlt) ;
      T3 = T10 + 0.5 * (T1 + T2) - Vfb ; /* Vgpa for sqrt calc. */
      T4 = T1 / T2 ;
      T5 = T0 * 2.0 * psia2_dlt / T2 ;
      T3_dVg = T10_dVg
             + 0.5 * (T1_dVg
                   + (T4 * T1_dVg + T5 * T10_dVg ) ) ;
      T3_dVd = T10_dVd
             + 0.5 * (T1_dVd
                   + (T4 * T1_dVd + T5 * T10_dVd ) ) ;
      T3_dVb = T10_dVb
             + 0.5 * (T1_dVb
                   + (T4 * T1_dVb + T5 * T10_dVb ) ) ;
      T3_dT = T10_dT
             + 0.5 * (T1_dT
                   + (T4 * T1_dT + T5 * T10_dT ) ) ;
      T4 = 4.0 / cnstCoxi * beta_inv * beta_inv ; 
      T8 = 4.0 / cnstCoxi ;
      T9 =  beta_inv * beta_inv ;
      T4_dT = - 4.0 * cnstCoxi_dT / ( cnstCoxi * cnstCoxi ) * T9
            + T8 * 2.0 * beta_inv * beta_inv_dT ;
      T5 = beta * T3 - 1.0 ;
      T5_dT = beta_dT * T3 + beta * T3_dT ;
      T6 = T5 / cnstCoxi ;
      T1 = 1.0 + T5 * T4 ;
      T2 = beta * T4 ;
      T6 = T6 * T4 ;
      T1_dVg = (T2 * T3_dVg - T6 * cnstCoxi_dVg ) ;
      T1_dVd = (T2 * T3_dVd - T6 * cnstCoxi_dVd ) ;
      T1_dVb = (T2 * T3_dVb - T6 * cnstCoxi_dVb ) ; 
      T1_dT = T5_dT * T4 + T5 * T4_dT ;
      Fn_SZ( T1 ,T1, psia_dlt, T7) ;
          T1 += epsm10 ;

      T1_dVg *= T7 ;
      T1_dVd *= T7 ;
      T1_dVb *= T7 ;
      T1_dT *= T7 ;

      T2 = sqrt (T1) ;
      T5 = 0.5 / T2 ;
      T2_dVg = T1_dVg * T5 ;
      T2_dVd = T1_dVd * T5 ;
      T2_dVb = T1_dVb * T5 ;
      T2_dT = T1_dT * T5 ;

      T4 = 0.5 * beta ;
      Psi_a = T3 + cnstCoxi * T4 * (1.0 - T2) ;
      T5 = T4 *  (1.0 - T2) ;
      T6 = T4 * cnstCoxi ;
      Psi_a_dVg = T3_dVg
                + (cnstCoxi_dVg * T5 - T6 * T2_dVg) ;
      Psi_a_dVd = T3_dVd
                + (cnstCoxi_dVd * T5 - T6 * T2_dVd) ;
      Psi_a_dVb = T3_dVb
                + (cnstCoxi_dVb * T5 - T6 * T2_dVb) ;
      Psi_a_dT = T3_dT
                + (cnstCoxi_dT * T5 - T6 * T2_dT) 
                + cnstCoxi * 0.5 * beta_dT * ( 1.0 - T2 ) ;

      Fn_SU( Pb20a , Psi_a, Pb20, delta0, T2) ;
      Pb20a_dVb = T2 * Psi_a_dVb ;
      Pb20a_dVd = T2 * Psi_a_dVd ;
      Pb20a_dVg = T2 * Psi_a_dVg ;
      Pb20a_dT = T2 * Psi_a_dT ;
    }

    T1 = pParam->HSMHV_pthrou ;
    Pb20b = Pb20 + T1 * (Pb20a - Pb20) ;
    Pb20b_dVb = T1 * Pb20a_dVb ;
    Pb20b_dVd = T1 * Pb20a_dVd ;
    Pb20b_dVg = T1 * Pb20a_dVg ;
    Pb20b_dT = T1 * Pb20a_dT ;
    
    T0 = 0.95 ;
    T1 = T0 * Pb20b - Vbsz2 - 1.0e-3 ;
    T1_dVb = T0 * Pb20b_dVb - Vbsz2_dVbs ;
    T1_dVd = T0 * Pb20b_dVd - Vbsz2_dVds ;
    T1_dVg = T0 * Pb20b_dVg - Vbsz2_dVgs ;
    T1_dT = T0 * Pb20b_dT - Vbsz2_dT ;
    T2 = sqrt (T1 * T1 + 4.0 * T0 * Pb20b * 1.0e-3) ;
    T3 = T0 * Pb20b - 0.5 * (T1 + T2) ;
    T4 = 2.0 * T0 * 1.0e-3 ;
    T5 = T1 / T2 ;
    T6 = T4 / T2 ;
    T3_dVb = T0 * Pb20b_dVb
           - 0.5 * (T1_dVb + (T1_dVb * T5 + T6 * Pb20b_dVb ) ) ;
    T3_dVd = T0 * Pb20b_dVd
           - 0.5 * (T1_dVd + (T1_dVd * T5 + T6 * Pb20b_dVd ) ) ;
    T3_dVg = T0 * Pb20b_dVg
           - 0.5 * (T1_dVg + (T1_dVg * T5 + T6 * Pb20b_dVg ) ) ;
    T3_dT = T0 * Pb20b_dT
           - 0.5 * (T1_dT + (T1_dT * T5 + T6 * Pb20b_dT ) ) ;
    Pbsum = Pb20b - T3 ;
    Pbsum_dVb = Pb20b_dVb - T3_dVb ;
    Pbsum_dVd = Pb20b_dVd - T3_dVd ;
    Pbsum_dVg = Pb20b_dVg - T3_dVg ;
    Pbsum_dT = Pb20b_dT - T3_dT ;

    sqrt_Pbsum = sqrt( Pbsum ) ;

    /*-------------------------------------------*
     * dVthLP : Short-channel effect induced by pocket.
     * - Vth0 : Vth without pocket.
     *-----------------*/
    if ( model->HSMHV_lp != 0.0 ) {
      T1 = here->HSMHV_2qnsub_esi ;
      T2 = model->HSMHV_bs2 - Vbsz2 ;
      T3 = T2 + small ;
      T4 = sqrt (T3 * T3 + 4.0 * vth_dlt) ;
      T5 = 0.5 * (T3 + T4) ;
      T6 = 0.5 * (1.0 + T3 / T4) ;
      T5_dVb = - Vbsz2_dVbs * T6 ;
      T5_dVd = - Vbsz2_dVds * T6 ;
      T5_dVg = - Vbsz2_dVgs * T6 ;
      T5_dT = - Vbsz2_dT * T6 ;
      T7 = 1.0 / T5 ;
      bs12 = model->HSMHV_bs1 * T7 ;
      T8 = - bs12 * T7 ;
      bs12_dVb = T8 * T5_dVb ;
      bs12_dVd = T8 * T5_dVd ;
      bs12_dVg = T8 * T5_dVg ;
      bs12_dT = T8 * T5_dT ;
      Fn_SU( T10 , Vbsz2 + bs12, 0.93 * Pb20, vth_dlt, T0) ;
      Qbmm = sqrt (T1 * (Pb20 - T10 )) ;
      T9 = T0 / Qbmm ;
      Qbmm_dVb = 0.5 * T1 * - (Vbsz2_dVbs + bs12_dVb) * T9 ;
      Qbmm_dVd = 0.5 * T1 * - (Vbsz2_dVds + bs12_dVd) * T9 ;
      Qbmm_dVg = 0.5 * T1 * - (Vbsz2_dVgs + bs12_dVg) * T9 ;
      Qbmm_dT = 0.5 * T1 * - (Vbsz2_dT + bs12_dT) * T9 ;

      dqb = (Qb0 - Qbmm) * Cox_inv ;
      dqb_dVb = Vthp_dVb - Qbmm_dVb * Cox_inv - Qbmm * Cox_inv_dVb ;
      dqb_dVd = Vthp_dVd - Qbmm_dVd * Cox_inv - Qbmm * Cox_inv_dVd ;
      dqb_dVg = Vthp_dVg - Qbmm_dVg * Cox_inv - Qbmm * Cox_inv_dVg ;
      dqb_dT  = Qb0_dT  * Cox_inv + Qb0 * Cox_inv_dT - Qbmm_dT * Cox_inv - Qbmm * Cox_inv_dT ;

      T1 = 2.0 * C_QE * here->HSMHV_nsubc * C_ESI ;
      T2 = sqrt( T1 * ( Pb2c - Vbsz2 ) ) ;
      Vth0 = Pb2c + Vfb + T2 * Cox_inv ;
      T3 = 0.5 * T1 / T2 * Cox_inv ;
      Vth0_dVb = T3 * ( - Vbsz2_dVbs ) + T2 * Cox_inv_dVb ;
      Vth0_dVd = T3 * ( - Vbsz2_dVds ) + T2 * Cox_inv_dVd ;
      Vth0_dVg = T3 * ( - Vbsz2_dVgs ) + T2 * Cox_inv_dVg ;
      Vth0_dT = T3 * ( - Vbsz2_dT ) + T2 * Cox_inv_dT ;

      T1 = C_ESI * Cox_inv ;
      T2 = here->HSMHV_wdplp ;
      T4 = 1.0e0 / ( model->HSMHV_lp * model->HSMHV_lp ) ;
      T5 = 2.0e0 * ( model->HSMHV_vbi - Pb20b ) * T1 * T2 * T4 ;
      dVth0 = T5 * sqrt_Pbsum ;
      T6 = 0.5 * T5 / sqrt_Pbsum ;
      T7 = 2.0e0 * ( model->HSMHV_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
      T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
      dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
      dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
      dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;
      dVth0_dT = T6 * Pbsum_dT + T7 * Cox_inv_dT + T8 * Pb20b_dT ;
      
      T1 = Vthp - Vth0 ;
      T1_dVb = Vthp_dVb - Vth0_dVb ;
      T2 = pParam->HSMHV_scp1 + pParam->HSMHV_scp3 * Pbsum / model->HSMHV_lp ;
      T2_dVb = pParam->HSMHV_scp3 * Pbsum_dVb / model->HSMHV_lp ;
      T3 = T2 + pParam->HSMHV_scp2 * Vdsz ;
      T3_dVb = T2_dVb + pParam->HSMHV_scp2 * Vdsz_dVbs ;


      Vdx = model->HSMHV_scp21 + Vdsz ;
      Vdx_dVbs = Vdsz_dVbs ;
/*      Vdx_dT = Vdsz_dT ;*/
      Vdx2 = Vdx * Vdx ;
      Vdx2_dVbs = 2 * Vdx_dVbs * Vdx ;
/*      Vdx2_dT = 2 * Vdx_dT * Vdx ;*/
      
      dVthLP = T1 * dVth0 * T3 + dqb - here->HSMHV_msc / Vdx2 ;
      dVthLP_dVb = T1_dVb * dVth0 * T3 + T1 * dVth0_dVb * T3 +  T1 * dVth0 * T3_dVb 
                   + dqb_dVb + here->HSMHV_msc / Vdx2 /Vdx2 *Vdx2_dVbs;
      T4 = T1 * dVth0 * pParam->HSMHV_scp3 / model->HSMHV_lp ;
      dVthLP_dVd = (Vthp_dVd - Vth0_dVd) * dVth0 * T3 + T1 * dVth0_dVd * T3 
                     + T4 * Pbsum_dVd
                     + T1 * dVth0 * pParam->HSMHV_scp2 * Vdsz_dVds
                     + dqb_dVd
                 + 2.0e0 * here->HSMHV_msc * Vdx * Vdsz_dVds / ( Vdx2 * Vdx2 ) ;
      dVthLP_dVg = (Vthp_dVg - Vth0_dVg) * dVth0 * T3 + T1 * dVth0_dVg * T3
                     + T4 * Pbsum_dVg + dqb_dVg ;
      dVthLP_dT = (Vthp_dT - Vth0_dT) * dVth0 * T3 + T1 * dVth0_dT * T3
                     + T4 * Pbsum_dT
                     + T1 * dVth0 * pParam->HSMHV_scp2 * Vdsz_dT
                     + dqb_dT
                 + 2.0e0 * here->HSMHV_msc * Vdx * Vdsz_dT / ( Vdx2 * Vdx2 );
    } else {
      dVthLP = 0.0e0 ;
      dVthLP_dVb = 0.0e0 ;
      dVthLP_dVd = 0.0e0 ;
      dVthLP_dVg = 0.0e0 ;
      dVthLP_dT = 0.0e0 ;
    }

    /*---------------------------------------------------*
     * dVthSC : Short-channel effect induced by Vds.
     *-----------------*/
    T1 = C_ESI * Cox_inv ;
    T2 = here->HSMHV_wdpl ;
    T3 = here->HSMHV_lgate - model->HSMHV_parl2 ;
    T4 = 1.0e0 / ( T3 * T3 ) ;
    T5 = 2.0e0 * ( model->HSMHV_vbi - Pb20b ) * T1 * T2 * T4 ;

    dVth0 = T5 * sqrt_Pbsum ;
    T6 = T5 / 2.0 / sqrt_Pbsum ;
    T7 = 2.0e0 * ( model->HSMHV_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
    T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
    dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
    dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
    dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;
    dVth0_dT = T6 * Pbsum_dT  + T7 * Cox_inv_dT + T8 * Pb20b_dT ;

    T1 = pParam->HSMHV_sc3 / here->HSMHV_lgate ;
    T4 = pParam->HSMHV_sc1 + T1 * Pbsum ;
    T4_dVb = T1 * Pbsum_dVb ;
    T4_dVd = T1 * Pbsum_dVd ;
    T4_dVg = T1 * Pbsum_dVg ;
    T4_dT = T1 * Pbsum_dT ;

    T5 = T4 + pParam->HSMHV_sc2 * Vdsz * ( 1.0 +  model->HSMHV_sc4 * Pbsum );
    T5_dVb = T4_dVb + pParam->HSMHV_sc2 * Vdsz * model->HSMHV_sc4 * Pbsum_dVb 
             + pParam->HSMHV_sc2 * Vdsz_dVbs * model->HSMHV_sc4 * Pbsum;
    T5_dVd = T4_dVd + pParam->HSMHV_sc2 * Vdsz_dVds * ( 1.0 + model->HSMHV_sc4 * Pbsum )
             + pParam->HSMHV_sc2 * Vdsz * model->HSMHV_sc4 * Pbsum_dVd;
    T5_dVg = T4_dVg + pParam->HSMHV_sc2 * Vdsz * model->HSMHV_sc4 * Pbsum_dVg;
    T5_dT  = T4_dT  + pParam->HSMHV_sc2 * Vdsz * model->HSMHV_sc4 * Pbsum_dT
            + pParam->HSMHV_sc2 * Vdsz_dT * model->HSMHV_sc4 * Pbsum;

    dVthSC = dVth0 * T5 ;
    dVthSC_dVb = dVth0_dVb * T5 + dVth0 * T5_dVb ;
    dVthSC_dVd = dVth0_dVd * T5 + dVth0 * T5_dVd ;
    dVthSC_dVg = dVth0_dVg * T5 + dVth0 * T5_dVg ;
    dVthSC_dT  = dVth0_dT  * T5 + dVth0 * T5_dT ;

    /*---------------------------------------------------*
     * dVthW : narrow-channel effect.
     *-----------------*/
    T1 = 1.0 / Cox ;
    T2 = T1 * T1 ;
    T3 = 1.0 / ( Cox +  pParam->HSMHV_wfc / Weff ) ;
    T4 = T3 * T3 ;
    T5 = T1 - T3 ;
    T6 = Qb0 * ( T2 - T4 ) ;

    dVthW = Qb0 * T5 + pParam->HSMHV_wvth0 / here->HSMHV_wg ;
    dVthW_dVb = Qb0_dVb * T5 - Cox_dVb * T6 ;
    dVthW_dVd = Qb0_dVd * T5 - Cox_dVd * T6 ;
    dVthW_dVg =              - Cox_dVg * T6 ;
    dVthW_dT = Qb0_dT * T5 - Cox_dT * T6 ;

    /*---------------------------------------------------*
     * dVth : Total variation. 
     * - Positive dVth means the decrease in Vth.
     *-----------------*/
    dVth = dVthSC + dVthLP + dVthW + here->HSMHV_dVthsm ;
    dVth_dVb = dVthSC_dVb + dVthLP_dVb + dVthW_dVb ;
    dVth_dVd = dVthSC_dVd + dVthLP_dVd + dVthW_dVd ;
    dVth_dVg = dVthSC_dVg + dVthLP_dVg + dVthW_dVg ;
    dVth_dT = dVthSC_dT + dVthLP_dT + dVthW_dT ;

    /*---------------------------------------------------*
     * Vth : Threshold voltage. 
     *-----------------*/
    Vth = Vthq - dVth ;

    /*-----------------------------------------------------------*
     * Constants in the equation of Ps0 . 
     *-----------------*/

    fac1 = cnst0 * Cox_inv ;
    fac1_dVbs = cnst0 * Cox_inv_dVb ;
    fac1_dVds = cnst0 * Cox_inv_dVd ;
    fac1_dVgs = cnst0 * Cox_inv_dVg ;

    fac1p2 = fac1 * fac1 ;
    fac1_dT = Cox_inv * cnst0_dT ;
    fac1p2_dT = 2.0 * fac1 * fac1_dT ;

    /*---------------------------------------------------*
     * Poly-Depletion Effect
     *-----------------*/

  if ( here->HSMHV_flg_pgd == 0 ) {
    dPpg = 0.0 ;
    dPpg_dVb = 0.0 ;
    dPpg_dVd = 0.0 ;
    dPpg_dVg = 0.0 ;
    dPpg_dT = 0.0 ;
  } else {
    T7 = Vgs ;
    T7_dVd = 0.0 ;
    T7_dVg = 1.0 ;

    T8 = Vds ;
    T8_dVd = 1.0 ;

    T0 = here->HSMHV_cnstpgd ;

    TX = pParam->HSMHV_pgd3 ;
    TY = FMDVDS * TX + ( 1.0 - FMDVDS ) * 0.5 ;
    T1 = TX - 0.5 ;
    TY_dVbs = T1 * FMDVDS_dVbs ;
    TY_dVds = T1 * FMDVDS_dVds ;
    TY_dVgs = T1 * FMDVDS_dVgs ;

    FMDVGS = 1.0 ;
    FMDVGS_dVgs = 0.0 ;
    if ( model->HSMHV_pgd2 > Vfb ) {
        T1 = model->HSMHV_pgd2 - Vfb ;
        T2 = ( Vgs - Vfb ) / T1 ;
        Fn_SZ( T3 , T2 , 1e-3 , T4 ) ;
        Fn_SU( T5 , T3 , 1.0 , 1e-3 , T6 ) ;
        T5_dVg = T4 * T6 / T1 ;
        FMDVGS = T5 * T5 ;
        FMDVGS_dVgs = 2 * T5 * T5_dVg ;
    }
    FMDPG = FMDVDS * FMDVGS ;
    FMDPG_dVbs = FMDVDS_dVbs * FMDVGS ;
    FMDPG_dVds = FMDVDS_dVds * FMDVGS ;
    FMDPG_dVgs = FMDVDS_dVgs * FMDVGS + FMDVDS * FMDVGS_dVgs ;
    FMDPG_dT = FMDVDS_dT * FMDVGS ;


    TX = pParam->HSMHV_pgd3 ;
    TY = FMDPG * TX + ( 1.0 - FMDPG ) * 0.5 ;
    T1 = TX - 0.5 ;
    TY_dVbs = T1 * FMDPG_dVbs ;
    TY_dVds = T1 * FMDPG_dVds ;
    TY_dVgs = T1 * FMDPG_dVgs ;
    TY_dT = T1 * FMDPG_dT ;
    if ( TX == 0.0 )  { TY =0.0 ; TY_dVbs =0.0 ; TY_dVds =0.0 ; TY_dVgs =0.0 ; TY_dT =0.0 ; }

    T3 = T7 - model->HSMHV_pgd2 - TY * T8 ;
    T3_dVb = - TY_dVbs * T8 ;
    T3_dVd = T7_dVd - ( TY_dVds * T8 + TY * T8_dVd ) ;
    T3_dVg = T7_dVg - ( TY_dVgs * T8 ) ;
    T3_dT = -  TY_dT * T8;

    Fn_ExpLim( dPpg , T3 , T6 ) ;
    dPpg *= T0 ;
    dPpg_dVb = T0 * T6 * T3_dVb ;
    dPpg_dVd = T0 * T6 * T3_dVd ;
    dPpg_dVg = T0 * T6 * T3_dVg ;
    dPpg_dT = T0 * T6 * T3_dT ;

    Fn_SU( dPpg , dPpg , pol_b , pol_dlt , T9 ) ;
    dPpg_dVb *= T9 ;
    dPpg_dVd *= T9 ;
    dPpg_dVg *= T9 ;
    dPpg_dT *= T9 ;

    /* damping in accumulation zone */

    T0 = Vfb + Vbsz ;
    T0_dVb = Vbsz_dVbs ;
    T1 = 0.6 * ( Vthq - T0 ) ;
    T1_dVb = 0.6 * ( Vthq_dVb - Vbsz_dVbs ) ;
    T1_dVd = 0.6 * ( Vthq_dVd - Vbsz_dVds ) ;
    Fn_SZ( T1 , T1 , 1e-2 , T2 ) ;
    T1_dVb *= T2 ;
    T1_dVd *= T2 ;
    T1 += T0 ;
    T1_dVb += Vbsz_dVbs ;
    T4 = 1.0 / ( T1 - T0 ) ;
    T5 = T4 * T4 ;
    T4_dVb = - ( T1_dVb - Vbsz_dVbs ) * T5 ;
    T4_dVd = - ( T1_dVd ) * T5 ;

    T6 = Vgsz - T0 ;
    T6_dVb = Vgsz_dVbs - T0_dVb ;
    dmpacc = T6 * T4 ;
    dmpacc_dVbs = T6 * T4_dVb + T6_dVb * T4 ;
    dmpacc_dVds = T6 * T4_dVd + ( Vgsz_dVds - Vbsz_dVds ) * T4 ;
    dmpacc_dVgs = Vgsz_dVgs * T4 ;
    Fn_SZ( dmpacc , dmpacc , 0.3 , T1 ) ;
    dmpacc_dVbs *= T1 ;
    dmpacc_dVds *= T1 ;
    dmpacc_dVgs *= T1 ;
    Fn_SU( dmpacc , dmpacc ,1.0 , 0.1 , T1 ) ;
    dmpacc_dVbs *= T1 ;
    dmpacc_dVds *= T1 ;
    dmpacc_dVgs *= T1 ;

  }

    /*---------------------------------------------------*
     * Vgp : Effective gate bias with SCE & RSCE & flatband. 
     *-----------------*/
    Vgp = Vgs - Vfb + dVth - dPpg ;
    Vgp_dVbs = dVth_dVb - dPpg_dVb ;
    Vgp_dVds = dVth_dVd - dPpg_dVd ;
    Vgp_dVgs = 1.0e0 + dVth_dVg - dPpg_dVg ;
    Vgp_dT = dVth_dT - dPpg_dT ;
   
 
    /*---------------------------------------------------*
     * Vgs_fb : Actual flatband voltage taking account Vbscl. 
     * - note: if Vgs == Vgs_fb then Vgp == Ps0 == Vbscl . 
     *------------------*/
    Vgs_fb = Vfb - dVth + dPpg + Vbscl ;

    
    /*---------------------------------------------------*
     * Vfbsft : Vfb shift (trial for Vbscl >> 0)
     *-----------------*/
    Vfbsft = 0.0 ;
    Vfbsft_dVbs = 0.0 ;
    Vfbsft_dVds = 0.0 ;
    Vfbsft_dVgs = 0.0 ;

    if ( Vbscl > 0.0 ) {
      /* values at D2/D3 boundary + beta */
      /* Ps0 */
      T1 = Vbscl + ( znbd5 + 1 ) * beta_inv ; 
      T1_dT = Vbscl_dT + ( znbd5 + 1 ) * beta_inv_dT ;
      /* Qb0 */
      /* T2 = cnst0 * sqrt( znbd5 ) */ 
      T2 = cnst0 * 2.23606797749979 ;
      T2_dT = cnst0_dT * 2.23606797749979 ;

      /* Vgp assuming Qn0=0 */
      T3 = T2 * Cox_inv + T1 ; 
      T3_dT = T2_dT * Cox_inv + T1_dT ;

      /* Vgp difference */
      TX = T3 - Vgp ;
      TX_dVbs = T2 * Cox_inv_dVb + Vbscl_dVbs - Vgp_dVbs ;
      TX_dVds = T2 * Cox_inv_dVd - Vgp_dVds ;
      TX_dVgs = T2 * Cox_inv_dVg - Vgp_dVgs ;
      TX_dT = T3_dT - Vgp_dT ;

      /* set lower limit to 0 */
      Fn_SZ( TX , TX , 0.1 , T4 ) ;
      TX_dVbs *= T4 ;
      TX_dVds *= T4 ;
      TX_dVgs *= T4 ;
      TX_dT *=   T4 ;

      /* TY: damping factor */
      T1 = 0.5 ;
      T5 = Vbscl / T1 ;
      T5_dVb = Vbscl_dVbs / T1 ;
      T5_dT = Vbscl_dT / T1 ;
      T0 = T5 * T5 ;
      T6 = T0 * T0 ; 
      T6_dVb = 4 * T0 * T5 * T5_dVb ;
      T6_dT = 4 * T0 * T5 * T5_dT ;
      T7 = 1.0 / ( 1.0 + T6 ) ;
      T8 = T7 * T7 ;
      TY = 1.0 - T7 ;
      TY_dVbs = T8 * T6_dVb ;
      TY_dT =   T8 * T6_dT ;

      TX = TY = 0.0 ;
      Vfbsft = TX * TY ;
      Vfbsft_dVbs = TX_dVbs * TY + TX * TY_dVbs ;
      Vfbsft_dVds = TX_dVds * TY ;
      Vfbsft_dVgs = TX_dVgs * TY ;
      Vfbsft_dT =   TX_dT   * TY + TX * TY_dT ;

      Vgs_fb -= Vfbsft ;

      Vgp += Vfbsft ;
      Vgp_dVbs += Vfbsft_dVbs ;
      Vgp_dVds += Vfbsft_dVds ;
      Vgp_dVgs += Vfbsft_dVgs ;
      Vgp_dT += Vfbsft_dT ;
   
    }


    /*-----------------------------------------------------------*
     * Accumulation zone. (zone-A)
     * - evaluate basic characteristics and exit from this part.
     *-----------------*/
    if ( Vgs < Vgs_fb ) { 

      flg_zone = -1 ;
 
      /*---------------------------------------------------*
       * Evaluation of Ps0.
       * - Psa : Analytical solution of 
       *             Cox( Vgp - Psa ) = cnst0 * Qacc
       *         where Qacc is the 3-degree series of (fdep)^{1/2}.
       *         The unkown is transformed to Chi=beta(Ps0-Vbs).
       * - Ps0_min : |Ps0_min| when Vbs=0.
       *-----------------*/
      Ps0_min = here->HSMHV_eg - Pb2 ;
      Ps0_min_dT = Eg_dT - Pb2_dT ;

      TX = beta * ( Vgp - Vbscl ) ;
      TX_dVbs = beta * ( Vgp_dVbs - Vbscl_dVbs ) ;
      TX_dVds = beta * Vgp_dVds ;
      TX_dVgs = beta * Vgp_dVgs ;
      TX_dT   = beta_dT * ( Vgp - Vbscl ) + beta * ( Vgp_dT - Vbscl_dT );

      T1 = 1.0 / ( beta * cnst0 ) ;
      TY = T1 * Cox ;
      TY_dVbs = T1 * Cox_dVb ;
      TY_dVds = T1 * Cox_dVd ;
      TY_dVgs = T1 * Cox_dVg ;
      T1_dT = - T1 / ( beta * cnst0 ) * ( beta_dT * cnst0 + beta * cnst0_dT ) ;
      TY_dT = T1_dT * Cox ;

      Ac41 = 2.0 + 3.0 * C_SQRT_2 * TY ;
      
      Ac4 = 8.0 * Ac41 * Ac41 * Ac41 ;
      T1 = 72.0 * Ac41 * Ac41 * C_SQRT_2 ;
      Ac4_dVbs = T1 * TY_dVbs ;
      Ac4_dVds = T1 * TY_dVds ;
      Ac4_dVgs = T1 * TY_dVgs ;
      Ac4_dT = T1 * TY_dT ;

      T4 = ( TX - 2.0 ) ;
      T5 = 9.0 * TY * T4 ;
      T5_dVb = 9.0 * ( TY_dVbs * T4 + TY * TX_dVbs ) ;
      T5_dVd = 9.0 * ( TY_dVds * T4 + TY * TX_dVds ) ;
      T5_dVg = 9.0 * ( TY_dVgs * T4 + TY * TX_dVgs ) ;
      T5_dT = 9.0 * ( TY_dT * T4 + TY * TX_dT ) ;


      Ac31 = 7.0 * C_SQRT_2 - T5 ;
      Ac31_dVbs = -T5_dVb ;
      Ac31_dVds = -T5_dVd ;
      Ac31_dVgs = -T5_dVg ;
      Ac31_dT = -T5_dT ;

      Ac3 = Ac31 * Ac31 ;
      T1 = 2.0 * Ac31 ;
      Ac3_dVbs = T1 * Ac31_dVbs ;
      Ac3_dVds = T1 * Ac31_dVds ;
      Ac3_dVgs = T1 * Ac31_dVgs ;
      Ac3_dT   = T1 * Ac31_dT ;

      Ac2 = sqrt( Ac4 + Ac3 ) ;
      T1 = 0.5 / Ac2 ;
      Ac2_dVbs = T1 * ( Ac4_dVbs + Ac3_dVbs ) ;
      Ac2_dVds = T1 * ( Ac4_dVds + Ac3_dVds ) ;
      Ac2_dVgs = T1 * ( Ac4_dVgs + Ac3_dVgs ) ;
      Ac2_dT = T1 * ( Ac4_dT + Ac3_dT ) ;
      

      Ac1 = -7.0 * C_SQRT_2 + Ac2 + T5 ;
      Ac1_dVbs = Ac2_dVbs + T5_dVb ;
      Ac1_dVds = Ac2_dVds + T5_dVd ;
      Ac1_dVgs = Ac2_dVgs + T5_dVg ;
      Ac1_dT   = Ac2_dT + T5_dT ;

      Acd = Fn_Pow( Ac1 , C_1o3 ) ;
      T1 = C_1o3 / ( Acd * Acd ) ;
      Acd_dVbs = Ac1_dVbs * T1 ;
      Acd_dVds = Ac1_dVds * T1 ;
      Acd_dVgs = Ac1_dVgs * T1 ;
      Acd_dT = Ac1_dT * T1 ;

      Acn = -4.0 * C_SQRT_2 - 12.0 * TY + 2.0 * Acd + C_SQRT_2 * Acd * Acd ;
      T1 = 2.0 + 2.0 * C_SQRT_2 * Acd ;
      Acn_dVbs = - 12.0 * TY_dVbs + T1 * Acd_dVbs ;
      Acn_dVds = - 12.0 * TY_dVds + T1 * Acd_dVds ;
      Acn_dVgs = - 12.0 * TY_dVgs + T1 * Acd_dVgs ;
      Acn_dT = - 12.0 * TY_dT + T1 * Acd_dT ;
     

      T1 = 1.0 / Acd ;
      Chi = Acn * T1 ;
      Chi_dVbs = ( Acn_dVbs - Chi * Acd_dVbs ) * T1 ;
      Chi_dVds = ( Acn_dVds - Chi * Acd_dVds ) * T1 ;
      Chi_dVgs = ( Acn_dVgs - Chi * Acd_dVgs ) * T1 ;
      Chi_dT = ( Acn_dT - Chi * Acd_dT ) * T1 ;


      Psa = Chi * beta_inv + Vbscl ;
      Psa_dVbs = Chi_dVbs * beta_inv + Vbscl_dVbs ;
      Psa_dVds = Chi_dVds * beta_inv ;
      Psa_dVgs = Chi_dVgs * beta_inv ;
      Psa_dT = Chi_dT * beta_inv + Chi * beta_inv_dT + Vbscl_dT;
      
      T1 = Psa - Vbscl ;
      T1_dT = Psa_dT - Vbscl_dT ;
      T2 = T1 / Ps0_min ;
      T2_dT = ( T1_dT * Ps0_min - T1 * Ps0_min_dT ) / ( Ps0_min * Ps0_min ) ;
      T3 = sqrt( 1.0 + ( T2 * T2 ) ) ;
      
      T3_dT = 1.0 / T3 * T2 * T2_dT ;

      T9 = T2 / T3 / Ps0_min ;
      T3_dVb = T9 * ( Psa_dVbs - Vbscl_dVbs ) ;
      T3_dVd = T9 * ( Psa_dVds ) ;
      T3_dVg = T9 * ( Psa_dVgs ) ;
      
      Ps0 = T1 / T3 + Vbscl ;
      T9 = 1.0 / ( T3 * T3 ) ;
      Ps0_dVbs = T9 * ( ( Psa_dVbs - Vbscl_dVbs ) * T3 - T1 * T3_dVb ) + Vbscl_dVbs ;
      Ps0_dVds = T9 * ( Psa_dVds * T3 - T1 * T3_dVd ) ;
      Ps0_dVgs = T9 * ( Psa_dVgs * T3 - T1 * T3_dVg ) ;
      Ps0_dT   = T9 * ( ( Psa_dT - Vbscl_dT )* T3 - T1 * T3_dT ) + Vbscl_dT;

      /*---------------------------------------------------*
       * Characteristics. 
       *-----------------*/
      Psl = Ps0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVds = Ps0_dVds ;
      Psl_dVgs = Ps0_dVgs ;
      Psl_dT   = Ps0_dT ;
      
      /** (reminder)
      Psdl = Psl ;
      Psdl_dVbs = Psl_dVbs ;
      Psdl_dVds = Psl_dVds ;      
      Psdl_dVgs = Psl_dVgs ;
      **/
    
      T2 = ( Vgp - Ps0 ) ;
      T2_dT = Vgp_dT - Ps0_dT ;
      Qbu = Cox * T2 ;
      Qbu_dVbs = Cox * ( Vgp_dVbs - Ps0_dVbs ) + Cox_dVb * T2 ;
      Qbu_dVds = Cox * ( Vgp_dVds - Ps0_dVds ) + Cox_dVd * T2 ;
      Qbu_dVgs = Cox * ( Vgp_dVgs - Ps0_dVgs ) + Cox_dVg * T2 ;
      Qbu_dT   = Cox * T2_dT ;

      Qiu = 0.0e0 ;
      Qiu_dVbs = 0.0e0 ;
      Qiu_dVds = 0.0e0 ;
      Qiu_dVgs = 0.0e0 ;
      Qiu_dT = 0.0e0 ;

      Qdrat = 0.0e0 ;
      Qdrat_dVbs = 0.0e0 ;
      Qdrat_dVds = 0.0e0 ;
      Qdrat_dVgs = 0.0e0 ;
      Qdrat_dT = 0.0 ;

      Lred = 0.0e0 ;
      Lred_dVbs = 0.0e0 ;
      Lred_dVds = 0.0e0 ;
      Lred_dVgs = 0.0e0 ;
      Lred_dT = 0.0e0 ;

      Ids = 0.0e0 ;
      Ids_dVbs = 0.0e0 ;
      Ids_dVds = 0.0e0 ;
      Ids_dVgs = 0.0e0 ;
      Ids_dT =  0.0e0 ;

      VgVt = 0.0 ;
      
      flg_noqi = 1 ;

      goto end_of_part_1 ;
    } 

    
    /*-----------------------------------------------------------*
     * Initial guess for Ps0. 
     *-----------------*/

    /*---------------------------------------------------*
     * Ps0_iniA: solution of subthreshold equation assuming zone-D1/D2.
     *-----------------*/
    TX = 1.0e0 + 4.0e0 
      * ( beta * ( Vgp - Vbscl ) - 1.0e0 ) / ( fac1p2 * beta2 ) ;
    TX = Fn_Max( TX , epsm10 ) ;
    Ps0_iniA = Vgp + fac1p2 * beta * 0.5 * ( 1.0e0 - sqrt( TX ) ) ;

    /* use analytical value in subthreshold region. */
    if ( Vgs < ( Vfb + Vth ) * 0.5 ) {
        flg_pprv = 0 ;
    }


    if ( flg_pprv >= 1 ) {
      /*---------------------------------------------------*
       * Use previous value.
       *-----------------*/

      T1  = Ps0_dVbs * dVbs + Ps0_dVds * dVds  + Ps0_dVgs * dVgs ;
      Ps0_ini  = Ps0 + T1 ;

        T2 = here->HSMHV_ps0_dtemp_prv * Temp_dif ;
        if ( fabs( T1 + T2 ) < dP_max ) { Ps0_ini += T2 ; } 

      if ( flg_pprv == 2 ) {
        /* TX_dVxs = d^2 Ps0 / d Vxs^2 here */
        if ( Vbsc_dif2 > epsm10 ) {
          TX_dVbs = ( here->HSMHV_ps0_dvbs_prv - here->HSMHV_ps0_dvbs_prv2 )
                  / Vbsc_dif2 ;
        } else {
          TX_dVbs = 0.0 ;
        }
        if ( Vdsc_dif2 > epsm10 ) {
          TX_dVds = ( here->HSMHV_ps0_dvds_prv - here->HSMHV_ps0_dvds_prv2 )
                  / Vdsc_dif2 ;
        } else {
          TX_dVds = 0.0 ;
        }
        if ( Vgsc_dif2 > epsm10 ) {
          TX_dVgs = ( here->HSMHV_ps0_dvgs_prv - here->HSMHV_ps0_dvgs_prv2 )
                  / Vgsc_dif2 ;
        } else {
          TX_dVgs = 0.0 ;
        }
        T2 = ( dVbs * dVbs ) / 2 * TX_dVbs
           + ( dVds * dVds ) / 2 * TX_dVds
           + ( dVgs * dVgs ) / 2 * TX_dVgs ;

        if ( fabs( T2 ) < fabs( 0.5 * T1 ) ) {
          Ps0_ini += T2 ;
        } else {
          flg_pprv = 1 ;
        }
      }

      T1 = Ps0_ini - Ps0 ;
      if ( T1 < - dP_max || T1 > dP_max ) {
        flg_pprv = 0 ; /* flag changes to analytical */
      } else {
        Ps0_iniA = Fn_Max( Ps0_ini , Ps0_iniA ) ;
      }
    } /* end of (flg_pprv >=1) if-block */

    if ( flg_pprv == 0 ) {

      /*---------------------------------------------------*
       * Analytical initial guess.
       *-----------------*/
      /*-------------------------------------------*
       * Common part.
       *-----------------*/
      Chi = beta * ( Ps0_iniA - Vbscl ) ;
 
      if ( Chi < znbd3 ) { 
        /*-----------------------------------*
         * zone-D1/D2
         * - Ps0_ini is the analytical solution of Qs=Qb0 with
         *   Qb0 being approximated to 3-degree polynomial.
         *-----------------*/
        TY = beta * ( Vgp - Vbscl ) ;
        T1 = 1.0e0 / ( cn_nc3 * beta * fac1 ) ;
        T2 = 81.0 + 3.0 * T1 ;
        T3 = -2916.0 - 81.0 * T1 + 27.0 * T1 * TY ;
        T4 = 1458.0 - 81.0 * ( 54.0 + T1 ) + 27.0 * T1 * TY ;
        T4 = T4 * T4 ;
        T5 = Fn_Pow( T3 + sqrt( 4 * T2 * T2 * T2 + T4 ) , C_1o3 ) ;
        TX = 3.0 - ( C_2p_1o3 * T2 ) / ( 3.0 * T5 )
           + 1 / ( 3.0 * C_2p_1o3 ) * T5 ;
        
        Ps0_iniA = TX * beta_inv + Vbscl ;
        Ps0_ini = Ps0_iniA ;
       
      } else if ( Vgs <= Vth ) { 
       /*-----------------------------------*
        * Weak inversion zone.  
        *-----------------*/
        Ps0_ini = Ps0_iniA ;
        
      } else { 
       /*-----------------------------------*
        * Strong inversion zone.  
        * - Ps0_iniB : upper bound.
        *-----------------*/
        T1 = 1.0 / cnst1 / cnstCoxi ;
        T2 = T1 * Vgp * Vgp ;
        T3 = beta + 2.0 / Vgp ;
        
        Ps0_iniB = log( T2 ) / T3 ;
        
        Fn_SU( Ps0_ini , Ps0_iniA, Ps0_iniB, c_ps0ini_2, T1) ;
      } 
    }
  
    TX = Vbscl + ps_conv / 2 ;
    if ( Ps0_ini < TX ) Ps0_ini = TX ;


    /*---------------------------------------------------*
     * Assign initial guess.
     *-----------------*/
    Ps0 = Ps0_ini ;
    Psl_lim = Ps0_iniA ;

    /*---------------------------------------------------*
     * Calculation of Ps0. (beginning of Newton loop) 
     * - Fs0 : Fs0 = 0 is the equation to be solved. 
     * - dPs0 : correction value. 
     *-----------------*/
    exp_bVbs = exp( beta * Vbscl ) ;
    cfs1 = cnst1 * exp_bVbs ;
    
    flg_conv = 0 ;
    for ( lp_s0 = 1 ; lp_s0 <= lp_s0_max + 1 ; lp_s0 ++ ) { 
      
      Chi = beta * ( Ps0 - Vbscl ) ;
    
      if ( Chi < znbd5 ) { 
        /*-------------------------------------------*
         * zone-D1/D2.  (Ps0)
         * - Qb0 is approximated to 5-degree polynomial.
         *-----------------*/
        fi = Chi * Chi * Chi 
          * ( cn_im53 + Chi * ( cn_im54 + Chi * cn_im55 ) ) ;
        fi_dChi = Chi * Chi 
          * ( 3 * cn_im53 + Chi * ( 4 * cn_im54 + Chi * 5 * cn_im55 ) ) ;
      
        fs01 = cfs1 * fi * fi ;
        fs01_dPs0 = cfs1 * beta * 2 * fi * fi_dChi ;

        fb = Chi * ( cn_nc51 
           + Chi * ( cn_nc52 
           + Chi * ( cn_nc53 
           + Chi * ( cn_nc54 + Chi * cn_nc55 ) ) ) ) ;
        fb_dChi = cn_nc51 
           + Chi * ( 2 * cn_nc52 
           + Chi * ( 3 * cn_nc53
           + Chi * ( 4 * cn_nc54 + Chi * 5 * cn_nc55 ) ) ) ;
        
        fs02 = sqrt( fb * fb + fs01 ) ;
        fs02_dPs0 = ( beta * fb_dChi * 2 * fb + fs01_dPs0 ) / ( fs02 + fs02 ) ;
     
      } else { 
        /*-------------------------------------------*
         * zone-D3.  (Ps0)
         *-----------------*/
        if ( Chi < large_arg ) { /* avoid exp_Chi to become extremely large */
           exp_Chi = exp( Chi ) ;
           fs01 = cfs1 * ( exp_Chi - 1.0e0 ) ;
           fs01_dPs0 = cfs1 * beta * ( exp_Chi ) ;
        } else {
          exp_bPs0 = exp( beta*Ps0 ) ;
          fs01     = cnst1 * ( exp_bPs0 - exp_bVbs ) ;
          fs01_dPs0 = cnst1 * beta * exp_bPs0 ;
        }
        fs02 = sqrt( Chi - 1.0 + fs01 ) ;
        fs02_dPs0 = ( beta + fs01_dPs0 ) / ( fs02 + fs02 ) ;
      } /* end of if ( Chi ... ) else block */

      Fs0 = Vgp - Ps0 - fac1 * fs02 ;
      Fs0_dPs0 = - 1.0e0 - fac1 * fs02_dPs0 ;

      if ( flg_conv == 1 ) break ;

      dPs0 = - Fs0 / Fs0_dPs0 ;


      /*-------------------------------------------*
       * Update Ps0 . 
       * - clamped to Vbscl if Ps0 < Vbscl .
       *-----------------*/
      dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Ps0))) ;
      if ( fabs( dPs0 ) > dPlim ) dPs0 = dPlim * Fn_Sgn( dPs0 ) ;
      
      Ps0 = Ps0 + dPs0 ;

      TX = Vbscl + ps_conv / 2 ;
      if ( Ps0 < TX ) Ps0 = TX ;
      
      /*-------------------------------------------*
       * Check convergence. 
       * NOTE: This condition may be too rigid. 
       *-----------------*/
      if ( fabs( dPs0 ) <= ps_conv && fabs( Fs0 ) <= gs_conv ) {
        flg_conv = 1 ;
      }
      
    } /* end of Ps0 Newton loop */

    /* Reduce loop count to exclude the sweep for derivative calculation */
    lp_s0 -- ;

    /*-------------------------------------------*
     * Procedure for diverged case.
     *-----------------*/
    if ( flg_conv == 0 ) { 
      fprintf( stderr , 
               "*** warning(HiSIM_HV): Went Over Iteration Maximum (Ps0)\n" ) ;
      fprintf( stderr , 
               " Vbse   = %7.3f Vdse = %7.3f Vgse = %7.3f\n" , 
               Vbse , Vdse , Vgse ) ;
      if ( flg_info >= 2 ) {
        printf( "*** warning(HiSIM_HV): Went Over Iteration Maximum (Ps0)\n" ) ;
      }
   } 

    /*---------------------------------------------------*
     * Evaluate derivatives of  Ps0. 
     * - note: Here, fs01_dVbs and fs02_dVbs are derivatives 
     *   w.r.t. explicit Vbs. So, Ps0 in the fs01 and fs02
     *   expressions is regarded as a constant.
     *-----------------*/
    /* self heating */
    Chi_dT = beta_dT *( Ps0 - Vbscl ) - beta * Vbscl_dT ;
    exp_bVbs_dT = ( beta_dT * Vbscl + beta * Vbscl_dT ) * exp_bVbs ;
    cfs1_dT     = exp_bVbs * cnst1_dT + exp_bVbs_dT * cnst1 ;
  
    /* derivatives of fs0* w.r.t. explicit Vbs */
    if ( Chi < znbd5 ) { 
      fs01_dVbs = cfs1 * beta * fi * ( fi - 2 * fi_dChi ) * Vbscl_dVbs ;
      fs01_dT = cfs1 * 2 * fi * fi_dChi * Chi_dT + fi * fi * cfs1_dT ;
      T2 = 1.0e0 / ( fs02 + fs02 ) ;
      fs02_dVbs = ( - beta * Vbscl_dVbs * fb_dChi * 2 * fb + fs01_dVbs ) * T2 ;
      fs02_dT = ( 2 * fb * fb_dChi * Chi_dT + fs01_dT ) * T2 ;
    } else {
      if ( Chi < large_arg ) {
        fs01_dVbs   = - cfs1 * beta * Vbscl_dVbs ;
        exp_Chi_dT  = exp_Chi * Chi_dT ;
        fs01_dT     = ( exp_Chi - 1.0e0 ) * cfs1_dT + cfs1 * exp_Chi_dT ;
      } else {
        fs01_dVbs   = - cfs1 * beta * Vbscl_dVbs ;
        exp_bPs0_dT = exp_bPs0 * Ps0 * beta_dT ;
        fs01_dT     = cnst1_dT*(exp_bPs0-exp_bVbs) + cnst1*(exp_bPs0_dT-exp_bVbs_dT) ;
      }
      T2 = 0.5e0 / fs02 ;
      fs02_dVbs = ( - beta * Vbscl_dVbs + fs01_dVbs ) * T2 ;
      fs02_dT = T2 * ( Chi_dT + fs01_dT ) ;
    }

    T1 = 1.0 / Fs0_dPs0 ;
    Ps0_dVbs = - ( Vgp_dVbs - ( fac1 * fs02_dVbs + fac1_dVbs * fs02 ) ) * T1 ;
    Ps0_dVds = - ( Vgp_dVds -                      fac1_dVds * fs02   ) * T1 ;
    Ps0_dVgs = - ( Vgp_dVgs -                      fac1_dVgs * fs02   ) * T1 ;
    Ps0_dT =   - ( Vgp_dT   - ( fac1 * fs02_dT   + fac1_dT   * fs02 ) ) * T1 ;

    Chi_dT = beta_dT *( Ps0 - Vbscl ) + beta * ( Ps0_dT - Vbscl_dT ) ;

    if ( Chi < znbd5 ) { 
      /*-------------------------------------------*
       * zone-D1/D2. (Ps0)
       * Xi0 := fdep0^2 = fb * fb  [D1,D2]
       *-----------------*/
      Xi0 = fb * fb + epsm10 ;
      T1 = 2 * fb * fb_dChi * beta ;
      Xi0_dVbs = T1 * ( Ps0_dVbs - Vbscl_dVbs ) ;
      Xi0_dVds = T1 * Ps0_dVds ;
      Xi0_dVgs = T1 * Ps0_dVgs ;
      Xi0_dT = 2 * fb * fb_dChi * Chi_dT ;

      Xi0p12 = fb + epsm10 ;
      T1 = fb_dChi * beta ;
      Xi0p12_dVbs = T1 * ( Ps0_dVbs - Vbscl_dVbs ) ;
      Xi0p12_dVds = T1 * Ps0_dVds ;
      Xi0p12_dVgs = T1 * Ps0_dVgs ;
      Xi0p12_dT = fb_dChi * Chi_dT ;

      Xi0p32 = fb * fb * fb + epsm10 ;
      T1 = 3 * fb * fb * fb_dChi * beta ;
      Xi0p32_dVbs = T1 * ( Ps0_dVbs - Vbscl_dVbs ) ;
      Xi0p32_dVds = T1 * Ps0_dVds ;
      Xi0p32_dVgs = T1 * Ps0_dVgs ;
      Xi0p32_dT = 3 * fb * fb * fb_dChi * Chi_dT ;

      fs01_dT = cfs1 * 2 * fi * fi_dChi * Chi_dT + fi * fi * cfs1_dT ;
      fs02_dT = ( 2 * fb * fb_dChi * Chi_dT + fs01_dT ) * T2 ;
      
    } else { 
      /*-------------------------------------------*
       * zone-D3. (Ps0)
       *-----------------*/
      flg_zone = 3 ;
      flg_noqi = 0 ;

      /*-----------------------------------*
       * Xi0 := fdep0^2 = Chi - 1 = beta * ( Ps0 - Vbscl ) - 1 [D3]
       *-----------------*/
      Xi0 = Chi - 1.0e0 ;
      Xi0_dVbs = beta * ( Ps0_dVbs - Vbscl_dVbs ) ;
      Xi0_dVds = beta * Ps0_dVds ;
      Xi0_dVgs = beta * Ps0_dVgs ;
      Xi0_dT = Chi_dT ;
 
      Xi0p12 = sqrt( Xi0 ) ;
      T1 = 0.5e0 / Xi0p12 ;
      Xi0p12_dVbs = T1 * Xi0_dVbs ;
      Xi0p12_dVds = T1 * Xi0_dVds ;
      Xi0p12_dVgs = T1 * Xi0_dVgs ;
      Xi0p12_dT = T1 * Xi0_dT ;
 
      Xi0p32 = Xi0 * Xi0p12 ;
      T1 = 1.5e0 * Xi0p12 ;
      Xi0p32_dVbs = T1 * Xi0_dVbs ;
      Xi0p32_dVds = T1 * Xi0_dVds ;
      Xi0p32_dVgs = T1 * Xi0_dVgs ;
      Xi0p32_dT = T1 * Xi0_dT ;

      if ( Chi < large_arg ) {
        exp_Chi_dT = exp_Chi * Chi_dT ;
        fs01_dT = ( exp_Chi - 1.0e0 ) * cfs1_dT + cfs1 * exp_Chi_dT ;
      } else {
        exp_bPs0_dT = exp_bPs0 * (beta_dT * Ps0 + beta * Ps0_dT) ;
        fs01_dT     = cnst1_dT*(exp_bPs0-exp_bVbs) + cnst1*(exp_bPs0_dT-exp_bVbs_dT) ;
      }
      fs02_dT = T2 * ( Chi_dT + fs01_dT ) ;

    } /* end of if ( Chi  ... ) block */
    
     /*-----------------------------------------------------------*
     * - Recalculate the derivatives of fs01 and fs02.
     * note: fs01  = cnst1 * exp( Vbs ) * ( exp( Chi ) - Chi - 1.0e0 ) ;
     *       fs02  = sqrt( Xi0 + fs01 ) ;
     *-----------------*/
    fs01_dVbs = Ps0_dVbs * fs01_dPs0 + fs01_dVbs ;
    fs01_dVds = Ps0_dVds * fs01_dPs0 ;
    fs01_dVgs = Ps0_dVgs * fs01_dPs0 ;
    fs02_dVbs = Ps0_dVbs * fs02_dPs0 + fs02_dVbs ;
    fs02_dVds = Ps0_dVds * fs02_dPs0 ;
    fs02_dVgs = Ps0_dVgs * fs02_dPs0 ;

    /*-----------------------------------------------------------*
     * Qb0 : Qb at source side.
     * Qn0 : Qi at source side.
     *-----------------*/

    Qb0 = cnst0 * Xi0p12 ;
    Qb0_dVb = cnst0 * Xi0p12_dVbs ;
    Qb0_dVd = cnst0 * Xi0p12_dVds ;
    Qb0_dVg = cnst0 * Xi0p12_dVgs ;
    Qb0_dT = cnst0 * Xi0p12_dT + cnst0_dT * Xi0p12 ;

    T1 = 1.0 / ( fs02 + Xi0p12 ) ;
    Qn0 = cnst0 * fs01 * T1 ;
    T1_dT = - T1 * T1 * ( fs02_dT + Xi0p12_dT ) ; 
    Qn0_dT = cnst0 * ( fs01 * T1_dT + T1 * fs01_dT ) + fs01 * T1 * cnst0_dT ;

    T2 = 1.0 / ( fs01 + epsm10 ) ;
 
    Qn0_dVbs = Qn0 * ( fs01_dVbs * T2 - ( fs02_dVbs + Xi0p12_dVbs ) * T1 ) ;
    Qn0_dVds = Qn0 * ( fs01_dVds * T2 - ( fs02_dVds + Xi0p12_dVds ) * T1 ) ;
    Qn0_dVgs = Qn0 * ( fs01_dVgs * T2 - ( fs02_dVgs + Xi0p12_dVgs ) * T1 ) ;
  

    /*-----------------------------------------------------------*
     * zone-D1 and D2
     *-----------------*/
    if ( Chi < znbd5 ) {
      if ( Chi < znbd3 ) {
        /*-------------------------------------------*
         * zone-D1. (Ps0)
         *-----------------*/
        flg_zone = 1 ;
        flg_noqi = 1 ; /** !! to be revisited !! **/

        Qiu = Qn0 ;
        Qiu_dVbs = Qn0_dVbs ;
        Qiu_dVds = Qn0_dVds ;
        Qiu_dVgs = Qn0_dVgs ;
        Qiu_dT = Qn0_dT;

        Qbu = Qb0 ;
        Qbu_dVbs = Qb0_dVb ;
        Qbu_dVds = Qb0_dVd ;
        Qbu_dVgs = Qb0_dVg ;
        Qbu_dT = Qb0_dT ;

        Qdrat = 0.5 ;
        Qdrat_dVbs = 0.0 ;
        Qdrat_dVds = 0.0 ;
        Qdrat_dVgs = 0.0 ;
        Qdrat_dT = 0.0;
 
        Lred = 0.0e0 ;
        Lred_dVbs = 0.0e0 ;
        Lred_dVds = 0.0e0 ;
        Lred_dVgs = 0.0e0 ;
        Lred_dT = 0.0e0 ;

        /** (reminder)
        *Psdl = Psl ;
        *Psdl_dVbs = Psl_dVbs ;
        *Psdl_dVds = Psl_dVds ;      
        *Psdl_dVgs = Psl_dVgs ;
        **/
 
      } else {
        /*-------------------------------------------*
         * zone-D2 (Ps0)
         *-----------------*/
        flg_zone = 2 ;
        flg_noqi = 0 ;
        /*-----------------------------------------------------------*
         * FD2 : connecting function for zone-D2.
         * - Qiu, Qbu, Qdrat and Lred should be interpolated later.
         *-----------------*/
        T1 = 1.0 / ( znbd5 - znbd3 ) ;
        TX = T1 * ( Chi - znbd3 ) ;
        TX_dVbs = beta * T1 * ( Ps0_dVbs - Vbscl_dVbs ) ;
        TX_dVds = beta * T1 * Ps0_dVds ;
        TX_dVgs = beta * T1 * Ps0_dVgs ;
        TX_dT = T1 * Chi_dT ;

        FD2 = TX * TX * TX * ( 10.0 + TX * ( -15.0 + TX * 6.0 ) ) ;
        T4 = TX * TX * ( 30.0 + TX * ( -60.0 + TX * 30.0 ) ) ;

        FD2_dVbs = T4 * TX_dVbs ;
        FD2_dVds = T4 * TX_dVds ;
        FD2_dVgs = T4 * TX_dVgs ;
        FD2_dT = T4 * TX_dT ;

      } /* end of zone-D2 */
    }


    /*---------------------------------------------------*
     * VgVt : Vgp - Vth_qi. ( Vth_qi is Vth for Qi evaluation. ) 
     *-----------------*/
    VgVt = Qn0 * Cox_inv ;
    VgVt_dVbs = Qn0_dVbs * Cox_inv + Qn0 * Cox_inv_dVb ;
    VgVt_dVds = Qn0_dVds * Cox_inv + Qn0 * Cox_inv_dVd ;
    VgVt_dVgs = Qn0_dVgs * Cox_inv + Qn0 * Cox_inv_dVg ;
    VgVt_dT = Qn0_dT * Cox_inv ;

    /*-----------------------------------------------------------*
     * make Qi=Qd=Ids=0 if VgVt <= VgVt_small 
     *-----------------*/
    if ( VgVt <= VgVt_small ) {
      flg_zone = 4 ;
      flg_noqi = 1 ;
        
      Psl = Ps0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVds = Ps0_dVds ;
      Psl_dVgs = Ps0_dVgs ;
      Psl_dT = Ps0_dT ;

      /** (reminder)
      *Psdl = Psl ;
      *Psdl_dVbs = Psl_dVbs ;
      *Psdl_dVds = Psl_dVds ;      
      *Psdl_dVgs = Psl_dVgs ;
      **/
      
      Pds = 0.0 ;
      Pds_dVbs = 0.0 ;
      Pds_dVds = 0.0 ;
      Pds_dVgs = 0.0 ;
      Pds_dT = 0.0 ;

      Qbu = Qb0 ;
      Qbu_dVbs = Qb0_dVb ;
      Qbu_dVds = Qb0_dVd ;
      Qbu_dVgs = Qb0_dVg ;
      Qbu_dT = Qb0_dT ;

      Qiu = 0.0 ;
      Qiu_dVbs = 0.0 ;
      Qiu_dVds = 0.0 ;
      Qiu_dVgs = 0.0 ;
      Qiu_dT = 0.0 ;
    
      Qdrat = 0.5 ;
      Qdrat_dVbs = 0.0 ;
      Qdrat_dVds = 0.0 ;
      Qdrat_dVgs = 0.0 ;
      Qdrat_dT = 0.0 ;

      Lred = 0.0 ;
      Lred_dVbs = 0.0 ;
      Lred_dVds = 0.0 ;
      Lred_dVgs = 0.0 ;
      Lred_dT = 0.0 ;
    
      Ids = 0.0e0 ;
      Ids_dVbs = 0.0e0 ;
      Ids_dVds = 0.0e0 ;
      Ids_dVgs = 0.0e0 ;
      Ids_dT =  0.0e0 ;
      
      goto end_of_part_1 ;
    }


    /*-----------------------------------------------------------*
     * Start point of Psl (= Ps0 + Pds) calculation. (label)
     *-----------------*/
/*  start_of_Psl: */


    /* Vdseff (begin) */
    Vdsorg = Vds ;

    T2 = here->HSMHV_qnsub_esi / ( Cox * Cox ) ;
    T4 = - 2.0e0 * T2 / Cox ;
    T2_dVb = T4 * Cox_dVb ;
    T2_dVd = T4 * Cox_dVd ;
    T2_dVg = T4 * Cox_dVg ;
    T2_dT  = T4 * Cox_dT  ;

    T0 = Vgp - beta_inv - Vbsz ;
    T0_dT = Vgp_dT - beta_inv_dT - Vbsz_dT ;
    Fn_SZ( T9, 1.0e0 + 2.0e0 / T2 * T0, 1e-3, TX ) ;
    T3 = sqrt( T9 ) ;
    T4 = 0.5e0 / T3 ;
    T5 = 1.0e0 / ( T2 * T2 ) ;
    T6 = T4 * 2.0e0 * T5 * TX ;
    T7 = T6 * T0 ;
    T8 = T6 * T2 ;
    T3_dVb = - T2_dVb * T7 + T8 * ( Vgp_dVbs - Vbsz_dVbs ) ;
    T3_dVd = - T2_dVd * T7 + T8 * ( Vgp_dVds - Vbsz_dVds ) ;
    T3_dVg = - T2_dVg * T7 + T8 * Vgp_dVgs ;
    T3_dT  = - T2_dT  * T7 + T8 * T0_dT ;

    T10 = Vgp + T2 * ( 1.0e0 - T3 ) ; 
    T10_dVb = Vgp_dVbs + T2_dVb * ( 1.0e0 - T3 ) - T2 * T3_dVb ;
    T10_dVd = Vgp_dVds + T2_dVd * ( 1.0e0 - T3 ) - T2 * T3_dVd ;
    T10_dVg = Vgp_dVgs + T2_dVg * ( 1.0e0 - T3 ) - T2 * T3_dVg ;
    T10_dT  = Vgp_dT   + T2_dT  * ( 1.0e0 - T3 ) - T2 * T3_dT ;
    Fn_SZ( T10 , T10 , 0.01 , T0 ) ;
    T10 += epsm10 ;
    T10_dVb *= T0 ;
    T10_dVd *= T0 ;
    T10_dVg *= T0 ;
    T10_dT *=  T0 ;

    T1 = Vds / T10 ;
    T2 = Fn_Pow( T1 , here->HSMHV_ddlt - 1.0e0 ) ;
    T7 = T2 * T1 ;
    T0 = here->HSMHV_ddlt * T2 / ( T10 * T10 ) ;
    T7_dVb = T0 * ( - Vds * T10_dVb ) ;
    T7_dVd = T0 * ( T10 - Vds * T10_dVd ) ;
    T7_dVg = T0 * ( - Vds * T10_dVg ) ;
    T7_dT =  T0 * ( - Vds * T10_dT ) ;

    T3 = 1.0 + T7 ;
    T4 = Fn_Pow( T3 , 1.0 / here->HSMHV_ddlt - 1.0 ) ;
    T6 = T4 * T3 ;
    T0 = T4 / here->HSMHV_ddlt ;
    T6_dVb = T0 * T7_dVb ;
    T6_dVd = T0 * T7_dVd ;
    T6_dVg = T0 * T7_dVg ;
    T6_dT =  T0 * T7_dT ;

    Vdseff = Vds / T6 ;
    T0 = 1.0 / ( T6 * T6 ) ;
    Vdseff_dVbs =      - Vds * T6_dVb   * T0 ;
    Vdseff_dVds = ( T6 - Vds * T6_dVd ) * T0 ;
    Vdseff_dVgs =      - Vds * T6_dVg   * T0 ;
    Vdseff_dT =        - Vds * T6_dT    * T0 ;

    Vds = Vdseff ;
    /* Vdseff (end) */


    exp_bVbsVds = exp( beta * ( Vbscl - Vds ) ) ;
    exp_bVbsVds_dT = ( beta_dT * ( Vbscl - Vds ) + beta * (Vbscl_dT - Vdseff_dT) ) * exp_bVbsVds ;

  
    /*---------------------------------------------------*
     * Skip Psl calculation when Vds is very small.
     *-----------------*/
    if ( Vds <= 0.0 ) {
      Pds = 0.0 ;
      Psl = Ps0 ;
//    flg_conv = 1 ;
      goto start_of_loopl ;
    }

    /*-----------------------------------------------------------*
     * Initial guess for Pds ( = Psl - Ps0 ). 
     *-----------------*/
    if ( flg_pprv >= 1  ) {
      /*---------------------------------------------------*
       * Use previous value.
       *-----------------*/

      T1  = Pds_dVbs * dVbs + Pds_dVds * dVds  + Pds_dVgs * dVgs ;
      Pds_ini  = Pds + T1 ;

      /* self heating */
        T2 = here->HSMHV_pds_dtemp_prv * Temp_dif ;
        if ( fabs( T1 + T2 ) < dP_max ) { Pds_ini += T2 ; }

      if ( flg_pprv == 2 ) {
        /* TX_dVxs = d^2 Pds / d Vxs^2 here */
        if ( Vbsc_dif2 > epsm10 ) {
          TX_dVbs = ( here->HSMHV_pds_dvbs_prv - here->HSMHV_pds_dvbs_prv2 )
                  / Vbsc_dif2 ;
        } else {
          TX_dVbs = 0.0 ;
        }
        if ( Vdsc_dif2 > epsm10 ) {
          TX_dVds = ( here->HSMHV_pds_dvds_prv - here->HSMHV_pds_dvds_prv2 )
                  / Vdsc_dif2 ;
        } else {
          TX_dVds = 0.0 ;
        }
        if ( Vgsc_dif2 > epsm10 ) {
          TX_dVgs = ( here->HSMHV_pds_dvgs_prv - here->HSMHV_pds_dvgs_prv2 )
                  / Vgsc_dif2 ;
        } else {
          TX_dVgs = 0.0 ;
        }
        T2 = ( dVbs * dVbs ) / 2 * TX_dVbs
           + ( dVds * dVds ) / 2 * TX_dVds
           + ( dVgs * dVgs ) / 2 * TX_dVgs ;

        if ( fabs( T2 ) < fabs( 0.5 * T1 ) ) {
          Pds_ini += T2 ;
        } else {
          flg_pprv = 1 ;
        }
      }
      
      T1 = Pds_ini - Pds ;
      if ( T1 < - dP_max || T1 > dP_max ) flg_pprv = 0 ; /* flag changes */

    } /* end of (flg_pprv>=1) if-block */
    
    if ( flg_pprv == 0 ) {
      /*---------------------------------------------------*
       * Analytical initial guess.
       *-----------------*/
      Pds_max = Fn_Max( Psl_lim - Ps0 , 0.0e0 ) ;
      
      Fn_SU( Pds_ini , Vds, (1.0e0 + c_pslini_1) * Pds_max, c_pslini_2, T1 ) ;
      Pds_ini = Fn_Min( Pds_ini , Pds_max ) ;
    }

    if ( Pds_ini < 0.0 ) Pds_ini = 0.0 ;
    else if ( Pds_ini > Vds ) Pds_ini = Vds ;


    /*---------------------------------------------------*
     * Assign initial guess.
     *-----------------*/
    Pds = Pds_ini ;
    Psl = Ps0 + Pds ;
    TX = Vbscl + ps_conv / 2 ;
    if ( Psl < TX ) Psl = TX ;
 
    /*---------------------------------------------------*
     * Calculation of Psl by solving Poisson eqn.
     * (beginning of Newton loop)
     * - Fsl : Fsl = 0 is the equation to be solved.
     * - dPsl : correction value.
     *-----------------*/
    flg_conv = 0 ;
    
    /*---------------------------------------------------*
     * start of Psl calculation. (label)
     *-----------------*/
start_of_loopl: 

    for ( lp_sl = 1 ; lp_sl <= lp_sl_max + 1 ; lp_sl ++ ) { 
      
      Chi = beta * ( Psl - Vbscl ) ;

      if ( Chi  < znbd5 ) { 
        /*-------------------------------------------*
         * zone-D2.  (Psl)
         * - Qb0 is approximated to 5-degree polynomial.
         *-----------------*/
        
        fi = Chi * Chi * Chi 
          * ( cn_im53 + Chi * ( cn_im54 + Chi * cn_im55 ) ) ;
        fi_dChi = Chi * Chi 
          * ( 3 * cn_im53 + Chi * ( 4 * cn_im54 + Chi * 5 * cn_im55 ) ) ;

        cfs1 = cnst1 * exp_bVbsVds ;

        fsl1 = cfs1 * fi * fi ;
        fsl1_dPsl = cfs1 * beta * 2 * fi * fi_dChi ;

        fb = Chi * ( cn_nc51 
           + Chi * ( cn_nc52 
           + Chi * ( cn_nc53 
           + Chi * ( cn_nc54 + Chi * cn_nc55 ) ) ) ) ;
        fb_dChi = cn_nc51 
           + Chi * ( 2 * cn_nc52 
           + Chi * ( 3 * cn_nc53
           + Chi * ( 4 * cn_nc54 + Chi * 5 * cn_nc55 ) ) ) ;

        fsl2 = sqrt( fb * fb + fsl1 ) ;
        fsl2_dPsl = ( beta * fb_dChi * 2 * fb + fsl1_dPsl ) / ( fsl2 + fsl2 ) ;
        
      } else { 
        /*-------------------------------------------*
         * zone-D3.  (Psl)
         *-----------------*/
        Rho = beta * ( Psl - Vds ) ;
        exp_Rho = exp( Rho ) ;

        fsl1 = cnst1 * ( exp_Rho - exp_bVbsVds ) ;
        fsl1_dPsl = cnst1 * beta * ( exp_Rho ) ;
        Xil = Chi - 1.0e0 ;
        fsl2 = sqrt( Xil + fsl1 ) ;
        fsl2_dPsl = ( beta + fsl1_dPsl ) / ( fsl2 + fsl2 ) ;
      }
    
      Fsl = Vgp - Psl - fac1 * fsl2 ;
      Fsl_dPsl = - 1.0e0 - fac1 * fsl2_dPsl ;

      if ( flg_conv == 1 ) break ;

      dPsl = - Fsl / Fsl_dPsl ;

      /*-------------------------------------------*
       * Update Psl .
       * - clamped to Vbscl if Psl < Vbscl .
       *-----------------*/
      dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Psl))) ;
      if ( fabs( dPsl ) > dPlim ) dPsl = dPlim * Fn_Sgn( dPsl ) ;
      if (Psl + dPsl < Ps0 ) {
         dPsl = Ps0 - Psl; 
         Psl = Ps0 ;
      } else {
         Psl = Psl + dPsl ;
      }

      TX = Vbscl + ps_conv / 2 ;
      if ( Psl < TX ) Psl = TX ;

      /*-------------------------------------------*
       * Check convergence.
       * NOTE: This condition may be too rigid.
       *-----------------*/
      if ( fabs( dPsl ) <= ps_conv && fabs( Fsl ) <= gs_conv ) {
        flg_conv = 1 ;
      }
    } /* end of Psl Newton loop */

    /* Reduce loop count to exclude derivative calculation sweep */
    lp_sl -- ;

    /*-------------------------------------------*
     * Procedure for diverged case.
     *-----------------*/
    if ( flg_conv == 0 ) {
      fprintf( stderr ,
               "*** warning(HiSIM_HV): Went Over Iteration Maximum (Psl)\n" ) ;
      fprintf( stderr ,
               " Vbse   = %7.3f Vdse = %7.3f Vgse = %7.3f\n" ,
               Vbse , Vdse , Vgse ) ;
      if ( flg_info >= 2 ) {
        printf("*** warning(HiSIM_HV): Went Over Iteration Maximum (Psl)\n" ) ;
      }
    }


    /*---------------------------------------------------*
     * Evaluate derivatives of  Psl.
     * - note: Here, fsl1_dVbs and fsl2_dVbs are derivatives 
     *   w.r.t. explicit Vbscl. So, Psl in the fsl1 and fsl2
     *   expressions is regarded as a constant.
     *-----------------*/
    Chi_dT = ( Psl - Vbscl ) * beta_dT - Vbscl_dT * beta ;

    if ( Chi < znbd5 ) { 
      T1 = cfs1 * beta * fi ;
      fsl1_dVbs = T1 * ( ( Vbscl_dVbs - Vdseff_dVbs ) * fi - 2.0 * fi_dChi * Vbscl_dVbs ) ;
      fsl1_dVds = - T1 * fi * Vdseff_dVds ;
      fsl1_dVgs = - T1 * fi * Vdseff_dVgs ;
      cfs1_dT = exp_bVbsVds * cnst1_dT + cnst1 * exp_bVbsVds_dT ;
      fsl1_dT = fi * fi * cfs1_dT + 2 * cfs1 * fi * fi_dChi * Chi_dT ;
      T2 =  0.5 / fsl2 ;
      fsl2_dVbs = ( - beta * fb_dChi * 2 * fb * Vbscl_dVbs + fsl1_dVbs ) * T2 ;
      fsl2_dVds = fsl1_dVds * T2 ;
      fsl2_dVgs = fsl1_dVgs * T2 ; 
      fsl2_dT = ( 2 * fb * fb_dChi * Chi_dT + fsl1_dT ) * T2 ;
    } else {
      Rho_dT = beta_dT * ( Psl - Vds ) - beta * Vdseff_dT ;
      exp_Rho_dT = Rho_dT * exp_Rho ;

      T1 = cnst1 * beta ;
      fsl1_dVbs = - T1 * ( exp_Rho * Vdseff_dVbs
                         + ( Vbscl_dVbs - Vdseff_dVbs ) * exp_bVbsVds );
      fsl1_dVds = - T1 * Vdseff_dVds * ( exp_Rho - exp_bVbsVds );
      fsl1_dVgs =   T1 * Vdseff_dVgs * ( - exp_Rho + exp_bVbsVds );
      fsl1_dT = cnst1 * ( exp_Rho_dT - exp_bVbsVds_dT ) + cnst1_dT * ( exp_Rho - exp_bVbsVds ) ;
      T2 = 0.5e0 / fsl2 ;
      fsl2_dVbs = ( - beta * Vbscl_dVbs + fsl1_dVbs ) * T2 ;
      fsl2_dVds = ( fsl1_dVds ) * T2 ;
      fsl2_dVgs = ( fsl1_dVgs ) * T2 ;
      fsl2_dT = ( Chi_dT + fsl1_dT ) * T2 ;
    }

    T1 = 1.0 / Fsl_dPsl ;
    Psl_dVbs = - ( Vgp_dVbs - ( fac1 * fsl2_dVbs + fac1_dVbs * fsl2 ) ) * T1 ;
    Psl_dVds = - ( Vgp_dVds - ( fac1 * fsl2_dVds + fac1_dVds * fsl2 ) ) * T1 ;
    Psl_dVgs = - ( Vgp_dVgs - ( fac1 * fsl2_dVgs + fac1_dVgs * fsl2 ) ) * T1 ;
    Psl_dT =   - ( Vgp_dT   - ( fac1 * fsl2_dT   + fac1_dT   * fsl2 ) ) * T1 ;

    Chi_dT = ( Psl - Vbscl ) * beta_dT + beta * ( Psl_dT - Vbscl_dT );
    exp_Chi_dT = exp_Chi * Chi_dT ;


    if ( Chi < znbd5 ) { 
      /*-------------------------------------------*
       * zone-D1/D2. (Psl)
       *-----------------*/
      Xil = fb * fb + epsm10 ;
      T1 = 2 * fb * fb_dChi * beta ;
      Xil_dVbs = T1 * ( Psl_dVbs - Vbscl_dVbs ) ;
      Xil_dVds = T1 * Psl_dVds ;
      Xil_dVgs = T1 * Psl_dVgs ;
      Xil_dT = 2 * fb * fb_dChi * Chi_dT ;

      Xilp12 = fb + epsm10 ;
      T1 = fb_dChi * beta ;
      Xilp12_dVbs = T1 * ( Psl_dVbs - Vbscl_dVbs ) ;
      Xilp12_dVds = T1 * Psl_dVds ;
      Xilp12_dVgs = T1 * Psl_dVgs ;
      Xilp12_dT = fb_dChi * Chi_dT ;
    
      Xilp32 = fb * fb * fb + epsm10 ;
      T1 = 3 * fb * fb * fb_dChi * beta ;
      Xilp32_dVbs = T1 * ( Psl_dVbs - Vbscl_dVbs ) ;
      Xilp32_dVds = T1 * Psl_dVds ;
      Xilp32_dVgs = T1 * Psl_dVgs ;
      Xilp32_dT = 3 * fb * fb * fb_dChi * Chi_dT ;
    
    } else { 
      /*-------------------------------------------*
       * zone-D3. (Psl)
       *-----------------*/

      Xil = Chi - 1.0e0 ;
      Xil_dVbs = beta * ( Psl_dVbs - Vbscl_dVbs ) ;
      Xil_dVds = beta * Psl_dVds ;
      Xil_dVgs = beta * Psl_dVgs ;
      Xil_dT = Chi_dT ;
 
      Xilp12 = sqrt( Xil ) ;
      T1 = 0.5e0 / Xilp12 ;
      Xilp12_dVbs = T1 * Xil_dVbs ;
      Xilp12_dVds = T1 * Xil_dVds ;
      Xilp12_dVgs = T1 * Xil_dVgs ;
      Xilp12_dT = T1 * Xil_dT ;

      Xilp32 = Xil * Xilp12 ;
      T1 = 1.5e0 * Xilp12 ;
      Xilp32_dVbs = T1 * Xil_dVbs ;
      Xilp32_dVds = T1 * Xil_dVds ;
      Xilp32_dVgs = T1 * Xil_dVgs ;
      Xilp32_dT = T1 * Xil_dT ;

    }

    /*---------------------------------------------------*
     * Assign Pds.
     *-----------------*/
    Pds = Psl - Ps0 ;

 /* if ( Pds < ps_conv ) { */
    if ( Pds < 0.0 ) { /* take care of numerical noise */
      Pds = 0.0 ;
      Psl = Ps0 ;
    }

    Pds_dVbs = Psl_dVbs - Ps0_dVbs ;
    Pds_dVds = Psl_dVds - Ps0_dVds ;
    Pds_dVgs = Psl_dVgs - Ps0_dVgs ;
    Pds_dT = Psl_dT - Ps0_dT ;
  
    /* if ( Pds < ps_conv ) { */
    if ( Pds < 0.0 ) {
      Pds_dVbs = 0.0 ;
      Pds_dVgs = 0.0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVgs = Ps0_dVgs ;
      Pds_dT = 0.0 ;
      Psl_dT = Ps0_dT ;
    }

    /* Vdseff */
    Vds = Vdsorg;

    /*-----------------------------------------------------------*
     * Evaluate Idd. 
     * - Eta : substantial variable of QB'/Pds and Idd/Pds. 
     * - note: Eta   = 4 * GAMMA_{hisim_0} 
     *-----------------*/
    T1 = beta / Xi0 ;
    Eta = T1 * Pds ;
    T2 = Eta * beta_inv ;
    Eta_dVbs = T1 * ( Pds_dVbs - Xi0_dVbs * T2 ) ;
    Eta_dVds = T1 * ( Pds_dVds - Xi0_dVds * T2 ) ;
    Eta_dVgs = T1 * ( Pds_dVgs - Xi0_dVgs * T2 ) ;
    T1_dT = ( beta_dT * Xi0 - beta * Xi0_dT ) / Xi0 / Xi0 ;
    Eta_dT = T1_dT * Pds + T1 * Pds_dT ;

    /* ( Eta + 1 )^n */
    Eta1 = Eta + 1.0e0 ;
    Eta1p12 = sqrt( Eta1 ) ;
    Eta1p32 = Eta1p12 * Eta1 ;
    Eta1p52 = Eta1p32 * Eta1 ;
    Eta1_dT = Eta_dT ;
    Eta1p12_dT = 0.5e0 / Eta1p12 * Eta1_dT ;
    Eta1p32_dT = Eta1p12_dT * Eta1 + Eta1p12 * Eta1_dT ;
    Eta1p52_dT = Eta1p32_dT * Eta1 + Eta1p32 * Eta1_dT ;
 
    /* 1 / ( ( Eta + 1 )^n + 1 ) */
    Zeta12 = 1.0e0 / ( Eta1p12 + 1.0e0 ) ;
    Zeta32 = 1.0e0 / ( Eta1p32 + 1.0e0 ) ;
    Zeta52 = 1.0e0 / ( Eta1p52 + 1.0e0 ) ;
    Zeta12_dT = - 1.0e0 / ( Eta1p12 + 1.0e0 ) / ( Eta1p12 + 1.0e0 ) * Eta1p12_dT ;
    Zeta32_dT = - 1.0e0 / ( Eta1p32 + 1.0e0 ) / ( Eta1p32 + 1.0e0 ) * Eta1p32_dT ;
    Zeta52_dT = - 1.0e0 / ( Eta1p52 + 1.0e0 ) / ( Eta1p52 + 1.0e0 ) * Eta1p52_dT ;    

    /*---------------------------------------------------*
     * F00 := PS00/Pds (n=1/2) 
     *-----------------*/
    F00 = Zeta12 / Xi0p12 ;
    T3 = - 1 / Xi0  ;
    T4 = - 0.5e0 / Eta1p12 * F00 ;
    T5 = Zeta12 * T3 ;
    T6 = Zeta12 * T4 ;
    F00_dVbs = ( Xi0p12_dVbs * T5 + Eta_dVbs * T6 ) ;
    F00_dVds = ( Xi0p12_dVds * T5 + Eta_dVds * T6 ) ;
    F00_dVgs = ( Xi0p12_dVgs * T5 + Eta_dVgs * T6 ) ;
    F00_dT = ( Zeta12_dT * Xi0p12 - Zeta12 * Xi0p12_dT ) / Xi0p12 / Xi0p12 ;

    /*---------------------------------------------------*
     * F10 := PS10/Pds (n=3/2) 
     *-----------------*/
    T1 = 3.0e0 + Eta * ( 3.0e0 + Eta ) ;
    F10 = C_2o3 * Xi0p12 * Zeta32 * T1 ;
    T2 = 3.0e0 + Eta * 2.0e0 ;
    T3 = C_2o3 * T1 ;
    T4 = - 1.5e0 * Eta1p12 * F10 + C_2o3 * Xi0p12 * T2 ;
    T5 = Zeta32 * T3 ;
    T6 = Zeta32 * T4 ;
    F10_dVbs = ( Xi0p12_dVbs * T5 + Eta_dVbs * T6 ) ;
    F10_dVds = ( Xi0p12_dVds * T5 + Eta_dVds * T6 ) ;
    F10_dVgs = ( Xi0p12_dVgs * T5 + Eta_dVgs * T6 ) ;
    T1_dT = Eta_dT * ( 3.0e0 + Eta ) + Eta * Eta_dT ;
    F10_dT = C_2o3 * Xi0p12 * Zeta32 * T1_dT
      + C_2o3 * T1 * ( Xi0p12 * Zeta32_dT + Zeta32 * Xi0p12_dT ) ;

    /*---------------------------------------------------*
     * F30 := PS30/Pds (n=5/2) 
     *-----------------*/
    T1 = 5e0 + Eta * ( 10e0 + Eta * ( 10e0 + Eta * ( 5e0 + Eta ) ) ) ;
    F30 = 4e0 / ( 15e0 * beta ) * Xi0p32 * Zeta52 * T1 ;
    T2 = 10e0 + Eta * ( 20e0 + Eta * ( 15e0 + Eta * 4e0 ) ) ;
    T3 = 4e0 / ( 15e0 * beta ) * T1 ;
    T4 = - ( 5e0 / 2e0 ) * Eta1p32 * F30 + 4e0 / ( 15e0 * beta ) * Xi0p32 * T2 ;
    T5 = Zeta52 * T3 ;
    T6 = Zeta52 * T4 ;
    F30_dVbs = ( Xi0p32_dVbs * T5 + Eta_dVbs * T6 ) ;
    F30_dVds = ( Xi0p32_dVds * T5 + Eta_dVds * T6 ) ;
    F30_dVgs = ( Xi0p32_dVgs * T5 + Eta_dVgs * T6 ) ;
    T1_dT = ( 10e0 + 20e0 * Eta + 15e0 * Eta * Eta + 4e0 * Eta * Eta * Eta ) * Eta_dT ;
    F30_dT = 4e0 / 15e0 * beta_inv_dT * ( Xi0p32 * Zeta52 * T1 )
      + 4e0 / 15e0 * beta_inv * ( Xi0p32_dT * Zeta52 * T1 + Xi0p32 * Zeta52_dT * T1 + Xi0p32 * Zeta52 * T1_dT ) ;

    /*---------------------------------------------------*
     * F11 := PS11/Pds. 
     *-----------------*/
    F11 = Ps0 * F10 + C_2o3 * beta_inv * Xilp32 - F30 ;
    T1 = C_2o3 * beta_inv ;
    F11_dVbs = Ps0_dVbs * F10 + Ps0 * F10_dVbs 
      + T1 * Xilp32_dVbs - F30_dVbs ;
    F11_dVds = Ps0_dVds * F10 + Ps0 * F10_dVds 
      + T1 * Xilp32_dVds - F30_dVds ;
    F11_dVgs = Ps0_dVgs * F10 + Ps0 * F10_dVgs 
      + T1 * Xilp32_dVgs - F30_dVgs ;
    F11_dT = Ps0_dT * F10 + Ps0 * F10_dT
      + C_2o3 *( beta_inv_dT * Xilp32 + beta_inv * Xilp32_dT ) - F30_dT ;

    /*---------------------------------------------------*
     * Fdd := Idd/Pds. 
     *-----------------*/
    T1 = Vgp + beta_inv - 0.5e0 * ( 2.0e0 * Ps0 + Pds ) ;
    T2 = - F10 + F00 ;
    T3 = beta * Cox ;
    T4 = beta * cnst0 ;
    Fdd = T3 * T1 + T4 * T2 ;
    Fdd_dVbs = T3 * ( Vgp_dVbs - Ps0_dVbs - 0.5e0 * Pds_dVbs ) 
      + beta * Cox_dVb * T1 + T4 * ( - F10_dVbs + F00_dVbs ) ;
    Fdd_dVds = T3 * ( Vgp_dVds - Ps0_dVds - 0.5e0 * Pds_dVds ) 
      + beta * Cox_dVd * T1 + T4 * ( - F10_dVds + F00_dVds ) ;
    Fdd_dVgs = T3 * ( Vgp_dVgs - Ps0_dVgs - 0.5e0 * Pds_dVgs ) 
      + beta * Cox_dVg * T1 + T4 * ( - F10_dVgs + F00_dVgs ) ;
    T1_dT = Vgp_dT + beta_inv_dT - 0.5e0 * ( 2.0e0 * Ps0_dT + Pds_dT ) ;
    T2_dT = -F10_dT + F00_dT ;
    T3_dT = Cox * beta_dT ;
    T4_dT = beta * cnst0_dT + cnst0 * beta_dT ;
    Fdd_dT = T1 * T3_dT + T3 * T1_dT + T2 * T4_dT + T4 * T2_dT ;
  
      
    /*---------------------------------------------------*
     *  Idd: 
     *-----------------*/
    Idd = Pds * Fdd ;
    Idd_dVbs = Pds_dVbs * Fdd + Pds * Fdd_dVbs ;
    Idd_dVds = Pds_dVds * Fdd + Pds * Fdd_dVds ;
    Idd_dVgs = Pds_dVgs * Fdd + Pds * Fdd_dVgs ;
    Idd_dT = Fdd * Pds_dT + Pds * Fdd_dT ;

    /*-----------------------------------------------------------*
     * Skip CLM and integrated charges if zone==D1
     *-----------------*/
    if( flg_zone == 1 ) {
      goto start_of_mobility ;
    }

    /*-----------------------------------------------------------*
     * Channel Length Modulation. Lred: \Delta L
     *-----------------*/
    if( pParam->HSMHV_clm2 < epsm10 && pParam->HSMHV_clm3 < epsm10 ) {
      Lred = 0.0e0 ;
      Lred_dVbs = 0.0e0 ;
      Lred_dVds = 0.0e0 ;
      Lred_dVgs = 0.0e0 ;
      Lred_dT = 0.0e0 ;

      Psdl = Psl ;
      Psdl_dVbs = Psl_dVbs ;
      Psdl_dVds = Psl_dVds ;
      Psdl_dVgs = Psl_dVgs ;
      Psdl_dT = Psl_dT ;

      if ( Psdl > Ps0 + Vds - epsm10 ) {
        Psdl = Ps0 + Vds - epsm10 ;
        Psdl_dVbs = Ps0_dVbs ;
        Psdl_dVds = Ps0_dVds + 1.0 ;
        Psdl_dVgs = Ps0_dVgs ;
        Psdl_dT = Ps0_dT ;
      }

    } else {
      T1 = here->HSMHV_wdpl ;
      T8 = sqrt (Psl - Vbscl) ;
      Wd = T1 * T8 ;
      T9 = 0.5 * T1 / T8 ;
      Wd_dVbs = T9 * (Psl_dVbs - Vbscl_dVbs) ;
      Wd_dVds = T9 * Psl_dVds ;
      Wd_dVgs = T9 * Psl_dVgs ;
      Wd_dT = T9 * (Psl_dT - Vbscl_dT) ;

      T0 = 1.0 / Wd ;
      T1 = Qn0 * T0 ;
      T2 = pParam->HSMHV_clm3 * T1 ;
      T3 = pParam->HSMHV_clm3 * T0 ;
      T2_dVb = T3 * (Qn0_dVbs - T1 * Wd_dVbs) ;
      T2_dVd = T3 * (Qn0_dVds - T1 * Wd_dVds) ;
      T2_dVg = T3 * (Qn0_dVgs - T1 * Wd_dVgs) ;
      T2_dT = T3 * (Qn0_dT - T1 * Wd_dT) ;

      T5 = pParam->HSMHV_clm2 * q_Nsub + T2 ;
      T1 = 1.0 / T5 ;
      T4 = C_ESI * T1 ;
      T4_dVb = - T4 * T2_dVb * T1 ;
      T4_dVd = - T4 * T2_dVd * T1 ;
      T4_dVg = - T4 * T2_dVg * T1 ;
      T4_dT = -T4 * T2_dT * T1 ;

      T1 = (1.0e0 - pParam->HSMHV_clm1) ;
      Psdl = pParam->HSMHV_clm1 * (Vds + Ps0) + T1 * Psl ;
      Psdl_dVbs = pParam->HSMHV_clm1 * Ps0_dVbs + T1 * Psl_dVbs ;
      Psdl_dVds = pParam->HSMHV_clm1 * (1.0 + Ps0_dVds) + T1 * Psl_dVds ;
      Psdl_dVgs = pParam->HSMHV_clm1 * Ps0_dVgs + T1 * Psl_dVgs ;
      Psdl_dT = pParam->HSMHV_clm1 * Ps0_dT + T1 * Psl_dT ;

      if ( Psdl > Ps0 + Vds - epsm10 ) {
        Psdl = Ps0 + Vds - epsm10 ;
        Psdl_dVbs = Ps0_dVbs ;
        Psdl_dVds = Ps0_dVds + 1.0 ;
        Psdl_dVgs = Ps0_dVgs ;
        Psdl_dT = Ps0_dT ;
      }
      T6 = Psdl - Psl ; 
      T6_dVb = Psdl_dVbs - Psl_dVbs ;
      T6_dVd = Psdl_dVds - Psl_dVds ;
      T6_dVg = Psdl_dVgs - Psl_dVgs ;
      T6_dT = Psdl_dT - Psl_dT ;

      T3 = beta * Qn0 ;
      T1 = 1.0 / T3 ;
      T5 = Idd * T1 ;
      T3_dT = beta * Qn0_dT + beta_dT * Qn0 ;
      T1_dT = - T1 * T1 * T3_dT ;
      T5_dT = Idd_dT * T1 + Idd * T1_dT ;
      T2 = T5 * beta ;
      T5_dVb = (Idd_dVbs - T2 * Qn0_dVbs) * T1 ;
      T5_dVd = (Idd_dVds - T2 * Qn0_dVds) * T1 ;
      T5_dVg = (Idd_dVgs - T2 * Qn0_dVgs) * T1 ;
      
      T10 = q_Nsub / C_ESI ;
      T1 = 1.0e5 ;
      T2 = 1.0 / Leff ;
      T11 = (2.0 * T5 + 2.0 * T10 * T6 * T4 + T1 * T4) * T2 ;
      T3 = T2 * T4 ;
      T7 = T11 * T4 ;
      T7_dVb = (2.0 * T5_dVb + 2.0 * T10 * (T6_dVb * T4 + T6 * T4_dVb) + T1 * T4_dVb) * T3 + T11 * T4_dVb ;
      T7_dVd = (2.0 * T5_dVd + 2.0 * T10 * (T6_dVd * T4 + T6 * T4_dVd) + T1 * T4_dVd) * T3 + T11 * T4_dVd ;
      T7_dVg = (2.0 * T5_dVg + 2.0 * T10 * (T6_dVg * T4 + T6 * T4_dVg) + T1 * T4_dVg) * T3 + T11 * T4_dVg ;
      T7_dT = (2.0 * T5_dT + 2.0 * T10 * ( T6_dT * T4 + T6 * T4_dT ) + T1 * T4_dT ) * T3 + T11 * T4_dT ;

      T11 = 4.0 * (2.0 * T10 * T6 + T1) ;
      T1 = 8.0 * T10 * T4 * T4 ;
      T2 = 2.0 * T11 * T4 ;
      T8 = T11 * T4 * T4 ;
      T8_dVb = ( T1 * T6_dVb + T2 * T4_dVb) ;
      T8_dVd = ( T1 * T6_dVd + T2 * T4_dVd) ;
      T8_dVg = ( T1 * T6_dVg + T2 * T4_dVg) ;
      T8_dT = ( T1 * T6_dT + T2 * T4_dT) ;

      T9 = sqrt (T7 * T7 + T8);
      T1 = 1.0 / T9 ;
      T2 = T7 * T1 ;
      T3 = 0.5 * T1 ;
      T9_dVb = (T2 * T7_dVb + T3 * T8_dVb) ;
      T9_dVd = (T2 * T7_dVd + T3 * T8_dVd) ;
      T9_dVg = (T2 * T7_dVg + T3 * T8_dVg) ;
      T9_dT = (T2 * T7_dT + T3 * T8_dT) ;

      Lred = 0.5 * (- T7 + T9) ;
      Lred_dVbs = 0.5 * (- T7_dVb + T9_dVb) ;
      Lred_dVds = 0.5 * (- T7_dVd + T9_dVd) ;
      Lred_dVgs = 0.5 * (- T7_dVg + T9_dVg) ;
      Lred_dT = 0.5 * (- T7_dT + T9_dT ) ;
      /*---------------------------------------------------*
       * Modify Lred for symmetry.
       *-----------------*/
      T1 = Lred ;
      Lred = FMDVDS * T1 ;
      Lred_dVbs = FMDVDS_dVbs * T1 + FMDVDS * Lred_dVbs ;
      Lred_dVds = FMDVDS_dVds * T1 + FMDVDS * Lred_dVds ;
      Lred_dVgs = FMDVDS_dVgs * T1 + FMDVDS * Lred_dVgs ;
      Lred_dT = FMDVDS_dT * T1 + FMDVDS * Lred_dT ;
    }

    /* CLM5 & CLM6 */
    Lred *= here->HSMHV_clmmod ;
    Lred_dVbs *= here->HSMHV_clmmod ;
    Lred_dVds *= here->HSMHV_clmmod ;
    Lred_dVgs *= here->HSMHV_clmmod ;
    Lred_dT *= here->HSMHV_clmmod ;

    /*---------------------------------------------------*
     * Qbu : -Qb in unit area.
     *-----------------*/
    T1 = Vgp + beta_inv ;
    T2 = T1 * F10 - F11 ;
    T1_dT = Vgp_dT + beta_inv_dT ;
    T2_dT = T1_dT * F10 + T1 * F10_dT - F11_dT ;

    
    Qbnm = cnst0 * ( cnst0 * ( 1.5e0 - ( Xi0 + 1.0e0 ) - 0.5e0 * beta * Pds )
                   + Cox * T2 ) ;
    Qbnm_dVbs = cnst0 * ( cnst0 * ( - Xi0_dVbs - 0.5e0 * beta * Pds_dVbs ) 
                          + Cox * ( Vgp_dVbs * F10 + T1 * F10_dVbs - F11_dVbs )
                          + Cox_dVb * T2 ) ;
    Qbnm_dVds = cnst0 * ( cnst0 * ( - Xi0_dVds - 0.5e0 * beta * Pds_dVds ) 
                          + Cox * ( Vgp_dVds * F10 + T1 * F10_dVds - F11_dVds )
                          + Cox_dVd * T2 ) ;
    Qbnm_dVgs = cnst0 * ( cnst0 * ( - Xi0_dVgs - 0.5e0 * beta * Pds_dVgs ) 
                          + Cox * ( Vgp_dVgs * F10 + T1 * F10_dVgs - F11_dVgs )
                          + Cox_dVg * T2 ) ;
    Qbnm_dT = cnst0_dT * ( cnst0 * ( 1.5e0 - ( Xi0 + 1.0e0 ) - 0.5e0 * beta * Pds )
                           + Cox * T2 ) 
      + cnst0 * ( cnst0_dT * ( 1.5e0 - ( Xi0 + 1.0e0 ) - 0.5e0 * beta * Pds )
                  + cnst0 * ( - Xi0_dT - 0.5 * beta_dT * Pds - 0.5 * beta * Pds_dT )
                  + Cox * T2_dT );    

    T1 = beta ;
    Qbu = T1 * Qbnm / Fdd ;
    T2 = T1 / ( Fdd * Fdd ) ;
    Qbu_dVbs = T2 * ( Fdd * Qbnm_dVbs - Qbnm * Fdd_dVbs ) ;
    Qbu_dVds = T2 * ( Fdd * Qbnm_dVds - Qbnm * Fdd_dVds ) ;
    Qbu_dVgs = T2 * ( Fdd * Qbnm_dVgs - Qbnm * Fdd_dVgs ) ;
    T1_dT = beta_dT ;
    Qbu_dT = ( Fdd * ( T1_dT * Qbnm + T1 * Qbnm_dT ) - T1 * Qbnm * Fdd_dT ) / ( Fdd * Fdd ) ;

    /*---------------------------------------------------*
     * preparation for Qi and Qd. 
     * - DtPds: Delta * Pds ;
     * - Achi: (1+Delta) * Pds ;
     *-----------------*/
    T1 = 2.0e0 * fac1 ;
    DtPds = T1 * ( F10 - Xi0p12 ) ;
    T2 = 2.0 * ( F10 - Xi0p12 ) ;
    DtPds_dVbs  = T1 * ( F10_dVbs - Xi0p12_dVbs )
                + T2 * fac1_dVbs ;
    DtPds_dVds  = T1 * ( F10_dVds - Xi0p12_dVds ) 
                + T2 * fac1_dVds ;
    DtPds_dVgs  = T1 * ( F10_dVgs - Xi0p12_dVgs ) 
                + T2 * fac1_dVgs ;
    T1_dT = 2.0e0 * fac1_dT ;
    DtPds_dT = T1_dT * ( F10 - Xi0p12 ) + T1 * ( F10_dT -Xi0p12_dT ) ;

    Achi = Pds + DtPds ;
    Achi_dVbs = Pds_dVbs + DtPds_dVbs ;
    Achi_dVds = Pds_dVds + DtPds_dVds ;
    Achi_dVgs = Pds_dVgs + DtPds_dVgs ;
    Achi_dT = Pds_dT + DtPds_dT ;
    
    /*-----------------------------------------------------------*
     * Alpha : parameter to evaluate charges. 
     * - Achi: (1+Delta) * Pds ;
     * - clamped to 0 if Alpha < 0. 
     *-----------------*/
    T1 = 1.0 / VgVt ;
    T2 = Achi * T1 ;
    T3 = 1.0e0 - T2 ;
    TX = 1.0 - T3 ;
    Fn_CP( TY , TX , 1.0 , 4 , T4 ) ;
    Alpha = 1.0 - TY ;
    T5 = T1 * T4 ;
    Alpha_dVbs = - ( Achi_dVbs - T2 * VgVt_dVbs ) * T5 ;
    Alpha_dVds = - ( Achi_dVds - T2 * VgVt_dVds ) * T5 ;
    Alpha_dVgs = - ( Achi_dVgs - T2 * VgVt_dVgs ) * T5 ;
    Alpha_dT = - ( Achi_dT - T2 * VgVt_dT ) * T5 ;


    /*-----------------------------------------------------------*
     * Qiu : -Qi in unit area.
     *-----------------*/

    Qinm = 1.0e0 + Alpha * ( 1.0e0 + Alpha ) ;
    T1 = 1.0e0 + Alpha + Alpha ;
    Qinm_dVbs = Alpha_dVbs * T1 ;
    Qinm_dVds = Alpha_dVds * T1 ;
    Qinm_dVgs = Alpha_dVgs * T1 ;
    Qinm_dT = Alpha_dT * T1 ;
 
    Qidn = Fn_Max( 1.0e0 + Alpha , epsm10 ) ;
    Qidn_dVbs = Alpha_dVbs ;
    Qidn_dVds = Alpha_dVds ;
    Qidn_dVgs = Alpha_dVgs ;
    Qidn_dT = Alpha_dT ;
    
    T1 = C_2o3 * VgVt * Qinm / Qidn ;
    Qiu = T1 * Cox ;
    T2 = 1.0 / VgVt ;
    T3 = 1.0 / Qinm ;
    T4 = 1.0 / Qidn ;
    Qiu_dVbs = Qiu * ( VgVt_dVbs * T2 + Qinm_dVbs * T3 - Qidn_dVbs * T4 )
             + T1 * Cox_dVb ;
    Qiu_dVds = Qiu * ( VgVt_dVds * T2 + Qinm_dVds * T3 - Qidn_dVds * T4) 
             + T1 * Cox_dVd ;
    Qiu_dVgs = Qiu * ( VgVt_dVgs * T2 + Qinm_dVgs * T3 - Qidn_dVgs * T4)
             + T1 * Cox_dVg ;
    T1_dT = C_2o3 * ( Qidn * ( VgVt_dT * Qinm + VgVt * Qinm_dT ) - Qidn_dT * VgVt * Qinm ) 
      / ( Qidn * Qidn ) ;
    Qiu_dT = Cox * T1_dT ;
    
  
    /*-----------------------------------------------------------*
     * Qdrat : Qd/Qi
     *-----------------*/
    Qdnm = 0.5e0 + Alpha ;
    Qdnm_dVbs = Alpha_dVbs ;
    Qdnm_dVds = Alpha_dVds ;
    Qdnm_dVgs = Alpha_dVgs ;
    Qdnm_dT = Alpha_dT ;
    
    Qddn = Qidn * Qinm ;
    Qddn_dVbs = Qidn_dVbs * Qinm + Qidn * Qinm_dVbs ;
    Qddn_dVds = Qidn_dVds * Qinm + Qidn * Qinm_dVds ;
    Qddn_dVgs = Qidn_dVgs * Qinm + Qidn * Qinm_dVgs ;
    Qddn_dT = Qidn_dT * Qinm + Qidn * Qinm_dT ;
    
    Quot = 0.4e0 * Qdnm  / Qddn ;
    Qdrat = 0.6e0 - Quot ;
 
    if ( Qdrat <= 0.5e0 ) { 
      T1 = 1.0 / Qddn ;
      T2 = 1.0 / Qdnm ;
      Qdrat_dVbs = Quot * ( Qddn_dVbs * T1 - Qdnm_dVbs * T2 ) ;
      Qdrat_dVds = Quot * ( Qddn_dVds * T1 - Qdnm_dVds * T2 ) ;
      Qdrat_dVgs = Quot * ( Qddn_dVgs * T1 - Qdnm_dVgs * T2 ) ;
      Qdrat_dT = Quot * ( Qddn_dT * T1 - Qdnm_dT * T2 ) ;
    } else { 
      Qdrat = 0.5e0 ;
      Qdrat_dVbs = 0.0e0 ;
      Qdrat_dVds = 0.0e0 ;
      Qdrat_dVgs = 0.0e0 ;
      Qdrat_dT = 0.0e0 ;

    } 
 

    /*-----------------------------------------------------------*
     * Interpolate charges and CLM for zone-D2.
     *-----------------*/

    if ( flg_zone == 2 ) {
      T1 = Qbu ;
      Qbu = FD2 * Qbu + ( 1.0 - FD2 ) * Qb0 ;
      Qbu_dVbs = FD2 * Qbu_dVbs + FD2_dVbs * T1 
        + ( 1.0 - FD2 ) * Qb0_dVb - FD2_dVbs * Qb0 ;
      Qbu_dVds = FD2 * Qbu_dVds + FD2_dVds * T1
        + ( 1.0 - FD2 ) * Qb0_dVd - FD2_dVds * Qb0 ;
      Qbu_dVgs = FD2 * Qbu_dVgs + FD2_dVgs * T1 
        + ( 1.0 - FD2 ) * Qb0_dVg - FD2_dVgs * Qb0 ;
      Qbu_dT = FD2 * Qbu_dT + FD2_dT * T1
        + ( 1.0 - FD2 ) * Qb0_dT - FD2_dT * Qb0 ;

      if ( Qbu < 0.0 ) {
        Qbu = 0.0 ;
        Qbu_dVbs = 0.0 ;
        Qbu_dVds = 0.0 ;
        Qbu_dVgs = 0.0 ;
        Qbu_dT = 0.0 ;
      }
      T1 = Qiu ;
      Qiu = FD2 * Qiu + ( 1.0 - FD2 ) * Qn0 ;
      Qiu_dVbs = FD2 * Qiu_dVbs + FD2_dVbs * T1 
        + ( 1.0 - FD2 ) * Qn0_dVbs - FD2_dVbs * Qn0 ;
      Qiu_dVds = FD2 * Qiu_dVds + FD2_dVds * T1 
        + ( 1.0 - FD2 ) * Qn0_dVds - FD2_dVds * Qn0 ;
      Qiu_dVgs = FD2 * Qiu_dVgs + FD2_dVgs * T1 
        + ( 1.0 - FD2 ) * Qn0_dVgs - FD2_dVgs * Qn0 ;
      Qiu_dT = FD2 * Qiu_dT + FD2_dT * T1 
        + ( 1.0 - FD2 ) * Qn0_dT - FD2_dT * Qn0 ;

      if ( Qiu < 0.0 ) {
        Qiu = 0.0 ;
        Qiu_dVbs = 0.0 ;
        Qiu_dVds = 0.0 ;
        Qiu_dVgs = 0.0 ;
        Qiu_dT = 0.0 ;
      }

      T1 = Qdrat ;
      Qdrat = FD2 * Qdrat + ( 1.0 - FD2 ) * 0.5e0 ;
      Qdrat_dVbs = FD2 * Qdrat_dVbs + FD2_dVbs * T1 - FD2_dVbs * 0.5e0 ;
      Qdrat_dVds = FD2 * Qdrat_dVds + FD2_dVds * T1 - FD2_dVds * 0.5e0 ;
      Qdrat_dVgs = FD2 * Qdrat_dVgs + FD2_dVgs * T1 - FD2_dVgs * 0.5e0 ;
      Qdrat_dT = FD2 * Qdrat_dT + FD2_dT * T1 - FD2_dT * 0.5e0 ;

      /* note: Lred=0 in zone-D1 */
      T1 = Lred ;
      Lred = FD2 * Lred ;
      Lred_dVbs = FD2 * Lred_dVbs + FD2_dVbs * T1 ;
      Lred_dVds = FD2 * Lred_dVds + FD2_dVds * T1 ;
      Lred_dVgs = FD2 * Lred_dVgs + FD2_dVgs * T1 ;
      Lred_dT = FD2 * Lred_dT + FD2_dT * T1 ;

    } /* end of flg_zone==2 if-block */


start_of_mobility:

    Lch = Leff - Lred ;
    if ( Lch < 1.0e-9 ) {
      fprintf ( stderr , "*** warning(HiSIM_HV): actual channel length is too small. (Lch=%e[m])\n" , Lch ) ;
      fprintf ( stderr , "                    CLM5 and/or CLM6 might be too large.\n" ) ;
      Lch = 1.0e-9 ; Lch_dVbs = Lch_dVds = Lch_dVgs = 0.0 ;
      Lch_dT = 0.0 ;
    } else { Lch_dVbs = - Lred_dVbs ; Lch_dVds = - Lred_dVds ; Lch_dVgs = - Lred_dVgs ; 
      Lch_dT = - Lred_dT ;
    }

    /*-----------------------------------------------------------*
     * Muun : universal mobility.  (CGS unit)
     *-----------------*/

    T1 = here->HSMHV_ndep_o_esi ;
    T2 = here->HSMHV_ninv_o_esi ;

    T0 = here->HSMHV_ninvd ;
    T4 = 1.0 + ( Psl - Ps0 ) * T0 ;
    T4_dVb = ( Psl_dVbs - Ps0_dVbs ) * T0 ;
    T4_dVd = ( Psl_dVds - Ps0_dVds ) * T0 ;
    T4_dVg = ( Psl_dVgs - Ps0_dVgs ) * T0 ;
    T4_dT =  ( Psl_dT - Ps0_dT )     * T0 + ( Psl - Ps0 ) * ninvd_dT ;
       
    T5     = T1 * Qbu      + T2 * Qiu ;
    T5_dVb = T1 * Qbu_dVbs + T2 * Qiu_dVbs ;
    T5_dVd = T1 * Qbu_dVds + T2 * Qiu_dVds ;
    T5_dVg = T1 * Qbu_dVgs + T2 * Qiu_dVgs ;
    T5_dT  = T1 * Qbu_dT   + T2 * Qiu_dT   ;

    T3     = T5 / T4 ;
    T3_dVb = ( - T4_dVb * T5 + T4 * T5_dVb ) / T4 / T4 ;
    T3_dVd = ( - T4_dVd * T5 + T4 * T5_dVd ) / T4 / T4 ;
    T3_dVg = ( - T4_dVg * T5 + T4 * T5_dVg ) / T4 / T4 ;
    T3_dT  = ( - T4_dT  * T5 + T4 * T5_dT  ) / T4 / T4 ;

    Eeff = T3 ;
    Eeff_dVbs = T3_dVb ;
    Eeff_dVds = T3_dVd ;
    Eeff_dVgs = T3_dVg ;
    Eeff_dT  = T3_dT ;

    T5 = Fn_Pow( Eeff , model->HSMHV_mueph0 - 1.0e0 ) ;
    T8 = T5 * Eeff ;
    T7 = Fn_Pow( Eeff , here->HSMHV_muesr - 1.0e0 ) ;
    T6 = T7 * Eeff ;
    T8_dT = model->HSMHV_mueph0 * T5 * Eeff_dT ;
    T6_dT = here->HSMHV_muesr * T7 * Eeff_dT ;


    T9 = C_QE * C_m2cm_p2 ;
    Rns = Qiu / T9 ;
    Rns_dT = Qiu_dT / T9 ;

    T1 = 1.0e0 / ( pParam->HSMHV_muecb0 + pParam->HSMHV_muecb1 * Rns / 1.0e11 ) 
      + here->HSMHV_mphn0 * T8 + T6 / pParam->HSMHV_muesr1 ;


    T1_dT =  - 1.0e0 / ( pParam->HSMHV_muecb0 + pParam->HSMHV_muecb1 * Rns / 1.0e11 )
      / ( pParam->HSMHV_muecb0 + pParam->HSMHV_muecb1 * Rns / 1.0e11 )
      * pParam->HSMHV_muecb1 * Rns_dT / 1.0e11 
      + here->HSMHV_mphn0 * T8_dT + mphn0_dT * T8 + T6_dT / pParam->HSMHV_muesr1 ;

    Muun = 1.0e0 / T1 ;
    Muun_dT = - Muun / T1 * T1_dT ;

    T1 = 1.0e0 / ( T1 * T1 ) ;
    T2 = pParam->HSMHV_muecb0 + pParam->HSMHV_muecb1 * Rns / 1.0e11 ;
    T2 = 1.0e0 / ( T2 * T2 ) ;
    T3 = here->HSMHV_mphn1 * T5 ;
    T4 = here->HSMHV_muesr * T7 / pParam->HSMHV_muesr1 ;
    T5 = - 1.0e-11 * pParam->HSMHV_muecb1 / C_QE * T2 / C_m2cm_p2 ;
    Muun_dVbs = - ( T5 * Qiu_dVbs 
                    + Eeff_dVbs * T3 + Eeff_dVbs * T4 ) * T1 ;
    Muun_dVds = - ( T5 * Qiu_dVds  
                    + Eeff_dVds * T3 + Eeff_dVds * T4 ) * T1 ;
    Muun_dVgs = - ( T5 * Qiu_dVgs  
                    + Eeff_dVgs * T3 + Eeff_dVgs * T4 ) * T1 ;

    /*  Change to MKS unit */
    Muun      /= C_m2cm_p2 ;
    Muun_dT   /= C_m2cm_p2 ;
    Muun_dVbs /= C_m2cm_p2 ;
    Muun_dVds /= C_m2cm_p2 ;
    Muun_dVgs /= C_m2cm_p2 ;

    /*-----------------------------------------------------------*
     * Mu : mobility 
     *-----------------*/
    T2 = beta * (Qn0 + small) * Lch ;


    T1 = 1.0e0 / T2 ;
    T3 = T1 * T1 ;
    T4 =  - beta * T3 ;
    T5 = T4 * Lch ;
    T6 = T4 * (Qn0 + small) ;
    T1_dVb = ( T5 * Qn0_dVbs + T6 * Lch_dVbs) ;
    T1_dVd = ( T5 * Qn0_dVds + T6 * Lch_dVds) ;
    T1_dVg = ( T5 * Qn0_dVgs + T6 * Lch_dVgs) ;
    T2_dT = beta_dT * (Qn0 + small) * Lch + beta * Qn0_dT * Lch + beta * (Qn0 + small) * Lch_dT ;
    T1_dT = - T1 / T2 * T2_dT ;

    TY = Idd * T1 ;
    TY_dVbs = Idd_dVbs * T1 + Idd * T1_dVb ;
    TY_dVds = Idd_dVds * T1 + Idd * T1_dVd ;
    TY_dVgs = Idd_dVgs * T1 + Idd * T1_dVg ;
    TY_dT = Idd_dT * T1 + Idd * T1_dT ;

    T2 = 0.2 * Vmax /  Muun ;
    T3 = - T2 / Muun ;
    T2_dVb = T3 * Muun_dVbs ;
    T2_dVd = T3 * Muun_dVds ;
    T2_dVg = T3 * Muun_dVgs ;
    T2_dT = 0.2 * ( Vmax_dT * Muun - Muun_dT * Vmax )/ ( Muun * Muun ) ;

    Ey = sqrt( TY * TY + T2 * T2 ) ;
    T4 = 1.0 / Ey ;
    Ey_dVbs = T4 * ( TY * TY_dVbs + T2 * T2_dVb ) ;
    Ey_dVds = T4 * ( TY * TY_dVds + T2 * T2_dVd ) ;
    Ey_dVgs = T4 * ( TY * TY_dVgs + T2 * T2_dVg ) ;
    Ey_dT = T4 * ( TY * TY_dT + T2 * T2_dT ) ;

    Em = Muun * Ey ;
    Em_dVbs = Muun_dVbs * Ey + Muun * Ey_dVbs ;
    Em_dVds = Muun_dVds * Ey + Muun * Ey_dVds ;
    Em_dVgs = Muun_dVgs * Ey + Muun * Ey_dVgs ;
    Em_dT = Ey * Muun_dT + Ey_dT * Muun ;
    
    T1  = Em / Vmax ;
    T1_dT = ( Em_dT * Vmax - Vmax_dT * Em ) / ( Vmax * Vmax );

    /* note: model->HSMHV_bb = 2 (electron) ;1 (hole) */
    if ( 1.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 1.0e0 + epsm10 ) {
      T3 = 1.0e0 ;
      T3_dT = 0.0e0 ;
    } else if ( 2.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 2.0e0 + epsm10 ) {
      T3 = T1 ;
      T3_dT = T1_dT ;
    } else {
      T3 = Fn_Pow( T1 , model->HSMHV_bb - 1.0e0 ) ;
      T3_dT = ( model->HSMHV_bb - 1.0e0 )* Fn_Pow( T1 , model->HSMHV_bb - 2.0e0 ) * T1_dT ;
    }
    T2 = T1 * T3 ;
    T4 = 1.0e0 + T2 ;
    T2_dT = T1 * T3_dT + T3 * T1_dT ;
    T4_dT = T2_dT ;

    if ( 1.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 1.0e0 + epsm10 ) {
      T5 = 1.0 / T4 ;
          T6 = T5 / T4 ;
      T5_dT = - T5 * T5 * T4_dT ; 
      T6_dT = T5 * T5 * ( T5_dT * T4 - T5 * T4_dT ) ;
    } else if ( 2.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 2.0e0 + epsm10 ) {
      T5 = 1.0 / sqrt( T4 ) ;
          T6 = T5 / T4 ;
      T5_dT = - 0.5e0 / ( T4 * sqrt(T4)) * T4_dT ;
      T6_dT = ( T5_dT * T4 - T5 * T4_dT ) / T4 / T4 ;
    } else {
      T6 = Fn_Pow( T4 , ( - 1.0e0 / model->HSMHV_bb - 1.0e0 ) ) ;
      T5 = T4 * T6 ;
      T6_dT =( - 1.0e0 / model->HSMHV_bb - 1.0e0 ) * Fn_Pow( T4 , ( - 1.0e0 / model->HSMHV_bb - 2.0e0 ) ) * T4_dT ;
      T5_dT = T4_dT * T6 + T4 * T6_dT ;
    }

    T7 = Muun / Vmax * T6 * T3 ;

    Mu = Muun * T5 ;
    Mu_dVbs = Muun_dVbs * T5 - T7 * Em_dVbs ;
    Mu_dVds = Muun_dVds * T5 - T7 * Em_dVds ;
    Mu_dVgs = Muun_dVgs * T5 - T7 * Em_dVgs ;
    Mu_dT = Muun_dT * T5 + Muun * T5_dT ;

/*  end_of_mobility : */

    /*-----------------------------------------------------------*
     * Ids: channel current.
     *-----------------*/
    betaWL = here->HSMHV_weff_nf * beta_inv / Lch ;
    T1 = - betaWL / Lch ;
    betaWL_dVbs = T1 * Lch_dVbs ;
    betaWL_dVds = T1 * Lch_dVds ;
    betaWL_dVgs = T1 * Lch_dVgs ;
    betaWL_dT = here->HSMHV_weff_nf * ( beta_inv_dT * Lch - beta_inv * Lch_dT ) / ( Lch * Lch ) ;

    Ids0 = betaWL * Idd * Mu ;
    T1 = betaWL * Idd ;
    T2 = Idd * Mu ;
    T3 = Mu * betaWL ;
    Ids0_dVbs = T3 * Idd_dVbs + T1 * Mu_dVbs + T2 * betaWL_dVbs ;
    Ids0_dVds = T3 * Idd_dVds + T1 * Mu_dVds + T2 * betaWL_dVds ;
    Ids0_dVgs = T3 * Idd_dVgs + T1 * Mu_dVgs + T2 * betaWL_dVgs ;
    Ids0_dT = T3 * Idd_dT + T1 * Mu_dT + T2 * betaWL_dT ;

    /* note: rpock procedure was removed. */
    if( flg_rsrd == 2 || flg_rsrd == 3 ){
      if( model->HSMHV_rd20 > 0.0 ){
        T4 = here->HSMHV_rd23 ;
        T1 = pParam->HSMHV_rd24 * ( Vgse - model->HSMHV_rd25 ) ;
        T1_dVg = pParam->HSMHV_rd24 ;       
        
        Fn_SL( T2 , T1 , T4 , delta_rd , T0 ) ;
        T2_dVg = T1_dVg * T0 ;
        T3 = T4 * ( model->HSMHV_rd20 + 1.0 ) ;
        Fn_SU( T7 , T2 , T3 , delta_rd , T0 ) ;
        T7_dVg = T2_dVg * T0 ;       

      }else{
        T7 = here->HSMHV_rd23;
        T7_dVg = 0.0e0 ;       
      }

      /* after testing we can remove Vdse_eff_dVbs, Vdse_eff_dVds, Vdse_eff_dVgs
         and Vdse_eff_dVbse, Vdse_eff_dVgse                                      */
      if (Vdse >= 0.0) {
        Vdse_eff = Vdse ;
        /* Vdse_eff_dVbs  = 0.0 ; */
        /* Vdse_eff_dVds  = 0.0 ; */
        /* Vdse_eff_dVgs  = 0.0 ; */
        /* Vdse_eff_dVbse = 0.0 ; */
        Vdse_eff_dVdse = 1.0 ;
        /* Vdse_eff_dVgse = 0.0 ; */
      } else {
        Vdse_eff = 0.0 ;
        /* Vdse_eff_dVbs  = 0.0 ; */
        /* Vdse_eff_dVds  = 0.0 ; */
        /* Vdse_eff_dVgs  = 0.0 ; */
        /* Vdse_eff_dVbse = 0.0 ; */
        Vdse_eff_dVdse = 0.0 ;
        /* Vdse_eff_dVgse = 0.0 ; */
      }
       
      /* smoothing of Ra for Vdse_eff close to zero */
      /* ... smoothing parameter is Ra_N            */
      if (Vdse_eff < Ra_N * small2) {
        Ra_alpha = pow( Ra_N+1.0 , model->HSMHV_rd21-1.0 )
                   * (Ra_N+1.0-0.5*model->HSMHV_rd21*Ra_N)
                   * pow( small2,model->HSMHV_rd21 );
        Ra_beta = 0.5*model->HSMHV_rd21
                  * pow( Ra_N+1.0 , model->HSMHV_rd21-1.0 ) / Ra_N
                  * pow( small2, model->HSMHV_rd21-2.0 );
        T1 = Ra_alpha + Ra_beta*Vdse_eff*Vdse_eff;
        T1_dVdse_eff = 2.0 * Ra_beta * Vdse_eff;
      } else {
        T1           = pow( Vdse_eff + small2 , model->HSMHV_rd21 ) ;
        T1_dVdse_eff = model->HSMHV_rd21 * pow( Vdse_eff + small2 , model->HSMHV_rd21 - 1.0 ) ;
      }

      T9           = pow( Vdse_eff + small2 , model->HSMHV_rd22d ) ;
      T9_dVdse_eff = model->HSMHV_rd22d * pow( Vdse_eff + small2 , model->HSMHV_rd22d - 1.0 ) ;

      Ra           = ( T7 * T1 + Vbse * pParam->HSMHV_rd22 * T9 ) / here->HSMHV_weff_nf ;
      Ra_dVdse_eff = ( T7 * T1_dVdse_eff + Vbse * pParam->HSMHV_rd22 * T9_dVdse_eff ) / here->HSMHV_weff_nf ;
      Ra_dVbs      =  Ra_dVdse_eff * Vdse_eff_dVbs ;
      Ra_dVds      =  Ra_dVdse_eff * Vdse_eff_dVds ;
      Ra_dVgs      =  Ra_dVdse_eff * Vdse_eff_dVgs  + T7_dVg * T1 / here->HSMHV_weff_nf ;
      Ra_dVbse     =  Ra_dVdse_eff * Vdse_eff_dVbse + pParam->HSMHV_rd22 * T9 / here->HSMHV_weff_nf ;
      Ra_dVdse     =  Ra_dVdse_eff * Vdse_eff_dVdse ;
      Ra_dVgse     =  Ra_dVdse_eff * Vdse_eff_dVgse ;

      T0 = Ra * Ids0 ;
      T0_dVb = Ra_dVbs * Ids0 + Ra * Ids0_dVbs ;
      T0_dVd = Ra_dVds * Ids0 + Ra * Ids0_dVds ;
      T0_dVg = Ra_dVgs * Ids0 + Ra * Ids0_dVgs ;
      T0_dT  =                  Ra * Ids0_dT ;

      T1 = Vds + small2 ;
      T2 = 1.0 / T1 ;
      T3 = 1.0 + T0 * T2 ;
      T3_dVb = T0_dVb * T2 ;
      T3_dVd = ( T0_dVd * T1 - T0 ) * T2 * T2 ;
      T3_dVg = T0_dVg  * T2 ;
      T3_dT = T0_dT * T2 ;

      T4 = 1.0 / T3 ;
      Ids = Ids0 * T4 ;
      T5 = T4 * T4 ;
      Ids_dVbs = ( Ids0_dVbs * T3 - Ids0 * T3_dVb ) * T5 ;
      Ids_dVds = ( Ids0_dVds * T3 - Ids0 * T3_dVd ) * T5 ;
      Ids_dVgs = ( Ids0_dVgs * T3 - Ids0 * T3_dVg ) * T5 ;
      Ids_dT = ( Ids0_dT * T3 - Ids0 * T3_dT ) * T5 ;
      Ids_dRa = - Ids * Ids / ( Vds + small ) ;

    } else {
      Ids = Ids0 ;
      Ids_dVbs = Ids0_dVbs ;
      Ids_dVds = Ids0_dVds ;
      Ids_dVgs = Ids0_dVgs ;
      Ids_dT = Ids0_dT ;
      Ra = 0.0 ;
      Ra_dVbs = Ra_dVds = Ra_dVgs = 0.0 ;
      Ra_dVbse = Ra_dVdse = Ra_dVgse = 0.0 ;
      Ids_dRa = 0.0 ;
    }
    /* just for testing  -- can be removed */
    /* if (!(ckt->CKTmode & MODEINITPRED))
      printf("rrb %e %e %e %e %e %e\n",ckt->CKTtime,here->HSMHV_mode*Vdse,Ra,Ra_dVdse,
                                       Vdse_eff,Vdse_eff_dVdse) ; */
    
    /* if ( Pds < ps_conv ) { */
    if ( Pds < 0.0 ) {
      Ids_dVbs = 0.0 ;
      Ids_dVgs = 0.0 ;
      Ids_dT = 0.0 ;
    }
  
    Ids += Gdsmin * Vds ;
    Ids_dVds += Gdsmin ;


    /*-----------------------------------------------------------*
     * STI
     *-----------------*/
    if ( model->HSMHV_coisti != 0 ) {
      /*---------------------------------------------------*
       * dVthSCSTI : Short-channel effect induced by Vds (STI).
       *-----------------*/      
      T1 = C_ESI * Cox_inv ;
      T2 = here->HSMHV_wdpl ;
      T3 =  here->HSMHV_lgatesm - model->HSMHV_parl2 ;
      T4 = 1.0 / (T3 * T3) ;
      T5 = 2.0 * (model->HSMHV_vbi - Pb20b) * T1 * T2 * T4 ;
      
      dVth0 = T5 * sqrt_Pbsum ;
      T6 = T5 * 0.5 / sqrt_Pbsum ;
      T7 = 2.0 * (model->HSMHV_vbi - Pb20b) * C_ESI * T2 * T4 * sqrt_Pbsum ;
      T8 = - 2.0 * T1 * T2 * T4 * sqrt_Pbsum ;
      dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
      dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
      dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;
      dVth0_dT = T6 * Pbsum_dT + T8 * Pb20b_dT ;

      T4 = pParam->HSMHV_scsti1 ;
      T6 = pParam->HSMHV_scsti2 ;
      T1  = T4 + T6 * Vdsz ;
      dVthSCSTI = dVth0 * T1 ;
      dVthSCSTI_dVb = dVth0_dVb * T1 + dVth0 * T6 * Vdsz_dVbs ;
      dVthSCSTI_dVd = dVth0_dVd * T1 + dVth0 * T6 * Vdsz_dVds ;
      dVthSCSTI_dVg = dVth0_dVg * T1 ;
      dVthSCSTI_dT  = dVth0_dT * T1  + dVth0 * T6 * Vdsz_dT ;

      T1 = pParam->HSMHV_vthsti - model->HSMHV_vdsti * Vds ;
      T1_dVd = - model->HSMHV_vdsti ;

      Vgssti = Vgsz - Vfb + T1 + dVthSCSTI ;
      Vgssti_dVbs = Vgsz_dVbs + dVthSCSTI_dVb ;
      Vgssti_dVds = Vgsz_dVds + T1_dVd + dVthSCSTI_dVd ;
      Vgssti_dVgs = Vgsz_dVgs + dVthSCSTI_dVg ;
      Vgssti_dT   = Vgsz_dT   + dVthSCSTI_dT ;
      
      costi0 = here->HSMHV_costi0 ;
      costi1 = here->HSMHV_costi1 ;

      costi3 = here->HSMHV_costi0_p2 * Cox_inv * Cox_inv ;
      T1 = 2.0 * here->HSMHV_costi0_p2 * Cox_inv ;
      costi3_dVb = T1 * Cox_inv_dVb ;
      costi3_dVd = T1 * Cox_inv_dVd ;
      costi3_dVg = T1 * Cox_inv_dVg ;
      costi3_dT = 2 * here->HSMHV_costi0 * here->HSMHV_costi00 * 0.5 / sqrt(here->HSMHV_beta_inv) * beta_inv_dT * Cox_inv * Cox_inv ; 
      T2 = 1.0 / costi3 ;
      costi3_dVb_c3 = costi3_dVb * T2 ;
      costi3_dVd_c3 = costi3_dVd * T2 ;
      costi3_dVg_c3 = costi3_dVg * T2 ;
/*      costi3_dT_c3 = costi3_dT * T2 ;*/

      costi4 = costi3 * beta * 0.5 ;
      costi4_dT = ( costi3_dT * beta + costi3 * beta_dT ) * 0.5 ;
      costi5 = costi4 * beta * 2.0 ;
      costi5_dT = ( costi4_dT * beta + costi4 * beta_dT ) * 2.0 ;

      T11 = beta * 0.25 ;
      T11_dT = beta_dT * 0.25 ;
      T10 = beta_inv - costi3 * T11 + Vfb - pParam->HSMHV_vthsti - dVthSCSTI + small ;
      T10_dVb = - T11 * costi3_dVb - dVthSCSTI_dVb ;
      T10_dVd = - T11 * costi3_dVd - dVthSCSTI_dVd ;
      T10_dVg = - T11 * costi3_dVg - dVthSCSTI_dVg ;
      T10_dT = beta_inv_dT - ( costi3_dT * T11 + costi3 * T11_dT ) - dVthSCSTI_dT ;

      T1 = Vgsz - T10 - psisti_dlt ;
      T1_dVb = Vgsz_dVbs - T10_dVb ;
      T1_dVd = Vgsz_dVds - T10_dVd ;
      T1_dVg = Vgsz_dVgs - T10_dVg ;
      T1_dT = - T10_dT ;
      T0 = Fn_Sgn(T10) ;
      T2 = sqrt (T1 * T1 + T0 * 4.0 * T10 * psisti_dlt) ;
      T3 = T10 + 0.5 * (T1 + T2) - Vfb + pParam->HSMHV_vthsti + dVthSCSTI - Vbsz ;
      T3_dVb = T10_dVb + 0.5 * (T1_dVb + (T1 * T1_dVb + T0 * 2.0 * T10_dVb * psisti_dlt) / T2) 
        + dVthSCSTI_dVb - Vbsz_dVbs ;
      T3_dVd = T10_dVd + 0.5 * (T1_dVd + (T1 * T1_dVd + T0 * 2.0 * T10_dVd * psisti_dlt) / T2) 
        + dVthSCSTI_dVd - Vbsz_dVds ;
      T3_dVg = T10_dVg + 0.5 * (T1_dVg + (T1 * T1_dVg + T0 * 2.0 * T10_dVg * psisti_dlt) / T2) 
        + dVthSCSTI_dVg ;
      T3_dT = T10_dT + 0.5 * (T1_dT + (T1 * T1_dT + T0 * 2.0 * T10_dT * psisti_dlt) / T2) 
        + dVthSCSTI_dT  - Vbsz_dT ;

      T4 = beta * T3 - 1.0 ;
      T4_dT = beta_dT * T3 + beta * T3_dT ; 
      T5 = 4.0 / costi5 ;
      T5_dT = - 4.0 / ( costi5 * costi5 ) * costi5_dT ;
      T1 = 1.0 + T4 * T5 ;
      T6 = beta * T5 ;
      T7 = T4 * T5 ;
      T1_dVb = (T6 * T3_dVb - T7 * costi3_dVb_c3) ;
      T1_dVd = (T6 * T3_dVd - T7 * costi3_dVd_c3) ;
      T1_dVg = (T6 * T3_dVg - T7 * costi3_dVg_c3) ;
      T1_dT = T4_dT * T5 + T4 * T5_dT ; 
      Fn_SZ( T1 , T1, 1.0e-2, T2) ;
      T1_dVb *= T2 ;
      T1_dVd *= T2 ;
      T1_dVg *= T2 ;
      T1_dT *= T2 ; 
      costi6 = sqrt(T1) ;
      costi6_dT = 0.5 / sqrt(T1) * T1_dT ;
      T0 = costi4 * (1.0 - costi6) ;
      T0_dT = costi4_dT * (1.0 - costi6) + costi4 * ( - costi6_dT ) ;
      Psasti = Vgssti + T0 ;
      T2 = 0.5 * costi4 / costi6 ;
      Psasti_dVbs = Vgssti_dVbs + costi3_dVb_c3 * T0 - T2 * T1_dVb ;
      Psasti_dVds = Vgssti_dVds + costi3_dVd_c3 * T0 - T2 * T1_dVd ;
      Psasti_dVgs = Vgssti_dVgs + costi3_dVg_c3 * T0 - T2 * T1_dVg ;
      Psasti_dT = Vgssti_dT + T0_dT ;

      T0 = 1.0 / (beta + 2.0 / (Vgssti + small)) ;
      T0_dT = - 1.0 / ((beta + 2.0 / (Vgssti + small)) * (beta + 2.0 / (Vgssti + small))) * ( beta_dT - 2 / ((Vgssti + small) * (Vgssti + small)) * Vgssti_dT );
      Psbsti = log (1.0 / costi1 / costi3 * (Vgssti * Vgssti)) * T0 ;
      T1 = 1 / costi1 / costi3 * (Vgssti * Vgssti) ;

      costi1_dT = 2 * here->HSMHV_nin * Nin_dT * here->HSMHV_nsti_p2 ;

      T1_dT = ( - 1 / costi1 / costi1 * costi1_dT / costi3 - 1 / costi3 / costi3 * costi3_dT / costi1 ) * Vgssti * Vgssti + 1 / costi1 / costi3 * 2 * Vgssti * Vgssti_dT ;
      T2 = 2.0 * T0 / (Vgssti + small) ;
      T3 = Psbsti / (Vgssti + small) ;
      Psbsti_dVbs = T2 * (Vgssti_dVbs - 0.5 * costi3_dVb_c3 * Vgssti 
                                + T3 * Vgssti_dVbs ) ;
      Psbsti_dVds = T2 * (Vgssti_dVds - 0.5 * costi3_dVd_c3 * Vgssti 
                                + T3 * Vgssti_dVds ) ;
      Psbsti_dVgs = T2 * (Vgssti_dVgs - 0.5 * costi3_dVg_c3 * Vgssti 
                                + T3 * Vgssti_dVgs ) ;
      Psbsti_dT = 1 / T1 * T1_dT * T0 + log( T1 ) * T0_dT ;
      
      Psab = Psbsti - Psasti - sti2_dlt ;
      Psab_dVbs = Psbsti_dVbs - Psasti_dVbs ;
      Psab_dVds = Psbsti_dVds - Psasti_dVds ;
      Psab_dVgs = Psbsti_dVgs - Psasti_dVgs ;
      Psab_dT = Psbsti_dT - Psasti_dT ; 
      T0 = sqrt (Psab * Psab + 4.0 * sti2_dlt * Psbsti) ;
      Psti = Psbsti - 0.5 * (Psab + T0) ;
      T1 = 1.0 / T0 ;
      Psti_dVbs = Psbsti_dVbs 
        - 0.5 * ( Psab_dVbs
                + ( Psab * Psab_dVbs + 2.0 * sti2_dlt * Psbsti_dVbs ) * T1 ) ;
      Psti_dVds = Psbsti_dVds 
        - 0.5 * ( Psab_dVds
                + ( Psab * Psab_dVds + 2.0 * sti2_dlt * Psbsti_dVds ) * T1 ) ;
      Psti_dVgs = Psbsti_dVgs 
        - 0.5 * ( Psab_dVgs
                + ( Psab * Psab_dVgs + 2.0 * sti2_dlt * Psbsti_dVgs ) * T1 ) ;
      Psti_dT = Psbsti_dT 
        - 0.5 * ( Psab_dT
                + ( Psab * Psab_dT + 2.0 * sti2_dlt * Psbsti_dT ) * T1 ) ;

      T0 = costi1 * exp (beta * Psti) ;
      T0_dT = costi1_dT * exp(beta * Psti) + costi1 * exp(beta * Psti) * ( beta_dT * Psti + beta * Psti_dT ) ;
      T1 = beta * (Psti - Vbsz) - 1.0 + T0 ;
      T1_dVb = beta * ((Psti_dVbs - Vbsz_dVbs) + T0 * Psti_dVbs) ;
      T1_dVd = beta * ((Psti_dVds - Vbsz_dVds) + T0 * Psti_dVds) ;
      T1_dVg = beta * (Psti_dVgs + T0 * Psti_dVgs) ;
      T1_dT = beta_dT * (Psti - Vbsz) + beta * (Psti_dT - Vbsz_dT) + T0_dT ;
      Fn_SZ ( T1 , T1, 1.0e-2, T0) ;
          T1 += epsm10 ;
      T1_dVb *= T0 ;
      T1_dVd *= T0 ;
      T1_dVg *= T0 ;
      T1_dT *= T0 ;
      sq1sti = sqrt (T1);
      T2 = 0.5 / sq1sti ;
      sq1sti_dVbs = T2 * T1_dVb ;
      sq1sti_dVds = T2 * T1_dVd ;
      sq1sti_dVgs = T2 * T1_dVg ;
      sq1sti_dT = T2 * T1_dT ;

      T1 = beta * (Psti - Vbsz) - 1.0;
      T1_dVb = beta * (Psti_dVbs - Vbsz_dVbs) ;
      T1_dVd = beta * (Psti_dVds - Vbsz_dVds) ;
      T1_dVg = beta * Psti_dVgs ;
      T1_dT = beta_dT * ( Psti - Vbsz ) + beta * (Psti_dT - Vbsz_dT) ;
      Fn_SZ( T1 , T1, 1.0e-2, T0) ;
          T1 += epsm10 ;
      T1_dVb *= T0 ;
      T1_dVd *= T0 ;
      T1_dVg *= T0 ;
      T1_dT *= T0 ;
      sq2sti = sqrt (T1);
      T2 = 0.5 / sq2sti ;
      sq2sti_dVbs = T2 * T1_dVb ;
      sq2sti_dVds = T2 * T1_dVd ;
      sq2sti_dVgs = T2 * T1_dVg ;
      sq2sti_dT = T2 * T1_dT ; 

      Qn0sti = costi0 * (sq1sti - sq2sti) ;
      Qn0sti_dVbs = costi0 * (sq1sti_dVbs - sq2sti_dVbs) ;
      Qn0sti_dVds = costi0 * (sq1sti_dVds - sq2sti_dVds) ;
      Qn0sti_dVgs = costi0 * (sq1sti_dVgs - sq2sti_dVgs) ;
      Qn0sti_dT = costi0 * (sq1sti_dT - sq2sti_dT) + here->HSMHV_costi00 * 0.5 / sqrt( here->HSMHV_beta_inv ) * beta_inv_dT * (sq1sti - sq2sti) ;

      /* T1: Vdsatsti */
      T1 = Psasti - Psti ;
      T1_dVb = Psasti_dVbs - Psti_dVbs ;
      T1_dVd = Psasti_dVds - Psti_dVds ;
      T1_dVg = Psasti_dVgs - Psti_dVgs ;
      T1_dT = Psasti_dT - Psti_dT ;

      Fn_SZ( T1 , T1 , 1.0e-1 , T2 ) ;
      T1_dVb *= T2 ;
      T1_dVd *= T2 ;
      T1_dVg *= T2 ;
      T1_dT *= T2 ;

      TX = Vds / T1 ;
      T2 = 1.0 / ( T1 * T1 ) ;
      TX_dVbs = T2 * ( - Vds * T1_dVb ) ;
      TX_dVds = T2 * ( T1 - Vds * T1_dVd ) ;
      TX_dVgs = T2 * ( - Vds * T1_dVg ) ;
      TX_dT = T2 * ( - Vds * T1_dT ) ;

      Fn_CP( TY , TX , 1.0 , 4 , T2 ) ;
      TY_dVbs = T2 * TX_dVbs ;
      TY_dVds = T2 * TX_dVds ;
      TY_dVgs = T2 * TX_dVgs ;
      TY_dT = T2 * TX_dT ; 

      costi7 = 2.0 * here->HSMHV_wsti * here->HSMHV_nf * beta_inv ;
      costi7_dT = 2.0 * here->HSMHV_wsti * here->HSMHV_nf * beta_inv_dT ;
      T1 = Lch ;
      Idssti = costi7 * Mu * Qn0sti * TY / T1 ;
      T3 = 1.0 / T1 ;
      T4 = Mu * Qn0sti * TY / T1 / T1 ;
      T5 = Mu * Qn0sti ;
      Idssti_dVbs = costi7 * (((Mu_dVbs * Qn0sti + Mu * Qn0sti_dVbs) * TY
                               +  T5 * TY_dVbs ) * T3 
                              - Lch_dVbs * T4 ) ;
      Idssti_dVds = costi7 * (((Mu_dVds * Qn0sti + Mu * Qn0sti_dVds) * TY 
                               +  T5 * TY_dVds ) * T3 
                              - Lch_dVds * T4 ) ;
      Idssti_dVgs = costi7 * (((Mu_dVgs * Qn0sti + Mu * Qn0sti_dVgs) * TY
                               +  T5 * TY_dVgs ) * T3 
                              - Lch_dVgs * T4 ) ;
      Idssti_dT = costi7 * (((Mu_dT * Qn0sti + Mu * Qn0sti_dT) * TY
                               +  T5 * TY_dT ) * T3 
                              - Lch_dT * T4 ) + costi7_dT * Mu * Qn0sti * TY / T1 ;

      Ids = Ids + Idssti ;
      Ids_dVbs = Ids_dVbs + Idssti_dVbs ;
      Ids_dVds = Ids_dVds + Idssti_dVds ;
      Ids_dVgs = Ids_dVgs + Idssti_dVgs ;
      Ids_dT  = Ids_dT  + Idssti_dT ;
    
    }



  /*----------------------------------------------------------*
   * induced gate noise. ( Part 1/3 )
   *----------------------*/
  if ( model->HSMHV_coign != 0 && model->HSMHV_cothrml != 0 ) {
      kusai00 = VgVt * VgVt ;
      kusaidd = 2.0e0 * beta_inv * Cox_inv * Idd ;
      kusaiL = kusai00 - kusaidd ;
      Fn_SZ( kusai00 , kusai00 , 1.0e-3 , T0 ) ;
      Fn_SZ( kusaiL , kusaiL , 1.0e-3 , T0 ) ;
      kusai00L = kusai00 - kusaiL ;
      if ( Qn0 < epsm10 || kusai00L < epsm10 ) flg_ign = 0 ;
      else flg_ign = 1 ;
  }

  /*-----------------------------------------------------------*
   * End of PART-1. (label) 
   *-----------------*/
 end_of_part_1: 

  /*----------------------------------------------------------*
   * Evaluate integrated charges in unit [C].
   *----------------------*/
    
    T1 = - here->HSMHV_weffcv_nf * Leff ;

    Qb = T1 * Qbu ;
    Qb_dVbs = T1 * Qbu_dVbs ;
    Qb_dVds = T1 * Qbu_dVds ;
    Qb_dVgs = T1 * Qbu_dVgs ;
    Qb_dT  = T1 * Qbu_dT ;
    
    Qi = T1 * Qiu ;
    Qi_dVbs = T1 * Qiu_dVbs ;
    Qi_dVds = T1 * Qiu_dVds ;
    Qi_dVgs = T1 * Qiu_dVgs ;
    Qi_dT  = T1 * Qiu_dT ;
 
    Qd = Qi * Qdrat ;
    Qd_dVbs = Qi_dVbs * Qdrat + Qi * Qdrat_dVbs ;
    Qd_dVds = Qi_dVds * Qdrat + Qi * Qdrat_dVds ;
    Qd_dVgs = Qi_dVgs * Qdrat + Qi * Qdrat_dVgs ;
    Qd_dT  = Qi_dT * Qdrat + Qi * Qdrat_dT ;

    
    /*-----------------------------------------------------------*
     * Modified potential for symmetry. 
     *-----------------*/
    T1 =  ( Vds - Pds ) / 2 ;
    Fn_SymAdd( Pzadd , T1 , model->HSMHV_pzadd0 , T2 ) ;
    T2 /= 2 ;
    Pzadd_dVbs = T2 * ( - Pds_dVbs ) ;
    Pzadd_dVds = T2 * ( 1.0 - Pds_dVds ) ;
    Pzadd_dVgs = T2 * ( - Pds_dVgs ) ;
    Pzadd_dT = T2 * ( -Pds_dT );

    
    if ( Pzadd  < epsm10 ) {
      Pzadd = epsm10 ;
      Pzadd_dVbs = 0.0 ;
      Pzadd_dVds = 0.0 ;
      Pzadd_dVgs = 0.0 ;
      Pzadd_dT = 0.0 ;
    }

    Ps0z = Ps0 + Pzadd ;
    Ps0z_dVbs = Ps0_dVbs + Pzadd_dVbs ;
    Ps0z_dVds = Ps0_dVds + Pzadd_dVds ;
    Ps0z_dVgs = Ps0_dVgs + Pzadd_dVgs ;
    Ps0z_dT = Ps0_dT + Pzadd_dT ;


  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-2: Substrate / gate / leak currents
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 
  /*-----------------------------------------------------------*
   * Isub : substrate current induced by impact ionization.
   *-----------------*/

  if ( flg_noqi == 1 || model->HSMHV_coisub == 0 ) {
   /* Accumulation zone or nonconductive case, in which Ids==0. */
    Isub = 0.0e0 ;
    Isub_dVbs = Isub_dVds = Isub_dVgs = 0.0e0 ;
    Isub_dT = 0.0;
  } else {
   /*-------------------------------------------*
    * Conductive case. 
    *-----------------*/
    if ( pParam->HSMHV_sub1 > 0.0e0 && pParam->HSMHV_vmax > 0.0e0 ) {
      T0 = here->HSMHV_vg2const ;
      T1 = T0 * Vgp ; 
      T1_dVd = T0 * Vgp_dVds ;
      T1_dVg = T0 * Vgp_dVgs ;
      T1_dVb = T0 * Vgp_dVbs ;
      T1_dT  = T0 * Vgp_dT   ;

      T7 = Cox0 * Cox0 ;
      T8 = here->HSMHV_qnsub_esi ;
      T3 = T8 / T7 ;

      T9 = 2.0 / T8 ;
      T4 = 1.0e0 + T9 * T7 ;

      T2 = here->HSMHV_xvbs ;
      T5 = T1 - beta_inv - T2 * Vbsz ;
      T5_dVd = T1_dVd - T2 * Vbsz_dVds;
      T5_dVg = T1_dVg ;
      T5_dVb = T1_dVb - T2 * Vbsz_dVbs;
      T5_dT  = - beta_inv_dT + T1_dT  - T2 * Vbsz_dT ;

      T6 = T4 * T5 ;
      T6_dVd = T4 * T5_dVd ;
      T6_dVg = T4 * T5_dVg ;
      T6_dVb = T4 * T5_dVb ;
      T6_dT = T4 * T5_dT ;
      Fn_SZ( T6 , T6, 1.0e-3, T9) ;
      T6 += small ;
      T6_dVd *= T9 ;
      T6_dVg *= T9 ;
      T6_dVb *= T9 ;
      T6_dT *= T9 ;
      T6 = sqrt( T6 ) ;
      T9 = 0.5 / T6 ;
      T6_dVd = T9 * T6_dVd ;
      T6_dVg = T9 * T6_dVg ;
      T6_dVb = T9 * T6_dVb ;
      T6_dT = T9 * T6_dT ;

      Psislsat = T1 + T3 * ( 1.0 - T6 ) ;
      Psislsat_dVd = T1_dVd - T3 * T6_dVd ;
      Psislsat_dVg = T1_dVg - T3 * T6_dVg ;
      Psislsat_dVb = T1_dVb - T3 * T6_dVb ;
      Psislsat_dT  = T1_dT  - T3 * T6_dT ;

      T2 = here->HSMHV_lgate / (here->HSMHV_xgate + here->HSMHV_lgate) ;

      Psisubsat = pParam->HSMHV_svds * Vdsz + Ps0z - T2 * Psislsat ;      
      Psisubsat_dVd = pParam->HSMHV_svds * Vdsz_dVds + Ps0z_dVds - T2 * Psislsat_dVd ; 
      Psisubsat_dVg = Ps0z_dVgs - T2 * Psislsat_dVg ; 
      Psisubsat_dVb =  pParam->HSMHV_svds * Vdsz_dVbs + Ps0z_dVbs - T2 * Psislsat_dVb ;
      Psisubsat_dT  =  pParam->HSMHV_svds * Vdsz_dT   + Ps0z_dT   - T2 * Psislsat_dT  ;
      Fn_SZ( Psisubsat , Psisubsat, 1.0e-3, T9 ) ; 
      Psisubsat += small ;
      Psisubsat_dVd *= T9 ;
      Psisubsat_dVg *= T9 ;
      Psisubsat_dVb *= T9 ;       
      Psisubsat_dT *= T9 ;

      T5 = here->HSMHV_xsub1 ;
      T6 = here->HSMHV_xsub2 ;
      T2 = exp( - T6 / Psisubsat ) ;
      T3 = T2 * T6 / ( Psisubsat * Psisubsat ) ;
      T2_dVd = T3 * Psisubsat_dVd ;
      T2_dVg = T3 * Psisubsat_dVg ;
      T2_dVb = T3 * Psisubsat_dVb ;
      T2_dT = T3 * Psisubsat_dT ;

      Isub = T5 * Psisubsat * Ids * T2 ;
      Isub_dVds =  T5 * ( Psisubsat_dVd * Ids * T2
                        + Psisubsat * Ids_dVds * T2 
                        + Psisubsat * Ids * T2_dVd ) ;
      Isub_dVgs =  T5 * ( Psisubsat_dVg * Ids * T2 
                        + Psisubsat * Ids_dVgs * T2
                        + Psisubsat * Ids * T2_dVg ) ;
      Isub_dVbs =  T5 * ( Psisubsat_dVb * Ids * T2
                        + Psisubsat * Ids_dVbs * T2 
                        + Psisubsat * Ids * T2_dVb ) ;
      Isub_dT =  T5 * ( Psisubsat_dT * Ids * T2
                        + Psisubsat * Ids_dT * T2 
                        + Psisubsat * Ids * T2_dT ) ;

    } else {
      Isub = 0.0e0 ;
      Isub_dVbs = Isub_dVds = Isub_dVgs = 0.0e0 ;
      Isub_dT  = 0.0e0 ;
    } /* end of if ( pParam->HSMHV_sub1 ... ) else block. */

    /*---------------------------------------------------*
     * Impact-Ionization Current in the Drift Region
     *-----------------*/
    T8 = here->HSMHV_subld1 ;
    if ( T8 != 0.0 ) {
      T0 = ( Vdse + Ps0 - Psdl  ) ;
      T0_dVb = Ps0_dVbs - Psdl_dVbs ;
      T0_dVd = Ps0_dVds - Psdl_dVds ;
      T0_dVg = Ps0_dVgs - Psdl_dVgs ;
      T0_dT =  Ps0_dT   - Psdl_dT   ;
 
      Fn_SZ( T0, T0, 1e-6, T1 ) ;
      T0_dVb *= T1 ;
      T0_dVd *= T1 ;
      T0_dVg *= T1 ;
      T0_dT  *= T1 ;
               
      T1 = sqrt( VgVt + small  ) ;
      T3 = 1.0 / ( 2.0 * T1 ) ;
      T1_dVb = VgVt_dVbs * T3  ;
      T1_dVd = VgVt_dVds * T3  ;
      T1_dVg = VgVt_dVgs * T3  ;
      T1_dT  = VgVt_dT * T3  ;

      T4 = 1.0 / ( T0 * T1 ) ;
      T7 = Ldrift * hereMKS->HSMHV_subld2 ;
      T2= exp( - T7 * T4 );
      T6 = T7 * T2 * T4 * T4 ;
      T2_dVb = ( T0_dVb * T1 + T0 * T1_dVb ) * T6  ;
      T2_dVd = ( T0_dVd * T1 + T0 * T1_dVd ) * T6  ;
      T2_dVg = ( T0_dVg * T1 + T0 * T1_dVg ) * T6  ;
      T2_dT = ( T0_dT * T1 + T0 * T1_dT ) * T6  ;
      T2_dVdse = T1 * T6 ;

      T5 = T8 * Ids * T0 * T2 ;
      T5_dVb = T8 * ( Ids_dVbs * T0 * T2 + Ids * T0_dVb * T2 + Ids * T0 * T2_dVb ) ;
      T5_dVd = T8 * ( Ids_dVds * T0 * T2 + Ids * T0_dVd * T2 + Ids * T0 * T2_dVd ) ;
      T5_dVg = T8 * ( Ids_dVgs * T0 * T2 + Ids * T0_dVg * T2 + Ids * T0 * T2_dVg ) ;
      T5_dT = T8 * ( Ids_dT * T0 * T2 + Ids * T0_dT * T2 + Ids * T0 * T2_dT ) ;
      T5_dVdse = T8 * ( Ids_dRa * Ra_dVdse * T0 * T2 + Ids * T2 + Ids * T0 * T2_dVdse ) ;

      Isub += T5 ;
      Isub_dVgs += T5_dVg ;
      Isub_dVds += T5_dVd ;
      Isub_dVbs += T5_dVb ;
      Isub_dT += T5_dT ; 
      Isub_dVdse += T5_dVdse ;
    }

  }


    /*---------------------------------------------------*
     * Impact-Ionization Induced Bulk Potential Change (IBPC)
     *-----------------*/
  if ( flg_noqi == 0 && Isub > 0e0 && pParam->HSMHV_ibpc1 != 0e0 ) {

    /* delta Vbs */
    T0 = 1e0 + pParam->HSMHV_ibpc2 * dVth ;
    dVbsIBPC = pParam->HSMHV_ibpc1 * T0 * Isub ;
    dVbsIBPC_dVbs = pParam->HSMHV_ibpc1 * ( pParam->HSMHV_ibpc2 * dVth_dVb * Isub + T0 * Isub_dVbs ) ;
    dVbsIBPC_dVds = pParam->HSMHV_ibpc1 * ( pParam->HSMHV_ibpc2 * dVth_dVd * Isub + T0 * Isub_dVds ) ;
    dVbsIBPC_dVgs = pParam->HSMHV_ibpc1 * ( pParam->HSMHV_ibpc2 * dVth_dVg * Isub + T0 * Isub_dVgs ) ;
    dVbsIBPC_dT = pParam->HSMHV_ibpc1 * ( pParam->HSMHV_ibpc2 * dVth_dT * Isub + T0 * Isub_dT ) ;

    /* dG3 & dG4 */
    T10 = 1e0 / Xi0 ;
    T10_dT = - T10 / Xi0 * Xi0_dT ;
    T1 = beta * dVbsIBPC * T10 ;
    T11 = T10 * T10;
    T1_dVb = beta * ( dVbsIBPC_dVbs * Xi0 - dVbsIBPC * Xi0_dVbs ) * T11 ;
    T1_dVd = beta * ( dVbsIBPC_dVds * Xi0 - dVbsIBPC * Xi0_dVds ) * T11 ;
    T1_dVg = beta * ( dVbsIBPC_dVgs * Xi0 - dVbsIBPC * Xi0_dVgs ) * T11 ;
    T1_dT = beta_dT * dVbsIBPC * T10 + beta * dVbsIBPC_dT * T10 + beta * dVbsIBPC * T10_dT ;


    T10 = 1e0 / Xil ;
    T10_dT = - T10 / Xil * Xil_dT ;
    T2 = beta * dVbsIBPC * T10 ;
    T11 = T10 * T10;
    T2_dVb = beta * ( dVbsIBPC_dVbs * Xil - dVbsIBPC * Xil_dVbs ) * T11 ;
    T2_dVd = beta * ( dVbsIBPC_dVds * Xil - dVbsIBPC * Xil_dVds ) * T11 ;
    T2_dVg = beta * ( dVbsIBPC_dVgs * Xil - dVbsIBPC * Xil_dVgs ) * T11 ;
    T2_dT = beta_dT * dVbsIBPC * T10 + beta * dVbsIBPC_dT * T10 + beta * dVbsIBPC * T10_dT ;


    dG3 = cnst0 * ( Xilp32 * T2 - Xi0p32 * T1 ) ;
    dG3_dVbs = cnst0 * ( Xilp32_dVbs * T2 + Xilp32 * T2_dVb - Xi0p32_dVbs * T1 - Xi0p32 * T1_dVb ) ;
    dG3_dVds = cnst0 * ( Xilp32_dVds * T2 + Xilp32 * T2_dVd - Xi0p32_dVds * T1 - Xi0p32 * T1_dVd ) ;
    dG3_dVgs = cnst0 * ( Xilp32_dVgs * T2 + Xilp32 * T2_dVg - Xi0p32_dVgs * T1 - Xi0p32 * T1_dVg ) ;
    dG3_dT =  cnst0 * ( Xilp32_dT * T2 + Xilp32 * T2_dT - Xi0p32_dT * T1 - Xi0p32 * T1_dT ) 
      + cnst0_dT * ( Xilp32 * T2 - Xi0p32 * T1 ) ;

    dG4 = cnst0 * 0.5 * ( - Xilp12 * T2 + Xi0p12 * T1 ) ;
    dG4_dVbs = cnst0 * 0.5 * ( - Xilp12_dVbs * T2 - Xilp12 * T2_dVb + Xi0p12_dVbs * T1 + Xi0p12 * T1_dVb ) ;
    dG4_dVds = cnst0 * 0.5 * ( - Xilp12_dVds * T2 - Xilp12 * T2_dVd + Xi0p12_dVds * T1 + Xi0p12 * T1_dVd ) ;
    dG4_dVgs = cnst0 * 0.5 * ( - Xilp12_dVgs * T2 - Xilp12 * T2_dVg + Xi0p12_dVgs * T1 + Xi0p12 * T1_dVg ) ;
    dG4_dT   = cnst0 * 0.5 * ( - Xilp12_dT   * T2 - Xilp12 * T2_dT  + Xi0p12_dT   * T1 + Xi0p12 * T1_dT  )
      + cnst0_dT *0.5 *  ( - Xilp12 * T2 + Xi0p12 * T1 ) ;

    /* Add IBPC current into Ids */
    dIdd = dG3 + dG4 ;
    dIdd_dVbs = dG3_dVbs + dG4_dVbs ;
    dIdd_dVds = dG3_dVds + dG4_dVds ;
    dIdd_dVgs = dG3_dVgs + dG4_dVgs ;
    dIdd_dT = dG3_dT + dG4_dT ;

    IdsIBPC = betaWL * dIdd * Mu ;
    IdsIBPC_dVbs = betaWL * ( Mu * dIdd_dVbs + dIdd * Mu_dVbs ) + betaWL_dVbs * Mu * dIdd ;
    IdsIBPC_dVds = betaWL * ( Mu * dIdd_dVds + dIdd * Mu_dVds ) + betaWL_dVds * Mu * dIdd ;
    IdsIBPC_dVgs = betaWL * ( Mu * dIdd_dVgs + dIdd * Mu_dVgs ) + betaWL_dVgs * Mu * dIdd ;
    IdsIBPC_dT  = betaWL * ( Mu * dIdd_dT   + dIdd * Mu_dT   ) + betaWL_dT   * Mu * dIdd ;


  } /* End if (IBPC) */

  T3 = 1 / TTEMP ;
  T0 =- model->HSMHV_igtemp2 * T3 * T3   
          - 2 * model->HSMHV_igtemp3 * T3 * T3 * T3   ;
  Egp12_dT = 0.5 * T0 / Egp12; 
  Egp32_dT = 1.5 * T0 * Egp12;

  /*-----------------------------------------------------------*
   * Igate : Gate current induced by tunneling.
   *-----------------*/
  if ( model->HSMHV_coiigs != 0 ) {
    /* Igate */
    if ( flg_noqi == 0 ) {
      Psdlz = Ps0z + Vdsz - epsm10 ;
      Psdlz_dVbs = Ps0z_dVbs + Vdsz_dVbs ;
      Psdlz_dVds = Ps0z_dVds + Vdsz_dVds ;
      Psdlz_dVgs = Ps0z_dVgs ;
      Psdlz_dT   = Ps0z_dT ;

      T1 = Vgsz - Vfb + modelMKS->HSMHV_gleak4 * (dVth - dPpg) * Leff - Psdlz * pParam->HSMHV_gleak3 ;
      T3 = modelMKS->HSMHV_gleak4 * Leff ;
      T1_dVg = Vgsz_dVgs + T3 * (dVth_dVg - dPpg_dVg) - Psdlz_dVgs * pParam->HSMHV_gleak3 ;
      T1_dVd = Vgsz_dVds + T3 * (dVth_dVd - dPpg_dVd)  - Psdlz_dVds * pParam->HSMHV_gleak3 ;
      T1_dVb = Vgsz_dVbs + T3 * ( dVth_dVb - dPpg_dVb )  - Psdlz_dVbs * pParam->HSMHV_gleak3 ;
      T1_dT  = Vgsz_dT   + T3 * ( dVth_dT - dPpg_dT )  - Psdlz_dT * pParam->HSMHV_gleak3 ;

      T3 = 2.0 * T1 ;
      T1_dVg = T3 * T1_dVg ;
      T1_dVd = T3 * T1_dVd ;
      T1_dVb = T3 * T1_dVb ;
      T1_dT = T3 * T1_dT ;
      T1 *= T1 ;
      
      T3 = 1.0 / Tox0 ;
      T2 = T1 * T3 ;
      T2_dVg = (T1_dVg ) * T3 ;
      T2_dVd = (T1_dVd ) * T3 ;
      T2_dVb = (T1_dVb ) * T3 ;
      T2_dT = T1_dT * T3 ;
     
      T3 = 1.0 / modelMKS->HSMHV_gleak5 ;
      T7 = 1.0 + Ey * T3 ; 
      T7_dVg = Ey_dVgs * T3 ;
      T7_dVd = Ey_dVds * T3 ;
      T7_dVb = Ey_dVbs * T3 ;
      T7_dT = Ey_dT * T3 ;

      Etun = T2 * T7 ;
      Etun_dVgs = T2_dVg * T7 + T7_dVg * T2 ;
      Etun_dVds = T2_dVd * T7 + T7_dVd * T2 ;
      Etun_dVbs = T2_dVb * T7 + T7_dVb * T2 ;
      Etun_dT = T2_dT * T7 + T7_dT * T2 ;

      Fn_SZ( Etun , Etun , igate_dlt , T5 ) ;
      Etun_dVgs *= T5 ;
      Etun_dVds *= T5 ;
      Etun_dVbs *= T5 ;
      Etun_dT *= T5 ;

      Fn_SZ( T3 , Vgsz , 1.0e-3 , T4 ) ;
      T3 -= model->HSMHV_vzadd0 ;
      T3_dVb = 0.5 * (Vgsz_dVbs + Vgsz * Vgsz_dVbs/TMF2);
      TX = T3 / cclmmdf ; 
      TX_dVbs = T3_dVb / cclmmdf ; 
      T2 = 1.0 +  TX * TX ;
      T2_dVb = 2 * TX_dVbs * TX ;
      T1 = 1.0 - 1.0 / T2 ;
      T1_dVb = T2_dVb / T2 / T2 ;
      T1_dVg = 2.0 * TX * T4 / ( T2 * T2 * cclmmdf ) ;
      T1_dVd = T1_dVg * Vgsz_dVds ;
      Etun_dVgs = T1 * Etun_dVgs + Etun * T1_dVg ;
      Etun_dVds = T1 * Etun_dVds + Etun * T1_dVd ;
      Etun_dVbs = Etun_dVbs * T1 + Etun * T1_dVb ;
      Etun_dT *= T1 ;
      Etun *= T1 ;

      T0 = Leff * here->HSMHV_weff_nf ;
      T7 = modelMKS->HSMHV_gleak7 / (modelMKS->HSMHV_gleak7 + T0) ;
      
      T6 = pParam->HSMHV_gleak6 ;
      T9 = T6 / (T6 + Vdsz) ;
      T9_dVb = - T9 / (T6 + Vdsz) * Vdsz_dVbs ;
      T9_dVd = - T9 / (T6 + Vdsz) * Vdsz_dVds ;
      
      T4 = 1 / (Etun + small ) ;
      T1 = - pParam->HSMHV_gleak2 * Egp32 * T4 ;
      T3 =  pParam->HSMHV_gleak2 * T4 * T4;
      T1_dT =  T3 * (Egp32 * Etun_dT  - Egp32_dT * (Etun + small ))  ;

      if ( T1 < - EXP_THR ) {
        Igate = 0.0 ;
        Igate_dVbs = Igate_dVds = Igate_dVgs = Igate_dT = 0.0 ;
      } else {
        T2 = exp ( T1 ) ;
        T2_dT = T1_dT * T2 ;

        T3 = pParam->HSMHV_gleak1 / Egp12 * C_QE * T0 ;
        T3_dT = - Egp12_dT * pParam->HSMHV_gleak1 / Egp12 / Egp12 * C_QE * T0 ;
  
	T5 = 1 / cnst0 ;
	T6 =  sqrt ((Qiu + Cox0 * VgVt_small )* T5 ) ;
	T6_dT = ( ( ( cnst0 * Qiu_dT - cnst0_dT * ( Qiu + Cox0 * VgVt_small ) ) * T5 * T5 ) ) / T6 * 0.5 ;
        T4 =  T2 * T3 * T6 ;
        T4_dT =  T2_dT * T3 * T6 + T2 * T3_dT * T6 +  T2 * T3 * T6_dT;
        T5 = T4 * Etun ;
        T6 = 0.5 * Etun / (Qiu + Cox0 * VgVt_small ) ;
        T10 = T5 * Etun ;
        T10_dVb = T5 * (2.0 * Etun_dVbs - T1 * Etun_dVbs + T6 * Qiu_dVbs) ;
        T10_dVd = T5 * (2.0 * Etun_dVds - T1 * Etun_dVds + T6 * Qiu_dVds) ;
        T10_dVg = T5 * (2.0 * Etun_dVgs - T1 * Etun_dVgs + T6 * Qiu_dVgs) ;
	T10_dT = 2 * T5 * Etun_dT + T4_dT * Etun * Etun ;
          
        Igate = T7 * T9 * T10 ;
        Igate_dVbs = T7 * (T9 * T10_dVb + T9_dVb * T10) ;
        Igate_dVds = T7 * (T9_dVd * T10 + T9 * T10_dVd) ;
        Igate_dVgs = T7 * T9 * T10_dVg ;
        Igate_dT = T7 * T9 * T10_dT ;
      }
    }

    /* Igs */
      T0 = - pParam->HSMHV_glksd2 * Vgs + modelMKS->HSMHV_glksd3 ;
      T2 = exp (Tox0 * T0);
      T2_dVg = (- Tox0 * pParam->HSMHV_glksd2) * T2;
        
      T0 = Vgs / Tox0 / Tox0 ;
      T3 = Vgs * T0 ;
      T3_dVg = 2.0 * T0 * (1.0 ) ;
      T4 = pParam->HSMHV_glksd1 / 1.0e6 * here->HSMHV_weff_nf ;
      Igs = T4 * T2 * T3 ;
      Igs_dVgs = T4 * (T2_dVg * T3 + T2 * T3_dVg) ;
      Igs_dVds = 0.0 ;
      Igs_dVbs = 0.0 ;
      Igs_dT =   0.0 ;
        
      if ( Vgs >= 0.0e0 ){
        Igs *= -1.0 ;
        Igs_dVgs *= -1.0 ;
        Igs_dVds *= -1.0 ; 
        Igs_dVbs *= -1.0 ;
      }


    /* Igd */
      T1 = Vgs - Vds ;
      T0 = - pParam->HSMHV_glksd2 * T1 + modelMKS->HSMHV_glksd3 ;
      T2 = exp (Tox0 * T0);
      T2_dVg = (- Tox0 * pParam->HSMHV_glksd2) * T2;
      T2_dVd = (+ Tox0 * pParam->HSMHV_glksd2) * T2;
      T2_dVb = 0.0 ;
        
      T0 = T1 / Tox0 / Tox0 ;
      T3 = T1 * T0 ;
      T3_dVg = 2.0 * T0 ;
      T3_dVd = - 2.0 * T0 ;
      T3_dVb = 0.0 ;
      T4 = pParam->HSMHV_glksd1 / 1.0e6 * here->HSMHV_weff_nf ;
      Igd = T4 * T2 * T3 ;
      Igd_dVgs = T4 * (T2_dVg * T3 + T2 * T3_dVg) ;
      Igd_dVds = T4 * (T2_dVd * T3 + T2 * T3_dVd) ;
      Igd_dVbs = 0.0 ;
      Igd_dT =   0.0 ;

      if( T1 >= 0.0e0 ){
        Igd  *= -1.0 ;
        Igd_dVgs *= -1.0 ;
        Igd_dVds *= -1.0 ;
        Igd_dVbs *= -1.0 ;
      }


    /* Igb */
      Etun = ( - ( Vgs - Vbs ) + Vfb + model->HSMHV_glkb3 ) / Tox0 ;
      Etun_dVgs = - 1.0 / Tox0 ;
      Etun_dVds = 0.0 ;
      Etun_dVbs =   1.0 / Tox0 ;

      Fn_SZ( Etun , Etun, igate_dlt, T5) ;
      Etun += small ;
      Etun_dVgs *= T5 ;
      Etun_dVbs *= T5 ;

      T1 = - pParam->HSMHV_glkb2 / Etun ;
      if ( T1 < - EXP_THR ) {
        Igb = 0.0 ;
        Igb_dVgs = Igb_dVds = Igb_dVbs = Igb_dT = 0.0 ;
      } else {
        T2 = exp ( T1 );
        T3 =  pParam->HSMHV_glkb2 / ( Etun * Etun ) * T2 ;
        T2_dVg = T3 * Etun_dVgs ;
	T2_dVb = T3 * Etun_dVbs ;
          
        T3 = pParam->HSMHV_glkb1 * here->HSMHV_weff_nf * Leff ;
        Igb = T3 * Etun * Etun * T2 ;
        Igb_dVgs = T3 * (2.0 * Etun * Etun_dVgs * T2 + Etun * Etun * T2_dVg);
        Igb_dVds = 0.0 ;
	Igb_dVbs = T3 * (2.0 * Etun * Etun_dVbs * T2 + Etun * Etun * T2_dVb);
	Igb_dT = 0.0;
      }

      /* Ifn: Fowler-Nordheim tunneling current */
      Eg12 = here->HSMHV_sqrt_eg ;
      Eg32 = here->HSMHV_eg * Eg12 ;
      T2 = - ( pParam->HSMHV_fvbs * Vbsz - Vgsz + dVthSC + dVthLP - pParam->HSMHV_fn3 ) / Tox0 ;
      T2_dVd = - ( pParam->HSMHV_fvbs * Vbsz_dVds - Vgsz_dVds + dVthSC_dVd + dVthLP_dVd 
                   ) / Tox0 ;
      T2_dVg = - ( - Vgsz_dVgs + dVthSC_dVg + dVthLP_dVg 
                   ) / Tox0 ;
      T2_dVb = - ( pParam->HSMHV_fvbs * Vbsz_dVbs -Vgsz_dVbs + dVthSC_dVb + dVthLP_dVb 
                   ) / Tox0 ;
      T2_dT = - (  pParam->HSMHV_fvbs * Vbsz_dT   -Vgsz_dT   + dVthSC_dT + dVthLP_dT
                   ) / Tox0 ;

      T0 = T2 * T2 ;
      T1 = pParam->HSMHV_fn2 * Eg32 ;
      T1_dT = 1.5 * Eg_dT * pParam->HSMHV_fn2 * Eg12 ;
      T3 = - T1 / T2  ;
      if ( T3 < - EXP_THR ) {
        T5 = 0.0 ;
        T5_dVd = T5_dVg = T5_dVb = T5_dT = 0.0 ;

      } else {
        T5 = exp( T3 ) ;
        T5_dVd = T5 * T1 * T2_dVd / T0 ;
        T5_dVg = T5 * T1 * T2_dVg / T0 ;
        T5_dVb = T5 * T1 * T2_dVb / T0 ;
        T5_dT = T5 * T1 * T2_dT / T0 ;
      }

      T4 = C_QE * pParam->HSMHV_fn1 * here->HSMHV_weff_nf * here->HSMHV_lgate / Eg12 ;  
      T4_dT = (- 0.5) * Eg_dT * T4 / here->HSMHV_eg ;  
      if ( 2e0 * T2 + T1 < 0e0 ){
        Ifn = 0.25e0 * T4 * T1 * T1 * c_exp_2 ; /* minimum value */
        Ifn_dVd = 0e0 ;
        Ifn_dVg = 0e0 ;
        Ifn_dVb = 0e0 ;
        Ifn_dT = 0.25e0 * T4_dT * T1 * T1 * c_exp_2 ;
      } else {
        Ifn = T4 * T0 * T5 ;
        Ifn_dVd = T4 * ( 2.0 * T2 * T2_dVd * T5 + T0 * T5_dVd ) ; 
        Ifn_dVg = T4 * ( 2.0 * T2 * T2_dVg * T5 + T0 * T5_dVg ) ; 
        Ifn_dVb = T4 * ( 2.0 * T2 * T2_dVb * T5 + T0 * T5_dVb ) ; 
        Ifn_dT = T4 * ( 2.0 * T2 * T2_dT * T5 + T0 * T5_dT ) +T4_dT * T0 * T5; 
      }
      Igb -= Ifn ;
      Igb_dVbs -= Ifn_dVb ;
      Igb_dVds -= Ifn_dVd ;
      Igb_dVgs -= Ifn_dVg ;
      Igb_dT -= Ifn_dT ;
  } /* if ( model->HSMHV_coiigs == 0 ) */

  
  /*-----------------------------------------------------------*
   * Vdsp : Vds modification for GIDL/GISL
   *-----------------*/
  if ( model->HSMHV_cogidl != 0 ) {
    T1 = Vds * (1.0 - gidla * Vds) - gidlvds_dlt ;
    T2 = sqrt (T1 * T1 + 4.0 * gidlvds_dlt * Vds) ;
    Vdsp = Vds - 0.5 * (T1 + T2) ;
    T3 = 1.0 - 2.0 * gidla * Vds ;
    Vdsp_dVd = 1.0 - 0.5 * (T3 + (T1 * T3 + 2.0 * gidlvds_dlt) / T2) ;
  }

  /*-----------------------------------------------------------*
   * Igidl : GIDL 
   *-----------------*/
  if( model->HSMHV_cogidl == 0 ){
    Igidl = 0.0e0 ;
    Igidl_dVbs = 0.0e0 ;
    Igidl_dVds = 0.0e0 ;
    Igidl_dVgs = 0.0e0 ;
    Igidl_dT = 0.0e0 ;
  } else {
    T1 = model->HSMHV_gidl3 * (Vdsp + model->HSMHV_gidl4) - Vgs + (dVthSC + dVthLP) * model->HSMHV_gidl5 ;
    T1_dT = (dVthSC_dT + dVthLP_dT) * model->HSMHV_gidl5 ;
    T2 = 1.0 / Tox0 ;
    E1 = T1 * T2 ;
    E1_dVb = ((model->HSMHV_gidl5 * (dVthSC_dVb + dVthLP_dVb)) ) * T2 ;
    E1_dVd = ((model->HSMHV_gidl3 * Vdsp_dVd) + model->HSMHV_gidl5 * (dVthSC_dVd + dVthLP_dVd)) * T2 ;
    E1_dVg = (-1.0 + model->HSMHV_gidl5 * (dVthSC_dVg + dVthLP_dVg) ) * T2 ;
    E1_dT = T1_dT * T2 ;

    Fn_SZ( Egidl , E1, eef_dlt, T5) ;
    Egidl_dVb = T5 * E1_dVb ;
    Egidl_dVd = T5 * E1_dVd ;
    Egidl_dVg = T5 * E1_dVg ;
    Egidl_dT = T5 * E1_dT ;

    T3 = 1 / (Egidl + small) ;
    T0 = - pParam->HSMHV_gidl2 * Egp32 * T3 ;
    T0_dT = - pParam->HSMHV_gidl2 * T3 *( Egp32_dT - Egidl_dT * T3 * Egp32 )  ;
    if ( T0 < - EXP_THR ) {
      Igidl = 0.0 ;
      Igidl_dVbs = Igidl_dVds = Igidl_dVgs = Igidl_dT = 0.0 ;
    } else {
      T1 = exp ( T0 ) ;
      T1_dT = T0_dT * T1 ;
      T2 = pParam->HSMHV_gidl1 / Egp12 * C_QE * here->HSMHV_weff_nf ;
      T2_dT = - Egp12_dT * pParam->HSMHV_gidl1 / Egp12 / Egp12 * C_QE * here->HSMHV_weff_nf ;
      Igidl = T2 * Egidl * Egidl * T1 ;
      T3 = T2 * T1 * Egidl * (2.0 + pParam->HSMHV_gidl2 * Egp32 * Egidl / (Egidl + small) / (Egidl + small)) ;
      Igidl_dVbs = T3 * Egidl_dVb ;
      Igidl_dVds = T3 * Egidl_dVd ;
      Igidl_dVgs = T3 * Egidl_dVg ;
      Igidl_dT = T2 * T1 * Egidl * 2.0 * Egidl_dT  +  T2 * Egidl * Egidl * T1_dT + T2_dT * Egidl * Egidl * T1;
    }
    
    /* bug-fix */
    Vdb = Vds - Vbs ;
    if ( Vdb > 0.0 ) {
      T2 = Vdb * Vdb ;
      T4 = T2 * Vdb ;
      T0 = T4 + C_gidl_delta ;
      T5 = T4 / T0 ;
      T7 = ( 3.0 * T2 * T0 - T4 * 3.0 * T2 ) / ( T0 * T0 ) ; /* == T5_dVdb */
      Igidl_dVbs = Igidl_dVbs * T5 + Igidl * T7 * ( - 1.0 ) ; /* Vdb_dVbs = -1 */
      Igidl_dVds = Igidl_dVds * T5 + Igidl * T7 * ( + 1.0 ) ; /* Vdb_dVds = +1 */
      Igidl_dVgs = Igidl_dVgs * T5 ; /* Vdb_dVgs = 0 */
      Igidl_dT =   Igidl_dT   * T5 ; /* Vdb_dT   = 0 */
      Igidl *= T5 ;
    } else {
      Igidl = 0.0 ;
      Igidl_dVbs = Igidl_dVds = Igidl_dVgs = Igidl_dT = 0.0 ;
    }
  }


  /*-----------------------------------------------------------*
   * Igisl : GISL
   *-----------------*/
  if( model->HSMHV_cogidl == 0){
    Igisl = 0.0e0 ;
    Igisl_dVbs = 0.0e0 ;
    Igisl_dVds = 0.0e0 ;
    Igisl_dVgs = 0.0e0 ;
    Igisl_dT = 0.0e0 ;
  } else {
    T1 = model->HSMHV_gidl3 * ( - Vdsp + model->HSMHV_gidl4 )
      - ( Vgs - Vdsp ) + ( dVthSC + dVthLP ) * model->HSMHV_gidl5 ;

    T1_dT = ( dVthSC_dT + dVthLP_dT ) * model->HSMHV_gidl5 ;
    T2 = 1.0 / Tox0 ;
    E1 = T1 * T2 ;
    E1_dVb = ((model->HSMHV_gidl5 * (dVthSC_dVb + dVthLP_dVb)) ) * T2 ;
    E1_dVd = (((1.0-model->HSMHV_gidl3 ) * Vdsp_dVd) + model->HSMHV_gidl5 * (dVthSC_dVd + dVthLP_dVd)) * T2 ;
    E1_dVg = (-1.0 + model->HSMHV_gidl5 * (dVthSC_dVg + dVthLP_dVg) ) * T2 ;
    E1_dT = T1_dT * T2 ;

    Fn_SZ( Egisl , E1, eef_dlt, T5) ;
    Egisl_dVb = T5 * E1_dVb ;
    Egisl_dVd = T5 * E1_dVd ;
    Egisl_dVg = T5 * E1_dVg ;
    Egisl_dT = T5 * E1_dT ;

    T3 =  1  / (Egisl + small) ;
    T0 = - pParam->HSMHV_gidl2 * Egp32 * T3 ;
    T0_dT = - pParam->HSMHV_gidl2 * T3 * ( Egp32_dT - Egisl_dT * T3 * Egp32 )  ;
    if ( T0 < - EXP_THR ) {
      Igisl = 0.0 ;
      Igisl_dVbs = Igisl_dVds = Igisl_dVgs = Igisl_dT = 0.0 ;
    } else {
      T1 = exp ( T0 ) ;
      T1_dT = T0_dT * T1 ;
      T3 = 1 / Egp12 ;
      T2 = pParam->HSMHV_gidl1 * T3 * C_QE * here->HSMHV_weff_nf ;
      T2_dT = - pParam->HSMHV_gidl1 * Egp12_dT * T3 * T3 * C_QE * here->HSMHV_weff_nf ;
      Igisl = T2 * Egisl * Egisl * T1 ;
      T3 = T2 * T1 * Egisl * (2.0 + pParam->HSMHV_gidl2 * Egp32 * Egisl / (Egisl + small) / (Egisl + small)) ;
      Igisl_dVbs = T3 * Egisl_dVb ;
      Igisl_dVds = T3 * Egisl_dVd ;
      Igisl_dVgs = T3 * Egisl_dVg ;
      Igisl_dT = T2 * T1 * Egisl * 2.0 * Egisl_dT + T2_dT * Egisl * Egisl * T1 + T2 * Egisl * Egisl * T1_dT ;
    }

    /* bug-fix */
    Vsb = - Vbs ;
    if ( Vsb > 0.0 ) {
      T2 = Vsb * Vsb ;
      T4 = T2 * Vsb ;
      T0 = T4 + C_gidl_delta ;
      T5 = T4 / T0 ;
      T7 = ( 3.0 * T2 * T0 - T4 * 3.0 * T2 ) / ( T0 * T0 ) ; /* == T5_dVsb */
      Igisl_dVbs = Igisl_dVbs * T5 + Igisl * T7 * ( - 1.0 ) ; /* Vsb_dVbs = -1 */
      Igisl_dVds = Igisl_dVds * T5 ; /* Vsb_dVds = 0 */
      Igisl_dVgs = Igisl_dVgs * T5 ; /* Vsb_dVgs = 0 */
      Igisl_dT =   Igisl_dT   * T5 ; /* Vsb_dT   = 0 */
      Igisl *= T5 ;
    } else {
      Igisl = 0.0 ;
      Igisl_dVbs = Igisl_dVds = Igisl_dVgs = Igisl_dT = 0.0 ;
    }
  }


  /*-----------------------------------------------------------*
   * End of PART-2. (label) 
   *-----------------*/
/* end_of_part_2: */


  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-3: Overlap charge
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
  Aclm = pParam->HSMHV_clm1 ;
  if ( flg_noqi != 0 ) {
    /*-------------------------------------------*
     * Calculation of Psdl for cases of flg_noqi==1.
     *-----------------*/
    Psdl = Aclm * (Vds + Ps0) + (1.0e0 - Aclm) * Psl ;
    Psdl_dVbs = Aclm * Ps0_dVbs + (1.0e0 - Aclm) * Psl_dVbs ;
    Psdl_dVds = Aclm * (1.0e0 + Ps0_dVds) + (1.0e0 - Aclm) * Psl_dVds ;
    Psdl_dVgs = Aclm * Ps0_dVgs + (1.0e0 - Aclm) * Psl_dVgs ;
    if ( Psdl > Ps0 + Vds - epsm10 ) {
      Psdl = Ps0 + Vds - epsm10 ;
      Psdl_dVbs = Ps0_dVbs ;
      Psdl_dVds = Ps0_dVds + 1.0 ;
      Psdl_dVgs = Ps0_dVgs ;
    }
    
    if (model->HSMHV_xqy !=0) {
      Ec = 0.0e0 ;
      Ec_dVbs =0.0e0 ;
      Ec_dVds =0.0e0 ;
      Ec_dVgs =0.0e0 ;
      Ec_dT =0.0e0 ; 
    }
  } else {
    /* Ec is removed from Lred calc. part */
    if (model->HSMHV_xqy !=0) {
      if ( Idd < C_IDD_MIN ) {
        Ec = 0.0e0 ;
        Ec_dVbs =0.0e0 ;
        Ec_dVds =0.0e0 ;
        Ec_dVgs =0.0e0 ;
      } else {
      T1 =  beta_inv / Leff ;
      T1_dT = beta_inv_dT / Leff ; 
      T2 = 1.0 / Qn0 ;
      T3 = T2 * T2 ;
      Ec = Idd * T1 * T2 ;
      Ec_dVbs = T1 * (Idd_dVbs * T2 - Idd * Qn0_dVbs * T3 ) ;
      Ec_dVds = T1 * (Idd_dVds * T2 - Idd * Qn0_dVds * T3 ) ;
      Ec_dVgs = T1 * (Idd_dVgs * T2 - Idd * Qn0_dVgs * T3 ) ;
      Ec_dT = T1 * (Idd_dT * T2 - Idd * Qn0_dT * T3 ) 
            + T1_dT * Idd * T2 ; 
      }
    }
  }

  /*-----------------------------------------------------------*
   * Overlap charges
   *-----------------*/
  Coovlps = (int)ModeNML * model->HSMHV_coovlps + (int)ModeRVS * model->HSMHV_coovlp ;
  Coovlpd = (int)ModeRVS * model->HSMHV_coovlps + (int)ModeNML * model->HSMHV_coovlp ;
  Lovers = ModeNML * here->HSMHV_lovers + ModeRVS * here->HSMHV_loverld ;
  Loverd = ModeRVS * here->HSMHV_lovers + ModeNML * here->HSMHV_loverld ;
  Novers = ModeNML * pParam->HSMHV_novers + ModeRVS * pParam->HSMHV_nover ;
  Noverd = ModeRVS * pParam->HSMHV_novers + ModeNML * pParam->HSMHV_nover ;
  CVDSOVER = pParam->HSMHV_cvdsover ;

  /*---------------------------------------------------*
   * Source side (mode-dependent)
   *-----------------*/
  /*-----------------------------------*
  * Constant capacitance model
  *-----------------*/
  if ( Coovlps == 0 ) {
    flg_overgiven = ( (int)ModeNML * model->HSMHV_cgso_Given 
                    + (int)ModeRVS * model->HSMHV_cgdo_Given  ) ;
    if ( flg_overgiven ) {
      Cgso = ModeNML * pParam->HSMHV_cgso + ModeRVS * pParam->HSMHV_cgdo ;
      Cgso *= - here->HSMHV_weffcv_nf ;
    } else {
      Cgso = - Cox0 * Lovers * here->HSMHV_weffcv_nf ;
    }

    Qgso = - Cgso * Vgse ;
    Qgso_dVbse = 0.0 ;
    Qgso_dVdse = 0.0 ;
    Qgso_dVgse = - Cgso ;

  /*-----------------------------------*
   * Simplified model
   *-----------------*/
  } else { /* Coovlps != 0 begin */
    if ( Lovers > 0.0 && Novers == 0.0 ){
      cov_slp = modelMKS->HSMHV_ovslp ;
      cov_mag = model->HSMHV_ovmag ;
      T1 = Cox0 * here->HSMHV_weffcv_nf ;
      T4 = cov_slp * T1 * ( cov_mag + Vgs ) ;
      T4_dVg = cov_slp * T1 ;
      T4_dVd = 0.0 ;
      T5 = Lovers * T1 ;
      T9 = 1.2e0 - Ps0 ;
      Qgos = Vgs * T5 - T4 * T9 ;
      Qgos_dVbs =      T4 * Ps0_dVbs ;
      Qgos_dVds =      T4 * Ps0_dVds - T9 * T4_dVd ;
      Qgos_dVgs = T5 + T4 * Ps0_dVgs - T9 * T4_dVg ;
      Qgos_dT  = T4 * Ps0_dT;

  /*-----------------------------------*
   * Surface potential model
   *------------------------*/
  } else if ( Lovers > 0.0 && Novers >= 0.0 ) {

    Vgbgmt = Vgs - Vbs ;
    Vgbgmt_dVbs = -1.0 ;
    Vgbgmt_dVds = 0.0 ;
    Vgbgmt_dVgs = 1.0 ;
    Vxbgmt = - Vbs ;
    Vxbgmt_dVbs = -1.0 ;
    Vxbgmt_dVds = 0.0 ;
    Vxbgmt_dVgs = 0.0 ;

    Nover_func = Novers ;
    cnst0over_func =    ModeNML * here->HSMHV_cnst0overs + ModeRVS * here->HSMHV_cnst0over ;
    cnst0over_func_dT = ModeNML * cnst0overs_dT         + ModeRVS * cnst0over_dT ;
/*    ps0ldinib_func =    ModeNML * here->HSMHV_ps0ldinibs + ModeRVS * here->HSMHV_ps0ldinib ;*/
/*    ps0ldinib_func_dT = ModeNML * ps0ldinibs_dT         + ModeRVS * ps0ldinib_dT ;*/
#include "hsmhveval_qover.h"

    T4 = here->HSMHV_weffcv_nf * Lovers * ( 1 - CVDSOVER ) ;

    Qovs =  T4 * QsuLD ;
    Qovs_dVds = T4 * QsuLD_dVds ;
    Qovs_dVgs = T4 * QsuLD_dVgs ;
    Qovs_dVbs = T4 * QsuLD_dVbs ;
    Qovs_dT  = T4 * QsuLD_dT ;
  
/*    QisLD = T4 * QiuLD ;*/
    QisLD_dVbs = T4 * QiuLD_dVbs ;
    QisLD_dVds = T4 * QiuLD_dVds ;
    QisLD_dVgs = T4 * QiuLD_dVgs ;
    QisLD_dT = T4 * QiuLD_dT ; 
  
    QbsLD = T4 * QbuLD ;
/*    QbsLD_dVbs = T4 * QbuLD_dVbs ;
    QbsLD_dVds = T4 * QbuLD_dVds ;
    QbsLD_dVgs = T4 * QbuLD_dVgs ;
    QbsLD_dT = T4 * QbuLD_dT ; 
*/


    if ( CVDSOVER != 0.0 ) { /* Qovsext begin */
      Vgbgmt = Vgse - Vbse ;
      Vgbgmt_dVbs = -1.0 ;
      Vgbgmt_dVds = 0.0 ;
      Vgbgmt_dVgs = 1.0 ;
      Vxbgmt = - Vbse ;
      Vxbgmt_dVbs = -1.0 ;
      Vxbgmt_dVds = 0.0 ;
      Vxbgmt_dVgs = 0.0 ;
      
#include "hsmhveval_qover.h"
      
      T4 = here->HSMHV_weffcv_nf * Lovers * CVDSOVER ;
      Qovsext =  T4 * QsuLD ;
      Qovsext_dVdse = T4 * QsuLD_dVds ;
      Qovsext_dVgse = T4 * QsuLD_dVgs ;
      Qovsext_dVbse = T4 * QsuLD_dVbs ;
      Qovsext_dT   = T4 * QsuLD_dT ;
      
/*      QisLDext = T4 * QiuLD ;*/
      QisLDext_dVbse = T4 * QiuLD_dVbs ;
      QisLDext_dVdse = T4 * QiuLD_dVds ;
      QisLDext_dVgse = T4 * QiuLD_dVgs ;
      QisLDext_dT = T4 * QiuLD_dT ; 
  
      QbsLDext = T4 * QbuLD ;
/*      QbsLDext_dVbse = T4 * QbuLD_dVbs ;
      QbsLDext_dVdse = T4 * QbuLD_dVds ;
      QbsLDext_dVgse = T4 * QbuLD_dVgs ;
      QbsLDext_dT = T4 * QbuLD_dT ; 
*/      
    } /* Qovsext end */

  } 

    /*-----------------------------------*
     * Additional constant capacitance model
     *-----------------*/
    flg_overgiven = ( (int)ModeNML * model->HSMHV_cgso_Given 
	             + (int)ModeRVS * model->HSMHV_cgdo_Given  ) ;
    if ( flg_overgiven ) {
      Cgso  = ModeNML * pParam->HSMHV_cgso + ModeRVS * pParam->HSMHV_cgdo ;
      Cgso *= - here->HSMHV_weffcv_nf ;
    }
    Qgso = - Cgso * Vgse ;
    Qgso_dVbse = 0.0 ;
    Qgso_dVdse = 0.0 ;
    Qgso_dVgse = - Cgso ;
  } /* Coovlps != 0 end */

  /*---------------------------------------------------*
   * Drain side (mode-dependent)
   *-----------------*/
  /*-----------------------------------*
  * Constant capacitance model
  *-----------------*/
  if ( Coovlpd == 0 ) {
    flg_overgiven = ( (int)ModeRVS * model->HSMHV_cgso_Given 
                    + (int)ModeNML * model->HSMHV_cgdo_Given  ) ;
    if ( flg_overgiven ) {
      Cgdo = ModeRVS * pParam->HSMHV_cgso + ModeNML * pParam->HSMHV_cgdo ;
      Cgdo *= - here->HSMHV_weffcv_nf ;
    } else {
      Cgdo = - Cox0 * Loverd * here->HSMHV_weffcv_nf ;
    }

    Qgdo = - Cgdo * (Vgse - Vdse) ;
    Qgdo_dVbse = 0.0 ;
    Qgdo_dVdse = Cgdo ;
    Qgdo_dVgse = - Cgdo ;

  /*-----------------------------------*
   * Simplified model
   *-----------------*/
  } else { /* Coovlpd != 0 begin */
    if ( Loverd > 0.0 && Noverd == 0.0 ){
      cov_slp = modelMKS->HSMHV_ovslp ;
      cov_mag = model->HSMHV_ovmag ;
      T1 = Cox0 * here->HSMHV_weffcv_nf ;
      T4 = cov_slp * T1 * ( cov_mag + Vgs - Vds ) ;
	  T4_dVg = cov_slp * T1 ;
	  T4_dVd = - cov_slp * T1 ;
	  T5 = Loverd * T1 ;
	  T9 = 1.2e0 + Vds - Psl ;
	  Qgod = ( Vgs - Vds ) * T5 - T4 * T9 ;
          Qgod_dVbs =      + T4 * Psl_dVbs ;
          Qgod_dVds = - T5 + T4 * ( -1.0 + Psl_dVds ) - T9 * T4_dVd ;
          Qgod_dVgs = + T5 + T4 * Psl_dVgs - T9 * T4_dVg ;
          Qgod_dT  = T4 * Psl_dT;


  /*-----------------------------------*
   * Surface potential model 
   *------------------------*/
  } else if ( Loverd > 0.0 && Noverd >= 0.0 ) {

    Vgbgmt = Vgs - Vbs ;
    Vgbgmt_dVbs = -1.0 ;
    Vgbgmt_dVds = 0.0 ;
    Vgbgmt_dVgs = 1.0 ;
    Vxbgmt = Vds - Vbs ;
    Vxbgmt_dVbs = -1.0 ;
    Vxbgmt_dVds = 1.0 ;
    Vxbgmt_dVgs = 0.0 ;

    Nover_func = Noverd ;
    cnst0over_func =    ModeNML * here->HSMHV_cnst0over + ModeRVS * here->HSMHV_cnst0overs ;
    cnst0over_func_dT = ModeNML * cnst0over_dT         + ModeRVS * cnst0overs_dT ;
/*    ps0ldinib_func =    ModeNML * here->HSMHV_ps0ldinib + ModeRVS * here->HSMHV_ps0ldinibs ;*/
/*    ps0ldinib_func_dT = ModeNML * ps0ldinib_dT         + ModeRVS * ps0ldinibs_dT ;*/
#include "hsmhveval_qover.h"

    T4 = here->HSMHV_weffcv_nf * Loverd * ( 1 - CVDSOVER ) ;
    Qovd =  T4 * QsuLD ;
    Qovd_dVds = T4 * QsuLD_dVds ;
    Qovd_dVgs = T4 * QsuLD_dVgs ;
    Qovd_dVbs = T4 * QsuLD_dVbs ;
    Qovd_dT = T4 * QsuLD_dT ;
  
/*    QidLD = T4 * QiuLD ;*/
    QidLD_dVbs = T4 * QiuLD_dVbs ;
    QidLD_dVds = T4 * QiuLD_dVds ;
    QidLD_dVgs = T4 * QiuLD_dVgs ;
    QidLD_dT = T4 * QiuLD_dT ; 
  
    QbdLD = T4 * QbuLD ;
    QbdLD_dVbs = T4 * QbuLD_dVbs ;
    QbdLD_dVds = T4 * QbuLD_dVds ;
    QbdLD_dVgs = T4 * QbuLD_dVgs ;
    QbdLD_dT = T4 * QbuLD_dT ; 


    if ( CVDSOVER != 0.0 ) { /* Qovdext begin */
      Vgbgmt = Vgse - Vbse ;
      Vgbgmt_dVbs = -1.0 ;
      Vgbgmt_dVds = 0.0 ;
      Vgbgmt_dVgs = 1.0 ;
      Vxbgmt = Vdse - Vbse ;
      Vxbgmt_dVbs = -1.0 ;
      Vxbgmt_dVds = 1.0 ;
      Vxbgmt_dVgs = 0.0 ;

#include "hsmhveval_qover.h"

      T4 = here->HSMHV_weffcv_nf * Loverd * CVDSOVER ;
      Qovdext =  T4 * QsuLD ;
      Qovdext_dVdse = T4 * QsuLD_dVds ;
      Qovdext_dVgse = T4 * QsuLD_dVgs ;
      Qovdext_dVbse = T4 * QsuLD_dVbs ;
      Qovdext_dT   = T4 * QsuLD_dT ;

/*      QidLDext = T4 * QiuLD ;*/
      QidLDext_dVbse = T4 * QiuLD_dVbs ;
      QidLDext_dVdse = T4 * QiuLD_dVds ;
      QidLDext_dVgse = T4 * QiuLD_dVgs ;
      QidLDext_dT = T4 * QiuLD_dT ; 
  
      QbdLDext = T4 * QbuLD ;
      QbdLDext_dVbse = T4 * QbuLD_dVbs ;
      QbdLDext_dVdse= T4 * QbuLD_dVds ;
      QbdLDext_dVgse= T4 * QbuLD_dVgs ;
      QbdLDext_dT = T4 * QbuLD_dT ; 

    } /* Qovdext end */

  }
    /*-----------------------------------*
     * Additional constant capacitance model
     *-----------------*/
    flg_overgiven = ( (int)ModeRVS * model->HSMHV_cgso_Given
                    + (int)ModeNML * model->HSMHV_cgdo_Given  ) ;
    if ( flg_overgiven ) {
      Cgdo  = ModeRVS * pParam->HSMHV_cgso + ModeNML * pParam->HSMHV_cgdo ;
      Cgdo *= - here->HSMHV_weffcv_nf ;
    }
    Qgdo = - Cgdo * (Vgse - Vdse) ;
    Qgdo_dVbse = 0.0 ;
    Qgdo_dVdse = Cgdo ;
    Qgdo_dVgse = - Cgdo ;
  } /* Coovlpd != 0 end */

  /*-------------------------------------------*
   * Gate/Bulk overlap charge: Qgbo
   *-----------------*/
  Cgbo_loc = - pParam->HSMHV_cgbo * here->HSMHV_lgate ;
  Qgbo = - Cgbo_loc * (Vgs -Vbs) ;
  Qgbo_dVgs = - Cgbo_loc ;
  Qgbo_dVbs =   Cgbo_loc ;
  Qgbo_dVds = 0.0 ;

  /*---------------------------------------------------*
   * Lateral-field-induced capacitance.
   *-----------------*/
  if ( model->HSMHV_xqy == 0 ){
    Qy = 0.0e0 ;
    Qy_dVds = 0.0e0 ;
    Qy_dVgs = 0.0e0 ;
    Qy_dVbs = 0.0e0 ;
    Qy_dT = 0.0e0 ;
  } else {
    Pslk = Ec * Leff + Ps0 ;
    Pslk_dVbs = Ec_dVbs * Leff + Ps0_dVbs;
    Pslk_dVds = Ec_dVds * Leff + Ps0_dVds;
    Pslk_dVgs = Ec_dVgs * Leff + Ps0_dVgs;
    Pslk_dT = Ec_dT * Leff + Ps0_dT; 
      
    T1 = Aclm * ( Vds + Ps0 ) + ( 1.0e0 - Aclm ) * Pslk ;
    T1_dVb = Aclm * (  Ps0_dVbs ) + ( 1.0e0 - Aclm ) * Pslk_dVbs ;
    T1_dVd = Aclm * ( 1.0 + Ps0_dVds ) + ( 1.0e0 - Aclm ) * Pslk_dVds ;
    T1_dVg = Aclm * (  Ps0_dVgs ) + ( 1.0e0 - Aclm ) * Pslk_dVgs ;
    T1_dT = Aclm * ( Ps0_dT ) + ( 1.0e0 - Aclm ) * Pslk_dT ; 
    T10 = here->HSMHV_wdpl ;
    T3 = T10 * 1.3 ;
    T2 = C_ESI * here->HSMHV_weffcv_nf * T3 ;
    Qy =  - ( ( Ps0 + Vds  - T1 ) / model->HSMHV_xqy ) * T2 ;
    Qy_dVds =  - ( ( Ps0_dVds + 1.0e0 - T1_dVd ) / model->HSMHV_xqy ) * T2 ;
    Qy_dVgs =  - ( ( Ps0_dVgs  - T1_dVg ) / model->HSMHV_xqy ) * T2 ;
    Qy_dVbs =  - ( ( Ps0_dVbs  - T1_dVb ) / model->HSMHV_xqy ) * T2 ;
    Qy_dT =  - ( ( Ps0_dT - T1_dT ) / model->HSMHV_xqy ) * T2 ; 
  }

  if ( model->HSMHV_xqy1 != 0.0 ){
    Qy += here->HSMHV_cqyb0 * Vbs ;
    Qy_dVbs += here->HSMHV_cqyb0 ;
  }

  /*---------------------------------------------------* 
   * Fringing capacitance.
   *-----------------*/ 
  Cfd = here->HSMHV_cfrng ;
  Cfs = here->HSMHV_cfrng ;
  Qfd = Cfd * ( Vgse - Vdse ) ;
  Qfs = Cfs * Vgse ;

  /*-----------------------------------------------------------* 
   * End of PART-3. (label) 
   *-----------------*/ 

/* end_of_part_3: */

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-4: Substrate-source/drain junction diode.
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
    
   
  /*-----------------------------------------------------------*
   * Cbsj, Cbdj: node-base S/D biases.
   *-----------------*/

  T10 = model->HSMHV_cvb * here->HSMHV_jd_nvtm_inv ;
  T11 = model->HSMHV_cvbk * here->HSMHV_jd_nvtm_inv ;

  T9 = model->HSMHV_cisb * here->HSMHV_exptemp ;
  T0 = here->HSMHV_isbd2 * T9 ;

  T2 = exp (- vbd_jct * T10 );
  T2_dVb = - T2 * T10 ;

  T3 = exp (- vbd_jct * T11 );
  T3_dVb = - T3 * T11 ;

  /* self heating */
  T10_dT = model->HSMHV_cvb * beta_dT / pParam->HSMHV_nj ;
  T11_dT = model->HSMHV_cvbk * beta_dT / pParam->HSMHV_nj ;
  T9_dT = model->HSMHV_cisb * exptemp_dT ;
  T0_dT = here->HSMHV_isbd2 * T9_dT + isbd2_dT * T9 ;
  T2_dT = -vbd_jct * T10 * T2 * beta_dT * beta_inv ;
  T3_dT = -vbd_jct * T11 * T3 * beta_dT * beta_inv ;


  /* ibd */
  if ( vbd_jct < here->HSMHV_vbdt ) {
    TX = vbd_jct * here->HSMHV_jd_nvtm_inv ;


    if ( TX < - EXP_THR ) {
      T1 = 0.0 ;
      T1_dVb = 0.0 ;
      T1_dT =  0.0 ;
    } else {
      T1 = exp ( TX ) ;
      T1_dVb = T1 * here->HSMHV_jd_nvtm_inv ;
      T1_dT = T1 * TX * beta_dT * beta_inv ;

    }

    Ibd = here->HSMHV_isbd * (T1 - 1.0) 
      + T0 * (T2 - 1.0) 
      + pParam->HSMHV_cisbk * (T3 - 1.0);   
    Gbd = here->HSMHV_isbd * T1_dVb 
      + T0 * T2_dVb 
      + pParam->HSMHV_cisbk * T3_dVb ;
    Ibd_dT = here->HSMHV_isbd * T1_dT + isbd_dT * ( T1 - 1.0 )
      + T0 * T2_dT + T0_dT * ( T2 - 1.0 )
      + pParam->HSMHV_cisbk * T3_dT ;

  } else {
    T1 = here->HSMHV_jd_expcd ;

    T4 = here->HSMHV_isbd * here->HSMHV_jd_nvtm_inv  * T1 ;

    Ibd = here->HSMHV_isbd * (T1 - 1.0) 
      + T4 * (vbd_jct - here->HSMHV_vbdt) 
      + T0 * (T2 - 1.0)
      + pParam->HSMHV_cisbk * (T3 - 1.0) ;
    Gbd = T4 
      + T0 * T2_dVb 
      + pParam->HSMHV_cisbk * T3_dVb ;

    T1_dT = jd_expcd_dT ;
    T4_dT = isbd_dT * here->HSMHV_jd_nvtm_inv * T1
      + here->HSMHV_isbd * jd_nvtm_inv_dT * T1
      + here->HSMHV_isbd * here->HSMHV_jd_nvtm_inv * T1_dT ;
    Ibd_dT = isbd_dT * ( T1 - 1.0 ) + here->HSMHV_isbd * T1_dT
      + T4_dT * ( vbd_jct - here->HSMHV_vbdt ) - T4 * vbdt_dT
      + T0_dT * ( T2 - 1.0 ) + T0 * T2_dT
      + pParam->HSMHV_cisbk * T3_dT ;
  }  
  T12 = model->HSMHV_divx * here->HSMHV_isbd2 ;
  Ibd += T12 * vbd_jct ;
  Gbd += T12 ;

  T12_dT = model->HSMHV_divx * isbd2_dT ;
  Ibd_dT += T12_dT * vbd_jct ;

  /* ibs */
  T0 = here->HSMHV_isbs2 * T9 ;
  T0_dT = here->HSMHV_isbs2 * T9_dT + isbs2_dT * T9 ;

  TX = - vbs_jct * T10 ;
  if ( TX < - EXP_THR ) {
    T2 = 0.0 ;
    T2_dVb = 0.0 ;
    T2_dT =  0.0 ;
  } else {
    T2 = exp ( TX );
    T2_dVb = - T2 * T10 ;
    T2_dT = T2 * TX * beta_dT * beta_inv ;
  }

  TX = - vbs_jct * T11 ;
  if ( TX < - EXP_THR ) {
    T3 = 0.0 ;
    T3_dVb = 0.0 ;
    T3_dT =  0.0 ;
  } else {
    T3 = exp ( TX );
    T3_dVb = - T3 * T11 ;
    T3_dT = T3 * TX * beta_dT * beta_inv ;
  }

  if ( vbs_jct < here->HSMHV_vbst ) {
    TX = vbs_jct * here->HSMHV_jd_nvtm_inv ;
    if ( TX < - EXP_THR ) {
      T1 = 0.0 ;
      T1_dVb = 0.0 ;
      T1_dT =  0.0 ;
    } else {
      T1 = exp ( TX ) ;
      T1_dVb = T1 * here->HSMHV_jd_nvtm_inv ;
      T1_dT = T1 * TX * beta_dT * beta_inv ;
    }
    Ibs = here->HSMHV_isbs * (T1 - 1.0) 
      + T0 * (T2 - 1.0) 
      + pParam->HSMHV_cisbk * (T3 - 1.0);
    Gbs = here->HSMHV_isbs * T1_dVb 
      + T0 * T2_dVb
      + pParam->HSMHV_cisbk * T3_dVb ;
    Ibs_dT = here->HSMHV_isbs * T1_dT + isbs_dT * ( T1 - 1.0 )
      + T0 * T2_dT + T0_dT * ( T2 - 1.0 )
      + pParam->HSMHV_cisbk * T3_dT ;
  } else {
    T1 = here->HSMHV_jd_expcs ;

    T4 = here->HSMHV_isbs * here->HSMHV_jd_nvtm_inv  * T1 ;

    Ibs = here->HSMHV_isbs * (T1 - 1.0)
      + T4 * (vbs_jct - here->HSMHV_vbst)
      + T0 * (T2 - 1.0) 
      + pParam->HSMHV_cisbk * (T3 - 1.0) ;
    Gbs = T4 
      + T0 * T2_dVb 
      + pParam->HSMHV_cisbk * T3_dVb ;

    T1_dT = jd_expcs_dT ;
    T4_dT = isbs_dT * here->HSMHV_jd_nvtm_inv * T1
      + here->HSMHV_isbs * jd_nvtm_inv_dT * T1
      + here->HSMHV_isbs * here->HSMHV_jd_nvtm_inv * T1_dT ;
    Ibs_dT = isbs_dT * ( T1 - 1.0 ) + here->HSMHV_isbs * T1_dT
      + T4_dT * ( vbs_jct - here->HSMHV_vbst) - T4 * vbst_dT
      + T0_dT * ( T2 - 1.0 ) + T0 * T2_dT
      + pParam->HSMHV_cisbk * T3_dT ;
  }
  T12 = model->HSMHV_divx * here->HSMHV_isbs2 ;
  Ibs += T12 * vbs_jct ;
  Gbs += T12 ;

  T12_dT = model->HSMHV_divx * isbs2_dT ;
  Ibs_dT += T12_dT * vbs_jct ;


  /*-----------------------------------------------------------*
   * Charges and Capacitances.
   *-----------------*/
  /*  charge storage elements
   *  bulk-drain and bulk-source depletion capacitances
   *  czbd : zero bias drain junction capacitance
   *  czbs : zero bias source junction capacitance
   *  czbdsw:zero bias drain junction sidewall capacitance
   *  czbssw:zero bias source junction sidewall capacitance
   */
  tcjbd = model->HSMHV_tcjbd ;
  tcjbs = model->HSMHV_tcjbs ;
  tcjbdsw = model->HSMHV_tcjbdsw ;
  tcjbssw = model->HSMHV_tcjbssw ;
  tcjbdswg = model->HSMHV_tcjbdswg ;
  tcjbsswg = model->HSMHV_tcjbsswg ;

  czbs = model->HSMHV_cj * here->HSMHV_as ;
  czbs = czbs * ( 1.0 + tcjbs * ( TTEMP - model->HSMHV_ktnom )) ;
  czbs_dT = ( model->HSMHV_cj * here->HSMHV_as ) * tcjbs ;

  czbd = model->HSMHV_cj * here->HSMHV_ad ;
  czbd = czbd * ( 1.0 + tcjbd * ( TTEMP - model->HSMHV_ktnom )) ;
  czbd_dT = ( model->HSMHV_cj * here->HSMHV_ad ) * tcjbd ;

  /* Source Bulk Junction */
  if (here->HSMHV_ps > here->HSMHV_weff_nf) {
    czbssw = model->HSMHV_cjsw * ( here->HSMHV_ps - here->HSMHV_weff_nf ) ;
    czbssw = czbssw * ( 1.0 + tcjbssw * ( TTEMP - model->HSMHV_ktnom )) ;
    czbssw_dT = ( model->HSMHV_cjsw * ( here->HSMHV_ps - here->HSMHV_weff_nf )) * tcjbssw ;

    czbsswg = model->HSMHV_cjswg * here->HSMHV_weff_nf ;
    czbsswg = czbsswg * ( 1.0 + tcjbsswg * ( TTEMP - model->HSMHV_ktnom )) ;
    czbsswg_dT = ( model->HSMHV_cjswg * here->HSMHV_weff_nf ) * tcjbsswg ;

//  if (vbs_jct == 0.0) {  
    if (0) {  
      Qbs = 0.0 ;
      Qbs_dT = 0.0 ;
      Capbs = czbs + czbssw + czbsswg ;
    } else if (vbs_jct < 0.0) { 
      if (czbs > 0.0) {
        arg = 1.0 - vbs_jct / model->HSMHV_pb ;
        if (model->HSMHV_mj == 0.5) 
          sarg = 1.0 / sqrt(arg) ;
        else 
          sarg = Fn_Pow( arg , -model->HSMHV_mj ) ;
        Qbs = model->HSMHV_pb * czbs * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Qbs_dT = model->HSMHV_pb * czbs_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Capbs = czbs * sarg ;
      } else {
        Qbs = 0.0 ;
        Qbs_dT = 0.0 ;
        Capbs = 0.0 ;
      }
      if (czbssw > 0.0) {
        arg = 1.0 - vbs_jct / model->HSMHV_pbsw ;
        if (model->HSMHV_mjsw == 0.5) 
          sarg = 1.0 / sqrt(arg) ;
        else 
          sarg = Fn_Pow( arg , -model->HSMHV_mjsw ) ;
        Qbs += model->HSMHV_pbsw * czbssw * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjsw) ;
        Qbs_dT += model->HSMHV_pbsw * czbssw_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjsw) ;
        Capbs += czbssw * sarg ;
      }
      if (czbsswg > 0.0) {
        arg = 1.0 - vbs_jct / model->HSMHV_pbswg ;
        if (model->HSMHV_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mjswg ) ;
        Qbs += model->HSMHV_pbswg * czbsswg * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Qbs_dT += model->HSMHV_pbswg * czbsswg_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Capbs += czbsswg * sarg ;
      }
    } else {
      T1 = czbs + czbssw + czbsswg ;
      T1_dT = czbs_dT + czbssw_dT + czbsswg_dT ;
      T2 = czbs * model->HSMHV_mj / model->HSMHV_pb 
        + czbssw * model->HSMHV_mjsw / model->HSMHV_pbsw 
        + czbsswg * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      T2_dT = czbs_dT * model->HSMHV_mj / model->HSMHV_pb 
        + czbssw_dT * model->HSMHV_mjsw / model->HSMHV_pbsw 
        + czbsswg_dT * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      Qbs = vbs_jct * (T1 + vbs_jct * 0.5 * T2) ;
      Qbs_dT = vbs_jct * (T1_dT + vbs_jct * 0.5 * T2_dT) ;
      Capbs = T1 + vbs_jct * T2 ;
    }
  } else {
    czbsswg = model->HSMHV_cjswg * here->HSMHV_ps ;
    czbsswg = czbsswg * ( 1.0 + tcjbsswg * ( TTEMP - model->HSMHV_ktnom )) ;
    czbsswg_dT = ( model->HSMHV_cjswg * here->HSMHV_ps ) * tcjbsswg ;
//  if (vbs_jct == 0.0) {
    if (0) {
      Qbs = 0.0 ;
      Qbs_dT = 0.0 ;
      Capbs = czbs + czbsswg ;
    } else if (vbs_jct < 0.0) {
      if (czbs > 0.0) {
        arg = 1.0 - vbs_jct / model->HSMHV_pb ;
        if (model->HSMHV_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mj ) ;
        Qbs = model->HSMHV_pb * czbs * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Qbs_dT = model->HSMHV_pb * czbs_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Capbs = czbs * sarg ;
      } else {
        Qbs = 0.0 ;
        Qbs_dT = 0.0 ;
        Capbs = 0.0 ;
      }
      if (czbsswg > 0.0) {
        arg = 1.0 - vbs_jct / model->HSMHV_pbswg ;
        if (model->HSMHV_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mjswg ) ;
        Qbs += model->HSMHV_pbswg * czbsswg * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Qbs_dT += model->HSMHV_pbswg * czbsswg_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Capbs += czbsswg * sarg ;
      }
    } else {
      T1 = czbs + czbsswg ;
      T1_dT = czbs_dT + czbsswg_dT ;
      T2 = czbs * model->HSMHV_mj / model->HSMHV_pb
        + czbsswg * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      T2_dT = czbs_dT * model->HSMHV_mj / model->HSMHV_pb
        + czbsswg_dT * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      Qbs = vbs_jct * (T1 + vbs_jct * 0.5 * T2) ;
      Qbs_dT = vbs_jct * (T1_dT + vbs_jct * 0.5 * T2_dT) ;
      Capbs = T1 + vbs_jct * T2 ;
    }
  }    
    
  /* Drain Bulk Junction */
  if (here->HSMHV_pd > here->HSMHV_weff_nf) {

    czbdsw = model->HSMHV_cjsw * ( here->HSMHV_pd - here->HSMHV_weff_nf ) ;
    czbdsw = czbdsw * ( 1.0 + tcjbdsw * ( TTEMP - model->HSMHV_ktnom )) ;
    czbdsw_dT = ( model->HSMHV_cjsw * ( here->HSMHV_pd - here->HSMHV_weff_nf )) * tcjbdsw ;

    czbdswg = model->HSMHV_cjswg * here->HSMHV_weff_nf ;
    czbdswg = czbdswg * ( 1.0 + tcjbdswg * ( TTEMP - model->HSMHV_ktnom )) ;
    czbdswg_dT = ( model->HSMHV_cjswg * here->HSMHV_weff_nf ) * tcjbdswg ;
//  if (vbd_jct == 0.0) {
    if (0) {
      Qbd = 0.0 ;
      Qbd_dT = 0.0 ;
      Capbd = czbd + czbdsw + czbdswg ;
    } else if (vbd_jct < 0.0) {
      if (czbd > 0.0) {
        arg = 1.0 - vbd_jct / model->HSMHV_pb ;
        if (model->HSMHV_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mj ) ;
        Qbd = model->HSMHV_pb * czbd * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Qbd_dT = model->HSMHV_pb * czbd_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Capbd = czbd * sarg ;
      } else {
        Qbd = 0.0 ;
        Qbd_dT = 0.0 ;
        Capbd = 0.0 ;
      }
      if (czbdsw > 0.0) {
        arg = 1.0 - vbd_jct / model->HSMHV_pbsw ;
        if (model->HSMHV_mjsw == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mjsw ) ;
        Qbd += model->HSMHV_pbsw * czbdsw * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjsw) ;
        Qbd_dT += model->HSMHV_pbsw * czbdsw_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjsw) ;
        Capbd += czbdsw * sarg ;
      }
      if (czbdswg > 0.0) {
        arg = 1.0 - vbd_jct / model->HSMHV_pbswg ;
        if (model->HSMHV_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mjswg ) ;
        Qbd += model->HSMHV_pbswg * czbdswg * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Qbd_dT += model->HSMHV_pbswg * czbdswg_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Capbd += czbdswg * sarg ;
      
      }
    } else {
      T1 = czbd + czbdsw + czbdswg ;
      T1_dT = czbd_dT + czbdsw_dT + czbdswg_dT ;
      T2 = czbd * model->HSMHV_mj / model->HSMHV_pb
        + czbdsw * model->HSMHV_mjsw / model->HSMHV_pbsw 
        + czbdswg * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      T2_dT = czbd_dT * model->HSMHV_mj / model->HSMHV_pb
        + czbdsw_dT * model->HSMHV_mjsw / model->HSMHV_pbsw 
        + czbdswg_dT * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      Qbd = vbd_jct * (T1 + vbd_jct * 0.5 * T2) ;
      Qbd_dT = vbd_jct * (T1_dT + vbd_jct * 0.5 * T2_dT) ;
      Capbd = T1 + vbd_jct * T2 ;
    }
    
  } else {
    czbdswg = model->HSMHV_cjswg * here->HSMHV_pd ;
    czbdswg = czbdswg * ( 1.0 + tcjbdswg * ( TTEMP - model->HSMHV_ktnom )) ;
    czbdswg_dT = ( model->HSMHV_cjswg * here->HSMHV_pd ) * tcjbdswg ;

//  if (vbd_jct == 0.0) {   
    if (0) {   
      Qbd = 0.0 ;
      Qbd_dT = 0.0 ;
      Capbd = czbd + czbdswg ;
    } else if (vbd_jct < 0.0) {
      if (czbd > 0.0) {
        arg = 1.0 - vbd_jct / model->HSMHV_pb ;
        if (model->HSMHV_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
	else
          sarg = Fn_Pow( arg , -model->HSMHV_mj ) ;
        Qbd = model->HSMHV_pb * czbd * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Qbd_dT = model->HSMHV_pb * czbd_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mj) ;
        Capbd = czbd * sarg ;
      } else {
        Qbd = 0.0 ;
        Qbd_dT = 0.0 ;
        Capbd = 0.0 ;
      }
      if (czbdswg > 0.0) {
        arg = 1.0 - vbd_jct / model->HSMHV_pbswg ;
        if (model->HSMHV_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSMHV_mjswg ) ;
        Qbd += model->HSMHV_pbswg * czbdswg * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Qbd_dT += model->HSMHV_pbswg * czbdswg_dT * (1.0 - arg * sarg) / (1.0 - model->HSMHV_mjswg) ;
        Capbd += czbdswg * sarg ;
      }
    } else {
      T1 = czbd + czbdswg ;
      T1_dT = czbd_dT + czbdswg_dT ;
      T2 = czbd * model->HSMHV_mj / model->HSMHV_pb
        + czbdswg * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      T2_dT = czbd_dT * model->HSMHV_mj / model->HSMHV_pb
        + czbdswg_dT * model->HSMHV_mjswg / model->HSMHV_pbswg ;
      Qbd = vbd_jct * (T1 + vbd_jct * 0.5 * T2) ;
      Qbd_dT = vbd_jct * (T1_dT + vbd_jct * 0.5 * T2_dT) ;
      Capbd = T1 + vbd_jct * T2 ;
    }
  }

  /*-----------------------------------------------------------* 
   * End of PART-4. (label) 
   *-----------------*/ 

/* end_of_part_4: */

  

  /*-----------------------------------------------------------* 
   * PART-5: NQS. (label) 
   *-----------------*/
  if (flg_nqs) {
    if(ckt->CKTmode & MODETRAN){
      if( ckt->CKTmode & MODEINITTRAN ){

        tau  = tau_dVds  = tau_dVgs  = tau_dVbs  = tau_dT  = 0.0 ;
        taub = taub_dVds = taub_dVgs = taub_dVbs = taub_dT = 0.0 ;

      } else {
        /* tau for inversion charge */
        if (flg_noqi == 0) {
          T12 = model->HSMHV_dly1;
          T10 = model->HSMHV_dly2;
            
          T3 = Lch ;
          T1 = T10 * T12 * T3 * T3 ;
          T2 = Mu * VgVt * T12 + T10 * T3 * T3 + small ;
          tau = T1 / T2 ;
            
          T1_dVg = T10 * T12 * 2.0 * T3 * Lch_dVgs ;
          T1_dVd = T10 * T12 * 2.0 * T3 * Lch_dVds ;
          T1_dVb = T10 * T12 * 2.0 * T3 * Lch_dVbs ;
          T1_dT = T10 * T12 * 2.0 * T3 * Lch_dT ; 

          T2_dVg = T12 * Mu_dVgs * VgVt 
            + T12 * Mu * VgVt_dVgs + T10 * 2.0 * T3 * Lch_dVgs ;
          T2_dVd = T12 * Mu_dVds * VgVt 
            + T12 * Mu * VgVt_dVds + T10 * 2.0 * T3 * Lch_dVds ;
          T2_dVb = T12 * Mu_dVbs * VgVt 
            + T12 * Mu * VgVt_dVbs + T10 * 2.0 * T3 * Lch_dVbs ;
          T2_dT = T12 * Mu_dT * VgVt 
            + T12 * Mu * VgVt_dT + T10 * 2.0 * T3 * Lch_dT ; 
          
          T4 = 1.0 / T2 ;
          tau_dVgs = ( T1_dVg - tau * T2_dVg ) * T4 ;
          tau_dVds = ( T1_dVd - tau * T2_dVd ) * T4 ;
          tau_dVbs = ( T1_dVb - tau * T2_dVb ) * T4 ;
          tau_dT =   ( T1_dT  - tau * T2_dT  ) * T4 ; 
       } else {
          tau = model->HSMHV_dly1 ;
          tau_dVgs = tau_dVds = tau_dVbs = tau_dT = 0.0 ;
        }
          
        T1 = ckt->CKTdelta ;

        /* tau for bulk charge */
        T2 = modelMKS->HSMHV_dly3 ;
        taub = T2 * Cox ;
        taub_dVgs = T2 * Cox_dVg ;
        taub_dVds = T2 * Cox_dVd ;
        taub_dVbs = T2 * Cox_dVb ;
        taub_dT   = 0.0          ;

      }
    } else { /* !(CKT_mode & MODETRAN) */
  
      tau  = tau_dVds  = tau_dVgs  = tau_dVbs  = tau_dT  = 0.0 ;
      taub = taub_dVds = taub_dVgs = taub_dVbs = taub_dT = 0.0 ;
    }
  }
    
  if ( flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG)) ) { /* ACNQS */

    if (flg_noqi == 0) {
      T12 = model->HSMHV_dly1 ;
      T10 = model->HSMHV_dly2 ;

      T3 = Lch ;
      T1 = T12 * T10 * T3 * T3 ;
      T2 = Mu * VgVt * T12 + T10 * T3 * T3 + small ;
      tau = T1 / T2 ;

      T1_dVg = T10 * T12 * 2.0 * T3 * Lch_dVgs ;
      T1_dVd = T10 * T12 * 2.0 * T3 * Lch_dVds ;
      T1_dVb = T10 * T12 * 2.0 * T3 * Lch_dVbs ;
      T1_dT = T10 * T12 * 2.0 * T3 * Lch_dT ;
 
      T2_dVg = T12 * Mu_dVgs * VgVt + T12 * Mu * VgVt_dVgs 
        + T10 * 2.0 * T3 * Lch_dVgs ;
      T2_dVd = T12 * Mu_dVds * VgVt + T12 * Mu * VgVt_dVds 
        + T10 * 2.0 * T3 * Lch_dVds ;
      T2_dVb = T12 * Mu_dVbs * VgVt + T12 * Mu * VgVt_dVbs 
        + T10 * 2.0 * T3 * Lch_dVbs ;
      T2_dT = T12 * Mu_dT * VgVt + T12 * Mu * VgVt_dT 
        + T10 * 2.0 * T3 * Lch_dT ; 
     
      T4 = 1.0 / T2 ;
      tau_dVgs = (T1_dVg - tau * T2_dVg) * T4 ;
      tau_dVds = (T1_dVd - tau * T2_dVd) * T4 ;
      tau_dVbs = (T1_dVb - tau * T2_dVb) * T4 ;
      tau_dT = (T1_dT - tau * T2_dT) * T4 ; 
    } else { 
      tau = model->HSMHV_dly1 ;
      tau_dVgs = tau_dVds = tau_dVbs = tau_dT = 0.0 ;
    }
   
    T2 = modelMKS->HSMHV_dly3 ;
    taub = T2 * Cox; 
    taub_dVgs = T2 * Cox_dVg ;
    taub_dVds = T2 * Cox_dVd ;
    taub_dVbs = T2 * Cox_dVb ;
    taub_dT   = 0.0 ;
  }

  /*-----------------------------------------------------------* 
   * End of PART-5. (label) 
   *-----------------*/ 
/* end_of_part_5: */

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-6: Noise Calculation.
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

  /*-----------------------------------------------------------*
   * 1/f noise.
   *-----------------*/
  if ( model->HSMHV_coflick != 0 && !flg_noqi ) {
    
    NFalp = pParam->HSMHV_nfalp ;
    NFtrp = pParam->HSMHV_nftrp ;
    Cit = modelMKS->HSMHV_cit ;

    T1 = Qn0 / C_QE ;
    T2 = ( Cox + Qn0 / ( Ps0 - Vbscl ) + Cit ) * beta_inv / C_QE ;
    T3 = -2.0E0 * Qi / C_QE / Lch / here->HSMHV_weff_nf - T1 ;
    if ( T3 != T1 ) {
      T4 = 1.0E0 / ( T1 + T2 ) / ( T3 + T2 ) + 2.0E0 * NFalp * Ey * Mu / ( T3 - T1 )
        * log( ( T3 + T2 ) / ( T1 + T2 ) ) + NFalp * Ey * Mu * NFalp * Ey * Mu ;
    }  else {
      T4 = 1.0 / ( T1 + T2 ) / ( T3 + T2 ) + 2.0 * NFalp * Ey * Mu / ( T1 + T2 )
        + NFalp * Ey * Mu * NFalp * Ey * Mu;
    }
    Nflic = Ids * Ids * NFtrp / ( Lch * beta * here->HSMHV_weff_nf ) * T4 ;
  } else {
    Nflic = 0.0 ;
  }

  /*-----------------------------------------------------------*
   * thermal noise.
   *-----------------*/
  if ( model->HSMHV_cothrml != 0 && !flg_noqi ) {

    Eyd = ( Psdl - Ps0 ) / Lch + small ;
    T12 = Muun * Eyd / 1.0e7 ;
    /* note: model->HSMHV_bb = 2 (electron) ;1 (hole) */
    if ( 1.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 1.0e0 + epsm10 ) {
      T7  = 1.0e0 ;
    } else if ( 2.0e0 - epsm10 <= model->HSMHV_bb && model->HSMHV_bb <= 2.0e0 + epsm10 ) {
      T7  = T12 ; 
    } else {
      T7  = Fn_Pow( Eyd, model->HSMHV_bb - 1.0e0 ) ;
    }
    T8 = T12 * T7 ;
    T9 = 1.0e0 + T8 ;
    T10 = Fn_Pow( T9, ( - 1.0e0 / model->HSMHV_bb - 1.0e0 ) ) ;
    T11 = T9 * T10 ;
    Mud_hoso = Muun * T11 ;
    Mu_Ave = ( Mu + Mud_hoso ) / 2.0 ;
    
    /* Sid_h = GAMMA * 4.0 * C_KB * model->HSMHV_temp * gds0_h2; */
    T0 = Alpha * Alpha ;
    Nthrml  = here->HSMHV_weff_nf * Cox * VgVt * Mu
      * ( ( 1e0 + 3e0 * Alpha + 6e0 * T0 ) * Mud_hoso * Mud_hoso
          + ( 3e0 + 4e0 * Alpha + 3e0 * T0 ) * Mud_hoso * Mu
          + ( 6e0 + 3e0 * Alpha + T0 ) * Mu * Mu )
      / ( 15e0 * Lch * ( 1e0 + Alpha ) * Mu_Ave * Mu_Ave ) ;
  } else {
    Nthrml = 0e0 ;
  }


  /*----------------------------------------------------------*
   * induced gate noise. ( Part 2/3 )
   *----------------------*/
  if ( model->HSMHV_coign != 0 && model->HSMHV_cothrml != 0 && flg_ign == 1 && !flg_noqi ) {
    sqrtkusaiL = sqrt( kusaiL ) ;
    T2 = VgVt + sqrtkusaiL ;
    T3 = kusai00 * kusai00 ;
    T4 = kusaiL * kusaiL ;
    T5 = 42.0e0 * kusai00 * kusaiL ;
    T5 += 4.0e0 * ( T3 + T4 ) ;
    T5 += 20.0e0 * sqrtkusaiL * VgVt * ( kusai00 + kusaiL ) ;
    T10 = T2 * T2 ;
    T10 *= T10 ;
    kusai_ig = T5 / ( T10 * T2 ) ; /* Induced Gate Noise parameter */
    gds0_ign = here->HSMHV_weff_nf / Lch * Mu * Cox ;
    gds0_h2 = gds0_ign * VgVt ;
    GAMMA = Nthrml / gds0_h2 ;
    T7 = kusai00 + 4.0e0 * VgVt * sqrtkusaiL + kusaiL ;
    /* cross-correlation coefficient (= Sigid/sqrt(Sig*Sid) ) */
    crl_f = c_sqrt_15 * kusai00L * T7
      / ( 6.0e0 * T2 * sqrt( GAMMA * T2 * VgVt * T5 ) ) ;
  }


  /*-----------------------------------------------------------* 
   * End of PART-6. (label) 
   *-----------------*/ 
/* end_of_part_6: */


  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-7: Evaluation of outputs. 
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
    
  /*-----------------------------------------------------------* 
   * Implicit quantities related to Alpha. 
   *-----------------*/ 
  if ( flg_noqi == 0 && VgVt > VgVt_small ) {
    Delta = fac1 * beta / ( 2 * Xi0p12 ) ;
    Pslsat = VgVt / ( 1.0 + Delta ) + Ps0 ;
  } else {
    Pslsat = 0.0 ;
  }
  Vdsat = Pslsat - Pb2 ;
  if ( Vdsat < 0.0 ) {
    Vdsat = 0.0 ;
  }

  /*-----------------------------------------------------------* 
   * Evaluate the derivatives w.r.t. external biases.
   * -> not to be done for the flat and Schur version!
   *-----------------*/ 

  Ids += IdsIBPC ;
  Ids_dVbs  += IdsIBPC_dVbs ;
  Ids_dVds  += IdsIBPC_dVds ;
  Ids_dVgs  += IdsIBPC_dVgs ;
  Ids_dT   += IdsIBPC_dT ;

  /*---------------------------------------------------* 
   * Derivatives of junction diode currents and charges.
   * - NOTE: These quantities are regarded as functions of
   *         external biases.
   * - NOTE: node-base S/D 
   *-----------------*/ 
  Gbse = Gbs ;
  Gbde = Gbd ;
  Capbse = Capbs ;
  Capbde = Capbd ;

  /*---------------------------------------------------*
   * Extrapolate quantities if external biases are out of bounds.
   * -> not necessary here
   *-----------------*/

  /*-----------------------------------------------------------* 
   * Warn negative conductance.
   * - T1 ( = d Ids / d Vds ) is the derivative w.r.t. circuit bias.
   *-----------------*/ 
  T1 = ModeNML * Ids_dVds + ModeRVS * ( Ids_dVbs + Ids_dVds + Ids_dVgs ) ;

  if ( flg_info >= 1 && 
       (Ids_dVbs < 0.0 || T1 < 0.0 || Ids_dVgs < 0.0) ) {
    printf( "*** warning(HiSIM_HV): Negative Conductance\n" ) ;
    printf( " type = %d  mode = %d\n" , model->HSMHV_type , here->HSMHV_mode ) ;
    printf( " Vbs = %12.5e Vds = %12.5e Vgse= %12.5e\n" , 
            Vbs , Vds , Vgs ) ;
    printf( " Ids_dVbs   = %12.5e\n" , Ids_dVbs ) ;
    printf( " Ids_dVds   = %12.5e\n" , T1 ) ;
    printf( " Ids_dVgs   = %12.5e\n" , Ids_dVgs ) ;
  }

 
  /*-----------------------------------------------------------* 
   * Assign outputs.
   *-----------------*/


  /*---------------------------------------------------*
   * Channel current and conductances. 
   *-----------------*/
  here->HSMHV_ids  = Mfactor * Ids ;
  here->HSMHV_dIds_dVdsi = Mfactor * Ids_dVds ;
  here->HSMHV_dIds_dVgsi = Mfactor * Ids_dVgs ;
  here->HSMHV_dIds_dVbsi = Mfactor * Ids_dVbs ;
  here->HSMHV_dIds_dTi   = Mfactor * Ids_dT  ;


  /* -------------------------------------*
   * Intrinsic charges / capacitances.
   *-----------------*/
  if (flg_nqs) { /* for flat handling of NQS: the NQS charges are added in hsmhvld */
    
    here->HSMHV_qg = 0.0 ;
    here->HSMHV_qd = 0.0 ;
    here->HSMHV_qs = 0.0 ;

    here->HSMHV_qdp = 0.0 ;
    here->HSMHV_qsp = 0.0 ;
    
    here->HSMHV_dqdp_dVdse = 0.0 ;
    here->HSMHV_dqdp_dVgse = 0.0 ;
    here->HSMHV_dqdp_dVbse = 0.0 ;
    here->HSMHV_dqdp_dTi   = 0.0 ;
    here->HSMHV_dqsp_dVdse = 0.0 ;
    here->HSMHV_dqsp_dVgse = 0.0 ;
    here->HSMHV_dqsp_dVbse = 0.0 ;
    here->HSMHV_dqsp_dTi   = 0.0 ;

    here->HSMHV_dQdi_dVdsi = 0.0 ;
    here->HSMHV_dQdi_dVgsi = 0.0 ;
    here->HSMHV_dQdi_dVbsi = 0.0 ;
    here->HSMHV_dQdi_dTi   = 0.0 ;
    here->HSMHV_dQg_dVdsi  = 0.0 ;
    here->HSMHV_dQg_dVgsi  = 0.0 ;
    here->HSMHV_dQg_dVbsi  = 0.0 ;
    here->HSMHV_dQg_dTi    = 0.0 ;
    here->HSMHV_dQb_dVdsi  = 0.0 ;
    here->HSMHV_dQb_dVgsi  = 0.0 ;
    here->HSMHV_dQb_dVbsi  = 0.0 ;
    here->HSMHV_dQb_dTi    = 0.0 ;

    here->HSMHV_qgext = 0.0 ;
    here->HSMHV_qdext = 0.0 ;
    here->HSMHV_qsext = 0.0 ;
    
    here->HSMHV_dQdext_dVdse = 0.0 ;
    here->HSMHV_dQdext_dVgse = 0.0 ;
    here->HSMHV_dQdext_dVbse = 0.0 ;
    here->HSMHV_dQdext_dTi   = 0.0 ;
    here->HSMHV_dQgext_dVdse = 0.0 ;
    here->HSMHV_dQgext_dVgse = 0.0 ;
    here->HSMHV_dQgext_dVbse = 0.0 ;
    here->HSMHV_dQgext_dTi   = 0.0 ;
    here->HSMHV_dQbext_dVdse = 0.0 ;
    here->HSMHV_dQbext_dVgse = 0.0 ;
    here->HSMHV_dQbext_dVbse = 0.0 ;
    here->HSMHV_dQbext_dTi   = 0.0 ;
    here->HSMHV_tau       = tau ;
    here->HSMHV_tau_dVgsi = tau_dVgs ;
    here->HSMHV_tau_dVdsi = tau_dVds ;
    here->HSMHV_tau_dVbsi = tau_dVbs ;
    here->HSMHV_tau_dTi   = tau_dT   ;
      
    here->HSMHV_taub       = taub ;
    here->HSMHV_taub_dVgsi = taub_dVgs ;
    here->HSMHV_taub_dVdsi = taub_dVds ;
    here->HSMHV_taub_dVbsi = taub_dVbs ;
    here->HSMHV_taub_dTi   = taub_dT   ;
      
    here->HSMHV_Xd       = Qdrat;
    here->HSMHV_Xd_dVgsi = Qdrat_dVgs ;
    here->HSMHV_Xd_dVdsi = Qdrat_dVds ;
    here->HSMHV_Xd_dVbsi = Qdrat_dVbs ;
    here->HSMHV_Xd_dTi   = Qdrat_dT   ;

    here->HSMHV_Qbulk       = Mfactor * Qb ;
    here->HSMHV_Qbulk_dVgsi = Mfactor * Qb_dVgs ;
    here->HSMHV_Qbulk_dVdsi = Mfactor * Qb_dVds ;
    here->HSMHV_Qbulk_dVbsi = Mfactor * Qb_dVbs ;
    here->HSMHV_Qbulk_dTi   = Mfactor * Qb_dT  ;
    
    here->HSMHV_Qi       = Mfactor * Qi ;
    here->HSMHV_Qi_dVgsi = Mfactor * Qi_dVgs ;
    here->HSMHV_Qi_dVdsi = Mfactor * Qi_dVds ;
    here->HSMHV_Qi_dVbsi = Mfactor * Qi_dVbs ;
    here->HSMHV_Qi_dTi   = Mfactor * Qi_dT  ;

  } else  { /* QS */

    here->HSMHV_qg = Mfactor * - (Qb + Qi) ;
    here->HSMHV_qd = Mfactor * Qd ;
    here->HSMHV_qs = Mfactor * ( Qi - Qd ) ;
    here->HSMHV_qdp = 0.0 ;
    here->HSMHV_qsp = 0.0 ;
    
    here->HSMHV_dqdp_dVdse = 0.0 ;
    here->HSMHV_dqdp_dVgse = 0.0 ;
    here->HSMHV_dqdp_dVbse = 0.0 ;
    here->HSMHV_dqdp_dTi   = 0.0 ;
    here->HSMHV_dqsp_dVdse = 0.0 ;
    here->HSMHV_dqsp_dVgse = 0.0 ;
    here->HSMHV_dqsp_dVbse = 0.0 ;
    here->HSMHV_dqsp_dTi   = 0.0 ;

    here->HSMHV_qgext = 0.0 ;
    here->HSMHV_qdext = 0.0 ;
    here->HSMHV_qsext = 0.0 ;
    
    here->HSMHV_dQdext_dVdse = 0.0 ;
    here->HSMHV_dQdext_dVgse = 0.0 ;
    here->HSMHV_dQdext_dVbse = 0.0 ;
    here->HSMHV_dQdext_dTi   = 0.0 ;
    here->HSMHV_dQgext_dVdse = 0.0 ;
    here->HSMHV_dQgext_dVgse = 0.0 ;
    here->HSMHV_dQgext_dVbse = 0.0 ;
    here->HSMHV_dQgext_dTi   = 0.0 ;
    here->HSMHV_dQbext_dVdse = 0.0 ;
    here->HSMHV_dQbext_dVgse = 0.0 ;
    here->HSMHV_dQbext_dVbse = 0.0 ;
    here->HSMHV_dQbext_dTi   = 0.0 ;
    here->HSMHV_dQdi_dVdsi = Mfactor * Qd_dVds ;
    here->HSMHV_dQdi_dVgsi = Mfactor * Qd_dVgs ;
    here->HSMHV_dQdi_dVbsi = Mfactor * Qd_dVbs ;
    here->HSMHV_dQdi_dTi   = Mfactor * Qd_dT  ;
    here->HSMHV_dQg_dVdsi  = Mfactor * ( - Qb_dVds - Qi_dVds ) ;
    here->HSMHV_dQg_dVgsi  = Mfactor * ( - Qb_dVgs - Qi_dVgs ) ;
    here->HSMHV_dQg_dVbsi  = Mfactor * ( - Qb_dVbs - Qi_dVbs ) ;
    here->HSMHV_dQg_dTi    = Mfactor * ( - Qb_dT  - Qi_dT  ) ;
    here->HSMHV_dQb_dVdsi  = Mfactor * Qb_dVds ;
    here->HSMHV_dQb_dVgsi  = Mfactor * Qb_dVgs ;
    here->HSMHV_dQb_dVbsi  = Mfactor * Qb_dVbs ;
    here->HSMHV_dQb_dTi    = Mfactor * Qb_dT  ;
  }

  /*---------------------------------------------------*
   * Add S/D overlap charges/capacitances to intrinsic ones.
   * - NOTE: This function depends on coadov, a control option.
   *-----------------*/
  if ( model->HSMHV_coadov == 1 ) {
      here->HSMHV_qg += Mfactor * ( Qgod + Qgos + Qgbo + Qy - Qovd - Qovs ) ;
      here->HSMHV_qd += Mfactor * ( - Qgod - Qy + QbdLD ) ;
      here->HSMHV_qs += Mfactor * ( - Qgos      + QbsLD ) ;

      here->HSMHV_qdp += Mfactor * ( - Qfd - Qgdo ) ;
      here->HSMHV_qsp += Mfactor * ( - Qfs - Qgso ) ;

      here->HSMHV_cddo        = Mfactor * ( - Qgod_dVds - Qy_dVds + QbdLD_dVds ) ;
      here->HSMHV_dQdi_dVdsi +=     here->HSMHV_cddo ;
      here->HSMHV_cdgo        = Mfactor * ( - Qgod_dVgs - Qy_dVgs + QbdLD_dVgs ) ;
      here->HSMHV_dQdi_dVgsi +=     here->HSMHV_cdgo ;
      here->HSMHV_cdbo        = Mfactor * ( - Qgod_dVbs - Qy_dVbs + QbdLD_dVbs ) ;
      here->HSMHV_dQdi_dVbsi +=     here->HSMHV_cdbo ;
      here->HSMHV_dQdi_dTi   += Mfactor * ( - Qgod_dT  - Qy_dT  + QbdLD_dT  ) ;
      here->HSMHV_cgdo        = Mfactor * (   Qgod_dVds + Qgos_dVds + Qgbo_dVds + Qy_dVds - Qovd_dVds - Qovs_dVds  ) ;
      here->HSMHV_dQg_dVdsi  +=     here->HSMHV_cgdo ;
      here->HSMHV_cggo        = Mfactor * (   Qgod_dVgs + Qgos_dVgs + Qgbo_dVgs + Qy_dVgs - Qovd_dVgs - Qovs_dVgs  ) ;
      here->HSMHV_dQg_dVgsi  +=     here->HSMHV_cggo ;
      here->HSMHV_cgbo        = Mfactor * (   Qgod_dVbs + Qgos_dVbs + Qgbo_dVbs + Qy_dVbs - Qovd_dVbs - Qovs_dVbs  ) ;
      here->HSMHV_dQg_dVbsi  +=     here->HSMHV_cgbo ;
      here->HSMHV_dQg_dTi    += Mfactor * (   Qgod_dT  + Qgos_dT  + Qgbo_dT  + Qy_dT  - Qovd_dT  - Qovs_dT   ) ;
      here->HSMHV_cbdo        = Mfactor * ( - Qgbo_dVds + QidLD_dVds + QisLD_dVds ) ;
      here->HSMHV_dQb_dVdsi  +=     here->HSMHV_cbdo ;
      here->HSMHV_cbgo        = Mfactor * ( - Qgbo_dVgs + QidLD_dVgs + QisLD_dVgs ) ;
      here->HSMHV_dQb_dVgsi  +=     here->HSMHV_cbgo ;
      here->HSMHV_cbbo        = Mfactor * ( - Qgbo_dVbs + QidLD_dVbs + QisLD_dVbs ) ;
      here->HSMHV_dQb_dVbsi  +=     here->HSMHV_cbbo ;
      here->HSMHV_dQb_dTi    += Mfactor * ( - Qgbo_dT  + QidLD_dT  + QisLD_dT  ) ;

      /* for fringing capacitances */
      here->HSMHV_dqdp_dVdse += Mfactor * (   Cfd - Qgdo_dVdse ) ;
      here->HSMHV_dqdp_dVgse += Mfactor * ( - Cfd - Qgdo_dVgse ) ;
      here->HSMHV_dqdp_dVbse += Mfactor * (       - Qgdo_dVbse ) ;
      here->HSMHV_dqdp_dTi   += 0.0 ;
      here->HSMHV_dqsp_dVdse += Mfactor * (       - Qgso_dVdse ) ;
      here->HSMHV_dqsp_dVgse += Mfactor * ( - Cfs - Qgso_dVgse ) ;
      here->HSMHV_dqsp_dVbse += Mfactor * (       - Qgso_dVbse ) ;
      here->HSMHV_dqsp_dTi   += 0.0 ;

      here->HSMHV_qgext += Mfactor * ( - Qovdext - Qovsext ) ;
      here->HSMHV_qdext += Mfactor * QbdLDext ;
      here->HSMHV_qsext += Mfactor * QbsLDext ;
    
      here->HSMHV_dQdext_dVdse += Mfactor * ( QbdLDext_dVdse ) ;
      here->HSMHV_dQdext_dVgse += Mfactor * ( QbdLDext_dVgse ) ;
      here->HSMHV_dQdext_dVbse += Mfactor * ( QbdLDext_dVbse ) ;
      here->HSMHV_dQdext_dTi   += Mfactor * ( QbdLDext_dT  ) ;
      here->HSMHV_dQgext_dVdse += Mfactor * ( - Qovdext_dVdse - Qovsext_dVdse  ) ;
      here->HSMHV_dQgext_dVgse += Mfactor * ( - Qovdext_dVgse - Qovsext_dVgse  ) ;
      here->HSMHV_dQgext_dVbse += Mfactor * ( - Qovdext_dVbse - Qovsext_dVbse  ) ;
      here->HSMHV_dQgext_dTi   += Mfactor * ( - Qovdext_dT    - Qovsext_dT ) ;
      here->HSMHV_dQbext_dVdse += Mfactor * ( QidLDext_dVdse + QisLDext_dVdse ) ;
      here->HSMHV_dQbext_dVgse += Mfactor * ( QidLDext_dVgse + QisLDext_dVgse ) ;
      here->HSMHV_dQbext_dVbse += Mfactor * ( QidLDext_dVbse + QisLDext_dVbse ) ;
      here->HSMHV_dQbext_dTi   += Mfactor * ( QidLDext_dT    + QisLDext_dT ) ;

  }
  here->HSMHV_dQsi_dVdsi = - (here->HSMHV_dQdi_dVdsi + here->HSMHV_dQg_dVdsi + here->HSMHV_dQb_dVdsi) ;
  here->HSMHV_dQsi_dVgsi = - (here->HSMHV_dQdi_dVgsi + here->HSMHV_dQg_dVgsi + here->HSMHV_dQb_dVgsi) ;
  here->HSMHV_dQsi_dVbsi = - (here->HSMHV_dQdi_dVbsi + here->HSMHV_dQg_dVbsi + here->HSMHV_dQb_dVbsi) ;
  here->HSMHV_dQsi_dTi   = - (here->HSMHV_dQdi_dTi   + here->HSMHV_dQg_dTi   + here->HSMHV_dQb_dTi  ) ;

  here->HSMHV_dQsext_dVdse = - (here->HSMHV_dQdext_dVdse + here->HSMHV_dQgext_dVdse + here->HSMHV_dQbext_dVdse) ;
  here->HSMHV_dQsext_dVgse = - (here->HSMHV_dQdext_dVgse + here->HSMHV_dQgext_dVgse + here->HSMHV_dQbext_dVgse) ;
  here->HSMHV_dQsext_dVbse = - (here->HSMHV_dQdext_dVbse + here->HSMHV_dQgext_dVbse + here->HSMHV_dQbext_dVbse) ;
  here->HSMHV_dQsext_dTi   = - (here->HSMHV_dQdext_dTi   + here->HSMHV_dQgext_dTi   + here->HSMHV_dQbext_dTi  ) ;

  /*---------------------------------------------------* 
   * Substrate/gate/leak currents.
   *-----------------*/ 

  here->HSMHV_isub = Mfactor * Isub ;
  here->HSMHV_dIsub_dVdsi = Mfactor * Isub_dVds ;
  here->HSMHV_dIsub_dVgsi = Mfactor * Isub_dVgs ;
  here->HSMHV_dIsub_dVbsi = Mfactor * Isub_dVbs ;
  here->HSMHV_dIsub_dTi   = Mfactor * Isub_dT  ;
  here->HSMHV_dIsub_dVdse = Mfactor * Isub_dVdse ; 
  
  here->HSMHV_igb   = Mfactor * -Igb ;
  here->HSMHV_dIgb_dVdsi = - Mfactor * Igb_dVds ;
  here->HSMHV_dIgb_dVgsi = - Mfactor * Igb_dVgs ;
  here->HSMHV_dIgb_dVbsi = - Mfactor * Igb_dVbs ;
  here->HSMHV_dIgb_dTi   = - Mfactor * Igb_dT  ;

  if (here->HSMHV_mode == HiSIM_NORMAL_MODE) {
    here->HSMHV_igd   = Mfactor * ( model->HSMHV_glpart1 * Igate - Igd ) ;
  } else {
    here->HSMHV_igd   = Mfactor * ( (1.0e0 - model->HSMHV_glpart1 ) * Igate - Igs ) ;
  }
  
  if (here->HSMHV_mode == HiSIM_NORMAL_MODE) {
    here->HSMHV_igs   = Mfactor * ( (1.0e0  - model->HSMHV_glpart1) * Igate - Igs ) ;
  } else {
    here->HSMHV_igs   = Mfactor * ( model->HSMHV_glpart1 * Igate - Igd ) ;
  }

  /* note: here->HSMHV_igd and here->HSMHV_igs are already subjected to mode handling,
     while the following derivatives here->HSMHV_dIgd_dVdsi, ... are not! */
  here->HSMHV_dIgd_dVdsi = Mfactor * ( model->HSMHV_glpart1 * Igate_dVds - Igd_dVds ) ;
  here->HSMHV_dIgd_dVgsi = Mfactor * ( model->HSMHV_glpart1 * Igate_dVgs - Igd_dVgs ) ;
  here->HSMHV_dIgd_dVbsi = Mfactor * ( model->HSMHV_glpart1 * Igate_dVbs - Igd_dVbs ) ;
  here->HSMHV_dIgd_dTi   = Mfactor * ( model->HSMHV_glpart1 * Igate_dT  - Igd_dT  ) ;
  here->HSMHV_dIgs_dVdsi = Mfactor * ( (1.0 - model->HSMHV_glpart1) * Igate_dVds - Igs_dVds ) ;
  here->HSMHV_dIgs_dVgsi = Mfactor * ( (1.0 - model->HSMHV_glpart1) * Igate_dVgs - Igs_dVgs ) ;
  here->HSMHV_dIgs_dVbsi = Mfactor * ( (1.0 - model->HSMHV_glpart1) * Igate_dVbs - Igs_dVbs ) ;
  here->HSMHV_dIgs_dTi   = Mfactor * ( (1.0 - model->HSMHV_glpart1) * Igate_dT  - Igs_dT  ) ;
  
  here->HSMHV_igidl    = Mfactor * Igidl ;
  here->HSMHV_dIgidl_dVdsi = Mfactor * Igidl_dVds ;
  here->HSMHV_dIgidl_dVgsi = Mfactor * Igidl_dVgs ;
  here->HSMHV_dIgidl_dVbsi = Mfactor * Igidl_dVbs ;
  here->HSMHV_dIgidl_dTi   = Mfactor * Igidl_dT  ;
  
  here->HSMHV_igisl    = Mfactor * Igisl ;
  here->HSMHV_dIgisl_dVdsi = Mfactor * Igisl_dVds ;
  here->HSMHV_dIgisl_dVgsi = Mfactor * Igisl_dVgs ;
  here->HSMHV_dIgisl_dVbsi = Mfactor * Igisl_dVbs ;
  here->HSMHV_dIgisl_dTi   = Mfactor * Igisl_dT  ;
  
  /*---------------------------------------------------* 
   * Von, Vdsat.
   *-----------------*/ 
  here->HSMHV_von = Vth ;
  here->HSMHV_vdsat = Vdsat ;

  

  /*---------------------------------------------------* 
   * Junction diode.
   *-----------------*/ 
  here->HSMHV_ibs = Mfactor * Ibs ;
  here->HSMHV_ibd = Mfactor * Ibd ;
  here->HSMHV_gbs = Mfactor * Gbse ;
  here->HSMHV_gbd = Mfactor * Gbde ;
  *(ckt->CKTstate0 + here->HSMHVqbs) = Mfactor * Qbs ;
  *(ckt->CKTstate0 + here->HSMHVqbd) = Mfactor * Qbd ;
  here->HSMHV_capbs = Mfactor * Capbse ;
  here->HSMHV_capbd = Mfactor * Capbde ;

  here->HSMHV_gbdT = Mfactor * Ibd_dT ;
  here->HSMHV_gbsT = Mfactor * Ibs_dT ;
  here->HSMHV_gcbdT = Mfactor * Qbd_dT ;
  here->HSMHV_gcbsT = Mfactor * Qbs_dT ;

  /*---------------------------------------------------*
   * Add Gjmin (gmin).
   *-----------------*/
  here->HSMHV_ibs += Mfactor * Gjmin * vbs_jct ;
  here->HSMHV_ibd += Mfactor * Gjmin * vbd_jct ;
  here->HSMHV_gbs += Mfactor * Gjmin ;
  here->HSMHV_gbd += Mfactor * Gjmin ;
  
  /*-----------------------------------------------------------* 
   * Warn floating-point exceptions.
   * - Function finite() in libm is called.
   * - Go to start with info==5.
   *-----------------*/
  T1 = here->HSMHV_ids + here->HSMHV_dIds_dVdsi + here->HSMHV_dIds_dVgsi + here->HSMHV_dIds_dVbsi ;
  T1 = T1 + here->HSMHV_qd - (here->HSMHV_dQdi_dVdsi + here->HSMHV_dQdi_dVgsi + here->HSMHV_dQdi_dVbsi) ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf (stderr ,
             "*** warning(HiSIM_HV): FP-exception (PART-1)\n" ) ;
    if ( flg_info >= 1 ) {
      printf ( "*** warning(HiSIM_HV): FP-exception\n") ;
      printf ( "Ids   = %e\n" , here->HSMHV_ids ) ;
      printf ( "Gmbs  = %e\n" , here->HSMHV_dIds_dVbsi ) ;
      printf ( "Gds   = %e\n" , here->HSMHV_dIds_dVdsi ) ;
      printf ( "Gm    = %e\n" , here->HSMHV_dIds_dVgsi ) ;
      printf ( "Qd    = %e\n" , here->HSMHV_qd  ) ;
      printf ( "Cds   = %e\n" , -(here->HSMHV_dQdi_dVdsi
                                                 + here->HSMHV_dQdi_dVgsi
                                                 + here->HSMHV_dQdi_dVbsi) ) ;
    }
  }
  
  T1 = here->HSMHV_isub + here->HSMHV_dIsub_dVbsi + here->HSMHV_dIsub_dVdsi + here->HSMHV_dIsub_dVgsi ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf (stderr ,
             "*** warning(HiSIM_HV): FP-exception (PART-2)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM_HV): FP-exception\n") ;
    }
  }
  
  T1 = here->HSMHV_dQg_dVdsi + here->HSMHV_dQg_dVgsi + here->HSMHV_dQg_dVbsi ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf(stderr ,
            "*** warning(HiSIM_HV): FP-exception (PART-3)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM_HV): FP-exception\n") ;
    }
  }
  
  T1 = here->HSMHV_ibs + here->HSMHV_ibd + here->HSMHV_gbs + here->HSMHV_gbd ;
  T1 = T1 + *(ckt->CKTstate0 + here->HSMHVqbs) + *(ckt->CKTstate0 + here->HSMHVqbd) + here->HSMHV_capbs + here->HSMHV_capbd ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf(stderr ,
            "*** warning(HiSIM_HV): FP-exception (PART-4)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM_HV): FP-exception\n") ;
    }
  }

  /*-----------------------------------------------------------* 
   * Exit for error case.
   *-----------------*/ 
  if ( flg_err != 0 ) {
    fprintf (stderr , "----- bias information (HiSIM_HV)\n" ) ;
    fprintf (stderr , "name: %s\n" , here->HSMHVname ) ;
    fprintf (stderr , "states: %d\n" , here->HSMHVstates ) ;
    fprintf (stderr , "Vdse= %.3e Vgse=%.3e Vbse=%.3e\n"
            , Vdse , Vgse , Vbse ) ;
    fprintf (stderr , "Vdsi= %.3e Vgsi=%.3e Vbsi=%.3e\n"
            , Vds , Vgs , Vbs ) ;
    fprintf (stderr , "vbs_jct= %12.5e vbd_jct= %12.5e\n"
            , vbs_jct , vbd_jct ) ;
    fprintf (stderr , "vd= %.3e vs= %.3e vdp= %.3e vgp= %.3e vbp= %.3e vsp= %.3e\n" 
            , *( ckt->CKTrhsOld + here->HSMHVdNode ) 
            , *( ckt->CKTrhsOld + here->HSMHVsNode ) 
            , *( ckt->CKTrhsOld + here->HSMHVdNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVgNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVbNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSMHVsNodePrime )  ) ;
    fprintf (stderr , "----- bias information (end)\n" ) ;
    return ( HiSIM_ERROR ) ;
  }


  /*-----------------------------------------------------------* 
   * Noise.
   *-----------------*/ 
  here->HSMHV_noiflick = Mfactor * Nflic ;
  here->HSMHV_noithrml = Mfactor * Nthrml ;

  /*----------------------------------------------------------*
   * induced gate noise. ( Part 3/3 )
   *----------------------*/
  if ( model->HSMHV_coign != 0 && model->HSMHV_cothrml != 0 && flg_ign == 1 && !flg_noqi ) {
    T0 = Cox_small * Cox * here->HSMHV_weff_nf * Leff ;
    T1 = -( here->HSMHV_dQg_dVdsi + here->HSMHV_dQg_dVgsi + here->HSMHV_dQg_dVbsi ) / Mfactor ; /* NQS case is not supported. */
    if( - T1 > T0 ){
      Nign0 = c_16o135 * C_QE * beta_inv * T1 * T1 / gds0_ign ;
      if ( kusai00L > epsm10 && Vds > epsm10 ) {
    MuModA = Muun / Mu ;
    MuModB = ( Muun / Mud_hoso - MuModA ) / Vds ;
    correct_w1 = MuModA + C_2o3 * MuModB
      * ( kusai00 + VgVt * sqrtkusaiL + kusaiL )
      / ( VgVt + sqrtkusaiL ) ;
      } else {
    correct_w1 = Muun / Mud_hoso ;
      }
      here->HSMHV_noiigate = Mfactor * Nign0 * kusai_ig * correct_w1 ;
      here->HSMHV_noicross = crl_f ;
      if ( here->HSMHV_noiigate < 0.0 ) here->HSMHV_noiigate = 0.0e0 ;
    }else{
      here->HSMHV_noiigate = 0.0e0 ;
      here->HSMHV_noicross = 0.0e0 ;
    }
  }else{
    here->HSMHV_noiigate = 0.0e0 ;
    here->HSMHV_noicross = 0.0e0 ;
  }
  
  /*-----------------------------------------------------------* 
   * Store values for next calculation.
   *-----------------*/ 

  /* Internal biases */
  if ( here->HSMHV_called >= 1 ) {
    here->HSMHV_vbsc_prv2 = here->HSMHV_vbsc_prv ;
    here->HSMHV_vdsc_prv2 = here->HSMHV_vdsc_prv ;
    here->HSMHV_vgsc_prv2 = here->HSMHV_vgsc_prv ;
    here->HSMHV_mode_prv2 = here->HSMHV_mode_prv ;
  }
  here->HSMHV_vbsc_prv = Vbs ;
  here->HSMHV_vdsc_prv = Vds ;
  here->HSMHV_vgsc_prv = Vgs ;
  here->HSMHV_mode_prv = here->HSMHV_mode ;
  here->HSMHV_temp_prv = TTEMP ;
  
  /* Surface potentials and derivatives w.r.t. internal biases */
  if ( here->HSMHV_called >= 1 ) {
    here->HSMHV_ps0_prv2 = here->HSMHV_ps0_prv ;
    here->HSMHV_ps0_dvbs_prv2 = here->HSMHV_ps0_dvbs_prv ;
    here->HSMHV_ps0_dvds_prv2 = here->HSMHV_ps0_dvds_prv ;
    here->HSMHV_ps0_dvgs_prv2 = here->HSMHV_ps0_dvgs_prv ;
    here->HSMHV_pds_prv2 = here->HSMHV_pds_prv ;
    here->HSMHV_pds_dvbs_prv2 = here->HSMHV_pds_dvbs_prv ;
    here->HSMHV_pds_dvds_prv2 = here->HSMHV_pds_dvds_prv ;
    here->HSMHV_pds_dvgs_prv2 = here->HSMHV_pds_dvgs_prv ;
  }

  here->HSMHV_ps0_prv = Ps0 ;
  here->HSMHV_ps0_dvbs_prv = Ps0_dVbs ;
  here->HSMHV_ps0_dvds_prv = Ps0_dVds ;
  here->HSMHV_ps0_dvgs_prv = Ps0_dVgs ;
  here->HSMHV_ps0_dtemp_prv = Ps0_dT ;
  here->HSMHV_pds_prv = Pds ;
  here->HSMHV_pds_dvbs_prv = Pds_dVbs ;
  here->HSMHV_pds_dvds_prv = Pds_dVds ;
  here->HSMHV_pds_dvgs_prv = Pds_dVgs ;
  here->HSMHV_pds_dtemp_prv = Pds_dT ;


  /* derivatives of channel current w.r.t. external bias (only due to Ra-dependencies!) */
  here->HSMHV_dIds_dVdse = Ids_dRa * Ra_dVdse * Mfactor ;
  here->HSMHV_dIds_dVgse = Ids_dRa * Ra_dVgse * Mfactor ;
  here->HSMHV_dIds_dVbse = Ids_dRa * Ra_dVbse * Mfactor ;

  if ( VdseModeNML > 0.0 ) {
    here->HSMHV_Rd = Rd / Mfactor ;
    here->HSMHV_dRd_dVdse = Rd_dVdse / Mfactor ;
    here->HSMHV_dRd_dVgse = Rd_dVgse / Mfactor ;
    here->HSMHV_dRd_dVbse = Rd_dVbse / Mfactor ;
    here->HSMHV_dRd_dVsubs = Rd_dVsubs / Mfactor ;
    here->HSMHV_dRd_dTi   = Rd_dT / Mfactor ;
    here->HSMHV_Rs = Rs / Mfactor ;
    here->HSMHV_dRs_dVdse = Rs_dVdse / Mfactor ;
    here->HSMHV_dRs_dVgse = Rs_dVgse / Mfactor ;
    here->HSMHV_dRs_dVbse = Rs_dVbse / Mfactor ;
    here->HSMHV_dRs_dVsubs = Rs_dVsubs / Mfactor ;
    here->HSMHV_dRs_dTi   = Rs_dT / Mfactor ;
  } else {
    here->HSMHV_Rd = Rs / Mfactor ;
    here->HSMHV_dRd_dVdse = - ( Rs_dVdse + Rs_dVgse + Rs_dVbse + Rs_dVsubs ) / Mfactor ;
    here->HSMHV_dRd_dVgse = Rs_dVgse / Mfactor ;
    here->HSMHV_dRd_dVbse = Rs_dVbse / Mfactor ;
    here->HSMHV_dRd_dVsubs = Rs_dVsubs / Mfactor ;
    here->HSMHV_dRd_dTi   = Rs_dT / Mfactor ;
    here->HSMHV_Rs = Rd / Mfactor ;
    here->HSMHV_dRs_dVdse = - ( Rd_dVdse + Rd_dVgse + Rd_dVbse + Rd_dVsubs ) / Mfactor ;
    here->HSMHV_dRs_dVgse = Rd_dVgse / Mfactor ;
    here->HSMHV_dRs_dVbse = Rd_dVbse / Mfactor ;
    here->HSMHV_dRs_dVsubs = Rd_dVsubs / Mfactor ;
    here->HSMHV_dRs_dTi   = Rd_dT / Mfactor ;
  }
  /* Clamping to Res_min */
  if(here->HSMHV_Rd < Res_min) {
     here->HSMHV_Rd = Res_min ;
     here->HSMHV_dRd_dVdse  = 0.0 ;
     here->HSMHV_dRd_dVgse  = 0.0 ;
     here->HSMHV_dRd_dVbse  = 0.0 ;
     here->HSMHV_dRd_dVsubs = 0.0 ;
     here->HSMHV_dRd_dTi    = 0.0 ;
  }
  if(here->HSMHV_Rs < Res_min) {
     here->HSMHV_Rs = Res_min ;
     here->HSMHV_dRs_dVdse  = 0.0 ;
     here->HSMHV_dRs_dVgse  = 0.0 ;
     here->HSMHV_dRs_dVbse  = 0.0 ;
     here->HSMHV_dRs_dVsubs = 0.0 ;
     here->HSMHV_dRs_dTi    = 0.0 ;
  }

  /*-----------------------------------------------------------*
   * End of PART-7. (label) 
   *-----------------*/ 
/* end_of_part_7: */

  /*-----------------------------------------------------------* 
   * Bottom of hsmhveval. 
   *-----------------*/ 

  return ( HiSIM_OK ) ;
  
} /* end of hsmhveval */
