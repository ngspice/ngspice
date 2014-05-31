/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) 
 
 FILE : hsm2eval.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

/**********************************************************************

The following source code, and all copyrights, trade secrets or other
intellectual property rights in and to the source code in its entirety,
is owned by the Hiroshima University and the STARC organization.

All users need to follow the "HiSIM2 Distribution Statement and
Copyright Notice" attached to HiSIM2 model.

-----HiSIM2 Distribution Statement and Copyright Notice--------------

Software is distributed as is, completely without warranty or service
support. Hiroshima University or STARC and its employees are not liable
for the condition or performance of the software.

Hiroshima University and STARC own the copyright and grant users a perpetual,
irrevocable, worldwide, non-exclusive, royalty-free license with respect 
to the software as set forth below.   

Hiroshima University and STARC hereby disclaim all implied warranties.

Hiroshima University and STARC grant the users the right to modify, copy,
and redistribute the software and documentation, both within the user's
organization and externally, subject to the following restrictions

1. The users agree not to charge for Hiroshima University and STARC code
itself but may charge for additions, extensions, or support.

2. In any product based on the software, the users agree to acknowledge
Hiroshima University and STARC that developed the software. This
acknowledgment shall appear in the product documentation.

3. The users agree to reproduce any copyright notice which appears on
the software on any copy or modification of such made available
to others."


*************************************************************************/

/*********************************************************************
* Memorandum on programming
* 
* (1) Bias (x: b|d|g)
*     . vxs : Input argument.
*     . Vxse: External bias taking account device type (pMOS->nMOS).
*     . Vxsc: Confined bias within a specified region. 
*     . Vxs : Internal bias.
*     . Y_dVxs denotes the partial derivative of Y w.r.t. Vxs.
* 
* (2) Device Mode
*     . Normal mode (Vds>0 for nMOS) is assumed.
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
#include "ngspice/cktdefs.h"

/*-----------------------------------*
* HiSIM macros
*-----------------*/
#include "hisim2.h"
#include "hsm2evalenv.h"

/*-----------------------------------*
* HiSIM constants
*-----------------*/
#define C_sce_dlt (1.0e-2)
#define C_gidl_delta 0.5
#define C_PSLK_DELTA 1e-3     /* delta for Pslk smoothing */
#define C_PSLK_SHIFT 1.0      /* constant value for temporary shift */
#define C_IDD_MIN    1.0e-15

/* local variables used in macro functions */
/*===========================================================*
* pow
*=================*/
#ifdef POW_TO_EXP_AND_LOG
#define Fn_Pow( x , y )  exp( (y) * log( x )  ) 
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
* Ceiling, smoothing functions.
*=================*/
/*---------------------------------------------------*
* smoothUpper: ceiling.
*      y = xmax - 0.5 ( arg + sqrt( arg^2 + 4 xmax delta ) )
*    arg = xmax - x - delta
*-----------------*/

#define Fn_SU( y , x , xmax , delta , dx ) { \
    double TMF1, TMF2; \
    TMF1 = ( xmax ) - ( x ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmax ) * ( delta) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    y = ( xmax ) - 0.5 * ( TMF1 + TMF2 ) ; \
  }

#define Fn_SU2( y , x , xmax , delta , dy_dx , dy_dxmax ) { \
    double TMF1, TMF2; \
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
    double TMF1, TMF2; \
    TMF1 = ( x ) - ( xmin ) - ( delta ) ; \
    TMF2 = 4.0 * ( xmin ) * ( delta ) ; \
    TMF2 = TMF2 > 0.0 ?  TMF2 : - ( TMF2 ) ; \
    TMF2 = sqrt ( TMF1 * TMF1 + TMF2 ) ; \
    dx = 0.5 * ( 1.0 + TMF1 / TMF2 ) ; \
    y = ( xmin ) + 0.5 * ( TMF1 + TMF2 ) ; \
  }

/*---------------------------------------------------*
* smoothZero: flooring to zero.
*      y = 0.5 ( x + sqrt( x^2 + 4 delta^2 ) )
*-----------------*/

#define Fn_SZ( y , x , delta , dx ) { \
    double TMF2; \
    TMF2 = sqrt ( ( x ) *  ( x ) + 4.0 * ( delta ) * ( delta) ) ; \
    dx = 0.5 * ( 1.0 + ( x ) / TMF2 ) ; \
    y = 0.5 * ( ( x ) + TMF2 ) ; \
    if( y < 0.0 ) { y=0.0; dx=0.0; } \
  }
#if 0
/*---------------------------------------------------*
* smoothZero: flooring to zero.
*      y = 0.5 ( x + sqrt( x^2 + 4 delta^2 ) )
*-----------------*/
static double smoothZero
(
 double x,
 double delta,
 double *dx
 )
{
  double sqr = sqrt ( x * x + 4.0 * delta * delta) ;
  if (dx) *dx = 0.5 * ( 1.0 + x / sqr ) ;
  return 0.5 * ( x + sqr ) ;
}
/*---------------------------------------------------*
* CeilingPow: ceiling for positive x, flooring for negative x.
*      y = x * xmax / ( x^{2m} + xmax^{2m} )^{1/(2m)}
* note:
*   - xmax has to be positive.
*   - -xmax < y < xmax.
*   - dy/dx|_{x=0} = 1.
*-----------------*/
static double CeilingPow
(
 double x,
 double xmax,
 int    pw,
 double *dx
 )
{
  double x2 = x * x ;
  double xmax2 = xmax * xmax ;
  double xp = 1.0 , xmp = 1.0 ;
  int   m , mm ;
  double arg , dnm ;
  double result ;

  for ( m = 0 ; m < pw ; m ++ ) {
    xp *= x2 ;
    xmp *= xmax2 ;
  }
  arg = xp + xmp ;
  dnm = arg ;
  if ( pw == 1 || pw == 2 || pw == 4 || pw == 8 ) {
    if ( pw == 1 ) {
        mm = 1 ;
    } else if ( pw == 2 ) {
        mm = 2 ;
    } else if ( pw == 4 ) {
        mm = 3 ;
    } else if ( pw == 8 ) {
        mm = 4 ;
    }
    for ( m = 0 ; m < mm ; m ++ ) {
        dnm = sqrt( dnm ) ;
    }
  } else {
        dnm = pow( dnm , 1.0 / ( 2.0 * pw ) ) ;
  }
  dnm = 1.0 / dnm ;
  result = x * xmax * dnm ;
  (*dx) = xmax * xmp * dnm / arg ;
  return result ;
}
#endif
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
  int   m , mm ; \
  double arg , dnm ; \
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


/*===========================================================*
* Functions for symmetry.
*=================*/

/*---------------------------------------------------*
* Declining function using a polynomial.
*-----------------*/

#define Fn_DclPoly4( y , x , dx ) { \
  double TMF2, TMF3, TMF4; \
  TMF2 = (x) * (x) ; \
  TMF3 = TMF2 * (x) ; \
  TMF4 = TMF2 * TMF2 ; \
  y = 1.0 / ( 1.0 + (x) + TMF2 + TMF3 + TMF4 ) ; \
  dx = - ( 1.0 + 2.0 * (x) + 3.0 * TMF2 + 4.0 * TMF3 )  * y * y  ; \
} 

/*---------------------------------------------------*
* "smoothUpper" uasing a polynomial
*-----------------*/

#define Fn_SUPoly4( y , x , xmax , dx ) { \
 double TMF1; \
 TMF1 = (x) / xmax ; \
 Fn_DclPoly4( y , TMF1 , dx ) ; \
 y = xmax * ( 1.0 - y ) ; \
 dx = - dx ; \
}
 
/*---------------------------------------------------*
* SymAdd: evaluate additional term for symmetry.
*-----------------*/

#define Fn_SymAdd( y , x , add0 , dx ) \
{ \
    double TMF1, TMF2, TMF3; \
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
* Function hsm2evaluate.
*=================*/
int HSM2evaluate
(
 double        vds,
 double        vgs,
 double        vbs,
 double        vbs_jct,
 double        vbd_jct,
 HSM2instance *here,
 HSM2model    *model,
 CKTcircuit   *ckt
 ) 
{
  HSM2binningParam *pParam = &here->pParam ;
  HSM2modelMKSParam *modelMKS = &model->modelMKS ;
/*  HSM2hereMKSParam *hereMKS = &here->hereMKS ;*/
  /*-----------------------------------*
   * Constants for Smoothing functions
   *---------------*/
  const double vth_dlt = 1.0e-3 ;
  /*  const double cclmmdf = 1.0e-2 ;*/
  const double cclmmdf = 1.0e-1 ;
  const double qme_dlt = 1.0e-4 ;
  const double eef_dlt = 1.0e-2 * C_m2cm ;
  const double sti2_dlt = 2.0e-3 ;
  const double pol_dlt = 5.0e-2 ; 
  const double psisti_dlt = 5.0e-3 ;

  /*---------------------------------------------------*
   * Local variables. 
   *-----------------*/
  /* Constants ----------------------- */
  const int lp_s0_max = 20 ;
  const int lp_sl_max = 20 ;
  int lp_bs_max = 10 ;
  const double Ids_tol = 1.0e-10 ;
  const double Ids_maxvar = 1.0e-1 ;
  const double dP_max  = 0.1e0 ;
  const double ps_conv = 5.0e-13 ;
  /* double  ps_conv = 1.0e-13 ;*/
  const double gs_conv = 1.0e-8 ;
  const double mini_current = 1.0e-15 ;
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
  const double Vbs_min = -10.5e0 ;
  const double epsm10 = 10.0e0 * C_EPS_M ;
  const double small = 1.0e-50 ;
  const double small2 = 1.0e-12 ;  /* for Qover */

  double Vbs_max = 0.8e0 ;
  double Vbs_bnd = 0.4e0 ; /* start point of positive Vbs bending */
  double Gdsmin = 0.0 ;
  double Gjmin = ckt->CKTgmin ;

  /* Internal flags  --------------------*/
  int flg_err = 0 ;  /* error level */
  int flg_rsrd = 0 ; /* Flag for bias loop accounting Rs and Rd */
  int flg_iprv = 0 ; /* Flag for initial guess of Ids */
  int flg_pprv = 0 ; /* Flag for initial guesses of Ps0 and Pds */
  int flg_noqi =0;     /* Flag for the cases regarding Qi=Qd=0 */
  int flg_vbsc = 0 ; /* Flag for Vbs confining */
  int flg_info = 0 ; 
  int flg_conv = 0 ; /* Flag for Poisson loop convergence */
  int flg_qme = 0 ; /* Flag for QME */

  /* flag for NQS calculation */
  int flg_nqs=0 ;
  
  /* Important Variables in HiSIM -------*/
  /* external bias */
  double Vbse =0.0, Vdse =0.0, Vgse =0.0 ;
  /* confine bias */
  double Vbsc =0.0, Vdsc =0.0, Vgsc =0.0 ;
  double Vbsc_dVbse = 1.0 ;
  /* internal bias */
  double Vbs =0.0, Vds =0.0, Vgs =0.0, Vdb =0.0, Vsb =0.0 ;
  double Vbs_dVbse = 1.0 , Vbs_dVdse = 0.0 , Vbs_dVgse = 0.0 ;
  double Vds_dVbse = 0.0 , Vds_dVdse = 1.0 , Vds_dVgse = 0.0 ;
  double Vgs_dVbse = 0.0 , Vgs_dVdse = 0.0 , Vgs_dVgse = 1.0 ;
  double Vgp =0.0 ;
  double Vgp_dVbs =0.0, Vgp_dVds =0.0, Vgp_dVgs =0.0 ;
  double Vgs_fb =0.0 ;
  /* Ps0 : surface potential at the source side */
  double Ps0 =0.0 ;
  double Ps0_dVbs =0.0, Ps0_dVds =0.0, Ps0_dVgs =0.0 ;
  double Ps0_ini =0.0, Ps0_iniA =0.0, Ps0_iniB =0.0 ;
  /* Psl : surface potential at the drain side */
  double Psl =0.0 ;
  double Psl_dVbs =0.0, Psl_dVds =0.0, Psl_dVgs =0.0 ;
  double Psl_lim =0.0, dPlim =0.0 ;
  /* Pds := Psl - Ps0 */
  double Pds = 0.0 ;
  double Pds_dVbs = 0.0, Pds_dVds = 0.0 , Pds_dVgs  = 0.0 ;
  double Pds_ini =0.0 ;
  double Pds_max =0.0 ;
  /* iteration numbers of Ps0 and Psl equations. */
  int lp_s0 = 0 , lp_sl = 0 ;
  /* Xi0 := beta * ( Ps0 - Vbs ) - 1. */
  double Xi0 =0.0 ;
  double Xi0_dVbs =0.0, Xi0_dVds =0.0, Xi0_dVgs =0.0 ;
  double Xi0p12 =0.0 ;
  double Xi0p12_dVbs =0.0, Xi0p12_dVds =0.0, Xi0p12_dVgs =0.0 ;
  double Xi0p32 =0.0 ;
  /* Xil := beta * ( Psl - Vbs ) - 1. */
  double Xilp12 =0.0 ;
  double Xilp32 =0.0 ;
  double Xil =0.0 ;
  /* modified bias and potential for sym.*/
  double Vbsz =0.0, Vdsz =0.0, Vgsz =0.0 ;
  double Vbsz_dVbs =0.0, Vbsz_dVds =0.0 ;
  double Vdsz_dVds =0.0 ;
  double Vgsz_dVgs =0.0, Vgsz_dVds =0.0 ;
  double Vzadd =0.0, Vzadd_dVds =0.0 ;
  double Ps0z =0.0, Ps0z_dVbs =0.0, Ps0z_dVds =0.0, Ps0z_dVgs =0.0 ;
  double Pzadd =0.0, Pzadd_dVbs =0.0, Pzadd_dVds =0.0, Pzadd_dVgs =0.0 ;
  double Vgpz , Vgpz_dVbs , Vgpz_dVds , Vgpz_dVgs ; /* (tmp) */

  /* IBPC */
  double dVbsIBPC =0.0, dVbsIBPC_dVbs =0.0, dVbsIBPC_dVds =0.0, dVbsIBPC_dVgs =0.0 ;
  double betaWL =0.0, betaWL_dVbs =0.0, betaWL_dVds =0.0, betaWL_dVgs =0.0 ;
  double Xi0p32_dVbs =0.0, Xi0p32_dVds =0.0, Xi0p32_dVgs =0.0 ;
  double Xil_dVbs =0.0, Xil_dVds =0.0, Xil_dVgs =0.0 ;
  double Xilp12_dVbs =0.0, Xilp12_dVds =0.0, Xilp12_dVgs =0.0 ;
  double Xilp32_dVbs =0.0, Xilp32_dVds =0.0, Xilp32_dVgs =0.0 ;
  double dG3 =0.0, dG3_dVbs =0.0, dG3_dVds =0.0, dG3_dVgs =0.0 ;
  double dG4 =0.0, dG4_dVbs =0.0, dG4_dVds =0.0, dG4_dVgs =0.0 ;
  double dIdd =0.0, dIdd_dVbs =0.0, dIdd_dVds =0.0, dIdd_dVgs =0.0 ;

  /* Chi := beta * ( Ps{0/l} - Vbs ) */
  double Chi =0.0 ;
  double Chi_dVbs =0.0, Chi_dVds =0.0, Chi_dVgs =0.0 ;
  /* Rho := beta * ( Psl - Vds ) */
  double Rho =0.0 ;
  /* threshold voltage */
  double Vth =0.0 ;
  double Vth0 =0.0 ;
  double Vth0_dVb =0.0, Vth0_dVd =0.0, Vth0_dVg =0.0 ;
  /* variation of threshold voltage */
  double dVth =0.0 ;
  double dVth_dVb =0.0, dVth_dVd =0.0, dVth_dVg =0.0 ;
  double dVth0 =0.0 ;
  double dVth0_dVb =0.0, dVth0_dVd =0.0, dVth0_dVg =0.0 ;
  double dVthSC =0.0 ;
  double dVthSC_dVb =0.0, dVthSC_dVd =0.0, dVthSC_dVg =0.0 ;
  double Pb20b =0.0 ;
  double Pb20b_dVg =0.0, Pb20b_dVb =0.0, Pb20b_dVd =0.0 ;
  double dVthW =0.0 ;
  double dVthW_dVb =0.0, dVthW_dVd =0.0, dVthW_dVg =0.0 ;
  /* Alpha and related parameters */
  double Alpha =0.0 ;
  double Alpha_dVbs =0.0, Alpha_dVds =0.0, Alpha_dVgs =0.0 ;
  double Achi =0.0 ;
  double Achi_dVbs =0.0, Achi_dVds =0.0, Achi_dVgs =0.0 ;
  double VgVt = 0.0 ;
  double VgVt_dVbs = 0.0, VgVt_dVds = 0.0, VgVt_dVgs = 0.0 ;
  double Pslsat = 0.0 ;
  double Vdsat = 0.0 ;
  double VdsatS = 0.0 ;
  double VdsatS_dVbs = 0.0, VdsatS_dVds = 0.0, VdsatS_dVgs = 0.0 ;
  double Delta =0.0 ;
  /* Q_B and capacitances */
  double Qb =0.0, Qb_dVbs =0.0, Qb_dVds =0.0, Qb_dVgs =0.0 ;
  double Qb_dVbse =0.0, Qb_dVdse =0.0, Qb_dVgse =0.0 ;
  double Qbu = 0.0 , Qbu_dVbs = 0.0 , Qbu_dVds = 0.0 , Qbu_dVgs = 0.0 ;
  /* Q_I and capacitances */
  double Qi =0.0, Qi_dVbs =0.0, Qi_dVds =0.0, Qi_dVgs =0.0 ;
  double Qi_dVbse =0.0, Qi_dVdse =0.0, Qi_dVgse =0.0 ;
  double Qiu = 0.0 , Qiu_dVbs = 0.0 , Qiu_dVds = 0.0 , Qiu_dVgs = 0.0 ;
  /* Q_D and capacitances */
  double Qd =0.0, Qd_dVbs =0.0, Qd_dVds =0.0, Qd_dVgs =0.0 ;
  double Qd_dVbse =0.0, Qd_dVdse =0.0, Qd_dVgse =0.0 ;
  double qd_dVgse=0.0, qd_dVdse=0.0, qd_dVbse=0.0, qd_dVsse =0.0 ;
  /* channel current */
  double Ids =0.0 ;
  double Ids_dVbs =0.0, Ids_dVds =0.0, Ids_dVgs =0.0 ;
  double Ids_dVbse =0.0, Ids_dVdse =0.0, Ids_dVgse =0.0 ;
  double Ids0 =0.0 ;
  double Ids0_dVbs =0.0, Ids0_dVds =0.0, Ids0_dVgs =0.0 ;
  /* STI */
  double dVthSCSTI =0.0 ;
  double dVthSCSTI_dVg =0.0, dVthSCSTI_dVd =0.0, dVthSCSTI_dVb =0.0 ;
  double Vgssti =0.0 ;
  double Vgssti_dVbs =0.0, Vgssti_dVds =0.0, Vgssti_dVgs =0.0 ;
  double costi0 =0.0, costi1 =0.0, costi3 =0.0 ;
  double costi4 =0.0, costi5 =0.0, costi6 =0.0, costi7 =0.0 ;
  double costi3_dVb =0.0, costi3_dVd=0.0, costi3_dVg =0.0 ;
  double costi3_dVb_c3 =0.0, costi3_dVd_c3=0.0, costi3_dVg_c3 =0.0 ;
  double Psasti =0.0 ;
  double Psasti_dVbs =0.0, Psasti_dVds =0.0, Psasti_dVgs =0.0 ;
  double Psbsti =0.0 ;
  double Psbsti_dVbs =0.0, Psbsti_dVds =0.0, Psbsti_dVgs =0.0 ;
  double Psab =0.0 ;
  double Psab_dVbs =0.0, Psab_dVds =0.0, Psab_dVgs =0.0 ;
  double Psti =0.0 ;
  double Psti_dVbs =0.0, Psti_dVds =0.0, Psti_dVgs =0.0 ;
  double sq1sti =0.0 ;
  double sq1sti_dVbs =0.0, sq1sti_dVds =0.0, sq1sti_dVgs =0.0 ;
  double sq2sti =0.0 ;
  double sq2sti_dVbs =0.0, sq2sti_dVds =0.0, sq2sti_dVgs =0.0 ;
  double Qn0sti =0.0 ;
  double Qn0sti_dVbs =0.0, Qn0sti_dVds =0.0, Qn0sti_dVgs =0.0 ;
  double Idssti =0.0 ;
  double Idssti_dVbs =0.0, Idssti_dVds =0.0, Idssti_dVgs =0.0 ;
  /* constants */
  double beta =0.0, beta_inv =0.0 ;
  double beta2 =0.0 ;
  double Pb2 =0.0 ;
  double Pb20 =0.0 ;
  double Pb2c =0.0 ;
  double Vfb =0.0 ;
  double c_eox =0.0 ;
  double Leff=0.0, Weff =0.0 ;
  double q_Nsub =0.0 ;
  /* PART-1 */
  /* Accumulation zone */
  double Psa =0.0 ;
  double Psa_dVbs =0.0, Psa_dVds =0.0, Psa_dVgs =0.0 ;
  /* CLM*/
  double Psdl =0.0, Psdl_dVbs =0.0, Psdl_dVds =0.0, Psdl_dVgs =0.0 ;
  double Lred =0.0, Lred_dVbs =0.0, Lred_dVds =0.0, Lred_dVgs =0.0 ;
  double Lch =0.0, Lch_dVbs =0.0, Lch_dVds =0.0, Lch_dVgs =0.0 ;
  double Wd =0.0, Wd_dVbs =0.0, Wd_dVds =0.0, Wd_dVgs =0.0 ;
  double Aclm =0.0 ;
  /* Pocket Implant */
  double Vthp=0.0, Vthp_dVb=0.0, Vthp_dVd=0.0, Vthp_dVg =0.0 ;
  double dVthLP=0.0, dVthLP_dVb=0.0, dVthLP_dVd=0.0, dVthLP_dVg =0.0 ;
  double bs12=0.0, bs12_dVb=0.0, bs12_dVd =0.0, bs12_dVg =0.0 ;
  double Qbmm=0.0, Qbmm_dVb=0.0, Qbmm_dVd =0.0, Qbmm_dVg =0.0 ;
  double dqb=0.0, dqb_dVb=0.0, dqb_dVg=0.0, dqb_dVd =0.0 ;
  double Vdx=0.0, Vdx2 =0.0 ;
  double Pbsum=0.0, sqrt_Pbsum =0.0 ;
  double Pbsum_dVb=0.0, Pbsum_dVd=0.0, Pbsum_dVg =0.0 ;
  /* Poly-Depletion Effect */
  const double pol_b = 1.0 ;
  double dPpg =0.0, dPpg_dVb =0.0, dPpg_dVd =0.0, dPpg_dVg =0.0; 
  /* Quantum Effect */
  double Tox =0.0, Tox_dVb =0.0, Tox_dVd =0.0, Tox_dVg =0.0 ;
  double dTox =0.0, dTox_dVb =0.0, dTox_dVd =0.0, dTox_dVg =0.0 ;
  double Cox =0.0, Cox_dVb =0.0, Cox_dVd =0.0, Cox_dVg =0.0 ;
  double Cox_inv =0.0, Cox_inv_dVb =0.0, Cox_inv_dVd =0.0, Cox_inv_dVg =0.0 ;
  double Tox0 =0.0, Cox0 =0.0, Cox0_inv =0.0 ;
  double Vthq=0.0, Vthq_dVb =0.0, Vthq_dVd =0.0 ;
  /* Igate , Igidl , Igisl */
  const double igate_dlt = 1.0e-2 ;
  double Psdlz =0.0, Psdlz_dVbs =0.0, Psdlz_dVds =0.0, Psdlz_dVgs =0.0 ;
  double Egp12 =0.0, Egp32 =0.0 ;
  double E1 =0.0, E1_dVb =0.0, E1_dVd =0.0, E1_dVg =0.0 ;
  double Qb0Cox =0.0, Qb0Cox_dVb =0.0, Qb0Cox_dVd =0.0, Qb0Cox_dVg =0.0 ;
  double Etun =0.0, Etun_dVbs =0.0, Etun_dVds =0.0, Etun_dVgs =0.0 ;
  double Egidl =0.0, Egidl_dVb =0.0, Egidl_dVd =0.0, Egidl_dVg =0.0 ;
  double Egisl =0.0, Egisl_dVb =0.0, Egisl_dVd =0.0, Egisl_dVg =0.0 ;
  double Igate =0.0, Igate_dVbs =0.0, Igate_dVds =0.0, Igate_dVgs =0.0 ;
  double Igate_dVbse =0.0, Igate_dVdse =0.0, Igate_dVgse =0.0 ;
  double Igs =0.0, Igd =0.0, Igb =0.0 ;
  double Igs_dVbs =0.0, Igs_dVds =0.0, Igs_dVgs =0.0 ;
  double Igs_dVbse =0.0, Igs_dVdse =0.0, Igs_dVgse =0.0 ;
  double Igd_dVbs =0.0, Igd_dVds =0.0, Igd_dVgs =0.0 ;
  double Igd_dVbse =0.0, Igd_dVdse =0.0, Igd_dVgse =0.0 ;
  double Igb_dVbs =0.0, Igb_dVds =0.0, Igb_dVgs =0.0 ;
  double Igb_dVbse =0.0, Igb_dVdse =0.0, Igb_dVgse =0.0 ;
  double Igidl =0.0, Igidl_dVbs =0.0, Igidl_dVds =0.0, Igidl_dVgs =0.0 ;
  double Igidl_dVbse =0.0, Igidl_dVdse =0.0, Igidl_dVgse =0.0 ;
  double Igisl =0.0, Igisl_dVbs =0.0, Igisl_dVds =0.0, Igisl_dVgs =0.0 ;
  double Igisl_dVbse =0.0, Igisl_dVdse =0.0, Igisl_dVgse =0.0 ;
  /* connecting function */
  double FD2 =0.0, FD2_dVbs =0.0, FD2_dVds =0.0, FD2_dVgs =0.0 ;
  double FMDVDS =0.0, FMDVDS_dVbs =0.0, FMDVDS_dVds =0.0, FMDVDS_dVgs =0.0 ;

  double cnst0 =0.0, cnst1 =0.0 ;
  double cnstCoxi =0.0 , cnstCoxi_dVg =0.0 , cnstCoxi_dVd =0.0 , cnstCoxi_dVb =0.0 ;
  double fac1 =0.0 ;
  double fac1_dVbs =0.0, fac1_dVds =0.0, fac1_dVgs =0.0 ;
  double fac1p2 =0.0 ;
  double fs01 =0.0 ;
  double fs01_dPs0 =0.0 ;
  double fs01_dVbs =0.0, fs01_dVds =0.0, fs01_dVgs =0.0 ;
  double fs02 =0.0 ;
  double fs02_dPs0 =0.0 ;
  double fs02_dVbs =0.0, fs02_dVds =0.0, fs02_dVgs =0.0 ;
  double fsl1 =0.0 ;
  double fsl1_dPsl =0.0 ;
  double fsl1_dVbs =0.0, fsl1_dVds =0.0, fsl1_dVgs =0.0; /* Vdseff */
  double fsl2 =0.0 ;
  double fsl2_dPsl =0.0 ;
  double fsl2_dVbs =0.0, fsl2_dVds =0.0, fsl2_dVgs =0.0; /* Vdseff */
  double cfs1 =0.0 ;
  double fb =0.0, fb_dChi =0.0 ;
  double fi =0.0, fi_dChi =0.0 ;
  double exp_Chi =0.0, exp_Rho =0.0, exp_bVbs =0.0, exp_bVbsVds =0.0 ;
  double Fs0=0.0, Fsl =0.0 ;
  double Fs0_dPs0 =0.0, Fsl_dPsl =0.0 ;
  double dPs0 =0.0, dPsl =0.0 ;
  double Qn0 = 0.0e0 ;
  double Qn0_dVbs =0.0, Qn0_dVds =0.0, Qn0_dVgs =0.0 ;
  double Qb0 =0.0 ;
  double Qb0_dVb =0.0, Qb0_dVd =0.0, Qb0_dVg =0.0 ;
  double Qbnm =0.0 ;
  double Qbnm_dVbs =0.0, Qbnm_dVds =0.0, Qbnm_dVgs =0.0 ;
  double DtPds =0.0 ;
  double DtPds_dVbs =0.0, DtPds_dVds =0.0, DtPds_dVgs =0.0 ;
  double Qinm =0.0 ;
  double Qinm_dVbs =0.0, Qinm_dVds =0.0, Qinm_dVgs =0.0 ;
  double Qidn =0.0 ;
  double Qidn_dVbs =0.0, Qidn_dVds =0.0, Qidn_dVgs =0.0 ;
  double Qdnm =0.0 ;
  double Qdnm_dVbs =0.0, Qdnm_dVds =0.0, Qdnm_dVgs =0.0 ;
  double Qddn =0.0 ;
  double Qddn_dVbs =0.0, Qddn_dVds =0.0, Qddn_dVgs =0.0 ;
  double Quot =0.0 ;
  double Qdrat = 0.5 ;
  double Qdrat_dVbs = 0.0 , Qdrat_dVds = 0.0, Qdrat_dVgs = 0.0 ;
  double Qdrat_dVbse =0.0, Qdrat_dVdse =0.0, Qdrat_dVgse =0.0 ;
  double Idd =0.0 ;
  double Idd_dVbs =0.0, Idd_dVds =0.0, Idd_dVgs =0.0 ;
  double Fdd =0.0 ;
  double Fdd_dVbs =0.0, Fdd_dVds =0.0, Fdd_dVgs =0.0 ;
  double Eeff =0.0 ;
  double Eeff_dVbs =0.0, Eeff_dVds =0.0, Eeff_dVgs =0.0 ;
  double Rns =0.0 ;
  double Mu = 0.0 ;
  double Mu_dVbs =0.0, Mu_dVds =0.0, Mu_dVgs =0.0 ;
  double Muun =0.0, Muun_dVbs =0.0, Muun_dVds =0.0, Muun_dVgs =0.0 ;
  double Ey = 0e0 ;
  double Ey_dVbs =0.0, Ey_dVds =0.0, Ey_dVgs =0.0 ;
  double Em =0.0 ;
  double Em_dVbs =0.0, Em_dVds =0.0, Em_dVgs =0.0 ;
  double Vmax =0.0 ;
  double Eta =0.0 ;
  double Eta_dVbs =0.0, Eta_dVds =0.0, Eta_dVgs =0.0 ;
  double Eta1 =0.0, Eta1p12 =0.0, Eta1p32 =0.0, Eta1p52 =0.0 ;
  double Zeta12 =0.0, Zeta32 =0.0, Zeta52 =0.0 ;
  double F00 =0.0 ;
  double F00_dVbs =0.0, F00_dVds =0.0, F00_dVgs =0.0 ;
  double F10 =0.0 ;
  double F10_dVbs =0.0, F10_dVds =0.0, F10_dVgs =0.0 ;
  double F30 =0.0 ;
  double F30_dVbs =0.0, F30_dVds =0.0, F30_dVgs =0.0 ;
  double F11 =0.0 ;
  double F11_dVbs =0.0, F11_dVds =0.0, F11_dVgs =0.0 ;
  double Ps0_min =0.0 ;
  double Ps0_min_dVbs =0.0, Ps0_min_dVds =0.0, Ps0_min_dVgs =0.0 ;
  double Acn =0.0, Acd =0.0, Ac1 =0.0, Ac2 =0.0, Ac3 =0.0, Ac4 =0.0, Ac31 =0.0, Ac41 =0.0 ;
  double Acn_dVbs =0.0, Acn_dVds =0.0, Acn_dVgs =0.0 ;
  double Acd_dVbs =0.0, Acd_dVds =0.0, Acd_dVgs =0.0 ;
  double Ac1_dVbs =0.0, Ac1_dVds =0.0, Ac1_dVgs =0.0 ;
  double Ac2_dVbs =0.0, Ac2_dVds =0.0, Ac2_dVgs =0.0 ;
  double Ac3_dVbs =0.0, Ac3_dVds =0.0, Ac3_dVgs =0.0 ;
  double Ac4_dVbs =0.0, Ac4_dVds =0.0, Ac4_dVgs =0.0 ;
  double Ac31_dVbs =0.0, Ac31_dVds =0.0, Ac31_dVgs =0.0 ;
  /* PART-2 (Isub) */
  double Isub =0.0 ;
  double Isub_dVbs =0.0, Isub_dVds =0.0, Isub_dVgs =0.0 ;
  double Isub_dVbse =0.0, Isub_dVdse =0.0, Isub_dVgse =0.0 ;
  double Psislsat=0.0, Psisubsat =0.0 ;
  double Psislsat_dVd=0.0, Psislsat_dVg=0.0, Psislsat_dVb =0.0 ;
  double Psisubsat_dVd=0.0, Psisubsat_dVg=0.0, Psisubsat_dVb =0.0 ;
  /* PART-3 (overlap) */
  double cov_slp =0.0, cov_mag =0.0, covvg =0.0, covvg_dVgs =0.0 ;
  double Lov =0.0 ;
  double Qgos = 0.0, Qgos_dVbs = 0.0, Qgos_dVds = 0.0, Qgos_dVgs = 0.0 ;
  double Qgos_dVbse =0.0, Qgos_dVdse =0.0, Qgos_dVgse =0.0 ;
  double Qgod = 0.0, Qgod_dVbs = 0.0, Qgod_dVds = 0.0, Qgod_dVgs = 0.0 ;
  double Qgod_dVbse =0.0, Qgod_dVdse =0.0, Qgod_dVgse =0.0 ;

  int flg_overgiven =0 ;

  double Qgbo =0.0, Qgbo_dVbs =0.0, Qgbo_dVds =0.0, Qgbo_dVgs =0.0 ;
  double Qgbo_dVbse =0.0, Qgbo_dVdse =0.0, Qgbo_dVgse =0.0 ;
  double Cggo = 0.0 , Cgdo = 0.0 , Cgso = 0.0 , Cgbo = 0.0 , Cgbo_loc=0.0 ;
  /* fringing capacitance */
  double Cf =0.0 ;
  double Qfd =0.0, Qfs =0.0 ;
  /* Cqy */
  double Ec =0.0, Ec_dVbs =0.0, Ec_dVds =0.0, Ec_dVgs =0.0 ;
  double Pslk =0.0, Pslk_dVbs =0.0, Pslk_dVds =0.0, Pslk_dVgs =0.0 ;
  double Qy =0.0 ;
  double Cqyd=0.0, Cqyg=0.0, Cqys=0.0, Cqyb =0.0 ;
  double Qy_dVbs =0.0, Qy_dVds =0.0, Qy_dVgs=0.0 ;
  double Qy_dVbse =0.0, Qy_dVdse=0.0, Qy_dVgse=0.0 ;
  double Qys=0.0, Qys_dVbse =0.0, Qys_dVdse=0.0, Qys_dVgse=0.0 ;
  /* PART-4 (junction diode) */
  double Ibs =0.0, Ibd =0.0, Gbs =0.0, Gbd =0.0, Gbse =0.0, Gbde =0.0 ;
  /* junction capacitance */
  double Qbs =0.0, Qbd =0.0, Capbs =0.0, Capbd =0.0, Capbse =0.0, Capbde =0.0 ;
  double czbd =0.0, czbdsw =0.0, czbdswg =0.0, czbs =0.0, czbssw =0.0, czbsswg =0.0 ;
  double arg =0.0, sarg =0.0 ;
  /* PART-5 (NQS) */
  double tau=0.0, Qi_prev =0.0; 
  double tau_dVgs=0.0, tau_dVds=0.0, tau_dVbs =0.0 ;
  double tau_dVgse=0.0, tau_dVdse=0.0, tau_dVbse =0.0 ;
  double Qi_nqs =0.0 ;
  double Qi_dVbs_nqs=0.0, Qi_dVds_nqs=0.0, Qi_dVgs_nqs =0.0 ;
  double Qi_dVbse_nqs=0.0, Qi_dVdse_nqs=0.0, Qi_dVgse_nqs =0.0 ;
  double taub=0.0, Qb_prev =0.0; 
  double taub_dVgs=0.0, taub_dVds=0.0, taub_dVbs =0.0 ;
  double taub_dVgse=0.0, taub_dVdse=0.0, taub_dVbse =0.0 ;
  double Qb_nqs =0.0 ;
  double Qb_dVbs_nqs=0.0, Qb_dVds_nqs=0.0, Qb_dVgs_nqs =0.0 ;
  double Qb_dVbse_nqs=0.0, Qb_dVdse_nqs=0.0, Qb_dVgse_nqs =0.0 ;
  /* PART-6 (noise) */
  /* 1/f */
  double NFalp =0.0, NFtrp =0.0, Cit =0.0, Nflic =0.0 ;
  /* thermal */
  double Eyd=0.0, Mu_Ave=0.0, Nthrml=0.0, Mud_hoso =0.0 ;
  /* induced gate noise ( Part 0/3 ) */
  double kusai00 =0.0, kusaidd =0.0, kusaiL =0.0, kusai00L =0.0 ;
  int flg_ign = 0 ;
  double sqrtkusaiL =0.0, kusai_ig =0.0, gds0_ign =0.0, gds0_h2 =0.0, GAMMA =0.0, crl_f =0.0 ;
  const double c_sqrt_15 = 3.872983346207417e0 ; /* sqrt(15) */
  const double Cox_small = 1.0e-6 ;
  const double c_16o135 = 1.185185185185185e-1 ; /* 16/135 */
  double Nign0 =0.0, MuModA =0.0, MuModB =0.0, correct_w1 =0.0 ;

  /* Bias iteration accounting Rs/Rd */
  int lp_bs =0 ;
  double Ids_last =0.0 ;
  double vtol_iprv = 2.0e-1 ;
  double vtol_pprv = 1.01e-1 ;
  double Vbsc_dif =0.0, Vdsc_dif =0.0, Vgsc_dif =0.0, sum_vdif =0.0 ;
  double Vbsc_dif2 =0.0, Vdsc_dif2 =0.0, Vgsc_dif2 =0.0, sum_vdif2 =0.0 ;
  double Rs =0.0, Rd =0.0 ;
  double Fbs =0.0, Fds =0.0, Fgs =0.0 ;
  double DJ =0.0, DJI =0.0 ;
  double JI11 =0.0, JI12 =0.0, JI13 =0.0, JI21 =0.0, JI22 =0.0, JI23 =0.0, JI31 =0.0, JI32 =0.0, JI33 =0.0 ;
  double dVbs =0.0, dVds =0.0, dVgs =0.0 ;
  double dV_sum =0.0 ;
  /* temporary vars. */
  double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12 ;
  double TX =0.0, TX_dVbs =0.0, TX_dVds =0.0, TX_dVgs =0.0 ;
  double TY =0.0, TY_dVbs =0.0, TY_dVds =0.0, TY_dVgs =0.0 ;
  double T1_dVb =0.0, T1_dVd =0.0, T1_dVg =0.0 ;
  double T2_dVb =0.0, T2_dVd =0.0, T2_dVg =0.0 ;
  double T3_dVb =0.0, T3_dVd =0.0, T3_dVg =0.0 ;
  double T4_dVb =0.0, T4_dVd =0.0, T4_dVg =0.0 ;
  double T5_dVb =0.0, T5_dVd =0.0, T5_dVg =0.0 ;
  double T6_dVb =0.0, T6_dVd =0.0, T6_dVg =0.0 ;
  double T7_dVb =0.0, T7_dVd =0.0, T7_dVg =0.0 ;
  double T8_dVb =0.0, T8_dVd =0.0, T8_dVg =0.0 ;
  double T9_dVb =0.0, T9_dVd =0.0, T9_dVg =0.0 ;
  double T10_dVb =0.0, T10_dVd =0.0, T10_dVg =0.0 ;

  int   flg_zone = 0 ;
  double Vfbsft = 0.0 , Vfbsft_dVbs =0.0, Vfbsft_dVds =0.0, Vfbsft_dVgs =0.0 ;

  /* Vdseff */
  double Vdseff =0.0, Vdsorg =0.0 ;
  double Vdseff_dVbs =0.0, Vdseff_dVds =0.0, Vdseff_dVgs =0.0 ;

  /* G/S and G/D Overlap Charges: Qovs/Qovd */
  double Qovd = 0.0, Qovd_dVbse = 0.0, Qovd_dVdse = 0.0, Qovd_dVgse = 0.0 ;
  double Qovd_dVbs = 0.0, Qovd_dVds = 0.0, Qovd_dVgs = 0.0 ;
  double Qovs = 0.0, Qovs_dVbse = 0.0, Qovs_dVdse = 0.0, Qovs_dVgse = 0.0 ;
  double Qovs_dVbs = 0.0, Qovs_dVds = 0.0, Qovs_dVgs = 0.0 ;
  int    lcover = 0, flg_ovloops = 0, flg_ovloopd = 0 ;
  int    flg_overs = 0, flg_overd = 0 ;
  double VgpLD =0.0 ;
  double QbdLD = 0.0 , QbdLD_dVbs = 0.0 , QbdLD_dVds = 0.0 , QbdLD_dVgs = 0.0 ;
  double QidLD = 0.0 , QidLD_dVbs = 0.0 , QidLD_dVds = 0.0 , QidLD_dVgs = 0.0 ;
  double QbsLD = 0.0 , QbsLD_dVbs = 0.0 , QbsLD_dVds = 0.0 , QbsLD_dVgs = 0.0 ;
  double QisLD = 0.0 , QisLD_dVbs = 0.0 , QisLD_dVds = 0.0 , QisLD_dVgs = 0.0 ;
  double QbdLD_dVbse = 0.0 , QbdLD_dVdse = 0.0 , QbdLD_dVgse = 0.0 ;
  double QidLD_dVbse = 0.0 , QidLD_dVdse = 0.0 , QidLD_dVgse = 0.0 ;
  double QbsLD_dVbse = 0.0 , QbsLD_dVdse = 0.0 , QbsLD_dVgse = 0.0 ;
  double QisLD_dVbse = 0.0 , QisLD_dVdse = 0.0 , QisLD_dVgse = 0.0 ;
  double QbuLD = 0.0 , QbuLD_dVbs = 0.0 , QbuLD_dVds = 0.0 , QbuLD_dVgs = 0.0 ;
  double QsuLD = 0.0 , QsuLD_dVbs = 0.0 , QsuLD_dVds = 0.0 , QsuLD_dVgs = 0.0 ;
  double QiuLD = 0.0 , QiuLD_dVbs = 0.0 , QiuLD_dVds = 0.0 , QiuLD_dVgs = 0.0 ;
  double Ps0LD = 0.0 , Ps0LD_dVds = 0.0 ;
  double QbuLD_dVxb = 0.0 , QbuLD_dVgb = 0.0 ;
  double QsuLD_dVxb = 0.0 , QsuLD_dVgb = 0.0 ;
  int   flg_ovzone = 0 ;

  /* Vgsz for SCE and PGD */
  double Vbsz2 =0.0, Vbsz2_dVbs =0.0, Vbsz2_dVds =0.0, Vbsz2_dVgs =0.0 ;

  /* Multiplication factor of a MOSFET instance */
  double M = 1.0 ;

  /* Mode flag ( = 0 | 1 )  */
  double ModeNML =0.0, ModeRVS =0.0 ;

  /* For Gate Leak Current Partitioning */
  double GLPART1 ; 
  double GLPART1_dVgs=0.0, GLPART1_dVds=0.0, GLPART1_dVbs =0.0; 
  double GLPART1_dVgse=0.0, GLPART1_dVdse=0.0, GLPART1_dVbse =0.0; 

  /* IBPC */
  double IdsIBPC = 0.0 ;
  double IdsIBPC_dVbs = 0.0 , IdsIBPC_dVds = 0.0 , IdsIBPC_dVgs = 0.0 ;
  double IdsIBPC_dVbse = 0.0 , IdsIBPC_dVdse = 0.0 , IdsIBPC_dVgse = 0.0 ;

  /* Overlap Charge: Qover */
  double Vbsgmt =0.0, Vdsgmt =0.0, Vgsgmt =0.0, Vdbgmt =0.0, Vgbgmt =0.0, Vsbgmt =0.0, Vxbgmt =0.0 ;
  double Vxbgmtcl = 0.0, Vxbgmtcl_dVxbgmt = 0.0 ;

  double Pb2over =0.0 ;

  /* Qover Iterative and Analytical Model */
  const double large_arg = 80 ;
  int lp_ld =0 ;
  double T1_dVxb=0.0, T1_dVgb=0.0, T5_dVxb=0.0, T5_dVgb =0.0 ;
  double Vgb_fb_LD=0.0,  VgpLD_dVgb =0.0 ;
  double VgpLD_shift=0.0, VgpLD_shift_dVxb=0.0, VgpLD_shift_dVgb =0.0 ;
  double TX_dVxb=0.0, TX_dVgb=0.0, TY_dVxb=0.0, TY_dVgb =0.0 ;
  double Ac1_dVxb=0.0, Ac1_dVgb=0.0, Ac2_dVxb=0.0, Ac2_dVgb =0.0 ;
  double Ac3_dVxb=0.0, Ac3_dVgb=0.0, Ac31_dVxb=0.0, Ac31_dVgb =0.0; 
  double Acd_dVxb=0.0, Acd_dVgb=0.0, Acn_dVxb=0.0, Acn_dVgb =0.0;  
  double Ta = 9.3868e-3, Tb = -0.1047839 ;
  double Tc=0.0, Tp =0.0 ;
  double Td=0.0, Td_dVxb=0.0, Td_dVgb =0.0 ;
  double Tq=0.0, Tq_dVxb=0.0, Tq_dVgb =0.0 ;
  double Tu=0.0, Tu_dVxb=0.0, Tu_dVgb =0.0 ;
  double Tv=0.0, Tv_dVxb=0.0, Tv_dVgb =0.0 ;
  double exp_bVbs_dVxb=0.0, exp_bPs0 =0.0 ;
  double cnst1over =0.0; 
  double gamma=0.0, gamma_dVxb =0.0; 
  double Chi_dVxb=0.0, Chi_dVgb =0.0 ;
  double Chi_A=0.0, Chi_A_dVxb=0.0, Chi_A_dVgb =0.0 ;
  double Chi_B=0.0, Chi_B_dVxb=0.0, Chi_B_dVgb =0.0 ;
  double Chi_1=0.0, Chi_1_dVxb=0.0, Chi_1_dVgb =0.0 ;
/*  double psi_B=0.0, arg_B =0.0 ;*/
  double psi=0.0, psi_dVgb=0.0, psi_dVxb =0.0 ;
  double Ps0_iniA_dVxb=0.0, Ps0_iniA_dVgb =0.0 ;
/*  double Ps0_iniB_dVxb=0.0 , Ps0_iniB_dVgb =0.0 ;*/
  double Psa_dVxb=0.0, Psa_dVgb=0.0, Ps0LD_dVxb=0.0, Ps0LD_dVgb =0.0 ;
  double /*fs02_dVxb=0.0,*/ fs02_dVgb =0.0 ;

  /* SCE LOOP */
  double A =0.0, A_dVgs=0.0, A_dVds=0.0, A_dVbs =0.0 ;
  int    NNN =0 ;
/*  double PS0_SCE=0 ,  PS0_SCE_dVgs = 0 ,  PS0_SCE_dVds = 0 ,  PS0_SCE_dVbs = 0 ;*/
  double PS0Z_SCE=0 , PS0Z_SCE_dVgs = 0 , PS0Z_SCE_dVds = 0 , PS0Z_SCE_dVbs = 0 ;
   /* double arg0 = 0.01 , arg1 = 0.04 ; */
   double arg0 = 0.01 ;
   double arg2 = here->HSM2_2qnsub_esi * 1.0e-4 ;
   int MAX_LOOP_SCE =0 ;

  int codqb = 0 ;
  int corecip = model->HSM2_corecip ;


  /* modify Qy in accumulation region */
  double Aclm_eff=0.0, Aclm_eff_dVds=0.0, Aclm_eff_dVgs=0.0, Aclm_eff_dVbs =0.0 ;

  double Idd1 =0.0, Idd1_dVbs =0.0,  Idd1_dVgs =0.0,  Idd1_dVds =0.0 ;

  double tcjbs=0.0, tcjbssw=0.0, tcjbsswg=0.0,
         tcjbd=0.0, tcjbdsw=0.0, tcjbdswg=0.0 ;
  double TTEMP =0.0 ;

  double PS0_SCE_tol = 4.0e-7 ;
  double PS0_SCE_deriv_tol = 1.0e-8 ;
  double Ps0_ini_dVds  =0.0,  Ps0_ini_dVgs =0.0, Ps0_ini_dVbs  =0.0 ;
  double Ps0_iniA_dVds =0.0, Ps0_iniA_dVgs =0.0, Ps0_iniA_dVbs =0.0 ;
  double Ps0_iniB_dVds =0.0, Ps0_iniB_dVgs =0.0, Ps0_iniB_dVbs =0.0 ;

  double A_dPS0Z = 0.0,      dqb_dPS0Z = 0.0,    dVth_dPS0Z = 0.0,   dVth0_dPS0Z = 0.0,
         dVthLP_dPS0Z = 0.0, dVthSC_dPS0Z = 0.0, Qbmm_dPS0Z = 0.0,   Vfbsft_dPS0Z = 0.0,
         Vgp_dPS0Z = 0.0,    Vgpz_dPS0Z = 0.0,                       Vthp_dPS0Z = 0.0,
         Vth0_dPS0Z = 0.0 ;
  double T1_dPS0Z=0.0,           T3_dPS0Z=0.0,           T4_dPS0Z=0.0,           T5_dPS0Z=0.0,
         T6_dPS0Z=0.0,           T7_dPS0Z=0.0,           T8_dPS0Z=0.0,           T9_dPS0Z=0.0,
         T10_dPS0Z=0.0,          TX_dPS0Z =0.0 ;
  double Ac1_dPS0Z=0.0,          Ac2_dPS0Z=0.0,          Ac3_dPS0Z=0.0,          Ac31_dPS0Z=0.0,
         Acd_dPS0Z=0.0,          Acn_dPS0Z=0.0,          Chi_dPS0Z=0.0,          Psa_dPS0Z =0.0 ;
  double Fs0_dPS0Z=0.0,          Fsl_dPS0Z=0.0,          Ps0_dPS0Z=0.0,          Psl_dPS0Z=0.0,
         Pds_dPS0Z=0.0,          Pzadd_dPS0Z=0.0,        Ps0z_dPS0Z =0.0 ;
  double G=0.0,                  delta_PS0Z_SCE=0.0,     delta_PS0Z_SCE_dVds=0.0,
         delta_PS0Z_SCE_dVgs=0.0, delta_PS0Z_SCE_dVbs =0.0 ;


  double Vgs_min =0.0 ;

  /*================ Start of executable code.=================*/

  if (here->HSM2_mode == HiSIM_NORMAL_MODE) {
    ModeNML = 1.0 ;
    ModeRVS = 0.0 ;
  } else {
    ModeNML = 0.0 ;
    ModeRVS = 1.0 ;
  }
  
  T1 = vbs + vds + vgs + vbd_jct + vbs_jct ;
  if ( ! finite (T1) ) {
    fprintf (stderr ,
       "*** warning(HiSIM): Unacceptable Bias(es).\n" ) ;
    fprintf (stderr , "----- bias information (HiSIM)\n" ) ;
    fprintf (stderr , "name: %s\n" , here->HSM2name ) ;
    fprintf (stderr , "states: %d\n" , here->HSM2states ) ;
    fprintf (stderr , "vds= %12.5e vgs=%12.5e vbs=%12.5e\n"
            , vds , vgs , vbs ) ;
    fprintf (stderr , "vbs_jct= %12.5e vbd_jct= %12.5e\n"
            , vbs_jct , vbd_jct ) ;
    fprintf (stderr , "vd= %12.5e vg= %12.5e vb= %12.5e vs= %12.5e\n" 
            , *( ckt->CKTrhsOld + here->HSM2dNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2gNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2bNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2sNodePrime )  ) ;
    if ( here->HSM2_called >= 1 ) {
      fprintf (stderr , "vdsc_prv= %12.5e vgsc_prv=%12.5e vbsc_prv=%12.5e\n"
              , here->HSM2_vdsc_prv , here->HSM2_vgsc_prv
              , here->HSM2_vbsc_prv ) ;
    }
    fprintf (stderr , "----- bias information (end)\n" ) ;
    return ( HiSIM_ERROR ) ;
  }

  flg_info = model->HSM2_info ;
  flg_nqs = model->HSM2_conqs ;
  
  /*-----------------------------------------------------------*
   * Start of the routine. (label)
   *-----------------*/
/*start_of_routine:*/

  TTEMP = ckt->CKTtemp ;
  if ( here->HSM2_dtemp_Given ) { TTEMP = TTEMP + here->HSM2_dtemp ; }

  beta = here->HSM2_beta ;

  /* Inverse of the thermal voltage */
  beta_inv = here->HSM2_beta_inv ;
  beta2 = here->HSM2_beta2 ;

  /* Bandgap */
  Egp12 = here->HSM2_egp12 ;
  Egp32 = here->HSM2_egp32 ;

  /* Metallurgical channel geometry */
  Leff = here->HSM2_leff ;
  Weff = here->HSM2_weff ;

  /* Flat band voltage */
  Vfb = pParam->HSM2_vfbc ;

  /* Surface impurity profile */
  q_Nsub = here->HSM2_qnsub ;
  
  /* Velocity Temperature Dependence */
  Vmax = here->HSM2_vmax ;
   
  /* 2 phi_B */
  Pb2 = here->HSM2_pb2 ;
  Pb20 = here->HSM2_pb20 ; 
  Pb2c = here->HSM2_pb2c ;

  /* Coefficient of the F function for bulk charge */
  cnst0 = here->HSM2_cnst0 ;

  /* cnst1: n_{p0} / p_{p0} */
  cnst1 = here->HSM2_cnst1 ;

  /* c_eox: Permitivity in ox  */
  c_eox = here->HSM2_cecox ;

  /* Tox and Cox without QME */
   Tox0 = model->HSM2_tox ;
   Cox0 = c_eox / Tox0 ;
   Cox0_inv = 1.0 / Cox0 ;

  /* for calculation of Ps0_min */
  Vgs_min = model->HSM2_type * model->HSM2_Vgsmin ;

  /*-----------------------------------------------------------*
   * Exchange bias conditions according to MOS type.
   * - Vxse are external biases for HiSIM. ( type=NMOS , Vds >= 0
   *   are assumed.) 
   *-----------------*/
  Vbse = vbs ;
  Vdse = vds ;
  Vgse = vgs ;

  /*---------------------------------------------------*
   * Clamp too large biases. 
   * -note: Quantities are extrapolated in PART-5.
   *-----------------*/

  if ( Pb2 - model->HSM2_vzadd0 < Vbs_max ) {
    Vbs_max = Pb2 - model->HSM2_vzadd0 ;
  }
  if ( Pb20 - model->HSM2_vzadd0 < Vbs_max ) {
    Vbs_max = Pb20 - model->HSM2_vzadd0 ;
  }
  if ( Pb2c - model->HSM2_vzadd0 < Vbs_max ) {
    Vbs_max = Pb2c - model->HSM2_vzadd0 ;
  }

  if ( Vbs_bnd > Vbs_max * 0.5 ) {
    Vbs_bnd = 0.5 * Vbs_max ;
  }

  if ( Vbse > Vbs_bnd ) {
    flg_vbsc = 1 ;
    T1 = Vbse - Vbs_bnd ;
    T2 = Vbs_max - Vbs_bnd ;
    Fn_SUPoly4( TY , T1 , T2 , Vbsc_dVbse ) ; 
    Vbsc = Vbs_bnd + TY ;
  }  else if ( Vbse < Vbs_min ) {
    flg_vbsc = -1 ;
    Vbsc = Vbs_min ;
    Vbsc_dVbse = 1.0 ;
  }  else {
    flg_vbsc =  0 ;
    Vbsc = Vbse ;
    Vbsc_dVbse = 1.0 ;
  }

  Vdsc = Vdse ;
  Vgsc = Vgse ;

  if (here->HSM2_rs > 0.0 || here->HSM2_rd > 0.0) {
    if ( model->HSM2_corsrd == 1 ) flg_rsrd  = 1 ;
    if ( model->HSM2_corsrd == 2 ) flg_rsrd  = 2 ;
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * PART-1: Basic device characteristics. 
   *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
   * Prepare for potential initial guesses using previous values
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  flg_pprv = 0 ;

  if ( here->HSM2_called >= 1 ) {

    Vbsc_dif = Vbsc - here->HSM2_vbsc_prv ;
    Vdsc_dif = Vdsc - here->HSM2_vdsc_prv ;
    Vgsc_dif = Vgsc - here->HSM2_vgsc_prv ;

    sum_vdif  = fabs( Vbsc_dif ) + fabs( Vdsc_dif ) 
              + fabs( Vgsc_dif ) ;

    if ( model->HSM2_copprv >= 1 && sum_vdif <= vtol_pprv  &&
         here->HSM2_mode * here->HSM2_mode_prv > 0 ) { flg_pprv = 1 ;}

    if ( here->HSM2_called >= 2 && flg_pprv == 1 ) {
      Vbsc_dif2 = here->HSM2_vbsc_prv - here->HSM2_vbsc_prv2 ;
      Vdsc_dif2 = here->HSM2_vdsc_prv - here->HSM2_vdsc_prv2 ;
      Vgsc_dif2 = here->HSM2_vgsc_prv - here->HSM2_vgsc_prv2 ;
      sum_vdif2  = fabs( Vbsc_dif2 ) + fabs( Vdsc_dif2 ) 
                + fabs( Vgsc_dif2 ) ;
      if ( epsm10 < sum_vdif2 && sum_vdif2 <= vtol_pprv &&
           here->HSM2_mode_prv * here->HSM2_mode_prv2 > 0 ) { flg_pprv = 2 ; }
    }
  } else {
    Vbsc_dif = 0.0 ;
    Vdsc_dif = 0.0 ;
    Vgsc_dif = 0.0 ;
    sum_vdif = 0.0 ;
    Vbsc_dif2 = 0.0 ;
    Vdsc_dif2 = 0.0 ;
    Vgsc_dif2 = 0.0 ;
    sum_vdif2 = 0.0 ;
    flg_iprv = 0 ;
    flg_pprv = 0 ;
  }

  dVbs = Vbsc_dif ;
  dVds = Vdsc_dif ;
  dVgs = Vgsc_dif ;

  if ( flg_pprv >= 1 ) {
    Ps0 = here->HSM2_ps0_prv ;
    Ps0_dVbs = here->HSM2_ps0_dvbs_prv ;
    Ps0_dVds = here->HSM2_ps0_dvds_prv ;
    Ps0_dVgs = here->HSM2_ps0_dvgs_prv ;
  
    Pds = here->HSM2_pds_prv ;
    Pds_dVbs = here->HSM2_pds_dvbs_prv ;
    Pds_dVds = here->HSM2_pds_dvds_prv ;
    Pds_dVgs = here->HSM2_pds_dvgs_prv ;
  }


  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
   * Bias loop: iteration to solve the system of equations of 
   *            the small circuit taking into account Rs and Rd.
   * - Vxs are internal (or effective) biases.
   * - Equations:
   *     Vbs = Vbsc - Rs * Ids
   *     Vds = Vdsc - ( Rs + Rd ) * Ids
   *     Vgs = Vgsc - Rs * Ids
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  if ( flg_rsrd == 1 ) {

    if (here->HSM2_mode == HiSIM_NORMAL_MODE) {
      Rs  = here->HSM2_rs ;
      Rd  = here->HSM2_rd ;
    } else {
      Rs  = here->HSM2_rd ;
      Rd  = here->HSM2_rs ;
    }
    
    if ( here->HSM2_called >= 1 ) {
      
      if ( model->HSM2_coiprv >= 1 && 
           0.0 < sum_vdif && sum_vdif <= vtol_iprv ) { flg_iprv = 1 ;}
    }

    /*-----------------------------------------------------------*
     * Initial guesses using the previous values.
     *-----------------*/
    if ( flg_iprv == 1 ) {
      here->HSM2_ids_dvbs_prv = Fn_Max( 0.0 , here->HSM2_ids_dvbs_prv ) ;
      here->HSM2_ids_dvds_prv = Fn_Max( 0.0 , here->HSM2_ids_dvds_prv ) ;
      here->HSM2_ids_dvgs_prv = Fn_Max( 0.0 , here->HSM2_ids_dvgs_prv ) ;
    
      dVbs = Vbsc_dif * ( 1.0 - 1.0 / ( 1.0 + Rs * here->HSM2_ids_dvbs_prv ) ) ;
      dVds = Vdsc_dif * ( 1.0 - 1.0 / ( 1.0 + ( Rs + Rd ) * here->HSM2_ids_dvds_prv ) ) ;
      dVgs = Vgsc_dif * ( 1.0 - 1.0 / ( 1.0 + Rs * here->HSM2_ids_dvgs_prv ) ) ;

      Ids = here->HSM2_ids_prv 
          + here->HSM2_ids_dvbs_prv * dVbs 
          + here->HSM2_ids_dvds_prv * dVds 
          + here->HSM2_ids_dvgs_prv * dVgs ;
      
      T1  = ( Ids - here->HSM2_ids_prv ) ;
      T2  = fabs( T1 ) ;
      if ( Ids_maxvar * here->HSM2_ids_prv < T2 ) {
            Ids = here->HSM2_ids_prv * ( 1.0 + Fn_Sgn( T1 ) * Ids_maxvar ) ;
      }
      if ( Ids < 0 ) Ids = 0.0 ;
      
    } else {
      Ids = 0.0 ;
      if ( flg_pprv >= 1 ) {
        dVbs = Vbsc_dif ;
        dVds = Vdsc_dif ;
        dVgs = Vgsc_dif ;
      }
    } /* end of flg_iprv if-blocks */
    
    Vbs = Vbsc - Ids * Rs ;
    Vds = Vdsc - Ids * ( Rs + Rd ) ;
    if ( Vds * Vdsc <= 0.0 ) { Vds = 0.0 ; } 
    
    Vgs = Vgsc - Ids * Rs ;

  } else {
    lp_bs_max = 1 ;
    Ids = 0.0 ;
    Vbs = Vbsc ;
    Vds = Vdsc ;
    Vgs = Vgsc ;
  } /* end of flg_rsrd if-blocks */

  /*-----------------------------------------------------------*
   * start of the loop.
   *-----------------*/
  for ( lp_bs = 1 ; lp_bs <= lp_bs_max ; lp_bs ++ ) {

    Ids_last = Ids ;
    /* Initialization of counters is needed for restart. */
    lp_s0   = 0 ;
    lp_sl   = 0 ;

    /*-----------------------------------------------------------*
     * Vxsz: Modified bias introduced to realize symmetry at Vds=0.
     *-----------------*/

    T1 = Vbsc_dVbse * Vds / 2 ;
    Fn_SymAdd(  Vzadd , T1 , model->HSM2_vzadd0 , T2 ) ;
    T2 *= Vbsc_dVbse / 2 ;
    Vzadd_dVds = T2 ;

    if ( Vzadd < ps_conv ) {
      Vzadd = ps_conv ;
      Vzadd_dVds = 0.0 ;
    }

    Vbsz = Vbs + Vzadd ;
    Vbsz_dVbs = 1.0 ;
    Vbsz_dVds = Vzadd_dVds ;
  
    Vdsz = Vds + 2.0 * Vzadd ;
    Vdsz_dVds = 1.0 + 2.0 * Vzadd_dVds ;

    Vgsz = Vgs + Vzadd ;
    Vgsz_dVgs = 1.0 ;
    Vgsz_dVds = Vzadd_dVds ;

    /*---------------------------------------------------*
     * Factor of modification for symmetry.
     *-----------------*/

    T1 = here->HSM2_qnsub_esi * Cox0_inv * Cox0_inv ;
    T2 = Vgs - Vfb ;
    T3 = 1 + 2.0 / T1 * ( T2 - beta_inv - Vbs ) ;
    Fn_SZ( T4 , T3 , 1e-3 , T5 ) ;
    TX = sqrt( T4 ) ;
    Pslsat = T2 + T1 * ( 1.0 - TX ) ; 
    VdsatS = Pslsat - Pb2 ;
    Fn_SL( VdsatS , VdsatS , 0.1 , 5e-2 , T6 ) ;

    VdsatS_dVbs = ( TX ? (T6 * T5 / TX ) : 0.0 ) ;
    VdsatS_dVds = 0.0 ;
    VdsatS_dVgs = ( TX ? (T6 * ( 1.0 - T5 / TX )) : 0.0 ) ;


    T1 = Vds / VdsatS ;
    Fn_SUPoly4( TX , T1 , 1.0 , T0 ) ; 
    FMDVDS = TX * TX ;
    T2 = 2 * TX * T0 ;
    T3 = T2 / ( VdsatS * VdsatS ) ;
    FMDVDS_dVbs = T3 * ( - Vds * VdsatS_dVbs ) ;
    FMDVDS_dVds = T3 * ( 1.0 * VdsatS - Vds * VdsatS_dVds ) ;
    FMDVDS_dVgs = T3 * ( - Vds * VdsatS_dVgs ) ;

    /*-----------------------------------------------------------*
     * Quantum effect
     *-----------------*/
    if ( model->HSM2_flg_qme == 0 ) {
      flg_qme = 0 ;
    } else {
      flg_qme = 1 ;
    }


      T1 = here->HSM2_2qnsub_esi ;
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
      } else {

	T5     = Vgs  - Vbs - Vthq + model->HSM2_qme2 ;
	T5_dVb = -1.0 - Vthq_dVb ;
	T5_dVd =      - Vthq_dVd ;
	T5_dVg =  1.0 ;
  	Fn_SZ( T2 , T5 , qme_dlt, T3) ;
  	T2 = T2 + small ;
  	T2_dVb = T3 * T5_dVb ; 
  	T2_dVd = T3 * T5_dVd ; 
  	T2_dVg = T3 * T5_dVg ; 

	T3 = 1.0 /  T2 ;
	T7 = -1.0 / ( T2 * T2 ) ;
  	T3_dVb = T7 * T2_dVb ; 
  	T3_dVd = T7 * T2_dVd ; 
  	T3_dVg = T7 * T2_dVg ; 

        T4 = 2.0 * fabs(Vthq) ;
        T6 = Vfb - Vthq + model->HSM2_qme2 ;
        if(T6 > T4) { T4 = T6; }

  	Fn_SU( T2 , T3 , 1.0 / T4  , qme_dlt, T6 ) ; 
  	T2_dVb = T6 * T3_dVb ;
  	T2_dVd = T6 * T3_dVd ;
  	T2_dVg = T6 * T3_dVg ;

	dTox = model->HSM2_qme1 * T2 + model->HSM2_qme3 ;
        T7   = model->HSM2_qme1 ;
	dTox_dVb = T7 * T2_dVb ;
	dTox_dVd = T7 * T2_dVd ;
	dTox_dVg = T7 * T2_dVg ;

	if ( dTox * 1.0e12 < Tox0 ) {
	  dTox = 0.0 ;
	  dTox_dVb = 0.0 ;
	  dTox_dVd = 0.0 ;
	  dTox_dVg = 0.0 ;
	  flg_qme = 0 ;
	}

	Tox = Tox0 + dTox ;
	Tox_dVb = dTox_dVb ;
	Tox_dVd = dTox_dVd ;
	Tox_dVg = dTox_dVg ;

	Cox = c_eox / Tox ;
	T1  = - c_eox / ( Tox * Tox ) ;
	Cox_dVb = T1 * Tox_dVb ;
	Cox_dVd = T1 * Tox_dVd ;
	Cox_dVg = T1 * Tox_dVg ;
  
	Cox_inv  = Tox / c_eox ; 
	T1  = 1.0 / c_eox ;
	Cox_inv_dVb = T1 * Tox_dVb ;
	Cox_inv_dVd = T1 * Tox_dVd ;
	Cox_inv_dVg = T1 * Tox_dVg ;
  
	T0 = cnst0 * cnst0 * Cox_inv ;
	cnstCoxi = T0 * Cox_inv ;
	T1 = 2.0 * T0 ;
	cnstCoxi_dVb = T1 * Cox_inv_dVb ;
	cnstCoxi_dVd = T1 * Cox_inv_dVd ;
	cnstCoxi_dVg = T1 * Cox_inv_dVg ;
      }


    fac1 = cnst0 * Cox_inv ;
    fac1_dVds = cnst0 * Cox_inv_dVd ;
    fac1_dVgs = cnst0 * Cox_inv_dVg ;
    fac1_dVbs = cnst0 * Cox_inv_dVb ;
    fac1p2 = fac1 * fac1 ;

    /* Ps0_min: approx. solution of Poisson equation at Vgs_min */
    /*          ( easy to improve, if necessary  )              */
    Ps0_min = 2.0 * beta_inv * log(-Vgs_min/fac1) ;
    Ps0_min_dVds = -2.0 * beta_inv * fac1_dVds / fac1 ;
    Ps0_min_dVgs = -2.0 * beta_inv * fac1_dVgs / fac1 ;
    Ps0_min_dVbs = -2.0 * beta_inv * fac1_dVbs / fac1 ;


    /*---------------------------------------------------*
     * Vbsz2 : Vbs for dVth
     *-----------------*/
    Vbsz2 = Vbsz ;
    Vbsz2_dVbs = Vbsz_dVbs ;
    Vbsz2_dVds = Vbsz_dVds ;
    Vbsz2_dVgs = 0.0 ;

    if ( corecip ) {

      /* ************************** */
      /* Initial value for PS0Z_SCE */
      /* ************************** */

      T1 = dP_max + dP_max ;

      if ( flg_pprv >= 1 ) {
        /* -------------------------- *
         * Extrapolate previous value *
         * -------------------------- */

        T1  =   here->HSM2_PS0Z_SCE_dvbs_prv * dVbs
              + here->HSM2_PS0Z_SCE_dvds_prv * dVds
              + here->HSM2_PS0Z_SCE_dvgs_prv * dVgs ;

        if ( fabs(T1) <= dP_max ) {
          Ps0_ini = here->HSM2_PS0Z_SCE_prv + T1 ; /* take extrapolated value */
          Ps0_ini_dVds = here->HSM2_PS0Z_SCE_dvds_prv ;
          Ps0_ini_dVgs = here->HSM2_PS0Z_SCE_dvgs_prv ;
          Ps0_ini_dVbs = here->HSM2_PS0Z_SCE_dvbs_prv ;
        }
      } /* end of (flg_pprv >=1) if-block */

      if ( fabs(T1) > dP_max) {
	/* ------------------------------------- *
	 * Analytical initial value for PS0Z_SCE *
         * ------------------------------------- */


	T1 = here->HSM2_2qnsub_esi ;
	T2 = sqrt( T1 * ( Pb20 - Vbsz ) ) ;
	Vthq = Pb20 + Vfb + T2 * Cox0_inv ;
	Vth = Vthq ;

	TX     = 4.0e0 * ( beta * ( Vgs - Vbs ) - 1.0e0 ) / ( fac1p2 * beta2 ) ;
        TX_dVds = - 2.0 * TX / fac1 * fac1_dVds ;
        TX_dVgs = - 2.0 * TX / fac1 * fac1_dVgs + 4.0 * beta / ( fac1p2 * beta2 ) ;
        TX_dVbs = - 2.0 * TX / fac1 * fac1_dVbs - 4.0 * beta / ( fac1p2 * beta2 ) ;
        TX    += 1.0 ; 
        if ( TX > epsm10 ) {
          T3     = sqrt( TX ) ;
          T3_dVd = 0.5 * TX_dVds / T3 ;
          T3_dVg = 0.5 * TX_dVgs / T3 ;
          T3_dVb = 0.5 * TX_dVbs / T3 ;
        } else {
          T3 = sqrt( epsm10 ) ;
          T3_dVd = T3_dVg = T3_dVb = 0.0 ;
        }
	Ps0_iniA      = Vgs + fac1p2 * beta * 0.5 * ( 1.0e0 - T3 ) ;
        Ps0_iniA_dVds =       fac1 * beta * ( 1.0 - T3 ) * fac1_dVds - fac1p2 * beta * 0.5 * T3_dVd ;
        Ps0_iniA_dVgs = 1.0 + fac1 * beta * ( 1.0 - T3 ) * fac1_dVgs - fac1p2 * beta * 0.5 * T3_dVg ;
        Ps0_iniA_dVbs =       fac1 * beta * ( 1.0 - T3 ) * fac1_dVbs - fac1p2 * beta * 0.5 * T3_dVb ;

	Chi = beta * ( Ps0_iniA - Vbs ) ;
 
	if ( Chi < znbd3 ) { 
	  /*-----------------------------------*
	   * zone-D1/D2
	   * - Ps0_ini is the analytical solution of Qs=Qb0 with
	   *   Qb0 being approximated to 3-degree polynomial.
	   *-----------------*/
	  TY = beta * ( Vgs - Vbs ) ;
          TY_dVds = 0.0 ;
          TY_dVgs = beta ;
          TY_dVbs = - beta ;
	  T1 = 1.0e0 / ( cn_nc3 * beta * fac1 ) ;
          T1_dVd = - T1 / fac1 * fac1_dVds ;
          T1_dVg = - T1 / fac1 * fac1_dVgs ;
          T1_dVb = - T1 / fac1 * fac1_dVbs ;
	  T2 = 81.0 + 3.0 * T1 ;
          T2_dVd = 3.0 * T1_dVd ;
          T2_dVg = 3.0 * T1_dVg ;
          T2_dVb = 3.0 * T1_dVb ;
	  T3 = -2916.0 - 81.0 * T1 + 27.0 * T1 * TY ;
          T3_dVd = ( - 81.0 + 27.0 * TY ) * T1_dVd + 27.0 * T1 * TY_dVds ;
          T3_dVg = ( - 81.0 + 27.0 * TY ) * T1_dVg + 27.0 * T1 * TY_dVgs ;
          T3_dVb = ( - 81.0 + 27.0 * TY ) * T1_dVb + 27.0 * T1 * TY_dVbs ;
          T4     = T3 ;
          T4_dVd = T3_dVd ;
          T4_dVg = T3_dVg ;
          T4_dVb = T3_dVb ;
          T6     = sqrt( 4 * T2 * T2 * T2 + T4 * T4 ) ;
          T6_dVd = ( 6.0 * T2 * T2 * T2_dVd + T4 * T4_dVd ) / T6 ;
          T6_dVg = ( 6.0 * T2 * T2 * T2_dVg + T4 * T4_dVg ) / T6 ;
          T6_dVb = ( 6.0 * T2 * T2 * T2_dVb + T4 * T4_dVb ) / T6 ;
	  T5     = Fn_Pow( T3 + T6 , C_1o3 ) ;
          T5_dVd = ( T3_dVd + T6_dVd ) / ( 3.0 * T5 * T5 ) ;
          T5_dVg = ( T3_dVg + T6_dVg ) / ( 3.0 * T5 * T5 ) ;
          T5_dVb = ( T3_dVb + T6_dVb ) / ( 3.0 * T5 * T5 ) ;

	  TX = 3.0 - ( C_2p_1o3 * T2 ) / ( 3.0 * T5 )
	    + 1 / ( 3.0 * C_2p_1o3 ) * T5 ;
          TX_dVds = - C_2p_1o3 / (3.0 * T5) * T2_dVd + ( C_2p_1o3 * T2 / (3.0 * T5 * T5) + 1.0 / (3.0 * C_2p_1o3) ) * T5_dVd ;
          TX_dVgs = - C_2p_1o3 / (3.0 * T5) * T2_dVg + ( C_2p_1o3 * T2 / (3.0 * T5 * T5) + 1.0 / (3.0 * C_2p_1o3) ) * T5_dVg ;
          TX_dVbs = - C_2p_1o3 / (3.0 * T5) * T2_dVb + ( C_2p_1o3 * T2 / (3.0 * T5 * T5) + 1.0 / (3.0 * C_2p_1o3) ) * T5_dVb ;
        
	  Ps0_iniA = TX * beta_inv + Vbs ;
          Ps0_iniA_dVds = TX_dVds * beta_inv ;
          Ps0_iniA_dVgs = TX_dVgs * beta_inv ;
          Ps0_iniA_dVbs = TX_dVbs * beta_inv + 1.0 ;
	  Ps0_ini = Ps0_iniA ;
          Ps0_ini_dVds  = Ps0_iniA_dVds ;
          Ps0_ini_dVgs  = Ps0_iniA_dVgs ;
          Ps0_ini_dVbs  = Ps0_iniA_dVbs ;
        
	} else if ( Vgs <= Vth ) { 
	  /*-----------------------------------*
	   * Weak inversion zone.  
	   *-----------------*/
	  Ps0_ini = Ps0_iniA ;
          Ps0_ini_dVds = Ps0_iniA_dVds ;
          Ps0_ini_dVgs = Ps0_iniA_dVgs ;
          Ps0_ini_dVbs = Ps0_iniA_dVbs ;
        
	} else { 
	  /*-----------------------------------*
	   * Strong inversion zone.  
	   * - Ps0_iniB : upper bound.
	   *-----------------*/
	  T1 = 1.0 / cnst1 / cnstCoxi ;
          T1_dVd = - T1 / cnstCoxi * cnstCoxi_dVd ;
          T1_dVg = - T1 / cnstCoxi * cnstCoxi_dVg ;
          T1_dVb = - T1 / cnstCoxi * cnstCoxi_dVb ;
	  T0 = Vgs - Vfb ;
	  T2 = T1 * T0 * T0 ;
          T2_dVd = T1_dVd * T0 * T0 ;
          T2_dVg = T1_dVg * T0 * T0 + 2.0 * T1 * T0 ;
          T2_dVb = T1_dVb * T0 * T0 ;
	  T3 = beta + 2.0 / T0 ;
          T3_dVg = -2.0 / (T0 * T0) ;  
        
	  Ps0_iniB = log( T2 + small ) / T3 ;
          Ps0_iniB_dVds = T2_dVd / (T2 * T3) ;
          Ps0_iniB_dVgs = T2_dVg / (T2 * T3) - T3_dVg * Ps0_iniB / T3 ;
          Ps0_iniB_dVbs = T2_dVb / (T2 * T3) ;
        
	  Fn_SU2( Ps0_ini , Ps0_iniA, Ps0_iniB, c_ps0ini_2, T1,T2) ;
          Ps0_ini_dVds = Ps0_iniA_dVds * T1 + Ps0_iniB_dVds * T2 ;
          Ps0_ini_dVgs = Ps0_iniA_dVgs * T1 + Ps0_iniB_dVgs * T2 ;
          Ps0_ini_dVbs = Ps0_iniA_dVbs * T1 + Ps0_iniB_dVbs * T2 ;
	} 
      } /* end of initial value calulation */
      
      /**************************/

	/* initial value for SCE LOOP */
/*	PS0_SCE = Ps0_ini ;   
      PS0_SCE_dVds  = Ps0_ini_dVds ;
      PS0_SCE_dVgs  = Ps0_ini_dVgs ;
      PS0_SCE_dVbs  = Ps0_ini_dVbs ; 
*/
	PS0Z_SCE = Ps0_ini ;
      PS0Z_SCE_dVds = Ps0_ini_dVds ;
      PS0Z_SCE_dVgs = Ps0_ini_dVgs ;
      PS0Z_SCE_dVbs = Ps0_ini_dVbs ;
    } /* end of corecip=1 case (initial value calculation) */

	MAX_LOOP_SCE = 5 ;
	NNN = 0 ;

    /* ************************************************************************* */
   
    START_OF_SCE_LOOP : /* outer loop of multi level Newton framework            */

    /* ************************************************************************* */
   
                        /* for multi level Newton method we need the derivatives */
                        /* with respect to PS0Z_SCE                              */
                        /* naming convention: ..._dPS0Z means d.../dPS0Z_SCE     */
	

  if ( flg_qme == 1 ) {
    /*---------------------------------------------------*
     * Vthp : Vth with pocket.
     *-----------------*/

    if( corecip ){

      T1 = here->HSM2_2qnsub_esi ;
      Qb0 = sqrt ( T1 ) ;
      Qb0_dVb = 0.0 ;
      Qb0_dVd = 0.0 ;
      Qb0_dVg = 0.0 ;

      Vthp     = PS0Z_SCE + Vfb + Qb0     * Cox_inv                     + here->HSM2_ptovr ;
      Vthp_dVb = PS0Z_SCE_dVbs  + Qb0_dVb * Cox_inv + Qb0 * Cox_inv_dVb ;
      Vthp_dVd = PS0Z_SCE_dVds  + Qb0_dVd * Cox_inv + Qb0 * Cox_inv_dVd ;
      Vthp_dVg = PS0Z_SCE_dVgs  + Qb0_dVg * Cox_inv + Qb0 * Cox_inv_dVg ;
      Vthp_dPS0Z = 1.0 ;

    }else{  /* original */
      T1 = here->HSM2_2qnsub_esi ;
      Qb0 = sqrt (T1 * (Pb20 - Vbsz2)) ;
      T2 = 0.5 * T1 / Qb0 ;
      Qb0_dVb = T2 * (- Vbsz2_dVbs) ;
      Qb0_dVd = T2 * (- Vbsz2_dVds) ;
      Qb0_dVg = T2 * (- Vbsz2_dVgs) ;

      Vthp = Pb20 + Vfb + Qb0 * Cox_inv + here->HSM2_ptovr;
      Vthp_dVb = Qb0_dVb * Cox_inv + Qb0 * Cox_inv_dVb ;
      Vthp_dVd = Qb0_dVd * Cox_inv + Qb0 * Cox_inv_dVd ;
      Vthp_dVg = Qb0_dVg * Cox_inv + Qb0 * Cox_inv_dVg ;
    }


    Pb20b = Pb20 ;
    Pb20b_dVb = 0.0 ;
    Pb20b_dVd = 0.0 ;
    Pb20b_dVg = 0.0 ;
    
    T0 = 0.95 ;
    T1 = T0 * Pb20b - Vbsz2 - 1.0e-3 ;
    T1_dVb = T0 * Pb20b_dVb - Vbsz2_dVbs ;
    T1_dVd = T0 * Pb20b_dVd - Vbsz2_dVds ;
    T1_dVg = T0 * Pb20b_dVg - Vbsz2_dVgs ;
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
    Pbsum = Pb20b - T3 ;
    Pbsum_dVb = Pb20b_dVb - T3_dVb ;
    Pbsum_dVd = Pb20b_dVd - T3_dVd ;
    Pbsum_dVg = Pb20b_dVg - T3_dVg ;

    sqrt_Pbsum = sqrt( Pbsum ) ;

    /*-------------------------------------------*
     * dVthLP : Short-channel effect induced by pocket.
     * - Vth0 : Vth without pocket.
     *-----------------*/
    if ( model->HSM2_lp != 0.0 ) {

      if( corecip ){

	T1 = here->HSM2_2qnsub_esi ;
	T2 = model->HSM2_bs2 - Vbsz2 ;
	T3 = T2 + small ;
	T4 = sqrt (T3 * T3 + 4.0 * vth_dlt) ;
	T5 = 0.5 * (T3 + T4) ;
	T6 = 0.5 * (1.0 + T3 / T4) ;
	T5_dVb = - Vbsz2_dVbs * T6 ;
	T5_dVd = - Vbsz2_dVds * T6 ;
	T5_dVg = - Vbsz2_dVgs * T6 ;
	T7 = 1.0 / T5 ;
	bs12 = model->HSM2_bs1 * T7 ;
	T8 = - bs12 * T7 ;
	bs12_dVb = T8 * T5_dVb ;
	bs12_dVd = T8 * T5_dVd ;
	bs12_dVg = T8 * T5_dVg ;

	T1     = 0.93 * ( PS0Z_SCE + Ps0_min - Vbsz2 );
	T1_dVb = 0.93 * ( PS0Z_SCE_dVbs + Ps0_min_dVbs - Vbsz2_dVbs );  
	T1_dVd = 0.93 * ( PS0Z_SCE_dVds + Ps0_min_dVds - Vbsz2_dVds ); 
	T1_dVg = 0.93 * ( PS0Z_SCE_dVgs + Ps0_min_dVgs - Vbsz2_dVgs );
        T1_dPS0Z = 0.93 ;

	T2 = bs12 ;
	T2_dVb = bs12_dVb ;
	T2_dVd = bs12_dVd ;
	T2_dVg = bs12_dVg ;

	Fn_SU2( T10 , T2 , T1 , vth_dlt, T0, T3 ) ;
	T10_dVb = T2_dVb * T0 + T1_dVb * T3 ;
	T10_dVd = T2_dVd * T0 + T1_dVd * T3 ;
	T10_dVg = T2_dVg * T0 + T1_dVg * T3 ;
        T10_dPS0Z = T1_dPS0Z * T3 ;

	T4     = here->HSM2_2qnsub_esi * ( PS0Z_SCE + Ps0_min - Vbsz2 - T10 ) ; 
	T4_dVb = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVbs + Ps0_min_dVbs - Vbsz2_dVbs - T10_dVb ) ;
	T4_dVd = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVds + Ps0_min_dVds - Vbsz2_dVds - T10_dVd ) ;
	T4_dVg = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVgs + Ps0_min_dVgs - Vbsz2_dVgs - T10_dVg ) ;
        T4_dPS0Z = here->HSM2_2qnsub_esi * ( 1.0 - T10_dPS0Z ) ;

        if (T4 > arg2){
	  Qbmm = sqrt ( T4 ) ;
	  Qbmm_dVb = 0.5 / Qbmm * T4_dVb ;
	  Qbmm_dVd = 0.5 / Qbmm * T4_dVd ;
	  Qbmm_dVg = 0.5 / Qbmm * T4_dVg ;
          Qbmm_dPS0Z = 0.5 / Qbmm * T4_dPS0Z ;
        } else {
          Qbmm = sqrt(arg2) + 0.5 / sqrt(arg2) * ( T4 - arg2) ;
          Qbmm_dVb = 0.5 / sqrt(arg2) * T4_dVb ;
          Qbmm_dVd = 0.5 / sqrt(arg2) * T4_dVd ;
          Qbmm_dVg = 0.5 / sqrt(arg2) * T4_dVg ;
          Qbmm_dPS0Z = 0.5 / sqrt(arg2) * T4_dPS0Z ;
        }

	dqb = ( Qb0 - Qbmm ) * Cox_inv ;
	dqb_dVb = ( Qb0_dVb - Qbmm_dVb ) * Cox_inv + ( Qb0 - Qbmm ) * Cox_inv_dVb ;
	dqb_dVd = ( Qb0_dVd - Qbmm_dVd ) * Cox_inv + ( Qb0 - Qbmm ) * Cox_inv_dVd ;
	dqb_dVg = ( Qb0_dVg - Qbmm_dVg ) * Cox_inv + ( Qb0 - Qbmm ) * Cox_inv_dVg ;
        dqb_dPS0Z =      - Qbmm_dPS0Z * Cox_inv ;

	if( codqb == 0 ){
	  dqb = 0.0 ;
	  dqb_dVb = 0.0 ;
	  dqb_dVd = 0.0 ;
	  dqb_dVg = 0.0 ;
          dqb_dPS0Z = 0.0 ;
	}

	T1 = 2.0 * C_QE * here->HSM2_nsubc * C_ESI ;
	T2 = sqrt( T1 ) ;
	T2_dVb = 0.0 ;
	T2_dVd = 0.0 ;
	T2_dVg = 0.0 ;
	
	Vth0     = PS0Z_SCE + Vfb +    T2 * Cox_inv ;
	Vth0_dVb = PS0Z_SCE_dVbs + T2_dVb * Cox_inv + T2 * Cox_inv_dVb ;
	Vth0_dVd = PS0Z_SCE_dVds + T2_dVd * Cox_inv + T2 * Cox_inv_dVd ;
	Vth0_dVg = PS0Z_SCE_dVgs + T2_dVg * Cox_inv + T2 * Cox_inv_dVg ;
        Vth0_dPS0Z = 1.0 ;
	
	T1     = C_ESI * Cox_inv ;
	T1_dVb = C_ESI * Cox_inv_dVb ;
	T1_dVd = C_ESI * Cox_inv_dVd ;
	T1_dVg = C_ESI * Cox_inv_dVg ;
	T2 = here->HSM2_wdplp ;

	T4 = 1.0e0 / ( model->HSM2_lp * model->HSM2_lp ) ;
	T3 = 2.0 * ( model->HSM2_vbi - Pb20 ) * T2 * T4 ;

	T5 = T1 * T3 ;
	T5_dVb = T1_dVb * T3 ;
	T5_dVd = T1_dVd * T3 ;
	T5_dVg = T1_dVg * T3 ;

	T6     = PS0Z_SCE      - Vbsz ;
	T6_dVb = PS0Z_SCE_dVbs - Vbsz_dVbs ;
	T6_dVd = PS0Z_SCE_dVds - Vbsz_dVds ;
	T6_dVg = PS0Z_SCE_dVgs ;
        T6_dPS0Z = 1.0 ;

  	Fn_SZ( T6, T6, C_sce_dlt, T0 );
	T6_dVb *= T0 ;
	T6_dVd *= T0 ;
	T6_dVg *= T0 ;
        T6_dPS0Z *= T0 ;

	dVth0     = T5 * sqrt( T6 ) ;
	dVth0_dVb = T5 * 0.5 / sqrt( T6 ) * T6_dVb + T5_dVb * sqrt( T6 );
	dVth0_dVd = T5 * 0.5 / sqrt( T6 ) * T6_dVd + T5_dVd * sqrt( T6 ) ;
	dVth0_dVg = T5 * 0.5 / sqrt( T6 ) * T6_dVg + T5_dVg * sqrt( T6 ) ;
	 dVth0_dPS0Z = T5 * 0.5 / sqrt( T6 ) * T6_dPS0Z ;

	T1     = Vthp     - Vth0 ;
	T1_dVb = Vthp_dVb - Vth0_dVb ;
	T1_dVd = Vthp_dVd - Vth0_dVd ;
	T1_dVg = Vthp_dVg - Vth0_dVg ;
        T1_dPS0Z = Vthp_dPS0Z - Vth0_dPS0Z ;

	T9     = PS0Z_SCE      - Vbsz2 ;
	T9_dVb = PS0Z_SCE_dVbs - Vbsz2_dVbs ;
	T9_dVd = PS0Z_SCE_dVds - Vbsz2_dVds ;
	T9_dVg = PS0Z_SCE_dVgs - Vbsz2_dVgs ;
        T9_dPS0Z = 1.0 ;

	T3     = pParam->HSM2_scp1 + pParam->HSM2_scp3 * T9 / model->HSM2_lp + pParam->HSM2_scp2 * Vdsz ;
	T3_dVb = pParam->HSM2_scp3 * T9_dVb / model->HSM2_lp ;
	T3_dVd = pParam->HSM2_scp3 * T9_dVd / model->HSM2_lp + pParam->HSM2_scp2 * Vdsz_dVds ; 
	T3_dVg = pParam->HSM2_scp3 * T9_dVg / model->HSM2_lp ;
        T3_dPS0Z = pParam->HSM2_scp3 * T9_dPS0Z / model->HSM2_lp ;

	Vdx = model->HSM2_scp21 + Vdsz ;
	Vdx2 = Vdx * Vdx + small ;
	T4 = Vdx * Vdx + small ;
	T4_dVb = 0.0 ;
	T4_dVd = 2.0 * Vdx * Vdsz_dVds ;
	T4_dVg = 0.0 ;

	T5 = 1.0 / T4 ;
	T5_dVb = - T4_dVb / T4 / T4 ;
	T5_dVd = - T4_dVd / T4 / T4 ;
	T5_dVg = - T4_dVg / T4 / T4 ;

	dVthLP     = T1 * dVth0 * T3 + dqb - here->HSM2_msc * T5  ;
	dVthLP_dVb = T1_dVb * dVth0 * T3 + T1 * dVth0_dVb * T3 + T1 * dVth0 * T3_dVb + dqb_dVb - here->HSM2_msc * T5_dVb  ;
	dVthLP_dVd = T1_dVd * dVth0 * T3 + T1 * dVth0_dVd * T3 + T1 * dVth0 * T3_dVd + dqb_dVd - here->HSM2_msc * T5_dVd  ;
	dVthLP_dVg = T1_dVg * dVth0 * T3 + T1 * dVth0_dVg * T3 + T1 * dVth0 * T3_dVg + dqb_dVg - here->HSM2_msc * T5_dVg  ;
        dVthLP_dPS0Z = T1_dPS0Z * dVth0 * T3 + T1 * dVth0_dPS0Z * T3 + T1 * dVth0 * T3_dPS0Z + dqb_dPS0Z   ;

      }else{  /* original */

	T1 = here->HSM2_2qnsub_esi ;
	T2 = model->HSM2_bs2 - Vbsz2 ;
	T3 = T2 + small ;
	T4 = sqrt (T3 * T3 + 4.0 * vth_dlt) ;
	T5 = 0.5 * (T3 + T4) ;
	T6 = 0.5 * (1.0 + T3 / T4) ;
	T5_dVb = - Vbsz2_dVbs * T6 ;
	T5_dVd = - Vbsz2_dVds * T6 ;
	T5_dVg = - Vbsz2_dVgs * T6 ;
	T7 = 1.0 / T5 ;
	bs12 = model->HSM2_bs1 * T7 ;
	T8 = - bs12 * T7 ;
	bs12_dVb = T8 * T5_dVb ;
	bs12_dVd = T8 * T5_dVd ;
	bs12_dVg = T8 * T5_dVg ;
	Fn_SU( T10 , Vbsz2 + bs12, 0.93 * Pb20, vth_dlt, T0) ;
	Qbmm = sqrt (T1 * (Pb20 - T10 )) ;
	T9 = T0 / Qbmm ;
	Qbmm_dVb = 0.5 * T1 * - (Vbsz2_dVbs + bs12_dVb) * T9 ;
	Qbmm_dVd = 0.5 * T1 * - (Vbsz2_dVds + bs12_dVd) * T9 ;
	Qbmm_dVg = 0.5 * T1 * - (Vbsz2_dVgs + bs12_dVg) * T9 ;

	dqb = (Qb0 - Qbmm) * Cox_inv ;
	dqb_dVb = Vthp_dVb - Qbmm_dVb * Cox_inv - Qbmm * Cox_inv_dVb ;
	dqb_dVd = Vthp_dVd - Qbmm_dVd * Cox_inv - Qbmm * Cox_inv_dVd ;
	dqb_dVg = Vthp_dVg - Qbmm_dVg * Cox_inv - Qbmm * Cox_inv_dVg ;

	T1 = 2.0 * C_QE * here->HSM2_nsubc * C_ESI ;
	T2 = sqrt( T1 * ( Pb2c - Vbsz2 ) ) ;
	Vth0 = Pb2c + Vfb + T2 * Cox_inv ;
	T3 = 0.5 * T1 / T2 * Cox_inv ;
	Vth0_dVb = T3 * ( - Vbsz2_dVbs ) + T2 * Cox_inv_dVb ;
	Vth0_dVd = T3 * ( - Vbsz2_dVds ) + T2 * Cox_inv_dVd ;
	Vth0_dVg = T3 * ( - Vbsz2_dVgs ) + T2 * Cox_inv_dVg ;

	T1 = C_ESI * Cox_inv ;
	T2 = here->HSM2_wdplp ;
	T4 = 1.0e0 / ( model->HSM2_lp * model->HSM2_lp ) ;
	T5 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * T1 * T2 * T4 ;
	dVth0 = T5 * sqrt_Pbsum ;
	T6 = 0.5 * T5 / sqrt_Pbsum ;
	T7 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
	T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
	dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
	dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
	dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;
      
	T1 = Vthp - Vth0 ;
	T2 = pParam->HSM2_scp1 + pParam->HSM2_scp3 * Pbsum / model->HSM2_lp ;
	T3 = T2 + pParam->HSM2_scp2 * Vdsz ;
    
	Vdx = model->HSM2_scp21 + Vdsz ;
	Vdx2 = Vdx * Vdx + small ;
      
	dVthLP = T1 * dVth0 * T3 + dqb - here->HSM2_msc / Vdx2 ;
	T4 = T1 * dVth0 * pParam->HSM2_scp3 / model->HSM2_lp ;
	dVthLP_dVb = (Vthp_dVb - Vth0_dVb) * dVth0 * T3 + T1 * dVth0_dVb * T3 
	  + T4 * Pbsum_dVb + dqb_dVb ;
	dVthLP_dVd = (Vthp_dVd - Vth0_dVd) * dVth0 * T3 + T1 * dVth0_dVd * T3 
	  + T4 * Pbsum_dVd
	  + T1 * dVth0 * pParam->HSM2_scp2 * Vdsz_dVds
	  + dqb_dVd
	  + 2.0e0 * here->HSM2_msc * Vdx * Vdsz_dVds / ( Vdx2 * Vdx2 ) ;
	dVthLP_dVg = (Vthp_dVg - Vth0_dVg) * dVth0 * T3 + T1 * dVth0_dVg * T3 
	  + T4 * Pbsum_dVg + dqb_dVg ;
      }

    } else {
      dVthLP = 0.0e0 ;
      dVthLP_dVb = 0.0e0 ;
      dVthLP_dVd = 0.0e0 ;
      dVthLP_dVg = 0.0e0 ;
    }

    /*---------------------------------------------------*
     * dVthSC : Short-channel effect induced by Vds.
     *-----------------*/

    if( corecip ){ 

      T3 = here->HSM2_lgate - model->HSM2_parl2 ;
      T4 = 1.0e0 / ( T3 * T3 ) ;
      T5 = pParam->HSM2_sc3 / here->HSM2_lgate ;

      T6     = pParam->HSM2_sc1 + T5 * ( PS0Z_SCE - Vbsz ) ;
      T6_dVb = T5 * ( PS0Z_SCE_dVbs - Vbsz_dVbs );
      T6_dVd = T5 * ( PS0Z_SCE_dVds - Vbsz_dVds );
      T6_dVg = T5 *   PS0Z_SCE_dVgs ;
       T6_dPS0Z = T5 ;

      /* QME:1  CORECIP:1 */
      if( pParam->HSM2_sc4 != 0 ){  
	T8     = pParam->HSM2_sc4 * Vdsz * ( PS0Z_SCE - Vbsz ) ;
	T8_dVd = pParam->HSM2_sc4 * ( Vdsz_dVds * ( PS0Z_SCE - Vbsz ) + Vdsz * ( PS0Z_SCE_dVds - Vbsz_dVds ) ) ;
	T8_dVb = pParam->HSM2_sc4 *                                     Vdsz * ( PS0Z_SCE_dVbs - Vbsz_dVbs )   ;
	T8_dVg = pParam->HSM2_sc4 *                                     Vdsz *   PS0Z_SCE_dVgs                 ;
        T8_dPS0Z = pParam->HSM2_sc4 * Vdsz ;

	T1     = T6     + pParam->HSM2_sc2 * Vdsz + T8 ;
	T1_dVd = T6_dVd + pParam->HSM2_sc2 * Vdsz_dVds + T8_dVd ;
	T1_dVb = T6_dVb + T8_dVb ;
	T1_dVg = T6_dVg + T8_dVg ;
        T1_dPS0Z = T6_dPS0Z + T8_dPS0Z ;
      }else{
	T1     = T6     + pParam->HSM2_sc2 * Vdsz      ;
	T1_dVb = T6_dVb                                ;
	T1_dVd = T6_dVd + pParam->HSM2_sc2 * Vdsz_dVds ;
	T1_dVg = T6_dVg                                ;
        T1_dPS0Z = T6_dPS0Z                      ;
      }

      T0 = C_ESI * here->HSM2_wdpl * 2.0e0 * ( model->HSM2_vbi - Pb20 ) * T4 ;

      T2     = T0 * Cox_inv ;
      T2_dVb = T0 * Cox_inv_dVb ;
      T2_dVd = T0 * Cox_inv_dVd ;
      T2_dVg = T0 * Cox_inv_dVg ;

      A      = T2 * T1 ;
      A_dVbs = T2 * T1_dVb + T1 * T2_dVb ;
      A_dVds = T2 * T1_dVd + T1 * T2_dVd ;
      A_dVgs = T2 * T1_dVg + T1 * T2_dVg ;
       A_dPS0Z = T2 * T1_dPS0Z      ;

      T9     = PS0Z_SCE      - Vbsz + Ps0_min ;
      T9_dVb = PS0Z_SCE_dVbs - Vbsz_dVbs + Ps0_min_dVbs ;
      T9_dVd = PS0Z_SCE_dVds - Vbsz_dVds + Ps0_min_dVds ;
      T9_dVg = PS0Z_SCE_dVgs             + Ps0_min_dVgs ;
       T9_dPS0Z = 1.0 ;

      if ( T9 > arg0 ) {
        T8 = sqrt( T9 ) ;
        T8_dVb = 0.5 * T9_dVb / T8 ;
        T8_dVd = 0.5 * T9_dVd / T8 ;
        T8_dVg = 0.5 * T9_dVg / T8 ;
        T8_dPS0Z = 0.5 * T9_dPS0Z / T8 ;
      } else {
        T8 = sqrt(arg0) + 0.5 / sqrt(arg0) * ( T9 - arg0) ;
        T8_dVb = 0.5 / sqrt(arg0) * T9_dVb ;
        T8_dVd = 0.5 / sqrt(arg0) * T9_dVd ;
        T8_dVg = 0.5 / sqrt(arg0) * T9_dVg ;
        T8_dPS0Z = 0.5 / sqrt(arg0) * T9_dPS0Z ;
      }

      dVthSC = A * T8 ;
      dVthSC_dVb = A * T8_dVb + A_dVbs * T8;
      dVthSC_dVd = A * T8_dVd + A_dVds * T8;
      dVthSC_dVg = A * T8_dVg + A_dVgs * T8;
      dVthSC_dPS0Z = A * T8_dPS0Z + A_dPS0Z * T8;


    }else{  /* original */

      T1 = C_ESI * Cox_inv ;
      T2 = here->HSM2_wdpl ;
      T3 = here->HSM2_lgate - model->HSM2_parl2 ;
      T4 = 1.0e0 / ( T3 * T3 ) ;
      T5 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * T1 * T2 * T4 ;

      dVth0 = T5 * sqrt_Pbsum ;
      T6 = T5 / 2.0 / sqrt_Pbsum ;
      T7 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
      T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
      dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
      dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
      dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;

      T1 = pParam->HSM2_sc3 / here->HSM2_lgate ;
      T4 = pParam->HSM2_sc1 + T1 * Pbsum ;
      T4_dVb = T1 * Pbsum_dVb ;
      T4_dVd = T1 * Pbsum_dVd ;
      T4_dVg = T1 * Pbsum_dVg ;

      /* QME:1  CORECIP:0 */
      if( pParam->HSM2_sc4 != 0 ){  
	T8     = pParam->HSM2_sc4 *   Vdsz * Pbsum ;
	T8_dVd = pParam->HSM2_sc4 * ( Vdsz_dVds * Pbsum + Vdsz * Pbsum_dVd ) ;
	T8_dVb = pParam->HSM2_sc4 *                       Vdsz * Pbsum_dVb   ;
	T8_dVg = pParam->HSM2_sc4 *                       Vdsz * Pbsum_dVg   ;

	T5     = T4     + pParam->HSM2_sc2 * Vdsz      + T8     ;
	T5_dVd = T4_dVd + pParam->HSM2_sc2 * Vdsz_dVds + T8_dVd ;
	T5_dVb = T4_dVb                                + T8_dVb ;
	T5_dVg = T4_dVg                                + T8_dVg ;
      }else{
	T5     = T4     + pParam->HSM2_sc2 * Vdsz      ;
	T5_dVb = T4_dVb                                ;
	T5_dVd = T4_dVd + pParam->HSM2_sc2 * Vdsz_dVds ;
	T5_dVg = T4_dVg                                ;
      }

      dVthSC = dVth0 * T5 ;
      dVthSC_dVb = dVth0_dVb * T5 + dVth0 * T5_dVb ;
      dVthSC_dVd = dVth0_dVd * T5 + dVth0 * T5_dVd ;
      dVthSC_dVg = dVth0_dVg * T5 + dVth0 * T5_dVg ;

    }


    /*---------------------------------------------------*
     * dVthW : narrow-channel effect.
     *-----------------*/
    T1 = 1.0 / Cox ;
    T2 = T1 * T1 ;
    T3 = 1.0 / ( Cox +  pParam->HSM2_wfc / Weff ) ;
    T4 = T3 * T3 ;
    T5 = T1 - T3 ;
    T6 = Qb0 * ( T2 - T4 ) ;

    if( corecip ){ 
      dVthW = Qb0 * T5 + pParam->HSM2_wvth0 / here->HSM2_wg ;
      dVthW_dVb = Qb0_dVb * T5 - Cox_dVb * T6 ;
      dVthW_dVd = Qb0_dVd * T5 - Cox_dVd * T6 ;
      dVthW_dVg = Qb0_dVg * T5 - Cox_dVg * T6 ;
    }else{ /* original */
      dVthW = Qb0 * T5 + pParam->HSM2_wvth0 / here->HSM2_wg ;
      dVthW_dVb = Qb0_dVb * T5 - Cox_dVb * T6 ;
      dVthW_dVd = Qb0_dVd * T5 - Cox_dVd * T6 ;
      dVthW_dVg =              - Cox_dVg * T6 ;
    }

    /* end of case flg_qme = 1 */

  } else {

    /* now case flg_qme = 0    */

    /*---------------------------------------------------*
     * Vthp : Vth with pocket.
     *-----------------*/

    if( corecip ){
      T1 = here->HSM2_2qnsub_esi ;
      Qb0 = sqrt ( T1 ) ;
      Qb0_dVb = 0.0 ;
      Qb0_dVd = 0.0 ;
      Qb0_dVg = 0.0 ;

      Vthp     = PS0Z_SCE + Vfb + Qb0 * Cox_inv + here->HSM2_ptovr;
      Vthp_dVb = PS0Z_SCE_dVbs + Qb0_dVb * Cox_inv ;
      Vthp_dVd = PS0Z_SCE_dVds + Qb0_dVd * Cox_inv ;
      Vthp_dVg = PS0Z_SCE_dVgs + Qb0_dVg * Cox_inv ;
       Vthp_dPS0Z = 1.0 ;
    }else{  /* original */
      T1 = here->HSM2_2qnsub_esi ;
      Qb0 = sqrt (T1 * (Pb20 - Vbsz2)) ;
      T2 = 0.5 * T1 / Qb0 ;
      Qb0_dVb = T2 * (- Vbsz2_dVbs) ;
      Qb0_dVd = T2 * (- Vbsz2_dVds) ;
      Qb0_dVg = T2 * (- Vbsz2_dVgs) ;

      Vthp = Pb20 + Vfb + Qb0 * Cox_inv + here->HSM2_ptovr;
      Vthp_dVb = Qb0_dVb * Cox_inv ;
      Vthp_dVd = Qb0_dVd * Cox_inv ;
      Vthp_dVg = Qb0_dVg * Cox_inv ;
    }
    
    Pb20b = Pb20 ;
    Pb20b_dVb = 0.0 ;
    Pb20b_dVd = 0.0 ;
    Pb20b_dVg = 0.0 ;
    
    T0 = 0.95 ;
    T1 = T0 * Pb20b - Vbsz2 - 1.0e-3 ;
    T1_dVb = T0 * Pb20b_dVb - Vbsz2_dVbs ;
    T1_dVd = T0 * Pb20b_dVd - Vbsz2_dVds ;
    T1_dVg = T0 * Pb20b_dVg - Vbsz2_dVgs ;
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
    Pbsum = Pb20b - T3 ;
    Pbsum_dVb = Pb20b_dVb - T3_dVb ;
    Pbsum_dVd = Pb20b_dVd - T3_dVd ;
    Pbsum_dVg = Pb20b_dVg - T3_dVg ;

    sqrt_Pbsum = sqrt( Pbsum ) ;

    /*-------------------------------------------*
     * dVthLP : Short-channel effect induced by pocket.
     * - Vth0 : Vth without pocket.
     *-----------------*/
    if ( model->HSM2_lp != 0.0 ) {

      if( corecip ){
	T1 = here->HSM2_2qnsub_esi ;
	T2 = model->HSM2_bs2 - Vbsz2 ;
	T3 = T2 + small ;
	T4 = sqrt (T3 * T3 + 4.0 * vth_dlt) ;
	T5 = 0.5 * (T3 + T4) ;
	T6 = 0.5 * (1.0 + T3 / T4) ;
	T5_dVb = - Vbsz2_dVbs * T6 ;
	T5_dVd = - Vbsz2_dVds * T6 ;
	T5_dVg = - Vbsz2_dVgs * T6 ;
	T7 = 1.0 / T5 ;
	bs12 = model->HSM2_bs1 * T7 ;
	T8 = - bs12 * T7 ;
	bs12_dVb = T8 * T5_dVb ;
	bs12_dVd = T8 * T5_dVd ;
	bs12_dVg = T8 * T5_dVg ;

	T1     = 0.93 * ( PS0Z_SCE + Ps0_min - Vbsz2 );
	T1_dVb = 0.93 * ( PS0Z_SCE_dVbs + Ps0_min_dVbs - Vbsz2_dVbs );  
	T1_dVd = 0.93 * ( PS0Z_SCE_dVds + Ps0_min_dVds - Vbsz2_dVds ); 
	T1_dVg = 0.93 * ( PS0Z_SCE_dVgs + Ps0_min_dVgs - Vbsz2_dVgs );
        T1_dPS0Z = 0.93 ;

	T2 = bs12 ;
	T2_dVb = bs12_dVb ;
	T2_dVd = bs12_dVd ;
	T2_dVg = bs12_dVg ;

	Fn_SU2( T10 , T2 , T1 , vth_dlt, T0, T3 ) ;
	T10_dVb = T2_dVb * T0 + T1_dVb * T3 ;
	T10_dVd = T2_dVd * T0 + T1_dVd * T3 ;
	T10_dVg = T2_dVg * T0 + T1_dVg * T3 ;
        T10_dPS0Z = T1_dPS0Z * T3 ;

	T4     = here->HSM2_2qnsub_esi * ( PS0Z_SCE + Ps0_min - Vbsz2 - T10 ) ; 
	T4_dVb = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVbs + Ps0_min_dVbs - Vbsz2_dVbs - T10_dVb ) ;
	T4_dVd = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVds + Ps0_min_dVds - Vbsz2_dVds - T10_dVd ) ;
	T4_dVg = here->HSM2_2qnsub_esi * ( PS0Z_SCE_dVgs + Ps0_min_dVgs - Vbsz2_dVgs - T10_dVg ) ;
        T4_dPS0Z = here->HSM2_2qnsub_esi * ( 1.0 - T10_dPS0Z ) ;

        if (T4 > arg2){
	  Qbmm = sqrt ( T4 ) ;
	  Qbmm_dVb = 0.5 / Qbmm * T4_dVb ;
	  Qbmm_dVd = 0.5 / Qbmm * T4_dVd ;
	  Qbmm_dVg = 0.5 / Qbmm * T4_dVg ;
          Qbmm_dPS0Z = 0.5 / Qbmm * T4_dPS0Z ;
        } else {
          Qbmm = sqrt(arg2) + 0.5 / sqrt(arg2) * ( T4 - arg2) ;
          Qbmm_dVb = 0.5 / sqrt(arg2) * T4_dVb ;
          Qbmm_dVd = 0.5 / sqrt(arg2) * T4_dVd ;
          Qbmm_dVg = 0.5 / sqrt(arg2) * T4_dVg ;
          Qbmm_dPS0Z = 0.5 / sqrt(arg2) * T4_dPS0Z ;
        }

	dqb = ( Qb0 - Qbmm ) * Cox_inv ;
	dqb_dVb = ( Qb0_dVb - Qbmm_dVb ) * Cox_inv ;
	dqb_dVd = ( Qb0_dVd - Qbmm_dVd ) * Cox_inv ;
	dqb_dVg = ( Qb0_dVg - Qbmm_dVg ) * Cox_inv ;
        dqb_dPS0Z = ( - Qbmm_dPS0Z ) * Cox_inv ;

        /*  W/O QME PART  */
	if( codqb == 0 ){
	  dqb = 0 ;
	  dqb_dVb = 0 ;
	  dqb_dVd = 0 ;
	  dqb_dVg = 0 ;
          dqb_dPS0Z = 0.0 ;
	}

	T1 = 2.0 * C_QE * here->HSM2_nsubc * C_ESI ;
	T2 = sqrt( T1 ) ;
	T2_dVb = 0.0 ;
	T2_dVd = 0.0 ;
	T2_dVg = 0.0 ;

	Vth0     = PS0Z_SCE + Vfb +    T2 * Cox_inv ;
	Vth0_dVb = PS0Z_SCE_dVbs + T2_dVb * Cox_inv ;
	Vth0_dVd = PS0Z_SCE_dVds + T2_dVd * Cox_inv ;
	Vth0_dVg = PS0Z_SCE_dVgs + T2_dVg * Cox_inv ;
        Vth0_dPS0Z = 1.0 ;

	T1 = C_ESI * Cox_inv ;
	T2 = here->HSM2_wdplp ;

	T4 = 1.0e0 / ( model->HSM2_lp * model->HSM2_lp ) ;
	T5 = 2.0e0 * ( model->HSM2_vbi - Pb20 ) * T1 * T2 * T4 ;
	T5_dVb = 0.0 ;
	T5_dVd = 0.0 ;
	T5_dVg = 0.0 ;

	T6     = PS0Z_SCE      - Vbsz      ;
	T6_dVb = PS0Z_SCE_dVbs - Vbsz_dVbs ;
	T6_dVd = PS0Z_SCE_dVds - Vbsz_dVds ;
	T6_dVg = PS0Z_SCE_dVgs                         ;
        T6_dPS0Z = 1.0 ;   

	Fn_SZ(T6, T6, C_sce_dlt, T0 );
        T6 += small ;
	T6_dVb *= T0 ;
	T6_dVd *= T0 ;
	T6_dVg *= T0 ;
        T6_dPS0Z *= T0 ;

	dVth0 = T5 * sqrt( T6 ) ;
	dVth0_dVb = T5 * 0.5 / sqrt( T6 ) * T6_dVb + T5_dVb * sqrt( T6 );
	dVth0_dVd = T5 * 0.5 / sqrt( T6 ) * T6_dVd + T5_dVd * sqrt( T6 ) ;
	dVth0_dVg = T5 * 0.5 / sqrt( T6 ) * T6_dVg + T5_dVg * sqrt( T6 ) ;
        dVth0_dPS0Z = T5 * 0.5 / sqrt( T6 ) * T6_dPS0Z ;

	T1 = Vthp - Vth0 ;
	T1_dVb = Vthp_dVb - Vth0_dVb ;
	T1_dVd = Vthp_dVd - Vth0_dVd ;
	T1_dVg = Vthp_dVg - Vth0_dVg ;
        T1_dPS0Z = Vthp_dPS0Z - Vth0_dPS0Z ;

	T9 = PS0Z_SCE - Vbsz2 ;
	T9_dVb = PS0Z_SCE_dVbs - Vbsz2_dVbs ;
	T9_dVd = PS0Z_SCE_dVds - Vbsz2_dVds ;
	T9_dVg = PS0Z_SCE_dVgs - Vbsz2_dVgs ;
        T9_dPS0Z = 1.0 ;

	T3 = pParam->HSM2_scp1 + pParam->HSM2_scp3 * T9 / model->HSM2_lp + pParam->HSM2_scp2 * Vdsz ;
	T3_dVb = pParam->HSM2_scp3 * T9_dVb / model->HSM2_lp ;
	T3_dVd = pParam->HSM2_scp3 * T9_dVd / model->HSM2_lp + pParam->HSM2_scp2 * Vdsz_dVds ; 
	T3_dVg = pParam->HSM2_scp3 * T9_dVg / model->HSM2_lp ;
        T3_dPS0Z = pParam->HSM2_scp3 * T9_dPS0Z / model->HSM2_lp ;

	Vdx = model->HSM2_scp21 + Vdsz ;
	Vdx2 = Vdx * Vdx + small ;
	T4 = Vdx * Vdx + small ;
	T4_dVb = 0.0 ;
	T4_dVd = 2.0 * Vdx * Vdsz_dVds ;
	T4_dVg = 0.0 ;

	T5 = 1.0 / T4 ;
	T5_dVb = - T4_dVb / T4 / T4 ;
	T5_dVd = - T4_dVd / T4 / T4 ;
	T5_dVg = - T4_dVg / T4 / T4 ;

	dVthLP     = T1 * dVth0 * T3 + dqb - here->HSM2_msc * T5  ;
	dVthLP_dVb = T1_dVb * dVth0 * T3 + T1 * dVth0_dVb * T3 + T1 * dVth0 * T3_dVb + dqb_dVb - here->HSM2_msc * T5_dVb  ;
	dVthLP_dVd = T1_dVd * dVth0 * T3 + T1 * dVth0_dVd * T3 + T1 * dVth0 * T3_dVd + dqb_dVd - here->HSM2_msc * T5_dVd  ;
	dVthLP_dVg = T1_dVg * dVth0 * T3 + T1 * dVth0_dVg * T3 + T1 * dVth0 * T3_dVg + dqb_dVg - here->HSM2_msc * T5_dVg  ;
        dVthLP_dPS0Z = T1_dPS0Z * dVth0 * T3 + T1 * dVth0_dPS0Z * T3 + T1 * dVth0 * T3_dPS0Z + dqb_dPS0Z   ;

      }else{  /* Original */

	T1 = here->HSM2_2qnsub_esi ;
	T2 = model->HSM2_bs2 - Vbsz2 ;
	T3 = T2 + small ;
	T4 = sqrt (T3 * T3 + 4.0 * vth_dlt) ;
	T5 = 0.5 * (T3 + T4) ;
	T6 = 0.5 * (1.0 + T3 / T4) ;
	T5_dVb = - Vbsz2_dVbs * T6 ;
	T5_dVd = - Vbsz2_dVds * T6 ;
	T5_dVg = - Vbsz2_dVgs * T6 ;
	T7 = 1.0 / T5 ;
	bs12 = model->HSM2_bs1 * T7 ;
	T8 = - bs12 * T7 ;
	bs12_dVb = T8 * T5_dVb ;
	bs12_dVd = T8 * T5_dVd ;
	bs12_dVg = T8 * T5_dVg ;
	Fn_SU( T10 , Vbsz2 + bs12, 0.93 * Pb20, vth_dlt, T0) ;
	Qbmm = sqrt (T1 * (Pb20 - T10 )) ;

	T9 = T0 / Qbmm ;
	Qbmm_dVb = 0.5 * T1 * - (Vbsz2_dVbs + bs12_dVb) * T9 ;
	Qbmm_dVd = 0.5 * T1 * - (Vbsz2_dVds + bs12_dVd) * T9 ;
	Qbmm_dVg = 0.5 * T1 * - (Vbsz2_dVgs + bs12_dVg) * T9 ;

	dqb = (Qb0 - Qbmm) * Cox_inv ;
	dqb_dVb = Vthp_dVb - Qbmm_dVb * Cox_inv ;
	dqb_dVd = Vthp_dVd - Qbmm_dVd * Cox_inv ;
	dqb_dVg = Vthp_dVg - Qbmm_dVg * Cox_inv ;

	T1 = 2.0 * C_QE * here->HSM2_nsubc * C_ESI ;
	T2 = sqrt( T1 * ( Pb2c - Vbsz2 ) ) ;
	Vth0 = Pb2c + Vfb + T2 * Cox_inv ;
	T3 = 0.5 * T1 / T2 * Cox_inv ;
	Vth0_dVb = T3 * ( - Vbsz2_dVbs ) ;
	Vth0_dVd = T3 * ( - Vbsz2_dVds ) ;
	Vth0_dVg = T3 * ( - Vbsz2_dVgs ) ;

	T1 = C_ESI * Cox_inv ;
	T2 = here->HSM2_wdplp ;

	T4 = 1.0e0 / ( model->HSM2_lp * model->HSM2_lp ) ;
	T5 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * T1 * T2 * T4 ;
	dVth0 = T5 * sqrt_Pbsum ;
	T6 = 0.5 * T5 / sqrt_Pbsum ;
	T7 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
	T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
	dVth0_dVb = T6 * Pbsum_dVb + T8 * Pb20b_dVb ;
	dVth0_dVd = T6 * Pbsum_dVd + T8 * Pb20b_dVd ;
	dVth0_dVg = T6 * Pbsum_dVg + T8 * Pb20b_dVg ;
      
	T1 = Vthp - Vth0 ;
	T2 = pParam->HSM2_scp1 + pParam->HSM2_scp3 * Pbsum / model->HSM2_lp ;
	T3 = T2 + pParam->HSM2_scp2 * Vdsz ;
    
	Vdx = model->HSM2_scp21 + Vdsz ;
	Vdx2 = Vdx * Vdx + small ;
      
	dVthLP = T1 * dVth0 * T3 + dqb - here->HSM2_msc / Vdx2 ;
	T4 = T1 * dVth0 * pParam->HSM2_scp3 / model->HSM2_lp ;
	dVthLP_dVb = (Vthp_dVb - Vth0_dVb) * dVth0 * T3 + T1 * dVth0_dVb * T3 
	  + T4 * Pbsum_dVb + dqb_dVb ;
	dVthLP_dVd = (Vthp_dVd - Vth0_dVd) * dVth0 * T3 + T1 * dVth0_dVd * T3 
	  + T4 * Pbsum_dVd
	  + T1 * dVth0 * pParam->HSM2_scp2 * Vdsz_dVds
	  + dqb_dVd
	  + 2.0e0 * here->HSM2_msc * Vdx * Vdsz_dVds / ( Vdx2 * Vdx2 ) ;
	dVthLP_dVg = (Vthp_dVg - Vth0_dVg) * dVth0 * T3 + T1 * dVth0_dVg * T3 
	  + T4 * Pbsum_dVg + dqb_dVg ;

      }

    } else {
      dVthLP = 0.0e0 ;
      dVthLP_dVb = 0.0e0 ;
      dVthLP_dVd = 0.0e0 ;
      dVthLP_dVg = 0.0e0 ;
    }

    /*---------------------------------------------------*
     * dVthSC : Short-channel effect induced by Vds.
     *-----------------*/

    if( corecip ){
      T3 = here->HSM2_lgate - model->HSM2_parl2 ;
      T4 = 1.0e0 / ( T3 * T3 ) ;
      T5 = pParam->HSM2_sc3 / here->HSM2_lgate ;

      T6     = pParam->HSM2_sc1 + T5 * ( PS0Z_SCE - Vbsz ) ;
      T6_dVb = T5 * ( PS0Z_SCE_dVbs - Vbsz_dVbs );
      T6_dVd = T5 * ( PS0Z_SCE_dVds - Vbsz_dVds );
      T6_dVg = T5 *   PS0Z_SCE_dVgs ;
      T6_dPS0Z = T5 ;

      /* QME:0  CORECIP:1 */
      if( pParam->HSM2_sc4 != 0 ){
	T8     = pParam->HSM2_sc4 * Vdsz * ( PS0Z_SCE - Vbsz ) ;
	T8_dVd = pParam->HSM2_sc4 * ( Vdsz_dVds * ( PS0Z_SCE - Vbsz ) + Vdsz * ( PS0Z_SCE_dVds - Vbsz_dVds ) ) ;
	T8_dVb = pParam->HSM2_sc4 *                                     Vdsz * ( PS0Z_SCE_dVbs - Vbsz_dVbs )   ;
	T8_dVg = pParam->HSM2_sc4 *                                     Vdsz *   PS0Z_SCE_dVgs                 ;
        T8_dPS0Z = pParam->HSM2_sc4 * Vdsz ;

	T1     = T6     + pParam->HSM2_sc2 * Vdsz + T8 ;
	T1_dVd = T6_dVd + pParam->HSM2_sc2 * Vdsz_dVds + T8_dVd ;
	T1_dVb = T6_dVb + T8_dVb ;
	T1_dVg = T6_dVg + T8_dVg ;
        T1_dPS0Z = T6_dPS0Z + T8_dPS0Z ;
      }else{
	T1     = T6     + pParam->HSM2_sc2 * Vdsz      ;
	T1_dVb = T6_dVb                                ;
	T1_dVd = T6_dVd + pParam->HSM2_sc2 * Vdsz_dVds ;
	T1_dVg = T6_dVg                                ;
        T1_dPS0Z = T6_dPS0Z                      ;
      }

      T2 = C_ESI * Cox_inv * here->HSM2_wdpl * 2.0e0 * ( model->HSM2_vbi - Pb20 ) * T4 ;

      A = T2 * T1 ;
      A_dVbs = T2 * T1_dVb ;
      A_dVds = T2 * T1_dVd ;
      A_dVgs = T2 * T1_dVg ;
       A_dPS0Z = T2 * T1_dPS0Z ;

      T7     = PS0Z_SCE      - Vbsz + Ps0_min ;
      T7_dVb = PS0Z_SCE_dVbs - Vbsz_dVbs + Ps0_min_dVbs;
      T7_dVd = PS0Z_SCE_dVds - Vbsz_dVds + Ps0_min_dVds ;
      T7_dVg = PS0Z_SCE_dVgs             + Ps0_min_dVgs ;
      T7_dPS0Z = 1.0 ;

      if ( T7 > arg0 ) {
        T8 = sqrt( T7 ) ;
        T8_dVb = 0.5 * T7_dVb / T8 ;
        T8_dVd = 0.5 * T7_dVd / T8 ;
        T8_dVg = 0.5 * T7_dVg / T8 ;
        T8_dPS0Z  = 0.5 * T7_dPS0Z  / T8 ;
      } else {
        T8 = sqrt(arg0) + 0.5 / sqrt(arg0) * ( T7 - arg0) ;
        T8_dVb = 0.5 / sqrt(arg0) * T7_dVb ;
        T8_dVd = 0.5 / sqrt(arg0) * T7_dVd ;
        T8_dVg = 0.5 / sqrt(arg0) * T7_dVg ;
        T8_dPS0Z = 0.5 / sqrt(arg0) * T7_dPS0Z ;
      }

      dVthSC = A * T8 ;
      dVthSC_dVb = A * T8_dVb + A_dVbs * T8;
      dVthSC_dVd = A * T8_dVd + A_dVds * T8;
      dVthSC_dVg = A * T8_dVg + A_dVgs * T8;
      dVthSC_dPS0Z = A * T8_dPS0Z + A_dPS0Z * T8;

    }else{ /* original */

      T1 = C_ESI * Cox_inv ;
      T2 = here->HSM2_wdpl ;
      T3 = here->HSM2_lgate - model->HSM2_parl2 ;
      T4 = 1.0e0 / ( T3 * T3 ) ;
      T5 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * T1 * T2 * T4 ;

      dVth0 = T5 * sqrt_Pbsum ;
      T6 = T5 / 2.0 / sqrt_Pbsum ;
      T7 = 2.0e0 * ( model->HSM2_vbi - Pb20b ) * C_ESI * T2 * T4 * sqrt_Pbsum ;
      T8 = - 2.0e0 * T1 * T2 * T4 * sqrt_Pbsum ;
      dVth0_dVb = T6 * Pbsum_dVb + T8 * Pb20b_dVb ;
      dVth0_dVd = T6 * Pbsum_dVd + T8 * Pb20b_dVd ;
      dVth0_dVg = T6 * Pbsum_dVg + T8 * Pb20b_dVg ;

      T1 = pParam->HSM2_sc3 / here->HSM2_lgate ;
      T4 = pParam->HSM2_sc1 + T1 * Pbsum ;
      T4_dVb = T1 * Pbsum_dVb ;
      T4_dVd = T1 * Pbsum_dVd ;
      T4_dVg = T1 * Pbsum_dVg ;

      /* QME:0  CORECIP:0 */
      if( pParam->HSM2_sc4 != 0 ){
	T8     = pParam->HSM2_sc4 *   Vdsz * Pbsum ;
	T8_dVd = pParam->HSM2_sc4 * ( Vdsz_dVds * Pbsum + Vdsz * Pbsum_dVd ) ;
	T8_dVb = pParam->HSM2_sc4 *                       Vdsz * Pbsum_dVb   ;
	T8_dVg = pParam->HSM2_sc4 *                       Vdsz * Pbsum_dVg   ;

	T5     = T4     + pParam->HSM2_sc2 * Vdsz      + T8     ;
	T5_dVd = T4_dVd + pParam->HSM2_sc2 * Vdsz_dVds + T8_dVd ;
	T5_dVb = T4_dVb                                + T8_dVb ;
	T5_dVg = T4_dVg                                + T8_dVg ;
      }else{
	T5     = T4     + pParam->HSM2_sc2 * Vdsz      ;
	T5_dVb = T4_dVb                                ;
	T5_dVd = T4_dVd + pParam->HSM2_sc2 * Vdsz_dVds ;
	T5_dVg = T4_dVg                                ;
      }

      dVthSC = dVth0 * T5 ;
      dVthSC_dVb = dVth0_dVb * T5 + dVth0 * T5_dVb ;
      dVthSC_dVd = dVth0_dVd * T5 + dVth0 * T5_dVd ;
      dVthSC_dVg = dVth0_dVg * T5 + dVth0 * T5_dVg ;

    }

    /*---------------------------------------------------*
     * dVthW : narrow-channel effect.
     *-----------------*/
    T1 = 1.0 / Cox ;
    T3 = 1.0 / ( Cox +  pParam->HSM2_wfc / Weff ) ;
    T5 = T1 - T3 ;

    if( corecip ){
      dVthW = Qb0 * T5 + pParam->HSM2_wvth0 / here->HSM2_wg ;
      dVthW_dVb = Qb0_dVb * T5 ;
      dVthW_dVd = Qb0_dVd * T5 ;
      dVthW_dVg = Qb0_dVg * T5 ;
    }else{ /* original */
      dVthW = Qb0 * T5 + pParam->HSM2_wvth0 / here->HSM2_wg ;
      dVthW_dVb = Qb0_dVb * T5 ;
      dVthW_dVd = Qb0_dVd * T5 ;
      dVthW_dVg = 0.0 ;
    }



  } /* end of flg_qme if-blocks */

    /*---------------------------------------------------*
     * dVth : Total variation. 
     * - Positive dVth means the decrease in Vth.
     *-----------------*/
    dVth = dVthSC + dVthLP + dVthW + here->HSM2_dVthsm ;
    dVth_dVb = dVthSC_dVb + dVthLP_dVb + dVthW_dVb ;
    dVth_dVd = dVthSC_dVd + dVthLP_dVd + dVthW_dVd ;
    dVth_dVg = dVthSC_dVg + dVthLP_dVg + dVthW_dVg ;
    dVth_dPS0Z = dVthSC_dPS0Z + dVthLP_dPS0Z ;

    /*---------------------------------------------------*
     * Vth : Threshold voltagei for OP. 
     *-----------------*/
    T2 = sqrt( here->HSM2_2qnsub_esi * Pb2 ) ;
    Vth = Pb2 + Vfb + T2 * Cox0_inv - dVth ;

    /*-----------------------------------------------------------*
     * Constants in the equation of Ps0 . 
     *-----------------*/

    fac1 = cnst0 * Cox_inv ;
    fac1_dVbs = cnst0 * Cox_inv_dVb ;
    fac1_dVds = cnst0 * Cox_inv_dVd ;
    fac1_dVgs = cnst0 * Cox_inv_dVg ;

    fac1p2 = fac1 * fac1 ;
    /*---------------------------------------------------*
     * Poly-Depletion Effect
     *-----------------*/

  if ( here->HSM2_flg_pgd == 0 ) {
    dPpg = 0.0 ;
    dPpg_dVb = 0.0 ;
    dPpg_dVd = 0.0 ;
    dPpg_dVg = 0.0 ;
  } else {

    T7 = Vgsz ;
    T7_dVb =  0.0 ;
    T7_dVd =  Vgsz_dVds ;
    T7_dVg =  Vgsz_dVgs ;

    T0 = here->HSM2_cnstpgd ;

    T3 = T7 - model->HSM2_pgd2 ;
    T3_dVb = T7_dVb ;
    T3_dVd = T7_dVd ;
    T3_dVg = T7_dVg ;

    Fn_ExpLim( dPpg , T3 , T6 ) ;
    dPpg_dVb = T6 * T3_dVb ;
    dPpg_dVd = T6 * T3_dVd ;
    dPpg_dVg = T6 * T3_dVg ;

    Fn_SZ( dPpg , dPpg - 1.0 , 0.1 , T6 ) ;
    dPpg_dVb *= T6 ;
    dPpg_dVd *= T6 ;
    dPpg_dVg *= T6 ;

    dPpg *= T0 ;
    dPpg_dVb *= T0 ;
    dPpg_dVd *= T0 ;
    dPpg_dVg *= T0 ;

    Fn_SU( dPpg , dPpg , pol_b , pol_dlt , T9 ) ;
    dPpg_dVb *= T9 ;
    dPpg_dVd *= T9 ;
    dPpg_dVg *= T9 ;

  }

    /*---------------------------------------------------*
     * Vgp : Effective gate bias with SCE & RSCE & flatband. 
     *-----------------*/
    Vgp = Vgs - Vfb + dVth - dPpg ;
    Vgp_dVbs = dVth_dVb - dPpg_dVb ;
    Vgp_dVds = dVth_dVd - dPpg_dVd ;
    Vgp_dVgs = 1.0e0 + dVth_dVg - dPpg_dVg ;
    Vgp_dPS0Z = dVth_dPS0Z ;
   
    Vgpz =  Vgsz - Vfb + dVth - dPpg ; /* (tmp) */
    Vgpz_dVbs = dVth_dVb ;
    Vgpz_dVds = Vgsz_dVds + dVth_dVd - dPpg_dVd ;
    Vgpz_dVgs = Vgsz_dVgs + dVth_dVg - dPpg_dVg ;
    Vgpz_dPS0Z = dVth_dPS0Z ;
 
    /*---------------------------------------------------*
     * Vgs_fb : Actual flatband voltage taking account Vbs. 
     * - note: if Vgs == Vgs_fb then Vgp == Ps0 == Vbs . 
     *------------------*/
    Vgs_fb = Vfb - dVth + dPpg + Vbs ;

    
    /*---------------------------------------------------*
     * Vfbsft : Vfb shift 
     *-----------------*/
    Vfbsft = 0.0 ;
    Vfbsft_dVbs = 0.0 ;
    Vfbsft_dVds = 0.0 ;
    Vfbsft_dVgs = 0.0 ;
    Vfbsft_dPS0Z = 0.0 ;

    if ( Vbs > 0.0 ) {
      /* values at D2/D3 boundary + beta */
      /* Ps0 */
      T1 = Vbs + ( znbd5 + 1 ) * beta_inv ; 
      /* Qb0 */
      /* T2 = cnst0 * sqrt( znbd5 ) */ 
      T2 = cnst0 * 2.23606797749979 ;

      /* Vgp assuming Qn0=0 */
      T3 = T2 * Cox_inv + T1 ; 

      /* Vgp difference */
      TX = T3 - Vgp ;
      TX_dVbs = 1.0 - Vgp_dVbs ;
      TX_dVds = - Vgp_dVds ;
      TX_dVgs = - Vgp_dVgs ;
      TX_dPS0Z = - Vgp_dPS0Z ;

      /* set lower limit to 0 */
      Fn_SZ( TX , TX , 0.1 , T4 ) ;
      TX_dVbs *= T4 ;
      TX_dVds *= T4 ;
      TX_dVgs *= T4 ;
      TX_dPS0Z *= T4 ;

      /* TY: damping factor */
      T1 = 0.5 ;
      T5 = Vbs / T1 ;
      T5_dVb = 1.0 / T1 ;
      T6 = T5 * T5 ;
      T6 *= T6 ; 
      T6_dVb = 4 * T5 * T5 * T5 * T5_dVb ;
      T7 = 1.0 / ( 1.0 + T6 ) ;
      T8 = T7 * T7 ;
      TY = 1.0 - T7 ;
      TY_dVbs = T8 * T6_dVb ;

      TX = TY = 0.0 ;
      TX_dVbs = TX_dVds = TX_dVgs = TX_dPS0Z = TY_dVbs = 0.0 ;

      Vfbsft = TX * TY ;
      Vfbsft_dVbs = TX_dVbs * TY + TX * TY_dVbs ;
      Vfbsft_dVds = TX_dVds * TY ;
      Vfbsft_dVgs = TX_dVgs * TY ;
      Vfbsft_dPS0Z = TX_dPS0Z * TY ;
      Vgs_fb -= Vfbsft ;

      Vgp += Vfbsft ;
      Vgp_dVbs += Vfbsft_dVbs ;
      Vgp_dVds += Vfbsft_dVds ;
      Vgp_dVgs += Vfbsft_dVgs ;
      Vgp_dPS0Z += Vfbsft_dPS0Z ;
   
      Vgpz += Vfbsft ; 
      Vgpz_dVbs += Vfbsft_dVbs ;
      Vgpz_dVds += Vfbsft_dVds ;
      Vgpz_dVgs += Vfbsft_dVgs ;
      Vgpz_dPS0Z += Vfbsft_dPS0Z ;
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

      /* Ps0_min = here->HSM2_eg - Pb2 ;                                    */
      /* -> replaced by approx. solving Poisson equation at Vgs=Vgs_min     */
      /* Ps0_min = 2.0 * beta_inv * log(-Vgs_min/fac1) ; already done above */
  
      TX = beta * ( Vgp - Vbs ) ;
      TX_dVbs = beta * ( Vgp_dVbs - 1.0 ) ;
      TX_dVds = beta * Vgp_dVds ;
      TX_dVgs = beta * Vgp_dVgs ;
      TX_dPS0Z = beta * Vgp_dPS0Z ;
      
      T1 = 1.0 / ( beta * cnst0 ) ;
      TY = T1 * Cox ;
      TY_dVbs = T1 * Cox_dVb ;
      TY_dVds = T1 * Cox_dVd ;
      TY_dVgs = T1 * Cox_dVg ;

      Ac41 = 2.0 + 3.0 * C_SQRT_2 * TY ;
      
      Ac4 = 8.0 * Ac41 * Ac41 * Ac41 ;
      T1 = 72.0 * Ac41 * Ac41 * C_SQRT_2 ;
      Ac4_dVbs = T1 * TY_dVbs ;
      Ac4_dVds = T1 * TY_dVds ;
      Ac4_dVgs = T1 * TY_dVgs ;

      T4 = ( TX - 2.0 ) ;
      T5 = 9.0 * TY * T4 ;
      T5_dVb = 9.0 * ( TY_dVbs * T4 + TY * TX_dVbs ) ;
      T5_dVd = 9.0 * ( TY_dVds * T4 + TY * TX_dVds ) ;
      T5_dVg = 9.0 * ( TY_dVgs * T4 + TY * TX_dVgs ) ;
      T5_dPS0Z = 9.0 * (              TY * TX_dPS0Z ) ;

      Ac31 = 7.0 * C_SQRT_2 - T5 ;
      Ac31_dVbs = -T5_dVb ;
      Ac31_dVds = -T5_dVd ;
      Ac31_dVgs = -T5_dVg ;
      Ac31_dPS0Z = -T5_dPS0Z ;

      Ac3 = Ac31 * Ac31 ;
      T1 = 2.0 * Ac31 ;
      Ac3_dVbs = T1 * Ac31_dVbs ;
      Ac3_dVds = T1 * Ac31_dVds ;
      Ac3_dVgs = T1 * Ac31_dVgs ;
      Ac3_dPS0Z = T1 * Ac31_dPS0Z ;

      Ac2 = sqrt( Ac4 + Ac3 ) ;
      T1 = 0.5 / Ac2 ;
      Ac2_dVbs = T1 * ( Ac4_dVbs + Ac3_dVbs ) ;
      Ac2_dVds = T1 * ( Ac4_dVds + Ac3_dVds ) ;
      Ac2_dVgs = T1 * ( Ac4_dVgs + Ac3_dVgs ) ;
      Ac2_dPS0Z = T1 * (           Ac3_dPS0Z ) ;
      
      Ac1 = -7.0 * C_SQRT_2 + Ac2 + T5 ;
      Ac1_dVbs = Ac2_dVbs + T5_dVb ;
      Ac1_dVds = Ac2_dVds + T5_dVd ;
      Ac1_dVgs = Ac2_dVgs + T5_dVg ;
       Ac1_dPS0Z = Ac2_dPS0Z + T5_dPS0Z ;

      Acd = Fn_Pow( Ac1 , C_1o3 ) ;
      T1 = C_1o3 / ( Acd * Acd ) ;
      Acd_dVbs = Ac1_dVbs * T1 ;
      Acd_dVds = Ac1_dVds * T1 ;
      Acd_dVgs = Ac1_dVgs * T1 ;
      Acd_dPS0Z = Ac1_dPS0Z * T1 ;

      Acn = -4.0 * C_SQRT_2 - 12.0 * TY + 2.0 * Acd + C_SQRT_2 * Acd * Acd ;
      T1 = 2.0 + 2.0 * C_SQRT_2 * Acd ;
      Acn_dVbs = - 12.0 * TY_dVbs + T1 * Acd_dVbs ;
      Acn_dVds = - 12.0 * TY_dVds + T1 * Acd_dVds ;
      Acn_dVgs = - 12.0 * TY_dVgs + T1 * Acd_dVgs ;
      Acn_dPS0Z =                   T1 * Acd_dPS0Z ;

      T1 = 1.0 / Acd ;
      Chi = Acn * T1 ;
      Chi_dVbs = ( Acn_dVbs - Chi * Acd_dVbs ) * T1 ;
      Chi_dVds = ( Acn_dVds - Chi * Acd_dVds ) * T1 ;
      Chi_dVgs = ( Acn_dVgs - Chi * Acd_dVgs ) * T1 ;
      Chi_dPS0Z = ( Acn_dPS0Z - Chi * Acd_dPS0Z ) * T1 ;

      Psa = Chi * beta_inv + Vbs ;
      Psa_dVbs = Chi_dVbs * beta_inv + 1.0 ;
      Psa_dVds = Chi_dVds * beta_inv ;
      Psa_dVgs = Chi_dVgs * beta_inv ;
      Psa_dPS0Z = Chi_dPS0Z * beta_inv ;
      
      T1 = Psa - Vbs ;
      T2 = T1 / Ps0_min ;
      T3 = sqrt( 1.0 + ( T2 * T2 ) ) ;

      T9 = T2 / T3 / Ps0_min ;
      T3_dVb = T9 * ( Psa_dVbs - 1.0 ) ;
      T3_dVd = T9 * ( Psa_dVds ) ;
      T3_dVg = T9 * ( Psa_dVgs ) ;
      T3_dPS0Z = T9 * ( Psa_dPS0Z ) ;
      
      Ps0 = T1 / T3 + Vbs ;
      T9 = 1.0 / ( T3 * T3 ) ;
      Ps0_dVbs = T9 * ( ( Psa_dVbs - 1.0 ) * T3 - T1 * T3_dVb ) + 1.0 ;
      Ps0_dVds = T9 * ( Psa_dVds * T3 - T1 * T3_dVd ) ;
      Ps0_dVgs = T9 * ( Psa_dVgs * T3 - T1 * T3_dVg ) ;
      Ps0_dPS0Z = T9 * ( Psa_dPS0Z * T3 - T1 * T3_dPS0Z ) ;

      /*---------------------------------------------------*
       * Characteristics. 
       *-----------------*/
      Psl = Ps0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVds = Ps0_dVds ;
      Psl_dVgs = Ps0_dVgs ;
      
      T2 = ( Vgp - Ps0 ) ;
      Qbu = Cox * T2 ;
      Qbu_dVbs = Cox * ( Vgp_dVbs - Ps0_dVbs ) + Cox_dVb * T2 ;
      Qbu_dVds = Cox * ( Vgp_dVds - Ps0_dVds ) + Cox_dVd * T2 ;
      Qbu_dVgs = Cox * ( Vgp_dVgs - Ps0_dVgs ) + Cox_dVg * T2 ;
    
      Qiu = 0.0e0 ;
      Qiu_dVbs = 0.0e0 ;
      Qiu_dVds = 0.0e0 ;
      Qiu_dVgs = 0.0e0 ;
 
      Qdrat = 0.5e0 ;
      Qdrat_dVbs = 0.0e0 ;
      Qdrat_dVds = 0.0e0 ;
      Qdrat_dVgs = 0.0e0 ;
 
      Lred = 0.0e0 ;
      Lred_dVbs = 0.0e0 ;
      Lred_dVds = 0.0e0 ;
      Lred_dVgs = 0.0e0 ;
 
      Ids = 0.0e0 ;
      Ids_dVbs = 0.0e0 ;
      Ids_dVds = 0.0e0 ;
      Ids_dVgs = 0.0e0 ;

      VgVt = 0.0 ;
      
      flg_noqi = 1 ;

      DJI = 1.0 ;

      if( corecip ){
	/*!! This is for Accumulation Region !!*/

	T1 = Vds * 0.5  ;
	Fn_SymAdd( Pzadd , T1 , model->HSM2_pzadd0 , T2 ) ;
	T2 /= 2 ;
	Pzadd_dVbs = 0.0 ;
	Pzadd_dVds = T2 * ( 1.0 ) ;
	Pzadd_dVgs = 0.0 ;
        Pzadd_dPS0Z = 0.0 ;
	
	if ( Pzadd < epsm10 ) {
	  Pzadd = epsm10 ;
	  Pzadd_dVbs = 0.0 ;
	  Pzadd_dVds = 0.0 ;
	  Pzadd_dVgs = 0.0 ;
          Pzadd_dPS0Z = 0.0 ;
	}
	
	Ps0z = Ps0 + Pzadd ;
	Ps0z_dVbs = Ps0_dVbs + Pzadd_dVbs ;
	Ps0z_dVds = Ps0_dVds + Pzadd_dVds ;
	Ps0z_dVgs = Ps0_dVgs + Pzadd_dVgs ;
        Ps0z_dPS0Z = Ps0_dPS0Z + Pzadd_dPS0Z ;

        /* calculate Newton correction: */
        G = PS0Z_SCE - Ps0z ;
        delta_PS0Z_SCE = - G / (1.0 - Ps0z_dPS0Z) ;

        delta_PS0Z_SCE_dVbs = Ps0z_dVbs - PS0Z_SCE_dVbs;
        delta_PS0Z_SCE_dVds = Ps0z_dVds - PS0Z_SCE_dVds;
        delta_PS0Z_SCE_dVgs = Ps0z_dVgs - PS0Z_SCE_dVgs;
        PS0Z_SCE += delta_PS0Z_SCE ;
	PS0Z_SCE_dVbs = Ps0z_dVbs ;
        PS0Z_SCE_dVds = Ps0z_dVds ;
        PS0Z_SCE_dVgs = Ps0z_dVgs ;
/*
        PS0_SCE = PS0Z_SCE - Pzadd ;
        PS0_SCE_dVbs = Ps0_dVbs ;
        PS0_SCE_dVds = Ps0_dVds ;
        PS0_SCE_dVgs = Ps0_dVgs ;
*/
	NNN += 1 ;

	if( (    fabs(delta_PS0Z_SCE) > PS0_SCE_tol
              || fabs(delta_PS0Z_SCE_dVbs) > PS0_SCE_deriv_tol
              || fabs(delta_PS0Z_SCE_dVds) > PS0_SCE_deriv_tol
              || fabs(delta_PS0Z_SCE_dVgs) > PS0_SCE_deriv_tol                        
            ) && (NNN < MAX_LOOP_SCE)  ){
	  goto START_OF_SCE_LOOP; 
	}
      }

      goto end_of_part_1 ;
    } 

    
    /*-----------------------------------------------------------*
     * Initial guess for Ps0. 
     *-----------------*/

    /*---------------------------------------------------*
     * Ps0_iniA: solution of subthreshold equation assuming zone-D1/D2.
     *-----------------*/
    TX = 1.0e0 + 4.0e0 
      * ( beta * ( Vgp - Vbs ) - 1.0e0 ) / ( fac1p2 * beta2 ) ;
    TX = Fn_Max( TX , epsm10 ) ;
    Ps0_iniA = Vgp + fac1p2 * beta * 0.5 * ( 1.0e0 - sqrt( TX ) ) ;

    /* use analitical value in subthreshold region. */
    if ( Vgs < ( Vfb + Vth ) * 0.5 ) {
        flg_pprv = 0 ;
    }

    if ( flg_pprv >= 1 ) {
      /*---------------------------------------------------*
       * Use previous value.
       *-----------------*/

      T1  = Ps0_dVbs * dVbs + Ps0_dVds * dVds  + Ps0_dVgs * dVgs ;
      Ps0_ini  = Ps0 + T1 ;

      if ( flg_pprv == 2 ) {
        /* TX_dVxs = d^2 Ps0 / d Vxs^2 here */
        if ( Vbsc_dif2 > epsm10 ) {
          TX_dVbs = ( here->HSM2_ps0_dvbs_prv - here->HSM2_ps0_dvbs_prv2 )
                  / Vbsc_dif2 ;
        } else {
          TX_dVbs = 0.0 ;
        }
        if ( Vdsc_dif2 > epsm10 ) {
          TX_dVds = ( here->HSM2_ps0_dvds_prv - here->HSM2_ps0_dvds_prv2 )
                  / Vdsc_dif2 ;
        } else {
          TX_dVds = 0.0 ;
        }
        if ( Vgsc_dif2 > epsm10 ) {
          TX_dVgs = ( here->HSM2_ps0_dvgs_prv - here->HSM2_ps0_dvgs_prv2 )
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
        flg_pprv = 0 ; /* flag changes to analitical */
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
      Chi = beta * ( Ps0_iniA - Vbs ) ;
 
      if ( Chi < znbd3 ) { 
        /*-----------------------------------*
         * zone-D1/D2
         * - Ps0_ini is the analytical solution of Qs=Qb0 with
         *   Qb0 being approximated to 3-degree polynomial.
         *-----------------*/
        TY = beta * ( Vgp - Vbs ) ;
        T1 = 1.0e0 / ( cn_nc3 * beta * fac1 ) ;
        T2 = 81.0 + 3.0 * T1 ;
        T3 = -2916.0 - 81.0 * T1 + 27.0 * T1 * TY ;
        T4 = 1458.0 - 81.0 * ( 54.0 + T1 ) + 27.0 * T1 * TY ;
        T4 = T4 * T4 ;
        T5 = Fn_Pow( T3 + sqrt( 4 * T2 * T2 * T2 + T4 ) , C_1o3 ) ;
        TX = 3.0 - ( C_2p_1o3 * T2 ) / ( 3.0 * T5 )
           + 1 / ( 3.0 * C_2p_1o3 ) * T5 ;
        
        Ps0_iniA = TX * beta_inv + Vbs ;
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
        
        Ps0_iniB = log( T2 + small ) / T3 ;
        
        Fn_SU( Ps0_ini , Ps0_iniA, Ps0_iniB, c_ps0ini_2, T1) ;
      } 
    }
  
    TX = Vbs + ps_conv / 2 ;
    if ( Ps0_ini < TX ) Ps0_ini = TX ;


    /*---------------------------------------------------*
     * Assign initial guess.
     *-----------------*/
    if ( NNN == 0 ) {
       Ps0 = Ps0_ini ;
    }
    Psl_lim = Ps0_iniA ;

    /*---------------------------------------------------*
     * Calculation of Ps0. (beginning of Newton loop) 
     * - Fs0 : Fs0 = 0 is the equation to be solved. 
     * - dPs0 : correction value. 
     *-----------------*/
    exp_bVbs = exp( beta * Vbs ) ;
    cfs1 = cnst1 * exp_bVbs ;
    
    flg_conv = 0 ;
    for ( lp_s0 = 1 ; lp_s0 <= lp_s0_max + 1 ; lp_s0 ++ ) { 
      
      Chi = beta * ( Ps0 - Vbs ) ;
    
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
            fs01 = cnst1 * ( exp_bPs0 - exp_bVbs ) ;
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
       * - clamped to Vbs if Ps0 < Vbs .
       *-----------------*/
      dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Ps0))) ;
      if ( fabs( dPs0 ) > dPlim ) dPs0 = dPlim * Fn_Sgn( dPs0 ) ;

      Ps0 = Ps0 + dPs0 ;

      TX = Vbs + ps_conv / 2 ;
      if ( Ps0 < TX ) Ps0 = TX ;
      
      /*-------------------------------------------*
       * Check convergence. 
       * NOTE: This condition may be too rigid. 
       *-----------------*/
      if ( fabs( dPs0 ) <= ps_conv && fabs( Fs0 ) <= gs_conv ) {
        flg_conv = 1 ;
      }
      
    } /* end of Ps0 Newton loop */

    /* Eliminate loop count to take account of the derivative loop */
    lp_s0 -- ;

    /*-------------------------------------------*
     * Procedure for diverged case.
     *-----------------*/
    if ( flg_conv == 0 ) { 
      fprintf( stderr , 
               "*** warning(HiSIM): Went Over Iteration Maximum (Ps0)\n" ) ;
      fprintf( stderr , 
               " Vbse   = %7.3f Vdse = %7.3f Vgse = %7.3f\n" , 
               Vbse , Vdse , Vgse ) ;
      if ( flg_info >= 2 ) {
        printf( "*** warning(HiSIM): Went Over Iteration Maximum (Ps0)\n" ) ;
      }
    } 

    /*---------------------------------------------------*
     * Evaluate derivatives of  Ps0. 
     * - note: Here, fs01_dVbs and fs02_dVbs are derivatives 
     *   w.r.t. explicit Vbs. So, Ps0 in the fs01 and fs02
     *   expressions is regarded as a constant.
     *-----------------*/
  
    /* derivatives of fs0* w.r.t. explicit Vbs */
    if ( Chi < znbd5 ) { 
      fs01_dVbs = cfs1 * beta * fi * ( fi - 2 * fi_dChi ) ;
      T2 = 1.0e0 / ( fs02 + fs02 ) ;
      fs02_dVbs = ( - beta * fb_dChi * 2 * fb + fs01_dVbs ) * T2 ;
    } else {
      fs01_dVbs = - cfs1 * beta ;
      T2 = 0.5e0 / fs02 ;
      fs02_dVbs = ( - beta + fs01_dVbs ) * T2 ;
    }

    T1 = 1.0 / Fs0_dPs0 ;
    Ps0_dVbs = - ( Vgp_dVbs - ( fac1 * fs02_dVbs + fac1_dVbs * fs02 ) ) * T1 ;
    Ps0_dVds = - ( Vgp_dVds - fac1_dVds * fs02 ) * T1 ;
    Ps0_dVgs = - ( Vgp_dVgs - fac1_dVgs * fs02 ) * T1 ;

    if ( Chi < znbd5 ) { 
      /*-------------------------------------------*
       * zone-D1/D2. (Ps0)
       * Xi0 := fdep0^2 = fb * fb  [D1,D2]
       *-----------------*/
      Xi0 = fb * fb + epsm10 ;
      T1 = 2 * fb * fb_dChi * beta ;
      Xi0_dVbs = T1 * ( Ps0_dVbs - 1.0 ) ;
      Xi0_dVds = T1 * Ps0_dVds ;
      Xi0_dVgs = T1 * Ps0_dVgs ;

      Xi0p12 = fb + epsm10 ;
      T1 = fb_dChi * beta ;
      Xi0p12_dVbs = T1 * ( Ps0_dVbs - 1.0 ) ;
      Xi0p12_dVds = T1 * Ps0_dVds ;
      Xi0p12_dVgs = T1 * Ps0_dVgs ;

      Xi0p32 = fb * fb * fb + epsm10 ;
      T1 = 3 * fb * fb * fb_dChi * beta ;
      Xi0p32_dVbs = T1 * ( Ps0_dVbs - 1.0 ) ;
      Xi0p32_dVds = T1 * Ps0_dVds ;
      Xi0p32_dVgs = T1 * Ps0_dVgs ;
      
    } else { 
      /*-------------------------------------------*
       * zone-D3. (Ps0)
       *-----------------*/
      flg_zone = 3 ;
      flg_noqi = 0 ;

      /*-----------------------------------*
       * Xi0 := fdep0^2 = Chi - 1 = beta * ( Ps0 - Vbs ) - 1 [D3]
       *-----------------*/
      Xi0 = Chi - 1.0e0 ;
      Xi0_dVbs = beta * ( Ps0_dVbs - 1.0e0 ) ;
      Xi0_dVds = beta * Ps0_dVds ;
      Xi0_dVgs = beta * Ps0_dVgs ;
 
      Xi0p12 = sqrt( Xi0 ) ;
      T1 = 0.5e0 / Xi0p12 ;
      Xi0p12_dVbs = T1 * Xi0_dVbs ;
      Xi0p12_dVds = T1 * Xi0_dVds ;
      Xi0p12_dVgs = T1 * Xi0_dVgs ;
 
      Xi0p32 = Xi0 * Xi0p12 ;
      T1 = 1.5e0 * Xi0p12 ;
      Xi0p32_dVbs = T1 * Xi0_dVbs ;
      Xi0p32_dVds = T1 * Xi0_dVds ;
      Xi0p32_dVgs = T1 * Xi0_dVgs ;

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

    T1 = 1.0 / ( fs02 + Xi0p12 ) ;
    Qn0 = cnst0 * fs01 * T1 ;

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

        Qbu = Qb0 ;
        Qbu_dVbs = Qb0_dVb ;
        Qbu_dVds = Qb0_dVd ;
        Qbu_dVgs = Qb0_dVg ;

        Qdrat = 0.5 ;
        Qdrat_dVbs = 0.0 ;
        Qdrat_dVds = 0.0 ;
        Qdrat_dVgs = 0.0 ;
 
        Lred = 0.0e0 ;
        Lred_dVbs = 0.0e0 ;
        Lred_dVds = 0.0e0 ;
        Lred_dVgs = 0.0e0 ;

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
        TX_dVbs = beta * T1 * ( Ps0_dVbs - 1.0 ) ;
        TX_dVds = beta * T1 * Ps0_dVds ;
        TX_dVgs = beta * T1 * Ps0_dVgs ;

        FD2 = TX * TX * TX * ( 10.0 + TX * ( -15.0 + TX * 6.0 ) ) ;
        T4 = TX * TX * ( 30.0 + TX * ( -60.0 + TX * 30.0 ) ) ;

        FD2_dVbs = T4 * TX_dVbs ;
        FD2_dVds = T4 * TX_dVds ;
        FD2_dVgs = T4 * TX_dVgs ;
      } /* end of zone-D2 */
    }


    /*---------------------------------------------------*
     * VgVt : Vgp - Vth_qi. ( Vth_qi is Vth for Qi evaluation. ) 
     *-----------------*/
    VgVt = Qn0 * Cox_inv ;
    VgVt_dVbs = Qn0_dVbs * Cox_inv + Qn0 * Cox_inv_dVb ;
    VgVt_dVds = Qn0_dVds * Cox_inv + Qn0 * Cox_inv_dVd ;
    VgVt_dVgs = Qn0_dVgs * Cox_inv + Qn0 * Cox_inv_dVg ;

    /*-----------------------------------------------------------*
     * make Qi=Qd=Ids=0 if VgVt <= VgVt_small 
     *-----------------*/
    if ( VgVt <= VgVt_small && model->HSM2_bypass_enable ) {
      flg_zone = 4 ;
      flg_noqi = 1 ;
        
      Psl = Ps0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVds = Ps0_dVds ;
      Psl_dVgs = Ps0_dVgs ;

      
      Pds = 0.0 ;
      Pds_dVbs = 0.0 ;
      Pds_dVds = 0.0 ;
      Pds_dVgs = 0.0 ;

      Qbu = Qb0 ;
      Qbu_dVbs = Qb0_dVb ;
      Qbu_dVds = Qb0_dVd ;
      Qbu_dVgs = Qb0_dVg ;

      Qiu = 0.0 ;
      Qiu_dVbs = 0.0 ;
      Qiu_dVds = 0.0 ;
      Qiu_dVgs = 0.0 ;
    
      Qdrat = 0.5 ;
      Qdrat_dVbs = 0.0 ;
      Qdrat_dVds = 0.0 ;
      Qdrat_dVgs = 0.0 ;

      Lred = 0.0 ;
      Lred_dVbs = 0.0 ;
      Lred_dVds = 0.0 ;
      Lred_dVgs = 0.0 ;
    
      Ids = 0.0e0 ;
      Ids_dVbs = 0.0e0 ;
      Ids_dVds = 0.0e0 ;
      Ids_dVgs = 0.0e0 ;
      
      DJI = 1.0 ;
            

      if( corecip ){

        Fs0_dPS0Z = Vgp_dPS0Z ;
        Ps0_dPS0Z = - Fs0_dPS0Z / Fs0_dPs0 ;

	T1 = Vds * 0.5  ;
	Fn_SymAdd( Pzadd , T1 , model->HSM2_pzadd0 , T2 ) ;
	T2 /= 2 ;
	Pzadd_dVbs = 0.0 ;
	Pzadd_dVds = T2 * ( 1.0 ) ;
	Pzadd_dVgs = 0.0 ;
        Pzadd_dPS0Z = 0.0 ;
	if ( Pzadd < epsm10 ) {
	  Pzadd = epsm10 ;
	  Pzadd_dVbs = 0.0 ;
	  Pzadd_dVds = 0.0 ;
	  Pzadd_dVgs = 0.0 ;
          Pzadd_dPS0Z = 0.0 ;
	}
	
	Ps0z = Ps0 + Pzadd ;
	Ps0z_dVbs = Ps0_dVbs + Pzadd_dVbs ;
	Ps0z_dVds = Ps0_dVds + Pzadd_dVds ;
	Ps0z_dVgs = Ps0_dVgs + Pzadd_dVgs ;
        Ps0z_dPS0Z = Ps0_dPS0Z + Pzadd_dPS0Z ;

        /* calculate Newton correction: */
        G = PS0Z_SCE - Ps0z ;
        delta_PS0Z_SCE = - G / (1.0 - Ps0z_dPS0Z) ;
  
        delta_PS0Z_SCE_dVbs = Ps0z_dVbs - PS0Z_SCE_dVbs;
        delta_PS0Z_SCE_dVds = Ps0z_dVds - PS0Z_SCE_dVds;
        delta_PS0Z_SCE_dVgs = Ps0z_dVgs - PS0Z_SCE_dVgs;
        PS0Z_SCE += delta_PS0Z_SCE ;
	PS0Z_SCE_dVbs = Ps0z_dVbs ;
        PS0Z_SCE_dVds = Ps0z_dVds ;
        PS0Z_SCE_dVgs = Ps0z_dVgs ;
/*
        PS0_SCE = PS0Z_SCE - Pzadd ;
        PS0_SCE_dVbs = Ps0_dVbs ;
        PS0_SCE_dVds = Ps0_dVds ;
        PS0_SCE_dVgs = Ps0_dVgs ;
*/

	NNN += 1 ;

	if( (    fabs(delta_PS0Z_SCE) > PS0_SCE_tol
              || fabs(delta_PS0Z_SCE_dVbs) > PS0_SCE_deriv_tol
              || fabs(delta_PS0Z_SCE_dVds) > PS0_SCE_deriv_tol
              || fabs(delta_PS0Z_SCE_dVgs) > PS0_SCE_deriv_tol                        
            ) && (NNN < MAX_LOOP_SCE)  ){
	  goto START_OF_SCE_LOOP; 
	}


      }

      goto end_of_part_1 ;
    }


    /*-----------------------------------------------------------*
     * Start point of Psl (= Ps0 + Pds) calculation. (label)
     *-----------------*/
/*  start_of_Psl:*/


    /* Vdseff (begin) */
    Vdsorg = Vds ;

    T2 = here->HSM2_qnsub_esi / ( Cox * Cox ) ;
    T4 = - 2.0e0 * T2 / Cox ;
    T2_dVb = T4 * Cox_dVb ;
    T2_dVd = T4 * Cox_dVd ;
    T2_dVg = T4 * Cox_dVg ;

    T5     = Vgpz - beta_inv - Vbsz ;
    T5_dVb = Vgpz_dVbs       - Vbsz_dVbs ;
    T5_dVd = Vgpz_dVds       - Vbsz_dVds ;
    T5_dVg = Vgpz_dVgs ;

    T0 = 2.0 / T2 ;
    T1     = 1.0e0 + T0 * T5 ;
    T1_dVb = ( T5_dVb * T2 - T5 * T2_dVb ) * T0 / T2 ;
    T1_dVd = ( T5_dVd * T2 - T5 * T2_dVd ) * T0 / T2 ;
    T1_dVg = ( T5_dVg * T2 - T5 * T2_dVg ) * T0 / T2 ;

    Fn_SZ( T1 , T1 , 0.05 , T9 ) ;
    T1 += small ;
    T1_dVb *= T9 ;
    T1_dVd *= T9 ;
    T1_dVg *= T9 ;

    T3 = sqrt( T1 ) ; 
    T3_dVb = 0.5 / T3 * T1_dVb ;
    T3_dVd = 0.5 / T3 * T1_dVd ;
    T3_dVg = 0.5 / T3 * T1_dVg ;

    T10 = Vgpz + T2 * ( 1.0e0 - T3 ) ; 
    T10_dVb = Vgpz_dVbs + T2_dVb * ( 1.0e0 - T3 ) - T2 * T3_dVb ;
    T10_dVd = Vgpz_dVds + T2_dVd * ( 1.0e0 - T3 ) - T2 * T3_dVd ;
    T10_dVg = Vgpz_dVgs + T2_dVg * ( 1.0e0 - T3 ) - T2 * T3_dVg ;
    Fn_SZ( T10 , T10 , 0.01 , T0 ) ;
    T10 += epsm10 ;
    T10_dVb *= T0 ;
    T10_dVd *= T0 ;
    T10_dVg *= T0 ;


    T1 = Vds / T10 ;
    T2 = Fn_Pow( T1 , here->HSM2_ddlt - 1.0e0 ) ;
    T7 = T2 * T1 ;
    T0 = here->HSM2_ddlt * T2 / ( T10 * T10 ) ;
    T7_dVb = T0 * ( - Vds * T10_dVb ) ;
    T7_dVd = T0 * ( T10 - Vds * T10_dVd ) ;
    T7_dVg = T0 * ( - Vds * T10_dVg ) ;

    T3 = 1.0 + T7 ;
    T4 = Fn_Pow( T3 , 1.0 / here->HSM2_ddlt - 1.0 ) ;
    T6 = T4 * T3 ;
    T0 = T4 / here->HSM2_ddlt ;
    T6_dVb = T0 * T7_dVb ;
    T6_dVd = T0 * T7_dVd ;
    T6_dVg = T0 * T7_dVg ;

    Vdseff = Vds / T6 ;
    T0 = 1.0 / ( T6 * T6 ) ;
    Vdseff_dVbs =      - Vds * T6_dVb   * T0 ;
    Vdseff_dVds = ( T6 - Vds * T6_dVd ) * T0 ;
    Vdseff_dVgs =      - Vds * T6_dVg   * T0 ;

    Vds = Vdseff ;
    /* Vdseff (end) */


    exp_bVbsVds = exp( beta * ( Vbs - Vds ) ) ;
  
  
    /*---------------------------------------------------*
     * Skip Psl calculation when Vds is very small.
     *-----------------*/
    if ( Vds <= 0.0 ) {
      Pds = 0.0 ;
      Psl = Ps0 ;
      flg_conv = 1 ;
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

      if ( flg_pprv == 2 ) {
        /* TX_dVxs = d^2 Pds / d Vxs^2 here */
        if ( Vbsc_dif2 > epsm10 ) {
          TX_dVbs = ( here->HSM2_pds_dvbs_prv - here->HSM2_pds_dvbs_prv2 )
                  / Vbsc_dif2 ;
        } else {
          TX_dVbs = 0.0 ;
        }
        if ( Vdsc_dif2 > epsm10 ) {
          TX_dVds = ( here->HSM2_pds_dvds_prv - here->HSM2_pds_dvds_prv2 )
                  / Vdsc_dif2 ;
        } else {
          TX_dVds = 0.0 ;
        }
        if ( Vgsc_dif2 > epsm10 ) {
          TX_dVgs = ( here->HSM2_pds_dvgs_prv - here->HSM2_pds_dvgs_prv2 )
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
    if ( NNN == 0 ) {
    Pds = Pds_ini ;
    Psl = Ps0 + Pds ;
    } else { 
      /* take solution from previous PS0_SCE_loop as initial value */
      Pds = Psl - Ps0 ;
    }
    TX = Vbs + ps_conv / 2 ;
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
      
      Chi = beta * ( Psl - Vbs ) ;
    
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
       * - clamped to Vbs if Psl < Vbs .
       *-----------------*/
      dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Psl))) ;
      if ( fabs( dPsl ) > dPlim ) dPsl = dPlim * Fn_Sgn( dPsl ) ;

      Psl = Psl + dPsl ;

      TX = Vbs + ps_conv / 2 ;
      if ( Psl < TX ) Psl = TX ;

      /*-------------------------------------------*
       * Check convergence.
       * NOTE: This condition may be too rigid.
       *-----------------*/
      if ( fabs( dPsl ) <= ps_conv && fabs( Fsl ) <= gs_conv ) {
        flg_conv = 1 ;
      }
    } /* end of Psl Newton loop */

    /* Eliminate loop count to take account of the derivative loop */
    lp_sl -- ;

    /*-------------------------------------------*
     * Procedure for diverged case.
     *-----------------*/
    if ( flg_conv == 0 ) {
      fprintf( stderr ,
               "*** warning(HiSIM): Went Over Iteration Maximum (Psl)\n" ) ;
      fprintf( stderr ,
               " Vbse   = %7.3f Vdse = %7.3f Vgse = %7.3f\n" ,
               Vbse , Vdse , Vgse ) ;
      if ( flg_info >= 2 ) {
        printf("*** warning(HiSIM): Went Over Iteration Maximum (Psl)\n" ) ;
      }
    }


    /*---------------------------------------------------*
     * Evaluate derivatives of  Psl.
     * - note: Here, fsl1_dVbs and fsl2_dVbs are derivatives 
     *   w.r.t. explicit Vbs. So, Psl in the fsl1 and fsl2
     *   expressions is regarded as a constant.
     *-----------------*/
    if ( Chi < znbd5 ) { 
      T1 = cfs1 * beta * fi ;
      fsl1_dVbs = T1 * ( ( 1.0 - Vdseff_dVbs ) * fi - 2.0 * fi_dChi ) ;
      fsl1_dVds = - T1 * fi * Vdseff_dVds ;
      fsl1_dVgs = - T1 * fi * Vdseff_dVgs ;
      T2 =  0.5 / fsl2 ;
      fsl2_dVbs = ( - beta * fb_dChi * 2 * fb + fsl1_dVbs ) * T2 ;
      fsl2_dVds = fsl1_dVds * T2 ;
      fsl2_dVgs = fsl1_dVgs * T2 ; 
    } else {
      T1 = cnst1 * beta ;
      fsl1_dVbs = - T1 * ( exp_Rho * Vdseff_dVbs
                         + ( 1.0 - Vdseff_dVbs ) * exp_bVbsVds );
      fsl1_dVds = - T1 * Vdseff_dVds * ( exp_Rho - exp_bVbsVds );
      fsl1_dVgs =   T1 * Vdseff_dVgs * ( - exp_Rho + exp_bVbsVds );
      T2 = 0.5e0 / fsl2 ;
      fsl2_dVbs = ( - beta + fsl1_dVbs ) * T2 ;
      fsl2_dVds = ( fsl1_dVds ) * T2 ;
      fsl2_dVgs = ( fsl1_dVgs ) * T2 ;
    }

    T1 = 1.0 / Fsl_dPsl ;
    Psl_dVbs = - ( Vgp_dVbs - ( fac1 * fsl2_dVbs + fac1_dVbs * fsl2 ) ) * T1 ;
    Psl_dVds = - ( Vgp_dVds - ( fac1 * fsl2_dVds + fac1_dVds * fsl2 ) ) * T1 ;
    Psl_dVgs = - ( Vgp_dVgs - ( fac1 * fsl2_dVgs + fac1_dVgs * fsl2 ) ) * T1 ;

    if ( Chi < znbd5 ) { 
      /*-------------------------------------------*
       * zone-D1/D2. (Psl)
       *-----------------*/
      Xil = fb * fb + epsm10 ;
      T1 = 2 * fb * fb_dChi * beta ;
      Xil_dVbs = T1 * ( Psl_dVbs - 1.0 ) ;
      Xil_dVds = T1 * Psl_dVds ;
      Xil_dVgs = T1 * Psl_dVgs ;

      Xilp12 = fb + epsm10 ;
      T1 = fb_dChi * beta ;
      Xilp12_dVbs = T1 * ( Psl_dVbs - 1.0 ) ;
      Xilp12_dVds = T1 * Psl_dVds ;
      Xilp12_dVgs = T1 * Psl_dVgs ;
    
      Xilp32 = fb * fb * fb + epsm10 ;
      T1 = 3 * fb * fb * fb_dChi * beta ;
      Xilp32_dVbs = T1 * ( Psl_dVbs - 1.0 ) ;
      Xilp32_dVds = T1 * Psl_dVds ;
      Xilp32_dVgs = T1 * Psl_dVgs ;
    
    } else { 
      /*-------------------------------------------*
       * zone-D3. (Psl)
       *-----------------*/

      Xil = Chi - 1.0e0 ;
      Xil_dVbs = beta * ( Psl_dVbs - 1.0e0 ) ;
      Xil_dVds = beta * Psl_dVds ;
      Xil_dVgs = beta * Psl_dVgs ;
 
      Xilp12 = sqrt( Xil ) ;
      T1 = 0.5e0 / Xilp12 ;
      Xilp12_dVbs = T1 * Xil_dVbs ;
      Xilp12_dVds = T1 * Xil_dVds ;
      Xilp12_dVgs = T1 * Xil_dVgs ;

      Xilp32 = Xil * Xilp12 ;
      T1 = 1.5e0 * Xilp12 ;
      Xilp32_dVbs = T1 * Xil_dVbs ;
      Xilp32_dVds = T1 * Xil_dVds ;
      Xilp32_dVgs = T1 * Xil_dVgs ;

    }

    /*---------------------------------------------------*
     * Assign Pds.
     *-----------------*/
    Pds = Psl - Ps0 ;

    if ( Pds < 0.0 ) {
      Pds = 0.0 ;
      Psl = Ps0 ;
    }

    Pds_dVbs = Psl_dVbs - Ps0_dVbs ;
    Pds_dVds = Psl_dVds - Ps0_dVds ;
    Pds_dVgs = Psl_dVgs - Ps0_dVgs ;
  
    if ( Pds < ps_conv ) {
      Pds_dVbs = 0.0 ;
      Pds_dVgs = 0.0 ;
      Psl_dVbs = Ps0_dVbs ;
      Psl_dVgs = Ps0_dVgs ;
    }

    /* Vdseff */
    Vds = Vdsorg;


    if( corecip ){

      /* For multi level Newton method to calculate PS0Z_SCE        */
      /* remember that for the inner Newton:                        */                                        
      /* Fs0(PS0Z_SCE,Ps0) = Vgp(PSOZ_SCE) - Ps0 - fac1 * fs02(Ps0) */
      /* -> dPs0/dPS0Z_SCE = - (dFs0/dPS0Z_SCE) / (dFs0/dPs0)       */
      /* and                                                        */
      /* Fsl(PS0Z_SCE,Psl) = Vgp(PS0Z_SCE) - Psl - fac1 * fsl1(Psl) */
      /* -> dPsl/dPS0Z_SCE = - (dFsl/dPS0Z_SCE) / (dFsl/dPsl)       */

      /* Outer Newton:                                                                       */
      /* PS0Z_SCE = Ps0 + Pzadd(Ps0,Psl)                                                     */
      /* -> G(PS0Z_SCE) := PS0Z_SCE - Ps0(PS0Z_SCE) - Pzadd(Ps0(PS0Z_SCE),Psl(Ps0Z_SCE)) = 0 */
      /* -> Newton correction delta_PS0Z_SCE from:                                           */
      /*         (1.0 - dPs0/dPS0Z_SCE - dPzadd/dPS0Z_SCE) * delta_PS0Z_SCE = - G            */

      Fs0_dPS0Z = Vgp_dPS0Z ;
      Ps0_dPS0Z = - Fs0_dPS0Z / Fs0_dPs0 ;

      Fsl_dPS0Z = Vgp_dPS0Z ;
      Psl_dPS0Z = - Fsl_dPS0Z / Fsl_dPsl ;

      if ( Pds < ps_conv ) {
        Pds_dPS0Z = 0.0 ;
      } else {
        Pds_dPS0Z = Psl_dPS0Z - Ps0_dPS0Z ;
      }

      T1 =  ( Vds - Pds ) / 2 ;
      Fn_SymAdd( Pzadd , T1 , model->HSM2_pzadd0 , T2 ) ;
      T2 /= 2 ;
      Pzadd_dVbs = T2 * ( - Pds_dVbs ) ;
      Pzadd_dVds = T2 * ( 1.0 - Pds_dVds ) ;
      Pzadd_dVgs = T2 * ( - Pds_dVgs ) ;
      Pzadd_dPS0Z = T2 * ( - Pds_dPS0Z) ;
  
      if ( Pzadd < epsm10 ) {
	Pzadd = epsm10 ;
	Pzadd_dVbs = 0.0 ;
	Pzadd_dVds = 0.0 ;
	Pzadd_dVgs = 0.0 ;
        Pzadd_dPS0Z = 0.0 ;
      }

      Ps0z = Ps0 + Pzadd ;
      Ps0z_dVbs = Ps0_dVbs + Pzadd_dVbs ;
      Ps0z_dVds = Ps0_dVds + Pzadd_dVds ;
      Ps0z_dVgs = Ps0_dVgs + Pzadd_dVgs ;
      Ps0z_dPS0Z = Ps0_dPS0Z + Pzadd_dPS0Z ;


      /* calculate Newton correction: */
      G = PS0Z_SCE - Ps0z ;
      delta_PS0Z_SCE = - G / (1.0 - Ps0z_dPS0Z) ;

      delta_PS0Z_SCE_dVbs = Ps0z_dVbs - PS0Z_SCE_dVbs;
      delta_PS0Z_SCE_dVds = Ps0z_dVds - PS0Z_SCE_dVds;
      delta_PS0Z_SCE_dVgs = Ps0z_dVgs - PS0Z_SCE_dVgs;
      PS0Z_SCE += delta_PS0Z_SCE ;
      PS0Z_SCE_dVbs = Ps0z_dVbs ;
      PS0Z_SCE_dVds = Ps0z_dVds ;
      PS0Z_SCE_dVgs = Ps0z_dVgs ;
/*
      PS0_SCE = PS0Z_SCE - Pzadd ;
      PS0_SCE_dVbs = Ps0_dVbs ;
      PS0_SCE_dVds = Ps0_dVds ;
      PS0_SCE_dVgs = Ps0_dVgs ;
*/
      NNN += 1 ;

      if( (    fabs(delta_PS0Z_SCE) > PS0_SCE_tol
            || fabs(delta_PS0Z_SCE_dVbs) > PS0_SCE_deriv_tol
            || fabs(delta_PS0Z_SCE_dVds) > PS0_SCE_deriv_tol
            || fabs(delta_PS0Z_SCE_dVgs) > PS0_SCE_deriv_tol                        
          ) && (NNN < MAX_LOOP_SCE)  ){
        goto START_OF_SCE_LOOP;
      }


    }


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

    /* ( Eta + 1 )^n */
    Eta1 = Eta + 1.0e0 ;
    Eta1p12 = sqrt( Eta1 ) ;
    Eta1p32 = Eta1p12 * Eta1 ;
    Eta1p52 = Eta1p32 * Eta1 ;
 
    /* 1 / ( ( Eta + 1 )^n + 1 ) */
    Zeta12 = 1.0e0 / ( Eta1p12 + 1.0e0 ) ;
    Zeta32 = 1.0e0 / ( Eta1p32 + 1.0e0 ) ;
    Zeta52 = 1.0e0 / ( Eta1p52 + 1.0e0 ) ;

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
  
      
    /*---------------------------------------------------*
     *  Idd: 
     *-----------------*/
    Idd = Pds * Fdd ;
    Idd_dVbs = Pds_dVbs * Fdd + Pds * Fdd_dVbs ;
    Idd_dVds = Pds_dVds * Fdd + Pds * Fdd_dVds ;
    Idd_dVgs = Pds_dVgs * Fdd + Pds * Fdd_dVgs ;

    /*-----------------------------------------------------------*
     * Skip CLM and integrated charges if zone==D1
     *-----------------*/
    if( flg_zone == 1 ) {
      goto start_of_mobility ;
    }

    /*-----------------------------------------------------------*
     * Channel Length Modulation. Lred: \Delta L
     *-----------------*/
    if( pParam->HSM2_clm2 < epsm10 && pParam->HSM2_clm3 < epsm10 ) {
      Lred = 0.0e0 ;
      Lred_dVbs = 0.0e0 ;
      Lred_dVds = 0.0e0 ;
      Lred_dVgs = 0.0e0 ;

      Psdl = Psl ;
      Psdl_dVbs = Psl_dVbs ;
      Psdl_dVds = Psl_dVds ;
      Psdl_dVgs = Psl_dVgs ;
      if ( Psdl > Ps0 + Vds - epsm10 ) {
        Psdl = Ps0 + Vds - epsm10 ;
        Psdl_dVbs = Ps0_dVbs ;
        Psdl_dVds = Ps0_dVds + 1.0 ;
        Psdl_dVgs = Ps0_dVgs ;
      }

    } else {
      T1 = here->HSM2_wdpl ;
      T8 = sqrt (Psl - Vbs) ;
      Wd = T1 * T8 ;
      T9 = 0.5 * T1 / T8 ;
      Wd_dVbs = T9 * (Psl_dVbs - 1.0) ;
      Wd_dVds = T9 * Psl_dVds ;
      Wd_dVgs = T9 * Psl_dVgs ;

      T0 = 1.0 / Wd ;
      T1 = Qn0 * T0 ;
      T2 = pParam->HSM2_clm3 * T1 ;
      T3 = pParam->HSM2_clm3 * T0 ;
      T2_dVb = T3 * (Qn0_dVbs - T1 * Wd_dVbs) ;
      T2_dVd = T3 * (Qn0_dVds - T1 * Wd_dVds) ;
      T2_dVg = T3 * (Qn0_dVgs - T1 * Wd_dVgs) ;

      T5 = pParam->HSM2_clm2 * q_Nsub + T2 ;
      T1 = 1.0 / T5 ;
      T4 = C_ESI * T1 ;
      T4_dVb = - T4 * T2_dVb * T1 ;
      T4_dVd = - T4 * T2_dVd * T1 ;
      T4_dVg = - T4 * T2_dVg * T1 ;

      T1 = (1.0e0 - pParam->HSM2_clm1) ;
      Psdl = pParam->HSM2_clm1 * (Vds + Ps0) + T1 * Psl ;
      Psdl_dVbs = pParam->HSM2_clm1 * Ps0_dVbs + T1 * Psl_dVbs ;
      Psdl_dVds = pParam->HSM2_clm1 * (1.0 + Ps0_dVds) + T1 * Psl_dVds ;
      Psdl_dVgs = pParam->HSM2_clm1 * Ps0_dVgs + T1 * Psl_dVgs ;
      if ( Psdl > Ps0 + Vds - epsm10 ) {
        Psdl = Ps0 + Vds - epsm10 ;
        Psdl_dVbs = Ps0_dVbs ;
        Psdl_dVds = Ps0_dVds + 1.0 ;
        Psdl_dVgs = Ps0_dVgs ;
      }
      T6 = Psdl - Psl ; 
      T6_dVb = Psdl_dVbs - Psl_dVbs ;
      T6_dVd = Psdl_dVds - Psl_dVds ;
      T6_dVg = Psdl_dVgs - Psl_dVgs ;

      T3 = beta * Qn0 ;
      T1 = 1.0 / T3 ;
      T5 = Idd * T1 ;
      T2 = T5 * beta ;
      T5_dVb = (Idd_dVbs - T2 * Qn0_dVbs) * T1 ;
      T5_dVd = (Idd_dVds - T2 * Qn0_dVds) * T1 ;
      T5_dVg = (Idd_dVgs - T2 * Qn0_dVgs) * T1 ;
      
      T10 = q_Nsub / C_ESI ;
      T1 = C_E0_p2 ;  // E0^2
      T2 = 1.0 / Leff ;
      T11 = (2.0 * T5 + 2.0 * T10 * T6 * T4 + T1 * T4) * T2 ;
      T3 = T2 * T4 ;
      T7 = T11 * T4 ;
      T7_dVb = (2.0 * T5_dVb + 2.0 * T10 * (T6_dVb * T4 + T6 * T4_dVb) + T1 * T4_dVb) * T3 + T11 * T4_dVb ;
      T7_dVd = (2.0 * T5_dVd + 2.0 * T10 * (T6_dVd * T4 + T6 * T4_dVd) + T1 * T4_dVd) * T3 + T11 * T4_dVd ;
      T7_dVg = (2.0 * T5_dVg + 2.0 * T10 * (T6_dVg * T4 + T6 * T4_dVg) + T1 * T4_dVg) * T3 + T11 * T4_dVg ;

      T11 = 4.0 * (2.0 * T10 * T6 + T1) ;
      T1 = 8.0 * T10 * T4 * T4 ;
      T2 = 2.0 * T11 * T4 ;
      T8 = T11 * T4 * T4 ;
      T8_dVb = ( T1 * T6_dVb + T2 * T4_dVb) ;
      T8_dVd = ( T1 * T6_dVd + T2 * T4_dVd) ;
      T8_dVg = ( T1 * T6_dVg + T2 * T4_dVg) ;

      T9 = sqrt (T7 * T7 + T8);
      T1 = 1.0 / T9 ;
      T2 = T7 * T1 ;
      T3 = 0.5 * T1 ;
      T9_dVb = (T2 * T7_dVb + T3 * T8_dVb) ;
      T9_dVd = (T2 * T7_dVd + T3 * T8_dVd) ;
      T9_dVg = (T2 * T7_dVg + T3 * T8_dVg) ;
      
      Lred = 0.5 * (- T7 + T9) ;
      Lred_dVbs = 0.5 * (- T7_dVb + T9_dVb) ;
      Lred_dVds = 0.5 * (- T7_dVd + T9_dVd) ;
      Lred_dVgs = 0.5 * (- T7_dVg + T9_dVg) ;

      /*---------------------------------------------------*
       * Modify Lred for symmetry.
       *-----------------*/
      T1 = Lred ;
      Lred = FMDVDS * T1 ;
      Lred_dVbs = FMDVDS_dVbs * T1 + FMDVDS * Lred_dVbs ;
      Lred_dVds = FMDVDS_dVds * T1 + FMDVDS * Lred_dVds ;
      Lred_dVgs = FMDVDS_dVgs * T1 + FMDVDS * Lred_dVgs ;
    }

    /* CLM5 & CLM6 */
    Lred *= here->HSM2_clmmod ;
    Lred_dVbs *= here->HSM2_clmmod ;
    Lred_dVds *= here->HSM2_clmmod ;
    Lred_dVgs *= here->HSM2_clmmod ;


    /*---------------------------------------------------*
     * Qbu : -Qb in unit area.
     *-----------------*/
    T1 = Vgp + beta_inv ;
    T2 = T1 * F10 - F11 ;
    
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

    T1 = beta ;
    Qbu = T1 * Qbnm / Fdd ;
    T2 = T1 / ( Fdd * Fdd ) ;
    Qbu_dVbs = T2 * ( Fdd * Qbnm_dVbs - Qbnm * Fdd_dVbs ) ;
    Qbu_dVds = T2 * ( Fdd * Qbnm_dVds - Qbnm * Fdd_dVds ) ;
    Qbu_dVgs = T2 * ( Fdd * Qbnm_dVgs - Qbnm * Fdd_dVgs ) ;
    
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

    Achi = Pds + DtPds ;
    Achi_dVbs = Pds_dVbs + DtPds_dVbs ;
    Achi_dVds = Pds_dVds + DtPds_dVds ;
    Achi_dVgs = Pds_dVgs + DtPds_dVgs ;
    
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


    /*-----------------------------------------------------------*
     * Qiu : -Qi in unit area.
     *-----------------*/

    Qinm = 1.0e0 + Alpha * ( 1.0e0 + Alpha ) ;
    T1 = 1.0e0 + Alpha + Alpha ;
    Qinm_dVbs = Alpha_dVbs * T1 ;
    Qinm_dVds = Alpha_dVds * T1 ;
    Qinm_dVgs = Alpha_dVgs * T1 ;
 
    Qidn = Fn_Max( 1.0e0 + Alpha , epsm10 ) ;
    Qidn_dVbs = Alpha_dVbs ;
    Qidn_dVds = Alpha_dVds ;
    Qidn_dVgs = Alpha_dVgs ;
    
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
  
    /*-----------------------------------------------------------*
     * Qdrat : Qd/Qi
     *-----------------*/
    Qdnm = 0.5e0 + Alpha ;
    Qdnm_dVbs = Alpha_dVbs ;
    Qdnm_dVds = Alpha_dVds ;
    Qdnm_dVgs = Alpha_dVgs ;
    
    Qddn = Qidn * Qinm ;
    Qddn_dVbs = Qidn_dVbs * Qinm + Qidn * Qinm_dVbs ;
    Qddn_dVds = Qidn_dVds * Qinm + Qidn * Qinm_dVds ;
    Qddn_dVgs = Qidn_dVgs * Qinm + Qidn * Qinm_dVgs ;
    
    Quot = 0.4e0 * Qdnm  / Qddn ;
    Qdrat = 0.6e0 - Quot ;
 
    if ( Qdrat <= 0.5e0 ) { 
      T1 = 1.0 / Qddn ;
      T2 = 1.0 / Qdnm ;
      Qdrat_dVbs = Quot * ( Qddn_dVbs * T1 - Qdnm_dVbs * T2 ) ;
      Qdrat_dVds = Quot * ( Qddn_dVds * T1 - Qdnm_dVds * T2 ) ;
      Qdrat_dVgs = Quot * ( Qddn_dVgs * T1 - Qdnm_dVgs * T2 ) ;
    } else { 
      Qdrat = 0.5e0 ;
      Qdrat_dVbs = 0.0e0 ;
      Qdrat_dVds = 0.0e0 ;
      Qdrat_dVgs = 0.0e0 ;
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

      if ( Qbu < 0.0 ) {
        Qbu = 0.0 ;
        Qbu_dVbs = 0.0 ;
        Qbu_dVds = 0.0 ;
        Qbu_dVgs = 0.0 ;
      }
      T1 = Qiu ;
      Qiu = FD2 * Qiu + ( 1.0 - FD2 ) * Qn0 ;
      Qiu_dVbs = FD2 * Qiu_dVbs + FD2_dVbs * T1 
        + ( 1.0 - FD2 ) * Qn0_dVbs - FD2_dVbs * Qn0 ;
      Qiu_dVds = FD2 * Qiu_dVds + FD2_dVds * T1 
        + ( 1.0 - FD2 ) * Qn0_dVds - FD2_dVds * Qn0 ;
      Qiu_dVgs = FD2 * Qiu_dVgs + FD2_dVgs * T1 
        + ( 1.0 - FD2 ) * Qn0_dVgs - FD2_dVgs * Qn0 ;
      if ( Qiu < 0.0 ) {
        Qiu = 0.0 ;
        Qiu_dVbs = 0.0 ;
        Qiu_dVds = 0.0 ;
        Qiu_dVgs = 0.0 ;
      }

      T1 = Qdrat ;
      Qdrat = FD2 * Qdrat + ( 1.0 - FD2 ) * 0.5e0 ;
      Qdrat_dVbs = FD2 * Qdrat_dVbs + FD2_dVbs * T1 - FD2_dVbs * 0.5e0 ;
      Qdrat_dVds = FD2 * Qdrat_dVds + FD2_dVds * T1 - FD2_dVds * 0.5e0 ;
      Qdrat_dVgs = FD2 * Qdrat_dVgs + FD2_dVgs * T1 - FD2_dVgs * 0.5e0 ;

      /* note: Lred=0 in zone-D1 */
      T1 = Lred ;
      Lred = FD2 * Lred ;
      Lred_dVbs = FD2 * Lred_dVbs + FD2_dVbs * T1 ;
      Lred_dVds = FD2 * Lred_dVds + FD2_dVds * T1 ;
      Lred_dVgs = FD2 * Lred_dVgs + FD2_dVgs * T1 ;
    } /* end of flg_zone==2 if-block */


start_of_mobility:

    Lch = Leff - Lred ;
    if ( Lch < 1.0e-9 ) {
         Lch = 1.0e-9 ; Lch_dVbs = Lch_dVds = Lch_dVgs = 0.0 ;
    } else { Lch_dVbs = - Lred_dVbs ; Lch_dVds = - Lred_dVds ; Lch_dVgs = - Lred_dVgs ; }


    /*-----------------------------------------------------------*
     * Modified potential for symmetry.
     *-----------------*/
    T1 =  ( Vds - Pds ) / 2 ;
    Fn_SymAdd( Pzadd , T1 , model->HSM2_pzadd0 , T2 ) ;
    T2 /= 2 ;
    Pzadd_dVbs = T2 * ( - Pds_dVbs ) ;
    Pzadd_dVds = T2 * ( 1.0 - Pds_dVds ) ;
    Pzadd_dVgs = T2 * ( - Pds_dVgs ) ;  
  
    if ( Pzadd < epsm10 ) {
      Pzadd = epsm10 ;
      Pzadd_dVbs = 0.0 ;
      Pzadd_dVds = 0.0 ;
      Pzadd_dVgs = 0.0 ;
    }
  
    Ps0z = Ps0 + Pzadd ;
    Ps0z_dVbs = Ps0_dVbs + Pzadd_dVbs ;
    Ps0z_dVds = Ps0_dVds + Pzadd_dVds ;
    Ps0z_dVgs = Ps0_dVgs + Pzadd_dVgs ;


    /*-----------------------------------------------------------*
     * Muun : universal mobility. (CGS unit)
     *-----------------*/

    T1 = here->HSM2_ndep_o_esi ;
    T2 = here->HSM2_ninv_o_esi ;

    T0 = model->HSM2_ninvd ;
    T4 = 1.0 + ( Psl - Ps0 ) * T0 ;
    T4_dVb = ( Psl_dVbs - Ps0_dVbs ) * T0 ;
    T4_dVd = ( Psl_dVds - Ps0_dVds ) * T0 ;
    T4_dVg = ( Psl_dVgs - Ps0_dVgs ) * T0 ;

    T5     = T1 * Qbu      + T2 * Qiu ;
    T5_dVb = T1 * Qbu_dVbs + T2 * Qiu_dVbs ;
    T5_dVd = T1 * Qbu_dVds + T2 * Qiu_dVds ;
    T5_dVg = T1 * Qbu_dVgs + T2 * Qiu_dVgs ;

    T3     = T5 / T4 ;
    T3_dVb = ( - T4_dVb * T5 + T4 * T5_dVb ) / T4 / T4 ;
    T3_dVd = ( - T4_dVd * T5 + T4 * T5_dVd ) / T4 / T4 ;
    T3_dVg = ( - T4_dVg * T5 + T4 * T5_dVg ) / T4 / T4 ;

    Eeff = T3 ;
    Eeff_dVbs = T3_dVb ;
    Eeff_dVds = T3_dVd ;
    Eeff_dVgs = T3_dVg ;

    T5 = Fn_Pow( Eeff , model->HSM2_mueph0 - 1.0e0 ) ;
    T8 = T5 * Eeff ;
    T7 = Fn_Pow( Eeff , here->HSM2_muesr - 1.0e0 ) ;
    T6 = T7 * Eeff ;


    T9 = C_QE * C_m2cm_p2 ;
    Rns = Qiu / T9 ;
    T1 = 1.0e0 / ( here->HSM2_muecb0 + here->HSM2_muecb1 * Rns / 1.0e11 ) 
      + here->HSM2_mphn0 * T8 + T6 / pParam->HSM2_muesr1 ;
    Muun = 1.0e0 / T1 ;
 
    T1 = 1.0e0 / ( T1 * T1 ) ;
    T2 = here->HSM2_muecb0 + here->HSM2_muecb1 * Rns / 1.0e11 ;
    T2 = 1.0e0 / ( T2 * T2 ) ;
    T3 = here->HSM2_mphn1 * T5 ;
    T4 = here->HSM2_muesr * T7 / pParam->HSM2_muesr1 ;
    T5 = - 1.0e-11 * here->HSM2_muecb1 / C_QE * T2 / C_m2cm_p2 ;
    Muun_dVbs = - ( T5 * Qiu_dVbs 
                    + Eeff_dVbs * T3 + Eeff_dVbs * T4 ) * T1 ;
    Muun_dVds = - ( T5 * Qiu_dVds  
                    + Eeff_dVds * T3 + Eeff_dVds * T4 ) * T1 ;
    Muun_dVgs = - ( T5 * Qiu_dVgs  
                    + Eeff_dVgs * T3 + Eeff_dVgs * T4 ) * T1 ;

    /*  Change to MKS unit */
    Muun      /= C_m2cm_p2 ;
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

    TY = Idd * T1 ;
    TY_dVbs = Idd_dVbs * T1 + Idd * T1_dVb ;
    TY_dVds = Idd_dVds * T1 + Idd * T1_dVd ;
    TY_dVgs = Idd_dVgs * T1 + Idd * T1_dVg ;

    T2 = 0.2 * Vmax /  Muun ;
    T3 = - T2 / Muun ;
    T2_dVb = T3 * Muun_dVbs ;
    T2_dVd = T3 * Muun_dVds ;
    T2_dVg = T3 * Muun_dVgs ;

    Ey = sqrt( TY * TY + T2 * T2 ) ;
    T4 = 1.0 / Ey ;
    Ey_dVbs = T4 * ( TY * TY_dVbs + T2 * T2_dVb ) ;
    Ey_dVds = T4 * ( TY * TY_dVds + T2 * T2_dVd ) ;
    Ey_dVgs = T4 * ( TY * TY_dVgs + T2 * T2_dVg ) ;

    Em = Muun * Ey ;
    Em_dVbs = Muun_dVbs * Ey + Muun * Ey_dVbs ;
    Em_dVds = Muun_dVds * Ey + Muun * Ey_dVds ;
    Em_dVgs = Muun_dVgs * Ey + Muun * Ey_dVgs ;
    
    T1  = Em / Vmax ;

    /* note: model->HSM2_bb = 2 (electron) ;1 (hole) */
    if ( 1.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 1.0e0 + epsm10 ) {
      T3 = 1.0e0 ;
    } else if ( 2.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 2.0e0 + epsm10 ) {
      T3 = T1 ;
    } else {
      T3 = Fn_Pow( T1 , model->HSM2_bb - 1.0e0 ) ;
    }
    T2 = T1 * T3 ;
    T4 = 1.0e0 + T2 ;

    if ( 1.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 1.0e0 + epsm10 ) {
      T5 = 1.0 / T4 ;
          T6 = T5 / T4 ;
    } else if ( 2.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 2.0e0 + epsm10 ) {
      T5 = 1.0 / sqrt( T4 ) ;
          T6 = T5 / T4 ;
    } else {
      T6 = Fn_Pow( T4 , ( - 1.0e0 / model->HSM2_bb - 1.0e0 ) ) ;
      T5 = T4 * T6 ;
    }

    T7 = Muun / Vmax * T6 * T3 ;

    Mu = Muun * T5 ;
    Mu_dVbs = Muun_dVbs * T5 - T7 * Em_dVbs ;
    Mu_dVds = Muun_dVds * T5 - T7 * Em_dVds ;
    Mu_dVgs = Muun_dVgs * T5 - T7 * Em_dVgs ;

/*  end_of_mobility : */

    /*-----------------------------------------------------------*
     * Ids: channel current.
     *-----------------*/
    betaWL = here->HSM2_weff_nf * beta_inv / Lch ;
    T1 = - betaWL / Lch ;
    betaWL_dVbs = T1 * Lch_dVbs ;
    betaWL_dVds = T1 * Lch_dVds ;
    betaWL_dVgs = T1 * Lch_dVgs ;

    Ids0 = betaWL * Idd * Mu ;
    T1 = betaWL * Idd ;
    T2 = Idd * Mu ;
    T3 = Mu * betaWL ;
    Ids0_dVbs = T3 * Idd_dVbs + T1 * Mu_dVbs + T2 * betaWL_dVbs ;
    Ids0_dVds = T3 * Idd_dVds + T1 * Mu_dVds + T2 * betaWL_dVds ;
    Ids0_dVgs = T3 * Idd_dVgs + T1 * Mu_dVgs + T2 * betaWL_dVgs ;

    /*-----------------------------------------------------------*
     * Adding parasitic components to the channel current.
     *-----------------*/
    if( model->HSM2_ptl != 0 ){
      T1 =  0.5 * ( Vds - Pds ) ;
      Fn_SymAdd( T6 , T1 , 0.01 , T2 ) ;
      T2 *= 0.5 ;
      T6_dVb = T2 * ( - Pds_dVbs ) ;
      T6_dVd = T2 * ( 1.0 - Pds_dVds ) ;
      T6_dVg = T2 * ( - Pds_dVgs ) ;

      T1     = 1.1 - ( Ps0 + T6 );
      T1_dVb =     - ( Ps0_dVbs + T6_dVb );
      T1_dVd =     - ( Ps0_dVds + T6_dVd );
      T1_dVg =     - ( Ps0_dVgs + T6_dVg );

      Fn_SZ( T2 , T1 , 0.05 , T0 ) ;
      T2 += small ;
      T2_dVb = T1_dVb * T0 ;
      T2_dVd = T1_dVd * T0 ;
      T2_dVg = T1_dVg * T0 ;

      T0 = beta * here->HSM2_ptl0 ;
      T3 = Cox * T0 ;
      T3_dVb = Cox_dVb * T0 ;
      T3_dVd = Cox_dVd * T0 ;
      T3_dVg = Cox_dVg * T0 ;
      T0 = pow( T2 , model->HSM2_ptp ) ;
      T9     = T3 * T0 ;
      T9_dVb = T3 * model->HSM2_ptp * T0 / T2 * T2_dVb + T3_dVb * T0 ;
      T9_dVd = T3 * model->HSM2_ptp * T0 / T2 * T2_dVd + T3_dVd * T0 ;
      T9_dVg = T3 * model->HSM2_ptp * T0 / T2 * T2_dVg + T3_dVg * T0 ;


      T4 = 1.0 + Vdsz * model->HSM2_pt2 ;
      T4_dVb = 0.0 ;
      T4_dVd = Vdsz_dVds * model->HSM2_pt2 ;
      T4_dVg = 0.0 ;

      T0 = here->HSM2_pt40 ;
      T5 = Ps0 + T6 - Vbsz ;
      T5_dVb = Ps0_dVbs + T6_dVb - Vbsz_dVbs ;
      T5_dVd = Ps0_dVds + T6_dVd - Vbsz_dVds ;
      T5_dVg = Ps0_dVgs + T6_dVg ;
      T4 += Vdsz * T0 * T5 ;
      T4_dVb += Vdsz * T0 * T5_dVb ;
      T4_dVd += Vdsz * T0 * T5_dVd + Vdsz_dVds * T0 * T5 ;
      T4_dVg += Vdsz * T0 * T5_dVg ;

      T6     = T9     * T4 ;
      T9_dVb = T9_dVb * T4 + T9 * T4_dVb ;
      T9_dVd = T9_dVd * T4 + T9 * T4_dVd ;
      T9_dVg = T9_dVg * T4 + T9 * T4_dVg ;
      T9     = T6 ;

    }else{
      T9 = 0.0 ;
      T9_dVb = 0.0 ;
      T9_dVd = 0.0 ;
      T9_dVg = 0.0 ;
    }


    if( model->HSM2_gdl != 0 ){
      T1 = beta * here->HSM2_gdl0 ;
      T2 = Cox * T1 ;
      T2_dVb = Cox_dVb * T1 ;
      T2_dVd = Cox_dVd * T1 ;
      T2_dVg = Cox_dVg * T1 ;
      T8     = T2 * Vdsz ;
      T8_dVb = T2_dVb * Vdsz ;
      T8_dVd = T2_dVd * Vdsz + T2 * Vdsz_dVds ;
      T8_dVg = T2_dVg * Vdsz ;
    }else{
      T8 = 0.0 ;
      T8_dVb = 0.0 ;
      T8_dVd = 0.0 ;
      T8_dVg = 0.0 ;
    }


    if ( ( T9 + T8 ) > 0.0 ) {
      Idd1 = Pds * ( T9 + T8 ) ;
      Idd1_dVbs = Pds_dVbs * ( T9 + T8 ) + Pds * ( T9_dVb + T8_dVb ) ;
      Idd1_dVds = Pds_dVds * ( T9 + T8 ) + Pds * ( T9_dVd + T8_dVd ) ;
      Idd1_dVgs = Pds_dVgs * ( T9 + T8 ) + Pds * ( T9_dVg + T8_dVg ) ;

      Ids0 += betaWL * Idd1 * Mu ;
      T1 = betaWL * Idd1 ;
      T2 = Idd1 * Mu ;
      T3 = Mu * betaWL ;
      Ids0_dVbs += T3 * Idd1_dVbs + T1 * Mu_dVbs + T2 * betaWL_dVbs ;
      Ids0_dVds += T3 * Idd1_dVds + T1 * Mu_dVds + T2 * betaWL_dVds ;
      Ids0_dVgs += T3 * Idd1_dVgs + T1 * Mu_dVgs + T2 * betaWL_dVgs ;
    }


    /* (note: rpock procedure was removed. (2006.04.20) */
    if ( flg_rsrd == 2 ) {
      Rd  = here->HSM2_rd ;
      T0 = Rd * Ids0 ;
      T1 = Vds + small ;
      T2 = 1.0 / T1 ;
      T3 = 1.0 + T0 * T2 ;
      T3_dVb =   Rd * Ids0_dVbs             * T2 ;
      T3_dVd = ( Rd * Ids0_dVds * T1 - T0 ) * T2 * T2 ;
      T3_dVg =   Rd * Ids0_dVgs             * T2 ;
      T4 = 1.0 / T3 ;
      Ids = Ids0 * T4 ;
      T5 = T4 * T4 ;
      Ids_dVbs = ( Ids0_dVbs * T3 - Ids0 * T3_dVb ) * T5 ;
      Ids_dVds = ( Ids0_dVds * T3 - Ids0 * T3_dVd ) * T5 ;
      Ids_dVgs = ( Ids0_dVgs * T3 - Ids0 * T3_dVg ) * T5 ;
    } else {
      Ids = Ids0 ;
      Ids_dVbs = Ids0_dVbs ;
      Ids_dVds = Ids0_dVds ;
      Ids_dVgs = Ids0_dVgs ;
    }
    
    if ( Pds < ps_conv ) {
      Ids_dVbs = 0.0 ;
      Ids_dVgs = 0.0 ;
    }
  
    Ids += Gdsmin * Vds ;
    Ids_dVds += Gdsmin ;


    /*-----------------------------------------------------------*
     * STI
     *-----------------*/
    if ( model->HSM2_coisti != 0 ) {
      /*---------------------------------------------------*
       * dVthSCSTI : Short-channel effect induced by Vds (STI).
       *-----------------*/      
      T1 = C_ESI * Cox_inv ;
      T2 = here->HSM2_wdpl ;
      T3 =  here->HSM2_lgatesm - model->HSM2_parl2 ;
      T4 = 1.0 / (T3 * T3) ;
      T5 = 2.0 * (model->HSM2_vbi - Pb20b) * T1 * T2 * T4 ;
      
      dVth0 = T5 * sqrt_Pbsum ;
      T6 = T5 * 0.5 / sqrt_Pbsum ;
      T7 = 2.0 * (model->HSM2_vbi - Pb20b) * C_ESI * T2 * T4 * sqrt_Pbsum ;
      T8 = - 2.0 * T1 * T2 * T4 * sqrt_Pbsum ;
      dVth0_dVb = T6 * Pbsum_dVb + T7 * Cox_inv_dVb + T8 * Pb20b_dVb ;
      dVth0_dVd = T6 * Pbsum_dVd + T7 * Cox_inv_dVd + T8 * Pb20b_dVd ;
      dVth0_dVg = T6 * Pbsum_dVg + T7 * Cox_inv_dVg + T8 * Pb20b_dVg ;

      T4 = pParam->HSM2_scsti1 ;
      T6 = pParam->HSM2_scsti2 ;
      T1  = T4 + T6 * Vdsz ;
      dVthSCSTI = dVth0 * T1 ;
      dVthSCSTI_dVb = dVth0_dVb * T1 ;
      dVthSCSTI_dVd = dVth0_dVd * T1 + dVth0 * T6 * Vdsz_dVds ;
      dVthSCSTI_dVg = dVth0_dVg * T1 ;

      T1 = pParam->HSM2_vthsti - model->HSM2_vdsti * Vds ;
      T1_dVd = - model->HSM2_vdsti ;

      Vgssti = Vgsz - Vfb + T1 + dVthSCSTI ;
      Vgssti_dVbs = dVthSCSTI_dVb ;
      Vgssti_dVds = Vgsz_dVds + T1_dVd + dVthSCSTI_dVd ;
      Vgssti_dVgs = Vgsz_dVgs + dVthSCSTI_dVg ;
      
      costi0 = here->HSM2_costi0 ;
      costi1 = here->HSM2_costi1 ;

      costi3 = here->HSM2_costi0_p2 * Cox_inv * Cox_inv ;
      T1 = 2.0 * here->HSM2_costi0_p2 * Cox_inv ;
      costi3_dVb = T1 * Cox_inv_dVb ;
      costi3_dVd = T1 * Cox_inv_dVd ;
      costi3_dVg = T1 * Cox_inv_dVg ;
      T2 = 1.0 / costi3 ;
      costi3_dVb_c3 = costi3_dVb * T2 ;
      costi3_dVd_c3 = costi3_dVd * T2 ;
      costi3_dVg_c3 = costi3_dVg * T2 ;

      costi4 = costi3 * beta * 0.5 ;
      costi5 = costi4 * beta * 2.0 ;

      T11 = beta * 0.25 ;
      T10 = beta_inv - costi3 * T11 + Vfb - pParam->HSM2_vthsti - dVthSCSTI + small ;
      T10_dVb = - T11 * costi3_dVb - dVthSCSTI_dVb ;
      T10_dVd = - T11 * costi3_dVd - dVthSCSTI_dVd ;
      T10_dVg = - T11 * costi3_dVg - dVthSCSTI_dVg ;
      T1 = Vgsz - T10 - psisti_dlt ;
      T1_dVb = - T10_dVb ;
      T1_dVd = Vgsz_dVds - T10_dVd ;
      T1_dVg = Vgsz_dVgs - T10_dVg ;
      T0 = Fn_Sgn(T10) ;
      T2 = sqrt (T1 * T1 + T0 * 4.0 * T10 * psisti_dlt) ;
      T3 = T10 + 0.5 * (T1 + T2) - Vfb + pParam->HSM2_vthsti + dVthSCSTI - Vbsz ;
      T3_dVb = T10_dVb + 0.5 * (T1_dVb + (T1 * T1_dVb + T0 * 2.0 * T10_dVb * psisti_dlt) / T2) 
        + dVthSCSTI_dVb - Vbsz_dVbs ;
      T3_dVd = T10_dVd + 0.5 * (T1_dVd + (T1 * T1_dVd + T0 * 2.0 * T10_dVd * psisti_dlt) / T2) 
        + dVthSCSTI_dVd - Vbsz_dVds ;
      T3_dVg = T10_dVg + 0.5 * (T1_dVg + (T1 * T1_dVg + T0 * 2.0 * T10_dVg * psisti_dlt) / T2) 
        + dVthSCSTI_dVg ;
      T4 = beta * T3 - 1.0 ;
      T5 = 4.0 / costi5 ;
      T1 = 1.0 + T4 * T5 ;
      T6 = beta * T5 ;
      T7 = T4 * T5 ;
      T1_dVb = (T6 * T3_dVb - T7 * costi3_dVb_c3) ;
      T1_dVd = (T6 * T3_dVd - T7 * costi3_dVd_c3) ;
      T1_dVg = (T6 * T3_dVg - T7 * costi3_dVg_c3) ;
      Fn_SZ( T1 , T1, 1.0e-2, T2) ;
      T1 += small ;
      T1_dVb *= T2 ;
      T1_dVd *= T2 ;
      T1_dVg *= T2 ;
      costi6 = sqrt (T1) ;
      T0 = costi4 * (1.0 - costi6) ;
      Psasti = Vgssti + T0 ;
      T2 = 0.5 * costi4 / costi6 ;
      Psasti_dVbs = Vgssti_dVbs + costi3_dVb_c3 * T0 - T2 * T1_dVb ;
      Psasti_dVds = Vgssti_dVds + costi3_dVd_c3 * T0 - T2 * T1_dVd ;
      Psasti_dVgs = Vgssti_dVgs + costi3_dVg_c3 * T0 - T2 * T1_dVg ;

      T0 = 1.0 / (beta + 2.0 / (Vgssti + small)) ;
      Psbsti = log (1.0 / costi1 / costi3 * (Vgssti * Vgssti)) * T0 ;
      T2 = 2.0 * T0 / (Vgssti + small) ;
      T3 = Psbsti / (Vgssti + small) ;
      Psbsti_dVbs = T2 * (Vgssti_dVbs - 0.5 * costi3_dVb_c3 * Vgssti 
                                + T3 * Vgssti_dVbs ) ;
      Psbsti_dVds = T2 * (Vgssti_dVds - 0.5 * costi3_dVd_c3 * Vgssti 
                                + T3 * Vgssti_dVds ) ;
      Psbsti_dVgs = T2 * (Vgssti_dVgs - 0.5 * costi3_dVg_c3 * Vgssti 
                                + T3 * Vgssti_dVgs ) ;
      
      Psab = Psbsti - Psasti - sti2_dlt ;
      Psab_dVbs = Psbsti_dVbs - Psasti_dVbs ;
      Psab_dVds = Psbsti_dVds - Psasti_dVds ;
      Psab_dVgs = Psbsti_dVgs - Psasti_dVgs ;
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

      T0 = costi1 * exp (beta * Psti) ;
      T1 = beta * (Psti - Vbsz) - 1.0 + T0 ;
      T1_dVb = beta * ((Psti_dVbs - Vbsz_dVbs) + T0 * Psti_dVbs) ;
      T1_dVd = beta * ((Psti_dVds - Vbsz_dVds) + T0 * Psti_dVds) ;
      T1_dVg = beta * (Psti_dVgs + T0 * Psti_dVgs) ;
      Fn_SZ ( T1 , T1, 1.0e-2, T0) ;
          T1 += epsm10 ;
      T1_dVb *= T0 ;
      T1_dVd *= T0 ;
      T1_dVg *= T0 ;
      sq1sti = sqrt (T1);
      T2 = 0.5 / sq1sti ;
      sq1sti_dVbs = T2 * T1_dVb ;
      sq1sti_dVds = T2 * T1_dVd ;
      sq1sti_dVgs = T2 * T1_dVg ;

      T1 = beta * (Psti - Vbsz) - 1.0;
      T1_dVb = beta * (Psti_dVbs - Vbsz_dVbs) ;
      T1_dVd = beta * (Psti_dVds - Vbsz_dVds) ;
      T1_dVg = beta * Psti_dVgs ;
      Fn_SZ( T1 , T1, 1.0e-2, T0) ;
          T1 += epsm10 ;
      T1_dVb *= T0 ;
      T1_dVd *= T0 ;
      T1_dVg *= T0 ;
      sq2sti = sqrt (T1);
      T2 = 0.5 / sq2sti ;
      sq2sti_dVbs = T2 * T1_dVb ;
      sq2sti_dVds = T2 * T1_dVd ;
      sq2sti_dVgs = T2 * T1_dVg ;

      Qn0sti = costi0 * (sq1sti - sq2sti) ;     
      Qn0sti_dVbs = costi0 * (sq1sti_dVbs - sq2sti_dVbs) ;
      Qn0sti_dVds = costi0 * (sq1sti_dVds - sq2sti_dVds) ;
      Qn0sti_dVgs = costi0 * (sq1sti_dVgs - sq2sti_dVgs) ;

      /* T1: Vdsatsti */
      T1 = Psasti - Psti ;
      T1_dVb = Psasti_dVbs - Psti_dVbs ;
      T1_dVd = Psasti_dVds - Psti_dVds ;
      T1_dVg = Psasti_dVgs - Psti_dVgs ;

      Fn_SZ( T1 , T1 , 1.0e-1 , T2 ) ;
      T1 += epsm10 ;
      T1_dVb *= T2 ;
      T1_dVd *= T2 ;
      T1_dVg *= T2 ;


      TX = Vds / T1 ;
      T2 = 1.0 / ( T1 * T1 ) ;
      TX_dVbs = T2 * ( - Vds * T1_dVb ) ;
      TX_dVds = T2 * ( T1 - Vds * T1_dVd ) ;
      TX_dVgs = T2 * ( - Vds * T1_dVg ) ;

      Fn_CP( TY , TX , 1.0 , 4 , T2 ) ;
      TY_dVbs = T2 * TX_dVbs ;
      TY_dVds = T2 * TX_dVds ;
      TY_dVgs = T2 * TX_dVgs ;

      costi7 = 2.0 * here->HSM2_wsti * here->HSM2_nf * beta_inv ;
      T1 = Lch ;
      Idssti = costi7 * Mu * Qn0sti * TY / T1;
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

      Ids = Ids + Idssti ;
      Ids_dVbs = Ids_dVbs + Idssti_dVbs ;
      Ids_dVds = Ids_dVds + Idssti_dVds ;
      Ids_dVgs = Ids_dVgs + Idssti_dVgs ;
    
    }

    /*-----------------------------------------------------------*
     * Break point for the case of Rs=Rd=0.
     *-----------------*/
    if ( flg_rsrd == 0 ) {
      DJI = 1.0 ;
      break ;
    }

    /*-----------------------------------------------------------*
     * calculate corrections of biases.
     * - Fbs = 0, etc. are the small ciucuit equations.
     * - DJ, Jacobian of the small circuit matrix, is g.t. 1 
     *   provided Rs, Rd and conductances are positive.
     *-----------------*/
    Fbs = Vbs - Vbsc + Ids * Rs ;
    Fds = Vds - Vdsc + Ids * ( Rs + Rd ) ;
    Fgs = Vgs - Vgsc + Ids * Rs ;

    DJ = 1.0 + Rs * Ids_dVbs + ( Rs + Rd ) * Ids_dVds + Rs * Ids_dVgs ;
    DJI = 1.0 / DJ ;

    JI11 = 1 + ( Rs + Rd ) * Ids_dVds + Rs * Ids_dVgs ;
    JI12 = - Rs * Ids_dVds ;
    JI13 = - Rs * Ids_dVgs ;
    JI21 = - ( Rs + Rd ) * Ids_dVbs ;
    JI22 = 1 + Rs * Ids_dVbs + Rs * Ids_dVgs ;
    JI23 = - ( Rs + Rd ) * Ids_dVgs ;
    JI31 = - Rs * Ids_dVbs ;
    JI32 = - Rs * Ids_dVds ;
    JI33 = 1 + Rs * Ids_dVbs + ( Rs + Rd ) * Ids_dVds ;

    dVbs = - DJI * ( JI11 * Fbs + JI12 * Fds + JI13 * Fgs ) ;
    dVds = - DJI * ( JI21 * Fbs + JI22 * Fds + JI23 * Fgs ) ;
    dVgs = - DJI * ( JI31 * Fbs + JI32 * Fds + JI33 * Fgs ) ;

    dV_sum = fabs( dVbs ) + fabs( dVds ) + fabs( dVgs ) ; 


    /*-----------------------------------------------------------*
     * Break point for converged case.
     * - Exit from the bias loop.
     * - NOTE: Update of internal biases is avoided.
     *-----------------*/
    if ( Ids_last * Ids_tol >= fabs( Ids_last - Ids ) || mini_current >= fabs( Ids_last - Ids )
            || dV_sum < ps_conv ) break ;

    /*-----------------------------------------------------------*
     * Update the internal biases.
     *-----------------*/
    Vbs  = Vbs + dVbs ;
    Vds  = Vds + dVds ;
    Vgs  = Vgs + dVgs ;
    if ( Vds < 0.0 ) { 
      Vds  = 0.0 ; 
      dVds = 0.0 ; 
    } 

    /*-----------------------------------------------------------*
     * Bottom of bias loop. (label) 
     *-----------------*/
/*  bottom_of_bias_loop :*/

    /*-----------------------------------------------------------*
     * Make initial guess flag of potential ON.
     * - This effects for the 2nd and later iterations of bias loop.
     *-----------------*/
    flg_pprv = 1 ;

  } /*++ End of the bias loop +++++++++++++++++++++++++++++*/

  if ( lp_bs > lp_bs_max ) { lp_bs -- ; }



  /*----------------------------------------------------------*
   * induced gate noise. ( Part 1/3 )
   *----------------------*/
  if ( model->HSM2_coign != 0 && model->HSM2_cothrml != 0 ) {
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
   * Evaluate integrated chages in unit [C].
   *----------------------*/

  T1 = - here->HSM2_weff_nf * Leff ;

  Qb = T1 * Qbu ;
  Qb_dVbs = T1 * Qbu_dVbs ;
  Qb_dVds = T1 * Qbu_dVds ;
  Qb_dVgs = T1 * Qbu_dVgs ;
    
  Qi = T1 * Qiu ;
  Qi_dVbs = T1 * Qiu_dVbs ;
  Qi_dVds = T1 * Qiu_dVds ;
  Qi_dVgs = T1 * Qiu_dVgs ;
 
  Qd = Qi * Qdrat ;
  Qd_dVbs = Qi_dVbs * Qdrat + Qi * Qdrat_dVbs ;
  Qd_dVds = Qi_dVds * Qdrat + Qi * Qdrat_dVds ;
  Qd_dVgs = Qi_dVgs * Qdrat + Qi * Qdrat_dVgs ;


  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-2: Substrate / gate / leak currents
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 
  /*-----------------------------------------------------------*
   * Isub : substrate current induced by impact ionization.
   *-----------------*/

  if ( flg_noqi == 1 || model->HSM2_coisub == 0 ) {
   /* Accumulation zone or nonconductive case, in which Ids==0. */
    Isub = 0.0e0 ;
    Isub_dVbs = Isub_dVds = Isub_dVgs = 0.0e0 ;
  } else {
   /*-------------------------------------------*
    * Conductive case. 
    *-----------------*/
    if ( pParam->HSM2_sub1 > 0.0e0 && pParam->HSM2_vmax > 0.0e0 ) {
      T0 = here->HSM2_vg2const ;
      T1 = T0 * Vgp ; 
      T1_dVd = T0 * Vgp_dVds ;
      T1_dVg = T0 * Vgp_dVgs ;
      T1_dVb = T0 * Vgp_dVbs ;

      T7 = Cox0 * Cox0 ;
      T8 = here->HSM2_qnsub_esi ;
      T3 = T8 / T7 ;

      T9 = 2.0 / T8 ;
      T4 = 1.0e0 + T9 * T7 ;

      T2 = here->HSM2_xvbs ;
      T5 = T1 - beta_inv - T2 * Vbsz ;
      T5_dVd = T1_dVd - T2 * Vbsz_dVds;
      T5_dVg = T1_dVg ;
      T5_dVb = T1_dVb - T2 * Vbsz_dVbs;

      T6 = T4 * T5 ;
      T6_dVd = T4 * T5_dVd ;
      T6_dVg = T4 * T5_dVg ;
      T6_dVb = T4 * T5_dVb ;
      Fn_SZ( T6 , T6, 1.0e-3, T9) ;
      T6 += small ;
      T6_dVd *= T9 ;
      T6_dVg *= T9 ;
      T6_dVb *= T9 ;
      T6 = sqrt( T6 ) ;
      T9 = 0.5 / T6 ;
      T6_dVd = T9 * T6_dVd ;
      T6_dVg = T9 * T6_dVg ;
      T6_dVb = T9 * T6_dVb ;
  
      Psislsat = T1 + T3 * ( 1.0 - T6 ) ;
      Psislsat_dVd = T1_dVd - T3 * T6_dVd ;
      Psislsat_dVg = T1_dVg - T3 * T6_dVg ;
      Psislsat_dVb = T1_dVb - T3 * T6_dVb ;

      T2 = here->HSM2_lgate / (here->HSM2_xgate + here->HSM2_lgate) ;

      Psisubsat = pParam->HSM2_svds * Vdsz + Ps0z - T2 * Psislsat ;      
      Psisubsat_dVd = pParam->HSM2_svds * Vdsz_dVds + Ps0z_dVds - T2 * Psislsat_dVd ; 
      Psisubsat_dVg = Ps0z_dVgs - T2 * Psislsat_dVg ; 
      Psisubsat_dVb = Ps0z_dVbs - T2 * Psislsat_dVb ; 
      Fn_SZ( Psisubsat , Psisubsat, 1.0e-3, T9 ) ; 
      Psisubsat += small ;
      Psisubsat_dVd *= T9 ;
      Psisubsat_dVg *= T9 ;
      Psisubsat_dVb *= T9 ;       
 
      T5 = here->HSM2_xsub1 ;
      T6 = here->HSM2_xsub2 ;
      T2 = exp( - T6 / Psisubsat ) ;
      T3 = T2 * T6 / ( Psisubsat * Psisubsat ) ;
      T2_dVd = T3 * Psisubsat_dVd ;
      T2_dVg = T3 * Psisubsat_dVg ;
      T2_dVb = T3 * Psisubsat_dVb ;

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
    } else {
      Isub = 0.0e0 ;
      Isub_dVbs = Isub_dVds = Isub_dVgs = 0.0e0 ;
    } /* end of if ( pParam->HSM2_sub1 ... ) else block. */
  }


    /*---------------------------------------------------*
     * Impact-Ionization Induced Bulk Potential Change (IBPC)
     *-----------------*/
  if ( flg_noqi == 0 && Isub > 0e0 && pParam->HSM2_ibpc1 != 0e0 ) {

    /* delta Vbs */
    T0 = 1e0 + pParam->HSM2_ibpc2 * dVth ;
    dVbsIBPC = pParam->HSM2_ibpc1 * T0 * Isub ;
    dVbsIBPC_dVbs = pParam->HSM2_ibpc1 * ( pParam->HSM2_ibpc2 * dVth_dVb * Isub + T0 * Isub_dVbs ) ;
    dVbsIBPC_dVds = pParam->HSM2_ibpc1 * ( pParam->HSM2_ibpc2 * dVth_dVd * Isub + T0 * Isub_dVds ) ;
    dVbsIBPC_dVgs = pParam->HSM2_ibpc1 * ( pParam->HSM2_ibpc2 * dVth_dVg * Isub + T0 * Isub_dVgs ) ;

    /* dG3 & dG4 */
    T10 = 1e0 / Xi0 ;
    T1 = beta * dVbsIBPC * T10 ;
    T10 *= T10 ;
    T1_dVb = beta * ( dVbsIBPC_dVbs * Xi0 - dVbsIBPC * Xi0_dVbs ) * T10 ;
    T1_dVd = beta * ( dVbsIBPC_dVds * Xi0 - dVbsIBPC * Xi0_dVds ) * T10 ;
    T1_dVg = beta * ( dVbsIBPC_dVgs * Xi0 - dVbsIBPC * Xi0_dVgs ) * T10 ;

    T10 = 1e0 / Xil ;
    T2 = beta * dVbsIBPC * T10 ;
    T10 *= T10 ;
    T2_dVb = beta * ( dVbsIBPC_dVbs * Xil - dVbsIBPC * Xil_dVbs ) * T10 ;
    T2_dVd = beta * ( dVbsIBPC_dVds * Xil - dVbsIBPC * Xil_dVds ) * T10 ;
    T2_dVg = beta * ( dVbsIBPC_dVgs * Xil - dVbsIBPC * Xil_dVgs ) * T10 ;

    dG3 = cnst0 * ( Xilp32 * T2 - Xi0p32 * T1 ) ;
    dG3_dVbs = cnst0 * ( Xilp32_dVbs * T2 + Xilp32 * T2_dVb - Xi0p32_dVbs * T1 - Xi0p32 * T1_dVb ) ;
    dG3_dVds = cnst0 * ( Xilp32_dVds * T2 + Xilp32 * T2_dVd - Xi0p32_dVds * T1 - Xi0p32 * T1_dVd ) ;
    dG3_dVgs = cnst0 * ( Xilp32_dVgs * T2 + Xilp32 * T2_dVg - Xi0p32_dVgs * T1 - Xi0p32 * T1_dVg ) ;
    dG4 = cnst0 * 0.5 * ( - Xilp12 * T2 + Xi0p12 * T1 ) ;
    dG4_dVbs = cnst0 * 0.5 * ( - Xilp12_dVbs * T2 - Xilp12 * T2_dVb + Xi0p12_dVbs * T1 + Xi0p12 * T1_dVb ) ;
    dG4_dVds = cnst0 * 0.5 * ( - Xilp12_dVds * T2 - Xilp12 * T2_dVd + Xi0p12_dVds * T1 + Xi0p12 * T1_dVd ) ;
    dG4_dVgs = cnst0 * 0.5 * ( - Xilp12_dVgs * T2 - Xilp12 * T2_dVg + Xi0p12_dVgs * T1 + Xi0p12 * T1_dVg ) ;

    /* Add IBPC current into Ids */
    dIdd = dG3 + dG4 ;
    dIdd_dVbs = dG3_dVbs + dG4_dVbs ;
    dIdd_dVds = dG3_dVds + dG4_dVds ;
    dIdd_dVgs = dG3_dVgs + dG4_dVgs ;
    IdsIBPC = betaWL * dIdd * Mu ;
    IdsIBPC_dVbs = betaWL * ( Mu * dIdd_dVbs + dIdd * Mu_dVbs ) + betaWL_dVbs * Mu * dIdd ;
    IdsIBPC_dVds = betaWL * ( Mu * dIdd_dVds + dIdd * Mu_dVds ) + betaWL_dVds * Mu * dIdd ;
    IdsIBPC_dVgs = betaWL * ( Mu * dIdd_dVgs + dIdd * Mu_dVgs ) + betaWL_dVgs * Mu * dIdd ;

  } /* End if (IBPC) */

  /*-----------------------------------------------------------*
   * Igate : Gate current induced by tunneling.
   *-----------------*/
  if ( model->HSM2_coiigs == 0 ) {
    Igate = 0.0 ;
    Igate_dVbs = Igate_dVds = Igate_dVgs = 0.0 ;
    Igs = 0.0 ;
    Igs_dVbs = Igs_dVds = Igs_dVgs = 0.0 ;
    Igd = 0.0 ;
    Igd_dVbs = Igd_dVds = Igd_dVgs = 0.0 ;
    Igb = 0.0 ;
    Igb_dVbs = Igb_dVds = Igb_dVgs = 0.0 ;
    GLPART1 = 0.0 ;
    GLPART1_dVgs = GLPART1_dVds = GLPART1_dVbs = 0.0 ;
  } else {


    /* Igate */
    if ( flg_noqi == 0 ) {
      Psdlz = Ps0z + Vdsz - epsm10 ;
      Psdlz_dVbs = Ps0z_dVbs ;
      Psdlz_dVds = Ps0z_dVds + Vdsz_dVds ;
      Psdlz_dVgs = Ps0z_dVgs ;

      T1 = Vgsz - Vfb + modelMKS->HSM2_gleak4 * (dVth - dPpg) * Leff - Psdlz * pParam->HSM2_gleak3 ;
      T1_dVg = Vgsz_dVgs + modelMKS->HSM2_gleak4 * (dVth_dVg - dPpg_dVg) * Leff - Psdlz_dVgs * pParam->HSM2_gleak3 ;
      T1_dVd = Vgsz_dVds + modelMKS->HSM2_gleak4 * (dVth_dVd - dPpg_dVd) * Leff - Psdlz_dVds * pParam->HSM2_gleak3 ;
      T1_dVb = modelMKS->HSM2_gleak4 * ( dVth_dVb - dPpg_dVb ) * Leff - Psdlz_dVbs * pParam->HSM2_gleak3 ;

      T3 = 1.0 / Tox0 ;
      T2 = T1 * T3 ;
      T2_dVg = (T1_dVg ) * T3 ;
      T2_dVd = (T1_dVd ) * T3 ;
      T2_dVb = (T1_dVb ) * T3 ;
      
      T3 = 1.0 / modelMKS->HSM2_gleak5 ;

      if ( VgVt <= VgVt_small ) {
        Ey = 0.0 ;
        Ey_dVgs = 0.0 ;
        Ey_dVds = 0.0 ;
        Ey_dVbs = 0.0 ;
      }
      T7 = 1.0 + Ey * T3 ; 
      T7_dVg = Ey_dVgs * T3 ;
      T7_dVd = Ey_dVds * T3 ;
      T7_dVb = Ey_dVbs * T3 ;

      Etun = T2 * T7 ;
      Etun_dVgs = T2_dVg * T7 + T7_dVg * T2 ;
      Etun_dVds = T2_dVd * T7 + T7_dVd * T2 ;
      Etun_dVbs = T2_dVb * T7 + T7_dVb * T2 ;

      Fn_SZ( Etun , Etun , igate_dlt , T5 ) ;
      Etun_dVgs *= T5 ;
      Etun_dVds *= T5 ;
      Etun_dVbs *= T5 ;

      Fn_SZ( T3 , Vgsz , 1.0e-3 , T4 ) ;
      T3 -= model->HSM2_vzadd0 ;
      TX = T3 / cclmmdf ; 
      T2 = 1.0 +  TX * TX ;
      T1 = 1.0 - 1.0 / T2 ;
      T1_dVg = 2.0 * TX * T4 / ( T2 * T2 * cclmmdf ) ;
      T1_dVd = T1_dVg * Vgsz_dVds ;
      Etun_dVgs = T1 * Etun_dVgs + Etun * T1_dVg ;
      Etun_dVds = T1 * Etun_dVds + Etun * T1_dVd ;
      Etun_dVbs *= T1 ;
      Etun *= T1 ;

      T0 = Leff * here->HSM2_weff_nf ;
      T7 = modelMKS->HSM2_gleak7 / (modelMKS->HSM2_gleak7 + T0) ;
      
      T6 = pParam->HSM2_gleak6 ;
      T9 = T6 / (T6 + Vdsz) ;
      T9_dVd = - T9 / (T6 + Vdsz) * Vdsz_dVds ;
      
      T1 = - pParam->HSM2_gleak2 * Egp32 / (Etun + small) ;

      if ( T1 < - EXP_THR ) {
        Igate = 0.0 ;
        Igate_dVbs = Igate_dVds = Igate_dVgs = 0.0 ;
      } else {
        T2 = exp ( T1 ) ;
        T3 = pParam->HSM2_gleak1 / Egp12 * C_QE * T0 ;
  
        T4 =  T2 * T3 * sqrt ((Qiu + Cox0 * VgVt_small ) / cnst0) ;
        T5 = T4 * Etun ;
        T6 = 0.5 * Etun / (Qiu + Cox0 * VgVt_small ) ;
        T10 = T5 * Etun ;
        T10_dVb = T5 * (2.0 * Etun_dVbs - T1 * Etun_dVbs + T6 * Qiu_dVbs) ;
        T10_dVd = T5 * (2.0 * Etun_dVds - T1 * Etun_dVds + T6 * Qiu_dVds) ;
        T10_dVg = T5 * (2.0 * Etun_dVgs - T1 * Etun_dVgs + T6 * Qiu_dVgs) ;
          
        Igate = T7 * T9 * T10 ;
        Igate_dVbs = T7 * T9 * T10_dVb ;
        Igate_dVds = T7 * (T9_dVd * T10 + T9 * T10_dVd) ;
        Igate_dVgs = T7 * T9 * T10_dVg ;
      }
    } else {
      Igate = 0.0 ;
      Igate_dVbs = Igate_dVds = Igate_dVgs = 0.0 ;
    }

    /* Igs */
      T0 = - pParam->HSM2_glksd2 * Vgs + modelMKS->HSM2_glksd3 ;
      T2 = exp (Tox0 * T0);
      T2_dVg = (- Tox0 * pParam->HSM2_glksd2) * T2;
        
      T0 = Vgs / Tox0 / Tox0 ;
      T3 = Vgs * T0 ;
      T3_dVg = 2.0 * T0 * (1.0 ) ;
      T4 = pParam->HSM2_glksd1 / 1.0e6 * here->HSM2_weff_nf ;
      Igs = T4 * T2 * T3 ;
      Igs_dVgs = T4 * (T2_dVg * T3 + T2 * T3_dVg) ;
      Igs_dVds = 0.0 ;
      Igs_dVbs = 0.0 ;
        
      if ( Vgs >= 0.0e0 ){
        Igs *= -1.0 ;
        Igs_dVgs *= -1.0 ;
        Igs_dVds *= -1.0 ; 
        Igs_dVbs *= -1.0 ;
      }


    /* Igd */
      T1 = Vgs - Vds ;
      T0 = - pParam->HSM2_glksd2 * T1 + modelMKS->HSM2_glksd3 ;
      T2 = exp (Tox0 * T0);
      T2_dVg = (- Tox0 * pParam->HSM2_glksd2) * T2;
      T2_dVd = (+ Tox0 * pParam->HSM2_glksd2) * T2;
      T2_dVb = 0.0 ;
        
      T0 = T1 / Tox0 / Tox0 ;
      T3 = T1 * T0 ;
      T3_dVg = 2.0 * T0 ;
      T3_dVd = - 2.0 * T0 ;
      T3_dVb = 0.0 ;
      T4 = pParam->HSM2_glksd1 / 1.0e6 * here->HSM2_weff_nf ;
      Igd = T4 * T2 * T3 ;
      Igd_dVgs = T4 * (T2_dVg * T3 + T2 * T3_dVg) ;
      Igd_dVds = T4 * (T2_dVd * T3 + T2 * T3_dVd) ;
      Igd_dVbs = 0.0 ;

      if( T1 >= 0.0e0 ){
        Igd  *= -1.0 ;
        Igd_dVgs *= -1.0 ;
        Igd_dVds *= -1.0 ;
        Igd_dVbs *= -1.0 ;
      }


    /* Igb */
      Etun = (- Vgs + Vbs + Vfb + model->HSM2_glkb3 ) / Tox0 ;

      Etun_dVgs = - 1.0 / Tox0 ;
      Etun_dVds = 0.0 ;
      Etun_dVbs = 1.0 / Tox0 ;

      Fn_SZ( Etun , Etun, igate_dlt, T5) ;
      Etun += small ;
      Etun_dVgs *= T5 ;
      Etun_dVbs *= T5 ;

      T1 = - pParam->HSM2_glkb2 / Etun ;
      if ( T1 < - EXP_THR ) {
        Igb = 0.0 ;
        Igb_dVgs = Igb_dVds = Igb_dVbs = 0.0 ;
      } else {
        T2 = exp ( T1 );
        T3 =  pParam->HSM2_glkb2 / ( Etun * Etun ) * T2 ;
        T2_dVg = T3 * Etun_dVgs ;
        T2_dVb = T3 * Etun_dVbs ;
          
        T3 = pParam->HSM2_glkb1 * here->HSM2_weff_nf * Leff ;
        Igb = T3 * Etun * Etun * T2 ;
        Igb_dVgs = T3 * (2.0 * Etun * Etun_dVgs * T2 + Etun * Etun * T2_dVg);
        Igb_dVds = 0.0 ;
        Igb_dVbs = T3 * (2.0 * Etun * Etun_dVbs * T2 + Etun * Etun * T2_dVb);
      }

        GLPART1 = 0.5 ;
	GLPART1_dVgs = 0.0 ;
	GLPART1_dVds = 0.0 ;
	GLPART1_dVbs = 0.0 ;

  } /* if ( model->HSM2_coiigs == 0 ) */

  
  

  /*-----------------------------------------------------------*
   * Igidl : GIDL 
   *-----------------*/
  if( model->HSM2_cogidl == 0 ){
    Igidl = 0.0e0 ;
    Igidl_dVbs = 0.0e0 ;
    Igidl_dVds = 0.0e0 ;
    Igidl_dVgs = 0.0e0 ;
  } else {
    T3 = here->HSM2_2qnsub_esi ;
    Qb0 = sqrt (T3 * (Pb20 - Vbsz2)) ;
    T4 = 0.5 * T3 / Qb0 ;
    Qb0_dVb = T4 * (- Vbsz2_dVbs) ;
    Qb0_dVd = T4 * (- Vbsz2_dVds) ;
    Qb0_dVg = T4 * (- Vbsz2_dVgs) ;

    Qb0Cox = model->HSM2_gidl6 * Qb0 * Cox_inv ;
    Qb0Cox_dVb = model->HSM2_gidl6 * ( Qb0_dVb * Cox_inv + Qb0 * Cox_inv_dVb ) ;
    Qb0Cox_dVd = model->HSM2_gidl6 * ( Qb0_dVd * Cox_inv + Qb0 * Cox_inv_dVd ) ;
    Qb0Cox_dVg = model->HSM2_gidl6 * ( Qb0_dVg * Cox_inv + Qb0 * Cox_inv_dVg ) ;

    T1 = model->HSM2_gidl3 * (Vds  + model->HSM2_gidl4) - Vgs + (dVthSC + dVthLP) * model->HSM2_gidl5 - Qb0Cox ;

    T2 = 1.0 / Tox0 ;
    E1 = T1 * T2 ;

    E1_dVb = ( model->HSM2_gidl5 * (dVthSC_dVb + dVthLP_dVb) - Qb0Cox_dVb ) * T2 ;
    E1_dVd = ( model->HSM2_gidl3 + model->HSM2_gidl5 * (dVthSC_dVd + dVthLP_dVd) - Qb0Cox_dVd ) * T2 ;
    E1_dVg = ( -1.0 + model->HSM2_gidl5 * (dVthSC_dVg + dVthLP_dVg) - Qb0Cox_dVg ) * T2 ;

    Fn_SZ( Egidl , E1, eef_dlt, T5) ;
    Egidl_dVb = T5 * E1_dVb ;
    Egidl_dVd = T5 * E1_dVd ;
    Egidl_dVg = T5 * E1_dVg ;
    Egidl    += small ;

    T6 = pow(Egidl,model->HSM2_gidl7) ;
    T1 = model->HSM2_gidl7 * pow(Egidl,model->HSM2_gidl7 -1.0) ;
    T6_dVb = Egidl_dVb * T1 ;
    T6_dVd = Egidl_dVd * T1 ;
    T6_dVg = Egidl_dVg * T1 ;

    T0 = - pParam->HSM2_gidl2 * Egp32 / T6 ;
    if ( T0 < - EXP_THR ) {
      Igidl = 0.0 ;
      Igidl_dVbs = Igidl_dVds = Igidl_dVgs = 0.0 ;
    } else {
      T1 = exp ( T0 ) ;    
      T2 = pParam->HSM2_gidl1 / Egp12 * C_QE * here->HSM2_weff_nf ;
      Igidl = T2 * Egidl * Egidl * T1 ;
      T3 = T2 * T1 * Egidl * (2.0 + pParam->HSM2_gidl2 * Egp32 * Egidl / T6 / T6 ) ;
      Igidl_dVbs = T3 * Egidl_dVb ;
      Igidl_dVds = T3 * Egidl_dVd ;
      Igidl_dVgs = T3 * Egidl_dVg ;
    }

    Vdb = Vds - Vbs ;
    if ( Vdb > 0.0 ) {
      T2 = Vdb * Vdb ;
      T4 = T2 * Vdb ;
      T0 = T4 + C_gidl_delta ;
      T5 = T4 / T0 ;
      T7 = ( 3.0 * T2 * T0 - T4 * 3.0 * T2 ) / ( T0 * T0 ) ;  /* == T5_dVdb */
      Igidl_dVbs = Igidl_dVbs * T5 + Igidl * T7 * ( - 1.0 ) ; /* Vdb_dVbs = -1 */
      Igidl_dVds = Igidl_dVds * T5 + Igidl * T7 * ( + 1.0 ) ; /* Vdb_dVds = +1 */
      Igidl_dVgs = Igidl_dVgs * T5 ; /* Vdb_dVgs = 0 */
      Igidl *= T5 ;
    } else {
      Igidl = 0.0 ;
      Igidl_dVbs = Igidl_dVds = Igidl_dVgs = 0.0 ;
    }

  }


  /*-----------------------------------------------------------*
   * Igisl : GISL
   *-----------------*/
  if( model->HSM2_cogidl == 0){
    Igisl = 0.0e0 ;
    Igisl_dVbs = 0.0e0 ;
    Igisl_dVds = 0.0e0 ;
    Igisl_dVgs = 0.0e0 ;
  } else {
    T1 = model->HSM2_gidl3 * ( - Vds + model->HSM2_gidl4 )
      - ( Vgs - Vds ) + ( dVthSC + dVthLP ) * model->HSM2_gidl5 - Qb0Cox ;

    T2 = 1.0 / Tox0 ;
    E1 = T1 * T2 ;
    E1_dVb = ( model->HSM2_gidl5 * (dVthSC_dVb + dVthLP_dVb) - Qb0Cox_dVb ) * T2 ;
    E1_dVd = ( -1.0 * model->HSM2_gidl3 + 1.0 + model->HSM2_gidl5 * (dVthSC_dVd + dVthLP_dVd) - Qb0Cox_dVd ) * T2 ;
    E1_dVg = ( -1.0 + model->HSM2_gidl5 * (dVthSC_dVg + dVthLP_dVg) - Qb0Cox_dVg ) * T2 ;

    Fn_SZ( Egisl , E1, eef_dlt, T5) ;
    Egisl_dVb = T5 * E1_dVb ;
    Egisl_dVd = T5 * E1_dVd ;
    Egisl_dVg = T5 * E1_dVg ;
    Egisl    += small ;

    T6 = pow(Egisl,model->HSM2_gidl7) ;
    T1 = model->HSM2_gidl7 * pow(Egisl,model->HSM2_gidl7 -1.0) ;
    T6_dVb = Egisl_dVb * T1 ;
    T6_dVd = Egisl_dVd * T1 ;
    T6_dVg = Egisl_dVg * T1 ;

    T0 = - pParam->HSM2_gidl2 * Egp32 / T6 ;
    if ( T0 < - EXP_THR ) {
      Igisl = 0.0 ;
      Igisl_dVbs = Igisl_dVds = Igisl_dVgs = 0.0 ;
    } else {
      T1 = exp ( T0 ) ;
      T2 = pParam->HSM2_gidl1 / Egp12 * C_QE * here->HSM2_weff_nf ;
      Igisl = T2 * Egisl * Egisl * T1 ;
      T3 = T2 * T1 * Egisl * (2.0 + pParam->HSM2_gidl2 * Egp32 * Egisl / T6 / T6) ;
      Igisl_dVbs = T3 * Egisl_dVb ;
      Igisl_dVds = T3 * Egisl_dVd ;
      Igisl_dVgs = T3 * Egisl_dVg ;
    }

    Vsb = - Vbs ;
    if ( Vsb > 0.0 ) {
      T2 = Vsb * Vsb ;
      T4 = T2 * Vsb ;
      T0 = T4 + C_gidl_delta ;
      T5 = T4 / T0 ;
      T7 = ( 3.0 * T2 * T0 - T4 * 3.0 * T2 ) / ( T0 * T0 ) ;  /* == T5_dVsb */
      Igisl_dVbs = Igisl_dVbs * T5 + Igisl * T7 * ( - 1.0 ) ; /* Vsb_dVbs = -1 */
      Igisl_dVds = Igisl_dVds * T5 ; /* Vsb_dVds = 0 */
      Igisl_dVgs = Igisl_dVgs * T5 ; /* Vsb_dVgs = 0 */
      Igisl *= T5 ;
    } else {
      Igisl = 0.0 ;
      Igisl_dVbs = Igisl_dVds = Igisl_dVgs = 0.0 ;
    }

  }


  /*-----------------------------------------------------------*
   * End of PART-2. (label) 
   *-----------------*/
/* end_of_part_2: */

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-3: Overlap charge
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
  Aclm = pParam->HSM2_clm1 ;
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
    
    if (model->HSM2_xqy !=0) {
      Ec = 0.0e0 ;
      Ec_dVbs =0.0e0 ;
      Ec_dVds =0.0e0 ;
      Ec_dVgs =0.0e0 ;
    }
  } else {
    /* Ec is removed from Lred calc. part */
    if (model->HSM2_xqy !=0) {
      if ( Idd < C_IDD_MIN ) {
        Ec = 0.0e0 ;
        Ec_dVbs =0.0e0 ;
        Ec_dVds =0.0e0 ;
        Ec_dVgs =0.0e0 ;
      } else {
      T1 =  beta_inv / Leff ;
      T2 = 1.0 / Qn0 ;
      T3 = T2 * T2 ;
      Ec = Idd * T1 * T2 ;
      Ec_dVbs = T1 * (Idd_dVbs * T2 - Idd * Qn0_dVbs * T3 ) ;
      Ec_dVds = T1 * (Idd_dVds * T2 - Idd * Qn0_dVds * T3 ) ;
      Ec_dVgs = T1 * (Idd_dVgs * T2 - Idd * Qn0_dVgs * T3 ) ;
      }
    }
  }


  /*-------------------------------------------*
   * Overlap charges: Qgod, Qgos, and Qover 
   *-----------------*/
  if ( model->HSM2_coovlp >= 1 && pParam->HSM2_lover > 0.0 ){
    cov_slp = modelMKS->HSM2_ovslp ;
    cov_mag = model->HSM2_ovmag ;

    covvg = Vgs ;
    covvg_dVgs = 1.0 ;

    T1 = Cox0 * here->HSM2_weff_nf ;
    Lov = pParam->HSM2_lover ;

    if ( pParam->HSM2_nover == 0.0 ){
      T4 = cov_slp * T1 * ( cov_mag + covvg ) ;
      T4_dVg = cov_slp * T1 * covvg_dVgs ;
      T4_dVd = 0.0 ;
      T5 = Lov * T1 ;
     
      TX = Ps0 ;
      TX_dVbs = Ps0_dVbs ;
      TX_dVds = Ps0_dVds ;
      TX_dVgs = Ps0_dVgs ;

      T9 = 1.2e0 - TX ;
      Qgos = Vgs * T5 - T9 * T4 ;
      Qgos_dVbs =      T4 * TX_dVbs ;
      Qgos_dVds =      T4 * TX_dVds - T9 * T4_dVd ;
      Qgos_dVgs = T5 + T4 * TX_dVgs - T9 * T4_dVg ;

      T4 = cov_slp * T1 * ( cov_mag + covvg - Vds ) ;
      T4_dVg = cov_slp * T1 * covvg_dVgs ;
      T4_dVd = - cov_slp * T1 ;
      TX = Psl - Vds ;
      TX_dVbs = Psl_dVbs ;
      TX_dVds = Psl_dVds - 1.0 ;
      TX_dVgs = Psl_dVgs ;
      T9 = 1.2e0 - TX ;
      Qgod = (Vgs - Vds) * T5 - T4 * T9 ;
      Qgod_dVbs =        T4 * TX_dVbs ;
      Qgod_dVds = - T5 + T4 * TX_dVds - T9 * T4_dVd ;
      Qgod_dVgs =   T5 + T4 * TX_dVgs - T9 * T4_dVg ;
    } else {
      for ( lcover = -1 ; lcover <= +1 ; lcover += 2 ) {
        flg_ovloops = ( 1 - lcover ) / 2 ; /* 1 in Source overlap calc. */
        flg_ovloopd = ( 1 + lcover ) / 2 ; /* 1 in Drain overlap calc. */
        /*-------------------------------------------*
         * Qover (G/D overlap charge)  |  note: _dVxs means _dVxse 
         *------------------------*/
        Vbsgmt = ModeNML * Vbse + ModeRVS * ( Vbse - Vdse ) ;
        Vdsgmt = ModeNML * Vdse + ModeRVS * ( - Vdse ) ; 
        Vgsgmt = ModeNML * Vgse + ModeRVS * ( Vgse - Vdse ) ;
        Vdbgmt = Vdsgmt - Vbsgmt ;
        Vgbgmt = Vgsgmt - Vbsgmt ;
        Vsbgmt = - Vbsgmt ;
        flg_overs = flg_ovloops * (int)ModeNML + flg_ovloopd * (int)ModeRVS ; /* geometrical source */
        flg_overd = flg_ovloops * (int)ModeRVS + flg_ovloopd * (int)ModeNML ; /* geometrical drain */
        Vxbgmt = flg_overs * Vsbgmt + flg_overd * Vdbgmt + epsm10 ;

    
        /*---------------------------------------------------*
         * Clamp -Vxbgmt.
         *-----------------*/
        T0 = - Vxbgmt;
        if ( T0 > Vbs_bnd ) {
          T1 =    T0   - Vbs_bnd;
          T2 =    Vbs_max    - Vbs_bnd;
    
          Fn_SUPoly4( TY, T1, T2, T11 );
    
          T10 = Vbs_bnd + TY ;
        }  else {
          T10 = T0 ;
          T11 = 1.0 ;
        }
        Vxbgmtcl = - T10 - small2 ;
        Vxbgmtcl_dVxbgmt = T11;
    
        fac1 = here->HSM2_cnst0over * Cox0_inv ;
        fac1_dVbs = 0.0; fac1_dVds = 0.0; fac1_dVgs = 0.0;
    
        fac1p2 = fac1 * fac1 ;
      
        VgpLD = - Vgbgmt + pParam->HSM2_vfbover;  
        VgpLD_dVgb = - 1.0e0 ;
    
        T0 = pParam->HSM2_nover / here->HSM2_nin ;
        Pb2over = 2.0 / beta * log( T0 ) ;
    
        Vgb_fb_LD =  - Vxbgmtcl ;
    
        /*-----------------------------------*
         * QsuLD: total charge = Accumulation | Depletion+inversion
         *-----------------*/
        if (   VgpLD  < Vgb_fb_LD ){   
          /*---------------------------*
           * Accumulation
           *-----------------*/
          flg_ovzone = -1 ; 
          T1 = 1.0 / ( beta * here->HSM2_cnst0over ) ;
          TY = T1 * Cox0 ;
          Ac41 = 2.0 + 3.0 * C_SQRT_2 * TY ;
          Ac4 = 8.0 * Ac41 * Ac41 * Ac41 ;
      
          Ps0_min = here->HSM2_eg - Pb2over ;
      
          TX = beta * ( VgpLD + Vxbgmtcl ) ;
          TX_dVxb = beta * Vxbgmtcl_dVxbgmt ;
          TX_dVgb = beta * VgpLD_dVgb ;
      
          Ac31 = 7.0 * C_SQRT_2 - 9.0 * TY * ( TX - 2.0 ) ;
          Ac31_dVxb = - 9.0 * TY * TX_dVxb ;
          Ac31_dVgb = - 9.0 * TY * TX_dVgb ;
      
          Ac3 = Ac31 * Ac31 ;
          T1 = 2.0 * Ac31 ;
          Ac3_dVxb = T1 * Ac31_dVxb ;
          Ac3_dVgb = T1 * Ac31_dVgb ;
      
          Ac2 = sqrt( Ac4 + Ac3 ) ;
          T1 = 0.5 / Ac2 ;
          Ac2_dVxb = T1 *  Ac3_dVxb ;
          Ac2_dVgb = T1 *  Ac3_dVgb ;
        
          Ac1 = -7.0 * C_SQRT_2 + Ac2 + 9.0 * TY * ( TX - 2.0 ) ;
          Ac1_dVxb = Ac2_dVxb + 9.0 * TY * TX_dVxb ;
          Ac1_dVgb = Ac2_dVgb + 9.0 * TY * TX_dVgb ;
      
          Acd = pow( Ac1 , C_1o3 ) ;
          T1 = C_1o3 / ( Acd * Acd ) ;
          Acd_dVxb = Ac1_dVxb * T1 ;
          Acd_dVgb = Ac1_dVgb * T1 ;
      
          Acn = -4.0 * C_SQRT_2 - 12.0 * TY + 2.0 * Acd + C_SQRT_2 * Acd * Acd ;
          T1 = 2.0 + 2.0 * C_SQRT_2 * Acd ;
          Acn_dVxb = T1 * Acd_dVxb ;
          Acn_dVgb = T1 * Acd_dVgb ;
       
          Chi = Acn / Acd ;
          T1 = 1.0 / ( Acd * Acd ) ;
          Chi_dVxb = ( Acn_dVxb * Acd - Acn * Acd_dVxb ) * T1 ;
          Chi_dVgb = ( Acn_dVgb * Acd - Acn * Acd_dVgb ) * T1 ;
      
          Psa = Chi * beta_inv - Vxbgmtcl ;
          Psa_dVxb = Chi_dVxb * beta_inv - Vxbgmtcl_dVxbgmt ;
          Psa_dVgb = Chi_dVgb * beta_inv ;
      
          T1 = Psa + Vxbgmtcl ;
          T2 = T1 / Ps0_min ;
          T3 = sqrt( 1.0 + ( T2 * T2 ) ) ;
      
          T9 = T2 / T3 / Ps0_min ;
          T3_dVd = T9 * ( Psa_dVxb + Vxbgmtcl_dVxbgmt ) ;
          T3_dVg = T9 * Psa_dVgb ;
    
          Ps0LD = T1 / T3 - Vxbgmtcl ;
          T9 = 1.0 / ( T3 * T3 ) ;
          Ps0LD_dVxb = T9 * ( ( Psa_dVxb + Vxbgmtcl_dVxbgmt ) * T3 - T1 * T3_dVd ) - Vxbgmtcl_dVxbgmt ;
          Ps0LD_dVgb = T9 * ( Psa_dVgb * T3 - T1 * T3_dVg );
         
          T2 = ( VgpLD - Ps0LD ) ;
          QsuLD = Cox0 * T2 ;
          QsuLD_dVxb = - Cox0 * Ps0LD_dVxb ;
          QsuLD_dVgb = Cox0 * ( VgpLD_dVgb - Ps0LD_dVgb ) ;
      
          QbuLD = QsuLD ;
          QbuLD_dVxb = QsuLD_dVxb ;
          QbuLD_dVgb = QsuLD_dVgb ;
      
        } else {
      
          /*---------------------------*
           * Depletion and inversion
           *-----------------*/
    
          /* initial value for a few fixpoint iterations
             to get Ps0_iniA from simplified Poisson equation: */
           flg_ovzone = 2 ;
           Chi = znbd3 ;
           Chi_dVxb = 0.0 ; Chi_dVgb = 0.0 ;

           Ps0_iniA= Chi/beta - Vxbgmtcl ;
           Ps0_iniA_dVxb = Chi_dVxb/beta - Vxbgmtcl_dVxbgmt ;
           Ps0_iniA_dVgb = Chi_dVgb/beta ;
          
          /* 1 .. 2 relaxation steps should be sufficient */
          for ( lp_ld = 1; lp_ld <= 2; lp_ld ++ ) {
            TY = exp(-Chi);
            TY_dVxb = -Chi_dVxb * TY;
            TY_dVgb = -Chi_dVgb * TY;
            TX = 1.0e0 + 4.0e0 
               * ( beta * ( VgpLD + Vxbgmtcl ) - 1.0e0 + TY ) / ( fac1p2 * beta2 ) ;
            TX_dVxb = 4.0e0 * ( beta * ( Vxbgmtcl_dVxbgmt ) + TY_dVxb ) / ( fac1p2 * beta2 );
            TX_dVgb = 4.0e0 * ( beta * ( VgpLD_dVgb       ) + TY_dVgb ) / ( fac1p2 * beta2 );
            T1 = ( beta * ( VgpLD + Vxbgmtcl ) - 1.0e0 + TY );
            T3 = fac1p2 * beta2 ;
            if ( TX < epsm10) {
              TX = epsm10; TX_dVxb = 0.0; TX_dVgb = 0.0;
            }
    
            Ps0_iniA = VgpLD + fac1p2 * beta / 2.0e0 * ( 1.0e0 - sqrt( TX ) ) ;
            Ps0_iniA_dVxb =            - fac1p2 * beta / 2.0e0 * TX_dVxb * 0.5 / sqrt( TX );
            Ps0_iniA_dVgb = VgpLD_dVgb - fac1p2 * beta / 2.0e0 * TX_dVgb * 0.5 / sqrt( TX );
            T1 = fac1p2 * beta ;
            T2 = 1.0 - sqrt( TX );
      
            Chi = beta * ( Ps0_iniA + Vxbgmtcl ) ;
            Chi_dVxb = beta * ( Ps0_iniA_dVxb + Vxbgmtcl_dVxbgmt ) ;
            Chi_dVgb = beta * ( Ps0_iniA_dVgb ) ;
          } /* End of iteration */
    
          if ( Chi < znbd3 ) { 
    
            flg_ovzone = 1 ; 
    
            /*-----------------------------------*
             * zone-D1
             * - Ps0_iniA is the analytical solution of QovLD=Qb0 with
             *   Qb0 being approximated by 3-degree polynomial.
             *
             *   new: Inclusion of exp(-Chi) term at right border
             *-----------------*/
            Ta =  1.0/(9.0*sqrt(2.0)) - (5.0+7.0*exp(-3.0)) / (54.0*sqrt(2.0+exp(-3.0)));
            Tb = (1.0+exp(-3.0)) / (2.0*sqrt(2.0+exp(-3.0))) - sqrt(2.0) / 3.0;
            Tc =  1.0/sqrt(2.0) + 1.0/(beta*fac1);
            Td = - (VgpLD + Vxbgmtcl) / fac1;
            Td_dVxb = - Vxbgmtcl_dVxbgmt / fac1;
            Td_dVgb = - VgpLD_dVgb / fac1;
            Tq = Tb*Tb*Tb / (27.0*Ta*Ta*Ta) - Tb*Tc/(6.0*Ta*Ta) + Td/(2.0*Ta);
            Tq_dVxb = Td_dVxb/(2.0*Ta);
            Tq_dVgb = Td_dVgb / (2.0*Ta);
            Tp = (3.0*Ta*Tc-Tb*Tb)/(9.0*Ta*Ta);
            T5      = sqrt(Tq*Tq + Tp*Tp*Tp);
            T5_dVxb = 2.0*Tq*Tq_dVxb / (2.0*T5);
            T5_dVgb = 2.0*Tq*Tq_dVgb / (2.0*T5);
            Tu = pow(-Tq + T5,C_1o3);
            Tu_dVxb = Tu / (3.0 * (-Tq + T5)) * (-Tq_dVxb + T5_dVxb);
            Tu_dVgb = Tu / (3.0 * (-Tq + T5)) * (-Tq_dVgb + T5_dVgb);
            Tv = -pow(Tq + T5,C_1o3);
            Tv_dVxb = Tv / (3.0 * (-Tq - T5)) * (-Tq_dVxb - T5_dVxb);
            Tv_dVgb = Tv / (3.0 * (-Tq - T5)) * (-Tq_dVgb - T5_dVgb);
            TX      = Tu + Tv - Tb/(3.0*Ta);
            TX_dVxb = Tu_dVxb + Tv_dVxb;
            TX_dVgb = Tu_dVgb + Tv_dVgb;
            
            Ps0_iniA = TX * beta_inv - Vxbgmtcl ;
            Ps0_iniA_dVxb = TX_dVxb * beta_inv - Vxbgmtcl_dVxbgmt;
            Ps0_iniA_dVgb = TX_dVgb * beta_inv;
    
            Chi = beta * ( Ps0_iniA + Vxbgmtcl ) ;
            Chi_dVxb = beta * ( Ps0_iniA_dVxb + Vxbgmtcl_dVxbgmt ) ;
            Chi_dVgb = beta * ( Ps0_iniA_dVgb ) ;
          }
    
          if ( model->HSM2_coqovsm > 0 ) {
    	  /*-----------------------------------*
    	   * - Ps0_iniB : upper bound.
    	   *-----------------*/
            flg_ovzone += 2;
    
            VgpLD_shift = VgpLD + Vxbgmtcl + 0.1;
            VgpLD_shift_dVgb = VgpLD_dVgb;
            VgpLD_shift_dVxb = Vxbgmtcl_dVxbgmt;
            exp_bVbs = exp( beta * - Vxbgmtcl ) + small ;
            exp_bVbs_dVxb = - exp_bVbs * beta * Vxbgmtcl_dVxbgmt;
            T0 = here->HSM2_nin / pParam->HSM2_nover;
            cnst1over = T0 * T0;
            gamma = cnst1over * exp_bVbs ;
            gamma_dVxb = cnst1over * exp_bVbs_dVxb;
         
            T0    = beta2 * fac1p2;
    
            psi = beta*VgpLD_shift;
            psi_dVgb = beta*VgpLD_shift_dVgb;
            psi_dVxb = beta*VgpLD_shift_dVxb;
            Chi_1      = log(gamma*T0 + psi*psi) - log(cnst1over*T0) + beta*Vxbgmtcl;
            Chi_1_dVgb = 2.0*psi*psi_dVgb/ (gamma*T0 + psi*psi);
            Chi_1_dVxb = (gamma_dVxb*T0+2.0*psi*psi_dVxb)/(gamma*T0+psi*psi)
                                + beta*Vxbgmtcl_dVxbgmt;    
    
            Fn_SU2( Chi_1, Chi_1, psi, 1.0, T1, T2 );
            Chi_1_dVgb = Chi_1_dVgb*T1 + psi_dVgb*T2;
            Chi_1_dVxb = Chi_1_dVxb*T1 + psi_dVxb*T2;
    
         /* 1 fixpoint step for getting more accurate Chi_B */
            psi      -= Chi_1 ;
            psi_dVgb -= Chi_1_dVgb ;
            psi_dVxb -= Chi_1_dVxb ;
         
            psi      += beta*0.1 ;
/*    
            psi_B = psi;
            arg_B = psi*psi/(gamma*T0);*/
            Chi_B = log(gamma*T0 + psi*psi) - log(cnst1over*T0) + beta*Vxbgmtcl;
            Chi_B_dVgb = 2.0*psi*psi_dVgb/ (gamma*T0 + psi*psi);
            Chi_B_dVxb = (gamma_dVxb*T0+2.0*psi*psi_dVxb)/(gamma*T0+psi*psi)
                                + beta*Vxbgmtcl_dVxbgmt;    
            Ps0_iniB      = Chi_B/beta - Vxbgmtcl ;
/*            Ps0_iniB_dVgb = Chi_B_dVgb/beta;
            Ps0_iniB_dVxb = Chi_B_dVxb/beta- Vxbgmtcl_dVxbgmt;
*/    
            
            /* construction of Ps0LD by taking Ps0_iniB as an upper limit of Ps0_iniA
             *
             * Limiting is done for Chi rather than for Ps0LD, to avoid shifting
             * for Fn_SU2 */
    
            Chi_A = Chi;
            Chi_A_dVxb = Chi_dVxb;
            Chi_A_dVgb = Chi_dVgb;
    
            Fn_SU2( Chi, Chi_A, Chi_B, c_ps0ini_2*75.00, T1, T2 ); /* org: 50 */
            Chi_dVgb = Chi_A_dVgb * T1 + Chi_B_dVgb * T2;  
            Chi_dVxb = Chi_A_dVxb * T1 + Chi_B_dVxb * T2;
    
          }
    
            /* updating Ps0LD */
            Ps0LD= Chi/beta - Vxbgmtcl ;
            Ps0LD_dVgb = Chi_dVgb/beta;
            Ps0LD_dVxb = Chi_dVxb/beta- Vxbgmtcl_dVxbgmt;
    
          T1      = Chi - 1.0 + exp(-Chi);
          T1_dVxb = (1.0 - exp(-Chi)) * Chi_dVxb ;
          T1_dVgb = (1.0 - exp(-Chi)) * Chi_dVgb ;
          if (T1 < epsm10) {
             T1 = epsm10 ;
             T1_dVxb = 0.0 ;
             T1_dVgb = 0.0 ;
          }
          T2 = sqrt(T1);
          QbuLD = here->HSM2_cnst0over * T2 ;
          T3 = here->HSM2_cnst0over * 0.5 / T2 ;
          QbuLD_dVxb = T3 * T1_dVxb ;
          QbuLD_dVgb = T3 * T1_dVgb ;
         
          /*-----------------------------------------------------------*
           * QsuLD : Qovs or Qovd in unit area.
           * note: QsuLD = Qdep+Qinv. 
           *-----------------*/
          QsuLD = Cox0 * ( VgpLD - Ps0LD ) ;
          QsuLD_dVxb = Cox0 * ( - Ps0LD_dVxb ) ;
          QsuLD_dVgb = Cox0 * ( VgpLD_dVgb - Ps0LD_dVgb ) ;
    
          if ( model->HSM2_coqovsm == 1 ) { /* take initial values from analytical model */ 
    
       
            /*---------------------------------------------------*
             * Calculation of Ps0LD. (beginning of Newton loop) 
             * - Fs0 : Fs0 = 0 is the equation to be solved. 
             * - dPs0 : correction value. 
             *-----------------*/
    
            /* initial value too close to flat band should not be used */
            exp_bVbs = exp( beta * - Vxbgmtcl ) ;
            T0 = here->HSM2_nin / pParam->HSM2_nover;
            cnst1over = T0 * T0;
            cfs1 = cnst1over * exp_bVbs ;
        
            flg_conv = 0 ;
            for ( lp_s0 = 1 ; lp_s0 <= 2*lp_s0_max + 1 ; lp_s0 ++ ) { 
    
                Chi = beta * ( Ps0LD + Vxbgmtcl ) ;
       
                if ( Chi < znbd5 ) { 
                  /*-------------------------------------------*
                   * zone-D1/D2. (Ps0LD)
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
    
                  fs02 = sqrt( fb * fb + fs01 + small ) ;
                  fs02_dPs0 = ( beta * fb_dChi * 2 * fb + fs01_dPs0 ) / ( fs02 + fs02 ) ;
    
                } else {
                 /*-------------------------------------------*
                  * zone-D3. (Ps0LD)
                  *-----------------*/
                 if ( Chi < large_arg ) { /* avoid exp_Chi to become extremely large */
    	        exp_Chi = exp( Chi ) ;
    	        fs01 = cfs1 * ( exp_Chi - 1.0e0 ) ;
    	        fs01_dPs0 = cfs1 * beta * ( exp_Chi ) ;
                 } else {
                    exp_bPs0 = exp( beta*Ps0LD ) ;
                    fs01     = cnst1over * ( exp_bPs0 - exp_bVbs ) ;
                    fs01_dPs0 = cnst1over * beta * exp_bPs0 ;
                 }
                 fs02 = sqrt( Chi - 1.0 + fs01 ) ;
                 fs02_dPs0 = ( beta + fs01_dPs0 ) / fs02 * 0.5 ;
       
                } /* end of if ( Chi  ... ) block */
                /*-----------------------------------------------------------*
                 * Fs0
                 *-----------------*/
                Fs0 = VgpLD - Ps0LD - fac1 * fs02 ;
                Fs0_dPs0 = - 1.0e0 - fac1 * fs02_dPs0 ;
    
                if ( flg_conv == 1 ) break ;
    
                dPs0 = - Fs0 / Fs0_dPs0 ;
    
                /*-------------------------------------------*
                 * Update Ps0LD .
                 *-----------------*/
                dPlim = 0.5*dP_max*(1.0 + Fn_Max(1.e0,fabs(Ps0LD))) ;
                if ( fabs( dPs0 ) > dPlim ) dPs0 = dPlim * Fn_Sgn( dPs0 ) ;
    
                Ps0LD = Ps0LD + dPs0 ;
    
                TX = -Vxbgmtcl + ps_conv / 2 ;
                if ( Ps0LD < TX ) Ps0LD = TX ;
          
                /*-------------------------------------------*
                 * Check convergence. 
                 *-----------------*/
                if ( fabs( dPs0 ) <= ps_conv && fabs( Fs0 ) <= gs_conv ) {
                  flg_conv = 1 ;
                }
          
            } /* end of Ps0LD Newton loop */
    
            /*-------------------------------------------*
             * Procedure for diverged case.
             *-----------------*/
            if ( flg_conv == 0 ) { 
              fprintf( stderr , 
                       "*** warning(HiSIM_HV): Went Over Iteration Maximum (Ps0LD)\n" ) ;
              fprintf( stderr , " -Vxbgmtcl = %e   Vgbgmt = %e\n" , -Vxbgmtcl , Vgbgmt ) ;
            } 
    
            /*---------------------------------------------------*
             * Evaluate derivatives of Ps0LD. 
             *-----------------*/
    
            if ( Chi < znbd5 ) { 
              fs01_dVbs = cfs1 * beta * fi * ( - fi + 2 * fi_dChi ) ; /* fs01_dVxbgmtcl */
              T2 = 1.0e0 / ( fs02 + fs02 ) ;
              fs02_dVbs = ( + beta * fb_dChi * 2 * fb + fs01_dVbs ) * T2 ; /* fs02_dVxbgmtcl */
            } else {
              if ( Chi < large_arg ) {
                fs01_dVbs = + cfs1 * beta ; /* fs01_dVxbgmtcl */
              } else {
                fs01_dVbs   = + cfs1 * beta ;
              }
              T2 = 0.5e0 / fs02 ;
              fs02_dVbs = ( + beta + fs01_dVbs ) * T2 ; /* fs02_dVxbgmtcl */
            }
    
            T1 = 1.0 / Fs0_dPs0 ;
            Ps0LD_dVxb = - ( - fac1 * fs02_dVbs ) * T1 ;
            Ps0LD_dVds = 0.0 ;
            Ps0LD_dVgb = - ( VgpLD_dVgb - fac1_dVgs * fs02 ) * T1 ;
    
    
            if ( Chi < znbd5 ) { 
              /*-------------------------------------------*
               * zone-D1/D2. (Ps0LD)
               *-----------------*/
              if ( Chi < znbd3 ) { flg_ovzone = 1; }
                            else { flg_ovzone = 2; }
    
              Xi0 = fb * fb + epsm10 ;
              T1 = 2 * fb * fb_dChi * beta ;
              Xi0_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0_dVxbgmtcl */
              Xi0_dVds = T1 * Ps0LD_dVds ;
              Xi0_dVgs = T1 * Ps0LD_dVgb ;
    
              Xi0p12 = fb + epsm10 ;
              T1 = fb_dChi * beta ;
              Xi0p12_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0p12_dVxbgmtcl */
              Xi0p12_dVds = T1 * Ps0LD_dVds ;
              Xi0p12_dVgs = T1 * Ps0LD_dVgb ;
    
              Xi0p32 = fb * fb * fb + epsm10 ;
              T1 = 3 * fb * fb * fb_dChi * beta ;
              Xi0p32_dVbs = T1 * ( Ps0LD_dVxb + 1.0 ) ; /* Xi0p32_dVxbgmtcl */
              Xi0p32_dVds = T1 * Ps0LD_dVds ;
              Xi0p32_dVgs = T1 * Ps0LD_dVgb ;
     
            } else { 
              /*-------------------------------------------*
               * zone-D3. (Ps0LD)
               *-----------------*/
              flg_ovzone = 3 ;
    
              Xi0 = Chi - 1.0e0 ;
              Xi0_dVbs = beta * ( Ps0LD_dVxb + 1.0e0 ) ; /* Xi0_dVxbgmtcl */
              Xi0_dVds = beta * Ps0LD_dVds ;
              Xi0_dVgs = beta * Ps0LD_dVgb ;
     
              Xi0p12 = sqrt( Xi0 ) ;
              T1 = 0.5e0 / Xi0p12 ;
              Xi0p12_dVbs = T1 * Xi0_dVbs ;
              Xi0p12_dVds = T1 * Xi0_dVds ;
              Xi0p12_dVgs = T1 * Xi0_dVgs ;
    
              Xi0p32 = Xi0 * Xi0p12 ;
              T1 = 1.5e0 * Xi0p12 ;
              Xi0p32_dVbs = T1 * Xi0_dVbs ;
              Xi0p32_dVds = T1 * Xi0_dVds ;
              Xi0p32_dVgs = T1 * Xi0_dVgs ;
    
              if ( Chi < large_arg ) {
              } else {
              }
            } /* end of if ( Chi  ... ) block */
        
            /*-----------------------------------------------------------*
             * - Recalculate the derivatives of fs01 and fs02.
             *-----------------*/
            fs01_dVbs = Ps0LD_dVxb * fs01_dPs0 + fs01_dVbs ;
            fs01_dVds = Ps0LD_dVds * fs01_dPs0 ;
            fs01_dVgs = Ps0LD_dVgb * fs01_dPs0 ;
            fs02_dVbs = Ps0LD_dVxb * fs02_dPs0 + fs02_dVbs ;
/*            fs02_dVxb = Ps0LD_dVds * fs02_dPs0 ;*/
            fs02_dVgb = Ps0LD_dVgb * fs02_dPs0 ;
    
            /*-----------------------------------------------------------*
             * QbuLD and QiuLD
             *-----------------*/
            QbuLD = here->HSM2_cnst0over * Xi0p12 ;
            QbuLD_dVxb = here->HSM2_cnst0over * Xi0p12_dVbs ;
            QbuLD_dVgb = here->HSM2_cnst0over * Xi0p12_dVgs ;
    
            T1 = 1.0 / ( fs02 + Xi0p12 ) ;
            QiuLD = here->HSM2_cnst0over * fs01 * T1 ;
            T2 = 1.0 / ( fs01 + epsm10 ) ;
            QiuLD_dVbs = QiuLD * ( fs01_dVbs * T2 - ( fs02_dVbs + Xi0p12_dVbs ) * T1 ) ;
            QiuLD_dVgs = QiuLD * ( fs01_dVgs * T2 - ( fs02_dVgb + Xi0p12_dVgs ) * T1 ) ;
    
            /*-----------------------------------------------------------*
             * Extrapolation: X_dVxbgmt = X_dVxbgmtcl * Vxbgmtcl_dVxbgmt
             *-----------------*/
            QbuLD_dVxb *= Vxbgmtcl_dVxbgmt ;
            QiuLD_dVbs *= Vxbgmtcl_dVxbgmt ;
    
            /*-----------------------------------------------------------*
             * Total overlap charge
             *-----------------*/
            QsuLD = QbuLD + QiuLD;
            QsuLD_dVxb = QbuLD_dVxb + QiuLD_dVbs;
            QsuLD_dVgb = QbuLD_dVgb + QiuLD_dVgs;
    
          } /* end of COQOVSM branches */
    
        } /* end of Vgbgmt region blocks */
      

        /* convert to source ref. */
        QsuLD_dVbs = - (QsuLD_dVxb + QsuLD_dVgb) ;
        QsuLD_dVds =    QsuLD_dVxb * flg_overd   ;
        QsuLD_dVgs =                 QsuLD_dVgb  ;

        QbuLD_dVbs = - (QbuLD_dVxb + QbuLD_dVgb) ;
        QbuLD_dVds =    QbuLD_dVxb * flg_overd   ;
        QbuLD_dVgs =                 QbuLD_dVgb  ;

        /* inversion charge = total - depletion */
        QiuLD = QsuLD - QbuLD  ;
        QiuLD_dVbs = QsuLD_dVbs - QbuLD_dVbs ;
        QiuLD_dVds = QsuLD_dVds - QbuLD_dVds ;
        QiuLD_dVgs = QsuLD_dVgs - QbuLD_dVgs ;

        /* assign final outputs of Qover model */
        /* note: Qovs and Qovd are exchanged in reverse mode */
        T4 = here->HSM2_weff_nf * Lov ;

        if(flg_ovloops) {
          Qovs =  T4 * QsuLD ;
          Qovs_dVbs = T4 * QsuLD_dVbs ;
          Qovs_dVds = T4 * QsuLD_dVds ;
          Qovs_dVgs = T4 * QsuLD_dVgs ;
          QisLD = T4 * QiuLD ;
          QisLD_dVbs = T4 * QiuLD_dVbs ;
          QisLD_dVds = T4 * QiuLD_dVds ;
          QisLD_dVgs = T4 * QiuLD_dVgs ;
          QbsLD = T4 * QbuLD ;
          QbsLD_dVbs = T4 * QbuLD_dVbs ;
          QbsLD_dVds = T4 * QbuLD_dVds ;
          QbsLD_dVgs = T4 * QbuLD_dVgs ;
        }

        if(flg_ovloopd) {
          Qovd =  T4 * QsuLD ;
          Qovd_dVbs = T4 * QsuLD_dVbs ;
          Qovd_dVds = T4 * QsuLD_dVds ;
          Qovd_dVgs = T4 * QsuLD_dVgs ;
          QidLD = T4 * QiuLD ;
          QidLD_dVbs = T4 * QiuLD_dVbs ;
          QidLD_dVds = T4 * QiuLD_dVds ;
          QidLD_dVgs = T4 * QiuLD_dVgs ;
          QbdLD = T4 * QbuLD ;
          QbdLD_dVbs = T4 * QbuLD_dVbs ;
          QbdLD_dVds = T4 * QbuLD_dVds ;
          QbdLD_dVgs = T4 * QbuLD_dVgs ;
        }
  
      } /* end of lcover loop */

      /* convert to the derivatives w.r.t. mode-dependent biases */
      Qovs_dVds = ModeNML * Qovs_dVds
       - ModeRVS * ( Qovs_dVds + Qovs_dVgs + Qovs_dVbs ) ;
      QisLD_dVds = ModeNML * QisLD_dVds
       - ModeRVS * ( QisLD_dVds + QisLD_dVgs + QisLD_dVbs ) ;
      QbsLD_dVds = ModeNML * QbsLD_dVds
       - ModeRVS * ( QbsLD_dVds + QbsLD_dVgs + QbsLD_dVbs ) ;
      Qovd_dVds = ModeNML * Qovd_dVds
       - ModeRVS * ( Qovd_dVds + Qovd_dVgs + Qovd_dVbs );
      QidLD_dVds = ModeNML * QidLD_dVds
       - ModeRVS * ( QidLD_dVds + QidLD_dVgs + QidLD_dVbs ) ;
      QbdLD_dVds = ModeNML * QbdLD_dVds
       - ModeRVS * ( QbdLD_dVds + QbdLD_dVgs + QbdLD_dVbs ) ;

    } /* end of if ( pParam->HSM2_nover == 0.0 ) */

    /*-----------------------------------*
     * Additional constant capacitance model
     *-----------------*/
    flg_overgiven = ( (int)ModeRVS * model->HSM2_cgso_Given
                    + (int)ModeNML * model->HSM2_cgdo_Given  ) ;
    if ( flg_overgiven ) {
      Cgdo  = ModeRVS * pParam->HSM2_cgso + ModeNML * pParam->HSM2_cgdo ;
      Cgdo *= - here->HSM2_weff_nf ;

      Qgod += - Cgdo * (Vgs - Vds) ;
      Qgod_dVds +=  Cgdo ;
      Qgod_dVgs += -Cgdo ;
    }

    flg_overgiven = ( (int)ModeNML * model->HSM2_cgso_Given
                    + (int)ModeRVS * model->HSM2_cgdo_Given ) ;
    if(flg_overgiven) {
      Cgso  = ModeNML * pParam->HSM2_cgso + ModeRVS * pParam->HSM2_cgdo ;
      Cgso *= - here->HSM2_weff_nf ;

      Qgos += - Cgso * Vgs ;
      Qgos_dVgs += -Cgso ;
    }

  } else { /* else case of if ( model->HSM2_coovlp >= 1 ) */
    if ( (here->HSM2_mode == HiSIM_NORMAL_MODE && !model->HSM2_cgdo_Given) ||
	 (here->HSM2_mode != HiSIM_NORMAL_MODE && !model->HSM2_cgso_Given) ) {
      Cgdo = - Cox0 * pParam->HSM2_lover * here->HSM2_weff_nf ;
    } else {
      Cgdo  = ModeRVS * pParam->HSM2_cgso + ModeNML * pParam->HSM2_cgdo ;
      Cgdo *= - here->HSM2_weff_nf ;
    }
    Qgod = - Cgdo * (Vgs - Vds) ;
    Qgod_dVbs = 0.0 ;
    Qgod_dVds = Cgdo ;
    Qgod_dVgs = - Cgdo ;
    
    if ( (here->HSM2_mode == HiSIM_NORMAL_MODE && !model->HSM2_cgso_Given) ||
         (here->HSM2_mode != HiSIM_NORMAL_MODE && !model->HSM2_cgdo_Given) ) {
      Cgso = - Cox0 * pParam->HSM2_lover * here->HSM2_weff_nf ;
    } else {
      Cgso  = ModeNML * pParam->HSM2_cgso + ModeRVS * pParam->HSM2_cgdo ;
      Cgso *= - here->HSM2_weff_nf ;
    }
    Qgos = - Cgso * Vgs ;
    Qgos_dVbs = 0.0 ;
    Qgos_dVds = 0.0 ;
    Qgos_dVgs = - Cgso ;
  } /* end of if ( model->HSM2_coovlp >= 1 ) */

  /*-------------------------------------------*
   * Gate/Bulk overlap charge: Qgbo
   *-----------------*/
  Cgbo_loc = - model->HSM2_cgbo * here->HSM2_lgate ;
  Qgbo = - Cgbo_loc * (Vgs -Vbs) ;
  Qgbo_dVgs = - Cgbo_loc ;
  Qgbo_dVbs =   Cgbo_loc ;
  Qgbo_dVds = 0.0 ;

  /*---------------------------------------------------*
   * Lateral-field-induced capacitance.
   *-----------------*/
  if ( model->HSM2_coqy == 0 || model->HSM2_xqy == 0 ){
    Qy = 0.0e0 ;
    Qy_dVds = 0.0e0 ;
    Qy_dVgs = 0.0e0 ;
    Qy_dVbs = 0.0e0 ;
  } else {
    Pslk = Ec * Leff + Ps0 ;
    Pslk_dVbs = Ec_dVbs * Leff + Ps0_dVbs;
    Pslk_dVds = Ec_dVds * Leff + Ps0_dVds;
    Pslk_dVgs = Ec_dVgs * Leff + Ps0_dVgs;
    Fn_SU2( T10, (Pslk + C_PSLK_SHIFT), (Psdl + C_PSLK_SHIFT), C_PSLK_DELTA, T1, T2 );
    Pslk_dVbs = Pslk_dVbs * T1 + Psdl_dVbs * T2;
    Pslk_dVds = Pslk_dVds * T1 + Psdl_dVds * T2;
    Pslk_dVgs = Pslk_dVgs * T1 + Psdl_dVgs * T2;
    Pslk = T10 - C_PSLK_SHIFT;

    /* suppress Qy in accumulation region */
    /*
    Aclm_eff = 1.0 - Pds / (eps_qy + Pds) * (1.0 - Aclm) ;
    Aclm_eff_dVds = eps_qy * Pds_dVds / ((eps_qy + Pds)*(eps_qy + Pds)) ;
    Aclm_eff_dVgs = eps_qy * Pds_dVgs / ((eps_qy + Pds)*(eps_qy + Pds)) ;
    Aclm_eff_dVbs = eps_qy * Pds_dVbs / ((eps_qy + Pds)*(eps_qy + Pds)) ;
    */

    Aclm_eff = Aclm ; Aclm_eff_dVds = Aclm_eff_dVgs = Aclm_eff_dVbs = 0.0 ;

    T1 = Aclm_eff * ( Vds + Ps0 ) + ( 1.0e0 - Aclm_eff ) * Pslk ;
    T1_dVb = Aclm_eff * (  Ps0_dVbs ) + ( 1.0e0 - Aclm_eff ) * Pslk_dVbs
             + Aclm_eff_dVbs * ( Vds + Ps0 - Pslk ) ;
    T1_dVd = Aclm_eff * ( 1.0 + Ps0_dVds ) + ( 1.0e0 - Aclm_eff ) * Pslk_dVds
             + Aclm_eff_dVds * ( Vds + Ps0 - Pslk ) ;
    T1_dVg = Aclm_eff * (  Ps0_dVgs ) + ( 1.0e0 - Aclm_eff ) * Pslk_dVgs
             + Aclm_eff_dVgs * ( Vds + Ps0 - Pslk ) ;
    T10 = here->HSM2_wdpl ;
    T3 = T10 * 1.3 ;
    T2 = C_ESI * here->HSM2_weff_nf * T3 ;
    T7 = 1.0e-9 ; /* 1nm */
    T0 = Fn_Max( model->HSM2_xqy , T7 ) ;
    T4 = T2 / T0 ;
    Qy =       - ( Ps0      + Vds   - T1     ) * T4 ;
    Qy_dVds =  - ( Ps0_dVds + 1.0e0 - T1_dVd ) * T4 ;
    Qy_dVgs =  - ( Ps0_dVgs         - T1_dVg ) * T4 ;
    Qy_dVbs =  - ( Ps0_dVbs         - T1_dVb ) * T4 ;


  }

  if ( model->HSM2_xqy1 != 0.0 ){
    Qy += here->HSM2_cqyb0 * Vbs ;
    Qy_dVbs += here->HSM2_cqyb0 ;
  }

  Qy      = Qy * FMDVDS ;
  Qy_dVbs = Qy_dVbs * FMDVDS + Qy * FMDVDS_dVbs ;
  Qy_dVds = Qy_dVds * FMDVDS + Qy * FMDVDS_dVds ;
  Qy_dVgs = Qy_dVgs * FMDVDS + Qy * FMDVDS_dVgs ;

  /*---------------------------------------------------* 
   * Fringing capacitance.
   *-----------------*/ 
  Cf = here->HSM2_cfrng ;
  Qfd = Cf * ( Vgs - Vds ) ;
  Qfs = Cf * Vgs ;

  /*-----------------------------------------------------------* 
   * End of PART-3. (label) 
   *-----------------*/ 

/* end_of_part_3:*/

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   * PART-4: Substrate-source/drain junction diode.
   *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
    
  /*-----------------------------------------------------------*
   * Cbsj, Cbdj: node-base S/D biases.
   *-----------------*/

  T10 = model->HSM2_cvb * here->HSM2_jd_nvtm_inv ;
  T11 = model->HSM2_cvbk * here->HSM2_jd_nvtm_inv ;

  T9 = model->HSM2_cisb * here->HSM2_exptemp ;
  T0 = here->HSM2_isbd2 * T9 ;

  T2 = exp (- vbd_jct * T10 );
  T2_dVb = - T2 * T10 ;

  T3 = exp (- vbd_jct * T11 );
  T3_dVb = - T3 * T11 ;


  /* ibd */
  if ( vbd_jct < here->HSM2_vbdt ) {
    TX = vbd_jct * here->HSM2_jd_nvtm_inv ;
    if ( TX < - EXP_THR ) {
      T1 = 0.0 ;
      T1_dVb = 0.0 ;
    } else {
      T1 = exp ( TX ) ;
      T1_dVb = T1 * here->HSM2_jd_nvtm_inv ;
    }

    Ibd = here->HSM2_isbd * (T1 - 1.0) 
      + T0 * (T2 - 1.0) 
      + pParam->HSM2_cisbk * (T3 - 1.0);   
    Gbd = here->HSM2_isbd * T1_dVb 
      + T0 * T2_dVb 
      + pParam->HSM2_cisbk * T3_dVb ;

  } else {
    T1 = here->HSM2_jd_expcd ;

    T4 = here->HSM2_isbd * here->HSM2_jd_nvtm_inv  * T1 ;

    Ibd = here->HSM2_isbd * (T1 - 1.0) 
      + T4 * (vbd_jct - here->HSM2_vbdt) 
      + T0 * (T2 - 1.0)
      + pParam->HSM2_cisbk * (T3 - 1.0) ;
    Gbd = T4 
      + T0 * T2_dVb 
      + pParam->HSM2_cisbk * T3_dVb ;
  }  
  T12 = model->HSM2_divx * here->HSM2_isbd2 ;
  Ibd += T12 * vbd_jct ;
  Gbd += T12 ;

  /* ibs */
  T0 = here->HSM2_isbs2 * T9 ;

  TX = - vbs_jct * T10 ;
  if ( TX < - EXP_THR ) {
    T2 = 0.0 ;
    T2_dVb = 0.0 ;
  } else {
    T2 = exp ( TX );
    T2_dVb = - T2 * T10 ;
  }

  TX = - vbs_jct * T11 ;
  if ( TX < - EXP_THR ) {
    T3 = 0.0 ;
    T3_dVb = 0.0 ;
  } else {
    T3 = exp ( TX );
    T3_dVb = - T3 * T11 ;
  }

  if ( vbs_jct < here->HSM2_vbst ) {
    TX = vbs_jct * here->HSM2_jd_nvtm_inv ;
    if ( TX < - EXP_THR ) {
      T1 = 0.0 ;
      T1_dVb = 0.0 ;
    } else {
      T1 = exp ( TX ) ;
      T1_dVb = T1 * here->HSM2_jd_nvtm_inv ;
    }
    Ibs = here->HSM2_isbs * (T1 - 1.0) 
      + T0 * (T2 - 1.0) 
      + pParam->HSM2_cisbk * (T3 - 1.0);
    Gbs = here->HSM2_isbs * T1_dVb 
      + T0 * T2_dVb
      + pParam->HSM2_cisbk * T3_dVb ;
  } else {
    T1 = here->HSM2_jd_expcs ;

    T4 = here->HSM2_isbs * here->HSM2_jd_nvtm_inv  * T1 ;

    Ibs = here->HSM2_isbs * (T1 - 1.0)
      + T4 * (vbs_jct - here->HSM2_vbst)
      + T0 * (T2 - 1.0) 
      + pParam->HSM2_cisbk * (T3 - 1.0) ;
    Gbs = T4 
      + T0 * T2_dVb 
      + pParam->HSM2_cisbk * T3_dVb ;
  }
  T12 = model->HSM2_divx * here->HSM2_isbs2 ;
  Ibs += T12 * vbs_jct ;
  Gbs += T12 ;

  /*---------------------------------------------------*
   * Add Gjmin.
   *-----------------*/
  Ibd += Gjmin * vbd_jct ;
  Ibs += Gjmin * vbs_jct ;
  
  Gbd += Gjmin ;
  Gbs += Gjmin ;

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
  /*  add new parameters 
      tcjbs	: temperature dependence of czbs
      tcjbd	: temperature dependence of czbd
      tcjbssw	: temperature dependence of czbssw
      tcjbdsw	: temperature dependence of czbdsw
      tcjbsswg	: temperature dependence of czbsswg
      tcjbdswg	: temperature dependence of czbdswg
  */
  tcjbd=model->HSM2_tcjbd;
  tcjbs=model->HSM2_tcjbs;
  tcjbdsw=model->HSM2_tcjbdsw;
  tcjbssw=model->HSM2_tcjbssw;
  tcjbdswg=model->HSM2_tcjbdswg;
  tcjbsswg=model->HSM2_tcjbsswg;

  czbs = model->HSM2_cj * here->HSM2_as ;
  czbs = czbs * ( 1.0 + tcjbs * ( TTEMP - model->HSM2_ktnom )) ;

  czbd = model->HSM2_cj * here->HSM2_ad ;
  czbd = czbd * ( 1.0 + tcjbd * ( TTEMP - model->HSM2_ktnom )) ;

  /* Source Bulk Junction */
  if (here->HSM2_ps > here->HSM2_weff_nf) {
    czbssw = model->HSM2_cjsw * ( here->HSM2_ps - here->HSM2_weff_nf ) ;
    czbssw = czbssw * ( 1.0 + tcjbssw * ( TTEMP - model->HSM2_ktnom )) ;

    czbsswg = model->HSM2_cjswg * here->HSM2_weff_nf ;
    czbsswg = czbsswg * ( 1.0 + tcjbsswg * ( TTEMP - model->HSM2_ktnom )) ;

    if (vbs_jct == 0.0) {  
      Qbs = 0.0 ;
      Capbs = czbs + czbssw + czbsswg ;
    } else if (vbs_jct < 0.0) { 
      if (czbs > 0.0) {
        arg = 1.0 - vbs_jct / model->HSM2_pb ;
        if (model->HSM2_mj == 0.5) 
          sarg = 1.0 / sqrt(arg) ;
        else 
          sarg = Fn_Pow( arg , -model->HSM2_mj ) ;
        Qbs = model->HSM2_pb * czbs * (1.0 - arg * sarg) / (1.0 - model->HSM2_mj) ;
        Capbs = czbs * sarg ;
      } else {
        Qbs = 0.0 ;
        Capbs = 0.0 ;
      }
      if (czbssw > 0.0) {
        arg = 1.0 - vbs_jct / model->HSM2_pbsw ;
        if (model->HSM2_mjsw == 0.5) 
          sarg = 1.0 / sqrt(arg) ;
        else 
          sarg = Fn_Pow( arg , -model->HSM2_mjsw ) ;
        Qbs += model->HSM2_pbsw * czbssw * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjsw) ;
        Capbs += czbssw * sarg ;
      }
      if (czbsswg > 0.0) {
        arg = 1.0 - vbs_jct / model->HSM2_pbswg ;
        if (model->HSM2_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mjswg ) ;
        Qbs += model->HSM2_pbswg * czbsswg * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjswg) ;
        Capbs += czbsswg * sarg ;
      }
    } else {
      T1 = czbs + czbssw + czbsswg ;
      T2 = czbs * model->HSM2_mj / model->HSM2_pb 
        + czbssw * model->HSM2_mjsw / model->HSM2_pbsw 
        + czbsswg * model->HSM2_mjswg / model->HSM2_pbswg ;
      Qbs = vbs_jct * (T1 + vbs_jct * 0.5 * T2) ;
      Capbs = T1 + vbs_jct * T2 ;
    }
  } else {
    czbsswg = model->HSM2_cjswg * here->HSM2_ps ;
    czbsswg = czbsswg * ( 1.0 + tcjbsswg * ( TTEMP - model->HSM2_ktnom )) ;
    if (vbs_jct == 0.0) {
      Qbs = 0.0 ;
      Capbs = czbs + czbsswg ;
    } else if (vbs_jct < 0.0) {
      if (czbs > 0.0) {
        arg = 1.0 - vbs_jct / model->HSM2_pb ;
        if (model->HSM2_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mj ) ;
          Qbs = model->HSM2_pb * czbs * (1.0 - arg * sarg) / (1.0 - model->HSM2_mj) ;
          Capbs = czbs * sarg ;
      } else {
        Qbs = 0.0 ;
        Capbs = 0.0 ;
      }
      if (czbsswg > 0.0) {
        arg = 1.0 - vbs_jct / model->HSM2_pbswg ;
        if (model->HSM2_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mjswg ) ;
        Qbs += model->HSM2_pbswg * czbsswg * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjswg) ;
        Capbs += czbsswg * sarg ;
      }
    } else {
      T1 = czbs + czbsswg ;
      T2 = czbs * model->HSM2_mj / model->HSM2_pb
        + czbsswg * model->HSM2_mjswg / model->HSM2_pbswg ;
      Qbs = vbs_jct * (T1 + vbs_jct * 0.5 * T2) ;
      Capbs = T1 + vbs_jct * T2 ;
    }
  }    
    
  /* Drain Bulk Junction */
  if (here->HSM2_pd > here->HSM2_weff_nf) {
    czbdsw = model->HSM2_cjsw * ( here->HSM2_pd - here->HSM2_weff_nf ) ;
    czbdsw = czbdsw * ( 1.0 + tcjbdsw * ( TTEMP - model->HSM2_ktnom )) ; 
    czbdswg = model->HSM2_cjswg * here->HSM2_weff_nf ;
    czbdswg = czbdswg * ( 1.0 + tcjbdswg * ( TTEMP - model->HSM2_ktnom )) ;
    if (vbd_jct == 0.0) {
      Qbd = 0.0 ;
      Capbd = czbd + czbdsw + czbdswg ;
    } else if (vbd_jct < 0.0) {
      if (czbd > 0.0) {
        arg = 1.0 - vbd_jct / model->HSM2_pb ;
        if (model->HSM2_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mj ) ;
        Qbd = model->HSM2_pb * czbd * (1.0 - arg * sarg) / (1.0 - model->HSM2_mj) ;
        Capbd = czbd * sarg ;
      } else {
        Qbd = 0.0 ;
        Capbd = 0.0 ;
      }
      if (czbdsw > 0.0) {
        arg = 1.0 - vbd_jct / model->HSM2_pbsw ;
        if (model->HSM2_mjsw == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mjsw ) ;
        Qbd += model->HSM2_pbsw * czbdsw * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjsw) ;
        Capbd += czbdsw * sarg ;
      }
      if (czbdswg > 0.0) {
        arg = 1.0 - vbd_jct / model->HSM2_pbswg ;
        if (model->HSM2_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mjswg ) ;
        Qbd += model->HSM2_pbswg * czbdswg * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjswg) ;
        Capbd += czbdswg * sarg ;
      }
    } else {
      T1 = czbd + czbdsw + czbdswg ;
      T2 = czbd * model->HSM2_mj / model->HSM2_pb
        + czbdsw * model->HSM2_mjsw / model->HSM2_pbsw 
        + czbdswg * model->HSM2_mjswg / model->HSM2_pbswg ;
      Qbd = vbd_jct * (T1 + vbd_jct * 0.5 * T2) ;
      Capbd = T1 + vbd_jct * T2 ;
    }
    
  } else {
    czbdswg = model->HSM2_cjswg * here->HSM2_pd ;
    czbdswg = czbdswg * ( 1.0 + tcjbdswg * ( TTEMP - model->HSM2_ktnom )) ;
    if (vbd_jct == 0.0) {   
      Qbd = 0.0 ;
      Capbd = czbd + czbdswg ;
    } else if (vbd_jct < 0.0) {
      if (czbd > 0.0) {
        arg = 1.0 - vbd_jct / model->HSM2_pb ;
        if (model->HSM2_mj == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mj ) ;
        Qbd = model->HSM2_pb * czbd * (1.0 - arg * sarg) / (1.0 - model->HSM2_mj) ;
        Capbd = czbd * sarg ;
      } else {
        Qbd = 0.0 ;
        Capbd = 0.0 ;
      }
      if (czbdswg > 0.0) {
        arg = 1.0 - vbd_jct / model->HSM2_pbswg ;
        if (model->HSM2_mjswg == 0.5)
          sarg = 1.0 / sqrt(arg) ;
        else
          sarg = Fn_Pow( arg , -model->HSM2_mjswg ) ;
        Qbd += model->HSM2_pbswg * czbdswg * (1.0 - arg * sarg) / (1.0 - model->HSM2_mjswg) ;
        Capbd += czbdswg * sarg ;
      }
    } else {
      T1 = czbd + czbdswg ;
      T2 = czbd * model->HSM2_mj / model->HSM2_pb
        + czbdswg * model->HSM2_mjswg / model->HSM2_pbswg ;
      Qbd = vbd_jct * (T1 + vbd_jct * 0.5 * T2) ;
      Capbd = T1 + vbd_jct * T2 ;
    }
  }
  
  /*-----------------------------------------------------------* 
   * End of PART-4. (label) 
   *-----------------*/ 

/* end_of_part_4:*/

  

  /*-----------------------------------------------------------* 
   * PART-5: NQS. (label) 
   *-----------------*/
  if (flg_nqs) {

    if(ckt->CKTmode & MODETRAN){
      if( ckt->CKTmode & MODEINITTRAN ){
        Qi_nqs = Qi ;
        Qi_dVgs_nqs = Qi_dVgs ;
        Qi_dVds_nqs = Qi_dVds ;
        Qi_dVbs_nqs = Qi_dVbs ;
          
        Qb_nqs = Qb ;
        Qb_dVgs_nqs = Qb_dVgs ;
        Qb_dVds_nqs = Qb_dVds ;
        Qb_dVbs_nqs = Qb_dVbs ;
      } else {
        /* tau for inversion charge */
        if (flg_noqi == 0) {
          T12 = model->HSM2_dly1;
          T10 = model->HSM2_dly2;
            
          T3 = Lch ;
          T1 = T10 * T12 * T3 * T3 ;
          T2 = Mu * VgVt * T12 + T10 * T3 * T3 + small ;
          tau = T1 / T2 ;

          T1_dVg = T10 * T12 * 2.0 * T3 * Lch_dVgs ;
          T1_dVd = T10 * T12 * 2.0 * T3 * Lch_dVds ;
          T1_dVb = T10 * T12 * 2.0 * T3 * Lch_dVbs ;

          T2_dVg = T12 * Mu_dVgs * VgVt 
            + T12 * Mu * VgVt_dVgs + T10 * 2.0 * T3 * Lch_dVgs ;
          T2_dVd = T12 * Mu_dVds * VgVt 
            + T12 * Mu * VgVt_dVds + T10 * 2.0 * T3 * Lch_dVds ;
          T2_dVb = T12 * Mu_dVbs * VgVt 
            + T12 * Mu * VgVt_dVbs + T10 * 2.0 * T3 * Lch_dVbs ;
      
          T4 = 1.0 / ( T2 * T2 ) ;
          tau_dVgs = ( T2 * T1_dVg - T1 * T2_dVg ) * T4 ;
          tau_dVds = ( T2 * T1_dVd - T1 * T2_dVd ) * T4 ;
          tau_dVbs = ( T2 * T1_dVb - T1 * T2_dVb ) * T4 ;
        } else {
          tau = model->HSM2_dly1 + small ;
          tau_dVgs = tau_dVds = tau_dVbs = 0.0 ;
        }
          
        T1 = ckt->CKTdelta ;

        /* Calculation of Qi */
        Qi_prev = *(ckt->CKTstate1 + here->HSM2qi_nqs) ;
        T2 = T1 + tau ;
        T0 = Qi - Qi_prev ;
        Qi_nqs = Qi_prev + T1 / T2 * T0;
        T3 = T1 / T2 ;
        T4 = T0 / T2 ;
        Qi_dVgs_nqs = T3 * (Qi_dVgs - T4 * tau_dVgs);
        Qi_dVds_nqs = T3 * (Qi_dVds - T4 * tau_dVds);
        Qi_dVbs_nqs = T3 * (Qi_dVbs - T4 * tau_dVbs);

        /* tau for bulk charge */
        T2 = modelMKS->HSM2_dly3 ;
        taub = T2 * Cox ;
        taub_dVgs = T2 * Cox_dVg ;
        taub_dVds = T2 * Cox_dVd ;
        taub_dVbs = T2 * Cox_dVb ;
        /* Calculation of Qb */
        Qb_prev = *(ckt->CKTstate1 + here->HSM2qb_nqs) ;
        T2 = T1 + taub ;
        T0 = Qb - Qb_prev ;
        Qb_nqs = Qb_prev + T1 / T2 * T0 ;
        T3 = T1 / T2 ;
        T4 = T0 / T2 ;
        Qb_dVgs_nqs = T3 * (Qb_dVgs - T4 * taub_dVgs) ;
        Qb_dVds_nqs = T3 * (Qb_dVds - T4 * taub_dVds) ;
        Qb_dVbs_nqs = T3 * (Qb_dVbs - T4 * taub_dVbs) ;
      }
    } else { /* !(CKT_mode & MODETRAN) */
      Qi_nqs = Qi ;
      Qi_dVgs_nqs = Qi_dVgs ;
      Qi_dVds_nqs = Qi_dVds ;
      Qi_dVbs_nqs = Qi_dVbs ;
      
      Qb_nqs = Qb ;
      Qb_dVgs_nqs = Qb_dVgs ;
      Qb_dVds_nqs = Qb_dVds ;
      Qb_dVbs_nqs = Qb_dVbs ;
    }
  }
    
  if ( flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG)) ) { /* ACNQS */

    if (flg_noqi == 0) {
      T10 = model->HSM2_dly1 ;
      T11 = model->HSM2_dly2 ;
      T12 = Lch ;

      T1 = T10 * T11 * T12 * T12 ;
      T2 = Mu * VgVt * T10 + T11 * T12 * T12 + small ;
      tau = T1 / T2 ;

      T1_dVg = T10 * T11 * 2.0 * T12 * Lch_dVgs ;
      T1_dVd = T10 * T11 * 2.0 * T12 * Lch_dVds ;
      T1_dVb = T10 * T11 * 2.0 * T12 * Lch_dVbs ;
      T2_dVg = T10 * Mu_dVgs * VgVt + T10 * Mu * VgVt_dVgs 
        + T11 * 2.0 * T12 * Lch_dVgs ;
      T2_dVd = T10 * Mu_dVds * VgVt + T10 * Mu * VgVt_dVds 
        + T11 * 2.0 * T12 * Lch_dVds ;
      T2_dVb = T10 * Mu_dVbs * VgVt + T10 * Mu * VgVt_dVbs 
        + T11 * 2.0 * T12 * Lch_dVbs ;
     
      T3 = 1.0 / T2 ;
      tau_dVgs = (T1_dVg - tau * T2_dVg) * T3 ;
      tau_dVds = (T1_dVd - tau * T2_dVd) * T3 ;
      tau_dVbs = (T1_dVb - tau * T2_dVb) * T3 ;
    } else {
      tau = model->HSM2_dly1 + small ;
      tau_dVgs = tau_dVds = tau_dVbs = 0.0 ;
    }

    T1 = modelMKS->HSM2_dly3 ;
    taub = T1 * Cox ;
    taub_dVgs = T1 * Cox_dVg ;
    taub_dVds = T1 * Cox_dVd ;
    taub_dVbs = T1 * Cox_dVb ;
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
  if ( model->HSM2_coflick != 0 && !flg_noqi ) {
    
    NFalp = pParam->HSM2_nfalp ;
    NFtrp = pParam->HSM2_nftrp ;
    Cit = modelMKS->HSM2_cit ;

    T1 = Qn0 / C_QE ;
    T2 = ( Cox + Qn0 / ( Ps0 - Vbs ) + Cit ) * beta_inv / C_QE ;
    T3 = -2.0E0 * Qi / C_QE / Lch / here->HSM2_weff_nf - T1 ;
    if ( T3 != T1 ) {
      T4 = 1.0E0 / ( T1 + T2 ) / ( T3 + T2 ) + 2.0E0 * NFalp * Ey * Mu / ( T3 - T1 )
        * log( ( T3 + T2 ) / ( T1 + T2 ) ) + NFalp * Ey * Mu * NFalp * Ey * Mu ;
    }  else {
      T4 = 1.0 / ( T1 + T2 ) / ( T3 + T2 ) + 2.0 * NFalp * Ey * Mu / ( T1 + T2 )
        + NFalp * Ey * Mu * NFalp * Ey * Mu;
    }
    Nflic = Ids * Ids * NFtrp / ( Lch * beta * here->HSM2_weff_nf ) * T4 ;
  } else {
    Nflic = 0.0 ;
  }
  
  /*-----------------------------------------------------------*
   * thermal noise.
   *-----------------*/
  if ( model->HSM2_cothrml != 0 && !flg_noqi ) {

    Eyd = ( Psdl - Ps0 ) / Lch ;
    T12 = Muun * Eyd / C_vmax ;
    /* note: model->HSM2_bb = 2 (electron) ;1 (hole) */
    if ( 1.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 1.0e0 + epsm10 ) {
      T7  = 1.0e0 ;
    } else if ( 2.0e0 - epsm10 <= model->HSM2_bb && model->HSM2_bb <= 2.0e0 + epsm10 ) {
      T7  = T12 ; 
    } else {
      T7  = Fn_Pow( Eyd, model->HSM2_bb - 1.0e0 ) ;
    }
    T8 = T12 * T7 ;
    T9 = 1.0e0 + T8 ;
    T10 = Fn_Pow( T9, ( - 1.0e0 / model->HSM2_bb - 1.0e0 ) ) ;
    T11 = T9 * T10 ;
    Mud_hoso = Muun * T11 ;
    Mu_Ave = ( Mu + Mud_hoso ) / 2.0 ;
    
    /* Sid_h = GAMMA * 4.0 * C_KB * model->HSM2_temp * gds0_h2; */
    T0 = Alpha * Alpha ;
    Nthrml  = here->HSM2_weff_nf * Cox * VgVt * Mu
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
  if ( model->HSM2_coign != 0 && model->HSM2_cothrml != 0 && flg_ign == 1 && !flg_noqi ) {
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
    gds0_ign = here->HSM2_weff_nf / Lch * Mu * Cox ;
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
  /* note: T1 = 1 + Delta */
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
   * - All derivatives that influence outputs must be modified here.
   *-----------------*/ 

  /*---------------------------------------------------* 
   * Case ignoring Rs/Rd.
   *-----------------*/ 
  if ( flg_rsrd != 1 || Ids < 0.0 ) {
    Ids_dVbse = Ids_dVbs ; Ids_dVdse = Ids_dVds ; Ids_dVgse = Ids_dVgs ;
    Qb_dVbse = Qb_dVbs ; Qb_dVdse = Qb_dVds ; Qb_dVgse = Qb_dVgs ;
    Qi_dVbse = Qi_dVbs ; Qi_dVdse = Qi_dVds ; Qi_dVgse = Qi_dVgs ;
    Qd_dVbse = Qd_dVbs ; Qd_dVdse = Qd_dVds ; Qd_dVgse = Qd_dVgs ;
    Isub_dVbse = Isub_dVbs ; Isub_dVdse = Isub_dVds ; Isub_dVgse = Isub_dVgs ;
    IdsIBPC_dVbse = IdsIBPC_dVbs ; IdsIBPC_dVdse = IdsIBPC_dVds ; IdsIBPC_dVgse = IdsIBPC_dVgs ;
    Igate_dVbse = Igate_dVbs ; Igate_dVdse = Igate_dVds ; Igate_dVgse = Igate_dVgs ;
    Igd_dVbse = Igd_dVbs ; Igd_dVdse = Igd_dVds ; Igd_dVgse = Igd_dVgs ;
    Igs_dVbse = Igs_dVbs ; Igs_dVdse = Igs_dVds ; Igs_dVgse = Igs_dVgs ;
    Igb_dVbse = Igb_dVbs ; Igb_dVdse = Igb_dVds ; Igb_dVgse = Igb_dVgs ;    
    Igidl_dVbse = Igidl_dVbs ; Igidl_dVdse = Igidl_dVds ; Igidl_dVgse = Igidl_dVgs ;
    Igisl_dVbse = Igisl_dVbs ; Igisl_dVdse = Igisl_dVds ; Igisl_dVgse = Igisl_dVgs ;
    Qgos_dVbse = Qgos_dVbs ; Qgos_dVdse = Qgos_dVds ; Qgos_dVgse = Qgos_dVgs ;
    Qgod_dVbse = Qgod_dVbs ; Qgod_dVdse = Qgod_dVds ; Qgod_dVgse = Qgod_dVgs ;
    Qgbo_dVbse = Qgbo_dVbs ; Qgbo_dVdse = Qgbo_dVds ; Qgbo_dVgse = Qgbo_dVgs ;
    Qovd_dVbse = Qovd_dVbs ; Qovd_dVdse = Qovd_dVds ; Qovd_dVgse = Qovd_dVgs ;
    QidLD_dVbse = QidLD_dVbs ; QidLD_dVdse = QidLD_dVds ; QidLD_dVgse = QidLD_dVgs ;
    QbdLD_dVbse = QbdLD_dVbs ; QbdLD_dVdse = QbdLD_dVds ; QbdLD_dVgse = QbdLD_dVgs ;
    Qovs_dVbse = Qovs_dVbs ; Qovs_dVdse = Qovs_dVds ; Qovs_dVgse = Qovs_dVgs ;
    QisLD_dVbse = QisLD_dVbs ; QisLD_dVdse = QisLD_dVds ; QisLD_dVgse = QisLD_dVgs ;
    QbsLD_dVbse = QbsLD_dVbs ; QbsLD_dVdse = QbsLD_dVds ; QbsLD_dVgse = QbsLD_dVgs ;
    Qy_dVbse = Qy_dVbs ; Qy_dVdse = Qy_dVds ; Qy_dVgse = Qy_dVgs ;
    Qdrat_dVbse = Qdrat_dVbs; Qdrat_dVdse = Qdrat_dVds; Qdrat_dVgse = Qdrat_dVgs;
    GLPART1_dVbse = GLPART1_dVbs ;GLPART1_dVdse = GLPART1_dVds ;GLPART1_dVgse = GLPART1_dVgs ;

    if (flg_nqs) { /* NQS */
      Qi_dVgse_nqs = Qi_dVgs_nqs; Qi_dVdse_nqs = Qi_dVds_nqs; Qi_dVbse_nqs = Qi_dVbs_nqs;
      Qb_dVgse_nqs = Qb_dVgs_nqs; Qb_dVdse_nqs = Qb_dVds_nqs; Qb_dVbse_nqs = Qb_dVbs_nqs;
    }
    if (flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG))) { /* ACNQS */
      tau_dVgse = tau_dVgs ; tau_dVdse = tau_dVds ; tau_dVbse = tau_dVbs ;
      taub_dVgse = taub_dVgs ; taub_dVdse = taub_dVds ; taub_dVbse = taub_dVbs ;
    }

  } else {
    /*---------------------------------------------------* 
     * Case Rs>0 or Rd>0
     *-----------------*/ 

    /*-------------------------------------------* 
     * Conductances w.r.t. confined biases.
     *-----------------*/ 
    Ids_dVbse = Ids_dVbs * DJI ;
    Ids_dVdse = Ids_dVds * DJI ;
    Ids_dVgse = Ids_dVgs * DJI ;
      
    /*-------------------------------------------* 
     * Derivatives of internal biases  w.r.t. external biases.
     *-----------------*/ 
    Vbs_dVbse = ( 1.0 - Rs * Ids_dVbse ) ;
    Vbs_dVdse = - Rs * Ids_dVdse ;
    Vbs_dVgse = - Rs * Ids_dVgse ;
    
    Vds_dVbse = - ( Rs + Rd ) * Ids_dVbse ;
    Vds_dVdse = ( 1.0 - ( Rs + Rd ) * Ids_dVdse ) ;
    Vds_dVgse = - ( Rs + Rd ) * Ids_dVgse ;
      
    Vgs_dVbse = - Rs * Ids_dVbse ;
    Vgs_dVdse = - Rs * Ids_dVdse ;
    Vgs_dVgse = ( 1.0 - Rs * Ids_dVgse ) ;

    /*-------------------------------------------* 
     * Derivatives of charges.
     *-----------------*/ 
    Qb_dVbse = Qb_dVbs * Vbs_dVbse + Qb_dVds * Vds_dVbse + Qb_dVgs * Vgs_dVbse ;
    Qb_dVdse = Qb_dVbs * Vbs_dVdse + Qb_dVds * Vds_dVdse + Qb_dVgs * Vgs_dVdse ;
    Qb_dVgse = Qb_dVbs * Vbs_dVgse + Qb_dVds * Vds_dVgse + Qb_dVgs * Vgs_dVgse ;
    Qi_dVbse = Qi_dVbs * Vbs_dVbse + Qi_dVds * Vds_dVbse + Qi_dVgs * Vgs_dVbse ;
    Qi_dVdse = Qi_dVbs * Vbs_dVdse + Qi_dVds * Vds_dVdse + Qi_dVgs * Vgs_dVdse ;
    Qi_dVgse = Qi_dVbs * Vbs_dVgse + Qi_dVds * Vds_dVgse + Qi_dVgs * Vgs_dVgse ;
    Qd_dVbse = Qd_dVbs * Vbs_dVbse + Qd_dVds * Vds_dVbse + Qd_dVgs * Vgs_dVbse ;
    Qd_dVdse = Qd_dVbs * Vbs_dVdse + Qd_dVds * Vds_dVdse + Qd_dVgs * Vgs_dVdse ;
    Qd_dVgse = Qd_dVbs * Vbs_dVgse + Qd_dVds * Vds_dVgse + Qd_dVgs * Vgs_dVgse ;

    /*-------------------------------------------* 
     * Substrate/gate/leak conductances.
     *-----------------*/ 
    Isub_dVbse = Isub_dVbs * Vbs_dVbse + Isub_dVds * Vds_dVbse + Isub_dVgs * Vgs_dVbse ;
    Isub_dVdse = Isub_dVbs * Vbs_dVdse + Isub_dVds * Vds_dVdse + Isub_dVgs * Vgs_dVdse ;
    Isub_dVgse = Isub_dVbs * Vbs_dVgse + Isub_dVds * Vds_dVgse + Isub_dVgs * Vgs_dVgse ;
    IdsIBPC_dVbse = IdsIBPC_dVbs * Vbs_dVbse + IdsIBPC_dVds * Vds_dVbse + IdsIBPC_dVgs * Vgs_dVbse ;
    IdsIBPC_dVdse = IdsIBPC_dVbs * Vbs_dVdse + IdsIBPC_dVds * Vds_dVdse + IdsIBPC_dVgs * Vgs_dVdse ;
    IdsIBPC_dVgse = IdsIBPC_dVbs * Vbs_dVgse + IdsIBPC_dVds * Vds_dVgse + IdsIBPC_dVgs * Vgs_dVgse ;
    Igate_dVbse = Igate_dVbs * Vbs_dVbse + Igate_dVds * Vds_dVbse + Igate_dVgs * Vgs_dVbse ;
    Igate_dVdse = Igate_dVbs * Vbs_dVdse + Igate_dVds * Vds_dVdse + Igate_dVgs * Vgs_dVdse ;
    Igate_dVgse = Igate_dVbs * Vbs_dVgse + Igate_dVds * Vds_dVgse + Igate_dVgs * Vgs_dVgse ;
    Igb_dVbse = Igb_dVbs * Vbs_dVbse + Igb_dVds * Vds_dVbse + Igb_dVgs * Vgs_dVbse ;
    Igb_dVdse = Igb_dVbs * Vbs_dVdse + Igb_dVds * Vds_dVdse + Igb_dVgs * Vgs_dVdse ;
    Igb_dVgse = Igb_dVbs * Vbs_dVgse + Igb_dVds * Vds_dVgse + Igb_dVgs * Vgs_dVgse ;
    Igd_dVbse = Igd_dVbs * Vbs_dVbse + Igd_dVds * Vds_dVbse + Igd_dVgs * Vgs_dVbse ;
    Igd_dVdse = Igd_dVbs * Vbs_dVdse + Igd_dVds * Vds_dVdse + Igd_dVgs * Vgs_dVdse ;
    Igd_dVgse = Igd_dVbs * Vbs_dVgse + Igd_dVds * Vds_dVgse + Igd_dVgs * Vgs_dVgse ;
    Igs_dVbse = Igs_dVbs * Vbs_dVbse + Igs_dVds * Vds_dVbse + Igs_dVgs * Vgs_dVbse ;
    Igs_dVdse = Igs_dVbs * Vbs_dVdse + Igs_dVds * Vds_dVdse + Igs_dVgs * Vgs_dVdse ;
    Igs_dVgse = Igs_dVbs * Vbs_dVgse + Igs_dVds * Vds_dVgse + Igs_dVgs * Vgs_dVgse ;
    Igidl_dVbse = Igidl_dVbs * Vbs_dVbse + Igidl_dVds * Vds_dVbse + Igidl_dVgs * Vgs_dVbse ;
    Igidl_dVdse = Igidl_dVbs * Vbs_dVdse + Igidl_dVds * Vds_dVdse + Igidl_dVgs * Vgs_dVdse ;
    Igidl_dVgse = Igidl_dVbs * Vbs_dVgse + Igidl_dVds * Vds_dVgse + Igidl_dVgs * Vgs_dVgse ;
    Igisl_dVbse = Igisl_dVbs * Vbs_dVbse + Igisl_dVds * Vds_dVbse + Igisl_dVgs * Vgs_dVbse ;
    Igisl_dVdse = Igisl_dVbs * Vbs_dVdse + Igisl_dVds * Vds_dVdse + Igisl_dVgs * Vgs_dVdse ;
    Igisl_dVgse = Igisl_dVbs * Vbs_dVgse + Igisl_dVds * Vds_dVgse + Igisl_dVgs * Vgs_dVgse ;
      
    GLPART1_dVbse = GLPART1_dVbs * Vbs_dVbse + GLPART1_dVds * Vds_dVbse + GLPART1_dVgs * Vgs_dVbse ;
    GLPART1_dVdse = GLPART1_dVbs * Vbs_dVdse + GLPART1_dVds * Vds_dVdse + GLPART1_dVgs * Vgs_dVdse ;
    GLPART1_dVgse = GLPART1_dVbs * Vbs_dVgse + GLPART1_dVds * Vds_dVgse + GLPART1_dVgs * Vgs_dVgse ;

    /*---------------------------------------------------* 
     * Derivatives of overlap charges.
     *-----------------*/ 
    Qgos_dVbse = Qgos_dVbs * Vbs_dVbse + Qgos_dVds * Vds_dVbse + Qgos_dVgs * Vgs_dVbse ;
    Qgos_dVdse = Qgos_dVbs * Vbs_dVdse + Qgos_dVds * Vds_dVdse + Qgos_dVgs * Vgs_dVdse ;
    Qgos_dVgse = Qgos_dVbs * Vbs_dVgse + Qgos_dVds * Vds_dVgse + Qgos_dVgs * Vgs_dVgse ;
    Qgod_dVbse = Qgod_dVbs * Vbs_dVbse + Qgod_dVds * Vds_dVbse + Qgod_dVgs * Vgs_dVbse ;
    Qgod_dVdse = Qgod_dVbs * Vbs_dVdse + Qgod_dVds * Vds_dVdse + Qgod_dVgs * Vgs_dVdse ;
    Qgod_dVgse = Qgod_dVbs * Vbs_dVgse + Qgod_dVds * Vds_dVgse + Qgod_dVgs * Vgs_dVgse ;
    Qgbo_dVbse = Qgbo_dVbs * Vbs_dVbse + Qgbo_dVds * Vds_dVbse + Qgbo_dVgs * Vgs_dVbse ;
    Qgbo_dVdse = Qgbo_dVbs * Vbs_dVdse + Qgbo_dVds * Vds_dVdse + Qgbo_dVgs * Vgs_dVdse ;
    Qgbo_dVgse = Qgbo_dVbs * Vbs_dVgse + Qgbo_dVds * Vds_dVgse + Qgbo_dVgs * Vgs_dVgse ;
    Qovd_dVbse = Qovd_dVbs * Vbs_dVbse + Qovd_dVds * Vds_dVbse + Qovd_dVgs * Vgs_dVbse ;
    Qovd_dVdse = Qovd_dVbs * Vbs_dVdse + Qovd_dVds * Vds_dVdse + Qovd_dVgs * Vgs_dVdse ;
    Qovd_dVgse = Qovd_dVbs * Vbs_dVgse + Qovd_dVds * Vds_dVgse + Qovd_dVgs * Vgs_dVgse ;
    QidLD_dVbse = QidLD_dVbs * Vbs_dVbse + QidLD_dVds * Vds_dVbse + QidLD_dVgs * Vgs_dVbse ;
    QidLD_dVdse = QidLD_dVbs * Vbs_dVdse + QidLD_dVds * Vds_dVdse + QidLD_dVgs * Vgs_dVdse ;
    QidLD_dVgse = QidLD_dVbs * Vbs_dVgse + QidLD_dVds * Vds_dVgse + QidLD_dVgs * Vgs_dVgse ;
    QbdLD_dVbse = QbdLD_dVbs * Vbs_dVbse + QbdLD_dVds * Vds_dVbse + QbdLD_dVgs * Vgs_dVbse ;
    QbdLD_dVdse = QbdLD_dVbs * Vbs_dVdse + QbdLD_dVds * Vds_dVdse + QbdLD_dVgs * Vgs_dVdse ;
    QbdLD_dVgse = QbdLD_dVbs * Vbs_dVgse + QbdLD_dVds * Vds_dVgse + QbdLD_dVgs * Vgs_dVgse ;
    Qovs_dVbse = Qovs_dVbs * Vbs_dVbse + Qovs_dVds * Vds_dVbse + Qovs_dVgs * Vgs_dVbse ;
    Qovs_dVdse = Qovs_dVbs * Vbs_dVdse + Qovs_dVds * Vds_dVdse + Qovs_dVgs * Vgs_dVdse ;
    Qovs_dVgse = Qovs_dVbs * Vbs_dVgse + Qovs_dVds * Vds_dVgse + Qovs_dVgs * Vgs_dVgse ;
    QisLD_dVbse = QisLD_dVbs * Vbs_dVbse + QisLD_dVds * Vds_dVbse + QisLD_dVgs * Vgs_dVbse ;
    QisLD_dVdse = QisLD_dVbs * Vbs_dVdse + QisLD_dVds * Vds_dVdse + QisLD_dVgs * Vgs_dVdse ;
    QisLD_dVgse = QisLD_dVbs * Vbs_dVgse + QisLD_dVds * Vds_dVgse + QisLD_dVgs * Vgs_dVgse ;
    QbsLD_dVbse = QbsLD_dVbs * Vbs_dVbse + QbsLD_dVds * Vds_dVbse + QbsLD_dVgs * Vgs_dVbse ;
    QbsLD_dVdse = QbsLD_dVbs * Vbs_dVdse + QbsLD_dVds * Vds_dVdse + QbsLD_dVgs * Vgs_dVdse ;
    QbsLD_dVgse = QbsLD_dVbs * Vbs_dVgse + QbsLD_dVds * Vds_dVgse + QbsLD_dVgs * Vgs_dVgse ;
    Qy_dVbse = Qy_dVbs * Vbs_dVbse + Qy_dVds * Vds_dVbse + Qy_dVgs * Vgs_dVbse ;
    Qy_dVdse = Qy_dVbs * Vbs_dVdse + Qy_dVds * Vds_dVdse + Qy_dVgs * Vgs_dVdse ;
    Qy_dVgse = Qy_dVbs * Vbs_dVgse + Qy_dVds * Vds_dVgse + Qy_dVgs * Vgs_dVgse ;
    Qdrat_dVbse = Qdrat_dVbs * Vbs_dVbse + Qdrat_dVds * Vds_dVbse + Qdrat_dVgs * Vgs_dVbse ;
    Qdrat_dVdse = Qdrat_dVbs * Vbs_dVdse + Qdrat_dVds * Vds_dVdse + Qdrat_dVgs * Vgs_dVdse ;
    Qdrat_dVgse = Qdrat_dVbs * Vbs_dVgse + Qdrat_dVds * Vds_dVgse + Qdrat_dVgs * Vgs_dVgse ;

    if (flg_nqs) { /* NQS */
      Qi_dVgse_nqs = Qi_dVgs_nqs * Vgs_dVgse + Qi_dVds_nqs * Vds_dVgse + Qi_dVbs_nqs * Vbs_dVgse;
      Qi_dVdse_nqs = Qi_dVgs_nqs * Vgs_dVdse + Qi_dVds_nqs * Vds_dVdse + Qi_dVbs_nqs * Vbs_dVdse;
      Qi_dVbse_nqs = Qi_dVgs_nqs * Vgs_dVbse + Qi_dVds_nqs * Vds_dVbse + Qi_dVbs_nqs * Vbs_dVbse;
      Qb_dVgse_nqs = Qb_dVgs_nqs * Vgs_dVgse + Qb_dVds_nqs * Vds_dVgse + Qb_dVbs_nqs * Vbs_dVgse;
      Qb_dVdse_nqs = Qb_dVgs_nqs * Vgs_dVdse + Qb_dVds_nqs * Vds_dVdse + Qb_dVbs_nqs * Vbs_dVdse;
      Qb_dVbse_nqs = Qb_dVgs_nqs * Vgs_dVbse + Qb_dVds_nqs * Vds_dVbse + Qb_dVbs_nqs * Vbs_dVbse;
    }
    if (flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG))) { /* ACNQS */
      tau_dVgse = tau_dVgs * Vgs_dVgse + tau_dVds * Vds_dVgse + tau_dVbs * Vbs_dVgse;
      tau_dVdse = tau_dVgs * Vgs_dVdse + tau_dVds * Vds_dVdse + tau_dVbs * Vbs_dVdse;
      tau_dVbse = tau_dVgs * Vgs_dVbse + tau_dVds * Vds_dVbse + tau_dVbs * Vbs_dVbse;
      taub_dVgse = taub_dVgs * Vgs_dVgse + taub_dVds * Vds_dVgse + taub_dVbs * Vbs_dVgse;
      taub_dVdse = taub_dVgs * Vgs_dVdse + taub_dVds * Vds_dVdse + taub_dVbs * Vbs_dVdse;
      taub_dVbse = taub_dVgs * Vgs_dVbse + taub_dVds * Vds_dVbse + taub_dVbs * Vbs_dVbse;
    }
  } /* end of if ( flg_rsrd == 0 ) blocks */

  /*-------------------------------------------*
   * Add IdsIBPC to Ids.
   *-----------------*/ 
  Ids += IdsIBPC ;
  Ids_dVbse += IdsIBPC_dVbse ;
  Ids_dVdse += IdsIBPC_dVdse ;
  Ids_dVgse += IdsIBPC_dVgse ;

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
   *-----------------*/
  if ( flg_vbsc == 1 ) {
    Ids_dVbse   *= Vbsc_dVbse ;
    Qb_dVbse    *= Vbsc_dVbse ;
    Qi_dVbse    *= Vbsc_dVbse ;
    Qd_dVbse    *= Vbsc_dVbse ;
    Isub_dVbse  *= Vbsc_dVbse ;
    Igate_dVbse *= Vbsc_dVbse ;
    Igs_dVbse   *= Vbsc_dVbse ;
    Igd_dVbse   *= Vbsc_dVbse ;
    Igb_dVbse   *= Vbsc_dVbse ;
    Igidl_dVbse *= Vbsc_dVbse ;
    Igisl_dVbse *= Vbsc_dVbse ;
    Qgos_dVbse  *= Vbsc_dVbse ;
    Qgod_dVbse  *= Vbsc_dVbse ;
    Qgbo_dVbse  *= Vbsc_dVbse ;
    Qy_dVbse    *= Vbsc_dVbse ;
    if (flg_nqs) { 
      Qi_dVbse_nqs *= Vbsc_dVbse ;
      Qb_dVbse_nqs *= Vbsc_dVbse ;
    }
    if (flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG))) { /* ACNQS */
      tau_dVbse *= Vbsc_dVbse ;
      taub_dVbse *= Vbsc_dVbse ;
    }    
  } else if ( flg_vbsc == -1 ) {
    T1 = Vbse - Vbsc ;

    TX = Ids + T1 * Ids_dVbse ;
    if ( TX * Ids >= 0.0 ) {
      Ids = TX ;
    } else {
      Ids_dVbse = 0.0 ;
      Ids_dVdse = 0.0 ;
      Ids_dVgse = 0.0 ;
      Ids = 0.0 ;
    }

    TX = Qb  + T1 * Qb_dVbse ;
    /*note: The sign of Qb can be changed.*/
    Qb = TX ;

    TX = Qd  + T1 * Qd_dVbse ;
    if ( TX * Qd >= 0.0 ) {
      Qd = TX ;
    } else {
      Qd_dVbse = 0.0 ;
      Qd_dVdse = 0.0 ;
      Qd_dVgse = 0.0 ;
      Qd = 0.0 ;
    }

    TX = Qi  + T1 * Qi_dVbse ;
    if ( TX * Qi >= 0.0 ) {
      Qi = TX ;
    } else {
      Qi_dVbse = 0.0 ;
      Qi_dVdse = 0.0 ;
      Qi_dVgse = 0.0 ;
      Qi = 0.0 ;
    }

    TX = Isub + T1 * Isub_dVbse ;
    if ( TX * Isub >= 0.0 ) {
      Isub = TX ;
    } else {
      Isub_dVbse = 0.0 ;
      Isub_dVdse = 0.0 ;
      Isub_dVgse = 0.0 ;
      Isub = 0.0 ;
    }

    TX = Igate + T1 * Igate_dVbse ;
    if ( TX * Igate >= 0.0 ) {
      Igate = TX ;
    } else {
      Igate_dVbse = 0.0 ;
      Igate_dVdse = 0.0 ;
      Igate_dVgse = 0.0 ;
      Igate = 0.0 ;
    }

    TX = Igs + T1 * Igs_dVbse ;
    if ( TX * Igs >= 0.0 ) {
      Igs = TX ;
    } else {
      Igs_dVbse = 0.0 ;
      Igs_dVdse = 0.0 ;
      Igs_dVgse = 0.0 ;
      Igs = 0.0 ;
    }

    TX = Igd + T1 * Igd_dVbse ;
    if ( TX * Igd >= 0.0 ) {
      Igd = TX ;
    } else {
      Igd_dVbse = 0.0 ;
      Igd_dVdse = 0.0 ;
      Igd_dVgse = 0.0 ;
      Igd = 0.0 ;
    }

    TX = Igb + T1 * Igb_dVbse ;
    if ( TX * Igb >= 0.0 ) {
      Igb = TX ;
    } else {
      Igb_dVbse = 0.0 ;
      Igb_dVdse = 0.0 ;
      Igb_dVgse = 0.0 ;
      Igb = 0.0 ;
    }

    TX = Igidl + T1 * Igidl_dVbse ;
    if ( TX * Igidl >= 0.0 ) {
      Igidl = TX ;
    } else {
      Igidl_dVbse = 0.0 ;
      Igidl_dVdse = 0.0 ;
      Igidl_dVgse = 0.0 ;
      Igidl = 0.0 ;
    }

    TX = Igisl + T1 * Igisl_dVbse ;
    if ( TX * Igisl >= 0.0 ) {
      Igisl = TX ;
    } else {
      Igisl_dVbse = 0.0 ;
      Igisl_dVdse = 0.0 ;
      Igisl_dVgse = 0.0 ;
      Igisl = 0.0 ;
    }

    TX = GLPART1 + T1 * GLPART1_dVbse ;
    if ( TX * GLPART1 >= 0.0 ) {
      GLPART1 = TX ;
    } else{
      GLPART1_dVbse = 0.0 ;
      GLPART1_dVdse = 0.0 ;
      GLPART1_dVgse = 0.0 ;
      GLPART1 = 0.0 ;
    }

    TX = Qgod + T1 * Qgod_dVbse ;
    if ( TX * Qgod >= 0.0 ) {
      Qgod = TX ;
    } else {
      Qgod_dVbse = 0.0 ;
      Qgod_dVdse = 0.0 ;
      Qgod_dVgse = 0.0 ;
      Qgod = 0.0 ;
    }

    TX = Qgos + T1 * Qgos_dVbse ;
    if ( TX * Qgos >= 0.0 ) {
      Qgos = TX ;
    } else {
      Qgos_dVbse = 0.0 ;
      Qgos_dVdse = 0.0 ;
      Qgos_dVgse = 0.0 ;
      Qgos = 0.0 ;
    }

    TX = Qgbo + T1 * Qgbo_dVbse ;
    if ( TX * Qgbo >= 0.0 ) {
      Qgbo = TX ;
    } else {
      Qgbo_dVbse = 0.0 ;
      Qgbo_dVdse = 0.0 ;
      Qgbo_dVgse = 0.0 ;
      Qgbo = 0.0 ;
    }
    
    TX = Qy + T1 * Qy_dVbse ;
    if ( TX * Qy >= 0.0 ) {
      Qy = TX ;
    } else {
      Qy_dVbse = 0.0 ;
      Qy_dVdse = 0.0 ;
      Qy_dVgse = 0.0 ;
      Qy = 0.0 ;
    }

    TX = Qdrat  + T1 * Qdrat_dVbse ;
    if ( TX * Qdrat >= 0.0 ) {
      Qdrat = TX ;
    } else{
      Qdrat_dVbse = 0.0 ;
      Qdrat_dVdse = 0.0 ;
      Qdrat_dVgse = 0.0 ;
      Qdrat = 0.0 ;
    }
    
    TX = Qovd + T1 * Qovd_dVbse ;
    if ( TX * Qovd >= 0.0 ) {
      Qovd = TX ;
    } else{
      Qovd_dVbse = 0.0 ;
      Qovd_dVdse = 0.0 ;
      Qovd_dVgse = 0.0 ;
      Qovd = 0.0 ;
    }

    TX = QidLD + T1 * QidLD_dVbse ;
    if ( TX * QidLD >= 0.0 ) {
      QidLD = TX ;
    } else{
      QidLD_dVbse = 0.0 ;
      QidLD_dVdse = 0.0 ;
      QidLD_dVgse = 0.0 ;
      QidLD = 0.0 ;
    }

    TX = QbdLD + T1 * QbdLD_dVbse ;
    if ( TX * QbdLD >= 0.0 ) {
      QbdLD = TX ;
    } else{
      QbdLD_dVbse = 0.0 ;
      QbdLD_dVdse = 0.0 ;
      QbdLD_dVgse = 0.0 ;
      QbdLD = 0.0 ;
    }

    TX = Qovs + T1 * Qovs_dVbse ;
    if ( TX * Qovs >= 0.0 ) {
      Qovs = TX ;
    } else{
      T7 = Qovs / ( Qovs - TX ) ;
      Qovs_dVbse *= T7 ;
      Qovs_dVdse *= T7 ;
      Qovs_dVgse *= T7 ;
      Qovs = 0.0 ;
    }

    TX = QisLD + T1 * QisLD_dVbse ;
    if ( TX * QisLD >= 0.0 ) {
      QisLD = TX ;
    } else{
      QisLD_dVbse = 0.0 ;
      QisLD_dVdse = 0.0 ;
      QisLD_dVgse = 0.0 ;
      QisLD = 0.0 ;
    }

    TX = QbsLD + T1 * QbsLD_dVbse ;
    if ( TX * QbsLD >= 0.0 ) {
      QbsLD = TX ;
    } else{
      QbsLD_dVbse = 0.0 ;
      QbsLD_dVdse = 0.0 ;
      QbsLD_dVgse = 0.0 ;
      QbsLD = 0.0 ;
    }

    if (flg_nqs) { /* for NQS charge */
      TX = Qi_nqs  + T1 * Qi_dVbse_nqs ;
      if ( TX * Qi_nqs >= 0.0 ) {
        Qi_nqs = TX ;
      } else {
        Qi_dVbse_nqs = 0.0 ;
        Qi_dVdse_nqs = 0.0 ;
        Qi_dVgse_nqs = 0.0 ;
        Qi_nqs = 0.0 ;
      }
      
      TX = Qb_nqs  + T1 * Qb_dVbse_nqs ;
      Qb_nqs = TX ;
    }

    if (flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG))) { /* ACNQS */
      TX = tau  + T1 * tau_dVbse ;
      if ( TX * tau >= 0.0 ) {
        tau = TX ;
      } else {
        tau_dVbse = 0.0 ;
        tau_dVdse = 0.0 ;
        tau_dVgse = 0.0 ;
        tau = 0.0 ;
      }

      TX = taub  + T1 * taub_dVbse ;
      if ( TX * taub >= 0.0 ) {
        taub = TX ;
      } else {
        taub_dVbse = 0.0 ;
        taub_dVdse = 0.0 ;
        taub_dVgse = 0.0 ;
        taub = 0.0 ;
      }
    }
  }

  /*-----------------------------------------------------------* 
   * Warn negative conductance.
   * - T1 ( = d Ids / d Vds ) is the derivative w.r.t. circuit bias.
   *-----------------*/ 
  if ( here->HSM2_mode == HiSIM_NORMAL_MODE ) {
    T1 = Ids_dVdse ;
  } else {
    T1 = Ids_dVbse + Ids_dVdse + Ids_dVgse ; /* Ids_dVss * -1 */
  }
  
  if ( flg_info >= 1 && 
       (Ids_dVbse < 0.0 || T1 < 0.0 || Ids_dVgse < 0.0) ) {
    printf( "*** warning(HiSIM): Negative Conductance\n" ) ;
    printf( " type = %d  mode = %d\n" , model->HSM2_type , here->HSM2_mode ) ;
    printf( " Vbse = %12.5e Vdse = %12.5e Vgse = %12.5e\n" , 
            Vbse , Vdse , Vgse ) ;
    printf( " Ids_dVbse   = %12.5e\n" , Ids_dVbse ) ;
    printf( " Ids_dVdse   = %12.5e\n" , T1 ) ;
    printf( " Ids_dVgse   = %12.5e\n" , Ids_dVgse ) ;
  }

  /*-----------------------------------------------------------* 
   * Redefine overlap charges/capacitances.
   *-----------------*/
  /*---------------------------------------------------* 
   * Overlap capacitance.
   *-----------------*/ 
  Cggo = Qgos_dVgse + Qgod_dVgse + Qgbo_dVgse ;
  Cgdo = Qgos_dVdse + Qgod_dVdse ;
  Cgso = - (Qgos_dVbse + Qgod_dVbse + Qgos_dVdse + Qgod_dVdse
            + Qgos_dVgse + Qgod_dVgse) ;
  Cgbo = Qgos_dVbse + Qgod_dVbse + Qgbo_dVbse ;
  
  /*---------------------------------------------------* 
   * Add fringing charge/capacitance to overlap.
   *-----------------*/ 
  Qgod += Qfd ;
  Qgos += Qfs ;
  
  Cggo += 2.0 * Cf ;
  Cgdo += - Cf ;
  Cgso += - Cf ;
  /*-----------------------------------------------------------* 
   * Assign outputs.
   *-----------------*/

  /*---------------------------------------------------*
   * Multiplication factor of a MOSFET instance.
   *-----------------*/
  M = here->HSM2_m ;

  /*---------------------------------------------------*
   * Channel current and conductances. 
   *-----------------*/
  here->HSM2_ids  = M * Ids ;
  here->HSM2_gmbs = M * Ids_dVbse  ;
  here->HSM2_gds  = M * Ids_dVdse ;
  here->HSM2_gm   = M * Ids_dVgse ;

  /*---------------------------------------------------*
   * Overlap capacitances.
   *-----------------*/
  /* Q_dVsx */
  T2 = - ( Qovd_dVbse + Qovd_dVdse + Qovd_dVgse ) ;
  T6 = - ( Qovs_dVbse + Qovs_dVdse + Qovs_dVgse ) ;
  T5 = - ( QbdLD_dVbse + QbdLD_dVdse + QbdLD_dVgse ) ;
  T7 = - ( QbsLD_dVbse + QbsLD_dVdse + QbsLD_dVgse ) ;

  here->HSM2_cgdo = M * (  Cgdo - Qovd_dVdse - Qovs_dVdse ) ;
  here->HSM2_cgso = M * (  Cgso - T2 - T6 ) ;
  here->HSM2_cgbo = M * (  Cgbo - Qovd_dVbse - Qovs_dVbse ) ;

  here->HSM2_cdgo = M * (  - Qgod_dVgse - Cf + QbdLD_dVgse ) ;
  here->HSM2_cddo = M * (  - Qgod_dVdse + Cf + QbdLD_dVdse ) ;
  here->HSM2_cdso = M * (  Qgod_dVbse + Qgod_dVdse + Qgod_dVgse + T5 ) ;

  here->HSM2_csgo = M * (  - Qgos_dVgse - Cf + QbsLD_dVgse ) ;
  here->HSM2_csdo = M * (  - Qgos_dVdse      + QbsLD_dVdse ) ;
  here->HSM2_csso = M * (  Qgos_dVbse + Qgos_dVdse + Qgos_dVgse + Cf + T7 ) ;
  
  /*---------------------------------------------------*
   * Lateral-field-induced capacitance.
   *-----------------*/
  T0 = model->HSM2_qyrat ;
  T1 = 1.0 - T0 ;
  Qys = Qy * T1 ;
  Qys_dVdse = Qy_dVdse * T1 ;
  Qys_dVgse = Qy_dVgse * T1 ;
  Qys_dVbse = Qy_dVbse * T1 ;
  Qy = Qy * T0 ;
  Qy_dVdse = Qy_dVdse * T0 ;
  Qy_dVgse = Qy_dVgse * T0 ;
  Qy_dVbse = Qy_dVbse * T0 ;

  Cqyd = Qy_dVdse ;
  Cqyg = Qy_dVgse ;
  Cqyb = Qy_dVbse ;
  Cqys = - ( Cqyb +  Cqyd + Cqyg ) ;
  here->HSM2_cqyd = M * Cqyd ;
  here->HSM2_cqyg = M * Cqyg ;
  here->HSM2_cqyb = M * Cqyb ;

  /* -------------------------------------*
   * Intrinsic charges / capacitances.
   *-----------------*/
  if ( flg_nqs && ((ckt->CKTmode & MODETRAN) ||
                   (ckt->CKTmode & MODEINITFIX)) ) { /* NQS (tran. analysis) */
    
    *(ckt->CKTstate0 + here->HSM2qi_nqs) = Qi_nqs ;
    *(ckt->CKTstate0 + here->HSM2qb_nqs) = Qb_nqs ;

    here->HSM2_qg = M * - (Qb_nqs + Qi_nqs) ;
    here->HSM2_qd = M * Qi_nqs * Qdrat ;
    here->HSM2_qs = M * Qi_nqs * (1.0 - Qdrat) ;

    here->HSM2_cbgb = M * Qb_dVgse_nqs ;
    here->HSM2_cbdb = M * Qb_dVdse_nqs ;
    here->HSM2_cbsb = M * - (Qb_dVbse_nqs + Qb_dVdse_nqs + Qb_dVgse_nqs) ;

    here->HSM2_cggb = M * ( - Qb_dVgse_nqs - Qi_dVgse_nqs ) ;
    here->HSM2_cgdb = M * ( - Qb_dVdse_nqs - Qi_dVdse_nqs ) ;
    here->HSM2_cgsb = M * ( Qb_dVbse_nqs + Qb_dVdse_nqs + Qb_dVgse_nqs
			    + Qi_dVbse_nqs + Qi_dVdse_nqs + Qi_dVgse_nqs ) ;

    qd_dVgse = Qi_dVgse_nqs * Qdrat + Qdrat_dVgse * Qi_nqs ;
    qd_dVdse = Qi_dVdse_nqs * Qdrat + Qdrat_dVdse * Qi_nqs ;
    qd_dVbse = Qi_dVbse_nqs * Qdrat + Qdrat_dVbse * Qi_nqs ;
    qd_dVsse = - ( qd_dVgse + qd_dVdse + qd_dVbse ) ;
    here->HSM2_cdgb = M * qd_dVgse ;
    here->HSM2_cddb = M * qd_dVdse ;
    here->HSM2_cdsb = M * qd_dVsse ;
  } else  { /* QS or NQS (ac dc analysis) */

    here->HSM2_qg = M * - (Qb + Qi) ;
    here->HSM2_qd = M * Qd ;
    here->HSM2_qs = M * ( Qi - Qd ) ;

    here->HSM2_cbgb = M * Qb_dVgse ;
    here->HSM2_cbdb = M * Qb_dVdse ;
    here->HSM2_cbsb = M * - (Qb_dVbse + Qb_dVdse + Qb_dVgse) ;
      
    here->HSM2_cggb = M * ( - Qb_dVgse - Qi_dVgse ) ;
    here->HSM2_cgdb = M * ( - Qb_dVdse - Qi_dVdse ) ;
    here->HSM2_cgsb = M * ( Qb_dVbse + Qb_dVdse + Qb_dVgse
			    + Qi_dVbse + Qi_dVdse + Qi_dVgse ) ;

    here->HSM2_cdgb = M * Qd_dVgse ;
    here->HSM2_cddb = M * Qd_dVdse ;
    here->HSM2_cdsb = M * - (Qd_dVgse + Qd_dVdse + Qd_dVbse) ;
  }


  /*---------------------------------------------------*
   * Add lateral-field-induced charges/capacitances to intrinsic ones.
   * - NOTE: This function depends on coqy, a control option.
   *-----------------*/
  if ( model->HSM2_coqy == 1 ) {
      here->HSM2_qg += M * (   Qy + Qys ) ;
      here->HSM2_qd += M * ( - Qy       ) ;
      here->HSM2_qs += M * (      - Qys ) ;

      T8 = - ( Qys_dVbse + Qys_dVdse + Qys_dVgse ) ;
      here->HSM2_cggb += M * (   Cqyg + Qys_dVgse ) ;
      here->HSM2_cgdb += M * (   Cqyd + Qys_dVdse ) ;
      here->HSM2_cgsb += M * (   Cqys + T8        ) ;

      here->HSM2_cdgb += M * ( - Cqyg ) ;
      here->HSM2_cddb += M * ( - Cqyd ) ;
      here->HSM2_cdsb += M * ( - Cqys ) ;
  }

  /*---------------------------------------------------*
   * Add S/D overlap charges/capacitances to intrinsic ones.
   * - NOTE: This function depends on coadov, a control option.
   *-----------------*/
  if ( model->HSM2_coadov == 1 ) {
    /* Q_dVsb */
    T0 = - ( Qgbo_dVbse + Qgbo_dVdse + Qgbo_dVgse ) ;
    T1 = - ( Qovd_dVbse + Qovd_dVdse + Qovd_dVgse ) ;
    T3 = - ( Qovs_dVbse + Qovs_dVdse + Qovs_dVgse ) ;
    T4 = - ( QidLD_dVbse + QidLD_dVdse + QidLD_dVgse + QisLD_dVbse + QisLD_dVdse + QisLD_dVgse ) ;
    T5 = - ( QbdLD_dVbse + QbdLD_dVdse + QbdLD_dVgse ) ;
    T7 = - ( Qgod_dVbse + Qgod_dVdse + Qgod_dVgse ) ;

      here->HSM2_qg += M * ( Qgod + Qgos + Qgbo - Qovd - Qovs ) ;
      here->HSM2_qd += M * ( - Qgod + QbdLD ) ;
      here->HSM2_qs += M * ( - Qgos + QbsLD ) ;

      here->HSM2_cbgb += M * ( - Qgbo_dVgse + QidLD_dVgse + QisLD_dVgse ) ;
      here->HSM2_cbdb += M * ( - Qgbo_dVdse + QidLD_dVdse + QisLD_dVdse ) ;
      here->HSM2_cbsb += M * ( - T0         + T4                        ) ;

      here->HSM2_cggb += M * ( Cggo - Qovd_dVgse - Qovs_dVgse ) ;
      here->HSM2_cgdb += M * ( Cgdo - Qovd_dVdse - Qovs_dVdse ) ;
      here->HSM2_cgsb += M * ( Cgso - T1         - T3         ) ;

      here->HSM2_cdgb += M * ( - Qgod_dVgse - Cf + QbdLD_dVgse ) ;
      here->HSM2_cddb += M * ( - Qgod_dVdse + Cf + QbdLD_dVdse ) ;
      here->HSM2_cdsb += M * ( - T7              + T5         ) ;
  }



  /*---------------------------------------------------*
   * tau (channel/bulk charge) for ACNQS.
   *-----------------*/
  if (flg_nqs && (ckt->CKTmode & (MODEDCOP | MODEINITSMSIG))) { 
    here->HSM2_tau = tau ;
    here->HSM2_tau_dVgs = tau_dVgse ;
    here->HSM2_tau_dVds = tau_dVdse ;
    here->HSM2_tau_dVbs = tau_dVbse ;
      
    here->HSM2_taub = taub ;
    here->HSM2_taub_dVgs = taub_dVgse ;
    here->HSM2_taub_dVds = taub_dVdse ;
    here->HSM2_taub_dVbs = taub_dVbse ;
      
    here->HSM2_Xd = Qdrat;
    here->HSM2_Xd_dVgs = Qdrat_dVgse ;
    here->HSM2_Xd_dVds = Qdrat_dVdse ;
    here->HSM2_Xd_dVbs = Qdrat_dVbse ;
   
    here->HSM2_Qb      = M * Qb ;
    here->HSM2_Qb_dVgs = M * Qb_dVgse ;
    here->HSM2_Qb_dVds = M * Qb_dVdse ;
    here->HSM2_Qb_dVbs = M * Qb_dVbse ;

    here->HSM2_Qi      = M * Qi ;
    here->HSM2_Qi_dVgs = M * Qi_dVgse ;
    here->HSM2_Qi_dVds = M * Qi_dVdse ;
    here->HSM2_Qi_dVbs = M * Qi_dVbse ;

    here->HSM2_alpha = Alpha ;
  }

  /*---------------------------------------------------* 
   * Substrate/gate/leak currents.
   *-----------------*/ 

  here->HSM2_isub = M * Isub ;
  here->HSM2_gbbs = M * Isub_dVbse ;
  here->HSM2_gbds = M * Isub_dVdse ;
  here->HSM2_gbgs = M * Isub_dVgse ;
    
  here->HSM2_igb   = M * -Igb ;
  here->HSM2_gigbb = M * -Igb_dVbse ;
  here->HSM2_gigbg = M * -Igb_dVgse ;
  if (here->HSM2_mode == HiSIM_NORMAL_MODE) {
    here->HSM2_gigbd = M * -Igb_dVdse ;
    here->HSM2_gigbs = M * ( Igb_dVbse + Igb_dVdse + Igb_dVgse ) ;
  } else {
    here->HSM2_gigbd = M * ( Igb_dVbse + Igb_dVdse + Igb_dVgse ) ;
    here->HSM2_gigbs = M * -Igb_dVdse ;
  }

  if (here->HSM2_mode == HiSIM_NORMAL_MODE) {

    here->HSM2_igd   = M * ( GLPART1 * Igate - Igd ) ;
    here->HSM2_gigdb = M * ( GLPART1 * Igate_dVbse + GLPART1_dVbse * Igate - Igd_dVbse ) ;
    here->HSM2_gigdd = M * ( GLPART1 * Igate_dVdse + GLPART1_dVdse * Igate - Igd_dVdse ) ;
    here->HSM2_gigdg = M * ( GLPART1 * Igate_dVgse + GLPART1_dVgse * Igate - Igd_dVgse ) ;

  } else {

    T1 = 1.0 - GLPART1 ;
    T1_dVb = - GLPART1_dVbse ;
    T1_dVd = - GLPART1_dVdse ;
    T1_dVg = - GLPART1_dVgse ;
    here->HSM2_igd = M * ( T1 * Igate - Igs ) ;
    here->HSM2_gigdb = M * ( T1 * Igate_dVbse + T1_dVb * Igate - Igs_dVbse ) ;
    here->HSM2_gigdd = M * ( T1 * - ( Igate_dVgse + Igate_dVbse + Igate_dVdse ) 
                     + ( - T1_dVb - T1_dVg - T1_dVd ) * Igate + ( Igs_dVgse + Igs_dVbse + Igs_dVdse ) ) ;
    here->HSM2_gigdg = M * ( T1 * Igate_dVgse + T1_dVg * Igate - Igs_dVgse ) ;
  }
  here->HSM2_gigds   = -(here->HSM2_gigdb + here->HSM2_gigdd + here->HSM2_gigdg) ;

  if (here->HSM2_mode == HiSIM_NORMAL_MODE) {

    T1 = 1.0 - GLPART1 ;
    T1_dVb = - GLPART1_dVbse ;
    T1_dVd = - GLPART1_dVdse ;
    T1_dVg = - GLPART1_dVgse ;
    here->HSM2_igs   = M * ( T1 * Igate - Igs ) ;
    here->HSM2_gigsb = M * ( T1 * Igate_dVbse + T1_dVb * Igate - Igs_dVbse ) ;
    here->HSM2_gigsd = M * ( T1 * Igate_dVdse + T1_dVd * Igate - Igs_dVdse ) ;
    here->HSM2_gigsg = M * ( T1 * Igate_dVgse + T1_dVg * Igate - Igs_dVgse ) ;

  } else {

    here->HSM2_igs   = M * ( GLPART1 * Igate - Igd ) ;
    here->HSM2_gigsb = M * ( GLPART1 * Igate_dVbse + GLPART1_dVbse * Igate - Igd_dVbse ) ;
    here->HSM2_gigsd = M * ( GLPART1 * -(Igate_dVgse + Igate_dVbse + Igate_dVdse) 
                 - Igate * ( GLPART1_dVbse + GLPART1_dVdse + GLPART1_dVgse ) + (Igs_dVgse + Igs_dVbse + Igs_dVdse) ) ;
    here->HSM2_gigsg = M * ( GLPART1 * Igate_dVgse + GLPART1_dVgse * Igate - Igd_dVgse ) ;

  }
  here->HSM2_gigss   = -(here->HSM2_gigsb + here->HSM2_gigsd + here->HSM2_gigsg) ;
  
  here->HSM2_igidl    = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igidl       : M * Igisl ;
  here->HSM2_gigidlbs = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igidl_dVbse : M * Igisl_dVbse ;
  here->HSM2_gigidlds = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igidl_dVdse : M * ( - Igisl_dVbse - Igisl_dVdse - Igisl_dVgse ) ;
  here->HSM2_gigidlgs = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igidl_dVgse : M * Igisl_dVgse ;
  
  here->HSM2_igisl    = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igisl       : M * Igidl ;
  here->HSM2_gigislbd = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igisl_dVbse : M * Igidl_dVbse ;
  here->HSM2_gigislsd = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igisl_dVdse : M * ( - Igidl_dVbse - Igidl_dVdse - Igidl_dVgse ) ;
  here->HSM2_gigislgd = (here->HSM2_mode == HiSIM_NORMAL_MODE) ? M * Igisl_dVgse : M * Igidl_dVgse ;

  /*---------------------------------------------------* 
   * Von, Vdsat.
   *-----------------*/ 
  here->HSM2_von = Vth ;
  here->HSM2_vdsat = Vdsat ;

  /*---------------------------------------------------* 
   * Junction diode.
   *-----------------*/ 
  here->HSM2_ibs = M * Ibs ;
  here->HSM2_ibd = M * Ibd ;
  here->HSM2_gbs = M * Gbse ;
  here->HSM2_gbd = M * Gbde ;
  *(ckt->CKTstate0 + here->HSM2qbs) = M * Qbs ;
  *(ckt->CKTstate0 + here->HSM2qbd) = M * Qbd ;
  here->HSM2_capbs = M * Capbse ;
  here->HSM2_capbd = M * Capbde ;
  
  /*-----------------------------------------------------------* 
   * Warn floating-point exceptions.
   * - Function finite() in libm is called.
   * - Go to start with info==5.
   *-----------------*/ 
  T1 = here->HSM2_ids + here->HSM2_gmbs + here->HSM2_gds + here->HSM2_gm ;
  T1 = T1 + here->HSM2_qd + here->HSM2_cdsb ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf (stderr ,
             "*** warning(HiSIM): FP-exception (PART-1)\n" ) ;
    if ( flg_info >= 1 ) {
      printf ( "*** warning(HiSIM): FP-exception\n") ;
      printf ( "here->HSM2_ids   = %12.5e\n" , here->HSM2_ids ) ;
      printf ( "here->HSM2_gmbs  = %12.5e\n" , here->HSM2_gmbs) ;
      printf ( "here->HSM2_gds   = %12.5e\n" , here->HSM2_gds ) ;
      printf ( "here->HSM2_gm    = %12.5e\n" , here->HSM2_gm  ) ;
      printf ( "here->HSM2_qd    = %12.5e\n" , here->HSM2_qd  ) ;
      printf ( "here->HSM2_cdsb  = %12.5e\n" , here->HSM2_cdsb) ;
    }
  }
  
  T1 = here->HSM2_isub + here->HSM2_gbbs + here->HSM2_gbds + here->HSM2_gbgs ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf (stderr ,
             "*** warning(HiSIM): FP-exception (PART-2)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM): FP-exception\n") ;
    }
  }
  
  T1 = here->HSM2_cgbo + Cgdo + Cgso + Cggo ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf(stderr ,
            "*** warning(HiSIM): FP-exception (PART-3)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM): FP-exception\n") ;
    }
  }
  
  T1 = here->HSM2_ibs + here->HSM2_ibd + here->HSM2_gbs + here->HSM2_gbd ;
  T1 = T1 + *(ckt->CKTstate0 + here->HSM2qbs) + *(ckt->CKTstate0 + here->HSM2qbd) 
          + here->HSM2_capbs + here->HSM2_capbd ;
  if ( ! finite (T1) ) {
    flg_err = 1 ;
    fprintf(stderr ,
            "*** warning(HiSIM): FP-exception (PART-4)\n") ;
    if ( flg_info >= 1 ) {
      printf ("*** warning(HiSIM): FP-exception\n") ;
    }
  }

  /*-----------------------------------------------------------* 
   * Exit for error case.
   *-----------------*/ 
  if ( flg_err != 0 ) {
    fprintf (stderr , "----- bias information (HiSIM)\n" ) ;
    fprintf (stderr , "name: %s\n" , here->HSM2name ) ;
    fprintf (stderr , "stetes: %d\n" , here->HSM2states ) ;
    fprintf (stderr , "vds= %12.5e vgs=%12.5e vbs=%12.5e\n"
            , vds , vgs , vbs ) ;
    fprintf (stderr , "vbs_jct= %12.5e vbd_jct= %12.5e\n"
            , vbs_jct , vbd_jct ) ;
    fprintf (stderr , "vd= %12.5e vg= %12.5e vb= %12.5e vs= %12.5e\n" 
            , *( ckt->CKTrhsOld + here->HSM2dNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2gNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2bNodePrime ) 
            , *( ckt->CKTrhsOld + here->HSM2sNodePrime )  ) ;
    if ( here->HSM2_called >= 1 ) {
      fprintf (stderr , "vdsc_prv= %12.5e vgsc_prv=%12.5e vbsc_prv=%12.5e\n"
              , here->HSM2_vdsc_prv , here->HSM2_vgsc_prv
              , here->HSM2_vbsc_prv ) ;
    }
    fprintf (stderr , "----- bias information (end)\n" ) ;
  }

  if ( flg_err != 0 ) return ( HiSIM_ERROR ) ;

  /*-----------------------------------------------------------* 
   * Noise.
   *-----------------*/ 
  here->HSM2_noiflick = M * Nflic ;
  here->HSM2_noithrml = M * Nthrml ;

  /*----------------------------------------------------------*
   * induced gate noise. ( Part 3/3 )
   *----------------------*/
  if ( model->HSM2_coign != 0 && model->HSM2_cothrml != 0 && flg_ign == 1 && !flg_noqi ) {
    T0 = Cox_small * Cox * here->HSM2_weff_nf * Leff ;
    T1 = here->HSM2_cgsb / M ;
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
      here->HSM2_noiigate = M * Nign0 * kusai_ig * correct_w1 ;
      here->HSM2_noicross = crl_f ;
      if ( here->HSM2_noiigate < 0.0 ) here->HSM2_noiigate = 0.0e0 ;
    }else{
    here->HSM2_noiigate = 0.0e0 ;
    here->HSM2_noicross = 0.0e0 ;
  }
  here->HSM2_Qdrat = Qdrat ; /* needed for calculating induced gate noise */
  }else{
    here->HSM2_noiigate = 0.0e0 ;
    here->HSM2_noicross = 0.0e0 ;
  }
  
  /*-----------------------------------------------------------* 
   * Restore values for next calculation.
   *-----------------*/ 

  /* Confined biases */
  if ( here->HSM2_called >= 1 ) {
    here->HSM2_vbsc_prv2 = here->HSM2_vbsc_prv ;
    here->HSM2_vdsc_prv2 = here->HSM2_vdsc_prv ;
    here->HSM2_vgsc_prv2 = here->HSM2_vgsc_prv ;
    here->HSM2_mode_prv2 = here->HSM2_mode_prv ;
  }
  here->HSM2_vbsc_prv = Vbsc ;
  here->HSM2_vdsc_prv = Vdsc ;
  here->HSM2_vgsc_prv = Vgsc ;
  here->HSM2_mode_prv = here->HSM2_mode ;
  
  /* Surface potentials and derivatives w.r.t. internal biases */
  if ( here->HSM2_called >= 1 ) {
    here->HSM2_ps0_prv2 = here->HSM2_ps0_prv ;
    here->HSM2_ps0_dvbs_prv2 = here->HSM2_ps0_dvbs_prv ;
    here->HSM2_ps0_dvds_prv2 = here->HSM2_ps0_dvds_prv ;
    here->HSM2_ps0_dvgs_prv2 = here->HSM2_ps0_dvgs_prv ;
    here->HSM2_pds_prv2 = here->HSM2_pds_prv ;
    here->HSM2_pds_dvbs_prv2 = here->HSM2_pds_dvbs_prv ;
    here->HSM2_pds_dvds_prv2 = here->HSM2_pds_dvds_prv ;
    here->HSM2_pds_dvgs_prv2 = here->HSM2_pds_dvgs_prv ;
  }

  here->HSM2_ps0_prv = Ps0 ;
  here->HSM2_ps0_dvbs_prv = Ps0_dVbs ;
  here->HSM2_ps0_dvds_prv = Ps0_dVds ;
  here->HSM2_ps0_dvgs_prv = Ps0_dVgs ;
  here->HSM2_pds_prv = Pds ;
  here->HSM2_pds_dvbs_prv = Pds_dVbs ;
  here->HSM2_pds_dvds_prv = Pds_dVds ;
  here->HSM2_pds_dvgs_prv = Pds_dVgs ;

  /* Derivatives of channel current w.r.t. internal biases */
  here->HSM2_ids_prv = Ids ;
  here->HSM2_ids_dvbs_prv = Ids_dVbs ;
  here->HSM2_ids_dvds_prv = Ids_dVds ;
  here->HSM2_ids_dvgs_prv = Ids_dVgs ;
      
  /* For CORECIP = 1 */
  if ( corecip ) {
    here->HSM2_PS0Z_SCE_prv = PS0Z_SCE ;
    here->HSM2_PS0Z_SCE_dvds_prv = PS0Z_SCE_dVds ;
    here->HSM2_PS0Z_SCE_dvgs_prv = PS0Z_SCE_dVgs ;
    here->HSM2_PS0Z_SCE_dvbs_prv = PS0Z_SCE_dVbs ;
    /* here->HSM2_nnn = NNN ; */
  }
      
  /*-----------------------------------------------------------* 
   * End of PART-7. (label) 
   *-----------------*/ 
/* end_of_part_7: */
   
  /*-----------------------------------------------------------* 
   * Bottom of hsm2eval. 
   *-----------------*/ 


  return ( HiSIM_OK ) ;
  
} /* end of hsm2eval */

