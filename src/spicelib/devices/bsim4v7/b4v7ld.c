/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/
/**** OpenMP support ngspice 06/28/2010 ****/
/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4ld.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Authors: 2008- Wenwei Yang, Ali Niknejad, Chenming Hu 
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 11/15/2002.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, 03/04/2004.
 * Modified by Xuemei Xi, Mohan Dunga, 07/29/2005.
 * Modified by Mohan Dunga, 12/13/2006.
 * Modified by Mohan Dunga, Wenwei Yang, 05/18/2007.
 * Modified by Wenwei Yang, 07/31/2008.
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
 **********/


#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h."
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

#define MAX_EXPL 2.688117142e+43
#define MIN_EXPL 3.720075976e-44
#define EXPL_THRESHOLD 100.0

#define MAX_EXP 5.834617425e14
#define MIN_EXP 1.713908431e-15
#define EXP_THRESHOLD 34.0
#define EPS0 8.85418e-12
#define EPSSI 1.03594e-10
#define Charge_q 1.60219e-19
#define DELTA_1 0.02
#define DELTA_2 0.02
#define DELTA_3 0.02
#define DELTA_4 0.02
#define MM  3  /* smooth coeff */
#define DEXP(A,B,C) {                                                         \
        if (A > EXP_THRESHOLD) {                                              \
            B = MAX_EXP*(1.0+(A)-EXP_THRESHOLD);                              \
            C = MAX_EXP;                                                      \
        } else if (A < -EXP_THRESHOLD)  {                                     \
            B = MIN_EXP;                                                      \
            C = 0;                                                            \
        } else   {                                                            \
            B = exp(A);                                                       \
            C = B;                                                            \
        }                                                                     \
    }

#ifdef USE_OMP
int BSIM4v7LoadOMP(BSIM4v7instance *here, CKTcircuit *ckt);
void BSIM4v7LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt);
#endif

int BSIM4v7polyDepletion(double phi, double ngate,double epsgate, double coxe, double Vgs, double *Vgs_eff, double *dVgs_eff_dVg);

int
BSIM4v7load(
GENmodel *inModel,
CKTcircuit *ckt)
{
#ifdef USE_OMP
    int idx;
    BSIM4v7model *model = (BSIM4v7model*)inModel;
    int good = 0;
    BSIM4v7instance *here;
    BSIM4v7instance **InstArray;
    InstArray = model->BSIM4v7InstanceArray;

#pragma omp parallel for private(here)
    for (idx = 0; idx < model->BSIM4v7InstCount; idx++) {
        here = InstArray[idx];
        good = BSIM4v7LoadOMP(here, ckt);
    }

    BSIM4v7LoadRhsMat(inModel, ckt);
    
    return good;
}


int BSIM4v7LoadOMP(BSIM4v7instance *here, CKTcircuit *ckt) {
BSIM4v7model *model;
#else
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;
#endif
double ceqgstot, dgstot_dvd, dgstot_dvg, dgstot_dvs, dgstot_dvb;
double ceqgdtot, dgdtot_dvd, dgdtot_dvg, dgdtot_dvs, dgdtot_dvb;
double gstot, gstotd, gstotg, gstots, gstotb, gspr, Rs, Rd;
double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
double vgs_eff, vgd_eff, dvgs_eff_dvg, dvgd_eff_dvg;
double dRs_dvg, dRd_dvg, dRs_dvb, dRd_dvb;
double dT0_dvg, dT1_dvb, dT3_dvg, dT3_dvb;
double vses, vdes, vdedo, delvses, delvded, delvdes;
double Isestot, cseshat, Idedtot, cdedhat;
#ifndef NEWCONV
double tol0, tol1, tol2, tol3, tol4, tol5, tol6;
#endif

double geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb, ceqgcrg;
double vges, vgms, vgedo, vgmdo, vged, vgmd, delvged, delvgmd;
double delvges, delvgms, vgmb;
double gcgmgmb=0.0, gcgmdb=0.0, gcgmsb=0.0, gcdgmb, gcsgmb;
double gcgmbb=0.0, gcbgmb, qgmb, qgmid=0.0, ceqqgmid;

double vbd, vbs, vds, vgb, vgd, vgs, vgdo;
#ifndef PREDICTOR
double xfact;
#endif
double vdbs, vdbd, vsbs, vsbdo, vsbd;
double delvdbs, delvdbd, delvsbs;
double delvbd_jct, delvbs_jct, vbs_jct, vbd_jct;

double SourceSatCurrent, DrainSatCurrent;
double ag0, qgb, von, cbhat, VgstNVt, ExpVgst;
double ceqqb, ceqqd, ceqqg, ceqqjd=0.0, ceqqjs=0.0, ceq, geq;
double cdrain, cdhat, ceqdrn, ceqbd, ceqbs, ceqjd, ceqjs, gjbd, gjbs;
double czbd, czbdsw, czbdswg, czbs, czbssw, czbsswg, evbd, evbs, arg, sarg;
double delvbd, delvbs, delvds, delvgd, delvgs;
double Vfbeff, dVfbeff_dVg, dVfbeff_dVb, V3, V4;
double gcbdb, gcbgb, gcbsb, gcddb, gcdgb, gcdsb, gcgdb, gcggb, gcgsb, gcsdb;
double gcgbb, gcdbb, gcsbb, gcbbb;
double gcdbdb, gcsbsb;
double gcsgb, gcssb, MJD, MJSWD, MJSWGD, MJS, MJSWS, MJSWGS;
double qgate=0.0, qbulk=0.0, qdrn=0.0, qsrc, cqgate, cqbody, cqdrn;
double Vdb, Vds, Vgs, Vbs, Gmbs, FwdSum, RevSum;
double Igidl, Ggidld, Ggidlg, Ggidlb;
double Voxacc=0.0, dVoxacc_dVg=0.0, dVoxacc_dVb=0.0;
double Voxdepinv=0.0, dVoxdepinv_dVg=0.0, dVoxdepinv_dVd=0.0, dVoxdepinv_dVb=0.0;
double VxNVt=0.0, ExpVxNVt, Vaux=0.0, dVaux_dVg=0.0, dVaux_dVd=0.0, dVaux_dVb=0.0;
double Igc, dIgc_dVg, dIgc_dVd, dIgc_dVb;
double Igcs, dIgcs_dVg, dIgcs_dVd, dIgcs_dVb;
double Igcd, dIgcd_dVg, dIgcd_dVd, dIgcd_dVb;
double Igs, dIgs_dVg, dIgs_dVs, Igd, dIgd_dVg, dIgd_dVd;
double Igbacc, dIgbacc_dVg, dIgbacc_dVb;
double Igbinv, dIgbinv_dVg, dIgbinv_dVd, dIgbinv_dVb;
double Pigcd, dPigcd_dVg, dPigcd_dVd, dPigcd_dVb;
double Istoteq, gIstotg, gIstotd, gIstots, gIstotb;
double Idtoteq, gIdtotg, gIdtotd, gIdtots, gIdtotb;
double Ibtoteq, gIbtotg, gIbtotd, gIbtots, gIbtotb;
double Igtoteq, gIgtotg, gIgtotd, gIgtots, gIgtotb;
double Igstot, cgshat, Igdtot, cgdhat, Igbtot, cgbhat;
double Vgs_eff, Vfb=0.0, Vth_NarrowW;
/* double Vgd_eff, dVgd_eff_dVg;          v4.7.0 */
double Phis, dPhis_dVb, sqrtPhis, dsqrtPhis_dVb, Vth, dVth_dVb, dVth_dVd;
double Vgst, dVgst_dVg, dVgst_dVb, dVgs_eff_dVg, Nvtms, Nvtmd;
double Vtm, Vtm0;
double n, dn_dVb, dn_dVd, voffcv, noff, dnoff_dVd, dnoff_dVb;
double V0, CoxWLcen, QovCox, LINK;
double DeltaPhi, dDeltaPhi_dVg, VgDP, dVgDP_dVg;
double Cox, Tox, Tcen, dTcen_dVg, dTcen_dVd, dTcen_dVb;
double Ccen, Coxeff, dCoxeff_dVd, dCoxeff_dVg, dCoxeff_dVb;
double Denomi, dDenomi_dVg, dDenomi_dVd, dDenomi_dVb;
double ueff, dueff_dVg, dueff_dVd, dueff_dVb; 
double Esat, Vdsat;
double EsatL, dEsatL_dVg, dEsatL_dVd, dEsatL_dVb;
double dVdsat_dVg, dVdsat_dVb, dVdsat_dVd, Vasat, dAlphaz_dVg, dAlphaz_dVb; 
double dVasat_dVg, dVasat_dVb, dVasat_dVd, Va, dVa_dVd, dVa_dVg, dVa_dVb; 
double Vbseff, dVbseff_dVb, VbseffCV, dVbseffCV_dVb; 
double VgsteffVth, dT11_dVg;
double Arg1, One_Third_CoxWL, Two_Third_CoxWL, Alphaz, CoxWL; 
double T0=0.0, dT0_dVg, dT0_dVd, dT0_dVb;
double T1, dT1_dVg, dT1_dVd, dT1_dVb;
double T2, dT2_dVg, dT2_dVd, dT2_dVb;
double T3, dT3_dVg, dT3_dVd, dT3_dVb;
double T4, dT4_dVd, dT4_dVb;
double T5, dT5_dVg, dT5_dVd, dT5_dVb;
double T6, dT6_dVg, dT6_dVd, dT6_dVb;
double T7, dT7_dVg, dT7_dVd, dT7_dVb;
double T8, dT8_dVg, dT8_dVd, dT8_dVb;
double T9, dT9_dVg, dT9_dVd, dT9_dVb;
double T10, dT10_dVg, dT10_dVb, dT10_dVd; 
double T11, T12, T13, T14;
double tmp, Abulk, dAbulk_dVb, Abulk0, dAbulk0_dVb;
double Cclm, dCclm_dVg, dCclm_dVd, dCclm_dVb;
double FP, dFP_dVg, PvagTerm, dPvagTerm_dVg, dPvagTerm_dVd, dPvagTerm_dVb;
double VADITS, dVADITS_dVg, dVADITS_dVd;
double Lpe_Vb, dDITS_Sft_dVb, dDITS_Sft_dVd;
double DITS_Sft2, dDITS_Sft2_dVd;        /* v4.7 New DITS */
double VACLM, dVACLM_dVg, dVACLM_dVd, dVACLM_dVb;
double VADIBL, dVADIBL_dVg, dVADIBL_dVd, dVADIBL_dVb;
double Xdep, dXdep_dVb, lt1, dlt1_dVb, ltw, dltw_dVb, Delt_vth, dDelt_vth_dVb;
double Theta0, dTheta0_dVb;
double TempRatio, tmp1, tmp2, tmp3, tmp4;
double DIBL_Sft, dDIBL_Sft_dVd, Lambda, dLambda_dVg;
double Idtot, Ibtot, a1, ScalingFactor;

double Vgsteff, dVgsteff_dVg, dVgsteff_dVd, dVgsteff_dVb; 
double Vdseff, dVdseff_dVg, dVdseff_dVd, dVdseff_dVb; 
double VdseffCV, dVdseffCV_dVg, dVdseffCV_dVd, dVdseffCV_dVb; 
double diffVds, dAbulk_dVg;
double beta, dbeta_dVg, dbeta_dVd, dbeta_dVb;
double gche, dgche_dVg, dgche_dVd, dgche_dVb;
double fgche1, dfgche1_dVg, dfgche1_dVd, dfgche1_dVb;
double fgche2, dfgche2_dVg, dfgche2_dVd, dfgche2_dVb;
double Idl, dIdl_dVg, dIdl_dVd, dIdl_dVb;
double Idsa, dIdsa_dVg, dIdsa_dVd, dIdsa_dVb;
double Ids, Gm, Gds, Gmb, devbs_dvb, devbd_dvb;
double Isub, Gbd, Gbg, Gbb;
double VASCBE, dVASCBE_dVg, dVASCBE_dVd, dVASCBE_dVb;
double CoxeffWovL;
double Rds, dRds_dVg, dRds_dVb, WVCox, WVCoxRds;
double Vgst2Vtm, VdsatCV;
double Leff, Weff, dWeff_dVg, dWeff_dVb;
double AbulkCV, dAbulkCV_dVb;
double qcheq, qdef, gqdef=0.0, cqdef=0.0, cqcheq=0.0;
double gcqdb=0.0, gcqsb=0.0, gcqgb=0.0, gcqbb=0.0;
double dxpart, sxpart, ggtg, ggtd, ggts, ggtb;
double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
double gbspsp, gbbdp, gbbsp, gbspg, gbspb, gbspdp; 
double gbdpdp, gbdpg, gbdpb, gbdpsp; 
double qgdo, qgso, cgdo, cgso;
double Cgg, Cgd, Cgb, Cdg, Cdd, Cds;
double Csg, Csd, Css, Csb, Cbg, Cbd, Cbb;
double Cgg1, Cgd1, Cgb1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;
double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;
double Igisl, Ggislg, Ggislb, Ggisls;
double Nvtmrss, Nvtmrssws, Nvtmrsswgs;
double Nvtmrsd, Nvtmrsswd, Nvtmrsswgd;

double vs, Fsevl, dvs_dVg, dvs_dVd, dvs_dVb, dFsevl_dVg, dFsevl_dVd, dFsevl_dVb;
double vgdx, vgsx, epssub, toxe, epsrox;
struct bsim4SizeDependParam *pParam;
int ByPass, ChargeComputationNeeded, error, Check, Check1, Check2;

double m;

#ifdef USE_OMP
model = here->BSIM4v7modPtr;
#endif

ScalingFactor = 1.0e-9;
ChargeComputationNeeded =  
                 ((ckt->CKTmode & (MODEDCTRANCURVE | MODEAC | MODETRAN | MODEINITSMSIG)) ||
                 ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)))
                 ? 1 : 0;

#ifndef USE_OMP
for (; model != NULL; model = model->BSIM4v7nextModel)
{    for (here = model->BSIM4v7instances; here != NULL; 
          here = here->BSIM4v7nextInstance)
     {
#endif

          Check = Check1 = Check2 = 1;
          ByPass = 0;
          pParam = here->pParam;

          if ((ckt->CKTmode & MODEINITSMSIG))
          {   vds = *(ckt->CKTstate0 + here->BSIM4v7vds);
              vgs = *(ckt->CKTstate0 + here->BSIM4v7vgs);
              vbs = *(ckt->CKTstate0 + here->BSIM4v7vbs);
              vges = *(ckt->CKTstate0 + here->BSIM4v7vges);
              vgms = *(ckt->CKTstate0 + here->BSIM4v7vgms);
              vdbs = *(ckt->CKTstate0 + here->BSIM4v7vdbs);
              vsbs = *(ckt->CKTstate0 + here->BSIM4v7vsbs);
              vses = *(ckt->CKTstate0 + here->BSIM4v7vses);
              vdes = *(ckt->CKTstate0 + here->BSIM4v7vdes);

              qdef = *(ckt->CKTstate0 + here->BSIM4v7qdef);
          }
          else if ((ckt->CKTmode & MODEINITTRAN))
          {   vds = *(ckt->CKTstate1 + here->BSIM4v7vds);
              vgs = *(ckt->CKTstate1 + here->BSIM4v7vgs);
              vbs = *(ckt->CKTstate1 + here->BSIM4v7vbs);
              vges = *(ckt->CKTstate1 + here->BSIM4v7vges);
              vgms = *(ckt->CKTstate1 + here->BSIM4v7vgms);
              vdbs = *(ckt->CKTstate1 + here->BSIM4v7vdbs);
              vsbs = *(ckt->CKTstate1 + here->BSIM4v7vsbs);
              vses = *(ckt->CKTstate1 + here->BSIM4v7vses);
              vdes = *(ckt->CKTstate1 + here->BSIM4v7vdes);

              qdef = *(ckt->CKTstate1 + here->BSIM4v7qdef);
          }
          else if ((ckt->CKTmode & MODEINITJCT) && !here->BSIM4v7off)
          {   vds = model->BSIM4v7type * here->BSIM4v7icVDS;
              vgs = vges = vgms = model->BSIM4v7type * here->BSIM4v7icVGS;
              vbs = vdbs = vsbs = model->BSIM4v7type * here->BSIM4v7icVBS;
              if (vds > 0.0)
              {   vdes = vds + 0.01;
                  vses = -0.01;
              }
              else if (vds < 0.0)
              {   vdes = vds - 0.01;
                  vses = 0.01;
              }
              else
                  vdes = vses = 0.0;

              qdef = 0.0;

              if ((vds == 0.0) && (vgs == 0.0) && (vbs == 0.0) &&
                  ((ckt->CKTmode & (MODETRAN | MODEAC|MODEDCOP |
                   MODEDCTRANCURVE)) || (!(ckt->CKTmode & MODEUIC))))
              {   vds = 0.1;
                  vdes = 0.11;
                    vses = -0.01;
                  vgs = vges = vgms = model->BSIM4v7type 
                                    * here->BSIM4v7vth0 + 0.1;
                  vbs = vdbs = vsbs = 0.0;
              }
          }
          else if ((ckt->CKTmode & (MODEINITJCT | MODEINITFIX)) && 
                  (here->BSIM4v7off)) 
          {   vds = vgs = vbs = vges = vgms = 0.0;
              vdbs = vsbs = vdes = vses = qdef = 0.0;
          }
          else
          {
#ifndef PREDICTOR
               if ((ckt->CKTmode & MODEINITPRED))
               {   xfact = ckt->CKTdelta / ckt->CKTdeltaOld[1];
                   *(ckt->CKTstate0 + here->BSIM4v7vds) = 
                         *(ckt->CKTstate1 + here->BSIM4v7vds);
                   vds = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vds))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vds)));
                   *(ckt->CKTstate0 + here->BSIM4v7vgs) = 
                         *(ckt->CKTstate1 + here->BSIM4v7vgs);
                   vgs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vgs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vgs)));
                   *(ckt->CKTstate0 + here->BSIM4v7vges) =
                         *(ckt->CKTstate1 + here->BSIM4v7vges);
                   vges = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vges))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vges)));
                   *(ckt->CKTstate0 + here->BSIM4v7vgms) =
                         *(ckt->CKTstate1 + here->BSIM4v7vgms);
                   vgms = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vgms))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vgms)));
                   *(ckt->CKTstate0 + here->BSIM4v7vbs) = 
                         *(ckt->CKTstate1 + here->BSIM4v7vbs);
                   vbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vbs)));
                   *(ckt->CKTstate0 + here->BSIM4v7vbd) = 
                         *(ckt->CKTstate0 + here->BSIM4v7vbs)
                         - *(ckt->CKTstate0 + here->BSIM4v7vds);
                   *(ckt->CKTstate0 + here->BSIM4v7vdbs) =
                         *(ckt->CKTstate1 + here->BSIM4v7vdbs);
                   vdbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vdbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vdbs)));
                   *(ckt->CKTstate0 + here->BSIM4v7vdbd) =
                         *(ckt->CKTstate0 + here->BSIM4v7vdbs)
                         - *(ckt->CKTstate0 + here->BSIM4v7vds);
                   *(ckt->CKTstate0 + here->BSIM4v7vsbs) =
                         *(ckt->CKTstate1 + here->BSIM4v7vsbs);
                   vsbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vsbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vsbs)));
                   *(ckt->CKTstate0 + here->BSIM4v7vses) =
                         *(ckt->CKTstate1 + here->BSIM4v7vses);
                   vses = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vses))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vses)));
                   *(ckt->CKTstate0 + here->BSIM4v7vdes) =
                         *(ckt->CKTstate1 + here->BSIM4v7vdes);
                   vdes = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7vdes))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v7vdes)));

                   *(ckt->CKTstate0 + here->BSIM4v7qdef) =
                         *(ckt->CKTstate1 + here->BSIM4v7qdef);
                   qdef = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v7qdef))
                        -(xfact * (*(ckt->CKTstate2 + here->BSIM4v7qdef)));
               }
               else
               {
#endif /* PREDICTOR */
                   vds = model->BSIM4v7type
                       * (*(ckt->CKTrhsOld + here->BSIM4v7dNodePrime)
                       - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vgs = model->BSIM4v7type
                       * (*(ckt->CKTrhsOld + here->BSIM4v7gNodePrime) 
                       - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vbs = model->BSIM4v7type
                       * (*(ckt->CKTrhsOld + here->BSIM4v7bNodePrime)
                       - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vges = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7gNodeExt)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vgms = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7gNodeMid)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vdbs = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7dbNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vsbs = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7sbNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vses = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7sNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   vdes = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7dNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
                   qdef = model->BSIM4v7type
                        * (*(ckt->CKTrhsOld + here->BSIM4v7qNode));
#ifndef PREDICTOR
               }
#endif /* PREDICTOR */

               vgdo = *(ckt->CKTstate0 + here->BSIM4v7vgs)
                    - *(ckt->CKTstate0 + here->BSIM4v7vds);
               vgedo = *(ckt->CKTstate0 + here->BSIM4v7vges)
                     - *(ckt->CKTstate0 + here->BSIM4v7vds);
               vgmdo = *(ckt->CKTstate0 + here->BSIM4v7vgms)
                     - *(ckt->CKTstate0 + here->BSIM4v7vds);

               vbd = vbs - vds;
               vdbd = vdbs - vds;
               vgd = vgs - vds;
                 vged = vges - vds;
                 vgmd = vgms - vds;

               delvbd = vbd - *(ckt->CKTstate0 + here->BSIM4v7vbd);
               delvdbd = vdbd - *(ckt->CKTstate0 + here->BSIM4v7vdbd);
               delvgd = vgd - vgdo;
               delvged = vged - vgedo;
               delvgmd = vgmd - vgmdo;

               delvds = vds - *(ckt->CKTstate0 + here->BSIM4v7vds);
               delvgs = vgs - *(ckt->CKTstate0 + here->BSIM4v7vgs);
               delvges = vges - *(ckt->CKTstate0 + here->BSIM4v7vges);
               delvgms = vgms - *(ckt->CKTstate0 + here->BSIM4v7vgms);
               delvbs = vbs - *(ckt->CKTstate0 + here->BSIM4v7vbs);
               delvdbs = vdbs - *(ckt->CKTstate0 + here->BSIM4v7vdbs);
               delvsbs = vsbs - *(ckt->CKTstate0 + here->BSIM4v7vsbs);

               delvses = vses - (*(ckt->CKTstate0 + here->BSIM4v7vses));
               vdedo = *(ckt->CKTstate0 + here->BSIM4v7vdes)
                     - *(ckt->CKTstate0 + here->BSIM4v7vds);
               delvdes = vdes - *(ckt->CKTstate0 + here->BSIM4v7vdes);
               delvded = vdes - vds - vdedo;

               delvbd_jct = (!here->BSIM4v7rbodyMod) ? delvbd : delvdbd;
               delvbs_jct = (!here->BSIM4v7rbodyMod) ? delvbs : delvsbs;
               if (here->BSIM4v7mode >= 0)
               {   Idtot = here->BSIM4v7cd + here->BSIM4v7csub - here->BSIM4v7cbd
                         + here->BSIM4v7Igidl;
                   cdhat = Idtot - here->BSIM4v7gbd * delvbd_jct
                         + (here->BSIM4v7gmbs + here->BSIM4v7gbbs + here->BSIM4v7ggidlb) * delvbs
                         + (here->BSIM4v7gm + here->BSIM4v7gbgs + here->BSIM4v7ggidlg) * delvgs 
                         + (here->BSIM4v7gds + here->BSIM4v7gbds + here->BSIM4v7ggidld) * delvds;
                   Ibtot = here->BSIM4v7cbs + here->BSIM4v7cbd 
                         - here->BSIM4v7Igidl - here->BSIM4v7Igisl - here->BSIM4v7csub;
                   cbhat = Ibtot + here->BSIM4v7gbd * delvbd_jct
                         + here->BSIM4v7gbs * delvbs_jct - (here->BSIM4v7gbbs + here->BSIM4v7ggidlb)
                         * delvbs - (here->BSIM4v7gbgs + here->BSIM4v7ggidlg) * delvgs
                         - (here->BSIM4v7gbds + here->BSIM4v7ggidld - here->BSIM4v7ggisls) * delvds 
                         - here->BSIM4v7ggislg * delvgd - here->BSIM4v7ggislb* delvbd;

                   Igstot = here->BSIM4v7Igs + here->BSIM4v7Igcs;
                   cgshat = Igstot + (here->BSIM4v7gIgsg + here->BSIM4v7gIgcsg) * delvgs
                          + here->BSIM4v7gIgcsd * delvds + here->BSIM4v7gIgcsb * delvbs;

                   Igdtot = here->BSIM4v7Igd + here->BSIM4v7Igcd;
                   cgdhat = Igdtot + here->BSIM4v7gIgdg * delvgd + here->BSIM4v7gIgcdg * delvgs
                          + here->BSIM4v7gIgcdd * delvds + here->BSIM4v7gIgcdb * delvbs;

                    Igbtot = here->BSIM4v7Igb;
                   cgbhat = here->BSIM4v7Igb + here->BSIM4v7gIgbg * delvgs + here->BSIM4v7gIgbd
                          * delvds + here->BSIM4v7gIgbb * delvbs;
               }
               else
               {   Idtot = here->BSIM4v7cd + here->BSIM4v7cbd - here->BSIM4v7Igidl; /* bugfix */
                   cdhat = Idtot + here->BSIM4v7gbd * delvbd_jct + here->BSIM4v7gmbs 
                         * delvbd + here->BSIM4v7gm * delvgd 
                         - (here->BSIM4v7gds + here->BSIM4v7ggidls) * delvds 
                         - here->BSIM4v7ggidlg * delvgs - here->BSIM4v7ggidlb * delvbs;
                   Ibtot = here->BSIM4v7cbs + here->BSIM4v7cbd 
                         - here->BSIM4v7Igidl - here->BSIM4v7Igisl - here->BSIM4v7csub;
                   cbhat = Ibtot + here->BSIM4v7gbs * delvbs_jct + here->BSIM4v7gbd 
                         * delvbd_jct - (here->BSIM4v7gbbs + here->BSIM4v7ggislb) * delvbd
                         - (here->BSIM4v7gbgs + here->BSIM4v7ggislg) * delvgd
                         + (here->BSIM4v7gbds + here->BSIM4v7ggisld - here->BSIM4v7ggidls) * delvds
                         - here->BSIM4v7ggidlg * delvgs - here->BSIM4v7ggidlb * delvbs; 

                   Igstot = here->BSIM4v7Igs + here->BSIM4v7Igcd;
                   cgshat = Igstot + here->BSIM4v7gIgsg * delvgs + here->BSIM4v7gIgcdg * delvgd
                          - here->BSIM4v7gIgcdd * delvds + here->BSIM4v7gIgcdb * delvbd;

                   Igdtot = here->BSIM4v7Igd + here->BSIM4v7Igcs;
                   cgdhat = Igdtot + (here->BSIM4v7gIgdg + here->BSIM4v7gIgcsg) * delvgd
                          - here->BSIM4v7gIgcsd * delvds + here->BSIM4v7gIgcsb * delvbd;

                   Igbtot = here->BSIM4v7Igb;
                   cgbhat = here->BSIM4v7Igb + here->BSIM4v7gIgbg * delvgd - here->BSIM4v7gIgbd
                          * delvds + here->BSIM4v7gIgbb * delvbd;
               }

               Isestot = here->BSIM4v7gstot * (*(ckt->CKTstate0 + here->BSIM4v7vses));
               cseshat = Isestot + here->BSIM4v7gstot * delvses
                       + here->BSIM4v7gstotd * delvds + here->BSIM4v7gstotg * delvgs
                       + here->BSIM4v7gstotb * delvbs;

               Idedtot = here->BSIM4v7gdtot * vdedo;
               cdedhat = Idedtot + here->BSIM4v7gdtot * delvded
                       + here->BSIM4v7gdtotd * delvds + here->BSIM4v7gdtotg * delvgs
                       + here->BSIM4v7gdtotb * delvbs;


#ifndef NOBYPASS
               /* Following should be one IF statement, but some C compilers 
                * can't handle that all at once, so we split it into several
                * successive IF's */

               if ((!(ckt->CKTmode & MODEINITPRED)) && (ckt->CKTbypass))
               if ((fabs(delvds) < (ckt->CKTreltol * MAX(fabs(vds),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v7vds))) + ckt->CKTvoltTol)))
               if ((fabs(delvgs) < (ckt->CKTreltol * MAX(fabs(vgs),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v7vgs))) + ckt->CKTvoltTol)))
               if ((fabs(delvbs) < (ckt->CKTreltol * MAX(fabs(vbs),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v7vbs))) + ckt->CKTvoltTol)))
               if ((fabs(delvbd) < (ckt->CKTreltol * MAX(fabs(vbd),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v7vbd))) + ckt->CKTvoltTol)))
               if ((here->BSIM4v7rgateMod == 0) || (here->BSIM4v7rgateMod == 1) 
                         || (fabs(delvges) < (ckt->CKTreltol * MAX(fabs(vges),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v7vges))) + ckt->CKTvoltTol)))
               if ((here->BSIM4v7rgateMod != 3) || (fabs(delvgms) < (ckt->CKTreltol
                   * MAX(fabs(vgms), fabs(*(ckt->CKTstate0 + here->BSIM4v7vgms)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v7rbodyMod) || (fabs(delvdbs) < (ckt->CKTreltol
                   * MAX(fabs(vdbs), fabs(*(ckt->CKTstate0 + here->BSIM4v7vdbs)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v7rbodyMod) || (fabs(delvdbd) < (ckt->CKTreltol
                   * MAX(fabs(vdbd), fabs(*(ckt->CKTstate0 + here->BSIM4v7vdbd)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v7rbodyMod) || (fabs(delvsbs) < (ckt->CKTreltol
                   * MAX(fabs(vsbs), fabs(*(ckt->CKTstate0 + here->BSIM4v7vsbs)))
                   + ckt->CKTvoltTol)))
               if ((!model->BSIM4v7rdsMod) || (fabs(delvses) < (ckt->CKTreltol
                   * MAX(fabs(vses), fabs(*(ckt->CKTstate0 + here->BSIM4v7vses)))
                   + ckt->CKTvoltTol)))
               if ((!model->BSIM4v7rdsMod) || (fabs(delvdes) < (ckt->CKTreltol
                   * MAX(fabs(vdes), fabs(*(ckt->CKTstate0 + here->BSIM4v7vdes)))
                   + ckt->CKTvoltTol)))
               if ((fabs(cdhat - Idtot) < ckt->CKTreltol
                   * MAX(fabs(cdhat), fabs(Idtot)) + ckt->CKTabstol))
               if ((fabs(cbhat - Ibtot) < ckt->CKTreltol
                   * MAX(fabs(cbhat), fabs(Ibtot)) + ckt->CKTabstol))
               if ((!model->BSIM4v7igcMod) || ((fabs(cgshat - Igstot) < ckt->CKTreltol
                   * MAX(fabs(cgshat), fabs(Igstot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v7igcMod) || ((fabs(cgdhat - Igdtot) < ckt->CKTreltol
                   * MAX(fabs(cgdhat), fabs(Igdtot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v7igbMod) || ((fabs(cgbhat - Igbtot) < ckt->CKTreltol
                   * MAX(fabs(cgbhat), fabs(Igbtot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v7rdsMod) || ((fabs(cseshat - Isestot) < ckt->CKTreltol
                   * MAX(fabs(cseshat), fabs(Isestot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v7rdsMod) || ((fabs(cdedhat - Idedtot) < ckt->CKTreltol
                   * MAX(fabs(cdedhat), fabs(Idedtot)) + ckt->CKTabstol)))
               {   vds = *(ckt->CKTstate0 + here->BSIM4v7vds);
                   vgs = *(ckt->CKTstate0 + here->BSIM4v7vgs);
                   vbs = *(ckt->CKTstate0 + here->BSIM4v7vbs);
                   vges = *(ckt->CKTstate0 + here->BSIM4v7vges);
                   vgms = *(ckt->CKTstate0 + here->BSIM4v7vgms);

                   vbd = *(ckt->CKTstate0 + here->BSIM4v7vbd);
                   vdbs = *(ckt->CKTstate0 + here->BSIM4v7vdbs);
                   vdbd = *(ckt->CKTstate0 + here->BSIM4v7vdbd);
                   vsbs = *(ckt->CKTstate0 + here->BSIM4v7vsbs);
                   vses = *(ckt->CKTstate0 + here->BSIM4v7vses);
                   vdes = *(ckt->CKTstate0 + here->BSIM4v7vdes);

                   vgd = vgs - vds;
                   vgb = vgs - vbs;
                   vged = vges - vds;
                   vgmd = vgms - vds;
                   vgmb = vgms - vbs;

                   vbs_jct = (!here->BSIM4v7rbodyMod) ? vbs : vsbs;
                   vbd_jct = (!here->BSIM4v7rbodyMod) ? vbd : vdbd;

/*** qdef should not be kept fixed even if vgs, vds & vbs has converged 
****               qdef = *(ckt->CKTstate0 + here->BSIM4v7qdef);  
***/
                   cdrain = here->BSIM4v7cd;

                   if ((ckt->CKTmode & (MODETRAN | MODEAC)) || 
                       ((ckt->CKTmode & MODETRANOP) && 
                       (ckt->CKTmode & MODEUIC)))
                   {   ByPass = 1;

                       qgate = here->BSIM4v7qgate;
                       qbulk = here->BSIM4v7qbulk;
                       qdrn = here->BSIM4v7qdrn;
                       cgdo = here->BSIM4v7cgdo;
                       qgdo = here->BSIM4v7qgdo;
                       cgso = here->BSIM4v7cgso;
                       qgso = here->BSIM4v7qgso;

                       goto line755;
                    }
                    else
                       goto line850;
               }
#endif /*NOBYPASS*/

               von = here->BSIM4v7von;
               if (*(ckt->CKTstate0 + here->BSIM4v7vds) >= 0.0)
               {   vgs = DEVfetlim(vgs, *(ckt->CKTstate0 + here->BSIM4v7vgs), von);
                   vds = vgs - vgd;
                   vds = DEVlimvds(vds, *(ckt->CKTstate0 + here->BSIM4v7vds));
                   vgd = vgs - vds;
                   if (here->BSIM4v7rgateMod == 3)
                   {   vges = DEVfetlim(vges, *(ckt->CKTstate0 + here->BSIM4v7vges), von);
                       vgms = DEVfetlim(vgms, *(ckt->CKTstate0 + here->BSIM4v7vgms), von);
                       vged = vges - vds;
                       vgmd = vgms - vds;
                   }
                   else if ((here->BSIM4v7rgateMod == 1) || (here->BSIM4v7rgateMod == 2))
                   {   vges = DEVfetlim(vges, *(ckt->CKTstate0 + here->BSIM4v7vges), von);
                       vged = vges - vds;
                   }

                   if (model->BSIM4v7rdsMod)
                   {   vdes = DEVlimvds(vdes, *(ckt->CKTstate0 + here->BSIM4v7vdes));
                       vses = -DEVlimvds(-vses, -(*(ckt->CKTstate0 + here->BSIM4v7vses)));
                   }

               }
               else
               {   vgd = DEVfetlim(vgd, vgdo, von);
                   vds = vgs - vgd;
                   vds = -DEVlimvds(-vds, -(*(ckt->CKTstate0 + here->BSIM4v7vds)));
                   vgs = vgd + vds;

                   if (here->BSIM4v7rgateMod == 3)
                   {   vged = DEVfetlim(vged, vgedo, von);
                       vges = vged + vds;
                       vgmd = DEVfetlim(vgmd, vgmdo, von);
                       vgms = vgmd + vds;
                   }
                   if ((here->BSIM4v7rgateMod == 1) || (here->BSIM4v7rgateMod == 2))
                   {   vged = DEVfetlim(vged, vgedo, von);
                       vges = vged + vds;
                   }

                   if (model->BSIM4v7rdsMod)
                   {   vdes = -DEVlimvds(-vdes, -(*(ckt->CKTstate0 + here->BSIM4v7vdes)));
                       vses = DEVlimvds(vses, *(ckt->CKTstate0 + here->BSIM4v7vses));
                   }
               }

               if (vds >= 0.0)
               {   vbs = DEVpnjlim(vbs, *(ckt->CKTstate0 + here->BSIM4v7vbs),
                                   CONSTvt0, model->BSIM4v7vcrit, &Check);
                   vbd = vbs - vds;
                   if (here->BSIM4v7rbodyMod)
                   {   vdbs = DEVpnjlim(vdbs, *(ckt->CKTstate0 + here->BSIM4v7vdbs),
                                        CONSTvt0, model->BSIM4v7vcrit, &Check1);
                       vdbd = vdbs - vds;
                       vsbs = DEVpnjlim(vsbs, *(ckt->CKTstate0 + here->BSIM4v7vsbs),
                                        CONSTvt0, model->BSIM4v7vcrit, &Check2);
                       if ((Check1 == 0) && (Check2 == 0))
                           Check = 0;
                       else 
                           Check = 1;
                   }
               }
               else
               {   vbd = DEVpnjlim(vbd, *(ckt->CKTstate0 + here->BSIM4v7vbd),
                                   CONSTvt0, model->BSIM4v7vcrit, &Check); 
                   vbs = vbd + vds;
                   if (here->BSIM4v7rbodyMod)
                   {   vdbd = DEVpnjlim(vdbd, *(ckt->CKTstate0 + here->BSIM4v7vdbd),
                                        CONSTvt0, model->BSIM4v7vcrit, &Check1);
                       vdbs = vdbd + vds;
                       vsbdo = *(ckt->CKTstate0 + here->BSIM4v7vsbs)
                             - *(ckt->CKTstate0 + here->BSIM4v7vds);
                       vsbd = vsbs - vds;
                       vsbd = DEVpnjlim(vsbd, vsbdo, CONSTvt0, model->BSIM4v7vcrit, &Check2);
                       vsbs = vsbd + vds;
                       if ((Check1 == 0) && (Check2 == 0))
                           Check = 0;
                       else
                           Check = 1;
                   }
               }
          }

          /* Calculate DC currents and their derivatives */
          vbd = vbs - vds;
          vgd = vgs - vds;
          vgb = vgs - vbs;
          vged = vges - vds;
          vgmd = vgms - vds;
          vgmb = vgms - vbs;
          vdbd = vdbs - vds;

          vbs_jct = (!here->BSIM4v7rbodyMod) ? vbs : vsbs;
          vbd_jct = (!here->BSIM4v7rbodyMod) ? vbd : vdbd;

          /* Source/drain junction diode DC model begins */
          Nvtms = model->BSIM4v7vtm * model->BSIM4v7SjctEmissionCoeff;
/*          if ((here->BSIM4v7Aseff <= 0.0) && (here->BSIM4v7Pseff <= 0.0))
          {   SourceSatCurrent = 1.0e-14;
          } v4.7 */
          if ((here->BSIM4v7Aseff <= 0.0) && (here->BSIM4v7Pseff <= 0.0))
          {   SourceSatCurrent = 0.0;
          }
          else
          {   SourceSatCurrent = here->BSIM4v7Aseff * model->BSIM4v7SjctTempSatCurDensity
                               + here->BSIM4v7Pseff * model->BSIM4v7SjctSidewallTempSatCurDensity
                               + pParam->BSIM4v7weffCJ * here->BSIM4v7nf
                               * model->BSIM4v7SjctGateSidewallTempSatCurDensity;
          }

          if (SourceSatCurrent <= 0.0)
          {   here->BSIM4v7gbs = ckt->CKTgmin;
              here->BSIM4v7cbs = here->BSIM4v7gbs * vbs_jct;
          }
          else
          {   switch(model->BSIM4v7dioMod)
              {   case 0:
                      evbs = exp(vbs_jct / Nvtms);
                      T1 = model->BSIM4v7xjbvs * exp(-(model->BSIM4v7bvs + vbs_jct) / Nvtms);
                      /* WDLiu: Magic T1 in this form; different from BSIM4v7 beta. */
                      here->BSIM4v7gbs = SourceSatCurrent * (evbs + T1) / Nvtms + ckt->CKTgmin;
                      here->BSIM4v7cbs = SourceSatCurrent * (evbs + here->BSIM4v7XExpBVS
                                     - T1 - 1.0) + ckt->CKTgmin * vbs_jct;
                      break;
                  case 1:
                      T2 = vbs_jct / Nvtms;
                      if (T2 < -EXP_THRESHOLD)
                      {   here->BSIM4v7gbs = ckt->CKTgmin;
                          here->BSIM4v7cbs = SourceSatCurrent * (MIN_EXP - 1.0)
                                         + ckt->CKTgmin * vbs_jct;
                      }
                      else if (vbs_jct <= here->BSIM4v7vjsmFwd)
                      {   evbs = exp(T2);
                          here->BSIM4v7gbs = SourceSatCurrent * evbs / Nvtms + ckt->CKTgmin;
                          here->BSIM4v7cbs = SourceSatCurrent * (evbs - 1.0)
                                         + ckt->CKTgmin * vbs_jct;
                      }
                      else
                      {   T0 = here->BSIM4v7IVjsmFwd / Nvtms;
                          here->BSIM4v7gbs = T0 + ckt->CKTgmin;
                          here->BSIM4v7cbs = here->BSIM4v7IVjsmFwd - SourceSatCurrent + T0 
                                         * (vbs_jct - here->BSIM4v7vjsmFwd) + ckt->CKTgmin * vbs_jct;
                      }        
                      break;
                  case 2:
                      if (vbs_jct < here->BSIM4v7vjsmRev)
                      {   T0 = vbs_jct / Nvtms;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbs = MIN_EXP;
                               devbs_dvb = 0.0;
                          }
                          else
                          {    evbs = exp(T0);
                               devbs_dvb = evbs / Nvtms;
                          }

                          T1 = evbs - 1.0;
                          T2 = here->BSIM4v7IVjsmRev + here->BSIM4v7SslpRev
                             * (vbs_jct - here->BSIM4v7vjsmRev);
                          here->BSIM4v7gbs = devbs_dvb * T2 + T1 * here->BSIM4v7SslpRev + ckt->CKTgmin;
                          here->BSIM4v7cbs = T1 * T2 + ckt->CKTgmin * vbs_jct;
                      }         
                      else if (vbs_jct <= here->BSIM4v7vjsmFwd)
                      {   T0 = vbs_jct / Nvtms;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbs = MIN_EXP;
                               devbs_dvb = 0.0;
                          }
                          else
                          {    evbs = exp(T0);
                               devbs_dvb = evbs / Nvtms;
                          }

                          T1 = (model->BSIM4v7bvs + vbs_jct) / Nvtms;
                          if (T1 > EXP_THRESHOLD)
                          {   T2 = MIN_EXP;
                              T3 = 0.0;
                          }
                          else
                          {   T2 = exp(-T1);
                              T3 = -T2 /Nvtms;
                          }
                          here->BSIM4v7gbs = SourceSatCurrent * (devbs_dvb - model->BSIM4v7xjbvs * T3)
                                         + ckt->CKTgmin;
                          here->BSIM4v7cbs = SourceSatCurrent * (evbs + here->BSIM4v7XExpBVS - 1.0
                                         - model->BSIM4v7xjbvs * T2) + ckt->CKTgmin * vbs_jct;
                      }
                      else
                      {   here->BSIM4v7gbs = here->BSIM4v7SslpFwd + ckt->CKTgmin;
                          here->BSIM4v7cbs = here->BSIM4v7IVjsmFwd + here->BSIM4v7SslpFwd * (vbs_jct
                                         - here->BSIM4v7vjsmFwd) + ckt->CKTgmin * vbs_jct;
                      }
                      break;
                  default: break;
              }
          }

          Nvtmd = model->BSIM4v7vtm * model->BSIM4v7DjctEmissionCoeff;
/*          if ((here->BSIM4v7Adeff <= 0.0) && (here->BSIM4v7Pdeff <= 0.0))
          {   DrainSatCurrent = 1.0e-14;
          } v4.7 */
          if ((here->BSIM4v7Adeff <= 0.0) && (here->BSIM4v7Pdeff <= 0.0))
          {   DrainSatCurrent = 0.0;
          }
          else
          {   DrainSatCurrent = here->BSIM4v7Adeff * model->BSIM4v7DjctTempSatCurDensity
                              + here->BSIM4v7Pdeff * model->BSIM4v7DjctSidewallTempSatCurDensity
                              + pParam->BSIM4v7weffCJ * here->BSIM4v7nf
                              * model->BSIM4v7DjctGateSidewallTempSatCurDensity;
          }

          if (DrainSatCurrent <= 0.0)
          {   here->BSIM4v7gbd = ckt->CKTgmin;
              here->BSIM4v7cbd = here->BSIM4v7gbd * vbd_jct;
          }
          else
          {   switch(model->BSIM4v7dioMod)
              {   case 0:
                      evbd = exp(vbd_jct / Nvtmd);
                      T1 = model->BSIM4v7xjbvd * exp(-(model->BSIM4v7bvd + vbd_jct) / Nvtmd);
                      /* WDLiu: Magic T1 in this form; different from BSIM4v7 beta. */
                      here->BSIM4v7gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + ckt->CKTgmin;
                      here->BSIM4v7cbd = DrainSatCurrent * (evbd + here->BSIM4v7XExpBVD
                                     - T1 - 1.0) + ckt->CKTgmin * vbd_jct;
                      break;
                  case 1:
                      T2 = vbd_jct / Nvtmd;
                      if (T2 < -EXP_THRESHOLD)
                      {   here->BSIM4v7gbd = ckt->CKTgmin;
                          here->BSIM4v7cbd = DrainSatCurrent * (MIN_EXP - 1.0)
                                         + ckt->CKTgmin * vbd_jct;
                      }
                      else if (vbd_jct <= here->BSIM4v7vjdmFwd)
                      {   evbd = exp(T2);
                          here->BSIM4v7gbd = DrainSatCurrent * evbd / Nvtmd + ckt->CKTgmin;
                          here->BSIM4v7cbd = DrainSatCurrent * (evbd - 1.0)
                                         + ckt->CKTgmin * vbd_jct;
                      }
                      else
                      {   T0 = here->BSIM4v7IVjdmFwd / Nvtmd;
                          here->BSIM4v7gbd = T0 + ckt->CKTgmin;
                          here->BSIM4v7cbd = here->BSIM4v7IVjdmFwd - DrainSatCurrent + T0
                                         * (vbd_jct - here->BSIM4v7vjdmFwd) + ckt->CKTgmin * vbd_jct;
                      }
                      break;
                  case 2:
                      if (vbd_jct < here->BSIM4v7vjdmRev)
                      {   T0 = vbd_jct / Nvtmd;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbd = MIN_EXP;
                               devbd_dvb = 0.0;
                          }
                          else
                          {    evbd = exp(T0);
                               devbd_dvb = evbd / Nvtmd;
                          }

                          T1 = evbd - 1.0;
                          T2 = here->BSIM4v7IVjdmRev + here->BSIM4v7DslpRev
                             * (vbd_jct - here->BSIM4v7vjdmRev);
                          here->BSIM4v7gbd = devbd_dvb * T2 + T1 * here->BSIM4v7DslpRev + ckt->CKTgmin;
                          here->BSIM4v7cbd = T1 * T2 + ckt->CKTgmin * vbd_jct;
                      }
                      else if (vbd_jct <= here->BSIM4v7vjdmFwd)
                      {   T0 = vbd_jct / Nvtmd;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbd = MIN_EXP;
                               devbd_dvb = 0.0;
                          }
                          else
                          {    evbd = exp(T0);
                               devbd_dvb = evbd / Nvtmd;
                          }

                          T1 = (model->BSIM4v7bvd + vbd_jct) / Nvtmd;
                          if (T1 > EXP_THRESHOLD)
                          {   T2 = MIN_EXP;
                              T3 = 0.0;
                          }
                          else
                          {   T2 = exp(-T1);
                              T3 = -T2 /Nvtmd;
                          }     
                          here->BSIM4v7gbd = DrainSatCurrent * (devbd_dvb - model->BSIM4v7xjbvd * T3)
                                         + ckt->CKTgmin;
                          here->BSIM4v7cbd = DrainSatCurrent * (evbd + here->BSIM4v7XExpBVD - 1.0
                                         - model->BSIM4v7xjbvd * T2) + ckt->CKTgmin * vbd_jct;
                      }
                      else
                      {   here->BSIM4v7gbd = here->BSIM4v7DslpFwd + ckt->CKTgmin;
                          here->BSIM4v7cbd = here->BSIM4v7IVjdmFwd + here->BSIM4v7DslpFwd * (vbd_jct
                                         - here->BSIM4v7vjdmFwd) + ckt->CKTgmin * vbd_jct;
                      }
                      break;
                  default: break;
              }
          } 

           /* trap-assisted tunneling and recombination current for reverse bias  */
          Nvtmrssws = model->BSIM4v7vtm0 * model->BSIM4v7njtsswstemp;
          Nvtmrsswgs = model->BSIM4v7vtm0 * model->BSIM4v7njtsswgstemp;
          Nvtmrss = model->BSIM4v7vtm0 * model->BSIM4v7njtsstemp;
          Nvtmrsswd = model->BSIM4v7vtm0 * model->BSIM4v7njtsswdtemp;
          Nvtmrsswgd = model->BSIM4v7vtm0 * model->BSIM4v7njtsswgdtemp;
          Nvtmrsd = model->BSIM4v7vtm0 * model->BSIM4v7njtsdtemp;

        if ((model->BSIM4v7vtss - vbs_jct) < (model->BSIM4v7vtss * 1e-3))
        { T9 = 1.0e3; 
          T0 = - vbs_jct / Nvtmrss * T9;
          DEXP(T0, T1, T10);
          dT1_dVb = T10 / Nvtmrss * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtss - vbs_jct);
          T0 = -vbs_jct / Nvtmrss * model->BSIM4v7vtss * T9;
          dT0_dVb = model->BSIM4v7vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T1, T10);
          dT1_dVb = T10 * dT0_dVb;
        }

       if ((model->BSIM4v7vtsd - vbd_jct) < (model->BSIM4v7vtsd * 1e-3) )
        { T9 = 1.0e3;
          T0 = -vbd_jct / Nvtmrsd * T9;
          DEXP(T0, T2, T10);
          dT2_dVb = T10 / Nvtmrsd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtsd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsd * model->BSIM4v7vtsd * T9;
          dT0_dVb = model->BSIM4v7vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T2, T10);
          dT2_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v7vtssws - vbs_jct) < (model->BSIM4v7vtssws * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbs_jct / Nvtmrssws * T9;
          DEXP(T0, T3, T10);
          dT3_dVb = T10 / Nvtmrssws * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtssws - vbs_jct);
          T0 = -vbs_jct / Nvtmrssws * model->BSIM4v7vtssws * T9;
          dT0_dVb = model->BSIM4v7vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T3, T10);
          dT3_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v7vtsswd - vbd_jct) < (model->BSIM4v7vtsswd * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbd_jct / Nvtmrsswd * T9;
          DEXP(T0, T4, T10);
          dT4_dVb = T10 / Nvtmrsswd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtsswd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsswd * model->BSIM4v7vtsswd * T9;
          dT0_dVb = model->BSIM4v7vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T4, T10);
          dT4_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v7vtsswgs - vbs_jct) < (model->BSIM4v7vtsswgs * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbs_jct / Nvtmrsswgs * T9;
          DEXP(T0, T5, T10);
          dT5_dVb = T10 / Nvtmrsswgs * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtsswgs - vbs_jct);
          T0 = -vbs_jct / Nvtmrsswgs * model->BSIM4v7vtsswgs * T9;
          dT0_dVb = model->BSIM4v7vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T5, T10);
          dT5_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v7vtsswgd - vbd_jct) < (model->BSIM4v7vtsswgd * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbd_jct / Nvtmrsswgd * T9;
          DEXP(T0, T6, T10);
          dT6_dVb = T10 / Nvtmrsswgd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v7vtsswgd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsswgd * model->BSIM4v7vtsswgd * T9;
          dT0_dVb = model->BSIM4v7vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T6, T10);
          dT6_dVb = T10 * dT0_dVb;
        }

          here->BSIM4v7gbs += here->BSIM4v7SjctTempRevSatCur * dT1_dVb
                                  + here->BSIM4v7SswTempRevSatCur * dT3_dVb
                                  + here->BSIM4v7SswgTempRevSatCur * dT5_dVb; 
          here->BSIM4v7cbs -= here->BSIM4v7SjctTempRevSatCur * (T1 - 1.0)
                                  + here->BSIM4v7SswTempRevSatCur * (T3 - 1.0)
                                  + here->BSIM4v7SswgTempRevSatCur * (T5 - 1.0); 
          here->BSIM4v7gbd += here->BSIM4v7DjctTempRevSatCur * dT2_dVb
                                  + here->BSIM4v7DswTempRevSatCur * dT4_dVb
                                  + here->BSIM4v7DswgTempRevSatCur * dT6_dVb; 
          here->BSIM4v7cbd -= here->BSIM4v7DjctTempRevSatCur * (T2 - 1.0) 
                                  + here->BSIM4v7DswTempRevSatCur * (T4 - 1.0)
                                  + here->BSIM4v7DswgTempRevSatCur * (T6 - 1.0); 

          /* End of diode DC model */

          if (vds >= 0.0)
          {   here->BSIM4v7mode = 1;
              Vds = vds;
              Vgs = vgs;
              Vbs = vbs;
              Vdb = vds - vbs;  /* WDLiu: for GIDL */

          }
          else
          {   here->BSIM4v7mode = -1;
              Vds = -vds;
              Vgs = vgd;
              Vbs = vbd;
              Vdb = -vbs;
          }


         /* dunga */
         if(model->BSIM4v7mtrlMod)
           {
             epsrox = 3.9;
             toxe = model->BSIM4v7eot;
             epssub = EPS0 * model->BSIM4v7epsrsub;
           }
         else
           {
             epsrox = model->BSIM4v7epsrox;
             toxe = model->BSIM4v7toxe;
             epssub = EPSSI;
           }


          T0 = Vbs - here->BSIM4v7vbsc - 0.001;
          T1 = sqrt(T0 * T0 - 0.004 * here->BSIM4v7vbsc);
          if (T0 >= 0.0)
          {   Vbseff = here->BSIM4v7vbsc + 0.5 * (T0 + T1);
              dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
          }
          else
          {   T2 = -0.002 / (T1 - T0);
              Vbseff = here->BSIM4v7vbsc * (1.0 + T2);
              dVbseff_dVb = T2 * here->BSIM4v7vbsc / T1;
          }

        /* JX: Correction to forward body bias  */
          T9 = 0.95 * pParam->BSIM4v7phi;
          T0 = T9 - Vbseff - 0.001;
          T1 = sqrt(T0 * T0 + 0.004 * T9);
          Vbseff = T9 - 0.5 * (T0 + T1);
          dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);
          Phis = pParam->BSIM4v7phi - Vbseff;
          dPhis_dVb = -1.0;
          sqrtPhis = sqrt(Phis);
          dsqrtPhis_dVb = -0.5 / sqrtPhis; 

          Xdep = pParam->BSIM4v7Xdep0 * sqrtPhis / pParam->BSIM4v7sqrtPhi;
          dXdep_dVb = (pParam->BSIM4v7Xdep0 / pParam->BSIM4v7sqrtPhi)
                    * dsqrtPhis_dVb;

          Leff = pParam->BSIM4v7leff;
          Vtm = model->BSIM4v7vtm;
          Vtm0 = model->BSIM4v7vtm0;

          /* Vth Calculation */
          T3 = sqrt(Xdep);
          V0 = pParam->BSIM4v7vbi - pParam->BSIM4v7phi;

          T0 = pParam->BSIM4v7dvt2 * Vbseff;
          if (T0 >= - 0.5)
          {   T1 = 1.0 + T0;
              T2 = pParam->BSIM4v7dvt2;
          }
          else
          {   T4 = 1.0 / (3.0 + 8.0 * T0);
              T1 = (1.0 + 3.0 * T0) * T4; 
              T2 = pParam->BSIM4v7dvt2 * T4 * T4;
          }
          lt1 = model->BSIM4v7factor1 * T3 * T1;
          dlt1_dVb = model->BSIM4v7factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

          T0 = pParam->BSIM4v7dvt2w * Vbseff;
          if (T0 >= - 0.5)
          {   T1 = 1.0 + T0;
              T2 = pParam->BSIM4v7dvt2w;
          }
          else
          {   T4 = 1.0 / (3.0 + 8.0 * T0);
              T1 = (1.0 + 3.0 * T0) * T4; 
              T2 = pParam->BSIM4v7dvt2w * T4 * T4;
          }
          ltw = model->BSIM4v7factor1 * T3 * T1;
          dltw_dVb = model->BSIM4v7factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

          T0 = pParam->BSIM4v7dvt1 * Leff / lt1;
          if (T0 < EXP_THRESHOLD)
          {   T1 = exp(T0);
              T2 = T1 - 1.0;
              T3 = T2 * T2;
              T4 = T3 + 2.0 * T1 * MIN_EXP;
              Theta0 = T1 / T4;
              dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
              dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
          }
          else
          {   Theta0 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
              dTheta0_dVb = 0.0;
          }
          here->BSIM4v7thetavth = pParam->BSIM4v7dvt0 * Theta0;
          Delt_vth = here->BSIM4v7thetavth * V0;
          dDelt_vth_dVb = pParam->BSIM4v7dvt0 * dTheta0_dVb * V0;

          T0 = pParam->BSIM4v7dvt1w * pParam->BSIM4v7weff * Leff / ltw;
          if (T0 < EXP_THRESHOLD)
          {   T1 = exp(T0);
              T2 = T1 - 1.0;
              T3 = T2 * T2;
              T4 = T3 + 2.0 * T1 * MIN_EXP;
              T5 = T1 / T4;
              dT1_dVb = -T0 * T1 * dltw_dVb / ltw; 
              dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
          }
          else
          {   T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
              dT5_dVb = 0.0;
          }
          T0 = pParam->BSIM4v7dvt0w * T5;
          T2 = T0 * V0;
          dT2_dVb = pParam->BSIM4v7dvt0w * dT5_dVb * V0;

          TempRatio =  ckt->CKTtemp / model->BSIM4v7tnom - 1.0;
          T0 = sqrt(1.0 + pParam->BSIM4v7lpe0 / Leff);
          T1 = pParam->BSIM4v7k1ox * (T0 - 1.0) * pParam->BSIM4v7sqrtPhi
             + (pParam->BSIM4v7kt1 + pParam->BSIM4v7kt1l / Leff
             + pParam->BSIM4v7kt2 * Vbseff) * TempRatio;
          Vth_NarrowW = toxe * pParam->BSIM4v7phi
                      / (pParam->BSIM4v7weff + pParam->BSIM4v7w0);

          T3 = here->BSIM4v7eta0 + pParam->BSIM4v7etab * Vbseff;
          if (T3 < 1.0e-4)
          {   T9 = 1.0 / (3.0 - 2.0e4 * T3);
              T3 = (2.0e-4 - T3) * T9;
              T4 = T9 * T9;
          }
          else
          {   T4 = 1.0;
          }
          dDIBL_Sft_dVd = T3 * pParam->BSIM4v7theta0vb0;
          DIBL_Sft = dDIBL_Sft_dVd * Vds;

           Lpe_Vb = sqrt(1.0 + pParam->BSIM4v7lpeb / Leff);

          Vth = model->BSIM4v7type * here->BSIM4v7vth0 + (pParam->BSIM4v7k1ox * sqrtPhis
              - pParam->BSIM4v7k1 * pParam->BSIM4v7sqrtPhi) * Lpe_Vb
              - here->BSIM4v7k2ox * Vbseff - Delt_vth - T2 + (pParam->BSIM4v7k3
              + pParam->BSIM4v7k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

          dVth_dVb = Lpe_Vb * pParam->BSIM4v7k1ox * dsqrtPhis_dVb - here->BSIM4v7k2ox
                   - dDelt_vth_dVb - dT2_dVb + pParam->BSIM4v7k3b * Vth_NarrowW
                   - pParam->BSIM4v7etab * Vds * pParam->BSIM4v7theta0vb0 * T4
                   + pParam->BSIM4v7kt2 * TempRatio;
          dVth_dVd = -dDIBL_Sft_dVd;


          /* Calculate n */
          tmp1 = epssub / Xdep;
          here->BSIM4v7nstar = model->BSIM4v7vtm / Charge_q * (model->BSIM4v7coxe
                           + tmp1 + pParam->BSIM4v7cit);  
          tmp2 = pParam->BSIM4v7nfactor * tmp1;
          tmp3 = pParam->BSIM4v7cdsc + pParam->BSIM4v7cdscb * Vbseff
               + pParam->BSIM4v7cdscd * Vds;
          tmp4 = (tmp2 + tmp3 * Theta0 + pParam->BSIM4v7cit) / model->BSIM4v7coxe;
          if (tmp4 >= -0.5)
          {   n = 1.0 + tmp4;
              dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                     + pParam->BSIM4v7cdscb * Theta0) / model->BSIM4v7coxe;
              dn_dVd = pParam->BSIM4v7cdscd * Theta0 / model->BSIM4v7coxe;
          }
          else
          {   T0 = 1.0 / (3.0 + 8.0 * tmp4);
              n = (1.0 + 3.0 * tmp4) * T0;
              T0 *= T0;
              dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                     + pParam->BSIM4v7cdscb * Theta0) / model->BSIM4v7coxe * T0;
              dn_dVd = pParam->BSIM4v7cdscd * Theta0 / model->BSIM4v7coxe * T0;
          }


          /* Vth correction for Pocket implant */
           if (pParam->BSIM4v7dvtp0 > 0.0)
          {   T0 = -pParam->BSIM4v7dvtp1 * Vds;
              if (T0 < -EXP_THRESHOLD)
              {   T2 = MIN_EXP;
                  dT2_dVd = 0.0;
              }
              else
              {   T2 = exp(T0);
                  dT2_dVd = -pParam->BSIM4v7dvtp1 * T2;
              }

              T3 = Leff + pParam->BSIM4v7dvtp0 * (1.0 + T2);
              dT3_dVd = pParam->BSIM4v7dvtp0 * dT2_dVd;
              if (model->BSIM4v7tempMod < 2)
              {
                T4 = Vtm * log(Leff / T3);
                dT4_dVd = -Vtm * dT3_dVd / T3;
              }
              else
              {
                T4 = model->BSIM4v7vtm0 * log(Leff / T3);
                dT4_dVd = -model->BSIM4v7vtm0 * dT3_dVd / T3;
              }
              dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
              dDITS_Sft_dVb = T4 * dn_dVb;

              Vth -= n * T4;
              dVth_dVd -= dDITS_Sft_dVd;
              dVth_dVb -= dDITS_Sft_dVb;
          }
        
        /* v4.7 DITS_SFT2  */
        if ((pParam->BSIM4v7dvtp4  == 0.0) || (pParam->BSIM4v7dvtp2factor == 0.0)) {
          T0 = 0.0;
            DITS_Sft2 = 0.0;   
        }
        else 
        {
          //T0 = exp(2.0 * pParam->BSIM4v7dvtp4 * Vds);   /* beta code */
          T1 = 2.0 * pParam->BSIM4v7dvtp4 * Vds;
          DEXP(T1, T0, T10);
            DITS_Sft2 = pParam->BSIM4v7dvtp2factor * (T0-1) / (T0+1);   
          //dDITS_Sft2_dVd = pParam->BSIM4v7dvtp2factor * pParam->BSIM4v7dvtp4 * 4.0 * T0 / ((T0+1) * (T0+1));   /* beta code */
          dDITS_Sft2_dVd = pParam->BSIM4v7dvtp2factor * pParam->BSIM4v7dvtp4 * 4.0 * T10 / ((T0+1) * (T0+1));
          Vth -= DITS_Sft2;
          dVth_dVd -= dDITS_Sft2_dVd;
        }
        
        

          here->BSIM4v7von = Vth;

          
          /* Poly Gate Si Depletion Effect */
          T0 = here->BSIM4v7vfb + pParam->BSIM4v7phi;
          if(model->BSIM4v7mtrlMod == 0)                 
            T1 = EPSSI;
          else
            T1 = model->BSIM4v7epsrgate * EPS0;


              BSIM4v7polyDepletion(T0, pParam->BSIM4v7ngate, T1, model->BSIM4v7coxe, vgs, &vgs_eff, &dvgs_eff_dvg);

              BSIM4v7polyDepletion(T0, pParam->BSIM4v7ngate, T1, model->BSIM4v7coxe, vgd, &vgd_eff, &dvgd_eff_dvg);
              
              if(here->BSIM4v7mode>0) {
                      Vgs_eff = vgs_eff;
                      dVgs_eff_dVg = dvgs_eff_dvg;
              } else {
                      Vgs_eff = vgd_eff;
                      dVgs_eff_dVg = dvgd_eff_dvg;
              }
              here->BSIM4v7vgs_eff = vgs_eff;
              here->BSIM4v7vgd_eff = vgd_eff;
              here->BSIM4v7dvgs_eff_dvg = dvgs_eff_dvg;
              here->BSIM4v7dvgd_eff_dvg = dvgd_eff_dvg;


          Vgst = Vgs_eff - Vth;

          /* Calculate Vgsteff */
          T0 = n * Vtm;
          T1 = pParam->BSIM4v7mstar * Vgst;
          T2 = T1 / T0;
          if (T2 > EXP_THRESHOLD)
          {   T10 = T1;
              dT10_dVg = pParam->BSIM4v7mstar * dVgs_eff_dVg;
              dT10_dVd = -dVth_dVd * pParam->BSIM4v7mstar;
              dT10_dVb = -dVth_dVb * pParam->BSIM4v7mstar;
          }
          else if (T2 < -EXP_THRESHOLD)
          {   T10 = Vtm * log(1.0 + MIN_EXP);
              dT10_dVg = 0.0;
              dT10_dVd = T10 * dn_dVd;
              dT10_dVb = T10 * dn_dVb;
              T10 *= n;
          }
          else
          {   ExpVgst = exp(T2);
              T3 = Vtm * log(1.0 + ExpVgst);
              T10 = n * T3;
              dT10_dVg = pParam->BSIM4v7mstar * ExpVgst / (1.0 + ExpVgst);
              dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
              dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
              dT10_dVg *= dVgs_eff_dVg;
          }

          T1 = pParam->BSIM4v7voffcbn - (1.0 - pParam->BSIM4v7mstar) * Vgst;
          T2 = T1 / T0;
          if (T2 < -EXP_THRESHOLD)
          {   T3 = model->BSIM4v7coxe * MIN_EXP / pParam->BSIM4v7cdep0;
              T9 = pParam->BSIM4v7mstar + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
          }
          else if (T2 > EXP_THRESHOLD)
          {   T3 = model->BSIM4v7coxe * MAX_EXP / pParam->BSIM4v7cdep0;
              T9 = pParam->BSIM4v7mstar + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
          }
          else
          {   ExpVgst = exp(T2);
              T3 = model->BSIM4v7coxe / pParam->BSIM4v7cdep0;
              T4 = T3 * ExpVgst;
              T5 = T1 * T4 / T0;
              T9 = pParam->BSIM4v7mstar + n * T4;
              dT9_dVg = T3 * (pParam->BSIM4v7mstar - 1.0) * ExpVgst / Vtm;
              dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
              dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
              dT9_dVg *= dVgs_eff_dVg;
          }
          here->BSIM4v7Vgsteff = Vgsteff = T10 / T9;
          T11 = T9 * T9;
          dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
          dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
          dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

          /* Calculate Effective Channel Geometry */
          T9 = sqrtPhis - pParam->BSIM4v7sqrtPhi;
          Weff = pParam->BSIM4v7weff - 2.0 * (pParam->BSIM4v7dwg * Vgsteff 
               + pParam->BSIM4v7dwb * T9); 
          dWeff_dVg = -2.0 * pParam->BSIM4v7dwg;
          dWeff_dVb = -2.0 * pParam->BSIM4v7dwb * dsqrtPhis_dVb;

          if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
          {   T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
              Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
              T0 *= T0 * 4.0e-16;
              dWeff_dVg *= T0;
              dWeff_dVb *= T0;
          }

          if (model->BSIM4v7rdsMod == 1)
              Rds = dRds_dVg = dRds_dVb = 0.0;
          else
          {   T0 = 1.0 + pParam->BSIM4v7prwg * Vgsteff;
              dT0_dVg = -pParam->BSIM4v7prwg / T0 / T0;
              T1 = pParam->BSIM4v7prwb * T9;
              dT1_dVb = pParam->BSIM4v7prwb * dsqrtPhis_dVb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01); /* 0.01 = 4.0 * 0.05 * 0.05 */
              dT3_dVg = 1.0 + T2 / (T3 - T2);
              dT3_dVb = dT3_dVg * dT1_dVb;
              dT3_dVg *= dT0_dVg;

              T4 = pParam->BSIM4v7rds0 * 0.5;
              Rds = pParam->BSIM4v7rdswmin + T3 * T4;
              dRds_dVg = T4 * dT3_dVg;
              dRds_dVb = T4 * dT3_dVb;

              if (Rds > 0.0)
                  here->BSIM4v7grdsw = 1.0 / Rds* here->BSIM4v7nf; /*4.6.2*/
              else
                  here->BSIM4v7grdsw = 0.0;
          }
          
          /* Calculate Abulk */
          T9 = 0.5 * pParam->BSIM4v7k1ox * Lpe_Vb / sqrtPhis;
          T1 = T9 + here->BSIM4v7k2ox - pParam->BSIM4v7k3b * Vth_NarrowW;
          dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

          T9 = sqrt(pParam->BSIM4v7xj * Xdep);
          tmp1 = Leff + 2.0 * T9;
          T5 = Leff / tmp1; 
          tmp2 = pParam->BSIM4v7a0 * T5;
          tmp3 = pParam->BSIM4v7weff + pParam->BSIM4v7b1; 
          tmp4 = pParam->BSIM4v7b0 / tmp3;
          T2 = tmp2 + tmp4;
          dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
          T6 = T5 * T5;
          T7 = T5 * T6;

          Abulk0 = 1.0 + T1 * T2; 
          dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

          T8 = pParam->BSIM4v7ags * pParam->BSIM4v7a0 * T7;
          dAbulk_dVg = -T1 * T8;
          Abulk = Abulk0 + dAbulk_dVg * Vgsteff; 
          dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb
                     + 3.0 * T1 * dT2_dVb);

          if (Abulk0 < 0.1) /* added to avoid the problems caused by Abulk0 */
          {   T9 = 1.0 / (3.0 - 20.0 * Abulk0);
              Abulk0 = (0.2 - Abulk0) * T9;
              dAbulk0_dVb *= T9 * T9;
          }

          if (Abulk < 0.1)
          {   T9 = 1.0 / (3.0 - 20.0 * Abulk);
              Abulk = (0.2 - Abulk) * T9;
              T10 = T9 * T9;
              dAbulk_dVb *= T10;
              dAbulk_dVg *= T10;
          }
          here->BSIM4v7Abulk = Abulk;

          T2 = pParam->BSIM4v7keta * Vbseff;
          if (T2 >= -0.9)
          {   T0 = 1.0 / (1.0 + T2);
              dT0_dVb = -pParam->BSIM4v7keta * T0 * T0;
          }
          else
          {   T1 = 1.0 / (0.8 + T2);
              T0 = (17.0 + 20.0 * T2) * T1;
              dT0_dVb = -pParam->BSIM4v7keta * T1 * T1;
          }
          dAbulk_dVg *= T0;
          dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
          dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
          Abulk *= T0;
          Abulk0 *= T0;

          /* Mobility calculation */
          if (model->BSIM4v7mtrlMod && model->BSIM4v7mtrlCompatMod == 0)
            T14 = 2.0 * model->BSIM4v7type *(model->BSIM4v7phig - model->BSIM4v7easub - 0.5*model->BSIM4v7Eg0 + 0.45);   
          else
            T14 = 0.0;

          if (model->BSIM4v7mobMod == 0)
            { T0 = Vgsteff + Vth + Vth - T14;
              T2 = pParam->BSIM4v7ua + pParam->BSIM4v7uc * Vbseff;
              T3 = T0 / toxe;
              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v7ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T3 * (T2 + pParam->BSIM4v7ub * T3) + T6;
              T7 = - 2.0 * T6 * T9;
              T11 = T7 * Vth/T12;
              dDenomi_dVg = (T2 + 2.0 * pParam->BSIM4v7ub * T3) / toxe;
              T13 = 2.0 * (dDenomi_dVg + T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + pParam->BSIM4v7uc * T3;
              dDenomi_dVg+= T7;
          }
          else if (model->BSIM4v7mobMod == 1)
          {   T0 = Vgsteff + Vth + Vth - T14;
              T2 = 1.0 + pParam->BSIM4v7uc * Vbseff;
              T3 = T0 / toxe;
              T4 = T3 * (pParam->BSIM4v7ua + pParam->BSIM4v7ub * T3);
              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v7ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T4 * T2 + T6;
              T7 = - 2.0 * T6 * T9;
              T11 = T7 * Vth/T12;
              dDenomi_dVg = (pParam->BSIM4v7ua + 2.0 * pParam->BSIM4v7ub * T3) * T2 / toxe;
              T13 = 2.0 * (dDenomi_dVg + T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + pParam->BSIM4v7uc * T4;
              dDenomi_dVg+= T7;
          }
          else if (model->BSIM4v7mobMod == 2)
          {   T0 = (Vgsteff + here->BSIM4v7vtfbphi1) / toxe;
              T1 = exp(pParam->BSIM4v7eu * log(T0));
              dT1_dVg = T1 * pParam->BSIM4v7eu / T0 / toxe;
              T2 = pParam->BSIM4v7ua + pParam->BSIM4v7uc * Vbseff;
              T3 = T0 / toxe; /*Do we need it?*/

              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v7ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T1 * T2 + T6;
              T7 = - 2.0 * T6 * T9;
              T11 = T7 * Vth/T12;
              dDenomi_dVg = T2 * dT1_dVg + T7;
              T13 = 2.0 * (T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + T1 * pParam->BSIM4v7uc;
          }
          /*high K mobility*/
         else 
      {                                 
                
                             
                   /*univsersal mobility*/
                 T0 = (Vgsteff + here->BSIM4v7vtfbphi1)* 1.0e-8 / toxe/6.0;
             T1 = exp(pParam->BSIM4v7eu * log(T0)); 
                 dT1_dVg = T1 * pParam->BSIM4v7eu * 1.0e-8/ T0 / toxe/6.0;
             T2 = pParam->BSIM4v7ua + pParam->BSIM4v7uc * Vbseff;
                 
                  /*Coulombic*/
                 VgsteffVth = pParam->BSIM4v7VgsteffVth;
                 
                 T10 = exp(pParam->BSIM4v7ucs * log(0.5 + 0.5 * Vgsteff/VgsteffVth));
                          T11 =  pParam->BSIM4v7ud/T10;
                 dT11_dVg = - 0.5 * pParam->BSIM4v7ucs * T11 /(0.5 + 0.5*Vgsteff/VgsteffVth)/VgsteffVth;
                 
                dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = T1 * pParam->BSIM4v7uc;
                
                T5 = T1 * T2 + T11;
          }

                 
         
      
        
          if (T5 >= -0.8)
          {   Denomi = 1.0 + T5;
          }
          else
          {   T9 = 1.0 / (7.0 + 10.0 * T5);
              Denomi = (0.6 + T5) * T9;
              T9 *= T9;
              dDenomi_dVg *= T9;
              dDenomi_dVd *= T9;
              dDenomi_dVb *= T9;
          }
          
        
          here->BSIM4v7ueff = ueff = here->BSIM4v7u0temp / Denomi;
          T9 = -ueff / Denomi;
          dueff_dVg = T9 * dDenomi_dVg;
          dueff_dVd = T9 * dDenomi_dVd;
          dueff_dVb = T9 * dDenomi_dVb;

          /* Saturation Drain Voltage  Vdsat */
          WVCox = Weff * here->BSIM4v7vsattemp * model->BSIM4v7coxe;
          WVCoxRds = WVCox * Rds; 

          Esat = 2.0 * here->BSIM4v7vsattemp / ueff;
          here->BSIM4v7EsatL = EsatL = Esat * Leff;
          T0 = -EsatL /ueff;
          dEsatL_dVg = T0 * dueff_dVg;
          dEsatL_dVd = T0 * dueff_dVd;
          dEsatL_dVb = T0 * dueff_dVb;
  
          /* Sqrt() */
          a1 = pParam->BSIM4v7a1;
          if (a1 == 0.0)
          {   Lambda = pParam->BSIM4v7a2;
              dLambda_dVg = 0.0;
          }
          else if (a1 > 0.0)
          {   T0 = 1.0 - pParam->BSIM4v7a2;
              T1 = T0 - pParam->BSIM4v7a1 * Vgsteff - 0.0001;
              T2 = sqrt(T1 * T1 + 0.0004 * T0);
              Lambda = pParam->BSIM4v7a2 + T0 - 0.5 * (T1 + T2);
              dLambda_dVg = 0.5 * pParam->BSIM4v7a1 * (1.0 + T1 / T2);
          }
          else
          {   T1 = pParam->BSIM4v7a2 + pParam->BSIM4v7a1 * Vgsteff - 0.0001;
              T2 = sqrt(T1 * T1 + 0.0004 * pParam->BSIM4v7a2);
              Lambda = 0.5 * (T1 + T2);
              dLambda_dVg = 0.5 * pParam->BSIM4v7a1 * (1.0 + T1 / T2);
          }

          Vgst2Vtm = Vgsteff + 2.0 * Vtm;
          if (Rds > 0)
          {   tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
              tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
          }
          else
          {   tmp2 = dWeff_dVg / Weff;
              tmp3 = dWeff_dVb / Weff;
          }
          if ((Rds == 0.0) && (Lambda == 1.0))
          {   T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
              tmp1 = 0.0;
              T1 = T0 * T0;
              T2 = Vgst2Vtm * T0;
              T3 = EsatL * Vgst2Vtm;
              Vdsat = T3 * T0;
                           
              dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
              dT0_dVd = -(Abulk * dEsatL_dVd) * T1; 
              dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;   

              dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
              dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
              dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;   
          }
          else
          {   tmp1 = dLambda_dVg / (Lambda * Lambda);
              T9 = Abulk * WVCoxRds;
              T8 = Abulk * T9;
              T7 = Vgst2Vtm * T9;
              T6 = Vgst2Vtm * WVCoxRds;
              T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda); 
              dT0_dVg = 2.0 * (T8 * tmp2 - Abulk * tmp1
                      + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);
             
              dT0_dVb = 2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3)
                      + (1.0 / Lambda - 1.0) * dAbulk_dVb);
              dT0_dVd = 0.0; 
              T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;
             
              dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
                      + Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 3.0 * (T9
                      + T7 * tmp2 + T6 * dAbulk_dVg);
              dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb
                      + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
              dT1_dVd = Abulk * dEsatL_dVd;

              T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
              dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg
                      + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
              dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
              dT2_dVd = Vgst2Vtm * dEsatL_dVd;

              T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
              Vdsat = (T1 - T3) / T0;

              dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg))
                      / T3;
              dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd))
                      / T3;
              dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb))
                      / T3;

              dVdsat_dVg = (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2
                         - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg) / T0;
              dVdsat_dVb = (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2
                         - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb) / T0;
              dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
          }
          here->BSIM4v7vdsat = Vdsat;

          /* Calculate Vdseff */
          T1 = Vdsat - Vds - pParam->BSIM4v7delta;
          dT1_dVg = dVdsat_dVg;
          dT1_dVd = dVdsat_dVd - 1.0;
          dT1_dVb = dVdsat_dVb;

          T2 = sqrt(T1 * T1 + 4.0 * pParam->BSIM4v7delta * Vdsat);
          T0 = T1 / T2;
             T9 = 2.0 * pParam->BSIM4v7delta;
          T3 = T9 / T2;
          dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
          dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
          dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

          if (T1 >= 0.0)
          {   Vdseff = Vdsat - 0.5 * (T1 + T2);
              dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
              dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
              dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
          }
          else
          {   T4 = T9 / (T2 - T1);
              T5 = 1.0 - T4;
              T6 = Vdsat * T4 / (T2 - T1);
              Vdseff = Vdsat * T5;
              dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
              dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
              dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
          }

          if (Vds == 0.0)
          {  Vdseff = 0.0;
             dVdseff_dVg = 0.0;
             dVdseff_dVb = 0.0; 
          }

          if (Vdseff > Vds)
              Vdseff = Vds;
          diffVds = Vds - Vdseff;
          here->BSIM4v7Vdseff = Vdseff;
          
          /* Velocity Overshoot */
        if((model->BSIM4v7lambdaGiven) && (model->BSIM4v7lambda > 0.0) )
        {  
          T1 =  Leff * ueff;
          T2 = pParam->BSIM4v7lambda / T1;
          T3 = -T2 / T1 * Leff;
          dT2_dVd = T3 * dueff_dVd;
          dT2_dVg = T3 * dueff_dVg;
          dT2_dVb = T3 * dueff_dVb;
          T5 = 1.0 / (Esat * pParam->BSIM4v7litl);
          T4 = -T5 / EsatL;
          dT5_dVg = dEsatL_dVg * T4;
          dT5_dVd = dEsatL_dVd * T4; 
          dT5_dVb = dEsatL_dVb * T4; 
          T6 = 1.0 + diffVds  * T5;
          dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
          dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
          dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
          T7 = 2.0 / (T6 * T6 + 1.0);
          T8 = 1.0 - T7;
          T9 = T6 * T7 * T7;
          dT8_dVg = T9 * dT6_dVg;
          dT8_dVd = T9 * dT6_dVd;
          dT8_dVb = T9 * dT6_dVb;
          T10 = 1.0 + T2 * T8;
          dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
          dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
          dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
          if(T10 == 1.0)
                dT10_dVg = dT10_dVd = dT10_dVb = 0.0;

          dEsatL_dVg *= T10;
          dEsatL_dVg += EsatL * dT10_dVg;
          dEsatL_dVd *= T10;
          dEsatL_dVd += EsatL * dT10_dVd;
          dEsatL_dVb *= T10;
          dEsatL_dVb += EsatL * dT10_dVb;
          EsatL *= T10;
          Esat = EsatL / Leff;  /* bugfix by Wenwei Yang (4.6.4) */
          here->BSIM4v7EsatL = EsatL;
        }

          /* Calculate Vasat */
          tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
          T9 = WVCoxRds * Vgsteff;
          T8 = T9 / Vgst2Vtm;
          T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;
         
          T7 = 2.0 * WVCoxRds * tmp4;
          dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
                  - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm
                  + Vdsat * dAbulk_dVg);   
                  
          dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
                  - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
          dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

          T9 = WVCoxRds * Abulk; 
          T1 = 2.0 / Lambda - 1.0 + T9; 
          dT1_dVg = -2.0 * tmp1 +  WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
          dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

          Vasat = T0 / T1;
          dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
          dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
          dVasat_dVd = dT0_dVd / T1;

          /* Calculate Idl first */
              
          tmp1 = here->BSIM4v7vtfbphi2;
          tmp2 = 2.0e8 * model->BSIM4v7toxp;
          dT0_dVg = 1.0 / tmp2;
          T0 = (Vgsteff + tmp1) * dT0_dVg;

          tmp3 = exp(model->BSIM4v7bdos * 0.7 * log(T0));
          T1 = 1.0 + tmp3;
          T2 = model->BSIM4v7bdos * 0.7 * tmp3 / T0;
          Tcen = model->BSIM4v7ados * 1.9e-9 / T1;
          dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

          Coxeff = epssub * model->BSIM4v7coxp
                 / (epssub + model->BSIM4v7coxp * Tcen);
          here->BSIM4v7Coxeff = Coxeff;
          dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

          CoxeffWovL = Coxeff * Weff / Leff;
          beta = ueff * CoxeffWovL;
          T3 = ueff / Leff;
          dbeta_dVg = CoxeffWovL * dueff_dVg + T3
                    * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
          dbeta_dVd = CoxeffWovL * dueff_dVd;
          dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

          here->BSIM4v7AbovVgst2Vtm = Abulk / Vgst2Vtm;
          T0 = 1.0 - 0.5 * Vdseff * here->BSIM4v7AbovVgst2Vtm;
          dT0_dVg = -0.5 * (Abulk * dVdseff_dVg
                  - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg) / Vgst2Vtm;
          dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
          dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff)
                  / Vgst2Vtm;

          fgche1 = Vgsteff * T0;
          dfgche1_dVg = Vgsteff * dT0_dVg + T0;
          dfgche1_dVd = Vgsteff * dT0_dVd;
          dfgche1_dVb = Vgsteff * dT0_dVb;

          T9 = Vdseff / EsatL;
          fgche2 = 1.0 + T9;
          dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
          dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
          dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

          gche = beta * fgche1 / fgche2;
          dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg
                    - gche * dfgche2_dVg) / fgche2;
          dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd
                    - gche * dfgche2_dVd) / fgche2;
          dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb
                    - gche * dfgche2_dVb) / fgche2;

          T0 = 1.0 + gche * Rds;
          Idl = gche / T0;
          T1 = (1.0 - Idl * Rds) / T0;
          T2 = Idl * Idl;
          dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
          dIdl_dVd = T1 * dgche_dVd;
          dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

          /* Calculate degradation factor due to pocket implant */

          if (pParam->BSIM4v7fprout <= 0.0)
          {   FP = 1.0;
              dFP_dVg = 0.0;
          }
          else
          {   T9 = pParam->BSIM4v7fprout * sqrt(Leff) / Vgst2Vtm;
              FP = 1.0 / (1.0 + T9);
              dFP_dVg = FP * FP * T9 / Vgst2Vtm;
          }

          /* Calculate VACLM */
          T8 = pParam->BSIM4v7pvag / EsatL;
          T9 = T8 * Vgsteff;
          if (T9 > -0.9)
          {   PvagTerm = 1.0 + T9;
              dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
              dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
              dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
          }
          else
          {   T4 = 1.0 / (17.0 + 20.0 * T9);
              PvagTerm = (0.8 + T9) * T4;
              T4 *= T4;
              dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
              T9 *= T4 / EsatL;
              dPvagTerm_dVb = -T9 * dEsatL_dVb;
              dPvagTerm_dVd = -T9 * dEsatL_dVd;
          }

          if ((pParam->BSIM4v7pclm > MIN_EXP) && (diffVds > 1.0e-10))
          {   T0 = 1.0 + Rds * Idl;
              dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
              dT0_dVd = Rds * dIdl_dVd;
              dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

              T2 = Vdsat / Esat;
              T1 = Leff + T2;
              dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
              dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
              dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

              Cclm = FP * PvagTerm * T0 * T1 / (pParam->BSIM4v7pclm * pParam->BSIM4v7litl);
              dCclm_dVg = Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm
                        + dT0_dVg / T0 + dT1_dVg / T1);
              dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0
                        + dT1_dVb / T1);
              dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0
                        + dT1_dVd / T1);
              VACLM = Cclm * diffVds;

              dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
              dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
              dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
          }
          else
          {   VACLM = Cclm = MAX_EXP;
              dVACLM_dVd = dVACLM_dVg = dVACLM_dVb = 0.0;
              dCclm_dVd = dCclm_dVg = dCclm_dVb = 0.0;
          }

          /* Calculate VADIBL */
          if (pParam->BSIM4v7thetaRout > MIN_EXP)
          {   T8 = Abulk * Vdsat;
              T0 = Vgst2Vtm * T8;
              dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8
                      + Vgst2Vtm * Vdsat * dAbulk_dVg;
              dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
              dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

              T1 = Vgst2Vtm + T8;
              dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
              dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
              dT1_dVd = Abulk * dVdsat_dVd;

              T9 = T1 * T1;
              T2 = pParam->BSIM4v7thetaRout;
              VADIBL = (Vgst2Vtm - T0 / T1) / T2;
              dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
              dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
              dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

              T7 = pParam->BSIM4v7pdiblb * Vbseff;
              if (T7 >= -0.9)
              {   T3 = 1.0 / (1.0 + T7);
                  VADIBL *= T3;
                  dVADIBL_dVg *= T3;
                  dVADIBL_dVb = (dVADIBL_dVb - VADIBL * pParam->BSIM4v7pdiblb)
                              * T3;
                  dVADIBL_dVd *= T3;
              }
              else
              {   T4 = 1.0 / (0.8 + T7);
                  T3 = (17.0 + 20.0 * T7) * T4;
                  dVADIBL_dVg *= T3;
                  dVADIBL_dVb = dVADIBL_dVb * T3
                              - VADIBL * pParam->BSIM4v7pdiblb * T4 * T4;
                  dVADIBL_dVd *= T3;
                  VADIBL *= T3;
              }

              dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
              dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
              dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
              VADIBL *= PvagTerm;
          }
          else
          {   VADIBL = MAX_EXP;
              dVADIBL_dVd = dVADIBL_dVg = dVADIBL_dVb = 0.0;
          }

          /* Calculate Va */
          Va = Vasat + VACLM;
          dVa_dVg = dVasat_dVg + dVACLM_dVg;
          dVa_dVb = dVasat_dVb + dVACLM_dVb;
          dVa_dVd = dVasat_dVd + dVACLM_dVd;

          /* Calculate VADITS */
          T0 = pParam->BSIM4v7pditsd * Vds;
          if (T0 > EXP_THRESHOLD)
          {   T1 = MAX_EXP;
              dT1_dVd = 0;
          }
          else
          {   T1 = exp(T0);
              dT1_dVd = T1 * pParam->BSIM4v7pditsd;
          }

          if (pParam->BSIM4v7pdits > MIN_EXP)
          {   T2 = 1.0 + model->BSIM4v7pditsl * Leff;
              VADITS = (1.0 + T2 * T1) / pParam->BSIM4v7pdits;
              dVADITS_dVg = VADITS * dFP_dVg;
              dVADITS_dVd = FP * T2 * dT1_dVd / pParam->BSIM4v7pdits;
              VADITS *= FP;
          }
          else
          {   VADITS = MAX_EXP;
              dVADITS_dVg = dVADITS_dVd = 0;
          }

          /* Calculate VASCBE */
          if ((pParam->BSIM4v7pscbe2 > 0.0)&&(pParam->BSIM4v7pscbe1>=0.0)) /*4.6.2*/
          {   if (diffVds > pParam->BSIM4v7pscbe1 * pParam->BSIM4v7litl
                  / EXP_THRESHOLD)
              {   T0 =  pParam->BSIM4v7pscbe1 * pParam->BSIM4v7litl / diffVds;
                  VASCBE = Leff * exp(T0) / pParam->BSIM4v7pscbe2;
                  T1 = T0 * VASCBE / diffVds;
                  dVASCBE_dVg = T1 * dVdseff_dVg;
                  dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                  dVASCBE_dVb = T1 * dVdseff_dVb;
              }
              else
              {   VASCBE = MAX_EXP * Leff/pParam->BSIM4v7pscbe2;
                  dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
              }
          }
          else
          {   VASCBE = MAX_EXP;
              dVASCBE_dVg = dVASCBE_dVd = dVASCBE_dVb = 0.0;
          }

          /* Add DIBL to Ids */
          T9 = diffVds / VADIBL;
          T0 = 1.0 + T9;
          Idsa = Idl * T0;
          dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
          dIdsa_dVd = T0 * dIdl_dVd + Idl
                    * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
          dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

          /* Add DITS to Ids */
          T9 = diffVds / VADITS;
          T0 = 1.0 + T9;
          dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
          dIdsa_dVd = T0 * dIdsa_dVd + Idsa 
                    * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
          dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
          Idsa *= T0;

          /* Add CLM to Ids */
          T0 = log(Va / Vasat);
          dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
          dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
          dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
          T1 = T0 / Cclm;
          T9 = 1.0 + T1;
          dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
          dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
          dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

          dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
          dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
          dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
          Idsa *= T9;

          /* Substrate current begins */
          tmp = pParam->BSIM4v7alpha0 + pParam->BSIM4v7alpha1 * Leff;
          if ((tmp <= 0.0) || (pParam->BSIM4v7beta0 <= 0.0))
          {   Isub = Gbd = Gbb = Gbg = 0.0;
          }
          else
          {   T2 = tmp / Leff;
              if (diffVds > pParam->BSIM4v7beta0 / EXP_THRESHOLD)
              {   T0 = -pParam->BSIM4v7beta0 / diffVds;
                  T1 = T2 * diffVds * exp(T0);
                  T3 = T1 / diffVds * (T0 - 1.0);
                  dT1_dVg = T3 * dVdseff_dVg;
                  dT1_dVd = T3 * (dVdseff_dVd - 1.0);
                  dT1_dVb = T3 * dVdseff_dVb;
              }
              else
              {   T3 = T2 * MIN_EXP;
                  T1 = T3 * diffVds;
                  dT1_dVg = -T3 * dVdseff_dVg;
                  dT1_dVd = T3 * (1.0 - dVdseff_dVd);
                  dT1_dVb = -T3 * dVdseff_dVb;
              }
              T4 = Idsa * Vdseff;
              Isub = T1 * T4;
              Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg)
                  + T4 * dT1_dVg;
              Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd)
                  + T4 * dT1_dVd;
              Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb)
                  + T4 * dT1_dVb;

              Gbd += Gbg * dVgsteff_dVd;
              Gbb += Gbg * dVgsteff_dVb;
              Gbg *= dVgsteff_dVg;
              Gbb *= dVbseff_dVb;
          }
          here->BSIM4v7csub = Isub;
          here->BSIM4v7gbbs = Gbb;
          here->BSIM4v7gbgs = Gbg;
          here->BSIM4v7gbds = Gbd;

          /* Add SCBE to Ids */
          T9 = diffVds / VASCBE;
          T0 = 1.0 + T9;
          Ids = Idsa * T0;

          Gm = T0 * dIdsa_dVg - Idsa 
             * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
          Gds = T0 * dIdsa_dVd + Idsa 
              * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
          Gmb = T0 * dIdsa_dVb - Idsa
              * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;


          tmp1 = Gds + Gm * dVgsteff_dVd;
          tmp2 = Gmb + Gm * dVgsteff_dVb;
          tmp3 = Gm;

          Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
          Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd)
              + Vdseff * tmp1;
          Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb)
              + Vdseff * tmp2) * dVbseff_dVb;

          cdrain = Ids * Vdseff;

          /* Source End Velocity Limit  */
        if((model->BSIM4v7vtlGiven) && (model->BSIM4v7vtl > 0.0) ) {
          T12 = 1.0 / Leff / CoxeffWovL;
          T11 = T12 / Vgsteff;
          T10 = -T11 / Vgsteff;
          vs = cdrain * T11; /* vs */
          dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
          dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
          dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
          T0 = 2 * MM;
          T1 = vs / (pParam->BSIM4v7vtl * pParam->BSIM4v7tfactor);
          if(T1 > 0.0)  
          {        T2 = 1.0 + exp(T0 * log(T1));
                  T3 = (T2 - 1.0) * T0 / vs; 
                  Fsevl = 1.0 / exp(log(T2)/ T0);
                  dT2_dVg = T3 * dvs_dVg;
                  dT2_dVd = T3 * dvs_dVd;
                  dT2_dVb = T3 * dvs_dVb;
                  T4 = -1.0 / T0 * Fsevl / T2;
                  dFsevl_dVg = T4 * dT2_dVg;
                  dFsevl_dVd = T4 * dT2_dVd;
                  dFsevl_dVb = T4 * dT2_dVb;
          } else {
                  Fsevl = 1.0;
                  dFsevl_dVg = 0.0;
                  dFsevl_dVd = 0.0;
                  dFsevl_dVb = 0.0;
          }
          Gm *=Fsevl;
          Gm += cdrain * dFsevl_dVg;
          Gmb *=Fsevl;
          Gmb += cdrain * dFsevl_dVb;
          Gds *=Fsevl;
          Gds += cdrain * dFsevl_dVd;

          cdrain *= Fsevl; 
        } 

          here->BSIM4v7gds = Gds;
          here->BSIM4v7gm = Gm;
          here->BSIM4v7gmbs = Gmb;
          here->BSIM4v7IdovVds = Ids;
          if( here->BSIM4v7IdovVds <= 1.0e-9) here->BSIM4v7IdovVds = 1.0e-9;

          /* Calculate Rg */
          if ((here->BSIM4v7rgateMod > 1) ||
              (here->BSIM4v7trnqsMod != 0) || (here->BSIM4v7acnqsMod != 0))
          {   T9 = pParam->BSIM4v7xrcrg2 * model->BSIM4v7vtm;
              T0 = T9 * beta;
              dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
              dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
              dT0_dVg = dbeta_dVg * T9;

              here->BSIM4v7gcrg = pParam->BSIM4v7xrcrg1 * ( T0 + Ids);
              here->BSIM4v7gcrgd = pParam->BSIM4v7xrcrg1 * (dT0_dVd + tmp1);
              here->BSIM4v7gcrgb = pParam->BSIM4v7xrcrg1 * (dT0_dVb + tmp2)
                                * dVbseff_dVb;        
              here->BSIM4v7gcrgg = pParam->BSIM4v7xrcrg1 * (dT0_dVg + tmp3)
                               * dVgsteff_dVg;

              if (here->BSIM4v7nf != 1.0)
              {   here->BSIM4v7gcrg *= here->BSIM4v7nf; 
                  here->BSIM4v7gcrgg *= here->BSIM4v7nf;
                  here->BSIM4v7gcrgd *= here->BSIM4v7nf;
                  here->BSIM4v7gcrgb *= here->BSIM4v7nf;
              }

              if (here->BSIM4v7rgateMod == 2)
              {   T10 = here->BSIM4v7grgeltd * here->BSIM4v7grgeltd;
                  T11 = here->BSIM4v7grgeltd + here->BSIM4v7gcrg;
                  here->BSIM4v7gcrg = here->BSIM4v7grgeltd * here->BSIM4v7gcrg / T11;
                  T12 = T10 / T11 / T11;
                  here->BSIM4v7gcrgg *= T12;
                  here->BSIM4v7gcrgd *= T12;
                  here->BSIM4v7gcrgb *= T12;
              }
              here->BSIM4v7gcrgs = -(here->BSIM4v7gcrgg + here->BSIM4v7gcrgd
                               + here->BSIM4v7gcrgb);
          }


          /* Calculate bias-dependent external S/D resistance */
          if (model->BSIM4v7rdsMod)
          {   /* Rs(V) */
              T0 = vgs - pParam->BSIM4v7vfbsd;
              T1 = sqrt(T0 * T0 + 1.0e-4);
              vgs_eff = 0.5 * (T0 + T1);
              dvgs_eff_dvg = vgs_eff / T1;

              T0 = 1.0 + pParam->BSIM4v7prwg * vgs_eff;
              dT0_dvg = -pParam->BSIM4v7prwg / T0 / T0 * dvgs_eff_dvg;
              T1 = -pParam->BSIM4v7prwb * vbs;
              dT1_dvb = -pParam->BSIM4v7prwb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01);
              dT3_dvg = T3 / (T3 - T2);
              dT3_dvb = dT3_dvg * dT1_dvb;
              dT3_dvg *= dT0_dvg;

              T4 = pParam->BSIM4v7rs0 * 0.5;
              Rs = pParam->BSIM4v7rswmin + T3 * T4;
              dRs_dvg = T4 * dT3_dvg;
              dRs_dvb = T4 * dT3_dvb;

              T0 = 1.0 + here->BSIM4v7sourceConductance * Rs;
              here->BSIM4v7gstot = here->BSIM4v7sourceConductance / T0;
              T0 = -here->BSIM4v7gstot * here->BSIM4v7gstot;
              dgstot_dvd = 0.0; /* place holder */
              dgstot_dvg = T0 * dRs_dvg;
              dgstot_dvb = T0 * dRs_dvb;
              dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

              /* Rd(V) */
              T0 = vgd - pParam->BSIM4v7vfbsd;
              T1 = sqrt(T0 * T0 + 1.0e-4);
              vgd_eff = 0.5 * (T0 + T1);
              dvgd_eff_dvg = vgd_eff / T1;

              T0 = 1.0 + pParam->BSIM4v7prwg * vgd_eff;
              dT0_dvg = -pParam->BSIM4v7prwg / T0 / T0 * dvgd_eff_dvg;
              T1 = -pParam->BSIM4v7prwb * vbd;
              dT1_dvb = -pParam->BSIM4v7prwb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01);
              dT3_dvg = T3 / (T3 - T2);
              dT3_dvb = dT3_dvg * dT1_dvb;
              dT3_dvg *= dT0_dvg;

              T4 = pParam->BSIM4v7rd0 * 0.5;
              Rd = pParam->BSIM4v7rdwmin + T3 * T4;
              dRd_dvg = T4 * dT3_dvg;
              dRd_dvb = T4 * dT3_dvb;

              T0 = 1.0 + here->BSIM4v7drainConductance * Rd;
              here->BSIM4v7gdtot = here->BSIM4v7drainConductance / T0;
              T0 = -here->BSIM4v7gdtot * here->BSIM4v7gdtot;
              dgdtot_dvs = 0.0;
              dgdtot_dvg = T0 * dRd_dvg;
              dgdtot_dvb = T0 * dRd_dvb;
              dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

              here->BSIM4v7gstotd = vses * dgstot_dvd;
              here->BSIM4v7gstotg = vses * dgstot_dvg;
              here->BSIM4v7gstots = vses * dgstot_dvs;
              here->BSIM4v7gstotb = vses * dgstot_dvb;

              T2 = vdes - vds;
              here->BSIM4v7gdtotd = T2 * dgdtot_dvd;
              here->BSIM4v7gdtotg = T2 * dgdtot_dvg;
              here->BSIM4v7gdtots = T2 * dgdtot_dvs;
              here->BSIM4v7gdtotb = T2 * dgdtot_dvb;
          }
          else /* WDLiu: for bypass */
          {   here->BSIM4v7gstot = here->BSIM4v7gstotd = here->BSIM4v7gstotg = 0.0;
              here->BSIM4v7gstots = here->BSIM4v7gstotb = 0.0;
              here->BSIM4v7gdtot = here->BSIM4v7gdtotd = here->BSIM4v7gdtotg = 0.0;
              here->BSIM4v7gdtots = here->BSIM4v7gdtotb = 0.0;
          }

          /* GIDL/GISL Models */

          if(model->BSIM4v7mtrlMod == 0)
            T0 = 3.0 * toxe;
          else
            T0 = model->BSIM4v7epsrsub * toxe / epsrox;
          
          /* Calculate GIDL current */

          vgs_eff = here->BSIM4v7vgs_eff;
          dvgs_eff_dvg = here->BSIM4v7dvgs_eff_dvg;
          vgd_eff = here->BSIM4v7vgd_eff;
          dvgd_eff_dvg = here->BSIM4v7dvgd_eff_dvg;

          if (model->BSIM4v7gidlMod==0){
          
          if(model->BSIM4v7mtrlMod ==0)
            T1 = (vds - vgs_eff - pParam->BSIM4v7egidl ) / T0;
          else
            T1 = (vds - vgs_eff - pParam->BSIM4v7egidl + pParam->BSIM4v7vfbsd) / T0;
          
          if ((pParam->BSIM4v7agidl <= 0.0) || (pParam->BSIM4v7bgidl <= 0.0)
              || (T1 <= 0.0) || (pParam->BSIM4v7cgidl <= 0.0) || (vbd > 0.0))
              Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
          else {
              dT1_dVd = 1.0 / T0;
              dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
              T2 = pParam->BSIM4v7bgidl / T1;
              if (T2 < 100.0)
              {   Igidl = pParam->BSIM4v7agidl * pParam->BSIM4v7weffCJ * T1 * exp(-T2);
                  T3 = Igidl * (1.0 + T2) / T1;
                  Ggidld = T3 * dT1_dVd;
                  Ggidlg = T3 * dT1_dVg;
              }
              else
              {   Igidl = pParam->BSIM4v7agidl * pParam->BSIM4v7weffCJ * 3.720075976e-44;
                  Ggidld = Igidl * dT1_dVd;
                  Ggidlg = Igidl * dT1_dVg;
                  Igidl *= T1;
              }
        
              T4 = vbd * vbd;
              T5 = -vbd * T4;
              T6 = pParam->BSIM4v7cgidl + T5;
              T7 = T5 / T6;
              T8 = 3.0 * pParam->BSIM4v7cgidl * T4 / T6 / T6;
              Ggidld = Ggidld * T7 + Igidl * T8;
              Ggidlg = Ggidlg * T7;
              Ggidlb = -Igidl * T8;
              Igidl *= T7;
          }
          here->BSIM4v7Igidl = Igidl;
          here->BSIM4v7ggidld = Ggidld;
          here->BSIM4v7ggidlg = Ggidlg;
          here->BSIM4v7ggidlb = Ggidlb;
        /* Calculate GISL current  */
          
          if(model->BSIM4v7mtrlMod ==0)
          T1 = (-vds - vgd_eff - pParam->BSIM4v7egisl ) / T0;
          else
          T1 = (-vds - vgd_eff - pParam->BSIM4v7egisl + pParam->BSIM4v7vfbsd ) / T0;

          if ((pParam->BSIM4v7agisl <= 0.0) || (pParam->BSIM4v7bgisl <= 0.0)
              || (T1 <= 0.0) || (pParam->BSIM4v7cgisl <= 0.0) || (vbs > 0.0))
              Igisl = Ggisls = Ggislg = Ggislb = 0.0;
          else {
              dT1_dVd = 1.0 / T0;
              dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
              T2 = pParam->BSIM4v7bgisl / T1;
              if (T2 < 100.0) 
              {   Igisl = pParam->BSIM4v7agisl * pParam->BSIM4v7weffCJ * T1 * exp(-T2);
                  T3 = Igisl * (1.0 + T2) / T1;
                  Ggisls = T3 * dT1_dVd;
                  Ggislg = T3 * dT1_dVg;
              }
              else 
              {   Igisl = pParam->BSIM4v7agisl * pParam->BSIM4v7weffCJ * 3.720075976e-44;
                  Ggisls = Igisl * dT1_dVd;
                  Ggislg = Igisl * dT1_dVg;
                  Igisl *= T1;
              }
        
              T4 = vbs * vbs;
              T5 = -vbs * T4;
              T6 = pParam->BSIM4v7cgisl + T5;
              T7 = T5 / T6;
              T8 = 3.0 * pParam->BSIM4v7cgisl * T4 / T6 / T6;
              Ggisls = Ggisls * T7 + Igisl * T8;
              Ggislg = Ggislg * T7;
              Ggislb = -Igisl * T8;
              Igisl *= T7;
          }
          here->BSIM4v7Igisl = Igisl;
          here->BSIM4v7ggisls = Ggisls;
          here->BSIM4v7ggislg = Ggislg;
          here->BSIM4v7ggislb = Ggislb;
          }
          else{ 
          /* v4.7 New Gidl/GISL model */
          
                    /* GISL */
                    if (model->BSIM4v7mtrlMod==0) 
                       T1 = (-vds - pParam->BSIM4v7rgisl * vgd_eff - pParam->BSIM4v7egisl) / T0;     
                    else 
                       T1 = (-vds - pParam->BSIM4v7rgisl * vgd_eff - pParam->BSIM4v7egisl + pParam->BSIM4v7vfbsd) / T0; 
                    
                    if ((pParam->BSIM4v7agisl <= 0.0) ||
                            (pParam->BSIM4v7bgisl <= 0.0) || (T1 <= 0.0) ||
                            (pParam->BSIM4v7cgisl < 0.0)  )
                        Igisl = Ggisls = Ggislg = Ggislb = 0.0; 
                    else
                    {
                        dT1_dVd = 1 / T0;                      
                        dT1_dVg = - pParam->BSIM4v7rgisl * dT1_dVd * dvgd_eff_dvg;
                        T2 = pParam->BSIM4v7bgisl / T1;
                        if (T2 < EXPL_THRESHOLD)
                        {
                            Igisl = pParam->BSIM4v7weffCJ * pParam->BSIM4v7agisl * T1 * exp(-T2);
                            T3 = Igisl / T1 * (T2 + 1);
                            Ggisls = T3 * dT1_dVd;
                            Ggislg = T3 * dT1_dVg;
                        } 
                        else
                        {
                            T3 = pParam->BSIM4v7weffCJ * pParam->BSIM4v7agisl * MIN_EXPL;
                            Igisl = T3 * T1 ;
                            Ggisls  = T3 * dT1_dVd;
                            Ggislg  = T3 * dT1_dVg;
                        
                        }
                        T4 = vbs - pParam->BSIM4v7fgisl;
                        
                        if (T4==0)
                            T5 = EXPL_THRESHOLD;
                        else
                            T5 = pParam->BSIM4v7kgisl / T4;
                        if (T5<EXPL_THRESHOLD)
                        {T6 = exp(T5);
                            Ggislb = -Igisl * T6 * T5 / T4;
                        }
                        else
                        {T6 = MAX_EXPL;
                            Ggislb=0.0;
                        }
                        Ggisls*=T6;
                        Ggislg*=T6;
                        Igisl*=T6;

                    }
                    here->BSIM4v7Igisl = Igisl;
                      here->BSIM4v7ggisls = Ggisls;
                      here->BSIM4v7ggislg = Ggislg;
                       here->BSIM4v7ggislb = Ggislb;
                    /* End of GISL */
                    
                    /* GIDL */
                    if (model->BSIM4v7mtrlMod==0) 
                        T1 = (vds - pParam->BSIM4v7rgidl * vgs_eff - pParam->BSIM4v7egidl) /  T0;                                           
                    else 
                        T1 = (vds - pParam->BSIM4v7rgidl * vgs_eff - pParam->BSIM4v7egidl + pParam->BSIM4v7vfbsd) / T0;
                                       
                    
                
                    if ((pParam->BSIM4v7agidl <= 0.0) ||
                            (pParam->BSIM4v7bgidl <= 0.0) || (T1 <= 0.0) ||
                            (pParam->BSIM4v7cgidl < 0.0)  )
                        Igidl = Ggidld = Ggidlg = Ggidlb = 0.0; 
                    else
                    {
                        dT1_dVd = 1 / T0;
                        dT1_dVg = - pParam->BSIM4v7rgidl * dT1_dVd * dvgs_eff_dvg;
                        T2 = pParam->BSIM4v7bgidl / T1;
                        if (T2 < EXPL_THRESHOLD)
                        {
                            Igidl = pParam->BSIM4v7weffCJ * pParam->BSIM4v7agidl * T1 * exp(-T2);
                            T3 = Igidl / T1 * (T2 + 1);
                            Ggidld = T3 * dT1_dVd;
                            Ggidlg = T3 * dT1_dVg;
                            
                        } else
                        {
                            T3 = pParam->BSIM4v7weffCJ * pParam->BSIM4v7agidl * MIN_EXPL;
                            Igidl = T3 * T1 ;
                            Ggidld  = T3 * dT1_dVd;
                            Ggidlg  = T3 * dT1_dVg;
                        }
                        T4 = vbd - pParam->BSIM4v7fgidl;
                        if (T4==0)
                            T5 = EXPL_THRESHOLD;
                        else
                            T5 = pParam->BSIM4v7kgidl / T4;
                        if (T5<EXPL_THRESHOLD)
                        {T6 = exp(T5);
                            Ggidlb = -Igidl * T6 * T5 / T4;
                        }
                        else
                        {T6 = MAX_EXPL;
                            Ggidlb=0.0;
                        }
                        Ggidld *= T6;
                        Ggidlg *= T6;
                        Igidl *= T6;
                    }                                        
                here->BSIM4v7Igidl = Igidl;
                     here->BSIM4v7ggidld = Ggidld;
                  here->BSIM4v7ggidlg = Ggidlg;
                  here->BSIM4v7ggidlb = Ggidlb;
                /* End of New GIDL */
                }
                   /*End of Gidl*/
          


          /* Calculate gate tunneling current */
          if ((model->BSIM4v7igcMod != 0) || (model->BSIM4v7igbMod != 0))
          {   Vfb = here->BSIM4v7vfbzb;
              V3 = Vfb - Vgs_eff + Vbseff - DELTA_3;
              if (Vfb <= 0.0)
                  T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
              else
                  T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
              T1 = 0.5 * (1.0 + V3 / T0);
              Vfbeff = Vfb - 0.5 * (V3 + T0);
              dVfbeff_dVg = T1 * dVgs_eff_dVg;
              dVfbeff_dVb = -T1; /* WDLiu: -No surprise? No. -Good! */

              Voxacc = Vfb - Vfbeff;
              dVoxacc_dVg = -dVfbeff_dVg;
              dVoxacc_dVb = -dVfbeff_dVb;
              if (Voxacc < 0.0) /* WDLiu: Avoiding numerical instability. */
                  Voxacc = dVoxacc_dVg = dVoxacc_dVb = 0.0;

              T0 = 0.5 * pParam->BSIM4v7k1ox;
              T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
              if (pParam->BSIM4v7k1ox == 0.0)
                  Voxdepinv = dVoxdepinv_dVg = dVoxdepinv_dVd
                            = dVoxdepinv_dVb = 0.0;
              else if (T3 < 0.0)
              {   Voxdepinv = -T3;
                  dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg
                                 + dVgsteff_dVg;
                  dVoxdepinv_dVd = dVgsteff_dVd;
                  dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
              }
              else
              {   T1 = sqrt(T0 * T0 + T3);
                  T2 = T0 / T1;
                  Voxdepinv = pParam->BSIM4v7k1ox * (T1 - T0);
                  dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg
                                 - dVgsteff_dVg);
                  dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
                  dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
              }

              Voxdepinv += Vgsteff;
              dVoxdepinv_dVg += dVgsteff_dVg;
              dVoxdepinv_dVd += dVgsteff_dVd;
              dVoxdepinv_dVb += dVgsteff_dVb;
          }

          if(model->BSIM4v7tempMod < 2)
                  tmp = Vtm;
          else /* model->BSIM4v7tempMod = 2 , 3*/
                  tmp = Vtm0;
          if (model->BSIM4v7igcMod)
          {   T0 = tmp * pParam->BSIM4v7nigc;
              if(model->BSIM4v7igcMod == 1) {
                      VxNVt = (Vgs_eff - model->BSIM4v7type * here->BSIM4v7vth0) / T0;
                      if (VxNVt > EXP_THRESHOLD)
                      {   Vaux = Vgs_eff - model->BSIM4v7type * here->BSIM4v7vth0;
                          dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                      }
              } else if (model->BSIM4v7igcMod == 2) {
                VxNVt = (Vgs_eff - here->BSIM4v7von) / T0;
                if (VxNVt > EXP_THRESHOLD)
                {   Vaux = Vgs_eff - here->BSIM4v7von;
                    dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = -dVth_dVd;
                    dVaux_dVb = -dVth_dVb;
                }
              } 
              if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
              }
              else if ((VxNVt >= -EXP_THRESHOLD) && (VxNVt <= EXP_THRESHOLD))
              {   ExpVxNVt = exp(VxNVt);
                  Vaux = T0 * log(1.0 + ExpVxNVt);
                  dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                  if(model->BSIM4v7igcMod == 1) {
                        dVaux_dVd = 0.0;
                          dVaux_dVb = 0.0;
                  } else if (model->BSIM4v7igcMod == 2) {
                        dVaux_dVd = -dVgs_eff_dVg * dVth_dVd;
                        dVaux_dVb = -dVgs_eff_dVg * dVth_dVb;
                  }
                  dVaux_dVg *= dVgs_eff_dVg;
              }

              T2 = Vgs_eff * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
              dT2_dVd = Vgs_eff * dVaux_dVd;
              dT2_dVb = Vgs_eff * dVaux_dVb;

              T11 = pParam->BSIM4v7Aechvb;
              T12 = pParam->BSIM4v7Bechvb;
              T3 = pParam->BSIM4v7aigc * pParam->BSIM4v7cigc
                 - pParam->BSIM4v7bigc;
              T4 = pParam->BSIM4v7bigc * pParam->BSIM4v7cigc;
              T5 = T12 * (pParam->BSIM4v7aigc + T3 * Voxdepinv
                 - T4 * Voxdepinv * Voxdepinv);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else
              {   T6 = exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                  dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                  dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                  dT6_dVg *= dVoxdepinv_dVg;
              }

              Igc = T11 * T2 * T6;
              dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
              dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

              if (model->BSIM4v7pigcdGiven)
              {   Pigcd = pParam->BSIM4v7pigcd;
                  dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
              }
              else
              {  /* T11 = pParam->BSIM4v7Bechvb * toxe; v4.7 */
                  T11 = -pParam->BSIM4v7Bechvb;
                  T12 = Vgsteff + 1.0e-20;
                  T13 = T11 / T12 / T12;
                  T14 = -T13 / T12;
                  Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12); 
                  dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg
                              - 3.0 * Vdseff / T12));
                  dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
                  dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
              }

              T7 = -Pigcd * Vdseff; /* bugfix */
              dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
              dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
              dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
              dT7_dVg *= dVgsteff_dVg;
              dT7_dVb *= dVbseff_dVb;
              T8 = T7 * T7 + 2.0e-4;
              dT8_dVg = 2.0 * T7;
              dT8_dVd = dT8_dVg * dT7_dVd;
              dT8_dVb = dT8_dVg * dT7_dVb;
              dT8_dVg *= dT7_dVg;

              if (T7 > EXP_THRESHOLD)
              {   T9 = MAX_EXP;
                  dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
              }
              else if (T7 < -EXP_THRESHOLD)
              {   T9 = MIN_EXP;
                  dT9_dVg = dT9_dVd = dT9_dVb = 0.0;
              }
              else
              {   T9 = exp(T7);
                  dT9_dVg = T9 * dT7_dVg;
                  dT9_dVd = T9 * dT7_dVd;
                  dT9_dVb = T9 * dT7_dVb;
              }

              T0 = T8 * T8;
              T1 = T9 - 1.0 + 1.0e-4;
              T10 = (T1 - T7) / T8;
              dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
              dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
              dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

              Igcs = Igc * T10;
              dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
              dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
              dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

              T1 = T9 - 1.0 - 1.0e-4;
              T10 = (T7 * T9 - T1) / T8;
              dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg
                       - T10 * dT8_dVg) / T8;
              dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd
                       - T10 * dT8_dVd) / T8;
              dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb
                       - T10 * dT8_dVb) / T8;
              Igcd = Igc * T10;
              dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
              dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
              dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

              here->BSIM4v7Igcs = Igcs;
              here->BSIM4v7gIgcsg = dIgcs_dVg;
              here->BSIM4v7gIgcsd = dIgcs_dVd;
              here->BSIM4v7gIgcsb =  dIgcs_dVb * dVbseff_dVb;
              here->BSIM4v7Igcd = Igcd;
              here->BSIM4v7gIgcdg = dIgcd_dVg;
              here->BSIM4v7gIgcdd = dIgcd_dVd;
              here->BSIM4v7gIgcdb = dIgcd_dVb * dVbseff_dVb;

              T0 = vgs - (pParam->BSIM4v7vfbsd + pParam->BSIM4v7vfbsdoff);
              vgs_eff = sqrt(T0 * T0 + 1.0e-4);
              dvgs_eff_dvg = T0 / vgs_eff;

              T2 = vgs * vgs_eff;
              dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
              T11 = pParam->BSIM4v7AechvbEdgeS;
              T12 = pParam->BSIM4v7BechvbEdge;
              T3 = pParam->BSIM4v7aigs * pParam->BSIM4v7cigs
                 - pParam->BSIM4v7bigs;
              T4 = pParam->BSIM4v7bigs * pParam->BSIM4v7cigs;
              T5 = T12 * (pParam->BSIM4v7aigs + T3 * vgs_eff
                 - T4 * vgs_eff * vgs_eff);
              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = 0.0;
              }
              else
              {   T6 = exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff)
                          * dvgs_eff_dvg;
              }
              Igs = T11 * T2 * T6;
              dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgs_dVs = -dIgs_dVg;


              T0 = vgd - (pParam->BSIM4v7vfbsd + pParam->BSIM4v7vfbsdoff);
              vgd_eff = sqrt(T0 * T0 + 1.0e-4);
              dvgd_eff_dvg = T0 / vgd_eff;

              T2 = vgd * vgd_eff;
              dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
              T11 = pParam->BSIM4v7AechvbEdgeD;
              T3 = pParam->BSIM4v7aigd * pParam->BSIM4v7cigd
                 - pParam->BSIM4v7bigd;
              T4 = pParam->BSIM4v7bigd * pParam->BSIM4v7cigd;
              T5 = T12 * (pParam->BSIM4v7aigd + T3 * vgd_eff
                 - T4 * vgd_eff * vgd_eff);
              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = 0.0;
              }
              else
              {   T6 = exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff)
                          * dvgd_eff_dvg;
              }
              Igd = T11 * T2 * T6;
              dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgd_dVd = -dIgd_dVg;

              here->BSIM4v7Igs = Igs;
              here->BSIM4v7gIgsg = dIgs_dVg;
              here->BSIM4v7gIgss = dIgs_dVs;
              here->BSIM4v7Igd = Igd;
              here->BSIM4v7gIgdg = dIgd_dVg;
              here->BSIM4v7gIgdd = dIgd_dVd;
          }
          else
          {   here->BSIM4v7Igcs = here->BSIM4v7gIgcsg = here->BSIM4v7gIgcsd
                               = here->BSIM4v7gIgcsb = 0.0;
              here->BSIM4v7Igcd = here->BSIM4v7gIgcdg = here->BSIM4v7gIgcdd
                                    = here->BSIM4v7gIgcdb = 0.0;
              here->BSIM4v7Igs = here->BSIM4v7gIgsg = here->BSIM4v7gIgss = 0.0;
              here->BSIM4v7Igd = here->BSIM4v7gIgdg = here->BSIM4v7gIgdd = 0.0;
          }

          if (model->BSIM4v7igbMod)
          {   T0 = tmp * pParam->BSIM4v7nigbacc;
              T1 = -Vgs_eff + Vbseff + Vfb;
              VxNVt = T1 / T0;
              if (VxNVt > EXP_THRESHOLD)
              {   Vaux = T1;
                  dVaux_dVg = -dVgs_eff_dVg;
                  dVaux_dVb = 1.0;
              }
              else if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVb = 0.0;
              }
              else
              {   ExpVxNVt = exp(VxNVt);
                  Vaux = T0 * log(1.0 + ExpVxNVt);
                  dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt); 
                  dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
              }

              T2 = (Vgs_eff - Vbseff) * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
              dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

              T11 = 4.97232e-7 * pParam->BSIM4v7weff
                  * pParam->BSIM4v7leff * pParam->BSIM4v7ToxRatio;
              T12 = -7.45669e11 * toxe;
              T3 = pParam->BSIM4v7aigbacc * pParam->BSIM4v7cigbacc
                 - pParam->BSIM4v7bigbacc;
              T4 = pParam->BSIM4v7bigbacc * pParam->BSIM4v7cigbacc;
              T5 = T12 * (pParam->BSIM4v7aigbacc + T3 * Voxacc
                 - T4 * Voxacc * Voxacc);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVb = 0.0;
              }
              else
              {   T6 = exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
                  dT6_dVb = dT6_dVg * dVoxacc_dVb;
                  dT6_dVg *= dVoxacc_dVg;
              }

              Igbacc = T11 * T2 * T6;
              dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);


              T0 = tmp * pParam->BSIM4v7nigbinv;
              T1 = Voxdepinv - pParam->BSIM4v7eigbinv;
              VxNVt = T1 / T0;
              if (VxNVt > EXP_THRESHOLD)
              {   Vaux = T1;
                  dVaux_dVg = dVoxdepinv_dVg;
                  dVaux_dVd = dVoxdepinv_dVd;
                  dVaux_dVb = dVoxdepinv_dVb;
              }
              else if (VxNVt < -EXP_THRESHOLD)
              {   Vaux = T0 * log(1.0 + MIN_EXP);
                  dVaux_dVg = dVaux_dVd = dVaux_dVb = 0.0;
              }
              else
              {   ExpVxNVt = exp(VxNVt);
                  Vaux = T0 * log(1.0 + ExpVxNVt);
                  dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                  dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
                  dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
                  dVaux_dVg *= dVoxdepinv_dVg;
              }

              T2 = (Vgs_eff - Vbseff) * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
              dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
              dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

              T11 *= 0.75610;
              T12 *= 1.31724;
              T3 = pParam->BSIM4v7aigbinv * pParam->BSIM4v7cigbinv
                 - pParam->BSIM4v7bigbinv;
              T4 = pParam->BSIM4v7bigbinv * pParam->BSIM4v7cigbinv;
              T5 = T12 * (pParam->BSIM4v7aigbinv + T3 * Voxdepinv
                 - T4 * Voxdepinv * Voxdepinv);

              if (T5 > EXP_THRESHOLD)
              {   T6 = MAX_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else if (T5 < -EXP_THRESHOLD)
              {   T6 = MIN_EXP;
                  dT6_dVg = dT6_dVd = dT6_dVb = 0.0;
              }
              else
              {   T6 = exp(T5);
                  dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                  dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                  dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                  dT6_dVg *= dVoxdepinv_dVg;
              }

              Igbinv = T11 * T2 * T6;
              dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
              dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
              dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

              here->BSIM4v7Igb = Igbinv + Igbacc;
              here->BSIM4v7gIgbg = dIgbinv_dVg + dIgbacc_dVg;
              here->BSIM4v7gIgbd = dIgbinv_dVd;
              here->BSIM4v7gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
          }
          else
          {  here->BSIM4v7Igb = here->BSIM4v7gIgbg = here->BSIM4v7gIgbd
                            = here->BSIM4v7gIgbs = here->BSIM4v7gIgbb = 0.0;   
          } /* End of Gate current */

          if (here->BSIM4v7nf != 1.0)
          {   cdrain *= here->BSIM4v7nf;
              here->BSIM4v7gds *= here->BSIM4v7nf;
              here->BSIM4v7gm *= here->BSIM4v7nf;
              here->BSIM4v7gmbs *= here->BSIM4v7nf;
              here->BSIM4v7IdovVds *= here->BSIM4v7nf;
                   
              here->BSIM4v7gbbs *= here->BSIM4v7nf;
              here->BSIM4v7gbgs *= here->BSIM4v7nf;
              here->BSIM4v7gbds *= here->BSIM4v7nf;
              here->BSIM4v7csub *= here->BSIM4v7nf;

              here->BSIM4v7Igidl *= here->BSIM4v7nf;
              here->BSIM4v7ggidld *= here->BSIM4v7nf;
              here->BSIM4v7ggidlg *= here->BSIM4v7nf;
              here->BSIM4v7ggidlb *= here->BSIM4v7nf;
        
              here->BSIM4v7Igisl *= here->BSIM4v7nf;
              here->BSIM4v7ggisls *= here->BSIM4v7nf;
              here->BSIM4v7ggislg *= here->BSIM4v7nf;
              here->BSIM4v7ggislb *= here->BSIM4v7nf;

              here->BSIM4v7Igcs *= here->BSIM4v7nf;
              here->BSIM4v7gIgcsg *= here->BSIM4v7nf;
              here->BSIM4v7gIgcsd *= here->BSIM4v7nf;
              here->BSIM4v7gIgcsb *= here->BSIM4v7nf;
              here->BSIM4v7Igcd *= here->BSIM4v7nf;
              here->BSIM4v7gIgcdg *= here->BSIM4v7nf;
              here->BSIM4v7gIgcdd *= here->BSIM4v7nf;
              here->BSIM4v7gIgcdb *= here->BSIM4v7nf;

              here->BSIM4v7Igs *= here->BSIM4v7nf;
              here->BSIM4v7gIgsg *= here->BSIM4v7nf;
              here->BSIM4v7gIgss *= here->BSIM4v7nf;
              here->BSIM4v7Igd *= here->BSIM4v7nf;
              here->BSIM4v7gIgdg *= here->BSIM4v7nf;
              here->BSIM4v7gIgdd *= here->BSIM4v7nf;

              here->BSIM4v7Igb *= here->BSIM4v7nf;
              here->BSIM4v7gIgbg *= here->BSIM4v7nf;
              here->BSIM4v7gIgbd *= here->BSIM4v7nf;
              here->BSIM4v7gIgbb *= here->BSIM4v7nf;
          }

          here->BSIM4v7ggidls = -(here->BSIM4v7ggidld + here->BSIM4v7ggidlg
                            + here->BSIM4v7ggidlb);
          here->BSIM4v7ggisld = -(here->BSIM4v7ggisls + here->BSIM4v7ggislg
                            + here->BSIM4v7ggislb);
          here->BSIM4v7gIgbs = -(here->BSIM4v7gIgbg + here->BSIM4v7gIgbd
                           + here->BSIM4v7gIgbb);
          here->BSIM4v7gIgcss = -(here->BSIM4v7gIgcsg + here->BSIM4v7gIgcsd
                            + here->BSIM4v7gIgcsb);
          here->BSIM4v7gIgcds = -(here->BSIM4v7gIgcdg + here->BSIM4v7gIgcdd
                            + here->BSIM4v7gIgcdb);
          here->BSIM4v7cd = cdrain;


          /* Calculations for noise analysis */

          if (model->BSIM4v7tnoiMod == 0)
          {   Abulk = Abulk0 * pParam->BSIM4v7abulkCVfactor;
              Vdsat = Vgsteff / Abulk;
              T0 = Vdsat - Vds - DELTA_4;
              T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * Vdsat);
              if (T0 >= 0.0)
                  Vdseff = Vdsat - 0.5 * (T0 + T1);
              else
              {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                  T4 = 1.0 - T3;
                  T5 = Vdsat * T3 / (T1 - T0);
                  Vdseff = Vdsat * T4;
              }
              if (Vds == 0.0)
                  Vdseff = 0.0;

              T0 = Abulk * Vdseff;
              T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
              T2 = Vdseff / T1;
              T3 = T0 * T2;
              here->BSIM4v7qinv = Coxeff * pParam->BSIM4v7weffCV * here->BSIM4v7nf
                              * pParam->BSIM4v7leffCV
                              * (Vgsteff - 0.5 * T0 + Abulk * T3);
          }
          else if(model->BSIM4v7tnoiMod == 2)
          {
              here->BSIM4v7noiGd0 = here->BSIM4v7nf * beta * Vgsteff / (1.0 + gche * Rds);
          }

          /*
           *  BSIM4v7 C-V begins
           */

          if ((model->BSIM4v7xpart < 0) || (!ChargeComputationNeeded))
          {   qgate  = qdrn = qsrc = qbulk = 0.0;
              here->BSIM4v7cggb = here->BSIM4v7cgsb = here->BSIM4v7cgdb = 0.0;
              here->BSIM4v7cdgb = here->BSIM4v7cdsb = here->BSIM4v7cddb = 0.0;
              here->BSIM4v7cbgb = here->BSIM4v7cbsb = here->BSIM4v7cbdb = 0.0;
              here->BSIM4v7csgb = here->BSIM4v7cssb = here->BSIM4v7csdb = 0.0;
              here->BSIM4v7cgbb = here->BSIM4v7csbb = here->BSIM4v7cdbb = here->BSIM4v7cbbb = 0.0;
              here->BSIM4v7cqdb = here->BSIM4v7cqsb = here->BSIM4v7cqgb 
                              = here->BSIM4v7cqbb = 0.0;
              here->BSIM4v7gtau = 0.0;
              goto finished;
          }
          else if (model->BSIM4v7capMod == 0)
          {
              if (Vbseff < 0.0)
              {   VbseffCV = Vbs; /*4.6.2*/
                  dVbseffCV_dVb = 1.0;
              }
              else
              {   VbseffCV = pParam->BSIM4v7phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb; /*4.6.2*/
              }

              Vfb = pParam->BSIM4v7vfbcv;
              Vth = Vfb + pParam->BSIM4v7phi + pParam->BSIM4v7k1ox * sqrtPhis; 
              Vgst = Vgs_eff - Vth;
              dVth_dVb = pParam->BSIM4v7k1ox * dsqrtPhis_dVb *dVbseff_dVb; /*4.6.2*/
              dVgst_dVb = -dVth_dVb;
              dVgst_dVg = dVgs_eff_dVg; 

              CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV
                    * pParam->BSIM4v7leffCV * here->BSIM4v7nf;
              Arg1 = Vgs_eff - VbseffCV - Vfb;

              if (Arg1 <= 0.0)
              {   qgate = CoxWL * Arg1;
                  qbulk = -qgate;
                  qdrn = 0.0;

                  here->BSIM4v7cggb = CoxWL * dVgs_eff_dVg;
                  here->BSIM4v7cgdb = 0.0;
                  here->BSIM4v7cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

                  here->BSIM4v7cdgb = 0.0;
                  here->BSIM4v7cddb = 0.0;
                  here->BSIM4v7cdsb = 0.0;

                  here->BSIM4v7cbgb = -CoxWL * dVgs_eff_dVg;
                  here->BSIM4v7cbdb = 0.0;
                  here->BSIM4v7cbsb = -here->BSIM4v7cgsb;
              } /* Arg1 <= 0.0, end of accumulation */
              else if (Vgst <= 0.0)
              {   T1 = 0.5 * pParam->BSIM4v7k1ox;
                  T2 = sqrt(T1 * T1 + Arg1);
                  qgate = CoxWL * pParam->BSIM4v7k1ox * (T2 - T1);
                  qbulk = -qgate;
                  qdrn = 0.0;

                  T0 = CoxWL * T1 / T2;
                  here->BSIM4v7cggb = T0 * dVgs_eff_dVg;
                  here->BSIM4v7cgdb = 0.0;
                  here->BSIM4v7cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);
   
                  here->BSIM4v7cdgb = 0.0;
                  here->BSIM4v7cddb = 0.0;
                  here->BSIM4v7cdsb = 0.0;

                  here->BSIM4v7cbgb = -here->BSIM4v7cggb;
                  here->BSIM4v7cbdb = 0.0;
                  here->BSIM4v7cbsb = -here->BSIM4v7cgsb;
              } /* Vgst <= 0.0, end of depletion */
              else
              {   One_Third_CoxWL = CoxWL / 3.0;
                  Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                  AbulkCV = Abulk0 * pParam->BSIM4v7abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v7abulkCVfactor * dAbulk0_dVb*dVbseff_dVb;
                  
                  dVdsat_dVg = 1.0 / AbulkCV;  /*4.6.2*/
                          Vdsat = Vgst * dVdsat_dVg;
                  dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)* dVdsat_dVg; 

                  if (model->BSIM4v7xpart > 0.5)
                  {   /* 0/100 Charge partition model */
                      if (Vdsat <= Vds)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v7phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.0;

                          here->BSIM4v7cggb = One_Third_CoxWL * (3.0
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb + T2);
                          here->BSIM4v7cgdb = 0.0;
       
                          here->BSIM4v7cdgb = 0.0;
                          here->BSIM4v7cddb = 0.0;
                          here->BSIM4v7cdsb = 0.0;

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v7cbsb = -(here->BSIM4v7cbgb + T3);
                          here->BSIM4v7cbdb = 0.0;
                      }
                      else
                      {   /* linear region */
                          Alphaz = Vgst / Vdsat;
                          T1 = 2.0 * Vdsat - Vds;
                          T2 = Vds / (3.0 * T1);
                          T3 = T2 * Vds;
                          T9 = 0.25 * CoxWL;
                          T4 = T9 * Alphaz;
                          T7 = 2.0 * Vds - T1 - 3.0 * T3;
                          T8 = T3 - T1 - 2.0 * Vds;
                          qgate = CoxWL * (Vgs_eff - Vfb 
                                - pParam->BSIM4v7phi - 0.5 * (Vds - T3));
                          T10 = T4 * T8;
                          qdrn = T4 * T7;
                          qbulk = -(qgate + qdrn + T10);
  
                          T5 = T3 / T1;
                          here->BSIM4v7cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          T11 = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v7cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb + T11
                                          + here->BSIM4v7cgdb);
                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                          T7 = T9 * T7;
                          T8 = T9 * T8;
                          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                          here->BSIM4v7cdgb = (T7 * dAlphaz_dVg - T9
                                          * dVdsat_dVg) * dVgs_eff_dVg;
                          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          here->BSIM4v7cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                          here->BSIM4v7cdsb = -(here->BSIM4v7cdgb + T12
                                          + here->BSIM4v7cddb);

                          T9 = 2.0 * T4 * (1.0 + T5);
                          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
                              * dVgs_eff_dVg;
                          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          T12 = T4 * (2.0 * T2 + T5 - 1.0); 
                          T0 = -(T10 + T11 + T12);

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb
                                          + here->BSIM4v7cdgb + T10);
                          here->BSIM4v7cbdb = -(here->BSIM4v7cgdb 
                                          + here->BSIM4v7cddb + T12);
                          here->BSIM4v7cbsb = -(here->BSIM4v7cgsb
                                          + here->BSIM4v7cdsb + T0);
                      }
                  }
                  else if (model->BSIM4v7xpart < 0.5)
                  {   /* 40/60 Charge partition model */
                      if (Vds >= Vdsat)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v7phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.4 * T2;

                          here->BSIM4v7cggb = One_Third_CoxWL * (3.0 
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb + T2);
                          here->BSIM4v7cgdb = 0.0;
       
                          T3 = 0.4 * Two_Third_CoxWL;
                          here->BSIM4v7cdgb = -T3 * dVgs_eff_dVg;
                          here->BSIM4v7cddb = 0.0;
                          T4 = T3 * dVth_dVb;
                          here->BSIM4v7cdsb = -(T4 + here->BSIM4v7cdgb);

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb 
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v7cbsb = -(here->BSIM4v7cbgb + T3);
                          here->BSIM4v7cbdb = 0.0;
                      }
                      else
                      {   /* linear region  */
                          Alphaz = Vgst / Vdsat;
                          T1 = 2.0 * Vdsat - Vds;
                          T2 = Vds / (3.0 * T1);
                          T3 = T2 * Vds;
                          T9 = 0.25 * CoxWL;
                          T4 = T9 * Alphaz;
                          qgate = CoxWL * (Vgs_eff - Vfb - pParam->BSIM4v7phi
                                - 0.5 * (Vds - T3));

                          T5 = T3 / T1;
                          here->BSIM4v7cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v7cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb 
                                          + here->BSIM4v7cgdb + tmp);

                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                          T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds
                             + 1.2 * Vds * Vds;
                          T8 = T2 / T1;
                          T7 = Vds - T1 - T8 * T6;
                          qdrn = T4 * T7;
                          T7 *= T9;
                          tmp = T8 / T1;
                          tmp1 = T4 * (2.0 - 4.0 * tmp * T6
                               + T8 * (16.0 * Vdsat - 6.0 * Vds));

                          here->BSIM4v7cdgb = (T7 * dAlphaz_dVg - tmp1
                                          * dVdsat_dVg) * dVgs_eff_dVg;
                          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                          here->BSIM4v7cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                          * T1) + 2.0 * tmp) * T6 + T8
                                          * (6.0 * Vdsat - 2.4 * Vds));
                          here->BSIM4v7cdsb = -(here->BSIM4v7cdgb 
                                          + T10 + here->BSIM4v7cddb);

                          T7 = 2.0 * (T1 + T3);
                          qbulk = -(qgate - T4 * T7);
                          T7 *= T9;
                          T0 = 4.0 * T4 * (1.0 - T5);
                          T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg
                                  - here->BSIM4v7cdgb;  /*4.6.2*/
                          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5) 
                              - here->BSIM4v7cddb;
                          tmp = -(T10 + T11 + T12);

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb 
                                          + here->BSIM4v7cdgb + T12);
                          here->BSIM4v7cbdb = -(here->BSIM4v7cgdb
                                          + here->BSIM4v7cddb + T10);  
                          here->BSIM4v7cbsb = -(here->BSIM4v7cgsb
                                          + here->BSIM4v7cdsb + tmp);
                      }
                  }
                  else
                  {   /* 50/50 partitioning */
                      if (Vds >= Vdsat)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v7phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.5 * T2;

                          here->BSIM4v7cggb = One_Third_CoxWL * (3.0
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb + T2);
                          here->BSIM4v7cgdb = 0.0;
       
                          here->BSIM4v7cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                          here->BSIM4v7cddb = 0.0;
                          T4 = One_Third_CoxWL * dVth_dVb;
                          here->BSIM4v7cdsb = -(T4 + here->BSIM4v7cdgb);

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb 
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v7cbsb = -(here->BSIM4v7cbgb + T3);
                          here->BSIM4v7cbdb = 0.0;
                      }
                      else
                      {   /* linear region */
                          Alphaz = Vgst / Vdsat;
                          T1 = 2.0 * Vdsat - Vds;
                          T2 = Vds / (3.0 * T1);
                          T3 = T2 * Vds;
                          T9 = 0.25 * CoxWL;
                          T4 = T9 * Alphaz;
                          qgate = CoxWL * (Vgs_eff - Vfb - pParam->BSIM4v7phi
                                - 0.5 * (Vds - T3));

                          T5 = T3 / T1;
                          here->BSIM4v7cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v7cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v7cgsb = -(here->BSIM4v7cggb 
                                          + here->BSIM4v7cgdb + tmp);

                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                          T7 = T1 + T3;
                          qdrn = -T4 * T7;
                          qbulk = - (qgate + qdrn + qdrn);
                          T7 *= T9;
                          T0 = T4 * (2.0 * T5 - 2.0);

                          here->BSIM4v7cdgb = (T0 * dVdsat_dVg - T7
                                          * dAlphaz_dVg) * dVgs_eff_dVg;
                          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                          here->BSIM4v7cddb = T4 * (1.0 - 2.0 * T2 - T5);
                          here->BSIM4v7cdsb = -(here->BSIM4v7cdgb + T12
                                          + here->BSIM4v7cddb);

                          here->BSIM4v7cbgb = -(here->BSIM4v7cggb
                                          + 2.0 * here->BSIM4v7cdgb);
                          here->BSIM4v7cbdb = -(here->BSIM4v7cgdb
                                          + 2.0 * here->BSIM4v7cddb);
                          here->BSIM4v7cbsb = -(here->BSIM4v7cgsb
                                          + 2.0 * here->BSIM4v7cdsb);
                      } /* end of linear region */
                  } /* end of 50/50 partition */
              } /* end of inversion */
          } /* end of capMod=0 */ 
          else
          {   if (Vbseff < 0.0)
              {   VbseffCV = Vbseff;
                  dVbseffCV_dVb = 1.0;
              }
              else
              {   VbseffCV = pParam->BSIM4v7phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb;
              }

              CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV
                    * pParam->BSIM4v7leffCV * here->BSIM4v7nf;

              if(model->BSIM4v7cvchargeMod == 0)
                {
                  /* Seperate VgsteffCV with noff and voffcv */
                  noff = n * pParam->BSIM4v7noff;
                  dnoff_dVd = pParam->BSIM4v7noff * dn_dVd;
                  dnoff_dVb = pParam->BSIM4v7noff * dn_dVb;
                  T0 = Vtm * noff;
                  voffcv = pParam->BSIM4v7voffcv;
                  VgstNVt = (Vgst - voffcv) / T0;
                  
                  if (VgstNVt > EXP_THRESHOLD)
                    {   
                      Vgsteff = Vgst - voffcv;
                      dVgsteff_dVg = dVgs_eff_dVg;
                      dVgsteff_dVd = -dVth_dVd;
                      dVgsteff_dVb = -dVth_dVb;
                    }
                  else if (VgstNVt < -EXP_THRESHOLD)
                    {   
                      Vgsteff = T0 * log(1.0 + MIN_EXP);
                      dVgsteff_dVg = 0.0;
                      dVgsteff_dVd = Vgsteff / noff;
                      dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
                      dVgsteff_dVd *= dnoff_dVd;
                    }
                  else
                    {   
                      ExpVgst = exp(VgstNVt);
                      Vgsteff = T0 * log(1.0 + ExpVgst);
                      dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
                      dVgsteff_dVd = -dVgsteff_dVg * (dVth_dVd + (Vgst - voffcv)
                                   / noff * dnoff_dVd) + Vgsteff / noff * dnoff_dVd;
                      dVgsteff_dVb = -dVgsteff_dVg * (dVth_dVb + (Vgst - voffcv)
                                   / noff * dnoff_dVb) + Vgsteff / noff * dnoff_dVb;
                      dVgsteff_dVg *= dVgs_eff_dVg;
                    } 
                  /* End of VgsteffCV for cvchargeMod = 0 */
                }
              else
                {                
                  T0 = n * Vtm;
                  T1 = pParam->BSIM4v7mstarcv * Vgst;
                  T2 = T1 / T0;
                  if (T2 > EXP_THRESHOLD)
                    {   
                      T10 = T1;
                      dT10_dVg = pParam->BSIM4v7mstarcv * dVgs_eff_dVg;
                      dT10_dVd = -dVth_dVd * pParam->BSIM4v7mstarcv;
                      dT10_dVb = -dVth_dVb * pParam->BSIM4v7mstarcv;
                    }
                  else if (T2 < -EXP_THRESHOLD)
                    {   
                      T10 = Vtm * log(1.0 + MIN_EXP);
                      dT10_dVg = 0.0;
                      dT10_dVd = T10 * dn_dVd;
                      dT10_dVb = T10 * dn_dVb;
                      T10 *= n;
                    }
                  else
                    {   
                      ExpVgst = exp(T2);
                      T3 = Vtm * log(1.0 + ExpVgst);
                      T10 = n * T3;
                      dT10_dVg = pParam->BSIM4v7mstarcv * ExpVgst / (1.0 + ExpVgst);
                      dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
                      dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
                      dT10_dVg *= dVgs_eff_dVg;
                    }
                  
                  T1 = pParam->BSIM4v7voffcbncv - (1.0 - pParam->BSIM4v7mstarcv) * Vgst;
                  T2 = T1 / T0;
                  if (T2 < -EXP_THRESHOLD)
                    {   
                      T3 = model->BSIM4v7coxe * MIN_EXP / pParam->BSIM4v7cdep0;
                      T9 = pParam->BSIM4v7mstarcv + T3 * n;
                      dT9_dVg = 0.0;
                      dT9_dVd = dn_dVd * T3;
                      dT9_dVb = dn_dVb * T3;
                    }
                  else if (T2 > EXP_THRESHOLD)
                    {   
                      T3 = model->BSIM4v7coxe * MAX_EXP / pParam->BSIM4v7cdep0;
                      T9 = pParam->BSIM4v7mstarcv + T3 * n;
                      dT9_dVg = 0.0;
                      dT9_dVd = dn_dVd * T3;
                      dT9_dVb = dn_dVb * T3;
                    }
                  else
                    {   
                      ExpVgst = exp(T2);
                      T3 = model->BSIM4v7coxe / pParam->BSIM4v7cdep0;
                      T4 = T3 * ExpVgst;
                      T5 = T1 * T4 / T0;
                      T9 = pParam->BSIM4v7mstarcv + n * T4;
                      dT9_dVg = T3 * (pParam->BSIM4v7mstarcv - 1.0) * ExpVgst / Vtm;
                      dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
                      dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
                      dT9_dVg *= dVgs_eff_dVg;
                    }
                  
                  Vgsteff = T10 / T9;
                  T11 = T9 * T9;
                  dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
                  dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
                  dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
                  /* End of VgsteffCV for cvchargeMod = 1 */
                }
          

              if (model->BSIM4v7capMod == 1)
              {   Vfb = here->BSIM4v7vfbzb;
                  V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
                  if (Vfb <= 0.0)
                      T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                  else
                      T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);

                  T1 = 0.5 * (1.0 + V3 / T0);
                  Vfbeff = Vfb - 0.5 * (V3 + T0);
                  dVfbeff_dVg = T1 * dVgs_eff_dVg;
                  dVfbeff_dVb = -T1 * dVbseffCV_dVb;
                  Qac0 = CoxWL * (Vfbeff - Vfb);
                  dQac0_dVg = CoxWL * dVfbeff_dVg;
                  dQac0_dVb = CoxWL * dVfbeff_dVb;

                  T0 = 0.5 * pParam->BSIM4v7k1ox;
                  T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam->BSIM4v7k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
                  else if (T3 < 0.0)
                  {   T1 = T0 + T3 / pParam->BSIM4v7k1ox;
                      T2 = CoxWL;
                  }
                  else
                  {   T1 = sqrt(T0 * T0 + T3);
                      T2 = CoxWL * T0 / T1;
                  }

                  Qsub0 = CoxWL * pParam->BSIM4v7k1ox * (T1 - T0);

                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb 
                             + dVgsteff_dVb);

                  AbulkCV = Abulk0 * pParam->BSIM4v7abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v7abulkCVfactor * dAbulk0_dVb;
                  VdsatCV = Vgsteff / AbulkCV;

                   T0 = VdsatCV - Vds - DELTA_4;
                  dT0_dVg = 1.0 / AbulkCV;
                  dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV; 
                  T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                  dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                  dT1_dVd = -T0 / T1;
                  dT1_dVb = dT1_dVg * dT0_dVb;
                  dT1_dVg *= dT0_dVg;
                  if (T0 >= 0.0)
                  {   VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                      dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                      dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                      dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                  }
                  else
                  {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                      T4 = 1.0 - T3;
                      T5 = VdsatCV * T3 / (T1 - T0);
                      VdseffCV = VdsatCV * T4;
                      dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                      dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                      dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                  }

                  if (Vds == 0.0)
                  {  VdseffCV = 0.0;
                     dVdseffCV_dVg = 0.0;
                     dVdseffCV_dVb = 0.0;
                  }

                  T0 = AbulkCV * VdseffCV;
                  T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
                  T2 = VdseffCV / T1;
                  T3 = T0 * T2;
                  
                  T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
                  T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
                  T6 = 12.0 * T2 * T2 * Vgsteff;
                  
                  qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
                  Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                  Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
                  Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                    + Cgg1 * dVgsteff_dVb;
                  Cgg1 *= dVgsteff_dVg;
                  
                  T7 = 1.0 - AbulkCV;
                  qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
                  T4 = -T7 * (T4 - 1.0);
                  T5 = -T7 * T5;
                  T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
                  Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                  Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
                  Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                    + Cbg1 * dVgsteff_dVb;
                  Cbg1 *= dVgsteff_dVg;
                  
                  if (model->BSIM4v7xpart > 0.5)
                    {   /* 0/100 Charge petition model */
                      T1 = T1 + T1;
                      qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0
                                       - T0 * T0 / T1);
                      T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
                      T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
                      T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
                      T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
                      Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                      Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                      Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                        + Csg * dVgsteff_dVb;
                      Csg *= dVgsteff_dVg;
                    }
                  else if (model->BSIM4v7xpart < 0.5)
                    {   /* 40/60 Charge petition model */
                      T1 = T1 / 12.0;
                      T2 = 0.5 * CoxWL / (T1 * T1);
                      T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff
                                      * (Vgsteff - 4.0 * T0 / 3.0))
                        - 2.0 * T0 * T0 * T0 / 15.0;
                      qsrc = -T2 * T3;
                      T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0)
                        + 0.4 * T0 * T0;
                      T4 = -2.0 * qsrc / T1 - T2 * (Vgsteff * (3.0
                                                               * Vgsteff - 8.0 * T0 / 3.0)
                                                    + 2.0 * T0 * T0 / 3.0);
                      T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
                      T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
                      Csg = (T4 + T5 * dVdseffCV_dVg);
                      Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                      Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                        + Csg * dVgsteff_dVb;
                      Csg *= dVgsteff_dVg;
                    }
                  else
                    {   /* 50/50 Charge petition model */
                      qsrc = -0.5 * (qgate + qbulk);
                      Csg = -0.5 * (Cgg1 + Cbg1);
                      Csb = -0.5 * (Cgb1 + Cbb1); 
                      Csd = -0.5 * (Cgd1 + Cbd1); 
                    }
                  
                  qgate += Qac0 + Qsub0;
                  qbulk -= (Qac0 + Qsub0);
                  qdrn = -(qgate + qbulk + qsrc);
                  
                  Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
                  Cgd = dQsub0_dVd + Cgd1;
                  Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;
                  
                  Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                  Cbd = Cbd1 - dQsub0_dVd;
                  Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;
                  
                  Cgb *= dVbseff_dVb;
                  Cbb *= dVbseff_dVb;
                  Csb *= dVbseff_dVb;
                  
                  here->BSIM4v7cggb = Cgg;
                  here->BSIM4v7cgsb = -(Cgg + Cgd + Cgb);
                  here->BSIM4v7cgdb = Cgd;
                  here->BSIM4v7cdgb = -(Cgg + Cbg + Csg);
                  here->BSIM4v7cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                                     + Csg + Csd + Csb);
                  here->BSIM4v7cddb = -(Cgd + Cbd + Csd);
                  here->BSIM4v7cbgb = Cbg;
                  here->BSIM4v7cbsb = -(Cbg + Cbd + Cbb);
                  here->BSIM4v7cbdb = Cbd;
              }
              
              /* Charge-Thickness capMod (CTM) begins */
              else if (model->BSIM4v7capMod == 2)
              {   V3 = here->BSIM4v7vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                  if (here->BSIM4v7vfbzb <= 0.0)
                      T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * here->BSIM4v7vfbzb);
                  else
                      T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * here->BSIM4v7vfbzb);

                  T1 = 0.5 * (1.0 + V3 / T0);
                  Vfbeff = here->BSIM4v7vfbzb - 0.5 * (V3 + T0);
                  dVfbeff_dVg = T1 * dVgs_eff_dVg;
                  dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                  Cox = model->BSIM4v7coxp;
                  Tox = 1.0e8 * model->BSIM4v7toxp;
                  T0 = (Vgs_eff - VbseffCV - here->BSIM4v7vfbzb) / Tox;
                  dT0_dVg = dVgs_eff_dVg / Tox;
                  dT0_dVb = -dVbseffCV_dVb / Tox;

                  tmp = T0 * pParam->BSIM4v7acde;
                  if ((-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD))
                  {   Tcen = pParam->BSIM4v7ldeb * exp(tmp);
                      dTcen_dVg = pParam->BSIM4v7acde * Tcen;
                      dTcen_dVb = dTcen_dVg * dT0_dVb;
                      dTcen_dVg *= dT0_dVg;
                  }
                  else if (tmp <= -EXP_THRESHOLD)
                  {   Tcen = pParam->BSIM4v7ldeb * MIN_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }
                  else
                  {   Tcen = pParam->BSIM4v7ldeb * MAX_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }

                  LINK = 1.0e-3 * model->BSIM4v7toxp;
                  V3 = pParam->BSIM4v7ldeb - Tcen - LINK;
                  V4 = sqrt(V3 * V3 + 4.0 * LINK * pParam->BSIM4v7ldeb);
                  Tcen = pParam->BSIM4v7ldeb - 0.5 * (V3 + V4);
                  T1 = 0.5 * (1.0 + V3 / V4);
                  dTcen_dVg *= T1;
                  dTcen_dVb *= T1;

                  Ccen = epssub / Tcen;
                  T2 = Cox / (Cox + Ccen);
                  Coxeff = T2 * Ccen;
                  T3 = -Ccen / Tcen;
                  dCoxeff_dVg = T2 * T2 * T3;
                  dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                  dCoxeff_dVg *= dTcen_dVg;
                  CoxWLcen = CoxWL * Coxeff / model->BSIM4v7coxe;

                  Qac0 = CoxWLcen * (Vfbeff - here->BSIM4v7vfbzb);
                  QovCox = Qac0 / Coxeff;
                  dQac0_dVg = CoxWLcen * dVfbeff_dVg
                            + QovCox * dCoxeff_dVg;
                  dQac0_dVb = CoxWLcen * dVfbeff_dVb 
                            + QovCox * dCoxeff_dVb;

                  T0 = 0.5 * pParam->BSIM4v7k1ox;
                  T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam->BSIM4v7k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
                  else if (T3 < 0.0)
                  {   T1 = T0 + T3 / pParam->BSIM4v7k1ox;
                      T2 = CoxWLcen;
                  }
                  else
                  {   T1 = sqrt(T0 * T0 + T3);
                      T2 = CoxWLcen * T0 / T1;
                  }

                  Qsub0 = CoxWLcen * pParam->BSIM4v7k1ox * (T1 - T0);
                  QovCox = Qsub0 / Coxeff;
                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                             + QovCox * dCoxeff_dVg;
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                             + QovCox * dCoxeff_dVb;

                  /* Gate-bias dependent delta Phis begins */
                  if (pParam->BSIM4v7k1ox <= 0.0)
                  {   Denomi = 0.25 * pParam->BSIM4v7moin * Vtm;
                      T0 = 0.5 * pParam->BSIM4v7sqrtPhi;
                  }
                  else
                  {   Denomi = pParam->BSIM4v7moin * Vtm 
                             * pParam->BSIM4v7k1ox * pParam->BSIM4v7k1ox;
                      T0 = pParam->BSIM4v7k1ox * pParam->BSIM4v7sqrtPhi;
                  }
                  T1 = 2.0 * T0 + Vgsteff;

                  DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
                  dDeltaPhi_dVg = 2.0 * Vtm * (T1 -T0) / (Denomi + T1 * Vgsteff);
                  /* End of delta Phis */

                  /* VgDP = Vgsteff - DeltaPhi */
                  T0 = Vgsteff - DeltaPhi - 0.001;
                  dT0_dVg = 1.0 - dDeltaPhi_dVg;
                  T1 = sqrt(T0 * T0 + Vgsteff * 0.004);
                  VgDP = 0.5 * (T0 + T1);
                  dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);                  
                  
                  Tox += Tox; /* WDLiu: Tcen reevaluated below due to different Vgsteff */
                  T0 = (Vgsteff + here->BSIM4v7vtfbphi2) / Tox;
                  tmp = exp(model->BSIM4v7bdos * 0.7 * log(T0));
                  T1 = 1.0 + tmp;
                  T2 = model->BSIM4v7bdos * 0.7 * tmp / (T0 * Tox);
                  Tcen = model->BSIM4v7ados * 1.9e-9 / T1;
                  dTcen_dVg = -Tcen * T2 / T1;
                  dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
                  dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
                  dTcen_dVg *= dVgsteff_dVg;

                  Ccen = epssub / Tcen;
                  T0 = Cox / (Cox + Ccen);
                  Coxeff = T0 * Ccen;
                  T1 = -Ccen / Tcen;
                  dCoxeff_dVg = T0 * T0 * T1;
                  dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
                  dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                  dCoxeff_dVg *= dTcen_dVg;
                  CoxWLcen = CoxWL * Coxeff / model->BSIM4v7coxe;

                  AbulkCV = Abulk0 * pParam->BSIM4v7abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v7abulkCVfactor * dAbulk0_dVb;
                  VdsatCV = VgDP / AbulkCV;

                  T0 = VdsatCV - Vds - DELTA_4;
                  dT0_dVg = dVgDP_dVg / AbulkCV;
                  dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                  T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                  dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                  dT1_dVd = -T0 / T1;
                  dT1_dVb = dT1_dVg * dT0_dVb;
                  dT1_dVg *= dT0_dVg;
                  if (T0 >= 0.0)
                  {   VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                      dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                      dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                      dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                  }
                  else
                  {   T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                      T4 = 1.0 - T3;
                      T5 = VdsatCV * T3 / (T1 - T0);
                      VdseffCV = VdsatCV * T4;
                      dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                      dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                      dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                  }

                  if (Vds == 0.0)
                  {  VdseffCV = 0.0;
                     dVdseffCV_dVg = 0.0;
                     dVdseffCV_dVb = 0.0;
                  }

                  T0 = AbulkCV * VdseffCV;
                  T1 = VgDP;
                  T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
                  T3 = T0 / T2;
                  T4 = 1.0 - 12.0 * T3 * T3;
                  T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
                  T6 = T5 * VdseffCV / AbulkCV;

                  qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
                  QovCox = qgate / Coxeff;
                  Cgg1 = CoxWLcen * (T4 * dVgDP_dVg 
                       + T5 * dVdseffCV_dVg);
                  Cgd1 = CoxWLcen * T5 * dVdseffCV_dVd + Cgg1 
                       * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                  Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb) 
                       + Cgg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                  Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;


                  T7 = 1.0 - AbulkCV;
                  T8 = T2 * T2;
                  T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
                  T10 = T9 * dVgDP_dVg;
                  T11 = -T7 * T5 / AbulkCV;
                  T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

                  qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
                  QovCox = qbulk / Coxeff;
                  Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
                  Cbd1 = CoxWLcen * T11 * dVdseffCV_dVd + Cbg1
                       * dVgsteff_dVd + QovCox * dCoxeff_dVd; 
                  Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
                       + Cbg1 * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                  Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                  if (model->BSIM4v7xpart > 0.5)
                  {   /* 0/100 partition */
                      qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 
                           - 0.5 * T0 * T0 / T2);
                      QovCox = qsrc / Coxeff;
                      T2 += T2;
                      T3 = T2 * T2;
                      T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
                      T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
                      T5 = T7 * AbulkCV;
                      T6 = T7 * VdseffCV;

                      Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
                      Csd = CoxWLcen * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd
                          + QovCox * dCoxeff_dVd;
                      Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                          + Csg * dVgsteff_dVb + QovCox * dCoxeff_dVb;
                      Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                  }
                  else if (model->BSIM4v7xpart < 0.5)
                  {   /* 40/60 partition */
                      T2 = T2 / 12.0;
                      T3 = 0.5 * CoxWLcen / (T2 * T2);
                      T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 
                         * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
                      qsrc = -T3 * T4;
                      QovCox = qsrc / Coxeff;
                      T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
                      T5 = -2.0 * qsrc / T2 - T3 * (T1 * (3.0 * T1 - 8.0 
                         * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
                      T6 = AbulkCV * (qsrc / T2 + T3 * T8);
                      T7 = T6 * VdseffCV / AbulkCV; 

                      Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg; 
                      Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd 
                          + QovCox * dCoxeff_dVd;
                      Csb = Csg * dVgsteff_dVb + T6 * dVdseffCV_dVb 
                          + T7 * dAbulkCV_dVb + QovCox * dCoxeff_dVb; 
                      Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                  }
                  else
                  {   /* 50/50 partition */
                      qsrc = -0.5 * qgate;
                      Csg = -0.5 * Cgg1;
                      Csd = -0.5 * Cgd1; 
                      Csb = -0.5 * Cgb1; 
                  }

                  qgate += Qac0 + Qsub0 - qbulk;
                  qbulk -= (Qac0 + Qsub0);
                  qdrn = -(qgate + qbulk + qsrc);

                  Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                  Cbd = Cbd1 - dQsub0_dVd;
                  Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                  Cgg = Cgg1 - Cbg;
                  Cgd = Cgd1 - Cbd;
                  Cgb = Cgb1 - Cbb;

                  Cgb *= dVbseff_dVb;
                  Cbb *= dVbseff_dVb;
                  Csb *= dVbseff_dVb;

                  here->BSIM4v7cggb = Cgg;
                  here->BSIM4v7cgsb = -(Cgg + Cgd + Cgb);
                  here->BSIM4v7cgdb = Cgd;
                  here->BSIM4v7cdgb = -(Cgg + Cbg + Csg);
                  here->BSIM4v7cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                                  + Csg + Csd + Csb);
                  here->BSIM4v7cddb = -(Cgd + Cbd + Csd);
                  here->BSIM4v7cbgb = Cbg;
                  here->BSIM4v7cbsb = -(Cbg + Cbd + Cbb);
                  here->BSIM4v7cbdb = Cbd;
              }  /* End of CTM */
          }

          here->BSIM4v7csgb = - here->BSIM4v7cggb - here->BSIM4v7cdgb - here->BSIM4v7cbgb;
          here->BSIM4v7csdb = - here->BSIM4v7cgdb - here->BSIM4v7cddb - here->BSIM4v7cbdb;
          here->BSIM4v7cssb = - here->BSIM4v7cgsb - here->BSIM4v7cdsb - here->BSIM4v7cbsb;
          here->BSIM4v7cgbb = - here->BSIM4v7cgdb - here->BSIM4v7cggb - here->BSIM4v7cgsb;
          here->BSIM4v7cdbb = - here->BSIM4v7cddb - here->BSIM4v7cdgb - here->BSIM4v7cdsb;
          here->BSIM4v7cbbb = - here->BSIM4v7cbgb - here->BSIM4v7cbdb - here->BSIM4v7cbsb;
          here->BSIM4v7csbb = - here->BSIM4v7cgbb - here->BSIM4v7cdbb - here->BSIM4v7cbbb;
          here->BSIM4v7qgate = qgate;
          here->BSIM4v7qbulk = qbulk;
          here->BSIM4v7qdrn = qdrn;
          here->BSIM4v7qsrc = -(qgate + qbulk + qdrn);

          /* NQS begins */
          if ((here->BSIM4v7trnqsMod) || (here->BSIM4v7acnqsMod))
          {   here->BSIM4v7qchqs = qcheq = -(qbulk + qgate);
              here->BSIM4v7cqgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb);
              here->BSIM4v7cqdb = -(here->BSIM4v7cgdb + here->BSIM4v7cbdb);
              here->BSIM4v7cqsb = -(here->BSIM4v7cgsb + here->BSIM4v7cbsb);
              here->BSIM4v7cqbb = -(here->BSIM4v7cqgb + here->BSIM4v7cqdb
                              + here->BSIM4v7cqsb);

              CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV * here->BSIM4v7nf
                    * pParam->BSIM4v7leffCV;
              T1 = here->BSIM4v7gcrg / CoxWL; /* 1 / tau */
              here->BSIM4v7gtau = T1 * ScalingFactor;

              if (here->BSIM4v7acnqsMod)
                  here->BSIM4v7taunet = 1.0 / T1;

              *(ckt->CKTstate0 + here->BSIM4v7qcheq) = qcheq;
              if (ckt->CKTmode & MODEINITTRAN)
                  *(ckt->CKTstate1 + here->BSIM4v7qcheq) =
                                   *(ckt->CKTstate0 + here->BSIM4v7qcheq);
              if (here->BSIM4v7trnqsMod)
              {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qcheq);
                  if (error)
                      return(error);
              }
          }


finished: 

          /* Calculate junction C-V */
          if (ChargeComputationNeeded)
          {   czbd = model->BSIM4v7DunitAreaTempJctCap * here->BSIM4v7Adeff; /* bug fix */
              czbs = model->BSIM4v7SunitAreaTempJctCap * here->BSIM4v7Aseff;
              czbdsw = model->BSIM4v7DunitLengthSidewallTempJctCap * here->BSIM4v7Pdeff;
              czbdswg = model->BSIM4v7DunitLengthGateSidewallTempJctCap
                      * pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
              czbssw = model->BSIM4v7SunitLengthSidewallTempJctCap * here->BSIM4v7Pseff;
              czbsswg = model->BSIM4v7SunitLengthGateSidewallTempJctCap
                      * pParam->BSIM4v7weffCJ * here->BSIM4v7nf;

              MJS = model->BSIM4v7SbulkJctBotGradingCoeff;
              MJSWS = model->BSIM4v7SbulkJctSideGradingCoeff;
              MJSWGS = model->BSIM4v7SbulkJctGateSideGradingCoeff;

              MJD = model->BSIM4v7DbulkJctBotGradingCoeff;
              MJSWD = model->BSIM4v7DbulkJctSideGradingCoeff;
              MJSWGD = model->BSIM4v7DbulkJctGateSideGradingCoeff;

              /* Source Bulk Junction */
              if (vbs_jct == 0.0)
              {   *(ckt->CKTstate0 + here->BSIM4v7qbs) = 0.0;
                  here->BSIM4v7capbs = czbs + czbssw + czbsswg;
              }
              else if (vbs_jct < 0.0)
              {   if (czbs > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v7PhiBS;
                      if (MJS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbs) = model->BSIM4v7PhiBS * czbs 
                                       * (1.0 - arg * sarg) / (1.0 - MJS);
                      here->BSIM4v7capbs = czbs * sarg;
                  }
                  else
                  {   *(ckt->CKTstate0 + here->BSIM4v7qbs) = 0.0;
                      here->BSIM4v7capbs = 0.0;
                  }
                  if (czbssw > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v7PhiBSWS;
                      if (MJSWS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbs) += model->BSIM4v7PhiBSWS * czbssw
                                       * (1.0 - arg * sarg) / (1.0 - MJSWS);
                      here->BSIM4v7capbs += czbssw * sarg;
                  }
                  if (czbsswg > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v7PhiBSWGS;
                      if (MJSWGS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWGS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbs) += model->BSIM4v7PhiBSWGS * czbsswg
                                       * (1.0 - arg * sarg) / (1.0 - MJSWGS);
                      here->BSIM4v7capbs += czbsswg * sarg;
                  }

              }
              else
              {   T0 = czbs + czbssw + czbsswg;
                  T1 = vbs_jct * (czbs * MJS / model->BSIM4v7PhiBS + czbssw * MJSWS 
                     / model->BSIM4v7PhiBSWS + czbsswg * MJSWGS / model->BSIM4v7PhiBSWGS);    
                  *(ckt->CKTstate0 + here->BSIM4v7qbs) = vbs_jct * (T0 + 0.5 * T1);
                  here->BSIM4v7capbs = T0 + T1;
              }

              /* Drain Bulk Junction */
              if (vbd_jct == 0.0)
              {   *(ckt->CKTstate0 + here->BSIM4v7qbd) = 0.0;
                  here->BSIM4v7capbd = czbd + czbdsw + czbdswg;
              }
              else if (vbd_jct < 0.0)
              {   if (czbd > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v7PhiBD;
                      if (MJD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbd) = model->BSIM4v7PhiBD* czbd 
                                       * (1.0 - arg * sarg) / (1.0 - MJD);
                      here->BSIM4v7capbd = czbd * sarg;
                  }
                  else
                  {   *(ckt->CKTstate0 + here->BSIM4v7qbd) = 0.0;
                      here->BSIM4v7capbd = 0.0;
                  }
                  if (czbdsw > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v7PhiBSWD;
                      if (MJSWD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbd) += model->BSIM4v7PhiBSWD * czbdsw 
                                       * (1.0 - arg * sarg) / (1.0 - MJSWD);
                      here->BSIM4v7capbd += czbdsw * sarg;
                  }
                  if (czbdswg > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v7PhiBSWGD;
                      if (MJSWGD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWGD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v7qbd) += model->BSIM4v7PhiBSWGD * czbdswg
                                       * (1.0 - arg * sarg) / (1.0 - MJSWGD);
                      here->BSIM4v7capbd += czbdswg * sarg;
                  }
              }
              else
              {   T0 = czbd + czbdsw + czbdswg;
                  T1 = vbd_jct * (czbd * MJD / model->BSIM4v7PhiBD + czbdsw * MJSWD
                     / model->BSIM4v7PhiBSWD + czbdswg * MJSWGD / model->BSIM4v7PhiBSWGD);
                  *(ckt->CKTstate0 + here->BSIM4v7qbd) = vbd_jct * (T0 + 0.5 * T1);
                  here->BSIM4v7capbd = T0 + T1; 
              }
          }


          /*
           *  check convergence
           */

          if ((here->BSIM4v7off == 0) || (!(ckt->CKTmode & MODEINITFIX)))
          {   if (Check == 1)
              {   ckt->CKTnoncon++;
#ifndef NEWCONV
              } 
              else
              {   if (here->BSIM4v7mode >= 0)
                  {   Idtot = here->BSIM4v7cd + here->BSIM4v7csub
                            + here->BSIM4v7Igidl - here->BSIM4v7cbd;
                  }
                  else
                  {   Idtot = here->BSIM4v7cd + here->BSIM4v7cbd - here->BSIM4v7Igidl; /* bugfix */
                  }
                  tol0 = ckt->CKTreltol * MAX(fabs(cdhat), fabs(Idtot))
                       + ckt->CKTabstol;
                  tol1 = ckt->CKTreltol * MAX(fabs(cseshat), fabs(Isestot))
                       + ckt->CKTabstol;
                  tol2 = ckt->CKTreltol * MAX(fabs(cdedhat), fabs(Idedtot))
                       + ckt->CKTabstol;
                  tol3 = ckt->CKTreltol * MAX(fabs(cgshat), fabs(Igstot))
                       + ckt->CKTabstol;
                  tol4 = ckt->CKTreltol * MAX(fabs(cgdhat), fabs(Igdtot))
                       + ckt->CKTabstol;
                  tol5 = ckt->CKTreltol * MAX(fabs(cgbhat), fabs(Igbtot))
                       + ckt->CKTabstol;
                  if ((fabs(cdhat - Idtot) >= tol0) || (fabs(cseshat - Isestot) >= tol1)
                      || (fabs(cdedhat - Idedtot) >= tol2))
                  {   ckt->CKTnoncon++;
                  }
                  else if ((fabs(cgshat - Igstot) >= tol3) || (fabs(cgdhat - Igdtot) >= tol4)
                      || (fabs(cgbhat - Igbtot) >= tol5))
                  {   ckt->CKTnoncon++;
                  }
                  else
                  {   Ibtot = here->BSIM4v7cbs + here->BSIM4v7cbd
                            - here->BSIM4v7Igidl - here->BSIM4v7Igisl - here->BSIM4v7csub;
                      tol6 = ckt->CKTreltol * MAX(fabs(cbhat), fabs(Ibtot))
                          + ckt->CKTabstol;
                      if (fabs(cbhat - Ibtot) > tol6)
                      {   ckt->CKTnoncon++;
                      }
                  }
#endif /* NEWCONV */
              }
          }
          *(ckt->CKTstate0 + here->BSIM4v7vds) = vds;
          *(ckt->CKTstate0 + here->BSIM4v7vgs) = vgs;
          *(ckt->CKTstate0 + here->BSIM4v7vbs) = vbs;
          *(ckt->CKTstate0 + here->BSIM4v7vbd) = vbd;
          *(ckt->CKTstate0 + here->BSIM4v7vges) = vges;
          *(ckt->CKTstate0 + here->BSIM4v7vgms) = vgms;
          *(ckt->CKTstate0 + here->BSIM4v7vdbs) = vdbs;
          *(ckt->CKTstate0 + here->BSIM4v7vdbd) = vdbd;
          *(ckt->CKTstate0 + here->BSIM4v7vsbs) = vsbs;
          *(ckt->CKTstate0 + here->BSIM4v7vses) = vses;
          *(ckt->CKTstate0 + here->BSIM4v7vdes) = vdes;
          *(ckt->CKTstate0 + here->BSIM4v7qdef) = qdef;


          if (!ChargeComputationNeeded)
              goto line850; 

          if (here->BSIM4v7rgateMod == 3) 
          {   
                  vgdx = vgmd; 
                  vgsx = vgms;
          }  
          else  /* For rgateMod == 0, 1 and 2 */
          {
                  vgdx = vgd;
                  vgsx = vgs;
          }
          if (model->BSIM4v7capMod == 0) 
          {  
                  cgdo = pParam->BSIM4v7cgdo; 
                  qgdo = pParam->BSIM4v7cgdo * vgdx;
                  cgso = pParam->BSIM4v7cgso;
                  qgso = pParam->BSIM4v7cgso * vgsx;
          }
          else /* For both capMod == 1 and 2 */
          {   T0 = vgdx + DELTA_1;
              T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
              T2 = 0.5 * (T0 - T1);

              T3 = pParam->BSIM4v7weffCV * pParam->BSIM4v7cgdl;
              T4 = sqrt(1.0 - 4.0 * T2 / pParam->BSIM4v7ckappad); 
              cgdo = pParam->BSIM4v7cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                   * (0.5 - 0.5 * T0 / T1);
              qgdo = (pParam->BSIM4v7cgdo + T3) * vgdx - T3 * (T2
                   + 0.5 * pParam->BSIM4v7ckappad * (T4 - 1.0));

              T0 = vgsx + DELTA_1;
              T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
              T2 = 0.5 * (T0 - T1);
              T3 = pParam->BSIM4v7weffCV * pParam->BSIM4v7cgsl;
              T4 = sqrt(1.0 - 4.0 * T2 / pParam->BSIM4v7ckappas);
              cgso = pParam->BSIM4v7cgso + T3 - T3 * (1.0 - 1.0 / T4)
                   * (0.5 - 0.5 * T0 / T1);
              qgso = (pParam->BSIM4v7cgso + T3) * vgsx - T3 * (T2
                   + 0.5 * pParam->BSIM4v7ckappas * (T4 - 1.0));
          }

          if (here->BSIM4v7nf != 1.0)
          {   cgdo *= here->BSIM4v7nf;
              cgso *= here->BSIM4v7nf;
              qgdo *= here->BSIM4v7nf;
              qgso *= here->BSIM4v7nf;
          }        
          here->BSIM4v7cgdo = cgdo;
          here->BSIM4v7qgdo = qgdo;
          here->BSIM4v7cgso = cgso;
          here->BSIM4v7qgso = qgso;

#ifndef NOBYPASS
line755:
#endif
          ag0 = ckt->CKTag[0];
          if (here->BSIM4v7mode > 0)
          {   if (here->BSIM4v7trnqsMod == 0)
              {   qdrn -= qgdo;
                  if (here->BSIM4v7rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcggb = here->BSIM4v7cggb * ag0;
                      gcgdb = here->BSIM4v7cgdb * ag0;
                      gcgsb = here->BSIM4v7cgsb * ag0;   
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = here->BSIM4v7cdgb * ag0;
                      gcsgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                            + here->BSIM4v7cdgb) * ag0;
                      gcbgb = here->BSIM4v7cbgb * ag0;

                      qgmb = pParam->BSIM4v7cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qbulk -= qgmb;
                      qsrc = -(qgate + qgmid + qbulk + qdrn);
                  }
                    else
                  {   gcggb = (here->BSIM4v7cggb + cgdo + cgso
                            + pParam->BSIM4v7cgbo ) * ag0;
                      gcgdb = (here->BSIM4v7cgdb - cgdo) * ag0;
                      gcgsb = (here->BSIM4v7cgsb - cgso) * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = (here->BSIM4v7cdgb - cgdo) * ag0;
                      gcsgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                            + here->BSIM4v7cdgb + cgso) * ag0;
                      gcbgb = (here->BSIM4v7cbgb - pParam->BSIM4v7cgbo) * ag0;

                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v7cgbo * vgb;
                      qgate += qgdo + qgso + qgb;
                      qbulk -= qgb;
                      qsrc = -(qgate + qbulk + qdrn);
                  }
                  gcddb = (here->BSIM4v7cddb + here->BSIM4v7capbd + cgdo) * ag0;
                  gcdsb = here->BSIM4v7cdsb * ag0;

                  gcsdb = -(here->BSIM4v7cgdb + here->BSIM4v7cbdb
                        + here->BSIM4v7cddb) * ag0;
                  gcssb = (here->BSIM4v7capbs + cgso - (here->BSIM4v7cgsb
                        + here->BSIM4v7cbsb + here->BSIM4v7cdsb)) * ag0;

                  if (!here->BSIM4v7rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                      gcbdb = (here->BSIM4v7cbdb - here->BSIM4v7capbd) * ag0;
                      gcbsb = (here->BSIM4v7cbsb - here->BSIM4v7capbs) * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb  = -(here->BSIM4v7cddb + here->BSIM4v7cdgb 
                             + here->BSIM4v7cdsb) * ag0;
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb)
                            + here->BSIM4v7capbs * ag0;
                      gcbdb = here->BSIM4v7cbdb * ag0;
                      gcbsb = here->BSIM4v7cbsb * ag0;

                      gcdbdb = -here->BSIM4v7capbd * ag0;
                      gcsbsb = -here->BSIM4v7capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);

                  ggtg = ggtd = ggtb = ggts = 0.0;
                  sxpart = 0.6;
                  dxpart = 0.4;
                  ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                  dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
              }
              else
              {   qcheq = here->BSIM4v7qchqs;
                  CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV * here->BSIM4v7nf
                        * pParam->BSIM4v7leffCV;
                  T0 = qdef * ScalingFactor / CoxWL;

                  ggtg = here->BSIM4v7gtg = T0 * here->BSIM4v7gcrgg;
                  ggtd = here->BSIM4v7gtd = T0 * here->BSIM4v7gcrgd;
                  ggts = here->BSIM4v7gts = T0 * here->BSIM4v7gcrgs;
                  ggtb = here->BSIM4v7gtb = T0 * here->BSIM4v7gcrgb;
                  gqdef = ScalingFactor * ag0;

                  gcqgb = here->BSIM4v7cqgb * ag0;
                  gcqdb = here->BSIM4v7cqdb * ag0;
                  gcqsb = here->BSIM4v7cqsb * ag0;
                  gcqbb = here->BSIM4v7cqbb * ag0;

                  if (fabs(qcheq) <= 1.0e-5 * CoxWL)
                  {   if (model->BSIM4v7xpart < 0.5)
                      {   dxpart = 0.4;
                      }
                      else if (model->BSIM4v7xpart > 0.5)
                      {   dxpart = 0.0;
                      }
                      else
                      {   dxpart = 0.5;
                      }
                      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                  = ddxpart_dVs = 0.0;
                  }
                  else
                  {   dxpart = qdrn / qcheq;
                      Cdd = here->BSIM4v7cddb;
                      Csd = -(here->BSIM4v7cgdb + here->BSIM4v7cddb
                          + here->BSIM4v7cbdb);
                      ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                      Cdg = here->BSIM4v7cdgb;
                      Csg = -(here->BSIM4v7cggb + here->BSIM4v7cdgb
                          + here->BSIM4v7cbgb);
                      ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                      Cds = here->BSIM4v7cdsb;
                      Css = -(here->BSIM4v7cgsb + here->BSIM4v7cdsb
                          + here->BSIM4v7cbsb);
                      ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                  }
                  sxpart = 1.0 - dxpart;
                  dsxpart_dVd = -ddxpart_dVd;
                  dsxpart_dVg = -ddxpart_dVg;
                  dsxpart_dVs = -ddxpart_dVs;
                  dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                  if (here->BSIM4v7rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcdgb = gcsgb = gcbgb = 0.0;
                      gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                      qgmb = pParam->BSIM4v7cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qgate = 0.0;
                      qbulk = -qgmb;
                      qdrn = -qgdo;
                      qsrc = -(qgmid + qbulk + qdrn);
                  }
                  else
                  {   gcggb = (cgdo + cgso + pParam->BSIM4v7cgbo ) * ag0;
                      gcgdb = -cgdo * ag0;
                      gcgsb = -cgso * ag0;
                      gcgbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgb = gcgdb;
                      gcsgb = gcgsb;
                      gcbgb = gcgbb;
                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v7cgbo * vgb;
                      qgate = qgdo + qgso + qgb;
                      qbulk = -qgb;
                      qdrn = -qgdo;
                      qsrc = -(qgate + qbulk + qdrn);
                  }

                  gcddb = (here->BSIM4v7capbd + cgdo) * ag0;
                  gcdsb = gcsdb = 0.0;
                  gcssb = (here->BSIM4v7capbs + cgso) * ag0;

                  if (!here->BSIM4v7rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                      gcsbb = -(gcsgb + gcssb + gcsgmb);
                      gcbdb = -here->BSIM4v7capbd * ag0;
                      gcbsb = -here->BSIM4v7capbs * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                      gcdbdb = -here->BSIM4v7capbd * ag0;
                      gcsbsb = -here->BSIM4v7capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
              }
          }
          else
          {   if (here->BSIM4v7trnqsMod == 0)
              {   qsrc = qdrn - qgso;
                  if (here->BSIM4v7rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                          gcgmsb = -cgso * ag0;
                          gcgmbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcggb = here->BSIM4v7cggb * ag0;
                      gcgdb = here->BSIM4v7cgsb * ag0;
                      gcgsb = here->BSIM4v7cgdb * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                            + here->BSIM4v7cdgb) * ag0;
                      gcsgb = here->BSIM4v7cdgb * ag0;
                      gcbgb = here->BSIM4v7cbgb * ag0;

                      qgmb = pParam->BSIM4v7cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qbulk -= qgmb;
                      qdrn = -(qgate + qgmid + qbulk + qsrc);
                  }
                  else
                  {   gcggb = (here->BSIM4v7cggb + cgdo + cgso
                            + pParam->BSIM4v7cgbo ) * ag0;
                      gcgdb = (here->BSIM4v7cgsb - cgdo) * ag0;
                      gcgsb = (here->BSIM4v7cgdb - cgso) * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                            + here->BSIM4v7cdgb + cgdo) * ag0;
                      gcsgb = (here->BSIM4v7cdgb - cgso) * ag0;
                      gcbgb = (here->BSIM4v7cbgb - pParam->BSIM4v7cgbo) * ag0;

                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v7cgbo * vgb;
                      qgate += qgdo + qgso + qgb;
                      qbulk -= qgb;
                      qdrn = -(qgate + qbulk + qsrc);
                   }
                  gcddb = (here->BSIM4v7capbd + cgdo - (here->BSIM4v7cgsb
                        + here->BSIM4v7cbsb + here->BSIM4v7cdsb)) * ag0;
                  gcdsb = -(here->BSIM4v7cgdb + here->BSIM4v7cbdb
                        + here->BSIM4v7cddb) * ag0;

                  gcsdb = here->BSIM4v7cdsb * ag0;
                  gcssb = (here->BSIM4v7cddb + here->BSIM4v7capbs + cgso) * ag0;

                  if (!here->BSIM4v7rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                      gcbdb = (here->BSIM4v7cbsb - here->BSIM4v7capbd) * ag0;
                      gcbsb = (here->BSIM4v7cbdb - here->BSIM4v7capbs) * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb)
                            + here->BSIM4v7capbd * ag0;
                      gcsbb = -(here->BSIM4v7cddb + here->BSIM4v7cdgb
                            + here->BSIM4v7cdsb) * ag0;
                      gcbdb = here->BSIM4v7cbsb * ag0;
                      gcbsb = here->BSIM4v7cbdb * ag0;
                      gcdbdb = -here->BSIM4v7capbd * ag0;
                      gcsbsb = -here->BSIM4v7capbs * ag0;
                  }
                  gcbbb = -(gcbgb + gcbdb + gcbsb + gcbgmb);

                  ggtg = ggtd = ggtb = ggts = 0.0;
                  sxpart = 0.4;
                  dxpart = 0.6;
                  ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                  dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
              }
              else
              {   qcheq = here->BSIM4v7qchqs;
                  CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV * here->BSIM4v7nf
                        * pParam->BSIM4v7leffCV;
                  T0 = qdef * ScalingFactor / CoxWL;
                  ggtg = here->BSIM4v7gtg = T0 * here->BSIM4v7gcrgg;
                  ggts = here->BSIM4v7gts = T0 * here->BSIM4v7gcrgd;
                  ggtd = here->BSIM4v7gtd = T0 * here->BSIM4v7gcrgs;
                  ggtb = here->BSIM4v7gtb = T0 * here->BSIM4v7gcrgb;
                  gqdef = ScalingFactor * ag0;

                  gcqgb = here->BSIM4v7cqgb * ag0;
                  gcqdb = here->BSIM4v7cqsb * ag0;
                  gcqsb = here->BSIM4v7cqdb * ag0;
                  gcqbb = here->BSIM4v7cqbb * ag0;

                  if (fabs(qcheq) <= 1.0e-5 * CoxWL)
                  {   if (model->BSIM4v7xpart < 0.5)
                      {   sxpart = 0.4;
                      }
                      else if (model->BSIM4v7xpart > 0.5)
                      {   sxpart = 0.0;
                      }
                      else
                      {   sxpart = 0.5;
                      }
                      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                  = dsxpart_dVs = 0.0;
                  }
                  else
                  {   sxpart = qdrn / qcheq;
                      Css = here->BSIM4v7cddb;
                      Cds = -(here->BSIM4v7cgdb + here->BSIM4v7cddb
                          + here->BSIM4v7cbdb);
                      dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                      Csg = here->BSIM4v7cdgb;
                      Cdg = -(here->BSIM4v7cggb + here->BSIM4v7cdgb
                          + here->BSIM4v7cbgb);
                      dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                      Csd = here->BSIM4v7cdsb;
                      Cdd = -(here->BSIM4v7cgsb + here->BSIM4v7cdsb
                          + here->BSIM4v7cbsb);
                      dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                  }
                  dxpart = 1.0 - sxpart;
                  ddxpart_dVd = -dsxpart_dVd;
                  ddxpart_dVg = -dsxpart_dVg;
                  ddxpart_dVs = -dsxpart_dVs;
                  ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                  if (here->BSIM4v7rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcdgb = gcsgb = gcbgb = 0.0;
                      gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                      qgmb = pParam->BSIM4v7cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qgate = 0.0;
                      qbulk = -qgmb;
                      qdrn = -qgdo;
                      qsrc = -qgso;
                  }
                  else
                  {   gcggb = (cgdo + cgso + pParam->BSIM4v7cgbo ) * ag0;
                      gcgdb = -cgdo * ag0;
                      gcgsb = -cgso * ag0;
                      gcgbb = -pParam->BSIM4v7cgbo * ag0;

                      gcdgb = gcgdb;
                      gcsgb = gcgsb;
                      gcbgb = gcgbb;
                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v7cgbo * vgb;
                      qgate = qgdo + qgso + qgb;
                      qbulk = -qgb;
                      qdrn = -qgdo;
                      qsrc = -qgso;
                  }

                  gcddb = (here->BSIM4v7capbd + cgdo) * ag0;
                  gcdsb = gcsdb = 0.0;
                  gcssb = (here->BSIM4v7capbs + cgso) * ag0;
                  if (!here->BSIM4v7rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                      gcsbb = -(gcsgb + gcssb + gcsgmb);
                      gcbdb = -here->BSIM4v7capbd * ag0;
                      gcbsb = -here->BSIM4v7capbs * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                      gcdbdb = -here->BSIM4v7capbd * ag0;
                      gcsbsb = -here->BSIM4v7capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
              }
          }


          if (here->BSIM4v7trnqsMod)
          {   *(ckt->CKTstate0 + here->BSIM4v7qcdump) = qdef * ScalingFactor;
              if (ckt->CKTmode & MODEINITTRAN)
                  *(ckt->CKTstate1 + here->BSIM4v7qcdump) =
                                   *(ckt->CKTstate0 + here->BSIM4v7qcdump);
              error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qcdump);
              if (error)
                  return(error);
          }

          if (ByPass) goto line860;

          *(ckt->CKTstate0 + here->BSIM4v7qg) = qgate;
          *(ckt->CKTstate0 + here->BSIM4v7qd) = qdrn
                           - *(ckt->CKTstate0 + here->BSIM4v7qbd);
          *(ckt->CKTstate0 + here->BSIM4v7qs) = qsrc
                           - *(ckt->CKTstate0 + here->BSIM4v7qbs);
          if (here->BSIM4v7rgateMod == 3)
              *(ckt->CKTstate0 + here->BSIM4v7qgmid) = qgmid;

          if (!here->BSIM4v7rbodyMod)
          {   *(ckt->CKTstate0 + here->BSIM4v7qb) = qbulk
                               + *(ckt->CKTstate0 + here->BSIM4v7qbd)
                               + *(ckt->CKTstate0 + here->BSIM4v7qbs);
          }
          else
              *(ckt->CKTstate0 + here->BSIM4v7qb) = qbulk;


          /* Store small signal parameters */
          if (ckt->CKTmode & MODEINITSMSIG)
          {   goto line1000;
          }

          if (!ChargeComputationNeeded)
              goto line850;

          if (ckt->CKTmode & MODEINITTRAN)
          {   *(ckt->CKTstate1 + here->BSIM4v7qb) =
                    *(ckt->CKTstate0 + here->BSIM4v7qb);
              *(ckt->CKTstate1 + here->BSIM4v7qg) =
                    *(ckt->CKTstate0 + here->BSIM4v7qg);
              *(ckt->CKTstate1 + here->BSIM4v7qd) =
                    *(ckt->CKTstate0 + here->BSIM4v7qd);
              if (here->BSIM4v7rgateMod == 3)
                  *(ckt->CKTstate1 + here->BSIM4v7qgmid) =
                        *(ckt->CKTstate0 + here->BSIM4v7qgmid);
              if (here->BSIM4v7rbodyMod)
              {   *(ckt->CKTstate1 + here->BSIM4v7qbs) =
                                   *(ckt->CKTstate0 + here->BSIM4v7qbs);
                  *(ckt->CKTstate1 + here->BSIM4v7qbd) =
                                   *(ckt->CKTstate0 + here->BSIM4v7qbd);
              }
          }

          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qb);
          if (error) 
              return(error);
          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qg);
          if (error) 
              return(error);
          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qd);
          if (error) 
              return(error);

          if (here->BSIM4v7rgateMod == 3)
          {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qgmid);
              if (error) return(error);
          }

          if (here->BSIM4v7rbodyMod)
          {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qbs);
              if (error) 
                  return(error);
              error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v7qbd);
              if (error) 
                  return(error);
          }

          goto line860;


line850:
          /* Zero gcap and ceqcap if (!ChargeComputationNeeded) */
          ceqqg = ceqqb = ceqqd = 0.0;
          ceqqjd = ceqqjs = 0.0;
          cqcheq = cqdef = 0.0;

          gcdgb = gcddb = gcdsb = gcdbb = 0.0;
          gcsgb = gcsdb = gcssb = gcsbb = 0.0;
          gcggb = gcgdb = gcgsb = gcgbb = 0.0;
          gcbdb = gcbgb = gcbsb = gcbbb = 0.0;

          gcgmgmb = gcgmdb = gcgmsb = gcgmbb = 0.0;
          gcdgmb = gcsgmb = gcbgmb = ceqqgmid = 0.0;
          gcdbdb = gcsbsb = 0.0;

          gqdef = gcqgb = gcqdb = gcqsb = gcqbb = 0.0;
          ggtg = ggtd = ggtb = ggts = 0.0;
          sxpart = (1.0 - (dxpart = (here->BSIM4v7mode > 0) ? 0.4 : 0.6));
          ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
          dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

          if (here->BSIM4v7trnqsMod)
          {   CoxWL = model->BSIM4v7coxe * pParam->BSIM4v7weffCV * here->BSIM4v7nf
                    * pParam->BSIM4v7leffCV;
              T1 = here->BSIM4v7gcrg / CoxWL;
              here->BSIM4v7gtau = T1 * ScalingFactor;
          }
          else
              here->BSIM4v7gtau = 0.0;

          goto line900;

            
line860:
          /* Calculate equivalent charge current */

          cqgate = *(ckt->CKTstate0 + here->BSIM4v7cqg);
          cqbody = *(ckt->CKTstate0 + here->BSIM4v7cqb);
          cqdrn = *(ckt->CKTstate0 + here->BSIM4v7cqd);

          ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
          ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb)
                * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
          ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb
                + gcbdb * vbd + gcbsb * vbs;


          if (here->BSIM4v7rgateMod == 3)
              ceqqgmid = *(ckt->CKTstate0 + here->BSIM4v7cqgmid)
                       + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;
          else
               ceqqgmid = 0.0;        

          if (here->BSIM4v7rbodyMod)
          {   ceqqjs = *(ckt->CKTstate0 + here->BSIM4v7cqbs) + gcsbsb * vbs_jct;
              ceqqjd = *(ckt->CKTstate0 + here->BSIM4v7cqbd) + gcdbdb * vbd_jct; 
          }

          if (here->BSIM4v7trnqsMod)
          {   T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
              ceqqg += T0;
              T1 = qdef * here->BSIM4v7gtau;
              ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                     * vbd - ddxpart_dVs * vbs);
              cqdef = *(ckt->CKTstate0 + here->BSIM4v7cqcdump) - gqdef * qdef;
              cqcheq = *(ckt->CKTstate0 + here->BSIM4v7cqcheq)
                     - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
          }

          if (ckt->CKTmode & MODEINITTRAN)
          {   *(ckt->CKTstate1 + here->BSIM4v7cqb) =
                               *(ckt->CKTstate0 + here->BSIM4v7cqb);
              *(ckt->CKTstate1 + here->BSIM4v7cqg) =
                               *(ckt->CKTstate0 + here->BSIM4v7cqg);
              *(ckt->CKTstate1 + here->BSIM4v7cqd) =
                               *(ckt->CKTstate0 + here->BSIM4v7cqd);

              if (here->BSIM4v7rgateMod == 3)
                  *(ckt->CKTstate1 + here->BSIM4v7cqgmid) =
                                   *(ckt->CKTstate0 + here->BSIM4v7cqgmid);

              if (here->BSIM4v7rbodyMod)
              {   *(ckt->CKTstate1 + here->BSIM4v7cqbs) =
                                   *(ckt->CKTstate0 + here->BSIM4v7cqbs);
                  *(ckt->CKTstate1 + here->BSIM4v7cqbd) =
                                   *(ckt->CKTstate0 + here->BSIM4v7cqbd);
              }
          }


          /*
           *  Load current vector
           */

line900:
          if (here->BSIM4v7mode >= 0)
          {   Gm = here->BSIM4v7gm;
              Gmbs = here->BSIM4v7gmbs;
              FwdSum = Gm + Gmbs;
              RevSum = 0.0;

              ceqdrn = model->BSIM4v7type * (cdrain - here->BSIM4v7gds * vds
                     - Gm * vgs - Gmbs * vbs);
              ceqbd = model->BSIM4v7type * (here->BSIM4v7csub + here->BSIM4v7Igidl
                    - (here->BSIM4v7gbds + here->BSIM4v7ggidld) * vds
                    - (here->BSIM4v7gbgs + here->BSIM4v7ggidlg) * vgs
                    - (here->BSIM4v7gbbs + here->BSIM4v7ggidlb) * vbs);
              ceqbs = model->BSIM4v7type * (here->BSIM4v7Igisl + here->BSIM4v7ggisls * vds 
                          - here->BSIM4v7ggislg * vgd - here->BSIM4v7ggislb * vbd);

              gbbdp = -(here->BSIM4v7gbds);
              gbbsp = here->BSIM4v7gbds + here->BSIM4v7gbgs + here->BSIM4v7gbbs; 
                    
              gbdpg = here->BSIM4v7gbgs;
              gbdpdp = here->BSIM4v7gbds;
              gbdpb = here->BSIM4v7gbbs;
              gbdpsp = -(gbdpg + gbdpdp + gbdpb);

              gbspg = 0.0;
              gbspdp = 0.0;
              gbspb = 0.0;
              gbspsp = 0.0;

              if (model->BSIM4v7igcMod)
              {   gIstotg = here->BSIM4v7gIgsg + here->BSIM4v7gIgcsg;
                  gIstotd = here->BSIM4v7gIgcsd;
                  gIstots = here->BSIM4v7gIgss + here->BSIM4v7gIgcss;
                  gIstotb = here->BSIM4v7gIgcsb;
                   Istoteq = model->BSIM4v7type * (here->BSIM4v7Igs + here->BSIM4v7Igcs
                            - gIstotg * vgs - here->BSIM4v7gIgcsd * vds
                          - here->BSIM4v7gIgcsb * vbs);

                  gIdtotg = here->BSIM4v7gIgdg + here->BSIM4v7gIgcdg;
                  gIdtotd = here->BSIM4v7gIgdd + here->BSIM4v7gIgcdd;
                  gIdtots = here->BSIM4v7gIgcds;
                  gIdtotb = here->BSIM4v7gIgcdb;
                  Idtoteq = model->BSIM4v7type * (here->BSIM4v7Igd + here->BSIM4v7Igcd
                          - here->BSIM4v7gIgdg * vgd - here->BSIM4v7gIgcdg * vgs
                          - here->BSIM4v7gIgcdd * vds - here->BSIM4v7gIgcdb * vbs);
              }
              else
              {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                  gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
              }

              if (model->BSIM4v7igbMod)
              {   gIbtotg = here->BSIM4v7gIgbg;
                  gIbtotd = here->BSIM4v7gIgbd;
                  gIbtots = here->BSIM4v7gIgbs;
                  gIbtotb = here->BSIM4v7gIgbb;
                  Ibtoteq = model->BSIM4v7type * (here->BSIM4v7Igb
                          - here->BSIM4v7gIgbg * vgs - here->BSIM4v7gIgbd * vds
                          - here->BSIM4v7gIgbb * vbs);
              }
              else
                  gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

              if ((model->BSIM4v7igcMod != 0) || (model->BSIM4v7igbMod != 0))
              {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                  gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                  gIgtots = gIstots + gIdtots + gIbtots;
                  gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  Igtoteq = Istoteq + Idtoteq + Ibtoteq; 
              }
              else
                  gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


              if (here->BSIM4v7rgateMod == 2)
                  T0 = vges - vgs;
              else if (here->BSIM4v7rgateMod == 3)
                  T0 = vgms - vgs;
              if (here->BSIM4v7rgateMod > 1)
              {   gcrgd = here->BSIM4v7gcrgd * T0;
                  gcrgg = here->BSIM4v7gcrgg * T0;
                  gcrgs = here->BSIM4v7gcrgs * T0;
                  gcrgb = here->BSIM4v7gcrgb * T0;
                  ceqgcrg = -(gcrgd * vds + gcrgg * vgs
                          + gcrgb * vbs);
                  gcrgg -= here->BSIM4v7gcrg;
                  gcrg = here->BSIM4v7gcrg;
              }
              else
                  ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
          }
          else
          {   Gm = -here->BSIM4v7gm;
              Gmbs = -here->BSIM4v7gmbs;
              FwdSum = 0.0;
              RevSum = -(Gm + Gmbs);

              ceqdrn = -model->BSIM4v7type * (cdrain + here->BSIM4v7gds * vds
                     + Gm * vgd + Gmbs * vbd);

              ceqbs = model->BSIM4v7type * (here->BSIM4v7csub + here->BSIM4v7Igisl 
                    + (here->BSIM4v7gbds + here->BSIM4v7ggisls) * vds
                    - (here->BSIM4v7gbgs + here->BSIM4v7ggislg) * vgd
                    - (here->BSIM4v7gbbs + here->BSIM4v7ggislb) * vbd);
              ceqbd = model->BSIM4v7type * (here->BSIM4v7Igidl - here->BSIM4v7ggidld * vds 
                              - here->BSIM4v7ggidlg * vgs - here->BSIM4v7ggidlb * vbs);

              gbbsp = -(here->BSIM4v7gbds);
              gbbdp = here->BSIM4v7gbds + here->BSIM4v7gbgs + here->BSIM4v7gbbs; 

              gbdpg = 0.0;
              gbdpsp = 0.0;
              gbdpb = 0.0;
              gbdpdp = 0.0;

              gbspg = here->BSIM4v7gbgs;
              gbspsp = here->BSIM4v7gbds;
              gbspb = here->BSIM4v7gbbs;
              gbspdp = -(gbspg + gbspsp + gbspb);

              if (model->BSIM4v7igcMod)
              {   gIstotg = here->BSIM4v7gIgsg + here->BSIM4v7gIgcdg;
                  gIstotd = here->BSIM4v7gIgcds;
                  gIstots = here->BSIM4v7gIgss + here->BSIM4v7gIgcdd;
                  gIstotb = here->BSIM4v7gIgcdb;
                  Istoteq = model->BSIM4v7type * (here->BSIM4v7Igs + here->BSIM4v7Igcd
                          - here->BSIM4v7gIgsg * vgs - here->BSIM4v7gIgcdg * vgd
                          + here->BSIM4v7gIgcdd * vds - here->BSIM4v7gIgcdb * vbd);

                  gIdtotg = here->BSIM4v7gIgdg + here->BSIM4v7gIgcsg;
                  gIdtotd = here->BSIM4v7gIgdd + here->BSIM4v7gIgcss;
                  gIdtots = here->BSIM4v7gIgcsd;
                  gIdtotb = here->BSIM4v7gIgcsb;
                  Idtoteq = model->BSIM4v7type * (here->BSIM4v7Igd + here->BSIM4v7Igcs
                          - (here->BSIM4v7gIgdg + here->BSIM4v7gIgcsg) * vgd
                          + here->BSIM4v7gIgcsd * vds - here->BSIM4v7gIgcsb * vbd);
              }
              else
              {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                  gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
              }

              if (model->BSIM4v7igbMod)
              {   gIbtotg = here->BSIM4v7gIgbg;
                  gIbtotd = here->BSIM4v7gIgbs;
                  gIbtots = here->BSIM4v7gIgbd;
                  gIbtotb = here->BSIM4v7gIgbb;
                  Ibtoteq = model->BSIM4v7type * (here->BSIM4v7Igb
                          - here->BSIM4v7gIgbg * vgd + here->BSIM4v7gIgbd * vds
                          - here->BSIM4v7gIgbb * vbd);
              }
              else
                  gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

              if ((model->BSIM4v7igcMod != 0) || (model->BSIM4v7igbMod != 0))
              {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                  gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                  gIgtots = gIstots + gIdtots + gIbtots;
                  gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  Igtoteq = Istoteq + Idtoteq + Ibtoteq;
              }
              else
                  gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


              if (here->BSIM4v7rgateMod == 2)
                  T0 = vges - vgs;
              else if (here->BSIM4v7rgateMod == 3)
                  T0 = vgms - vgs;
              if (here->BSIM4v7rgateMod > 1)
              {   gcrgd = here->BSIM4v7gcrgs * T0;
                  gcrgg = here->BSIM4v7gcrgg * T0;
                  gcrgs = here->BSIM4v7gcrgd * T0;
                  gcrgb = here->BSIM4v7gcrgb * T0;
                  ceqgcrg = -(gcrgg * vgd - gcrgs * vds
                          + gcrgb * vbd);
                  gcrgg -= here->BSIM4v7gcrg;
                  gcrg = here->BSIM4v7gcrg;
              }
              else
                  ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
          }

          if (model->BSIM4v7rdsMod == 1)
          {   ceqgstot = model->BSIM4v7type * (here->BSIM4v7gstotd * vds
                       + here->BSIM4v7gstotg * vgs + here->BSIM4v7gstotb * vbs);
              /* WDLiu: ceqgstot flowing away from sNodePrime */
              gstot = here->BSIM4v7gstot;
              gstotd = here->BSIM4v7gstotd;
              gstotg = here->BSIM4v7gstotg;
              gstots = here->BSIM4v7gstots - gstot;
              gstotb = here->BSIM4v7gstotb;

              ceqgdtot = -model->BSIM4v7type * (here->BSIM4v7gdtotd * vds
                       + here->BSIM4v7gdtotg * vgs + here->BSIM4v7gdtotb * vbs);
              /* WDLiu: ceqgdtot defined as flowing into dNodePrime */
              gdtot = here->BSIM4v7gdtot;
              gdtotd = here->BSIM4v7gdtotd - gdtot;
              gdtotg = here->BSIM4v7gdtotg;
              gdtots = here->BSIM4v7gdtots;
              gdtotb = here->BSIM4v7gdtotb;
          }
          else
          {   gstot = gstotd = gstotg = gstots = gstotb = ceqgstot = 0.0;
              gdtot = gdtotd = gdtotg = gdtots = gdtotb = ceqgdtot = 0.0;
          }

           if (model->BSIM4v7type > 0)
           {   ceqjs = (here->BSIM4v7cbs - here->BSIM4v7gbs * vbs_jct);
               ceqjd = (here->BSIM4v7cbd - here->BSIM4v7gbd * vbd_jct);
           }
           else
           {   ceqjs = -(here->BSIM4v7cbs - here->BSIM4v7gbs * vbs_jct); 
               ceqjd = -(here->BSIM4v7cbd - here->BSIM4v7gbd * vbd_jct);
               ceqqg = -ceqqg;
               ceqqd = -ceqqd;
               ceqqb = -ceqqb;
               ceqgcrg = -ceqgcrg;

               if (here->BSIM4v7trnqsMod)
               {   cqdef = -cqdef;
                   cqcheq = -cqcheq;
               }

               if (here->BSIM4v7rbodyMod)
               {   ceqqjs = -ceqqjs;
                   ceqqjd = -ceqqjd;
               }

               if (here->BSIM4v7rgateMod == 3)
                   ceqqgmid = -ceqqgmid;
           }


           /*
            *  Loading RHS
            */

              m = here->BSIM4v7m;

#ifdef USE_OMP
       here->BSIM4v7rhsdPrime = m * (ceqjd - ceqbd + ceqgdtot
                                                    - ceqdrn - ceqqd + Idtoteq);
       here->BSIM4v7rhsgPrime = m * (ceqqg - ceqgcrg + Igtoteq);

       if (here->BSIM4v7rgateMod == 2)
           here->BSIM4v7rhsgExt = m * ceqgcrg;
       else if (here->BSIM4v7rgateMod == 3)
               here->BSIM4v7grhsMid = m * (ceqqgmid + ceqgcrg);

       if (!here->BSIM4v7rbodyMod)
       {   here->BSIM4v7rhsbPrime = m * (ceqbd + ceqbs - ceqjd
                                                        - ceqjs - ceqqb + Ibtoteq);
           here->BSIM4v7rhssPrime = m * (ceqdrn - ceqbs + ceqjs 
                              + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq);
        }
        else
        {   here->BSIM4v7rhsdb = m * (ceqjd + ceqqjd);
            here->BSIM4v7rhsbPrime = m * (ceqbd + ceqbs - ceqqb + Ibtoteq);
            here->BSIM4v7rhssb = m * (ceqjs + ceqqjs);
            here->BSIM4v7rhssPrime = m * (ceqdrn - ceqbs + ceqjs + ceqqd 
                + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq);
        }

        if (model->BSIM4v7rdsMod)
        {   here->BSIM4v7rhsd = m * ceqgdtot; 
            here->BSIM4v7rhss = m * ceqgstot;
        }

        if (here->BSIM4v7trnqsMod)
           here->BSIM4v7rhsq = m * (cqcheq - cqdef);
#else
        (*(ckt->CKTrhs + here->BSIM4v7dNodePrime) += m * (ceqjd - ceqbd + ceqgdtot
                                                    - ceqdrn - ceqqd + Idtoteq));
        (*(ckt->CKTrhs + here->BSIM4v7gNodePrime) -= m * (ceqqg - ceqgcrg + Igtoteq));

        if (here->BSIM4v7rgateMod == 2)
            (*(ckt->CKTrhs + here->BSIM4v7gNodeExt) -= m * ceqgcrg);
        else if (here->BSIM4v7rgateMod == 3)
            (*(ckt->CKTrhs + here->BSIM4v7gNodeMid) -= m * (ceqqgmid + ceqgcrg));

        if (!here->BSIM4v7rbodyMod)
        {   (*(ckt->CKTrhs + here->BSIM4v7bNodePrime) += m * (ceqbd + ceqbs - ceqjd
                                                        - ceqjs - ceqqb + Ibtoteq));
            (*(ckt->CKTrhs + here->BSIM4v7sNodePrime) += m * (ceqdrn - ceqbs + ceqjs 
                              + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq));
        }

        else
        {   (*(ckt->CKTrhs + here->BSIM4v7dbNode) -= m * (ceqjd + ceqqjd));
            (*(ckt->CKTrhs + here->BSIM4v7bNodePrime) += m * (ceqbd + ceqbs - ceqqb + Ibtoteq));
            (*(ckt->CKTrhs + here->BSIM4v7sbNode) -= m * (ceqjs + ceqqjs));
            (*(ckt->CKTrhs + here->BSIM4v7sNodePrime) += m * (ceqdrn - ceqbs + ceqjs + ceqqd 
                + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq));
        }

        if (model->BSIM4v7rdsMod)
        {   (*(ckt->CKTrhs + here->BSIM4v7dNode) -= m * ceqgdtot); 
            (*(ckt->CKTrhs + here->BSIM4v7sNode) += m * ceqgstot);
        }

        if (here->BSIM4v7trnqsMod)
            *(ckt->CKTrhs + here->BSIM4v7qNode) += m * (cqcheq - cqdef);
#endif

        /*
         *  Loading matrix
         */

           if (!here->BSIM4v7rbodyMod)
           {   gjbd = here->BSIM4v7gbd;
               gjbs = here->BSIM4v7gbs;
           }
       else
          gjbd = gjbs = 0.0;

       if (!model->BSIM4v7rdsMod)
       {   gdpr = here->BSIM4v7drainConductance;
           gspr = here->BSIM4v7sourceConductance;
       }
       else
           gdpr = gspr = 0.0;

       geltd = here->BSIM4v7grgeltd;

       T1 = qdef * here->BSIM4v7gtau;
#ifdef USE_OMP
       if (here->BSIM4v7rgateMod == 1)
       {   here->BSIM4v7_1 = m * geltd;
           here->BSIM4v7_2 = m * geltd;
           here->BSIM4v7_3 = m * geltd;
           here->BSIM4v7_4 = m * (gcggb + geltd - ggtg + gIgtotg);
           here->BSIM4v7_5 = m * (gcgdb - ggtd + gIgtotd);
           here->BSIM4v7_6 = m * (gcgsb - ggts + gIgtots);
           here->BSIM4v7_7 = m * (gcgbb - ggtb + gIgtotb);
       } /* WDLiu: gcrg already subtracted from all gcrgg below */
       else if (here->BSIM4v7rgateMod == 2)        
       {   here->BSIM4v7_8 = m * gcrg;
           here->BSIM4v7_9 = m * gcrgg;
           here->BSIM4v7_10 = m * gcrgd;
           here->BSIM4v7_11 = m * gcrgs;
           here->BSIM4v7_12 = m * gcrgb;        

           here->BSIM4v7_13 = m * gcrg;
           here->BSIM4v7_14 = m * (gcggb  - gcrgg - ggtg + gIgtotg);
           here->BSIM4v7_15 = m * (gcgdb - gcrgd - ggtd + gIgtotd);
           here->BSIM4v7_16 = m * (gcgsb - gcrgs - ggts + gIgtots);
           here->BSIM4v7_17 = m * (gcgbb - gcrgb - ggtb + gIgtotb);
       }
       else if (here->BSIM4v7rgateMod == 3)
       {   here->BSIM4v7_18 = m * geltd;
           here->BSIM4v7_19 = m * geltd;
           here->BSIM4v7_20 = m * geltd;
           here->BSIM4v7_21 = m * (geltd + gcrg + gcgmgmb);

           here->BSIM4v7_22 = m * (gcrgd + gcgmdb);
           here->BSIM4v7_23 = m * gcrgg;
           here->BSIM4v7_24 = m * (gcrgs + gcgmsb);
           here->BSIM4v7_25 = m * (gcrgb + gcgmbb);

           here->BSIM4v7_26 = m * gcdgmb;
           here->BSIM4v7_26 = m * gcrg;
           here->BSIM4v7_28 = m * gcsgmb;
           here->BSIM4v7_29 = m * gcbgmb;

           here->BSIM4v7_30 = m * (gcggb - gcrgg - ggtg + gIgtotg);
           here->BSIM4v7_31 = m * (gcgdb - gcrgd - ggtd + gIgtotd);
           here->BSIM4v7_32 = m * (gcgsb - gcrgs - ggts + gIgtots);
           here->BSIM4v7_33 = m * (gcgbb - gcrgb - ggtb + gIgtotb);
       }
       else
       {   here->BSIM4v7_34 = m * (gcggb - ggtg + gIgtotg);
           here->BSIM4v7_35 = m * (gcgdb - ggtd + gIgtotd);
           here->BSIM4v7_36 = m * (gcgsb - ggts + gIgtots);
           here->BSIM4v7_37 = m * (gcgbb - ggtb + gIgtotb);
       }

       if (model->BSIM4v7rdsMod)
       {   here->BSIM4v7_38 = m * gdtotg;
           here->BSIM4v7_39 = m * gdtots;
           here->BSIM4v7_40 = m * gdtotb;
           here->BSIM4v7_41 = m * gstotd;
           here->BSIM4v7_42 = m * gstotg;
           here->BSIM4v7_43 = m * gstotb;
       }

       here->BSIM4v7_44 = m * (gdpr + here->BSIM4v7gds + here->BSIM4v7gbd + T1 * ddxpart_dVd
                                   - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd);
       here->BSIM4v7_45 = m * (gdpr + gdtot);
       here->BSIM4v7_46 = m * (Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                                   + dxpart * ggtg + T1 * ddxpart_dVg);
       here->BSIM4v7_47 = m * (here->BSIM4v7gds + gdtots - dxpart * ggts + gIdtots
                                   - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp);
       here->BSIM4v7_48 = m * (gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                                   - T1 * ddxpart_dVb - dxpart * ggtb);

       here->BSIM4v7_49 = m * (gdpr - gdtotd);
       here->BSIM4v7_50 = m * (gdpr + gdtot);

       here->BSIM4v7_51 = m * (here->BSIM4v7gds + gstotd + RevSum - gcsdb - gbspdp
                                   - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd);
       here->BSIM4v7_52 = m * (gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                                   + T1 * dsxpart_dVg - gIstotg);
       here->BSIM4v7_53 = m * (gspr + here->BSIM4v7gds + here->BSIM4v7gbs + T1 * dsxpart_dVs
                                   - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots);
       here->BSIM4v7_54 = m * (gspr + gstot);
       here->BSIM4v7_55 = m * (gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                                   - T1 * dsxpart_dVb + gIstotb);

       here->BSIM4v7_56 = m * (gspr - gstots);
       here->BSIM4v7_57 = m * (gspr + gstot);

       here->BSIM4v7_58 = m * (gcbdb - gjbd + gbbdp - gIbtotd);
       here->BSIM4v7_59 = m * (gcbgb - here->BSIM4v7gbgs - gIbtotg);
       here->BSIM4v7_60 = m * (gcbsb - gjbs + gbbsp - gIbtots);
       here->BSIM4v7_61 = m * (gjbd + gjbs + gcbbb - here->BSIM4v7gbbs - gIbtotb);

       ggidld = here->BSIM4v7ggidld;
       ggidlg = here->BSIM4v7ggidlg;
       ggidlb = here->BSIM4v7ggidlb;
       ggislg = here->BSIM4v7ggislg;
       ggisls = here->BSIM4v7ggisls;
       ggislb = here->BSIM4v7ggislb;

       /* stamp gidl */
       here->BSIM4v7_62 = m * ggidld;
       here->BSIM4v7_63 = m * ggidlg;
       here->BSIM4v7_64 = m * (ggidlg + ggidld + ggidlb);
       here->BSIM4v7_65 = m * ggidlb;
       here->BSIM4v7_66 = m * ggidld;
       here->BSIM4v7_67 = m * ggidlg;
       here->BSIM4v7_68 = m * (ggidlg + ggidld + ggidlb);
       here->BSIM4v7_69 = m * ggidlb;
       /* stamp gisl */
       here->BSIM4v7_70 = m * (ggisls + ggislg + ggislb);
       here->BSIM4v7_71 = m * ggislg;
       here->BSIM4v7_72 = m * ggisls;
       here->BSIM4v7_73 = m * ggislb;
       here->BSIM4v7_74 = m * (ggislg + ggisls + ggislb);
       here->BSIM4v7_75 = m * ggislg;
       here->BSIM4v7_76 = m * ggisls;
       here->BSIM4v7_77 = m * ggislb;

       if (here->BSIM4v7rbodyMod)
       {   here->BSIM4v7_78 = m * (gcdbdb - here->BSIM4v7gbd);
           here->BSIM4v7_79 = m * (here->BSIM4v7gbs - gcsbsb);

           here->BSIM4v7_80 = m * (gcdbdb - here->BSIM4v7gbd);
           here->BSIM4v7_81 = m * (here->BSIM4v7gbd - gcdbdb 
                          + here->BSIM4v7grbpd + here->BSIM4v7grbdb);
           here->BSIM4v7_82 = m * here->BSIM4v7grbpd;
           here->BSIM4v7_83 = m * here->BSIM4v7grbdb;

           here->BSIM4v7_84 = m * here->BSIM4v7grbpd;
           here->BSIM4v7_85 = m * here->BSIM4v7grbpb;
           here->BSIM4v7_86 = m * here->BSIM4v7grbps;
           here->BSIM4v7_87 = m * (here->BSIM4v7grbpd + here->BSIM4v7grbps 
                          + here->BSIM4v7grbpb);
           /* WDLiu: (gcbbb - here->BSIM4v7gbbs) already added to BPbpPtr */        

           here->BSIM4v7_88 = m * (gcsbsb - here->BSIM4v7gbs);
           here->BSIM4v7_89 = m * here->BSIM4v7grbps;
           here->BSIM4v7_90 = m * here->BSIM4v7grbsb;
           here->BSIM4v7_91 = m * (here->BSIM4v7gbs - gcsbsb 
                          + here->BSIM4v7grbps + here->BSIM4v7grbsb);

           here->BSIM4v7_92 = m * here->BSIM4v7grbdb;
           here->BSIM4v7_93 = m * here->BSIM4v7grbpb;
           here->BSIM4v7_94 = m * here->BSIM4v7grbsb;
           here->BSIM4v7_95 = m * (here->BSIM4v7grbsb + here->BSIM4v7grbdb
                           + here->BSIM4v7grbpb);
       }

           if (here->BSIM4v7trnqsMod)
           {   here->BSIM4v7_96 = m * (gqdef + here->BSIM4v7gtau);
               here->BSIM4v7_97 = m * (ggtg - gcqgb);
               here->BSIM4v7_98 = m * (ggtd - gcqdb);
               here->BSIM4v7_99 = m * (ggts - gcqsb);
               here->BSIM4v7_100 = m * (ggtb - gcqbb);

               here->BSIM4v7_101 = m * dxpart * here->BSIM4v7gtau;
               here->BSIM4v7_102 = m * sxpart * here->BSIM4v7gtau;
               here->BSIM4v7_103 = m * here->BSIM4v7gtau;
           }
#else
           if (here->BSIM4v7rgateMod == 1)
           {   (*(here->BSIM4v7GEgePtr) += m * geltd);
               (*(here->BSIM4v7GPgePtr) -= m * geltd);
               (*(here->BSIM4v7GEgpPtr) -= m * geltd);
               (*(here->BSIM4v7GPgpPtr) += m * (gcggb + geltd - ggtg + gIgtotg));
               (*(here->BSIM4v7GPdpPtr) += m * (gcgdb - ggtd + gIgtotd));
               (*(here->BSIM4v7GPspPtr) += m * (gcgsb - ggts + gIgtots));
               (*(here->BSIM4v7GPbpPtr) += m * (gcgbb - ggtb + gIgtotb));
           } /* WDLiu: gcrg already subtracted from all gcrgg below */
           else if (here->BSIM4v7rgateMod == 2)        
           {   (*(here->BSIM4v7GEgePtr) += m * gcrg);
               (*(here->BSIM4v7GEgpPtr) += m * gcrgg);
               (*(here->BSIM4v7GEdpPtr) += m * gcrgd);
               (*(here->BSIM4v7GEspPtr) += m * gcrgs);
               (*(here->BSIM4v7GEbpPtr) += m * gcrgb);        

               (*(here->BSIM4v7GPgePtr) -= m * gcrg);
               (*(here->BSIM4v7GPgpPtr) += m * (gcggb  - gcrgg - ggtg + gIgtotg));
               (*(here->BSIM4v7GPdpPtr) += m * (gcgdb - gcrgd - ggtd + gIgtotd));
               (*(here->BSIM4v7GPspPtr) += m * (gcgsb - gcrgs - ggts + gIgtots));
               (*(here->BSIM4v7GPbpPtr) += m * (gcgbb - gcrgb - ggtb + gIgtotb));
           }
           else if (here->BSIM4v7rgateMod == 3)
           {   (*(here->BSIM4v7GEgePtr) += m * geltd);
               (*(here->BSIM4v7GEgmPtr) -= m * geltd);
               (*(here->BSIM4v7GMgePtr) -= m * geltd);
               (*(here->BSIM4v7GMgmPtr) += m * (geltd + gcrg + gcgmgmb));

               (*(here->BSIM4v7GMdpPtr) += m * (gcrgd + gcgmdb));
               (*(here->BSIM4v7GMgpPtr) += m * gcrgg);
               (*(here->BSIM4v7GMspPtr) += m * (gcrgs + gcgmsb));
               (*(here->BSIM4v7GMbpPtr) += m * (gcrgb + gcgmbb));

               (*(here->BSIM4v7DPgmPtr) += m * gcdgmb);
               (*(here->BSIM4v7GPgmPtr) -= m * gcrg);
               (*(here->BSIM4v7SPgmPtr) += m * gcsgmb);
               (*(here->BSIM4v7BPgmPtr) += m * gcbgmb);

               (*(here->BSIM4v7GPgpPtr) += m * (gcggb - gcrgg - ggtg + gIgtotg));
               (*(here->BSIM4v7GPdpPtr) += m * (gcgdb - gcrgd - ggtd + gIgtotd));
               (*(here->BSIM4v7GPspPtr) += m * (gcgsb - gcrgs - ggts + gIgtots));
               (*(here->BSIM4v7GPbpPtr) += m * (gcgbb - gcrgb - ggtb + gIgtotb));
           }
            else
           {   (*(here->BSIM4v7GPgpPtr) += m * (gcggb - ggtg + gIgtotg));
               (*(here->BSIM4v7GPdpPtr) += m * (gcgdb - ggtd + gIgtotd));
               (*(here->BSIM4v7GPspPtr) += m * (gcgsb - ggts + gIgtots));
               (*(here->BSIM4v7GPbpPtr) += m * (gcgbb - ggtb + gIgtotb));
           }

           if (model->BSIM4v7rdsMod)
           {   (*(here->BSIM4v7DgpPtr) += m * gdtotg);
               (*(here->BSIM4v7DspPtr) += m * gdtots);
               (*(here->BSIM4v7DbpPtr) += m * gdtotb);
               (*(here->BSIM4v7SdpPtr) += m * gstotd);
               (*(here->BSIM4v7SgpPtr) += m * gstotg);
               (*(here->BSIM4v7SbpPtr) += m * gstotb);
           }

           (*(here->BSIM4v7DPdpPtr) += m * (gdpr + here->BSIM4v7gds + here->BSIM4v7gbd + T1 * ddxpart_dVd
                                   - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd));
           (*(here->BSIM4v7DPdPtr) -= m * (gdpr + gdtot));
           (*(here->BSIM4v7DPgpPtr) += m * (Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                                   + dxpart * ggtg + T1 * ddxpart_dVg));
           (*(here->BSIM4v7DPspPtr) -= m * (here->BSIM4v7gds + gdtots - dxpart * ggts + gIdtots
                                   - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp));
           (*(here->BSIM4v7DPbpPtr) -= m * (gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                                   - T1 * ddxpart_dVb - dxpart * ggtb));

           (*(here->BSIM4v7DdpPtr) -= m * (gdpr - gdtotd));
           (*(here->BSIM4v7DdPtr) += m * (gdpr + gdtot));

           (*(here->BSIM4v7SPdpPtr) -= m * (here->BSIM4v7gds + gstotd + RevSum - gcsdb - gbspdp
                                   - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd));
           (*(here->BSIM4v7SPgpPtr) += m * (gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                                   + T1 * dsxpart_dVg - gIstotg));
           (*(here->BSIM4v7SPspPtr) += m * (gspr + here->BSIM4v7gds + here->BSIM4v7gbs + T1 * dsxpart_dVs
                                   - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots));
           (*(here->BSIM4v7SPsPtr) -= m * (gspr + gstot));
           (*(here->BSIM4v7SPbpPtr) -= m * (gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                                   - T1 * dsxpart_dVb + gIstotb));

           (*(here->BSIM4v7SspPtr) -= m * (gspr - gstots));
           (*(here->BSIM4v7SsPtr) += m * (gspr + gstot));

           (*(here->BSIM4v7BPdpPtr) += m * (gcbdb - gjbd + gbbdp - gIbtotd));
           (*(here->BSIM4v7BPgpPtr) += m * (gcbgb - here->BSIM4v7gbgs - gIbtotg));
           (*(here->BSIM4v7BPspPtr) += m * (gcbsb - gjbs + gbbsp - gIbtots));
           (*(here->BSIM4v7BPbpPtr) += m * (gjbd + gjbs + gcbbb - here->BSIM4v7gbbs
                                   - gIbtotb));

           ggidld = here->BSIM4v7ggidld;
           ggidlg = here->BSIM4v7ggidlg;
           ggidlb = here->BSIM4v7ggidlb;
           ggislg = here->BSIM4v7ggislg;
           ggisls = here->BSIM4v7ggisls;
           ggislb = here->BSIM4v7ggislb;

           /* stamp gidl */
           (*(here->BSIM4v7DPdpPtr) += m * ggidld);
           (*(here->BSIM4v7DPgpPtr) += m * ggidlg);
           (*(here->BSIM4v7DPspPtr) -= m * (ggidlg + ggidld + ggidlb));
           (*(here->BSIM4v7DPbpPtr) += m * ggidlb);
           (*(here->BSIM4v7BPdpPtr) -= m * ggidld);
           (*(here->BSIM4v7BPgpPtr) -= m * ggidlg);
           (*(here->BSIM4v7BPspPtr) += m * (ggidlg + ggidld + ggidlb));
           (*(here->BSIM4v7BPbpPtr) -= m * ggidlb);
            /* stamp gisl */
           (*(here->BSIM4v7SPdpPtr) -= m * (ggisls + ggislg + ggislb));
           (*(here->BSIM4v7SPgpPtr) += m * ggislg);
           (*(here->BSIM4v7SPspPtr) += m * ggisls);
           (*(here->BSIM4v7SPbpPtr) += m * ggislb);
           (*(here->BSIM4v7BPdpPtr) += m * (ggislg + ggisls + ggislb));
           (*(here->BSIM4v7BPgpPtr) -= m * ggislg);
           (*(here->BSIM4v7BPspPtr) -= m * ggisls);
           (*(here->BSIM4v7BPbpPtr) -= m * ggislb);


           if (here->BSIM4v7rbodyMod)
           {   (*(here->BSIM4v7DPdbPtr) += m * (gcdbdb - here->BSIM4v7gbd));
               (*(here->BSIM4v7SPsbPtr) -= m * (here->BSIM4v7gbs - gcsbsb));

               (*(here->BSIM4v7DBdpPtr) += m * (gcdbdb - here->BSIM4v7gbd));
               (*(here->BSIM4v7DBdbPtr) += m * (here->BSIM4v7gbd - gcdbdb 
                                       + here->BSIM4v7grbpd + here->BSIM4v7grbdb));
               (*(here->BSIM4v7DBbpPtr) -= m * here->BSIM4v7grbpd);
               (*(here->BSIM4v7DBbPtr) -= m * here->BSIM4v7grbdb);

               (*(here->BSIM4v7BPdbPtr) -= m * here->BSIM4v7grbpd);
               (*(here->BSIM4v7BPbPtr) -= m * here->BSIM4v7grbpb);
               (*(here->BSIM4v7BPsbPtr) -= m * here->BSIM4v7grbps);
               (*(here->BSIM4v7BPbpPtr) += m * (here->BSIM4v7grbpd + here->BSIM4v7grbps 
                                       + here->BSIM4v7grbpb));
               /* WDLiu: (gcbbb - here->BSIM4v7gbbs) already added to BPbpPtr */        

               (*(here->BSIM4v7SBspPtr) += m * (gcsbsb - here->BSIM4v7gbs));
               (*(here->BSIM4v7SBbpPtr) -= m * here->BSIM4v7grbps);
               (*(here->BSIM4v7SBbPtr) -= m * here->BSIM4v7grbsb);
               (*(here->BSIM4v7SBsbPtr) += m * (here->BSIM4v7gbs - gcsbsb 
                                       + here->BSIM4v7grbps + here->BSIM4v7grbsb));

               (*(here->BSIM4v7BdbPtr) -= m * here->BSIM4v7grbdb);
               (*(here->BSIM4v7BbpPtr) -= m * here->BSIM4v7grbpb);
               (*(here->BSIM4v7BsbPtr) -= m * here->BSIM4v7grbsb);
               (*(here->BSIM4v7BbPtr) += m * (here->BSIM4v7grbsb + here->BSIM4v7grbdb
                                     + here->BSIM4v7grbpb));
           }

           if (here->BSIM4v7trnqsMod)
           {   (*(here->BSIM4v7QqPtr) += m * (gqdef + here->BSIM4v7gtau));
               (*(here->BSIM4v7QgpPtr) += m * (ggtg - gcqgb));
               (*(here->BSIM4v7QdpPtr) += m * (ggtd - gcqdb));
               (*(here->BSIM4v7QspPtr) += m * (ggts - gcqsb));
               (*(here->BSIM4v7QbpPtr) += m * (ggtb - gcqbb));

               (*(here->BSIM4v7DPqPtr) += m * dxpart * here->BSIM4v7gtau);
               (*(here->BSIM4v7SPqPtr) += m * sxpart * here->BSIM4v7gtau);
               (*(here->BSIM4v7GPqPtr) -= m * here->BSIM4v7gtau);
           }
#endif

line1000:  ;

#ifndef USE_OMP
     }  /* End of MOSFET Instance */
}   /* End of Model Instance */
#endif

return(OK);
}

/* function to compute poly depletion effect */
int BSIM4v7polyDepletion(
    double  phi,
    double  ngate,
    double  epsgate,
    double  coxe,
    double  Vgs,
    double *Vgs_eff,
    double *dVgs_eff_dVg)
{
    double T1, T2, T3, T4, T5, T6, T7, T8;

    /* Poly Gate Si Depletion Effect */
    if ((ngate > 1.0e18) && 
        (ngate < 1.0e25) && (Vgs > phi) && (epsgate!=0)
       ){
        T1 = 1.0e6 * CHARGE * epsgate * ngate / (coxe * coxe);
        T8 = Vgs - phi;
        T4 = sqrt(1.0 + 2.0 * T8 / T1);
        T2 = 2.0 * T8 / (T4 + 1.0);
        T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
        T7 = 1.12 - T3 - 0.05;
        T6 = sqrt(T7 * T7 + 0.224);
        T5 = 1.12 - 0.5 * (T7 + T6);
        *Vgs_eff = Vgs - T5;
        *dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6); 
    }
    else {
        *Vgs_eff = Vgs;
        *dVgs_eff_dVg = 1.0;
    }
    return(0);
}

#ifdef USE_OMP
void BSIM4v7LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt)
{
    unsigned int InstCount, idx;
    BSIM4v7instance **InstArray;
    BSIM4v7instance *here;
    BSIM4v7model *model = (BSIM4v7model*)inModel;

    InstArray = model->BSIM4v7InstanceArray;
    InstCount = model->BSIM4v7InstCount;

    for(idx = 0; idx < InstCount; idx++) {
       here = InstArray[idx];
        /* Update b for Ax = b */
           (*(ckt->CKTrhs + here->BSIM4v7dNodePrime) += here->BSIM4v7rhsdPrime);
           (*(ckt->CKTrhs + here->BSIM4v7gNodePrime) -= here->BSIM4v7rhsgPrime);

           if (here->BSIM4v7rgateMod == 2)
               (*(ckt->CKTrhs + here->BSIM4v7gNodeExt) -= here->BSIM4v7rhsgExt);
           else if (here->BSIM4v7rgateMod == 3)
               (*(ckt->CKTrhs + here->BSIM4v7gNodeMid) -= here->BSIM4v7grhsMid);

           if (!here->BSIM4v7rbodyMod)
           {   (*(ckt->CKTrhs + here->BSIM4v7bNodePrime) += here->BSIM4v7rhsbPrime);
               (*(ckt->CKTrhs + here->BSIM4v7sNodePrime) += here->BSIM4v7rhssPrime);
           }
           else
           {   (*(ckt->CKTrhs + here->BSIM4v7dbNode) -= here->BSIM4v7rhsdb);
               (*(ckt->CKTrhs + here->BSIM4v7bNodePrime) += here->BSIM4v7rhsbPrime);
               (*(ckt->CKTrhs + here->BSIM4v7sbNode) -= here->BSIM4v7rhssb);
               (*(ckt->CKTrhs + here->BSIM4v7sNodePrime) += here->BSIM4v7rhssPrime);
           }

           if (model->BSIM4v7rdsMod)
           {   (*(ckt->CKTrhs + here->BSIM4v7dNode) -= here->BSIM4v7rhsd); 
               (*(ckt->CKTrhs + here->BSIM4v7sNode) += here->BSIM4v7rhss);
           }

           if (here->BSIM4v7trnqsMod)
               *(ckt->CKTrhs + here->BSIM4v7qNode) += here->BSIM4v7rhsq;


        /* Update A for Ax = b */
           if (here->BSIM4v7rgateMod == 1)
           {   (*(here->BSIM4v7GEgePtr) += here->BSIM4v7_1);
               (*(here->BSIM4v7GPgePtr) -= here->BSIM4v7_2);
               (*(here->BSIM4v7GEgpPtr) -= here->BSIM4v7_3);
               (*(here->BSIM4v7GPgpPtr) += here->BSIM4v7_4);
               (*(here->BSIM4v7GPdpPtr) += here->BSIM4v7_5);
               (*(here->BSIM4v7GPspPtr) += here->BSIM4v7_6);
               (*(here->BSIM4v7GPbpPtr) += here->BSIM4v7_7);
           }
           else if (here->BSIM4v7rgateMod == 2)        
           {   (*(here->BSIM4v7GEgePtr) += here->BSIM4v7_8);
               (*(here->BSIM4v7GEgpPtr) += here->BSIM4v7_9);
               (*(here->BSIM4v7GEdpPtr) += here->BSIM4v7_10);
               (*(here->BSIM4v7GEspPtr) += here->BSIM4v7_11);
               (*(here->BSIM4v7GEbpPtr) += here->BSIM4v7_12);        

               (*(here->BSIM4v7GPgePtr) -= here->BSIM4v7_13);
               (*(here->BSIM4v7GPgpPtr) += here->BSIM4v7_14);
               (*(here->BSIM4v7GPdpPtr) += here->BSIM4v7_15);
               (*(here->BSIM4v7GPspPtr) += here->BSIM4v7_16);
               (*(here->BSIM4v7GPbpPtr) += here->BSIM4v7_17);
           }
           else if (here->BSIM4v7rgateMod == 3)
           {   (*(here->BSIM4v7GEgePtr) += here->BSIM4v7_18);
               (*(here->BSIM4v7GEgmPtr) -= here->BSIM4v7_19);
               (*(here->BSIM4v7GMgePtr) -= here->BSIM4v7_20);
               (*(here->BSIM4v7GMgmPtr) += here->BSIM4v7_21);

               (*(here->BSIM4v7GMdpPtr) += here->BSIM4v7_22);
               (*(here->BSIM4v7GMgpPtr) += here->BSIM4v7_23);
               (*(here->BSIM4v7GMspPtr) += here->BSIM4v7_24);
               (*(here->BSIM4v7GMbpPtr) += here->BSIM4v7_25);

               (*(here->BSIM4v7DPgmPtr) += here->BSIM4v7_26);
               (*(here->BSIM4v7GPgmPtr) -= here->BSIM4v7_27);
               (*(here->BSIM4v7SPgmPtr) += here->BSIM4v7_28);
               (*(here->BSIM4v7BPgmPtr) += here->BSIM4v7_29);

               (*(here->BSIM4v7GPgpPtr) += here->BSIM4v7_30);
               (*(here->BSIM4v7GPdpPtr) += here->BSIM4v7_31);
               (*(here->BSIM4v7GPspPtr) += here->BSIM4v7_32);
               (*(here->BSIM4v7GPbpPtr) += here->BSIM4v7_33);
           }


            else
           {   (*(here->BSIM4v7GPgpPtr) += here->BSIM4v7_34);
               (*(here->BSIM4v7GPdpPtr) += here->BSIM4v7_35);
               (*(here->BSIM4v7GPspPtr) += here->BSIM4v7_36);
               (*(here->BSIM4v7GPbpPtr) += here->BSIM4v7_37);
           }


           if (model->BSIM4v7rdsMod)
           {   (*(here->BSIM4v7DgpPtr) += here->BSIM4v7_38);
               (*(here->BSIM4v7DspPtr) += here->BSIM4v7_39);
               (*(here->BSIM4v7DbpPtr) += here->BSIM4v7_40);
               (*(here->BSIM4v7SdpPtr) += here->BSIM4v7_41);
               (*(here->BSIM4v7SgpPtr) += here->BSIM4v7_42);
               (*(here->BSIM4v7SbpPtr) += here->BSIM4v7_43);
           }

           (*(here->BSIM4v7DPdpPtr) += here->BSIM4v7_44);
           (*(here->BSIM4v7DPdPtr) -= here->BSIM4v7_45);
           (*(here->BSIM4v7DPgpPtr) += here->BSIM4v7_46);
           (*(here->BSIM4v7DPspPtr) -= here->BSIM4v7_47);
           (*(here->BSIM4v7DPbpPtr) -= here->BSIM4v7_48);

           (*(here->BSIM4v7DdpPtr) -= here->BSIM4v7_49);
           (*(here->BSIM4v7DdPtr) += here->BSIM4v7_50);

           (*(here->BSIM4v7SPdpPtr) -= here->BSIM4v7_51);
           (*(here->BSIM4v7SPgpPtr) += here->BSIM4v7_52);
           (*(here->BSIM4v7SPspPtr) += here->BSIM4v7_53);
           (*(here->BSIM4v7SPsPtr) -= here->BSIM4v7_54);
           (*(here->BSIM4v7SPbpPtr) -= here->BSIM4v7_55);

           (*(here->BSIM4v7SspPtr) -= here->BSIM4v7_56);
           (*(here->BSIM4v7SsPtr) += here->BSIM4v7_57);

           (*(here->BSIM4v7BPdpPtr) += here->BSIM4v7_58);
           (*(here->BSIM4v7BPgpPtr) += here->BSIM4v7_59);
           (*(here->BSIM4v7BPspPtr) += here->BSIM4v7_60);
           (*(here->BSIM4v7BPbpPtr) += here->BSIM4v7_61);

           /* stamp gidl */
           (*(here->BSIM4v7DPdpPtr) += here->BSIM4v7_62);
           (*(here->BSIM4v7DPgpPtr) += here->BSIM4v7_63);
           (*(here->BSIM4v7DPspPtr) -= here->BSIM4v7_64);
           (*(here->BSIM4v7DPbpPtr) += here->BSIM4v7_65);
           (*(here->BSIM4v7BPdpPtr) -= here->BSIM4v7_66);
           (*(here->BSIM4v7BPgpPtr) -= here->BSIM4v7_67);
           (*(here->BSIM4v7BPspPtr) += here->BSIM4v7_68);
           (*(here->BSIM4v7BPbpPtr) -= here->BSIM4v7_69);
            /* stamp gisl */
           (*(here->BSIM4v7SPdpPtr) -= here->BSIM4v7_70);
           (*(here->BSIM4v7SPgpPtr) += here->BSIM4v7_71);
           (*(here->BSIM4v7SPspPtr) += here->BSIM4v7_72);
           (*(here->BSIM4v7SPbpPtr) += here->BSIM4v7_73);
           (*(here->BSIM4v7BPdpPtr) += here->BSIM4v7_74);
           (*(here->BSIM4v7BPgpPtr) -= here->BSIM4v7_75);
           (*(here->BSIM4v7BPspPtr) -= here->BSIM4v7_76);
           (*(here->BSIM4v7BPbpPtr) -= here->BSIM4v7_77);


           if (here->BSIM4v7rbodyMod)
           {   (*(here->BSIM4v7DPdbPtr) += here->BSIM4v7_78);
               (*(here->BSIM4v7SPsbPtr) -= here->BSIM4v7_79);

               (*(here->BSIM4v7DBdpPtr) += here->BSIM4v7_80);
               (*(here->BSIM4v7DBdbPtr) += here->BSIM4v7_81);
               (*(here->BSIM4v7DBbpPtr) -= here->BSIM4v7_82);
               (*(here->BSIM4v7DBbPtr) -= here->BSIM4v7_83);

               (*(here->BSIM4v7BPdbPtr) -= here->BSIM4v7_84);
               (*(here->BSIM4v7BPbPtr) -= here->BSIM4v7_85);
               (*(here->BSIM4v7BPsbPtr) -= here->BSIM4v7_86);
               (*(here->BSIM4v7BPbpPtr) += here->BSIM4v7_87);

               (*(here->BSIM4v7SBspPtr) += here->BSIM4v7_88);
               (*(here->BSIM4v7SBbpPtr) -= here->BSIM4v7_89);
               (*(here->BSIM4v7SBbPtr) -= here->BSIM4v7_90);
               (*(here->BSIM4v7SBsbPtr) += here->BSIM4v7_91);

               (*(here->BSIM4v7BdbPtr) -= here->BSIM4v7_92);
               (*(here->BSIM4v7BbpPtr) -= here->BSIM4v7_93);
               (*(here->BSIM4v7BsbPtr) -= here->BSIM4v7_94);
               (*(here->BSIM4v7BbPtr) += here->BSIM4v7_95);
           }

           if (here->BSIM4v7trnqsMod)
           {   (*(here->BSIM4v7QqPtr) += here->BSIM4v7_96);
               (*(here->BSIM4v7QgpPtr) += here->BSIM4v7_97);
               (*(here->BSIM4v7QdpPtr) += here->BSIM4v7_98);
               (*(here->BSIM4v7QspPtr) += here->BSIM4v7_99);
               (*(here->BSIM4v7QbpPtr) += here->BSIM4v7_100);

               (*(here->BSIM4v7DPqPtr) += here->BSIM4v7_101);
               (*(here->BSIM4v7SPqPtr) += here->BSIM4v7_102);
               (*(here->BSIM4v7GPqPtr) -= here->BSIM4v7_103);
           }
    }
}

#endif
