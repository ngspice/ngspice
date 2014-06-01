/**** BSIM4.6.2 Released by Wenwei Yang 07/31/2008****/
/**** BSIM4.6.5 Update ngspice 09/22/2009 ****/
/**** OpenMP support ngspice 06/28/2010 ****/
/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4ld.c of BSIM4.6.2.
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
 **********/


#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v6def.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

/*
#define MAX_EXP 2.688117142e+43
#define MIN_EXP 3.720075976e-44
#define EXP_THRESHOLD 100.0
*/
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
int BSIM4v6LoadOMP(BSIM4v6instance *here, CKTcircuit *ckt);
void BSIM4v6LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt);
#endif

int BSIM4v6polyDepletion(double phi, double ngate,double epsgate, double coxe, double Vgs, double *Vgs_eff, double *dVgs_eff_dVg);

int
BSIM4v6load(
GENmodel *inModel,
CKTcircuit *ckt)
{
#ifdef USE_OMP
    int idx;
    BSIM4v6model *model = (BSIM4v6model*)inModel;
    int good = 0;
    BSIM4v6instance **InstArray;
    InstArray = model->BSIM4v6InstanceArray;

#pragma omp parallel for
    for (idx = 0; idx < model->BSIM4v6InstCount; idx++) {
        BSIM4v6instance *here = InstArray[idx];
        int local_good = BSIM4v6LoadOMP(here, ckt);
        if (local_good)
            good = local_good;
    }

    BSIM4v6LoadRhsMat(inModel, ckt);
    
    return good;
}


int BSIM4v6LoadOMP(BSIM4v6instance *here, CKTcircuit *ckt) {
BSIM4v6model *model;
#else
BSIM4v6model *model = (BSIM4v6model*)inModel;
BSIM4v6instance *here;
#endif
double ceqgstot, dgstot_dvd, dgstot_dvg, dgstot_dvs, dgstot_dvb;
double ceqgdtot, dgdtot_dvd, dgdtot_dvg, dgdtot_dvs, dgdtot_dvb;
double gstot, gstotd, gstotg, gstots, gstotb, gspr, Rs, Rd;
double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
double vgs_eff, vgd_eff, dvgs_eff_dvg, dvgd_eff_dvg;
double dRs_dvg, dRd_dvg, dRs_dvb, dRd_dvb;
double dT0_dvg, dT1_dvb, dT3_dvg, dT3_dvb;
double vses, vdes, vdedo, delvses, delvded, delvdes;
double Isestot=0.0, cseshat=0.0, Idedtot=0.0, cdedhat=0.0;
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
double ag0, qgb, von, cbhat=0.0, VgstNVt, ExpVgst;
double ceqqb, ceqqd, ceqqg, ceqqjd=0.0, ceqqjs=0.0, ceq, geq;
double cdrain, cdhat=0.0, ceqdrn, ceqbd, ceqbs, ceqjd, ceqjs, gjbd, gjbs;
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
double Igstot=0.0, cgshat=0.0, Igdtot=0.0, cgdhat=0.0, Igbtot=0.0, cgbhat=0.0;
double Vgs_eff, Vfb=0.0, Vth_NarrowW;
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
double Cgg1, Cgb1, Cgd1, Cbg1, Cbb1, Cbd1, Qac0, Qsub0;
double dQac0_dVg, dQac0_dVb, dQsub0_dVg, dQsub0_dVd, dQsub0_dVb;
double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;
double Igisl, Ggislg, Ggislb, Ggisls;
double Nvtmrss, Nvtmrssws, Nvtmrsswgs;
double Nvtmrsd, Nvtmrsswd, Nvtmrsswgd;

double vs, Fsevl, dvs_dVg, dvs_dVd, dvs_dVb, dFsevl_dVg, dFsevl_dVd, dFsevl_dVb;
double vgdx, vgsx, epssub, toxe, epsrox;
struct bsim4v6SizeDependParam *pParam;
int ByPass, ChargeComputationNeeded, error, Check, Check1, Check2;

double m;

#ifdef USE_OMP
model = here->BSIM4v6modPtr;
#endif

ScalingFactor = 1.0e-9;
ChargeComputationNeeded =  
                 ((ckt->CKTmode & (MODEAC | MODETRAN | MODEINITSMSIG)) ||
                 ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)))
                 ? 1 : 0;

#ifndef USE_OMP
for (; model != NULL; model = model->BSIM4v6nextModel)
{    for (here = model->BSIM4v6instances; here != NULL; 
          here = here->BSIM4v6nextInstance)
     {
#endif

          Check = Check1 = Check2 = 1;
          ByPass = 0;
          pParam = here->pParam;

          if ((ckt->CKTmode & MODEINITSMSIG))
          {   vds = *(ckt->CKTstate0 + here->BSIM4v6vds);
              vgs = *(ckt->CKTstate0 + here->BSIM4v6vgs);
              vbs = *(ckt->CKTstate0 + here->BSIM4v6vbs);
              vges = *(ckt->CKTstate0 + here->BSIM4v6vges);
              vgms = *(ckt->CKTstate0 + here->BSIM4v6vgms);
              vdbs = *(ckt->CKTstate0 + here->BSIM4v6vdbs);
              vsbs = *(ckt->CKTstate0 + here->BSIM4v6vsbs);
              vses = *(ckt->CKTstate0 + here->BSIM4v6vses);
              vdes = *(ckt->CKTstate0 + here->BSIM4v6vdes);

              qdef = *(ckt->CKTstate0 + here->BSIM4v6qdef);
          }
          else if ((ckt->CKTmode & MODEINITTRAN))
          {   vds = *(ckt->CKTstate1 + here->BSIM4v6vds);
              vgs = *(ckt->CKTstate1 + here->BSIM4v6vgs);
              vbs = *(ckt->CKTstate1 + here->BSIM4v6vbs);
              vges = *(ckt->CKTstate1 + here->BSIM4v6vges);
              vgms = *(ckt->CKTstate1 + here->BSIM4v6vgms);
              vdbs = *(ckt->CKTstate1 + here->BSIM4v6vdbs);
              vsbs = *(ckt->CKTstate1 + here->BSIM4v6vsbs);
              vses = *(ckt->CKTstate1 + here->BSIM4v6vses);
              vdes = *(ckt->CKTstate1 + here->BSIM4v6vdes);

              qdef = *(ckt->CKTstate1 + here->BSIM4v6qdef);
          }
          else if ((ckt->CKTmode & MODEINITJCT) && !here->BSIM4v6off)
          {   vds = model->BSIM4v6type * here->BSIM4v6icVDS;
              vgs = vges = vgms = model->BSIM4v6type * here->BSIM4v6icVGS;
              vbs = vdbs = vsbs = model->BSIM4v6type * here->BSIM4v6icVBS;
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
                  vgs = vges = vgms = model->BSIM4v6type 
                                    * here->BSIM4v6vth0 + 0.1;
                  vbs = vdbs = vsbs = 0.0;
              }
          }
          else if ((ckt->CKTmode & (MODEINITJCT | MODEINITFIX)) && 
                  (here->BSIM4v6off)) 
          {   vds = vgs = vbs = vges = vgms = 0.0;
              vdbs = vsbs = vdes = vses = qdef = 0.0;
          }
          else
          {
#ifndef PREDICTOR
               if ((ckt->CKTmode & MODEINITPRED))
               {   xfact = ckt->CKTdelta / ckt->CKTdeltaOld[1];
                   *(ckt->CKTstate0 + here->BSIM4v6vds) = 
                         *(ckt->CKTstate1 + here->BSIM4v6vds);
                   vds = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vds))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vds)));
                   *(ckt->CKTstate0 + here->BSIM4v6vgs) = 
                         *(ckt->CKTstate1 + here->BSIM4v6vgs);
                   vgs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vgs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vgs)));
                   *(ckt->CKTstate0 + here->BSIM4v6vges) =
                         *(ckt->CKTstate1 + here->BSIM4v6vges);
                   vges = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vges))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vges)));
                   *(ckt->CKTstate0 + here->BSIM4v6vgms) =
                         *(ckt->CKTstate1 + here->BSIM4v6vgms);
                   vgms = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vgms))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vgms)));
                   *(ckt->CKTstate0 + here->BSIM4v6vbs) = 
                         *(ckt->CKTstate1 + here->BSIM4v6vbs);
                   vbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vbs)));
                   *(ckt->CKTstate0 + here->BSIM4v6vbd) = 
                         *(ckt->CKTstate0 + here->BSIM4v6vbs)
                         - *(ckt->CKTstate0 + here->BSIM4v6vds);
                   *(ckt->CKTstate0 + here->BSIM4v6vdbs) =
                         *(ckt->CKTstate1 + here->BSIM4v6vdbs);
                   vdbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vdbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vdbs)));
                   *(ckt->CKTstate0 + here->BSIM4v6vdbd) =
                         *(ckt->CKTstate0 + here->BSIM4v6vdbs)
                         - *(ckt->CKTstate0 + here->BSIM4v6vds);
                   *(ckt->CKTstate0 + here->BSIM4v6vsbs) =
                         *(ckt->CKTstate1 + here->BSIM4v6vsbs);
                   vsbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vsbs))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vsbs)));
                   *(ckt->CKTstate0 + here->BSIM4v6vses) =
                         *(ckt->CKTstate1 + here->BSIM4v6vses);
                   vses = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vses))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vses)));
                   *(ckt->CKTstate0 + here->BSIM4v6vdes) =
                         *(ckt->CKTstate1 + here->BSIM4v6vdes);
                   vdes = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6vdes))
                         - (xfact * (*(ckt->CKTstate2 + here->BSIM4v6vdes)));

                   *(ckt->CKTstate0 + here->BSIM4v6qdef) =
                         *(ckt->CKTstate1 + here->BSIM4v6qdef);
                   qdef = (1.0 + xfact)* (*(ckt->CKTstate1 + here->BSIM4v6qdef))
                        -(xfact * (*(ckt->CKTstate2 + here->BSIM4v6qdef)));
               }
               else
               {
#endif /* PREDICTOR */
                   vds = model->BSIM4v6type
                       * (*(ckt->CKTrhsOld + here->BSIM4v6dNodePrime)
                       - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vgs = model->BSIM4v6type
                       * (*(ckt->CKTrhsOld + here->BSIM4v6gNodePrime) 
                       - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vbs = model->BSIM4v6type
                       * (*(ckt->CKTrhsOld + here->BSIM4v6bNodePrime)
                       - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vges = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6gNodeExt)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vgms = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6gNodeMid)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vdbs = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6dbNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vsbs = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6sbNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vses = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6sNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   vdes = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6dNode)
                        - *(ckt->CKTrhsOld + here->BSIM4v6sNodePrime));
                   qdef = model->BSIM4v6type
                        * (*(ckt->CKTrhsOld + here->BSIM4v6qNode));
#ifndef PREDICTOR
               }
#endif /* PREDICTOR */

               vgdo = *(ckt->CKTstate0 + here->BSIM4v6vgs)
                    - *(ckt->CKTstate0 + here->BSIM4v6vds);
               vgedo = *(ckt->CKTstate0 + here->BSIM4v6vges)
                     - *(ckt->CKTstate0 + here->BSIM4v6vds);
               vgmdo = *(ckt->CKTstate0 + here->BSIM4v6vgms)
                     - *(ckt->CKTstate0 + here->BSIM4v6vds);

               vbd = vbs - vds;
               vdbd = vdbs - vds;
               vgd = vgs - vds;
                 vged = vges - vds;
                 vgmd = vgms - vds;

               delvbd = vbd - *(ckt->CKTstate0 + here->BSIM4v6vbd);
               delvdbd = vdbd - *(ckt->CKTstate0 + here->BSIM4v6vdbd);
               delvgd = vgd - vgdo;
               delvged = vged - vgedo;
               delvgmd = vgmd - vgmdo;

               delvds = vds - *(ckt->CKTstate0 + here->BSIM4v6vds);
               delvgs = vgs - *(ckt->CKTstate0 + here->BSIM4v6vgs);
               delvges = vges - *(ckt->CKTstate0 + here->BSIM4v6vges);
               delvgms = vgms - *(ckt->CKTstate0 + here->BSIM4v6vgms);
               delvbs = vbs - *(ckt->CKTstate0 + here->BSIM4v6vbs);
               delvdbs = vdbs - *(ckt->CKTstate0 + here->BSIM4v6vdbs);
               delvsbs = vsbs - *(ckt->CKTstate0 + here->BSIM4v6vsbs);

               delvses = vses - (*(ckt->CKTstate0 + here->BSIM4v6vses));
               vdedo = *(ckt->CKTstate0 + here->BSIM4v6vdes)
                     - *(ckt->CKTstate0 + here->BSIM4v6vds);
               delvdes = vdes - *(ckt->CKTstate0 + here->BSIM4v6vdes);
               delvded = vdes - vds - vdedo;

               delvbd_jct = (!here->BSIM4v6rbodyMod) ? delvbd : delvdbd;
               delvbs_jct = (!here->BSIM4v6rbodyMod) ? delvbs : delvsbs;
               if (here->BSIM4v6mode >= 0)
               {   Idtot = here->BSIM4v6cd + here->BSIM4v6csub - here->BSIM4v6cbd
                         + here->BSIM4v6Igidl;
                   cdhat = Idtot - here->BSIM4v6gbd * delvbd_jct
                         + (here->BSIM4v6gmbs + here->BSIM4v6gbbs + here->BSIM4v6ggidlb) * delvbs
                         + (here->BSIM4v6gm + here->BSIM4v6gbgs + here->BSIM4v6ggidlg) * delvgs 
                         + (here->BSIM4v6gds + here->BSIM4v6gbds + here->BSIM4v6ggidld) * delvds;
                   Ibtot = here->BSIM4v6cbs + here->BSIM4v6cbd 
                         - here->BSIM4v6Igidl - here->BSIM4v6Igisl - here->BSIM4v6csub;
                   cbhat = Ibtot + here->BSIM4v6gbd * delvbd_jct
                         + here->BSIM4v6gbs * delvbs_jct - (here->BSIM4v6gbbs + here->BSIM4v6ggidlb)
                         * delvbs - (here->BSIM4v6gbgs + here->BSIM4v6ggidlg) * delvgs
                         - (here->BSIM4v6gbds + here->BSIM4v6ggidld - here->BSIM4v6ggisls) * delvds 
                         - here->BSIM4v6ggislg * delvgd - here->BSIM4v6ggislb* delvbd;

                   Igstot = here->BSIM4v6Igs + here->BSIM4v6Igcs;
                   cgshat = Igstot + (here->BSIM4v6gIgsg + here->BSIM4v6gIgcsg) * delvgs
                          + here->BSIM4v6gIgcsd * delvds + here->BSIM4v6gIgcsb * delvbs;

                   Igdtot = here->BSIM4v6Igd + here->BSIM4v6Igcd;
                   cgdhat = Igdtot + here->BSIM4v6gIgdg * delvgd + here->BSIM4v6gIgcdg * delvgs
                          + here->BSIM4v6gIgcdd * delvds + here->BSIM4v6gIgcdb * delvbs;

                    Igbtot = here->BSIM4v6Igb;
                   cgbhat = here->BSIM4v6Igb + here->BSIM4v6gIgbg * delvgs + here->BSIM4v6gIgbd
                          * delvds + here->BSIM4v6gIgbb * delvbs;
               }
               else
               {   Idtot = here->BSIM4v6cd + here->BSIM4v6cbd - here->BSIM4v6Igidl; /* bugfix */
                   cdhat = Idtot + here->BSIM4v6gbd * delvbd_jct + here->BSIM4v6gmbs 
                         * delvbd + here->BSIM4v6gm * delvgd 
                         - (here->BSIM4v6gds + here->BSIM4v6ggidls) * delvds 
                         - here->BSIM4v6ggidlg * delvgs - here->BSIM4v6ggidlb * delvbs;
                   Ibtot = here->BSIM4v6cbs + here->BSIM4v6cbd 
                         - here->BSIM4v6Igidl - here->BSIM4v6Igisl - here->BSIM4v6csub;
                   cbhat = Ibtot + here->BSIM4v6gbs * delvbs_jct + here->BSIM4v6gbd 
                         * delvbd_jct - (here->BSIM4v6gbbs + here->BSIM4v6ggislb) * delvbd
                         - (here->BSIM4v6gbgs + here->BSIM4v6ggislg) * delvgd
                         + (here->BSIM4v6gbds + here->BSIM4v6ggisld - here->BSIM4v6ggidls) * delvds
                         - here->BSIM4v6ggidlg * delvgs - here->BSIM4v6ggidlb * delvbs; 

                   Igstot = here->BSIM4v6Igs + here->BSIM4v6Igcd;
                   cgshat = Igstot + here->BSIM4v6gIgsg * delvgs + here->BSIM4v6gIgcdg * delvgd
                          - here->BSIM4v6gIgcdd * delvds + here->BSIM4v6gIgcdb * delvbd;

                   Igdtot = here->BSIM4v6Igd + here->BSIM4v6Igcs;
                   cgdhat = Igdtot + (here->BSIM4v6gIgdg + here->BSIM4v6gIgcsg) * delvgd
                          - here->BSIM4v6gIgcsd * delvds + here->BSIM4v6gIgcsb * delvbd;

                   Igbtot = here->BSIM4v6Igb;
                   cgbhat = here->BSIM4v6Igb + here->BSIM4v6gIgbg * delvgd - here->BSIM4v6gIgbd
                          * delvds + here->BSIM4v6gIgbb * delvbd;
               }

               Isestot = here->BSIM4v6gstot * (*(ckt->CKTstate0 + here->BSIM4v6vses));
               cseshat = Isestot + here->BSIM4v6gstot * delvses
                       + here->BSIM4v6gstotd * delvds + here->BSIM4v6gstotg * delvgs
                       + here->BSIM4v6gstotb * delvbs;

               Idedtot = here->BSIM4v6gdtot * vdedo;
               cdedhat = Idedtot + here->BSIM4v6gdtot * delvded
                       + here->BSIM4v6gdtotd * delvds + here->BSIM4v6gdtotg * delvgs
                       + here->BSIM4v6gdtotb * delvbs;


#ifndef NOBYPASS
               /* Following should be one IF statement, but some C compilers 
                * can't handle that all at once, so we split it into several
                * successive IF's */

               if ((!(ckt->CKTmode & MODEINITPRED)) && (ckt->CKTbypass))
               if ((fabs(delvds) < (ckt->CKTreltol * MAX(fabs(vds),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v6vds))) + ckt->CKTvoltTol)))
               if ((fabs(delvgs) < (ckt->CKTreltol * MAX(fabs(vgs),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v6vgs))) + ckt->CKTvoltTol)))
               if ((fabs(delvbs) < (ckt->CKTreltol * MAX(fabs(vbs),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v6vbs))) + ckt->CKTvoltTol)))
               if ((fabs(delvbd) < (ckt->CKTreltol * MAX(fabs(vbd),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v6vbd))) + ckt->CKTvoltTol)))
               if ((here->BSIM4v6rgateMod == 0) || (here->BSIM4v6rgateMod == 1) 
                         || (fabs(delvges) < (ckt->CKTreltol * MAX(fabs(vges),
                   fabs(*(ckt->CKTstate0 + here->BSIM4v6vges))) + ckt->CKTvoltTol)))
               if ((here->BSIM4v6rgateMod != 3) || (fabs(delvgms) < (ckt->CKTreltol
                   * MAX(fabs(vgms), fabs(*(ckt->CKTstate0 + here->BSIM4v6vgms)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v6rbodyMod) || (fabs(delvdbs) < (ckt->CKTreltol
                   * MAX(fabs(vdbs), fabs(*(ckt->CKTstate0 + here->BSIM4v6vdbs)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v6rbodyMod) || (fabs(delvdbd) < (ckt->CKTreltol
                   * MAX(fabs(vdbd), fabs(*(ckt->CKTstate0 + here->BSIM4v6vdbd)))
                   + ckt->CKTvoltTol)))
               if ((!here->BSIM4v6rbodyMod) || (fabs(delvsbs) < (ckt->CKTreltol
                   * MAX(fabs(vsbs), fabs(*(ckt->CKTstate0 + here->BSIM4v6vsbs)))
                   + ckt->CKTvoltTol)))
               if ((!model->BSIM4v6rdsMod) || (fabs(delvses) < (ckt->CKTreltol
                   * MAX(fabs(vses), fabs(*(ckt->CKTstate0 + here->BSIM4v6vses)))
                   + ckt->CKTvoltTol)))
               if ((!model->BSIM4v6rdsMod) || (fabs(delvdes) < (ckt->CKTreltol
                   * MAX(fabs(vdes), fabs(*(ckt->CKTstate0 + here->BSIM4v6vdes)))
                   + ckt->CKTvoltTol)))
               if ((fabs(cdhat - Idtot) < ckt->CKTreltol
                   * MAX(fabs(cdhat), fabs(Idtot)) + ckt->CKTabstol))
               if ((fabs(cbhat - Ibtot) < ckt->CKTreltol
                   * MAX(fabs(cbhat), fabs(Ibtot)) + ckt->CKTabstol))
               if ((!model->BSIM4v6igcMod) || ((fabs(cgshat - Igstot) < ckt->CKTreltol
                   * MAX(fabs(cgshat), fabs(Igstot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v6igcMod) || ((fabs(cgdhat - Igdtot) < ckt->CKTreltol
                   * MAX(fabs(cgdhat), fabs(Igdtot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v6igbMod) || ((fabs(cgbhat - Igbtot) < ckt->CKTreltol
                   * MAX(fabs(cgbhat), fabs(Igbtot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v6rdsMod) || ((fabs(cseshat - Isestot) < ckt->CKTreltol
                   * MAX(fabs(cseshat), fabs(Isestot)) + ckt->CKTabstol)))
               if ((!model->BSIM4v6rdsMod) || ((fabs(cdedhat - Idedtot) < ckt->CKTreltol
                   * MAX(fabs(cdedhat), fabs(Idedtot)) + ckt->CKTabstol)))
               {   vds = *(ckt->CKTstate0 + here->BSIM4v6vds);
                   vgs = *(ckt->CKTstate0 + here->BSIM4v6vgs);
                   vbs = *(ckt->CKTstate0 + here->BSIM4v6vbs);
                   vges = *(ckt->CKTstate0 + here->BSIM4v6vges);
                   vgms = *(ckt->CKTstate0 + here->BSIM4v6vgms);

                   vbd = *(ckt->CKTstate0 + here->BSIM4v6vbd);
                   vdbs = *(ckt->CKTstate0 + here->BSIM4v6vdbs);
                   vdbd = *(ckt->CKTstate0 + here->BSIM4v6vdbd);
                   vsbs = *(ckt->CKTstate0 + here->BSIM4v6vsbs);
                   vses = *(ckt->CKTstate0 + here->BSIM4v6vses);
                   vdes = *(ckt->CKTstate0 + here->BSIM4v6vdes);

                   vgd = vgs - vds;
                   vgb = vgs - vbs;
                   vged = vges - vds;
                   vgmd = vgms - vds;
                   vgmb = vgms - vbs;

                   vbs_jct = (!here->BSIM4v6rbodyMod) ? vbs : vsbs;
                   vbd_jct = (!here->BSIM4v6rbodyMod) ? vbd : vdbd;

/*** qdef should not be kept fixed even if vgs, vds & vbs has converged 
****               qdef = *(ckt->CKTstate0 + here->BSIM4v6qdef);  
***/
                   cdrain = here->BSIM4v6cd;

                   if ((ckt->CKTmode & (MODETRAN | MODEAC)) || 
                       ((ckt->CKTmode & MODETRANOP) && 
                       (ckt->CKTmode & MODEUIC)))
                   {   ByPass = 1;

                       qgate = here->BSIM4v6qgate;
                       qbulk = here->BSIM4v6qbulk;
                       qdrn = here->BSIM4v6qdrn;
                       cgdo = here->BSIM4v6cgdo;
                       qgdo = here->BSIM4v6qgdo;
                       cgso = here->BSIM4v6cgso;
                       qgso = here->BSIM4v6qgso;

                       goto line755;
                    }
                    else
                       goto line850;
               }
#endif /*NOBYPASS*/

               von = here->BSIM4v6von;
               if (*(ckt->CKTstate0 + here->BSIM4v6vds) >= 0.0)
               {   vgs = DEVfetlim(vgs, *(ckt->CKTstate0 + here->BSIM4v6vgs), von);
                   vds = vgs - vgd;
                   vds = DEVlimvds(vds, *(ckt->CKTstate0 + here->BSIM4v6vds));
                   vgd = vgs - vds;
                   if (here->BSIM4v6rgateMod == 3)
                   {   vges = DEVfetlim(vges, *(ckt->CKTstate0 + here->BSIM4v6vges), von);
                       vgms = DEVfetlim(vgms, *(ckt->CKTstate0 + here->BSIM4v6vgms), von);
                       vged = vges - vds;
                       vgmd = vgms - vds;
                   }
                   else if ((here->BSIM4v6rgateMod == 1) || (here->BSIM4v6rgateMod == 2))
                   {   vges = DEVfetlim(vges, *(ckt->CKTstate0 + here->BSIM4v6vges), von);
                       vged = vges - vds;
                   }

                   if (model->BSIM4v6rdsMod)
                   {   vdes = DEVlimvds(vdes, *(ckt->CKTstate0 + here->BSIM4v6vdes));
                       vses = -DEVlimvds(-vses, -(*(ckt->CKTstate0 + here->BSIM4v6vses)));
                   }

               }
               else
               {   vgd = DEVfetlim(vgd, vgdo, von);
                   vds = vgs - vgd;
                   vds = -DEVlimvds(-vds, -(*(ckt->CKTstate0 + here->BSIM4v6vds)));
                   vgs = vgd + vds;

                   if (here->BSIM4v6rgateMod == 3)
                   {   vged = DEVfetlim(vged, vgedo, von);
                       vges = vged + vds;
                       vgmd = DEVfetlim(vgmd, vgmdo, von);
                       vgms = vgmd + vds;
                   }
                   if ((here->BSIM4v6rgateMod == 1) || (here->BSIM4v6rgateMod == 2))
                   {   vged = DEVfetlim(vged, vgedo, von);
                       vges = vged + vds;
                   }

                   if (model->BSIM4v6rdsMod)
                   {   vdes = -DEVlimvds(-vdes, -(*(ckt->CKTstate0 + here->BSIM4v6vdes)));
                       vses = DEVlimvds(vses, *(ckt->CKTstate0 + here->BSIM4v6vses));
                   }
               }

               if (vds >= 0.0)
               {   vbs = DEVpnjlim(vbs, *(ckt->CKTstate0 + here->BSIM4v6vbs),
                                   CONSTvt0, model->BSIM4v6vcrit, &Check);
                   vbd = vbs - vds;
                   if (here->BSIM4v6rbodyMod)
                   {   vdbs = DEVpnjlim(vdbs, *(ckt->CKTstate0 + here->BSIM4v6vdbs),
                                        CONSTvt0, model->BSIM4v6vcrit, &Check1);
                       vdbd = vdbs - vds;
                       vsbs = DEVpnjlim(vsbs, *(ckt->CKTstate0 + here->BSIM4v6vsbs),
                                        CONSTvt0, model->BSIM4v6vcrit, &Check2);
                       if ((Check1 == 0) && (Check2 == 0))
                           Check = 0;
                       else 
                           Check = 1;
                   }
               }
               else
               {   vbd = DEVpnjlim(vbd, *(ckt->CKTstate0 + here->BSIM4v6vbd),
                                   CONSTvt0, model->BSIM4v6vcrit, &Check); 
                   vbs = vbd + vds;
                   if (here->BSIM4v6rbodyMod)
                   {   vdbd = DEVpnjlim(vdbd, *(ckt->CKTstate0 + here->BSIM4v6vdbd),
                                        CONSTvt0, model->BSIM4v6vcrit, &Check1);
                       vdbs = vdbd + vds;
                       vsbdo = *(ckt->CKTstate0 + here->BSIM4v6vsbs)
                             - *(ckt->CKTstate0 + here->BSIM4v6vds);
                       vsbd = vsbs - vds;
                       vsbd = DEVpnjlim(vsbd, vsbdo, CONSTvt0, model->BSIM4v6vcrit, &Check2);
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

          vbs_jct = (!here->BSIM4v6rbodyMod) ? vbs : vsbs;
          vbd_jct = (!here->BSIM4v6rbodyMod) ? vbd : vdbd;

          /* Source/drain junction diode DC model begins */
          Nvtms = model->BSIM4v6vtm * model->BSIM4v6SjctEmissionCoeff;
          if ((here->BSIM4v6Aseff <= 0.0) && (here->BSIM4v6Pseff <= 0.0))
          {   SourceSatCurrent = 0.0;
          }
          else
          {   SourceSatCurrent = here->BSIM4v6Aseff * model->BSIM4v6SjctTempSatCurDensity
                               + here->BSIM4v6Pseff * model->BSIM4v6SjctSidewallTempSatCurDensity
                               + pParam->BSIM4v6weffCJ * here->BSIM4v6nf
                               * model->BSIM4v6SjctGateSidewallTempSatCurDensity;
          }

          if (SourceSatCurrent <= 0.0)
          {   here->BSIM4v6gbs = ckt->CKTgmin;
              here->BSIM4v6cbs = here->BSIM4v6gbs * vbs_jct;
          }
          else
          {   switch(model->BSIM4v6dioMod)
              {   case 0:
                      evbs = exp(vbs_jct / Nvtms);
                      T1 = model->BSIM4v6xjbvs * exp(-(model->BSIM4v6bvs + vbs_jct) / Nvtms);
                      /* WDLiu: Magic T1 in this form; different from BSIM4v6 beta. */
                      here->BSIM4v6gbs = SourceSatCurrent * (evbs + T1) / Nvtms + ckt->CKTgmin;
                      here->BSIM4v6cbs = SourceSatCurrent * (evbs + here->BSIM4v6XExpBVS
                                     - T1 - 1.0) + ckt->CKTgmin * vbs_jct;
                      break;
                  case 1:
                      T2 = vbs_jct / Nvtms;
                      if (T2 < -EXP_THRESHOLD)
                      {   here->BSIM4v6gbs = ckt->CKTgmin;
                          here->BSIM4v6cbs = SourceSatCurrent * (MIN_EXP - 1.0)
                                         + ckt->CKTgmin * vbs_jct;
                      }
                      else if (vbs_jct <= here->BSIM4v6vjsmFwd)
                      {   evbs = exp(T2);
                          here->BSIM4v6gbs = SourceSatCurrent * evbs / Nvtms + ckt->CKTgmin;
                          here->BSIM4v6cbs = SourceSatCurrent * (evbs - 1.0)
                                         + ckt->CKTgmin * vbs_jct;
                      }
                      else
                      {   T0 = here->BSIM4v6IVjsmFwd / Nvtms;
                          here->BSIM4v6gbs = T0 + ckt->CKTgmin;
                          here->BSIM4v6cbs = here->BSIM4v6IVjsmFwd - SourceSatCurrent + T0 
                                         * (vbs_jct - here->BSIM4v6vjsmFwd) + ckt->CKTgmin * vbs_jct;
                      }        
                      break;
                  case 2:
                      if (vbs_jct < here->BSIM4v6vjsmRev)
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
                          T2 = here->BSIM4v6IVjsmRev + here->BSIM4v6SslpRev
                             * (vbs_jct - here->BSIM4v6vjsmRev);
                          here->BSIM4v6gbs = devbs_dvb * T2 + T1 * here->BSIM4v6SslpRev + ckt->CKTgmin;
                          here->BSIM4v6cbs = T1 * T2 + ckt->CKTgmin * vbs_jct;
                      }         
                      else if (vbs_jct <= here->BSIM4v6vjsmFwd)
                      {   T0 = vbs_jct / Nvtms;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbs = MIN_EXP;
                               devbs_dvb = 0.0;
                          }
                          else
                          {    evbs = exp(T0);
                               devbs_dvb = evbs / Nvtms;
                          }

                          T1 = (model->BSIM4v6bvs + vbs_jct) / Nvtms;
                          if (T1 > EXP_THRESHOLD)
                          {   T2 = MIN_EXP;
                              T3 = 0.0;
                          }
                          else
                          {   T2 = exp(-T1);
                              T3 = -T2 /Nvtms;
                          }
                          here->BSIM4v6gbs = SourceSatCurrent * (devbs_dvb - model->BSIM4v6xjbvs * T3)
                                         + ckt->CKTgmin;
                          here->BSIM4v6cbs = SourceSatCurrent * (evbs + here->BSIM4v6XExpBVS - 1.0
                                         - model->BSIM4v6xjbvs * T2) + ckt->CKTgmin * vbs_jct;
                      }
                      else
                      {   here->BSIM4v6gbs = here->BSIM4v6SslpFwd + ckt->CKTgmin;
                          here->BSIM4v6cbs = here->BSIM4v6IVjsmFwd + here->BSIM4v6SslpFwd * (vbs_jct
                                         - here->BSIM4v6vjsmFwd) + ckt->CKTgmin * vbs_jct;
                      }
                      break;
                  default: break;
              }
          }

          Nvtmd = model->BSIM4v6vtm * model->BSIM4v6DjctEmissionCoeff;
          if ((here->BSIM4v6Adeff <= 0.0) && (here->BSIM4v6Pdeff <= 0.0))
          {   DrainSatCurrent = 0.0;
          }
          else
          {   DrainSatCurrent = here->BSIM4v6Adeff * model->BSIM4v6DjctTempSatCurDensity
                              + here->BSIM4v6Pdeff * model->BSIM4v6DjctSidewallTempSatCurDensity
                              + pParam->BSIM4v6weffCJ * here->BSIM4v6nf
                              * model->BSIM4v6DjctGateSidewallTempSatCurDensity;
          }

          if (DrainSatCurrent <= 0.0)
          {   here->BSIM4v6gbd = ckt->CKTgmin;
              here->BSIM4v6cbd = here->BSIM4v6gbd * vbd_jct;
          }
          else
          {   switch(model->BSIM4v6dioMod)
              {   case 0:
                      evbd = exp(vbd_jct / Nvtmd);
                      T1 = model->BSIM4v6xjbvd * exp(-(model->BSIM4v6bvd + vbd_jct) / Nvtmd);
                      /* WDLiu: Magic T1 in this form; different from BSIM4v6 beta. */
                      here->BSIM4v6gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + ckt->CKTgmin;
                      here->BSIM4v6cbd = DrainSatCurrent * (evbd + here->BSIM4v6XExpBVD
                                     - T1 - 1.0) + ckt->CKTgmin * vbd_jct;
                      break;
                  case 1:
                      T2 = vbd_jct / Nvtmd;
                      if (T2 < -EXP_THRESHOLD)
                      {   here->BSIM4v6gbd = ckt->CKTgmin;
                          here->BSIM4v6cbd = DrainSatCurrent * (MIN_EXP - 1.0)
                                         + ckt->CKTgmin * vbd_jct;
                      }
                      else if (vbd_jct <= here->BSIM4v6vjdmFwd)
                      {   evbd = exp(T2);
                          here->BSIM4v6gbd = DrainSatCurrent * evbd / Nvtmd + ckt->CKTgmin;
                          here->BSIM4v6cbd = DrainSatCurrent * (evbd - 1.0)
                                         + ckt->CKTgmin * vbd_jct;
                      }
                      else
                      {   T0 = here->BSIM4v6IVjdmFwd / Nvtmd;
                          here->BSIM4v6gbd = T0 + ckt->CKTgmin;
                          here->BSIM4v6cbd = here->BSIM4v6IVjdmFwd - DrainSatCurrent + T0
                                         * (vbd_jct - here->BSIM4v6vjdmFwd) + ckt->CKTgmin * vbd_jct;
                      }
                      break;
                  case 2:
                      if (vbd_jct < here->BSIM4v6vjdmRev)
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
                          T2 = here->BSIM4v6IVjdmRev + here->BSIM4v6DslpRev
                             * (vbd_jct - here->BSIM4v6vjdmRev);
                          here->BSIM4v6gbd = devbd_dvb * T2 + T1 * here->BSIM4v6DslpRev + ckt->CKTgmin;
                          here->BSIM4v6cbd = T1 * T2 + ckt->CKTgmin * vbd_jct;
                      }
                      else if (vbd_jct <= here->BSIM4v6vjdmFwd)
                      {   T0 = vbd_jct / Nvtmd;
                          if (T0 < -EXP_THRESHOLD)
                          {    evbd = MIN_EXP;
                               devbd_dvb = 0.0;
                          }
                          else
                          {    evbd = exp(T0);
                               devbd_dvb = evbd / Nvtmd;
                          }

                          T1 = (model->BSIM4v6bvd + vbd_jct) / Nvtmd;
                          if (T1 > EXP_THRESHOLD)
                          {   T2 = MIN_EXP;
                              T3 = 0.0;
                          }
                          else
                          {   T2 = exp(-T1);
                              T3 = -T2 /Nvtmd;
                          }     
                          here->BSIM4v6gbd = DrainSatCurrent * (devbd_dvb - model->BSIM4v6xjbvd * T3)
                                         + ckt->CKTgmin;
                          here->BSIM4v6cbd = DrainSatCurrent * (evbd + here->BSIM4v6XExpBVD - 1.0
                                         - model->BSIM4v6xjbvd * T2) + ckt->CKTgmin * vbd_jct;
                      }
                      else
                      {   here->BSIM4v6gbd = here->BSIM4v6DslpFwd + ckt->CKTgmin;
                          here->BSIM4v6cbd = here->BSIM4v6IVjdmFwd + here->BSIM4v6DslpFwd * (vbd_jct
                                         - here->BSIM4v6vjdmFwd) + ckt->CKTgmin * vbd_jct;
                      }
                      break;
                  default: break;
              }
          } 

           /* trap-assisted tunneling and recombination current for reverse bias  */
          Nvtmrssws = model->BSIM4v6vtm0 * model->BSIM4v6njtsswstemp;
          Nvtmrsswgs = model->BSIM4v6vtm0 * model->BSIM4v6njtsswgstemp;
          Nvtmrss = model->BSIM4v6vtm0 * model->BSIM4v6njtsstemp;
          Nvtmrsswd = model->BSIM4v6vtm0 * model->BSIM4v6njtsswdtemp;
          Nvtmrsswgd = model->BSIM4v6vtm0 * model->BSIM4v6njtsswgdtemp;
          Nvtmrsd = model->BSIM4v6vtm0 * model->BSIM4v6njtsdtemp;

        if ((model->BSIM4v6vtss - vbs_jct) < (model->BSIM4v6vtss * 1e-3))
        { T9 = 1.0e3; 
          T0 = - vbs_jct / Nvtmrss * T9;
          DEXP(T0, T1, T10);
          dT1_dVb = T10 / Nvtmrss * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtss - vbs_jct);
          T0 = -vbs_jct / Nvtmrss * model->BSIM4v6vtss * T9;
          dT0_dVb = model->BSIM4v6vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T1, T10);
          dT1_dVb = T10 * dT0_dVb;
        }

       if ((model->BSIM4v6vtsd - vbd_jct) < (model->BSIM4v6vtsd * 1e-3) )
        { T9 = 1.0e3;
          T0 = -vbd_jct / Nvtmrsd * T9;
          DEXP(T0, T2, T10);
          dT2_dVb = T10 / Nvtmrsd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtsd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsd * model->BSIM4v6vtsd * T9;
          dT0_dVb = model->BSIM4v6vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T2, T10);
          dT2_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v6vtssws - vbs_jct) < (model->BSIM4v6vtssws * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbs_jct / Nvtmrssws * T9;
          DEXP(T0, T3, T10);
          dT3_dVb = T10 / Nvtmrssws * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtssws - vbs_jct);
          T0 = -vbs_jct / Nvtmrssws * model->BSIM4v6vtssws * T9;
          dT0_dVb = model->BSIM4v6vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T3, T10);
          dT3_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v6vtsswd - vbd_jct) < (model->BSIM4v6vtsswd * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbd_jct / Nvtmrsswd * T9;
          DEXP(T0, T4, T10);
          dT4_dVb = T10 / Nvtmrsswd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtsswd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsswd * model->BSIM4v6vtsswd * T9;
          dT0_dVb = model->BSIM4v6vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T4, T10);
          dT4_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v6vtsswgs - vbs_jct) < (model->BSIM4v6vtsswgs * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbs_jct / Nvtmrsswgs * T9;
          DEXP(T0, T5, T10);
          dT5_dVb = T10 / Nvtmrsswgs * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtsswgs - vbs_jct);
          T0 = -vbs_jct / Nvtmrsswgs * model->BSIM4v6vtsswgs * T9;
          dT0_dVb = model->BSIM4v6vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9) ;
          DEXP(T0, T5, T10);
          dT5_dVb = T10 * dT0_dVb;
        }

        if ((model->BSIM4v6vtsswgd - vbd_jct) < (model->BSIM4v6vtsswgd * 1e-3) )
        { T9 = 1.0e3; 
          T0 = -vbd_jct / Nvtmrsswgd * T9;
          DEXP(T0, T6, T10);
          dT6_dVb = T10 / Nvtmrsswgd * T9; 
        } else {
          T9 = 1.0 / (model->BSIM4v6vtsswgd - vbd_jct);
          T0 = -vbd_jct / Nvtmrsswgd * model->BSIM4v6vtsswgd * T9;
          dT0_dVb = model->BSIM4v6vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9) ;
          DEXP(T0, T6, T10);
          dT6_dVb = T10 * dT0_dVb;
        }

          here->BSIM4v6gbs += here->BSIM4v6SjctTempRevSatCur * dT1_dVb
                                  + here->BSIM4v6SswTempRevSatCur * dT3_dVb
                                  + here->BSIM4v6SswgTempRevSatCur * dT5_dVb; 
          here->BSIM4v6cbs -= here->BSIM4v6SjctTempRevSatCur * (T1 - 1.0)
                                  + here->BSIM4v6SswTempRevSatCur * (T3 - 1.0)
                                  + here->BSIM4v6SswgTempRevSatCur * (T5 - 1.0); 
          here->BSIM4v6gbd += here->BSIM4v6DjctTempRevSatCur * dT2_dVb
                                  + here->BSIM4v6DswTempRevSatCur * dT4_dVb
                                  + here->BSIM4v6DswgTempRevSatCur * dT6_dVb; 
          here->BSIM4v6cbd -= here->BSIM4v6DjctTempRevSatCur * (T2 - 1.0) 
                                  + here->BSIM4v6DswTempRevSatCur * (T4 - 1.0)
                                  + here->BSIM4v6DswgTempRevSatCur * (T6 - 1.0); 

          /* End of diode DC model */

          if (vds >= 0.0)
          {   here->BSIM4v6mode = 1;
              Vds = vds;
              Vgs = vgs;
              Vbs = vbs;
              Vdb = vds - vbs;  /* WDLiu: for GIDL */
          }
          else
          {   here->BSIM4v6mode = -1;
              Vds = -vds;
              Vgs = vgd;
              Vbs = vbd;
              Vdb = -vbs;
          }


         /* dunga */
         if(model->BSIM4v6mtrlMod)
           {
             epsrox = 3.9;
             toxe = model->BSIM4v6eot;
             epssub = EPS0 * model->BSIM4v6epsrsub;
           }
         else
           {
             epsrox = model->BSIM4v6epsrox;
             toxe = model->BSIM4v6toxe;
             epssub = EPSSI;
           }


          T0 = Vbs - here->BSIM4v6vbsc - 0.001;
          T1 = sqrt(T0 * T0 - 0.004 * here->BSIM4v6vbsc);
          if (T0 >= 0.0)
          {   Vbseff = here->BSIM4v6vbsc + 0.5 * (T0 + T1);
              dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
          }
          else
          {   T2 = -0.002 / (T1 - T0);
              Vbseff = here->BSIM4v6vbsc * (1.0 + T2);
              dVbseff_dVb = T2 * here->BSIM4v6vbsc / T1;
          }

        /* JX: Correction to forward body bias  */
          T9 = 0.95 * pParam->BSIM4v6phi;
          T0 = T9 - Vbseff - 0.001;
          T1 = sqrt(T0 * T0 + 0.004 * T9);
          Vbseff = T9 - 0.5 * (T0 + T1);
          dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);
          Phis = pParam->BSIM4v6phi - Vbseff;
          dPhis_dVb = -1.0;
          sqrtPhis = sqrt(Phis);
          dsqrtPhis_dVb = -0.5 / sqrtPhis; 

          Xdep = pParam->BSIM4v6Xdep0 * sqrtPhis / pParam->BSIM4v6sqrtPhi;
          dXdep_dVb = (pParam->BSIM4v6Xdep0 / pParam->BSIM4v6sqrtPhi)
                    * dsqrtPhis_dVb;

          Leff = pParam->BSIM4v6leff;
          Vtm = model->BSIM4v6vtm;
          Vtm0 = model->BSIM4v6vtm0;

          /* Vth Calculation */
          T3 = sqrt(Xdep);
          V0 = pParam->BSIM4v6vbi - pParam->BSIM4v6phi;

          T0 = pParam->BSIM4v6dvt2 * Vbseff;
          if (T0 >= - 0.5)
          {   T1 = 1.0 + T0;
              T2 = pParam->BSIM4v6dvt2;
          }
          else
          {   T4 = 1.0 / (3.0 + 8.0 * T0);
              T1 = (1.0 + 3.0 * T0) * T4; 
              T2 = pParam->BSIM4v6dvt2 * T4 * T4;
          }
          lt1 = model->BSIM4v6factor1 * T3 * T1;
          dlt1_dVb = model->BSIM4v6factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

          T0 = pParam->BSIM4v6dvt2w * Vbseff;
          if (T0 >= - 0.5)
          {   T1 = 1.0 + T0;
              T2 = pParam->BSIM4v6dvt2w;
          }
          else
          {   T4 = 1.0 / (3.0 + 8.0 * T0);
              T1 = (1.0 + 3.0 * T0) * T4; 
              T2 = pParam->BSIM4v6dvt2w * T4 * T4;
          }
          ltw = model->BSIM4v6factor1 * T3 * T1;
          dltw_dVb = model->BSIM4v6factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

          T0 = pParam->BSIM4v6dvt1 * Leff / lt1;
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
          here->BSIM4v6thetavth = pParam->BSIM4v6dvt0 * Theta0;
          Delt_vth = here->BSIM4v6thetavth * V0;
          dDelt_vth_dVb = pParam->BSIM4v6dvt0 * dTheta0_dVb * V0;

          T0 = pParam->BSIM4v6dvt1w * pParam->BSIM4v6weff * Leff / ltw;
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
          T0 = pParam->BSIM4v6dvt0w * T5;
          T2 = T0 * V0;
          dT2_dVb = pParam->BSIM4v6dvt0w * dT5_dVb * V0;

          TempRatio =  ckt->CKTtemp / model->BSIM4v6tnom - 1.0;
          T0 = sqrt(1.0 + pParam->BSIM4v6lpe0 / Leff);
          T1 = pParam->BSIM4v6k1ox * (T0 - 1.0) * pParam->BSIM4v6sqrtPhi
             + (pParam->BSIM4v6kt1 + pParam->BSIM4v6kt1l / Leff
             + pParam->BSIM4v6kt2 * Vbseff) * TempRatio;
          Vth_NarrowW = toxe * pParam->BSIM4v6phi
                      / (pParam->BSIM4v6weff + pParam->BSIM4v6w0);

          T3 = here->BSIM4v6eta0 + pParam->BSIM4v6etab * Vbseff;
          if (T3 < 1.0e-4)
          {   T9 = 1.0 / (3.0 - 2.0e4 * T3);
              T3 = (2.0e-4 - T3) * T9;
              T4 = T9 * T9;
          }
          else
          {   T4 = 1.0;
          }
          dDIBL_Sft_dVd = T3 * pParam->BSIM4v6theta0vb0;
          DIBL_Sft = dDIBL_Sft_dVd * Vds;

           Lpe_Vb = sqrt(1.0 + pParam->BSIM4v6lpeb / Leff);

          Vth = model->BSIM4v6type * here->BSIM4v6vth0 + (pParam->BSIM4v6k1ox * sqrtPhis
              - pParam->BSIM4v6k1 * pParam->BSIM4v6sqrtPhi) * Lpe_Vb
              - here->BSIM4v6k2ox * Vbseff - Delt_vth - T2 + (pParam->BSIM4v6k3
              + pParam->BSIM4v6k3b * Vbseff) * Vth_NarrowW + T1 - DIBL_Sft;

          dVth_dVb = Lpe_Vb * pParam->BSIM4v6k1ox * dsqrtPhis_dVb - here->BSIM4v6k2ox
                   - dDelt_vth_dVb - dT2_dVb + pParam->BSIM4v6k3b * Vth_NarrowW
                   - pParam->BSIM4v6etab * Vds * pParam->BSIM4v6theta0vb0 * T4
                   + pParam->BSIM4v6kt2 * TempRatio;
          dVth_dVd = -dDIBL_Sft_dVd;


          /* Calculate n */
          tmp1 = epssub / Xdep;
          here->BSIM4v6nstar = model->BSIM4v6vtm / Charge_q * (model->BSIM4v6coxe
                           + tmp1 + pParam->BSIM4v6cit);  
          tmp2 = pParam->BSIM4v6nfactor * tmp1;
          tmp3 = pParam->BSIM4v6cdsc + pParam->BSIM4v6cdscb * Vbseff
               + pParam->BSIM4v6cdscd * Vds;
          tmp4 = (tmp2 + tmp3 * Theta0 + pParam->BSIM4v6cit) / model->BSIM4v6coxe;
          if (tmp4 >= -0.5)
          {   n = 1.0 + tmp4;
              dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                     + pParam->BSIM4v6cdscb * Theta0) / model->BSIM4v6coxe;
              dn_dVd = pParam->BSIM4v6cdscd * Theta0 / model->BSIM4v6coxe;
          }
          else
          {   T0 = 1.0 / (3.0 + 8.0 * tmp4);
              n = (1.0 + 3.0 * tmp4) * T0;
              T0 *= T0;
              dn_dVb = (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb
                     + pParam->BSIM4v6cdscb * Theta0) / model->BSIM4v6coxe * T0;
              dn_dVd = pParam->BSIM4v6cdscd * Theta0 / model->BSIM4v6coxe * T0;
          }


          /* Vth correction for Pocket implant */
           if (pParam->BSIM4v6dvtp0 > 0.0)
          {   T0 = -pParam->BSIM4v6dvtp1 * Vds;
              if (T0 < -EXP_THRESHOLD)
              {   T2 = MIN_EXP;
                  dT2_dVd = 0.0;
              }
              else
              {   T2 = exp(T0);
                  dT2_dVd = -pParam->BSIM4v6dvtp1 * T2;
              }

              T3 = Leff + pParam->BSIM4v6dvtp0 * (1.0 + T2);
              dT3_dVd = pParam->BSIM4v6dvtp0 * dT2_dVd;
              if (model->BSIM4v6tempMod < 2)
              {
                T4 = Vtm * log(Leff / T3);
                dT4_dVd = -Vtm * dT3_dVd / T3;
              }
              else
              {
                T4 = model->BSIM4v6vtm0 * log(Leff / T3);
                dT4_dVd = -model->BSIM4v6vtm0 * dT3_dVd / T3;
              }
              dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
              dDITS_Sft_dVb = T4 * dn_dVb;

              Vth -= n * T4;
              dVth_dVd -= dDITS_Sft_dVd;
              dVth_dVb -= dDITS_Sft_dVb;
          }
          here->BSIM4v6von = Vth;

          
          /* Poly Gate Si Depletion Effect */
          T0 = here->BSIM4v6vfb + pParam->BSIM4v6phi;
          if(model->BSIM4v6mtrlMod == 0)                 
            T1 = EPSSI;
          else
            T1 = model->BSIM4v6epsrgate * EPS0;


              BSIM4v6polyDepletion(T0, pParam->BSIM4v6ngate, T1, model->BSIM4v6coxe, vgs, &vgs_eff, &dvgs_eff_dvg);

              BSIM4v6polyDepletion(T0, pParam->BSIM4v6ngate, T1, model->BSIM4v6coxe, vgd, &vgd_eff, &dvgd_eff_dvg);
              
              if(here->BSIM4v6mode>0) {
                      Vgs_eff = vgs_eff;
                      dVgs_eff_dVg = dvgs_eff_dvg;
              } else {
                      Vgs_eff = vgd_eff;
                      dVgs_eff_dVg = dvgd_eff_dvg;
              }
              here->BSIM4v6vgs_eff = vgs_eff;
              here->BSIM4v6vgd_eff = vgd_eff;
              here->BSIM4v6dvgs_eff_dvg = dvgs_eff_dvg;
              here->BSIM4v6dvgd_eff_dvg = dvgd_eff_dvg;


          Vgst = Vgs_eff - Vth;

          /* Calculate Vgsteff */
          T0 = n * Vtm;
          T1 = pParam->BSIM4v6mstar * Vgst;
          T2 = T1 / T0;
          if (T2 > EXP_THRESHOLD)
          {   T10 = T1;
              dT10_dVg = pParam->BSIM4v6mstar * dVgs_eff_dVg;
              dT10_dVd = -dVth_dVd * pParam->BSIM4v6mstar;
              dT10_dVb = -dVth_dVb * pParam->BSIM4v6mstar;
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
              dT10_dVg = pParam->BSIM4v6mstar * ExpVgst / (1.0 + ExpVgst);
              dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
              dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
              dT10_dVg *= dVgs_eff_dVg;
          }

          T1 = pParam->BSIM4v6voffcbn - (1.0 - pParam->BSIM4v6mstar) * Vgst;
          T2 = T1 / T0;
          if (T2 < -EXP_THRESHOLD)
          {   T3 = model->BSIM4v6coxe * MIN_EXP / pParam->BSIM4v6cdep0;
              T9 = pParam->BSIM4v6mstar + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
          }
          else if (T2 > EXP_THRESHOLD)
          {   T3 = model->BSIM4v6coxe * MAX_EXP / pParam->BSIM4v6cdep0;
              T9 = pParam->BSIM4v6mstar + T3 * n;
              dT9_dVg = 0.0;
              dT9_dVd = dn_dVd * T3;
              dT9_dVb = dn_dVb * T3;
          }
          else
          {   ExpVgst = exp(T2);
              T3 = model->BSIM4v6coxe / pParam->BSIM4v6cdep0;
              T4 = T3 * ExpVgst;
              T5 = T1 * T4 / T0;
              T9 = pParam->BSIM4v6mstar + n * T4;
              dT9_dVg = T3 * (pParam->BSIM4v6mstar - 1.0) * ExpVgst / Vtm;
              dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
              dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
              dT9_dVg *= dVgs_eff_dVg;
          }

          here->BSIM4v6Vgsteff = Vgsteff = T10 / T9;
          T11 = T9 * T9;
          dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
          dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
          dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

          /* Calculate Effective Channel Geometry */
          T9 = sqrtPhis - pParam->BSIM4v6sqrtPhi;
          Weff = pParam->BSIM4v6weff - 2.0 * (pParam->BSIM4v6dwg * Vgsteff 
               + pParam->BSIM4v6dwb * T9); 
          dWeff_dVg = -2.0 * pParam->BSIM4v6dwg;
          dWeff_dVb = -2.0 * pParam->BSIM4v6dwb * dsqrtPhis_dVb;

          if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
          {   T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
              Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
              T0 *= T0 * 4.0e-16;
              dWeff_dVg *= T0;
              dWeff_dVb *= T0;
          }

          if (model->BSIM4v6rdsMod == 1)
              Rds = dRds_dVg = dRds_dVb = 0.0;
          else
          {   T0 = 1.0 + pParam->BSIM4v6prwg * Vgsteff;
              dT0_dVg = -pParam->BSIM4v6prwg / T0 / T0;
              T1 = pParam->BSIM4v6prwb * T9;
              dT1_dVb = pParam->BSIM4v6prwb * dsqrtPhis_dVb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01); /* 0.01 = 4.0 * 0.05 * 0.05 */
              dT3_dVg = 1.0 + T2 / (T3 - T2);
              dT3_dVb = dT3_dVg * dT1_dVb;
              dT3_dVg *= dT0_dVg;

              T4 = pParam->BSIM4v6rds0 * 0.5;
              Rds = pParam->BSIM4v6rdswmin + T3 * T4;
              dRds_dVg = T4 * dT3_dVg;
              dRds_dVb = T4 * dT3_dVb;

              if (Rds > 0.0)
                  here->BSIM4v6grdsw = 1.0 / Rds* here->BSIM4v6nf; /*4.6.2*/
              else
                  here->BSIM4v6grdsw = 0.0;
          }
          
          /* Calculate Abulk */
          T9 = 0.5 * pParam->BSIM4v6k1ox * Lpe_Vb / sqrtPhis;
          T1 = T9 + here->BSIM4v6k2ox - pParam->BSIM4v6k3b * Vth_NarrowW;
          dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

          T9 = sqrt(pParam->BSIM4v6xj * Xdep);
          tmp1 = Leff + 2.0 * T9;
          T5 = Leff / tmp1; 
          tmp2 = pParam->BSIM4v6a0 * T5;
          tmp3 = pParam->BSIM4v6weff + pParam->BSIM4v6b1; 
          tmp4 = pParam->BSIM4v6b0 / tmp3;
          T2 = tmp2 + tmp4;
          dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
          T6 = T5 * T5;
          T7 = T5 * T6;

          Abulk0 = 1.0 + T1 * T2; 
          dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

          T8 = pParam->BSIM4v6ags * pParam->BSIM4v6a0 * T7;
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
          here->BSIM4v6Abulk = Abulk;

          T2 = pParam->BSIM4v6keta * Vbseff;
          if (T2 >= -0.9)
          {   T0 = 1.0 / (1.0 + T2);
              dT0_dVb = -pParam->BSIM4v6keta * T0 * T0;
          }
          else
          {   T1 = 1.0 / (0.8 + T2);
              T0 = (17.0 + 20.0 * T2) * T1;
              dT0_dVb = -pParam->BSIM4v6keta * T1 * T1;
          }
          dAbulk_dVg *= T0;
          dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
          dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
          Abulk *= T0;
          Abulk0 *= T0;

          /* Mobility calculation */
          if (model->BSIM4v6mtrlMod)
            T14 = 2.0 * model->BSIM4v6type *(model->BSIM4v6phig - model->BSIM4v6easub - 0.5*model->BSIM4v6Eg0 + 0.45);   
          else
            T14 = 0.0;

          if (model->BSIM4v6mobMod == 0)
            { T0 = Vgsteff + Vth + Vth - T14;
              T2 = pParam->BSIM4v6ua + pParam->BSIM4v6uc * Vbseff;
              T3 = T0 / toxe;
              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v6ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T3 * (T2 + pParam->BSIM4v6ub * T3) + T6;
              T7 = - 2.0 * T6 * T9;
               T11 = T7 * Vth/T12;
              dDenomi_dVg = (T2 + 2.0 * pParam->BSIM4v6ub * T3) / toxe;
              T13 = 2.0 * (dDenomi_dVg + T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + pParam->BSIM4v6uc * T3;
              dDenomi_dVg+= T7;
          }
          else if (model->BSIM4v6mobMod == 1)
          {   T0 = Vgsteff + Vth + Vth - T14;
              T2 = 1.0 + pParam->BSIM4v6uc * Vbseff;
              T3 = T0 / toxe;
              T4 = T3 * (pParam->BSIM4v6ua + pParam->BSIM4v6ub * T3);
              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v6ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T4 * T2 + T6;
              T7 = - 2.0 * T6 * T9;
               T11 = T7 * Vth/T12;
              dDenomi_dVg = (pParam->BSIM4v6ua + 2.0 * pParam->BSIM4v6ub * T3) * T2 / toxe;
              T13 = 2.0 * (dDenomi_dVg + T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + pParam->BSIM4v6uc * T4;
              dDenomi_dVg+= T7;
          }
          else if (model->BSIM4v6mobMod == 2)
          {   T0 = (Vgsteff + here->BSIM4v6vtfbphi1) / toxe;
              T1 = exp(pParam->BSIM4v6eu * log(T0));
              dT1_dVg = T1 * pParam->BSIM4v6eu / T0 / toxe;
              T2 = pParam->BSIM4v6ua + pParam->BSIM4v6uc * Vbseff;
              T3 = T0 / toxe; /*Do we need it?*/
                      
              T12 = sqrt(Vth * Vth + 0.0001);
              T9 = 1.0/(Vgsteff + 2*T12);
              T10 = T9*toxe;
              T8 = pParam->BSIM4v6ud * T10 * T10 * Vth;
              T6 = T8 * Vth;
              T5 = T1 * T2 + T6;
              T7 = - 2.0 * T6 * T9;
               T11 = T7 * Vth/T12;
               dDenomi_dVg = T2 * dT1_dVg + T7;
              T13 = 2.0 * (T11 + T8);
              dDenomi_dVd = T13 * dVth_dVd;
              dDenomi_dVb = T13 * dVth_dVb + T1 * pParam->BSIM4v6uc;
          }
          /*high K mobility*/
         else 
      {                                 
                
                             
                   /*univsersal mobility*/
                 T0 = (Vgsteff + here->BSIM4v6vtfbphi1)* 1.0e-8 / toxe/6.0;
             T1 = exp(pParam->BSIM4v6eu * log(T0)); 
                 dT1_dVg = T1 * pParam->BSIM4v6eu * 1.0e-8/ T0 / toxe/6.0;
             T2 = pParam->BSIM4v6ua + pParam->BSIM4v6uc * Vbseff;
                 
                  /*Coulombic*/
                 VgsteffVth = pParam->BSIM4v6VgsteffVth;
                 
                 T10 = exp(pParam->BSIM4v6ucs * log(0.5 + 0.5 * Vgsteff/VgsteffVth));
                          T11 =  pParam->BSIM4v6ud/T10;
                 dT11_dVg = - 0.5 * pParam->BSIM4v6ucs * T11 /(0.5 + 0.5*Vgsteff/VgsteffVth)/VgsteffVth;
                 
                dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
                dDenomi_dVd = 0.0;
                dDenomi_dVb = T1 * pParam->BSIM4v6uc;
                
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
          
        
          here->BSIM4v6ueff = ueff = here->BSIM4v6u0temp / Denomi;
          T9 = -ueff / Denomi;
          dueff_dVg = T9 * dDenomi_dVg;
          dueff_dVd = T9 * dDenomi_dVd;
          dueff_dVb = T9 * dDenomi_dVb;

          /* Saturation Drain Voltage  Vdsat */
          WVCox = Weff * here->BSIM4v6vsattemp * model->BSIM4v6coxe;
          WVCoxRds = WVCox * Rds; 

          Esat = 2.0 * here->BSIM4v6vsattemp / ueff;
          here->BSIM4v6EsatL = EsatL = Esat * Leff;
          T0 = -EsatL /ueff;
          dEsatL_dVg = T0 * dueff_dVg;
          dEsatL_dVd = T0 * dueff_dVd;
          dEsatL_dVb = T0 * dueff_dVb;
  
          /* Sqrt() */
          a1 = pParam->BSIM4v6a1;
          if (a1 == 0.0)
          {   Lambda = pParam->BSIM4v6a2;
              dLambda_dVg = 0.0;
          }
          else if (a1 > 0.0)
          {   T0 = 1.0 - pParam->BSIM4v6a2;
              T1 = T0 - pParam->BSIM4v6a1 * Vgsteff - 0.0001;
              T2 = sqrt(T1 * T1 + 0.0004 * T0);
              Lambda = pParam->BSIM4v6a2 + T0 - 0.5 * (T1 + T2);
              dLambda_dVg = 0.5 * pParam->BSIM4v6a1 * (1.0 + T1 / T2);
          }
          else
          {   T1 = pParam->BSIM4v6a2 + pParam->BSIM4v6a1 * Vgsteff - 0.0001;
              T2 = sqrt(T1 * T1 + 0.0004 * pParam->BSIM4v6a2);
              Lambda = 0.5 * (T1 + T2);
              dLambda_dVg = 0.5 * pParam->BSIM4v6a1 * (1.0 + T1 / T2);
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
          here->BSIM4v6vdsat = Vdsat;

          /* Calculate Vdseff */
          T1 = Vdsat - Vds - pParam->BSIM4v6delta;
          dT1_dVg = dVdsat_dVg;
          dT1_dVd = dVdsat_dVd - 1.0;
          dT1_dVb = dVdsat_dVb;

          T2 = sqrt(T1 * T1 + 4.0 * pParam->BSIM4v6delta * Vdsat);
          T0 = T1 / T2;
             T9 = 2.0 * pParam->BSIM4v6delta;
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
          here->BSIM4v6Vdseff = Vdseff;
          
          /* Velocity Overshoot */
        if((model->BSIM4v6lambdaGiven) && (model->BSIM4v6lambda > 0.0) )
        {  
          T1 =  Leff * ueff;
          T2 = pParam->BSIM4v6lambda / T1;
          T3 = -T2 / T1 * Leff;
          dT2_dVd = T3 * dueff_dVd;
          dT2_dVg = T3 * dueff_dVg;
          dT2_dVb = T3 * dueff_dVb;
          T5 = 1.0 / (Esat * pParam->BSIM4v6litl);
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
          here->BSIM4v6EsatL = EsatL;
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
              
          tmp1 = here->BSIM4v6vtfbphi2;
          tmp2 = 2.0e8 * model->BSIM4v6toxp;
          dT0_dVg = 1.0 / tmp2;
          T0 = (Vgsteff + tmp1) * dT0_dVg;

          tmp3 = exp(model->BSIM4v6bdos * 0.7 * log(T0));
          T1 = 1.0 + tmp3;
          T2 = model->BSIM4v6bdos * 0.7 * tmp3 / T0;
          Tcen = model->BSIM4v6ados * 1.9e-9 / T1;
          dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

          Coxeff = epssub * model->BSIM4v6coxp
                 / (epssub + model->BSIM4v6coxp * Tcen);
          dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

          CoxeffWovL = Coxeff * Weff / Leff;
          beta = ueff * CoxeffWovL;
          T3 = ueff / Leff;
          dbeta_dVg = CoxeffWovL * dueff_dVg + T3
                    * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
          dbeta_dVd = CoxeffWovL * dueff_dVd;
          dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

          here->BSIM4v6AbovVgst2Vtm = Abulk / Vgst2Vtm;
          T0 = 1.0 - 0.5 * Vdseff * here->BSIM4v6AbovVgst2Vtm;
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

          if (pParam->BSIM4v6fprout <= 0.0)
          {   FP = 1.0;
              dFP_dVg = 0.0;
          }
          else
          {   T9 = pParam->BSIM4v6fprout * sqrt(Leff) / Vgst2Vtm;
              FP = 1.0 / (1.0 + T9);
              dFP_dVg = FP * FP * T9 / Vgst2Vtm;
          }

          /* Calculate VACLM */
          T8 = pParam->BSIM4v6pvag / EsatL;
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

          if ((pParam->BSIM4v6pclm > MIN_EXP) && (diffVds > 1.0e-10))
          {   T0 = 1.0 + Rds * Idl;
              dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
              dT0_dVd = Rds * dIdl_dVd;
              dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

              T2 = Vdsat / Esat;
              T1 = Leff + T2;
              dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
              dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
              dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

              Cclm = FP * PvagTerm * T0 * T1 / (pParam->BSIM4v6pclm * pParam->BSIM4v6litl);
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
          if (pParam->BSIM4v6thetaRout > MIN_EXP)
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
              T2 = pParam->BSIM4v6thetaRout;
              VADIBL = (Vgst2Vtm - T0 / T1) / T2;
              dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
              dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
              dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

              T7 = pParam->BSIM4v6pdiblb * Vbseff;
              if (T7 >= -0.9)
              {   T3 = 1.0 / (1.0 + T7);
                  VADIBL *= T3;
                  dVADIBL_dVg *= T3;
                  dVADIBL_dVb = (dVADIBL_dVb - VADIBL * pParam->BSIM4v6pdiblb)
                              * T3;
                  dVADIBL_dVd *= T3;
              }
              else
              {   T4 = 1.0 / (0.8 + T7);
                  T3 = (17.0 + 20.0 * T7) * T4;
                  dVADIBL_dVg *= T3;
                  dVADIBL_dVb = dVADIBL_dVb * T3
                              - VADIBL * pParam->BSIM4v6pdiblb * T4 * T4;
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
          T0 = pParam->BSIM4v6pditsd * Vds;
          if (T0 > EXP_THRESHOLD)
          {   T1 = MAX_EXP;
              dT1_dVd = 0;
          }
          else
          {   T1 = exp(T0);
              dT1_dVd = T1 * pParam->BSIM4v6pditsd;
          }

          if (pParam->BSIM4v6pdits > MIN_EXP)
          {   T2 = 1.0 + model->BSIM4v6pditsl * Leff;
              VADITS = (1.0 + T2 * T1) / pParam->BSIM4v6pdits;
              dVADITS_dVg = VADITS * dFP_dVg;
              dVADITS_dVd = FP * T2 * dT1_dVd / pParam->BSIM4v6pdits;
              VADITS *= FP;
          }
          else
          {   VADITS = MAX_EXP;
              dVADITS_dVg = dVADITS_dVd = 0;
          }

          /* Calculate VASCBE */
          if ((pParam->BSIM4v6pscbe2 > 0.0)&&(pParam->BSIM4v6pscbe1>=0.0)) /*4.6.2*/
          {   if (diffVds > pParam->BSIM4v6pscbe1 * pParam->BSIM4v6litl
                  / EXP_THRESHOLD)
              {   T0 =  pParam->BSIM4v6pscbe1 * pParam->BSIM4v6litl / diffVds;
                  VASCBE = Leff * exp(T0) / pParam->BSIM4v6pscbe2;
                  T1 = T0 * VASCBE / diffVds;
                  dVASCBE_dVg = T1 * dVdseff_dVg;
                  dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                  dVASCBE_dVb = T1 * dVdseff_dVb;
              }
              else
              {   VASCBE = MAX_EXP * Leff/pParam->BSIM4v6pscbe2;
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
          tmp = pParam->BSIM4v6alpha0 + pParam->BSIM4v6alpha1 * Leff;
          if ((tmp <= 0.0) || (pParam->BSIM4v6beta0 <= 0.0))
          {   Isub = Gbd = Gbb = Gbg = 0.0;
          }
          else
          {   T2 = tmp / Leff;
              if (diffVds > pParam->BSIM4v6beta0 / EXP_THRESHOLD)
              {   T0 = -pParam->BSIM4v6beta0 / diffVds;
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
          here->BSIM4v6csub = Isub;
          here->BSIM4v6gbbs = Gbb;
          here->BSIM4v6gbgs = Gbg;
          here->BSIM4v6gbds = Gbd;

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
        if((model->BSIM4v6vtlGiven) && (model->BSIM4v6vtl > 0.0) ) {
          T12 = 1.0 / Leff / CoxeffWovL;
          T11 = T12 / Vgsteff;
          T10 = -T11 / Vgsteff;
          vs = cdrain * T11; /* vs */
          dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
          dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
          dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
          T0 = 2 * MM;
          T1 = vs / (pParam->BSIM4v6vtl * pParam->BSIM4v6tfactor);
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

          here->BSIM4v6gds = Gds;
          here->BSIM4v6gm = Gm;
          here->BSIM4v6gmbs = Gmb;
          here->BSIM4v6IdovVds = Ids;
          if( here->BSIM4v6IdovVds <= 1.0e-9) here->BSIM4v6IdovVds = 1.0e-9;

          /* Calculate Rg */
          if ((here->BSIM4v6rgateMod > 1) ||
              (here->BSIM4v6trnqsMod != 0) || (here->BSIM4v6acnqsMod != 0))
          {   T9 = pParam->BSIM4v6xrcrg2 * model->BSIM4v6vtm;
              T0 = T9 * beta;
              dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
              dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
              dT0_dVg = dbeta_dVg * T9;

              here->BSIM4v6gcrg = pParam->BSIM4v6xrcrg1 * ( T0 + Ids);
              here->BSIM4v6gcrgd = pParam->BSIM4v6xrcrg1 * (dT0_dVd + tmp1);
              here->BSIM4v6gcrgb = pParam->BSIM4v6xrcrg1 * (dT0_dVb + tmp2)
                                * dVbseff_dVb;        
              here->BSIM4v6gcrgg = pParam->BSIM4v6xrcrg1 * (dT0_dVg + tmp3)
                               * dVgsteff_dVg;

              if (here->BSIM4v6nf != 1.0)
              {   here->BSIM4v6gcrg *= here->BSIM4v6nf; 
                  here->BSIM4v6gcrgg *= here->BSIM4v6nf;
                  here->BSIM4v6gcrgd *= here->BSIM4v6nf;
                  here->BSIM4v6gcrgb *= here->BSIM4v6nf;
              }

              if (here->BSIM4v6rgateMod == 2)
              {   T10 = here->BSIM4v6grgeltd * here->BSIM4v6grgeltd;
                  T11 = here->BSIM4v6grgeltd + here->BSIM4v6gcrg;
                  here->BSIM4v6gcrg = here->BSIM4v6grgeltd * here->BSIM4v6gcrg / T11;
                  T12 = T10 / T11 / T11;
                  here->BSIM4v6gcrgg *= T12;
                  here->BSIM4v6gcrgd *= T12;
                  here->BSIM4v6gcrgb *= T12;
              }
              here->BSIM4v6gcrgs = -(here->BSIM4v6gcrgg + here->BSIM4v6gcrgd
                               + here->BSIM4v6gcrgb);
          }


          /* Calculate bias-dependent external S/D resistance */
          if (model->BSIM4v6rdsMod)
          {   /* Rs(V) */
              T0 = vgs - pParam->BSIM4v6vfbsd;
              T1 = sqrt(T0 * T0 + 1.0e-4);
              vgs_eff = 0.5 * (T0 + T1);
              dvgs_eff_dvg = vgs_eff / T1;

              T0 = 1.0 + pParam->BSIM4v6prwg * vgs_eff;
              dT0_dvg = -pParam->BSIM4v6prwg / T0 / T0 * dvgs_eff_dvg;
              T1 = -pParam->BSIM4v6prwb * vbs;
              dT1_dvb = -pParam->BSIM4v6prwb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01);
              dT3_dvg = T3 / (T3 - T2);
              dT3_dvb = dT3_dvg * dT1_dvb;
              dT3_dvg *= dT0_dvg;

              T4 = pParam->BSIM4v6rs0 * 0.5;
              Rs = pParam->BSIM4v6rswmin + T3 * T4;
              dRs_dvg = T4 * dT3_dvg;
              dRs_dvb = T4 * dT3_dvb;

              T0 = 1.0 + here->BSIM4v6sourceConductance * Rs;
              here->BSIM4v6gstot = here->BSIM4v6sourceConductance / T0;
              T0 = -here->BSIM4v6gstot * here->BSIM4v6gstot;
              dgstot_dvd = 0.0; /* place holder */
              dgstot_dvg = T0 * dRs_dvg;
              dgstot_dvb = T0 * dRs_dvb;
              dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

              /* Rd(V) */
              T0 = vgd - pParam->BSIM4v6vfbsd;
              T1 = sqrt(T0 * T0 + 1.0e-4);
              vgd_eff = 0.5 * (T0 + T1);
              dvgd_eff_dvg = vgd_eff / T1;

              T0 = 1.0 + pParam->BSIM4v6prwg * vgd_eff;
              dT0_dvg = -pParam->BSIM4v6prwg / T0 / T0 * dvgd_eff_dvg;
              T1 = -pParam->BSIM4v6prwb * vbd;
              dT1_dvb = -pParam->BSIM4v6prwb;

              T2 = 1.0 / T0 + T1;
              T3 = T2 + sqrt(T2 * T2 + 0.01);
              dT3_dvg = T3 / (T3 - T2);
              dT3_dvb = dT3_dvg * dT1_dvb;
              dT3_dvg *= dT0_dvg;

              T4 = pParam->BSIM4v6rd0 * 0.5;
              Rd = pParam->BSIM4v6rdwmin + T3 * T4;
              dRd_dvg = T4 * dT3_dvg;
              dRd_dvb = T4 * dT3_dvb;

              T0 = 1.0 + here->BSIM4v6drainConductance * Rd;
              here->BSIM4v6gdtot = here->BSIM4v6drainConductance / T0;
              T0 = -here->BSIM4v6gdtot * here->BSIM4v6gdtot;
              dgdtot_dvs = 0.0;
              dgdtot_dvg = T0 * dRd_dvg;
              dgdtot_dvb = T0 * dRd_dvb;
              dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

              here->BSIM4v6gstotd = vses * dgstot_dvd;
              here->BSIM4v6gstotg = vses * dgstot_dvg;
              here->BSIM4v6gstots = vses * dgstot_dvs;
              here->BSIM4v6gstotb = vses * dgstot_dvb;

              T2 = vdes - vds;
              here->BSIM4v6gdtotd = T2 * dgdtot_dvd;
              here->BSIM4v6gdtotg = T2 * dgdtot_dvg;
              here->BSIM4v6gdtots = T2 * dgdtot_dvs;
              here->BSIM4v6gdtotb = T2 * dgdtot_dvb;
          }
          else /* WDLiu: for bypass */
          {   here->BSIM4v6gstot = here->BSIM4v6gstotd = here->BSIM4v6gstotg = 0.0;
              here->BSIM4v6gstots = here->BSIM4v6gstotb = 0.0;
              here->BSIM4v6gdtot = here->BSIM4v6gdtotd = here->BSIM4v6gdtotg = 0.0;
              here->BSIM4v6gdtots = here->BSIM4v6gdtotb = 0.0;
          }

          /* Calculate GIDL current */
          vgs_eff = here->BSIM4v6vgs_eff;
          dvgs_eff_dvg = here->BSIM4v6dvgs_eff_dvg;
          if(model->BSIM4v6mtrlMod == 0)
            T0 = 3.0 * toxe;
          else
            T0 = model->BSIM4v6epsrsub * toxe / epsrox;

          if(model->BSIM4v6mtrlMod ==0)
            T1 = (vds - vgs_eff - pParam->BSIM4v6egidl ) / T0;
          else
            T1 = (vds - vgs_eff - pParam->BSIM4v6egidl + pParam->BSIM4v6vfbsd) / T0;
          if ((pParam->BSIM4v6agidl <= 0.0) || (pParam->BSIM4v6bgidl <= 0.0)
              || (T1 <= 0.0) || (pParam->BSIM4v6cgidl <= 0.0) || (vbd > 0.0))
              Igidl = Ggidld = Ggidlg = Ggidlb = 0.0;
          else {
              dT1_dVd = 1.0 / T0;
              dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
              T2 = pParam->BSIM4v6bgidl / T1;
              if (T2 < 100.0)
              {   Igidl = pParam->BSIM4v6agidl * pParam->BSIM4v6weffCJ * T1 * exp(-T2);
                  T3 = Igidl * (1.0 + T2) / T1;
                  Ggidld = T3 * dT1_dVd;
                  Ggidlg = T3 * dT1_dVg;
              }
              else
              {   Igidl = pParam->BSIM4v6agidl * pParam->BSIM4v6weffCJ * 3.720075976e-44;
                  Ggidld = Igidl * dT1_dVd;
                  Ggidlg = Igidl * dT1_dVg;
                  Igidl *= T1;
              }
        
              T4 = vbd * vbd;
              T5 = -vbd * T4;
              T6 = pParam->BSIM4v6cgidl + T5;
              T7 = T5 / T6;
              T8 = 3.0 * pParam->BSIM4v6cgidl * T4 / T6 / T6;
              Ggidld = Ggidld * T7 + Igidl * T8;
              Ggidlg = Ggidlg * T7;
              Ggidlb = -Igidl * T8;
              Igidl *= T7;
          }
          here->BSIM4v6Igidl = Igidl;
          here->BSIM4v6ggidld = Ggidld;
          here->BSIM4v6ggidlg = Ggidlg;
          here->BSIM4v6ggidlb = Ggidlb;

        /* Calculate GISL current  */
          vgd_eff = here->BSIM4v6vgd_eff;
          dvgd_eff_dvg = here->BSIM4v6dvgd_eff_dvg;

          if(model->BSIM4v6mtrlMod ==0)
          T1 = (-vds - vgd_eff - pParam->BSIM4v6egisl ) / T0;
          else
          T1 = (-vds - vgd_eff - pParam->BSIM4v6egisl + pParam->BSIM4v6vfbsd ) / T0;

          if ((pParam->BSIM4v6agisl <= 0.0) || (pParam->BSIM4v6bgisl <= 0.0)
              || (T1 <= 0.0) || (pParam->BSIM4v6cgisl <= 0.0) || (vbs > 0.0))
              Igisl = Ggisls = Ggislg = Ggislb = 0.0;
          else {
              dT1_dVd = 1.0 / T0;
              dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
              T2 = pParam->BSIM4v6bgisl / T1;
              if (T2 < 100.0) 
              {   Igisl = pParam->BSIM4v6agisl * pParam->BSIM4v6weffCJ * T1 * exp(-T2);
                  T3 = Igisl * (1.0 + T2) / T1;
                  Ggisls = T3 * dT1_dVd;
                  Ggislg = T3 * dT1_dVg;
              }
              else 
              {   Igisl = pParam->BSIM4v6agisl * pParam->BSIM4v6weffCJ * 3.720075976e-44;
                  Ggisls = Igisl * dT1_dVd;
                  Ggislg = Igisl * dT1_dVg;
                  Igisl *= T1;
              }
        
              T4 = vbs * vbs;
              T5 = -vbs * T4;
              T6 = pParam->BSIM4v6cgisl + T5;
              T7 = T5 / T6;
              T8 = 3.0 * pParam->BSIM4v6cgisl * T4 / T6 / T6;
              Ggisls = Ggisls * T7 + Igisl * T8;
              Ggislg = Ggislg * T7;
              Ggislb = -Igisl * T8;
              Igisl *= T7;
          }
          here->BSIM4v6Igisl = Igisl;
          here->BSIM4v6ggisls = Ggisls;
          here->BSIM4v6ggislg = Ggislg;
          here->BSIM4v6ggislb = Ggislb;


          /* Calculate gate tunneling current */
          if ((model->BSIM4v6igcMod != 0) || (model->BSIM4v6igbMod != 0))
          {   Vfb = here->BSIM4v6vfbzb;
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

              T0 = 0.5 * pParam->BSIM4v6k1ox;
              T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
              if (pParam->BSIM4v6k1ox == 0.0)
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
                  Voxdepinv = pParam->BSIM4v6k1ox * (T1 - T0);
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

          if(model->BSIM4v6tempMod < 2)
                  tmp = Vtm;
          else /* model->BSIM4v6tempMod = 2 , 3*/
                  tmp = Vtm0;
          if (model->BSIM4v6igcMod)
          {   T0 = tmp * pParam->BSIM4v6nigc;
              if(model->BSIM4v6igcMod == 1) {
                      VxNVt = (Vgs_eff - model->BSIM4v6type * here->BSIM4v6vth0) / T0;
                      if (VxNVt > EXP_THRESHOLD)
                      {   Vaux = Vgs_eff - model->BSIM4v6type * here->BSIM4v6vth0;
                          dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                      }
              } else if (model->BSIM4v6igcMod == 2) {
                VxNVt = (Vgs_eff - here->BSIM4v6von) / T0;
                if (VxNVt > EXP_THRESHOLD)
                {   Vaux = Vgs_eff - here->BSIM4v6von;
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
                  if(model->BSIM4v6igcMod == 1) {
                        dVaux_dVd = 0.0;
                          dVaux_dVb = 0.0;
                  } else if (model->BSIM4v6igcMod == 2) {
                        dVaux_dVd = -dVgs_eff_dVg * dVth_dVd;
                        dVaux_dVb = -dVgs_eff_dVg * dVth_dVb;
                  }
                  dVaux_dVg *= dVgs_eff_dVg;
              }

              T2 = Vgs_eff * Vaux;
              dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
              dT2_dVd = Vgs_eff * dVaux_dVd;
              dT2_dVb = Vgs_eff * dVaux_dVb;

              T11 = pParam->BSIM4v6Aechvb;
              T12 = pParam->BSIM4v6Bechvb;
              T3 = pParam->BSIM4v6aigc * pParam->BSIM4v6cigc
                 - pParam->BSIM4v6bigc;
              T4 = pParam->BSIM4v6bigc * pParam->BSIM4v6cigc;
              T5 = T12 * (pParam->BSIM4v6aigc + T3 * Voxdepinv
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

              if (model->BSIM4v6pigcdGiven)
              {   Pigcd = pParam->BSIM4v6pigcd;
                  dPigcd_dVg = dPigcd_dVd = dPigcd_dVb = 0.0;
              }
              else
              {   T11 = pParam->BSIM4v6Bechvb * toxe;
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

              here->BSIM4v6Igcs = Igcs;
              here->BSIM4v6gIgcsg = dIgcs_dVg;
              here->BSIM4v6gIgcsd = dIgcs_dVd;
              here->BSIM4v6gIgcsb =  dIgcs_dVb * dVbseff_dVb;
              here->BSIM4v6Igcd = Igcd;
              here->BSIM4v6gIgcdg = dIgcd_dVg;
              here->BSIM4v6gIgcdd = dIgcd_dVd;
              here->BSIM4v6gIgcdb = dIgcd_dVb * dVbseff_dVb;

              T0 = vgs - (pParam->BSIM4v6vfbsd + pParam->BSIM4v6vfbsdoff);
              vgs_eff = sqrt(T0 * T0 + 1.0e-4);
              dvgs_eff_dvg = T0 / vgs_eff;

              T2 = vgs * vgs_eff;
              dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
              T11 = pParam->BSIM4v6AechvbEdgeS;
              T12 = pParam->BSIM4v6BechvbEdge;
              T3 = pParam->BSIM4v6aigs * pParam->BSIM4v6cigs
                 - pParam->BSIM4v6bigs;
              T4 = pParam->BSIM4v6bigs * pParam->BSIM4v6cigs;
              T5 = T12 * (pParam->BSIM4v6aigs + T3 * vgs_eff
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


              T0 = vgd - (pParam->BSIM4v6vfbsd + pParam->BSIM4v6vfbsdoff);
              vgd_eff = sqrt(T0 * T0 + 1.0e-4);
              dvgd_eff_dvg = T0 / vgd_eff;

              T2 = vgd * vgd_eff;
              dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
              T11 = pParam->BSIM4v6AechvbEdgeD;
              T3 = pParam->BSIM4v6aigd * pParam->BSIM4v6cigd
                 - pParam->BSIM4v6bigd;
              T4 = pParam->BSIM4v6bigd * pParam->BSIM4v6cigd;
              T5 = T12 * (pParam->BSIM4v6aigd + T3 * vgd_eff
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

              here->BSIM4v6Igs = Igs;
              here->BSIM4v6gIgsg = dIgs_dVg;
              here->BSIM4v6gIgss = dIgs_dVs;
              here->BSIM4v6Igd = Igd;
              here->BSIM4v6gIgdg = dIgd_dVg;
              here->BSIM4v6gIgdd = dIgd_dVd;
          }
          else
          {   here->BSIM4v6Igcs = here->BSIM4v6gIgcsg = here->BSIM4v6gIgcsd
                               = here->BSIM4v6gIgcsb = 0.0;
              here->BSIM4v6Igcd = here->BSIM4v6gIgcdg = here->BSIM4v6gIgcdd
                                    = here->BSIM4v6gIgcdb = 0.0;
              here->BSIM4v6Igs = here->BSIM4v6gIgsg = here->BSIM4v6gIgss = 0.0;
              here->BSIM4v6Igd = here->BSIM4v6gIgdg = here->BSIM4v6gIgdd = 0.0;
          }

          if (model->BSIM4v6igbMod)
          {   T0 = tmp * pParam->BSIM4v6nigbacc;
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

              T11 = 4.97232e-7 * pParam->BSIM4v6weff
                  * pParam->BSIM4v6leff * pParam->BSIM4v6ToxRatio;
              T12 = -7.45669e11 * toxe;
              T3 = pParam->BSIM4v6aigbacc * pParam->BSIM4v6cigbacc
                 - pParam->BSIM4v6bigbacc;
              T4 = pParam->BSIM4v6bigbacc * pParam->BSIM4v6cigbacc;
              T5 = T12 * (pParam->BSIM4v6aigbacc + T3 * Voxacc
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


              T0 = tmp * pParam->BSIM4v6nigbinv;
              T1 = Voxdepinv - pParam->BSIM4v6eigbinv;
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
              T3 = pParam->BSIM4v6aigbinv * pParam->BSIM4v6cigbinv
                 - pParam->BSIM4v6bigbinv;
              T4 = pParam->BSIM4v6bigbinv * pParam->BSIM4v6cigbinv;
              T5 = T12 * (pParam->BSIM4v6aigbinv + T3 * Voxdepinv
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

              here->BSIM4v6Igb = Igbinv + Igbacc;
              here->BSIM4v6gIgbg = dIgbinv_dVg + dIgbacc_dVg;
              here->BSIM4v6gIgbd = dIgbinv_dVd;
              here->BSIM4v6gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
          }
          else
          {  here->BSIM4v6Igb = here->BSIM4v6gIgbg = here->BSIM4v6gIgbd
                            = here->BSIM4v6gIgbs = here->BSIM4v6gIgbb = 0.0;   
          } /* End of Gate current */

          if (here->BSIM4v6nf != 1.0)
          {   cdrain *= here->BSIM4v6nf;
              here->BSIM4v6gds *= here->BSIM4v6nf;
              here->BSIM4v6gm *= here->BSIM4v6nf;
              here->BSIM4v6gmbs *= here->BSIM4v6nf;
              here->BSIM4v6IdovVds *= here->BSIM4v6nf;
                   
              here->BSIM4v6gbbs *= here->BSIM4v6nf;
              here->BSIM4v6gbgs *= here->BSIM4v6nf;
              here->BSIM4v6gbds *= here->BSIM4v6nf;
              here->BSIM4v6csub *= here->BSIM4v6nf;

              here->BSIM4v6Igidl *= here->BSIM4v6nf;
              here->BSIM4v6ggidld *= here->BSIM4v6nf;
              here->BSIM4v6ggidlg *= here->BSIM4v6nf;
              here->BSIM4v6ggidlb *= here->BSIM4v6nf;
        
              here->BSIM4v6Igisl *= here->BSIM4v6nf;
              here->BSIM4v6ggisls *= here->BSIM4v6nf;
              here->BSIM4v6ggislg *= here->BSIM4v6nf;
              here->BSIM4v6ggislb *= here->BSIM4v6nf;

              here->BSIM4v6Igcs *= here->BSIM4v6nf;
              here->BSIM4v6gIgcsg *= here->BSIM4v6nf;
              here->BSIM4v6gIgcsd *= here->BSIM4v6nf;
              here->BSIM4v6gIgcsb *= here->BSIM4v6nf;
              here->BSIM4v6Igcd *= here->BSIM4v6nf;
              here->BSIM4v6gIgcdg *= here->BSIM4v6nf;
              here->BSIM4v6gIgcdd *= here->BSIM4v6nf;
              here->BSIM4v6gIgcdb *= here->BSIM4v6nf;

              here->BSIM4v6Igs *= here->BSIM4v6nf;
              here->BSIM4v6gIgsg *= here->BSIM4v6nf;
              here->BSIM4v6gIgss *= here->BSIM4v6nf;
              here->BSIM4v6Igd *= here->BSIM4v6nf;
              here->BSIM4v6gIgdg *= here->BSIM4v6nf;
              here->BSIM4v6gIgdd *= here->BSIM4v6nf;

              here->BSIM4v6Igb *= here->BSIM4v6nf;
              here->BSIM4v6gIgbg *= here->BSIM4v6nf;
              here->BSIM4v6gIgbd *= here->BSIM4v6nf;
              here->BSIM4v6gIgbb *= here->BSIM4v6nf;
          }

          here->BSIM4v6ggidls = -(here->BSIM4v6ggidld + here->BSIM4v6ggidlg
                            + here->BSIM4v6ggidlb);
          here->BSIM4v6ggisld = -(here->BSIM4v6ggisls + here->BSIM4v6ggislg
                            + here->BSIM4v6ggislb);
          here->BSIM4v6gIgbs = -(here->BSIM4v6gIgbg + here->BSIM4v6gIgbd
                           + here->BSIM4v6gIgbb);
          here->BSIM4v6gIgcss = -(here->BSIM4v6gIgcsg + here->BSIM4v6gIgcsd
                            + here->BSIM4v6gIgcsb);
          here->BSIM4v6gIgcds = -(here->BSIM4v6gIgcdg + here->BSIM4v6gIgcdd
                            + here->BSIM4v6gIgcdb);
          here->BSIM4v6cd = cdrain;


          if (model->BSIM4v6tnoiMod == 0)
          {   Abulk = Abulk0 * pParam->BSIM4v6abulkCVfactor;
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
              here->BSIM4v6qinv = Coxeff * pParam->BSIM4v6weffCV * here->BSIM4v6nf
                              * pParam->BSIM4v6leffCV
                              * (Vgsteff - 0.5 * T0 + Abulk * T3);
          }

          /*
           *  BSIM4v6 C-V begins
           */

          if ((model->BSIM4v6xpart < 0) || (!ChargeComputationNeeded))
          {   qgate  = qdrn = qsrc = qbulk = 0.0;
              here->BSIM4v6cggb = here->BSIM4v6cgsb = here->BSIM4v6cgdb = 0.0;
              here->BSIM4v6cdgb = here->BSIM4v6cdsb = here->BSIM4v6cddb = 0.0;
              here->BSIM4v6cbgb = here->BSIM4v6cbsb = here->BSIM4v6cbdb = 0.0;
              here->BSIM4v6csgb = here->BSIM4v6cssb = here->BSIM4v6csdb = 0.0;
              here->BSIM4v6cgbb = here->BSIM4v6csbb = here->BSIM4v6cdbb = here->BSIM4v6cbbb = 0.0;
              here->BSIM4v6cqdb = here->BSIM4v6cqsb = here->BSIM4v6cqgb 
                              = here->BSIM4v6cqbb = 0.0;
              here->BSIM4v6gtau = 0.0;
              goto finished;
          }
          else if (model->BSIM4v6capMod == 0)
          {
              if (Vbseff < 0.0)
              {   VbseffCV = Vbs; /*4.6.2*/
                  dVbseffCV_dVb = 1.0;
              }
              else
              {   VbseffCV = pParam->BSIM4v6phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb; /*4.6.2*/
              }

              Vfb = pParam->BSIM4v6vfbcv;
              Vth = Vfb + pParam->BSIM4v6phi + pParam->BSIM4v6k1ox * sqrtPhis; 
              Vgst = Vgs_eff - Vth;
              dVth_dVb = pParam->BSIM4v6k1ox * dsqrtPhis_dVb *dVbseff_dVb; /*4.6.2*/
              dVgst_dVb = -dVth_dVb;
              dVgst_dVg = dVgs_eff_dVg; 

              CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV
                    * pParam->BSIM4v6leffCV * here->BSIM4v6nf;
              Arg1 = Vgs_eff - VbseffCV - Vfb;

              if (Arg1 <= 0.0)
              {   qgate = CoxWL * Arg1;
                  qbulk = -qgate;
                  qdrn = 0.0;

                  here->BSIM4v6cggb = CoxWL * dVgs_eff_dVg;
                  here->BSIM4v6cgdb = 0.0;
                  here->BSIM4v6cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

                  here->BSIM4v6cdgb = 0.0;
                  here->BSIM4v6cddb = 0.0;
                  here->BSIM4v6cdsb = 0.0;

                  here->BSIM4v6cbgb = -CoxWL * dVgs_eff_dVg;
                  here->BSIM4v6cbdb = 0.0;
                  here->BSIM4v6cbsb = -here->BSIM4v6cgsb;
              } /* Arg1 <= 0.0, end of accumulation */
              else if (Vgst <= 0.0)
              {   T1 = 0.5 * pParam->BSIM4v6k1ox;
                  T2 = sqrt(T1 * T1 + Arg1);
                  qgate = CoxWL * pParam->BSIM4v6k1ox * (T2 - T1);
                  qbulk = -qgate;
                  qdrn = 0.0;

                  T0 = CoxWL * T1 / T2;
                  here->BSIM4v6cggb = T0 * dVgs_eff_dVg;
                  here->BSIM4v6cgdb = 0.0;
                  here->BSIM4v6cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);
   
                  here->BSIM4v6cdgb = 0.0;
                  here->BSIM4v6cddb = 0.0;
                  here->BSIM4v6cdsb = 0.0;

                  here->BSIM4v6cbgb = -here->BSIM4v6cggb;
                  here->BSIM4v6cbdb = 0.0;
                  here->BSIM4v6cbsb = -here->BSIM4v6cgsb;
              } /* Vgst <= 0.0, end of depletion */
              else
              {   One_Third_CoxWL = CoxWL / 3.0;
                  Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                  AbulkCV = Abulk0 * pParam->BSIM4v6abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v6abulkCVfactor * dAbulk0_dVb*dVbseff_dVb;
                  
                  dVdsat_dVg = 1.0 / AbulkCV;  /*4.6.2*/
                          Vdsat = Vgst * dVdsat_dVg;
                  dVdsat_dVb = - (Vdsat * dAbulkCV_dVb + dVth_dVb)* dVdsat_dVg; 

                  if (model->BSIM4v6xpart > 0.5)
                  {   /* 0/100 Charge partition model */
                      if (Vdsat <= Vds)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v6phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.0;

                          here->BSIM4v6cggb = One_Third_CoxWL * (3.0
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb + T2);
                          here->BSIM4v6cgdb = 0.0;
       
                          here->BSIM4v6cdgb = 0.0;
                          here->BSIM4v6cddb = 0.0;
                          here->BSIM4v6cdsb = 0.0;

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v6cbsb = -(here->BSIM4v6cbgb + T3);
                          here->BSIM4v6cbdb = 0.0;
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
                                - pParam->BSIM4v6phi - 0.5 * (Vds - T3));
                          T10 = T4 * T8;
                          qdrn = T4 * T7;
                          qbulk = -(qgate + qdrn + T10);
  
                          T5 = T3 / T1;
                          here->BSIM4v6cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          T11 = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v6cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb + T11
                                          + here->BSIM4v6cgdb);
                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                          T7 = T9 * T7;
                          T8 = T9 * T8;
                          T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                          here->BSIM4v6cdgb = (T7 * dAlphaz_dVg - T9
                                          * dVdsat_dVg) * dVgs_eff_dVg;
                          T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          here->BSIM4v6cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                          here->BSIM4v6cdsb = -(here->BSIM4v6cdgb + T12
                                          + here->BSIM4v6cddb);

                          T9 = 2.0 * T4 * (1.0 + T5);
                          T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg)
                              * dVgs_eff_dVg;
                          T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                          T12 = T4 * (2.0 * T2 + T5 - 1.0); 
                          T0 = -(T10 + T11 + T12);

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb
                                          + here->BSIM4v6cdgb + T10);
                          here->BSIM4v6cbdb = -(here->BSIM4v6cgdb 
                                          + here->BSIM4v6cddb + T12);
                          here->BSIM4v6cbsb = -(here->BSIM4v6cgsb
                                          + here->BSIM4v6cdsb + T0);
                      }
                  }
                  else if (model->BSIM4v6xpart < 0.5)
                  {   /* 40/60 Charge partition model */
                      if (Vds >= Vdsat)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v6phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.4 * T2;

                          here->BSIM4v6cggb = One_Third_CoxWL * (3.0 
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb + T2);
                          here->BSIM4v6cgdb = 0.0;
       
                          T3 = 0.4 * Two_Third_CoxWL;
                          here->BSIM4v6cdgb = -T3 * dVgs_eff_dVg;
                          here->BSIM4v6cddb = 0.0;
                          T4 = T3 * dVth_dVb;
                          here->BSIM4v6cdsb = -(T4 + here->BSIM4v6cdgb);

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb 
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v6cbsb = -(here->BSIM4v6cbgb + T3);
                          here->BSIM4v6cbdb = 0.0;
                      }
                      else
                      {   /* linear region  */
                          Alphaz = Vgst / Vdsat;
                          T1 = 2.0 * Vdsat - Vds;
                          T2 = Vds / (3.0 * T1);
                          T3 = T2 * Vds;
                          T9 = 0.25 * CoxWL;
                          T4 = T9 * Alphaz;
                          qgate = CoxWL * (Vgs_eff - Vfb - pParam->BSIM4v6phi
                                - 0.5 * (Vds - T3));

                          T5 = T3 / T1;
                          here->BSIM4v6cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v6cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb 
                                          + here->BSIM4v6cgdb + tmp);

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

                          here->BSIM4v6cdgb = (T7 * dAlphaz_dVg - tmp1
                                          * dVdsat_dVg) * dVgs_eff_dVg;
                          T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                          here->BSIM4v6cddb = T4 * (2.0 - (1.0 / (3.0 * T1
                                          * T1) + 2.0 * tmp) * T6 + T8
                                          * (6.0 * Vdsat - 2.4 * Vds));
                          here->BSIM4v6cdsb = -(here->BSIM4v6cdgb 
                                          + T10 + here->BSIM4v6cddb);

                          T7 = 2.0 * (T1 + T3);
                          qbulk = -(qgate - T4 * T7);
                          T7 *= T9;
                          T0 = 4.0 * T4 * (1.0 - T5);
                          T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg
                                  - here->BSIM4v6cdgb;  /*4.6.2*/
                          T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                          T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5) 
                              - here->BSIM4v6cddb;
                          tmp = -(T10 + T11 + T12);

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb 
                                          + here->BSIM4v6cdgb + T12);
                          here->BSIM4v6cbdb = -(here->BSIM4v6cgdb
                                          + here->BSIM4v6cddb + T10);  
                          here->BSIM4v6cbsb = -(here->BSIM4v6cgsb
                                          + here->BSIM4v6cdsb + tmp);
                      }
                  }
                  else
                  {   /* 50/50 partitioning */
                      if (Vds >= Vdsat)
                      {   /* saturation region */
                          T1 = Vdsat / 3.0;
                          qgate = CoxWL * (Vgs_eff - Vfb
                                - pParam->BSIM4v6phi - T1);
                          T2 = -Two_Third_CoxWL * Vgst;
                          qbulk = -(qgate + T2);
                          qdrn = 0.5 * T2;

                          here->BSIM4v6cggb = One_Third_CoxWL * (3.0
                                          - dVdsat_dVg) * dVgs_eff_dVg;
                          T2 = -One_Third_CoxWL * dVdsat_dVb;
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb + T2);
                          here->BSIM4v6cgdb = 0.0;
       
                          here->BSIM4v6cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                          here->BSIM4v6cddb = 0.0;
                          T4 = One_Third_CoxWL * dVth_dVb;
                          here->BSIM4v6cdsb = -(T4 + here->BSIM4v6cdgb);

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb 
                                          - Two_Third_CoxWL * dVgs_eff_dVg);
                          T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                          here->BSIM4v6cbsb = -(here->BSIM4v6cbgb + T3);
                          here->BSIM4v6cbdb = 0.0;
                      }
                      else
                      {   /* linear region */
                          Alphaz = Vgst / Vdsat;
                          T1 = 2.0 * Vdsat - Vds;
                          T2 = Vds / (3.0 * T1);
                          T3 = T2 * Vds;
                          T9 = 0.25 * CoxWL;
                          T4 = T9 * Alphaz;
                          qgate = CoxWL * (Vgs_eff - Vfb - pParam->BSIM4v6phi
                                - 0.5 * (Vds - T3));

                          T5 = T3 / T1;
                          here->BSIM4v6cggb = CoxWL * (1.0 - T5 * dVdsat_dVg)
                                          * dVgs_eff_dVg;
                          tmp = -CoxWL * T5 * dVdsat_dVb;
                          here->BSIM4v6cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                          here->BSIM4v6cgsb = -(here->BSIM4v6cggb 
                                          + here->BSIM4v6cgdb + tmp);

                          T6 = 1.0 / Vdsat;
                          dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                          dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                          T7 = T1 + T3;
                          qdrn = -T4 * T7;
                          qbulk = - (qgate + qdrn + qdrn);
                          T7 *= T9;
                          T0 = T4 * (2.0 * T5 - 2.0);

                          here->BSIM4v6cdgb = (T0 * dVdsat_dVg - T7
                                          * dAlphaz_dVg) * dVgs_eff_dVg;
                          T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                          here->BSIM4v6cddb = T4 * (1.0 - 2.0 * T2 - T5);
                          here->BSIM4v6cdsb = -(here->BSIM4v6cdgb + T12
                                          + here->BSIM4v6cddb);

                          here->BSIM4v6cbgb = -(here->BSIM4v6cggb
                                          + 2.0 * here->BSIM4v6cdgb);
                          here->BSIM4v6cbdb = -(here->BSIM4v6cgdb
                                          + 2.0 * here->BSIM4v6cddb);
                          here->BSIM4v6cbsb = -(here->BSIM4v6cgsb
                                          + 2.0 * here->BSIM4v6cdsb);
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
              {   VbseffCV = pParam->BSIM4v6phi - Phis;
                  dVbseffCV_dVb = -dPhis_dVb;
              }

              CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV
                    * pParam->BSIM4v6leffCV * here->BSIM4v6nf;

              if(model->BSIM4v6cvchargeMod == 0)
                {
                  /* Seperate VgsteffCV with noff and voffcv */
                  noff = n * pParam->BSIM4v6noff;
                  dnoff_dVd = pParam->BSIM4v6noff * dn_dVd;
                  dnoff_dVb = pParam->BSIM4v6noff * dn_dVb;
                  T0 = Vtm * noff;
                  voffcv = pParam->BSIM4v6voffcv;
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
                  T1 = pParam->BSIM4v6mstarcv * Vgst;
                  T2 = T1 / T0;
                  if (T2 > EXP_THRESHOLD)
                    {   
                      T10 = T1;
                      dT10_dVg = pParam->BSIM4v6mstarcv * dVgs_eff_dVg;
                      dT10_dVd = -dVth_dVd * pParam->BSIM4v6mstarcv;
                      dT10_dVb = -dVth_dVb * pParam->BSIM4v6mstarcv;
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
                      dT10_dVg = pParam->BSIM4v6mstarcv * ExpVgst / (1.0 + ExpVgst);
                      dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
                      dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
                      dT10_dVg *= dVgs_eff_dVg;
                    }
                  
                  T1 = pParam->BSIM4v6voffcbncv - (1.0 - pParam->BSIM4v6mstarcv) * Vgst;
                  T2 = T1 / T0;
                  if (T2 < -EXP_THRESHOLD)
                    {   
                      T3 = model->BSIM4v6coxe * MIN_EXP / pParam->BSIM4v6cdep0;
                      T9 = pParam->BSIM4v6mstarcv + T3 * n;
                      dT9_dVg = 0.0;
                      dT9_dVd = dn_dVd * T3;
                      dT9_dVb = dn_dVb * T3;
                    }
                  else if (T2 > EXP_THRESHOLD)
                    {   
                      T3 = model->BSIM4v6coxe * MAX_EXP / pParam->BSIM4v6cdep0;
                      T9 = pParam->BSIM4v6mstarcv + T3 * n;
                      dT9_dVg = 0.0;
                      dT9_dVd = dn_dVd * T3;
                      dT9_dVb = dn_dVb * T3;
                    }
                  else
                    {   
                      ExpVgst = exp(T2);
                      T3 = model->BSIM4v6coxe / pParam->BSIM4v6cdep0;
                      T4 = T3 * ExpVgst;
                      T5 = T1 * T4 / T0;
                      T9 = pParam->BSIM4v6mstarcv + n * T4;
                      dT9_dVg = T3 * (pParam->BSIM4v6mstarcv - 1.0) * ExpVgst / Vtm;
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
          

              if (model->BSIM4v6capMod == 1)
              {   Vfb = here->BSIM4v6vfbzb;
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

                  T0 = 0.5 * pParam->BSIM4v6k1ox;
                  T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam->BSIM4v6k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
                  else if (T3 < 0.0)
                  {   T1 = T0 + T3 / pParam->BSIM4v6k1ox;
                      T2 = CoxWL;
                  }
                  else
                  {   T1 = sqrt(T0 * T0 + T3);
                      T2 = CoxWL * T0 / T1;
                  }

                  Qsub0 = CoxWL * pParam->BSIM4v6k1ox * (T1 - T0);

                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb 
                             + dVgsteff_dVb);

                  AbulkCV = Abulk0 * pParam->BSIM4v6abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v6abulkCVfactor * dAbulk0_dVb;
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
                  
                  if (model->BSIM4v6xpart > 0.5)
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
                  else if (model->BSIM4v6xpart < 0.5)
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
                  
                  here->BSIM4v6cggb = Cgg;
                  here->BSIM4v6cgsb = -(Cgg + Cgd + Cgb);
                  here->BSIM4v6cgdb = Cgd;
                  here->BSIM4v6cdgb = -(Cgg + Cbg + Csg);
                  here->BSIM4v6cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                                     + Csg + Csd + Csb);
                  here->BSIM4v6cddb = -(Cgd + Cbd + Csd);
                  here->BSIM4v6cbgb = Cbg;
                  here->BSIM4v6cbsb = -(Cbg + Cbd + Cbb);
                  here->BSIM4v6cbdb = Cbd;
              }
              
              /* Charge-Thickness capMod (CTM) begins */
              else if (model->BSIM4v6capMod == 2)
              {   V3 = here->BSIM4v6vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                  if (here->BSIM4v6vfbzb <= 0.0)
                      T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * here->BSIM4v6vfbzb);
                  else
                      T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * here->BSIM4v6vfbzb);

                  T1 = 0.5 * (1.0 + V3 / T0);
                  Vfbeff = here->BSIM4v6vfbzb - 0.5 * (V3 + T0);
                  dVfbeff_dVg = T1 * dVgs_eff_dVg;
                  dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                  Cox = model->BSIM4v6coxp;
                  Tox = 1.0e8 * model->BSIM4v6toxp;
                  T0 = (Vgs_eff - VbseffCV - here->BSIM4v6vfbzb) / Tox;
                  dT0_dVg = dVgs_eff_dVg / Tox;
                  dT0_dVb = -dVbseffCV_dVb / Tox;

                  tmp = T0 * pParam->BSIM4v6acde;
                  if ((-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD))
                  {   Tcen = pParam->BSIM4v6ldeb * exp(tmp);
                      dTcen_dVg = pParam->BSIM4v6acde * Tcen;
                      dTcen_dVb = dTcen_dVg * dT0_dVb;
                      dTcen_dVg *= dT0_dVg;
                  }
                  else if (tmp <= -EXP_THRESHOLD)
                  {   Tcen = pParam->BSIM4v6ldeb * MIN_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }
                  else
                  {   Tcen = pParam->BSIM4v6ldeb * MAX_EXP;
                      dTcen_dVg = dTcen_dVb = 0.0;
                  }

                  LINK = 1.0e-3 * model->BSIM4v6toxp;
                  V3 = pParam->BSIM4v6ldeb - Tcen - LINK;
                  V4 = sqrt(V3 * V3 + 4.0 * LINK * pParam->BSIM4v6ldeb);
                  Tcen = pParam->BSIM4v6ldeb - 0.5 * (V3 + V4);
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
                  CoxWLcen = CoxWL * Coxeff / model->BSIM4v6coxe;

                  Qac0 = CoxWLcen * (Vfbeff - here->BSIM4v6vfbzb);
                  QovCox = Qac0 / Coxeff;
                  dQac0_dVg = CoxWLcen * dVfbeff_dVg
                            + QovCox * dCoxeff_dVg;
                  dQac0_dVb = CoxWLcen * dVfbeff_dVb 
                            + QovCox * dCoxeff_dVb;

                  T0 = 0.5 * pParam->BSIM4v6k1ox;
                  T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                  if (pParam->BSIM4v6k1ox == 0.0)
                  {   T1 = 0.0;
                      T2 = 0.0;
                  }
                  else if (T3 < 0.0)
                  {   T1 = T0 + T3 / pParam->BSIM4v6k1ox;
                      T2 = CoxWLcen;
                  }
                  else
                  {   T1 = sqrt(T0 * T0 + T3);
                      T2 = CoxWLcen * T0 / T1;
                  }

                  Qsub0 = CoxWLcen * pParam->BSIM4v6k1ox * (T1 - T0);
                  QovCox = Qsub0 / Coxeff;
                  dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg)
                             + QovCox * dCoxeff_dVg;
                  dQsub0_dVd = -T2 * dVgsteff_dVd;
                  dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb)
                             + QovCox * dCoxeff_dVb;

                  /* Gate-bias dependent delta Phis begins */
                  if (pParam->BSIM4v6k1ox <= 0.0)
                  {   Denomi = 0.25 * pParam->BSIM4v6moin * Vtm;
                      T0 = 0.5 * pParam->BSIM4v6sqrtPhi;
                  }
                  else
                  {   Denomi = pParam->BSIM4v6moin * Vtm 
                             * pParam->BSIM4v6k1ox * pParam->BSIM4v6k1ox;
                      T0 = pParam->BSIM4v6k1ox * pParam->BSIM4v6sqrtPhi;
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
                  T0 = (Vgsteff + here->BSIM4v6vtfbphi2) / Tox;
                  tmp = exp(model->BSIM4v6bdos * 0.7 * log(T0));
                  T1 = 1.0 + tmp;
                  T2 = model->BSIM4v6bdos * 0.7 * tmp / (T0 * Tox);
                  Tcen = model->BSIM4v6ados * 1.9e-9 / T1;
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
                  CoxWLcen = CoxWL * Coxeff / model->BSIM4v6coxe;

                  AbulkCV = Abulk0 * pParam->BSIM4v6abulkCVfactor;
                  dAbulkCV_dVb = pParam->BSIM4v6abulkCVfactor * dAbulk0_dVb;
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

                  if (model->BSIM4v6xpart > 0.5)
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
                  else if (model->BSIM4v6xpart < 0.5)
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

                  here->BSIM4v6cggb = Cgg;
                  here->BSIM4v6cgsb = -(Cgg + Cgd + Cgb);
                  here->BSIM4v6cgdb = Cgd;
                  here->BSIM4v6cdgb = -(Cgg + Cbg + Csg);
                  here->BSIM4v6cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb
                                  + Csg + Csd + Csb);
                  here->BSIM4v6cddb = -(Cgd + Cbd + Csd);
                  here->BSIM4v6cbgb = Cbg;
                  here->BSIM4v6cbsb = -(Cbg + Cbd + Cbb);
                  here->BSIM4v6cbdb = Cbd;
              }  /* End of CTM */
          }

          here->BSIM4v6csgb = - here->BSIM4v6cggb - here->BSIM4v6cdgb - here->BSIM4v6cbgb;
          here->BSIM4v6csdb = - here->BSIM4v6cgdb - here->BSIM4v6cddb - here->BSIM4v6cbdb;
          here->BSIM4v6cssb = - here->BSIM4v6cgsb - here->BSIM4v6cdsb - here->BSIM4v6cbsb;
          here->BSIM4v6cgbb = - here->BSIM4v6cgdb - here->BSIM4v6cggb - here->BSIM4v6cgsb;
          here->BSIM4v6cdbb = - here->BSIM4v6cddb - here->BSIM4v6cdgb - here->BSIM4v6cdsb;
          here->BSIM4v6cbbb = - here->BSIM4v6cbgb - here->BSIM4v6cbdb - here->BSIM4v6cbsb;
          here->BSIM4v6csbb = - here->BSIM4v6cgbb - here->BSIM4v6cdbb - here->BSIM4v6cbbb;
          here->BSIM4v6qgate = qgate;
          here->BSIM4v6qbulk = qbulk;
          here->BSIM4v6qdrn = qdrn;
          here->BSIM4v6qsrc = -(qgate + qbulk + qdrn);

          /* NQS begins */
          if ((here->BSIM4v6trnqsMod) || (here->BSIM4v6acnqsMod))
          {   here->BSIM4v6qchqs = qcheq = -(qbulk + qgate);
              here->BSIM4v6cqgb = -(here->BSIM4v6cggb + here->BSIM4v6cbgb);
              here->BSIM4v6cqdb = -(here->BSIM4v6cgdb + here->BSIM4v6cbdb);
              here->BSIM4v6cqsb = -(here->BSIM4v6cgsb + here->BSIM4v6cbsb);
              here->BSIM4v6cqbb = -(here->BSIM4v6cqgb + here->BSIM4v6cqdb
                              + here->BSIM4v6cqsb);

              CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV * here->BSIM4v6nf
                    * pParam->BSIM4v6leffCV;
              T1 = here->BSIM4v6gcrg / CoxWL; /* 1 / tau */
              here->BSIM4v6gtau = T1 * ScalingFactor;

              if (here->BSIM4v6acnqsMod)
                  here->BSIM4v6taunet = 1.0 / T1;

              *(ckt->CKTstate0 + here->BSIM4v6qcheq) = qcheq;
              if (ckt->CKTmode & MODEINITTRAN)
                  *(ckt->CKTstate1 + here->BSIM4v6qcheq) =
                                   *(ckt->CKTstate0 + here->BSIM4v6qcheq);
              if (here->BSIM4v6trnqsMod)
              {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qcheq);
                  if (error)
                      return(error);
              }
          }


finished: 

          /* Calculate junction C-V */
          if (ChargeComputationNeeded)
          {   czbd = model->BSIM4v6DunitAreaTempJctCap * here->BSIM4v6Adeff; /* bug fix */
              czbs = model->BSIM4v6SunitAreaTempJctCap * here->BSIM4v6Aseff;
              czbdsw = model->BSIM4v6DunitLengthSidewallTempJctCap * here->BSIM4v6Pdeff;
              czbdswg = model->BSIM4v6DunitLengthGateSidewallTempJctCap
                      * pParam->BSIM4v6weffCJ * here->BSIM4v6nf;
              czbssw = model->BSIM4v6SunitLengthSidewallTempJctCap * here->BSIM4v6Pseff;
              czbsswg = model->BSIM4v6SunitLengthGateSidewallTempJctCap
                      * pParam->BSIM4v6weffCJ * here->BSIM4v6nf;

              MJS = model->BSIM4v6SbulkJctBotGradingCoeff;
              MJSWS = model->BSIM4v6SbulkJctSideGradingCoeff;
              MJSWGS = model->BSIM4v6SbulkJctGateSideGradingCoeff;

              MJD = model->BSIM4v6DbulkJctBotGradingCoeff;
              MJSWD = model->BSIM4v6DbulkJctSideGradingCoeff;
              MJSWGD = model->BSIM4v6DbulkJctGateSideGradingCoeff;

              /* Source Bulk Junction */
              if (vbs_jct == 0.0)
              {   *(ckt->CKTstate0 + here->BSIM4v6qbs) = 0.0;
                  here->BSIM4v6capbs = czbs + czbssw + czbsswg;
              }
              else if (vbs_jct < 0.0)
              {   if (czbs > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v6PhiBS;
                      if (MJS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbs) = model->BSIM4v6PhiBS * czbs 
                                       * (1.0 - arg * sarg) / (1.0 - MJS);
                      here->BSIM4v6capbs = czbs * sarg;
                  }
                  else
                  {   *(ckt->CKTstate0 + here->BSIM4v6qbs) = 0.0;
                      here->BSIM4v6capbs = 0.0;
                  }
                  if (czbssw > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v6PhiBSWS;
                      if (MJSWS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbs) += model->BSIM4v6PhiBSWS * czbssw
                                       * (1.0 - arg * sarg) / (1.0 - MJSWS);
                      here->BSIM4v6capbs += czbssw * sarg;
                  }
                  if (czbsswg > 0.0)
                  {   arg = 1.0 - vbs_jct / model->BSIM4v6PhiBSWGS;
                      if (MJSWGS == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWGS * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbs) += model->BSIM4v6PhiBSWGS * czbsswg
                                       * (1.0 - arg * sarg) / (1.0 - MJSWGS);
                      here->BSIM4v6capbs += czbsswg * sarg;
                  }

              }
              else
              {   T0 = czbs + czbssw + czbsswg;
                  T1 = vbs_jct * (czbs * MJS / model->BSIM4v6PhiBS + czbssw * MJSWS 
                     / model->BSIM4v6PhiBSWS + czbsswg * MJSWGS / model->BSIM4v6PhiBSWGS);    
                  *(ckt->CKTstate0 + here->BSIM4v6qbs) = vbs_jct * (T0 + 0.5 * T1);
                  here->BSIM4v6capbs = T0 + T1;
              }

              /* Drain Bulk Junction */
              if (vbd_jct == 0.0)
              {   *(ckt->CKTstate0 + here->BSIM4v6qbd) = 0.0;
                  here->BSIM4v6capbd = czbd + czbdsw + czbdswg;
              }
              else if (vbd_jct < 0.0)
              {   if (czbd > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v6PhiBD;
                      if (MJD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbd) = model->BSIM4v6PhiBD* czbd 
                                       * (1.0 - arg * sarg) / (1.0 - MJD);
                      here->BSIM4v6capbd = czbd * sarg;
                  }
                  else
                  {   *(ckt->CKTstate0 + here->BSIM4v6qbd) = 0.0;
                      here->BSIM4v6capbd = 0.0;
                  }
                  if (czbdsw > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v6PhiBSWD;
                      if (MJSWD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbd) += model->BSIM4v6PhiBSWD * czbdsw 
                                       * (1.0 - arg * sarg) / (1.0 - MJSWD);
                      here->BSIM4v6capbd += czbdsw * sarg;
                  }
                  if (czbdswg > 0.0)
                  {   arg = 1.0 - vbd_jct / model->BSIM4v6PhiBSWGD;
                      if (MJSWGD == 0.5)
                          sarg = 1.0 / sqrt(arg);
                      else
                          sarg = exp(-MJSWGD * log(arg));
                      *(ckt->CKTstate0 + here->BSIM4v6qbd) += model->BSIM4v6PhiBSWGD * czbdswg
                                       * (1.0 - arg * sarg) / (1.0 - MJSWGD);
                      here->BSIM4v6capbd += czbdswg * sarg;
                  }
              }
              else
              {   T0 = czbd + czbdsw + czbdswg;
                  T1 = vbd_jct * (czbd * MJD / model->BSIM4v6PhiBD + czbdsw * MJSWD
                     / model->BSIM4v6PhiBSWD + czbdswg * MJSWGD / model->BSIM4v6PhiBSWGD);
                  *(ckt->CKTstate0 + here->BSIM4v6qbd) = vbd_jct * (T0 + 0.5 * T1);
                  here->BSIM4v6capbd = T0 + T1; 
              }
          }


          /*
           *  check convergence
           */

          if ((here->BSIM4v6off == 0) || (!(ckt->CKTmode & MODEINITFIX)))
          {   if (Check == 1)
              {   ckt->CKTnoncon++;
#ifndef NEWCONV
              } 
              else
              {   if (here->BSIM4v6mode >= 0)
                  {   Idtot = here->BSIM4v6cd + here->BSIM4v6csub
                            + here->BSIM4v6Igidl - here->BSIM4v6cbd;
                  }
                  else
                  {   Idtot = here->BSIM4v6cd + here->BSIM4v6cbd - here->BSIM4v6Igidl; /* bugfix */
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
                  {   Ibtot = here->BSIM4v6cbs + here->BSIM4v6cbd
                            - here->BSIM4v6Igidl - here->BSIM4v6Igisl - here->BSIM4v6csub;
                      tol6 = ckt->CKTreltol * MAX(fabs(cbhat), fabs(Ibtot))
                          + ckt->CKTabstol;
                      if (fabs(cbhat - Ibtot) > tol6)
                      {   ckt->CKTnoncon++;
                      }
                  }
#endif /* NEWCONV */
              }
          }
          *(ckt->CKTstate0 + here->BSIM4v6vds) = vds;
          *(ckt->CKTstate0 + here->BSIM4v6vgs) = vgs;
          *(ckt->CKTstate0 + here->BSIM4v6vbs) = vbs;
          *(ckt->CKTstate0 + here->BSIM4v6vbd) = vbd;
          *(ckt->CKTstate0 + here->BSIM4v6vges) = vges;
          *(ckt->CKTstate0 + here->BSIM4v6vgms) = vgms;
          *(ckt->CKTstate0 + here->BSIM4v6vdbs) = vdbs;
          *(ckt->CKTstate0 + here->BSIM4v6vdbd) = vdbd;
          *(ckt->CKTstate0 + here->BSIM4v6vsbs) = vsbs;
          *(ckt->CKTstate0 + here->BSIM4v6vses) = vses;
          *(ckt->CKTstate0 + here->BSIM4v6vdes) = vdes;
          *(ckt->CKTstate0 + here->BSIM4v6qdef) = qdef;


          if (!ChargeComputationNeeded)
              goto line850; 

          if (here->BSIM4v6rgateMod == 3) 
          {   
                  vgdx = vgmd; 
                  vgsx = vgms;
          }  
          else  /* For rgateMod == 0, 1 and 2 */
          {
                  vgdx = vgd;
                  vgsx = vgs;
          }
          if (model->BSIM4v6capMod == 0) 
          {  
                  cgdo = pParam->BSIM4v6cgdo; 
                  qgdo = pParam->BSIM4v6cgdo * vgdx;
                  cgso = pParam->BSIM4v6cgso;
                  qgso = pParam->BSIM4v6cgso * vgsx;
          }
          else /* For both capMod == 1 and 2 */
          {   T0 = vgdx + DELTA_1;
              T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
              T2 = 0.5 * (T0 - T1);

              T3 = pParam->BSIM4v6weffCV * pParam->BSIM4v6cgdl;
              T4 = sqrt(1.0 - 4.0 * T2 / pParam->BSIM4v6ckappad); 
              cgdo = pParam->BSIM4v6cgdo + T3 - T3 * (1.0 - 1.0 / T4)
                   * (0.5 - 0.5 * T0 / T1);
              qgdo = (pParam->BSIM4v6cgdo + T3) * vgdx - T3 * (T2
                   + 0.5 * pParam->BSIM4v6ckappad * (T4 - 1.0));

              T0 = vgsx + DELTA_1;
              T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
              T2 = 0.5 * (T0 - T1);
              T3 = pParam->BSIM4v6weffCV * pParam->BSIM4v6cgsl;
              T4 = sqrt(1.0 - 4.0 * T2 / pParam->BSIM4v6ckappas);
              cgso = pParam->BSIM4v6cgso + T3 - T3 * (1.0 - 1.0 / T4)
                   * (0.5 - 0.5 * T0 / T1);
              qgso = (pParam->BSIM4v6cgso + T3) * vgsx - T3 * (T2
                   + 0.5 * pParam->BSIM4v6ckappas * (T4 - 1.0));
          }

          if (here->BSIM4v6nf != 1.0)
          {   cgdo *= here->BSIM4v6nf;
              cgso *= here->BSIM4v6nf;
              qgdo *= here->BSIM4v6nf;
              qgso *= here->BSIM4v6nf;
          }        
          here->BSIM4v6cgdo = cgdo;
          here->BSIM4v6qgdo = qgdo;
          here->BSIM4v6cgso = cgso;
          here->BSIM4v6qgso = qgso;

#ifndef NOBYPASS
line755:
#endif
          ag0 = ckt->CKTag[0];
          if (here->BSIM4v6mode > 0)
          {   if (here->BSIM4v6trnqsMod == 0)
              {   qdrn -= qgdo;
                  if (here->BSIM4v6rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v6cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcggb = here->BSIM4v6cggb * ag0;
                      gcgdb = here->BSIM4v6cgdb * ag0;
                      gcgsb = here->BSIM4v6cgsb * ag0;   
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = here->BSIM4v6cdgb * ag0;
                      gcsgb = -(here->BSIM4v6cggb + here->BSIM4v6cbgb
                            + here->BSIM4v6cdgb) * ag0;
                      gcbgb = here->BSIM4v6cbgb * ag0;

                      qgmb = pParam->BSIM4v6cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qbulk -= qgmb;
                      qsrc = -(qgate + qgmid + qbulk + qdrn);
                  }
                    else
                  {   gcggb = (here->BSIM4v6cggb + cgdo + cgso
                            + pParam->BSIM4v6cgbo ) * ag0;
                      gcgdb = (here->BSIM4v6cgdb - cgdo) * ag0;
                      gcgsb = (here->BSIM4v6cgsb - cgso) * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = (here->BSIM4v6cdgb - cgdo) * ag0;
                      gcsgb = -(here->BSIM4v6cggb + here->BSIM4v6cbgb
                            + here->BSIM4v6cdgb + cgso) * ag0;
                      gcbgb = (here->BSIM4v6cbgb - pParam->BSIM4v6cgbo) * ag0;

                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v6cgbo * vgb;
                      qgate += qgdo + qgso + qgb;
                      qbulk -= qgb;
                      qsrc = -(qgate + qbulk + qdrn);
                  }
                  gcddb = (here->BSIM4v6cddb + here->BSIM4v6capbd + cgdo) * ag0;
                  gcdsb = here->BSIM4v6cdsb * ag0;

                  gcsdb = -(here->BSIM4v6cgdb + here->BSIM4v6cbdb
                        + here->BSIM4v6cddb) * ag0;
                  gcssb = (here->BSIM4v6capbs + cgso - (here->BSIM4v6cgsb
                        + here->BSIM4v6cbsb + here->BSIM4v6cdsb)) * ag0;

                  if (!here->BSIM4v6rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                      gcbdb = (here->BSIM4v6cbdb - here->BSIM4v6capbd) * ag0;
                      gcbsb = (here->BSIM4v6cbsb - here->BSIM4v6capbs) * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb  = -(here->BSIM4v6cddb + here->BSIM4v6cdgb 
                             + here->BSIM4v6cdsb) * ag0;
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb)
                            + here->BSIM4v6capbs * ag0;
                      gcbdb = here->BSIM4v6cbdb * ag0;
                      gcbsb = here->BSIM4v6cbsb * ag0;

                      gcdbdb = -here->BSIM4v6capbd * ag0;
                      gcsbsb = -here->BSIM4v6capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);

                  ggtg = ggtd = ggtb = ggts = 0.0;
                  sxpart = 0.6;
                  dxpart = 0.4;
                  ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                  dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
              }
              else
              {   qcheq = here->BSIM4v6qchqs;
                  CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV * here->BSIM4v6nf
                        * pParam->BSIM4v6leffCV;
                  T0 = qdef * ScalingFactor / CoxWL;

                  ggtg = here->BSIM4v6gtg = T0 * here->BSIM4v6gcrgg;
                  ggtd = here->BSIM4v6gtd = T0 * here->BSIM4v6gcrgd;
                  ggts = here->BSIM4v6gts = T0 * here->BSIM4v6gcrgs;
                  ggtb = here->BSIM4v6gtb = T0 * here->BSIM4v6gcrgb;
                  gqdef = ScalingFactor * ag0;

                  gcqgb = here->BSIM4v6cqgb * ag0;
                  gcqdb = here->BSIM4v6cqdb * ag0;
                  gcqsb = here->BSIM4v6cqsb * ag0;
                  gcqbb = here->BSIM4v6cqbb * ag0;

                  if (fabs(qcheq) <= 1.0e-5 * CoxWL)
                  {   if (model->BSIM4v6xpart < 0.5)
                      {   dxpart = 0.4;
                      }
                      else if (model->BSIM4v6xpart > 0.5)
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
                      Cdd = here->BSIM4v6cddb;
                      Csd = -(here->BSIM4v6cgdb + here->BSIM4v6cddb
                          + here->BSIM4v6cbdb);
                      ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                      Cdg = here->BSIM4v6cdgb;
                      Csg = -(here->BSIM4v6cggb + here->BSIM4v6cdgb
                          + here->BSIM4v6cbgb);
                      ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                      Cds = here->BSIM4v6cdsb;
                      Css = -(here->BSIM4v6cgsb + here->BSIM4v6cdsb
                          + here->BSIM4v6cbsb);
                      ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                  }
                  sxpart = 1.0 - dxpart;
                  dsxpart_dVd = -ddxpart_dVd;
                  dsxpart_dVg = -ddxpart_dVg;
                  dsxpart_dVs = -ddxpart_dVs;
                  dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                  if (here->BSIM4v6rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v6cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcdgb = gcsgb = gcbgb = 0.0;
                      gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                      qgmb = pParam->BSIM4v6cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qgate = 0.0;
                      qbulk = -qgmb;
                      qdrn = -qgdo;
                      qsrc = -(qgmid + qbulk + qdrn);
                  }
                  else
                  {   gcggb = (cgdo + cgso + pParam->BSIM4v6cgbo ) * ag0;
                      gcgdb = -cgdo * ag0;
                      gcgsb = -cgso * ag0;
                      gcgbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgb = gcgdb;
                      gcsgb = gcgsb;
                      gcbgb = gcgbb;
                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v6cgbo * vgb;
                      qgate = qgdo + qgso + qgb;
                      qbulk = -qgb;
                      qdrn = -qgdo;
                      qsrc = -(qgate + qbulk + qdrn);
                  }

                  gcddb = (here->BSIM4v6capbd + cgdo) * ag0;
                  gcdsb = gcsdb = 0.0;
                  gcssb = (here->BSIM4v6capbs + cgso) * ag0;

                  if (!here->BSIM4v6rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                      gcsbb = -(gcsgb + gcssb + gcsgmb);
                      gcbdb = -here->BSIM4v6capbd * ag0;
                      gcbsb = -here->BSIM4v6capbs * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                      gcdbdb = -here->BSIM4v6capbd * ag0;
                      gcsbsb = -here->BSIM4v6capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
              }
          }
          else
          {   if (here->BSIM4v6trnqsMod == 0)
              {   qsrc = qdrn - qgso;
                  if (here->BSIM4v6rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v6cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                          gcgmsb = -cgso * ag0;
                          gcgmbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcggb = here->BSIM4v6cggb * ag0;
                      gcgdb = here->BSIM4v6cgsb * ag0;
                      gcgsb = here->BSIM4v6cgdb * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = -(here->BSIM4v6cggb + here->BSIM4v6cbgb
                            + here->BSIM4v6cdgb) * ag0;
                      gcsgb = here->BSIM4v6cdgb * ag0;
                      gcbgb = here->BSIM4v6cbgb * ag0;

                      qgmb = pParam->BSIM4v6cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qbulk -= qgmb;
                      qdrn = -(qgate + qgmid + qbulk + qsrc);
                  }
                  else
                  {   gcggb = (here->BSIM4v6cggb + cgdo + cgso
                            + pParam->BSIM4v6cgbo ) * ag0;
                      gcgdb = (here->BSIM4v6cgsb - cgdo) * ag0;
                      gcgsb = (here->BSIM4v6cgdb - cgso) * ag0;
                      gcgbb = -(gcggb + gcgdb + gcgsb);

                      gcdgb = -(here->BSIM4v6cggb + here->BSIM4v6cbgb
                            + here->BSIM4v6cdgb + cgdo) * ag0;
                      gcsgb = (here->BSIM4v6cdgb - cgso) * ag0;
                      gcbgb = (here->BSIM4v6cbgb - pParam->BSIM4v6cgbo) * ag0;

                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v6cgbo * vgb;
                      qgate += qgdo + qgso + qgb;
                      qbulk -= qgb;
                      qdrn = -(qgate + qbulk + qsrc);
                   }
                  gcddb = (here->BSIM4v6capbd + cgdo - (here->BSIM4v6cgsb
                        + here->BSIM4v6cbsb + here->BSIM4v6cdsb)) * ag0;
                  gcdsb = -(here->BSIM4v6cgdb + here->BSIM4v6cbdb
                        + here->BSIM4v6cddb) * ag0;

                  gcsdb = here->BSIM4v6cdsb * ag0;
                  gcssb = (here->BSIM4v6cddb + here->BSIM4v6capbs + cgso) * ag0;

                  if (!here->BSIM4v6rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                      gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                      gcbdb = (here->BSIM4v6cbsb - here->BSIM4v6capbd) * ag0;
                      gcbsb = (here->BSIM4v6cbdb - here->BSIM4v6capbs) * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb)
                            + here->BSIM4v6capbd * ag0;
                      gcsbb = -(here->BSIM4v6cddb + here->BSIM4v6cdgb
                            + here->BSIM4v6cdsb) * ag0;
                      gcbdb = here->BSIM4v6cbsb * ag0;
                      gcbsb = here->BSIM4v6cbdb * ag0;
                      gcdbdb = -here->BSIM4v6capbd * ag0;
                      gcsbsb = -here->BSIM4v6capbs * ag0;
                  }
                  gcbbb = -(gcbgb + gcbdb + gcbsb + gcbgmb);

                  ggtg = ggtd = ggtb = ggts = 0.0;
                  sxpart = 0.4;
                  dxpart = 0.6;
                  ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
                  dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;
              }
              else
              {   qcheq = here->BSIM4v6qchqs;
                  CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV * here->BSIM4v6nf
                        * pParam->BSIM4v6leffCV;
                  T0 = qdef * ScalingFactor / CoxWL;
                  ggtg = here->BSIM4v6gtg = T0 * here->BSIM4v6gcrgg;
                  ggts = here->BSIM4v6gts = T0 * here->BSIM4v6gcrgd;
                  ggtd = here->BSIM4v6gtd = T0 * here->BSIM4v6gcrgs;
                  ggtb = here->BSIM4v6gtb = T0 * here->BSIM4v6gcrgb;
                  gqdef = ScalingFactor * ag0;

                  gcqgb = here->BSIM4v6cqgb * ag0;
                  gcqdb = here->BSIM4v6cqsb * ag0;
                  gcqsb = here->BSIM4v6cqdb * ag0;
                  gcqbb = here->BSIM4v6cqbb * ag0;

                  if (fabs(qcheq) <= 1.0e-5 * CoxWL)
                  {   if (model->BSIM4v6xpart < 0.5)
                      {   sxpart = 0.4;
                      }
                      else if (model->BSIM4v6xpart > 0.5)
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
                      Css = here->BSIM4v6cddb;
                      Cds = -(here->BSIM4v6cgdb + here->BSIM4v6cddb
                          + here->BSIM4v6cbdb);
                      dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                      Csg = here->BSIM4v6cdgb;
                      Cdg = -(here->BSIM4v6cggb + here->BSIM4v6cdgb
                          + here->BSIM4v6cbgb);
                      dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                      Csd = here->BSIM4v6cdsb;
                      Cdd = -(here->BSIM4v6cgsb + here->BSIM4v6cdsb
                          + here->BSIM4v6cbsb);
                      dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                  }
                  dxpart = 1.0 - sxpart;
                  ddxpart_dVd = -dsxpart_dVd;
                  ddxpart_dVg = -dsxpart_dVg;
                  ddxpart_dVs = -dsxpart_dVs;
                  ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                  if (here->BSIM4v6rgateMod == 3)
                  {   gcgmgmb = (cgdo + cgso + pParam->BSIM4v6cgbo) * ag0;
                      gcgmdb = -cgdo * ag0;
                      gcgmsb = -cgso * ag0;
                      gcgmbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgmb = gcgmdb;
                      gcsgmb = gcgmsb;
                      gcbgmb = gcgmbb;

                      gcdgb = gcsgb = gcbgb = 0.0;
                      gcggb = gcgdb = gcgsb = gcgbb = 0.0;

                      qgmb = pParam->BSIM4v6cgbo * vgmb;
                      qgmid = qgdo + qgso + qgmb;
                      qgate = 0.0;
                      qbulk = -qgmb;
                      qdrn = -qgdo;
                      qsrc = -qgso;
                  }
                  else
                  {   gcggb = (cgdo + cgso + pParam->BSIM4v6cgbo ) * ag0;
                      gcgdb = -cgdo * ag0;
                      gcgsb = -cgso * ag0;
                      gcgbb = -pParam->BSIM4v6cgbo * ag0;

                      gcdgb = gcgdb;
                      gcsgb = gcgsb;
                      gcbgb = gcgbb;
                      gcdgmb = gcsgmb = gcbgmb = 0.0;

                      qgb = pParam->BSIM4v6cgbo * vgb;
                      qgate = qgdo + qgso + qgb;
                      qbulk = -qgb;
                      qdrn = -qgdo;
                      qsrc = -qgso;
                  }

                  gcddb = (here->BSIM4v6capbd + cgdo) * ag0;
                  gcdsb = gcsdb = 0.0;
                  gcssb = (here->BSIM4v6capbs + cgso) * ag0;
                  if (!here->BSIM4v6rbodyMod)
                  {   gcdbb = -(gcdgb + gcddb + gcdgmb);
                      gcsbb = -(gcsgb + gcssb + gcsgmb);
                      gcbdb = -here->BSIM4v6capbd * ag0;
                      gcbsb = -here->BSIM4v6capbs * ag0;
                      gcdbdb = 0.0; gcsbsb = 0.0;
                  }
                  else
                  {   gcdbb = gcsbb = gcbdb = gcbsb = 0.0;
                      gcdbdb = -here->BSIM4v6capbd * ag0;
                      gcsbsb = -here->BSIM4v6capbs * ag0;
                  }
                  gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
              }
          }


          if (here->BSIM4v6trnqsMod)
          {   *(ckt->CKTstate0 + here->BSIM4v6qcdump) = qdef * ScalingFactor;
              if (ckt->CKTmode & MODEINITTRAN)
                  *(ckt->CKTstate1 + here->BSIM4v6qcdump) =
                                   *(ckt->CKTstate0 + here->BSIM4v6qcdump);
              error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qcdump);
              if (error)
                  return(error);
          }

          if (ByPass) goto line860;

          *(ckt->CKTstate0 + here->BSIM4v6qg) = qgate;
          *(ckt->CKTstate0 + here->BSIM4v6qd) = qdrn
                           - *(ckt->CKTstate0 + here->BSIM4v6qbd);
          *(ckt->CKTstate0 + here->BSIM4v6qs) = qsrc
                           - *(ckt->CKTstate0 + here->BSIM4v6qbs);
          if (here->BSIM4v6rgateMod == 3)
              *(ckt->CKTstate0 + here->BSIM4v6qgmid) = qgmid;

          if (!here->BSIM4v6rbodyMod)
          {   *(ckt->CKTstate0 + here->BSIM4v6qb) = qbulk
                               + *(ckt->CKTstate0 + here->BSIM4v6qbd)
                               + *(ckt->CKTstate0 + here->BSIM4v6qbs);
          }
          else
              *(ckt->CKTstate0 + here->BSIM4v6qb) = qbulk;


          /* Store small signal parameters */
          if (ckt->CKTmode & MODEINITSMSIG)
          {   goto line1000;
          }

          if (!ChargeComputationNeeded)
              goto line850;

          if (ckt->CKTmode & MODEINITTRAN)
          {   *(ckt->CKTstate1 + here->BSIM4v6qb) =
                    *(ckt->CKTstate0 + here->BSIM4v6qb);
              *(ckt->CKTstate1 + here->BSIM4v6qg) =
                    *(ckt->CKTstate0 + here->BSIM4v6qg);
              *(ckt->CKTstate1 + here->BSIM4v6qd) =
                    *(ckt->CKTstate0 + here->BSIM4v6qd);
              if (here->BSIM4v6rgateMod == 3)
                  *(ckt->CKTstate1 + here->BSIM4v6qgmid) =
                        *(ckt->CKTstate0 + here->BSIM4v6qgmid);
              if (here->BSIM4v6rbodyMod)
              {   *(ckt->CKTstate1 + here->BSIM4v6qbs) =
                                   *(ckt->CKTstate0 + here->BSIM4v6qbs);
                  *(ckt->CKTstate1 + here->BSIM4v6qbd) =
                                   *(ckt->CKTstate0 + here->BSIM4v6qbd);
              }
          }

          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qb);
          if (error) 
              return(error);
          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qg);
          if (error) 
              return(error);
          error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qd);
          if (error) 
              return(error);

          if (here->BSIM4v6rgateMod == 3)
          {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qgmid);
              if (error) return(error);
          }

          if (here->BSIM4v6rbodyMod)
          {   error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qbs);
              if (error) 
                  return(error);
              error = NIintegrate(ckt, &geq, &ceq, 0.0, here->BSIM4v6qbd);
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
          sxpart = (1.0 - (dxpart = (here->BSIM4v6mode > 0) ? 0.4 : 0.6));
          ddxpart_dVd = ddxpart_dVg = ddxpart_dVb = ddxpart_dVs = 0.0;
          dsxpart_dVd = dsxpart_dVg = dsxpart_dVb = dsxpart_dVs = 0.0;

          if (here->BSIM4v6trnqsMod)
          {   CoxWL = model->BSIM4v6coxe * pParam->BSIM4v6weffCV * here->BSIM4v6nf
                    * pParam->BSIM4v6leffCV;
              T1 = here->BSIM4v6gcrg / CoxWL;
              here->BSIM4v6gtau = T1 * ScalingFactor;
          }
          else
              here->BSIM4v6gtau = 0.0;

          goto line900;

            
line860:
          /* Calculate equivalent charge current */

          cqgate = *(ckt->CKTstate0 + here->BSIM4v6cqg);
          cqbody = *(ckt->CKTstate0 + here->BSIM4v6cqb);
          cqdrn = *(ckt->CKTstate0 + here->BSIM4v6cqd);

          ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
          ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb)
                * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
          ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb
                + gcbdb * vbd + gcbsb * vbs;


          if (here->BSIM4v6rgateMod == 3)
              ceqqgmid = *(ckt->CKTstate0 + here->BSIM4v6cqgmid)
                       + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;
          else
               ceqqgmid = 0.0;        

          if (here->BSIM4v6rbodyMod)
          {   ceqqjs = *(ckt->CKTstate0 + here->BSIM4v6cqbs) + gcsbsb * vbs_jct;
              ceqqjd = *(ckt->CKTstate0 + here->BSIM4v6cqbd) + gcdbdb * vbd_jct; 
          }

          if (here->BSIM4v6trnqsMod)
          {   T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
              ceqqg += T0;
              T1 = qdef * here->BSIM4v6gtau;
              ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd
                     * vbd - ddxpart_dVs * vbs);
              cqdef = *(ckt->CKTstate0 + here->BSIM4v6cqcdump) - gqdef * qdef;
              cqcheq = *(ckt->CKTstate0 + here->BSIM4v6cqcheq)
                     - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
          }

          if (ckt->CKTmode & MODEINITTRAN)
          {   *(ckt->CKTstate1 + here->BSIM4v6cqb) =
                               *(ckt->CKTstate0 + here->BSIM4v6cqb);
              *(ckt->CKTstate1 + here->BSIM4v6cqg) =
                               *(ckt->CKTstate0 + here->BSIM4v6cqg);
              *(ckt->CKTstate1 + here->BSIM4v6cqd) =
                               *(ckt->CKTstate0 + here->BSIM4v6cqd);

              if (here->BSIM4v6rgateMod == 3)
                  *(ckt->CKTstate1 + here->BSIM4v6cqgmid) =
                                   *(ckt->CKTstate0 + here->BSIM4v6cqgmid);

              if (here->BSIM4v6rbodyMod)
              {   *(ckt->CKTstate1 + here->BSIM4v6cqbs) =
                                   *(ckt->CKTstate0 + here->BSIM4v6cqbs);
                  *(ckt->CKTstate1 + here->BSIM4v6cqbd) =
                                   *(ckt->CKTstate0 + here->BSIM4v6cqbd);
              }
          }


          /*
           *  Load current vector
           */

line900:
          if (here->BSIM4v6mode >= 0)
          {   Gm = here->BSIM4v6gm;
              Gmbs = here->BSIM4v6gmbs;
              FwdSum = Gm + Gmbs;
              RevSum = 0.0;

              ceqdrn = model->BSIM4v6type * (cdrain - here->BSIM4v6gds * vds
                     - Gm * vgs - Gmbs * vbs);
              ceqbd = model->BSIM4v6type * (here->BSIM4v6csub + here->BSIM4v6Igidl
                    - (here->BSIM4v6gbds + here->BSIM4v6ggidld) * vds
                    - (here->BSIM4v6gbgs + here->BSIM4v6ggidlg) * vgs
                    - (here->BSIM4v6gbbs + here->BSIM4v6ggidlb) * vbs);
              ceqbs = model->BSIM4v6type * (here->BSIM4v6Igisl + here->BSIM4v6ggisls * vds 
                          - here->BSIM4v6ggislg * vgd - here->BSIM4v6ggislb * vbd);

              gbbdp = -(here->BSIM4v6gbds);
              gbbsp = here->BSIM4v6gbds + here->BSIM4v6gbgs + here->BSIM4v6gbbs; 
                    
              gbdpg = here->BSIM4v6gbgs;
              gbdpdp = here->BSIM4v6gbds;
              gbdpb = here->BSIM4v6gbbs;
              gbdpsp = -(gbdpg + gbdpdp + gbdpb);

              gbspg = 0.0;
              gbspdp = 0.0;
              gbspb = 0.0;
              gbspsp = 0.0;

              if (model->BSIM4v6igcMod)
              {   gIstotg = here->BSIM4v6gIgsg + here->BSIM4v6gIgcsg;
                  gIstotd = here->BSIM4v6gIgcsd;
                  gIstots = here->BSIM4v6gIgss + here->BSIM4v6gIgcss;
                  gIstotb = here->BSIM4v6gIgcsb;
                   Istoteq = model->BSIM4v6type * (here->BSIM4v6Igs + here->BSIM4v6Igcs
                            - gIstotg * vgs - here->BSIM4v6gIgcsd * vds
                          - here->BSIM4v6gIgcsb * vbs);

                  gIdtotg = here->BSIM4v6gIgdg + here->BSIM4v6gIgcdg;
                  gIdtotd = here->BSIM4v6gIgdd + here->BSIM4v6gIgcdd;
                  gIdtots = here->BSIM4v6gIgcds;
                  gIdtotb = here->BSIM4v6gIgcdb;
                  Idtoteq = model->BSIM4v6type * (here->BSIM4v6Igd + here->BSIM4v6Igcd
                          - here->BSIM4v6gIgdg * vgd - here->BSIM4v6gIgcdg * vgs
                          - here->BSIM4v6gIgcdd * vds - here->BSIM4v6gIgcdb * vbs);
              }
              else
              {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                  gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
              }

              if (model->BSIM4v6igbMod)
              {   gIbtotg = here->BSIM4v6gIgbg;
                  gIbtotd = here->BSIM4v6gIgbd;
                  gIbtots = here->BSIM4v6gIgbs;
                  gIbtotb = here->BSIM4v6gIgbb;
                  Ibtoteq = model->BSIM4v6type * (here->BSIM4v6Igb
                          - here->BSIM4v6gIgbg * vgs - here->BSIM4v6gIgbd * vds
                          - here->BSIM4v6gIgbb * vbs);
              }
              else
                  gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

              if ((model->BSIM4v6igcMod != 0) || (model->BSIM4v6igbMod != 0))
              {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                  gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                  gIgtots = gIstots + gIdtots + gIbtots;
                  gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  Igtoteq = Istoteq + Idtoteq + Ibtoteq; 
              }
              else
                  gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


              if (here->BSIM4v6rgateMod == 2)
                  T0 = vges - vgs;
              else if (here->BSIM4v6rgateMod == 3)
                  T0 = vgms - vgs;
              if (here->BSIM4v6rgateMod > 1)
              {   gcrgd = here->BSIM4v6gcrgd * T0;
                  gcrgg = here->BSIM4v6gcrgg * T0;
                  gcrgs = here->BSIM4v6gcrgs * T0;
                  gcrgb = here->BSIM4v6gcrgb * T0;
                  ceqgcrg = -(gcrgd * vds + gcrgg * vgs
                          + gcrgb * vbs);
                  gcrgg -= here->BSIM4v6gcrg;
                  gcrg = here->BSIM4v6gcrg;
              }
              else
                  ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
          }
          else
          {   Gm = -here->BSIM4v6gm;
              Gmbs = -here->BSIM4v6gmbs;
              FwdSum = 0.0;
              RevSum = -(Gm + Gmbs);

              ceqdrn = -model->BSIM4v6type * (cdrain + here->BSIM4v6gds * vds
                     + Gm * vgd + Gmbs * vbd);

              ceqbs = model->BSIM4v6type * (here->BSIM4v6csub + here->BSIM4v6Igisl 
                    + (here->BSIM4v6gbds + here->BSIM4v6ggisls) * vds
                    - (here->BSIM4v6gbgs + here->BSIM4v6ggislg) * vgd
                    - (here->BSIM4v6gbbs + here->BSIM4v6ggislb) * vbd);
              ceqbd = model->BSIM4v6type * (here->BSIM4v6Igidl - here->BSIM4v6ggidld * vds 
                              - here->BSIM4v6ggidlg * vgs - here->BSIM4v6ggidlb * vbs);

              gbbsp = -(here->BSIM4v6gbds);
              gbbdp = here->BSIM4v6gbds + here->BSIM4v6gbgs + here->BSIM4v6gbbs; 

              gbdpg = 0.0;
              gbdpsp = 0.0;
              gbdpb = 0.0;
              gbdpdp = 0.0;

              gbspg = here->BSIM4v6gbgs;
              gbspsp = here->BSIM4v6gbds;
              gbspb = here->BSIM4v6gbbs;
              gbspdp = -(gbspg + gbspsp + gbspb);

              if (model->BSIM4v6igcMod)
              {   gIstotg = here->BSIM4v6gIgsg + here->BSIM4v6gIgcdg;
                  gIstotd = here->BSIM4v6gIgcds;
                  gIstots = here->BSIM4v6gIgss + here->BSIM4v6gIgcdd;
                  gIstotb = here->BSIM4v6gIgcdb;
                  Istoteq = model->BSIM4v6type * (here->BSIM4v6Igs + here->BSIM4v6Igcd
                          - here->BSIM4v6gIgsg * vgs - here->BSIM4v6gIgcdg * vgd
                          + here->BSIM4v6gIgcdd * vds - here->BSIM4v6gIgcdb * vbd);

                  gIdtotg = here->BSIM4v6gIgdg + here->BSIM4v6gIgcsg;
                  gIdtotd = here->BSIM4v6gIgdd + here->BSIM4v6gIgcss;
                  gIdtots = here->BSIM4v6gIgcsd;
                  gIdtotb = here->BSIM4v6gIgcsb;
                  Idtoteq = model->BSIM4v6type * (here->BSIM4v6Igd + here->BSIM4v6Igcs
                          - (here->BSIM4v6gIgdg + here->BSIM4v6gIgcsg) * vgd
                          + here->BSIM4v6gIgcsd * vds - here->BSIM4v6gIgcsb * vbd);
              }
              else
              {   gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
                  gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
              }

              if (model->BSIM4v6igbMod)
              {   gIbtotg = here->BSIM4v6gIgbg;
                  gIbtotd = here->BSIM4v6gIgbs;
                  gIbtots = here->BSIM4v6gIgbd;
                  gIbtotb = here->BSIM4v6gIgbb;
                  Ibtoteq = model->BSIM4v6type * (here->BSIM4v6Igb
                          - here->BSIM4v6gIgbg * vgd + here->BSIM4v6gIgbd * vds
                          - here->BSIM4v6gIgbb * vbd);
              }
              else
                  gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;

              if ((model->BSIM4v6igcMod != 0) || (model->BSIM4v6igbMod != 0))
              {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                  gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                  gIgtots = gIstots + gIdtots + gIbtots;
                  gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  Igtoteq = Istoteq + Idtoteq + Ibtoteq;
              }
              else
                  gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;


              if (here->BSIM4v6rgateMod == 2)
                  T0 = vges - vgs;
              else if (here->BSIM4v6rgateMod == 3)
                  T0 = vgms - vgs;
              if (here->BSIM4v6rgateMod > 1)
              {   gcrgd = here->BSIM4v6gcrgs * T0;
                  gcrgg = here->BSIM4v6gcrgg * T0;
                  gcrgs = here->BSIM4v6gcrgd * T0;
                  gcrgb = here->BSIM4v6gcrgb * T0;
                  ceqgcrg = -(gcrgg * vgd - gcrgs * vds
                          + gcrgb * vbd);
                  gcrgg -= here->BSIM4v6gcrg;
                  gcrg = here->BSIM4v6gcrg;
              }
              else
                  ceqgcrg = gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;
          }

          if (model->BSIM4v6rdsMod == 1)
          {   ceqgstot = model->BSIM4v6type * (here->BSIM4v6gstotd * vds
                       + here->BSIM4v6gstotg * vgs + here->BSIM4v6gstotb * vbs);
              /* WDLiu: ceqgstot flowing away from sNodePrime */
              gstot = here->BSIM4v6gstot;
              gstotd = here->BSIM4v6gstotd;
              gstotg = here->BSIM4v6gstotg;
              gstots = here->BSIM4v6gstots - gstot;
              gstotb = here->BSIM4v6gstotb;

              ceqgdtot = -model->BSIM4v6type * (here->BSIM4v6gdtotd * vds
                       + here->BSIM4v6gdtotg * vgs + here->BSIM4v6gdtotb * vbs);
              /* WDLiu: ceqgdtot defined as flowing into dNodePrime */
              gdtot = here->BSIM4v6gdtot;
              gdtotd = here->BSIM4v6gdtotd - gdtot;
              gdtotg = here->BSIM4v6gdtotg;
              gdtots = here->BSIM4v6gdtots;
              gdtotb = here->BSIM4v6gdtotb;
          }
          else
          {   gstot = gstotd = gstotg = gstots = gstotb = ceqgstot = 0.0;
              gdtot = gdtotd = gdtotg = gdtots = gdtotb = ceqgdtot = 0.0;
          }

           if (model->BSIM4v6type > 0)
           {   ceqjs = (here->BSIM4v6cbs - here->BSIM4v6gbs * vbs_jct);
               ceqjd = (here->BSIM4v6cbd - here->BSIM4v6gbd * vbd_jct);
           }
           else
           {   ceqjs = -(here->BSIM4v6cbs - here->BSIM4v6gbs * vbs_jct); 
               ceqjd = -(here->BSIM4v6cbd - here->BSIM4v6gbd * vbd_jct);
               ceqqg = -ceqqg;
               ceqqd = -ceqqd;
               ceqqb = -ceqqb;
               ceqgcrg = -ceqgcrg;

               if (here->BSIM4v6trnqsMod)
               {   cqdef = -cqdef;
                   cqcheq = -cqcheq;
               }

               if (here->BSIM4v6rbodyMod)
               {   ceqqjs = -ceqqjs;
                   ceqqjd = -ceqqjd;
               }

               if (here->BSIM4v6rgateMod == 3)
                   ceqqgmid = -ceqqgmid;
           }


           /*
            *  Loading RHS
            */

              m = here->BSIM4v6m;

#ifdef USE_OMP
       here->BSIM4v6rhsdPrime = m * (ceqjd - ceqbd + ceqgdtot
                                                    - ceqdrn - ceqqd + Idtoteq);
       here->BSIM4v6rhsgPrime = m * (ceqqg - ceqgcrg + Igtoteq);

       if (here->BSIM4v6rgateMod == 2)
           here->BSIM4v6rhsgExt = m * ceqgcrg;
       else if (here->BSIM4v6rgateMod == 3)
               here->BSIM4v6grhsMid = m * (ceqqgmid + ceqgcrg);

       if (!here->BSIM4v6rbodyMod)
       {   here->BSIM4v6rhsbPrime = m * (ceqbd + ceqbs - ceqjd
                                                        - ceqjs - ceqqb + Ibtoteq);
           here->BSIM4v6rhssPrime = m * (ceqdrn - ceqbs + ceqjs 
                              + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq);
        }
        else
        {   here->BSIM4v6rhsdb = m * (ceqjd + ceqqjd);
            here->BSIM4v6rhsbPrime = m * (ceqbd + ceqbs - ceqqb + Ibtoteq);
            here->BSIM4v6rhssb = m * (ceqjs + ceqqjs);
            here->BSIM4v6rhssPrime = m * (ceqdrn - ceqbs + ceqjs + ceqqd 
                + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq);
        }

        if (model->BSIM4v6rdsMod)
        {   here->BSIM4v6rhsd = m * ceqgdtot; 
            here->BSIM4v6rhss = m * ceqgstot;
        }

        if (here->BSIM4v6trnqsMod)
           here->BSIM4v6rhsq = m * (cqcheq - cqdef);
#else
        (*(ckt->CKTrhs + here->BSIM4v6dNodePrime) += m * (ceqjd - ceqbd + ceqgdtot
                                                    - ceqdrn - ceqqd + Idtoteq));
        (*(ckt->CKTrhs + here->BSIM4v6gNodePrime) -= m * (ceqqg - ceqgcrg + Igtoteq));

        if (here->BSIM4v6rgateMod == 2)
            (*(ckt->CKTrhs + here->BSIM4v6gNodeExt) -= m * ceqgcrg);
        else if (here->BSIM4v6rgateMod == 3)
            (*(ckt->CKTrhs + here->BSIM4v6gNodeMid) -= m * (ceqqgmid + ceqgcrg));

        if (!here->BSIM4v6rbodyMod)
        {   (*(ckt->CKTrhs + here->BSIM4v6bNodePrime) += m * (ceqbd + ceqbs - ceqjd
                                                        - ceqjs - ceqqb + Ibtoteq));
            (*(ckt->CKTrhs + here->BSIM4v6sNodePrime) += m * (ceqdrn - ceqbs + ceqjs 
                              + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq));
        }

        else
        {   (*(ckt->CKTrhs + here->BSIM4v6dbNode) -= m * (ceqjd + ceqqjd));
            (*(ckt->CKTrhs + here->BSIM4v6bNodePrime) += m * (ceqbd + ceqbs - ceqqb + Ibtoteq));
            (*(ckt->CKTrhs + here->BSIM4v6sbNode) -= m * (ceqjs + ceqqjs));
            (*(ckt->CKTrhs + here->BSIM4v6sNodePrime) += m * (ceqdrn - ceqbs + ceqjs + ceqqd 
                + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq));
        }

        if (model->BSIM4v6rdsMod)
        {   (*(ckt->CKTrhs + here->BSIM4v6dNode) -= m * ceqgdtot); 
            (*(ckt->CKTrhs + here->BSIM4v6sNode) += m * ceqgstot);
        }

        if (here->BSIM4v6trnqsMod)
            *(ckt->CKTrhs + here->BSIM4v6qNode) += m * (cqcheq - cqdef);
#endif

        /*
         *  Loading matrix
         */

           if (!here->BSIM4v6rbodyMod)
           {   gjbd = here->BSIM4v6gbd;
               gjbs = here->BSIM4v6gbs;
           }
       else
          gjbd = gjbs = 0.0;

       if (!model->BSIM4v6rdsMod)
       {   gdpr = here->BSIM4v6drainConductance;
           gspr = here->BSIM4v6sourceConductance;
       }
       else
           gdpr = gspr = 0.0;

       geltd = here->BSIM4v6grgeltd;

       T1 = qdef * here->BSIM4v6gtau;
#ifdef USE_OMP
       if (here->BSIM4v6rgateMod == 1)
       {   here->BSIM4v6_1 = m * geltd;
           here->BSIM4v6_2 = m * geltd;
           here->BSIM4v6_3 = m * geltd;
           here->BSIM4v6_4 = m * (gcggb + geltd - ggtg + gIgtotg);
           here->BSIM4v6_5 = m * (gcgdb - ggtd + gIgtotd);
           here->BSIM4v6_6 = m * (gcgsb - ggts + gIgtots);
           here->BSIM4v6_7 = m * (gcgbb - ggtb + gIgtotb);
       } /* WDLiu: gcrg already subtracted from all gcrgg below */
       else if (here->BSIM4v6rgateMod == 2)        
       {   here->BSIM4v6_8 = m * gcrg;
           here->BSIM4v6_9 = m * gcrgg;
           here->BSIM4v6_10 = m * gcrgd;
           here->BSIM4v6_11 = m * gcrgs;
           here->BSIM4v6_12 = m * gcrgb;        

           here->BSIM4v6_13 = m * gcrg;
           here->BSIM4v6_14 = m * (gcggb  - gcrgg - ggtg + gIgtotg);
           here->BSIM4v6_15 = m * (gcgdb - gcrgd - ggtd + gIgtotd);
           here->BSIM4v6_16 = m * (gcgsb - gcrgs - ggts + gIgtots);
           here->BSIM4v6_17 = m * (gcgbb - gcrgb - ggtb + gIgtotb);
       }
       else if (here->BSIM4v6rgateMod == 3)
       {   here->BSIM4v6_18 = m * geltd;
           here->BSIM4v6_19 = m * geltd;
           here->BSIM4v6_20 = m * geltd;
           here->BSIM4v6_21 = m * (geltd + gcrg + gcgmgmb);

           here->BSIM4v6_22 = m * (gcrgd + gcgmdb);
           here->BSIM4v6_23 = m * gcrgg;
           here->BSIM4v6_24 = m * (gcrgs + gcgmsb);
           here->BSIM4v6_25 = m * (gcrgb + gcgmbb);

           here->BSIM4v6_26 = m * gcdgmb;
           here->BSIM4v6_26 = m * gcrg;
           here->BSIM4v6_28 = m * gcsgmb;
           here->BSIM4v6_29 = m * gcbgmb;

           here->BSIM4v6_30 = m * (gcggb - gcrgg - ggtg + gIgtotg);
           here->BSIM4v6_31 = m * (gcgdb - gcrgd - ggtd + gIgtotd);
           here->BSIM4v6_32 = m * (gcgsb - gcrgs - ggts + gIgtots);
           here->BSIM4v6_33 = m * (gcgbb - gcrgb - ggtb + gIgtotb);
       }
       else
       {   here->BSIM4v6_34 = m * (gcggb - ggtg + gIgtotg);
           here->BSIM4v6_35 = m * (gcgdb - ggtd + gIgtotd);
           here->BSIM4v6_36 = m * (gcgsb - ggts + gIgtots);
           here->BSIM4v6_37 = m * (gcgbb - ggtb + gIgtotb);
       }

       if (model->BSIM4v6rdsMod)
       {   here->BSIM4v6_38 = m * gdtotg;
           here->BSIM4v6_39 = m * gdtots;
           here->BSIM4v6_40 = m * gdtotb;
           here->BSIM4v6_41 = m * gstotd;
           here->BSIM4v6_42 = m * gstotg;
           here->BSIM4v6_43 = m * gstotb;
       }

       here->BSIM4v6_44 = m * (gdpr + here->BSIM4v6gds + here->BSIM4v6gbd + T1 * ddxpart_dVd
                                   - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd);
       here->BSIM4v6_45 = m * (gdpr + gdtot);
       here->BSIM4v6_46 = m * (Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                                   + dxpart * ggtg + T1 * ddxpart_dVg);
       here->BSIM4v6_47 = m * (here->BSIM4v6gds + gdtots - dxpart * ggts + gIdtots
                                   - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp);
       here->BSIM4v6_48 = m * (gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                                   - T1 * ddxpart_dVb - dxpart * ggtb);

       here->BSIM4v6_49 = m * (gdpr - gdtotd);
       here->BSIM4v6_50 = m * (gdpr + gdtot);

       here->BSIM4v6_51 = m * (here->BSIM4v6gds + gstotd + RevSum - gcsdb - gbspdp
                                   - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd);
       here->BSIM4v6_52 = m * (gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                                   + T1 * dsxpart_dVg - gIstotg);
       here->BSIM4v6_53 = m * (gspr + here->BSIM4v6gds + here->BSIM4v6gbs + T1 * dsxpart_dVs
                                   - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots);
       here->BSIM4v6_54 = m * (gspr + gstot);
       here->BSIM4v6_55 = m * (gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                                   - T1 * dsxpart_dVb + gIstotb);

       here->BSIM4v6_56 = m * (gspr - gstots);
       here->BSIM4v6_57 = m * (gspr + gstot);

       here->BSIM4v6_58 = m * (gcbdb - gjbd + gbbdp - gIbtotd);
       here->BSIM4v6_59 = m * (gcbgb - here->BSIM4v6gbgs - gIbtotg);
       here->BSIM4v6_60 = m * (gcbsb - gjbs + gbbsp - gIbtots);
       here->BSIM4v6_61 = m * (gjbd + gjbs + gcbbb - here->BSIM4v6gbbs - gIbtotb);

       ggidld = here->BSIM4v6ggidld;
       ggidlg = here->BSIM4v6ggidlg;
       ggidlb = here->BSIM4v6ggidlb;
       ggislg = here->BSIM4v6ggislg;
       ggisls = here->BSIM4v6ggisls;
       ggislb = here->BSIM4v6ggislb;

       /* stamp gidl */
       here->BSIM4v6_62 = m * ggidld;
       here->BSIM4v6_63 = m * ggidlg;
       here->BSIM4v6_64 = m * (ggidlg + ggidld + ggidlb);
       here->BSIM4v6_65 = m * ggidlb;
       here->BSIM4v6_66 = m * ggidld;
       here->BSIM4v6_67 = m * ggidlg;
       here->BSIM4v6_68 = m * (ggidlg + ggidld + ggidlb);
       here->BSIM4v6_69 = m * ggidlb;
       /* stamp gisl */
       here->BSIM4v6_70 = m * (ggisls + ggislg + ggislb);
       here->BSIM4v6_71 = m * ggislg;
       here->BSIM4v6_72 = m * ggisls;
       here->BSIM4v6_73 = m * ggislb;
       here->BSIM4v6_74 = m * (ggislg + ggisls + ggislb);
       here->BSIM4v6_75 = m * ggislg;
       here->BSIM4v6_76 = m * ggisls;
       here->BSIM4v6_77 = m * ggislb;

       if (here->BSIM4v6rbodyMod)
       {   here->BSIM4v6_78 = m * (gcdbdb - here->BSIM4v6gbd);
           here->BSIM4v6_79 = m * (here->BSIM4v6gbs - gcsbsb);

           here->BSIM4v6_80 = m * (gcdbdb - here->BSIM4v6gbd);
           here->BSIM4v6_81 = m * (here->BSIM4v6gbd - gcdbdb 
                          + here->BSIM4v6grbpd + here->BSIM4v6grbdb);
           here->BSIM4v6_82 = m * here->BSIM4v6grbpd;
           here->BSIM4v6_83 = m * here->BSIM4v6grbdb;

           here->BSIM4v6_84 = m * here->BSIM4v6grbpd;
           here->BSIM4v6_85 = m * here->BSIM4v6grbpb;
           here->BSIM4v6_86 = m * here->BSIM4v6grbps;
           here->BSIM4v6_87 = m * (here->BSIM4v6grbpd + here->BSIM4v6grbps 
                          + here->BSIM4v6grbpb);
           /* WDLiu: (gcbbb - here->BSIM4v6gbbs) already added to BPbpPtr */        

           here->BSIM4v6_88 = m * (gcsbsb - here->BSIM4v6gbs);
           here->BSIM4v6_89 = m * here->BSIM4v6grbps;
           here->BSIM4v6_90 = m * here->BSIM4v6grbsb;
           here->BSIM4v6_91 = m * (here->BSIM4v6gbs - gcsbsb 
                          + here->BSIM4v6grbps + here->BSIM4v6grbsb);

           here->BSIM4v6_92 = m * here->BSIM4v6grbdb;
           here->BSIM4v6_93 = m * here->BSIM4v6grbpb;
           here->BSIM4v6_94 = m * here->BSIM4v6grbsb;
           here->BSIM4v6_95 = m * (here->BSIM4v6grbsb + here->BSIM4v6grbdb
                           + here->BSIM4v6grbpb);
       }

           if (here->BSIM4v6trnqsMod)
           {   here->BSIM4v6_96 = m * (gqdef + here->BSIM4v6gtau);
               here->BSIM4v6_97 = m * (ggtg - gcqgb);
               here->BSIM4v6_98 = m * (ggtd - gcqdb);
               here->BSIM4v6_99 = m * (ggts - gcqsb);
               here->BSIM4v6_100 = m * (ggtb - gcqbb);

               here->BSIM4v6_101 = m * dxpart * here->BSIM4v6gtau;
               here->BSIM4v6_102 = m * sxpart * here->BSIM4v6gtau;
               here->BSIM4v6_103 = m * here->BSIM4v6gtau;
           }
#else
           if (here->BSIM4v6rgateMod == 1)
           {   (*(here->BSIM4v6GEgePtr) += m * geltd);
               (*(here->BSIM4v6GPgePtr) -= m * geltd);
               (*(here->BSIM4v6GEgpPtr) -= m * geltd);
               (*(here->BSIM4v6GPgpPtr) += m * (gcggb + geltd - ggtg + gIgtotg));
               (*(here->BSIM4v6GPdpPtr) += m * (gcgdb - ggtd + gIgtotd));
               (*(here->BSIM4v6GPspPtr) += m * (gcgsb - ggts + gIgtots));
               (*(here->BSIM4v6GPbpPtr) += m * (gcgbb - ggtb + gIgtotb));
           } /* WDLiu: gcrg already subtracted from all gcrgg below */
           else if (here->BSIM4v6rgateMod == 2)        
           {   (*(here->BSIM4v6GEgePtr) += m * gcrg);
               (*(here->BSIM4v6GEgpPtr) += m * gcrgg);
               (*(here->BSIM4v6GEdpPtr) += m * gcrgd);
               (*(here->BSIM4v6GEspPtr) += m * gcrgs);
               (*(here->BSIM4v6GEbpPtr) += m * gcrgb);        

               (*(here->BSIM4v6GPgePtr) -= m * gcrg);
               (*(here->BSIM4v6GPgpPtr) += m * (gcggb  - gcrgg - ggtg + gIgtotg));
               (*(here->BSIM4v6GPdpPtr) += m * (gcgdb - gcrgd - ggtd + gIgtotd));
               (*(here->BSIM4v6GPspPtr) += m * (gcgsb - gcrgs - ggts + gIgtots));
               (*(here->BSIM4v6GPbpPtr) += m * (gcgbb - gcrgb - ggtb + gIgtotb));
           }
           else if (here->BSIM4v6rgateMod == 3)
           {   (*(here->BSIM4v6GEgePtr) += m * geltd);
               (*(here->BSIM4v6GEgmPtr) -= m * geltd);
               (*(here->BSIM4v6GMgePtr) -= m * geltd);
               (*(here->BSIM4v6GMgmPtr) += m * (geltd + gcrg + gcgmgmb));

               (*(here->BSIM4v6GMdpPtr) += m * (gcrgd + gcgmdb));
               (*(here->BSIM4v6GMgpPtr) += m * gcrgg);
               (*(here->BSIM4v6GMspPtr) += m * (gcrgs + gcgmsb));
               (*(here->BSIM4v6GMbpPtr) += m * (gcrgb + gcgmbb));

               (*(here->BSIM4v6DPgmPtr) += m * gcdgmb);
               (*(here->BSIM4v6GPgmPtr) -= m * gcrg);
               (*(here->BSIM4v6SPgmPtr) += m * gcsgmb);
               (*(here->BSIM4v6BPgmPtr) += m * gcbgmb);

               (*(here->BSIM4v6GPgpPtr) += m * (gcggb - gcrgg - ggtg + gIgtotg));
               (*(here->BSIM4v6GPdpPtr) += m * (gcgdb - gcrgd - ggtd + gIgtotd));
               (*(here->BSIM4v6GPspPtr) += m * (gcgsb - gcrgs - ggts + gIgtots));
               (*(here->BSIM4v6GPbpPtr) += m * (gcgbb - gcrgb - ggtb + gIgtotb));
           }
            else
           {   (*(here->BSIM4v6GPgpPtr) += m * (gcggb - ggtg + gIgtotg));
               (*(here->BSIM4v6GPdpPtr) += m * (gcgdb - ggtd + gIgtotd));
               (*(here->BSIM4v6GPspPtr) += m * (gcgsb - ggts + gIgtots));
               (*(here->BSIM4v6GPbpPtr) += m * (gcgbb - ggtb + gIgtotb));
           }

           if (model->BSIM4v6rdsMod)
           {   (*(here->BSIM4v6DgpPtr) += m * gdtotg);
               (*(here->BSIM4v6DspPtr) += m * gdtots);
               (*(here->BSIM4v6DbpPtr) += m * gdtotb);
               (*(here->BSIM4v6SdpPtr) += m * gstotd);
               (*(here->BSIM4v6SgpPtr) += m * gstotg);
               (*(here->BSIM4v6SbpPtr) += m * gstotb);
           }

           (*(here->BSIM4v6DPdpPtr) += m * (gdpr + here->BSIM4v6gds + here->BSIM4v6gbd + T1 * ddxpart_dVd
                                   - gdtotd + RevSum + gcddb + gbdpdp + dxpart * ggtd - gIdtotd));
           (*(here->BSIM4v6DPdPtr) -= m * (gdpr + gdtot));
           (*(here->BSIM4v6DPgpPtr) += m * (Gm + gcdgb - gdtotg + gbdpg - gIdtotg
                                   + dxpart * ggtg + T1 * ddxpart_dVg));
           (*(here->BSIM4v6DPspPtr) -= m * (here->BSIM4v6gds + gdtots - dxpart * ggts + gIdtots
                                   - T1 * ddxpart_dVs + FwdSum - gcdsb - gbdpsp));
           (*(here->BSIM4v6DPbpPtr) -= m * (gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb
                                   - T1 * ddxpart_dVb - dxpart * ggtb));

           (*(here->BSIM4v6DdpPtr) -= m * (gdpr - gdtotd));
           (*(here->BSIM4v6DdPtr) += m * (gdpr + gdtot));

           (*(here->BSIM4v6SPdpPtr) -= m * (here->BSIM4v6gds + gstotd + RevSum - gcsdb - gbspdp
                                   - T1 * dsxpart_dVd - sxpart * ggtd + gIstotd));
           (*(here->BSIM4v6SPgpPtr) += m * (gcsgb - Gm - gstotg + gbspg + sxpart * ggtg
                                   + T1 * dsxpart_dVg - gIstotg));
           (*(here->BSIM4v6SPspPtr) += m * (gspr + here->BSIM4v6gds + here->BSIM4v6gbs + T1 * dsxpart_dVs
                                   - gstots + FwdSum + gcssb + gbspsp + sxpart * ggts - gIstots));
           (*(here->BSIM4v6SPsPtr) -= m * (gspr + gstot));
           (*(here->BSIM4v6SPbpPtr) -= m * (gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb
                                   - T1 * dsxpart_dVb + gIstotb));

           (*(here->BSIM4v6SspPtr) -= m * (gspr - gstots));
           (*(here->BSIM4v6SsPtr) += m * (gspr + gstot));

           (*(here->BSIM4v6BPdpPtr) += m * (gcbdb - gjbd + gbbdp - gIbtotd));
           (*(here->BSIM4v6BPgpPtr) += m * (gcbgb - here->BSIM4v6gbgs - gIbtotg));
           (*(here->BSIM4v6BPspPtr) += m * (gcbsb - gjbs + gbbsp - gIbtots));
           (*(here->BSIM4v6BPbpPtr) += m * (gjbd + gjbs + gcbbb - here->BSIM4v6gbbs
                                   - gIbtotb));

           ggidld = here->BSIM4v6ggidld;
           ggidlg = here->BSIM4v6ggidlg;
           ggidlb = here->BSIM4v6ggidlb;
           ggislg = here->BSIM4v6ggislg;
           ggisls = here->BSIM4v6ggisls;
           ggislb = here->BSIM4v6ggislb;

           /* stamp gidl */
           (*(here->BSIM4v6DPdpPtr) += m * ggidld);
           (*(here->BSIM4v6DPgpPtr) += m * ggidlg);
           (*(here->BSIM4v6DPspPtr) -= m * (ggidlg + ggidld + ggidlb));
           (*(here->BSIM4v6DPbpPtr) += m * ggidlb);
           (*(here->BSIM4v6BPdpPtr) -= m * ggidld);
           (*(here->BSIM4v6BPgpPtr) -= m * ggidlg);
           (*(here->BSIM4v6BPspPtr) += m * (ggidlg + ggidld + ggidlb));
           (*(here->BSIM4v6BPbpPtr) -= m * ggidlb);
            /* stamp gisl */
           (*(here->BSIM4v6SPdpPtr) -= m * (ggisls + ggislg + ggislb));
           (*(here->BSIM4v6SPgpPtr) += m * ggislg);
           (*(here->BSIM4v6SPspPtr) += m * ggisls);
           (*(here->BSIM4v6SPbpPtr) += m * ggislb);
           (*(here->BSIM4v6BPdpPtr) += m * (ggislg + ggisls + ggislb));
           (*(here->BSIM4v6BPgpPtr) -= m * ggislg);
           (*(here->BSIM4v6BPspPtr) -= m * ggisls);
           (*(here->BSIM4v6BPbpPtr) -= m * ggislb);


           if (here->BSIM4v6rbodyMod)
           {   (*(here->BSIM4v6DPdbPtr) += m * (gcdbdb - here->BSIM4v6gbd));
               (*(here->BSIM4v6SPsbPtr) -= m * (here->BSIM4v6gbs - gcsbsb));

               (*(here->BSIM4v6DBdpPtr) += m * (gcdbdb - here->BSIM4v6gbd));
               (*(here->BSIM4v6DBdbPtr) += m * (here->BSIM4v6gbd - gcdbdb 
                                       + here->BSIM4v6grbpd + here->BSIM4v6grbdb));
               (*(here->BSIM4v6DBbpPtr) -= m * here->BSIM4v6grbpd);
               (*(here->BSIM4v6DBbPtr) -= m * here->BSIM4v6grbdb);

               (*(here->BSIM4v6BPdbPtr) -= m * here->BSIM4v6grbpd);
               (*(here->BSIM4v6BPbPtr) -= m * here->BSIM4v6grbpb);
               (*(here->BSIM4v6BPsbPtr) -= m * here->BSIM4v6grbps);
               (*(here->BSIM4v6BPbpPtr) += m * (here->BSIM4v6grbpd + here->BSIM4v6grbps 
                                       + here->BSIM4v6grbpb));
               /* WDLiu: (gcbbb - here->BSIM4v6gbbs) already added to BPbpPtr */        

               (*(here->BSIM4v6SBspPtr) += m * (gcsbsb - here->BSIM4v6gbs));
               (*(here->BSIM4v6SBbpPtr) -= m * here->BSIM4v6grbps);
               (*(here->BSIM4v6SBbPtr) -= m * here->BSIM4v6grbsb);
               (*(here->BSIM4v6SBsbPtr) += m * (here->BSIM4v6gbs - gcsbsb 
                                       + here->BSIM4v6grbps + here->BSIM4v6grbsb));

               (*(here->BSIM4v6BdbPtr) -= m * here->BSIM4v6grbdb);
               (*(here->BSIM4v6BbpPtr) -= m * here->BSIM4v6grbpb);
               (*(here->BSIM4v6BsbPtr) -= m * here->BSIM4v6grbsb);
               (*(here->BSIM4v6BbPtr) += m * (here->BSIM4v6grbsb + here->BSIM4v6grbdb
                                     + here->BSIM4v6grbpb));
           }

           if (here->BSIM4v6trnqsMod)
           {   (*(here->BSIM4v6QqPtr) += m * (gqdef + here->BSIM4v6gtau));
               (*(here->BSIM4v6QgpPtr) += m * (ggtg - gcqgb));
               (*(here->BSIM4v6QdpPtr) += m * (ggtd - gcqdb));
               (*(here->BSIM4v6QspPtr) += m * (ggts - gcqsb));
               (*(here->BSIM4v6QbpPtr) += m * (ggtb - gcqbb));

               (*(here->BSIM4v6DPqPtr) += m * dxpart * here->BSIM4v6gtau);
               (*(here->BSIM4v6SPqPtr) += m * sxpart * here->BSIM4v6gtau);
               (*(here->BSIM4v6GPqPtr) -= m * here->BSIM4v6gtau);
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
int BSIM4v6polyDepletion(
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
void BSIM4v6LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt)
{
    unsigned int InstCount, idx;
    BSIM4v6instance **InstArray;
    BSIM4v6instance *here;
    BSIM4v6model *model = (BSIM4v6model*)inModel;

    InstArray = model->BSIM4v6InstanceArray;
    InstCount = model->BSIM4v6InstCount;

    for(idx = 0; idx < InstCount; idx++) {
       here = InstArray[idx];
        /* Update b for Ax = b */
           (*(ckt->CKTrhs + here->BSIM4v6dNodePrime) += here->BSIM4v6rhsdPrime);
           (*(ckt->CKTrhs + here->BSIM4v6gNodePrime) -= here->BSIM4v6rhsgPrime);

           if (here->BSIM4v6rgateMod == 2)
               (*(ckt->CKTrhs + here->BSIM4v6gNodeExt) -= here->BSIM4v6rhsgExt);
           else if (here->BSIM4v6rgateMod == 3)
               (*(ckt->CKTrhs + here->BSIM4v6gNodeMid) -= here->BSIM4v6grhsMid);

           if (!here->BSIM4v6rbodyMod)
           {   (*(ckt->CKTrhs + here->BSIM4v6bNodePrime) += here->BSIM4v6rhsbPrime);
               (*(ckt->CKTrhs + here->BSIM4v6sNodePrime) += here->BSIM4v6rhssPrime);
           }
           else
           {   (*(ckt->CKTrhs + here->BSIM4v6dbNode) -= here->BSIM4v6rhsdb);
               (*(ckt->CKTrhs + here->BSIM4v6bNodePrime) += here->BSIM4v6rhsbPrime);
               (*(ckt->CKTrhs + here->BSIM4v6sbNode) -= here->BSIM4v6rhssb);
               (*(ckt->CKTrhs + here->BSIM4v6sNodePrime) += here->BSIM4v6rhssPrime);
           }

           if (model->BSIM4v6rdsMod)
           {   (*(ckt->CKTrhs + here->BSIM4v6dNode) -= here->BSIM4v6rhsd); 
               (*(ckt->CKTrhs + here->BSIM4v6sNode) += here->BSIM4v6rhss);
           }

           if (here->BSIM4v6trnqsMod)
               *(ckt->CKTrhs + here->BSIM4v6qNode) += here->BSIM4v6rhsq;


        /* Update A for Ax = b */
           if (here->BSIM4v6rgateMod == 1)
           {   (*(here->BSIM4v6GEgePtr) += here->BSIM4v6_1);
               (*(here->BSIM4v6GPgePtr) -= here->BSIM4v6_2);
               (*(here->BSIM4v6GEgpPtr) -= here->BSIM4v6_3);
               (*(here->BSIM4v6GPgpPtr) += here->BSIM4v6_4);
               (*(here->BSIM4v6GPdpPtr) += here->BSIM4v6_5);
               (*(here->BSIM4v6GPspPtr) += here->BSIM4v6_6);
               (*(here->BSIM4v6GPbpPtr) += here->BSIM4v6_7);
           }
           else if (here->BSIM4v6rgateMod == 2)        
           {   (*(here->BSIM4v6GEgePtr) += here->BSIM4v6_8);
               (*(here->BSIM4v6GEgpPtr) += here->BSIM4v6_9);
               (*(here->BSIM4v6GEdpPtr) += here->BSIM4v6_10);
               (*(here->BSIM4v6GEspPtr) += here->BSIM4v6_11);
               (*(here->BSIM4v6GEbpPtr) += here->BSIM4v6_12);        

               (*(here->BSIM4v6GPgePtr) -= here->BSIM4v6_13);
               (*(here->BSIM4v6GPgpPtr) += here->BSIM4v6_14);
               (*(here->BSIM4v6GPdpPtr) += here->BSIM4v6_15);
               (*(here->BSIM4v6GPspPtr) += here->BSIM4v6_16);
               (*(here->BSIM4v6GPbpPtr) += here->BSIM4v6_17);
           }
           else if (here->BSIM4v6rgateMod == 3)
           {   (*(here->BSIM4v6GEgePtr) += here->BSIM4v6_18);
               (*(here->BSIM4v6GEgmPtr) -= here->BSIM4v6_19);
               (*(here->BSIM4v6GMgePtr) -= here->BSIM4v6_20);
               (*(here->BSIM4v6GMgmPtr) += here->BSIM4v6_21);

               (*(here->BSIM4v6GMdpPtr) += here->BSIM4v6_22);
               (*(here->BSIM4v6GMgpPtr) += here->BSIM4v6_23);
               (*(here->BSIM4v6GMspPtr) += here->BSIM4v6_24);
               (*(here->BSIM4v6GMbpPtr) += here->BSIM4v6_25);

               (*(here->BSIM4v6DPgmPtr) += here->BSIM4v6_26);
               (*(here->BSIM4v6GPgmPtr) -= here->BSIM4v6_27);
               (*(here->BSIM4v6SPgmPtr) += here->BSIM4v6_28);
               (*(here->BSIM4v6BPgmPtr) += here->BSIM4v6_29);

               (*(here->BSIM4v6GPgpPtr) += here->BSIM4v6_30);
               (*(here->BSIM4v6GPdpPtr) += here->BSIM4v6_31);
               (*(here->BSIM4v6GPspPtr) += here->BSIM4v6_32);
               (*(here->BSIM4v6GPbpPtr) += here->BSIM4v6_33);
           }


            else
           {   (*(here->BSIM4v6GPgpPtr) += here->BSIM4v6_34);
               (*(here->BSIM4v6GPdpPtr) += here->BSIM4v6_35);
               (*(here->BSIM4v6GPspPtr) += here->BSIM4v6_36);
               (*(here->BSIM4v6GPbpPtr) += here->BSIM4v6_37);
           }


           if (model->BSIM4v6rdsMod)
           {   (*(here->BSIM4v6DgpPtr) += here->BSIM4v6_38);
               (*(here->BSIM4v6DspPtr) += here->BSIM4v6_39);
               (*(here->BSIM4v6DbpPtr) += here->BSIM4v6_40);
               (*(here->BSIM4v6SdpPtr) += here->BSIM4v6_41);
               (*(here->BSIM4v6SgpPtr) += here->BSIM4v6_42);
               (*(here->BSIM4v6SbpPtr) += here->BSIM4v6_43);
           }

           (*(here->BSIM4v6DPdpPtr) += here->BSIM4v6_44);
           (*(here->BSIM4v6DPdPtr) -= here->BSIM4v6_45);
           (*(here->BSIM4v6DPgpPtr) += here->BSIM4v6_46);
           (*(here->BSIM4v6DPspPtr) -= here->BSIM4v6_47);
           (*(here->BSIM4v6DPbpPtr) -= here->BSIM4v6_48);

           (*(here->BSIM4v6DdpPtr) -= here->BSIM4v6_49);
           (*(here->BSIM4v6DdPtr) += here->BSIM4v6_50);

           (*(here->BSIM4v6SPdpPtr) -= here->BSIM4v6_51);
           (*(here->BSIM4v6SPgpPtr) += here->BSIM4v6_52);
           (*(here->BSIM4v6SPspPtr) += here->BSIM4v6_53);
           (*(here->BSIM4v6SPsPtr) -= here->BSIM4v6_54);
           (*(here->BSIM4v6SPbpPtr) -= here->BSIM4v6_55);

           (*(here->BSIM4v6SspPtr) -= here->BSIM4v6_56);
           (*(here->BSIM4v6SsPtr) += here->BSIM4v6_57);

           (*(here->BSIM4v6BPdpPtr) += here->BSIM4v6_58);
           (*(here->BSIM4v6BPgpPtr) += here->BSIM4v6_59);
           (*(here->BSIM4v6BPspPtr) += here->BSIM4v6_60);
           (*(here->BSIM4v6BPbpPtr) += here->BSIM4v6_61);

           /* stamp gidl */
           (*(here->BSIM4v6DPdpPtr) += here->BSIM4v6_62);
           (*(here->BSIM4v6DPgpPtr) += here->BSIM4v6_63);
           (*(here->BSIM4v6DPspPtr) -= here->BSIM4v6_64);
           (*(here->BSIM4v6DPbpPtr) += here->BSIM4v6_65);
           (*(here->BSIM4v6BPdpPtr) -= here->BSIM4v6_66);
           (*(here->BSIM4v6BPgpPtr) -= here->BSIM4v6_67);
           (*(here->BSIM4v6BPspPtr) += here->BSIM4v6_68);
           (*(here->BSIM4v6BPbpPtr) -= here->BSIM4v6_69);
            /* stamp gisl */
           (*(here->BSIM4v6SPdpPtr) -= here->BSIM4v6_70);
           (*(here->BSIM4v6SPgpPtr) += here->BSIM4v6_71);
           (*(here->BSIM4v6SPspPtr) += here->BSIM4v6_72);
           (*(here->BSIM4v6SPbpPtr) += here->BSIM4v6_73);
           (*(here->BSIM4v6BPdpPtr) += here->BSIM4v6_74);
           (*(here->BSIM4v6BPgpPtr) -= here->BSIM4v6_75);
           (*(here->BSIM4v6BPspPtr) -= here->BSIM4v6_76);
           (*(here->BSIM4v6BPbpPtr) -= here->BSIM4v6_77);


           if (here->BSIM4v6rbodyMod)
           {   (*(here->BSIM4v6DPdbPtr) += here->BSIM4v6_78);
               (*(here->BSIM4v6SPsbPtr) -= here->BSIM4v6_79);

               (*(here->BSIM4v6DBdpPtr) += here->BSIM4v6_80);
               (*(here->BSIM4v6DBdbPtr) += here->BSIM4v6_81);
               (*(here->BSIM4v6DBbpPtr) -= here->BSIM4v6_82);
               (*(here->BSIM4v6DBbPtr) -= here->BSIM4v6_83);

               (*(here->BSIM4v6BPdbPtr) -= here->BSIM4v6_84);
               (*(here->BSIM4v6BPbPtr) -= here->BSIM4v6_85);
               (*(here->BSIM4v6BPsbPtr) -= here->BSIM4v6_86);
               (*(here->BSIM4v6BPbpPtr) += here->BSIM4v6_87);

               (*(here->BSIM4v6SBspPtr) += here->BSIM4v6_88);
               (*(here->BSIM4v6SBbpPtr) -= here->BSIM4v6_89);
               (*(here->BSIM4v6SBbPtr) -= here->BSIM4v6_90);
               (*(here->BSIM4v6SBsbPtr) += here->BSIM4v6_91);

               (*(here->BSIM4v6BdbPtr) -= here->BSIM4v6_92);
               (*(here->BSIM4v6BbpPtr) -= here->BSIM4v6_93);
               (*(here->BSIM4v6BsbPtr) -= here->BSIM4v6_94);
               (*(here->BSIM4v6BbPtr) += here->BSIM4v6_95);
           }

           if (here->BSIM4v6trnqsMod)
           {   (*(here->BSIM4v6QqPtr) += here->BSIM4v6_96);
               (*(here->BSIM4v6QgpPtr) += here->BSIM4v6_97);
               (*(here->BSIM4v6QdpPtr) += here->BSIM4v6_98);
               (*(here->BSIM4v6QspPtr) += here->BSIM4v6_99);
               (*(here->BSIM4v6QbpPtr) += here->BSIM4v6_100);

               (*(here->BSIM4v6DPqPtr) += here->BSIM4v6_101);
               (*(here->BSIM4v6SPqPtr) += here->BSIM4v6_102);
               (*(here->BSIM4v6GPqPtr) -= here->BSIM4v6_103);
           }
    }
}

#endif
