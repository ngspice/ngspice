/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4temp.c of BSIM4v7.0.
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
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

#define Kb 1.3806226e-23
#define KboQ 8.617087e-5
#define EPS0 8.85418e-12
#define EPSSI 1.03594e-10
#define PI 3.141592654
#define MAX_EXP 5.834617425e14
#define MIN_EXP 1.713908431e-15
#define EXP_THRESHOLD 34.0
#define Charge_q 1.60219e-19
#define DELTA  1.0E-9
#define DEXP(A,B) {                                                        \
        if (A > EXP_THRESHOLD) {                                           \
            B = MAX_EXP*(1.0+(A)-EXP_THRESHOLD);                           \
        } else if (A < -EXP_THRESHOLD)  {                                  \
            B = MIN_EXP;                                                   \
        } else   {                                                         \
            B = exp(A);                                                    \
        }                                                                  \
    }

static int
BSIM4v7DioIjthVjmEval(
double Nvtm, double Ijth, double Isb, double XExpBV,
double *Vjm)
{
double Tb, Tc, EVjmovNv;

       Tc = XExpBV;
       Tb = 1.0 + Ijth / Isb - Tc;
       EVjmovNv = 0.5 * (Tb + sqrt(Tb * Tb + 4.0 * Tc));
       *Vjm = Nvtm * log(EVjmovNv);

return 0;
}


int
BSIM4v7temp(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4v7model *model = (BSIM4v7model*) inModel;
BSIM4v7instance *here;
struct bsim4SizeDependParam *pSizeDependParamKnot, *pLastKnot, *pParam=NULL;
double tmp, tmp1, tmp2, tmp3, Eg, Eg0, ni, epssub;
double T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Lnew=0.0, Wnew;
double delTemp, Temp, TRatio, Inv_L, Inv_W, Inv_LW, Vtm0, Tnom;
double dumPs, dumPd, dumAs, dumAd, PowWeffWr;
double DMCGeff, DMCIeff, DMDGeff;
double Nvtms, Nvtmd, SourceSatCurrent, DrainSatCurrent;
double T10, T11;
double Inv_saref, Inv_sbref, Inv_sa, Inv_sb, rho, Ldrn, dvth0_lod;
double W_tmp, Inv_ODeff, OD_offset, dk2_lod, deta0_lod;
double lnl, lnw, lnnf, rbpbx, rbpby, rbsbx, rbsby, rbdbx, rbdby,bodymode;
double kvsat, wlod, sceff, Wdrn;
double V0, lt1, ltw, Theta0, Delt_vth, Vth_NarrowW, Lpe_Vb, Vth;
double n,n0, Vgsteff, Vgs_eff, niter, toxpf, toxpi, Tcen, toxe, epsrox, vddeot;
double vtfbphi2eot, phieot, TempRatioeot, Vtm0eot, Vtmeot,vbieot;

int Size_Not_Found, i;

    /*  loop through all the BSIM4v7 device models */
    for (; model != NULL; model = model->BSIM4v7nextModel)
    {    Temp = ckt->CKTtemp;
         if (model->BSIM4v7SbulkJctPotential < 0.1)
         {   model->BSIM4v7SbulkJctPotential = 0.1;
             fprintf(stderr, "Given pbs is less than 0.1. Pbs is set to 0.1.\n");
         }
         if (model->BSIM4v7SsidewallJctPotential < 0.1)
         {   model->BSIM4v7SsidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbsws is less than 0.1. Pbsws is set to 0.1.\n");
         }
         if (model->BSIM4v7SGatesidewallJctPotential < 0.1)
         {   model->BSIM4v7SGatesidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswgs is less than 0.1. Pbswgs is set to 0.1.\n");
         }

         if (model->BSIM4v7DbulkJctPotential < 0.1)
         {   model->BSIM4v7DbulkJctPotential = 0.1;
             fprintf(stderr, "Given pbd is less than 0.1. Pbd is set to 0.1.\n");
         }
         if (model->BSIM4v7DsidewallJctPotential < 0.1)
         {   model->BSIM4v7DsidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswd is less than 0.1. Pbswd is set to 0.1.\n");
         }
         if (model->BSIM4v7DGatesidewallJctPotential < 0.1)
         {   model->BSIM4v7DGatesidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswgd is less than 0.1. Pbswgd is set to 0.1.\n");
         }

         if(model->BSIM4v7mtrlMod == 0)
         {
             if ((model->BSIM4v7toxeGiven) && (model->BSIM4v7toxpGiven) && (model->BSIM4v7dtoxGiven)
                 && (model->BSIM4v7toxe != (model->BSIM4v7toxp + model->BSIM4v7dtox)))
                 printf("Warning: toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored.\n");
             else if ((model->BSIM4v7toxeGiven) && (!model->BSIM4v7toxpGiven))
               model->BSIM4v7toxp = model->BSIM4v7toxe - model->BSIM4v7dtox;
             else if ((!model->BSIM4v7toxeGiven) && (model->BSIM4v7toxpGiven)){
               model->BSIM4v7toxe = model->BSIM4v7toxp + model->BSIM4v7dtox;
                 if (!model->BSIM4v7toxmGiven)                        /* v4.7 */
                     model->BSIM4v7toxm = model->BSIM4v7toxe;
             }
         }
         else if(model->BSIM4v7mtrlCompatMod != 0) /* v4.7 */
         {
             T0 = model->BSIM4v7epsrox / 3.9;
             if ((model->BSIM4v7eotGiven) && (model->BSIM4v7toxpGiven) && (model->BSIM4v7dtoxGiven)
                 && (ABS(model->BSIM4v7eot * T0 - (model->BSIM4v7toxp + model->BSIM4v7dtox)) > 1.0e-20))
             {
                 printf("Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.\n");
             }
             else if ((model->BSIM4v7eotGiven) && (!model->BSIM4v7toxpGiven))
               model->BSIM4v7toxp = T0 * model->BSIM4v7eot - model->BSIM4v7dtox;
             else if ((!model->BSIM4v7eotGiven) && (model->BSIM4v7toxpGiven)){
               model->BSIM4v7eot = (model->BSIM4v7toxp + model->BSIM4v7dtox) / T0;
                 if (!model->BSIM4v7toxmGiven)
                     model->BSIM4v7toxm = model->BSIM4v7eot;
             }
         }

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


         model->BSIM4v7coxe = epsrox * EPS0 / toxe;
         if(model->BSIM4v7mtrlMod == 0 || model->BSIM4v7mtrlCompatMod != 0)
           model->BSIM4v7coxp = model->BSIM4v7epsrox * EPS0 / model->BSIM4v7toxp;

         if (!model->BSIM4v7cgdoGiven)
         {   if (model->BSIM4v7dlcGiven && (model->BSIM4v7dlc > 0.0))
                 model->BSIM4v7cgdo = model->BSIM4v7dlc * model->BSIM4v7coxe
                                  - model->BSIM4v7cgdl ;
             else
                 model->BSIM4v7cgdo = 0.6 * model->BSIM4v7xj * model->BSIM4v7coxe;
         }
         if (!model->BSIM4v7cgsoGiven)
         {   if (model->BSIM4v7dlcGiven && (model->BSIM4v7dlc > 0.0))
                 model->BSIM4v7cgso = model->BSIM4v7dlc * model->BSIM4v7coxe
                                  - model->BSIM4v7cgsl ;
             else
                 model->BSIM4v7cgso = 0.6 * model->BSIM4v7xj * model->BSIM4v7coxe;
         }
         if (!model->BSIM4v7cgboGiven)
             model->BSIM4v7cgbo = 2.0 * model->BSIM4v7dwc * model->BSIM4v7coxe;
         model->pSizeDependParamKnot = NULL;
         pLastKnot = NULL;

         Tnom = model->BSIM4v7tnom;
         TRatio = Temp / Tnom;

         model->BSIM4v7vcrit = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
         model->BSIM4v7factor1 = sqrt(epssub / (epsrox * EPS0)* toxe);

         Vtm0 = model->BSIM4v7vtm0 = KboQ * Tnom;

         if(model->BSIM4v7mtrlMod==0)
         {
             Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
             ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15)
                 * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
         }
         else
         {
           Eg0 = model->BSIM4v7bg0sub - model->BSIM4v7tbgasub * Tnom * Tnom
                                      / (Tnom + model->BSIM4v7tbgbsub);
           T0 =  model->BSIM4v7bg0sub - model->BSIM4v7tbgasub * 90090.0225
                                      / (300.15 + model->BSIM4v7tbgbsub);
           ni = model->BSIM4v7ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15)
                 * exp((T0 - Eg0) / (2.0 * Vtm0));
         }

         model->BSIM4v7Eg0 = Eg0;
         model->BSIM4v7vtm = KboQ * Temp;
         if(model->BSIM4v7mtrlMod == 0)
           Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
         else
           Eg = model->BSIM4v7bg0sub - model->BSIM4v7tbgasub * Temp * Temp
                                      / (Temp + model->BSIM4v7tbgbsub);
         if (Temp != Tnom)
         {   T0 = Eg0 / Vtm0 - Eg / model->BSIM4v7vtm;
             T1 = log(Temp / Tnom);
             T2 = T0 + model->BSIM4v7SjctTempExponent * T1;
             T3 = exp(T2 / model->BSIM4v7SjctEmissionCoeff);
             model->BSIM4v7SjctTempSatCurDensity = model->BSIM4v7SjctSatCurDensity
                                               * T3;
             model->BSIM4v7SjctSidewallTempSatCurDensity
                         = model->BSIM4v7SjctSidewallSatCurDensity * T3;
             model->BSIM4v7SjctGateSidewallTempSatCurDensity
                         = model->BSIM4v7SjctGateSidewallSatCurDensity * T3;

             T2 = T0 + model->BSIM4v7DjctTempExponent * T1;
             T3 = exp(T2 / model->BSIM4v7DjctEmissionCoeff);
             model->BSIM4v7DjctTempSatCurDensity = model->BSIM4v7DjctSatCurDensity
                                               * T3;
             model->BSIM4v7DjctSidewallTempSatCurDensity
                         = model->BSIM4v7DjctSidewallSatCurDensity * T3;
             model->BSIM4v7DjctGateSidewallTempSatCurDensity
                         = model->BSIM4v7DjctGateSidewallSatCurDensity * T3;
         }
         else
         {   model->BSIM4v7SjctTempSatCurDensity = model->BSIM4v7SjctSatCurDensity;
             model->BSIM4v7SjctSidewallTempSatCurDensity
                        = model->BSIM4v7SjctSidewallSatCurDensity;
             model->BSIM4v7SjctGateSidewallTempSatCurDensity
                        = model->BSIM4v7SjctGateSidewallSatCurDensity;
             model->BSIM4v7DjctTempSatCurDensity = model->BSIM4v7DjctSatCurDensity;
             model->BSIM4v7DjctSidewallTempSatCurDensity
                        = model->BSIM4v7DjctSidewallSatCurDensity;
             model->BSIM4v7DjctGateSidewallTempSatCurDensity
                        = model->BSIM4v7DjctGateSidewallSatCurDensity;
         }

         if (model->BSIM4v7SjctTempSatCurDensity < 0.0)
             model->BSIM4v7SjctTempSatCurDensity = 0.0;
         if (model->BSIM4v7SjctSidewallTempSatCurDensity < 0.0)
             model->BSIM4v7SjctSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4v7SjctGateSidewallTempSatCurDensity < 0.0)
             model->BSIM4v7SjctGateSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4v7DjctTempSatCurDensity < 0.0)
             model->BSIM4v7DjctTempSatCurDensity = 0.0;
         if (model->BSIM4v7DjctSidewallTempSatCurDensity < 0.0)
             model->BSIM4v7DjctSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4v7DjctGateSidewallTempSatCurDensity < 0.0)
             model->BSIM4v7DjctGateSidewallTempSatCurDensity = 0.0;

         /* Temperature dependence of D/B and S/B diode capacitance begins */
         delTemp = ckt->CKTtemp - model->BSIM4v7tnom;
         T0 = model->BSIM4v7tcj * delTemp;
         if (T0 >= -1.0)
         {   model->BSIM4v7SunitAreaTempJctCap = model->BSIM4v7SunitAreaJctCap *(1.0 + T0); /*bug_fix -JX */
             model->BSIM4v7DunitAreaTempJctCap = model->BSIM4v7DunitAreaJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4v7SunitAreaJctCap > 0.0)
             {   model->BSIM4v7SunitAreaTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjs to be negative. Cjs is clamped to zero.\n");
             }
             if (model->BSIM4v7DunitAreaJctCap > 0.0)
             {   model->BSIM4v7DunitAreaTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjd to be negative. Cjd is clamped to zero.\n");
             }
         }
         T0 = model->BSIM4v7tcjsw * delTemp;
                   if (model->BSIM4v7SunitLengthSidewallJctCap < 0.0)/*4.6.2*/
                      {model->BSIM4v7SunitLengthSidewallJctCap = 0.0;
                           fprintf(stderr, "CJSWS is negative. Cjsws is clamped to zero.\n");}
                  if (model->BSIM4v7DunitLengthSidewallJctCap < 0.0)
                      {model->BSIM4v7DunitLengthSidewallJctCap = 0.0;
                           fprintf(stderr, "CJSWD is negative. Cjswd is clamped to zero.\n");}
         if (T0 >= -1.0)
         {   model->BSIM4v7SunitLengthSidewallTempJctCap = model->BSIM4v7SunitLengthSidewallJctCap *(1.0 + T0);
             model->BSIM4v7DunitLengthSidewallTempJctCap = model->BSIM4v7DunitLengthSidewallJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4v7SunitLengthSidewallJctCap > 0.0)
             {   model->BSIM4v7SunitLengthSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.\n");
             }
             if (model->BSIM4v7DunitLengthSidewallJctCap > 0.0)
             {   model->BSIM4v7DunitLengthSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.\n");
             }
         }
         T0 = model->BSIM4v7tcjswg * delTemp;
         if (T0 >= -1.0)
         {   model->BSIM4v7SunitLengthGateSidewallTempJctCap = model->BSIM4v7SunitLengthGateSidewallJctCap *(1.0 + T0);
             model->BSIM4v7DunitLengthGateSidewallTempJctCap = model->BSIM4v7DunitLengthGateSidewallJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4v7SunitLengthGateSidewallJctCap > 0.0)
             {   model->BSIM4v7SunitLengthGateSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.\n");
             }
             if (model->BSIM4v7DunitLengthGateSidewallJctCap > 0.0)
             {   model->BSIM4v7DunitLengthGateSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.\n");
             }
         }

         model->BSIM4v7PhiBS = model->BSIM4v7SbulkJctPotential
                           - model->BSIM4v7tpb * delTemp;
         if (model->BSIM4v7PhiBS < 0.01)
         {   model->BSIM4v7PhiBS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.\n");
         }
         model->BSIM4v7PhiBD = model->BSIM4v7DbulkJctPotential
                           - model->BSIM4v7tpb * delTemp;
         if (model->BSIM4v7PhiBD < 0.01)
         {   model->BSIM4v7PhiBD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.\n");
         }

         model->BSIM4v7PhiBSWS = model->BSIM4v7SsidewallJctPotential
                             - model->BSIM4v7tpbsw * delTemp;
         if (model->BSIM4v7PhiBSWS <= 0.01)
         {   model->BSIM4v7PhiBSWS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.\n");
         }
         model->BSIM4v7PhiBSWD = model->BSIM4v7DsidewallJctPotential
                             - model->BSIM4v7tpbsw * delTemp;
         if (model->BSIM4v7PhiBSWD <= 0.01)
         {   model->BSIM4v7PhiBSWD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.\n");
         }

         model->BSIM4v7PhiBSWGS = model->BSIM4v7SGatesidewallJctPotential
                              - model->BSIM4v7tpbswg * delTemp;
         if (model->BSIM4v7PhiBSWGS <= 0.01)
         {   model->BSIM4v7PhiBSWGS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.\n");
         }
         model->BSIM4v7PhiBSWGD = model->BSIM4v7DGatesidewallJctPotential
                              - model->BSIM4v7tpbswg * delTemp;
         if (model->BSIM4v7PhiBSWGD <= 0.01)
         {   model->BSIM4v7PhiBSWGD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.\n");
         } /* End of junction capacitance */


         if (model->BSIM4v7ijthdfwd <= 0.0)
         {   model->BSIM4v7ijthdfwd = 0.0;
             fprintf(stderr, "Ijthdfwd reset to %g.\n", model->BSIM4v7ijthdfwd);
         }
         if (model->BSIM4v7ijthsfwd <= 0.0)
         {   model->BSIM4v7ijthsfwd = 0.0;
             fprintf(stderr, "Ijthsfwd reset to %g.\n", model->BSIM4v7ijthsfwd);
         }
         if (model->BSIM4v7ijthdrev <= 0.0)
         {   model->BSIM4v7ijthdrev = 0.0;
             fprintf(stderr, "Ijthdrev reset to %g.\n", model->BSIM4v7ijthdrev);
         }
         if (model->BSIM4v7ijthsrev <= 0.0)
         {   model->BSIM4v7ijthsrev = 0.0;
             fprintf(stderr, "Ijthsrev reset to %g.\n", model->BSIM4v7ijthsrev);
         }

         if ((model->BSIM4v7xjbvd <= 0.0) && (model->BSIM4v7dioMod == 2))
         {   model->BSIM4v7xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model->BSIM4v7xjbvd);
         }
         else if ((model->BSIM4v7xjbvd < 0.0) && (model->BSIM4v7dioMod == 0))
         {   model->BSIM4v7xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model->BSIM4v7xjbvd);
         }

         if (model->BSIM4v7bvd <= 0.0)   /*4.6.2*/
         {   model->BSIM4v7bvd = 0.0;
             fprintf(stderr, "BVD reset to %g.\n", model->BSIM4v7bvd);
         }

         if ((model->BSIM4v7xjbvs <= 0.0) && (model->BSIM4v7dioMod == 2))
         {   model->BSIM4v7xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model->BSIM4v7xjbvs);
         }
         else if ((model->BSIM4v7xjbvs < 0.0) && (model->BSIM4v7dioMod == 0))
         {   model->BSIM4v7xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model->BSIM4v7xjbvs);
         }

         if (model->BSIM4v7bvs <= 0.0)
         {   model->BSIM4v7bvs = 0.0;
             fprintf(stderr, "BVS reset to %g.\n", model->BSIM4v7bvs);
         }


         /* loop through all the instances of the model */
         for (here = model->BSIM4v7instances; here != NULL;
              here = here->BSIM4v7nextInstance)
         {
              pSizeDependParamKnot = model->pSizeDependParamKnot;
              Size_Not_Found = 1;
              while ((pSizeDependParamKnot != NULL) && Size_Not_Found)
              {   if ((here->BSIM4v7l == pSizeDependParamKnot->Length)
                      && (here->BSIM4v7w == pSizeDependParamKnot->Width)
                      && (here->BSIM4v7nf == pSizeDependParamKnot->NFinger))
                  {   Size_Not_Found = 0;
                      here->pParam = pSizeDependParamKnot;
                      pParam = here->pParam; /*bug-fix  */
                  }
                  else
                  {   pLastKnot = pSizeDependParamKnot;
                      pSizeDependParamKnot = pSizeDependParamKnot->pNext;
                  }
              }

              /* stress effect */
              Ldrn = here->BSIM4v7l;
              Wdrn = here->BSIM4v7w / here->BSIM4v7nf;

              if (Size_Not_Found)
              {   pParam = TMALLOC(struct bsim4SizeDependParam, 1);
                  if (pLastKnot == NULL)
                      model->pSizeDependParamKnot = pParam;
                  else
                      pLastKnot->pNext = pParam;
                  pParam->pNext = NULL;
                  here->pParam = pParam;

                  pParam->Length = here->BSIM4v7l;
                  pParam->Width = here->BSIM4v7w;
                  pParam->NFinger = here->BSIM4v7nf;
                  Lnew = here->BSIM4v7l  + model->BSIM4v7xl ;
                  Wnew = here->BSIM4v7w / here->BSIM4v7nf + model->BSIM4v7xw;

                  T0 = pow(Lnew, model->BSIM4v7Lln);
                  T1 = pow(Wnew, model->BSIM4v7Lwn);
                  tmp1 = model->BSIM4v7Ll / T0 + model->BSIM4v7Lw / T1
                       + model->BSIM4v7Lwl / (T0 * T1);
                  pParam->BSIM4v7dl = model->BSIM4v7Lint + tmp1;
                  tmp2 = model->BSIM4v7Llc / T0 + model->BSIM4v7Lwc / T1
                       + model->BSIM4v7Lwlc / (T0 * T1);
                  pParam->BSIM4v7dlc = model->BSIM4v7dlc + tmp2;

                  T2 = pow(Lnew, model->BSIM4v7Wln);
                  T3 = pow(Wnew, model->BSIM4v7Wwn);
                  tmp1 = model->BSIM4v7Wl / T2 + model->BSIM4v7Ww / T3
                       + model->BSIM4v7Wwl / (T2 * T3);
                  pParam->BSIM4v7dw = model->BSIM4v7Wint + tmp1;
                  tmp2 = model->BSIM4v7Wlc / T2 + model->BSIM4v7Wwc / T3
                       + model->BSIM4v7Wwlc / (T2 * T3);
                  pParam->BSIM4v7dwc = model->BSIM4v7dwc + tmp2;
                  pParam->BSIM4v7dwj = model->BSIM4v7dwj + tmp2;

                  pParam->BSIM4v7leff = Lnew - 2.0 * pParam->BSIM4v7dl;
                  if (pParam->BSIM4v7leff <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4v7: mosfet %s, model %s: Effective channel length <= 0",
                       model->BSIM4v7modName, here->BSIM4v7name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4v7weff = Wnew - 2.0 * pParam->BSIM4v7dw;
                  if (pParam->BSIM4v7weff <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4v7: mosfet %s, model %s: Effective channel width <= 0",
                       model->BSIM4v7modName, here->BSIM4v7name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4v7leffCV = Lnew - 2.0 * pParam->BSIM4v7dlc;
                  if (pParam->BSIM4v7leffCV <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4v7: mosfet %s, model %s: Effective channel length for C-V <= 0",
                       model->BSIM4v7modName, here->BSIM4v7name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4v7weffCV = Wnew - 2.0 * pParam->BSIM4v7dwc;
                  if (pParam->BSIM4v7weffCV <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4v7: mosfet %s, model %s: Effective channel width for C-V <= 0",
                       model->BSIM4v7modName, here->BSIM4v7name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4v7weffCJ = Wnew - 2.0 * pParam->BSIM4v7dwj;
                  if (pParam->BSIM4v7weffCJ <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4v7: mosfet %s, model %s: Effective channel width for S/D junctions <= 0",
                       model->BSIM4v7modName, here->BSIM4v7name);
                      return(E_BADPARM);
                  }


                  if (model->BSIM4v7binUnit == 1)
                  {   Inv_L = 1.0e-6 / pParam->BSIM4v7leff;
                      Inv_W = 1.0e-6 / pParam->BSIM4v7weff;
                      Inv_LW = 1.0e-12 / (pParam->BSIM4v7leff
                             * pParam->BSIM4v7weff);
                  }
                  else
                  {   Inv_L = 1.0 / pParam->BSIM4v7leff;
                      Inv_W = 1.0 / pParam->BSIM4v7weff;
                      Inv_LW = 1.0 / (pParam->BSIM4v7leff
                             * pParam->BSIM4v7weff);
                  }
                  pParam->BSIM4v7cdsc = model->BSIM4v7cdsc
                                    + model->BSIM4v7lcdsc * Inv_L
                                    + model->BSIM4v7wcdsc * Inv_W
                                    + model->BSIM4v7pcdsc * Inv_LW;
                  pParam->BSIM4v7cdscb = model->BSIM4v7cdscb
                                     + model->BSIM4v7lcdscb * Inv_L
                                     + model->BSIM4v7wcdscb * Inv_W
                                     + model->BSIM4v7pcdscb * Inv_LW;

                      pParam->BSIM4v7cdscd = model->BSIM4v7cdscd
                                     + model->BSIM4v7lcdscd * Inv_L
                                     + model->BSIM4v7wcdscd * Inv_W
                                     + model->BSIM4v7pcdscd * Inv_LW;

                  pParam->BSIM4v7cit = model->BSIM4v7cit
                                   + model->BSIM4v7lcit * Inv_L
                                   + model->BSIM4v7wcit * Inv_W
                                   + model->BSIM4v7pcit * Inv_LW;
                  pParam->BSIM4v7nfactor = model->BSIM4v7nfactor
                                       + model->BSIM4v7lnfactor * Inv_L
                                       + model->BSIM4v7wnfactor * Inv_W
                                       + model->BSIM4v7pnfactor * Inv_LW;
                  pParam->BSIM4v7tnfactor = model->BSIM4v7tnfactor                        /* v4.7 */
                                       + model->BSIM4v7ltnfactor * Inv_L
                                       + model->BSIM4v7wtnfactor * Inv_W
                                       + model->BSIM4v7ptnfactor * Inv_LW;
                  pParam->BSIM4v7xj = model->BSIM4v7xj
                                  + model->BSIM4v7lxj * Inv_L
                                  + model->BSIM4v7wxj * Inv_W
                                  + model->BSIM4v7pxj * Inv_LW;
                  pParam->BSIM4v7vsat = model->BSIM4v7vsat
                                    + model->BSIM4v7lvsat * Inv_L
                                    + model->BSIM4v7wvsat * Inv_W
                                    + model->BSIM4v7pvsat * Inv_LW;
                  pParam->BSIM4v7at = model->BSIM4v7at
                                  + model->BSIM4v7lat * Inv_L
                                  + model->BSIM4v7wat * Inv_W
                                  + model->BSIM4v7pat * Inv_LW;
                  pParam->BSIM4v7a0 = model->BSIM4v7a0
                                  + model->BSIM4v7la0 * Inv_L
                                  + model->BSIM4v7wa0 * Inv_W
                                  + model->BSIM4v7pa0 * Inv_LW;

                  pParam->BSIM4v7ags = model->BSIM4v7ags
                                  + model->BSIM4v7lags * Inv_L
                                  + model->BSIM4v7wags * Inv_W
                                  + model->BSIM4v7pags * Inv_LW;

                  pParam->BSIM4v7a1 = model->BSIM4v7a1
                                  + model->BSIM4v7la1 * Inv_L
                                  + model->BSIM4v7wa1 * Inv_W
                                  + model->BSIM4v7pa1 * Inv_LW;
                  pParam->BSIM4v7a2 = model->BSIM4v7a2
                                  + model->BSIM4v7la2 * Inv_L
                                  + model->BSIM4v7wa2 * Inv_W
                                  + model->BSIM4v7pa2 * Inv_LW;
                  pParam->BSIM4v7keta = model->BSIM4v7keta
                                    + model->BSIM4v7lketa * Inv_L
                                    + model->BSIM4v7wketa * Inv_W
                                    + model->BSIM4v7pketa * Inv_LW;
                  pParam->BSIM4v7nsub = model->BSIM4v7nsub
                                    + model->BSIM4v7lnsub * Inv_L
                                    + model->BSIM4v7wnsub * Inv_W
                                    + model->BSIM4v7pnsub * Inv_LW;
                  pParam->BSIM4v7ndep = model->BSIM4v7ndep
                                    + model->BSIM4v7lndep * Inv_L
                                    + model->BSIM4v7wndep * Inv_W
                                    + model->BSIM4v7pndep * Inv_LW;
                  pParam->BSIM4v7nsd = model->BSIM4v7nsd
                                   + model->BSIM4v7lnsd * Inv_L
                                   + model->BSIM4v7wnsd * Inv_W
                                   + model->BSIM4v7pnsd * Inv_LW;
                  pParam->BSIM4v7phin = model->BSIM4v7phin
                                    + model->BSIM4v7lphin * Inv_L
                                    + model->BSIM4v7wphin * Inv_W
                                    + model->BSIM4v7pphin * Inv_LW;
                  pParam->BSIM4v7ngate = model->BSIM4v7ngate
                                     + model->BSIM4v7lngate * Inv_L
                                     + model->BSIM4v7wngate * Inv_W
                                     + model->BSIM4v7pngate * Inv_LW;
                  pParam->BSIM4v7gamma1 = model->BSIM4v7gamma1
                                      + model->BSIM4v7lgamma1 * Inv_L
                                      + model->BSIM4v7wgamma1 * Inv_W
                                      + model->BSIM4v7pgamma1 * Inv_LW;
                  pParam->BSIM4v7gamma2 = model->BSIM4v7gamma2
                                      + model->BSIM4v7lgamma2 * Inv_L
                                      + model->BSIM4v7wgamma2 * Inv_W
                                      + model->BSIM4v7pgamma2 * Inv_LW;
                  pParam->BSIM4v7vbx = model->BSIM4v7vbx
                                   + model->BSIM4v7lvbx * Inv_L
                                   + model->BSIM4v7wvbx * Inv_W
                                   + model->BSIM4v7pvbx * Inv_LW;
                  pParam->BSIM4v7vbm = model->BSIM4v7vbm
                                   + model->BSIM4v7lvbm * Inv_L
                                   + model->BSIM4v7wvbm * Inv_W
                                   + model->BSIM4v7pvbm * Inv_LW;
                  pParam->BSIM4v7xt = model->BSIM4v7xt
                                   + model->BSIM4v7lxt * Inv_L
                                   + model->BSIM4v7wxt * Inv_W
                                   + model->BSIM4v7pxt * Inv_LW;
                  pParam->BSIM4v7vfb = model->BSIM4v7vfb
                                   + model->BSIM4v7lvfb * Inv_L
                                   + model->BSIM4v7wvfb * Inv_W
                                   + model->BSIM4v7pvfb * Inv_LW;
                  pParam->BSIM4v7k1 = model->BSIM4v7k1
                                  + model->BSIM4v7lk1 * Inv_L
                                  + model->BSIM4v7wk1 * Inv_W
                                  + model->BSIM4v7pk1 * Inv_LW;
                  pParam->BSIM4v7kt1 = model->BSIM4v7kt1
                                   + model->BSIM4v7lkt1 * Inv_L
                                   + model->BSIM4v7wkt1 * Inv_W
                                   + model->BSIM4v7pkt1 * Inv_LW;
                  pParam->BSIM4v7kt1l = model->BSIM4v7kt1l
                                    + model->BSIM4v7lkt1l * Inv_L
                                    + model->BSIM4v7wkt1l * Inv_W
                                    + model->BSIM4v7pkt1l * Inv_LW;
                  pParam->BSIM4v7k2 = model->BSIM4v7k2
                                  + model->BSIM4v7lk2 * Inv_L
                                  + model->BSIM4v7wk2 * Inv_W
                                  + model->BSIM4v7pk2 * Inv_LW;
                  pParam->BSIM4v7kt2 = model->BSIM4v7kt2
                                   + model->BSIM4v7lkt2 * Inv_L
                                   + model->BSIM4v7wkt2 * Inv_W
                                   + model->BSIM4v7pkt2 * Inv_LW;
                  pParam->BSIM4v7k3 = model->BSIM4v7k3
                                  + model->BSIM4v7lk3 * Inv_L
                                  + model->BSIM4v7wk3 * Inv_W
                                  + model->BSIM4v7pk3 * Inv_LW;
                  pParam->BSIM4v7k3b = model->BSIM4v7k3b
                                   + model->BSIM4v7lk3b * Inv_L
                                   + model->BSIM4v7wk3b * Inv_W
                                   + model->BSIM4v7pk3b * Inv_LW;
                  pParam->BSIM4v7w0 = model->BSIM4v7w0
                                  + model->BSIM4v7lw0 * Inv_L
                                  + model->BSIM4v7ww0 * Inv_W
                                  + model->BSIM4v7pw0 * Inv_LW;
                  pParam->BSIM4v7lpe0 = model->BSIM4v7lpe0
                                    + model->BSIM4v7llpe0 * Inv_L
                                     + model->BSIM4v7wlpe0 * Inv_W
                                    + model->BSIM4v7plpe0 * Inv_LW;
                  pParam->BSIM4v7lpeb = model->BSIM4v7lpeb
                                    + model->BSIM4v7llpeb * Inv_L
                                    + model->BSIM4v7wlpeb * Inv_W
                                    + model->BSIM4v7plpeb * Inv_LW;
                  pParam->BSIM4v7dvtp0 = model->BSIM4v7dvtp0
                                     + model->BSIM4v7ldvtp0 * Inv_L
                                     + model->BSIM4v7wdvtp0 * Inv_W
                                     + model->BSIM4v7pdvtp0 * Inv_LW;
                  pParam->BSIM4v7dvtp1 = model->BSIM4v7dvtp1
                                     + model->BSIM4v7ldvtp1 * Inv_L
                                     + model->BSIM4v7wdvtp1 * Inv_W
                                     + model->BSIM4v7pdvtp1 * Inv_LW;
                  pParam->BSIM4v7dvtp2 = model->BSIM4v7dvtp2                 /* v4.7  */
                                     + model->BSIM4v7ldvtp2 * Inv_L
                                     + model->BSIM4v7wdvtp2 * Inv_W
                                     + model->BSIM4v7pdvtp2 * Inv_LW;
                  pParam->BSIM4v7dvtp3 = model->BSIM4v7dvtp3                 /* v4.7  */
                                     + model->BSIM4v7ldvtp3 * Inv_L
                                     + model->BSIM4v7wdvtp3 * Inv_W
                                     + model->BSIM4v7pdvtp3 * Inv_LW;
                  pParam->BSIM4v7dvtp4 = model->BSIM4v7dvtp4                 /* v4.7  */
                                     + model->BSIM4v7ldvtp4 * Inv_L
                                     + model->BSIM4v7wdvtp4 * Inv_W
                                     + model->BSIM4v7pdvtp4 * Inv_LW;
                  pParam->BSIM4v7dvtp5 = model->BSIM4v7dvtp5                 /* v4.7  */
                                     + model->BSIM4v7ldvtp5 * Inv_L
                                     + model->BSIM4v7wdvtp5 * Inv_W
                                     + model->BSIM4v7pdvtp5 * Inv_LW;
                  pParam->BSIM4v7dvt0 = model->BSIM4v7dvt0
                                    + model->BSIM4v7ldvt0 * Inv_L
                                    + model->BSIM4v7wdvt0 * Inv_W
                                    + model->BSIM4v7pdvt0 * Inv_LW;
                  pParam->BSIM4v7dvt1 = model->BSIM4v7dvt1
                                    + model->BSIM4v7ldvt1 * Inv_L
                                    + model->BSIM4v7wdvt1 * Inv_W
                                    + model->BSIM4v7pdvt1 * Inv_LW;
                  pParam->BSIM4v7dvt2 = model->BSIM4v7dvt2
                                    + model->BSIM4v7ldvt2 * Inv_L
                                    + model->BSIM4v7wdvt2 * Inv_W
                                    + model->BSIM4v7pdvt2 * Inv_LW;
                  pParam->BSIM4v7dvt0w = model->BSIM4v7dvt0w
                                    + model->BSIM4v7ldvt0w * Inv_L
                                    + model->BSIM4v7wdvt0w * Inv_W
                                    + model->BSIM4v7pdvt0w * Inv_LW;
                  pParam->BSIM4v7dvt1w = model->BSIM4v7dvt1w
                                    + model->BSIM4v7ldvt1w * Inv_L
                                    + model->BSIM4v7wdvt1w * Inv_W
                                    + model->BSIM4v7pdvt1w * Inv_LW;
                  pParam->BSIM4v7dvt2w = model->BSIM4v7dvt2w
                                    + model->BSIM4v7ldvt2w * Inv_L
                                    + model->BSIM4v7wdvt2w * Inv_W
                                    + model->BSIM4v7pdvt2w * Inv_LW;
                  pParam->BSIM4v7drout = model->BSIM4v7drout
                                     + model->BSIM4v7ldrout * Inv_L
                                     + model->BSIM4v7wdrout * Inv_W
                                     + model->BSIM4v7pdrout * Inv_LW;
                  pParam->BSIM4v7dsub = model->BSIM4v7dsub
                                    + model->BSIM4v7ldsub * Inv_L
                                    + model->BSIM4v7wdsub * Inv_W
                                    + model->BSIM4v7pdsub * Inv_LW;
                  pParam->BSIM4v7vth0 = model->BSIM4v7vth0
                                    + model->BSIM4v7lvth0 * Inv_L
                                    + model->BSIM4v7wvth0 * Inv_W
                                    + model->BSIM4v7pvth0 * Inv_LW;
                  pParam->BSIM4v7ua = model->BSIM4v7ua
                                  + model->BSIM4v7lua * Inv_L
                                  + model->BSIM4v7wua * Inv_W
                                  + model->BSIM4v7pua * Inv_LW;
                  pParam->BSIM4v7ua1 = model->BSIM4v7ua1
                                   + model->BSIM4v7lua1 * Inv_L
                                   + model->BSIM4v7wua1 * Inv_W
                                   + model->BSIM4v7pua1 * Inv_LW;
                  pParam->BSIM4v7ub = model->BSIM4v7ub
                                  + model->BSIM4v7lub * Inv_L
                                  + model->BSIM4v7wub * Inv_W
                                  + model->BSIM4v7pub * Inv_LW;
                  pParam->BSIM4v7ub1 = model->BSIM4v7ub1
                                   + model->BSIM4v7lub1 * Inv_L
                                   + model->BSIM4v7wub1 * Inv_W
                                   + model->BSIM4v7pub1 * Inv_LW;
                  pParam->BSIM4v7uc = model->BSIM4v7uc
                                  + model->BSIM4v7luc * Inv_L
                                  + model->BSIM4v7wuc * Inv_W
                                  + model->BSIM4v7puc * Inv_LW;
                  pParam->BSIM4v7uc1 = model->BSIM4v7uc1
                                   + model->BSIM4v7luc1 * Inv_L
                                   + model->BSIM4v7wuc1 * Inv_W
                                   + model->BSIM4v7puc1 * Inv_LW;
                  pParam->BSIM4v7ud = model->BSIM4v7ud
                                  + model->BSIM4v7lud * Inv_L
                                  + model->BSIM4v7wud * Inv_W
                                  + model->BSIM4v7pud * Inv_LW;
                  pParam->BSIM4v7ud1 = model->BSIM4v7ud1
                                  + model->BSIM4v7lud1 * Inv_L
                                  + model->BSIM4v7wud1 * Inv_W
                                  + model->BSIM4v7pud1 * Inv_LW;
                  pParam->BSIM4v7up = model->BSIM4v7up
                                  + model->BSIM4v7lup * Inv_L
                                  + model->BSIM4v7wup * Inv_W
                                  + model->BSIM4v7pup * Inv_LW;
                  pParam->BSIM4v7lp = model->BSIM4v7lp
                                  + model->BSIM4v7llp * Inv_L
                                  + model->BSIM4v7wlp * Inv_W
                                  + model->BSIM4v7plp * Inv_LW;
                  pParam->BSIM4v7eu = model->BSIM4v7eu
                                  + model->BSIM4v7leu * Inv_L
                                  + model->BSIM4v7weu * Inv_W
                                  + model->BSIM4v7peu * Inv_LW;
                  pParam->BSIM4v7u0 = model->BSIM4v7u0
                                  + model->BSIM4v7lu0 * Inv_L
                                  + model->BSIM4v7wu0 * Inv_W
                                  + model->BSIM4v7pu0 * Inv_LW;
                  pParam->BSIM4v7ute = model->BSIM4v7ute
                                   + model->BSIM4v7lute * Inv_L
                                   + model->BSIM4v7wute * Inv_W
                                   + model->BSIM4v7pute * Inv_LW;
                /*high k mobility*/
                 pParam->BSIM4v7ucs = model->BSIM4v7ucs
                                  + model->BSIM4v7lucs * Inv_L
                                  + model->BSIM4v7wucs * Inv_W
                                  + model->BSIM4v7pucs * Inv_LW;
                  pParam->BSIM4v7ucste = model->BSIM4v7ucste
                           + model->BSIM4v7lucste * Inv_L
                                   + model->BSIM4v7wucste * Inv_W
                                   + model->BSIM4v7pucste * Inv_LW;

                  pParam->BSIM4v7voff = model->BSIM4v7voff
                                    + model->BSIM4v7lvoff * Inv_L
                                    + model->BSIM4v7wvoff * Inv_W
                                    + model->BSIM4v7pvoff * Inv_LW;
                  pParam->BSIM4v7tvoff = model->BSIM4v7tvoff
                                    + model->BSIM4v7ltvoff * Inv_L
                                    + model->BSIM4v7wtvoff * Inv_W
                                    + model->BSIM4v7ptvoff * Inv_LW;
                  pParam->BSIM4v7minv = model->BSIM4v7minv
                                    + model->BSIM4v7lminv * Inv_L
                                    + model->BSIM4v7wminv * Inv_W
                                    + model->BSIM4v7pminv * Inv_LW;
                  pParam->BSIM4v7minvcv = model->BSIM4v7minvcv
                                    + model->BSIM4v7lminvcv * Inv_L
                                    + model->BSIM4v7wminvcv * Inv_W
                                    + model->BSIM4v7pminvcv * Inv_LW;
                  pParam->BSIM4v7fprout = model->BSIM4v7fprout
                                     + model->BSIM4v7lfprout * Inv_L
                                     + model->BSIM4v7wfprout * Inv_W
                                     + model->BSIM4v7pfprout * Inv_LW;
                  pParam->BSIM4v7pdits = model->BSIM4v7pdits
                                     + model->BSIM4v7lpdits * Inv_L
                                     + model->BSIM4v7wpdits * Inv_W
                                     + model->BSIM4v7ppdits * Inv_LW;
                  pParam->BSIM4v7pditsd = model->BSIM4v7pditsd
                                      + model->BSIM4v7lpditsd * Inv_L
                                      + model->BSIM4v7wpditsd * Inv_W
                                      + model->BSIM4v7ppditsd * Inv_LW;
                  pParam->BSIM4v7delta = model->BSIM4v7delta
                                     + model->BSIM4v7ldelta * Inv_L
                                     + model->BSIM4v7wdelta * Inv_W
                                     + model->BSIM4v7pdelta * Inv_LW;
                  pParam->BSIM4v7rdsw = model->BSIM4v7rdsw
                                    + model->BSIM4v7lrdsw * Inv_L
                                    + model->BSIM4v7wrdsw * Inv_W
                                    + model->BSIM4v7prdsw * Inv_LW;
                  pParam->BSIM4v7rdw = model->BSIM4v7rdw
                                    + model->BSIM4v7lrdw * Inv_L
                                    + model->BSIM4v7wrdw * Inv_W
                                    + model->BSIM4v7prdw * Inv_LW;
                  pParam->BSIM4v7rsw = model->BSIM4v7rsw
                                    + model->BSIM4v7lrsw * Inv_L
                                    + model->BSIM4v7wrsw * Inv_W
                                    + model->BSIM4v7prsw * Inv_LW;
                  pParam->BSIM4v7prwg = model->BSIM4v7prwg
                                    + model->BSIM4v7lprwg * Inv_L
                                    + model->BSIM4v7wprwg * Inv_W
                                    + model->BSIM4v7pprwg * Inv_LW;
                  pParam->BSIM4v7prwb = model->BSIM4v7prwb
                                    + model->BSIM4v7lprwb * Inv_L
                                    + model->BSIM4v7wprwb * Inv_W
                                    + model->BSIM4v7pprwb * Inv_LW;
                  pParam->BSIM4v7prt = model->BSIM4v7prt
                                    + model->BSIM4v7lprt * Inv_L
                                    + model->BSIM4v7wprt * Inv_W
                                    + model->BSIM4v7pprt * Inv_LW;
                  pParam->BSIM4v7eta0 = model->BSIM4v7eta0
                                    + model->BSIM4v7leta0 * Inv_L
                                    + model->BSIM4v7weta0 * Inv_W
                                    + model->BSIM4v7peta0 * Inv_LW;
                  pParam->BSIM4v7teta0 = model->BSIM4v7teta0                 /* v4.7  */
                                    + model->BSIM4v7lteta0 * Inv_L
                                    + model->BSIM4v7wteta0 * Inv_W
                                    + model->BSIM4v7pteta0 * Inv_LW;
                  pParam->BSIM4v7etab = model->BSIM4v7etab
                                    + model->BSIM4v7letab * Inv_L
                                    + model->BSIM4v7wetab * Inv_W
                                    + model->BSIM4v7petab * Inv_LW;
                  pParam->BSIM4v7pclm = model->BSIM4v7pclm
                                    + model->BSIM4v7lpclm * Inv_L
                                    + model->BSIM4v7wpclm * Inv_W
                                    + model->BSIM4v7ppclm * Inv_LW;
                  pParam->BSIM4v7pdibl1 = model->BSIM4v7pdibl1
                                      + model->BSIM4v7lpdibl1 * Inv_L
                                      + model->BSIM4v7wpdibl1 * Inv_W
                                      + model->BSIM4v7ppdibl1 * Inv_LW;
                  pParam->BSIM4v7pdibl2 = model->BSIM4v7pdibl2
                                      + model->BSIM4v7lpdibl2 * Inv_L
                                      + model->BSIM4v7wpdibl2 * Inv_W
                                      + model->BSIM4v7ppdibl2 * Inv_LW;
                  pParam->BSIM4v7pdiblb = model->BSIM4v7pdiblb
                                      + model->BSIM4v7lpdiblb * Inv_L
                                      + model->BSIM4v7wpdiblb * Inv_W
                                      + model->BSIM4v7ppdiblb * Inv_LW;
                  pParam->BSIM4v7pscbe1 = model->BSIM4v7pscbe1
                                      + model->BSIM4v7lpscbe1 * Inv_L
                                      + model->BSIM4v7wpscbe1 * Inv_W
                                      + model->BSIM4v7ppscbe1 * Inv_LW;
                  pParam->BSIM4v7pscbe2 = model->BSIM4v7pscbe2
                                      + model->BSIM4v7lpscbe2 * Inv_L
                                      + model->BSIM4v7wpscbe2 * Inv_W
                                      + model->BSIM4v7ppscbe2 * Inv_LW;
                  pParam->BSIM4v7pvag = model->BSIM4v7pvag
                                    + model->BSIM4v7lpvag * Inv_L
                                    + model->BSIM4v7wpvag * Inv_W
                                    + model->BSIM4v7ppvag * Inv_LW;
                  pParam->BSIM4v7wr = model->BSIM4v7wr
                                  + model->BSIM4v7lwr * Inv_L
                                  + model->BSIM4v7wwr * Inv_W
                                  + model->BSIM4v7pwr * Inv_LW;
                  pParam->BSIM4v7dwg = model->BSIM4v7dwg
                                   + model->BSIM4v7ldwg * Inv_L
                                   + model->BSIM4v7wdwg * Inv_W
                                   + model->BSIM4v7pdwg * Inv_LW;
                  pParam->BSIM4v7dwb = model->BSIM4v7dwb
                                   + model->BSIM4v7ldwb * Inv_L
                                   + model->BSIM4v7wdwb * Inv_W
                                   + model->BSIM4v7pdwb * Inv_LW;
                  pParam->BSIM4v7b0 = model->BSIM4v7b0
                                  + model->BSIM4v7lb0 * Inv_L
                                  + model->BSIM4v7wb0 * Inv_W
                                  + model->BSIM4v7pb0 * Inv_LW;
                  pParam->BSIM4v7b1 = model->BSIM4v7b1
                                  + model->BSIM4v7lb1 * Inv_L
                                  + model->BSIM4v7wb1 * Inv_W
                                  + model->BSIM4v7pb1 * Inv_LW;
                  pParam->BSIM4v7alpha0 = model->BSIM4v7alpha0
                                      + model->BSIM4v7lalpha0 * Inv_L
                                      + model->BSIM4v7walpha0 * Inv_W
                                      + model->BSIM4v7palpha0 * Inv_LW;
                  pParam->BSIM4v7alpha1 = model->BSIM4v7alpha1
                                      + model->BSIM4v7lalpha1 * Inv_L
                                      + model->BSIM4v7walpha1 * Inv_W
                                      + model->BSIM4v7palpha1 * Inv_LW;
                  pParam->BSIM4v7beta0 = model->BSIM4v7beta0
                                     + model->BSIM4v7lbeta0 * Inv_L
                                     + model->BSIM4v7wbeta0 * Inv_W
                                     + model->BSIM4v7pbeta0 * Inv_LW;
                  pParam->BSIM4v7agidl = model->BSIM4v7agidl
                                     + model->BSIM4v7lagidl * Inv_L
                                     + model->BSIM4v7wagidl * Inv_W
                                     + model->BSIM4v7pagidl * Inv_LW;
                  pParam->BSIM4v7bgidl = model->BSIM4v7bgidl
                                     + model->BSIM4v7lbgidl * Inv_L
                                     + model->BSIM4v7wbgidl * Inv_W
                                     + model->BSIM4v7pbgidl * Inv_LW;
                  pParam->BSIM4v7cgidl = model->BSIM4v7cgidl
                                     + model->BSIM4v7lcgidl * Inv_L
                                     + model->BSIM4v7wcgidl * Inv_W
                                     + model->BSIM4v7pcgidl * Inv_LW;
                  pParam->BSIM4v7egidl = model->BSIM4v7egidl
                                     + model->BSIM4v7legidl * Inv_L
                                     + model->BSIM4v7wegidl * Inv_W
                                     + model->BSIM4v7pegidl * Inv_LW;
                  pParam->BSIM4v7rgidl = model->BSIM4v7rgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lrgidl * Inv_L
                                     + model->BSIM4v7wrgidl * Inv_W
                                     + model->BSIM4v7prgidl * Inv_LW;
                  pParam->BSIM4v7kgidl = model->BSIM4v7kgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lkgidl * Inv_L
                                     + model->BSIM4v7wkgidl * Inv_W
                                     + model->BSIM4v7pkgidl * Inv_LW;
                  pParam->BSIM4v7fgidl = model->BSIM4v7fgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lfgidl * Inv_L
                                     + model->BSIM4v7wfgidl * Inv_W
                                     + model->BSIM4v7pfgidl * Inv_LW;
                  pParam->BSIM4v7agisl = model->BSIM4v7agisl
                                     + model->BSIM4v7lagisl * Inv_L
                                     + model->BSIM4v7wagisl * Inv_W
                                     + model->BSIM4v7pagisl * Inv_LW;
                  pParam->BSIM4v7bgisl = model->BSIM4v7bgisl
                                     + model->BSIM4v7lbgisl * Inv_L
                                     + model->BSIM4v7wbgisl * Inv_W
                                     + model->BSIM4v7pbgisl * Inv_LW;
                  pParam->BSIM4v7cgisl = model->BSIM4v7cgisl
                                     + model->BSIM4v7lcgisl * Inv_L
                                     + model->BSIM4v7wcgisl * Inv_W
                                     + model->BSIM4v7pcgisl * Inv_LW;
                  pParam->BSIM4v7egisl = model->BSIM4v7egisl
                                     + model->BSIM4v7legisl * Inv_L
                                     + model->BSIM4v7wegisl * Inv_W
                                     + model->BSIM4v7pegisl * Inv_LW;
                   pParam->BSIM4v7rgisl = model->BSIM4v7rgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lrgisl * Inv_L
                                     + model->BSIM4v7wrgisl * Inv_W
                                     + model->BSIM4v7prgisl * Inv_LW;
                  pParam->BSIM4v7kgisl = model->BSIM4v7kgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lkgisl * Inv_L
                                     + model->BSIM4v7wkgisl * Inv_W
                                     + model->BSIM4v7pkgisl * Inv_LW;
                  pParam->BSIM4v7fgisl = model->BSIM4v7fgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4v7lfgisl * Inv_L
                                     + model->BSIM4v7wfgisl * Inv_W
                                     + model->BSIM4v7pfgisl * Inv_LW;
                  pParam->BSIM4v7aigc = model->BSIM4v7aigc
                                     + model->BSIM4v7laigc * Inv_L
                                     + model->BSIM4v7waigc * Inv_W
                                     + model->BSIM4v7paigc * Inv_LW;
                  pParam->BSIM4v7bigc = model->BSIM4v7bigc
                                     + model->BSIM4v7lbigc * Inv_L
                                     + model->BSIM4v7wbigc * Inv_W
                                     + model->BSIM4v7pbigc * Inv_LW;
                  pParam->BSIM4v7cigc = model->BSIM4v7cigc
                                     + model->BSIM4v7lcigc * Inv_L
                                     + model->BSIM4v7wcigc * Inv_W
                                     + model->BSIM4v7pcigc * Inv_LW;
                  pParam->BSIM4v7aigs = model->BSIM4v7aigs
                                     + model->BSIM4v7laigs * Inv_L
                                     + model->BSIM4v7waigs * Inv_W
                                     + model->BSIM4v7paigs * Inv_LW;
                  pParam->BSIM4v7bigs = model->BSIM4v7bigs
                                     + model->BSIM4v7lbigs * Inv_L
                                     + model->BSIM4v7wbigs * Inv_W
                                     + model->BSIM4v7pbigs * Inv_LW;
                  pParam->BSIM4v7cigs = model->BSIM4v7cigs
                                     + model->BSIM4v7lcigs * Inv_L
                                     + model->BSIM4v7wcigs * Inv_W
                                     + model->BSIM4v7pcigs * Inv_LW;
                  pParam->BSIM4v7aigd = model->BSIM4v7aigd
                                     + model->BSIM4v7laigd * Inv_L
                                     + model->BSIM4v7waigd * Inv_W
                                     + model->BSIM4v7paigd * Inv_LW;
                  pParam->BSIM4v7bigd = model->BSIM4v7bigd
                                     + model->BSIM4v7lbigd * Inv_L
                                     + model->BSIM4v7wbigd * Inv_W
                                     + model->BSIM4v7pbigd * Inv_LW;
                  pParam->BSIM4v7cigd = model->BSIM4v7cigd
                                     + model->BSIM4v7lcigd * Inv_L
                                     + model->BSIM4v7wcigd * Inv_W
                                     + model->BSIM4v7pcigd * Inv_LW;
                  pParam->BSIM4v7aigbacc = model->BSIM4v7aigbacc
                                       + model->BSIM4v7laigbacc * Inv_L
                                       + model->BSIM4v7waigbacc * Inv_W
                                       + model->BSIM4v7paigbacc * Inv_LW;
                  pParam->BSIM4v7bigbacc = model->BSIM4v7bigbacc
                                       + model->BSIM4v7lbigbacc * Inv_L
                                       + model->BSIM4v7wbigbacc * Inv_W
                                       + model->BSIM4v7pbigbacc * Inv_LW;
                  pParam->BSIM4v7cigbacc = model->BSIM4v7cigbacc
                                       + model->BSIM4v7lcigbacc * Inv_L
                                       + model->BSIM4v7wcigbacc * Inv_W
                                       + model->BSIM4v7pcigbacc * Inv_LW;
                  pParam->BSIM4v7aigbinv = model->BSIM4v7aigbinv
                                       + model->BSIM4v7laigbinv * Inv_L
                                       + model->BSIM4v7waigbinv * Inv_W
                                       + model->BSIM4v7paigbinv * Inv_LW;
                  pParam->BSIM4v7bigbinv = model->BSIM4v7bigbinv
                                       + model->BSIM4v7lbigbinv * Inv_L
                                       + model->BSIM4v7wbigbinv * Inv_W
                                       + model->BSIM4v7pbigbinv * Inv_LW;
                  pParam->BSIM4v7cigbinv = model->BSIM4v7cigbinv
                                       + model->BSIM4v7lcigbinv * Inv_L
                                       + model->BSIM4v7wcigbinv * Inv_W
                                       + model->BSIM4v7pcigbinv * Inv_LW;
                  pParam->BSIM4v7nigc = model->BSIM4v7nigc
                                       + model->BSIM4v7lnigc * Inv_L
                                       + model->BSIM4v7wnigc * Inv_W
                                       + model->BSIM4v7pnigc * Inv_LW;
                  pParam->BSIM4v7nigbacc = model->BSIM4v7nigbacc
                                       + model->BSIM4v7lnigbacc * Inv_L
                                       + model->BSIM4v7wnigbacc * Inv_W
                                       + model->BSIM4v7pnigbacc * Inv_LW;
                  pParam->BSIM4v7nigbinv = model->BSIM4v7nigbinv
                                       + model->BSIM4v7lnigbinv * Inv_L
                                       + model->BSIM4v7wnigbinv * Inv_W
                                       + model->BSIM4v7pnigbinv * Inv_LW;
                  pParam->BSIM4v7ntox = model->BSIM4v7ntox
                                    + model->BSIM4v7lntox * Inv_L
                                    + model->BSIM4v7wntox * Inv_W
                                    + model->BSIM4v7pntox * Inv_LW;
                  pParam->BSIM4v7eigbinv = model->BSIM4v7eigbinv
                                       + model->BSIM4v7leigbinv * Inv_L
                                       + model->BSIM4v7weigbinv * Inv_W
                                       + model->BSIM4v7peigbinv * Inv_LW;
                  pParam->BSIM4v7pigcd = model->BSIM4v7pigcd
                                     + model->BSIM4v7lpigcd * Inv_L
                                     + model->BSIM4v7wpigcd * Inv_W
                                     + model->BSIM4v7ppigcd * Inv_LW;
                  pParam->BSIM4v7poxedge = model->BSIM4v7poxedge
                                       + model->BSIM4v7lpoxedge * Inv_L
                                       + model->BSIM4v7wpoxedge * Inv_W
                                       + model->BSIM4v7ppoxedge * Inv_LW;
                  pParam->BSIM4v7xrcrg1 = model->BSIM4v7xrcrg1
                                      + model->BSIM4v7lxrcrg1 * Inv_L
                                      + model->BSIM4v7wxrcrg1 * Inv_W
                                      + model->BSIM4v7pxrcrg1 * Inv_LW;
                  pParam->BSIM4v7xrcrg2 = model->BSIM4v7xrcrg2
                                      + model->BSIM4v7lxrcrg2 * Inv_L
                                      + model->BSIM4v7wxrcrg2 * Inv_W
                                      + model->BSIM4v7pxrcrg2 * Inv_LW;
                  pParam->BSIM4v7lambda = model->BSIM4v7lambda
                                      + model->BSIM4v7llambda * Inv_L
                                      + model->BSIM4v7wlambda * Inv_W
                                      + model->BSIM4v7plambda * Inv_LW;
                  pParam->BSIM4v7vtl = model->BSIM4v7vtl
                                      + model->BSIM4v7lvtl * Inv_L
                                      + model->BSIM4v7wvtl * Inv_W
                                      + model->BSIM4v7pvtl * Inv_LW;
                  pParam->BSIM4v7xn = model->BSIM4v7xn
                                      + model->BSIM4v7lxn * Inv_L
                                      + model->BSIM4v7wxn * Inv_W
                                      + model->BSIM4v7pxn * Inv_LW;
                  pParam->BSIM4v7vfbsdoff = model->BSIM4v7vfbsdoff
                                      + model->BSIM4v7lvfbsdoff * Inv_L
                                      + model->BSIM4v7wvfbsdoff * Inv_W
                                      + model->BSIM4v7pvfbsdoff * Inv_LW;
                  pParam->BSIM4v7tvfbsdoff = model->BSIM4v7tvfbsdoff
                                      + model->BSIM4v7ltvfbsdoff * Inv_L
                                      + model->BSIM4v7wtvfbsdoff * Inv_W
                                      + model->BSIM4v7ptvfbsdoff * Inv_LW;

                  pParam->BSIM4v7cgsl = model->BSIM4v7cgsl
                                    + model->BSIM4v7lcgsl * Inv_L
                                    + model->BSIM4v7wcgsl * Inv_W
                                    + model->BSIM4v7pcgsl * Inv_LW;
                  pParam->BSIM4v7cgdl = model->BSIM4v7cgdl
                                    + model->BSIM4v7lcgdl * Inv_L
                                    + model->BSIM4v7wcgdl * Inv_W
                                    + model->BSIM4v7pcgdl * Inv_LW;
                  pParam->BSIM4v7ckappas = model->BSIM4v7ckappas
                                       + model->BSIM4v7lckappas * Inv_L
                                       + model->BSIM4v7wckappas * Inv_W
                                        + model->BSIM4v7pckappas * Inv_LW;
                  pParam->BSIM4v7ckappad = model->BSIM4v7ckappad
                                       + model->BSIM4v7lckappad * Inv_L
                                       + model->BSIM4v7wckappad * Inv_W
                                       + model->BSIM4v7pckappad * Inv_LW;
                  pParam->BSIM4v7cf = model->BSIM4v7cf
                                  + model->BSIM4v7lcf * Inv_L
                                  + model->BSIM4v7wcf * Inv_W
                                  + model->BSIM4v7pcf * Inv_LW;
                  pParam->BSIM4v7clc = model->BSIM4v7clc
                                   + model->BSIM4v7lclc * Inv_L
                                   + model->BSIM4v7wclc * Inv_W
                                   + model->BSIM4v7pclc * Inv_LW;
                  pParam->BSIM4v7cle = model->BSIM4v7cle
                                   + model->BSIM4v7lcle * Inv_L
                                   + model->BSIM4v7wcle * Inv_W
                                   + model->BSIM4v7pcle * Inv_LW;
                  pParam->BSIM4v7vfbcv = model->BSIM4v7vfbcv
                                     + model->BSIM4v7lvfbcv * Inv_L
                                     + model->BSIM4v7wvfbcv * Inv_W
                                     + model->BSIM4v7pvfbcv * Inv_LW;
                  pParam->BSIM4v7acde = model->BSIM4v7acde
                                    + model->BSIM4v7lacde * Inv_L
                                    + model->BSIM4v7wacde * Inv_W
                                    + model->BSIM4v7pacde * Inv_LW;
                  pParam->BSIM4v7moin = model->BSIM4v7moin
                                    + model->BSIM4v7lmoin * Inv_L
                                    + model->BSIM4v7wmoin * Inv_W
                                    + model->BSIM4v7pmoin * Inv_LW;
                  pParam->BSIM4v7noff = model->BSIM4v7noff
                                    + model->BSIM4v7lnoff * Inv_L
                                    + model->BSIM4v7wnoff * Inv_W
                                    + model->BSIM4v7pnoff * Inv_LW;
                  pParam->BSIM4v7voffcv = model->BSIM4v7voffcv
                                      + model->BSIM4v7lvoffcv * Inv_L
                                      + model->BSIM4v7wvoffcv * Inv_W
                                      + model->BSIM4v7pvoffcv * Inv_LW;
                  pParam->BSIM4v7kvth0we = model->BSIM4v7kvth0we
                                      + model->BSIM4v7lkvth0we * Inv_L
                                      + model->BSIM4v7wkvth0we * Inv_W
                                      + model->BSIM4v7pkvth0we * Inv_LW;
                  pParam->BSIM4v7k2we = model->BSIM4v7k2we
                                      + model->BSIM4v7lk2we * Inv_L
                                      + model->BSIM4v7wk2we * Inv_W
                                      + model->BSIM4v7pk2we * Inv_LW;
                  pParam->BSIM4v7ku0we = model->BSIM4v7ku0we
                                      + model->BSIM4v7lku0we * Inv_L
                                      + model->BSIM4v7wku0we * Inv_W
                                      + model->BSIM4v7pku0we * Inv_LW;

                  pParam->BSIM4v7abulkCVfactor = 1.0 + pow((pParam->BSIM4v7clc
                                             / pParam->BSIM4v7leffCV),
                                             pParam->BSIM4v7cle);

                  T0 = (TRatio - 1.0);

                  PowWeffWr = pow(pParam->BSIM4v7weffCJ * 1.0e6, pParam->BSIM4v7wr) * here->BSIM4v7nf;

                  T1 = T2 = T3 = T4 = 0.0;
                  pParam->BSIM4v7ucs = pParam->BSIM4v7ucs * pow(TRatio, pParam->BSIM4v7ucste);
                  if(model->BSIM4v7tempMod == 0)
                  {
                      pParam->BSIM4v7ua = pParam->BSIM4v7ua + pParam->BSIM4v7ua1 * T0;
                      pParam->BSIM4v7ub = pParam->BSIM4v7ub + pParam->BSIM4v7ub1 * T0;
                      pParam->BSIM4v7uc = pParam->BSIM4v7uc + pParam->BSIM4v7uc1 * T0;
                      pParam->BSIM4v7ud = pParam->BSIM4v7ud + pParam->BSIM4v7ud1 * T0;
                      pParam->BSIM4v7vsattemp = pParam->BSIM4v7vsat - pParam->BSIM4v7at * T0;
                      T10 = pParam->BSIM4v7prt * T0;
                      if(model->BSIM4v7rdsMod)
                      {
                            /* External Rd(V) */
                            T1 = pParam->BSIM4v7rdw + T10;
                            T2 = model->BSIM4v7rdwmin + T10;
                            /* External Rs(V) */
                            T3 = pParam->BSIM4v7rsw + T10;
                            T4 = model->BSIM4v7rswmin + T10;
                      }
                      /* Internal Rds(V) in IV */
                      pParam->BSIM4v7rds0 = (pParam->BSIM4v7rdsw + T10)
                                            * here->BSIM4v7nf / PowWeffWr;
                      pParam->BSIM4v7rdswmin = (model->BSIM4v7rdswmin + T10)
                                               * here->BSIM4v7nf / PowWeffWr;
                  }
                  else
                  {
                      if (model->BSIM4v7tempMod == 3)
                      {
                          pParam->BSIM4v7ua = pParam->BSIM4v7ua * pow(TRatio, pParam->BSIM4v7ua1) ;
                            pParam->BSIM4v7ub = pParam->BSIM4v7ub * pow(TRatio, pParam->BSIM4v7ub1);
                            pParam->BSIM4v7uc = pParam->BSIM4v7uc * pow(TRatio, pParam->BSIM4v7uc1);
                            pParam->BSIM4v7ud = pParam->BSIM4v7ud * pow(TRatio, pParam->BSIM4v7ud1);
                      }
                      else
                      {
                          /* tempMod = 1, 2 */
                            pParam->BSIM4v7ua = pParam->BSIM4v7ua * (1.0 + pParam->BSIM4v7ua1 * delTemp) ;
                            pParam->BSIM4v7ub = pParam->BSIM4v7ub * (1.0 + pParam->BSIM4v7ub1 * delTemp);
                            pParam->BSIM4v7uc = pParam->BSIM4v7uc * (1.0 + pParam->BSIM4v7uc1 * delTemp);
                            pParam->BSIM4v7ud = pParam->BSIM4v7ud * (1.0 + pParam->BSIM4v7ud1 * delTemp);
                      }
                      pParam->BSIM4v7vsattemp = pParam->BSIM4v7vsat * (1.0 - pParam->BSIM4v7at * delTemp);
                      T10 = 1.0 + pParam->BSIM4v7prt * delTemp;
                      if(model->BSIM4v7rdsMod)
                      {
                            /* External Rd(V) */
                            T1 = pParam->BSIM4v7rdw * T10;
                            T2 = model->BSIM4v7rdwmin * T10;

                            /* External Rs(V) */
                            T3 = pParam->BSIM4v7rsw * T10;
                            T4 = model->BSIM4v7rswmin * T10;
                      }
                      /* Internal Rds(V) in IV */
                      pParam->BSIM4v7rds0 = pParam->BSIM4v7rdsw * T10 * here->BSIM4v7nf / PowWeffWr;
                      pParam->BSIM4v7rdswmin = model->BSIM4v7rdswmin * T10 * here->BSIM4v7nf / PowWeffWr;
                  }

                  if (T1 < 0.0)
                  {   T1 = 0.0;
                      printf("Warning: Rdw at current temperature is negative; set to 0.\n");
                  }
                  if (T2 < 0.0)
                  {   T2 = 0.0;
                      printf("Warning: Rdwmin at current temperature is negative; set to 0.\n");
                  }
                  pParam->BSIM4v7rd0 = T1 / PowWeffWr;
                  pParam->BSIM4v7rdwmin = T2 / PowWeffWr;
                  if (T3 < 0.0)
                  {   T3 = 0.0;
                      printf("Warning: Rsw at current temperature is negative; set to 0.\n");
                  }
                  if (T4 < 0.0)
                  {   T4 = 0.0;
                      printf("Warning: Rswmin at current temperature is negative; set to 0.\n");
                  }
                  pParam->BSIM4v7rs0 = T3 / PowWeffWr;
                  pParam->BSIM4v7rswmin = T4 / PowWeffWr;

                  if (pParam->BSIM4v7u0 > 1.0)
                      pParam->BSIM4v7u0 = pParam->BSIM4v7u0 / 1.0e4;

                  /* mobility channel length dependence */
                  T5 = 1.0 - pParam->BSIM4v7up * exp( - pParam->BSIM4v7leff / pParam->BSIM4v7lp);
                  pParam->BSIM4v7u0temp = pParam->BSIM4v7u0 * T5
                                      * pow(TRatio, pParam->BSIM4v7ute);
                  if (pParam->BSIM4v7eu < 0.0)
                  {   pParam->BSIM4v7eu = 0.0;
                      printf("Warning: eu has been negative; reset to 0.0.\n");
                  }
                  if (pParam->BSIM4v7ucs < 0.0)
                  {   pParam->BSIM4v7ucs = 0.0;
                      printf("Warning: ucs has been negative; reset to 0.0.\n");
                  }

                 pParam->BSIM4v7vfbsdoff = pParam->BSIM4v7vfbsdoff * (1.0 + pParam->BSIM4v7tvfbsdoff * delTemp);
                 pParam->BSIM4v7voff = pParam->BSIM4v7voff * (1.0 + pParam->BSIM4v7tvoff * delTemp);

                 pParam->BSIM4v7nfactor = pParam->BSIM4v7nfactor + pParam->BSIM4v7tnfactor * delTemp / Tnom;  /* v4.7 temp dep of leakage currents */
                        pParam->BSIM4v7voffcv = pParam->BSIM4v7voffcv * (1.0 + pParam->BSIM4v7tvoffcv * delTemp);   /*         v4.7 temp dep of leakage currents */
                 pParam->BSIM4v7eta0 = pParam->BSIM4v7eta0 + pParam->BSIM4v7teta0 * delTemp / Tnom;   /*         v4.7 temp dep of leakage currents */

                /* Source End Velocity Limit  */
                        if((model->BSIM4v7vtlGiven) && (model->BSIM4v7vtl > 0.0) )
                      {
                     if(model->BSIM4v7lc < 0.0) pParam->BSIM4v7lc = 0.0;
                     else   pParam->BSIM4v7lc = model->BSIM4v7lc ;
                     T0 = pParam->BSIM4v7leff / (pParam->BSIM4v7xn * pParam->BSIM4v7leff + pParam->BSIM4v7lc);
                     pParam->BSIM4v7tfactor = (1.0 - T0) / (1.0 + T0 );
                       }

                  pParam->BSIM4v7cgdo = (model->BSIM4v7cgdo + pParam->BSIM4v7cf)
                                    * pParam->BSIM4v7weffCV;
                  pParam->BSIM4v7cgso = (model->BSIM4v7cgso + pParam->BSIM4v7cf)
                                    * pParam->BSIM4v7weffCV;
                  pParam->BSIM4v7cgbo = model->BSIM4v7cgbo * pParam->BSIM4v7leffCV * here->BSIM4v7nf;

                  if (!model->BSIM4v7ndepGiven && model->BSIM4v7gamma1Given)
                  {   T0 = pParam->BSIM4v7gamma1 * model->BSIM4v7coxe;
                      pParam->BSIM4v7ndep = 3.01248e22 * T0 * T0;
                  }

                  pParam->BSIM4v7phi = Vtm0 * log(pParam->BSIM4v7ndep / ni)
                                   + pParam->BSIM4v7phin + 0.4;

                  pParam->BSIM4v7sqrtPhi = sqrt(pParam->BSIM4v7phi);
                  pParam->BSIM4v7phis3 = pParam->BSIM4v7sqrtPhi * pParam->BSIM4v7phi;

                  pParam->BSIM4v7Xdep0 = sqrt(2.0 * epssub / (Charge_q
                                     * pParam->BSIM4v7ndep * 1.0e6))
                                     * pParam->BSIM4v7sqrtPhi;
                  pParam->BSIM4v7sqrtXdep0 = sqrt(pParam->BSIM4v7Xdep0);

                  if(model->BSIM4v7mtrlMod == 0)
                    pParam->BSIM4v7litl = sqrt(3.0 * 3.9 / epsrox * pParam->BSIM4v7xj * toxe);
                  else
                    pParam->BSIM4v7litl = sqrt(model->BSIM4v7epsrsub/epsrox * pParam->BSIM4v7xj * toxe);

                  pParam->BSIM4v7vbi = Vtm0 * log(pParam->BSIM4v7nsd
                                   * pParam->BSIM4v7ndep / (ni * ni));

                  if (model->BSIM4v7mtrlMod == 0)
                  {
                    if (pParam->BSIM4v7ngate > 0.0)
                    {   pParam->BSIM4v7vfbsd = Vtm0 * log(pParam->BSIM4v7ngate
                                         / pParam->BSIM4v7nsd);
                     }
                    else
                      pParam->BSIM4v7vfbsd = 0.0;
                  }
                  else
                  {
                    T0 = Vtm0 * log(pParam->BSIM4v7nsd/ni);
                    T1 = 0.5 * Eg0;
                    if(T0 > T1)
                      T0 = T1;
                    T2 = model->BSIM4v7easub + T1 - model->BSIM4v7type * T0;
                    pParam->BSIM4v7vfbsd = model->BSIM4v7phig - T2;
                  }

                  pParam->BSIM4v7cdep0 = sqrt(Charge_q * epssub
                                     * pParam->BSIM4v7ndep * 1.0e6 / 2.0
                                     / pParam->BSIM4v7phi);

                  pParam->BSIM4v7ToxRatio = exp(pParam->BSIM4v7ntox
                                        * log(model->BSIM4v7toxref / toxe))
                                        / toxe / toxe;
                  pParam->BSIM4v7ToxRatioEdge = exp(pParam->BSIM4v7ntox
                                            * log(model->BSIM4v7toxref
                                            / (toxe * pParam->BSIM4v7poxedge)))
                                            / toxe / toxe
                                            / pParam->BSIM4v7poxedge / pParam->BSIM4v7poxedge;
                  pParam->BSIM4v7Aechvb = (model->BSIM4v7type == NMOS) ? 4.97232e-7 : 3.42537e-7;
                  pParam->BSIM4v7Bechvb = (model->BSIM4v7type == NMOS) ? 7.45669e11 : 1.16645e12;
                  pParam->BSIM4v7AechvbEdgeS = pParam->BSIM4v7Aechvb * pParam->BSIM4v7weff
                                          * model->BSIM4v7dlcig * pParam->BSIM4v7ToxRatioEdge;
                  pParam->BSIM4v7AechvbEdgeD = pParam->BSIM4v7Aechvb * pParam->BSIM4v7weff
                                          * model->BSIM4v7dlcigd * pParam->BSIM4v7ToxRatioEdge;
                  pParam->BSIM4v7BechvbEdge = -pParam->BSIM4v7Bechvb
                                          * toxe * pParam->BSIM4v7poxedge;
                  pParam->BSIM4v7Aechvb *= pParam->BSIM4v7weff * pParam->BSIM4v7leff
                                       * pParam->BSIM4v7ToxRatio;
                  pParam->BSIM4v7Bechvb *= -toxe;


                  pParam->BSIM4v7mstar = 0.5 + atan(pParam->BSIM4v7minv) / PI;
                  pParam->BSIM4v7mstarcv = 0.5 + atan(pParam->BSIM4v7minvcv) / PI;
                  pParam->BSIM4v7voffcbn =  pParam->BSIM4v7voff + model->BSIM4v7voffl / pParam->BSIM4v7leff;
                  pParam->BSIM4v7voffcbncv =  pParam->BSIM4v7voffcv + model->BSIM4v7voffcvl / pParam->BSIM4v7leff;

                  pParam->BSIM4v7ldeb = sqrt(epssub * Vtm0 / (Charge_q
                                    * pParam->BSIM4v7ndep * 1.0e6)) / 3.0;
                  pParam->BSIM4v7acde *= pow((pParam->BSIM4v7ndep / 2.0e16), -0.25);


                  if (model->BSIM4v7k1Given || model->BSIM4v7k2Given)
                  {   if (!model->BSIM4v7k1Given)
                      {
                          if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) /* don't print in sensitivity */
                              fprintf(stdout, "Warning: k1 should be specified with k2.\n");
                          pParam->BSIM4v7k1 = 0.53;
                      }
                      if (!model->BSIM4v7k2Given)
                      {
                          if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) /* don't print in sensitivity */
                              fprintf(stdout, "Warning: k2 should be specified with k1.\n");
                          pParam->BSIM4v7k2 = -0.0186;
                      }
                      if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) { /* don't print in sensitivity */
                          if (model->BSIM4v7nsubGiven)
                              fprintf(stdout, "Warning: nsub is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4v7xtGiven)
                              fprintf(stdout, "Warning: xt is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4v7vbxGiven)
                              fprintf(stdout, "Warning: vbx is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4v7gamma1Given)
                              fprintf(stdout, "Warning: gamma1 is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4v7gamma2Given)
                              fprintf(stdout, "Warning: gamma2 is ignored because k1 or k2 is given.\n");
                      }
                  }
                  else
                  {   if (!model->BSIM4v7vbxGiven)
                          pParam->BSIM4v7vbx = pParam->BSIM4v7phi - 7.7348e-4
                                           * pParam->BSIM4v7ndep
                                           * pParam->BSIM4v7xt * pParam->BSIM4v7xt;
                      if (pParam->BSIM4v7vbx > 0.0)
                          pParam->BSIM4v7vbx = -pParam->BSIM4v7vbx;
                      if (pParam->BSIM4v7vbm > 0.0)
                          pParam->BSIM4v7vbm = -pParam->BSIM4v7vbm;

                      if (!model->BSIM4v7gamma1Given)
                          pParam->BSIM4v7gamma1 = 5.753e-12
                                              * sqrt(pParam->BSIM4v7ndep)
                                              / model->BSIM4v7coxe;
                      if (!model->BSIM4v7gamma2Given)
                          pParam->BSIM4v7gamma2 = 5.753e-12
                                              * sqrt(pParam->BSIM4v7nsub)
                                              / model->BSIM4v7coxe;

                      T0 = pParam->BSIM4v7gamma1 - pParam->BSIM4v7gamma2;
                      T1 = sqrt(pParam->BSIM4v7phi - pParam->BSIM4v7vbx)
                         - pParam->BSIM4v7sqrtPhi;
                      T2 = sqrt(pParam->BSIM4v7phi * (pParam->BSIM4v7phi
                         - pParam->BSIM4v7vbm)) - pParam->BSIM4v7phi;
                      pParam->BSIM4v7k2 = T0 * T1 / (2.0 * T2 + pParam->BSIM4v7vbm);
                      pParam->BSIM4v7k1 = pParam->BSIM4v7gamma2 - 2.0
                                      * pParam->BSIM4v7k2 * sqrt(pParam->BSIM4v7phi
                                      - pParam->BSIM4v7vbm);
                  }

                  if (!model->BSIM4v7vfbGiven)
                  {
                    if (model->BSIM4v7vth0Given)
                      {   pParam->BSIM4v7vfb = model->BSIM4v7type * pParam->BSIM4v7vth0
                                           - pParam->BSIM4v7phi - pParam->BSIM4v7k1
                                           * pParam->BSIM4v7sqrtPhi;
                      }
                      else
                      {
                        if ((model->BSIM4v7mtrlMod) && (model->BSIM4v7phigGiven) &&
                            (model->BSIM4v7nsubGiven))
                          {
                            T0 = Vtm0 * log(pParam->BSIM4v7nsub/ni);
                            T1 = 0.5 * Eg0;
                            if(T0 > T1)
                              T0 = T1;
                            T2 = model->BSIM4v7easub + T1 + model->BSIM4v7type * T0;
                            pParam->BSIM4v7vfb = model->BSIM4v7phig - T2;
                          }
                        else
                          {
                            pParam->BSIM4v7vfb = -1.0;
                          }
                      }
                  }
                   if (!model->BSIM4v7vth0Given)
                  {   pParam->BSIM4v7vth0 = model->BSIM4v7type * (pParam->BSIM4v7vfb
                                        + pParam->BSIM4v7phi + pParam->BSIM4v7k1
                                        * pParam->BSIM4v7sqrtPhi);
                  }

                  pParam->BSIM4v7k1ox = pParam->BSIM4v7k1 * toxe
                                    / model->BSIM4v7toxm;

                  tmp = sqrt(epssub / (epsrox * EPS0) * toxe * pParam->BSIM4v7Xdep0);
                    T0 = pParam->BSIM4v7dsub * pParam->BSIM4v7leff / tmp;
                  if (T0 < EXP_THRESHOLD)
                    {   T1 = exp(T0);
                            T2 = T1 - 1.0;
                            T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      pParam->BSIM4v7theta0vb0 = T1 / T4;
                  }
                  else
                      pParam->BSIM4v7theta0vb0 = 1.0 / (MAX_EXP - 2.0);

                   T0 = pParam->BSIM4v7drout * pParam->BSIM4v7leff / tmp;
                  if (T0 < EXP_THRESHOLD)
                         {   T1 = exp(T0);
                            T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T5 = T1 / T4;
                  }
                  else
                      T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                  pParam->BSIM4v7thetaRout = pParam->BSIM4v7pdibl1 * T5
                                         + pParam->BSIM4v7pdibl2;

                  tmp = sqrt(pParam->BSIM4v7Xdep0);
                  tmp1 = pParam->BSIM4v7vbi - pParam->BSIM4v7phi;
                  tmp2 = model->BSIM4v7factor1 * tmp;

                  T0 = pParam->BSIM4v7dvt1w * pParam->BSIM4v7weff
                     * pParam->BSIM4v7leff / tmp2;
                  if (T0 < EXP_THRESHOLD)
                  {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T8 = T1 / T4;
                  }
                  else
                      T8 = 1.0 / (MAX_EXP - 2.0);
                  T0 = pParam->BSIM4v7dvt0w * T8;
                  T8 = T0 * tmp1;

                  T0 = pParam->BSIM4v7dvt1 * pParam->BSIM4v7leff / tmp2;
                  if (T0 < EXP_THRESHOLD)
                  {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T9 = T1 / T4;
                  }
                  else
                      T9 = 1.0 / (MAX_EXP - 2.0);
                  T9 = pParam->BSIM4v7dvt0 * T9 * tmp1;

                  T4 = toxe * pParam->BSIM4v7phi
                     / (pParam->BSIM4v7weff + pParam->BSIM4v7w0);

                  T0 = sqrt(1.0 + pParam->BSIM4v7lpe0 / pParam->BSIM4v7leff);
                  if((model->BSIM4v7tempMod == 1) || (model->BSIM4v7tempMod == 0))
                          T3 = (pParam->BSIM4v7kt1 + pParam->BSIM4v7kt1l / pParam->BSIM4v7leff)
                                     * (TRatio - 1.0);
                  if((model->BSIM4v7tempMod == 2)||(model->BSIM4v7tempMod == 3))
                        T3 = - pParam->BSIM4v7kt1 * (TRatio - 1.0);

                  T5 = pParam->BSIM4v7k1ox * (T0 - 1.0) * pParam->BSIM4v7sqrtPhi
                     + T3;
                  pParam->BSIM4v7vfbzbfactor = - T8 - T9 + pParam->BSIM4v7k3 * T4 + T5
                                             - pParam->BSIM4v7phi - pParam->BSIM4v7k1 * pParam->BSIM4v7sqrtPhi;

                  /* stress effect */

                      wlod = model->BSIM4v7wlod;
                      if (model->BSIM4v7wlod < 0.0)
                  {   fprintf(stderr, "Warning: WLOD = %g is less than 0. 0.0 is used\n",model->BSIM4v7wlod);
                             wlod = 0.0;
                  }
                  T0 = pow(Lnew, model->BSIM4v7llodku0);
                  W_tmp = Wnew + wlod;
                  T1 = pow(W_tmp, model->BSIM4v7wlodku0);
                  tmp1 = model->BSIM4v7lku0 / T0 + model->BSIM4v7wku0 / T1
                         + model->BSIM4v7pku0 / (T0 * T1);
                  pParam->BSIM4v7ku0 = 1.0 + tmp1;

                  T0 = pow(Lnew, model->BSIM4v7llodvth);
                  T1 = pow(W_tmp, model->BSIM4v7wlodvth);
                  tmp1 = model->BSIM4v7lkvth0 / T0 + model->BSIM4v7wkvth0 / T1
                       + model->BSIM4v7pkvth0 / (T0 * T1);
                  pParam->BSIM4v7kvth0 = 1.0 + tmp1;
                  pParam->BSIM4v7kvth0 = sqrt(pParam->BSIM4v7kvth0*pParam->BSIM4v7kvth0 + DELTA);

                  T0 = (TRatio - 1.0);
                  pParam->BSIM4v7ku0temp = pParam->BSIM4v7ku0 * (1.0 + model->BSIM4v7tku0 *T0) + DELTA;

                  Inv_saref = 1.0/(model->BSIM4v7saref + 0.5*Ldrn);
                  Inv_sbref = 1.0/(model->BSIM4v7sbref + 0.5*Ldrn);
                  pParam->BSIM4v7inv_od_ref = Inv_saref + Inv_sbref;
                  pParam->BSIM4v7rho_ref = model->BSIM4v7ku0 / pParam->BSIM4v7ku0temp * pParam->BSIM4v7inv_od_ref;

              } /* End of SizeNotFound */

              /*  stress effect */
              if( (here->BSIM4v7sa > 0.0) && (here->BSIM4v7sb > 0.0) &&
                        ((here->BSIM4v7nf == 1.0) || ((here->BSIM4v7nf > 1.0) && (here->BSIM4v7sd > 0.0))) )
              {          Inv_sa = 0;
                        Inv_sb = 0;

                         kvsat = model->BSIM4v7kvsat;
                  if (model->BSIM4v7kvsat < -1.0 )
                  {   fprintf(stderr, "Warning: KVSAT = %g is too small; -1.0 is used.\n",model->BSIM4v7kvsat);
                             kvsat = -1.0;
                      }
                      if (model->BSIM4v7kvsat > 1.0)
                      {   fprintf(stderr, "Warning: KVSAT = %g is too big; 1.0 is used.\n",model->BSIM4v7kvsat);
                         kvsat = 1.0;
                      }

                        for(i = 0; i < here->BSIM4v7nf; i++){
                           T0 = 1.0 / here->BSIM4v7nf / (here->BSIM4v7sa + 0.5*Ldrn + i * (here->BSIM4v7sd +Ldrn));
                            T1 = 1.0 / here->BSIM4v7nf / (here->BSIM4v7sb + 0.5*Ldrn + i * (here->BSIM4v7sd +Ldrn));
                           Inv_sa += T0;
                            Inv_sb += T1;
                  }
                  Inv_ODeff = Inv_sa + Inv_sb;
                  rho = model->BSIM4v7ku0 / pParam->BSIM4v7ku0temp * Inv_ODeff;
                  T0 = (1.0 + rho)/(1.0 + pParam->BSIM4v7rho_ref);
                  here->BSIM4v7u0temp = pParam->BSIM4v7u0temp * T0;

                  T1 = (1.0 + kvsat * rho)/(1.0 + kvsat * pParam->BSIM4v7rho_ref);
                  here->BSIM4v7vsattemp = pParam->BSIM4v7vsattemp * T1;

                  OD_offset = Inv_ODeff - pParam->BSIM4v7inv_od_ref;
                  dvth0_lod = model->BSIM4v7kvth0 / pParam->BSIM4v7kvth0 * OD_offset;
                  dk2_lod = model->BSIM4v7stk2 / pow(pParam->BSIM4v7kvth0, model->BSIM4v7lodk2) *
                                   OD_offset;
                  deta0_lod = model->BSIM4v7steta0 / pow(pParam->BSIM4v7kvth0, model->BSIM4v7lodeta0) *
                                     OD_offset;
                  here->BSIM4v7vth0 = pParam->BSIM4v7vth0 + dvth0_lod;

                  here->BSIM4v7eta0 = pParam->BSIM4v7eta0 + deta0_lod;
                  here->BSIM4v7k2 = pParam->BSIM4v7k2 + dk2_lod;
               } else {
                      here->BSIM4v7u0temp = pParam->BSIM4v7u0temp;
                      here->BSIM4v7vth0 = pParam->BSIM4v7vth0;
                      here->BSIM4v7vsattemp = pParam->BSIM4v7vsattemp;
                      here->BSIM4v7eta0 = pParam->BSIM4v7eta0;
                      here->BSIM4v7k2 = pParam->BSIM4v7k2;
              }

              /*  Well Proximity Effect  */
              if (model->BSIM4v7wpemod)
              { if( (!here->BSIM4v7scaGiven) && (!here->BSIM4v7scbGiven) && (!here->BSIM4v7sccGiven) )
                {   if((here->BSIM4v7scGiven) && (here->BSIM4v7sc > 0.0) )
                          {   T1 = here->BSIM4v7sc + Wdrn;
                        T2 = 1.0 / model->BSIM4v7scref;
                        here->BSIM4v7sca = model->BSIM4v7scref * model->BSIM4v7scref
                                        / (here->BSIM4v7sc * T1);
                        here->BSIM4v7scb = ( (0.1 * here->BSIM4v7sc + 0.01 * model->BSIM4v7scref)
                                        * exp(-10.0 * here->BSIM4v7sc * T2)
                                        - (0.1 * T1 + 0.01 * model->BSIM4v7scref)
                                        * exp(-10.0 * T1 * T2) ) / Wdrn;
                        here->BSIM4v7scc = ( (0.05 * here->BSIM4v7sc + 0.0025 * model->BSIM4v7scref)
                                        * exp(-20.0 * here->BSIM4v7sc * T2)
                                        - (0.05 * T1 + 0.0025 * model->BSIM4v7scref)
                                        * exp(-20.0 * T1 * T2) ) / Wdrn;
                    } else {
                        fprintf(stderr, "Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.\n");
                    }
                }

                       if (here->BSIM4v7sca < 0.0)
                {
                    printf("Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4v7sca);
                    here->BSIM4v7sca = 0.0;
                }
                if (here->BSIM4v7scb < 0.0)
                {
                    printf("Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4v7scb);
                    here->BSIM4v7scb = 0.0;
                }
                if (here->BSIM4v7scc < 0.0)
                {
                    printf("Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4v7scc);
                    here->BSIM4v7scc = 0.0;
                }
                if (here->BSIM4v7sc < 0.0)
                {
                    printf("Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4v7sc);
                    here->BSIM4v7sc = 0.0;
                }
                                /*4.6.2*/
                sceff = here->BSIM4v7sca + model->BSIM4v7web * here->BSIM4v7scb
                      + model->BSIM4v7wec * here->BSIM4v7scc;
                here->BSIM4v7vth0 += pParam->BSIM4v7kvth0we * sceff;
                here->BSIM4v7k2 +=  pParam->BSIM4v7k2we * sceff;
                  T3 =  1.0 + pParam->BSIM4v7ku0we * sceff;
                if (T3 <= 0.0)
                {         T3 = 0.0;
                        fprintf(stderr, "Warning: ku0we = %g is negatively too high. Negative mobility! \n", pParam->BSIM4v7ku0we);
                }
                here->BSIM4v7u0temp *= T3;
              }

            /* adding delvto  */
            here->BSIM4v7vth0 += here->BSIM4v7delvto;
            here->BSIM4v7vfb = pParam->BSIM4v7vfb + model->BSIM4v7type * here->BSIM4v7delvto;

            /* Instance variables calculation  */
            T3 = model->BSIM4v7type * here->BSIM4v7vth0
               - here->BSIM4v7vfb - pParam->BSIM4v7phi;
            T4 = T3 + T3;
            T5 = 2.5 * T3;
            here->BSIM4v7vtfbphi1 = (model->BSIM4v7type == NMOS) ? T4 : T5;
            if (here->BSIM4v7vtfbphi1 < 0.0)
                here->BSIM4v7vtfbphi1 = 0.0;

            here->BSIM4v7vtfbphi2 = 4.0 * T3;
            if (here->BSIM4v7vtfbphi2 < 0.0)
                here->BSIM4v7vtfbphi2 = 0.0;

            if (here->BSIM4v7k2 < 0.0)
            {   T0 = 0.5 * pParam->BSIM4v7k1 / here->BSIM4v7k2;
                here->BSIM4v7vbsc = 0.9 * (pParam->BSIM4v7phi - T0 * T0);
                if (here->BSIM4v7vbsc > -3.0)
                    here->BSIM4v7vbsc = -3.0;
                else if (here->BSIM4v7vbsc < -30.0)
                    here->BSIM4v7vbsc = -30.0;
            }
            else
                here->BSIM4v7vbsc = -30.0;
            if (here->BSIM4v7vbsc > pParam->BSIM4v7vbm)
                here->BSIM4v7vbsc = pParam->BSIM4v7vbm;
            here->BSIM4v7k2ox = here->BSIM4v7k2 * toxe
                              / model->BSIM4v7toxm;

            here->BSIM4v7vfbzb = pParam->BSIM4v7vfbzbfactor
                                +  model->BSIM4v7type * here->BSIM4v7vth0 ;

              here->BSIM4v7cgso = pParam->BSIM4v7cgso;
              here->BSIM4v7cgdo = pParam->BSIM4v7cgdo;

              lnl = log(pParam->BSIM4v7leff * 1.0e6);
              lnw = log(pParam->BSIM4v7weff * 1.0e6);
              lnnf = log(here->BSIM4v7nf);

              bodymode = 5;
              if( ( !model->BSIM4v7rbps0Given) ||
                  ( !model->BSIM4v7rbpd0Given) )
                bodymode = 1;
              else
                if( (!model->BSIM4v7rbsbx0Given && !model->BSIM4v7rbsby0Given) ||
                      (!model->BSIM4v7rbdbx0Given && !model->BSIM4v7rbdby0Given) )
                  bodymode = 3;

              if(here->BSIM4v7rbodyMod == 2)
                {
                  if (bodymode == 5)
                    {
                      rbsbx =  exp( log(model->BSIM4v7rbsbx0) + model->BSIM4v7rbsdbxl * lnl +
                                    model->BSIM4v7rbsdbxw * lnw + model->BSIM4v7rbsdbxnf * lnnf );
                      rbsby =  exp( log(model->BSIM4v7rbsby0) + model->BSIM4v7rbsdbyl * lnl +
                                    model->BSIM4v7rbsdbyw * lnw + model->BSIM4v7rbsdbynf * lnnf );
                      here->BSIM4v7rbsb = rbsbx * rbsby / (rbsbx + rbsby);


                      rbdbx =  exp( log(model->BSIM4v7rbdbx0) + model->BSIM4v7rbsdbxl * lnl +
                                    model->BSIM4v7rbsdbxw * lnw + model->BSIM4v7rbsdbxnf * lnnf );
                      rbdby =  exp( log(model->BSIM4v7rbdby0) + model->BSIM4v7rbsdbyl * lnl +
                                    model->BSIM4v7rbsdbyw * lnw + model->BSIM4v7rbsdbynf * lnnf );
                      here->BSIM4v7rbdb = rbdbx * rbdby / (rbdbx + rbdby);
                    }

                  if ((bodymode == 3)|| (bodymode == 5))
                    {
                      here->BSIM4v7rbps = exp( log(model->BSIM4v7rbps0) + model->BSIM4v7rbpsl * lnl +
                                             model->BSIM4v7rbpsw * lnw + model->BSIM4v7rbpsnf * lnnf );
                      here->BSIM4v7rbpd = exp( log(model->BSIM4v7rbpd0) + model->BSIM4v7rbpdl * lnl +
                                             model->BSIM4v7rbpdw * lnw + model->BSIM4v7rbpdnf * lnnf );
                    }

                  rbpbx =  exp( log(model->BSIM4v7rbpbx0) + model->BSIM4v7rbpbxl * lnl +
                                model->BSIM4v7rbpbxw * lnw + model->BSIM4v7rbpbxnf * lnnf );
                  rbpby =  exp( log(model->BSIM4v7rbpby0) + model->BSIM4v7rbpbyl * lnl +
                                model->BSIM4v7rbpbyw * lnw + model->BSIM4v7rbpbynf * lnnf );
                  here->BSIM4v7rbpb = rbpbx*rbpby/(rbpbx + rbpby);
                }


              if ((here->BSIM4v7rbodyMod == 1 ) || ((here->BSIM4v7rbodyMod == 2 ) && (bodymode == 5)) )
              {   if (here->BSIM4v7rbdb < 1.0e-3)
                      here->BSIM4v7grbdb = 1.0e3; /* in mho */
                  else
                      here->BSIM4v7grbdb = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbdb;
                  if (here->BSIM4v7rbpb < 1.0e-3)
                      here->BSIM4v7grbpb = 1.0e3;
                  else
                      here->BSIM4v7grbpb = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbpb;
                  if (here->BSIM4v7rbps < 1.0e-3)
                      here->BSIM4v7grbps = 1.0e3;
                  else
                      here->BSIM4v7grbps = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbps;
                  if (here->BSIM4v7rbsb < 1.0e-3)
                      here->BSIM4v7grbsb = 1.0e3;
                  else
                      here->BSIM4v7grbsb = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbsb;
                  if (here->BSIM4v7rbpd < 1.0e-3)
                      here->BSIM4v7grbpd = 1.0e3;
                  else
                      here->BSIM4v7grbpd = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbpd;

              }

              if((here->BSIM4v7rbodyMod == 2) && (bodymode == 3))
              {
                      here->BSIM4v7grbdb = here->BSIM4v7grbsb = model->BSIM4v7gbmin;
                  if (here->BSIM4v7rbpb < 1.0e-3)
                      here->BSIM4v7grbpb = 1.0e3;
                  else
                      here->BSIM4v7grbpb = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbpb;
                  if (here->BSIM4v7rbps < 1.0e-3)
                      here->BSIM4v7grbps = 1.0e3;
                  else
                      here->BSIM4v7grbps = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbps;
                  if (here->BSIM4v7rbpd < 1.0e-3)
                      here->BSIM4v7grbpd = 1.0e3;
                  else
                      here->BSIM4v7grbpd = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbpd;
              }

              if((here->BSIM4v7rbodyMod == 2) && (bodymode == 1))
              {
                      here->BSIM4v7grbdb = here->BSIM4v7grbsb = model->BSIM4v7gbmin;
                      here->BSIM4v7grbps = here->BSIM4v7grbpd = 1.0e3;
                  if (here->BSIM4v7rbpb < 1.0e-3)
                      here->BSIM4v7grbpb = 1.0e3;
                  else
                      here->BSIM4v7grbpb = model->BSIM4v7gbmin + 1.0 / here->BSIM4v7rbpb;
              }


              /*
               * Process geomertry dependent parasitics
               */

              here->BSIM4v7grgeltd = model->BSIM4v7rshg * (here->BSIM4v7xgw
                      + pParam->BSIM4v7weffCJ / 3.0 / here->BSIM4v7ngcon) /
                      (here->BSIM4v7ngcon * here->BSIM4v7nf *
                      (Lnew - model->BSIM4v7xgl));
              if (here->BSIM4v7grgeltd > 0.0)
                  here->BSIM4v7grgeltd = 1.0 / here->BSIM4v7grgeltd;
              else
              {   here->BSIM4v7grgeltd = 1.0e3; /* mho */
                  if (here->BSIM4v7rgateMod != 0)
                  printf("Warning: The gate conductance reset to 1.0e3 mho.\n");
              }

              DMCGeff = model->BSIM4v7dmcg - model->BSIM4v7dmcgt;
              DMCIeff = model->BSIM4v7dmci;
              DMDGeff = model->BSIM4v7dmdg - model->BSIM4v7dmcgt;

/*              if (here->BSIM4v7sourcePerimeterGiven)
              {   if (model->BSIM4v7perMod == 0)
                      here->BSIM4v7Pseff = here->BSIM4v7sourcePerimeter;
                  else
                      here->BSIM4v7Pseff = here->BSIM4v7sourcePerimeter
                                       - pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
              }
              else
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &(here->BSIM4v7Pseff), &dumPd, &dumAs, &dumAd);
              if (here->BSIM4v7Pseff < 0.0) /4.6.2/
                      here->BSIM4v7Pseff = 0.0; */

        /* New Diode Model v4.7*/
              if (here->BSIM4v7sourcePerimeterGiven)
              {   /* given */
                  if (here->BSIM4v7sourcePerimeter == 0.0)
                          here->BSIM4v7Pseff = 0.0;
                  else if (here->BSIM4v7sourcePerimeter < 0.0)
                  {
                          printf("Warning: Source Perimeter is specified as negative, it is set to zero.\n");
                          here->BSIM4v7Pseff = 0.0;
                  } else
                  {
                          if (model->BSIM4v7perMod == 0)
                                  here->BSIM4v7Pseff = here->BSIM4v7sourcePerimeter;
                          else
                                  here->BSIM4v7Pseff = here->BSIM4v7sourcePerimeter
                                          - pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
                  }
              } else /* not given */
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &(here->BSIM4v7Pseff), &dumPd, &dumAs, &dumAd);

              if (here->BSIM4v7Pseff < 0.0){ /* v4.7 final check */
                              here->BSIM4v7Pseff = 0.0;
                              printf("Warning: Pseff is negative, it is set to zero.\n");
              }
              /*  if (here->BSIM4v7drainPerimeterGiven)
              {   if (model->BSIM4v7perMod == 0)
                      here->BSIM4v7Pdeff = here->BSIM4v7drainPerimeter;
                  else
                      here->BSIM4v7Pdeff = here->BSIM4v7drainPerimeter
                                       - pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
              }
              else
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &(here->BSIM4v7Pdeff), &dumAs, &dumAd);
               if (here->BSIM4v7Pdeff < 0.0) /4.6.2/
                      here->BSIM4v7Pdeff = 0.0; */

              if (here->BSIM4v7drainPerimeterGiven)
              {   /* given */
                  if (here->BSIM4v7drainPerimeter == 0.0)
                                here->BSIM4v7Pdeff = 0.0;
                  else if (here->BSIM4v7drainPerimeter < 0.0)
                  {
                                printf("Warning: Drain Perimeter is specified as negative, it is set to zero.\n");
                                here->BSIM4v7Pdeff = 0.0;
                  } else
                  {
                                if (model->BSIM4v7perMod == 0)
                                        here->BSIM4v7Pdeff = here->BSIM4v7drainPerimeter;
                                else
                                        here->BSIM4v7Pdeff = here->BSIM4v7drainPerimeter
                                                - pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
                  }
              } else /* not given */
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &(here->BSIM4v7Pdeff), &dumAs, &dumAd);

              if (here->BSIM4v7Pdeff < 0.0){
                      here->BSIM4v7Pdeff = 0.0; /*New Diode v4.7*/
                      printf("Warning: Pdeff is negative, it is set to zero.\n");
              }
              if (here->BSIM4v7sourceAreaGiven)
                  here->BSIM4v7Aseff = here->BSIM4v7sourceArea;
              else
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &dumPd, &(here->BSIM4v7Aseff), &dumAd);
              if (here->BSIM4v7Aseff < 0.0){
                  here->BSIM4v7Aseff = 0.0; /* v4.7 */
                  printf("Warning: Aseff is negative, it is set to zero.\n");
              }
              if (here->BSIM4v7drainAreaGiven)
                  here->BSIM4v7Adeff = here->BSIM4v7drainArea;
              else
                  BSIM4v7PAeffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod, here->BSIM4v7min,
                                    pParam->BSIM4v7weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &dumPd, &dumAs, &(here->BSIM4v7Adeff));
              if (here->BSIM4v7Adeff < 0.0){
                  here->BSIM4v7Adeff = 0.0; /* v4.7 */
                  printf("Warning: Adeff is negative, it is set to zero.\n");
              }
              /* Processing S/D resistance and conductance below */
              if(here->BSIM4v7sNodePrime != here->BSIM4v7sNode)
              {
                 here->BSIM4v7sourceConductance = 0.0;
                 if(here->BSIM4v7sourceSquaresGiven)
                 {
                    here->BSIM4v7sourceConductance = model->BSIM4v7sheetResistance
                                               * here->BSIM4v7sourceSquares;
                 } else if (here->BSIM4v7rgeoMod > 0)
                 {
                    BSIM4v7RdseffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod,
                      here->BSIM4v7rgeoMod, here->BSIM4v7min,
                      pParam->BSIM4v7weffCJ, model->BSIM4v7sheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 1, &(here->BSIM4v7sourceConductance));
                 } else
                 {
                    here->BSIM4v7sourceConductance = 0.0;
                 }

                 if (here->BSIM4v7sourceConductance > 0.0)
                     here->BSIM4v7sourceConductance = 1.0
                                            / here->BSIM4v7sourceConductance;
                 else
                 {
                     here->BSIM4v7sourceConductance = 1.0e3; /* mho */
                     printf ("Warning: Source conductance reset to 1.0e3 mho.\n");
                 }
              } else
              {
                  here->BSIM4v7sourceConductance = 0.0;
              }

              if(here->BSIM4v7dNodePrime != here->BSIM4v7dNode)
              {
                 here->BSIM4v7drainConductance = 0.0;
                 if(here->BSIM4v7drainSquaresGiven)
                 {
                    here->BSIM4v7drainConductance = model->BSIM4v7sheetResistance
                                              * here->BSIM4v7drainSquares;
                 } else if (here->BSIM4v7rgeoMod > 0)
                 {
                    BSIM4v7RdseffGeo(here->BSIM4v7nf, here->BSIM4v7geoMod,
                      here->BSIM4v7rgeoMod, here->BSIM4v7min,
                      pParam->BSIM4v7weffCJ, model->BSIM4v7sheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 0, &(here->BSIM4v7drainConductance));
                 } else
                 {
                    here->BSIM4v7drainConductance = 0.0;
                 }

                 if (here->BSIM4v7drainConductance > 0.0)
                     here->BSIM4v7drainConductance = 1.0
                                           / here->BSIM4v7drainConductance;
                 else
                 {
                     here->BSIM4v7drainConductance = 1.0e3; /* mho */
                     printf ("Warning: Drain conductance reset to 1.0e3 mho.\n");
                  }
              } else
              {
                  here->BSIM4v7drainConductance = 0.0;
              }

               /* End of Rsd processing */


              Nvtms = model->BSIM4v7vtm * model->BSIM4v7SjctEmissionCoeff;
              if ((here->BSIM4v7Aseff <= 0.0) && (here->BSIM4v7Pseff <= 0.0))
              {   SourceSatCurrent = 0.0; /* v4.7 */
                  /* SourceSatCurrent = 1.0e-14; */
              }
              else
              {   SourceSatCurrent = here->BSIM4v7Aseff * model->BSIM4v7SjctTempSatCurDensity
                                   + here->BSIM4v7Pseff * model->BSIM4v7SjctSidewallTempSatCurDensity
                                   + pParam->BSIM4v7weffCJ * here->BSIM4v7nf
                                   * model->BSIM4v7SjctGateSidewallTempSatCurDensity;
              }
              if (SourceSatCurrent > 0.0)
              {   switch(model->BSIM4v7dioMod)
                  {   case 0:
                          if ((model->BSIM4v7bvs / Nvtms) > EXP_THRESHOLD)
                              here->BSIM4v7XExpBVS = model->BSIM4v7xjbvs * MIN_EXP;
                          else
                              here->BSIM4v7XExpBVS = model->BSIM4v7xjbvs * exp(-model->BSIM4v7bvs / Nvtms);
                          break;
                      case 1:
                          BSIM4v7DioIjthVjmEval(Nvtms, model->BSIM4v7ijthsfwd, SourceSatCurrent,
                                              0.0, &(here->BSIM4v7vjsmFwd));
                          here->BSIM4v7IVjsmFwd = SourceSatCurrent * exp(here->BSIM4v7vjsmFwd / Nvtms);
                          break;
                      case 2:
                          if ((model->BSIM4v7bvs / Nvtms) > EXP_THRESHOLD)
                          {   here->BSIM4v7XExpBVS = model->BSIM4v7xjbvs * MIN_EXP;
                              tmp = MIN_EXP;
                          }
                          else
                          {   here->BSIM4v7XExpBVS = exp(-model->BSIM4v7bvs / Nvtms);
                              tmp = here->BSIM4v7XExpBVS;
                              here->BSIM4v7XExpBVS *= model->BSIM4v7xjbvs;
                          }

                          BSIM4v7DioIjthVjmEval(Nvtms, model->BSIM4v7ijthsfwd, SourceSatCurrent,
                                                     here->BSIM4v7XExpBVS, &(here->BSIM4v7vjsmFwd));
                          T0 = exp(here->BSIM4v7vjsmFwd / Nvtms);
                          here->BSIM4v7IVjsmFwd = SourceSatCurrent * (T0 - here->BSIM4v7XExpBVS / T0
                                                + here->BSIM4v7XExpBVS - 1.0);
                          here->BSIM4v7SslpFwd = SourceSatCurrent
                                               * (T0 + here->BSIM4v7XExpBVS / T0) / Nvtms;

                          T2 = model->BSIM4v7ijthsrev / SourceSatCurrent;
                          if (T2 < 1.0)
                          {   T2 = 10.0;
                              fprintf(stderr, "Warning: ijthsrev too small and set to 10 times IsbSat.\n");
                          }
                          here->BSIM4v7vjsmRev = -model->BSIM4v7bvs
                                             - Nvtms * log((T2 - 1.0) / model->BSIM4v7xjbvs);
                          T1 = model->BSIM4v7xjbvs * exp(-(model->BSIM4v7bvs
                             + here->BSIM4v7vjsmRev) / Nvtms);
                          here->BSIM4v7IVjsmRev = SourceSatCurrent * (1.0 + T1);
                          here->BSIM4v7SslpRev = -SourceSatCurrent * T1 / Nvtms;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model->BSIM4v7dioMod);
                  }
              }

              Nvtmd = model->BSIM4v7vtm * model->BSIM4v7DjctEmissionCoeff;
              if ((here->BSIM4v7Adeff <= 0.0) && (here->BSIM4v7Pdeff <= 0.0))
              {  /* DrainSatCurrent = 1.0e-14;         v4.7 */
                   DrainSatCurrent = 0.0;
              }
              else
              {   DrainSatCurrent = here->BSIM4v7Adeff * model->BSIM4v7DjctTempSatCurDensity
                                  + here->BSIM4v7Pdeff * model->BSIM4v7DjctSidewallTempSatCurDensity
                                  + pParam->BSIM4v7weffCJ * here->BSIM4v7nf
                                  * model->BSIM4v7DjctGateSidewallTempSatCurDensity;
              }
              if (DrainSatCurrent > 0.0)
              {   switch(model->BSIM4v7dioMod)
                  {   case 0:
                          if ((model->BSIM4v7bvd / Nvtmd) > EXP_THRESHOLD)
                          here->BSIM4v7XExpBVD = model->BSIM4v7xjbvd * MIN_EXP;
                          else
                          here->BSIM4v7XExpBVD = model->BSIM4v7xjbvd * exp(-model->BSIM4v7bvd / Nvtmd);
                          break;
                      case 1:
                          BSIM4v7DioIjthVjmEval(Nvtmd, model->BSIM4v7ijthdfwd, DrainSatCurrent,
                                              0.0, &(here->BSIM4v7vjdmFwd));
                          here->BSIM4v7IVjdmFwd = DrainSatCurrent * exp(here->BSIM4v7vjdmFwd / Nvtmd);
                          break;
                      case 2:
                          if ((model->BSIM4v7bvd / Nvtmd) > EXP_THRESHOLD)
                          {   here->BSIM4v7XExpBVD = model->BSIM4v7xjbvd * MIN_EXP;
                              tmp = MIN_EXP;
                          }
                          else
                          {   here->BSIM4v7XExpBVD = exp(-model->BSIM4v7bvd / Nvtmd);
                              tmp = here->BSIM4v7XExpBVD;
                              here->BSIM4v7XExpBVD *= model->BSIM4v7xjbvd;
                          }

                          BSIM4v7DioIjthVjmEval(Nvtmd, model->BSIM4v7ijthdfwd, DrainSatCurrent,
                                              here->BSIM4v7XExpBVD, &(here->BSIM4v7vjdmFwd));
                          T0 = exp(here->BSIM4v7vjdmFwd / Nvtmd);
                          here->BSIM4v7IVjdmFwd = DrainSatCurrent * (T0 - here->BSIM4v7XExpBVD / T0
                                              + here->BSIM4v7XExpBVD - 1.0);
                          here->BSIM4v7DslpFwd = DrainSatCurrent
                                               * (T0 + here->BSIM4v7XExpBVD / T0) / Nvtmd;

                          T2 = model->BSIM4v7ijthdrev / DrainSatCurrent;
                          if (T2 < 1.0)
                          {   T2 = 10.0;
                              fprintf(stderr, "Warning: ijthdrev too small and set to 10 times IdbSat.\n");
                          }
                          here->BSIM4v7vjdmRev = -model->BSIM4v7bvd
                                             - Nvtmd * log((T2 - 1.0) / model->BSIM4v7xjbvd); /* bugfix */
                          T1 = model->BSIM4v7xjbvd * exp(-(model->BSIM4v7bvd
                             + here->BSIM4v7vjdmRev) / Nvtmd);
                          here->BSIM4v7IVjdmRev = DrainSatCurrent * (1.0 + T1);
                          here->BSIM4v7DslpRev = -DrainSatCurrent * T1 / Nvtmd;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model->BSIM4v7dioMod);
                  }
              }

                /* GEDL current reverse bias */
                T0 = (TRatio - 1.0);
                model->BSIM4v7njtsstemp = model->BSIM4v7njts * (1.0 + model->BSIM4v7tnjts * T0);
                model->BSIM4v7njtsswstemp = model->BSIM4v7njtssw * (1.0 + model->BSIM4v7tnjtssw * T0);
                model->BSIM4v7njtsswgstemp = model->BSIM4v7njtsswg * (1.0 + model->BSIM4v7tnjtsswg * T0);
                model->BSIM4v7njtsdtemp = model->BSIM4v7njtsd * (1.0 + model->BSIM4v7tnjtsd * T0);
                model->BSIM4v7njtsswdtemp = model->BSIM4v7njtsswd * (1.0 + model->BSIM4v7tnjtsswd * T0);
                model->BSIM4v7njtsswgdtemp = model->BSIM4v7njtsswgd * (1.0 + model->BSIM4v7tnjtsswgd * T0);
                T7 = Eg0 / model->BSIM4v7vtm * T0;
                T9 = model->BSIM4v7xtss * T7;
                DEXP(T9, T1);
                T9 = model->BSIM4v7xtsd * T7;
                DEXP(T9, T2);
                T9 = model->BSIM4v7xtssws * T7;
                DEXP(T9, T3);
                T9 = model->BSIM4v7xtsswd * T7;
                DEXP(T9, T4);
                T9 = model->BSIM4v7xtsswgs * T7;
                DEXP(T9, T5);
                T9 = model->BSIM4v7xtsswgd * T7;
                DEXP(T9, T6);
                                /*IBM TAT*/
                                if(model->BSIM4v7jtweff < 0.0)
                              {   model->BSIM4v7jtweff = 0.0;
                          fprintf(stderr, "TAT width dependence effect is negative. Jtweff is clamped to zero.\n");
                      }
                                T11 = sqrt(model->BSIM4v7jtweff / pParam->BSIM4v7weffCJ) + 1.0;

                T10 = pParam->BSIM4v7weffCJ * here->BSIM4v7nf;
                here->BSIM4v7SjctTempRevSatCur = T1 * here->BSIM4v7Aseff * model->BSIM4v7jtss;
                here->BSIM4v7DjctTempRevSatCur = T2 * here->BSIM4v7Adeff * model->BSIM4v7jtsd;
                here->BSIM4v7SswTempRevSatCur = T3 * here->BSIM4v7Pseff * model->BSIM4v7jtssws;
                here->BSIM4v7DswTempRevSatCur = T4 * here->BSIM4v7Pdeff * model->BSIM4v7jtsswd;
                here->BSIM4v7SswgTempRevSatCur = T5 * T10 * T11 * model->BSIM4v7jtsswgs;
                here->BSIM4v7DswgTempRevSatCur = T6 * T10 * T11 * model->BSIM4v7jtsswgd;

                /*high k*/
                /*Calculate VgsteffVth for mobMod=3*/
                if(model->BSIM4v7mobMod==3)
                {        /*Calculate n @ Vbs=Vds=0*/
                    V0 = pParam->BSIM4v7vbi - pParam->BSIM4v7phi;
                    lt1 = model->BSIM4v7factor1* pParam->BSIM4v7sqrtXdep0;
                    ltw = lt1;
                    T0 = pParam->BSIM4v7dvt1 * pParam->BSIM4v7leff / lt1;
                    if (T0 < EXP_THRESHOLD)
                      {
                        T1 = exp(T0);
                        T2 = T1 - 1.0;
                        T3 = T2 * T2;
                        T4 = T3 + 2.0 * T1 * MIN_EXP;
                        Theta0 = T1 / T4;
                      }
                    else
                      Theta0 = 1.0 / (MAX_EXP - 2.0);

                     tmp1 = epssub / pParam->BSIM4v7Xdep0;
                    here->BSIM4v7nstar = model->BSIM4v7vtm / Charge_q *
                      (model->BSIM4v7coxe        + tmp1 + pParam->BSIM4v7cit);
                    tmp2 = pParam->BSIM4v7nfactor * tmp1;
                    tmp3 = (tmp2 + pParam->BSIM4v7cdsc * Theta0 + pParam->BSIM4v7cit) / model->BSIM4v7coxe;
                    if (tmp3 >= -0.5)
                      n0 = 1.0 + tmp3;
                    else
                      {
                        T0 = 1.0 / (3.0 + 8.0 * tmp3);
                        n0 = (1.0 + 3.0 * tmp3) * T0;
                      }

                  T0 = n0 * model->BSIM4v7vtm;
                  T1 = pParam->BSIM4v7voffcbn;
                  T2 = T1/T0;
                  if (T2 < -EXP_THRESHOLD)
                  {   T3 = model->BSIM4v7coxe * MIN_EXP / pParam->BSIM4v7cdep0;
                      T4 = pParam->BSIM4v7mstar + T3 * n0;
                  }
                  else if (T2 > EXP_THRESHOLD)
                  {   T3 = model->BSIM4v7coxe * MAX_EXP / pParam->BSIM4v7cdep0;
                      T4 = pParam->BSIM4v7mstar + T3 * n0;
                  }
                  else
                  {  T3 = exp(T2)* model->BSIM4v7coxe / pParam->BSIM4v7cdep0;
                               T4 = pParam->BSIM4v7mstar + T3 * n0;
                  }
                  pParam->BSIM4v7VgsteffVth = T0 * log(2.0)/T4;

                }

                /* New DITS term added in 4.7 */
                T0 = -pParam->BSIM4v7dvtp3 * log(pParam->BSIM4v7leff);
                DEXP(T0, T1);
                pParam->BSIM4v7dvtp2factor = pParam->BSIM4v7dvtp5 + pParam->BSIM4v7dvtp2 * T1;

                if(model->BSIM4v7mtrlMod != 0 && model->BSIM4v7mtrlCompatMod == 0)
                {
                    /* Calculate TOXP from EOT */
                    /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
                    Vtm0eot = KboQ * model->BSIM4v7tempeot;
                    Vtmeot  = Vtm0eot;
                    vbieot = Vtm0eot * log(pParam->BSIM4v7nsd
                                   * pParam->BSIM4v7ndep / (ni * ni));
                    phieot = Vtm0eot * log(pParam->BSIM4v7ndep / ni)
                                   + pParam->BSIM4v7phin + 0.4;
                    tmp2 = here->BSIM4v7vfb + phieot;
                    vddeot = model->BSIM4v7type * model->BSIM4v7vddeot;
                    T0 = model->BSIM4v7epsrgate * EPS0;
                    if ((pParam->BSIM4v7ngate > 1.0e18) && (pParam->BSIM4v7ngate < 1.0e25)
                        && (vddeot > tmp2) && (T0!=0))
                      {
                        T1 = 1.0e6 * CHARGE * T0 * pParam->BSIM4v7ngate /
                          (model->BSIM4v7coxe * model->BSIM4v7coxe);
                        T8 = vddeot - tmp2;
                        T4 = sqrt(1.0 + 2.0 * T8 / T1);
                        T2 = 2.0 * T8 / (T4 + 1.0);
                        T3 = 0.5 * T2 * T2 / T1;
                        T7 = 1.12 - T3 - 0.05;
                        T6 = sqrt(T7 * T7 + 0.224);
                        T5 = 1.12 - 0.5 * (T7 + T6);
                        Vgs_eff = vddeot - T5;
                      }
                    else
                      Vgs_eff = vddeot;

                    /* Calculate Vth @ Vds=Vbs=0 */

                    V0 = vbieot - phieot;
                    lt1 = model->BSIM4v7factor1* pParam->BSIM4v7sqrtXdep0;
                    ltw = lt1;
                    T0 = pParam->BSIM4v7dvt1 * model->BSIM4v7leffeot / lt1;
                    if (T0 < EXP_THRESHOLD)
                      {
                        T1 = exp(T0);
                        T2 = T1 - 1.0;
                        T3 = T2 * T2;
                        T4 = T3 + 2.0 * T1 * MIN_EXP;
                        Theta0 = T1 / T4;
                      }
                    else
                      Theta0 = 1.0 / (MAX_EXP - 2.0);
                    Delt_vth = pParam->BSIM4v7dvt0 * Theta0 * V0;
                    T0 = pParam->BSIM4v7dvt1w * model->BSIM4v7weffeot * model->BSIM4v7leffeot / ltw;
                    if (T0 < EXP_THRESHOLD)
                      {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T5 = T1 / T4;
                      }
                    else
                      T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                    T2 = pParam->BSIM4v7dvt0w * T5 * V0;
                    TempRatioeot =  model->BSIM4v7tempeot / model->BSIM4v7tnom - 1.0;
                    T0 = sqrt(1.0 + pParam->BSIM4v7lpe0 / model->BSIM4v7leffeot);
                    T1 = pParam->BSIM4v7k1ox * (T0 - 1.0) * sqrt(phieot)
                      + (pParam->BSIM4v7kt1 + pParam->BSIM4v7kt1l / model->BSIM4v7leffeot) * TempRatioeot;
                    Vth_NarrowW = toxe * phieot
                      / (model->BSIM4v7weffeot + pParam->BSIM4v7w0);
                    Lpe_Vb = sqrt(1.0 + pParam->BSIM4v7lpeb / model->BSIM4v7leffeot);
                    Vth = model->BSIM4v7type * here->BSIM4v7vth0 +
                      (pParam->BSIM4v7k1ox - pParam->BSIM4v7k1)*sqrt(phieot)*Lpe_Vb
                      - Delt_vth - T2 + pParam->BSIM4v7k3 * Vth_NarrowW + T1;

                    /* Calculate n */
                    tmp1 = epssub / pParam->BSIM4v7Xdep0;
                    here->BSIM4v7nstar = Vtmeot / Charge_q *
                      (model->BSIM4v7coxe        + tmp1 + pParam->BSIM4v7cit);
                    tmp2 = pParam->BSIM4v7nfactor * tmp1;
                    tmp3 = (tmp2 + pParam->BSIM4v7cdsc * Theta0 + pParam->BSIM4v7cit) / model->BSIM4v7coxe;
                    if (tmp3 >= -0.5)
                      n = 1.0 + tmp3;
                    else
                      {
                        T0 = 1.0 / (3.0 + 8.0 * tmp3);
                        n = (1.0 + 3.0 * tmp3) * T0;
                      }

                    /* Vth correction for Pocket implant */
                    if (pParam->BSIM4v7dvtp0 > 0.0)
                      {
                        T3 = model->BSIM4v7leffeot + pParam->BSIM4v7dvtp0 * 2.0;
                        if (model->BSIM4v7tempMod < 2)
                          T4 = Vtmeot * log(model->BSIM4v7leffeot / T3);
                        else
                          T4 = Vtm0eot * log(model->BSIM4v7leffeot / T3);
                        Vth -= n * T4;
                      }
                    Vgsteff = Vgs_eff-Vth;
                    /* calculating Toxp */
                        T3 = model->BSIM4v7type * here->BSIM4v7vth0
               - here->BSIM4v7vfb - phieot;
            T4 = T3 + T3;
            T5 = 2.5 * T3;

            vtfbphi2eot = 4.0 * T3;
            if (vtfbphi2eot < 0.0)
                vtfbphi2eot = 0.0;


                    niter = 0;
                    toxpf = toxe;
                    do
                      {
                        toxpi = toxpf;
                        tmp2 = 2.0e8 * toxpf;
                        T0 = (Vgsteff + vtfbphi2eot) / tmp2;
                        T1 = 1.0 + exp(model->BSIM4v7bdos * 0.7 * log(T0));
                        Tcen = model->BSIM4v7ados * 1.9e-9 / T1;
                        toxpf = toxe - epsrox/model->BSIM4v7epsrsub * Tcen;
                        niter++;
                      } while ((niter<=4)&&(ABS(toxpf-toxpi)>1e-12));
                      model->BSIM4v7toxp = toxpf;
                      model->BSIM4v7coxp = epsrox * EPS0 / model->BSIM4v7toxp;
                      }

              if (BSIM4v7checkModel(model, here, ckt))
              {
                  SPfrontEnd->IFerrorf (ERR_FATAL, "Fatal error(s) detected during BSIM4v7.6.0 parameter checking for %s in model %s", model->BSIM4v7modName, here->BSIM4v7name);
                  return(E_BADPARM);
              }
         } /* End instance */
    }
    return(OK);
}
