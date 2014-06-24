/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4temp.c of BSIM4.7.0.
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
 * Modified by Pankaj Kumar Thakur, 07/23/2012
 **********/


#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
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
BSIM4DioIjthVjmEval(
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
BSIM4temp(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4model *model = (BSIM4model*) inModel;
BSIM4instance *here;
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

    /*  loop through all the BSIM4 device models */
    for (; model != NULL; model = model->BSIM4nextModel)
    {    Temp = ckt->CKTtemp;
         if (model->BSIM4SbulkJctPotential < 0.1)
         {   model->BSIM4SbulkJctPotential = 0.1;
             fprintf(stderr, "Given pbs is less than 0.1. Pbs is set to 0.1.\n");
         }
         if (model->BSIM4SsidewallJctPotential < 0.1)
         {   model->BSIM4SsidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbsws is less than 0.1. Pbsws is set to 0.1.\n");
         }
         if (model->BSIM4SGatesidewallJctPotential < 0.1)
         {   model->BSIM4SGatesidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswgs is less than 0.1. Pbswgs is set to 0.1.\n");
         }

         if (model->BSIM4DbulkJctPotential < 0.1)
         {   model->BSIM4DbulkJctPotential = 0.1;
             fprintf(stderr, "Given pbd is less than 0.1. Pbd is set to 0.1.\n");
         }
         if (model->BSIM4DsidewallJctPotential < 0.1)
         {   model->BSIM4DsidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswd is less than 0.1. Pbswd is set to 0.1.\n");
         }
         if (model->BSIM4DGatesidewallJctPotential < 0.1)
         {   model->BSIM4DGatesidewallJctPotential = 0.1;
             fprintf(stderr, "Given pbswgd is less than 0.1. Pbswgd is set to 0.1.\n");
         }

         if(model->BSIM4mtrlMod == 0)
         {
             if ((model->BSIM4toxeGiven) && (model->BSIM4toxpGiven) && (model->BSIM4dtoxGiven)
                 && (model->BSIM4toxe != (model->BSIM4toxp + model->BSIM4dtox)))
                 printf("Warning: toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored.\n");
             else if ((model->BSIM4toxeGiven) && (!model->BSIM4toxpGiven))
               model->BSIM4toxp = model->BSIM4toxe - model->BSIM4dtox;
             else if ((!model->BSIM4toxeGiven) && (model->BSIM4toxpGiven)){
               model->BSIM4toxe = model->BSIM4toxp + model->BSIM4dtox;
                 if (!model->BSIM4toxmGiven)                        /* v4.7 */
                     model->BSIM4toxm = model->BSIM4toxe;
             }
         }
         else if(model->BSIM4mtrlCompatMod != 0) /* v4.7 */
         {
             T0 = model->BSIM4epsrox / 3.9;
             if ((model->BSIM4eotGiven) && (model->BSIM4toxpGiven) && (model->BSIM4dtoxGiven)
                 && (ABS(model->BSIM4eot * T0 - (model->BSIM4toxp + model->BSIM4dtox)) > 1.0e-20))
             {
                 printf("Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.\n");
             }
             else if ((model->BSIM4eotGiven) && (!model->BSIM4toxpGiven))
               model->BSIM4toxp = T0 * model->BSIM4eot - model->BSIM4dtox;
             else if ((!model->BSIM4eotGiven) && (model->BSIM4toxpGiven)){
               model->BSIM4eot = (model->BSIM4toxp + model->BSIM4dtox) / T0;
                 if (!model->BSIM4toxmGiven)
                     model->BSIM4toxm = model->BSIM4eot;
             }
         }

         if(model->BSIM4mtrlMod)
           {
             epsrox = 3.9;
             toxe = model->BSIM4eot;
             epssub = EPS0 * model->BSIM4epsrsub;
           }
         else
           {
             epsrox = model->BSIM4epsrox;
             toxe = model->BSIM4toxe;
             epssub = EPSSI;
           }


         model->BSIM4coxe = epsrox * EPS0 / toxe;
         if(model->BSIM4mtrlMod == 0 || model->BSIM4mtrlCompatMod != 0)
           model->BSIM4coxp = model->BSIM4epsrox * EPS0 / model->BSIM4toxp;

         if (!model->BSIM4cgdoGiven)
         {   if (model->BSIM4dlcGiven && (model->BSIM4dlc > 0.0))
                 model->BSIM4cgdo = model->BSIM4dlc * model->BSIM4coxe
                                  - model->BSIM4cgdl ;
             else
                 model->BSIM4cgdo = 0.6 * model->BSIM4xj * model->BSIM4coxe;
         }
         if (!model->BSIM4cgsoGiven)
         {   if (model->BSIM4dlcGiven && (model->BSIM4dlc > 0.0))
                 model->BSIM4cgso = model->BSIM4dlc * model->BSIM4coxe
                                  - model->BSIM4cgsl ;
             else
                 model->BSIM4cgso = 0.6 * model->BSIM4xj * model->BSIM4coxe;
         }
         if (!model->BSIM4cgboGiven)
             model->BSIM4cgbo = 2.0 * model->BSIM4dwc * model->BSIM4coxe;
         model->pSizeDependParamKnot = NULL;
         pLastKnot = NULL;

         Tnom = model->BSIM4tnom;
         TRatio = Temp / Tnom;

         model->BSIM4vcrit = CONSTvt0 * log(CONSTvt0 / (CONSTroot2 * 1.0e-14));
         model->BSIM4factor1 = sqrt(epssub / (epsrox * EPS0)* toxe);

         Vtm0 = model->BSIM4vtm0 = KboQ * Tnom;

         if(model->BSIM4mtrlMod==0)
         {
             Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
             ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15)
                 * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
         }
         else
         {
           Eg0 = model->BSIM4bg0sub - model->BSIM4tbgasub * Tnom * Tnom
                                      / (Tnom + model->BSIM4tbgbsub);
           T0 =  model->BSIM4bg0sub - model->BSIM4tbgasub * 90090.0225
                                      / (300.15 + model->BSIM4tbgbsub);
           ni = model->BSIM4ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15)
                 * exp((T0 - Eg0) / (2.0 * Vtm0));
         }

         model->BSIM4Eg0 = Eg0;
         model->BSIM4vtm = KboQ * Temp;
         if(model->BSIM4mtrlMod == 0)
           Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
         else
           Eg = model->BSIM4bg0sub - model->BSIM4tbgasub * Temp * Temp
                                      / (Temp + model->BSIM4tbgbsub);
         if (Temp != Tnom)
         {   T0 = Eg0 / Vtm0 - Eg / model->BSIM4vtm;
             T1 = log(Temp / Tnom);
             T2 = T0 + model->BSIM4SjctTempExponent * T1;
             T3 = exp(T2 / model->BSIM4SjctEmissionCoeff);
             model->BSIM4SjctTempSatCurDensity = model->BSIM4SjctSatCurDensity
                                               * T3;
             model->BSIM4SjctSidewallTempSatCurDensity
                         = model->BSIM4SjctSidewallSatCurDensity * T3;
             model->BSIM4SjctGateSidewallTempSatCurDensity
                         = model->BSIM4SjctGateSidewallSatCurDensity * T3;

             T2 = T0 + model->BSIM4DjctTempExponent * T1;
             T3 = exp(T2 / model->BSIM4DjctEmissionCoeff);
             model->BSIM4DjctTempSatCurDensity = model->BSIM4DjctSatCurDensity
                                               * T3;
             model->BSIM4DjctSidewallTempSatCurDensity
                         = model->BSIM4DjctSidewallSatCurDensity * T3;
             model->BSIM4DjctGateSidewallTempSatCurDensity
                         = model->BSIM4DjctGateSidewallSatCurDensity * T3;
         }
         else
         {   model->BSIM4SjctTempSatCurDensity = model->BSIM4SjctSatCurDensity;
             model->BSIM4SjctSidewallTempSatCurDensity
                        = model->BSIM4SjctSidewallSatCurDensity;
             model->BSIM4SjctGateSidewallTempSatCurDensity
                        = model->BSIM4SjctGateSidewallSatCurDensity;
             model->BSIM4DjctTempSatCurDensity = model->BSIM4DjctSatCurDensity;
             model->BSIM4DjctSidewallTempSatCurDensity
                        = model->BSIM4DjctSidewallSatCurDensity;
             model->BSIM4DjctGateSidewallTempSatCurDensity
                        = model->BSIM4DjctGateSidewallSatCurDensity;
         }

         if (model->BSIM4SjctTempSatCurDensity < 0.0)
             model->BSIM4SjctTempSatCurDensity = 0.0;
         if (model->BSIM4SjctSidewallTempSatCurDensity < 0.0)
             model->BSIM4SjctSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4SjctGateSidewallTempSatCurDensity < 0.0)
             model->BSIM4SjctGateSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4DjctTempSatCurDensity < 0.0)
             model->BSIM4DjctTempSatCurDensity = 0.0;
         if (model->BSIM4DjctSidewallTempSatCurDensity < 0.0)
             model->BSIM4DjctSidewallTempSatCurDensity = 0.0;
         if (model->BSIM4DjctGateSidewallTempSatCurDensity < 0.0)
             model->BSIM4DjctGateSidewallTempSatCurDensity = 0.0;

         /* Temperature dependence of D/B and S/B diode capacitance begins */
         delTemp = ckt->CKTtemp - model->BSIM4tnom;
         T0 = model->BSIM4tcj * delTemp;
         if (T0 >= -1.0)
         {   model->BSIM4SunitAreaTempJctCap = model->BSIM4SunitAreaJctCap *(1.0 + T0); /*bug_fix -JX */
             model->BSIM4DunitAreaTempJctCap = model->BSIM4DunitAreaJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4SunitAreaJctCap > 0.0)
             {   model->BSIM4SunitAreaTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjs to be negative. Cjs is clamped to zero.\n");
             }
             if (model->BSIM4DunitAreaJctCap > 0.0)
             {   model->BSIM4DunitAreaTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjd to be negative. Cjd is clamped to zero.\n");
             }
         }
         T0 = model->BSIM4tcjsw * delTemp;
                   if (model->BSIM4SunitLengthSidewallJctCap < 0.0)/*4.6.2*/
                      {model->BSIM4SunitLengthSidewallJctCap = 0.0;
                           fprintf(stderr, "CJSWS is negative. Cjsws is clamped to zero.\n");}
                  if (model->BSIM4DunitLengthSidewallJctCap < 0.0)
                      {model->BSIM4DunitLengthSidewallJctCap = 0.0;
                           fprintf(stderr, "CJSWD is negative. Cjswd is clamped to zero.\n");}
         if (T0 >= -1.0)
         {   model->BSIM4SunitLengthSidewallTempJctCap = model->BSIM4SunitLengthSidewallJctCap *(1.0 + T0);
             model->BSIM4DunitLengthSidewallTempJctCap = model->BSIM4DunitLengthSidewallJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4SunitLengthSidewallJctCap > 0.0)
             {   model->BSIM4SunitLengthSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.\n");
             }
             if (model->BSIM4DunitLengthSidewallJctCap > 0.0)
             {   model->BSIM4DunitLengthSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.\n");
             }
         }
         T0 = model->BSIM4tcjswg * delTemp;
         if (T0 >= -1.0)
         {   model->BSIM4SunitLengthGateSidewallTempJctCap = model->BSIM4SunitLengthGateSidewallJctCap *(1.0 + T0);
             model->BSIM4DunitLengthGateSidewallTempJctCap = model->BSIM4DunitLengthGateSidewallJctCap *(1.0 + T0);
         }
         else
         {   if (model->BSIM4SunitLengthGateSidewallJctCap > 0.0)
             {   model->BSIM4SunitLengthGateSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.\n");
             }
             if (model->BSIM4DunitLengthGateSidewallJctCap > 0.0)
             {   model->BSIM4DunitLengthGateSidewallTempJctCap = 0.0;
                 fprintf(stderr, "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.\n");
             }
         }

         model->BSIM4PhiBS = model->BSIM4SbulkJctPotential
                           - model->BSIM4tpb * delTemp;
         if (model->BSIM4PhiBS < 0.01)
         {   model->BSIM4PhiBS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.\n");
         }
         model->BSIM4PhiBD = model->BSIM4DbulkJctPotential
                           - model->BSIM4tpb * delTemp;
         if (model->BSIM4PhiBD < 0.01)
         {   model->BSIM4PhiBD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.\n");
         }

         model->BSIM4PhiBSWS = model->BSIM4SsidewallJctPotential
                             - model->BSIM4tpbsw * delTemp;
         if (model->BSIM4PhiBSWS <= 0.01)
         {   model->BSIM4PhiBSWS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.\n");
         }
         model->BSIM4PhiBSWD = model->BSIM4DsidewallJctPotential
                             - model->BSIM4tpbsw * delTemp;
         if (model->BSIM4PhiBSWD <= 0.01)
         {   model->BSIM4PhiBSWD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.\n");
         }

         model->BSIM4PhiBSWGS = model->BSIM4SGatesidewallJctPotential
                              - model->BSIM4tpbswg * delTemp;
         if (model->BSIM4PhiBSWGS <= 0.01)
         {   model->BSIM4PhiBSWGS = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.\n");
         }
         model->BSIM4PhiBSWGD = model->BSIM4DGatesidewallJctPotential
                              - model->BSIM4tpbswg * delTemp;
         if (model->BSIM4PhiBSWGD <= 0.01)
         {   model->BSIM4PhiBSWGD = 0.01;
             fprintf(stderr, "Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.\n");
         } /* End of junction capacitance */


         if (model->BSIM4ijthdfwd <= 0.0)
         {   model->BSIM4ijthdfwd = 0.0;
             fprintf(stderr, "Ijthdfwd reset to %g.\n", model->BSIM4ijthdfwd);
         }
         if (model->BSIM4ijthsfwd <= 0.0)
         {   model->BSIM4ijthsfwd = 0.0;
             fprintf(stderr, "Ijthsfwd reset to %g.\n", model->BSIM4ijthsfwd);
         }
         if (model->BSIM4ijthdrev <= 0.0)
         {   model->BSIM4ijthdrev = 0.0;
             fprintf(stderr, "Ijthdrev reset to %g.\n", model->BSIM4ijthdrev);
         }
         if (model->BSIM4ijthsrev <= 0.0)
         {   model->BSIM4ijthsrev = 0.0;
             fprintf(stderr, "Ijthsrev reset to %g.\n", model->BSIM4ijthsrev);
         }

         if ((model->BSIM4xjbvd <= 0.0) && (model->BSIM4dioMod == 2))
         {   model->BSIM4xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model->BSIM4xjbvd);
         }
         else if ((model->BSIM4xjbvd < 0.0) && (model->BSIM4dioMod == 0))
         {   model->BSIM4xjbvd = 0.0;
             fprintf(stderr, "Xjbvd reset to %g.\n", model->BSIM4xjbvd);
         }

         if (model->BSIM4bvd <= 0.0)   /*4.6.2*/
         {   model->BSIM4bvd = 0.0;
             fprintf(stderr, "BVD reset to %g.\n", model->BSIM4bvd);
         }

         if ((model->BSIM4xjbvs <= 0.0) && (model->BSIM4dioMod == 2))
         {   model->BSIM4xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model->BSIM4xjbvs);
         }
         else if ((model->BSIM4xjbvs < 0.0) && (model->BSIM4dioMod == 0))
         {   model->BSIM4xjbvs = 0.0;
             fprintf(stderr, "Xjbvs reset to %g.\n", model->BSIM4xjbvs);
         }

         if (model->BSIM4bvs <= 0.0)
         {   model->BSIM4bvs = 0.0;
             fprintf(stderr, "BVS reset to %g.\n", model->BSIM4bvs);
         }


         /* loop through all the instances of the model */
         for (here = model->BSIM4instances; here != NULL;
              here = here->BSIM4nextInstance)
         {
              pSizeDependParamKnot = model->pSizeDependParamKnot;
              Size_Not_Found = 1;
              while ((pSizeDependParamKnot != NULL) && Size_Not_Found)
              {   if ((here->BSIM4l == pSizeDependParamKnot->Length)
                      && (here->BSIM4w == pSizeDependParamKnot->Width)
                      && (here->BSIM4nf == pSizeDependParamKnot->NFinger))
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
              Ldrn = here->BSIM4l;
              Wdrn = here->BSIM4w / here->BSIM4nf;

              if (Size_Not_Found)
              {   pParam = TMALLOC(struct bsim4SizeDependParam, 1);
                  if (pLastKnot == NULL)
                      model->pSizeDependParamKnot = pParam;
                  else
                      pLastKnot->pNext = pParam;
                  pParam->pNext = NULL;
                  here->pParam = pParam;

                  pParam->Length = here->BSIM4l;
                  pParam->Width = here->BSIM4w;
                  pParam->NFinger = here->BSIM4nf;
                  Lnew = here->BSIM4l  + model->BSIM4xl ;
                  Wnew = here->BSIM4w / here->BSIM4nf + model->BSIM4xw;

                  T0 = pow(Lnew, model->BSIM4Lln);
                  T1 = pow(Wnew, model->BSIM4Lwn);
                  tmp1 = model->BSIM4Ll / T0 + model->BSIM4Lw / T1
                       + model->BSIM4Lwl / (T0 * T1);
                  pParam->BSIM4dl = model->BSIM4Lint + tmp1;
                  tmp2 = model->BSIM4Llc / T0 + model->BSIM4Lwc / T1
                       + model->BSIM4Lwlc / (T0 * T1);
                  pParam->BSIM4dlc = model->BSIM4dlc + tmp2;

                  T2 = pow(Lnew, model->BSIM4Wln);
                  T3 = pow(Wnew, model->BSIM4Wwn);
                  tmp1 = model->BSIM4Wl / T2 + model->BSIM4Ww / T3
                       + model->BSIM4Wwl / (T2 * T3);
                  pParam->BSIM4dw = model->BSIM4Wint + tmp1;
                  tmp2 = model->BSIM4Wlc / T2 + model->BSIM4Wwc / T3
                       + model->BSIM4Wwlc / (T2 * T3);
                  pParam->BSIM4dwc = model->BSIM4dwc + tmp2;
                  pParam->BSIM4dwj = model->BSIM4dwj + tmp2;

                  pParam->BSIM4leff = Lnew - 2.0 * pParam->BSIM4dl;
                  if (pParam->BSIM4leff <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4: mosfet %s, model %s: Effective channel length <= 0",
                       model->BSIM4modName, here->BSIM4name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4weff = Wnew - 2.0 * pParam->BSIM4dw;
                  if (pParam->BSIM4weff <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4: mosfet %s, model %s: Effective channel width <= 0",
                       model->BSIM4modName, here->BSIM4name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4leffCV = Lnew - 2.0 * pParam->BSIM4dlc;
                  if (pParam->BSIM4leffCV <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4: mosfet %s, model %s: Effective channel length for C-V <= 0",
                       model->BSIM4modName, here->BSIM4name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4weffCV = Wnew - 2.0 * pParam->BSIM4dwc;
                  if (pParam->BSIM4weffCV <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4: mosfet %s, model %s: Effective channel width for C-V <= 0",
                       model->BSIM4modName, here->BSIM4name);
                      return(E_BADPARM);
                  }

                  pParam->BSIM4weffCJ = Wnew - 2.0 * pParam->BSIM4dwj;
                  if (pParam->BSIM4weffCJ <= 0.0)
                  {
                      SPfrontEnd->IFerrorf(ERR_FATAL,
                      "BSIM4: mosfet %s, model %s: Effective channel width for S/D junctions <= 0",
                       model->BSIM4modName, here->BSIM4name);
                      return(E_BADPARM);
                  }


                  if (model->BSIM4binUnit == 1)
                  {   Inv_L = 1.0e-6 / pParam->BSIM4leff;
                      Inv_W = 1.0e-6 / pParam->BSIM4weff;
                      Inv_LW = 1.0e-12 / (pParam->BSIM4leff
                             * pParam->BSIM4weff);
                  }
                  else
                  {   Inv_L = 1.0 / pParam->BSIM4leff;
                      Inv_W = 1.0 / pParam->BSIM4weff;
                      Inv_LW = 1.0 / (pParam->BSIM4leff
                             * pParam->BSIM4weff);
                  }
                  pParam->BSIM4cdsc = model->BSIM4cdsc
                                    + model->BSIM4lcdsc * Inv_L
                                    + model->BSIM4wcdsc * Inv_W
                                    + model->BSIM4pcdsc * Inv_LW;
                  pParam->BSIM4cdscb = model->BSIM4cdscb
                                     + model->BSIM4lcdscb * Inv_L
                                     + model->BSIM4wcdscb * Inv_W
                                     + model->BSIM4pcdscb * Inv_LW;

                      pParam->BSIM4cdscd = model->BSIM4cdscd
                                     + model->BSIM4lcdscd * Inv_L
                                     + model->BSIM4wcdscd * Inv_W
                                     + model->BSIM4pcdscd * Inv_LW;

                  pParam->BSIM4cit = model->BSIM4cit
                                   + model->BSIM4lcit * Inv_L
                                   + model->BSIM4wcit * Inv_W
                                   + model->BSIM4pcit * Inv_LW;
                  pParam->BSIM4nfactor = model->BSIM4nfactor
                                       + model->BSIM4lnfactor * Inv_L
                                       + model->BSIM4wnfactor * Inv_W
                                       + model->BSIM4pnfactor * Inv_LW;
                  pParam->BSIM4tnfactor = model->BSIM4tnfactor                        /* v4.7 */
                                       + model->BSIM4ltnfactor * Inv_L
                                       + model->BSIM4wtnfactor * Inv_W
                                       + model->BSIM4ptnfactor * Inv_LW;
                  pParam->BSIM4xj = model->BSIM4xj
                                  + model->BSIM4lxj * Inv_L
                                  + model->BSIM4wxj * Inv_W
                                  + model->BSIM4pxj * Inv_LW;
                  pParam->BSIM4vsat = model->BSIM4vsat
                                    + model->BSIM4lvsat * Inv_L
                                    + model->BSIM4wvsat * Inv_W
                                    + model->BSIM4pvsat * Inv_LW;
                  pParam->BSIM4at = model->BSIM4at
                                  + model->BSIM4lat * Inv_L
                                  + model->BSIM4wat * Inv_W
                                  + model->BSIM4pat * Inv_LW;
                  pParam->BSIM4a0 = model->BSIM4a0
                                  + model->BSIM4la0 * Inv_L
                                  + model->BSIM4wa0 * Inv_W
                                  + model->BSIM4pa0 * Inv_LW;

                  pParam->BSIM4ags = model->BSIM4ags
                                  + model->BSIM4lags * Inv_L
                                  + model->BSIM4wags * Inv_W
                                  + model->BSIM4pags * Inv_LW;

                  pParam->BSIM4a1 = model->BSIM4a1
                                  + model->BSIM4la1 * Inv_L
                                  + model->BSIM4wa1 * Inv_W
                                  + model->BSIM4pa1 * Inv_LW;
                  pParam->BSIM4a2 = model->BSIM4a2
                                  + model->BSIM4la2 * Inv_L
                                  + model->BSIM4wa2 * Inv_W
                                  + model->BSIM4pa2 * Inv_LW;
                  pParam->BSIM4keta = model->BSIM4keta
                                    + model->BSIM4lketa * Inv_L
                                    + model->BSIM4wketa * Inv_W
                                    + model->BSIM4pketa * Inv_LW;
                  pParam->BSIM4nsub = model->BSIM4nsub
                                    + model->BSIM4lnsub * Inv_L
                                    + model->BSIM4wnsub * Inv_W
                                    + model->BSIM4pnsub * Inv_LW;
                  pParam->BSIM4ndep = model->BSIM4ndep
                                    + model->BSIM4lndep * Inv_L
                                    + model->BSIM4wndep * Inv_W
                                    + model->BSIM4pndep * Inv_LW;
                  pParam->BSIM4nsd = model->BSIM4nsd
                                   + model->BSIM4lnsd * Inv_L
                                   + model->BSIM4wnsd * Inv_W
                                   + model->BSIM4pnsd * Inv_LW;
                  pParam->BSIM4phin = model->BSIM4phin
                                    + model->BSIM4lphin * Inv_L
                                    + model->BSIM4wphin * Inv_W
                                    + model->BSIM4pphin * Inv_LW;
                  pParam->BSIM4ngate = model->BSIM4ngate
                                     + model->BSIM4lngate * Inv_L
                                     + model->BSIM4wngate * Inv_W
                                     + model->BSIM4pngate * Inv_LW;
                  pParam->BSIM4gamma1 = model->BSIM4gamma1
                                      + model->BSIM4lgamma1 * Inv_L
                                      + model->BSIM4wgamma1 * Inv_W
                                      + model->BSIM4pgamma1 * Inv_LW;
                  pParam->BSIM4gamma2 = model->BSIM4gamma2
                                      + model->BSIM4lgamma2 * Inv_L
                                      + model->BSIM4wgamma2 * Inv_W
                                      + model->BSIM4pgamma2 * Inv_LW;
                  pParam->BSIM4vbx = model->BSIM4vbx
                                   + model->BSIM4lvbx * Inv_L
                                   + model->BSIM4wvbx * Inv_W
                                   + model->BSIM4pvbx * Inv_LW;
                  pParam->BSIM4vbm = model->BSIM4vbm
                                   + model->BSIM4lvbm * Inv_L
                                   + model->BSIM4wvbm * Inv_W
                                   + model->BSIM4pvbm * Inv_LW;
                  pParam->BSIM4xt = model->BSIM4xt
                                   + model->BSIM4lxt * Inv_L
                                   + model->BSIM4wxt * Inv_W
                                   + model->BSIM4pxt * Inv_LW;
                  pParam->BSIM4vfb = model->BSIM4vfb
                                   + model->BSIM4lvfb * Inv_L
                                   + model->BSIM4wvfb * Inv_W
                                   + model->BSIM4pvfb * Inv_LW;
                  pParam->BSIM4k1 = model->BSIM4k1
                                  + model->BSIM4lk1 * Inv_L
                                  + model->BSIM4wk1 * Inv_W
                                  + model->BSIM4pk1 * Inv_LW;
                  pParam->BSIM4kt1 = model->BSIM4kt1
                                   + model->BSIM4lkt1 * Inv_L
                                   + model->BSIM4wkt1 * Inv_W
                                   + model->BSIM4pkt1 * Inv_LW;
                  pParam->BSIM4kt1l = model->BSIM4kt1l
                                    + model->BSIM4lkt1l * Inv_L
                                    + model->BSIM4wkt1l * Inv_W
                                    + model->BSIM4pkt1l * Inv_LW;
                  pParam->BSIM4k2 = model->BSIM4k2
                                  + model->BSIM4lk2 * Inv_L
                                  + model->BSIM4wk2 * Inv_W
                                  + model->BSIM4pk2 * Inv_LW;
                  pParam->BSIM4kt2 = model->BSIM4kt2
                                   + model->BSIM4lkt2 * Inv_L
                                   + model->BSIM4wkt2 * Inv_W
                                   + model->BSIM4pkt2 * Inv_LW;
                  pParam->BSIM4k3 = model->BSIM4k3
                                  + model->BSIM4lk3 * Inv_L
                                  + model->BSIM4wk3 * Inv_W
                                  + model->BSIM4pk3 * Inv_LW;
                  pParam->BSIM4k3b = model->BSIM4k3b
                                   + model->BSIM4lk3b * Inv_L
                                   + model->BSIM4wk3b * Inv_W
                                   + model->BSIM4pk3b * Inv_LW;
                  pParam->BSIM4w0 = model->BSIM4w0
                                  + model->BSIM4lw0 * Inv_L
                                  + model->BSIM4ww0 * Inv_W
                                  + model->BSIM4pw0 * Inv_LW;
                  pParam->BSIM4lpe0 = model->BSIM4lpe0
                                    + model->BSIM4llpe0 * Inv_L
                                     + model->BSIM4wlpe0 * Inv_W
                                    + model->BSIM4plpe0 * Inv_LW;
                  pParam->BSIM4lpeb = model->BSIM4lpeb
                                    + model->BSIM4llpeb * Inv_L
                                    + model->BSIM4wlpeb * Inv_W
                                    + model->BSIM4plpeb * Inv_LW;
                  pParam->BSIM4dvtp0 = model->BSIM4dvtp0
                                     + model->BSIM4ldvtp0 * Inv_L
                                     + model->BSIM4wdvtp0 * Inv_W
                                     + model->BSIM4pdvtp0 * Inv_LW;
                  pParam->BSIM4dvtp1 = model->BSIM4dvtp1
                                     + model->BSIM4ldvtp1 * Inv_L
                                     + model->BSIM4wdvtp1 * Inv_W
                                     + model->BSIM4pdvtp1 * Inv_LW;
                  pParam->BSIM4dvtp2 = model->BSIM4dvtp2                 /* v4.7  */
                                     + model->BSIM4ldvtp2 * Inv_L
                                     + model->BSIM4wdvtp2 * Inv_W
                                     + model->BSIM4pdvtp2 * Inv_LW;
                  pParam->BSIM4dvtp3 = model->BSIM4dvtp3                 /* v4.7  */
                                     + model->BSIM4ldvtp3 * Inv_L
                                     + model->BSIM4wdvtp3 * Inv_W
                                     + model->BSIM4pdvtp3 * Inv_LW;
                  pParam->BSIM4dvtp4 = model->BSIM4dvtp4                 /* v4.7  */
                                     + model->BSIM4ldvtp4 * Inv_L
                                     + model->BSIM4wdvtp4 * Inv_W
                                     + model->BSIM4pdvtp4 * Inv_LW;
                  pParam->BSIM4dvtp5 = model->BSIM4dvtp5                 /* v4.7  */
                                     + model->BSIM4ldvtp5 * Inv_L
                                     + model->BSIM4wdvtp5 * Inv_W
                                     + model->BSIM4pdvtp5 * Inv_LW;
                  pParam->BSIM4dvt0 = model->BSIM4dvt0
                                    + model->BSIM4ldvt0 * Inv_L
                                    + model->BSIM4wdvt0 * Inv_W
                                    + model->BSIM4pdvt0 * Inv_LW;
                  pParam->BSIM4dvt1 = model->BSIM4dvt1
                                    + model->BSIM4ldvt1 * Inv_L
                                    + model->BSIM4wdvt1 * Inv_W
                                    + model->BSIM4pdvt1 * Inv_LW;
                  pParam->BSIM4dvt2 = model->BSIM4dvt2
                                    + model->BSIM4ldvt2 * Inv_L
                                    + model->BSIM4wdvt2 * Inv_W
                                    + model->BSIM4pdvt2 * Inv_LW;
                  pParam->BSIM4dvt0w = model->BSIM4dvt0w
                                    + model->BSIM4ldvt0w * Inv_L
                                    + model->BSIM4wdvt0w * Inv_W
                                    + model->BSIM4pdvt0w * Inv_LW;
                  pParam->BSIM4dvt1w = model->BSIM4dvt1w
                                    + model->BSIM4ldvt1w * Inv_L
                                    + model->BSIM4wdvt1w * Inv_W
                                    + model->BSIM4pdvt1w * Inv_LW;
                  pParam->BSIM4dvt2w = model->BSIM4dvt2w
                                    + model->BSIM4ldvt2w * Inv_L
                                    + model->BSIM4wdvt2w * Inv_W
                                    + model->BSIM4pdvt2w * Inv_LW;
                  pParam->BSIM4drout = model->BSIM4drout
                                     + model->BSIM4ldrout * Inv_L
                                     + model->BSIM4wdrout * Inv_W
                                     + model->BSIM4pdrout * Inv_LW;
                  pParam->BSIM4dsub = model->BSIM4dsub
                                    + model->BSIM4ldsub * Inv_L
                                    + model->BSIM4wdsub * Inv_W
                                    + model->BSIM4pdsub * Inv_LW;
                  pParam->BSIM4vth0 = model->BSIM4vth0
                                    + model->BSIM4lvth0 * Inv_L
                                    + model->BSIM4wvth0 * Inv_W
                                    + model->BSIM4pvth0 * Inv_LW;
                  pParam->BSIM4ua = model->BSIM4ua
                                  + model->BSIM4lua * Inv_L
                                  + model->BSIM4wua * Inv_W
                                  + model->BSIM4pua * Inv_LW;
                  pParam->BSIM4ua1 = model->BSIM4ua1
                                   + model->BSIM4lua1 * Inv_L
                                   + model->BSIM4wua1 * Inv_W
                                   + model->BSIM4pua1 * Inv_LW;
                  pParam->BSIM4ub = model->BSIM4ub
                                  + model->BSIM4lub * Inv_L
                                  + model->BSIM4wub * Inv_W
                                  + model->BSIM4pub * Inv_LW;
                  pParam->BSIM4ub1 = model->BSIM4ub1
                                   + model->BSIM4lub1 * Inv_L
                                   + model->BSIM4wub1 * Inv_W
                                   + model->BSIM4pub1 * Inv_LW;
                  pParam->BSIM4uc = model->BSIM4uc
                                  + model->BSIM4luc * Inv_L
                                  + model->BSIM4wuc * Inv_W
                                  + model->BSIM4puc * Inv_LW;
                  pParam->BSIM4uc1 = model->BSIM4uc1
                                   + model->BSIM4luc1 * Inv_L
                                   + model->BSIM4wuc1 * Inv_W
                                   + model->BSIM4puc1 * Inv_LW;
                  pParam->BSIM4ud = model->BSIM4ud
                                  + model->BSIM4lud * Inv_L
                                  + model->BSIM4wud * Inv_W
                                  + model->BSIM4pud * Inv_LW;
                  pParam->BSIM4ud1 = model->BSIM4ud1
                                  + model->BSIM4lud1 * Inv_L
                                  + model->BSIM4wud1 * Inv_W
                                  + model->BSIM4pud1 * Inv_LW;
                  pParam->BSIM4up = model->BSIM4up
                                  + model->BSIM4lup * Inv_L
                                  + model->BSIM4wup * Inv_W
                                  + model->BSIM4pup * Inv_LW;
                  pParam->BSIM4lp = model->BSIM4lp
                                  + model->BSIM4llp * Inv_L
                                  + model->BSIM4wlp * Inv_W
                                  + model->BSIM4plp * Inv_LW;
                  pParam->BSIM4eu = model->BSIM4eu
                                  + model->BSIM4leu * Inv_L
                                  + model->BSIM4weu * Inv_W
                                  + model->BSIM4peu * Inv_LW;
                  pParam->BSIM4u0 = model->BSIM4u0
                                  + model->BSIM4lu0 * Inv_L
                                  + model->BSIM4wu0 * Inv_W
                                  + model->BSIM4pu0 * Inv_LW;
                  pParam->BSIM4ute = model->BSIM4ute
                                   + model->BSIM4lute * Inv_L
                                   + model->BSIM4wute * Inv_W
                                   + model->BSIM4pute * Inv_LW;
                /*high k mobility*/
                 pParam->BSIM4ucs = model->BSIM4ucs
                                  + model->BSIM4lucs * Inv_L
                                  + model->BSIM4wucs * Inv_W
                                  + model->BSIM4pucs * Inv_LW;
                  pParam->BSIM4ucste = model->BSIM4ucste
                           + model->BSIM4lucste * Inv_L
                                   + model->BSIM4wucste * Inv_W
                                   + model->BSIM4pucste * Inv_LW;

                  pParam->BSIM4voff = model->BSIM4voff
                                    + model->BSIM4lvoff * Inv_L
                                    + model->BSIM4wvoff * Inv_W
                                    + model->BSIM4pvoff * Inv_LW;
                  pParam->BSIM4tvoff = model->BSIM4tvoff
                                    + model->BSIM4ltvoff * Inv_L
                                    + model->BSIM4wtvoff * Inv_W
                                    + model->BSIM4ptvoff * Inv_LW;
                  pParam->BSIM4minv = model->BSIM4minv
                                    + model->BSIM4lminv * Inv_L
                                    + model->BSIM4wminv * Inv_W
                                    + model->BSIM4pminv * Inv_LW;
                  pParam->BSIM4minvcv = model->BSIM4minvcv
                                    + model->BSIM4lminvcv * Inv_L
                                    + model->BSIM4wminvcv * Inv_W
                                    + model->BSIM4pminvcv * Inv_LW;
                  pParam->BSIM4fprout = model->BSIM4fprout
                                     + model->BSIM4lfprout * Inv_L
                                     + model->BSIM4wfprout * Inv_W
                                     + model->BSIM4pfprout * Inv_LW;
                  pParam->BSIM4pdits = model->BSIM4pdits
                                     + model->BSIM4lpdits * Inv_L
                                     + model->BSIM4wpdits * Inv_W
                                     + model->BSIM4ppdits * Inv_LW;
                  pParam->BSIM4pditsd = model->BSIM4pditsd
                                      + model->BSIM4lpditsd * Inv_L
                                      + model->BSIM4wpditsd * Inv_W
                                      + model->BSIM4ppditsd * Inv_LW;
                  pParam->BSIM4delta = model->BSIM4delta
                                     + model->BSIM4ldelta * Inv_L
                                     + model->BSIM4wdelta * Inv_W
                                     + model->BSIM4pdelta * Inv_LW;
                  pParam->BSIM4rdsw = model->BSIM4rdsw
                                    + model->BSIM4lrdsw * Inv_L
                                    + model->BSIM4wrdsw * Inv_W
                                    + model->BSIM4prdsw * Inv_LW;
                  pParam->BSIM4rdw = model->BSIM4rdw
                                    + model->BSIM4lrdw * Inv_L
                                    + model->BSIM4wrdw * Inv_W
                                    + model->BSIM4prdw * Inv_LW;
                  pParam->BSIM4rsw = model->BSIM4rsw
                                    + model->BSIM4lrsw * Inv_L
                                    + model->BSIM4wrsw * Inv_W
                                    + model->BSIM4prsw * Inv_LW;
                  pParam->BSIM4prwg = model->BSIM4prwg
                                    + model->BSIM4lprwg * Inv_L
                                    + model->BSIM4wprwg * Inv_W
                                    + model->BSIM4pprwg * Inv_LW;
                  pParam->BSIM4prwb = model->BSIM4prwb
                                    + model->BSIM4lprwb * Inv_L
                                    + model->BSIM4wprwb * Inv_W
                                    + model->BSIM4pprwb * Inv_LW;
                  pParam->BSIM4prt = model->BSIM4prt
                                    + model->BSIM4lprt * Inv_L
                                    + model->BSIM4wprt * Inv_W
                                    + model->BSIM4pprt * Inv_LW;
                  pParam->BSIM4eta0 = model->BSIM4eta0
                                    + model->BSIM4leta0 * Inv_L
                                    + model->BSIM4weta0 * Inv_W
                                    + model->BSIM4peta0 * Inv_LW;
                  pParam->BSIM4teta0 = model->BSIM4teta0                 /* v4.7  */
                                    + model->BSIM4lteta0 * Inv_L
                                    + model->BSIM4wteta0 * Inv_W
                                    + model->BSIM4pteta0 * Inv_LW;
                  pParam->BSIM4tvoffcv = model->BSIM4tvoffcv /* v4.8.0  */
                                    + model->BSIM4ltvoffcv * Inv_L
                                    + model->BSIM4wtvoffcv * Inv_W
                                    + model->BSIM4ptvoffcv * Inv_LW;
                  pParam->BSIM4etab = model->BSIM4etab
                                    + model->BSIM4letab * Inv_L
                                    + model->BSIM4wetab * Inv_W
                                    + model->BSIM4petab * Inv_LW;
                  pParam->BSIM4pclm = model->BSIM4pclm
                                    + model->BSIM4lpclm * Inv_L
                                    + model->BSIM4wpclm * Inv_W
                                    + model->BSIM4ppclm * Inv_LW;
                  pParam->BSIM4pdibl1 = model->BSIM4pdibl1
                                      + model->BSIM4lpdibl1 * Inv_L
                                      + model->BSIM4wpdibl1 * Inv_W
                                      + model->BSIM4ppdibl1 * Inv_LW;
                  pParam->BSIM4pdibl2 = model->BSIM4pdibl2
                                      + model->BSIM4lpdibl2 * Inv_L
                                      + model->BSIM4wpdibl2 * Inv_W
                                      + model->BSIM4ppdibl2 * Inv_LW;
                  pParam->BSIM4pdiblb = model->BSIM4pdiblb
                                      + model->BSIM4lpdiblb * Inv_L
                                      + model->BSIM4wpdiblb * Inv_W
                                      + model->BSIM4ppdiblb * Inv_LW;
                  pParam->BSIM4pscbe1 = model->BSIM4pscbe1
                                      + model->BSIM4lpscbe1 * Inv_L
                                      + model->BSIM4wpscbe1 * Inv_W
                                      + model->BSIM4ppscbe1 * Inv_LW;
                  pParam->BSIM4pscbe2 = model->BSIM4pscbe2
                                      + model->BSIM4lpscbe2 * Inv_L
                                      + model->BSIM4wpscbe2 * Inv_W
                                      + model->BSIM4ppscbe2 * Inv_LW;
                  pParam->BSIM4pvag = model->BSIM4pvag
                                    + model->BSIM4lpvag * Inv_L
                                    + model->BSIM4wpvag * Inv_W
                                    + model->BSIM4ppvag * Inv_LW;
                  pParam->BSIM4wr = model->BSIM4wr
                                  + model->BSIM4lwr * Inv_L
                                  + model->BSIM4wwr * Inv_W
                                  + model->BSIM4pwr * Inv_LW;
                  pParam->BSIM4dwg = model->BSIM4dwg
                                   + model->BSIM4ldwg * Inv_L
                                   + model->BSIM4wdwg * Inv_W
                                   + model->BSIM4pdwg * Inv_LW;
                  pParam->BSIM4dwb = model->BSIM4dwb
                                   + model->BSIM4ldwb * Inv_L
                                   + model->BSIM4wdwb * Inv_W
                                   + model->BSIM4pdwb * Inv_LW;
                  pParam->BSIM4b0 = model->BSIM4b0
                                  + model->BSIM4lb0 * Inv_L
                                  + model->BSIM4wb0 * Inv_W
                                  + model->BSIM4pb0 * Inv_LW;
                  pParam->BSIM4b1 = model->BSIM4b1
                                  + model->BSIM4lb1 * Inv_L
                                  + model->BSIM4wb1 * Inv_W
                                  + model->BSIM4pb1 * Inv_LW;
                  pParam->BSIM4alpha0 = model->BSIM4alpha0
                                      + model->BSIM4lalpha0 * Inv_L
                                      + model->BSIM4walpha0 * Inv_W
                                      + model->BSIM4palpha0 * Inv_LW;
                  pParam->BSIM4alpha1 = model->BSIM4alpha1
                                      + model->BSIM4lalpha1 * Inv_L
                                      + model->BSIM4walpha1 * Inv_W
                                      + model->BSIM4palpha1 * Inv_LW;
                  pParam->BSIM4beta0 = model->BSIM4beta0
                                     + model->BSIM4lbeta0 * Inv_L
                                     + model->BSIM4wbeta0 * Inv_W
                                     + model->BSIM4pbeta0 * Inv_LW;
                  pParam->BSIM4agidl = model->BSIM4agidl
                                     + model->BSIM4lagidl * Inv_L
                                     + model->BSIM4wagidl * Inv_W
                                     + model->BSIM4pagidl * Inv_LW;
                  pParam->BSIM4bgidl = model->BSIM4bgidl
                                     + model->BSIM4lbgidl * Inv_L
                                     + model->BSIM4wbgidl * Inv_W
                                     + model->BSIM4pbgidl * Inv_LW;
                  pParam->BSIM4cgidl = model->BSIM4cgidl
                                     + model->BSIM4lcgidl * Inv_L
                                     + model->BSIM4wcgidl * Inv_W
                                     + model->BSIM4pcgidl * Inv_LW;
                  pParam->BSIM4egidl = model->BSIM4egidl
                                     + model->BSIM4legidl * Inv_L
                                     + model->BSIM4wegidl * Inv_W
                                     + model->BSIM4pegidl * Inv_LW;
                  pParam->BSIM4rgidl = model->BSIM4rgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lrgidl * Inv_L
                                     + model->BSIM4wrgidl * Inv_W
                                     + model->BSIM4prgidl * Inv_LW;
                  pParam->BSIM4kgidl = model->BSIM4kgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lkgidl * Inv_L
                                     + model->BSIM4wkgidl * Inv_W
                                     + model->BSIM4pkgidl * Inv_LW;
                  pParam->BSIM4fgidl = model->BSIM4fgidl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lfgidl * Inv_L
                                     + model->BSIM4wfgidl * Inv_W
                                     + model->BSIM4pfgidl * Inv_LW;
                  pParam->BSIM4agisl = model->BSIM4agisl
                                     + model->BSIM4lagisl * Inv_L
                                     + model->BSIM4wagisl * Inv_W
                                     + model->BSIM4pagisl * Inv_LW;
                  pParam->BSIM4bgisl = model->BSIM4bgisl
                                     + model->BSIM4lbgisl * Inv_L
                                     + model->BSIM4wbgisl * Inv_W
                                     + model->BSIM4pbgisl * Inv_LW;
                  pParam->BSIM4cgisl = model->BSIM4cgisl
                                     + model->BSIM4lcgisl * Inv_L
                                     + model->BSIM4wcgisl * Inv_W
                                     + model->BSIM4pcgisl * Inv_LW;
                  pParam->BSIM4egisl = model->BSIM4egisl
                                     + model->BSIM4legisl * Inv_L
                                     + model->BSIM4wegisl * Inv_W
                                     + model->BSIM4pegisl * Inv_LW;
                   pParam->BSIM4rgisl = model->BSIM4rgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lrgisl * Inv_L
                                     + model->BSIM4wrgisl * Inv_W
                                     + model->BSIM4prgisl * Inv_LW;
                  pParam->BSIM4kgisl = model->BSIM4kgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lkgisl * Inv_L
                                     + model->BSIM4wkgisl * Inv_W
                                     + model->BSIM4pkgisl * Inv_LW;
                  pParam->BSIM4fgisl = model->BSIM4fgisl                /* v4.7 New GIDL/GISL */
                                     + model->BSIM4lfgisl * Inv_L
                                     + model->BSIM4wfgisl * Inv_W
                                     + model->BSIM4pfgisl * Inv_LW;
                  pParam->BSIM4aigc = model->BSIM4aigc
                                     + model->BSIM4laigc * Inv_L
                                     + model->BSIM4waigc * Inv_W
                                     + model->BSIM4paigc * Inv_LW;
                  pParam->BSIM4bigc = model->BSIM4bigc
                                     + model->BSIM4lbigc * Inv_L
                                     + model->BSIM4wbigc * Inv_W
                                     + model->BSIM4pbigc * Inv_LW;
                  pParam->BSIM4cigc = model->BSIM4cigc
                                     + model->BSIM4lcigc * Inv_L
                                     + model->BSIM4wcigc * Inv_W
                                     + model->BSIM4pcigc * Inv_LW;
                  pParam->BSIM4aigsd = model->BSIM4aigsd
                                     + model->BSIM4laigsd * Inv_L
                                     + model->BSIM4waigsd * Inv_W
                                     + model->BSIM4paigsd * Inv_LW;
                  pParam->BSIM4bigsd = model->BSIM4bigsd
                                     + model->BSIM4lbigsd * Inv_L
                                     + model->BSIM4wbigsd * Inv_W
                                     + model->BSIM4pbigsd * Inv_LW;
                  pParam->BSIM4cigsd = model->BSIM4cigsd
                                     + model->BSIM4lcigsd * Inv_L
                                     + model->BSIM4wcigsd * Inv_W
                                     + model->BSIM4pcigsd * Inv_LW;
                  pParam->BSIM4aigs = model->BSIM4aigs
                                     + model->BSIM4laigs * Inv_L
                                     + model->BSIM4waigs * Inv_W
                                     + model->BSIM4paigs * Inv_LW;
                  pParam->BSIM4bigs = model->BSIM4bigs
                                     + model->BSIM4lbigs * Inv_L
                                     + model->BSIM4wbigs * Inv_W
                                     + model->BSIM4pbigs * Inv_LW;
                  pParam->BSIM4cigs = model->BSIM4cigs
                                     + model->BSIM4lcigs * Inv_L
                                     + model->BSIM4wcigs * Inv_W
                                     + model->BSIM4pcigs * Inv_LW;
                  pParam->BSIM4aigd = model->BSIM4aigd
                                     + model->BSIM4laigd * Inv_L
                                     + model->BSIM4waigd * Inv_W
                                     + model->BSIM4paigd * Inv_LW;
                  pParam->BSIM4bigd = model->BSIM4bigd
                                     + model->BSIM4lbigd * Inv_L
                                     + model->BSIM4wbigd * Inv_W
                                     + model->BSIM4pbigd * Inv_LW;
                  pParam->BSIM4cigd = model->BSIM4cigd
                                     + model->BSIM4lcigd * Inv_L
                                     + model->BSIM4wcigd * Inv_W
                                     + model->BSIM4pcigd * Inv_LW;
                  pParam->BSIM4aigbacc = model->BSIM4aigbacc
                                       + model->BSIM4laigbacc * Inv_L
                                       + model->BSIM4waigbacc * Inv_W
                                       + model->BSIM4paigbacc * Inv_LW;
                  pParam->BSIM4bigbacc = model->BSIM4bigbacc
                                       + model->BSIM4lbigbacc * Inv_L
                                       + model->BSIM4wbigbacc * Inv_W
                                       + model->BSIM4pbigbacc * Inv_LW;
                  pParam->BSIM4cigbacc = model->BSIM4cigbacc
                                       + model->BSIM4lcigbacc * Inv_L
                                       + model->BSIM4wcigbacc * Inv_W
                                       + model->BSIM4pcigbacc * Inv_LW;
                  pParam->BSIM4aigbinv = model->BSIM4aigbinv
                                       + model->BSIM4laigbinv * Inv_L
                                       + model->BSIM4waigbinv * Inv_W
                                       + model->BSIM4paigbinv * Inv_LW;
                  pParam->BSIM4bigbinv = model->BSIM4bigbinv
                                       + model->BSIM4lbigbinv * Inv_L
                                       + model->BSIM4wbigbinv * Inv_W
                                       + model->BSIM4pbigbinv * Inv_LW;
                  pParam->BSIM4cigbinv = model->BSIM4cigbinv
                                       + model->BSIM4lcigbinv * Inv_L
                                       + model->BSIM4wcigbinv * Inv_W
                                       + model->BSIM4pcigbinv * Inv_LW;
                  pParam->BSIM4nigc = model->BSIM4nigc
                                       + model->BSIM4lnigc * Inv_L
                                       + model->BSIM4wnigc * Inv_W
                                       + model->BSIM4pnigc * Inv_LW;
                  pParam->BSIM4nigbacc = model->BSIM4nigbacc
                                       + model->BSIM4lnigbacc * Inv_L
                                       + model->BSIM4wnigbacc * Inv_W
                                       + model->BSIM4pnigbacc * Inv_LW;
                  pParam->BSIM4nigbinv = model->BSIM4nigbinv
                                       + model->BSIM4lnigbinv * Inv_L
                                       + model->BSIM4wnigbinv * Inv_W
                                       + model->BSIM4pnigbinv * Inv_LW;
                  pParam->BSIM4ntox = model->BSIM4ntox
                                    + model->BSIM4lntox * Inv_L
                                    + model->BSIM4wntox * Inv_W
                                    + model->BSIM4pntox * Inv_LW;
                  pParam->BSIM4eigbinv = model->BSIM4eigbinv
                                       + model->BSIM4leigbinv * Inv_L
                                       + model->BSIM4weigbinv * Inv_W
                                       + model->BSIM4peigbinv * Inv_LW;
                  pParam->BSIM4pigcd = model->BSIM4pigcd
                                     + model->BSIM4lpigcd * Inv_L
                                     + model->BSIM4wpigcd * Inv_W
                                     + model->BSIM4ppigcd * Inv_LW;
                  pParam->BSIM4poxedge = model->BSIM4poxedge
                                       + model->BSIM4lpoxedge * Inv_L
                                       + model->BSIM4wpoxedge * Inv_W
                                       + model->BSIM4ppoxedge * Inv_LW;
                  pParam->BSIM4xrcrg1 = model->BSIM4xrcrg1
                                      + model->BSIM4lxrcrg1 * Inv_L
                                      + model->BSIM4wxrcrg1 * Inv_W
                                      + model->BSIM4pxrcrg1 * Inv_LW;
                  pParam->BSIM4xrcrg2 = model->BSIM4xrcrg2
                                      + model->BSIM4lxrcrg2 * Inv_L
                                      + model->BSIM4wxrcrg2 * Inv_W
                                      + model->BSIM4pxrcrg2 * Inv_LW;
                  pParam->BSIM4lambda = model->BSIM4lambda
                                      + model->BSIM4llambda * Inv_L
                                      + model->BSIM4wlambda * Inv_W
                                      + model->BSIM4plambda * Inv_LW;
                  pParam->BSIM4vtl = model->BSIM4vtl
                                      + model->BSIM4lvtl * Inv_L
                                      + model->BSIM4wvtl * Inv_W
                                      + model->BSIM4pvtl * Inv_LW;
                  pParam->BSIM4xn = model->BSIM4xn
                                      + model->BSIM4lxn * Inv_L
                                      + model->BSIM4wxn * Inv_W
                                      + model->BSIM4pxn * Inv_LW;
                  pParam->BSIM4vfbsdoff = model->BSIM4vfbsdoff
                                      + model->BSIM4lvfbsdoff * Inv_L
                                      + model->BSIM4wvfbsdoff * Inv_W
                                      + model->BSIM4pvfbsdoff * Inv_LW;
                  pParam->BSIM4tvfbsdoff = model->BSIM4tvfbsdoff
                                      + model->BSIM4ltvfbsdoff * Inv_L
                                      + model->BSIM4wtvfbsdoff * Inv_W
                                      + model->BSIM4ptvfbsdoff * Inv_LW;

                  pParam->BSIM4cgsl = model->BSIM4cgsl
                                    + model->BSIM4lcgsl * Inv_L
                                    + model->BSIM4wcgsl * Inv_W
                                    + model->BSIM4pcgsl * Inv_LW;
                  pParam->BSIM4cgdl = model->BSIM4cgdl
                                    + model->BSIM4lcgdl * Inv_L
                                    + model->BSIM4wcgdl * Inv_W
                                    + model->BSIM4pcgdl * Inv_LW;
                  pParam->BSIM4ckappas = model->BSIM4ckappas
                                       + model->BSIM4lckappas * Inv_L
                                       + model->BSIM4wckappas * Inv_W
                                        + model->BSIM4pckappas * Inv_LW;
                  pParam->BSIM4ckappad = model->BSIM4ckappad
                                       + model->BSIM4lckappad * Inv_L
                                       + model->BSIM4wckappad * Inv_W
                                       + model->BSIM4pckappad * Inv_LW;
                  pParam->BSIM4cf = model->BSIM4cf
                                  + model->BSIM4lcf * Inv_L
                                  + model->BSIM4wcf * Inv_W
                                  + model->BSIM4pcf * Inv_LW;
                  pParam->BSIM4clc = model->BSIM4clc
                                   + model->BSIM4lclc * Inv_L
                                   + model->BSIM4wclc * Inv_W
                                   + model->BSIM4pclc * Inv_LW;
                  pParam->BSIM4cle = model->BSIM4cle
                                   + model->BSIM4lcle * Inv_L
                                   + model->BSIM4wcle * Inv_W
                                   + model->BSIM4pcle * Inv_LW;
                  pParam->BSIM4vfbcv = model->BSIM4vfbcv
                                     + model->BSIM4lvfbcv * Inv_L
                                     + model->BSIM4wvfbcv * Inv_W
                                     + model->BSIM4pvfbcv * Inv_LW;
                  pParam->BSIM4acde = model->BSIM4acde
                                    + model->BSIM4lacde * Inv_L
                                    + model->BSIM4wacde * Inv_W
                                    + model->BSIM4pacde * Inv_LW;
                  pParam->BSIM4moin = model->BSIM4moin
                                    + model->BSIM4lmoin * Inv_L
                                    + model->BSIM4wmoin * Inv_W
                                    + model->BSIM4pmoin * Inv_LW;
                  pParam->BSIM4noff = model->BSIM4noff
                                    + model->BSIM4lnoff * Inv_L
                                    + model->BSIM4wnoff * Inv_W
                                    + model->BSIM4pnoff * Inv_LW;
                  pParam->BSIM4voffcv = model->BSIM4voffcv
                                      + model->BSIM4lvoffcv * Inv_L
                                      + model->BSIM4wvoffcv * Inv_W
                                      + model->BSIM4pvoffcv * Inv_LW;
                  pParam->BSIM4kvth0we = model->BSIM4kvth0we
                                      + model->BSIM4lkvth0we * Inv_L
                                      + model->BSIM4wkvth0we * Inv_W
                                      + model->BSIM4pkvth0we * Inv_LW;
                  pParam->BSIM4k2we = model->BSIM4k2we
                                      + model->BSIM4lk2we * Inv_L
                                      + model->BSIM4wk2we * Inv_W
                                      + model->BSIM4pk2we * Inv_LW;
                  pParam->BSIM4ku0we = model->BSIM4ku0we
                                      + model->BSIM4lku0we * Inv_L
                                      + model->BSIM4wku0we * Inv_W
                                      + model->BSIM4pku0we * Inv_LW;

                  pParam->BSIM4abulkCVfactor = 1.0 + pow((pParam->BSIM4clc
                                             / pParam->BSIM4leffCV),
                                             pParam->BSIM4cle);

                  T0 = (TRatio - 1.0);

                  PowWeffWr = pow(pParam->BSIM4weffCJ * 1.0e6, pParam->BSIM4wr) * here->BSIM4nf;

                  T1 = T2 = T3 = T4 = 0.0;
                  pParam->BSIM4ucs = pParam->BSIM4ucs * pow(TRatio, pParam->BSIM4ucste);
                  if(model->BSIM4tempMod == 0)
                  {
                      pParam->BSIM4ua = pParam->BSIM4ua + pParam->BSIM4ua1 * T0;
                      pParam->BSIM4ub = pParam->BSIM4ub + pParam->BSIM4ub1 * T0;
                      pParam->BSIM4uc = pParam->BSIM4uc + pParam->BSIM4uc1 * T0;
                      pParam->BSIM4ud = pParam->BSIM4ud + pParam->BSIM4ud1 * T0;
                      pParam->BSIM4vsattemp = pParam->BSIM4vsat - pParam->BSIM4at * T0;
                      T10 = pParam->BSIM4prt * T0;
                      if(model->BSIM4rdsMod)
                      {
                            /* External Rd(V) */
                            T1 = pParam->BSIM4rdw + T10;
                            T2 = model->BSIM4rdwmin + T10;
                            /* External Rs(V) */
                            T3 = pParam->BSIM4rsw + T10;
                            T4 = model->BSIM4rswmin + T10;
                      }
                      /* Internal Rds(V) in IV */
                      pParam->BSIM4rds0 = (pParam->BSIM4rdsw + T10)
                                            * here->BSIM4nf / PowWeffWr;
                      pParam->BSIM4rdswmin = (model->BSIM4rdswmin + T10)
                                               * here->BSIM4nf / PowWeffWr;
                  }
                  else
                  {
                      if (model->BSIM4tempMod == 3)
                      {
                          pParam->BSIM4ua = pParam->BSIM4ua * pow(TRatio, pParam->BSIM4ua1) ;
                            pParam->BSIM4ub = pParam->BSIM4ub * pow(TRatio, pParam->BSIM4ub1);
                            pParam->BSIM4uc = pParam->BSIM4uc * pow(TRatio, pParam->BSIM4uc1);
                            pParam->BSIM4ud = pParam->BSIM4ud * pow(TRatio, pParam->BSIM4ud1);
                      }
                      else
                      {
                          /* tempMod = 1, 2 */
                            pParam->BSIM4ua = pParam->BSIM4ua * (1.0 + pParam->BSIM4ua1 * delTemp) ;
                            pParam->BSIM4ub = pParam->BSIM4ub * (1.0 + pParam->BSIM4ub1 * delTemp);
                            pParam->BSIM4uc = pParam->BSIM4uc * (1.0 + pParam->BSIM4uc1 * delTemp);
                            pParam->BSIM4ud = pParam->BSIM4ud * (1.0 + pParam->BSIM4ud1 * delTemp);
                      }
                      pParam->BSIM4vsattemp = pParam->BSIM4vsat * (1.0 - pParam->BSIM4at * delTemp);
                      T10 = 1.0 + pParam->BSIM4prt * delTemp;
                      if(model->BSIM4rdsMod)
                      {
                            /* External Rd(V) */
                            T1 = pParam->BSIM4rdw * T10;
                            T2 = model->BSIM4rdwmin * T10;

                            /* External Rs(V) */
                            T3 = pParam->BSIM4rsw * T10;
                            T4 = model->BSIM4rswmin * T10;
                      }
                      /* Internal Rds(V) in IV */
                      pParam->BSIM4rds0 = pParam->BSIM4rdsw * T10 * here->BSIM4nf / PowWeffWr;
                      pParam->BSIM4rdswmin = model->BSIM4rdswmin * T10 * here->BSIM4nf / PowWeffWr;
                  }

                  if (T1 < 0.0)
                  {   T1 = 0.0;
                      printf("Warning: Rdw at current temperature is negative; set to 0.\n");
                  }
                  if (T2 < 0.0)
                  {   T2 = 0.0;
                      printf("Warning: Rdwmin at current temperature is negative; set to 0.\n");
                  }
                  pParam->BSIM4rd0 = T1 / PowWeffWr;
                  pParam->BSIM4rdwmin = T2 / PowWeffWr;
                  if (T3 < 0.0)
                  {   T3 = 0.0;
                      printf("Warning: Rsw at current temperature is negative; set to 0.\n");
                  }
                  if (T4 < 0.0)
                  {   T4 = 0.0;
                      printf("Warning: Rswmin at current temperature is negative; set to 0.\n");
                  }
                  pParam->BSIM4rs0 = T3 / PowWeffWr;
                  pParam->BSIM4rswmin = T4 / PowWeffWr;

                  if (pParam->BSIM4u0 > 1.0)
                      pParam->BSIM4u0 = pParam->BSIM4u0 / 1.0e4;

                  /* mobility channel length dependence */
                  T5 = 1.0 - pParam->BSIM4up * exp( - pParam->BSIM4leff / pParam->BSIM4lp);
                  pParam->BSIM4u0temp = pParam->BSIM4u0 * T5
                                      * pow(TRatio, pParam->BSIM4ute);
                  if (pParam->BSIM4eu < 0.0)
                  {   pParam->BSIM4eu = 0.0;
                      printf("Warning: eu has been negative; reset to 0.0.\n");
                  }
                  if (pParam->BSIM4ucs < 0.0)
                  {   pParam->BSIM4ucs = 0.0;
                      printf("Warning: ucs has been negative; reset to 0.0.\n");
                  }

                 pParam->BSIM4vfbsdoff = pParam->BSIM4vfbsdoff * (1.0 + pParam->BSIM4tvfbsdoff * delTemp);
                 pParam->BSIM4voff = pParam->BSIM4voff * (1.0 + pParam->BSIM4tvoff * delTemp);

                 pParam->BSIM4nfactor = pParam->BSIM4nfactor + pParam->BSIM4tnfactor * delTemp / Tnom;  /* v4.7 temp dep of leakage currents */
                        pParam->BSIM4voffcv = pParam->BSIM4voffcv * (1.0 + pParam->BSIM4tvoffcv * delTemp);   /*         v4.7 temp dep of leakage currents */
                 pParam->BSIM4eta0 = pParam->BSIM4eta0 + pParam->BSIM4teta0 * delTemp / Tnom;   /*         v4.7 temp dep of leakage currents */

                /* Source End Velocity Limit  */
                        if((model->BSIM4vtlGiven) && (model->BSIM4vtl > 0.0) )
                      {
                     if(model->BSIM4lc < 0.0) pParam->BSIM4lc = 0.0;
                     else   pParam->BSIM4lc = model->BSIM4lc ;
                     T0 = pParam->BSIM4leff / (pParam->BSIM4xn * pParam->BSIM4leff + pParam->BSIM4lc);
                     pParam->BSIM4tfactor = (1.0 - T0) / (1.0 + T0 );
                       }

                  pParam->BSIM4cgdo = (model->BSIM4cgdo + pParam->BSIM4cf)
                                    * pParam->BSIM4weffCV;
                  pParam->BSIM4cgso = (model->BSIM4cgso + pParam->BSIM4cf)
                                    * pParam->BSIM4weffCV;
                  pParam->BSIM4cgbo = model->BSIM4cgbo * pParam->BSIM4leffCV * here->BSIM4nf;

                  if (!model->BSIM4ndepGiven && model->BSIM4gamma1Given)
                  {   T0 = pParam->BSIM4gamma1 * model->BSIM4coxe;
                      pParam->BSIM4ndep = 3.01248e22 * T0 * T0;
                  }

                  pParam->BSIM4phi = Vtm0 * log(pParam->BSIM4ndep / ni)
                                   + pParam->BSIM4phin + 0.4;

                  pParam->BSIM4sqrtPhi = sqrt(pParam->BSIM4phi);
                  pParam->BSIM4phis3 = pParam->BSIM4sqrtPhi * pParam->BSIM4phi;

                  pParam->BSIM4Xdep0 = sqrt(2.0 * epssub / (Charge_q
                                     * pParam->BSIM4ndep * 1.0e6))
                                     * pParam->BSIM4sqrtPhi;
                  pParam->BSIM4sqrtXdep0 = sqrt(pParam->BSIM4Xdep0);

                  if(model->BSIM4mtrlMod == 0)
                    pParam->BSIM4litl = sqrt(3.0 * 3.9 / epsrox * pParam->BSIM4xj * toxe);
                  else
                    pParam->BSIM4litl = sqrt(model->BSIM4epsrsub/epsrox * pParam->BSIM4xj * toxe);

                  pParam->BSIM4vbi = Vtm0 * log(pParam->BSIM4nsd
                                   * pParam->BSIM4ndep / (ni * ni));

                  if (model->BSIM4mtrlMod == 0)
                  {
                    if (pParam->BSIM4ngate > 0.0)
                    {   pParam->BSIM4vfbsd = Vtm0 * log(pParam->BSIM4ngate
                                         / pParam->BSIM4nsd);
                     }
                    else
                      pParam->BSIM4vfbsd = 0.0;
                  }
                  else
                  {
                    T0 = Vtm0 * log(pParam->BSIM4nsd/ni);
                    T1 = 0.5 * Eg0;
                    if(T0 > T1)
                      T0 = T1;
                    T2 = model->BSIM4easub + T1 - model->BSIM4type * T0;
                    pParam->BSIM4vfbsd = model->BSIM4phig - T2;
                  }

                  pParam->BSIM4cdep0 = sqrt(Charge_q * epssub
                                     * pParam->BSIM4ndep * 1.0e6 / 2.0
                                     / pParam->BSIM4phi);

                  pParam->BSIM4ToxRatio = exp(pParam->BSIM4ntox
                                        * log(model->BSIM4toxref / toxe))
                                        / toxe / toxe;
                  pParam->BSIM4ToxRatioEdge = exp(pParam->BSIM4ntox
                                            * log(model->BSIM4toxref
                                            / (toxe * pParam->BSIM4poxedge)))
                                            / toxe / toxe
                                            / pParam->BSIM4poxedge / pParam->BSIM4poxedge;
                  pParam->BSIM4Aechvb = (model->BSIM4type == NMOS) ? 4.97232e-7 : 3.42537e-7;
                  pParam->BSIM4Bechvb = (model->BSIM4type == NMOS) ? 7.45669e11 : 1.16645e12;
                  pParam->BSIM4AechvbEdgeS = pParam->BSIM4Aechvb * pParam->BSIM4weff
                                          * model->BSIM4dlcig * pParam->BSIM4ToxRatioEdge;
                  pParam->BSIM4AechvbEdgeD = pParam->BSIM4Aechvb * pParam->BSIM4weff
                                          * model->BSIM4dlcigd * pParam->BSIM4ToxRatioEdge;
                  pParam->BSIM4BechvbEdge = -pParam->BSIM4Bechvb
                                          * toxe * pParam->BSIM4poxedge;
                  pParam->BSIM4Aechvb *= pParam->BSIM4weff * pParam->BSIM4leff
                                       * pParam->BSIM4ToxRatio;
                  pParam->BSIM4Bechvb *= -toxe;


                  pParam->BSIM4mstar = 0.5 + atan(pParam->BSIM4minv) / PI;
                  pParam->BSIM4mstarcv = 0.5 + atan(pParam->BSIM4minvcv) / PI;
                  pParam->BSIM4voffcbn =  pParam->BSIM4voff + model->BSIM4voffl / pParam->BSIM4leff;
                  pParam->BSIM4voffcbncv =  pParam->BSIM4voffcv + model->BSIM4voffcvl / pParam->BSIM4leff;

                  pParam->BSIM4ldeb = sqrt(epssub * Vtm0 / (Charge_q
                                    * pParam->BSIM4ndep * 1.0e6)) / 3.0;
                  pParam->BSIM4acde *= pow((pParam->BSIM4ndep / 2.0e16), -0.25);


                  if (model->BSIM4k1Given || model->BSIM4k2Given)
                  {   if (!model->BSIM4k1Given)
                      {
                          if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) /* don't print in sensitivity */
                              fprintf(stdout, "Warning: k1 should be specified with k2.\n");
                          pParam->BSIM4k1 = 0.53;
                      }
                      if (!model->BSIM4k2Given)
                      {
                          if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) /* don't print in sensitivity */
                              fprintf(stdout, "Warning: k2 should be specified with k1.\n");
                          pParam->BSIM4k2 = -0.0186;
                      }
                      if ((!ckt->CKTcurJob) || (ckt->CKTcurJob->JOBtype < 9)) { /* don't print in sensitivity */
                          if (model->BSIM4nsubGiven)
                              fprintf(stdout, "Warning: nsub is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4xtGiven)
                              fprintf(stdout, "Warning: xt is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4vbxGiven)
                              fprintf(stdout, "Warning: vbx is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4gamma1Given)
                              fprintf(stdout, "Warning: gamma1 is ignored because k1 or k2 is given.\n");
                          if (model->BSIM4gamma2Given)
                              fprintf(stdout, "Warning: gamma2 is ignored because k1 or k2 is given.\n");
                      }
                  }
                  else
                  {   if (!model->BSIM4vbxGiven)
                          pParam->BSIM4vbx = pParam->BSIM4phi - 7.7348e-4
                                           * pParam->BSIM4ndep
                                           * pParam->BSIM4xt * pParam->BSIM4xt;
                      if (pParam->BSIM4vbx > 0.0)
                          pParam->BSIM4vbx = -pParam->BSIM4vbx;
                      if (pParam->BSIM4vbm > 0.0)
                          pParam->BSIM4vbm = -pParam->BSIM4vbm;

                      if (!model->BSIM4gamma1Given)
                          pParam->BSIM4gamma1 = 5.753e-12
                                              * sqrt(pParam->BSIM4ndep)
                                              / model->BSIM4coxe;
                      if (!model->BSIM4gamma2Given)
                          pParam->BSIM4gamma2 = 5.753e-12
                                              * sqrt(pParam->BSIM4nsub)
                                              / model->BSIM4coxe;

                      T0 = pParam->BSIM4gamma1 - pParam->BSIM4gamma2;
                      T1 = sqrt(pParam->BSIM4phi - pParam->BSIM4vbx)
                         - pParam->BSIM4sqrtPhi;
                      T2 = sqrt(pParam->BSIM4phi * (pParam->BSIM4phi
                         - pParam->BSIM4vbm)) - pParam->BSIM4phi;
                      pParam->BSIM4k2 = T0 * T1 / (2.0 * T2 + pParam->BSIM4vbm);
                      pParam->BSIM4k1 = pParam->BSIM4gamma2 - 2.0
                                      * pParam->BSIM4k2 * sqrt(pParam->BSIM4phi
                                      - pParam->BSIM4vbm);
                  }

                  if (!model->BSIM4vfbGiven)
                  {
                    if (model->BSIM4vth0Given)
                      {   pParam->BSIM4vfb = model->BSIM4type * pParam->BSIM4vth0
                                           - pParam->BSIM4phi - pParam->BSIM4k1
                                           * pParam->BSIM4sqrtPhi;
                      }
                      else
                      {
                        if ((model->BSIM4mtrlMod) && (model->BSIM4phigGiven) &&
                            (model->BSIM4nsubGiven))
                          {
                            T0 = Vtm0 * log(pParam->BSIM4nsub/ni);
                            T1 = 0.5 * Eg0;
                            if(T0 > T1)
                              T0 = T1;
                            T2 = model->BSIM4easub + T1 + model->BSIM4type * T0;
                            pParam->BSIM4vfb = model->BSIM4phig - T2;
                          }
                        else
                          {
                            pParam->BSIM4vfb = -1.0;
                          }
                      }
                  }
                   if (!model->BSIM4vth0Given)
                  {   pParam->BSIM4vth0 = model->BSIM4type * (pParam->BSIM4vfb
                                        + pParam->BSIM4phi + pParam->BSIM4k1
                                        * pParam->BSIM4sqrtPhi);
                  }

                  pParam->BSIM4k1ox = pParam->BSIM4k1 * toxe
                                    / model->BSIM4toxm;

                  tmp = sqrt(epssub / (epsrox * EPS0) * toxe * pParam->BSIM4Xdep0);
                    T0 = pParam->BSIM4dsub * pParam->BSIM4leff / tmp;
                  if (T0 < EXP_THRESHOLD)
                    {   T1 = exp(T0);
                            T2 = T1 - 1.0;
                            T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      pParam->BSIM4theta0vb0 = T1 / T4;
                  }
                  else
                      pParam->BSIM4theta0vb0 = 1.0 / (MAX_EXP - 2.0);

                   T0 = pParam->BSIM4drout * pParam->BSIM4leff / tmp;
                  if (T0 < EXP_THRESHOLD)
                         {   T1 = exp(T0);
                            T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T5 = T1 / T4;
                  }
                  else
                      T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                  pParam->BSIM4thetaRout = pParam->BSIM4pdibl1 * T5
                                         + pParam->BSIM4pdibl2;

                  tmp = sqrt(pParam->BSIM4Xdep0);
                  tmp1 = pParam->BSIM4vbi - pParam->BSIM4phi;
                  tmp2 = model->BSIM4factor1 * tmp;

                  T0 = pParam->BSIM4dvt1w * pParam->BSIM4weff
                     * pParam->BSIM4leff / tmp2;
                  if (T0 < EXP_THRESHOLD)
                  {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T8 = T1 / T4;
                  }
                  else
                      T8 = 1.0 / (MAX_EXP - 2.0);
                  T0 = pParam->BSIM4dvt0w * T8;
                  T8 = T0 * tmp1;

                  T0 = pParam->BSIM4dvt1 * pParam->BSIM4leff / tmp2;
                  if (T0 < EXP_THRESHOLD)
                  {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T9 = T1 / T4;
                  }
                  else
                      T9 = 1.0 / (MAX_EXP - 2.0);
                  T9 = pParam->BSIM4dvt0 * T9 * tmp1;

                  T4 = toxe * pParam->BSIM4phi
                     / (pParam->BSIM4weff + pParam->BSIM4w0);

                  T0 = sqrt(1.0 + pParam->BSIM4lpe0 / pParam->BSIM4leff);
                  if((model->BSIM4tempMod == 1) || (model->BSIM4tempMod == 0))
                          T3 = (pParam->BSIM4kt1 + pParam->BSIM4kt1l / pParam->BSIM4leff)
                                     * (TRatio - 1.0);
                  if((model->BSIM4tempMod == 2)||(model->BSIM4tempMod == 3))
                        T3 = - pParam->BSIM4kt1 * (TRatio - 1.0);

                  T5 = pParam->BSIM4k1ox * (T0 - 1.0) * pParam->BSIM4sqrtPhi
                     + T3;
                  pParam->BSIM4vfbzbfactor = - T8 - T9 + pParam->BSIM4k3 * T4 + T5
                                             - pParam->BSIM4phi - pParam->BSIM4k1 * pParam->BSIM4sqrtPhi;

                  /* stress effect */

                      wlod = model->BSIM4wlod;
                      if (model->BSIM4wlod < 0.0)
                  {   fprintf(stderr, "Warning: WLOD = %g is less than 0. 0.0 is used\n",model->BSIM4wlod);
                             wlod = 0.0;
                  }
                  T0 = pow(Lnew, model->BSIM4llodku0);
                  W_tmp = Wnew + wlod;
                  T1 = pow(W_tmp, model->BSIM4wlodku0);
                  tmp1 = model->BSIM4lku0 / T0 + model->BSIM4wku0 / T1
                         + model->BSIM4pku0 / (T0 * T1);
                  pParam->BSIM4ku0 = 1.0 + tmp1;

                  T0 = pow(Lnew, model->BSIM4llodvth);
                  T1 = pow(W_tmp, model->BSIM4wlodvth);
                  tmp1 = model->BSIM4lkvth0 / T0 + model->BSIM4wkvth0 / T1
                       + model->BSIM4pkvth0 / (T0 * T1);
                  pParam->BSIM4kvth0 = 1.0 + tmp1;
                  pParam->BSIM4kvth0 = sqrt(pParam->BSIM4kvth0*pParam->BSIM4kvth0 + DELTA);

                  T0 = (TRatio - 1.0);
                  pParam->BSIM4ku0temp = pParam->BSIM4ku0 * (1.0 + model->BSIM4tku0 *T0) + DELTA;

                  Inv_saref = 1.0/(model->BSIM4saref + 0.5*Ldrn);
                  Inv_sbref = 1.0/(model->BSIM4sbref + 0.5*Ldrn);
                  pParam->BSIM4inv_od_ref = Inv_saref + Inv_sbref;
                  pParam->BSIM4rho_ref = model->BSIM4ku0 / pParam->BSIM4ku0temp * pParam->BSIM4inv_od_ref;

                  /*high k*/
                  /*Calculate VgsteffVth for mobMod=3*/
                  if(model->BSIM4mobMod==3)
                  {	/*Calculate n @ Vbs=Vds=0*/
                      lt1 = model->BSIM4factor1* pParam->BSIM4sqrtXdep0;
                      T0 = pParam->BSIM4dvt1 * pParam->BSIM4leff / lt1;
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

                      tmp1 = epssub / pParam->BSIM4Xdep0;
                      tmp2 = pParam->BSIM4nfactor * tmp1;
                      tmp3 = (tmp2 + pParam->BSIM4cdsc * Theta0 + pParam->BSIM4cit) / model->BSIM4coxe;
                      if (tmp3 >= -0.5)
                          n0 = 1.0 + tmp3;
                      else
                      {
                          T0 = 1.0 / (3.0 + 8.0 * tmp3);
                          n0 = (1.0 + 3.0 * tmp3) * T0;
                      }

                      T0 = n0 * model->BSIM4vtm;
                      T1 = pParam->BSIM4voffcbn;
                      T2 = T1/T0;
                      if (T2 < -EXP_THRESHOLD)
                      {   T3 = model->BSIM4coxe * MIN_EXP / pParam->BSIM4cdep0;
                          T4 = pParam->BSIM4mstar + T3 * n0;
                      }
                      else if (T2 > EXP_THRESHOLD)
                      {   T3 = model->BSIM4coxe * MAX_EXP / pParam->BSIM4cdep0;
                          T4 = pParam->BSIM4mstar + T3 * n0;
                      }
                      else
                      {  T3 = exp(T2)* model->BSIM4coxe / pParam->BSIM4cdep0;
                         T4 = pParam->BSIM4mstar + T3 * n0;
                      }
                      pParam->BSIM4VgsteffVth = T0 * log(2.0)/T4;
                  }

                  /* New DITS term added in 4.7 */
                  T0 = -pParam->BSIM4dvtp3 * log(pParam->BSIM4leff);
                  DEXP(T0, T1);
                  pParam->BSIM4dvtp2factor = pParam->BSIM4dvtp5 + pParam->BSIM4dvtp2 * T1;

              } /* End of SizeNotFound */

              /*  stress effect */
              if( (here->BSIM4sa > 0.0) && (here->BSIM4sb > 0.0) &&
                        ((here->BSIM4nf == 1.0) || ((here->BSIM4nf > 1.0) && (here->BSIM4sd > 0.0))) )
              {          Inv_sa = 0;
                        Inv_sb = 0;

                         kvsat = model->BSIM4kvsat;
                  if (model->BSIM4kvsat < -1.0 )
                  {   fprintf(stderr, "Warning: KVSAT = %g is too small; -1.0 is used.\n",model->BSIM4kvsat);
                             kvsat = -1.0;
                      }
                      if (model->BSIM4kvsat > 1.0)
                      {   fprintf(stderr, "Warning: KVSAT = %g is too big; 1.0 is used.\n",model->BSIM4kvsat);
                         kvsat = 1.0;
                      }

                        for(i = 0; i < here->BSIM4nf; i++){
                           T0 = 1.0 / here->BSIM4nf / (here->BSIM4sa + 0.5*Ldrn + i * (here->BSIM4sd +Ldrn));
                            T1 = 1.0 / here->BSIM4nf / (here->BSIM4sb + 0.5*Ldrn + i * (here->BSIM4sd +Ldrn));
                           Inv_sa += T0;
                            Inv_sb += T1;
                  }
                  Inv_ODeff = Inv_sa + Inv_sb;
                  rho = model->BSIM4ku0 / pParam->BSIM4ku0temp * Inv_ODeff;
                  T0 = (1.0 + rho)/(1.0 + pParam->BSIM4rho_ref);
                  here->BSIM4u0temp = pParam->BSIM4u0temp * T0;

                  T1 = (1.0 + kvsat * rho)/(1.0 + kvsat * pParam->BSIM4rho_ref);
                  here->BSIM4vsattemp = pParam->BSIM4vsattemp * T1;

                  OD_offset = Inv_ODeff - pParam->BSIM4inv_od_ref;
                  dvth0_lod = model->BSIM4kvth0 / pParam->BSIM4kvth0 * OD_offset;
                  dk2_lod = model->BSIM4stk2 / pow(pParam->BSIM4kvth0, model->BSIM4lodk2) *
                                   OD_offset;
                  deta0_lod = model->BSIM4steta0 / pow(pParam->BSIM4kvth0, model->BSIM4lodeta0) *
                                     OD_offset;
                  here->BSIM4vth0 = pParam->BSIM4vth0 + dvth0_lod;

                  here->BSIM4eta0 = pParam->BSIM4eta0 + deta0_lod;
                  here->BSIM4k2 = pParam->BSIM4k2 + dk2_lod;
               } else {
                      here->BSIM4u0temp = pParam->BSIM4u0temp;
                      here->BSIM4vth0 = pParam->BSIM4vth0;
                      here->BSIM4vsattemp = pParam->BSIM4vsattemp;
                      here->BSIM4eta0 = pParam->BSIM4eta0;
                      here->BSIM4k2 = pParam->BSIM4k2;
              }

              /*  Well Proximity Effect  */
              if (model->BSIM4wpemod)
              { if( (!here->BSIM4scaGiven) && (!here->BSIM4scbGiven) && (!here->BSIM4sccGiven) )
                {   if((here->BSIM4scGiven) && (here->BSIM4sc > 0.0) )
                          {   T1 = here->BSIM4sc + Wdrn;
                        T2 = 1.0 / model->BSIM4scref;
                        here->BSIM4sca = model->BSIM4scref * model->BSIM4scref
                                        / (here->BSIM4sc * T1);
                        here->BSIM4scb = ( (0.1 * here->BSIM4sc + 0.01 * model->BSIM4scref)
                                        * exp(-10.0 * here->BSIM4sc * T2)
                                        - (0.1 * T1 + 0.01 * model->BSIM4scref)
                                        * exp(-10.0 * T1 * T2) ) / Wdrn;
                        here->BSIM4scc = ( (0.05 * here->BSIM4sc + 0.0025 * model->BSIM4scref)
                                        * exp(-20.0 * here->BSIM4sc * T2)
                                        - (0.05 * T1 + 0.0025 * model->BSIM4scref)
                                        * exp(-20.0 * T1 * T2) ) / Wdrn;
                    } else {
                        fprintf(stderr, "Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.\n");
                    }
                }

                       if (here->BSIM4sca < 0.0)
                {
                    printf("Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4sca);
                    here->BSIM4sca = 0.0;
                }
                if (here->BSIM4scb < 0.0)
                {
                    printf("Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4scb);
                    here->BSIM4scb = 0.0;
                }
                if (here->BSIM4scc < 0.0)
                {
                    printf("Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4scc);
                    here->BSIM4scc = 0.0;
                }
                if (here->BSIM4sc < 0.0)
                {
                    printf("Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4sc);
                    here->BSIM4sc = 0.0;
                }
                                /*4.6.2*/
                sceff = here->BSIM4sca + model->BSIM4web * here->BSIM4scb
                      + model->BSIM4wec * here->BSIM4scc;
                here->BSIM4vth0 += pParam->BSIM4kvth0we * sceff;
                here->BSIM4k2 +=  pParam->BSIM4k2we * sceff;
                  T3 =  1.0 + pParam->BSIM4ku0we * sceff;
                if (T3 <= 0.0)
                {         T3 = 0.0;
                        fprintf(stderr, "Warning: ku0we = %g is negatively too high. Negative mobility! \n", pParam->BSIM4ku0we);
                }
                here->BSIM4u0temp *= T3;
              }

            /* adding delvto  */
            here->BSIM4vth0 += here->BSIM4delvto;
            here->BSIM4vfb = pParam->BSIM4vfb + model->BSIM4type * here->BSIM4delvto;

            /* Instance variables calculation  */
            T3 = model->BSIM4type * here->BSIM4vth0
               - here->BSIM4vfb - pParam->BSIM4phi;
            T4 = T3 + T3;
            T5 = 2.5 * T3;
            here->BSIM4vtfbphi1 = (model->BSIM4type == NMOS) ? T4 : T5;
            if (here->BSIM4vtfbphi1 < 0.0)
                here->BSIM4vtfbphi1 = 0.0;

            here->BSIM4vtfbphi2 = 4.0 * T3;
            if (here->BSIM4vtfbphi2 < 0.0)
                here->BSIM4vtfbphi2 = 0.0;

            if (here->BSIM4k2 < 0.0)
            {   T0 = 0.5 * pParam->BSIM4k1 / here->BSIM4k2;
                here->BSIM4vbsc = 0.9 * (pParam->BSIM4phi - T0 * T0);
                if (here->BSIM4vbsc > -3.0)
                    here->BSIM4vbsc = -3.0;
                else if (here->BSIM4vbsc < -30.0)
                    here->BSIM4vbsc = -30.0;
            }
            else
                here->BSIM4vbsc = -30.0;
            if (here->BSIM4vbsc > pParam->BSIM4vbm)
                here->BSIM4vbsc = pParam->BSIM4vbm;
            here->BSIM4k2ox = here->BSIM4k2 * toxe
                              / model->BSIM4toxm;

            here->BSIM4vfbzb = pParam->BSIM4vfbzbfactor
                                +  model->BSIM4type * here->BSIM4vth0 ;

              here->BSIM4cgso = pParam->BSIM4cgso;
              here->BSIM4cgdo = pParam->BSIM4cgdo;

              lnl = log(pParam->BSIM4leff * 1.0e6);
              lnw = log(pParam->BSIM4weff * 1.0e6);
              lnnf = log(here->BSIM4nf);

              bodymode = 5;
              if( ( !model->BSIM4rbps0Given) ||
                  ( !model->BSIM4rbpd0Given) )
                bodymode = 1;
              else
                if( (!model->BSIM4rbsbx0Given && !model->BSIM4rbsby0Given) ||
                      (!model->BSIM4rbdbx0Given && !model->BSIM4rbdby0Given) )
                  bodymode = 3;

              if(here->BSIM4rbodyMod == 2)
                {
                    if (bodymode == 5)
                      {
                        /*rbsbx =  exp( log(model->BSIM4rbsbx0) + model->BSIM4rbsdbxl * lnl +
                              model->BSIM4rbsdbxw * lnw + model->BSIM4rbsdbxnf * lnnf );
                        rbsby =  exp( log(model->BSIM4rbsby0) + model->BSIM4rbsdbyl * lnl +
                              model->BSIM4rbsdbyw * lnw + model->BSIM4rbsdbynf * lnnf );
                               */
                        rbsbx =  model->BSIM4rbsbx0 * exp( model->BSIM4rbsdbxl * lnl +
                              model->BSIM4rbsdbxw * lnw + model->BSIM4rbsdbxnf * lnnf );
                        rbsby =  model->BSIM4rbsby0 * exp( model->BSIM4rbsdbyl * lnl +
                              model->BSIM4rbsdbyw * lnw + model->BSIM4rbsdbynf * lnnf );
                        here->BSIM4rbsb = rbsbx * rbsby / (rbsbx + rbsby);


                        /*rbdbx =  exp( log(model->BSIM4rbdbx0) + model->BSIM4rbsdbxl * lnl +
                              model->BSIM4rbsdbxw * lnw + model->BSIM4rbsdbxnf * lnnf );
                        rbdby =  exp( log(model->BSIM4rbdby0) + model->BSIM4rbsdbyl * lnl +
                              model->BSIM4rbsdbyw * lnw + model->BSIM4rbsdbynf * lnnf );
                                */

                        rbdbx =  model->BSIM4rbdbx0 * exp( model->BSIM4rbsdbxl * lnl +
                              model->BSIM4rbsdbxw * lnw + model->BSIM4rbsdbxnf * lnnf );
                        rbdby =  model->BSIM4rbdby0 * exp(  model->BSIM4rbsdbyl * lnl +
                              model->BSIM4rbsdbyw * lnw + model->BSIM4rbsdbynf * lnnf );

                        here->BSIM4rbdb = rbdbx * rbdby / (rbdbx + rbdby);
                      }

                    if ((bodymode == 3)|| (bodymode == 5))
                      {
                        /*here->BSIM4rbps = exp( log(model->BSIM4rbps0) + model->BSIM4rbpsl * lnl +
                                   model->BSIM4rbpsw * lnw + model->BSIM4rbpsnf * lnnf );
                        here->BSIM4rbpd = exp( log(model->BSIM4rbpd0) + model->BSIM4rbpdl * lnl +
                                                 model->BSIM4rbpdw * lnw + model->BSIM4rbpdnf * lnnf );
                               */
                        here->BSIM4rbps = model->BSIM4rbps0 * exp( model->BSIM4rbpsl * lnl +
                                   model->BSIM4rbpsw * lnw + model->BSIM4rbpsnf * lnnf );
                        here->BSIM4rbpd = model->BSIM4rbpd0 * exp( model->BSIM4rbpdl * lnl +
                                                 model->BSIM4rbpdw * lnw + model->BSIM4rbpdnf * lnnf );

                      }

                    /*rbpbx =  exp( log(model->BSIM4rbpbx0) + model->BSIM4rbpbxl * lnl +
                      	model->BSIM4rbpbxw * lnw + model->BSIM4rbpbxnf * lnnf );
                    rbpby =  exp( log(model->BSIM4rbpby0) + model->BSIM4rbpbyl * lnl +
                      	model->BSIM4rbpbyw * lnw + model->BSIM4rbpbynf * lnnf );
                            */
                            rbpbx =  model->BSIM4rbpbx0 * exp(  model->BSIM4rbpbxl * lnl +
                      	model->BSIM4rbpbxw * lnw + model->BSIM4rbpbxnf * lnnf );
                    rbpby =  model->BSIM4rbpby0 * exp(  model->BSIM4rbpbyl * lnl +
                      	model->BSIM4rbpbyw * lnw + model->BSIM4rbpbynf * lnnf );

                    here->BSIM4rbpb = rbpbx*rbpby/(rbpbx + rbpby);
                  }

              if ((here->BSIM4rbodyMod == 1 ) || ((here->BSIM4rbodyMod == 2 ) && (bodymode == 5)) )
              {   if (here->BSIM4rbdb < 1.0e-3)
                      here->BSIM4grbdb = 1.0e3; /* in mho */
                  else
                      here->BSIM4grbdb = model->BSIM4gbmin + 1.0 / here->BSIM4rbdb;
                  if (here->BSIM4rbpb < 1.0e-3)
                      here->BSIM4grbpb = 1.0e3;
                  else
                      here->BSIM4grbpb = model->BSIM4gbmin + 1.0 / here->BSIM4rbpb;
                  if (here->BSIM4rbps < 1.0e-3)
                      here->BSIM4grbps = 1.0e3;
                  else
                      here->BSIM4grbps = model->BSIM4gbmin + 1.0 / here->BSIM4rbps;
                  if (here->BSIM4rbsb < 1.0e-3)
                      here->BSIM4grbsb = 1.0e3;
                  else
                      here->BSIM4grbsb = model->BSIM4gbmin + 1.0 / here->BSIM4rbsb;
                  if (here->BSIM4rbpd < 1.0e-3)
                      here->BSIM4grbpd = 1.0e3;
                  else
                      here->BSIM4grbpd = model->BSIM4gbmin + 1.0 / here->BSIM4rbpd;

              }

              if((here->BSIM4rbodyMod == 2) && (bodymode == 3))
              {
                      here->BSIM4grbdb = here->BSIM4grbsb = model->BSIM4gbmin;
                  if (here->BSIM4rbpb < 1.0e-3)
                      here->BSIM4grbpb = 1.0e3;
                  else
                      here->BSIM4grbpb = model->BSIM4gbmin + 1.0 / here->BSIM4rbpb;
                  if (here->BSIM4rbps < 1.0e-3)
                      here->BSIM4grbps = 1.0e3;
                  else
                      here->BSIM4grbps = model->BSIM4gbmin + 1.0 / here->BSIM4rbps;
                  if (here->BSIM4rbpd < 1.0e-3)
                      here->BSIM4grbpd = 1.0e3;
                  else
                      here->BSIM4grbpd = model->BSIM4gbmin + 1.0 / here->BSIM4rbpd;
              }

              if((here->BSIM4rbodyMod == 2) && (bodymode == 1))
              {
                      here->BSIM4grbdb = here->BSIM4grbsb = model->BSIM4gbmin;
                      here->BSIM4grbps = here->BSIM4grbpd = 1.0e3;
                  if (here->BSIM4rbpb < 1.0e-3)
                      here->BSIM4grbpb = 1.0e3;
                  else
                      here->BSIM4grbpb = model->BSIM4gbmin + 1.0 / here->BSIM4rbpb;
              }


              /*
               * Process geomertry dependent parasitics
               */

              here->BSIM4grgeltd = model->BSIM4rshg * (here->BSIM4xgw
                      + pParam->BSIM4weffCJ / 3.0 / here->BSIM4ngcon) /
                      (here->BSIM4ngcon * here->BSIM4nf *
                      (Lnew - model->BSIM4xgl));
              if (here->BSIM4grgeltd > 0.0)
                  here->BSIM4grgeltd = 1.0 / here->BSIM4grgeltd;
              else
              {   here->BSIM4grgeltd = 1.0e3; /* mho */
                  if (here->BSIM4rgateMod != 0)
                  printf("Warning: The gate conductance reset to 1.0e3 mho.\n");
              }

              DMCGeff = model->BSIM4dmcg - model->BSIM4dmcgt;
              DMCIeff = model->BSIM4dmci;
              DMDGeff = model->BSIM4dmdg - model->BSIM4dmcgt;

/*              if (here->BSIM4sourcePerimeterGiven)
              {   if (model->BSIM4perMod == 0)
                      here->BSIM4Pseff = here->BSIM4sourcePerimeter;
                  else
                      here->BSIM4Pseff = here->BSIM4sourcePerimeter
                                       - pParam->BSIM4weffCJ * here->BSIM4nf;
              }
              else
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &(here->BSIM4Pseff), &dumPd, &dumAs, &dumAd);
              if (here->BSIM4Pseff < 0.0) /4.6.2/
                      here->BSIM4Pseff = 0.0; */

        /* New Diode Model v4.7*/
              if (here->BSIM4sourcePerimeterGiven)
              {   /* given */
                  if (here->BSIM4sourcePerimeter == 0.0)
                          here->BSIM4Pseff = 0.0;
                  else if (here->BSIM4sourcePerimeter < 0.0)
                  {
                          printf("Warning: Source Perimeter is specified as negative, it is set to zero.\n");
                          here->BSIM4Pseff = 0.0;
                  } else
                  {
                          if (model->BSIM4perMod == 0)
                                  here->BSIM4Pseff = here->BSIM4sourcePerimeter;
                          else
                                  here->BSIM4Pseff = here->BSIM4sourcePerimeter
                                          - pParam->BSIM4weffCJ * here->BSIM4nf;
                  }
              } else /* not given */
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &(here->BSIM4Pseff), &dumPd, &dumAs, &dumAd);

              if (here->BSIM4Pseff < 0.0){ /* v4.7 final check */
                              here->BSIM4Pseff = 0.0;
                              printf("Warning: Pseff is negative, it is set to zero.\n");
              }
              /*  if (here->BSIM4drainPerimeterGiven)
              {   if (model->BSIM4perMod == 0)
                      here->BSIM4Pdeff = here->BSIM4drainPerimeter;
                  else
                      here->BSIM4Pdeff = here->BSIM4drainPerimeter
                                       - pParam->BSIM4weffCJ * here->BSIM4nf;
              }
              else
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &(here->BSIM4Pdeff), &dumAs, &dumAd);
               if (here->BSIM4Pdeff < 0.0) /4.6.2/
                      here->BSIM4Pdeff = 0.0; */

              if (here->BSIM4drainPerimeterGiven)
              {   /* given */
                  if (here->BSIM4drainPerimeter == 0.0)
                                here->BSIM4Pdeff = 0.0;
                  else if (here->BSIM4drainPerimeter < 0.0)
                  {
                                printf("Warning: Drain Perimeter is specified as negative, it is set to zero.\n");
                                here->BSIM4Pdeff = 0.0;
                  } else
                  {
                                if (model->BSIM4perMod == 0)
                                        here->BSIM4Pdeff = here->BSIM4drainPerimeter;
                                else
                                        here->BSIM4Pdeff = here->BSIM4drainPerimeter
                                                - pParam->BSIM4weffCJ * here->BSIM4nf;
                  }
              } else /* not given */
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &(here->BSIM4Pdeff), &dumAs, &dumAd);

              if (here->BSIM4Pdeff < 0.0){
                      here->BSIM4Pdeff = 0.0; /*New Diode v4.7*/
                      printf("Warning: Pdeff is negative, it is set to zero.\n");
              }
              if (here->BSIM4sourceAreaGiven)
                  here->BSIM4Aseff = here->BSIM4sourceArea;
              else
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &dumPd, &(here->BSIM4Aseff), &dumAd);
              if (here->BSIM4Aseff < 0.0){
                  here->BSIM4Aseff = 0.0; /* v4.7 */
                  printf("Warning: Aseff is negative, it is set to zero.\n");
              }
              if (here->BSIM4drainAreaGiven)
                  here->BSIM4Adeff = here->BSIM4drainArea;
              else
                  BSIM4PAeffGeo(here->BSIM4nf, here->BSIM4geoMod, here->BSIM4min,
                                    pParam->BSIM4weffCJ, DMCGeff, DMCIeff, DMDGeff,
                                    &dumPs, &dumPd, &dumAs, &(here->BSIM4Adeff));
              if (here->BSIM4Adeff < 0.0){
                  here->BSIM4Adeff = 0.0; /* v4.7 */
                  printf("Warning: Adeff is negative, it is set to zero.\n");
              }
              /* Processing S/D resistance and conductance below */
              if(here->BSIM4sNodePrime != here->BSIM4sNode)
              {
                 here->BSIM4sourceConductance = 0.0;
                 if(here->BSIM4sourceSquaresGiven)
                 {
                    here->BSIM4sourceConductance = model->BSIM4sheetResistance
                                               * here->BSIM4sourceSquares;
                 } else if (here->BSIM4rgeoMod > 0)
                 {
                    BSIM4RdseffGeo(here->BSIM4nf, here->BSIM4geoMod,
                      here->BSIM4rgeoMod, here->BSIM4min,
                      pParam->BSIM4weffCJ, model->BSIM4sheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 1, &(here->BSIM4sourceConductance));
                 } else
                 {
                    here->BSIM4sourceConductance = 0.0;
                 }

                 if (here->BSIM4sourceConductance > 0.0)
                     here->BSIM4sourceConductance = 1.0
                                            / here->BSIM4sourceConductance;
                 else
                 {
                     here->BSIM4sourceConductance = 1.0e3; /* mho */
                     printf ("Warning: Source conductance reset to 1.0e3 mho.\n");
                 }
              } else
              {
                  here->BSIM4sourceConductance = 0.0;
              }

              if(here->BSIM4dNodePrime != here->BSIM4dNode)
              {
                 here->BSIM4drainConductance = 0.0;
                 if(here->BSIM4drainSquaresGiven)
                 {
                    here->BSIM4drainConductance = model->BSIM4sheetResistance
                                              * here->BSIM4drainSquares;
                 } else if (here->BSIM4rgeoMod > 0)
                 {
                    BSIM4RdseffGeo(here->BSIM4nf, here->BSIM4geoMod,
                      here->BSIM4rgeoMod, here->BSIM4min,
                      pParam->BSIM4weffCJ, model->BSIM4sheetResistance,
                  DMCGeff, DMCIeff, DMDGeff, 0, &(here->BSIM4drainConductance));
                 } else
                 {
                    here->BSIM4drainConductance = 0.0;
                 }

                 if (here->BSIM4drainConductance > 0.0)
                     here->BSIM4drainConductance = 1.0
                                           / here->BSIM4drainConductance;
                 else
                 {
                     here->BSIM4drainConductance = 1.0e3; /* mho */
                     printf ("Warning: Drain conductance reset to 1.0e3 mho.\n");
                  }
              } else
              {
                  here->BSIM4drainConductance = 0.0;
              }

               /* End of Rsd processing */


              Nvtms = model->BSIM4vtm * model->BSIM4SjctEmissionCoeff;
              if ((here->BSIM4Aseff <= 0.0) && (here->BSIM4Pseff <= 0.0))
              {   SourceSatCurrent = 0.0; /* v4.7 */
                  /* SourceSatCurrent = 1.0e-14; */
              }
              else
              {   SourceSatCurrent = here->BSIM4Aseff * model->BSIM4SjctTempSatCurDensity
                                   + here->BSIM4Pseff * model->BSIM4SjctSidewallTempSatCurDensity
                                   + pParam->BSIM4weffCJ * here->BSIM4nf
                                   * model->BSIM4SjctGateSidewallTempSatCurDensity;
              }
              if (SourceSatCurrent > 0.0)
              {   switch(model->BSIM4dioMod)
                  {   case 0:
                          if ((model->BSIM4bvs / Nvtms) > EXP_THRESHOLD)
                              here->BSIM4XExpBVS = model->BSIM4xjbvs * MIN_EXP;
                          else
                              here->BSIM4XExpBVS = model->BSIM4xjbvs * exp(-model->BSIM4bvs / Nvtms);
                          break;
                      case 1:
                          BSIM4DioIjthVjmEval(Nvtms, model->BSIM4ijthsfwd, SourceSatCurrent,
                                              0.0, &(here->BSIM4vjsmFwd));
                          here->BSIM4IVjsmFwd = SourceSatCurrent * exp(here->BSIM4vjsmFwd / Nvtms);
                          break;
                      case 2:
                          if ((model->BSIM4bvs / Nvtms) > EXP_THRESHOLD)
                          {   here->BSIM4XExpBVS = model->BSIM4xjbvs * MIN_EXP;
                              tmp = MIN_EXP;
                          }
                          else
                          {   here->BSIM4XExpBVS = exp(-model->BSIM4bvs / Nvtms);
                              tmp = here->BSIM4XExpBVS;
                              here->BSIM4XExpBVS *= model->BSIM4xjbvs;
                          }

                          BSIM4DioIjthVjmEval(Nvtms, model->BSIM4ijthsfwd, SourceSatCurrent,
                                                     here->BSIM4XExpBVS, &(here->BSIM4vjsmFwd));
                          T0 = exp(here->BSIM4vjsmFwd / Nvtms);
                          here->BSIM4IVjsmFwd = SourceSatCurrent * (T0 - here->BSIM4XExpBVS / T0
                                                + here->BSIM4XExpBVS - 1.0);
                          here->BSIM4SslpFwd = SourceSatCurrent
                                               * (T0 + here->BSIM4XExpBVS / T0) / Nvtms;

                          T2 = model->BSIM4ijthsrev / SourceSatCurrent;
                          if (T2 < 1.0)
                          {   T2 = 10.0;
                              fprintf(stderr, "Warning: ijthsrev too small and set to 10 times IsbSat.\n");
                          }
                          here->BSIM4vjsmRev = -model->BSIM4bvs
                                             - Nvtms * log((T2 - 1.0) / model->BSIM4xjbvs);
                          T1 = model->BSIM4xjbvs * exp(-(model->BSIM4bvs
                             + here->BSIM4vjsmRev) / Nvtms);
                          here->BSIM4IVjsmRev = SourceSatCurrent * (1.0 + T1);
                          here->BSIM4SslpRev = -SourceSatCurrent * T1 / Nvtms;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model->BSIM4dioMod);
                  }
              }

              Nvtmd = model->BSIM4vtm * model->BSIM4DjctEmissionCoeff;
              if ((here->BSIM4Adeff <= 0.0) && (here->BSIM4Pdeff <= 0.0))
              {  /* DrainSatCurrent = 1.0e-14;         v4.7 */
                   DrainSatCurrent = 0.0;
              }
              else
              {   DrainSatCurrent = here->BSIM4Adeff * model->BSIM4DjctTempSatCurDensity
                                  + here->BSIM4Pdeff * model->BSIM4DjctSidewallTempSatCurDensity
                                  + pParam->BSIM4weffCJ * here->BSIM4nf
                                  * model->BSIM4DjctGateSidewallTempSatCurDensity;
              }
              if (DrainSatCurrent > 0.0)
              {   switch(model->BSIM4dioMod)
                  {   case 0:
                          if ((model->BSIM4bvd / Nvtmd) > EXP_THRESHOLD)
                          here->BSIM4XExpBVD = model->BSIM4xjbvd * MIN_EXP;
                          else
                          here->BSIM4XExpBVD = model->BSIM4xjbvd * exp(-model->BSIM4bvd / Nvtmd);
                          break;
                      case 1:
                          BSIM4DioIjthVjmEval(Nvtmd, model->BSIM4ijthdfwd, DrainSatCurrent,
                                              0.0, &(here->BSIM4vjdmFwd));
                          here->BSIM4IVjdmFwd = DrainSatCurrent * exp(here->BSIM4vjdmFwd / Nvtmd);
                          break;
                      case 2:
                          if ((model->BSIM4bvd / Nvtmd) > EXP_THRESHOLD)
                          {   here->BSIM4XExpBVD = model->BSIM4xjbvd * MIN_EXP;
                              tmp = MIN_EXP;
                          }
                          else
                          {   here->BSIM4XExpBVD = exp(-model->BSIM4bvd / Nvtmd);
                              tmp = here->BSIM4XExpBVD;
                              here->BSIM4XExpBVD *= model->BSIM4xjbvd;
                          }

                          BSIM4DioIjthVjmEval(Nvtmd, model->BSIM4ijthdfwd, DrainSatCurrent,
                                              here->BSIM4XExpBVD, &(here->BSIM4vjdmFwd));
                          T0 = exp(here->BSIM4vjdmFwd / Nvtmd);
                          here->BSIM4IVjdmFwd = DrainSatCurrent * (T0 - here->BSIM4XExpBVD / T0
                                              + here->BSIM4XExpBVD - 1.0);
                          here->BSIM4DslpFwd = DrainSatCurrent
                                               * (T0 + here->BSIM4XExpBVD / T0) / Nvtmd;

                          T2 = model->BSIM4ijthdrev / DrainSatCurrent;
                          if (T2 < 1.0)
                          {   T2 = 10.0;
                              fprintf(stderr, "Warning: ijthdrev too small and set to 10 times IdbSat.\n");
                          }
                          here->BSIM4vjdmRev = -model->BSIM4bvd
                                             - Nvtmd * log((T2 - 1.0) / model->BSIM4xjbvd); /* bugfix */
                          T1 = model->BSIM4xjbvd * exp(-(model->BSIM4bvd
                             + here->BSIM4vjdmRev) / Nvtmd);
                          here->BSIM4IVjdmRev = DrainSatCurrent * (1.0 + T1);
                          here->BSIM4DslpRev = -DrainSatCurrent * T1 / Nvtmd;
                          break;
                  default:
                          printf("Specified dioMod = %d not matched\n", model->BSIM4dioMod);
                  }
              }

                /* GEDL current reverse bias */
                T0 = (TRatio - 1.0);
                model->BSIM4njtsstemp = model->BSIM4njts * (1.0 + model->BSIM4tnjts * T0);
                model->BSIM4njtsswstemp = model->BSIM4njtssw * (1.0 + model->BSIM4tnjtssw * T0);
                model->BSIM4njtsswgstemp = model->BSIM4njtsswg * (1.0 + model->BSIM4tnjtsswg * T0);
                model->BSIM4njtsdtemp = model->BSIM4njtsd * (1.0 + model->BSIM4tnjtsd * T0);
                model->BSIM4njtsswdtemp = model->BSIM4njtsswd * (1.0 + model->BSIM4tnjtsswd * T0);
                model->BSIM4njtsswgdtemp = model->BSIM4njtsswgd * (1.0 + model->BSIM4tnjtsswgd * T0);
                T7 = Eg0 / model->BSIM4vtm * T0;
                T9 = model->BSIM4xtss * T7;
                DEXP(T9, T1);
                T9 = model->BSIM4xtsd * T7;
                DEXP(T9, T2);
                T9 = model->BSIM4xtssws * T7;
                DEXP(T9, T3);
                T9 = model->BSIM4xtsswd * T7;
                DEXP(T9, T4);
                T9 = model->BSIM4xtsswgs * T7;
                DEXP(T9, T5);
                T9 = model->BSIM4xtsswgd * T7;
                DEXP(T9, T6);
                                /*IBM TAT*/
                                if(model->BSIM4jtweff < 0.0)
                              {   model->BSIM4jtweff = 0.0;
                          fprintf(stderr, "TAT width dependence effect is negative. Jtweff is clamped to zero.\n");
                      }
                                T11 = sqrt(model->BSIM4jtweff / pParam->BSIM4weffCJ) + 1.0;

                T10 = pParam->BSIM4weffCJ * here->BSIM4nf;
                here->BSIM4SjctTempRevSatCur = T1 * here->BSIM4Aseff * model->BSIM4jtss;
                here->BSIM4DjctTempRevSatCur = T2 * here->BSIM4Adeff * model->BSIM4jtsd;
                here->BSIM4SswTempRevSatCur = T3 * here->BSIM4Pseff * model->BSIM4jtssws;
                here->BSIM4DswTempRevSatCur = T4 * here->BSIM4Pdeff * model->BSIM4jtsswd;
                here->BSIM4SswgTempRevSatCur = T5 * T10 * T11 * model->BSIM4jtsswgs;
                here->BSIM4DswgTempRevSatCur = T6 * T10 * T11 * model->BSIM4jtsswgd;

                if(model->BSIM4mtrlMod != 0 && model->BSIM4mtrlCompatMod == 0)
        {
            /* Calculate TOXP from EOT */
                    /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
                    Vtm0eot = KboQ * model->BSIM4tempeot;
                    Vtmeot  = Vtm0eot;
                    vbieot = Vtm0eot * log(pParam->BSIM4nsd
                                   * pParam->BSIM4ndep / (ni * ni));
                    phieot = Vtm0eot * log(pParam->BSIM4ndep / ni)
                                   + pParam->BSIM4phin + 0.4;
                    tmp2 = here->BSIM4vfb + phieot;
                    vddeot = model->BSIM4type * model->BSIM4vddeot;
                    T0 = model->BSIM4epsrgate * EPS0;
                    if ((pParam->BSIM4ngate > 1.0e18) && (pParam->BSIM4ngate < 1.0e25)
                        && (vddeot > tmp2) && (T0!=0))
                      {
                        T1 = 1.0e6 * CHARGE * T0 * pParam->BSIM4ngate /
                          (model->BSIM4coxe * model->BSIM4coxe);
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
                    lt1 = model->BSIM4factor1* pParam->BSIM4sqrtXdep0;
                    ltw = lt1;
                    T0 = pParam->BSIM4dvt1 * model->BSIM4leffeot / lt1;
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
                    Delt_vth = pParam->BSIM4dvt0 * Theta0 * V0;
                    T0 = pParam->BSIM4dvt1w * model->BSIM4weffeot * model->BSIM4leffeot / ltw;
                    if (T0 < EXP_THRESHOLD)
                      {   T1 = exp(T0);
                      T2 = T1 - 1.0;
                      T3 = T2 * T2;
                      T4 = T3 + 2.0 * T1 * MIN_EXP;
                      T5 = T1 / T4;
                      }
                    else
                      T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
                    T2 = pParam->BSIM4dvt0w * T5 * V0;
                    TempRatioeot =  model->BSIM4tempeot / model->BSIM4tnom - 1.0;
                    T0 = sqrt(1.0 + pParam->BSIM4lpe0 / model->BSIM4leffeot);
                    T1 = pParam->BSIM4k1ox * (T0 - 1.0) * sqrt(phieot)
                      + (pParam->BSIM4kt1 + pParam->BSIM4kt1l / model->BSIM4leffeot) * TempRatioeot;
                    Vth_NarrowW = toxe * phieot
                      / (model->BSIM4weffeot + pParam->BSIM4w0);
                    Lpe_Vb = sqrt(1.0 + pParam->BSIM4lpeb / model->BSIM4leffeot);
                    Vth = model->BSIM4type * here->BSIM4vth0 +
                      (pParam->BSIM4k1ox - pParam->BSIM4k1)*sqrt(phieot)*Lpe_Vb
                      - Delt_vth - T2 + pParam->BSIM4k3 * Vth_NarrowW + T1;

                    /* Calculate n */
                    tmp1 = epssub / pParam->BSIM4Xdep0;
                    here->BSIM4nstar = Vtmeot / Charge_q *
                      (model->BSIM4coxe        + tmp1 + pParam->BSIM4cit);
                    tmp2 = pParam->BSIM4nfactor * tmp1;
                    tmp3 = (tmp2 + pParam->BSIM4cdsc * Theta0 + pParam->BSIM4cit) / model->BSIM4coxe;
                    if (tmp3 >= -0.5)
                      n = 1.0 + tmp3;
                    else
                      {
                        T0 = 1.0 / (3.0 + 8.0 * tmp3);
                        n = (1.0 + 3.0 * tmp3) * T0;
                      }


                    /* Vth correction for Pocket implant */
                    if (pParam->BSIM4dvtp0 > 0.0)
                      {
                        T3 = model->BSIM4leffeot + pParam->BSIM4dvtp0 * 2.0;
                        if (model->BSIM4tempMod < 2)
                          T4 = Vtmeot * log(model->BSIM4leffeot / T3);
                        else
                          T4 = Vtm0eot * log(model->BSIM4leffeot / T3);
                        Vth -= n * T4;
                      }
                    Vgsteff = Vgs_eff-Vth;
                    /* calculating Toxp */
                        T3 = model->BSIM4type * here->BSIM4vth0
               - here->BSIM4vfb - phieot;
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
                        T1 = 1.0 + exp(model->BSIM4bdos * 0.7 * log(T0));
                        Tcen = model->BSIM4ados * 1.9e-9 / T1;
                        toxpf = toxe - epsrox/model->BSIM4epsrsub * Tcen;
                        niter++;
                      } while ((niter<=4)&&(ABS(toxpf-toxpi)>1e-12));
                      here->BSIM4toxp = toxpf;
                      here->BSIM4coxp = epsrox * EPS0 / model->BSIM4toxp;
                } else {
                    here->BSIM4toxp = model->BSIM4toxp;
                    here->BSIM4coxp = model->BSIM4coxp;
                }

              if (BSIM4checkModel(model, here, ckt))
              {
                  SPfrontEnd->IFerrorf (ERR_FATAL, "Fatal error(s) detected during BSIM4.6.0 parameter checking for %s in model %s", model->BSIM4modName, here->BSIM4name);
                  return(E_BADPARM);
              }
         } /* End instance */
    }
    return(OK);
}
