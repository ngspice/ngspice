/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4noi.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
  * Authors: 2008- Wenwei Yang,  Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, 03/04/2004.
 * Modified by Xuemei Xi, 07/29/2005.
 * Modified by Mohan Dunga, 12/13/2006
 * Modified by Wenwei Yang, 07/31/2008.
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
 * Modified by Pankaj Kumar Thakur, 07/23/2012
 **********/

#include "ngspice/ngspice.h"
#include "bsim4def.h"
#include "ngspice/cktdefs.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/noisedef.h"
#include "ngspice/suffix.h"
#include "ngspice/const.h"


/*
 * WDL: 1/f noise model has been smoothed out and enhanced with
 * bulk charge effect as well as physical N* equ. and necessary
 * conversion into the SI unit system.
 */

static double
Eval1ovFNoise(
double Vds,
BSIM4model *model,
BSIM4instance *here,
double freq, double temp)
{
struct bsim4SizeDependParam *pParam;
double cd, esat, DelClm, EffFreq, N0, Nl, Leff, Leffsq;
double T0=0.0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Ssi;

    pParam = here->pParam;
    cd = fabs(here->BSIM4cd);
    Leff = pParam->BSIM4leff - 2.0 * model->BSIM4lintnoi;
    Leffsq = Leff * Leff;
    esat = 2.0 * here->BSIM4vsattemp / here->BSIM4ueff;
    if(model->BSIM4em<=0.0) DelClm = 0.0; /* flicker noise modified -JX  */
    else {
            T0 = ((((Vds - here->BSIM4Vdseff) / pParam->BSIM4litl)
                       + model->BSIM4em) / esat);
            DelClm = pParam->BSIM4litl * log (MAX(T0, N_MINLOG));
            if (DelClm < 0.0)        DelClm = 0.0;  /* bugfix */
    }
    EffFreq = pow(freq, model->BSIM4ef);
    T1 = CHARGE * CHARGE * CONSTboltz * cd * temp * here->BSIM4ueff;
    T2 = 1.0e10 * EffFreq * here->BSIM4Abulk * model->BSIM4coxe * Leffsq;
    N0 = model->BSIM4coxe * here->BSIM4Vgsteff / CHARGE;
    Nl = model->BSIM4coxe * here->BSIM4Vgsteff
       * (1.0 - here->BSIM4AbovVgst2Vtm * here->BSIM4Vdseff) / CHARGE;

    T3 = model->BSIM4oxideTrapDensityA
       * log(MAX(((N0 + here->BSIM4nstar) / (Nl + here->BSIM4nstar)), N_MINLOG));
    T4 = model->BSIM4oxideTrapDensityB * (N0 - Nl);
    T5 = model->BSIM4oxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

    T6 = CONSTboltz * temp * cd * cd;
    T7 = 1.0e10 * EffFreq * Leffsq * pParam->BSIM4weff * here->BSIM4nf;
    T8 = model->BSIM4oxideTrapDensityA + model->BSIM4oxideTrapDensityB * Nl
       + model->BSIM4oxideTrapDensityC * Nl * Nl;
    T9 = (Nl + here->BSIM4nstar) * (Nl + here->BSIM4nstar);
    Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
    return Ssi;
}


int
BSIM4noise (
int mode, int operation,
GENmodel *inModel,
CKTcircuit *ckt,
Ndata *data,
double *OnDens)
{
NOISEAN *job = (NOISEAN *) ckt->CKTcurJob;

BSIM4model *model = (BSIM4model *)inModel;
BSIM4instance *here;
struct bsim4SizeDependParam *pParam;
char name[N_MXVLNTH];
double tempOnoise;
double tempInoise;
double noizDens[BSIM4NSRCS];
double lnNdens[BSIM4NSRCS];

double T0, T1, T2, T3, T4, T5, T6, T7, T8, T10, T11;
double Vds, Ssi, Swi;
double tmp=0.0, gdpr, gspr, npart_theta=0.0, npart_beta=0.0, igsquare, bodymode;

/* tnoiMod=2 (v4.7) */
double eta, Leff, Lvsat, gamma, delta, epsilon, GammaGd0=0.0;
double npart_c, sigrat=0.0, C0, omega, ctnoi=0.0;

int i;

double m;

    /* define the names of the noise sources */
    static char *BSIM4nNames[BSIM4NSRCS] =
    {   /* Note that we have to keep the order */
        ".rd",              /* noise due to rd */
        ".rs",              /* noise due to rs */
        ".rg",              /* noise due to rgeltd */
        ".rbps",            /* noise due to rbps */
        ".rbpd",            /* noise due to rbpd */
        ".rbpb",            /* noise due to rbpb */
        ".rbsb",            /* noise due to rbsb */
        ".rbdb",            /* noise due to rbdb */
        ".id",              /* noise due to id (for tnoiMod2: uncorrelated portion only) */
        ".1overf",          /* flicker (1/f) noise */
        ".igs",             /* shot noise due to IGS */
        ".igd",             /* shot noise due to IGD */
        ".igb",             /* shot noise due to IGB */
        ".corl",            /* contribution of correlated drain and induced gate noise */
        ""                  /* total transistor noise */
    };

    for (; model != NULL; model = model->BSIM4nextModel)
    {
         if(model->BSIM4tnoiMod != 2) {
             noizDens[BSIM4CORLNOIZ] = 0.0;
             lnNdens[BSIM4CORLNOIZ] = N_MINLOG;
         }
         for (here = model->BSIM4instances; here != NULL;
              here = here->BSIM4nextInstance)
         {    pParam = here->pParam;
              switch (operation)
              {  case N_OPEN:
                     /* see if we have to to produce a summary report */
                     /* if so, name all the noise generators */

                      if (job->NStpsSm != 0)
                      {   switch (mode)
                          {  case N_DENS:
                                  for (i = 0; i < BSIM4NSRCS; i++)
                                  {    (void) sprintf(name, "onoise.%s%s",
                                                      here->BSIM4name,
                                                      BSIM4nNames[i]);
                                       data->namelist = TREALLOC(IFuid,
                                             data->namelist,
                                             data->numPlots + 1);
                                       if (!data->namelist)
                                           return(E_NOMEM);
                                       SPfrontEnd->IFnewUid (ckt,
                                          &(data->namelist[data->numPlots++]),
                                          NULL, name, UID_OTHER, NULL);
                                       /* we've added one more plot */
                                  }
                                  break;
                             case INT_NOIZ:
                                  for (i = 0; i < BSIM4NSRCS; i++)
                                  {    (void) sprintf(name, "onoise_total.%s%s",
                                                      here->BSIM4name,
                                                      BSIM4nNames[i]);
                                       data->namelist = TREALLOC(IFuid,
                                             data->namelist,
                                             data->numPlots + 1);
                                       if (!data->namelist)
                                           return(E_NOMEM);
                                       SPfrontEnd->IFnewUid (ckt,
                                          &(data->namelist[data->numPlots++]),
                                          NULL, name, UID_OTHER, NULL);
                                       /* we've added one more plot */

                                       (void) sprintf(name, "inoise_total.%s%s",
                                                      here->BSIM4name,
                                                      BSIM4nNames[i]);
                                       data->namelist = TREALLOC(IFuid,
                                             data->namelist,
                                             data->numPlots + 1);
                                       if (!data->namelist)
                                           return(E_NOMEM);
                                       SPfrontEnd->IFnewUid (ckt,
                                          &(data->namelist[data->numPlots++]),
                                          NULL, name, UID_OTHER, NULL);
                                       /* we've added one more plot */
                                  }
                                  break;
                          }
                      }
                      break;
                 case N_CALC:
                      m = here->BSIM4m;
                      switch (mode)
                      {  case N_DENS:
                              if (model->BSIM4tnoiMod == 0)
                              {   if (model->BSIM4rdsMod == 0)
                                  {   gspr = here->BSIM4sourceConductance;
                                      gdpr = here->BSIM4drainConductance;
                                      if (here->BSIM4grdsw > 0.0)
                                          tmp = 1.0 / here->BSIM4grdsw; /* tmp used below */
                                      else
                                          tmp = 0.0;
                                  }
                                  else
                                  {   gspr = here->BSIM4gstot;
                                      gdpr = here->BSIM4gdtot;
                                      tmp = 0.0;
                                  }
                              }
                              else if(model->BSIM4tnoiMod == 1)
                              {   T5 = here->BSIM4Vgsteff / here->BSIM4EsatL;
                                  T5 *= T5;
                                  npart_beta = model->BSIM4rnoia * (1.0 + T5
                                             * model->BSIM4tnoia * pParam->BSIM4leff);
                                  npart_theta = model->BSIM4rnoib * (1.0 + T5
                                              * model->BSIM4tnoib * pParam->BSIM4leff);
                                  if(npart_theta > 0.9)
                                     npart_theta = 0.9;
                                  if(npart_theta > 0.9 * npart_beta)
                                     npart_theta = 0.9 * npart_beta; //4.6.2

                                  if (model->BSIM4rdsMod == 0)
                                  {   gspr = here->BSIM4sourceConductance;
                                      gdpr = here->BSIM4drainConductance;
                                  }
                                  else
                                  {   gspr = here->BSIM4gstot;
                                      gdpr = here->BSIM4gdtot;
                                  }

                                  if ((*(ckt->CKTstates[0] + here->BSIM4vds)) >= 0.0)
                                      gspr = gspr * (1.0 + npart_theta * npart_theta * gspr
                                           / here->BSIM4IdovVds);
                                  else
                                      gdpr = gdpr * (1.0 + npart_theta * npart_theta * gdpr
                                           / here->BSIM4IdovVds);
                              }
                              else
                              {   /* tnoiMod=2 (v4.7) */

                                  if (model->BSIM4rdsMod == 0)
                                  {   gspr = here->BSIM4sourceConductance;
                                      gdpr = here->BSIM4drainConductance;
                                  }
                                  else
                                  {   gspr = here->BSIM4gstot;
                                      gdpr = here->BSIM4gdtot;
                                  }

                              }

                              NevalSrc(&noizDens[BSIM4RDNOIZ],
                                       &lnNdens[BSIM4RDNOIZ], ckt, THERMNOISE,
                                       here->BSIM4dNodePrime, here->BSIM4dNode,
                                       gdpr * m);

                              NevalSrc(&noizDens[BSIM4RSNOIZ],
                                       &lnNdens[BSIM4RSNOIZ], ckt, THERMNOISE,
                                       here->BSIM4sNodePrime, here->BSIM4sNode,
                                       gspr * m);


                              if (here->BSIM4rgateMod == 1)
                              {   NevalSrc(&noizDens[BSIM4RGNOIZ],
                                       &lnNdens[BSIM4RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4gNodePrime, here->BSIM4gNodeExt,
                                       here->BSIM4grgeltd * m);
                              }
                              else if (here->BSIM4rgateMod == 2)
                              {
                                T0 = 1.0 + here->BSIM4grgeltd/here->BSIM4gcrg;
                                T1 = T0 * T0;
                                  NevalSrc(&noizDens[BSIM4RGNOIZ],
                                       &lnNdens[BSIM4RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4gNodePrime, here->BSIM4gNodeExt,
                                       here->BSIM4grgeltd * m / T1);
                              }
                              else if (here->BSIM4rgateMod == 3)
                              {   NevalSrc(&noizDens[BSIM4RGNOIZ],
                                       &lnNdens[BSIM4RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4gNodeMid, here->BSIM4gNodeExt,
                                       here->BSIM4grgeltd * m);
                              }
                              else
                              {    noizDens[BSIM4RGNOIZ] = 0.0;
                                   lnNdens[BSIM4RGNOIZ] =
                                          log(MAX(noizDens[BSIM4RGNOIZ], N_MINLOG));
                              }

                                    bodymode = 5;
                                    if (here->BSIM4rbodyMod == 2)
                                    {        if( ( !model->BSIM4rbps0Given) ||
                                      ( !model->BSIM4rbpd0Given) )
                                             bodymode = 1;
                                           else
                                     if( (!model->BSIM4rbsbx0Given && !model->BSIM4rbsby0Given) ||
                                          (!model->BSIM4rbdbx0Given && !model->BSIM4rbdby0Given) )
                                             bodymode = 3;
                                }

                              if (here->BSIM4rbodyMod)
                              {
                                if(bodymode == 5)
                                  {
                                    NevalSrc(&noizDens[BSIM4RBPSNOIZ],
                                             &lnNdens[BSIM4RBPSNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4sbNode,
                                             here->BSIM4grbps * m);
                                    NevalSrc(&noizDens[BSIM4RBPDNOIZ],
                                             &lnNdens[BSIM4RBPDNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4dbNode,
                                             here->BSIM4grbpd * m);
                                    NevalSrc(&noizDens[BSIM4RBPBNOIZ],
                                             &lnNdens[BSIM4RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4bNode,
                                             here->BSIM4grbpb * m);
                                    NevalSrc(&noizDens[BSIM4RBSBNOIZ],
                                             &lnNdens[BSIM4RBSBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNode, here->BSIM4sbNode,
                                             here->BSIM4grbsb * m);
                                    NevalSrc(&noizDens[BSIM4RBDBNOIZ],
                                             &lnNdens[BSIM4RBDBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNode, here->BSIM4dbNode,
                                             here->BSIM4grbdb * m);
                                  }
                                if(bodymode == 3)
                                  {
                                    NevalSrc(&noizDens[BSIM4RBPSNOIZ],
                                             &lnNdens[BSIM4RBPSNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4sbNode,
                                             here->BSIM4grbps * m);
                                    NevalSrc(&noizDens[BSIM4RBPDNOIZ],
                                             &lnNdens[BSIM4RBPDNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4dbNode,
                                             here->BSIM4grbpd * m);
                                    NevalSrc(&noizDens[BSIM4RBPBNOIZ],
                                             &lnNdens[BSIM4RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4bNode,
                                             here->BSIM4grbpb * m);
                                     noizDens[BSIM4RBSBNOIZ] = noizDens[BSIM4RBDBNOIZ] = 0.0;
                                     lnNdens[BSIM4RBSBNOIZ] =
                                       log(MAX(noizDens[BSIM4RBSBNOIZ], N_MINLOG));
                                     lnNdens[BSIM4RBDBNOIZ] =
                                       log(MAX(noizDens[BSIM4RBDBNOIZ], N_MINLOG));
                                  }
                                if(bodymode == 1)
                                  {
                                    NevalSrc(&noizDens[BSIM4RBPBNOIZ],
                                             &lnNdens[BSIM4RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4bNodePrime, here->BSIM4bNode,
                                             here->BSIM4grbpb * m);
                                    noizDens[BSIM4RBPSNOIZ] = noizDens[BSIM4RBPDNOIZ] = 0.0;
                                    noizDens[BSIM4RBSBNOIZ] = noizDens[BSIM4RBDBNOIZ] = 0.0;
                                    lnNdens[BSIM4RBPSNOIZ] =
                                      log(MAX(noizDens[BSIM4RBPSNOIZ], N_MINLOG));
                                    lnNdens[BSIM4RBPDNOIZ] =
                                      log(MAX(noizDens[BSIM4RBPDNOIZ], N_MINLOG));
                                    lnNdens[BSIM4RBSBNOIZ] =
                                      log(MAX(noizDens[BSIM4RBSBNOIZ], N_MINLOG));
                                    lnNdens[BSIM4RBDBNOIZ] =
                                      log(MAX(noizDens[BSIM4RBDBNOIZ], N_MINLOG));
                                  }
                              }
                              else
                              {   noizDens[BSIM4RBPSNOIZ] = noizDens[BSIM4RBPDNOIZ] = 0.0;
                                  noizDens[BSIM4RBPBNOIZ] = 0.0;
                                  noizDens[BSIM4RBSBNOIZ] = noizDens[BSIM4RBDBNOIZ] = 0.0;
                                  lnNdens[BSIM4RBPSNOIZ] =
                                          log(MAX(noizDens[BSIM4RBPSNOIZ], N_MINLOG));
                                  lnNdens[BSIM4RBPDNOIZ] =
                                          log(MAX(noizDens[BSIM4RBPDNOIZ], N_MINLOG));
                                  lnNdens[BSIM4RBPBNOIZ] =
                                          log(MAX(noizDens[BSIM4RBPBNOIZ], N_MINLOG));
                                  lnNdens[BSIM4RBSBNOIZ] =
                                          log(MAX(noizDens[BSIM4RBSBNOIZ], N_MINLOG));
                                  lnNdens[BSIM4RBDBNOIZ] =
                                          log(MAX(noizDens[BSIM4RBDBNOIZ], N_MINLOG));
                              }

                              if(model->BSIM4tnoiMod == 2)
                              {
                                  eta = 1.0 - here->BSIM4Vdseff * here->BSIM4AbovVgst2Vtm;
                                  T0 = 1.0 - eta;
                                  T1 = 1.0 + eta;
                                  T2 = T1 + 2.0 * here->BSIM4Abulk * model->BSIM4vtm / here->BSIM4Vgsteff;
                                  Leff = pParam->BSIM4leff;
                                  Lvsat = Leff * (1.0 + here->BSIM4Vdseff / here->BSIM4EsatL);
                                  T6 = Leff / Lvsat;
                                  /*Unwanted code for T5 commented*/
                                  /*T5 = here->BSIM4Vgsteff / here->BSIM4EsatL;
                                  T5 = T5 * T5;
                                  */
                                  gamma = T6 * (0.5 * T1 + T0 * T0 / (6.0 * T2));
                                  T3 = T2 * T2;
                                  T4 = T0 * T0;
                                  T5 = T3 * T3;
                                  delta = (T1 / T3 - (5.0 * T1 + T2) * T4 / (15.0 * T5) + T4 * T4 / (9.0 * T5 * T2)) / (6.0 * T6 * T6 * T6);
                                  T7 = T0 / T2;
                                  epsilon = (T7 - T7 * T7 * T7 / 3.0) / (6.0 * T6);

                                  T8 = here->BSIM4Vgsteff / here->BSIM4EsatL;
                                  T8 *= T8;
                                  npart_c = model->BSIM4rnoic * (1.0 + T8
                                          * model->BSIM4tnoic * Leff);
                                  ctnoi = epsilon / sqrt(gamma * delta)
                                      * (2.5316 * npart_c);

                                  npart_beta = model->BSIM4rnoia * (1.0 + T8
                                      * model->BSIM4tnoia * Leff);
                                  npart_theta = model->BSIM4rnoib * (1.0 + T8
                                      * model->BSIM4tnoib * Leff);
                                  gamma = gamma * (3.0 * npart_beta * npart_beta);
                                  delta = delta * (3.75 * npart_theta * npart_theta);

                                  GammaGd0 = gamma * here->BSIM4noiGd0;
                                  C0 = here->BSIM4Coxeff * pParam->BSIM4weffCV * here->BSIM4nf * pParam->BSIM4leffCV;
                                  T0 = C0 / here->BSIM4noiGd0;
                                  sigrat = T0 * sqrt(delta / gamma);
                              }
                              switch(model->BSIM4tnoiMod)
                              {  case 0:
                                      T0 = here->BSIM4ueff * fabs(here->BSIM4qinv);
                                      T1 = T0 * tmp + pParam->BSIM4leff
                                         * pParam->BSIM4leff;
                                      NevalSrc(&noizDens[BSIM4IDNOIZ],
                                               &lnNdens[BSIM4IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4dNodePrime,
                                               here->BSIM4sNodePrime,
                                               (T0 / T1) * model->BSIM4ntnoi * m);
                                      break;
                                 case 1:
                                      T0 = here->BSIM4gm + here->BSIM4gmbs + here->BSIM4gds;
                                      T0 *= T0;
                                      igsquare = npart_theta * npart_theta * T0 / here->BSIM4IdovVds;
                                      T1 = npart_beta * (here->BSIM4gm
                                         + here->BSIM4gmbs) + here->BSIM4gds;
                                      T2 = T1 * T1 / here->BSIM4IdovVds;
                                      NevalSrc(&noizDens[BSIM4IDNOIZ],
                                               &lnNdens[BSIM4IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4dNodePrime,
                                               here->BSIM4sNodePrime, (T2 - igsquare) * m);
                                      break;
                                  case 2:
                                      T2 = GammaGd0;
                                      T3 = ctnoi * ctnoi;
                                      T4 = 1.0 - T3;
                                      NevalSrc(&noizDens[BSIM4IDNOIZ],
                                               &lnNdens[BSIM4IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4dNodePrime,
                                               here->BSIM4sNodePrime, T2 * T4 * m);

                                     /* Evaluate output noise due to two correlated noise sources */
                                     omega = 2.0 * M_PI * data->freq;
                                     T5 = omega * sigrat;
                                     T6 = T5 * T5;
                                     T7 = T6 / (1.0 + T6);

                                     if (here->BSIM4mode >= 0)  {
                                         NevalSrc2(&noizDens[BSIM4CORLNOIZ],
                                               &lnNdens[BSIM4CORLNOIZ], ckt,
                                               THERMNOISE, here->BSIM4dNodePrime,
                                               here->BSIM4sNodePrime, T2 * T3 * m,
                                               here->BSIM4gNodePrime,
                                               here->BSIM4sNodePrime,
                                               T2 * T7 * m, 0.5 * M_PI);
                                     }
                                     else
                                     {
                                         NevalSrc2(&noizDens[BSIM4CORLNOIZ],
                                               &lnNdens[BSIM4CORLNOIZ], ckt,
                                               THERMNOISE, here->BSIM4sNodePrime,
                                               here->BSIM4dNodePrime, T2 * T3 * m,
                                               here->BSIM4gNodePrime,
                                               here->BSIM4dNodePrime,
                                               T2 * T7 * m, 0.5 * M_PI);
                                     }
                                     break;
                              }

                              NevalSrc(&noizDens[BSIM4FLNOIZ], (double*) NULL,
                                       ckt, N_GAIN, here->BSIM4dNodePrime,
                                       here->BSIM4sNodePrime, (double) 0.0);

                              switch(model->BSIM4fnoiMod)
                              {  case 0:
                                      noizDens[BSIM4FLNOIZ] *= m * model->BSIM4kf
                                            * exp(model->BSIM4af
                                            * log(MAX(fabs(here->BSIM4cd),
                                            N_MINLOG)))
                                            / (pow(data->freq, model->BSIM4ef)
                                            * pParam->BSIM4leff
                                            * pParam->BSIM4leff
                                            * model->BSIM4coxe);
                                      break;
                                 case 1:
                                      Vds = *(ckt->CKTstates[0] + here->BSIM4vds);
                                      if (Vds < 0.0)
                                          Vds = -Vds;

                                      Ssi = Eval1ovFNoise(Vds, model, here,
                                          data->freq, ckt->CKTtemp);
                                      T10 = model->BSIM4oxideTrapDensityA
                                          * CONSTboltz * ckt->CKTtemp;
                                      T11 = pParam->BSIM4weff * here->BSIM4nf * pParam->BSIM4leff
                                          * pow(data->freq, model->BSIM4ef) * 1.0e10
                                          * here->BSIM4nstar * here->BSIM4nstar;
                                      Swi = T10 / T11 * here->BSIM4cd
                                          * here->BSIM4cd;
                                      T1 = Swi + Ssi;
                                      if (T1 > 0.0)
                                          noizDens[BSIM4FLNOIZ] *= m * (Ssi * Swi) / T1;
                                      else
                                          noizDens[BSIM4FLNOIZ] *= 0.0;
                                      break;
                              }

                              lnNdens[BSIM4FLNOIZ] =
                                     log(MAX(noizDens[BSIM4FLNOIZ], N_MINLOG));


                        if(here->BSIM4mode >= 0) {  /* bugfix  */
                              NevalSrc(&noizDens[BSIM4IGSNOIZ],
                                   &lnNdens[BSIM4IGSNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4gNodePrime, here->BSIM4sNodePrime,
                                   m * (here->BSIM4Igs + here->BSIM4Igcs));
                              NevalSrc(&noizDens[BSIM4IGDNOIZ],
                                   &lnNdens[BSIM4IGDNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4gNodePrime, here->BSIM4dNodePrime,
                                   m * (here->BSIM4Igd + here->BSIM4Igcd));
                        } else {
                              NevalSrc(&noizDens[BSIM4IGSNOIZ],
                                   &lnNdens[BSIM4IGSNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4gNodePrime, here->BSIM4sNodePrime,
                                   m * (here->BSIM4Igs + here->BSIM4Igcd));
                              NevalSrc(&noizDens[BSIM4IGDNOIZ],
                                   &lnNdens[BSIM4IGDNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4gNodePrime, here->BSIM4dNodePrime,
                                   m * (here->BSIM4Igd + here->BSIM4Igcs));
                        }
                              NevalSrc(&noizDens[BSIM4IGBNOIZ],
                                   &lnNdens[BSIM4IGBNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4gNodePrime, here->BSIM4bNodePrime,
                                   m * here->BSIM4Igb);


                              noizDens[BSIM4TOTNOIZ] = noizDens[BSIM4RDNOIZ]
                                     + noizDens[BSIM4RSNOIZ] + noizDens[BSIM4RGNOIZ]
                                     + noizDens[BSIM4RBPSNOIZ] + noizDens[BSIM4RBPDNOIZ]
                                     + noizDens[BSIM4RBPBNOIZ]
                                     + noizDens[BSIM4RBSBNOIZ] + noizDens[BSIM4RBDBNOIZ]
                                     + noizDens[BSIM4IDNOIZ] + noizDens[BSIM4FLNOIZ]
                                     + noizDens[BSIM4IGSNOIZ] + noizDens[BSIM4IGDNOIZ]
                                     + noizDens[BSIM4IGBNOIZ] + noizDens[BSIM4CORLNOIZ];
                              lnNdens[BSIM4TOTNOIZ] =
                                     log(MAX(noizDens[BSIM4TOTNOIZ], N_MINLOG));

                              *OnDens += noizDens[BSIM4TOTNOIZ];

                              if (data->delFreq == 0.0)
                              {   /* if we haven't done any previous
                                     integration, we need to initialize our
                                     "history" variables.
                                    */

                                  for (i = 0; i < BSIM4NSRCS; i++)
                                  {    here->BSIM4nVar[LNLSTDENS][i] =
                                             lnNdens[i];
                                  }

                                  /* clear out our integration variables
                                     if it's the first pass
                                   */
                                  if (data->freq ==
                                      job->NstartFreq)
                                  {   for (i = 0; i < BSIM4NSRCS; i++)
                                      {    here->BSIM4nVar[OUTNOIZ][i] = 0.0;
                                           here->BSIM4nVar[INNOIZ][i] = 0.0;
                                      }
                                  }
                              }
                              else
                              {   /* data->delFreq != 0.0,
                                     we have to integrate.
                                   */
                                  for (i = 0; i < BSIM4NSRCS; i++)
                                  {    if (i != BSIM4TOTNOIZ)
                                       {   tempOnoise = Nintegrate(noizDens[i],
                                                lnNdens[i],
                                                here->BSIM4nVar[LNLSTDENS][i],
                                                data);
                                           tempInoise = Nintegrate(noizDens[i]
                                                * data->GainSqInv, lnNdens[i]
                                                + data->lnGainInv,
                                                here->BSIM4nVar[LNLSTDENS][i]
                                                + data->lnGainInv, data);
                                           here->BSIM4nVar[LNLSTDENS][i] =
                                                lnNdens[i];
                                           data->outNoiz += tempOnoise;
                                           data->inNoise += tempInoise;
                                           if (job->NStpsSm != 0)
                                           {   here->BSIM4nVar[OUTNOIZ][i]
                                                     += tempOnoise;
                                               here->BSIM4nVar[OUTNOIZ][BSIM4TOTNOIZ]
                                                     += tempOnoise;
                                               here->BSIM4nVar[INNOIZ][i]
                                                     += tempInoise;
                                               here->BSIM4nVar[INNOIZ][BSIM4TOTNOIZ]
                                                     += tempInoise;
                                           }
                                       }
                                  }
                              }
                              if (data->prtSummary)
                              {   for (i = 0; i < BSIM4NSRCS; i++)
                                  {    /* print a summary report */
                                       data->outpVector[data->outNumber++]
                                             = noizDens[i];
                                  }
                              }
                              break;
                         case INT_NOIZ:
                              /* already calculated, just output */
                              if (job->NStpsSm != 0)
                              {   for (i = 0; i < BSIM4NSRCS; i++)
                                  {    data->outpVector[data->outNumber++]
                                             = here->BSIM4nVar[OUTNOIZ][i];
                                       data->outpVector[data->outNumber++]
                                             = here->BSIM4nVar[INNOIZ][i];
                                  }
                              }
                              break;
                      }
                      break;
                 case N_CLOSE:
                      /* do nothing, the main calling routine will close */
                      return (OK);
                      break;   /* the plots */
              }       /* switch (operation) */
         }    /* for here */
    }    /* for model */

    return(OK);
}
