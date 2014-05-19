/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4noi.c of BSIM4v7.0.
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
 **********/

#include "ngspice/ngspice.h"
#include "bsim4v7def.h."
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
BSIM4v7model *model,
BSIM4v7instance *here,
double freq, double temp)
{
struct bsim4SizeDependParam *pParam;
double cd, esat, DelClm, EffFreq, N0, Nl, Leff, Leffsq;
double T0=0.0, T1, T2, T3, T4, T5, T6, T7, T8, T9, Ssi;

    pParam = here->pParam;
    cd = fabs(here->BSIM4v7cd);
    Leff = pParam->BSIM4v7leff - 2.0 * model->BSIM4v7lintnoi;
    Leffsq = Leff * Leff;
    esat = 2.0 * here->BSIM4v7vsattemp / here->BSIM4v7ueff;
    if(model->BSIM4v7em<=0.0) DelClm = 0.0; /* flicker noise modified -JX  */
    else {
            T0 = ((((Vds - here->BSIM4v7Vdseff) / pParam->BSIM4v7litl)
                       + model->BSIM4v7em) / esat);
            DelClm = pParam->BSIM4v7litl * log (MAX(T0, N_MINLOG));
            if (DelClm < 0.0)        DelClm = 0.0;  /* bugfix */
    }
    EffFreq = pow(freq, model->BSIM4v7ef);
    T1 = CHARGE * CHARGE * CONSTboltz * cd * temp * here->BSIM4v7ueff;
    T2 = 1.0e10 * EffFreq * here->BSIM4v7Abulk * model->BSIM4v7coxe * Leffsq;
    N0 = model->BSIM4v7coxe * here->BSIM4v7Vgsteff / CHARGE;
    Nl = model->BSIM4v7coxe * here->BSIM4v7Vgsteff
       * (1.0 - here->BSIM4v7AbovVgst2Vtm * here->BSIM4v7Vdseff) / CHARGE;

    T3 = model->BSIM4v7oxideTrapDensityA
       * log(MAX(((N0 + here->BSIM4v7nstar) / (Nl + here->BSIM4v7nstar)), N_MINLOG));
    T4 = model->BSIM4v7oxideTrapDensityB * (N0 - Nl);
    T5 = model->BSIM4v7oxideTrapDensityC * 0.5 * (N0 * N0 - Nl * Nl);

    T6 = CONSTboltz * temp * cd * cd;
    T7 = 1.0e10 * EffFreq * Leffsq * pParam->BSIM4v7weff * here->BSIM4v7nf;
    T8 = model->BSIM4v7oxideTrapDensityA + model->BSIM4v7oxideTrapDensityB * Nl
       + model->BSIM4v7oxideTrapDensityC * Nl * Nl;
    T9 = (Nl + here->BSIM4v7nstar) * (Nl + here->BSIM4v7nstar);
    Ssi = T1 / T2 * (T3 + T4 + T5) + T6 / T7 * DelClm * T8 / T9;
    return Ssi;
}


int
BSIM4v7noise (
int mode, int operation,
GENmodel *inModel,
CKTcircuit *ckt,
Ndata *data,
double *OnDens)
{
NOISEAN *job = (NOISEAN *) ckt->CKTcurJob;

BSIM4v7model *model = (BSIM4v7model *)inModel;
BSIM4v7instance *here;
struct bsim4SizeDependParam *pParam;
char name[N_MXVLNTH];
double tempOnoise;
double tempInoise;
double noizDens[BSIM4v7NSRCS];
double lnNdens[BSIM4v7NSRCS];

double T0, T1, T2, T3, T4, T5, T6, T7, T8, T10, T11;
double Vds, Ssi, Swi;
double tmp=0.0, gdpr, gspr, npart_theta=0.0, npart_beta=0.0, igsquare, bodymode;

/* tnoiMod=2 (v4.7) */
double eta, Leff, Lvsat, gamma, delta, epsilon, GammaGd0=0.0;
double npart_c, sigrat=0.0, C0, omega, ctnoi=0.0;

int i;

double m;

    /* define the names of the noise sources */
    static char *BSIM4v7nNames[BSIM4v7NSRCS] =
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

    for (; model != NULL; model = model->BSIM4v7nextModel)
    {
         if(model->BSIM4v7tnoiMod != 2) {
             noizDens[BSIM4v7CORLNOIZ] = 0.0;
             lnNdens[BSIM4v7CORLNOIZ] = N_MINLOG;
         }
         for (here = model->BSIM4v7instances; here != NULL;
              here = here->BSIM4v7nextInstance)
         {    pParam = here->pParam;
              switch (operation)
              {  case N_OPEN:
                     /* see if we have to to produce a summary report */
                     /* if so, name all the noise generators */

                      if (job->NStpsSm != 0)
                      {   switch (mode)
                          {  case N_DENS:
                                  for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    (void) sprintf(name, "onoise.%s%s",
                                                      here->BSIM4v7name,
                                                      BSIM4v7nNames[i]);
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
                                  for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    (void) sprintf(name, "onoise_total.%s%s",
                                                      here->BSIM4v7name,
                                                      BSIM4v7nNames[i]);
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
                                                      here->BSIM4v7name,
                                                      BSIM4v7nNames[i]);
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
                      m = here->BSIM4v7m;
                      switch (mode)
                      {  case N_DENS:
                              if (model->BSIM4v7tnoiMod == 0)
                              {   if (model->BSIM4v7rdsMod == 0)
                                  {   gspr = here->BSIM4v7sourceConductance;
                                      gdpr = here->BSIM4v7drainConductance;
                                      if (here->BSIM4v7grdsw > 0.0)
                                          tmp = 1.0 / here->BSIM4v7grdsw; /* tmp used below */
                                      else
                                          tmp = 0.0;
                                  }
                                  else
                                  {   gspr = here->BSIM4v7gstot;
                                      gdpr = here->BSIM4v7gdtot;
                                      tmp = 0.0;
                                  }
                              }
                              else if(model->BSIM4v7tnoiMod == 1)
                              {   T5 = here->BSIM4v7Vgsteff / here->BSIM4v7EsatL;
                                  T5 *= T5;
                                  npart_beta = model->BSIM4v7rnoia * (1.0 + T5
                                             * model->BSIM4v7tnoia * pParam->BSIM4v7leff);
                                  npart_theta = model->BSIM4v7rnoib * (1.0 + T5
                                              * model->BSIM4v7tnoib * pParam->BSIM4v7leff);
                                  if(npart_theta > 0.9)
                                     npart_theta = 0.9;
                                  if(npart_theta > 0.9 * npart_beta)
                                     npart_theta = 0.9 * npart_beta; //4.6.2

                                  if (model->BSIM4v7rdsMod == 0)
                                  {   gspr = here->BSIM4v7sourceConductance;
                                      gdpr = here->BSIM4v7drainConductance;
                                  }
                                  else
                                  {   gspr = here->BSIM4v7gstot;
                                      gdpr = here->BSIM4v7gdtot;
                                  }

                                  if ((*(ckt->CKTstates[0] + here->BSIM4v7vds)) >= 0.0)
                                      gspr = gspr * (1.0 + npart_theta * npart_theta * gspr
                                           / here->BSIM4v7IdovVds);
                                  else
                                      gdpr = gdpr * (1.0 + npart_theta * npart_theta * gdpr
                                           / here->BSIM4v7IdovVds);
                              }
                              else
                              {   /* tnoiMod=2 (v4.7) */

                                  if (model->BSIM4v7rdsMod == 0)
                                  {   gspr = here->BSIM4v7sourceConductance;
                                      gdpr = here->BSIM4v7drainConductance;
                                  }
                                  else
                                  {   gspr = here->BSIM4v7gstot;
                                      gdpr = here->BSIM4v7gdtot;
                                  }

                              }

                              NevalSrc(&noizDens[BSIM4v7RDNOIZ],
                                       &lnNdens[BSIM4v7RDNOIZ], ckt, THERMNOISE,
                                       here->BSIM4v7dNodePrime, here->BSIM4v7dNode,
                                       gdpr * m);

                              NevalSrc(&noizDens[BSIM4v7RSNOIZ],
                                       &lnNdens[BSIM4v7RSNOIZ], ckt, THERMNOISE,
                                       here->BSIM4v7sNodePrime, here->BSIM4v7sNode,
                                       gspr * m);


                              if (here->BSIM4v7rgateMod == 1)
                              {   NevalSrc(&noizDens[BSIM4v7RGNOIZ],
                                       &lnNdens[BSIM4v7RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4v7gNodePrime, here->BSIM4v7gNodeExt,
                                       here->BSIM4v7grgeltd * m);
                              }
                              else if (here->BSIM4v7rgateMod == 2)
                              {
                                T0 = 1.0 + here->BSIM4v7grgeltd/here->BSIM4v7gcrg;
                                T1 = T0 * T0;
                                  NevalSrc(&noizDens[BSIM4v7RGNOIZ],
                                       &lnNdens[BSIM4v7RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4v7gNodePrime, here->BSIM4v7gNodeExt,
                                       here->BSIM4v7grgeltd * m / T1);
                              }
                              else if (here->BSIM4v7rgateMod == 3)
                              {   NevalSrc(&noizDens[BSIM4v7RGNOIZ],
                                       &lnNdens[BSIM4v7RGNOIZ], ckt, THERMNOISE,
                                       here->BSIM4v7gNodeMid, here->BSIM4v7gNodeExt,
                                       here->BSIM4v7grgeltd * m);
                              }
                              else
                              {    noizDens[BSIM4v7RGNOIZ] = 0.0;
                                   lnNdens[BSIM4v7RGNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RGNOIZ], N_MINLOG));
                              }

                                    bodymode = 5;
                                    if (here->BSIM4v7rbodyMod == 2)
                                    {        if( ( !model->BSIM4v7rbps0Given) ||
                                      ( !model->BSIM4v7rbpd0Given) )
                                             bodymode = 1;
                                           else
                                     if( (!model->BSIM4v7rbsbx0Given && !model->BSIM4v7rbsby0Given) ||
                                          (!model->BSIM4v7rbdbx0Given && !model->BSIM4v7rbdby0Given) )
                                             bodymode = 3;
                                }

                              if (here->BSIM4v7rbodyMod)
                              {
                                if(bodymode == 5)
                                  {
                                    NevalSrc(&noizDens[BSIM4v7RBPSNOIZ],
                                             &lnNdens[BSIM4v7RBPSNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7sbNode,
                                             here->BSIM4v7grbps * m);
                                    NevalSrc(&noizDens[BSIM4v7RBPDNOIZ],
                                             &lnNdens[BSIM4v7RBPDNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7dbNode,
                                             here->BSIM4v7grbpd * m);
                                    NevalSrc(&noizDens[BSIM4v7RBPBNOIZ],
                                             &lnNdens[BSIM4v7RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7bNode,
                                             here->BSIM4v7grbpb * m);
                                    NevalSrc(&noizDens[BSIM4v7RBSBNOIZ],
                                             &lnNdens[BSIM4v7RBSBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNode, here->BSIM4v7sbNode,
                                             here->BSIM4v7grbsb * m);
                                    NevalSrc(&noizDens[BSIM4v7RBDBNOIZ],
                                             &lnNdens[BSIM4v7RBDBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNode, here->BSIM4v7dbNode,
                                             here->BSIM4v7grbdb * m);
                                  }
                                if(bodymode == 3)
                                  {
                                    NevalSrc(&noizDens[BSIM4v7RBPSNOIZ],
                                             &lnNdens[BSIM4v7RBPSNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7sbNode,
                                             here->BSIM4v7grbps * m);
                                    NevalSrc(&noizDens[BSIM4v7RBPDNOIZ],
                                             &lnNdens[BSIM4v7RBPDNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7dbNode,
                                             here->BSIM4v7grbpd * m);
                                    NevalSrc(&noizDens[BSIM4v7RBPBNOIZ],
                                             &lnNdens[BSIM4v7RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7bNode,
                                             here->BSIM4v7grbpb * m);
                                     noizDens[BSIM4v7RBSBNOIZ] = noizDens[BSIM4v7RBDBNOIZ] = 0.0;
                                     lnNdens[BSIM4v7RBSBNOIZ] =
                                       log(MAX(noizDens[BSIM4v7RBSBNOIZ], N_MINLOG));
                                     lnNdens[BSIM4v7RBDBNOIZ] =
                                       log(MAX(noizDens[BSIM4v7RBDBNOIZ], N_MINLOG));
                                  }
                                if(bodymode == 1)
                                  {
                                    NevalSrc(&noizDens[BSIM4v7RBPBNOIZ],
                                             &lnNdens[BSIM4v7RBPBNOIZ], ckt, THERMNOISE,
                                             here->BSIM4v7bNodePrime, here->BSIM4v7bNode,
                                             here->BSIM4v7grbpb * m);
                                    noizDens[BSIM4v7RBPSNOIZ] = noizDens[BSIM4v7RBPDNOIZ] = 0.0;
                                    noizDens[BSIM4v7RBSBNOIZ] = noizDens[BSIM4v7RBDBNOIZ] = 0.0;
                                    lnNdens[BSIM4v7RBPSNOIZ] =
                                      log(MAX(noizDens[BSIM4v7RBPSNOIZ], N_MINLOG));
                                    lnNdens[BSIM4v7RBPDNOIZ] =
                                      log(MAX(noizDens[BSIM4v7RBPDNOIZ], N_MINLOG));
                                    lnNdens[BSIM4v7RBSBNOIZ] =
                                      log(MAX(noizDens[BSIM4v7RBSBNOIZ], N_MINLOG));
                                    lnNdens[BSIM4v7RBDBNOIZ] =
                                      log(MAX(noizDens[BSIM4v7RBDBNOIZ], N_MINLOG));
                                  }
                              }
                              else
                              {   noizDens[BSIM4v7RBPSNOIZ] = noizDens[BSIM4v7RBPDNOIZ] = 0.0;
                                  noizDens[BSIM4v7RBPBNOIZ] = 0.0;
                                  noizDens[BSIM4v7RBSBNOIZ] = noizDens[BSIM4v7RBDBNOIZ] = 0.0;
                                  lnNdens[BSIM4v7RBPSNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RBPSNOIZ], N_MINLOG));
                                  lnNdens[BSIM4v7RBPDNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RBPDNOIZ], N_MINLOG));
                                  lnNdens[BSIM4v7RBPBNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RBPBNOIZ], N_MINLOG));
                                  lnNdens[BSIM4v7RBSBNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RBSBNOIZ], N_MINLOG));
                                  lnNdens[BSIM4v7RBDBNOIZ] =
                                          log(MAX(noizDens[BSIM4v7RBDBNOIZ], N_MINLOG));
                              }

                              if(model->BSIM4v7tnoiMod == 2)
                              {
                                  eta = 1.0 - here->BSIM4v7Vdseff * here->BSIM4v7AbovVgst2Vtm;
                                  T0 = 1.0 - eta;
                                  T1 = 1.0 + eta;
                                  T2 = T1 + 2.0 * here->BSIM4v7Abulk * model->BSIM4v7vtm / here->BSIM4v7Vgsteff;
                                  Leff = pParam->BSIM4v7leff;
                                  Lvsat = Leff * (1.0 + here->BSIM4v7Vdseff / here->BSIM4v7EsatL);
                                  T6 = Leff / Lvsat;

                                  T5 = here->BSIM4v7Vgsteff / here->BSIM4v7EsatL;
                                  T5 = T5 * T5;
                                  gamma = T6 * (0.5 * T1 + T0 * T0 / (6.0 * T2));
                                  T3 = T2 * T2;
                                  T4 = T0 * T0;
                                  T5 = T3 * T3;
                                  delta = (T1 / T3 - (5.0 * T1 + T2) * T4 / (15.0 * T5) + T4 * T4 / (9.0 * T5 * T2)) / (6.0 * T6 * T6 * T6);
                                  T7 = T0 / T2;
                                  epsilon = (T7 - T7 * T7 * T7 / 3.0) / (6.0 * T6);

                                  T8 = here->BSIM4v7Vgsteff / here->BSIM4v7EsatL;
                                  T8 *= T8;
                                  npart_c = model->BSIM4v7rnoic * (1.0 + T8
                                          * model->BSIM4v7tnoic * Leff);
                                  ctnoi = epsilon / sqrt(gamma * delta)
                                      * (2.5316 * npart_c);

                                  npart_beta = model->BSIM4v7rnoia * (1.0 + T8
                                      * model->BSIM4v7tnoia * Leff);
                                  npart_theta = model->BSIM4v7rnoib * (1.0 + T8
                                      * model->BSIM4v7tnoib * Leff);
                                  gamma = gamma * (3.0 * npart_beta * npart_beta);
                                  delta = delta * (3.75 * npart_theta * npart_theta);

                                  GammaGd0 = gamma * here->BSIM4v7noiGd0;
                                  C0 = here->BSIM4v7Coxeff * pParam->BSIM4v7weffCV * here->BSIM4v7nf * pParam->BSIM4v7leffCV;
                                  T0 = C0 / here->BSIM4v7noiGd0;
                                  sigrat = T0 * sqrt(delta / gamma);
                              }
                              switch(model->BSIM4v7tnoiMod)
                              {  case 0:
                                      T0 = here->BSIM4v7ueff * fabs(here->BSIM4v7qinv);
                                      T1 = T0 * tmp + pParam->BSIM4v7leff
                                         * pParam->BSIM4v7leff;
                                      NevalSrc(&noizDens[BSIM4v7IDNOIZ],
                                               &lnNdens[BSIM4v7IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4v7dNodePrime,
                                               here->BSIM4v7sNodePrime,
                                               (T0 / T1) * model->BSIM4v7ntnoi * m);
                                      break;
                                 case 1:
                                      T0 = here->BSIM4v7gm + here->BSIM4v7gmbs + here->BSIM4v7gds;
                                      T0 *= T0;
                                      igsquare = npart_theta * npart_theta * T0 / here->BSIM4v7IdovVds;
                                      T1 = npart_beta * (here->BSIM4v7gm
                                         + here->BSIM4v7gmbs) + here->BSIM4v7gds;
                                      T2 = T1 * T1 / here->BSIM4v7IdovVds;
                                      NevalSrc(&noizDens[BSIM4v7IDNOIZ],
                                               &lnNdens[BSIM4v7IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4v7dNodePrime,
                                               here->BSIM4v7sNodePrime, (T2 - igsquare) * m);
                                      break;
                                  case 2:
                                      T2 = GammaGd0;
                                      T3 = ctnoi * ctnoi;
                                      T4 = 1.0 - T3;
                                      NevalSrc(&noizDens[BSIM4v7IDNOIZ],
                                               &lnNdens[BSIM4v7IDNOIZ], ckt,
                                               THERMNOISE, here->BSIM4v7dNodePrime,
                                               here->BSIM4v7sNodePrime, T2 * T4 * m);

                                     /* Evaluate output noise due to two correlated noise sources */
                                     omega = 2.0 * M_PI * data->freq;
                                     T5 = omega * sigrat;
                                     T6 = T5 * T5;
                                     T7 = T6 / (1.0 + T6);

                                     if (here->BSIM4v7mode >= 0)  {
                                         NevalSrc2(&noizDens[BSIM4v7CORLNOIZ],
                                               &lnNdens[BSIM4v7CORLNOIZ], ckt,
                                               THERMNOISE, here->BSIM4v7dNodePrime,
                                               here->BSIM4v7sNodePrime, T2 * T3 * m,
                                               here->BSIM4v7gNodePrime,
                                               here->BSIM4v7sNodePrime,
                                               T2 * T7 * m, 0.5 * M_PI);
                                     }
                                     else
                                     {
                                         NevalSrc2(&noizDens[BSIM4v7CORLNOIZ],
                                               &lnNdens[BSIM4v7CORLNOIZ], ckt,
                                               THERMNOISE, here->BSIM4v7sNodePrime,
                                               here->BSIM4v7dNodePrime, T2 * T3 * m,
                                               here->BSIM4v7gNodePrime,
                                               here->BSIM4v7dNodePrime,
                                               T2 * T7 * m, 0.5 * M_PI);
                                     }
                                     break;
                              }

                              NevalSrc(&noizDens[BSIM4v7FLNOIZ], (double*) NULL,
                                       ckt, N_GAIN, here->BSIM4v7dNodePrime,
                                       here->BSIM4v7sNodePrime, (double) 0.0);

                              switch(model->BSIM4v7fnoiMod)
                              {  case 0:
                                      noizDens[BSIM4v7FLNOIZ] *= m * model->BSIM4v7kf
                                            * exp(model->BSIM4v7af
                                            * log(MAX(fabs(here->BSIM4v7cd),
                                            N_MINLOG)))
                                            / (pow(data->freq, model->BSIM4v7ef)
                                            * pParam->BSIM4v7leff
                                            * pParam->BSIM4v7leff
                                            * model->BSIM4v7coxe);
                                      break;
                                 case 1:
                                      Vds = *(ckt->CKTstates[0] + here->BSIM4v7vds);
                                      if (Vds < 0.0)
                                          Vds = -Vds;

                                      Ssi = Eval1ovFNoise(Vds, model, here,
                                          data->freq, ckt->CKTtemp);
                                      T10 = model->BSIM4v7oxideTrapDensityA
                                          * CONSTboltz * ckt->CKTtemp;
                                      T11 = pParam->BSIM4v7weff * here->BSIM4v7nf * pParam->BSIM4v7leff
                                          * pow(data->freq, model->BSIM4v7ef) * 1.0e10
                                          * here->BSIM4v7nstar * here->BSIM4v7nstar;
                                      Swi = T10 / T11 * here->BSIM4v7cd
                                          * here->BSIM4v7cd;
                                      T1 = Swi + Ssi;
                                      if (T1 > 0.0)
                                          noizDens[BSIM4v7FLNOIZ] *= m * (Ssi * Swi) / T1;
                                      else
                                          noizDens[BSIM4v7FLNOIZ] *= 0.0;
                                      break;
                              }

                              lnNdens[BSIM4v7FLNOIZ] =
                                     log(MAX(noizDens[BSIM4v7FLNOIZ], N_MINLOG));


                        if(here->BSIM4v7mode >= 0) {  /* bugfix  */
                              NevalSrc(&noizDens[BSIM4v7IGSNOIZ],
                                   &lnNdens[BSIM4v7IGSNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4v7gNodePrime, here->BSIM4v7sNodePrime,
                                   m * (here->BSIM4v7Igs + here->BSIM4v7Igcs));
                              NevalSrc(&noizDens[BSIM4v7IGDNOIZ],
                                   &lnNdens[BSIM4v7IGDNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4v7gNodePrime, here->BSIM4v7dNodePrime,
                                   m * (here->BSIM4v7Igd + here->BSIM4v7Igcd));
                        } else {
                              NevalSrc(&noizDens[BSIM4v7IGSNOIZ],
                                   &lnNdens[BSIM4v7IGSNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4v7gNodePrime, here->BSIM4v7sNodePrime,
                                   m * (here->BSIM4v7Igs + here->BSIM4v7Igcd));
                              NevalSrc(&noizDens[BSIM4v7IGDNOIZ],
                                   &lnNdens[BSIM4v7IGDNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4v7gNodePrime, here->BSIM4v7dNodePrime,
                                   m * (here->BSIM4v7Igd + here->BSIM4v7Igcs));
                        }
                              NevalSrc(&noizDens[BSIM4v7IGBNOIZ],
                                   &lnNdens[BSIM4v7IGBNOIZ], ckt, SHOTNOISE,
                                   here->BSIM4v7gNodePrime, here->BSIM4v7bNodePrime,
                                   m * here->BSIM4v7Igb);


                              noizDens[BSIM4v7TOTNOIZ] = noizDens[BSIM4v7RDNOIZ]
                                     + noizDens[BSIM4v7RSNOIZ] + noizDens[BSIM4v7RGNOIZ]
                                     + noizDens[BSIM4v7RBPSNOIZ] + noizDens[BSIM4v7RBPDNOIZ]
                                     + noizDens[BSIM4v7RBPBNOIZ]
                                     + noizDens[BSIM4v7RBSBNOIZ] + noizDens[BSIM4v7RBDBNOIZ]
                                     + noizDens[BSIM4v7IDNOIZ] + noizDens[BSIM4v7FLNOIZ]
                                     + noizDens[BSIM4v7IGSNOIZ] + noizDens[BSIM4v7IGDNOIZ]
                                     + noizDens[BSIM4v7IGBNOIZ] + noizDens[BSIM4v7CORLNOIZ];
                              lnNdens[BSIM4v7TOTNOIZ] =
                                     log(MAX(noizDens[BSIM4v7TOTNOIZ], N_MINLOG));

                              *OnDens += noizDens[BSIM4v7TOTNOIZ];

                              if (data->delFreq == 0.0)
                              {   /* if we haven't done any previous
                                     integration, we need to initialize our
                                     "history" variables.
                                    */

                                  for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    here->BSIM4v7nVar[LNLSTDENS][i] =
                                             lnNdens[i];
                                  }

                                  /* clear out our integration variables
                                     if it's the first pass
                                   */
                                  if (data->freq ==
                                      job->NstartFreq)
                                  {   for (i = 0; i < BSIM4v7NSRCS; i++)
                                      {    here->BSIM4v7nVar[OUTNOIZ][i] = 0.0;
                                           here->BSIM4v7nVar[INNOIZ][i] = 0.0;
                                      }
                                  }
                              }
                              else
                              {   /* data->delFreq != 0.0,
                                     we have to integrate.
                                   */
                                  for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    if (i != BSIM4v7TOTNOIZ)
                                       {   tempOnoise = Nintegrate(noizDens[i],
                                                lnNdens[i],
                                                here->BSIM4v7nVar[LNLSTDENS][i],
                                                data);
                                           tempInoise = Nintegrate(noizDens[i]
                                                * data->GainSqInv, lnNdens[i]
                                                + data->lnGainInv,
                                                here->BSIM4v7nVar[LNLSTDENS][i]
                                                + data->lnGainInv, data);
                                           here->BSIM4v7nVar[LNLSTDENS][i] =
                                                lnNdens[i];
                                           data->outNoiz += tempOnoise;
                                           data->inNoise += tempInoise;
                                           if (job->NStpsSm != 0)
                                           {   here->BSIM4v7nVar[OUTNOIZ][i]
                                                     += tempOnoise;
                                               here->BSIM4v7nVar[OUTNOIZ][BSIM4v7TOTNOIZ]
                                                     += tempOnoise;
                                               here->BSIM4v7nVar[INNOIZ][i]
                                                     += tempInoise;
                                               here->BSIM4v7nVar[INNOIZ][BSIM4v7TOTNOIZ]
                                                     += tempInoise;
                                           }
                                       }
                                  }
                              }
                              if (data->prtSummary)
                              {   for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    /* print a summary report */
                                       data->outpVector[data->outNumber++]
                                             = noizDens[i];
                                  }
                              }
                              break;
                         case INT_NOIZ:
                              /* already calculated, just output */
                              if (job->NStpsSm != 0)
                              {   for (i = 0; i < BSIM4v7NSRCS; i++)
                                  {    data->outpVector[data->outNumber++]
                                             = here->BSIM4v7nVar[OUTNOIZ][i];
                                       data->outpVector[data->outNumber++]
                                             = here->BSIM4v7nVar[INNOIZ][i];
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
