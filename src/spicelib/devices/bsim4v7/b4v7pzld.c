/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4pzld.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 10/05/2001.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/complex.h"
#include "ngspice/sperror.h"
#include "bsim4v7def.h"
#include "ngspice/suffix.h"

int
BSIM4v7pzLoad(
GENmodel *inModel,
CKTcircuit *ckt,
SPcomplex *s)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;

double gjbd, gjbs, geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb;
double xcggb, xcgdb, xcgsb, xcgbb, xcbgb, xcbdb, xcbsb, xcbbb;
double xcdgb, xcddb, xcdsb, xcdbb, xcsgb, xcsdb, xcssb, xcsbb;
double gds, capbd, capbs, FwdSum, RevSum, Gm, Gmbs;
double gstot, gstotd, gstotg, gstots, gstotb, gspr;
double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
double gIstotg, gIstotd, gIstots, gIstotb;
double gIdtotg, gIdtotd, gIdtots, gIdtotb;
double gIbtotg, gIbtotd, gIbtots, gIbtotb;
double gIgtotg, gIgtotd, gIgtots, gIgtotb;
double cgso, cgdo, cgbo;
double xcdbdb=0.0, xcsbsb=0.0, xcgmgmb=0.0, xcgmdb=0.0, xcgmsb=0.0, xcdgmb=0.0, xcsgmb=0.0;
double xcgmbb=0.0, xcbgmb=0.0;
double dxpart, sxpart, xgtg, xgtd, xgts, xgtb, xcqgb=0.0, xcqdb=0.0, xcqsb=0.0, xcqbb=0.0;
double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
double T0=0.0, T1, CoxWL, qcheq, Cdg, Cdd, Cds, Csg, Csd, Css;
double ScalingFactor = 1.0e-9;
struct bsim4SizeDependParam *pParam;
double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;

double m;

    for (; model != NULL; model = model->BSIM4v7nextModel)
    {    for (here = model->BSIM4v7instances; here!= NULL;
              here = here->BSIM4v7nextInstance)
         {
              pParam = here->pParam;
              capbd = here->BSIM4v7capbd;
              capbs = here->BSIM4v7capbs;
              cgso = here->BSIM4v7cgso;
              cgdo = here->BSIM4v7cgdo;
              cgbo = pParam->BSIM4v7cgbo;

              if (here->BSIM4v7mode >= 0)
              {   Gm = here->BSIM4v7gm;
                  Gmbs = here->BSIM4v7gmbs;
                  FwdSum = Gm + Gmbs;
                  RevSum = 0.0;

                  gbbdp = -(here->BSIM4v7gbds);
                  gbbsp = here->BSIM4v7gbds + here->BSIM4v7gbgs + here->BSIM4v7gbbs;
                  gbdpg = here->BSIM4v7gbgs;
                  gbdpdp = here->BSIM4v7gbds;
                  gbdpb = here->BSIM4v7gbbs;
                  gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                  gbspdp = 0.0;
                  gbspg = 0.0;
                  gbspb = 0.0;
                  gbspsp = 0.0;

                  if (model->BSIM4v7igcMod)
                  {   gIstotg = here->BSIM4v7gIgsg + here->BSIM4v7gIgcsg;
                      gIstotd = here->BSIM4v7gIgcsd;
                      gIstots = here->BSIM4v7gIgss + here->BSIM4v7gIgcss;
                      gIstotb = here->BSIM4v7gIgcsb;

                      gIdtotg = here->BSIM4v7gIgdg + here->BSIM4v7gIgcdg;
                      gIdtotd = here->BSIM4v7gIgdd + here->BSIM4v7gIgcdd;
                      gIdtots = here->BSIM4v7gIgcds;
                      gIdtotb = here->BSIM4v7gIgcdb;
                  }
                  else
                  {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                      gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
                  }

                  if (model->BSIM4v7igbMod)
                  {   gIbtotg = here->BSIM4v7gIgbg;
                      gIbtotd = here->BSIM4v7gIgbd;
                      gIbtots = here->BSIM4v7gIgbs;
                      gIbtotb = here->BSIM4v7gIgbb;
                  }
                  else
                      gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                  if ((model->BSIM4v7igcMod != 0) || (model->BSIM4v7igbMod != 0))
                  {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                      gIgtots = gIstots + gIdtots + gIbtots;
                      gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  }
                  else
                      gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                  if (here->BSIM4v7rgateMod == 2)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4v7vges)
                         - *(ckt->CKTstates[0] + here->BSIM4v7vgs);
                  else if (here->BSIM4v7rgateMod == 3)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4v7vgms)
                         - *(ckt->CKTstates[0] + here->BSIM4v7vgs);
                  if (here->BSIM4v7rgateMod > 1)
                  {   gcrgd = here->BSIM4v7gcrgd * T0;
                      gcrgg = here->BSIM4v7gcrgg * T0;
                      gcrgs = here->BSIM4v7gcrgs * T0;
                      gcrgb = here->BSIM4v7gcrgb * T0;
                      gcrgg -= here->BSIM4v7gcrg;
                      gcrg = here->BSIM4v7gcrg;
                  }
                  else
                      gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                  if (here->BSIM4v7acnqsMod == 0)
                  {   if (here->BSIM4v7rgateMod == 3)
                      {   xcgmgmb = cgdo + cgso + pParam->BSIM4v7cgbo;
                          xcgmdb = -cgdo;
                          xcgmsb = -cgso;
                          xcgmbb = -pParam->BSIM4v7cgbo;

                          xcdgmb = xcgmdb;
                          xcsgmb = xcgmsb;
                          xcbgmb = xcgmbb;

                          xcggb = here->BSIM4v7cggb;
                          xcgdb = here->BSIM4v7cgdb;
                          xcgsb = here->BSIM4v7cgsb;
                          xcgbb = -(xcggb + xcgdb + xcgsb);

                          xcdgb = here->BSIM4v7cdgb;
                          xcsgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                                + here->BSIM4v7cdgb);
                          xcbgb = here->BSIM4v7cbgb;
                      }
                      else
                      {   xcggb = here->BSIM4v7cggb + cgdo + cgso
                                + pParam->BSIM4v7cgbo;
                          xcgdb = here->BSIM4v7cgdb - cgdo;
                          xcgsb = here->BSIM4v7cgsb - cgso;
                          xcgbb = -(xcggb + xcgdb + xcgsb);

                          xcdgb = here->BSIM4v7cdgb - cgdo;
                          xcsgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                                + here->BSIM4v7cdgb + cgso);
                          xcbgb = here->BSIM4v7cbgb - pParam->BSIM4v7cgbo;

                          xcdgmb = xcsgmb = xcbgmb = 0.0;
                      }
                      xcddb = here->BSIM4v7cddb + here->BSIM4v7capbd + cgdo;
                      xcdsb = here->BSIM4v7cdsb;

                      xcsdb = -(here->BSIM4v7cgdb + here->BSIM4v7cbdb
                            + here->BSIM4v7cddb);
                      xcssb = here->BSIM4v7capbs + cgso - (here->BSIM4v7cgsb
                            + here->BSIM4v7cbsb + here->BSIM4v7cdsb);

                      if (!here->BSIM4v7rbodyMod)
                      {   xcdbb = -(xcdgb + xcddb + xcdsb + xcdgmb);
                          xcsbb = -(xcsgb + xcsdb + xcssb + xcsgmb);
                          xcbdb = here->BSIM4v7cbdb - here->BSIM4v7capbd;
                          xcbsb = here->BSIM4v7cbsb - here->BSIM4v7capbs;
                          xcdbdb = 0.0;
                      }
                      else
                      {   xcdbb  = -(here->BSIM4v7cddb + here->BSIM4v7cdgb
                                 + here->BSIM4v7cdsb);
                          xcsbb = -(xcsgb + xcsdb + xcssb + xcsgmb)
                                + here->BSIM4v7capbs;
                          xcbdb = here->BSIM4v7cbdb;
                          xcbsb = here->BSIM4v7cbsb;

                          xcdbdb = -here->BSIM4v7capbd;
                          xcsbsb = -here->BSIM4v7capbs;
                      }
                      xcbbb = -(xcbdb + xcbgb + xcbsb + xcbgmb);

                      xgtg = xgtd = xgts = xgtb = 0.0;
                      sxpart = 0.6;
                      dxpart = 0.4;
                      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                  = ddxpart_dVs = 0.0;
                      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                  = dsxpart_dVs = 0.0;
                  }
                  else
                  {   xcggb = xcgdb = xcgsb = xcgbb = 0.0;
                      xcbgb = xcbdb = xcbsb = xcbbb = 0.0;
                      xcdgb = xcddb = xcdsb = xcdbb = 0.0;
                      xcsgb = xcsdb = xcssb = xcsbb = 0.0;

                      xgtg = here->BSIM4v7gtg;
                      xgtd = here->BSIM4v7gtd;
                      xgts = here->BSIM4v7gts;
                      xgtb = here->BSIM4v7gtb;

                      xcqgb = here->BSIM4v7cqgb;
                      xcqdb = here->BSIM4v7cqdb;
                      xcqsb = here->BSIM4v7cqsb;
                      xcqbb = here->BSIM4v7cqbb;

                      CoxWL = model->BSIM4v7coxe * here->pParam->BSIM4v7weffCV
                            * here->BSIM4v7nf * here->pParam->BSIM4v7leffCV;
                      qcheq = -(here->BSIM4v7qgate + here->BSIM4v7qbulk);
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
                      {   dxpart = here->BSIM4v7qdrn / qcheq;
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

                          ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg
                                      + ddxpart_dVs);
                      }
                      sxpart = 1.0 - dxpart;
                      dsxpart_dVd = -ddxpart_dVd;
                      dsxpart_dVg = -ddxpart_dVg;
                      dsxpart_dVs = -ddxpart_dVs;
                      dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                  }
              }
              else
              {   Gm = -here->BSIM4v7gm;
                  Gmbs = -here->BSIM4v7gmbs;
                  FwdSum = 0.0;
                  RevSum = -(Gm + Gmbs);

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

                      gIdtotg = here->BSIM4v7gIgdg + here->BSIM4v7gIgcsg;
                      gIdtotd = here->BSIM4v7gIgdd + here->BSIM4v7gIgcss;
                      gIdtots = here->BSIM4v7gIgcsd;
                      gIdtotb = here->BSIM4v7gIgcsb;
                  }
                  else
                  {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                      gIdtotg = gIdtotd = gIdtots = gIdtotb  = 0.0;
                  }

                  if (model->BSIM4v7igbMod)
                  {   gIbtotg = here->BSIM4v7gIgbg;
                      gIbtotd = here->BSIM4v7gIgbs;
                      gIbtots = here->BSIM4v7gIgbd;
                      gIbtotb = here->BSIM4v7gIgbb;
                  }
                  else
                      gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                  if ((model->BSIM4v7igcMod != 0) || (model->BSIM4v7igbMod != 0))
                  {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                      gIgtots = gIstots + gIdtots + gIbtots;
                      gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  }
                  else
                      gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                  if (here->BSIM4v7rgateMod == 2)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4v7vges)
                         - *(ckt->CKTstates[0] + here->BSIM4v7vgs);
                  else if (here->BSIM4v7rgateMod == 3)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4v7vgms)
                         - *(ckt->CKTstates[0] + here->BSIM4v7vgs);
                  if (here->BSIM4v7rgateMod > 1)
                  {   gcrgd = here->BSIM4v7gcrgs * T0;
                      gcrgg = here->BSIM4v7gcrgg * T0;
                      gcrgs = here->BSIM4v7gcrgd * T0;
                      gcrgb = here->BSIM4v7gcrgb * T0;
                      gcrgg -= here->BSIM4v7gcrg;
                      gcrg = here->BSIM4v7gcrg;
                  }
                  else
                      gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                  if (here->BSIM4v7acnqsMod == 0)
                  {   if (here->BSIM4v7rgateMod == 3)
                      {   xcgmgmb = cgdo + cgso + pParam->BSIM4v7cgbo;
                          xcgmdb = -cgdo;
                          xcgmsb = -cgso;
                          xcgmbb = -pParam->BSIM4v7cgbo;

                          xcdgmb = xcgmdb;
                          xcsgmb = xcgmsb;
                          xcbgmb = xcgmbb;

                          xcggb = here->BSIM4v7cggb;
                          xcgdb = here->BSIM4v7cgsb;
                          xcgsb = here->BSIM4v7cgdb;
                          xcgbb = -(xcggb + xcgdb + xcgsb);

                          xcdgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                                + here->BSIM4v7cdgb);
                          xcsgb = here->BSIM4v7cdgb;
                          xcbgb = here->BSIM4v7cbgb;
                      }
                      else
                      {   xcggb = here->BSIM4v7cggb + cgdo + cgso
                                + pParam->BSIM4v7cgbo;
                          xcgdb = here->BSIM4v7cgsb - cgdo;
                          xcgsb = here->BSIM4v7cgdb - cgso;
                          xcgbb = -(xcggb + xcgdb + xcgsb);

                          xcdgb = -(here->BSIM4v7cggb + here->BSIM4v7cbgb
                                + here->BSIM4v7cdgb + cgdo);
                          xcsgb = here->BSIM4v7cdgb - cgso;
                          xcbgb = here->BSIM4v7cbgb - pParam->BSIM4v7cgbo;

                          xcdgmb = xcsgmb = xcbgmb = 0.0;
                      }
                      xcddb = here->BSIM4v7capbd + cgdo - (here->BSIM4v7cgsb
                            + here->BSIM4v7cbsb + here->BSIM4v7cdsb);
                      xcdsb = -(here->BSIM4v7cgdb + here->BSIM4v7cbdb
                            + here->BSIM4v7cddb);

                      xcsdb = here->BSIM4v7cdsb;
                      xcssb = here->BSIM4v7cddb + here->BSIM4v7capbs + cgso;

                      if (!here->BSIM4v7rbodyMod)
                      {   xcdbb = -(xcdgb + xcddb + xcdsb + xcdgmb);
                          xcsbb = -(xcsgb + xcsdb + xcssb + xcsgmb);
                          xcbdb = here->BSIM4v7cbsb - here->BSIM4v7capbd;
                          xcbsb = here->BSIM4v7cbdb - here->BSIM4v7capbs;
                          xcdbdb = 0.0;
                      }
                      else
                      {   xcdbb = -(xcdgb + xcddb + xcdsb + xcdgmb)
                                + here->BSIM4v7capbd;
                          xcsbb = -(here->BSIM4v7cddb + here->BSIM4v7cdgb
                                + here->BSIM4v7cdsb);
                          xcbdb = here->BSIM4v7cbsb;
                          xcbsb = here->BSIM4v7cbdb;
                          xcdbdb = -here->BSIM4v7capbd;
                          xcsbsb = -here->BSIM4v7capbs;
                      }
                      xcbbb = -(xcbgb + xcbdb + xcbsb + xcbgmb);

                      xgtg = xgtd = xgts = xgtb = 0.0;
                      sxpart = 0.4;
                      dxpart = 0.6;
                      ddxpart_dVd = ddxpart_dVg = ddxpart_dVb
                                  = ddxpart_dVs = 0.0;
                      dsxpart_dVd = dsxpart_dVg = dsxpart_dVb
                                  = dsxpart_dVs = 0.0;
                  }
                  else
                  {   xcggb = xcgdb = xcgsb = xcgbb = 0.0;
                      xcbgb = xcbdb = xcbsb = xcbbb = 0.0;
                      xcdgb = xcddb = xcdsb = xcdbb = 0.0;
                      xcsgb = xcsdb = xcssb = xcsbb = 0.0;

                      xgtg = here->BSIM4v7gtg;
                      xgtd = here->BSIM4v7gts;
                      xgts = here->BSIM4v7gtd;
                      xgtb = here->BSIM4v7gtb;

                      xcqgb = here->BSIM4v7cqgb;
                      xcqdb = here->BSIM4v7cqsb;
                      xcqsb = here->BSIM4v7cqdb;
                      xcqbb = here->BSIM4v7cqbb;

                      CoxWL = model->BSIM4v7coxe * here->pParam->BSIM4v7weffCV
                            * here->BSIM4v7nf * here->pParam->BSIM4v7leffCV;
                      qcheq = -(here->BSIM4v7qgate + here->BSIM4v7qbulk);
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
                      {   sxpart = here->BSIM4v7qdrn / qcheq;
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

                          dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg
                                      + dsxpart_dVs);
                      }
                      dxpart = 1.0 - sxpart;
                      ddxpart_dVd = -dsxpart_dVd;
                      ddxpart_dVg = -dsxpart_dVg;
                      ddxpart_dVs = -dsxpart_dVs;
                      ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                  }
              }

              if (model->BSIM4v7rdsMod == 1)
              {   gstot = here->BSIM4v7gstot;
                  gstotd = here->BSIM4v7gstotd;
                  gstotg = here->BSIM4v7gstotg;
                  gstots = here->BSIM4v7gstots - gstot;
                  gstotb = here->BSIM4v7gstotb;

                  gdtot = here->BSIM4v7gdtot;
                  gdtotd = here->BSIM4v7gdtotd - gdtot;
                  gdtotg = here->BSIM4v7gdtotg;
                  gdtots = here->BSIM4v7gdtots;
                  gdtotb = here->BSIM4v7gdtotb;
              }
              else
              {   gstot = gstotd = gstotg = gstots = gstotb = 0.0;
                  gdtot = gdtotd = gdtotg = gdtots = gdtotb = 0.0;
              }


              T1 = *(ckt->CKTstate0 + here->BSIM4v7qdef) * here->BSIM4v7gtau;
              gds = here->BSIM4v7gds;

              /*
               * Loading PZ matrix
               */
              m = here->BSIM4v7m;

              if (!model->BSIM4v7rdsMod)
              {   gdpr = here->BSIM4v7drainConductance;
                  gspr = here->BSIM4v7sourceConductance;
              }
              else
                  gdpr = gspr = 0.0;

              if (!here->BSIM4v7rbodyMod)
              {   gjbd = here->BSIM4v7gbd;
                  gjbs = here->BSIM4v7gbs;
              }
              else
                  gjbd = gjbs = 0.0;

              geltd = here->BSIM4v7grgeltd;

              if (here->BSIM4v7rgateMod == 1)
              {   *(here->BSIM4v7GEgePtr) += m * geltd;
                  *(here->BSIM4v7GPgePtr) -= m * geltd;
                  *(here->BSIM4v7GEgpPtr) -= m * geltd;

                  *(here->BSIM4v7GPgpPtr ) += m * xcggb * s->real;
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggb * s->imag;
                  *(here->BSIM4v7GPgpPtr) += m * (geltd - xgtg + gIgtotg);
                  *(here->BSIM4v7GPdpPtr ) += m * xcgdb * s->real;
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdb * s->imag;
                  *(here->BSIM4v7GPdpPtr) -= m * (xgtd - gIgtotd);
                  *(here->BSIM4v7GPspPtr ) += m * xcgsb * s->real;
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsb * s->imag;
                  *(here->BSIM4v7GPspPtr) -= m * (xgts - gIgtots);
                  *(here->BSIM4v7GPbpPtr ) += m * xcgbb * s->real;
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbb * s->imag;
                  *(here->BSIM4v7GPbpPtr) -= m * (xgtb - gIgtotb);
              }
              else if (here->BSIM4v7rgateMod == 2)
              {   *(here->BSIM4v7GEgePtr) += m * gcrg;
                  *(here->BSIM4v7GEgpPtr) += m * gcrgg;
                  *(here->BSIM4v7GEdpPtr) += m * gcrgd;
                  *(here->BSIM4v7GEspPtr) += m * gcrgs;
                  *(here->BSIM4v7GEbpPtr) += m * gcrgb;

                  *(here->BSIM4v7GPgePtr) -= m * gcrg;
                  *(here->BSIM4v7GPgpPtr ) += m * xcggb * s->real;
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggb * s->imag;
                  *(here->BSIM4v7GPgpPtr) -= m * (gcrgg + xgtg - gIgtotg);
                  *(here->BSIM4v7GPdpPtr ) += m * xcgdb * s->real;
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdb * s->imag;
                  *(here->BSIM4v7GPdpPtr) -= m * (gcrgd + xgtd - gIgtotd);
                  *(here->BSIM4v7GPspPtr ) += m * xcgsb * s->real;
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsb * s->imag;
                  *(here->BSIM4v7GPspPtr) -= m * (gcrgs + xgts - gIgtots);
                  *(here->BSIM4v7GPbpPtr ) += m * xcgbb * s->real;
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbb * s->imag;
                  *(here->BSIM4v7GPbpPtr) -= m * (gcrgb + xgtb - gIgtotb);
              }
              else if (here->BSIM4v7rgateMod == 3)
              {   *(here->BSIM4v7GEgePtr) += m * geltd;
                  *(here->BSIM4v7GEgmPtr) -= m * geltd;
                  *(here->BSIM4v7GMgePtr) -= m * geltd;
                  *(here->BSIM4v7GMgmPtr) += m * (geltd + gcrg);
                  *(here->BSIM4v7GMgmPtr ) += m * xcgmgmb * s->real;
                  *(here->BSIM4v7GMgmPtr +1) += m * xcgmgmb * s->imag;
  
                  *(here->BSIM4v7GMdpPtr) += m * gcrgd;
                  *(here->BSIM4v7GMdpPtr ) += m * xcgmdb * s->real;
                  *(here->BSIM4v7GMdpPtr +1) += m * xcgmdb * s->imag;
                  *(here->BSIM4v7GMgpPtr) += m * gcrgg;
                  *(here->BSIM4v7GMspPtr) += m * gcrgs;
                  *(here->BSIM4v7GMspPtr ) += m * xcgmsb * s->real;
                  *(here->BSIM4v7GMspPtr +1) += m * xcgmsb * s->imag;
                  *(here->BSIM4v7GMbpPtr) += m * gcrgb;
                  *(here->BSIM4v7GMbpPtr ) += m * xcgmbb * s->real;
                  *(here->BSIM4v7GMbpPtr +1) += m * xcgmbb * s->imag;
  
                  *(here->BSIM4v7DPgmPtr ) += m * xcdgmb * s->real;
                  *(here->BSIM4v7DPgmPtr +1) += m * xcdgmb * s->imag;
                  *(here->BSIM4v7GPgmPtr) -= m * gcrg;
                  *(here->BSIM4v7SPgmPtr ) += m * xcsgmb * s->real;
                  *(here->BSIM4v7SPgmPtr +1) += m * xcsgmb * s->imag;
                  *(here->BSIM4v7BPgmPtr ) += m * xcbgmb * s->real;
                  *(here->BSIM4v7BPgmPtr +1) += m * xcbgmb * s->imag;
  
                  *(here->BSIM4v7GPgpPtr) -= m * (gcrgg + xgtg - gIgtotg);
                  *(here->BSIM4v7GPgpPtr ) += m * xcggb * s->real;
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggb * s->imag;
                  *(here->BSIM4v7GPdpPtr) -= m * (gcrgd + xgtd - gIgtotd);
                  *(here->BSIM4v7GPdpPtr ) += m * xcgdb * s->real;
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdb * s->imag;
                  *(here->BSIM4v7GPspPtr) -= m * (gcrgs + xgts - gIgtots);
                  *(here->BSIM4v7GPspPtr ) += m * xcgsb * s->real;
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsb * s->imag;
                  *(here->BSIM4v7GPbpPtr) -= m * (gcrgb + xgtb - gIgtotb);
                  *(here->BSIM4v7GPbpPtr ) += m * xcgbb * s->real;
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbb * s->imag;
              }
              else
              {   *(here->BSIM4v7GPdpPtr ) += m * xcgdb * s->real;
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdb * s->imag;
                  *(here->BSIM4v7GPdpPtr) -= m * (xgtd - gIgtotd);
                  *(here->BSIM4v7GPgpPtr ) += m * xcggb * s->real;
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggb * s->imag;
                  *(here->BSIM4v7GPgpPtr) -= m * (xgtg - gIgtotg);
                  *(here->BSIM4v7GPspPtr ) += m * xcgsb * s->real;
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsb * s->imag;
                  *(here->BSIM4v7GPspPtr) -= m * (xgts - gIgtots);
                  *(here->BSIM4v7GPbpPtr ) += m * xcgbb * s->real;
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbb * s->imag;
                  *(here->BSIM4v7GPbpPtr) -= m * (xgtb - gIgtotb);
              }

              if (model->BSIM4v7rdsMod)
              {   (*(here->BSIM4v7DgpPtr) += m * gdtotg);
                  (*(here->BSIM4v7DspPtr) += m * gdtots);
                  (*(here->BSIM4v7DbpPtr) += m * gdtotb);
                  (*(here->BSIM4v7SdpPtr) += m * gstotd);
                  (*(here->BSIM4v7SgpPtr) += m * gstotg);
                  (*(here->BSIM4v7SbpPtr) += m * gstotb);
              }

              *(here->BSIM4v7DPdpPtr ) += m * xcddb * s->real;
              *(here->BSIM4v7DPdpPtr +1) += m * xcddb * s->imag;
              *(here->BSIM4v7DPdpPtr) += m * (gdpr + gds + here->BSIM4v7gbd
                                     - gdtotd + RevSum + gbdpdp - gIdtotd
                                     + dxpart * xgtd + T1 * ddxpart_dVd);
              *(here->BSIM4v7DPdPtr) -= m * (gdpr + gdtot);
              *(here->BSIM4v7DPgpPtr ) += m * xcdgb * s->real;
              *(here->BSIM4v7DPgpPtr +1) += m * xcdgb * s->imag;
              *(here->BSIM4v7DPgpPtr) += m * (Gm - gdtotg + gbdpg - gIdtotg
                                     + T1 * ddxpart_dVg + dxpart * xgtg);
              *(here->BSIM4v7DPspPtr ) += m * xcdsb * s->real;
              *(here->BSIM4v7DPspPtr +1) += m * xcdsb * s->imag;
              *(here->BSIM4v7DPspPtr) -= m * (gds + FwdSum + gdtots - gbdpsp + gIdtots
                                     - T1 * ddxpart_dVs - dxpart * xgts);
              *(here->BSIM4v7DPbpPtr ) += m * xcdbb * s->real;
              *(here->BSIM4v7DPbpPtr +1) += m * xcdbb * s->imag;
              *(here->BSIM4v7DPbpPtr) -= m * (gjbd + gdtotb - Gmbs - gbdpb + gIdtotb
                                     - T1 * ddxpart_dVb - dxpart * xgtb);

              *(here->BSIM4v7DdpPtr) -= m * (gdpr - gdtotd);
              *(here->BSIM4v7DdPtr) += m * (gdpr + gdtot);

              *(here->BSIM4v7SPdpPtr ) += m * xcsdb * s->real;
              *(here->BSIM4v7SPdpPtr +1) += m * xcsdb * s->imag;
              *(here->BSIM4v7SPdpPtr) -= m * (gds + gstotd + RevSum - gbspdp + gIstotd
                                     - T1 * dsxpart_dVd - sxpart * xgtd);
              *(here->BSIM4v7SPgpPtr ) += m * xcsgb * s->real;
              *(here->BSIM4v7SPgpPtr +1) += m * xcsgb * s->imag;
              *(here->BSIM4v7SPgpPtr) -= m * (Gm + gstotg - gbspg + gIstotg
                                     - T1 * dsxpart_dVg - sxpart * xgtg);
              *(here->BSIM4v7SPspPtr ) += m * xcssb * s->real;
              *(here->BSIM4v7SPspPtr +1) += m * xcssb * s->imag;
              *(here->BSIM4v7SPspPtr) += m * (gspr + gds + here->BSIM4v7gbs - gIstots
                                     - gstots + FwdSum + gbspsp
                                     + sxpart * xgts + T1 * dsxpart_dVs);
              *(here->BSIM4v7SPsPtr) -= m * (gspr + gstot);
              *(here->BSIM4v7SPbpPtr ) += m * xcsbb * s->real;
              *(here->BSIM4v7SPbpPtr +1) += m * xcsbb * s->imag;
              *(here->BSIM4v7SPbpPtr) -= m * (gjbs + gstotb + Gmbs - gbspb + gIstotb
                                     - T1 * dsxpart_dVb - sxpart * xgtb);

              *(here->BSIM4v7SspPtr) -= m * (gspr - gstots);
              *(here->BSIM4v7SsPtr) += m * (gspr + gstot);

              *(here->BSIM4v7BPdpPtr ) += m * xcbdb * s->real;
              *(here->BSIM4v7BPdpPtr +1) += m * xcbdb * s->imag;
              *(here->BSIM4v7BPdpPtr) -= m * (gjbd - gbbdp + gIbtotd);
              *(here->BSIM4v7BPgpPtr ) += m * xcbgb * s->real;
              *(here->BSIM4v7BPgpPtr +1) += m * xcbgb * s->imag;
              *(here->BSIM4v7BPgpPtr) -= m * (here->BSIM4v7gbgs + gIbtotg);
              *(here->BSIM4v7BPspPtr ) += m * xcbsb * s->real;
              *(here->BSIM4v7BPspPtr +1) += m * xcbsb * s->imag;
              *(here->BSIM4v7BPspPtr) -= m * (gjbs - gbbsp + gIbtots);
              *(here->BSIM4v7BPbpPtr ) += m * xcbbb * s->real;
              *(here->BSIM4v7BPbpPtr +1) += m * xcbbb * s->imag;
              *(here->BSIM4v7BPbpPtr) += m * (gjbd + gjbs - here->BSIM4v7gbbs
                                     - gIbtotb);
              ggidld = here->BSIM4v7ggidld;
              ggidlg = here->BSIM4v7ggidlg;
              ggidlb = here->BSIM4v7ggidlb;
              ggislg = here->BSIM4v7ggislg;
              ggisls = here->BSIM4v7ggisls;
              ggislb = here->BSIM4v7ggislb;

              /* stamp gidl */
              (*(here->BSIM4v7DPdpPtr) += m * ggidld);
              (*(here->BSIM4v7DPgpPtr) += m * ggidlg);
              (*(here->BSIM4v7DPspPtr) -= m * ((ggidlg + ggidld) + ggidlb));
              (*(here->BSIM4v7DPbpPtr) += m * ggidlb);
              (*(here->BSIM4v7BPdpPtr) -= m * ggidld);
              (*(here->BSIM4v7BPgpPtr) -= m * ggidlg);
              (*(here->BSIM4v7BPspPtr) += m * ((ggidlg + ggidld) + ggidlb));
              (*(here->BSIM4v7BPbpPtr) -= m * ggidlb);
               /* stamp gisl */
              (*(here->BSIM4v7SPdpPtr) -= m * ((ggisls + ggislg) + ggislb));
              (*(here->BSIM4v7SPgpPtr) += m * ggislg);
              (*(here->BSIM4v7SPspPtr) += m * ggisls);
              (*(here->BSIM4v7SPbpPtr) += m * ggislb);
              (*(here->BSIM4v7BPdpPtr) += m * ((ggislg + ggisls) + ggislb));
              (*(here->BSIM4v7BPgpPtr) -= m * ggislg);
              (*(here->BSIM4v7BPspPtr) -= m * ggisls);
              (*(here->BSIM4v7BPbpPtr) -= m * ggislb);

              if (here->BSIM4v7rbodyMod)
              {   (*(here->BSIM4v7DPdbPtr ) += m * xcdbdb * s->real);
                  (*(here->BSIM4v7DPdbPtr +1) += m * xcdbdb * s->imag);
                  (*(here->BSIM4v7DPdbPtr) -= m * here->BSIM4v7gbd);
                  (*(here->BSIM4v7SPsbPtr ) += m * xcsbsb * s->real);
                  (*(here->BSIM4v7SPsbPtr +1) += m * xcsbsb * s->imag);
                  (*(here->BSIM4v7SPsbPtr) -= m * here->BSIM4v7gbs);

                  (*(here->BSIM4v7DBdpPtr ) += m * xcdbdb * s->real);
                  (*(here->BSIM4v7DBdpPtr +1) += m * xcdbdb * s->imag);
                  (*(here->BSIM4v7DBdpPtr) -= m * here->BSIM4v7gbd);
                  (*(here->BSIM4v7DBdbPtr ) -= m * xcdbdb * s->real);
                  (*(here->BSIM4v7DBdbPtr +1) -= m * xcdbdb * s->imag);
                  (*(here->BSIM4v7DBdbPtr) += m * (here->BSIM4v7gbd + here->BSIM4v7grbpd
                                          + here->BSIM4v7grbdb));
                  (*(here->BSIM4v7DBbpPtr) -= m * here->BSIM4v7grbpd);
                  (*(here->BSIM4v7DBbPtr) -= m * here->BSIM4v7grbdb);

                  (*(here->BSIM4v7BPdbPtr) -= m * here->BSIM4v7grbpd);
                  (*(here->BSIM4v7BPbPtr) -= m * here->BSIM4v7grbpb);
                  (*(here->BSIM4v7BPsbPtr) -= m * here->BSIM4v7grbps);
                  (*(here->BSIM4v7BPbpPtr) += m * (here->BSIM4v7grbpd + here->BSIM4v7grbps
                                          + here->BSIM4v7grbpb));
                  /* WDL: (-here->BSIM4v7gbbs) already added to BPbpPtr */

                  (*(here->BSIM4v7SBspPtr ) += m * xcsbsb * s->real);
                  (*(here->BSIM4v7SBspPtr +1) += m * xcsbsb * s->imag);
                  (*(here->BSIM4v7SBspPtr) -= m * here->BSIM4v7gbs);
                  (*(here->BSIM4v7SBbpPtr) -= m * here->BSIM4v7grbps);
                  (*(here->BSIM4v7SBbPtr) -= m * here->BSIM4v7grbsb);
                  (*(here->BSIM4v7SBsbPtr ) -= m * xcsbsb * s->real);
                  (*(here->BSIM4v7SBsbPtr +1) -= m * xcsbsb * s->imag);
                  (*(here->BSIM4v7SBsbPtr) += m * (here->BSIM4v7gbs
                                          + here->BSIM4v7grbps + here->BSIM4v7grbsb));

                  (*(here->BSIM4v7BdbPtr) -= m * here->BSIM4v7grbdb);
                  (*(here->BSIM4v7BbpPtr) -= m * here->BSIM4v7grbpb);
                  (*(here->BSIM4v7BsbPtr) -= m * here->BSIM4v7grbsb);
                  (*(here->BSIM4v7BbPtr) += m * (here->BSIM4v7grbsb + here->BSIM4v7grbdb
                                        + here->BSIM4v7grbpb));
              }

              if (here->BSIM4v7acnqsMod)
              {   *(here->BSIM4v7QqPtr ) += m * s->real * ScalingFactor;
                  *(here->BSIM4v7QqPtr +1) += m * s->imag * ScalingFactor;
                  *(here->BSIM4v7QgpPtr ) -= m * xcqgb * s->real;
                  *(here->BSIM4v7QgpPtr +1) -= m * xcqgb * s->imag;
                  *(here->BSIM4v7QdpPtr ) -= m * xcqdb * s->real;
                  *(here->BSIM4v7QdpPtr +1) -= m * xcqdb * s->imag;
                  *(here->BSIM4v7QbpPtr ) -= m * xcqbb * s->real;
                  *(here->BSIM4v7QbpPtr +1) -= m * xcqbb * s->imag;
                  *(here->BSIM4v7QspPtr ) -= m * xcqsb * s->real;
                  *(here->BSIM4v7QspPtr +1) -= m * xcqsb * s->imag;

                  *(here->BSIM4v7GPqPtr) -= m * here->BSIM4v7gtau;
                  *(here->BSIM4v7DPqPtr) += m * dxpart * here->BSIM4v7gtau;
                  *(here->BSIM4v7SPqPtr) += m * sxpart * here->BSIM4v7gtau;

                  *(here->BSIM4v7QqPtr) += m * here->BSIM4v7gtau;
                  *(here->BSIM4v7QgpPtr) += m * xgtg;
                  *(here->BSIM4v7QdpPtr) += m * xgtd;
                  *(here->BSIM4v7QbpPtr) += m * xgtb;
                  *(here->BSIM4v7QspPtr) += m * xgts;
              }
         }
    }
    return(OK);
}
