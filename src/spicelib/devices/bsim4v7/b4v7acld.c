/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4acld.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 10/05/2001.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h."
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
BSIM4v7acLoad(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;

double gjbd, gjbs, geltd, gcrg, gcrgg, gcrgd, gcrgs, gcrgb;
double xcbgb, xcbdb, xcbsb, xcbbb;
double xcggbr, xcgdbr, xcgsbr, xcgbbr, xcggbi, xcgdbi, xcgsbi, xcgbbi;
double Cggr, Cgdr, Cgsr, Cgbr, Cggi, Cgdi, Cgsi, Cgbi;
double xcddbr, xcdgbr, xcdsbr, xcdbbr, xcsdbr, xcsgbr, xcssbr, xcsbbr;
double xcddbi, xcdgbi, xcdsbi, xcdbbi, xcsdbi, xcsgbi, xcssbi, xcsbbi;
double xcdbdb, xcsbsb=0.0, xcgmgmb=0.0, xcgmdb=0.0, xcgmsb=0.0, xcdgmb, xcsgmb;
double xcgmbb=0.0, xcbgmb;
double capbd, capbs, omega;
double gstot, gstotd, gstotg, gstots, gstotb, gspr;
double gdtot, gdtotd, gdtotg, gdtots, gdtotb, gdpr;
double gIstotg, gIstotd, gIstots, gIstotb;
double gIdtotg, gIdtotd, gIdtots, gIdtotb;
double gIbtotg, gIbtotd, gIbtots, gIbtotb;
double gIgtotg, gIgtotd, gIgtots, gIgtotb;
double cgso, cgdo, cgbo;
double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
double T0=0.0, T1, T2, T3;
double Csg, Csd, Css;
double Cdgr, Cddr, Cdsr, Cdbr, Csgr, Csdr, Cssr, Csbr;
double Cdgi, Cddi, Cdsi, Cdbi, Csgi, Csdi, Cssi, Csbi;
double gmr, gmi, gmbsr, gmbsi, gdsr, gdsi;
double FwdSumr, RevSumr, Gmr, Gmbsr;
double FwdSumi, RevSumi, Gmi, Gmbsi;
struct bsim4SizeDependParam *pParam;
double ggidld, ggidlg, ggidlb, ggislg, ggislb, ggisls;

double m;

    omega = ckt->CKTomega;
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

              Csd = -(here->BSIM4v7cddb + here->BSIM4v7cgdb + here->BSIM4v7cbdb);
              Csg = -(here->BSIM4v7cdgb + here->BSIM4v7cggb + here->BSIM4v7cbgb);
              Css = -(here->BSIM4v7cdsb + here->BSIM4v7cgsb + here->BSIM4v7cbsb);

              if (here->BSIM4v7acnqsMod)
              {   T0 = omega * here->BSIM4v7taunet;
                  T1 = T0 * T0;
                  T2 = 1.0 / (1.0 + T1);
                  T3 = T0 * T2;

                  gmr = here->BSIM4v7gm * T2;
                  gmbsr = here->BSIM4v7gmbs * T2;
                  gdsr = here->BSIM4v7gds * T2;

                  gmi = -here->BSIM4v7gm * T3;
                  gmbsi = -here->BSIM4v7gmbs * T3;
                  gdsi = -here->BSIM4v7gds * T3;

                  Cddr = here->BSIM4v7cddb * T2;
                  Cdgr = here->BSIM4v7cdgb * T2;
                  Cdsr = here->BSIM4v7cdsb * T2;
                  Cdbr = -(Cddr + Cdgr + Cdsr);

                  /* WDLiu: Cxyi mulitplied by jomega below, and actually to be of conductance */
                  Cddi = here->BSIM4v7cddb * T3 * omega;
                  Cdgi = here->BSIM4v7cdgb * T3 * omega;
                  Cdsi = here->BSIM4v7cdsb * T3 * omega;
                  Cdbi = -(Cddi + Cdgi + Cdsi);

                  Csdr = Csd * T2;
                  Csgr = Csg * T2;
                  Cssr = Css * T2;
                  Csbr = -(Csdr + Csgr + Cssr);

                  Csdi = Csd * T3 * omega;
                  Csgi = Csg * T3 * omega;
                  Cssi = Css * T3 * omega;
                  Csbi = -(Csdi + Csgi + Cssi);

                  Cgdr = -(Cddr + Csdr + here->BSIM4v7cbdb);
                  Cggr = -(Cdgr + Csgr + here->BSIM4v7cbgb);
                  Cgsr = -(Cdsr + Cssr + here->BSIM4v7cbsb);
                  Cgbr = -(Cgdr + Cggr + Cgsr);

                  Cgdi = -(Cddi + Csdi);
                  Cggi = -(Cdgi + Csgi);
                  Cgsi = -(Cdsi + Cssi);
                  Cgbi = -(Cgdi + Cggi + Cgsi);
              }
              else /* QS */
              {   gmr = here->BSIM4v7gm;
                  gmbsr = here->BSIM4v7gmbs;
                  gdsr = here->BSIM4v7gds;
                  gmi = gmbsi = gdsi = 0.0;

                  Cddr = here->BSIM4v7cddb;
                  Cdgr = here->BSIM4v7cdgb;
                  Cdsr = here->BSIM4v7cdsb;
                  Cdbr = -(Cddr + Cdgr + Cdsr);
                  Cddi = Cdgi = Cdsi = Cdbi = 0.0;

                  Csdr = Csd;
                  Csgr = Csg;
                  Cssr = Css;
                  Csbr = -(Csdr + Csgr + Cssr);
                  Csdi = Csgi = Cssi = Csbi = 0.0;

                  Cgdr = here->BSIM4v7cgdb;
                  Cggr = here->BSIM4v7cggb;
                  Cgsr = here->BSIM4v7cgsb;
                  Cgbr = -(Cgdr + Cggr + Cgsr);
                  Cgdi = Cggi = Cgsi = Cgbi = 0.0;
              }


              if (here->BSIM4v7mode >= 0) 
              {   Gmr = gmr;
                  Gmbsr = gmbsr;
                  FwdSumr = Gmr + Gmbsr;
                  RevSumr = 0.0;
                  Gmi = gmi;
                  Gmbsi = gmbsi;
                  FwdSumi = Gmi + Gmbsi;
                  RevSumi = 0.0;

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

                  if (here->BSIM4v7rgateMod == 3)
                  {   xcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * omega;
                      xcgmdb = -cgdo * omega;
                      xcgmsb = -cgso * omega;
                      xcgmbb = -pParam->BSIM4v7cgbo * omega;
    
                      xcdgmb = xcgmdb;
                      xcsgmb = xcgmsb;
                      xcbgmb = xcgmbb;
    
                      xcggbr = Cggr * omega;
                      xcgdbr = Cgdr * omega;
                      xcgsbr = Cgsr * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = Cdgr * omega;
                      xcsgbr = Csgr * omega;
                      xcbgb = here->BSIM4v7cbgb * omega;
                  }
                  else
                  {   xcggbr = (Cggr + cgdo + cgso + pParam->BSIM4v7cgbo ) * omega;
                      xcgdbr = (Cgdr - cgdo) * omega;
                      xcgsbr = (Cgsr - cgso) * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = (Cdgr - cgdo) * omega;
                      xcsgbr = (Csgr - cgso) * omega;
                      xcbgb = (here->BSIM4v7cbgb - pParam->BSIM4v7cgbo) * omega;
    
                      xcdgmb = xcsgmb = xcbgmb = 0.0;
                  }
                  xcddbr = (Cddr + here->BSIM4v7capbd + cgdo) * omega;
                  xcdsbr = Cdsr * omega;
                  xcsdbr = Csdr * omega;
                  xcssbr = (here->BSIM4v7capbs + cgso + Cssr) * omega;
    
                  if (!here->BSIM4v7rbodyMod)
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                      xcbdb = (here->BSIM4v7cbdb - here->BSIM4v7capbd) * omega;
                      xcbsb = (here->BSIM4v7cbsb - here->BSIM4v7capbs) * omega;
                      xcdbdb = 0.0;
                  }
                  else
                  {   xcdbbr = Cdbr * omega;
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb)
                             + here->BSIM4v7capbs * omega;

                      xcbdb = here->BSIM4v7cbdb * omega;
                      xcbsb = here->BSIM4v7cbsb * omega;
    
                      xcdbdb = -here->BSIM4v7capbd * omega;
                      xcsbsb = -here->BSIM4v7capbs * omega;
                  }
                  xcbbb = -(xcbdb + xcbgb + xcbsb + xcbgmb);

                  xcdgbi = Cdgi;
                  xcsgbi = Csgi;
                  xcddbi = Cddi;
                  xcdsbi = Cdsi;
                  xcsdbi = Csdi;
                  xcssbi = Cssi;
                  xcdbbi = Cdbi;
                  xcsbbi = Csbi;
                  xcggbi = Cggi;
                  xcgdbi = Cgdi;
                  xcgsbi = Cgsi;
                  xcgbbi = Cgbi;
              } 
              else /* Reverse mode */
              {   Gmr = -gmr;
                  Gmbsr = -gmbsr;
                  FwdSumr = 0.0;
                  RevSumr = -(Gmr + Gmbsr);
                  Gmi = -gmi;
                  Gmbsi = -gmbsi;
                  FwdSumi = 0.0;
                  RevSumi = -(Gmi + Gmbsi);

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

                  if (here->BSIM4v7rgateMod == 3)
                  {   xcgmgmb = (cgdo + cgso + pParam->BSIM4v7cgbo) * omega;
                      xcgmdb = -cgdo * omega;
                      xcgmsb = -cgso * omega;
                      xcgmbb = -pParam->BSIM4v7cgbo * omega;
    
                      xcdgmb = xcgmdb;
                      xcsgmb = xcgmsb;
                      xcbgmb = xcgmbb;
   
                      xcggbr = Cggr * omega;
                      xcgdbr = Cgsr * omega;
                      xcgsbr = Cgdr * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = Csgr * omega;
                      xcsgbr = Cdgr * omega;
                      xcbgb = here->BSIM4v7cbgb * omega;
                  }
                  else
                  {   xcggbr = (Cggr + cgdo + cgso + pParam->BSIM4v7cgbo ) * omega;
                      xcgdbr = (Cgsr - cgdo) * omega;
                      xcgsbr = (Cgdr - cgso) * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = (Csgr - cgdo) * omega;
                      xcsgbr = (Cdgr - cgso) * omega;
                      xcbgb = (here->BSIM4v7cbgb - pParam->BSIM4v7cgbo) * omega;
    
                      xcdgmb = xcsgmb = xcbgmb = 0.0;
                  }
                  xcddbr = (here->BSIM4v7capbd + cgdo + Cssr) * omega;
                  xcdsbr = Csdr * omega;
                  xcsdbr = Cdsr * omega;
                  xcssbr = (Cddr + here->BSIM4v7capbs + cgso) * omega;
    
                  if (!here->BSIM4v7rbodyMod)
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                      xcbdb = (here->BSIM4v7cbsb - here->BSIM4v7capbd) * omega;
                      xcbsb = (here->BSIM4v7cbdb - here->BSIM4v7capbs) * omega;
                      xcdbdb = 0.0;
                  }
                  else
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb)
                             + here->BSIM4v7capbd * omega;
                      xcsbbr = Cdbr * omega;

                      xcbdb = here->BSIM4v7cbsb * omega;
                      xcbsb = here->BSIM4v7cbdb * omega;
                      xcdbdb = -here->BSIM4v7capbd * omega;
                      xcsbsb = -here->BSIM4v7capbs * omega;
                  }
                  xcbbb = -(xcbgb + xcbdb + xcbsb + xcbgmb);

                  xcdgbi = Csgi;
                  xcsgbi = Cdgi;
                  xcddbi = Cssi;
                  xcdsbi = Csdi;
                  xcsdbi = Cdsi;
                  xcssbi = Cddi;
                  xcdbbi = Csbi;
                  xcsbbi = Cdbi;
                  xcggbi = Cggi;
                  xcgdbi = Cgsi;
                  xcgsbi = Cgdi;
                  xcgbbi = Cgbi;
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


              /*
               * Loading AC matrix
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

                  *(here->BSIM4v7GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4v7GPgpPtr) += m * (geltd + xcggbi + gIgtotg);
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4v7GPdpPtr) += m * (xcgdbi + gIgtotd);
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4v7GPspPtr) += m * (xcgsbi + gIgtots);
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4v7GPbpPtr) += m * (xcgbbi + gIgtotb);
              } /* WDLiu: gcrg already subtracted from all gcrgg below */
              else if (here->BSIM4v7rgateMod == 2)
              {   *(here->BSIM4v7GEgePtr) += m * gcrg;
                  *(here->BSIM4v7GEgpPtr) += m * gcrgg;
                  *(here->BSIM4v7GEdpPtr) += m * gcrgd;
                  *(here->BSIM4v7GEspPtr) += m * gcrgs;
                  *(here->BSIM4v7GEbpPtr) += m * gcrgb;

                  *(here->BSIM4v7GPgePtr) -= m * gcrg;
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4v7GPgpPtr) -= m * (gcrgg - xcggbi - gIgtotg);
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4v7GPdpPtr) -= m * (gcrgd - xcgdbi - gIgtotd);
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4v7GPspPtr) -= m * (gcrgs - xcgsbi - gIgtots);
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4v7GPbpPtr) -= m * (gcrgb - xcgbbi - gIgtotb);
              }
              else if (here->BSIM4v7rgateMod == 3)
              {   *(here->BSIM4v7GEgePtr) += m * geltd;
                  *(here->BSIM4v7GEgmPtr) -= m * geltd;
                  *(here->BSIM4v7GMgePtr) -= m * geltd;
                  *(here->BSIM4v7GMgmPtr) += m * (geltd + gcrg);
                  *(here->BSIM4v7GMgmPtr +1) += m * xcgmgmb;
   
                  *(here->BSIM4v7GMdpPtr) += m * gcrgd;
                  *(here->BSIM4v7GMdpPtr +1) += m * xcgmdb;
                  *(here->BSIM4v7GMgpPtr) += m * gcrgg;
                  *(here->BSIM4v7GMspPtr) += m * gcrgs;
                  *(here->BSIM4v7GMspPtr +1) += m * xcgmsb;
                  *(here->BSIM4v7GMbpPtr) += m * gcrgb;
                  *(here->BSIM4v7GMbpPtr +1) += m * xcgmbb;
   
                  *(here->BSIM4v7DPgmPtr +1) += m * xcdgmb;
                  *(here->BSIM4v7GPgmPtr) -= m * gcrg;
                  *(here->BSIM4v7SPgmPtr +1) += m * xcsgmb;
                  *(here->BSIM4v7BPgmPtr +1) += m * xcbgmb;
   
                  *(here->BSIM4v7GPgpPtr) -= m * (gcrgg - xcggbi - gIgtotg);
                  *(here->BSIM4v7GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4v7GPdpPtr) -= m * (gcrgd - xcgdbi - gIgtotd);
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4v7GPspPtr) -= m * (gcrgs - xcgsbi - gIgtots);
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4v7GPbpPtr) -= m * (gcrgb - xcgbbi - gIgtotb);
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbbr;
              }
              else
              {   *(here->BSIM4v7GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4v7GPgpPtr) += m * (xcggbi + gIgtotg);
                  *(here->BSIM4v7GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4v7GPdpPtr) += m * (xcgdbi + gIgtotd);
                  *(here->BSIM4v7GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4v7GPspPtr) += m * (xcgsbi + gIgtots);
                  *(here->BSIM4v7GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4v7GPbpPtr) += m * (xcgbbi + gIgtotb);
              }

              if (model->BSIM4v7rdsMod)
              {   (*(here->BSIM4v7DgpPtr) += m * gdtotg);
                  (*(here->BSIM4v7DspPtr) += m * gdtots);
                  (*(here->BSIM4v7DbpPtr) += m * gdtotb);
                  (*(here->BSIM4v7SdpPtr) += m * gstotd);
                  (*(here->BSIM4v7SgpPtr) += m * gstotg);
                  (*(here->BSIM4v7SbpPtr) += m * gstotb);
              }

              *(here->BSIM4v7DPdpPtr +1) += m * (xcddbr + gdsi + RevSumi);
              *(here->BSIM4v7DPdpPtr) += m * (gdpr + xcddbi + gdsr + here->BSIM4v7gbd 
                                     - gdtotd + RevSumr + gbdpdp - gIdtotd);
              *(here->BSIM4v7DPdPtr) -= m * (gdpr + gdtot);
              *(here->BSIM4v7DPgpPtr +1) += m * (xcdgbr + Gmi);
              *(here->BSIM4v7DPgpPtr) += m * (Gmr + xcdgbi - gdtotg + gbdpg - gIdtotg);
              *(here->BSIM4v7DPspPtr +1) += m * (xcdsbr - gdsi - FwdSumi);
              *(here->BSIM4v7DPspPtr) -= m * (gdsr - xcdsbi + FwdSumr + gdtots - gbdpsp + gIdtots);
              *(here->BSIM4v7DPbpPtr +1) += m * (xcdbbr + Gmbsi);
              *(here->BSIM4v7DPbpPtr) -= m * (gjbd + gdtotb - xcdbbi - Gmbsr - gbdpb + gIdtotb);

              *(here->BSIM4v7DdpPtr) -= m * (gdpr - gdtotd);
              *(here->BSIM4v7DdPtr) += m * (gdpr + gdtot);

              *(here->BSIM4v7SPdpPtr +1) += m * (xcsdbr - gdsi - RevSumi);
              *(here->BSIM4v7SPdpPtr) -= m * (gdsr - xcsdbi + gstotd + RevSumr - gbspdp + gIstotd);
              *(here->BSIM4v7SPgpPtr +1) += m * (xcsgbr - Gmi);
              *(here->BSIM4v7SPgpPtr) -= m * (Gmr - xcsgbi + gstotg - gbspg + gIstotg);
              *(here->BSIM4v7SPspPtr +1) += m * (xcssbr + gdsi + FwdSumi);
              *(here->BSIM4v7SPspPtr) += m * (gspr + xcssbi + gdsr + here->BSIM4v7gbs
                                     - gstots + FwdSumr + gbspsp - gIstots);
              *(here->BSIM4v7SPsPtr) -= m * (gspr + gstot);
              *(here->BSIM4v7SPbpPtr +1) += m * (xcsbbr - Gmbsi);
              *(here->BSIM4v7SPbpPtr) -= m * (gjbs + gstotb - xcsbbi + Gmbsr - gbspb + gIstotb);

              *(here->BSIM4v7SspPtr) -= m * (gspr - gstots);
              *(here->BSIM4v7SsPtr) += m * (gspr + gstot);

              *(here->BSIM4v7BPdpPtr +1) += m * xcbdb;
              *(here->BSIM4v7BPdpPtr) -= m * (gjbd - gbbdp + gIbtotd);
              *(here->BSIM4v7BPgpPtr +1) += m * xcbgb;
              *(here->BSIM4v7BPgpPtr) -= m * (here->BSIM4v7gbgs + gIbtotg);
              *(here->BSIM4v7BPspPtr +1) += m * xcbsb;
              *(here->BSIM4v7BPspPtr) -= m * (gjbs - gbbsp + gIbtots);
              *(here->BSIM4v7BPbpPtr +1) += m * xcbbb;
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
              {   (*(here->BSIM4v7DPdbPtr +1) += m * xcdbdb);
                  (*(here->BSIM4v7DPdbPtr) -= m * here->BSIM4v7gbd);
                  (*(here->BSIM4v7SPsbPtr +1) += m * xcsbsb);
                  (*(here->BSIM4v7SPsbPtr) -= m * here->BSIM4v7gbs);

                  (*(here->BSIM4v7DBdpPtr +1) += m * xcdbdb);
                  (*(here->BSIM4v7DBdpPtr) -= m * here->BSIM4v7gbd);
                  (*(here->BSIM4v7DBdbPtr +1) -= m * xcdbdb);
                  (*(here->BSIM4v7DBdbPtr) += m * (here->BSIM4v7gbd + here->BSIM4v7grbpd 
                                          + here->BSIM4v7grbdb));
                  (*(here->BSIM4v7DBbpPtr) -= m * here->BSIM4v7grbpd);
                  (*(here->BSIM4v7DBbPtr) -= m * here->BSIM4v7grbdb);

                  (*(here->BSIM4v7BPdbPtr) -= m * here->BSIM4v7grbpd);
                  (*(here->BSIM4v7BPbPtr) -= m * here->BSIM4v7grbpb);
                  (*(here->BSIM4v7BPsbPtr) -= m * here->BSIM4v7grbps);
                  (*(here->BSIM4v7BPbpPtr) += m * (here->BSIM4v7grbpd + here->BSIM4v7grbps 
                                          + here->BSIM4v7grbpb));
                  /* WDLiu: (-here->BSIM4v7gbbs) already added to BPbpPtr */

                  (*(here->BSIM4v7SBspPtr +1) += m * xcsbsb);
                  (*(here->BSIM4v7SBspPtr) -= m * here->BSIM4v7gbs);
                  (*(here->BSIM4v7SBbpPtr) -= m * here->BSIM4v7grbps);
                  (*(here->BSIM4v7SBbPtr) -= m * here->BSIM4v7grbsb);
                  (*(here->BSIM4v7SBsbPtr +1) -= m * xcsbsb);
                  (*(here->BSIM4v7SBsbPtr) += m * (here->BSIM4v7gbs
                                          + here->BSIM4v7grbps + here->BSIM4v7grbsb));

                  (*(here->BSIM4v7BdbPtr) -= m * here->BSIM4v7grbdb);
                  (*(here->BSIM4v7BbpPtr) -= m * here->BSIM4v7grbpb);
                  (*(here->BSIM4v7BsbPtr) -= m * here->BSIM4v7grbsb);
                  (*(here->BSIM4v7BbPtr) += m * (here->BSIM4v7grbsb + here->BSIM4v7grbdb
                                        + here->BSIM4v7grbpb));
              }


           /*
            * WDLiu: The internal charge node generated for transient NQS is not needed for
            *        AC NQS. The following is not doing a real job, but we have to keep it;
            *        otherwise a singular AC NQS matrix may occur if the transient NQS is on.
            *        The charge node is isolated from the instance.
            */
           if (here->BSIM4v7trnqsMod)
           {   (*(here->BSIM4v7QqPtr) += m * 1.0);
               (*(here->BSIM4v7QgpPtr) += 0.0);
               (*(here->BSIM4v7QdpPtr) += 0.0);
               (*(here->BSIM4v7QspPtr) += 0.0);
               (*(here->BSIM4v7QbpPtr) += 0.0);

               (*(here->BSIM4v7DPqPtr) += 0.0);
               (*(here->BSIM4v7SPqPtr) += 0.0);
               (*(here->BSIM4v7GPqPtr) += 0.0);
           }
         }
    }
    return(OK);
}
