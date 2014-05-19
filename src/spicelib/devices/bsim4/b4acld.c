/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4acld.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 10/05/2001.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
BSIM4acLoad(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4model *model = (BSIM4model*)inModel;
BSIM4instance *here;

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
    for (; model != NULL; model = model->BSIM4nextModel) 
    {    for (here = model->BSIM4instances; here!= NULL;
              here = here->BSIM4nextInstance) 
         {
              pParam = here->pParam;
              capbd = here->BSIM4capbd;
              capbs = here->BSIM4capbs;
              cgso = here->BSIM4cgso;
              cgdo = here->BSIM4cgdo;
              cgbo = pParam->BSIM4cgbo;

              Csd = -(here->BSIM4cddb + here->BSIM4cgdb + here->BSIM4cbdb);
              Csg = -(here->BSIM4cdgb + here->BSIM4cggb + here->BSIM4cbgb);
              Css = -(here->BSIM4cdsb + here->BSIM4cgsb + here->BSIM4cbsb);

              if (here->BSIM4acnqsMod)
              {   T0 = omega * here->BSIM4taunet;
                  T1 = T0 * T0;
                  T2 = 1.0 / (1.0 + T1);
                  T3 = T0 * T2;

                  gmr = here->BSIM4gm * T2;
                  gmbsr = here->BSIM4gmbs * T2;
                  gdsr = here->BSIM4gds * T2;

                  gmi = -here->BSIM4gm * T3;
                  gmbsi = -here->BSIM4gmbs * T3;
                  gdsi = -here->BSIM4gds * T3;

                  Cddr = here->BSIM4cddb * T2;
                  Cdgr = here->BSIM4cdgb * T2;
                  Cdsr = here->BSIM4cdsb * T2;
                  Cdbr = -(Cddr + Cdgr + Cdsr);

                  /* WDLiu: Cxyi mulitplied by jomega below, and actually to be of conductance */
                  Cddi = here->BSIM4cddb * T3 * omega;
                  Cdgi = here->BSIM4cdgb * T3 * omega;
                  Cdsi = here->BSIM4cdsb * T3 * omega;
                  Cdbi = -(Cddi + Cdgi + Cdsi);

                  Csdr = Csd * T2;
                  Csgr = Csg * T2;
                  Cssr = Css * T2;
                  Csbr = -(Csdr + Csgr + Cssr);

                  Csdi = Csd * T3 * omega;
                  Csgi = Csg * T3 * omega;
                  Cssi = Css * T3 * omega;
                  Csbi = -(Csdi + Csgi + Cssi);

                  Cgdr = -(Cddr + Csdr + here->BSIM4cbdb);
                  Cggr = -(Cdgr + Csgr + here->BSIM4cbgb);
                  Cgsr = -(Cdsr + Cssr + here->BSIM4cbsb);
                  Cgbr = -(Cgdr + Cggr + Cgsr);

                  Cgdi = -(Cddi + Csdi);
                  Cggi = -(Cdgi + Csgi);
                  Cgsi = -(Cdsi + Cssi);
                  Cgbi = -(Cgdi + Cggi + Cgsi);
              }
              else /* QS */
              {   gmr = here->BSIM4gm;
                  gmbsr = here->BSIM4gmbs;
                  gdsr = here->BSIM4gds;
                  gmi = gmbsi = gdsi = 0.0;

                  Cddr = here->BSIM4cddb;
                  Cdgr = here->BSIM4cdgb;
                  Cdsr = here->BSIM4cdsb;
                  Cdbr = -(Cddr + Cdgr + Cdsr);
                  Cddi = Cdgi = Cdsi = Cdbi = 0.0;

                  Csdr = Csd;
                  Csgr = Csg;
                  Cssr = Css;
                  Csbr = -(Csdr + Csgr + Cssr);
                  Csdi = Csgi = Cssi = Csbi = 0.0;

                  Cgdr = here->BSIM4cgdb;
                  Cggr = here->BSIM4cggb;
                  Cgsr = here->BSIM4cgsb;
                  Cgbr = -(Cgdr + Cggr + Cgsr);
                  Cgdi = Cggi = Cgsi = Cgbi = 0.0;
              }


              if (here->BSIM4mode >= 0) 
              {   Gmr = gmr;
                  Gmbsr = gmbsr;
                  FwdSumr = Gmr + Gmbsr;
                  RevSumr = 0.0;
                  Gmi = gmi;
                  Gmbsi = gmbsi;
                  FwdSumi = Gmi + Gmbsi;
                  RevSumi = 0.0;

                  gbbdp = -(here->BSIM4gbds);
                  gbbsp = here->BSIM4gbds + here->BSIM4gbgs + here->BSIM4gbbs;
                  gbdpg = here->BSIM4gbgs;
                  gbdpdp = here->BSIM4gbds;
                  gbdpb = here->BSIM4gbbs;
                  gbdpsp = -(gbdpg + gbdpdp + gbdpb);

                  gbspdp = 0.0;
                  gbspg = 0.0;
                  gbspb = 0.0;
                  gbspsp = 0.0;

                  if (model->BSIM4igcMod)
                  {   gIstotg = here->BSIM4gIgsg + here->BSIM4gIgcsg;
                      gIstotd = here->BSIM4gIgcsd;
                      gIstots = here->BSIM4gIgss + here->BSIM4gIgcss;
                      gIstotb = here->BSIM4gIgcsb;

                      gIdtotg = here->BSIM4gIgdg + here->BSIM4gIgcdg;
                      gIdtotd = here->BSIM4gIgdd + here->BSIM4gIgcdd;
                      gIdtots = here->BSIM4gIgcds;
                      gIdtotb = here->BSIM4gIgcdb;
                  }
                  else
                  {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                      gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
                  }

                  if (model->BSIM4igbMod)
                  {   gIbtotg = here->BSIM4gIgbg;
                      gIbtotd = here->BSIM4gIgbd;
                      gIbtots = here->BSIM4gIgbs;
                      gIbtotb = here->BSIM4gIgbb;
                  }
                  else
                      gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                  if ((model->BSIM4igcMod != 0) || (model->BSIM4igbMod != 0))
                  {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                      gIgtots = gIstots + gIdtots + gIbtots;
                      gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  }
                  else
                      gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                  if (here->BSIM4rgateMod == 2)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4vges)
                         - *(ckt->CKTstates[0] + here->BSIM4vgs);
                  else if (here->BSIM4rgateMod == 3)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4vgms)
                         - *(ckt->CKTstates[0] + here->BSIM4vgs);
                  if (here->BSIM4rgateMod > 1)
                  {   gcrgd = here->BSIM4gcrgd * T0;
                      gcrgg = here->BSIM4gcrgg * T0;
                      gcrgs = here->BSIM4gcrgs * T0;
                      gcrgb = here->BSIM4gcrgb * T0;
                      gcrgg -= here->BSIM4gcrg;
                      gcrg = here->BSIM4gcrg;
                  }
                  else
                      gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                  if (here->BSIM4rgateMod == 3)
                  {   xcgmgmb = (cgdo + cgso + pParam->BSIM4cgbo) * omega;
                      xcgmdb = -cgdo * omega;
                      xcgmsb = -cgso * omega;
                      xcgmbb = -pParam->BSIM4cgbo * omega;
    
                      xcdgmb = xcgmdb;
                      xcsgmb = xcgmsb;
                      xcbgmb = xcgmbb;
    
                      xcggbr = Cggr * omega;
                      xcgdbr = Cgdr * omega;
                      xcgsbr = Cgsr * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = Cdgr * omega;
                      xcsgbr = Csgr * omega;
                      xcbgb = here->BSIM4cbgb * omega;
                  }
                  else
                  {   xcggbr = (Cggr + cgdo + cgso + pParam->BSIM4cgbo ) * omega;
                      xcgdbr = (Cgdr - cgdo) * omega;
                      xcgsbr = (Cgsr - cgso) * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = (Cdgr - cgdo) * omega;
                      xcsgbr = (Csgr - cgso) * omega;
                      xcbgb = (here->BSIM4cbgb - pParam->BSIM4cgbo) * omega;
    
                      xcdgmb = xcsgmb = xcbgmb = 0.0;
                  }
                  xcddbr = (Cddr + here->BSIM4capbd + cgdo) * omega;
                  xcdsbr = Cdsr * omega;
                  xcsdbr = Csdr * omega;
                  xcssbr = (here->BSIM4capbs + cgso + Cssr) * omega;
    
                  if (!here->BSIM4rbodyMod)
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                      xcbdb = (here->BSIM4cbdb - here->BSIM4capbd) * omega;
                      xcbsb = (here->BSIM4cbsb - here->BSIM4capbs) * omega;
                      xcdbdb = 0.0;
                  }
                  else
                  {   xcdbbr = Cdbr * omega;
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb)
                             + here->BSIM4capbs * omega;

                      xcbdb = here->BSIM4cbdb * omega;
                      xcbsb = here->BSIM4cbsb * omega;
    
                      xcdbdb = -here->BSIM4capbd * omega;
                      xcsbsb = -here->BSIM4capbs * omega;
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

                  gbbsp = -(here->BSIM4gbds);
                  gbbdp = here->BSIM4gbds + here->BSIM4gbgs + here->BSIM4gbbs;

                  gbdpg = 0.0;
                  gbdpsp = 0.0;
                  gbdpb = 0.0;
                  gbdpdp = 0.0;

                  gbspg = here->BSIM4gbgs;
                  gbspsp = here->BSIM4gbds;
                  gbspb = here->BSIM4gbbs;
                  gbspdp = -(gbspg + gbspsp + gbspb);

                  if (model->BSIM4igcMod)
                  {   gIstotg = here->BSIM4gIgsg + here->BSIM4gIgcdg;
                      gIstotd = here->BSIM4gIgcds;
                      gIstots = here->BSIM4gIgss + here->BSIM4gIgcdd;
                      gIstotb = here->BSIM4gIgcdb;

                      gIdtotg = here->BSIM4gIgdg + here->BSIM4gIgcsg;
                      gIdtotd = here->BSIM4gIgdd + here->BSIM4gIgcss;
                      gIdtots = here->BSIM4gIgcsd;
                      gIdtotb = here->BSIM4gIgcsb;
                  }
                  else
                  {   gIstotg = gIstotd = gIstots = gIstotb = 0.0;
                      gIdtotg = gIdtotd = gIdtots = gIdtotb  = 0.0;
                  }

                  if (model->BSIM4igbMod)
                  {   gIbtotg = here->BSIM4gIgbg;
                      gIbtotd = here->BSIM4gIgbs;
                      gIbtots = here->BSIM4gIgbd;
                      gIbtotb = here->BSIM4gIgbb;
                  }
                  else
                      gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;

                  if ((model->BSIM4igcMod != 0) || (model->BSIM4igbMod != 0))
                  {   gIgtotg = gIstotg + gIdtotg + gIbtotg;
                      gIgtotd = gIstotd + gIdtotd + gIbtotd ;
                      gIgtots = gIstots + gIdtots + gIbtots;
                      gIgtotb = gIstotb + gIdtotb + gIbtotb;
                  }
                  else
                      gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

                  if (here->BSIM4rgateMod == 2)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4vges)
                         - *(ckt->CKTstates[0] + here->BSIM4vgs);
                  else if (here->BSIM4rgateMod == 3)
                      T0 = *(ckt->CKTstates[0] + here->BSIM4vgms)
                         - *(ckt->CKTstates[0] + here->BSIM4vgs);
                  if (here->BSIM4rgateMod > 1)
                  {   gcrgd = here->BSIM4gcrgs * T0;
                      gcrgg = here->BSIM4gcrgg * T0;
                      gcrgs = here->BSIM4gcrgd * T0;
                      gcrgb = here->BSIM4gcrgb * T0;
                      gcrgg -= here->BSIM4gcrg;
                      gcrg = here->BSIM4gcrg;
                  }
                  else
                      gcrg = gcrgd = gcrgg = gcrgs = gcrgb = 0.0;

                  if (here->BSIM4rgateMod == 3)
                  {   xcgmgmb = (cgdo + cgso + pParam->BSIM4cgbo) * omega;
                      xcgmdb = -cgdo * omega;
                      xcgmsb = -cgso * omega;
                      xcgmbb = -pParam->BSIM4cgbo * omega;
    
                      xcdgmb = xcgmdb;
                      xcsgmb = xcgmsb;
                      xcbgmb = xcgmbb;
   
                      xcggbr = Cggr * omega;
                      xcgdbr = Cgsr * omega;
                      xcgsbr = Cgdr * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = Csgr * omega;
                      xcsgbr = Cdgr * omega;
                      xcbgb = here->BSIM4cbgb * omega;
                  }
                  else
                  {   xcggbr = (Cggr + cgdo + cgso + pParam->BSIM4cgbo ) * omega;
                      xcgdbr = (Cgsr - cgdo) * omega;
                      xcgsbr = (Cgdr - cgso) * omega;
                      xcgbbr = -(xcggbr + xcgdbr + xcgsbr);
    
                      xcdgbr = (Csgr - cgdo) * omega;
                      xcsgbr = (Cdgr - cgso) * omega;
                      xcbgb = (here->BSIM4cbgb - pParam->BSIM4cgbo) * omega;
    
                      xcdgmb = xcsgmb = xcbgmb = 0.0;
                  }
                  xcddbr = (here->BSIM4capbd + cgdo + Cssr) * omega;
                  xcdsbr = Csdr * omega;
                  xcsdbr = Cdsr * omega;
                  xcssbr = (Cddr + here->BSIM4capbs + cgso) * omega;
    
                  if (!here->BSIM4rbodyMod)
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb);
                      xcsbbr = -(xcsgbr + xcsdbr + xcssbr + xcsgmb);

                      xcbdb = (here->BSIM4cbsb - here->BSIM4capbd) * omega;
                      xcbsb = (here->BSIM4cbdb - here->BSIM4capbs) * omega;
                      xcdbdb = 0.0;
                  }
                  else
                  {   xcdbbr = -(xcdgbr + xcddbr + xcdsbr + xcdgmb)
                             + here->BSIM4capbd * omega;
                      xcsbbr = Cdbr * omega;

                      xcbdb = here->BSIM4cbsb * omega;
                      xcbsb = here->BSIM4cbdb * omega;
                      xcdbdb = -here->BSIM4capbd * omega;
                      xcsbsb = -here->BSIM4capbs * omega;
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

              if (model->BSIM4rdsMod == 1)
              {   gstot = here->BSIM4gstot;
                  gstotd = here->BSIM4gstotd;
                  gstotg = here->BSIM4gstotg;
                  gstots = here->BSIM4gstots - gstot;
                  gstotb = here->BSIM4gstotb;

                  gdtot = here->BSIM4gdtot;
                  gdtotd = here->BSIM4gdtotd - gdtot;
                  gdtotg = here->BSIM4gdtotg;
                  gdtots = here->BSIM4gdtots;
                  gdtotb = here->BSIM4gdtotb;
              }
              else
              {   gstot = gstotd = gstotg = gstots = gstotb = 0.0;
                  gdtot = gdtotd = gdtotg = gdtots = gdtotb = 0.0;
              }


              /*
               * Loading AC matrix
               */
              m = here->BSIM4m;

              if (!model->BSIM4rdsMod)
              {   gdpr = here->BSIM4drainConductance;
                  gspr = here->BSIM4sourceConductance;
              }
              else
                  gdpr = gspr = 0.0;

              if (!here->BSIM4rbodyMod)
              {   gjbd = here->BSIM4gbd;
                  gjbs = here->BSIM4gbs;
              }
              else
                  gjbd = gjbs = 0.0;

              geltd = here->BSIM4grgeltd;

              if (here->BSIM4rgateMod == 1)
              {   *(here->BSIM4GEgePtr) += m * geltd;
                  *(here->BSIM4GPgePtr) -= m * geltd;
                  *(here->BSIM4GEgpPtr) -= m * geltd;

                  *(here->BSIM4GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4GPgpPtr) += m * (geltd + xcggbi + gIgtotg);
                  *(here->BSIM4GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4GPdpPtr) += m * (xcgdbi + gIgtotd);
                  *(here->BSIM4GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4GPspPtr) += m * (xcgsbi + gIgtots);
                  *(here->BSIM4GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4GPbpPtr) += m * (xcgbbi + gIgtotb);
              } /* WDLiu: gcrg already subtracted from all gcrgg below */
              else if (here->BSIM4rgateMod == 2)
              {   *(here->BSIM4GEgePtr) += m * gcrg;
                  *(here->BSIM4GEgpPtr) += m * gcrgg;
                  *(here->BSIM4GEdpPtr) += m * gcrgd;
                  *(here->BSIM4GEspPtr) += m * gcrgs;
                  *(here->BSIM4GEbpPtr) += m * gcrgb;

                  *(here->BSIM4GPgePtr) -= m * gcrg;
                  *(here->BSIM4GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4GPgpPtr) -= m * (gcrgg - xcggbi - gIgtotg);
                  *(here->BSIM4GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4GPdpPtr) -= m * (gcrgd - xcgdbi - gIgtotd);
                  *(here->BSIM4GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4GPspPtr) -= m * (gcrgs - xcgsbi - gIgtots);
                  *(here->BSIM4GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4GPbpPtr) -= m * (gcrgb - xcgbbi - gIgtotb);
              }
              else if (here->BSIM4rgateMod == 3)
              {   *(here->BSIM4GEgePtr) += m * geltd;
                  *(here->BSIM4GEgmPtr) -= m * geltd;
                  *(here->BSIM4GMgePtr) -= m * geltd;
                  *(here->BSIM4GMgmPtr) += m * (geltd + gcrg);
                  *(here->BSIM4GMgmPtr +1) += m * xcgmgmb;
   
                  *(here->BSIM4GMdpPtr) += m * gcrgd;
                  *(here->BSIM4GMdpPtr +1) += m * xcgmdb;
                  *(here->BSIM4GMgpPtr) += m * gcrgg;
                  *(here->BSIM4GMspPtr) += m * gcrgs;
                  *(here->BSIM4GMspPtr +1) += m * xcgmsb;
                  *(here->BSIM4GMbpPtr) += m * gcrgb;
                  *(here->BSIM4GMbpPtr +1) += m * xcgmbb;
   
                  *(here->BSIM4DPgmPtr +1) += m * xcdgmb;
                  *(here->BSIM4GPgmPtr) -= m * gcrg;
                  *(here->BSIM4SPgmPtr +1) += m * xcsgmb;
                  *(here->BSIM4BPgmPtr +1) += m * xcbgmb;
   
                  *(here->BSIM4GPgpPtr) -= m * (gcrgg - xcggbi - gIgtotg);
                  *(here->BSIM4GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4GPdpPtr) -= m * (gcrgd - xcgdbi - gIgtotd);
                  *(here->BSIM4GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4GPspPtr) -= m * (gcrgs - xcgsbi - gIgtots);
                  *(here->BSIM4GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4GPbpPtr) -= m * (gcrgb - xcgbbi - gIgtotb);
                  *(here->BSIM4GPbpPtr +1) += m * xcgbbr;
              }
              else
              {   *(here->BSIM4GPgpPtr +1) += m * xcggbr;
                  *(here->BSIM4GPgpPtr) += m * (xcggbi + gIgtotg);
                  *(here->BSIM4GPdpPtr +1) += m * xcgdbr;
                  *(here->BSIM4GPdpPtr) += m * (xcgdbi + gIgtotd);
                  *(here->BSIM4GPspPtr +1) += m * xcgsbr;
                  *(here->BSIM4GPspPtr) += m * (xcgsbi + gIgtots);
                  *(here->BSIM4GPbpPtr +1) += m * xcgbbr;
                  *(here->BSIM4GPbpPtr) += m * (xcgbbi + gIgtotb);
              }

              if (model->BSIM4rdsMod)
              {   (*(here->BSIM4DgpPtr) += m * gdtotg);
                  (*(here->BSIM4DspPtr) += m * gdtots);
                  (*(here->BSIM4DbpPtr) += m * gdtotb);
                  (*(here->BSIM4SdpPtr) += m * gstotd);
                  (*(here->BSIM4SgpPtr) += m * gstotg);
                  (*(here->BSIM4SbpPtr) += m * gstotb);
              }

              *(here->BSIM4DPdpPtr +1) += m * (xcddbr + gdsi + RevSumi);
              *(here->BSIM4DPdpPtr) += m * (gdpr + xcddbi + gdsr + here->BSIM4gbd 
                                     - gdtotd + RevSumr + gbdpdp - gIdtotd);
              *(here->BSIM4DPdPtr) -= m * (gdpr + gdtot);
              *(here->BSIM4DPgpPtr +1) += m * (xcdgbr + Gmi);
              *(here->BSIM4DPgpPtr) += m * (Gmr + xcdgbi - gdtotg + gbdpg - gIdtotg);
              *(here->BSIM4DPspPtr +1) += m * (xcdsbr - gdsi - FwdSumi);
              *(here->BSIM4DPspPtr) -= m * (gdsr - xcdsbi + FwdSumr + gdtots - gbdpsp + gIdtots);
              *(here->BSIM4DPbpPtr +1) += m * (xcdbbr + Gmbsi);
              *(here->BSIM4DPbpPtr) -= m * (gjbd + gdtotb - xcdbbi - Gmbsr - gbdpb + gIdtotb);

              *(here->BSIM4DdpPtr) -= m * (gdpr - gdtotd);
              *(here->BSIM4DdPtr) += m * (gdpr + gdtot);

              *(here->BSIM4SPdpPtr +1) += m * (xcsdbr - gdsi - RevSumi);
              *(here->BSIM4SPdpPtr) -= m * (gdsr - xcsdbi + gstotd + RevSumr - gbspdp + gIstotd);
              *(here->BSIM4SPgpPtr +1) += m * (xcsgbr - Gmi);
              *(here->BSIM4SPgpPtr) -= m * (Gmr - xcsgbi + gstotg - gbspg + gIstotg);
              *(here->BSIM4SPspPtr +1) += m * (xcssbr + gdsi + FwdSumi);
              *(here->BSIM4SPspPtr) += m * (gspr + xcssbi + gdsr + here->BSIM4gbs
                                     - gstots + FwdSumr + gbspsp - gIstots);
              *(here->BSIM4SPsPtr) -= m * (gspr + gstot);
              *(here->BSIM4SPbpPtr +1) += m * (xcsbbr - Gmbsi);
              *(here->BSIM4SPbpPtr) -= m * (gjbs + gstotb - xcsbbi + Gmbsr - gbspb + gIstotb);

              *(here->BSIM4SspPtr) -= m * (gspr - gstots);
              *(here->BSIM4SsPtr) += m * (gspr + gstot);

              *(here->BSIM4BPdpPtr +1) += m * xcbdb;
              *(here->BSIM4BPdpPtr) -= m * (gjbd - gbbdp + gIbtotd);
              *(here->BSIM4BPgpPtr +1) += m * xcbgb;
              *(here->BSIM4BPgpPtr) -= m * (here->BSIM4gbgs + gIbtotg);
              *(here->BSIM4BPspPtr +1) += m * xcbsb;
              *(here->BSIM4BPspPtr) -= m * (gjbs - gbbsp + gIbtots);
              *(here->BSIM4BPbpPtr +1) += m * xcbbb;
              *(here->BSIM4BPbpPtr) += m * (gjbd + gjbs - here->BSIM4gbbs
                                     - gIbtotb);
              ggidld = here->BSIM4ggidld;
              ggidlg = here->BSIM4ggidlg;
              ggidlb = here->BSIM4ggidlb;
              ggislg = here->BSIM4ggislg;
              ggisls = here->BSIM4ggisls;
              ggislb = here->BSIM4ggislb;

              /* stamp gidl */
              (*(here->BSIM4DPdpPtr) += m * ggidld);
              (*(here->BSIM4DPgpPtr) += m * ggidlg);
              (*(here->BSIM4DPspPtr) -= m * ((ggidlg + ggidld) + ggidlb));
              (*(here->BSIM4DPbpPtr) += m * ggidlb);
              (*(here->BSIM4BPdpPtr) -= m * ggidld);
              (*(here->BSIM4BPgpPtr) -= m * ggidlg);
              (*(here->BSIM4BPspPtr) += m * ((ggidlg + ggidld) + ggidlb));
              (*(here->BSIM4BPbpPtr) -= m * ggidlb);
               /* stamp gisl */
              (*(here->BSIM4SPdpPtr) -= m * ((ggisls + ggislg) + ggislb));
              (*(here->BSIM4SPgpPtr) += m * ggislg);
              (*(here->BSIM4SPspPtr) += m * ggisls);
              (*(here->BSIM4SPbpPtr) += m * ggislb);
              (*(here->BSIM4BPdpPtr) += m * ((ggislg + ggisls) + ggislb));
              (*(here->BSIM4BPgpPtr) -= m * ggislg);
              (*(here->BSIM4BPspPtr) -= m * ggisls);
              (*(here->BSIM4BPbpPtr) -= m * ggislb);

              if (here->BSIM4rbodyMod)
              {   (*(here->BSIM4DPdbPtr +1) += m * xcdbdb);
                  (*(here->BSIM4DPdbPtr) -= m * here->BSIM4gbd);
                  (*(here->BSIM4SPsbPtr +1) += m * xcsbsb);
                  (*(here->BSIM4SPsbPtr) -= m * here->BSIM4gbs);

                  (*(here->BSIM4DBdpPtr +1) += m * xcdbdb);
                  (*(here->BSIM4DBdpPtr) -= m * here->BSIM4gbd);
                  (*(here->BSIM4DBdbPtr +1) -= m * xcdbdb);
                  (*(here->BSIM4DBdbPtr) += m * (here->BSIM4gbd + here->BSIM4grbpd 
                                          + here->BSIM4grbdb));
                  (*(here->BSIM4DBbpPtr) -= m * here->BSIM4grbpd);
                  (*(here->BSIM4DBbPtr) -= m * here->BSIM4grbdb);

                  (*(here->BSIM4BPdbPtr) -= m * here->BSIM4grbpd);
                  (*(here->BSIM4BPbPtr) -= m * here->BSIM4grbpb);
                  (*(here->BSIM4BPsbPtr) -= m * here->BSIM4grbps);
                  (*(here->BSIM4BPbpPtr) += m * (here->BSIM4grbpd + here->BSIM4grbps 
                                          + here->BSIM4grbpb));
                  /* WDLiu: (-here->BSIM4gbbs) already added to BPbpPtr */

                  (*(here->BSIM4SBspPtr +1) += m * xcsbsb);
                  (*(here->BSIM4SBspPtr) -= m * here->BSIM4gbs);
                  (*(here->BSIM4SBbpPtr) -= m * here->BSIM4grbps);
                  (*(here->BSIM4SBbPtr) -= m * here->BSIM4grbsb);
                  (*(here->BSIM4SBsbPtr +1) -= m * xcsbsb);
                  (*(here->BSIM4SBsbPtr) += m * (here->BSIM4gbs
                                          + here->BSIM4grbps + here->BSIM4grbsb));

                  (*(here->BSIM4BdbPtr) -= m * here->BSIM4grbdb);
                  (*(here->BSIM4BbpPtr) -= m * here->BSIM4grbpb);
                  (*(here->BSIM4BsbPtr) -= m * here->BSIM4grbsb);
                  (*(here->BSIM4BbPtr) += m * (here->BSIM4grbsb + here->BSIM4grbdb
                                        + here->BSIM4grbpb));
              }


           /*
            * WDLiu: The internal charge node generated for transient NQS is not needed for
            *        AC NQS. The following is not doing a real job, but we have to keep it;
            *        otherwise a singular AC NQS matrix may occur if the transient NQS is on.
            *        The charge node is isolated from the instance.
            */
           if (here->BSIM4trnqsMod)
           {   (*(here->BSIM4QqPtr) += m * 1.0);
               (*(here->BSIM4QgpPtr) += 0.0);
               (*(here->BSIM4QdpPtr) += 0.0);
               (*(here->BSIM4QspPtr) += 0.0);
               (*(here->BSIM4QbpPtr) += 0.0);

               (*(here->BSIM4DPqPtr) += 0.0);
               (*(here->BSIM4SPqPtr) += 0.0);
               (*(here->BSIM4GPqPtr) += 0.0);
           }
         }
    }
    return(OK);
}
