/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4cvtest.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 05/09/2003.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h."
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
BSIM4v7convTest(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;
double delvbd, delvbs, delvds, delvgd, delvgs;
double delvdbd, delvsbs; 
double delvbd_jct, delvbs_jct;
double vds, vgs, vgd, vgdo, vbs, vbd;
double vdbd, vdbs, vsbs;
double cbhat, cdhat, Idtot, Ibtot;
double vses, vdes, vdedo, delvses, delvded, delvdes;
double Isestot, cseshat, Idedtot, cdedhat;
double Igstot, cgshat, Igdtot, cgdhat, Igbtot, cgbhat;
double tol0, tol1, tol2, tol3, tol4, tol5, tol6;

    for (; model != NULL; model = model->BSIM4v7nextModel)
    {    for (here = model->BSIM4v7instances; here != NULL ;
              here=here->BSIM4v7nextInstance) 
         {
              vds = model->BSIM4v7type
                  * (*(ckt->CKTrhsOld + here->BSIM4v7dNodePrime)
                  - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
              vgs = model->BSIM4v7type
                  * (*(ckt->CKTrhsOld + here->BSIM4v7gNodePrime)
                  - *(ckt->CKTrhsOld + here->BSIM4v7sNodePrime));
              vbs = model->BSIM4v7type
                  * (*(ckt->CKTrhsOld + here->BSIM4v7bNodePrime)
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
              
              vgdo = *(ckt->CKTstate0 + here->BSIM4v7vgs)
                    - *(ckt->CKTstate0 + here->BSIM4v7vds);
              vbd = vbs - vds;
              vdbd = vdbs - vds;
              vgd = vgs - vds;

              delvbd = vbd - *(ckt->CKTstate0 + here->BSIM4v7vbd);
              delvdbd = vdbd - *(ckt->CKTstate0 + here->BSIM4v7vdbd);
              delvgd = vgd - vgdo;

              delvds = vds - *(ckt->CKTstate0 + here->BSIM4v7vds);
              delvgs = vgs - *(ckt->CKTstate0 + here->BSIM4v7vgs);
              delvbs = vbs - *(ckt->CKTstate0 + here->BSIM4v7vbs);
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

              /*
               *  Check convergence
               */

              if ((here->BSIM4v7off == 0)  || (!(ckt->CKTmode & MODEINITFIX)))
              {   tol0 = ckt->CKTreltol * MAX(fabs(cdhat), fabs(Idtot))
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
                      return(OK);
                  } 

                  if ((fabs(cgshat - Igstot) >= tol3) || (fabs(cgdhat - Igdtot) >= tol4)
                      || (fabs(cgbhat - Igbtot) >= tol5))
                  {   ckt->CKTnoncon++;
                      return(OK);
                  }

                  Ibtot = here->BSIM4v7cbs + here->BSIM4v7cbd
                        - here->BSIM4v7Igidl - here->BSIM4v7Igisl - here->BSIM4v7csub;
                  if (here->BSIM4v7mode >= 0)
                  {   cbhat = Ibtot + here->BSIM4v7gbd * delvbd_jct
                            + here->BSIM4v7gbs * delvbs_jct - (here->BSIM4v7gbbs + here->BSIM4v7ggidlb)
                            * delvbs - (here->BSIM4v7gbgs + here->BSIM4v7ggidlg) * delvgs
                            - (here->BSIM4v7gbds + here->BSIM4v7ggidld) * delvds
                            - here->BSIM4v7ggislg * delvgd - here->BSIM4v7ggislb* delvbd + here->BSIM4v7ggisls * delvds ;
                  }
                  else
                  {   cbhat = Ibtot + here->BSIM4v7gbs * delvbs_jct + here->BSIM4v7gbd 
                         * delvbd_jct - (here->BSIM4v7gbbs + here->BSIM4v7ggislb) * delvbd
                         - (here->BSIM4v7gbgs + here->BSIM4v7ggislg) * delvgd
                         + (here->BSIM4v7gbds + here->BSIM4v7ggisld - here->BSIM4v7ggidls) * delvds
                         - here->BSIM4v7ggidlg * delvgs - here->BSIM4v7ggidlb * delvbs; 
                  }
                  tol6 = ckt->CKTreltol * MAX(fabs(cbhat), 
                        fabs(Ibtot)) + ckt->CKTabstol;
                  if (fabs(cbhat - Ibtot) > tol6) 
                  {   ckt->CKTnoncon++;
                      return(OK);
                  }
              }
         }
    }
    return(OK);
}
