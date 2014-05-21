/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2cvtest.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "hsm2def.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int HSM2convTest(
     GENmodel *inModel,
     CKTcircuit *ckt)
{
  HSM2model *model = (HSM2model*)inModel;
  HSM2instance *here;
  double delvbd, delvbs, delvds, delvgd, delvgs, vbd, vbs, vds;
  double cd, cdhat, tol0, tol1, tol2, tol3, tol4, vgd, vgdo, vgs;
  double Ibtot, cbhat, Igstot, cgshat, Igdtot, cgdhat, Igbtot, cgbhat;

  /*  loop through all the HSM2 device models */
  for ( ; model != NULL; model = model->HSM2nextModel ) {
    /* loop through all the instances of the model */
    for ( here = model->HSM2instances; here != NULL ;
	  here = here->HSM2nextInstance ) {
      vbs = model->HSM2_type * 
	(*(ckt->CKTrhsOld+here->HSM2bNode) - 
	 *(ckt->CKTrhsOld+here->HSM2sNodePrime));
      vgs = model->HSM2_type *
	(*(ckt->CKTrhsOld+here->HSM2gNode) - 
	 *(ckt->CKTrhsOld+here->HSM2sNodePrime));
      vds = model->HSM2_type * 
	(*(ckt->CKTrhsOld+here->HSM2dNodePrime) - 
	 *(ckt->CKTrhsOld+here->HSM2sNodePrime));
      vbd = vbs - vds;
      vgd = vgs - vds;
      vgdo = *(ckt->CKTstate0 + here->HSM2vgs) - 
	*(ckt->CKTstate0 + here->HSM2vds);
      delvbs = vbs - *(ckt->CKTstate0 + here->HSM2vbs);
      delvbd = vbd - *(ckt->CKTstate0 + here->HSM2vbd);
      delvgs = vgs - *(ckt->CKTstate0 + here->HSM2vgs);
      delvds = vds - *(ckt->CKTstate0 + here->HSM2vds);
      delvgd = vgd - vgdo;

      cd = here->HSM2_ids - here->HSM2_ibd;
      if ( here->HSM2_mode >= 0 ) {
	cd += here->HSM2_isub + here->HSM2_igidl;
	cdhat = cd - here->HSM2_gbd * delvbd 
	  + (here->HSM2_gmbs + here->HSM2_gbbs + here->HSM2_gigidlbs) * delvbs
	  + (here->HSM2_gm + here->HSM2_gbgs + here->HSM2_gigidlgs) * delvgs
	  + (here->HSM2_gds + here->HSM2_gbds + here->HSM2_gigidlds) * delvds;
	Ibtot = here->HSM2_ibs + here->HSM2_ibd - here->HSM2_isub
	  - here->HSM2_igidl - here->HSM2_igisl;
	cbhat = Ibtot + here->HSM2_gbd * delvbd
	  + (here->HSM2_gbs -  here->HSM2_gbbs - here->HSM2_gigidlbs) * delvbs
	  - (here->HSM2_gbgs + here->HSM2_gigidlgs) * delvgs
	  - (here->HSM2_gbds + here->HSM2_gigidlds) * delvds
	  - here->HSM2_gigislgd * delvgd - here->HSM2_gigislbd * delvbd
	  + here->HSM2_gigislsd * delvds;
	Igstot = here->HSM2_igs;
	cgshat = Igstot + here->HSM2_gigsg * delvgs + 
	  here->HSM2_gigsd * delvds + here->HSM2_gigsb * delvbs;
	Igdtot = here->HSM2_igd;
	cgdhat = Igdtot + here->HSM2_gigdg * delvgs + 
	  here->HSM2_gigdd * delvds + here->HSM2_gigdb * delvbs;
	Igbtot = here->HSM2_igb;
	cgbhat = Igbtot + here->HSM2_gigbg * delvgs + 
	  here->HSM2_gigbd * delvds + here->HSM2_gigbb * delvbs;
      }
      else {
	cd -= here->HSM2_igidl;
	cdhat = cd 
	  + (here->HSM2_gmbs + here->HSM2_gbd - here->HSM2_gigidlbs) * delvbd 
	  + (here->HSM2_gm - here->HSM2_gigidlgs) * delvgd 
	  + (- here->HSM2_gds + here->HSM2_gigidlds) * delvds;
	Ibtot = here->HSM2_ibs + here->HSM2_ibd - here->HSM2_isub
	  - here->HSM2_igidl - here->HSM2_igisl;
	cbhat = Ibtot + here->HSM2_gbs * delvbs
	  + (here->HSM2_gbd - here->HSM2_gbbs - here->HSM2_gigidlbs) * delvbd
	  - (here->HSM2_gbgs + here->HSM2_gigidlgs) * delvgd
	  + (here->HSM2_gbds + here->HSM2_gigidlds) * delvds
	  - here->HSM2_gigislgd * delvgd - here->HSM2_gigislbd * delvbd
	  + here->HSM2_gigislsd * delvds;
	Igbtot = here->HSM2_igb;
	cgbhat = Igbtot + here->HSM2_gigbg * delvgd 
	  - here->HSM2_gigbs * delvds + here->HSM2_gigbb * delvbd;
	Igstot = here->HSM2_igs;
	cgshat = Igstot + here->HSM2_gigsg * delvgd
	  - here->HSM2_gigss * delvds + here->HSM2_gigsb * delvbd;
	Igdtot = here->HSM2_igd;
	cgdhat = Igdtot + here->HSM2_gigdg * delvgd 
	  - here->HSM2_gigds * delvds + here->HSM2_gigdb * delvbd;
      }

      /*
       *  check convergence
       */
      if ( here->HSM2_off == 0  || !(ckt->CKTmode & MODEINITFIX) ) {
	tol0 = ckt->CKTreltol * MAX(fabs(cdhat), fabs(cd)) + ckt->CKTabstol;
	tol1 = ckt->CKTreltol * MAX(fabs(cgshat), fabs(Igstot)) + ckt->CKTabstol;
	tol2 = ckt->CKTreltol * MAX(fabs(cgdhat), fabs(Igdtot)) + ckt->CKTabstol;
	tol3 = ckt->CKTreltol * MAX(fabs(cgbhat), fabs(Igbtot)) + ckt->CKTabstol;
	tol4 = ckt->CKTreltol * MAX(fabs(cbhat), fabs(Ibtot)) + ckt->CKTabstol;

	if ( (fabs(cdhat - cd) >= tol0)
	     || (fabs(cgshat - Igstot) >= tol1) 
	     || (fabs(cgdhat - Igdtot) >= tol2)
	     || (fabs(cgbhat - Igbtot) >= tol3) 
	     || (fabs(cbhat - Ibtot) >= tol4) ) {
	  ckt->CKTnoncon++;
	  return(OK);
	}
      }
    }
  }
  return(OK);
}
