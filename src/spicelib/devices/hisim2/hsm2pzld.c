/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2pzld.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/complex.h"
#include "ngspice/sperror.h"
#include "hsm2def.h"
#include "ngspice/suffix.h"

int HSM2pzLoad(
     GENmodel *inModel,
     CKTcircuit *ckt,
     SPcomplex *s)
{
  HSM2model *model = (HSM2model*)inModel;
  HSM2instance *here;
  double xcggb, xcgdb, xcgsb, xcgbb, xcbgb, xcbdb = 0.0, xcbsb, xcbbb;
  double xcdgb, xcddb, xcdsb, xcdbb, xcsgb, xcsdb = 0.0, xcssb, xcsbb;
  double xcdbdb = 0.0, xcsbsb = 0.0;
  double gdpr, gspr, gds, gbd, gbs, capbd, capbs, FwdSum, RevSum, gm, gmbs;
  double gjbd, gjbs, grg;
  double cggb, cgdb, cgsb, cbgb, cbdb, cbsb, cddb, cdgb, cdsb;
  double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
  double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
  double gIbtotg, gIbtotd, gIbtots, gIbtotb;
  double gIgtotg, gIgtotd, gIgtots, gIgtotb;
  double gIdtotg, gIdtotd, gIdtots, gIdtotb;
  double gIstotg, gIstotd, gIstots, gIstotb;
  
  NG_IGNORE(ckt);

  for ( ;model != NULL ;model = model->HSM2nextModel ) {
    for ( here = model->HSM2instances ;here!= NULL ;
	  here = here->HSM2nextInstance ) {
      if ( here->HSM2_mode >= 0 ) {
	gm = here->HSM2_gm;
	gmbs = here->HSM2_gmbs;
	FwdSum = gm + gmbs;
	RevSum = 0.0;
	
	gbbdp = -here->HSM2_gbds;
	gbbsp = here->HSM2_gbds + here->HSM2_gbgs + here->HSM2_gbbs;
	
	gbdpg = here->HSM2_gbgs;
	gbdpdp = here->HSM2_gbds;
	gbdpb = here->HSM2_gbbs;
	gbdpsp = -(gbdpg + gbdpdp + gbdpb);
	
	gbspg = 0.0;
	gbspdp = 0.0;
	gbspb = 0.0;
	gbspsp = 0.0;

	if (model->HSM2_coiigs) {
	  gIbtotg = here->HSM2_gigbg;
	  gIbtotd = here->HSM2_gigbd;
	  gIbtots = here->HSM2_gigbs;
	  gIbtotb = here->HSM2_gigbb;

	  gIstotg = here->HSM2_gigsg;
	  gIstotd = here->HSM2_gigsd;
	  gIstots = here->HSM2_gigss;
	  gIstotb = here->HSM2_gigsb;

	  gIdtotg = here->HSM2_gigdg;
	  gIdtotd = here->HSM2_gigdd;
	  gIdtots = here->HSM2_gigds;
	  gIdtotb = here->HSM2_gigdb;

	}
	else {
	  gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;
	  gIstotg = gIstotd = gIstots = gIstotb = 0.0;
	  gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
	}

	if (model->HSM2_coiigs) {
	  gIgtotg = gIbtotg + gIstotg + gIdtotg;
	  gIgtotd = gIbtotd + gIstotd + gIdtotd;
	  gIgtots = gIbtots + gIstots + gIdtots;
	  gIgtotb = gIbtotb + gIstotb + gIdtotb;
	}
	else
	  gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;

	cggb = here->HSM2_cggb;
	cgsb = here->HSM2_cgsb;
	cgdb = here->HSM2_cgdb;
	
	cbgb = here->HSM2_cbgb;
	cbsb = here->HSM2_cbsb;
	cbdb = here->HSM2_cbdb;
	
	cdgb = here->HSM2_cdgb;
	cdsb = here->HSM2_cdsb;
	cddb = here->HSM2_cddb;
	
      }
      else {
	gm = -here->HSM2_gm;
	gmbs = -here->HSM2_gmbs;
	FwdSum = 0.0;
	RevSum = -(gm + gmbs);
	
	gbbsp = -here->HSM2_gbds;
	gbbdp = here->HSM2_gbds + here->HSM2_gbgs + here->HSM2_gbbs;
	
	gbdpg = 0.0;
	gbdpsp = 0.0;
	gbdpb = 0.0;
	gbdpdp = 0.0;

	if (model->HSM2_coiigs) {
	  gIbtotg = here->HSM2_gigbg;
	  gIbtotd = here->HSM2_gigbd;
	  gIbtots = here->HSM2_gigbs;
	  gIbtotb = here->HSM2_gigbb;

	  gIstotg = here->HSM2_gigsg;
	  gIstotd = here->HSM2_gigsd;
	  gIstots = here->HSM2_gigss;
	  gIstotb = here->HSM2_gigsb;

	  gIdtotg = here->HSM2_gigdg;
	  gIdtotd = here->HSM2_gigdd;
	  gIdtots = here->HSM2_gigds;
	  gIdtotb = here->HSM2_gigdb;
	}
	else {
	  gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;
	  gIstotg = gIstotd = gIstots = gIstotb = 0.0;
	  gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
	}
	
	if (model->HSM2_coiigs) {
	  gIgtotg = gIbtotg + gIstotg + gIdtotg;
	  gIgtotd = gIbtotd + gIstotd + gIdtotd;
	  gIgtots = gIbtots + gIstots + gIdtots;
	  gIgtotb = gIbtotb + gIstotb + gIdtotb;
	}
	else
	  gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;
	
	gbspg = here->HSM2_gbgs;
	gbspsp = here->HSM2_gbds;
	gbspb = here->HSM2_gbbs;
	gbspdp = -(gbspg + gbspsp + gbspb);
	
	cggb = here->HSM2_cggb;
	cgsb = here->HSM2_cgdb;
	cgdb = here->HSM2_cgsb;
	  
	cbgb = here->HSM2_cbgb;
	cbsb = here->HSM2_cbdb;
	cbdb = here->HSM2_cbsb;
	
	cdgb = -(here->HSM2_cdgb + cggb + cbgb);
	cdsb = -(here->HSM2_cddb + cgsb + cbsb);
	cddb = -(here->HSM2_cdsb + cgdb + cbdb);
      }
      
      gdpr = here->HSM2drainConductance;
      gspr = here->HSM2sourceConductance;
      gds = here->HSM2_gds;
      gbd = here->HSM2_gbd;
      gbs = here->HSM2_gbs;
      capbd = here->HSM2_capbd;
      capbs = here->HSM2_capbs;
      
      xcdgb = cdgb;
      xcddb = cddb + capbd;
      xcdsb = cdsb;
      xcdbb = -(xcdgb + xcddb + xcdsb);
      if (here->HSM2_corbnet == 1) xcdbb += capbd;

      xcsgb = -(cggb + cbgb + cdgb);
      xcsdb = -(cgdb + cbdb + cddb);
      xcssb = capbs - (cgsb + cbsb + cdsb);
      xcsbb = -(xcsgb + xcsdb + xcssb); 
      if (here->HSM2_corbnet == 1) xcsbb += capbs;

      xcggb = cggb;
      xcgdb = cgdb;
      xcgsb = cgsb;
      xcgbb = -(xcggb + xcgdb + xcgsb);

      xcbgb = cbgb;
      if (!here->HSM2_corbnet) {
	xcbdb = cbdb - capbd;
	xcbsb = cbsb - capbs;
      } else {
	xcbdb = cbdb;
	xcbsb = cbsb;
	xcdbdb = - capbd;
	xcsbsb = - capbs;
      }
      xcbbb = -(xcbgb + xcbdb + xcbsb);

      if (!here->HSM2_corbnet) {
	gjbd = gbd;
	gjbs = gbs;
      } else
	gjbd = gjbs = 0.0;

      if (here->HSM2_corg == 1) {
	grg = here->HSM2_grg;
	*(here->HSM2GgPtr) += grg;
	*(here->HSM2GPgPtr) -= grg;
	*(here->HSM2GgpPtr) -= grg;

	*(here->HSM2GPdpPtr ) += xcgdb * s->real;
	*(here->HSM2GPdpPtr +1) += xcgdb * s->imag;
	*(here->HSM2GPdpPtr) += grg + gIgtotd;
	*(here->HSM2GPgpPtr ) += xcggb * s->real;
	*(here->HSM2GPgpPtr +1) += xcggb * s->imag;
	*(here->HSM2GPgpPtr) += gIgtotg;
	*(here->HSM2GPspPtr ) += xcgsb * s->real;
	*(here->HSM2GPspPtr +1) += xcgsb * s->imag;
	*(here->HSM2GPspPtr) += gIgtots;
	*(here->HSM2GPbpPtr ) += xcgbb * s->real;
	*(here->HSM2GPbpPtr +1) += xcgbb * s->imag;
	*(here->HSM2GPbpPtr) += gIgtotb;

      } else {
	*(here->HSM2GPdpPtr ) += xcgdb * s->real;
	*(here->HSM2GPdpPtr +1) += xcgdb * s->imag;
	*(here->HSM2GPdpPtr) += gIgtotd;
	*(here->HSM2GPgpPtr ) += xcggb * s->real;
	*(here->HSM2GPgpPtr +1) += xcggb * s->imag;
	*(here->HSM2GPgpPtr) += gIgtotg;
	*(here->HSM2GPspPtr ) += xcgsb * s->real;
	*(here->HSM2GPspPtr +1) += xcgsb * s->imag;
	*(here->HSM2GPspPtr) += gIgtots;
	*(here->HSM2GPbpPtr ) += xcgbb * s->real;
	*(here->HSM2GPbpPtr +1) += xcgbb * s->imag;
	*(here->HSM2GPbpPtr) += gIgtotb;
      }

      *(here->HSM2DPdpPtr ) += xcddb * s->real;
      *(here->HSM2DPdpPtr +1) += xcddb * s->imag;
      *(here->HSM2DPdpPtr) += gdpr + gds + gbd + RevSum + gbdpdp - gIdtotd;
      *(here->HSM2DPdPtr) -= gdpr;
      *(here->HSM2DPgpPtr ) += xcdgb * s->real;
      *(here->HSM2DPgpPtr +1) += xcdgb * s->imag;
      *(here->HSM2DPgpPtr) += gm + gbdpg - gIdtotg;
      *(here->HSM2DPspPtr ) += xcdsb * s->real;
      *(here->HSM2DPspPtr +1) += xcdsb * s->imag;
      *(here->HSM2DPspPtr) -= gds + FwdSum - gbdpsp + gIdtots;
      *(here->HSM2DPbpPtr ) += xcdbb * s->real;
      *(here->HSM2DPbpPtr +1) += xcdbb * s->imag;
      *(here->HSM2DPbpPtr) -= gjbd - gmbs - gbdpb + gIdtotb;

      *(here->HSM2DdpPtr) -= gdpr;
      *(here->HSM2DdPtr) += gdpr;

      *(here->HSM2SPdpPtr ) += xcsdb * s->real;
      *(here->HSM2SPdpPtr +1) += xcsdb * s->imag;
      *(here->HSM2SPdpPtr) -= gds + RevSum - gbspdp + gIstotd;
      *(here->HSM2SPgpPtr ) += xcsgb * s->real;
      *(here->HSM2SPgpPtr +1) += xcsgb * s->imag;
      *(here->HSM2SPgpPtr) -= gm - gbspg + gIstotg;
      *(here->HSM2SPspPtr ) += xcssb * s->real;
      *(here->HSM2SPspPtr +1) += xcssb * s->imag;
      *(here->HSM2SPspPtr) += gspr + gds + gbs + FwdSum + gbspsp - gIstots;
      *(here->HSM2SPsPtr) -= gspr;
      *(here->HSM2SPbpPtr ) += xcsbb * s->real;
      *(here->HSM2SPbpPtr +1) += xcsbb * s->imag;
      *(here->HSM2SPbpPtr) -= gjbs + gmbs - gbspb + gIstotb;

      *(here->HSM2SspPtr) -= gspr;
      *(here->HSM2SsPtr) += gspr;

      *(here->HSM2BPdpPtr ) += xcbdb * s->real;
      *(here->HSM2BPdpPtr +1) += xcbdb * s->imag;
      *(here->HSM2BPdpPtr) -= gjbd - gbbdp + gIbtotd;
      *(here->HSM2BPgpPtr ) += xcbgb * s->real;
      *(here->HSM2BPgpPtr +1) += xcbgb * s->imag;
      *(here->HSM2BPgpPtr) -= here->HSM2_gbgs + gIbtotg;
      *(here->HSM2BPspPtr ) += xcbsb * s->real;
      *(here->HSM2BPspPtr +1) += xcbsb * s->imag;
      *(here->HSM2BPspPtr) -= gjbs - gbbsp + gIbtots;
      *(here->HSM2BPbpPtr ) += xcbbb * s->real;
      *(here->HSM2BPbpPtr +1) += xcbbb * s->imag;
      *(here->HSM2BPbpPtr) += gjbd + gjbs - here->HSM2_gbbs - gIbtotb;

      if (model->HSM2_cogidl) {
	/* stamp gidl */
	*(here->HSM2DPdpPtr) += here->HSM2_gigidlds;
	*(here->HSM2DPgpPtr) += here->HSM2_gigidlgs;
	*(here->HSM2DPspPtr) -= (here->HSM2_gigidlgs + 
				 here->HSM2_gigidlds + here->HSM2_gigidlbs);
	*(here->HSM2DPbpPtr) += here->HSM2_gigidlbs;
	*(here->HSM2BPdpPtr) -= here->HSM2_gigidlds;
	*(here->HSM2BPgpPtr) -= here->HSM2_gigidlgs;
	*(here->HSM2BPspPtr) += (here->HSM2_gigidlgs + 
				 here->HSM2_gigidlds + here->HSM2_gigidlbs);
	*(here->HSM2BPbpPtr) -= here->HSM2_gigidlbs;
	/* stamp gisl */
	*(here->HSM2SPdpPtr) -= (here->HSM2_gigislsd + 
				 here->HSM2_gigislgd + here->HSM2_gigislbd);
	*(here->HSM2SPgpPtr) += here->HSM2_gigislgd;
	*(here->HSM2SPspPtr) += here->HSM2_gigislsd;
	*(here->HSM2SPbpPtr) += here->HSM2_gigislbd;
	*(here->HSM2BPdpPtr) += (here->HSM2_gigislgd + 
				 here->HSM2_gigislsd + here->HSM2_gigislbd);
	*(here->HSM2BPgpPtr) -= here->HSM2_gigislgd;
	*(here->HSM2BPspPtr) -= here->HSM2_gigislsd;
	*(here->HSM2BPbpPtr) -= here->HSM2_gigislbd;
      }

      if (here->HSM2_corbnet) {
	*(here->HSM2DPdbPtr ) += xcdbdb * s->real;
	*(here->HSM2DPdbPtr +1) += xcdbdb * s->imag;
	*(here->HSM2DPdbPtr) -= gbd;
	*(here->HSM2SPsbPtr ) += xcsbsb * s->real;
	*(here->HSM2SPsbPtr +1) += xcsbsb * s->imag;
	*(here->HSM2SPsbPtr) -= gbs;

	*(here->HSM2DBdpPtr ) += xcdbdb * s->real;
	*(here->HSM2DBdpPtr +1) += xcdbdb * s->imag;
	*(here->HSM2DBdpPtr) -= gbd;
	*(here->HSM2DBdbPtr ) -= xcdbdb * s->real;
	*(here->HSM2DBdbPtr +1) -= xcdbdb * s->imag;
	*(here->HSM2DBdbPtr) += gbd + here->HSM2_grbpd + here->HSM2_grbdb;
	*(here->HSM2DBbpPtr) -= here->HSM2_grbpd;
	*(here->HSM2DBbPtr) -= here->HSM2_grbdb;

	*(here->HSM2BPdbPtr) -= here->HSM2_grbpd;
	*(here->HSM2BPbPtr) -= here->HSM2_grbpb;
	*(here->HSM2BPsbPtr) -= here->HSM2_grbps;
	*(here->HSM2BPbpPtr) += here->HSM2_grbpd + here->HSM2_grbps + here->HSM2_grbpb;

	*(here->HSM2SBspPtr ) += xcsbsb * s->real;
	*(here->HSM2SBspPtr +1) += xcsbsb * s->imag;
	*(here->HSM2SBspPtr) -= gbs;
	*(here->HSM2SBbpPtr) -= here->HSM2_grbps;
	*(here->HSM2SBbPtr) -= here->HSM2_grbsb;
	*(here->HSM2SBsbPtr ) -= xcsbsb * s->real;
	*(here->HSM2SBsbPtr +1) -= xcsbsb * s->imag;
	*(here->HSM2SBsbPtr) += gbs + here->HSM2_grbps + here->HSM2_grbsb;

	*(here->HSM2BdbPtr) -= here->HSM2_grbdb;
	*(here->HSM2BbpPtr) -= here->HSM2_grbpb;
	*(here->HSM2BsbPtr) -= here->HSM2_grbsb;
	*(here->HSM2BbPtr) += here->HSM2_grbsb + here->HSM2_grbdb + here->HSM2_grbpb;
      }

    }
  }
  return(OK);
}

