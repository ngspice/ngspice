/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2acld.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "hsm2def.h"


int HSM2acLoad(
     GENmodel *inModel,
     CKTcircuit *ckt)
{
  HSM2model *model = (HSM2model*)inModel;
  HSM2instance *here;
  double xcggb_r, xcgdb_r, xcgsb_r, xcgbb_r, xcggb_i, xcgdb_i, xcgsb_i, xcgbb_i;
  double xcbgb_r, xcbdb_r, xcbsb_r, xcbbb_r, xcbgb_i, xcbdb_i, xcbsb_i, xcbbb_i;
  double xcdgb_r, xcddb_r, xcdsb_r, xcdbb_r, xcdgb_i, xcddb_i, xcdsb_i, xcdbb_i;
  double xcsgb_r, xcsdb_r, xcssb_r, xcsbb_r, xcsgb_i, xcsdb_i, xcssb_i, xcsbb_i;
  double gdpr, gspr, gds, gbd, gbs, capbd, capbs, omega;
  double FwdSum, RevSum, gm, gmbs;
  double gbspsp, gbbdp, gbbsp, gbspg, gbspb;
  double gbspdp, gbdpdp, gbdpg, gbdpb, gbdpsp;
  double gIbtotg, gIbtotd, gIbtots, gIbtotb;
  double gIgtotg, gIgtotd, gIgtots, gIgtotb;
  double gIdtotg, gIdtotd, gIdtots, gIdtotb;
  double gIstotg, gIstotd, gIstots, gIstotb;
  double cggb_real = 0.0, cgsb_real = 0.0, cgdb_real = 0.0, cggb_imag = 0.0, cgsb_imag = 0.0, cgdb_imag = 0.0;
  double cdgb_real = 0.0, cdsb_real = 0.0, cddb_real = 0.0, cdgb_imag = 0.0, cdsb_imag = 0.0, cddb_imag = 0.0;
  double csgb_real = 0.0, cssb_real = 0.0, csdb_real = 0.0, csgb_imag = 0.0, cssb_imag = 0.0, csdb_imag = 0.0;
  double cbgb_real = 0.0, cbsb_real = 0.0, cbdb_real = 0.0, cbgb_imag = 0.0, cbsb_imag = 0.0, cbdb_imag = 0.0;
  double pyggb_r = 0.0, pygdb_r = 0.0, pygsb_r = 0.0, pygbb_r = 0.0;
  double pybgb_r = 0.0, pybdb_r = 0.0, pybsb_r = 0.0, pybbb_r = 0.0;
  double pydgb_r = 0.0, pyddb_r = 0.0, pydsb_r = 0.0, pydbb_r = 0.0;
  double pysgb_r = 0.0, pysdb_r = 0.0, pyssb_r = 0.0, pysbb_r = 0.0;
  double pyggb_i = 0.0, pygdb_i = 0.0, pygsb_i = 0.0, pygbb_i = 0.0;
  double pybgb_i = 0.0, pybdb_i = 0.0, pybsb_i = 0.0, pybbb_i = 0.0;
  double pydgb_i = 0.0, pyddb_i = 0.0, pydsb_i = 0.0, pydbb_i = 0.0;
  double pysgb_i = 0.0, pysdb_i = 0.0, pyssb_i = 0.0, pysbb_i = 0.0;
  double yggb_r, ygdb_r, ygsb_r, ygbb_r, yggb_i, ygdb_i, ygsb_i, ygbb_i;
  double ybgb_r, ybdb_r, ybsb_r, ybbb_r, ybgb_i, ybdb_i, ybsb_i, ybbb_i;
  double ydgb_r, yddb_r, ydsb_r, ydbb_r, ydgb_i, yddb_i, ydsb_i, ydbb_i;
  double ysgb_r, ysdb_r, yssb_r, ysbb_r, ysgb_i, ysdb_i, yssb_i, ysbb_i;
  double grg = 0.0, pxcbdb_i = 0.0, pxcbsb_i = 0.0;

  double Qi, Qi_dVgs, Qi_dVbs, Qi_dVds ;
#ifdef DEBUG_HISIM2CGG
  double Qb ;
#endif
  double Qb_dVgs, Qb_dVbs, Qb_dVds ;
  double tau ;
  double taub ;
  double Xd, Xd_dVgs, Xd_dVbs, Xd_dVds ;
  double T1, T2, T3, T4 ;
  double cdbs_real, cgbs_real, csbs_real, cbbs_real;
  double cdbs_imag, cgbs_imag, csbs_imag, cbbs_imag;

  omega = ckt->CKTomega;
  for ( ; model != NULL; model = model->HSM2nextModel ) {
    for ( here = model->HSM2instances; here!= NULL; here = here->HSM2nextInstance ) {

      gdpr = here->HSM2drainConductance;
      gspr = here->HSM2sourceConductance;
      gds = here->HSM2_gds;
      gbd = here->HSM2_gbd;
      gbs = here->HSM2_gbs;
      capbd = here->HSM2_capbd;
      capbs = here->HSM2_capbs;

      if (model->HSM2_conqs) { /* for nqs mode */

	tau = here->HSM2_tau ;

	taub = here->HSM2_taub ;

	Xd = here->HSM2_Xd;
	Xd_dVgs = here->HSM2_Xd_dVgs ;
	Xd_dVds = here->HSM2_Xd_dVds ;
	Xd_dVbs = here->HSM2_Xd_dVbs ;
 
	Qi = here->HSM2_Qi  ;
	Qi_dVgs = here->HSM2_Qi_dVgs ;
	Qi_dVds = here->HSM2_Qi_dVds ;
	Qi_dVbs = here->HSM2_Qi_dVbs ;

#ifdef DEBUG_HISIM2CGG
	Qb = here->HSM2_Qb  ;
#endif
	Qb_dVgs = here->HSM2_Qb_dVgs ;
	Qb_dVds = here->HSM2_Qb_dVds ;
	Qb_dVbs = here->HSM2_Qb_dVbs ;

	T1 = 1.0 + (tau * omega) * (tau * omega);
	T2 = tau * omega / T1;
	T3 = 1.0 + (taub * omega) * (taub * omega);
	T4 = taub * omega / T3;

        cddb_real = Xd_dVds*Qi + Xd/T1*Qi_dVds;
        cdgb_real = Xd_dVgs*Qi + Xd/T1*Qi_dVgs;
        cdbs_real = Xd_dVbs*Qi + Xd/T1*Qi_dVbs;
        cdsb_real = - (cddb_real + cdgb_real + cdbs_real);

        cddb_imag = - T2*Xd*Qi_dVds;
        cdgb_imag = - T2*Xd*Qi_dVgs;
        cdbs_imag = - T2*Xd*Qi_dVbs;
        cdsb_imag = - (cddb_imag + cdgb_imag + cdbs_imag);

        csdb_real = - Xd_dVds*Qi + (1.0-Xd)/T1*Qi_dVds;
        csgb_real = - Xd_dVgs*Qi + (1.0-Xd)/T1*Qi_dVgs;
        csbs_real = - Xd_dVbs*Qi + (1.0-Xd)/T1*Qi_dVbs;
        cssb_real = - (csdb_real + csgb_real + csbs_real);

        csdb_imag = - T2*(1.0-Xd)*Qi_dVds;
        csgb_imag = - T2*(1.0-Xd)*Qi_dVgs;
        csbs_imag = - T2*(1.0-Xd)*Qi_dVbs;
        cssb_imag = - (csdb_imag + csgb_imag + csbs_imag);

        cbdb_real = Qb_dVds/T3;
        cbgb_real = Qb_dVgs/T3;
        cbbs_real = Qb_dVbs/T3;
        cbsb_real = - (cbdb_real + cbgb_real + cbbs_real);

        cbdb_imag = - T4*Qb_dVds;
        cbgb_imag = - T4*Qb_dVgs;
        cbbs_imag = - T4*Qb_dVbs;
        cbsb_imag = - (cbdb_imag + cbgb_imag + cbbs_imag);

        cgdb_real = - Qi_dVds/T1 - Qb_dVds/T3;
        cggb_real = - Qi_dVgs/T1 - Qb_dVgs/T3;
        cgbs_real = - Qi_dVbs/T1 - Qb_dVbs/T3;
        cgsb_real = - (cgdb_real + cggb_real + cgbs_real);

        cgdb_imag = T2*Qi_dVds + T4*Qb_dVds;
        cggb_imag = T2*Qi_dVgs + T4*Qb_dVgs;
        cgbs_imag = T2*Qi_dVbs + T4*Qb_dVbs;
        cgsb_imag = - (cgdb_imag + cggb_imag + cgbs_imag); 



#ifdef DEBUG_HISIM2CGG
	printf("Freq. %e ", omega/(2*3.14159265358979) ) ;
	printf("mag[Cgg] %e ", sqrt( cggb_real * cggb_real + cggb_imag * cggb_imag ) ) ;
	printf("qi %e ", ( sqrt(T1) / T1 ) * Qi ) ;
	printf("qb %e ", ( sqrt(T3) / T3 ) * Qb ) ;
	printf("\n") ;
#endif

#ifdef DEBUG_HISIM2AC
	printf ("#1 cssb, cggb, cgdb, cgsb, cdgb, cddb, cdsb, csgb, csdb  %e %e %e %e %e %e %e %e %e %e %e %e\n",  cggb_real, cgdb_real, cgsb_real, cdgb_real, cddb_real, cdsb_real, csgb_real, csdb_real, cssb_real);
#endif
      }
	
      if ( here->HSM2_mode >= 0 ) {
	gm = here->HSM2_gm;
	gmbs = here->HSM2_gmbs;
	FwdSum = gm + gmbs;
	RevSum = 0.0;
	
	gbbdp = -here->HSM2_gbds;
	gbbsp = here->HSM2_gbds + here->HSM2_gbgs + here->HSM2_gbbs;
	
	gbdpg = here->HSM2_gbgs;
	gbdpb = here->HSM2_gbbs;
	gbdpdp = here->HSM2_gbds;
	gbdpsp = -(gbdpg + gbdpb + gbdpdp);
	
	gbspdp = 0.0;
	gbspg = 0.0;
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
	} else {
	  gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;
	  gIstotg = gIstotd = gIstots = gIstotb = 0.0;
	  gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
	}

	if (model->HSM2_coiigs) {
	  gIgtotg = gIbtotg + gIstotg + gIdtotg;
	  gIgtotd = gIbtotd + gIstotd + gIdtotd;
	  gIgtots = gIbtots + gIstots + gIdtots;
	  gIgtotb = gIbtotb + gIstotb + gIdtotb;
	} else {
	  gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;
	}

	if (model->HSM2_conqs) { /* for nqs mode */
          /*
	    cggb_real, cgsb_real, cgdb_real
	    cggb_imag, cgsb_imag, cgdb_imag
	    cdgb_real, cdsb_real, cddb_real
	    cdgb_imag, cdsb_imag, cddb_imag
	    csgb_real, cssb_real, csdb_real
	    csgb_imag, cssb_imag, csdb_imag
	    cbgb_real, cbsb_real, cbdb_real
	    cbgb_imag, cbsb_imag, cbdb_imag
	    have already obtained.
	  */
	  if (model->HSM2_coadov == 1) { /* add overlap caps to intrinsic caps */
	    pydgb_i =  (here->HSM2_cdgo - here->HSM2_cqyg) * omega ;
	    pyddb_i =  (here->HSM2_cddo - here->HSM2_cqyd) * omega ;
	    pydsb_i =  (here->HSM2_cdso + here->HSM2_cqyg + here->HSM2_cqyd + here->HSM2_cqyb) * omega ;
	    pydbb_i = -(pydgb_i + pyddb_i + pydsb_i) ;

	    pysgb_i =  here->HSM2_csgo * omega ;
	    pysdb_i =  here->HSM2_csdo * omega ;
	    pyssb_i =  here->HSM2_csso * omega ;
	    pysbb_i = -(pysgb_i + pysdb_i + pyssb_i) ;

	    pyggb_i = (-(here->HSM2_cgdo + here->HSM2_cgbo + here->HSM2_cgso) + here->HSM2_cqyg) * omega ;
	    pygdb_i =  (here->HSM2_cgdo + here->HSM2_cqyd) * omega ;
	    pygsb_i =  (here->HSM2_cgso - here->HSM2_cqyg - here->HSM2_cqyd - here->HSM2_cqyb) * omega ;
	    pygbb_i = -(pyggb_i + pygdb_i + pygsb_i) ;
	  } 

	} else { /* for qs mode */
	  /* if coadov = 1, coverlap caps have been arleady added to intrinsic caps (QS mode)*/
	  cggb_real = here->HSM2_cggb;
	  cgsb_real = here->HSM2_cgsb;
	  cgdb_real = here->HSM2_cgdb;
	  cggb_imag = cgsb_imag = cgdb_imag = 0.0;
	
	  cbgb_real = here->HSM2_cbgb;
	  cbsb_real = here->HSM2_cbsb;
	  cbdb_real = here->HSM2_cbdb;
	  cbgb_imag = cbsb_imag = cbdb_imag = 0.0;
	  
	  cdgb_real = here->HSM2_cdgb;
	  cdsb_real = here->HSM2_cdsb;
	  cddb_real = here->HSM2_cddb;
	  cdgb_imag = cdsb_imag = cddb_imag = 0.0;
	  
	  csgb_real = -(cdgb_real + cggb_real + cbgb_real);
         cssb_real = -(cdsb_real + cgsb_real + cbsb_real);
         csdb_real = -(cddb_real + cgdb_real + cbdb_real);
         csgb_imag = cssb_imag = csdb_imag = 0.0;

         pyggb_i = 0.0; pygdb_i = 0.0; pygsb_i = 0.0; pygbb_i = 0.0;
         pybgb_i = 0.0; pybdb_i = 0.0; pybsb_i = 0.0; pybbb_i = 0.0;
         pydgb_i = 0.0; pyddb_i = 0.0; pydsb_i = 0.0; pydbb_i = 0.0;
         pysgb_i = 0.0; pysdb_i = 0.0; pyssb_i = 0.0; pysbb_i = 0.0;
       }

      } else { /* reverse mode  */
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

	gbspg = here->HSM2_gbgs;
	gbspsp = here->HSM2_gbds;
	gbspb = here->HSM2_gbbs;
	gbspdp = -(gbspg + gbspsp + gbspb);

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
	} else {
	  gIbtotg = gIbtotd = gIbtots = gIbtotb = 0.0;
	  gIstotg = gIstotd = gIstots = gIstotb = 0.0;
	  gIdtotg = gIdtotd = gIdtots = gIdtotb = 0.0;
	}
	
	if (model->HSM2_coiigs) {
	  gIgtotg = gIbtotg + gIstotg + gIdtotg;
	  gIgtotd = gIbtotd + gIstotd + gIdtotd;
	  gIgtots = gIbtots + gIstots + gIdtots;
	  gIgtotb = gIbtotb + gIstotb + gIdtotb;
	} else {
	  gIgtotg = gIgtotd = gIgtots = gIgtotb = 0.0;
	}
	
	if (model->HSM2_conqs) { /* for nqs mode */
	  /* swap d with s, s with d */
          /* cggb_real has already obtained. */
	  T1 = cgsb_real;
	  cgsb_real = cgdb_real;
	  cgdb_real = T1;
	  /* cggb_imag has already obtained. */
	  T1 = cgsb_imag;
	  cgsb_imag = cgdb_imag;
	  cgdb_imag = T1;

	  T1 = cdgb_real;
	  cdgb_real = csgb_real;
	  csgb_real = T1;
	  T1 = cdsb_real;
	  cdsb_real = csdb_real;
	  csdb_real = T1;
	  T1 = cddb_real;
	  cddb_real = cssb_real;
	  cssb_real = T1;
	  T1 = cdgb_imag;
	  cdgb_imag = csgb_imag;
	  csgb_imag = T1;
	  T1 = cdsb_imag;
	  cdsb_imag = csdb_imag;
	  csdb_imag = T1;
	  T1 = cddb_imag;
	  cddb_imag = cssb_imag;
	  cssb_imag = T1;

	  /* cbgb_real has already obtained. */
	  T1 = cbsb_real;
	  cbsb_real = cbdb_real;
	  cbdb_real = T1;
	  /* cbgb_imag has already obtained. */
	  T1 = cbsb_imag;
	  cbsb_imag = cbdb_imag;
	  cbdb_imag = T1;

	  if (model->HSM2_coadov == 1) { /* add overlap caps to intrinsic caps */
	    pydgb_i = here->HSM2_csgo * omega ;
	    pyddb_i = here->HSM2_csso * omega ;
	    pydsb_i = here->HSM2_csdo * omega ;
	    pydbb_i = -(pydgb_i + pyddb_i + pydsb_i) ;

	    pysgb_i = (here->HSM2_cdgo - here->HSM2_cqyg) * omega ;
	    pysdb_i = (here->HSM2_cdso + here->HSM2_cqyg + here->HSM2_cqyd + here->HSM2_cqyb) * omega ;
	    pyssb_i = (here->HSM2_cddo - here->HSM2_cqyd) * omega ;
	    pysbb_i = -(pysgb_i + pysdb_i + pyssb_i) ;

	    pyggb_i = (-(here->HSM2_cgdo + here->HSM2_cgbo + here->HSM2_cgso) + here->HSM2_cqyg) * omega ;
	    pygdb_i = (here->HSM2_cgso - here->HSM2_cqyg - here->HSM2_cqyd - here->HSM2_cqyb) * omega ;
	    pygsb_i = (here->HSM2_cgdo + here->HSM2_cqyd) * omega ;
	    pygbb_i = -(pyggb_i + pygdb_i + pygsb_i) ;
	  } 

	} else { /* for qs mode */
	  /* if coadov = 1, coverlap caps have already been added to intrinsic caps (QS mode)*/
	  cggb_real = here->HSM2_cggb;
	  cgsb_real = here->HSM2_cgdb;
	  cgdb_real = here->HSM2_cgsb;
	  cggb_imag = cgsb_imag = cgdb_imag = 0.0;
	
	  cbgb_real = here->HSM2_cbgb;
	  cbsb_real = here->HSM2_cbdb;
	  cbdb_real = here->HSM2_cbsb;
	  cbgb_imag = cbsb_imag = cbdb_imag = 0.0;
	  
	  csgb_real = here->HSM2_cdgb ;
	  cssb_real = here->HSM2_cddb ;
	  csdb_real = here->HSM2_cdsb ;
	  csgb_imag = cssb_imag = csdb_imag = 0.0;
	  
	  cdgb_real = -(csgb_real + cggb_real + cbgb_real);
         cdsb_real = -(cssb_real + cgsb_real + cbsb_real);
         cddb_real = -(csdb_real + cgdb_real + cbdb_real);
         cdgb_imag = cdsb_imag = cddb_imag = 0.0;

         pyggb_i = 0.0; pygdb_i = 0.0; pygsb_i = 0.0; pygbb_i = 0.0;
         pybgb_i = 0.0; pybdb_i = 0.0; pybsb_i = 0.0; pybbb_i = 0.0;
         pydgb_i = 0.0; pyddb_i = 0.0; pydsb_i = 0.0; pydbb_i = 0.0;
         pysgb_i = 0.0; pysdb_i = 0.0; pyssb_i = 0.0; pysbb_i = 0.0;
       }

      }
#ifdef DEBUG_HISIM2AC
      printf ("#2 cssb, cggb, cgdb, cgsb, cdgb, cddb, cdsb, csgb, csdb  %e %e %e %e %e %e %e %e %e %e %e %e\n",  cggb_real, cgdb_real, cgsb_real, cdgb_real, cddb_real, cdsb_real, csgb_real, csdb_real, cssb_real);
#endif

      /* matrix elements for ac analysis (including real and imaginary parts) */
      xcdgb_r = cdgb_real * omega;
      xcddb_r = cddb_real * omega;
      xcdsb_r = cdsb_real * omega;
      xcdbb_r = -(xcdgb_r + xcddb_r + xcdsb_r);

      xcsgb_r = csgb_real * omega;
      xcsdb_r = csdb_real * omega;
      xcssb_r = cssb_real * omega;
      xcsbb_r = -(xcsgb_r + xcsdb_r + xcssb_r);

      xcggb_r = cggb_real * omega;
      xcgdb_r = cgdb_real * omega;
      xcgsb_r = cgsb_real * omega;
      xcgbb_r = -(xcggb_r + xcgdb_r + xcgsb_r);
      
      xcbgb_r = cbgb_real * omega;
      xcbdb_r = cbdb_real * omega;
      xcbsb_r = cbsb_real * omega;
      xcbbb_r = -(xcbgb_r + xcbdb_r + xcbsb_r);

      xcdgb_i = cdgb_imag * omega;
      xcddb_i = cddb_imag * omega;
      xcdsb_i = cdsb_imag * omega;
      xcdbb_i = -(xcddb_i + xcdgb_i + xcdsb_i);

      xcsgb_i = csgb_imag * omega;
      xcsdb_i = csdb_imag * omega;
      xcssb_i = cssb_imag * omega;
      xcsbb_i = -(xcsdb_i + xcsgb_i + xcssb_i);

      xcggb_i = cggb_imag * omega;
      xcgdb_i = cgdb_imag * omega;
      xcgsb_i = cgsb_imag * omega;
      xcgbb_i = -(xcggb_i + xcgdb_i + xcgsb_i);

      xcbgb_i = cbgb_imag * omega;
      xcbdb_i = cbdb_imag * omega;
      xcbsb_i = cbsb_imag * omega;
      xcbbb_i = -(xcbgb_i + xcbdb_i + xcbsb_i);

      /* stamp intrinsic y-parameters */
      yggb_r = - xcggb_i; yggb_i = xcggb_r;
      ygdb_r = - xcgdb_i; ygdb_i = xcgdb_r;
      ygsb_r = - xcgsb_i; ygsb_i = xcgsb_r;
      ygbb_r = - xcgbb_i; ygbb_i = xcgbb_r;

      ydgb_r = - xcdgb_i; ydgb_i = xcdgb_r;
      yddb_r = - xcddb_i; yddb_i = xcddb_r;
      ydsb_r = - xcdsb_i; ydsb_i = xcdsb_r;
      ydbb_r = - xcdbb_i; ydbb_i = xcdbb_r;
      ydgb_r += gm;
      yddb_r += gds + RevSum;
      ydsb_r += - gds - FwdSum;
      ydbb_r += gmbs;

      ysgb_r = - xcsgb_i; ysgb_i = xcsgb_r;
      ysdb_r = - xcsdb_i; ysdb_i = xcsdb_r;
      yssb_r = - xcssb_i; yssb_i = xcssb_r;
      ysbb_r = - xcsbb_i; ysbb_i = xcsbb_r;
      ysgb_r += - gm;
      ysdb_r += - gds - RevSum;
      yssb_r += gds + FwdSum;
      ysbb_r += - gmbs;

      ybgb_r = - xcbgb_i; ybgb_i = xcbgb_r;
      ybdb_r = - xcbdb_i; ybdb_i = xcbdb_r;
      ybsb_r = - xcbsb_i; ybsb_i = xcbsb_r;
      ybbb_r = - xcbbb_i; ybbb_i = xcbbb_r;
      
      /* Ibd, Ibs, Igate, Igd, Igs, Igb, Igidl, Igisl, Isub */
      pydgb_r = gbdpg - gIdtotg + here->HSM2_gigidlgs ;
      pyddb_r = gbd + gbdpdp - gIdtotd + here->HSM2_gigidlds ;
      pydsb_r = gbdpsp - gIdtots - (here->HSM2_gigidlgs + here->HSM2_gigidlds + here->HSM2_gigidlbs);
      pydbb_r = gbdpb - gIdtotb + here->HSM2_gigidlbs ;
      if (!here->HSM2_corbnet) pydbb_r += - gbd;

      pysgb_r = gbspg - gIstotg + here->HSM2_gigislgd ;
      pysdb_r = gbspdp - gIstotd - (here->HSM2_gigislsd + here->HSM2_gigislgd + here->HSM2_gigislbd);
      pyssb_r = gbs + gbspsp - gIstots + here->HSM2_gigislsd ;
      pysbb_r =gbspb - gIstotb + here->HSM2_gigislbd ;
      if (!here->HSM2_corbnet) pysbb_r += - gbs;
	
      pyggb_r = gIgtotg ;
      if (here->HSM2_corg == 1) {	
	grg = here->HSM2_grg;
	pyggb_r += grg;
      }
      pygdb_r = gIgtotd ;
      pygsb_r = gIgtots ;
      pygbb_r = gIgtotb ;

      pybgb_r = - here->HSM2_gbgs - gIbtotg - here->HSM2_gigidlgs - here->HSM2_gigislgd ;
      pybdb_r = gbbdp - gIbtotd 
	- here->HSM2_gigidlds + (here->HSM2_gigislgd + here->HSM2_gigislsd + here->HSM2_gigislbd) ;
      if (!here->HSM2_corbnet) pybdb_r += - gbd ;
      pybsb_r = gbbsp - gIbtots 
	+ (here->HSM2_gigidlgs + here->HSM2_gigidlds + here->HSM2_gigidlbs) - here->HSM2_gigislsd ;
      if (!here->HSM2_corbnet) pybsb_r += - gbs ;
      pybbb_r = - here->HSM2_gbbs - gIbtotb - here->HSM2_gigidlbs - here->HSM2_gigislbd ;
      if (!here->HSM2_corbnet) pybbb_r += gbd + gbs ;

      pybdb_i = -(pyddb_i + pygdb_i + pysdb_i); 
      pybgb_i = -(pydgb_i + pyggb_i + pysgb_i); 
      pybsb_i = -(pydsb_i + pygsb_i + pyssb_i); 
      pybbb_i = -(pydbb_i + pygbb_i + pysbb_i); 
      
      /* Cbd, Cbs */
      pyddb_i += capbd * omega ;
      pyssb_i += capbs * omega ;
      if (!here->HSM2_corbnet) {
	pydbb_i -= capbd * omega ; 
	pysbb_i -= capbs * omega ; 

 	pybdb_i -= capbd * omega ; 
 	pybsb_i -= capbs * omega ;  
 	pybbb_i += (capbd + capbs) * omega ; 
      } else {
	pxcbdb_i = - capbd * omega ;
	pxcbsb_i = - capbs * omega ;
      }

#ifdef DEBUG_HISIM2AC
      /* for representing y-parameters */ 
      printf("f ygg_r ygg_i %e %e %e\n",omega/(2.0*3.141592653589793),yggb_r+pyggb_r,yggb_i+pyggb_i);
      printf("f ygd_r ygd_i %e %e %e\n",omega/(2.0*3.141592653589793),ygdb_r+pygdb_r,ygdb_i+pygdb_i);
      printf("f ygs_r ygs_i %e %e %e\n",omega/(2.0*3.141592653589793),ygsb_r+pygsb_r,ygsb_i+pygsb_i);
      printf("f ygb_r ygb_i %e %e %e\n",omega/(2.0*3.141592653589793),ygbb_r+pygbb_r,ygbb_i+pygbb_i);

      printf("f ydg_r ydg_i %e %e %e\n",omega/(2.0*3.141592653589793),ydgb_r+pydgb_r,ydgb_i+pydgb_i);
      printf("f ydd_r ydd_i %e %e %e\n",omega/(2.0*3.141592653589793),yddb_r+pyddb_r,yddb_i+pyddb_i);
      printf("f yds_r yds_i %e %e %e\n",omega/(2.0*3.141592653589793),ydsb_r+pydsb_r,ydsb_i+pydsb_i);
      printf("f ydb_r ydb_i %e %e %e\n",omega/(2.0*3.141592653589793),ydbb_r+pydbb_r,ydbb_i+pydbb_i);

      printf("f ybg_r ybg_i %e %e %e\n",omega/(2.0*3.141592653589793),ybgb_r+pybgb_r,ybgb_i+pybgb_i);
      printf("f ybd_r ybd_i %e %e %e\n",omega/(2.0*3.141592653589793),ybdb_r+pybdb_r,ybdb_i+pybdb_i);
      printf("f ybs_r ybs_i %e %e %e\n",omega/(2.0*3.141592653589793),ybsb_r+pybsb_r,ybsb_i+pybsb_i);
      printf("f ybb_r ybb_i %e %e %e\n",omega/(2.0*3.141592653589793),ybbb_r+pybbb_r,ybbb_i+pybbb_i);

      printf("f ysg_r ysg_i %e %e %e\n",omega/(2.0*3.141592653589793),ysgb_r+pysgb_r,ysgb_i+pysgb_i);
      printf("f ysd_r ysd_i %e %e %e\n",omega/(2.0*3.141592653589793),ysdb_r+pysdb_r,ysdb_i+pysdb_i);
      printf("f yss_r yss_i %e %e %e\n",omega/(2.0*3.141592653589793),yssb_r+pyssb_r,yssb_i+pyssb_i);
      printf("f ysb_r ysb_i %e %e %e\n",omega/(2.0*3.141592653589793),ysbb_r+pysbb_r,ysbb_i+pysbb_i);

      printf("f y11r y11i y12r y12i y21r y21i y22r y22i %e %e %e %e %e %e %e %e %e\n",omega/(2.0*3.141592653589793),yggb_r+pyggb_r,yggb_i+pyggb_i, ygdb_r+pygdb_r,ygdb_i+pygdb_i, ydgb_r+pydgb_r,ydgb_i+pydgb_i, yddb_r+pyddb_r,yddb_i+pyddb_i);
#endif
      
      if (here->HSM2_corg == 1) {
	*(here->HSM2GgPtr) += grg;
	*(here->HSM2GPgPtr) -= grg;
	*(here->HSM2GgpPtr) -= grg;
      } 

      *(here->HSM2GPgpPtr +1)  += yggb_i + pyggb_i;
      *(here->HSM2GPgpPtr)     += yggb_r + pyggb_r;
      *(here->HSM2GPdpPtr +1)  += ygdb_i + pygdb_i;
      *(here->HSM2GPdpPtr)     += ygdb_r + pygdb_r;
      *(here->HSM2GPspPtr +1)  += ygsb_i + pygsb_i;
      *(here->HSM2GPspPtr)     += ygsb_r + pygsb_r;
      *(here->HSM2GPbpPtr +1)  += ygbb_i + pygbb_i;
      *(here->HSM2GPbpPtr)     += ygbb_r + pygbb_r;

      *(here->HSM2DPdpPtr +1) += yddb_i + pyddb_i;
      *(here->HSM2DPdpPtr)    += yddb_r + pyddb_r + gdpr;
      *(here->HSM2DPdPtr)     -= gdpr;
      *(here->HSM2DPgpPtr +1) += ydgb_i + pydgb_i;
      *(here->HSM2DPgpPtr)    += ydgb_r + pydgb_r;
      *(here->HSM2DPspPtr +1) += ydsb_i + pydsb_i;
      *(here->HSM2DPspPtr)    += ydsb_r + pydsb_r;
      *(here->HSM2DPbpPtr +1) += ydbb_i + pydbb_i;
      *(here->HSM2DPbpPtr)    += ydbb_r + pydbb_r;

      *(here->HSM2DdpPtr)     -= gdpr;
      *(here->HSM2DdPtr)      += gdpr;

      *(here->HSM2SPdpPtr +1) += ysdb_i + pysdb_i;
      *(here->HSM2SPdpPtr)    += ysdb_r + pysdb_r;
      *(here->HSM2SPgpPtr +1) += ysgb_i + pysgb_i;
      *(here->HSM2SPgpPtr)    += ysgb_r + pysgb_r;
      *(here->HSM2SPspPtr +1) += yssb_i + pyssb_i;
      *(here->HSM2SPspPtr)    += yssb_r + pyssb_r + gspr;
      *(here->HSM2SPsPtr)     -= gspr ;
      *(here->HSM2SPbpPtr +1) += ysbb_i + pysbb_i;
      *(here->HSM2SPbpPtr)    += ysbb_r + pysbb_r;

      *(here->HSM2SspPtr)     -= gspr;
      *(here->HSM2SsPtr)      += gspr;

      *(here->HSM2BPdpPtr +1)  += ybdb_i + pybdb_i;
      *(here->HSM2BPdpPtr)     += ybdb_r + pybdb_r;
      *(here->HSM2BPgpPtr +1)  += ybgb_i + pybgb_i;
      *(here->HSM2BPgpPtr)     += ybgb_r + pybgb_r;
      *(here->HSM2BPspPtr +1)  += ybsb_i + pybsb_i;
      *(here->HSM2BPspPtr)     += ybsb_r + pybsb_r;
      *(here->HSM2BPbpPtr +1)  += ybbb_i + pybbb_i;
      *(here->HSM2BPbpPtr)     += ybbb_r + pybbb_r;

      if (here->HSM2_corbnet == 1) {
	*(here->HSM2DPdbPtr +1) += pxcbdb_i;
	*(here->HSM2DPdbPtr) -= gbd;
	*(here->HSM2SPsbPtr +1) += pxcbsb_i;
	*(here->HSM2SPsbPtr) -= gbs;
	
	*(here->HSM2DBdpPtr +1) += pxcbdb_i;
	*(here->HSM2DBdpPtr) -= gbd;
	*(here->HSM2DBdbPtr +1) -= pxcbdb_i;
	*(here->HSM2DBdbPtr) += gbd + here->HSM2_grbpd + here->HSM2_grbdb;
	*(here->HSM2DBbpPtr) -= here->HSM2_grbpd;
	*(here->HSM2DBbPtr) -= here->HSM2_grbdb;

	*(here->HSM2BPdbPtr) -= here->HSM2_grbpd;
	*(here->HSM2BPbPtr) -= here->HSM2_grbpb;
	*(here->HSM2BPsbPtr) -= here->HSM2_grbps;
	*(here->HSM2BPbpPtr) += here->HSM2_grbpd + here->HSM2_grbps + here->HSM2_grbpb;

	*(here->HSM2SBspPtr +1) += pxcbsb_i;
	*(here->HSM2SBspPtr) -= gbs;
	*(here->HSM2SBbpPtr) -= here->HSM2_grbps;
	*(here->HSM2SBbPtr) -= here->HSM2_grbsb;
	*(here->HSM2SBsbPtr +1) -= pxcbsb_i;
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
