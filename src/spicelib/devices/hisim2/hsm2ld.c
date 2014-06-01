/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2ld.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "hsm2def.h"
#include "hisim2.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"

#define SHOW_EPS_QUANT 1.0e-15
#define BYP_TOL_FACTOR model->HSM2_byptol  

#ifdef MOS_MODEL_TIME
#ifdef USE_OMP
#error "MOS_MODEL_TIME is not supported when USE_OMP is active"
#endif
/** MOS Model Time **/
#include <sys/time.h>
extern char *mos_model_name ;
extern double mos_model_time ;
double gtodsecld(void) {
  struct timeval tv;
  const double sec2000 = 9.46e8 ;
  gettimeofday(&tv, NULL);
  return ( tv.tv_sec - sec2000 ) + (double)tv.tv_usec*1e-6;
}
double tm0 , tm1 ;
#ifdef PARAMOS_TIME
#include <time.h>
double vsum ;
static double vsum0 = 1.0e5 ;
#endif
#endif

#ifdef USE_OMP
int HSM2LoadOMP(HSM2instance *here, CKTcircuit *ckt);
void HSM2LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt);
#endif

static void ShowPhysVals
(
 HSM2instance *here,
 HSM2model *model,
 int isFirst,
 double vds,
 double vgs,
 double vbs,
 double vgd,
 double vbd,
 double vgb
)
{
  NG_IGNORE(vgd);
  NG_IGNORE(vbd);
  /* regard the epsilon-quantity as 0.0 */
  vds = (fabs(vds) < SHOW_EPS_QUANT) ? 0.0 : vds;
  vgs = (fabs(vgs) < SHOW_EPS_QUANT) ? 0.0 : vgs;
  vbs = (fabs(vbs) < SHOW_EPS_QUANT) ? 0.0 : vbs;
  vgb = (fabs(vgb) < SHOW_EPS_QUANT) ? 0.0 : vgb;
  switch (model->HSM2_show) {
  case 1:
    if (isFirst) printf("Vds        Ids\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_mode*here->HSM2_ids);
    break;
  case 2:
    if (isFirst) printf("Vgs        Ids\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_mode*here->HSM2_ids);
    break;
  case 3:
    if (isFirst) printf("Vgs        log10(|Ids|)\n");
    printf("%e %e\n", model->HSM2_type*vgs, log10(here->HSM2_ids));
    break;
  case 4:
    if (isFirst) printf("log10(|Ids|)    gm/|Ids|\n");
    if (here->HSM2_ids == 0.0)
      printf("I can't show gm/Ids - log10(Ids), because Ids = 0.\n");
    else
      printf("%e %e\n",  log10(here->HSM2_ids), here->HSM2_gm/here->HSM2_ids);
    break;
  case 5:
    if (isFirst) printf("Vds        gds\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_gds);
    break;
  case 6:
    if (isFirst) printf("Vgs        gm\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_gm);
    break;
  case 7:
    if (isFirst) printf("Vbs        gbs\n");
    printf("%e %e\n", model->HSM2_type*vbs, here->HSM2_gmbs);
    break;
  case 8:
    if (isFirst) printf("Vgs        Cgg\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_cggb);
    break;
  case 9:
    if (isFirst) printf("Vgs        Cgs\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_cgsb);
    break;
  case 10:
    if (isFirst) printf("Vgs        Cgd\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_cgdb);
    break;
  case 11:
    if (isFirst) printf("Vgs        Cgb\n");
    printf("%e %e\n", model->HSM2_type*vgs, -(here->HSM2_cggb+here->HSM2_cgsb+here->HSM2_cgdb));
    break;
  case 12:
    if (isFirst) printf("Vds        Csg\n");
    printf("%e %e\n", model->HSM2_type*vds, -(here->HSM2_cggb+here->HSM2_cbgb+here->HSM2_cdgb));
    break;
  case 13:
    if (isFirst) printf("Vds        Cdg\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_cdgb);
    break;
  case 14:
    if (isFirst) printf("Vds        Cbg\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_cbgb);
    break;
  case 15:
    if (isFirst) printf("Vds        Cgg\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_cggb);
    break;
  case 16:
    if (isFirst) printf("Vds        Cgs\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_cgsb);
    break;
  case 17:
    if (isFirst) printf("Vds        Cgd\n");
    printf("%e %e\n", model->HSM2_type*vds, here->HSM2_cgdb);
    break;
  case 18:
    if (isFirst) printf("Vds        Cgb\n");
    printf("%e %e\n", model->HSM2_type*vds, -(here->HSM2_cggb+here->HSM2_cgsb+here->HSM2_cgdb));
    break;
  case 19:
    if (isFirst) printf("Vgs        Csg\n");
    printf("%e %e\n", model->HSM2_type*vgs, -(here->HSM2_cggb+here->HSM2_cbgb+here->HSM2_cdgb));
    break;
  case 20:
    if (isFirst) printf("Vgs        Cdg\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_cdgb);
    break;
  case 21:
    if (isFirst) printf("Vgs        Cbg\n");
    printf("%e %e\n", model->HSM2_type*vgs, here->HSM2_cbgb);
    break;
  case 22:
    if (isFirst) printf("Vgb        Cgb\n");
    printf("%e %e\n", model->HSM2_type*vgb, -(here->HSM2_cggb+here->HSM2_cgsb+here->HSM2_cgdb));
    break;
  case 50:
    if (isFirst) printf("Vgs  Vds  Vbs  Vgb  Ids  log10(|Ids|)  gm/|Ids|  gm  gds  gbs  Cgg  Cgs  Cgb  Cgd  Csg  Cbg  Cdg\n");
    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", model->HSM2_type*vgs, model->HSM2_type*vds, model->HSM2_type*vbs, model->HSM2_type*vgb, here->HSM2_mode*here->HSM2_ids, log10(here->HSM2_ids), here->HSM2_gm/here->HSM2_ids, here->HSM2_gm, here->HSM2_gds, here->HSM2_gmbs, here->HSM2_cggb, here->HSM2_cgsb, -(here->HSM2_cggb+here->HSM2_cgsb+here->HSM2_cgdb), here->HSM2_cgdb, -(here->HSM2_cggb+here->HSM2_cbgb+here->HSM2_cdgb), here->HSM2_cbgb, here->HSM2_cdgb);
    break;
  default:
    /*
    printf("There is no physical value corrsponding to %d\n", flag);
    */
    break;
  }
}

int HSM2load(
     GENmodel *inModel,
     CKTcircuit *ckt)
     /* actually load the current value into the 
      * sparse matrix previously provided 
      */
{
#ifdef USE_OMP
    int idx;
    HSM2model *model = (HSM2model*)inModel;
    int good = 0;
    HSM2instance **InstArray;
    InstArray = model->HSM2InstanceArray;

#pragma omp parallel for
    for (idx = 0; idx < model->HSM2InstCount; idx++) {
        HSM2instance *here = InstArray[idx];
        int local_good = HSM2LoadOMP(here, ckt);
        if (local_good)
            good = local_good;
    }

    HSM2LoadRhsMat(inModel, ckt);

    return good;
}

int HSM2LoadOMP(HSM2instance *here, CKTcircuit *ckt)
{
  HSM2model *model;
#else
  HSM2model *model = (HSM2model*)inModel;
  HSM2instance *here;
#endif
/*  HSM2binningParam *pParam;*/
  double cbhat=0.0, cdrain=0.0, cdhat=0.0, cdreq=0.0, cgbhat=0.0, cgshat=0.0, cgdhat=0.0 ;
  double Ibtot=0.0, Idtot=0.0, Igbtot=0.0, Igstot=0.0, Igdtot=0.0 ;
  double ceq=0.0, ceqbd=0.0, ceqbs=0.0, ceqqb=0.0, ceqqd=0.0, ceqqg=0.0 ;
  double ceqjs=0.0, ceqjd=0.0, ceqqjs=0.0, ceqqjd=0.0 ;
  double delvbd=0.0, delvbs=0.0, delvds=0.0, delvgd=0.0, delvgs=0.0 ;
  double gcbdb=0.0, gcbgb=0.0, gcbsb=0.0, gcddb=0.0, gcdgb=0.0, gcdsb=0.0 ;
  double gcgdb=0.0, gcggb=0.0, gcgsb=0.0, gcgbb=0.0, gcsdb=0.0, gcsgb=0.0, gcssb=0.0 ;
  double geq=0.0, xfact=0.0 ;
  double vbd=0.0, vbs=0.0, vds=0.0, vgb=0.0, vgd=0.0, vgdo=0.0, vgs=0.0, von=0.0 ;
  double gbbdp=0.0, gbbsp=0.0, gbspg=0.0, gbspdp=0.0, gbspb=0.0, gbspsp=0.0 ;
  double qgate=0.0, qbulk=0.0, qdrn=0.0 ;
  double cqgate=0.0, cqbulk=0.0, cqdrn=0.0 ;
  double gbdpdp=0.0, gbdpg=0.0, gbdpb=0.0, gbdpsp=0.0; 
  double gm=0.0, gmbs=0.0, FwdSum=0.0, RevSum=0.0 ;
  double ag0=0.0 ;
  double Ibtoteq=0.0, gIbtotg=0.0, gIbtotd=0.0, gIbtots=0.0, gIbtotb=0.0 ;
  double Igtoteq=0.0, gIgtotg=0.0, gIgtotd=0.0, gIgtots=0.0, gIgtotb=0.0 ;
  double Idtoteq=0.0, gIdtotg=0.0, gIdtotd=0.0, gIdtots=0.0, gIdtotb=0.0 ;
  double Istoteq=0.0, gIstotg=0.0, gIstotd=0.0, gIstots=0.0, gIstotb=0.0 ;
  double ivds=0.0, ivgs=0.0, ivbs=0.0 ;
  double gjbs=0.0, gjbd=0.0, gcdbdb=0.0, gcsbsb=0.0, gcbbb=0.0, gcdbb=0.0, gcsbb=0.0, grg=0.0 ;
  double vdbs=0.0, vsbs=0.0, vdbd=0.0, delvdbs=0.0, delvsbs=0.0, delvdbd=0.0 ;
  double vges=0.0, vged=0.0, delvges=0.0,/* delvged=0.0,*/ vgedo=0.0 ;
  double vsbdo=0.0, vsbd=0.0; 
  double vbs_jct=0.0, vbd_jct=0.0, delvbs_jct=0.0, delvbd_jct=0.0 ;
  int ByPass=0, Check=0, Check1=0, Check2=0 ;
  int BYPASS_enable =0 ;
#ifndef NOBYPASS
  double tempv=0.0 ;
#endif /*NOBYPASS*/
#ifndef NEWCONV
  double tol=0.0, tol2=0.0, tol3=0.0, tol4=0.0 ;
#endif
  int ChargeComputationNeeded =  
    ((ckt->CKTmode & (MODEAC | MODETRAN | MODEINITSMSIG)) ||
     ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)))
    ? 1 : 0;
  int showPhysVal;
  int isConv;
  double vds_pre = 0.0;
  double reltol, abstol , voltTol ;
  
#ifdef MOS_MODEL_TIME
tm0 = gtodsecld() ;
#endif


#ifdef USE_OMP
    model = here->HSM2modPtr;
    reltol = ckt->CKTreltol * BYP_TOL_FACTOR ;
    abstol = ckt->CKTabstol * BYP_TOL_FACTOR ;
    voltTol= ckt->CKTvoltTol* BYP_TOL_FACTOR ;
    BYPASS_enable = (BYP_TOL_FACTOR > 0.0 && ckt->CKTbypass) ;
    model->HSM2_bypass_enable = BYPASS_enable ;
#else
  /*  loop through all the HSM2 device models */
  for ( ; model != NULL; model = model->HSM2nextModel ) {
    /* loop through all the instances of the model */

    reltol = ckt->CKTreltol * BYP_TOL_FACTOR ; 
    abstol = ckt->CKTabstol * BYP_TOL_FACTOR ;
    voltTol= ckt->CKTvoltTol* BYP_TOL_FACTOR ;
    BYPASS_enable = (BYP_TOL_FACTOR > 0.0 && ckt->CKTbypass) ;
    model->HSM2_bypass_enable = BYPASS_enable ;

    for (here = model->HSM2instances; here != NULL ;
	 here = here->HSM2nextInstance) {
#endif
/*      pParam = &here->pParam ;*/
      showPhysVal = 0;
      Check=1;
      ByPass = 0;

#ifdef DEBUG_HISIM2LD_VX
      printf("mode = %x\n", ckt->CKTmode);
      printf("Vd Vg Vs Vb %e %e %e %e\n", *(ckt->CKTrhsOld+here->HSM2dNodePrime),
	     *(ckt->CKTrhsOld+here->HSM2gNodePrime),
	     *(ckt->CKTrhsOld+here->HSM2sNodePrime),
	     *(ckt->CKTrhsOld+here->HSM2bNodePrime));
#endif

      if ( ckt->CKTmode & MODEINITSMSIG ) {
	vbs = *(ckt->CKTstate0 + here->HSM2vbs);
	vgs = *(ckt->CKTstate0 + here->HSM2vgs);
	vds = *(ckt->CKTstate0 + here->HSM2vds);

	vges = *(ckt->CKTstate0 + here->HSM2vges);
	vdbs = *(ckt->CKTstate0 + here->HSM2vdbs);
	vsbs = *(ckt->CKTstate0 + here->HSM2vsbs);
      } 
      else if ( ckt->CKTmode & MODEINITTRAN ) {
	vbs = *(ckt->CKTstate1 + here->HSM2vbs);
	vgs = *(ckt->CKTstate1 + here->HSM2vgs);
	vds = *(ckt->CKTstate1 + here->HSM2vds);

	vges = *(ckt->CKTstate1 + here->HSM2vges);
	vdbs = *(ckt->CKTstate1 + here->HSM2vdbs);
	vsbs = *(ckt->CKTstate1 + here->HSM2vsbs);
      } 
      else if ( (ckt->CKTmode & MODEINITJCT) && !here->HSM2_off ) {
	vds = model->HSM2_type * here->HSM2_icVDS;
	vgs = vges = model->HSM2_type * here->HSM2_icVGS;
	vbs = vdbs = vsbs = model->HSM2_type * here->HSM2_icVBS;
	if ( (vds == 0.0) && (vgs == 0.0) && (vbs == 0.0) && 
	     ( (ckt->CKTmode & (MODETRAN|MODEAC|MODEDCOP|MODEDCTRANCURVE)) ||
	       !(ckt->CKTmode & MODEUIC) ) ) { 
	  /* set biases for starting analysis */
	  vbs = vdbs = vsbs = 0.0;
	  /*
	  vgs = vges = model->HSM2_type * pParam->HSM2_vfbc + 0.1;
	  */
	  vgs = vges = 0.1;
	  vds = 0.1;
	}
      } 
      else if ( ( ckt->CKTmode & (MODEINITJCT | MODEINITFIX) ) && 
		here->HSM2_off ) {
	vbs = vgs = vds = 0.0; vges = 0.0; vdbs = vsbs = 0.0;
      } 
      else {
#ifndef PREDICTOR /* BSIM3 style */
	if (ckt->CKTmode & MODEINITPRED) {
	  xfact = ckt->CKTdelta / ckt->CKTdeltaOld[1];
	  *(ckt->CKTstate0 + here->HSM2vbs) = 
	    *(ckt->CKTstate1 + here->HSM2vbs);
	  vbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vbs))
	    -(xfact * (*(ckt->CKTstate2 + here->HSM2vbs)));
	  *(ckt->CKTstate0 + here->HSM2vgs) = 
	    *(ckt->CKTstate1 + here->HSM2vgs);
	  vgs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vgs))
	    -(xfact * (*(ckt->CKTstate2 + here->HSM2vgs)));
	  *(ckt->CKTstate0 + here->HSM2vds) = 
	    *(ckt->CKTstate1 + here->HSM2vds);
	  vds = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vds))
	    -(xfact * (*(ckt->CKTstate2 + here->HSM2vds)));
	  *(ckt->CKTstate0 + here->HSM2vbd) = 
	    *(ckt->CKTstate0 + here->HSM2vbs)-
	    *(ckt->CKTstate0 + here->HSM2vds);

	  *(ckt->CKTstate0 + here->HSM2vges) = 
	    *(ckt->CKTstate1 + here->HSM2vges);
	  vges = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vges))
	    -(xfact * (*(ckt->CKTstate2 + here->HSM2vges)));
	  *(ckt->CKTstate0 + here->HSM2vdbs) =
	    *(ckt->CKTstate1 + here->HSM2vdbs);
	  vdbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vdbs))
	    - (xfact * (*(ckt->CKTstate2 + here->HSM2vdbs)));
	  *(ckt->CKTstate0 + here->HSM2vdbd) =
	    *(ckt->CKTstate0 + here->HSM2vdbs)
	    - *(ckt->CKTstate0 + here->HSM2vds);
	  *(ckt->CKTstate0 + here->HSM2vsbs) =
	    *(ckt->CKTstate1 + here->HSM2vsbs);
	  vsbs = (1.0 + xfact)* (*(ckt->CKTstate1 + here->HSM2vsbs))
	    - (xfact * (*(ckt->CKTstate2 + here->HSM2vsbs)));
	} 
	else {
#endif /* PREDICTOR */
	  /* get biases from CKT */
	  vbs = model->HSM2_type * 
	    (*(ckt->CKTrhsOld+here->HSM2bNodePrime) -
	     *(ckt->CKTrhsOld+here->HSM2sNodePrime));
	  vgs = model->HSM2_type * 
	    (*(ckt->CKTrhsOld+here->HSM2gNodePrime) -
	     *(ckt->CKTrhsOld+here->HSM2sNodePrime));
	  vds = model->HSM2_type * 
	    (*(ckt->CKTrhsOld+here->HSM2dNodePrime) -
	     *(ckt->CKTrhsOld+here->HSM2sNodePrime));

	  vges = model->HSM2_type * 
	    (*(ckt->CKTrhsOld+here->HSM2gNode) -
	     *(ckt->CKTrhsOld+here->HSM2sNodePrime));
	  vdbs = model->HSM2_type
	    * (*(ckt->CKTrhsOld + here->HSM2dbNode)
	       - *(ckt->CKTrhsOld + here->HSM2sNodePrime));
	  vsbs = model->HSM2_type
	    * (*(ckt->CKTrhsOld + here->HSM2sbNode)
	       - *(ckt->CKTrhsOld + here->HSM2sNodePrime));
#ifndef PREDICTOR
	}
#endif /* PREDICTOR */

	vbd = vbs - vds;
	vgd = vgs - vds;
	vged = vges - vds;
	vdbd = vdbs - vds;
	vgdo = *(ckt->CKTstate0 + here->HSM2vgs) - *(ckt->CKTstate0 + here->HSM2vds);
	vgedo = *(ckt->CKTstate0 + here->HSM2vges) - *(ckt->CKTstate0 + here->HSM2vds);
	delvbs = vbs - *(ckt->CKTstate0 + here->HSM2vbs);
	delvbd = vbd - *(ckt->CKTstate0 + here->HSM2vbd);
	delvgs = vgs - *(ckt->CKTstate0 + here->HSM2vgs);
	delvges = vges - *(ckt->CKTstate0 + here->HSM2vges);
	delvds = vds - *(ckt->CKTstate0 + here->HSM2vds);
	delvdbs = vdbs - *(ckt->CKTstate0 + here->HSM2vdbs);
	delvsbs = vsbs - *(ckt->CKTstate0 + here->HSM2vsbs);
	delvdbd = vdbd - *(ckt->CKTstate0 + here->HSM2vdbd);
	delvgd = vgd - vgdo;
/*	delvged = vged - vgedo;*/

	delvbd_jct = (!here->HSM2_corbnet) ? delvbd : delvdbd;
	delvbs_jct = (!here->HSM2_corbnet) ? delvbs : delvsbs;
	if (here->HSM2_mode >= 0) {
	  Idtot = here->HSM2_ids + here->HSM2_isub - here->HSM2_ibd
	    + here->HSM2_igidl;
	  cdhat = Idtot - here->HSM2_gbd * delvbd_jct
	    + (here->HSM2_gmbs + here->HSM2_gbbs + here->HSM2_gigidlbs) * delvbs 
	    + (here->HSM2_gm + here->HSM2_gbgs + here->HSM2_gigidlgs) * delvgs
	    + (here->HSM2_gds + here->HSM2_gbds + here->HSM2_gigidlds) * delvds;
	  Ibtot = here->HSM2_ibs + here->HSM2_ibd - here->HSM2_isub
	    - here->HSM2_igidl - here->HSM2_igisl;
	  cbhat = Ibtot + here->HSM2_gbd * delvbd_jct
	    + here->HSM2_gbs * delvbs_jct - (here->HSM2_gbbs + here->HSM2_gigidlbs) * delvbs
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
	  Idtot = here->HSM2_ids + here->HSM2_ibd - here->HSM2_igidl;
	  cdhat = Idtot + here->HSM2_gbd * delvbd_jct + here->HSM2_gmbs * delvbd
	    + here->HSM2_gm * delvgd - here->HSM2_gds * delvds
	    - here->HSM2_gigidlgs * delvgd - here->HSM2_gigidlbs * delvbd 
	    + here->HSM2_gigidlds * delvds ;
	  Ibtot = here->HSM2_ibs + here->HSM2_ibd - here->HSM2_isub
	    - here->HSM2_igidl - here->HSM2_igisl;
	  cbhat = Ibtot + here->HSM2_gbs * delvbs_jct
	    + here->HSM2_gbd * delvbd_jct - (here->HSM2_gbbs + here->HSM2_gigidlbs) * delvbd
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

	vds_pre = vds;

#ifndef NOBYPASS /* BSIM3 style */
	/* now lets see if we can bypass (ugh) */
	
	/* following should be one big if connected by && all over
	 * the place, but some C compilers can't handle that, so
	 * we split it up here to let them digest it in stages
	 */
	if ( !(ckt->CKTmode & MODEINITPRED) && BYPASS_enable )
	  if ((!here->HSM2_corbnet) || 
	      (fabs(delvdbs) < 
	       (reltol
		* MAX(fabs(vdbs), fabs(*(ckt->CKTstate0 + here->HSM2vdbs)))
		+ voltTol)))
	    if ((!here->HSM2_corbnet) || 
		(fabs(delvdbd) < 
		 (reltol
		  * MAX(fabs(vdbd), fabs(*(ckt->CKTstate0 + here->HSM2vdbd)))
		  + voltTol)))
               if ((!here->HSM2_corbnet) || 
		   (fabs(delvsbs) < 
		    (reltol
		     * MAX(fabs(vsbs), fabs(*(ckt->CKTstate0 + here->HSM2vsbs)))
		     + voltTol)))
		 if ((here->HSM2_corg == 0) || (here->HSM2_corg == 1) || 
		     (fabs(delvges) < 
		      (reltol 
		       * MAX(fabs(vges), fabs(*(ckt->CKTstate0 + here->HSM2vges)))
		       + voltTol)))
	  if ( fabs(delvbs) < 
	       ( reltol * 
		 MAX(fabs(vbs), fabs(*(ckt->CKTstate0+here->HSM2vbs))) + 
		 voltTol ) )
	    if ( fabs(delvbd) < 
		 ( reltol * 
		   MAX(fabs(vbd), fabs(*(ckt->CKTstate0+here->HSM2vbd))) + 
		   voltTol ) )
	      if ( fabs(delvgs) < 
		   ( reltol * 
		     MAX(fabs(vgs), fabs(*(ckt->CKTstate0+here->HSM2vgs))) +
		     voltTol ) )
		if ( fabs(delvds) < 
		     ( reltol * 
		       MAX(fabs(vds), fabs(*(ckt->CKTstate0+here->HSM2vds))) + 
		       voltTol ) )
		  if ( fabs(cdhat - Idtot) < 
		       ( reltol * 
			 MAX(fabs(cdhat),fabs(Idtot)) + abstol ) ) 
		    if (!model->HSM2_coiigs || 
			(fabs(cgbhat - Igbtot) < reltol
			 * MAX(fabs(cgbhat), fabs(Igbtot)) + abstol))
		      if (!model->HSM2_coiigs || 
			  (fabs(cgshat - Igstot) < reltol
			   * MAX(fabs(cgshat), fabs(Igstot)) + abstol))
			if (!model->HSM2_coiigs || 
			    (fabs(cgdhat - Igdtot) < reltol
			     * MAX(fabs(cgdhat), fabs(Igdtot)) + abstol)){
			  tempv = MAX(fabs(cbhat),fabs(Ibtot)) + abstol;
			  if ((fabs(cbhat - Ibtot)) < reltol * tempv) {
			    /* bypass code */
			    vbs = *(ckt->CKTstate0 + here->HSM2vbs);
			    vbd = *(ckt->CKTstate0 + here->HSM2vbd);
			    vgs = *(ckt->CKTstate0 + here->HSM2vgs);
			    vds = *(ckt->CKTstate0 + here->HSM2vds);

			    vges = *(ckt->CKTstate0 + here->HSM2vges);
			    vdbs = *(ckt->CKTstate0 + here->HSM2vdbs);
			    vdbd = *(ckt->CKTstate0 + here->HSM2vdbd);
			    vsbs = *(ckt->CKTstate0 + here->HSM2vsbs);

			    vgd = vgs - vds;
			    vgb = vgs - vbs;
			    vged = vges - vds;

			    vbs_jct = (!here->HSM2_corbnet) ? vbs : vsbs;
			    vbd_jct = (!here->HSM2_corbnet) ? vbd : vdbd;

			    cdrain = here->HSM2_ids;
			    if ((ckt->CKTmode & (MODETRAN | MODEAC)) || 
				((ckt->CKTmode & MODETRANOP) && 
				 (ckt->CKTmode & MODEUIC))) {
			      ByPass = 1;
			      qgate = here->HSM2_qg;
			      qbulk = here->HSM2_qb;
			      qdrn = here->HSM2_qd;
			      goto line755;
			    }
			    else
			      goto line850;
			  }
			}
#endif /*NOBYPASS*/

#ifdef DEBUG_HISIM2LD_VX	
	printf( "vbd_p    = %12.5e\n" , vbd );
	printf( "vbs_p    = %12.5e\n" , vbs );
	printf( "vgs_p    = %12.5e\n" , vgs );
	printf( "vds_p    = %12.5e\n" , vds );
#endif

	von = here->HSM2_von; 
	if(*(ckt->CKTstate0 + here->HSM2vds) >= 0.0) {
	  vgs = DEVfetlim(vgs, *(ckt->CKTstate0 + here->HSM2vgs), von);
	  vds = vgs - vgd;
	  vds = DEVlimvds(vds, *(ckt->CKTstate0 + here->HSM2vds));
	  vgd = vgs - vds;

	  if (here->HSM2_corg == 1) {
	    vges = DEVfetlim(vges, *(ckt->CKTstate0 + here->HSM2vges), von);
	    vged = vges - vds;
	  }
	} 
	else {
	  vgd = DEVfetlim(vgd, vgdo, von);
	  vds = vgs - vgd;
	  vds = -DEVlimvds(-vds, -(*(ckt->CKTstate0 + here->HSM2vds)));
	  vgs = vgd + vds;

	  if (here->HSM2_corg == 1) {
	    vged = DEVfetlim(vged, vgedo, von);
	    vges = vged + vds;
	  }
	}
	if (vds >= 0.0) {
	  vbs = DEVpnjlim(vbs, *(ckt->CKTstate0 + here->HSM2vbs),
			  CONSTvt0, model->HSM2_vcrit, &Check);
	  vbd = vbs - vds;
	  if (here->HSM2_corbnet) {
	    vdbs = DEVpnjlim(vdbs, *(ckt->CKTstate0 + here->HSM2vdbs),
			     CONSTvt0, model->HSM2_vcrit, &Check1);
	    vdbd = vdbs - vds;
	    vsbs = DEVpnjlim(vsbs, *(ckt->CKTstate0 + here->HSM2vsbs),
			     CONSTvt0, model->HSM2_vcrit, &Check2);
	    if ((Check1 == 0) && (Check2 == 0)) Check = 0;
	    else Check = 1;
	  }
	} 
	else {
	  vbd = DEVpnjlim(vbd, *(ckt->CKTstate0 + here->HSM2vbd),
			  CONSTvt0, model->HSM2_vcrit, &Check);
	  vbs = vbd + vds;
	  if (here->HSM2_corbnet) {
	    vdbd = DEVpnjlim(vdbd, *(ckt->CKTstate0 + here->HSM2vdbd),
			     CONSTvt0, model->HSM2_vcrit, &Check1);
	    vdbs = vdbd + vds;
	    vsbdo = *(ckt->CKTstate0 + here->HSM2vsbs)
	      - *(ckt->CKTstate0 + here->HSM2vds);
	    vsbd = vsbs - vds;
	    vsbd = DEVpnjlim(vsbd, vsbdo, CONSTvt0, model->HSM2_vcrit, &Check2);
	    vsbs = vsbd + vds;
	    if ((Check1 == 0) && (Check2 == 0)) Check = 0;
	    else Check = 1;
	  }
	}
      }
      
      vbd = vbs - vds;
      vgd = vgs - vds;
      vgb = vgs - vbs;
      vged = vges - vds;
      vdbd = vdbs - vds;

      vbs_jct = (!here->HSM2_corbnet) ? vbs : vsbs;
      vbd_jct = (!here->HSM2_corbnet) ? vbd : vdbd;

#ifdef DEBUG_HISIM2LD_VX
      printf( "vbd    = %12.5e\n" , vbd );
      printf( "vbs    = %12.5e\n" , vbs );
      printf( "vgs    = %12.5e\n" , vgs );
      printf( "vds    = %12.5e\n" , vds );
#endif

      if (vds >= 0) { /* normal mode */
	here->HSM2_mode = 1;
	ivds = vds;
	ivgs = vgs;
	ivbs = vbs;
      } else { /* reverse mode */
	here->HSM2_mode = -1;
	ivds = -vds;
	ivgs = vgd;
	ivbs = vbd;
      }

      if ( model->HSM2_info >= 5 ) { /* mode, bias conditions ... */
	printf( "--- variables given to HSM2evaluate() ----\n" );
	printf( "type   = %s\n" , (model->HSM2_type>0) ? "NMOS" : "PMOS" );
	printf( "mode   = %s\n" , (here->HSM2_mode>0) ? "NORMAL" : "REVERSE" );
	
	printf( "vbs    = %12.5e\n" , ivbs );
	printf( "vds    = %12.5e\n" , ivds );
	printf( "vgs    = %12.5e\n" , ivgs );
      }
      if ( model->HSM2_info >= 6 ) { /* input flags */
	printf( "corsrd = %s\n" , (model->HSM2_corsrd) ? "true" : "false" ) ;
	printf( "coadov = %s\n" , (model->HSM2_coadov) ? "true" : "false" ) ;
	printf( "coisub = %s\n" , (model->HSM2_coisub) ? "true" : "false" ) ;
	printf( "coiigs = %s\n" , (model->HSM2_coiigs) ? "true" : "false" ) ;
	printf( "cogidl = %s\n" , (model->HSM2_cogidl) ? "true" : "false" ) ;
	printf( "coovlp = %s\n" , (model->HSM2_coovlp) ? "true" : "false" ) ;
	printf( "coflick = %s\n" , (model->HSM2_coflick) ? "true" : "false" ) ;
	printf( "coisti = %s\n" , (model->HSM2_coisti) ? "true" : "false" ) ;
	printf( "conqs  = %s\n" , (model->HSM2_conqs)  ? "true" : "false" ) ;
	printf( "cothrml = %s\n" , (model->HSM2_cothrml) ? "true" : "false" ) ;
	printf( "coign = %s\n" , (model->HSM2_coign) ? "true" : "false" ) ;
      }
      /* print inputs ------------AA */

#ifdef DEBUG_HISIM2CGG
      /* Print convergence flag */
      printf("isConv %d ", isConv );
      printf("CKTtime %e ", ckt->CKTtime );
      printf("Vb %1.3e ", (model->HSM2_type>0) ? vbs:-vbs );
      printf("Vd %1.3e ", (model->HSM2_type>0) ? vds:-vds );
      printf("Vg %1.3e ", (model->HSM2_type>0) ? vgs:-vgs );
#endif

      /* call model evaluation */
      if ( HSM2evaluate(ivds, ivgs, ivbs, vbs_jct, vbd_jct, here, model, ckt) == HiSIM_ERROR ) 
	return (HiSIM_ERROR);


#ifdef DEBUG_HISIM2CGG
      printf("HSM2_ids %e ", here->HSM2_ids ) ;
      printf("HSM2_cggb %e ", here->HSM2_cggb ) ;
      printf("\n") ;
#endif

      /* modified by T.Y. 2006.05.31 
      * if ( !here->HSM2_called ) here->HSM2_called = 1;
      */
      here->HSM2_called += 1;

      cdrain = here->HSM2_ids ; /* cdrain */
      qgate = here->HSM2_qg ; /* gate */
      qdrn = here->HSM2_qd ;  /* drain */
      qbulk = here->HSM2_qb = -1.0 * (here->HSM2_qg + here->HSM2_qd + here->HSM2_qs); /* bulk */

      /* print all outputs ------------VV */
      if ( model->HSM2_info >= 4 ) {
	printf( "--- variables returned from HSM2evaluate() ----\n" ) ;
	
	printf( "von    = %12.5e\n" , here->HSM2_von ) ;
	printf( "vdsat  = %12.5e\n" , here->HSM2_vdsat ) ;
	printf( "ids    = %12.5e\n" , here->HSM2_ids ) ;
	
	printf( "gds    = %12.5e\n" , here->HSM2_gds ) ;
	printf( "gm     = %12.5e\n" , here->HSM2_gm ) ;
	printf( "gmbs   = %12.5e\n" , here->HSM2_gmbs ) ;

	printf( "cggo   = %12.5e\n" , -(here->HSM2_cgdo + here->HSM2_cgso +here->HSM2_cgbo) ) ;	
	printf( "cgdo   = %12.5e\n" , here->HSM2_cgdo ) ;
	printf( "cgso   = %12.5e\n" , here->HSM2_cgso ) ;
	printf( "cdgo   = %12.5e\n" , here->HSM2_cdgo ) ;
	printf( "cddo   = %12.5e\n" , here->HSM2_cddo ) ;
	printf( "cdso   = %12.5e\n" , here->HSM2_cdso ) ;
	printf( "csgo   = %12.5e\n" , here->HSM2_csgo ) ;
	printf( "csdo   = %12.5e\n" , here->HSM2_csdo ) ;
	printf( "csso   = %12.5e\n" , here->HSM2_csso ) ;
	
	printf( "qg     = %12.5e\n" , here->HSM2_qg ) ;
	printf( "qd     = %12.5e\n" , here->HSM2_qd ) ;
	printf( "qs     = %12.5e\n" , here->HSM2_qs ) ;
	
	printf( "cggb   = %12.5e\n" , here->HSM2_cggb ) ;
	printf( "cgsb   = %12.5e\n" , here->HSM2_cgsb ) ;
	printf( "cgdb   = %12.5e\n" , here->HSM2_cgdb ) ;
	printf( "cbgb   = %12.5e\n" , here->HSM2_cbgb ) ;
	printf( "cbsb   = %12.5e\n" , here->HSM2_cbsb ) ;
	printf( "cbdb   = %12.5e\n" , here->HSM2_cbdb ) ;
	printf( "cdgb   = %12.5e\n" , here->HSM2_cdgb ) ;
	printf( "cdsb   = %12.5e\n" , here->HSM2_cdsb ) ;
	printf( "cddb   = %12.5e\n" , here->HSM2_cddb ) ;
      
	printf( "ibd    = %12.5e\n" , here->HSM2_ibd ) ;
	printf( "ibs    = %12.5e\n" , here->HSM2_ibs ) ;
	printf( "gbd    = %12.5e\n" , here->HSM2_gbd ) ;
	printf( "gbs    = %12.5e\n" , here->HSM2_gbs ) ;
	printf( "capbd  = %12.5e\n" , here->HSM2_capbd ) ;
	printf( "capbs  = %12.5e\n" , here->HSM2_capbs ) ;
	printf( "qbd    = %12.5e\n" , *(ckt->CKTstate0 + here->HSM2qbd) ) ;
	printf( "qbs    = %12.5e\n" , *(ckt->CKTstate0 + here->HSM2qbs) ) ;

	printf( "isub   = %12.5e\n" , here->HSM2_isub ) ;
	printf( "gbgs   = %12.5e\n" , here->HSM2_gbgs ) ;
	printf( "gbds   = %12.5e\n" , here->HSM2_gbds ) ;
	printf( "gbbs   = %12.5e\n" , here->HSM2_gbbs ) ;

 	printf( "S_flicker_noise * ( freq / gain ) = %.16e\n" , here->HSM2_noiflick ) ;
 	printf( "S_thermal_noise / ( gain * 4kT )  = %.16e\n" , here->HSM2_noithrml ) ;
 	printf( "S_induced_gate_noise / ( gain * freq^2 ) = %.16e\n" , here->HSM2_noiigate ) ;
 	printf( "cross-correlation coefficient (= Sigid/sqrt(Sig*Sid) ) = %.16e\n" , here->HSM2_noicross ) ;
	/* print Surface Potentials */
	printf( "ivds %e ivgs %e ivbs %e Ps0 %.16e Pds %.16e\n" ,
		ivds, ivgs, ivbs, here->HSM2_ps0_prv, here->HSM2_pds_prv ) ;
      }
      /* print all outputs ------------AA */

      if ( model->HSM2_info >= 3 ) { /* physical valiables vs bias */
	static int isFirst = 1;
	if (isFirst) {
	  printf("# vbs vds vgs cggb cgdb cgsb cbgb cbdb cbsb cdgb cddb cdsb\n");
#ifndef USE_OMP
	  isFirst = 0;
#endif
	}
	printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", 
	       vbs, vds, vgs , 
	       here->HSM2_cggb, here->HSM2_cgdb, here->HSM2_cgsb, 
	       here->HSM2_cbgb, here->HSM2_cbdb, here->HSM2_cbsb,
	       here->HSM2_cdgb, here->HSM2_cddb, here->HSM2_cdsb);
	  
      }
      /*
       *  check convergence
       */
      isConv = 1;
      if ( (here->HSM2_off == 0) || !(ckt->CKTmode & MODEINITFIX) ) {
	if (Check == 1) {
	  ckt->CKTnoncon++;
	  isConv = 0;
#ifndef NEWCONV
	} 
	else {
	  if (here->HSM2_mode >= 0) 
	    Idtot = here->HSM2_ids + here->HSM2_isub - here->HSM2_ibd + here->HSM2_igidl;
	  else
	    Idtot = here->HSM2_ids + here->HSM2_ibd - here->HSM2_igidl;
	  tol = ckt->CKTreltol * MAX(fabs(cdhat), fabs(Idtot)) + ckt->CKTabstol;
	  tol2 = ckt->CKTreltol * MAX(fabs(cgbhat), fabs(Igbtot)) + ckt->CKTabstol;
	  tol3 = ckt->CKTreltol * MAX(fabs(cgshat), fabs(Igstot)) + ckt->CKTabstol;
	  tol4 = ckt->CKTreltol * MAX(fabs(cgdhat), fabs(Igdtot)) + ckt->CKTabstol;
	  if (fabs(cdhat - Idtot) >= tol) {
	    ckt->CKTnoncon++;
	    isConv = 0;
	  }
	  else if (fabs(cgbhat - Igbtot) >= tol2 || 
		   fabs(cgshat - Igstot) >= tol3 ||
		   fabs(cgdhat - Igdtot) >= tol4) {
	    ckt->CKTnoncon++;
	    isConv = 0;
	  }
	  else {
	    Ibtot = here->HSM2_ibs + here->HSM2_ibd 
	      - here->HSM2_isub - here->HSM2_igidl - here->HSM2_igisl;
	    tol = ckt->CKTreltol * MAX(fabs(cbhat), fabs(Ibtot)) + ckt->CKTabstol;
	    if (fabs(cbhat - Ibtot) > tol) {
	      ckt->CKTnoncon++;
	      isConv = 0;
	    }
	  }
	}
#endif /* NEWCONV */
      }
    }
    *(ckt->CKTstate0 + here->HSM2vbs) = vbs;
    *(ckt->CKTstate0 + here->HSM2vbd) = vbd;
    *(ckt->CKTstate0 + here->HSM2vgs) = vgs;
    *(ckt->CKTstate0 + here->HSM2vds) = vds;
    *(ckt->CKTstate0 + here->HSM2vsbs) = vsbs;
    *(ckt->CKTstate0 + here->HSM2vdbs) = vdbs;
    *(ckt->CKTstate0 + here->HSM2vdbd) = vdbd;
    *(ckt->CKTstate0 + here->HSM2vges) = vges;

    if ((ckt->CKTmode & MODEDC) && 
	!(ckt->CKTmode & MODEINITFIX) && !(ckt->CKTmode & MODEINITJCT)) 
      showPhysVal = 1;
    if (model->HSM2_show_Given && showPhysVal && isConv) {
      static int isFirst = 1;
      if (vds != vds_pre) 
	ShowPhysVals(here, model, isFirst, vds_pre, vgs, vbs, vgd, vbd, vgb);
      else 
	ShowPhysVals(here, model, isFirst, vds, vgs, vbs, vgd, vbd, vgb);
#ifndef USE_OMP
      if (isFirst) isFirst = 0;
#endif
    }

    /* bulk and channel charge plus overlaps */
    
    if (!ChargeComputationNeeded) goto line850; 
    
 line755:

    ag0 = ckt->CKTag[0];
    if (here->HSM2_mode > 0) { /* NORMAL mode */
      gcggb = here->HSM2_cggb * ag0;
      gcgdb = here->HSM2_cgdb * ag0;
      gcgsb = here->HSM2_cgsb * ag0;
      gcgbb = -(gcggb + gcgdb + gcgsb);

      gcdgb = here->HSM2_cdgb * ag0;
      gcddb = (here->HSM2_cddb + here->HSM2_capbd) * ag0;
      gcdsb = here->HSM2_cdsb * ag0;

      gcsgb = -(here->HSM2_cggb + here->HSM2_cbgb + here->HSM2_cdgb) * ag0;
      gcsdb = -(here->HSM2_cgdb + here->HSM2_cbdb + here->HSM2_cddb) * ag0;
      gcssb = (here->HSM2_capbs 
	       - (here->HSM2_cgsb + here->HSM2_cbsb + here->HSM2_cdsb)) * ag0;

      gcbgb = here->HSM2_cbgb * ag0;

      if ( !here->HSM2_corbnet ) {
	gcdbb = -(gcdgb + gcddb + gcdsb);
	gcsbb = -(gcsgb + gcsdb + gcssb);
	gcbdb = (here->HSM2_cbdb - here->HSM2_capbd) * ag0;
	gcbsb = (here->HSM2_cbsb - here->HSM2_capbs) * ag0;
	gcdbdb = 0.0; gcsbsb = 0.0;
      } else {
	gcdbb = -(gcdgb + gcddb + gcdsb) + here->HSM2_capbd * ag0;
	gcsbb = -(gcsgb + gcsdb + gcssb) + here->HSM2_capbs * ag0;
	gcbdb = here->HSM2_cbdb * ag0;
	gcbsb = here->HSM2_cbsb * ag0;
	gcdbdb = - here->HSM2_capbd * ag0;
	gcsbsb = - here->HSM2_capbs * ag0;
      }
      gcbbb = -(gcbdb + gcbgb + gcbsb);

    }
    else { /* REVERSE mode */
      gcggb = here->HSM2_cggb * ag0;
      gcgdb = here->HSM2_cgsb * ag0;
      gcgsb = here->HSM2_cgdb * ag0;
      gcgbb = -(gcggb + gcgdb + gcgsb);
      
      gcdgb = -(here->HSM2_cggb + here->HSM2_cbgb + here->HSM2_cdgb) * ag0;
      gcddb = (here->HSM2_capbd 
	       - (here->HSM2_cgsb + here->HSM2_cbsb + here->HSM2_cdsb)) * ag0;
      gcdsb = -(here->HSM2_cgdb + here->HSM2_cbdb + here->HSM2_cddb) * ag0;
      
      gcsgb = here->HSM2_cdgb * ag0;
      gcsdb = here->HSM2_cdsb * ag0;
      gcssb = (here->HSM2_cddb + here->HSM2_capbs) * ag0;
      
      gcbgb = here->HSM2_cbgb * ag0;

      if ( !here->HSM2_corbnet ){
	gcdbb = -(gcdgb + gcddb + gcdsb);
	gcsbb = -(gcsgb + gcsdb + gcssb);
	gcbdb = (here->HSM2_cbsb - here->HSM2_capbd) * ag0;
	gcbsb = (here->HSM2_cbdb - here->HSM2_capbs) * ag0;
	gcdbdb = 0.0; gcsbsb = 0.0;
      } else {
	gcdbb = -(gcdgb + gcddb + gcdsb) + here->HSM2_capbd * ag0;
	gcsbb = -(gcsgb + gcsdb + gcssb) + here->HSM2_capbs * ag0;
	gcbdb = here->HSM2_cbsb * ag0;
	gcbsb = here->HSM2_cbdb * ag0;
	gcdbdb = - here->HSM2_capbd * ag0;
	gcsbsb = - here->HSM2_capbs * ag0;
      }
      gcbbb = -(gcbgb + gcbdb + gcbsb);

      qdrn = -(qgate + qbulk + qdrn);
    }

    if (ByPass) goto line860;
    
    *(ckt->CKTstate0 + here->HSM2qg) = qgate;
    *(ckt->CKTstate0 + here->HSM2qd) = qdrn - *(ckt->CKTstate0 + here->HSM2qbd);
    if ( !here->HSM2_corbnet ) {
      *(ckt->CKTstate0 + here->HSM2qb) = qbulk 
	+ *(ckt->CKTstate0 + here->HSM2qbd) + *(ckt->CKTstate0 + here->HSM2qbs);
    } else {
      *(ckt->CKTstate0 + here->HSM2qb) = qbulk;
    }

#ifdef DEBUG_HISIM2LD
    printf( "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq\n" ) ;
    printf( "HSM2qg   = %12.5e\n" , *(ckt->CKTstate0 + here->HSM2qg) ) ;    
    printf( "HSM2qd   = %12.5e\n" , *(ckt->CKTstate0 + here->HSM2qd) ) ;    
    printf( "HSM2qb   = %12.5e\n" , *(ckt->CKTstate0 + here->HSM2qb) ) ;    
    printf( "qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq\n" ) ;
#endif

    /* store small signal parameters */
    if (ckt->CKTmode & MODEINITSMSIG) goto line1000;
    if (!ChargeComputationNeeded) goto line850;

    if (ckt->CKTmode & MODEINITTRAN) {
      *(ckt->CKTstate1 + here->HSM2qb) = *(ckt->CKTstate0 + here->HSM2qb);
      *(ckt->CKTstate1 + here->HSM2qg) = *(ckt->CKTstate0 + here->HSM2qg);
      *(ckt->CKTstate1 + here->HSM2qd) = *(ckt->CKTstate0 + here->HSM2qd);
      if ( here->HSM2_corbnet ) {
	*(ckt->CKTstate1 + here->HSM2qbs) = *(ckt->CKTstate0 + here->HSM2qbs);
	*(ckt->CKTstate1 + here->HSM2qbd) = *(ckt->CKTstate0 + here->HSM2qbd);
      }
    }
    
    return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSM2qb));
    return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSM2qg));
    return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSM2qd));
    if ( here->HSM2_corbnet ) {
      return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSM2qbs));
      return_if_error (NIintegrate(ckt, &geq, &ceq, 0.0, here->HSM2qbd));
    }

    goto line860;
    
 line850:
    /* initialize to zero charge conductance and current */
    ceqqg = ceqqb = ceqqd = 0.0;
    ceqqjd = ceqqjs = 0.0;
	  
    gcdgb = gcddb = gcdsb = gcdbb = 0.0;
    gcsgb = gcsdb = gcssb = gcsbb = 0.0;
    gcggb = gcgdb = gcgsb = gcgbb = 0.0;
    gcbgb = gcbdb = gcbsb = gcbbb = 0.0;

    gcdbdb = gcsbsb = 0.0;

    goto line900;
    
 line860:
    /* evaluate equivalent charge current */
    
    cqgate = *(ckt->CKTstate0 + here->HSM2cqg);
    cqbulk = *(ckt->CKTstate0 + here->HSM2cqb);
    cqdrn = *(ckt->CKTstate0 + here->HSM2cqd);

#ifdef DEBUG_HISIM2LD
    printf( "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii\n" ) ;
    printf( "cqgate   = %12.5e\n" , cqgate ) ;
    printf( "cqbulk   = %12.5e\n" , cqbulk ) ;
    printf( "cqdrn   = %12.5e\n" , cqdrn ) ;
    printf( "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii\n" ) ;
#endif
    
    ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
    ceqqd = cqdrn - gcdgb * vgb + (gcddb + gcdbdb) * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
    ceqqb = cqbulk - gcbgb * vgb + gcbdb * vbd + gcbsb * vbs;

    if (here->HSM2_corbnet) {
      ceqqjs = *(ckt->CKTstate0 + here->HSM2cqbs) + gcsbsb * vbs_jct;
      ceqqjd = *(ckt->CKTstate0 + here->HSM2cqbd) + gcdbdb * vbd_jct; 
    }

    if (ckt->CKTmode & MODEINITTRAN) {
      *(ckt->CKTstate1 + here->HSM2cqb) = *(ckt->CKTstate0 + here->HSM2cqb);
      *(ckt->CKTstate1 + here->HSM2cqg) = *(ckt->CKTstate0 + here->HSM2cqg);
      *(ckt->CKTstate1 + here->HSM2cqd) = *(ckt->CKTstate0 + here->HSM2cqd);
      if (here->HSM2_corbnet) {
	*(ckt->CKTstate1 + here->HSM2cqbs) = *(ckt->CKTstate0 + here->HSM2cqbs);
	*(ckt->CKTstate1 + here->HSM2cqbd) = *(ckt->CKTstate0 + here->HSM2cqbd);
      }      
    }

    /*
     *  load current vector
     */
 line900:

    if (here->HSM2_mode >= 0) { /* NORMAL mode */
      gm = here->HSM2_gm;
      gmbs = here->HSM2_gmbs;
      FwdSum = gm + gmbs;
      RevSum = 0.0;

      cdreq = model->HSM2_type * 
	(cdrain - here->HSM2_gds * vds - gm * vgs - gmbs * vbs); 
      ceqbd = model->HSM2_type * (here->HSM2_isub + here->HSM2_igidl 
				  - (here->HSM2_gbds + here->HSM2_gigidlds) * vds 
				  - (here->HSM2_gbgs + here->HSM2_gigidlgs) * vgs 
				  - (here->HSM2_gbbs + here->HSM2_gigidlbs) * vbs);
      ceqbs = model->HSM2_type * (here->HSM2_igisl 
				  + here->HSM2_gigislsd * vds 
				  - here->HSM2_gigislgd * vgd 
				  - here->HSM2_gigislbd * vbd);
      
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
	Ibtoteq = model->HSM2_type * 
	  (here->HSM2_igb - here->HSM2_gigbg * vgs 
	   - here->HSM2_gigbd * vds - here->HSM2_gigbb * vbs);

	gIstotg = here->HSM2_gigsg;
	gIstotd = here->HSM2_gigsd;
	gIstots = here->HSM2_gigss;
	gIstotb = here->HSM2_gigsb;
	Istoteq = model->HSM2_type * 
	  (here->HSM2_igs - here->HSM2_gigsg * vgs
	   - here->HSM2_gigsd * vds - here->HSM2_gigsb * vbs);

	gIdtotg = here->HSM2_gigdg;
	gIdtotd = here->HSM2_gigdd;
	gIdtots = here->HSM2_gigds;
	gIdtotb = here->HSM2_gigdb;
	Idtoteq = model->HSM2_type * 
	  (here->HSM2_igd - here->HSM2_gigdg * vgs
	   - here->HSM2_gigdd * vds - here->HSM2_gigdb * vbs);
      }
      else {
	gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;
	gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
	gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
      }

      if (model->HSM2_coiigs) {
	gIgtotg = gIbtotg + gIstotg + gIdtotg;
	gIgtotd = gIbtotd + gIstotd + gIdtotd;
	gIgtots = gIbtots + gIstots + gIdtots;
	gIgtotb = gIbtotb + gIstotb + gIdtotb;
	Igtoteq = Ibtoteq + Istoteq + Idtoteq;
      }
      else
	gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;

    }
    else { /* REVERSE mode */
      gm = - here->HSM2_gm;
      gmbs = - here->HSM2_gmbs;
      FwdSum = 0.0;
      RevSum = -(gm + gmbs);
      cdreq = -model->HSM2_type * (cdrain + here->HSM2_gds * vds + gm * vgd + gmbs * vbd);

      ceqbs = model->HSM2_type * (here->HSM2_isub + here->HSM2_igisl 
				  + (here->HSM2_gbds + here->HSM2_gigislsd) * vds
				  - (here->HSM2_gbgs + here->HSM2_gigislgd) * vgd
				  - (here->HSM2_gbbs + here->HSM2_gigislbd) * vbd);
      ceqbd = model->HSM2_type * (here->HSM2_igidl 
				  - here->HSM2_gigidlds * vds 
				  - here->HSM2_gigidlgs * vgs 
				  - here->HSM2_gigidlbs * vbs);

      gbbsp = - here->HSM2_gbds;
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
	Ibtoteq = model->HSM2_type * 
	  (here->HSM2_igb - here->HSM2_gigbg * vgd 
	   + here->HSM2_gigbs * vds - here->HSM2_gigbb * vbd);

	gIstotg = here->HSM2_gigsg;
	gIstotd = here->HSM2_gigsd;
	gIstots = here->HSM2_gigss;
	gIstotb = here->HSM2_gigsb;
	Istoteq = model->HSM2_type * 
	  (here->HSM2_igs - here->HSM2_gigsg * vgd
	   + here->HSM2_gigss * vds - here->HSM2_gigsb * vbd);

	gIdtotg = here->HSM2_gigdg;
	gIdtotd = here->HSM2_gigdd;
	gIdtots = here->HSM2_gigds;
	gIdtotb = here->HSM2_gigdb;
	Idtoteq = model->HSM2_type * 
	  (here->HSM2_igd - here->HSM2_gigdg * vgd
	   + here->HSM2_gigds * vds - here->HSM2_gigdb * vbd);
      }
      else {
	gIbtotg = gIbtotd = gIbtots = gIbtotb = Ibtoteq = 0.0;
	gIstotg = gIstotd = gIstots = gIstotb = Istoteq = 0.0;
	gIdtotg = gIdtotd = gIdtots = gIdtotb = Idtoteq = 0.0;
      }

      if (model->HSM2_coiigs) {
	gIgtotg = gIbtotg + gIstotg + gIdtotg;
	gIgtotd = gIbtotd + gIstotd + gIdtotd;
	gIgtots = gIbtots + gIstots + gIdtots;
	gIgtotb = gIbtotb + gIstotb + gIdtotb;
	Igtoteq = Ibtoteq + Istoteq + Idtoteq;
      }
      else
	gIgtotg = gIgtotd = gIgtots = gIgtotb = Igtoteq = 0.0;

    }
    
    if (model->HSM2_type > 0) { 
      ceqjs = here->HSM2_ibs - here->HSM2_gbs * vbs_jct;
      ceqjd = here->HSM2_ibd - here->HSM2_gbd * vbd_jct;
    }
    else {
      ceqjs = -(here->HSM2_ibs - here->HSM2_gbs * vbs_jct); 
      ceqjd = -(here->HSM2_ibd - here->HSM2_gbd * vbd_jct);
      ceqqg = -ceqqg;
      ceqqb = -ceqqb;
      ceqqd = -ceqqd;

      if (here->HSM2_corbnet) {
	ceqqjs = -ceqqjs;
	ceqqjd = -ceqqjd;
      }
    }

#ifdef DEBUG_HISIM2LD
    printf( "----------------------------------------------------\n" ) ;
    printf( "ceqqg   = %12.5e\n" , ceqqg ) ;
    printf( "....................................................\n" ) ;
    printf( "ceqbs   = %12.5e\n" , ceqbs ) ;
    printf( "ceqbd   = %12.5e\n" , ceqbd ) ;
    printf( "ceqqb   = %12.5e\n" , ceqqb ) ;
    printf( "....................................................\n" ) ;
    printf( "ceqbd   = %12.5e\n" , ceqbd ) ;
    printf( "cdreq   = %12.5e\n" , cdreq ) ;
    printf( "ceqqd   = %12.5e\n" , ceqqd ) ;
    printf( "----------------------------------------------------\n" ) ;
#endif

#ifdef USE_OMP
    here->HSM2rhsdPrime = ceqjd - ceqbd - cdreq - ceqqd + Idtoteq;
    here->HSM2rhsgPrime = ceqqg + Igtoteq;

    if ( !here->HSM2_corbnet ) {
      here->HSM2rhsbPrime = ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq;
      here->HSM2rhssPrime = cdreq - ceqbs + ceqjs + ceqqg + ceqqb + ceqqd + Istoteq;
    } else {
      here->HSM2rhsdb = ceqjd + ceqqjd;
      here->HSM2rhsbPrime = ceqbd + ceqbs - ceqqb + Ibtoteq;
      here->HSM2rhssb = ceqjs + ceqqjs;
      here->HSM2rhssPrime = cdreq - ceqbs + ceqjs + ceqqd
                                + ceqqg + ceqqb + ceqqjd + ceqqjs + Istoteq;
    }
#else
    *(ckt->CKTrhs + here->HSM2dNodePrime) += ceqjd - ceqbd - cdreq - ceqqd + Idtoteq;
    *(ckt->CKTrhs + here->HSM2gNodePrime) -= ceqqg + Igtoteq;

    if ( !here->HSM2_corbnet ) {
      *(ckt->CKTrhs + here->HSM2bNodePrime) += ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq;
      *(ckt->CKTrhs + here->HSM2sNodePrime) += cdreq - ceqbs + ceqjs + ceqqg + ceqqb + ceqqd + Istoteq;
    } else {
      *(ckt->CKTrhs + here->HSM2dbNode) -= ceqjd + ceqqjd;
      *(ckt->CKTrhs + here->HSM2bNodePrime) += ceqbd + ceqbs - ceqqb + Ibtoteq;
      *(ckt->CKTrhs + here->HSM2sbNode) -= ceqjs + ceqqjs;
      *(ckt->CKTrhs + here->HSM2sNodePrime) += cdreq - ceqbs + ceqjs + ceqqd
	+ ceqqg + ceqqb + ceqqjd + ceqqjs + Istoteq;
    }
#endif

#ifdef DEBUG_HISIM2LD    
    printf ("id ig ib is %12.5e %12.5e %12.5e %12.5e\n", ceqjd - ceqbd - cdreq - ceqqd + Idtoteq,
	    -(ceqqg + Igtoteq), ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq, 
	    cdreq - ceqbs + ceqjs + ceqqg + ceqqb + ceqqd + Istoteq);
#endif
    
    /*
     *  load y matrix
     */

    if ( !here->HSM2_corbnet ){
      gjbd = here->HSM2_gbd;
      gjbs = here->HSM2_gbs;
    } else
      gjbd = gjbs = 0.0;

#ifdef USE_OMP
    if (here->HSM2_corg == 1) {
      grg = here->HSM2_grg;
      here->HSM2_1 = grg;
      here->HSM2_2 = grg;
      here->HSM2_3 = grg;
      here->HSM2_4 = gcggb + grg + gIgtotg;
      here->HSM2_5 = gcgdb + gIgtotd;
      here->HSM2_6 = gcgsb + gIgtots;
      here->HSM2_7 = gcgbb + gIgtotb;
    } else {
      here->HSM2_8 = gcggb + gIgtotg;
      here->HSM2_9 = gcgdb + gIgtotd;
      here->HSM2_10 = gcgsb + gIgtots;
      here->HSM2_11 = gcgbb + gIgtotb;
    }

    here->HSM2_12 = here->HSM2drainConductance
      + here->HSM2_gds + here->HSM2_gbd + RevSum + gcddb + gbdpdp - gIdtotd;
    here->HSM2_13 = here->HSM2drainConductance;
    here->HSM2_14 = gm + gcdgb + gbdpg - gIdtotg;
    here->HSM2_15 = here->HSM2_gds + FwdSum - gcdsb - gbdpsp + gIdtots;
    here->HSM2_16 = gjbd - gmbs - gcdbb - gbdpb + gIdtotb;

    here->HSM2_17 = here->HSM2drainConductance;
    here->HSM2_18 = here->HSM2drainConductance;

    here->HSM2_19 = here->HSM2_gds + RevSum - gcsdb - gbspdp + gIstotd;
    here->HSM2_20 = gcsgb - gm + gbspg - gIstotg;
    here->HSM2_21 = here->HSM2sourceConductance
      + here->HSM2_gds + here->HSM2_gbs + FwdSum + gcssb + gbspsp - gIstots;
    here->HSM2_22 = here->HSM2sourceConductance;
    here->HSM2_23 = gjbs + gmbs - gcsbb - gbspb + gIstotb;

    here->HSM2_24 = here->HSM2sourceConductance;
    here->HSM2_25 = here->HSM2sourceConductance;

    here->HSM2_26 = gcbdb - gjbd + gbbdp - gIbtotd;
    here->HSM2_27 = gcbgb - here->HSM2_gbgs - gIbtotg;
    here->HSM2_28 = gcbsb - gjbs + gbbsp - gIbtots;
    here->HSM2_29 = gjbd + gjbs + gcbbb - here->HSM2_gbbs - gIbtotb;

    if (model->HSM2_cogidl) {
      /* stamp GIDL */
      here->HSM2_30 = here->HSM2_gigidlds;
      here->HSM2_31 = here->HSM2_gigidlgs;
      here->HSM2_32 = (here->HSM2_gigidlgs +
                               here->HSM2_gigidlds + here->HSM2_gigidlbs);
      here->HSM2_33 = here->HSM2_gigidlbs;
      here->HSM2_34 = here->HSM2_gigidlds;
      here->HSM2_35 = here->HSM2_gigidlgs;
      here->HSM2_36 = (here->HSM2_gigidlgs +
                               here->HSM2_gigidlds + here->HSM2_gigidlbs);
      here->HSM2_37 = here->HSM2_gigidlbs;
      /* stamp GISL */
      here->HSM2_38 = (here->HSM2_gigislsd +
                               here->HSM2_gigislgd + here->HSM2_gigislbd);
      here->HSM2_39 = here->HSM2_gigislgd;
      here->HSM2_40 = here->HSM2_gigislsd;
      here->HSM2_41 = here->HSM2_gigislbd;
      here->HSM2_42 = (here->HSM2_gigislgd +
                               here->HSM2_gigislsd + here->HSM2_gigislbd);
      here->HSM2_43 = here->HSM2_gigislgd;
      here->HSM2_44 = here->HSM2_gigislsd;
      here->HSM2_45 = here->HSM2_gigislbd;
    }

    if (here->HSM2_corbnet) { /* body resistance network */
      here->HSM2_46 = gcdbdb - here->HSM2_gbd;
      here->HSM2_47 = here->HSM2_gbs - gcsbsb;

      here->HSM2_48 = gcdbdb - here->HSM2_gbd;
      here->HSM2_49 = here->HSM2_gbd - gcdbdb
                        + here->HSM2_grbpd + here->HSM2_grbdb;
      here->HSM2_50 = here->HSM2_grbpd;
      here->HSM2_51 = here->HSM2_grbdb;

      here->HSM2_52 = here->HSM2_grbpd;
      here->HSM2_53 = here->HSM2_grbpb;
      here->HSM2_54 = here->HSM2_grbps;
      here->HSM2_55 = here->HSM2_grbpd + here->HSM2_grbps + here->HSM2_grbpb;

      here->HSM2_56 = gcsbsb - here->HSM2_gbs;
      here->HSM2_57 = here->HSM2_grbps;
      here->HSM2_58 = here->HSM2_grbsb;
      here->HSM2_59 = here->HSM2_gbs - gcsbsb
                        + here->HSM2_grbps + here->HSM2_grbsb;

      here->HSM2_60 = here->HSM2_grbdb;
      here->HSM2_61 = here->HSM2_grbpb;
      here->HSM2_62 = here->HSM2_grbsb;
      here->HSM2_63 = here->HSM2_grbsb + here->HSM2_grbdb + here->HSM2_grbpb;
    }
#else
    if (here->HSM2_corg == 1) {
      grg = here->HSM2_grg;
      *(here->HSM2GgPtr) += grg;
      *(here->HSM2GPgPtr) -= grg;
      *(here->HSM2GgpPtr) -= grg;
      *(here->HSM2GPgpPtr) += gcggb + grg + gIgtotg;
      *(here->HSM2GPdpPtr) += gcgdb + gIgtotd;
      *(here->HSM2GPspPtr) += gcgsb + gIgtots;
      *(here->HSM2GPbpPtr) += gcgbb + gIgtotb;
    } else {
      *(here->HSM2GPgpPtr) += gcggb + gIgtotg;
      *(here->HSM2GPdpPtr) += gcgdb + gIgtotd;
      *(here->HSM2GPspPtr) += gcgsb + gIgtots;
      *(here->HSM2GPbpPtr) += gcgbb + gIgtotb;
    }

    *(here->HSM2DPdpPtr) += here->HSM2drainConductance
      + here->HSM2_gds + here->HSM2_gbd + RevSum + gcddb + gbdpdp - gIdtotd;
    *(here->HSM2DPdPtr) -= here->HSM2drainConductance;
    *(here->HSM2DPgpPtr) += gm + gcdgb + gbdpg - gIdtotg;
    *(here->HSM2DPspPtr) -= here->HSM2_gds + FwdSum - gcdsb - gbdpsp + gIdtots;
    *(here->HSM2DPbpPtr) -= gjbd - gmbs - gcdbb - gbdpb + gIdtotb;

    *(here->HSM2DdpPtr) -= here->HSM2drainConductance;
    *(here->HSM2DdPtr) += here->HSM2drainConductance;

    *(here->HSM2SPdpPtr) -= here->HSM2_gds + RevSum - gcsdb - gbspdp + gIstotd;
    *(here->HSM2SPgpPtr) += gcsgb - gm + gbspg - gIstotg;
    *(here->HSM2SPspPtr) += here->HSM2sourceConductance
      + here->HSM2_gds + here->HSM2_gbs + FwdSum + gcssb + gbspsp - gIstots;
    *(here->HSM2SPsPtr) -= here->HSM2sourceConductance;
    *(here->HSM2SPbpPtr) -= gjbs + gmbs - gcsbb - gbspb + gIstotb;

    *(here->HSM2SspPtr) -= here->HSM2sourceConductance;
    *(here->HSM2SsPtr) += here->HSM2sourceConductance;

    *(here->HSM2BPdpPtr) += gcbdb - gjbd + gbbdp - gIbtotd;
    *(here->HSM2BPgpPtr) += gcbgb - here->HSM2_gbgs - gIbtotg;
    *(here->HSM2BPspPtr) += gcbsb - gjbs + gbbsp - gIbtots;
    *(here->HSM2BPbpPtr) += gjbd + gjbs + gcbbb - here->HSM2_gbbs - gIbtotb;

    if (model->HSM2_cogidl) {
      /* stamp GIDL */
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
      /* stamp GISL */
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

    if (here->HSM2_corbnet) { /* body resistance network */
      *(here->HSM2DPdbPtr) += gcdbdb - here->HSM2_gbd;
      *(here->HSM2SPsbPtr) -= here->HSM2_gbs - gcsbsb;

      *(here->HSM2DBdpPtr) += gcdbdb - here->HSM2_gbd;
      *(here->HSM2DBdbPtr) += here->HSM2_gbd - gcdbdb 
	+ here->HSM2_grbpd + here->HSM2_grbdb;
      *(here->HSM2DBbpPtr) -= here->HSM2_grbpd;
      *(here->HSM2DBbPtr) -= here->HSM2_grbdb;

      *(here->HSM2BPdbPtr) -= here->HSM2_grbpd;
      *(here->HSM2BPbPtr) -= here->HSM2_grbpb;
      *(here->HSM2BPsbPtr) -= here->HSM2_grbps;
      *(here->HSM2BPbpPtr) += here->HSM2_grbpd + here->HSM2_grbps + here->HSM2_grbpb;
      
      *(here->HSM2SBspPtr) += gcsbsb - here->HSM2_gbs;
      *(here->HSM2SBbpPtr) -= here->HSM2_grbps;
      *(here->HSM2SBbPtr) -= here->HSM2_grbsb;
      *(here->HSM2SBsbPtr) += here->HSM2_gbs - gcsbsb 
	+ here->HSM2_grbps + here->HSM2_grbsb;

      *(here->HSM2BdbPtr) -= here->HSM2_grbdb;
      *(here->HSM2BbpPtr) -= here->HSM2_grbpb;
      *(here->HSM2BsbPtr) -= here->HSM2_grbsb;
      *(here->HSM2BbPtr) += here->HSM2_grbsb + here->HSM2_grbdb + here->HSM2_grbpb;
    }
#endif

  line1000:
    ;
    
#ifndef USE_OMP
   } /* End of MOSFET Instance */
  } /* End of Model Instance */
#endif
  
#ifdef MOS_MODEL_TIME
tm1 = gtodsecld() ;
mos_model_time += ( tm1 - tm0 ) ;
sprintf( mos_model_name , "HiSIM 240BSC1" ) ;
#ifdef PARAMOS_TIME
vsum = vbs + vds + vgs ;
if ( vsum < vsum0 - 1e-6 || vsum > vsum0 + 1e-6 ) {
printf( "PMVbs= %12.5e\n" , vbs ) ;
printf( "PMVds= %12.5e\n" , vds ) ;
printf( "PMVgs= %12.5e\n" , vgs ) ;
printf( "PMTime= %12.5e\n" , tm1 - tm0 ) ;
}
vsum0 = vsum ;
#endif
#endif

  return(OK);
}

#ifdef USE_OMP
void HSM2LoadRhsMat(GENmodel *inModel, CKTcircuit *ckt)
{
    unsigned int InstCount, idx;
    HSM2instance **InstArray;
    HSM2instance *here;
    HSM2model *model = (HSM2model*)inModel;

    InstArray = model->HSM2InstanceArray;
    InstCount = model->HSM2InstCount;

    for (idx = 0; idx < InstCount; idx++) {
       here = InstArray[idx];
        /* Update b for Ax = b */
        *(ckt->CKTrhs + here->HSM2dNodePrime) += here->HSM2rhsdPrime;
        *(ckt->CKTrhs + here->HSM2gNodePrime) -= here->HSM2rhsgPrime;

        if ( !here->HSM2_corbnet ) {
          *(ckt->CKTrhs + here->HSM2bNodePrime) += here->HSM2rhsbPrime;
          *(ckt->CKTrhs + here->HSM2sNodePrime) += here->HSM2rhssPrime;
        } else {
          *(ckt->CKTrhs + here->HSM2dbNode) -= here->HSM2rhsdb;
          *(ckt->CKTrhs + here->HSM2bNodePrime) += here->HSM2rhsbPrime;
          *(ckt->CKTrhs + here->HSM2sbNode) -= here->HSM2rhssb;
          *(ckt->CKTrhs + here->HSM2sNodePrime) += here->HSM2rhssPrime;
       }

       /* Update A for Ax = b */
       if (here->HSM2_corg == 1) {
         *(here->HSM2GgPtr) += here->HSM2_1;
         *(here->HSM2GPgPtr) -= here->HSM2_2;
         *(here->HSM2GgpPtr) -= here->HSM2_3;
         *(here->HSM2GPgpPtr) += here->HSM2_4;
         *(here->HSM2GPdpPtr) += here->HSM2_5;
         *(here->HSM2GPspPtr) += here->HSM2_6;
         *(here->HSM2GPbpPtr) += here->HSM2_7;
       } else {
         *(here->HSM2GPgpPtr) += here->HSM2_8;
         *(here->HSM2GPdpPtr) += here->HSM2_9;
         *(here->HSM2GPspPtr) += here->HSM2_10;
         *(here->HSM2GPbpPtr) += here->HSM2_11;
       }

       *(here->HSM2DPdpPtr) += here->HSM2_12;

       *(here->HSM2DPdPtr) -= here->HSM2_13;
       *(here->HSM2DPgpPtr) += here->HSM2_14;
       *(here->HSM2DPspPtr) -= here->HSM2_15;
       *(here->HSM2DPbpPtr) -= here->HSM2_16;

       *(here->HSM2DdpPtr) -= here->HSM2_17;
       *(here->HSM2DdPtr) += here->HSM2_18;

       *(here->HSM2SPdpPtr) -= here->HSM2_19;
       *(here->HSM2SPgpPtr) += here->HSM2_20;
       *(here->HSM2SPspPtr) += here->HSM2_21;

       *(here->HSM2SPsPtr) -= here->HSM2_22;
       *(here->HSM2SPbpPtr) -= here->HSM2_23;

       *(here->HSM2SspPtr) -= here->HSM2_24;
       *(here->HSM2SsPtr) += here->HSM2_25;

       *(here->HSM2BPdpPtr) += here->HSM2_26;
       *(here->HSM2BPgpPtr) += here->HSM2_27;
       *(here->HSM2BPspPtr) += here->HSM2_28;
       *(here->HSM2BPbpPtr) += here->HSM2_29;

       if (model->HSM2_cogidl) {
         /* stamp GIDL */
         *(here->HSM2DPdpPtr) += here->HSM2_30;
         *(here->HSM2DPgpPtr) += here->HSM2_31;
         *(here->HSM2DPspPtr) -= here->HSM2_32;

         *(here->HSM2DPbpPtr) += here->HSM2_33;
         *(here->HSM2BPdpPtr) -= here->HSM2_34;
         *(here->HSM2BPgpPtr) -= here->HSM2_35;
         *(here->HSM2BPspPtr) += here->HSM2_36;

         *(here->HSM2BPbpPtr) -= here->HSM2_37;
         /* stamp GISL */
         *(here->HSM2SPdpPtr) -= here->HSM2_38;

         *(here->HSM2SPgpPtr) += here->HSM2_39;
         *(here->HSM2SPspPtr) += here->HSM2_40;
         *(here->HSM2SPbpPtr) += here->HSM2_41;
         *(here->HSM2BPdpPtr) += here->HSM2_42;

         *(here->HSM2BPgpPtr) -= here->HSM2_43;
         *(here->HSM2BPspPtr) -= here->HSM2_44;
         *(here->HSM2BPbpPtr) -= here->HSM2_45;
       }

       if (here->HSM2_corbnet) { /* body resistance network */
         *(here->HSM2DPdbPtr) += here->HSM2_46;
         *(here->HSM2SPsbPtr) -= here->HSM2_47;

         *(here->HSM2DBdpPtr) += here->HSM2_48;
         *(here->HSM2DBdbPtr) += here->HSM2_49;

         *(here->HSM2DBbpPtr) -= here->HSM2_50;
         *(here->HSM2DBbPtr) -= here->HSM2_51;

         *(here->HSM2BPdbPtr) -= here->HSM2_52;
         *(here->HSM2BPbPtr) -= here->HSM2_53;
         *(here->HSM2BPsbPtr) -= here->HSM2_54;
         *(here->HSM2BPbpPtr) += here->HSM2_55;

         *(here->HSM2SBspPtr) += here->HSM2_56;
         *(here->HSM2SBbpPtr) -= here->HSM2_57;
         *(here->HSM2SBbPtr) -= here->HSM2_58;
         *(here->HSM2SBsbPtr) += here->HSM2_59;


         *(here->HSM2BdbPtr) -= here->HSM2_60;
         *(here->HSM2BbpPtr) -= here->HSM2_61;
         *(here->HSM2BsbPtr) -= here->HSM2_62;
         *(here->HSM2BbPtr) += here->HSM2_63;
       }
    }
}
#endif
