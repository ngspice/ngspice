/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM
 ( VERSION : 2  SUBVERSION : 7  REVISION : 0 ) Beta
 
 FILE : hsm2trunc.c

 Date : 2012.10.25

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "hsm2def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int HSM2trunc(
     GENmodel *inModel,
     CKTcircuit *ckt,
     double *timeStep)

{
  HSM2model *model = (HSM2model*)inModel;
  HSM2instance *here;
#ifdef STEPDEBUG
  double debugtemp;
#endif /* STEPDEBUG */
  
  for ( ;model != NULL ;model = model->HSM2nextModel ) {
    for ( here=model->HSM2instances ;here!=NULL ;
	  here = here->HSM2nextInstance ) {
#ifdef STEPDEBUG
      debugtemp = *timeStep;
#endif /* STEPDEBUG */
      CKTterr(here->HSM2qb,ckt,timeStep);
      CKTterr(here->HSM2qg,ckt,timeStep);
      CKTterr(here->HSM2qd,ckt,timeStep);
#ifdef STEPDEBUG
      if ( debugtemp != *timeStep ) 
	printf("device %s reduces step from %g to %g\n",
	       here->HSM2name, debugtemp, *timeStep);
#endif /* STEPDEBUG */
    }
  }
  return(OK);
}


