/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4trunc.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
BSIM4trunc(
GENmodel *inModel,
CKTcircuit *ckt,
double *timeStep)
{
BSIM4model *model = (BSIM4model*)inModel;
BSIM4instance *here;

#ifdef STEPDEBUG
    double debugtemp;
#endif /* STEPDEBUG */

    for (; model != NULL; model = model->BSIM4nextModel)
    {    for (here = model->BSIM4instances; here != NULL;
	      here = here->BSIM4nextInstance)
	      {
#ifdef STEPDEBUG
            debugtemp = *timeStep;
#endif /* STEPDEBUG */
            CKTterr(here->BSIM4qb,ckt,timeStep);
            CKTterr(here->BSIM4qg,ckt,timeStep);
            CKTterr(here->BSIM4qd,ckt,timeStep);
            if (here->BSIM4trnqsMod)
                CKTterr(here->BSIM4qcdump,ckt,timeStep);
            if (here->BSIM4rbodyMod)
            {   CKTterr(here->BSIM4qbs,ckt,timeStep);
                CKTterr(here->BSIM4qbd,ckt,timeStep);
	          }
	          if (here->BSIM4rgateMod == 3)
		        CKTterr(here->BSIM4qgmid,ckt,timeStep);
#ifdef STEPDEBUG
            if(debugtemp != *timeStep)
	          {  printf("device %s reduces step from %g to %g\n",
                       here->BSIM4name,debugtemp,*timeStep);
            }
#endif /* STEPDEBUG */
        }
    }
    return(OK);
}
