/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4trunc.c of BSIM4.7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


int
BSIM4v7trunc(
GENmodel *inModel,
CKTcircuit *ckt,
double *timeStep)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;

#ifdef STEPDEBUG
    double debugtemp;
#endif /* STEPDEBUG */

    for (; model != NULL; model = model->BSIM4v7nextModel)
    {    for (here = model->BSIM4v7instances; here != NULL;
	      here = here->BSIM4v7nextInstance)
	      {
#ifdef STEPDEBUG
            debugtemp = *timeStep;
#endif /* STEPDEBUG */
            CKTterr(here->BSIM4v7qb,ckt,timeStep);
            CKTterr(here->BSIM4v7qg,ckt,timeStep);
            CKTterr(here->BSIM4v7qd,ckt,timeStep);
            if (here->BSIM4v7trnqsMod)
                CKTterr(here->BSIM4v7qcdump,ckt,timeStep);
            if (here->BSIM4v7rbodyMod)
            {   CKTterr(here->BSIM4v7qbs,ckt,timeStep);
                CKTterr(here->BSIM4v7qbd,ckt,timeStep);
	          }
	          if (here->BSIM4v7rgateMod == 3)
		        CKTterr(here->BSIM4v7qgmid,ckt,timeStep);
#ifdef STEPDEBUG
            if(debugtemp != *timeStep)
	          {  printf("device %s reduces step from %g to %g\n",
                       here->BSIM4v7name,debugtemp,*timeStep);
            }
#endif /* STEPDEBUG */
        }
    }
    return(OK);
}
