/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4del.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/gendefs.h"
#include "ngspice/suffix.h"


int
BSIM4delete(
GENmodel *inModel,
IFuid name,
GENinstance **inInst)
{
BSIM4instance **fast = (BSIM4instance**)inInst;
BSIM4model *model = (BSIM4model*)inModel;
BSIM4instance **prev = NULL;
BSIM4instance *here;

    for (; model ; model = model->BSIM4nextModel) 
    {    prev = &(model->BSIM4instances);
         for (here = *prev; here ; here = *prev) 
	 {    if (here->BSIM4name == name || (fast && here==*fast))
	      {   *prev= here->BSIM4nextInstance;
                  FREE(here);
                  return(OK);
              }
              prev = &(here->BSIM4nextInstance);
         }
    }
    return(E_NODEV);
}
