/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4del.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "bsim4v7def.h"
#include "ngspice/sperror.h"
#include "ngspice/gendefs.h"
#include "ngspice/suffix.h"


int
BSIM4v7delete(
GENmodel *inModel,
IFuid name,
GENinstance **inInst)
{
BSIM4v7instance **fast = (BSIM4v7instance**)inInst;
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance **prev = NULL;
BSIM4v7instance *here;

    for (; model ; model = model->BSIM4v7nextModel) 
    {    prev = &(model->BSIM4v7instances);
         for (here = *prev; here ; here = *prev) 
	 {    if (here->BSIM4v7name == name || (fast && here==*fast))
	      {   *prev= here->BSIM4v7nextInstance;
                  FREE(here);
                  return(OK);
              }
              prev = &(here->BSIM4v7nextInstance);
         }
    }
    return(E_NODEV);
}
