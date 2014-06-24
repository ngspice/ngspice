/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4mdel.c of BSIM4.7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "bsim4v7def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4v7mDelete(
GENmodel **inModel,
IFuid modname,
GENmodel *kill)
{
BSIM4v7model **model = (BSIM4v7model**)inModel;
BSIM4v7model *modfast = (BSIM4v7model*)kill;
BSIM4v7instance *here;
BSIM4v7instance *prev = NULL;
BSIM4v7model **oldmod;

    oldmod = model;
    for (; *model ; model = &((*model)->BSIM4v7nextModel)) 
    {    if ((*model)->BSIM4v7modName == modname || 
             (modfast && *model == modfast))
	     goto delgot;
         oldmod = model;
    }
    return(E_NOMOD);

delgot:
    *oldmod = (*model)->BSIM4v7nextModel; /* cut deleted device out of list */
    for (here = (*model)->BSIM4v7instances; here; here = here->BSIM4v7nextInstance)
    {    if(prev) FREE(prev);
         prev = here;
    }
    if(prev) FREE(prev);
    FREE(*model);
    return(OK);
}
