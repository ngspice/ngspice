/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4dest.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "bsim4def.h"
#include "ngspice/suffix.h"

void
BSIM4destroy(
    GENmodel **inModel)
{
    BSIM4model **model = (BSIM4model**)inModel;
    BSIM4instance *here;
    BSIM4instance *prev = NULL;
    BSIM4model *mod = *model;
    BSIM4model *oldmod = NULL;

    for (; mod ; mod = mod->BSIM4nextModel) {
    /** added to get rid of link list pSizeDependParamKnot **/      
        struct bsim4SizeDependParam *pParam, *pParamOld=NULL;

        pParam = mod->pSizeDependParamKnot;

        for (; pParam ; pParam = pParam->pNext) {
            FREE(pParamOld);
            pParamOld = pParam;
        }
        FREE(pParamOld);
        pParam = NULL;
     /** end of extra code **/
        if(oldmod) {
            FREE(oldmod->BSIM4version);
            FREE(oldmod);
        }
        oldmod = mod;
        prev = (BSIM4instance *)NULL;
        for (here = mod->BSIM4instances; here; here = here->BSIM4nextInstance) {
            if(prev) FREE(prev);
            prev = here;
        }
        if(prev) FREE(prev);
    }
    if(oldmod) {
#ifdef USE_OMP
        /* free just once for all models */
        FREE(oldmod->BSIM4InstanceArray);
#endif
        /* oldmod->BSIM4modName to be freed in INPtabEnd() */
        FREE(oldmod->BSIM4version);
        FREE(oldmod);
    }
    *model = NULL;
    return;
}
