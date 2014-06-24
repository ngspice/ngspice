/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4dest.c of BSIM4.7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "bsim4v7def.h"
#include "ngspice/suffix.h"

void
BSIM4v7destroy(
    GENmodel **inModel)
{
    BSIM4v7model **model = (BSIM4v7model**)inModel;
    BSIM4v7instance *here;
    BSIM4v7instance *prev = NULL;
    BSIM4v7model *mod = *model;
    BSIM4v7model *oldmod = NULL;

    for (; mod ; mod = mod->BSIM4v7nextModel) {
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
            FREE(oldmod->BSIM4v7version);
            FREE(oldmod);
        }
        oldmod = mod;
        prev = (BSIM4v7instance *)NULL;
        for (here = mod->BSIM4v7instances; here; here = here->BSIM4v7nextInstance) {
            if(prev) FREE(prev);
            prev = here;
        }
        if(prev) FREE(prev);
    }
    if(oldmod) {
#ifdef USE_OMP
        /* free just once for all models */
        FREE(oldmod->BSIM4v7InstanceArray);
#endif
        /* oldmod->BSIM4v7modName to be freed in INPtabEnd() */
        FREE(oldmod->BSIM4v7version);
        FREE(oldmod);
    }
    *model = NULL;
    return;
}
