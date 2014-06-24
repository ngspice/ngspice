/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4getic.c of BSIM4.7.0.
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
BSIM4v7getic(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;

    for (; model ; model = model->BSIM4v7nextModel) 
    {    for (here = model->BSIM4v7instances; here; here = here->BSIM4v7nextInstance)
          {
              if (!here->BSIM4v7icVDSGiven) 
              {   here->BSIM4v7icVDS = *(ckt->CKTrhs + here->BSIM4v7dNode) 
                                   - *(ckt->CKTrhs + here->BSIM4v7sNode);
              }
              if (!here->BSIM4v7icVGSGiven) 
              {   here->BSIM4v7icVGS = *(ckt->CKTrhs + here->BSIM4v7gNodeExt) 
                                   - *(ckt->CKTrhs + here->BSIM4v7sNode);
              }
              if(!here->BSIM4v7icVBSGiven)
              {  here->BSIM4v7icVBS = *(ckt->CKTrhs + here->BSIM4v7bNode)
                                  - *(ckt->CKTrhs + here->BSIM4v7sNode);
              }
         }
    }
    return(OK);
}
