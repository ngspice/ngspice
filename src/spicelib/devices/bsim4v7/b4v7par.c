/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4par.c of BSIM4.7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 11/15/2002.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, Mohan Dunga, 07/29/2005.
**********/

#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "bsim4v7def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/fteext.h"

int
BSIM4v7param(
int param,
IFvalue *value,
GENinstance *inst,
IFvalue *select)
{
    double scale;

    BSIM4v7instance *here = (BSIM4v7instance*)inst;

    NG_IGNORE(select);

    if (!cp_getvar("scale", CP_REAL, &scale))
        scale = 1;

    switch(param) 
    {   case BSIM4v7_W:
            here->BSIM4v7w = value->rValue * scale;
            here->BSIM4v7wGiven = TRUE;
            break;
        case BSIM4v7_L:
            here->BSIM4v7l = value->rValue * scale;
            here->BSIM4v7lGiven = TRUE;
            break;
        case BSIM4v7_M:
            here->BSIM4v7m = value->rValue;
            here->BSIM4v7mGiven = TRUE;
            break;
        case BSIM4v7_NF:
            here->BSIM4v7nf = value->rValue;
            here->BSIM4v7nfGiven = TRUE;
            break;
        case BSIM4v7_MIN:
            here->BSIM4v7min = value->iValue;
            here->BSIM4v7minGiven = TRUE;
            break;
        case BSIM4v7_AS:
            here->BSIM4v7sourceArea = value->rValue * scale * scale;
            here->BSIM4v7sourceAreaGiven = TRUE;
            break;
        case BSIM4v7_AD:
            here->BSIM4v7drainArea = value->rValue * scale * scale;
            here->BSIM4v7drainAreaGiven = TRUE;
            break;
        case BSIM4v7_PS:
            here->BSIM4v7sourcePerimeter = value->rValue * scale;
            here->BSIM4v7sourcePerimeterGiven = TRUE;
            break;
        case BSIM4v7_PD:
            here->BSIM4v7drainPerimeter = value->rValue * scale;
            here->BSIM4v7drainPerimeterGiven = TRUE;
            break;
        case BSIM4v7_NRS:
            here->BSIM4v7sourceSquares = value->rValue;
            here->BSIM4v7sourceSquaresGiven = TRUE;
            break;
        case BSIM4v7_NRD:
            here->BSIM4v7drainSquares = value->rValue;
            here->BSIM4v7drainSquaresGiven = TRUE;
            break;
        case BSIM4v7_OFF:
            here->BSIM4v7off = value->iValue;
            break;
        case BSIM4v7_SA:
            here->BSIM4v7sa = value->rValue;
            here->BSIM4v7saGiven = TRUE;
            break;
        case BSIM4v7_SB:
            here->BSIM4v7sb = value->rValue;
            here->BSIM4v7sbGiven = TRUE;
            break;
        case BSIM4v7_SD:
            here->BSIM4v7sd = value->rValue;
            here->BSIM4v7sdGiven = TRUE;
            break;
        case BSIM4v7_SCA:
            here->BSIM4v7sca = value->rValue;
            here->BSIM4v7scaGiven = TRUE;
            break;
        case BSIM4v7_SCB:
            here->BSIM4v7scb = value->rValue;
            here->BSIM4v7scbGiven = TRUE;
            break;
        case BSIM4v7_SCC:
            here->BSIM4v7scc = value->rValue;
            here->BSIM4v7sccGiven = TRUE;
            break;
        case BSIM4v7_SC:
            here->BSIM4v7sc = value->rValue;
            here->BSIM4v7scGiven = TRUE;
            break;
        case BSIM4v7_RBSB:
            here->BSIM4v7rbsb = value->rValue;
            here->BSIM4v7rbsbGiven = TRUE;
            break;
        case BSIM4v7_RBDB:
            here->BSIM4v7rbdb = value->rValue;
            here->BSIM4v7rbdbGiven = TRUE;
            break;
        case BSIM4v7_RBPB:
            here->BSIM4v7rbpb = value->rValue;
            here->BSIM4v7rbpbGiven = TRUE;
            break;
        case BSIM4v7_RBPS:
            here->BSIM4v7rbps = value->rValue;
            here->BSIM4v7rbpsGiven = TRUE;
            break;
        case BSIM4v7_RBPD:
            here->BSIM4v7rbpd = value->rValue;
            here->BSIM4v7rbpdGiven = TRUE;
            break;
        case BSIM4v7_DELVTO:
            here->BSIM4v7delvto = value->rValue;
            here->BSIM4v7delvtoGiven = TRUE;
            break;
        case BSIM4v7_XGW:
            here->BSIM4v7xgw = value->rValue;
            here->BSIM4v7xgwGiven = TRUE;
            break;
        case BSIM4v7_NGCON:
            here->BSIM4v7ngcon = value->rValue;
            here->BSIM4v7ngconGiven = TRUE;
            break;
        case BSIM4v7_TRNQSMOD:
            here->BSIM4v7trnqsMod = value->iValue;
            here->BSIM4v7trnqsModGiven = TRUE;
            break;
        case BSIM4v7_ACNQSMOD:
            here->BSIM4v7acnqsMod = value->iValue;
            here->BSIM4v7acnqsModGiven = TRUE;
            break;
        case BSIM4v7_RBODYMOD:
            here->BSIM4v7rbodyMod = value->iValue;
            here->BSIM4v7rbodyModGiven = TRUE;
            break;
        case BSIM4v7_RGATEMOD:
            here->BSIM4v7rgateMod = value->iValue;
            here->BSIM4v7rgateModGiven = TRUE;
            break;
        case BSIM4v7_GEOMOD:
            here->BSIM4v7geoMod = value->iValue;
            here->BSIM4v7geoModGiven = TRUE;
            break;
        case BSIM4v7_RGEOMOD:
            here->BSIM4v7rgeoMod = value->iValue;
            here->BSIM4v7rgeoModGiven = TRUE;
            break;
        case BSIM4v7_IC_VDS:
            here->BSIM4v7icVDS = value->rValue;
            here->BSIM4v7icVDSGiven = TRUE;
            break;
        case BSIM4v7_IC_VGS:
            here->BSIM4v7icVGS = value->rValue;
            here->BSIM4v7icVGSGiven = TRUE;
            break;
        case BSIM4v7_IC_VBS:
            here->BSIM4v7icVBS = value->rValue;
            here->BSIM4v7icVBSGiven = TRUE;
            break;
        case BSIM4v7_IC:
            switch(value->v.numValue)
            {   case 3:
                    here->BSIM4v7icVBS = *(value->v.vec.rVec+2);
                    here->BSIM4v7icVBSGiven = TRUE;
                case 2:
                    here->BSIM4v7icVGS = *(value->v.vec.rVec+1);
                    here->BSIM4v7icVGSGiven = TRUE;
                case 1:
                    here->BSIM4v7icVDS = *(value->v.vec.rVec);
                    here->BSIM4v7icVDSGiven = TRUE;
                    break;
                default:
                    return(E_BADPARM);
            }
            break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}
