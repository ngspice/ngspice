/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4par.c of BSIM4.8.0.
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
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/fteext.h"

int
BSIM4param(
int param,
IFvalue *value,
GENinstance *inst,
IFvalue *select)
{
    double scale;

    BSIM4instance *here = (BSIM4instance*)inst;

    NG_IGNORE(select);

    if (!cp_getvar("scale", CP_REAL, &scale))
        scale = 1;

    switch(param) 
    {   case BSIM4_W:
            here->BSIM4w = value->rValue * scale;
            here->BSIM4wGiven = TRUE;
            break;
        case BSIM4_L:
            here->BSIM4l = value->rValue * scale;
            here->BSIM4lGiven = TRUE;
            break;
        case BSIM4_M:
            here->BSIM4m = value->rValue;
            here->BSIM4mGiven = TRUE;
            break;
        case BSIM4_NF:
            here->BSIM4nf = value->rValue;
            here->BSIM4nfGiven = TRUE;
            break;
        case BSIM4_MIN:
            here->BSIM4min = value->iValue;
            here->BSIM4minGiven = TRUE;
            break;
        case BSIM4_AS:
            here->BSIM4sourceArea = value->rValue * scale * scale;
            here->BSIM4sourceAreaGiven = TRUE;
            break;
        case BSIM4_AD:
            here->BSIM4drainArea = value->rValue * scale * scale;
            here->BSIM4drainAreaGiven = TRUE;
            break;
        case BSIM4_PS:
            here->BSIM4sourcePerimeter = value->rValue * scale;
            here->BSIM4sourcePerimeterGiven = TRUE;
            break;
        case BSIM4_PD:
            here->BSIM4drainPerimeter = value->rValue * scale;
            here->BSIM4drainPerimeterGiven = TRUE;
            break;
        case BSIM4_NRS:
            here->BSIM4sourceSquares = value->rValue;
            here->BSIM4sourceSquaresGiven = TRUE;
            break;
        case BSIM4_NRD:
            here->BSIM4drainSquares = value->rValue;
            here->BSIM4drainSquaresGiven = TRUE;
            break;
        case BSIM4_OFF:
            here->BSIM4off = value->iValue;
            break;
        case BSIM4_SA:
            here->BSIM4sa = value->rValue;
            here->BSIM4saGiven = TRUE;
            break;
        case BSIM4_SB:
            here->BSIM4sb = value->rValue;
            here->BSIM4sbGiven = TRUE;
            break;
        case BSIM4_SD:
            here->BSIM4sd = value->rValue;
            here->BSIM4sdGiven = TRUE;
            break;
        case BSIM4_SCA:
            here->BSIM4sca = value->rValue;
            here->BSIM4scaGiven = TRUE;
            break;
        case BSIM4_SCB:
            here->BSIM4scb = value->rValue;
            here->BSIM4scbGiven = TRUE;
            break;
        case BSIM4_SCC:
            here->BSIM4scc = value->rValue;
            here->BSIM4sccGiven = TRUE;
            break;
        case BSIM4_SC:
            here->BSIM4sc = value->rValue;
            here->BSIM4scGiven = TRUE;
            break;
        case BSIM4_RBSB:
            here->BSIM4rbsb = value->rValue;
            here->BSIM4rbsbGiven = TRUE;
            break;
        case BSIM4_RBDB:
            here->BSIM4rbdb = value->rValue;
            here->BSIM4rbdbGiven = TRUE;
            break;
        case BSIM4_RBPB:
            here->BSIM4rbpb = value->rValue;
            here->BSIM4rbpbGiven = TRUE;
            break;
        case BSIM4_RBPS:
            here->BSIM4rbps = value->rValue;
            here->BSIM4rbpsGiven = TRUE;
            break;
        case BSIM4_RBPD:
            here->BSIM4rbpd = value->rValue;
            here->BSIM4rbpdGiven = TRUE;
            break;
        case BSIM4_DELVTO:
            here->BSIM4delvto = value->rValue;
            here->BSIM4delvtoGiven = TRUE;
            break;
        case BSIM4_XGW:
            here->BSIM4xgw = value->rValue;
            here->BSIM4xgwGiven = TRUE;
            break;
        case BSIM4_NGCON:
            here->BSIM4ngcon = value->rValue;
            here->BSIM4ngconGiven = TRUE;
            break;
        case BSIM4_TRNQSMOD:
            here->BSIM4trnqsMod = value->iValue;
            here->BSIM4trnqsModGiven = TRUE;
            break;
        case BSIM4_ACNQSMOD:
            here->BSIM4acnqsMod = value->iValue;
            here->BSIM4acnqsModGiven = TRUE;
            break;
        case BSIM4_RBODYMOD:
            here->BSIM4rbodyMod = value->iValue;
            here->BSIM4rbodyModGiven = TRUE;
            break;
        case BSIM4_RGATEMOD:
            here->BSIM4rgateMod = value->iValue;
            here->BSIM4rgateModGiven = TRUE;
            break;
        case BSIM4_GEOMOD:
            here->BSIM4geoMod = value->iValue;
            here->BSIM4geoModGiven = TRUE;
            break;
        case BSIM4_RGEOMOD:
            here->BSIM4rgeoMod = value->iValue;
            here->BSIM4rgeoModGiven = TRUE;
            break;
        case BSIM4_IC_VDS:
            here->BSIM4icVDS = value->rValue;
            here->BSIM4icVDSGiven = TRUE;
            break;
        case BSIM4_IC_VGS:
            here->BSIM4icVGS = value->rValue;
            here->BSIM4icVGSGiven = TRUE;
            break;
        case BSIM4_IC_VBS:
            here->BSIM4icVBS = value->rValue;
            here->BSIM4icVBSGiven = TRUE;
            break;
        case BSIM4_IC:
            switch(value->v.numValue)
            {   case 3:
                    here->BSIM4icVBS = *(value->v.vec.rVec+2);
                    here->BSIM4icVBSGiven = TRUE;
                case 2:
                    here->BSIM4icVGS = *(value->v.vec.rVec+1);
                    here->BSIM4icVGSGiven = TRUE;
                case 1:
                    here->BSIM4icVDS = *(value->v.vec.rVec);
                    here->BSIM4icVDSGiven = TRUE;
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
