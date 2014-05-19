/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4ask.c of BSIM4.8.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, Mohan Dunga, 07/29/2005.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4ask(
CKTcircuit *ckt,
GENinstance *inst,
int which,
IFvalue *value,
IFvalue *select)
{
BSIM4instance *here = (BSIM4instance*)inst;

    NG_IGNORE(select);

    switch(which) 
    {   case BSIM4_L:
            value->rValue = here->BSIM4l;
            return(OK);
        case BSIM4_W:
            value->rValue = here->BSIM4w;
            return(OK);
        case BSIM4_M:
            value->rValue = here->BSIM4m;
            return(OK);
        case BSIM4_NF:
            value->rValue = here->BSIM4nf;
            return(OK);
        case BSIM4_MIN:
            value->iValue = here->BSIM4min;
            return(OK);
        case BSIM4_AS:
            value->rValue = here->BSIM4sourceArea;
            return(OK);
        case BSIM4_AD:
            value->rValue = here->BSIM4drainArea;
            return(OK);
        case BSIM4_PS:
            value->rValue = here->BSIM4sourcePerimeter;
            return(OK);
        case BSIM4_PD:
            value->rValue = here->BSIM4drainPerimeter;
            return(OK);
        case BSIM4_NRS:
            value->rValue = here->BSIM4sourceSquares;
            return(OK);
        case BSIM4_NRD:
            value->rValue = here->BSIM4drainSquares;
            return(OK);
        case BSIM4_OFF:
            value->rValue = here->BSIM4off;
            return(OK);
        case BSIM4_SA:
            value->rValue = here->BSIM4sa ;
            return(OK);
        case BSIM4_SB:
            value->rValue = here->BSIM4sb ;
            return(OK);
        case BSIM4_SD:
            value->rValue = here->BSIM4sd ;
            return(OK);
        case BSIM4_SCA:
            value->rValue = here->BSIM4sca ;
            return(OK);
        case BSIM4_SCB:
            value->rValue = here->BSIM4scb ;
            return(OK);
        case BSIM4_SCC:
            value->rValue = here->BSIM4scc ;
            return(OK);
        case BSIM4_SC:
            value->rValue = here->BSIM4sc ;
            return(OK);

        case BSIM4_RBSB:
            value->rValue = here->BSIM4rbsb;
            return(OK);
        case BSIM4_RBDB:
            value->rValue = here->BSIM4rbdb;
            return(OK);
        case BSIM4_RBPB:
            value->rValue = here->BSIM4rbpb;
            return(OK);
        case BSIM4_RBPS:
            value->rValue = here->BSIM4rbps;
            return(OK);
        case BSIM4_RBPD:
            value->rValue = here->BSIM4rbpd;
            return(OK);
        case BSIM4_DELVTO:
            value->rValue = here->BSIM4delvto;
            return(OK);
        case BSIM4_XGW:
            value->rValue = here->BSIM4xgw;
            return(OK);
        case BSIM4_NGCON:
            value->rValue = here->BSIM4ngcon;
            return(OK);
        case BSIM4_TRNQSMOD:
            value->iValue = here->BSIM4trnqsMod;
            return(OK);
        case BSIM4_ACNQSMOD:
            value->iValue = here->BSIM4acnqsMod;
            return(OK);
        case BSIM4_RBODYMOD:
            value->iValue = here->BSIM4rbodyMod;
            return(OK);
        case BSIM4_RGATEMOD:
            value->iValue = here->BSIM4rgateMod;
            return(OK);
        case BSIM4_GEOMOD:
            value->iValue = here->BSIM4geoMod;
            return(OK);
        case BSIM4_RGEOMOD:
            value->iValue = here->BSIM4rgeoMod;
            return(OK);
        case BSIM4_IC_VDS:
            value->rValue = here->BSIM4icVDS;
            return(OK);
        case BSIM4_IC_VGS:
            value->rValue = here->BSIM4icVGS;
            return(OK);
        case BSIM4_IC_VBS:
            value->rValue = here->BSIM4icVBS;
            return(OK);
        case BSIM4_DNODE:
            value->iValue = here->BSIM4dNode;
            return(OK);
        case BSIM4_GNODEEXT:
            value->iValue = here->BSIM4gNodeExt;
            return(OK);
        case BSIM4_SNODE:
            value->iValue = here->BSIM4sNode;
            return(OK);
        case BSIM4_BNODE:
            value->iValue = here->BSIM4bNode;
            return(OK);
        case BSIM4_DNODEPRIME:
            value->iValue = here->BSIM4dNodePrime;
            return(OK);
        case BSIM4_GNODEPRIME:
            value->iValue = here->BSIM4gNodePrime;
            return(OK);
        case BSIM4_GNODEMID:
            value->iValue = here->BSIM4gNodeMid;
            return(OK);
        case BSIM4_SNODEPRIME:
            value->iValue = here->BSIM4sNodePrime;
            return(OK);
        case BSIM4_DBNODE:
            value->iValue = here->BSIM4dbNode;
            return(OK);
        case BSIM4_BNODEPRIME:
            value->iValue = here->BSIM4bNodePrime;
            return(OK);
        case BSIM4_SBNODE:
            value->iValue = here->BSIM4sbNode;
            return(OK);
        case BSIM4_SOURCECONDUCT:
            value->rValue = here->BSIM4sourceConductance;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_DRAINCONDUCT:
            value->rValue = here->BSIM4drainConductance;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_VBD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4vbd);
            return(OK);
        case BSIM4_VBS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4vbs);
            return(OK);
        case BSIM4_VGS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4vgs);
            return(OK);
        case BSIM4_VDS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4vds);
            return(OK);
        case BSIM4_CD:
            value->rValue = here->BSIM4cd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBS:
            value->rValue = here->BSIM4cbs; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBD:
            value->rValue = here->BSIM4cbd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CSUB:
            value->rValue = here->BSIM4csub; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_QINV:
            value->rValue = here-> BSIM4qinv; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGIDL:
            value->rValue = here->BSIM4Igidl; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGISL:
            value->rValue = here->BSIM4Igisl; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGS:
            value->rValue = here->BSIM4Igs; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGD:
            value->rValue = here->BSIM4Igd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGB:
            value->rValue = here->BSIM4Igb; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGCS:
            value->rValue = here->BSIM4Igcs; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_IGCD:
            value->rValue = here->BSIM4Igcd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GM:
            value->rValue = here->BSIM4gm; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GDS:
            value->rValue = here->BSIM4gds; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GMBS:
            value->rValue = here->BSIM4gmbs; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GBD:
            value->rValue = here->BSIM4gbd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GBS:
            value->rValue = here->BSIM4gbs; 
            value->rValue *= here->BSIM4m;
            return(OK);
/*        case BSIM4_QB:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qb); 
            return(OK); */
        case BSIM4_CQB:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4cqb); 
            return(OK);
/*        case BSIM4_QG:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qg); 
            return(OK); */
        case BSIM4_CQG:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4cqg); 
            return(OK);
/*        case BSIM4_QD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qd); 
            return(OK); */
        case BSIM4_CQD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4cqd); 
            return(OK);
/*        case BSIM4_QS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qs); 
            return(OK); */
        case BSIM4_QB:
            value->rValue = here->BSIM4qbulk; 
            value->rValue *= here->BSIM4m;
            return(OK); 
        case BSIM4_QG:
            value->rValue = here->BSIM4qgate; 
            value->rValue *= here->BSIM4m;
            return(OK); 
        case BSIM4_QS:
            value->rValue = here->BSIM4qsrc; 
            value->rValue *= here->BSIM4m;
            return(OK); 
        case BSIM4_QD:
            value->rValue = here->BSIM4qdrn; 
            value->rValue *= here->BSIM4m;
            return(OK); 
        case BSIM4_QDEF:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qdef); 
            return(OK); 
        case BSIM4_GCRG:
            value->rValue = here->BSIM4gcrg;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_GTAU:
            value->rValue = here->BSIM4gtau;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CGGB:
            value->rValue = here->BSIM4cggb; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CGDB:
            value->rValue = here->BSIM4cgdb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CGSB:
            value->rValue = here->BSIM4cgsb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CDGB:
            value->rValue = here->BSIM4cdgb; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CDDB:
            value->rValue = here->BSIM4cddb; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CDSB:
            value->rValue = here->BSIM4cdsb; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBGB:
            value->rValue = here->BSIM4cbgb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBDB:
            value->rValue = here->BSIM4cbdb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBSB:
            value->rValue = here->BSIM4cbsb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CSGB:
            value->rValue = here->BSIM4csgb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CSDB:
            value->rValue = here->BSIM4csdb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CSSB:
            value->rValue = here->BSIM4cssb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CGBB:
            value->rValue = here->BSIM4cgbb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CDBB:
            value->rValue = here->BSIM4cdbb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CSBB:
            value->rValue = here->BSIM4csbb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CBBB:
            value->rValue = here->BSIM4cbbb;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CAPBD:
            value->rValue = here->BSIM4capbd; 
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_CAPBS:
            value->rValue = here->BSIM4capbs;
            value->rValue *= here->BSIM4m;
            return(OK);
        case BSIM4_VON:
            value->rValue = here->BSIM4von; 
            return(OK);
        case BSIM4_VDSAT:
            value->rValue = here->BSIM4vdsat; 
            return(OK);
        case BSIM4_QBS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qbs); 
            return(OK);
        case BSIM4_QBD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4qbd); 
            return(OK);
        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}

