/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4ask.c of BSIM4v7.0.
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
#include "bsim4v7def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4v7ask(
CKTcircuit *ckt,
GENinstance *inst,
int which,
IFvalue *value,
IFvalue *select)
{
BSIM4v7instance *here = (BSIM4v7instance*)inst;

    NG_IGNORE(select);

    switch(which) 
    {   case BSIM4v7_L:
            value->rValue = here->BSIM4v7l;
            return(OK);
        case BSIM4v7_W:
            value->rValue = here->BSIM4v7w;
            return(OK);
        case BSIM4v7_M:
            value->rValue = here->BSIM4v7m;
            return(OK);
        case BSIM4v7_NF:
            value->rValue = here->BSIM4v7nf;
            return(OK);
        case BSIM4v7_MIN:
            value->iValue = here->BSIM4v7min;
            return(OK);
        case BSIM4v7_AS:
            value->rValue = here->BSIM4v7sourceArea;
            return(OK);
        case BSIM4v7_AD:
            value->rValue = here->BSIM4v7drainArea;
            return(OK);
        case BSIM4v7_PS:
            value->rValue = here->BSIM4v7sourcePerimeter;
            return(OK);
        case BSIM4v7_PD:
            value->rValue = here->BSIM4v7drainPerimeter;
            return(OK);
        case BSIM4v7_NRS:
            value->rValue = here->BSIM4v7sourceSquares;
            return(OK);
        case BSIM4v7_NRD:
            value->rValue = here->BSIM4v7drainSquares;
            return(OK);
        case BSIM4v7_OFF:
            value->rValue = here->BSIM4v7off;
            return(OK);
        case BSIM4v7_SA:
            value->rValue = here->BSIM4v7sa ;
            return(OK);
        case BSIM4v7_SB:
            value->rValue = here->BSIM4v7sb ;
            return(OK);
        case BSIM4v7_SD:
            value->rValue = here->BSIM4v7sd ;
            return(OK);
        case BSIM4v7_SCA:
            value->rValue = here->BSIM4v7sca ;
            return(OK);
        case BSIM4v7_SCB:
            value->rValue = here->BSIM4v7scb ;
            return(OK);
        case BSIM4v7_SCC:
            value->rValue = here->BSIM4v7scc ;
            return(OK);
        case BSIM4v7_SC:
            value->rValue = here->BSIM4v7sc ;
            return(OK);

        case BSIM4v7_RBSB:
            value->rValue = here->BSIM4v7rbsb;
            return(OK);
        case BSIM4v7_RBDB:
            value->rValue = here->BSIM4v7rbdb;
            return(OK);
        case BSIM4v7_RBPB:
            value->rValue = here->BSIM4v7rbpb;
            return(OK);
        case BSIM4v7_RBPS:
            value->rValue = here->BSIM4v7rbps;
            return(OK);
        case BSIM4v7_RBPD:
            value->rValue = here->BSIM4v7rbpd;
            return(OK);
        case BSIM4v7_DELVTO:
            value->rValue = here->BSIM4v7delvto;
            return(OK);
        case BSIM4v7_XGW:
            value->rValue = here->BSIM4v7xgw;
            return(OK);
        case BSIM4v7_NGCON:
            value->rValue = here->BSIM4v7ngcon;
            return(OK);
        case BSIM4v7_TRNQSMOD:
            value->iValue = here->BSIM4v7trnqsMod;
            return(OK);
        case BSIM4v7_ACNQSMOD:
            value->iValue = here->BSIM4v7acnqsMod;
            return(OK);
        case BSIM4v7_RBODYMOD:
            value->iValue = here->BSIM4v7rbodyMod;
            return(OK);
        case BSIM4v7_RGATEMOD:
            value->iValue = here->BSIM4v7rgateMod;
            return(OK);
        case BSIM4v7_GEOMOD:
            value->iValue = here->BSIM4v7geoMod;
            return(OK);
        case BSIM4v7_RGEOMOD:
            value->iValue = here->BSIM4v7rgeoMod;
            return(OK);
        case BSIM4v7_IC_VDS:
            value->rValue = here->BSIM4v7icVDS;
            return(OK);
        case BSIM4v7_IC_VGS:
            value->rValue = here->BSIM4v7icVGS;
            return(OK);
        case BSIM4v7_IC_VBS:
            value->rValue = here->BSIM4v7icVBS;
            return(OK);
        case BSIM4v7_DNODE:
            value->iValue = here->BSIM4v7dNode;
            return(OK);
        case BSIM4v7_GNODEEXT:
            value->iValue = here->BSIM4v7gNodeExt;
            return(OK);
        case BSIM4v7_SNODE:
            value->iValue = here->BSIM4v7sNode;
            return(OK);
        case BSIM4v7_BNODE:
            value->iValue = here->BSIM4v7bNode;
            return(OK);
        case BSIM4v7_DNODEPRIME:
            value->iValue = here->BSIM4v7dNodePrime;
            return(OK);
        case BSIM4v7_GNODEPRIME:
            value->iValue = here->BSIM4v7gNodePrime;
            return(OK);
        case BSIM4v7_GNODEMID:
            value->iValue = here->BSIM4v7gNodeMid;
            return(OK);
        case BSIM4v7_SNODEPRIME:
            value->iValue = here->BSIM4v7sNodePrime;
            return(OK);
        case BSIM4v7_DBNODE:
            value->iValue = here->BSIM4v7dbNode;
            return(OK);
        case BSIM4v7_BNODEPRIME:
            value->iValue = here->BSIM4v7bNodePrime;
            return(OK);
        case BSIM4v7_SBNODE:
            value->iValue = here->BSIM4v7sbNode;
            return(OK);
        case BSIM4v7_SOURCECONDUCT:
            value->rValue = here->BSIM4v7sourceConductance;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_DRAINCONDUCT:
            value->rValue = here->BSIM4v7drainConductance;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_VBD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7vbd);
            return(OK);
        case BSIM4v7_VBS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7vbs);
            return(OK);
        case BSIM4v7_VGS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7vgs);
            return(OK);
        case BSIM4v7_VDS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7vds);
            return(OK);
        case BSIM4v7_CD:
            value->rValue = here->BSIM4v7cd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBS:
            value->rValue = here->BSIM4v7cbs; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBD:
            value->rValue = here->BSIM4v7cbd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CSUB:
            value->rValue = here->BSIM4v7csub; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_QINV:
            value->rValue = here-> BSIM4v7qinv; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGIDL:
            value->rValue = here->BSIM4v7Igidl; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGISL:
            value->rValue = here->BSIM4v7Igisl; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGS:
            value->rValue = here->BSIM4v7Igs; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGD:
            value->rValue = here->BSIM4v7Igd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGB:
            value->rValue = here->BSIM4v7Igb; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGCS:
            value->rValue = here->BSIM4v7Igcs; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_IGCD:
            value->rValue = here->BSIM4v7Igcd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GM:
            value->rValue = here->BSIM4v7gm; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GDS:
            value->rValue = here->BSIM4v7gds; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GMBS:
            value->rValue = here->BSIM4v7gmbs; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GBD:
            value->rValue = here->BSIM4v7gbd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GBS:
            value->rValue = here->BSIM4v7gbs; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
/*        case BSIM4v7_QB:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qb); 
            return(OK); */
        case BSIM4v7_CQB:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7cqb); 
            return(OK);
/*        case BSIM4v7_QG:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qg); 
            return(OK); */
        case BSIM4v7_CQG:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7cqg); 
            return(OK);
/*        case BSIM4v7_QD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qd); 
            return(OK); */
        case BSIM4v7_CQD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7cqd); 
            return(OK);
/*        case BSIM4v7_QS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qs); 
            return(OK); */
        case BSIM4v7_QB:
            value->rValue = here->BSIM4v7qbulk; 
            value->rValue *= here->BSIM4v7m;
            return(OK); 
        case BSIM4v7_QG:
            value->rValue = here->BSIM4v7qgate; 
            value->rValue *= here->BSIM4v7m;
            return(OK); 
        case BSIM4v7_QS:
            value->rValue = here->BSIM4v7qsrc; 
            value->rValue *= here->BSIM4v7m;
            return(OK); 
        case BSIM4v7_QD:
            value->rValue = here->BSIM4v7qdrn; 
            value->rValue *= here->BSIM4v7m;
            return(OK); 
        case BSIM4v7_QDEF:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qdef); 
            return(OK); 
        case BSIM4v7_GCRG:
            value->rValue = here->BSIM4v7gcrg;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_GTAU:
            value->rValue = here->BSIM4v7gtau;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CGGB:
            value->rValue = here->BSIM4v7cggb; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CGDB:
            value->rValue = here->BSIM4v7cgdb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CGSB:
            value->rValue = here->BSIM4v7cgsb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CDGB:
            value->rValue = here->BSIM4v7cdgb; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CDDB:
            value->rValue = here->BSIM4v7cddb; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CDSB:
            value->rValue = here->BSIM4v7cdsb; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBGB:
            value->rValue = here->BSIM4v7cbgb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBDB:
            value->rValue = here->BSIM4v7cbdb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBSB:
            value->rValue = here->BSIM4v7cbsb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CSGB:
            value->rValue = here->BSIM4v7csgb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CSDB:
            value->rValue = here->BSIM4v7csdb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CSSB:
            value->rValue = here->BSIM4v7cssb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CGBB:
            value->rValue = here->BSIM4v7cgbb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CDBB:
            value->rValue = here->BSIM4v7cdbb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CSBB:
            value->rValue = here->BSIM4v7csbb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CBBB:
            value->rValue = here->BSIM4v7cbbb;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CAPBD:
            value->rValue = here->BSIM4v7capbd; 
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_CAPBS:
            value->rValue = here->BSIM4v7capbs;
            value->rValue *= here->BSIM4v7m;
            return(OK);
        case BSIM4v7_VON:
            value->rValue = here->BSIM4v7von; 
            return(OK);
        case BSIM4v7_VDSAT:
            value->rValue = here->BSIM4v7vdsat; 
            return(OK);
        case BSIM4v7_QBS:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qbs); 
            return(OK);
        case BSIM4v7_QBD:
            value->rValue = *(ckt->CKTstate0 + here->BSIM4v7qbd); 
            return(OK);
        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}

