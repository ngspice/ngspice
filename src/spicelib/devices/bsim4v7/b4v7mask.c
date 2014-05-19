/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4mask.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 11/15/2002.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, Mohan Dunga, 07/29/2005.
 * Modified by Mohan Dunga, 12/13/2006
 * Modified by Mohan Dunga, Wenwei Yang, 07/31/2008.
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
 **********/


#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "bsim4v7def.h."
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4v7mAsk(
CKTcircuit *ckt,
GENmodel *inst,
int which,
IFvalue *value)
{
    BSIM4v7model *model = (BSIM4v7model *)inst;

    NG_IGNORE(ckt);

    switch(which) 
    {   case BSIM4v7_MOD_MOBMOD :
            value->iValue = model->BSIM4v7mobMod; 
            return(OK);
        case BSIM4v7_MOD_PARAMCHK :
            value->iValue = model->BSIM4v7paramChk; 
            return(OK);
        case BSIM4v7_MOD_BINUNIT :
            value->iValue = model->BSIM4v7binUnit; 
            return(OK);
        case BSIM4v7_MOD_CVCHARGEMOD :
            value->iValue = model->BSIM4v7cvchargeMod; 
            return(OK);
        case BSIM4v7_MOD_CAPMOD :
            value->iValue = model->BSIM4v7capMod; 
            return(OK);
        case BSIM4v7_MOD_DIOMOD :
            value->iValue = model->BSIM4v7dioMod;
            return(OK);
        case BSIM4v7_MOD_TRNQSMOD :
            value->iValue = model->BSIM4v7trnqsMod;
            return(OK);
        case BSIM4v7_MOD_ACNQSMOD :
            value->iValue = model->BSIM4v7acnqsMod;
            return(OK);
        case BSIM4v7_MOD_FNOIMOD :
            value->iValue = model->BSIM4v7fnoiMod; 
            return(OK);
        case BSIM4v7_MOD_TNOIMOD :
            value->iValue = model->BSIM4v7tnoiMod;
            return(OK);
        case BSIM4v7_MOD_RDSMOD :
            value->iValue = model->BSIM4v7rdsMod;
            return(OK);
        case BSIM4v7_MOD_RBODYMOD :
            value->iValue = model->BSIM4v7rbodyMod;
            return(OK);
        case BSIM4v7_MOD_RGATEMOD :
            value->iValue = model->BSIM4v7rgateMod;
            return(OK);
        case BSIM4v7_MOD_PERMOD :
            value->iValue = model->BSIM4v7perMod;
            return(OK);
        case BSIM4v7_MOD_GEOMOD :
            value->iValue = model->BSIM4v7geoMod;
            return(OK);
        case BSIM4v7_MOD_RGEOMOD :
            value->iValue = model->BSIM4v7rgeoMod;
            return(OK);
        case BSIM4v7_MOD_MTRLMOD :
            value->iValue = model->BSIM4v7mtrlMod;
            return(OK);
	case BSIM4v7_MOD_GIDLMOD :		/* v4.7 New GIDL/GISL*/
            value->iValue = model->BSIM4v7gidlMod;
            return(OK);
        case BSIM4v7_MOD_MTRLCOMPATMOD :
	    value->iValue = model->BSIM4v7mtrlCompatMod;
	    return(OK);
        case BSIM4v7_MOD_IGCMOD :
            value->iValue = model->BSIM4v7igcMod;
            return(OK);
        case BSIM4v7_MOD_IGBMOD :
            value->iValue = model->BSIM4v7igbMod;
            return(OK);
        case  BSIM4v7_MOD_TEMPMOD :
            value->iValue = model->BSIM4v7tempMod;
            return(OK);

        case  BSIM4v7_MOD_VERSION :
          value->sValue = model->BSIM4v7version;
            return(OK);
        case  BSIM4v7_MOD_TOXREF :
          value->rValue = model->BSIM4v7toxref;
          return(OK);
        case  BSIM4v7_MOD_EOT :
          value->rValue = model->BSIM4v7eot;
            return(OK);
        case  BSIM4v7_MOD_VDDEOT :
          value->rValue = model->BSIM4v7vddeot;
            return(OK);
		case  BSIM4v7_MOD_TEMPEOT :
          value->rValue = model->BSIM4v7tempeot;
            return(OK);
		case  BSIM4v7_MOD_LEFFEOT :
          value->rValue = model->BSIM4v7leffeot;
            return(OK);
		case  BSIM4v7_MOD_WEFFEOT :
          value->rValue = model->BSIM4v7weffeot;
            return(OK);
        case  BSIM4v7_MOD_ADOS :
          value->rValue = model->BSIM4v7ados;
            return(OK);
        case  BSIM4v7_MOD_BDOS :
          value->rValue = model->BSIM4v7bdos;
            return(OK);
        case  BSIM4v7_MOD_TOXE :
          value->rValue = model->BSIM4v7toxe;
            return(OK);
        case  BSIM4v7_MOD_TOXP :
          value->rValue = model->BSIM4v7toxp;
            return(OK);
        case  BSIM4v7_MOD_TOXM :
          value->rValue = model->BSIM4v7toxm;
            return(OK);
        case  BSIM4v7_MOD_DTOX :
          value->rValue = model->BSIM4v7dtox;
            return(OK);
        case  BSIM4v7_MOD_EPSROX :
          value->rValue = model->BSIM4v7epsrox;
            return(OK);
        case  BSIM4v7_MOD_CDSC :
          value->rValue = model->BSIM4v7cdsc;
            return(OK);
        case  BSIM4v7_MOD_CDSCB :
          value->rValue = model->BSIM4v7cdscb;
            return(OK);

        case  BSIM4v7_MOD_CDSCD :
          value->rValue = model->BSIM4v7cdscd;
            return(OK);

        case  BSIM4v7_MOD_CIT :
          value->rValue = model->BSIM4v7cit;
            return(OK);
        case  BSIM4v7_MOD_NFACTOR :
          value->rValue = model->BSIM4v7nfactor;
            return(OK);
        case BSIM4v7_MOD_XJ:
            value->rValue = model->BSIM4v7xj;
            return(OK);
        case BSIM4v7_MOD_VSAT:
            value->rValue = model->BSIM4v7vsat;
            return(OK);
        case BSIM4v7_MOD_VTL:
            value->rValue = model->BSIM4v7vtl;
            return(OK);
        case BSIM4v7_MOD_XN:
            value->rValue = model->BSIM4v7xn;
            return(OK);
        case BSIM4v7_MOD_LC:
            value->rValue = model->BSIM4v7lc;
            return(OK);
        case BSIM4v7_MOD_LAMBDA:
            value->rValue = model->BSIM4v7lambda;
            return(OK);
        case BSIM4v7_MOD_AT:
            value->rValue = model->BSIM4v7at;
            return(OK);
        case BSIM4v7_MOD_A0:
            value->rValue = model->BSIM4v7a0;
            return(OK);

        case BSIM4v7_MOD_AGS:
            value->rValue = model->BSIM4v7ags;
            return(OK);

        case BSIM4v7_MOD_A1:
            value->rValue = model->BSIM4v7a1;
            return(OK);
        case BSIM4v7_MOD_A2:
            value->rValue = model->BSIM4v7a2;
            return(OK);
        case BSIM4v7_MOD_KETA:
            value->rValue = model->BSIM4v7keta;
            return(OK);   
        case BSIM4v7_MOD_NSUB:
            value->rValue = model->BSIM4v7nsub;
            return(OK);
        case BSIM4v7_MOD_PHIG:
	    value->rValue = model->BSIM4v7phig;
	    return(OK);
        case BSIM4v7_MOD_EPSRGATE:
	    value->rValue = model->BSIM4v7epsrgate;
	    return(OK);
        case BSIM4v7_MOD_EASUB:
            value->rValue = model->BSIM4v7easub;
            return(OK);
        case BSIM4v7_MOD_EPSRSUB:
            value->rValue = model->BSIM4v7epsrsub;
            return(OK);
        case BSIM4v7_MOD_NI0SUB:
            value->rValue = model->BSIM4v7ni0sub;
            return(OK);
        case BSIM4v7_MOD_BG0SUB:
            value->rValue = model->BSIM4v7bg0sub;
            return(OK);
        case BSIM4v7_MOD_TBGASUB:
            value->rValue = model->BSIM4v7tbgasub;
            return(OK);
        case BSIM4v7_MOD_TBGBSUB:
            value->rValue = model->BSIM4v7tbgbsub;
            return(OK);
        case BSIM4v7_MOD_NDEP:
            value->rValue = model->BSIM4v7ndep;
            return(OK);
        case BSIM4v7_MOD_NSD:
            value->rValue = model->BSIM4v7nsd;
            return(OK);
        case BSIM4v7_MOD_NGATE:
            value->rValue = model->BSIM4v7ngate;
            return(OK);
        case BSIM4v7_MOD_GAMMA1:
            value->rValue = model->BSIM4v7gamma1;
            return(OK);
        case BSIM4v7_MOD_GAMMA2:
            value->rValue = model->BSIM4v7gamma2;
            return(OK);
        case BSIM4v7_MOD_VBX:
            value->rValue = model->BSIM4v7vbx;
            return(OK);
        case BSIM4v7_MOD_VBM:
            value->rValue = model->BSIM4v7vbm;
            return(OK);
        case BSIM4v7_MOD_XT:
            value->rValue = model->BSIM4v7xt;
            return(OK);
        case  BSIM4v7_MOD_K1:
          value->rValue = model->BSIM4v7k1;
            return(OK);
        case  BSIM4v7_MOD_KT1:
          value->rValue = model->BSIM4v7kt1;
            return(OK);
        case  BSIM4v7_MOD_KT1L:
          value->rValue = model->BSIM4v7kt1l;
            return(OK);
        case  BSIM4v7_MOD_KT2 :
          value->rValue = model->BSIM4v7kt2;
            return(OK);
        case  BSIM4v7_MOD_K2 :
          value->rValue = model->BSIM4v7k2;
            return(OK);
        case  BSIM4v7_MOD_K3:
          value->rValue = model->BSIM4v7k3;
            return(OK);
        case  BSIM4v7_MOD_K3B:
          value->rValue = model->BSIM4v7k3b;
            return(OK);
        case  BSIM4v7_MOD_W0:
          value->rValue = model->BSIM4v7w0;
            return(OK);
        case  BSIM4v7_MOD_LPE0:
          value->rValue = model->BSIM4v7lpe0;
            return(OK);
        case  BSIM4v7_MOD_LPEB:
          value->rValue = model->BSIM4v7lpeb;
            return(OK);
        case  BSIM4v7_MOD_DVTP0:
          value->rValue = model->BSIM4v7dvtp0;
            return(OK);
        case  BSIM4v7_MOD_DVTP1:
          value->rValue = model->BSIM4v7dvtp1;
            return(OK);
        case  BSIM4v7_MOD_DVTP2:
          value->rValue = model->BSIM4v7dvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4v7_MOD_DVTP3:
          value->rValue = model->BSIM4v7dvtp3;
            return(OK);
        case  BSIM4v7_MOD_DVTP4:
          value->rValue = model->BSIM4v7dvtp4;
            return(OK);
        case  BSIM4v7_MOD_DVTP5:
          value->rValue = model->BSIM4v7dvtp5;
            return(OK);
        case  BSIM4v7_MOD_DVT0 :                
          value->rValue = model->BSIM4v7dvt0;
            return(OK);
        case  BSIM4v7_MOD_DVT1 :             
          value->rValue = model->BSIM4v7dvt1;
            return(OK);
        case  BSIM4v7_MOD_DVT2 :             
          value->rValue = model->BSIM4v7dvt2;
            return(OK);
        case  BSIM4v7_MOD_DVT0W :                
          value->rValue = model->BSIM4v7dvt0w;
            return(OK);
        case  BSIM4v7_MOD_DVT1W :             
          value->rValue = model->BSIM4v7dvt1w;
            return(OK);
        case  BSIM4v7_MOD_DVT2W :             
          value->rValue = model->BSIM4v7dvt2w;
            return(OK);
        case  BSIM4v7_MOD_DROUT :           
          value->rValue = model->BSIM4v7drout;
            return(OK);
        case  BSIM4v7_MOD_DSUB :           
          value->rValue = model->BSIM4v7dsub;
            return(OK);
        case BSIM4v7_MOD_VTH0:
            value->rValue = model->BSIM4v7vth0; 
            return(OK);
        case BSIM4v7_MOD_EU:
            value->rValue = model->BSIM4v7eu;
            return(OK);
		 case BSIM4v7_MOD_UCS:
            value->rValue = model->BSIM4v7ucs;
            return(OK);
        case BSIM4v7_MOD_UA:
            value->rValue = model->BSIM4v7ua; 
            return(OK);
        case BSIM4v7_MOD_UA1:
            value->rValue = model->BSIM4v7ua1; 
            return(OK);
        case BSIM4v7_MOD_UB:
            value->rValue = model->BSIM4v7ub;  
            return(OK);
        case BSIM4v7_MOD_UB1:
            value->rValue = model->BSIM4v7ub1;  
            return(OK);
        case BSIM4v7_MOD_UC:
            value->rValue = model->BSIM4v7uc; 
            return(OK);
        case BSIM4v7_MOD_UC1:
            value->rValue = model->BSIM4v7uc1; 
            return(OK);
        case BSIM4v7_MOD_UD:
            value->rValue = model->BSIM4v7ud; 
            return(OK);
        case BSIM4v7_MOD_UD1:
            value->rValue = model->BSIM4v7ud1; 
            return(OK);
        case BSIM4v7_MOD_UP:
            value->rValue = model->BSIM4v7up; 
            return(OK);
        case BSIM4v7_MOD_LP:
            value->rValue = model->BSIM4v7lp; 
            return(OK);
        case BSIM4v7_MOD_U0:
            value->rValue = model->BSIM4v7u0;
            return(OK);
        case BSIM4v7_MOD_UTE:
            value->rValue = model->BSIM4v7ute;
            return(OK);
		 case BSIM4v7_MOD_UCSTE:
            value->rValue = model->BSIM4v7ucste;
            return(OK);
        case BSIM4v7_MOD_VOFF:
            value->rValue = model->BSIM4v7voff;
            return(OK);
        case BSIM4v7_MOD_TVOFF:
            value->rValue = model->BSIM4v7tvoff;
            return(OK);
        case BSIM4v7_MOD_TNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7tnfactor;
            return(OK);
        case BSIM4v7_MOD_TETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7teta0;
            return(OK);
        case BSIM4v7_MOD_TVOFFCV:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7tvoffcv;
            return(OK);
        case BSIM4v7_MOD_VFBSDOFF:
            value->rValue = model->BSIM4v7vfbsdoff;
            return(OK);
        case BSIM4v7_MOD_TVFBSDOFF:
            value->rValue = model->BSIM4v7tvfbsdoff;
            return(OK);
        case BSIM4v7_MOD_VOFFL:
            value->rValue = model->BSIM4v7voffl;
            return(OK);
        case BSIM4v7_MOD_VOFFCVL:
            value->rValue = model->BSIM4v7voffcvl;
            return(OK);
        case BSIM4v7_MOD_MINV:
            value->rValue = model->BSIM4v7minv;
            return(OK);
        case BSIM4v7_MOD_MINVCV:
            value->rValue = model->BSIM4v7minvcv;
            return(OK);
        case BSIM4v7_MOD_FPROUT:
            value->rValue = model->BSIM4v7fprout;
            return(OK);
        case BSIM4v7_MOD_PDITS:
            value->rValue = model->BSIM4v7pdits;
            return(OK);
        case BSIM4v7_MOD_PDITSD:
            value->rValue = model->BSIM4v7pditsd;
            return(OK);
        case BSIM4v7_MOD_PDITSL:
            value->rValue = model->BSIM4v7pditsl;
            return(OK);
        case BSIM4v7_MOD_DELTA:
            value->rValue = model->BSIM4v7delta;
            return(OK);
        case BSIM4v7_MOD_RDSW:
            value->rValue = model->BSIM4v7rdsw; 
            return(OK);
        case BSIM4v7_MOD_RDSWMIN:
            value->rValue = model->BSIM4v7rdswmin;
            return(OK);
        case BSIM4v7_MOD_RDWMIN:
            value->rValue = model->BSIM4v7rdwmin;
            return(OK);
        case BSIM4v7_MOD_RSWMIN:
            value->rValue = model->BSIM4v7rswmin;
            return(OK);
        case BSIM4v7_MOD_RDW:
            value->rValue = model->BSIM4v7rdw;
            return(OK);
        case BSIM4v7_MOD_RSW:
            value->rValue = model->BSIM4v7rsw;
            return(OK);
        case BSIM4v7_MOD_PRWG:
            value->rValue = model->BSIM4v7prwg; 
            return(OK);             
        case BSIM4v7_MOD_PRWB:
            value->rValue = model->BSIM4v7prwb; 
            return(OK);             
        case BSIM4v7_MOD_PRT:
            value->rValue = model->BSIM4v7prt; 
            return(OK);              
        case BSIM4v7_MOD_ETA0:
            value->rValue = model->BSIM4v7eta0; 
            return(OK);               
        case BSIM4v7_MOD_ETAB:
            value->rValue = model->BSIM4v7etab; 
            return(OK);               
        case BSIM4v7_MOD_PCLM:
            value->rValue = model->BSIM4v7pclm; 
            return(OK);               
        case BSIM4v7_MOD_PDIBL1:
            value->rValue = model->BSIM4v7pdibl1; 
            return(OK);               
        case BSIM4v7_MOD_PDIBL2:
            value->rValue = model->BSIM4v7pdibl2; 
            return(OK);               
        case BSIM4v7_MOD_PDIBLB:
            value->rValue = model->BSIM4v7pdiblb; 
            return(OK);               
        case BSIM4v7_MOD_PSCBE1:
            value->rValue = model->BSIM4v7pscbe1; 
            return(OK);               
        case BSIM4v7_MOD_PSCBE2:
            value->rValue = model->BSIM4v7pscbe2; 
            return(OK);               
        case BSIM4v7_MOD_PVAG:
            value->rValue = model->BSIM4v7pvag; 
            return(OK);               
        case BSIM4v7_MOD_WR:
            value->rValue = model->BSIM4v7wr;
            return(OK);
        case BSIM4v7_MOD_DWG:
            value->rValue = model->BSIM4v7dwg;
            return(OK);
        case BSIM4v7_MOD_DWB:
            value->rValue = model->BSIM4v7dwb;
            return(OK);
        case BSIM4v7_MOD_B0:
            value->rValue = model->BSIM4v7b0;
            return(OK);
        case BSIM4v7_MOD_B1:
            value->rValue = model->BSIM4v7b1;
            return(OK);
        case BSIM4v7_MOD_ALPHA0:
            value->rValue = model->BSIM4v7alpha0;
            return(OK);
        case BSIM4v7_MOD_ALPHA1:
            value->rValue = model->BSIM4v7alpha1;
            return(OK);
        case BSIM4v7_MOD_BETA0:
            value->rValue = model->BSIM4v7beta0;
            return(OK);
        case BSIM4v7_MOD_AGIDL:
            value->rValue = model->BSIM4v7agidl;
            return(OK);
        case BSIM4v7_MOD_BGIDL:
            value->rValue = model->BSIM4v7bgidl;
            return(OK);
        case BSIM4v7_MOD_CGIDL:
            value->rValue = model->BSIM4v7cgidl;
            return(OK);
        case BSIM4v7_MOD_EGIDL:
            value->rValue = model->BSIM4v7egidl;
            return(OK);
 	case BSIM4v7_MOD_FGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7fgidl;
            return(OK);
 	case BSIM4v7_MOD_KGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7kgidl;
            return(OK);
 	case BSIM4v7_MOD_RGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7rgidl;
            return(OK);
        case BSIM4v7_MOD_AGISL:
            value->rValue = model->BSIM4v7agisl;
            return(OK);
        case BSIM4v7_MOD_BGISL:
            value->rValue = model->BSIM4v7bgisl;
            return(OK);
        case BSIM4v7_MOD_CGISL:
            value->rValue = model->BSIM4v7cgisl;
            return(OK);
        case BSIM4v7_MOD_EGISL:
            value->rValue = model->BSIM4v7egisl;
            return(OK);
 	case BSIM4v7_MOD_FGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7fgisl;
            return(OK);
 	case BSIM4v7_MOD_KGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7kgisl;
            return(OK);
 	case BSIM4v7_MOD_RGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7rgisl;
            return(OK);
        case BSIM4v7_MOD_AIGC:
            value->rValue = model->BSIM4v7aigc;
            return(OK);
        case BSIM4v7_MOD_BIGC:
            value->rValue = model->BSIM4v7bigc;
            return(OK);
        case BSIM4v7_MOD_CIGC:
            value->rValue = model->BSIM4v7cigc;
            return(OK);
        case BSIM4v7_MOD_AIGSD:
            value->rValue = model->BSIM4v7aigsd;
            return(OK);
        case BSIM4v7_MOD_BIGSD:
            value->rValue = model->BSIM4v7bigsd;
            return(OK);
        case BSIM4v7_MOD_CIGSD:
            value->rValue = model->BSIM4v7cigsd;
            return(OK);
        case BSIM4v7_MOD_AIGS:
            value->rValue = model->BSIM4v7aigs;
            return(OK);
        case BSIM4v7_MOD_BIGS:
            value->rValue = model->BSIM4v7bigs;
            return(OK);
        case BSIM4v7_MOD_CIGS:
            value->rValue = model->BSIM4v7cigs;
            return(OK);
        case BSIM4v7_MOD_AIGD:
            value->rValue = model->BSIM4v7aigd;
            return(OK);
        case BSIM4v7_MOD_BIGD:
            value->rValue = model->BSIM4v7bigd;
            return(OK);
        case BSIM4v7_MOD_CIGD:
            value->rValue = model->BSIM4v7cigd;
            return(OK);
        case BSIM4v7_MOD_AIGBACC:
            value->rValue = model->BSIM4v7aigbacc;
            return(OK);
        case BSIM4v7_MOD_BIGBACC:
            value->rValue = model->BSIM4v7bigbacc;
            return(OK);
        case BSIM4v7_MOD_CIGBACC:
            value->rValue = model->BSIM4v7cigbacc;
            return(OK);
        case BSIM4v7_MOD_AIGBINV:
            value->rValue = model->BSIM4v7aigbinv;
            return(OK);
        case BSIM4v7_MOD_BIGBINV:
            value->rValue = model->BSIM4v7bigbinv;
            return(OK);
        case BSIM4v7_MOD_CIGBINV:
            value->rValue = model->BSIM4v7cigbinv;
            return(OK);
        case BSIM4v7_MOD_NIGC:
            value->rValue = model->BSIM4v7nigc;
            return(OK);
        case BSIM4v7_MOD_NIGBACC:
            value->rValue = model->BSIM4v7nigbacc;
            return(OK);
        case BSIM4v7_MOD_NIGBINV:
            value->rValue = model->BSIM4v7nigbinv;
            return(OK);
        case BSIM4v7_MOD_NTOX:
            value->rValue = model->BSIM4v7ntox;
            return(OK);
        case BSIM4v7_MOD_EIGBINV:
            value->rValue = model->BSIM4v7eigbinv;
            return(OK);
        case BSIM4v7_MOD_PIGCD:
            value->rValue = model->BSIM4v7pigcd;
            return(OK);
        case BSIM4v7_MOD_POXEDGE:
            value->rValue = model->BSIM4v7poxedge;
            return(OK);
        case BSIM4v7_MOD_PHIN:
            value->rValue = model->BSIM4v7phin;
            return(OK);
        case BSIM4v7_MOD_XRCRG1:
            value->rValue = model->BSIM4v7xrcrg1;
            return(OK);
        case BSIM4v7_MOD_XRCRG2:
            value->rValue = model->BSIM4v7xrcrg2;
            return(OK);
        case BSIM4v7_MOD_TNOIA:
            value->rValue = model->BSIM4v7tnoia;
            return(OK);
        case BSIM4v7_MOD_TNOIB:
            value->rValue = model->BSIM4v7tnoib;
            return(OK);
        case BSIM4v7_MOD_TNOIC:
            value->rValue = model->BSIM4v7tnoic;
            return(OK);
        case BSIM4v7_MOD_RNOIA:
            value->rValue = model->BSIM4v7rnoia;
            return(OK);
        case BSIM4v7_MOD_RNOIB:
            value->rValue = model->BSIM4v7rnoib;
            return(OK);
        case BSIM4v7_MOD_RNOIC:
            value->rValue = model->BSIM4v7rnoic;
            return(OK);
        case BSIM4v7_MOD_NTNOI:
            value->rValue = model->BSIM4v7ntnoi;
            return(OK);
        case BSIM4v7_MOD_IJTHDFWD:
            value->rValue = model->BSIM4v7ijthdfwd;
            return(OK);
        case BSIM4v7_MOD_IJTHSFWD:
            value->rValue = model->BSIM4v7ijthsfwd;
            return(OK);
        case BSIM4v7_MOD_IJTHDREV:
            value->rValue = model->BSIM4v7ijthdrev;
            return(OK);
        case BSIM4v7_MOD_IJTHSREV:
            value->rValue = model->BSIM4v7ijthsrev;
            return(OK);
        case BSIM4v7_MOD_XJBVD:
            value->rValue = model->BSIM4v7xjbvd;
            return(OK);
        case BSIM4v7_MOD_XJBVS:
            value->rValue = model->BSIM4v7xjbvs;
            return(OK);
        case BSIM4v7_MOD_BVD:
            value->rValue = model->BSIM4v7bvd;
            return(OK);
        case BSIM4v7_MOD_BVS:
            value->rValue = model->BSIM4v7bvs;
            return(OK);
        case BSIM4v7_MOD_VFB:
            value->rValue = model->BSIM4v7vfb;
            return(OK);

        case BSIM4v7_MOD_JTSS:
            value->rValue = model->BSIM4v7jtss;
            return(OK);
        case BSIM4v7_MOD_JTSD:
            value->rValue = model->BSIM4v7jtsd;
            return(OK);
        case BSIM4v7_MOD_JTSSWS:
            value->rValue = model->BSIM4v7jtssws;
            return(OK);
        case BSIM4v7_MOD_JTSSWD:
            value->rValue = model->BSIM4v7jtsswd;
            return(OK);
        case BSIM4v7_MOD_JTSSWGS:
            value->rValue = model->BSIM4v7jtsswgs;
            return(OK);
        case BSIM4v7_MOD_JTSSWGD:
            value->rValue = model->BSIM4v7jtsswgd;
            return(OK);
		case BSIM4v7_MOD_JTWEFF:
		    value->rValue = model->BSIM4v7jtweff;
			return(OK);
        case BSIM4v7_MOD_NJTS:
            value->rValue = model->BSIM4v7njts;
            return(OK);
        case BSIM4v7_MOD_NJTSSW:
            value->rValue = model->BSIM4v7njtssw;
            return(OK);
        case BSIM4v7_MOD_NJTSSWG:
            value->rValue = model->BSIM4v7njtsswg;
            return(OK);
        case BSIM4v7_MOD_NJTSD:
            value->rValue = model->BSIM4v7njtsd;
            return(OK);
        case BSIM4v7_MOD_NJTSSWD:
            value->rValue = model->BSIM4v7njtsswd;
            return(OK);
        case BSIM4v7_MOD_NJTSSWGD:
            value->rValue = model->BSIM4v7njtsswgd;
            return(OK);
        case BSIM4v7_MOD_XTSS:
            value->rValue = model->BSIM4v7xtss;
            return(OK);
        case BSIM4v7_MOD_XTSD:
            value->rValue = model->BSIM4v7xtsd;
            return(OK);
        case BSIM4v7_MOD_XTSSWS:
            value->rValue = model->BSIM4v7xtssws;
            return(OK);
        case BSIM4v7_MOD_XTSSWD:
            value->rValue = model->BSIM4v7xtsswd;
            return(OK);
        case BSIM4v7_MOD_XTSSWGS:
            value->rValue = model->BSIM4v7xtsswgs;
            return(OK);
        case BSIM4v7_MOD_XTSSWGD:
            value->rValue = model->BSIM4v7xtsswgd;
            return(OK);
        case BSIM4v7_MOD_TNJTS:
            value->rValue = model->BSIM4v7tnjts;
            return(OK);
        case BSIM4v7_MOD_TNJTSSW:
            value->rValue = model->BSIM4v7tnjtssw;
            return(OK);
        case BSIM4v7_MOD_TNJTSSWG:
            value->rValue = model->BSIM4v7tnjtsswg;
            return(OK);
        case BSIM4v7_MOD_TNJTSD:
            value->rValue = model->BSIM4v7tnjtsd;
            return(OK);
        case BSIM4v7_MOD_TNJTSSWD:
            value->rValue = model->BSIM4v7tnjtsswd;
            return(OK);
        case BSIM4v7_MOD_TNJTSSWGD:
            value->rValue = model->BSIM4v7tnjtsswgd;
            return(OK);
        case BSIM4v7_MOD_VTSS:
            value->rValue = model->BSIM4v7vtss;
            return(OK);
        case BSIM4v7_MOD_VTSD:
            value->rValue = model->BSIM4v7vtsd;
            return(OK);
        case BSIM4v7_MOD_VTSSWS:
            value->rValue = model->BSIM4v7vtssws;
            return(OK);
        case BSIM4v7_MOD_VTSSWD:
            value->rValue = model->BSIM4v7vtsswd;
            return(OK);
        case BSIM4v7_MOD_VTSSWGS:
            value->rValue = model->BSIM4v7vtsswgs;
            return(OK);
        case BSIM4v7_MOD_VTSSWGD:
            value->rValue = model->BSIM4v7vtsswgd;
            return(OK);

        case BSIM4v7_MOD_GBMIN:
            value->rValue = model->BSIM4v7gbmin;
            return(OK);
        case BSIM4v7_MOD_RBDB:
            value->rValue = model->BSIM4v7rbdb;
            return(OK);
        case BSIM4v7_MOD_RBPB:
            value->rValue = model->BSIM4v7rbpb;
            return(OK);
        case BSIM4v7_MOD_RBSB:
            value->rValue = model->BSIM4v7rbsb;
            return(OK);
        case BSIM4v7_MOD_RBPS:
            value->rValue = model->BSIM4v7rbps;
            return(OK);
        case BSIM4v7_MOD_RBPD:
            value->rValue = model->BSIM4v7rbpd;
            return(OK);

        case BSIM4v7_MOD_RBPS0:
            value->rValue = model->BSIM4v7rbps0;
            return(OK);
        case BSIM4v7_MOD_RBPSL:
            value->rValue = model->BSIM4v7rbpsl;
            return(OK);
        case BSIM4v7_MOD_RBPSW:
            value->rValue = model->BSIM4v7rbpsw;
            return(OK);
        case BSIM4v7_MOD_RBPSNF:
            value->rValue = model->BSIM4v7rbpsnf;
            return(OK);
        case BSIM4v7_MOD_RBPD0:
            value->rValue = model->BSIM4v7rbpd0;
            return(OK);
        case BSIM4v7_MOD_RBPDL:
            value->rValue = model->BSIM4v7rbpdl;
            return(OK);
        case BSIM4v7_MOD_RBPDW:
            value->rValue = model->BSIM4v7rbpdw;
            return(OK);
        case BSIM4v7_MOD_RBPDNF:
            value->rValue = model->BSIM4v7rbpdnf;
            return(OK);
        case BSIM4v7_MOD_RBPBX0:
            value->rValue = model->BSIM4v7rbpbx0;
            return(OK);
        case BSIM4v7_MOD_RBPBXL:
            value->rValue = model->BSIM4v7rbpbxl;
            return(OK);
        case BSIM4v7_MOD_RBPBXW:
            value->rValue = model->BSIM4v7rbpbxw;
            return(OK);
        case BSIM4v7_MOD_RBPBXNF:
            value->rValue = model->BSIM4v7rbpbxnf;
            return(OK);
        case BSIM4v7_MOD_RBPBY0:
            value->rValue = model->BSIM4v7rbpby0;
            return(OK);
        case BSIM4v7_MOD_RBPBYL:
            value->rValue = model->BSIM4v7rbpbyl;
            return(OK);
        case BSIM4v7_MOD_RBPBYW:
            value->rValue = model->BSIM4v7rbpbyw;
            return(OK);
        case BSIM4v7_MOD_RBPBYNF:
            value->rValue = model->BSIM4v7rbpbynf;
            return(OK);

        case BSIM4v7_MOD_RBSBX0:
            value->rValue = model->BSIM4v7rbsbx0;
            return(OK);
        case BSIM4v7_MOD_RBSBY0:
            value->rValue = model->BSIM4v7rbsby0;
            return(OK);
        case BSIM4v7_MOD_RBDBX0:
            value->rValue = model->BSIM4v7rbdbx0;
            return(OK);
        case BSIM4v7_MOD_RBDBY0:
            value->rValue = model->BSIM4v7rbdby0;
            return(OK);
        case BSIM4v7_MOD_RBSDBXL:
            value->rValue = model->BSIM4v7rbsdbxl;
            return(OK);
        case BSIM4v7_MOD_RBSDBXW:
            value->rValue = model->BSIM4v7rbsdbxw;
            return(OK);
        case BSIM4v7_MOD_RBSDBXNF:
            value->rValue = model->BSIM4v7rbsdbxnf;
            return(OK);
        case BSIM4v7_MOD_RBSDBYL:
            value->rValue = model->BSIM4v7rbsdbyl;
            return(OK);
        case BSIM4v7_MOD_RBSDBYW:
            value->rValue = model->BSIM4v7rbsdbyw;
            return(OK);
        case BSIM4v7_MOD_RBSDBYNF:
            value->rValue = model->BSIM4v7rbsdbynf;
            return(OK);


        case BSIM4v7_MOD_CGSL:
            value->rValue = model->BSIM4v7cgsl;
            return(OK);
        case BSIM4v7_MOD_CGDL:
            value->rValue = model->BSIM4v7cgdl;
            return(OK);
        case BSIM4v7_MOD_CKAPPAS:
            value->rValue = model->BSIM4v7ckappas;
            return(OK);
        case BSIM4v7_MOD_CKAPPAD:
            value->rValue = model->BSIM4v7ckappad;
            return(OK);
        case BSIM4v7_MOD_CF:
            value->rValue = model->BSIM4v7cf;
            return(OK);
        case BSIM4v7_MOD_CLC:
            value->rValue = model->BSIM4v7clc;
            return(OK);
        case BSIM4v7_MOD_CLE:
            value->rValue = model->BSIM4v7cle;
            return(OK);
        case BSIM4v7_MOD_DWC:
            value->rValue = model->BSIM4v7dwc;
            return(OK);
        case BSIM4v7_MOD_DLC:
            value->rValue = model->BSIM4v7dlc;
            return(OK);
        case BSIM4v7_MOD_XW:
            value->rValue = model->BSIM4v7xw;
            return(OK);
        case BSIM4v7_MOD_XL:
            value->rValue = model->BSIM4v7xl;
            return(OK);
        case BSIM4v7_MOD_DLCIG:
            value->rValue = model->BSIM4v7dlcig;
            return(OK);
        case BSIM4v7_MOD_DLCIGD:
            value->rValue = model->BSIM4v7dlcigd;
            return(OK);
        case BSIM4v7_MOD_DWJ:
            value->rValue = model->BSIM4v7dwj;
            return(OK);
        case BSIM4v7_MOD_VFBCV:
            value->rValue = model->BSIM4v7vfbcv; 
            return(OK);
        case BSIM4v7_MOD_ACDE:
            value->rValue = model->BSIM4v7acde;
            return(OK);
        case BSIM4v7_MOD_MOIN:
            value->rValue = model->BSIM4v7moin;
            return(OK);
        case BSIM4v7_MOD_NOFF:
            value->rValue = model->BSIM4v7noff;
            return(OK);
        case BSIM4v7_MOD_VOFFCV:
            value->rValue = model->BSIM4v7voffcv;
            return(OK);
        case BSIM4v7_MOD_DMCG:
            value->rValue = model->BSIM4v7dmcg;
            return(OK);
        case BSIM4v7_MOD_DMCI:
            value->rValue = model->BSIM4v7dmci;
            return(OK);
        case BSIM4v7_MOD_DMDG:
            value->rValue = model->BSIM4v7dmdg;
            return(OK);
        case BSIM4v7_MOD_DMCGT:
            value->rValue = model->BSIM4v7dmcgt;
            return(OK);
        case BSIM4v7_MOD_XGW:
            value->rValue = model->BSIM4v7xgw;
            return(OK);
        case BSIM4v7_MOD_XGL:
            value->rValue = model->BSIM4v7xgl;
            return(OK);
        case BSIM4v7_MOD_RSHG:
            value->rValue = model->BSIM4v7rshg;
            return(OK);
        case BSIM4v7_MOD_NGCON:
            value->rValue = model->BSIM4v7ngcon; 
            return(OK);
        case BSIM4v7_MOD_TCJ:
            value->rValue = model->BSIM4v7tcj;
            return(OK);
        case BSIM4v7_MOD_TPB:
            value->rValue = model->BSIM4v7tpb;
            return(OK);
        case BSIM4v7_MOD_TCJSW:
            value->rValue = model->BSIM4v7tcjsw;
            return(OK);
        case BSIM4v7_MOD_TPBSW:
            value->rValue = model->BSIM4v7tpbsw;
            return(OK);
        case BSIM4v7_MOD_TCJSWG:
            value->rValue = model->BSIM4v7tcjswg;
            return(OK);
        case BSIM4v7_MOD_TPBSWG:
            value->rValue = model->BSIM4v7tpbswg;
            return(OK);

	/* Length dependence */
        case  BSIM4v7_MOD_LCDSC :
          value->rValue = model->BSIM4v7lcdsc;
            return(OK);
        case  BSIM4v7_MOD_LCDSCB :
          value->rValue = model->BSIM4v7lcdscb;
            return(OK);
        case  BSIM4v7_MOD_LCDSCD :
          value->rValue = model->BSIM4v7lcdscd;
            return(OK);
        case  BSIM4v7_MOD_LCIT :
          value->rValue = model->BSIM4v7lcit;
            return(OK);
        case  BSIM4v7_MOD_LNFACTOR :
          value->rValue = model->BSIM4v7lnfactor;
            return(OK);
        case BSIM4v7_MOD_LXJ:
            value->rValue = model->BSIM4v7lxj;
            return(OK);
        case BSIM4v7_MOD_LVSAT:
            value->rValue = model->BSIM4v7lvsat;
            return(OK);
        case BSIM4v7_MOD_LAT:
            value->rValue = model->BSIM4v7lat;
            return(OK);
        case BSIM4v7_MOD_LA0:
            value->rValue = model->BSIM4v7la0;
            return(OK);
        case BSIM4v7_MOD_LAGS:
            value->rValue = model->BSIM4v7lags;
            return(OK);
        case BSIM4v7_MOD_LA1:
            value->rValue = model->BSIM4v7la1;
            return(OK);
        case BSIM4v7_MOD_LA2:
            value->rValue = model->BSIM4v7la2;
            return(OK);
        case BSIM4v7_MOD_LKETA:
            value->rValue = model->BSIM4v7lketa;
            return(OK);   
        case BSIM4v7_MOD_LNSUB:
            value->rValue = model->BSIM4v7lnsub;
            return(OK);
        case BSIM4v7_MOD_LNDEP:
            value->rValue = model->BSIM4v7lndep;
            return(OK);
        case BSIM4v7_MOD_LNSD:
            value->rValue = model->BSIM4v7lnsd;
            return(OK);
        case BSIM4v7_MOD_LNGATE:
            value->rValue = model->BSIM4v7lngate;
            return(OK);
        case BSIM4v7_MOD_LGAMMA1:
            value->rValue = model->BSIM4v7lgamma1;
            return(OK);
        case BSIM4v7_MOD_LGAMMA2:
            value->rValue = model->BSIM4v7lgamma2;
            return(OK);
        case BSIM4v7_MOD_LVBX:
            value->rValue = model->BSIM4v7lvbx;
            return(OK);
        case BSIM4v7_MOD_LVBM:
            value->rValue = model->BSIM4v7lvbm;
            return(OK);
        case BSIM4v7_MOD_LXT:
            value->rValue = model->BSIM4v7lxt;
            return(OK);
        case  BSIM4v7_MOD_LK1:
          value->rValue = model->BSIM4v7lk1;
            return(OK);
        case  BSIM4v7_MOD_LKT1:
          value->rValue = model->BSIM4v7lkt1;
            return(OK);
        case  BSIM4v7_MOD_LKT1L:
          value->rValue = model->BSIM4v7lkt1l;
            return(OK);
        case  BSIM4v7_MOD_LKT2 :
          value->rValue = model->BSIM4v7lkt2;
            return(OK);
        case  BSIM4v7_MOD_LK2 :
          value->rValue = model->BSIM4v7lk2;
            return(OK);
        case  BSIM4v7_MOD_LK3:
          value->rValue = model->BSIM4v7lk3;
            return(OK);
        case  BSIM4v7_MOD_LK3B:
          value->rValue = model->BSIM4v7lk3b;
            return(OK);
        case  BSIM4v7_MOD_LW0:
          value->rValue = model->BSIM4v7lw0;
            return(OK);
        case  BSIM4v7_MOD_LLPE0:
          value->rValue = model->BSIM4v7llpe0;
            return(OK);
        case  BSIM4v7_MOD_LLPEB:
          value->rValue = model->BSIM4v7llpeb;
            return(OK);
        case  BSIM4v7_MOD_LDVTP0:
          value->rValue = model->BSIM4v7ldvtp0;
            return(OK);
        case  BSIM4v7_MOD_LDVTP1:
          value->rValue = model->BSIM4v7ldvtp1;
            return(OK);
	case  BSIM4v7_MOD_LDVTP2:
          value->rValue = model->BSIM4v7ldvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4v7_MOD_LDVTP3:
          value->rValue = model->BSIM4v7ldvtp3;
            return(OK);
        case  BSIM4v7_MOD_LDVTP4:
          value->rValue = model->BSIM4v7ldvtp4;
            return(OK);
        case  BSIM4v7_MOD_LDVTP5:
          value->rValue = model->BSIM4v7ldvtp5;
            return(OK);
        case  BSIM4v7_MOD_LDVT0:                
          value->rValue = model->BSIM4v7ldvt0;
            return(OK);
        case  BSIM4v7_MOD_LDVT1 :             
          value->rValue = model->BSIM4v7ldvt1;
            return(OK);
        case  BSIM4v7_MOD_LDVT2 :             
          value->rValue = model->BSIM4v7ldvt2;
            return(OK);
        case  BSIM4v7_MOD_LDVT0W :                
          value->rValue = model->BSIM4v7ldvt0w;
            return(OK);
        case  BSIM4v7_MOD_LDVT1W :             
          value->rValue = model->BSIM4v7ldvt1w;
            return(OK);
        case  BSIM4v7_MOD_LDVT2W :             
          value->rValue = model->BSIM4v7ldvt2w;
            return(OK);
        case  BSIM4v7_MOD_LDROUT :           
          value->rValue = model->BSIM4v7ldrout;
            return(OK);
        case  BSIM4v7_MOD_LDSUB :           
          value->rValue = model->BSIM4v7ldsub;
            return(OK);
        case BSIM4v7_MOD_LVTH0:
            value->rValue = model->BSIM4v7lvth0; 
            return(OK);
        case BSIM4v7_MOD_LUA:
            value->rValue = model->BSIM4v7lua; 
            return(OK);
        case BSIM4v7_MOD_LUA1:
            value->rValue = model->BSIM4v7lua1; 
            return(OK);
        case BSIM4v7_MOD_LUB:
            value->rValue = model->BSIM4v7lub;  
            return(OK);
        case BSIM4v7_MOD_LUB1:
            value->rValue = model->BSIM4v7lub1;  
            return(OK);
        case BSIM4v7_MOD_LUC:
            value->rValue = model->BSIM4v7luc; 
            return(OK);
        case BSIM4v7_MOD_LUC1:
            value->rValue = model->BSIM4v7luc1; 
            return(OK);
        case BSIM4v7_MOD_LUD:
            value->rValue = model->BSIM4v7lud; 
            return(OK);
        case BSIM4v7_MOD_LUD1:
            value->rValue = model->BSIM4v7lud1; 
            return(OK);
        case BSIM4v7_MOD_LUP:
            value->rValue = model->BSIM4v7lup; 
            return(OK);
        case BSIM4v7_MOD_LLP:
            value->rValue = model->BSIM4v7llp; 
            return(OK);
        case BSIM4v7_MOD_LU0:
            value->rValue = model->BSIM4v7lu0;
            return(OK);
        case BSIM4v7_MOD_LUTE:
            value->rValue = model->BSIM4v7lute;
            return(OK);
		case BSIM4v7_MOD_LUCSTE:
            value->rValue = model->BSIM4v7lucste;
            return(OK);
        case BSIM4v7_MOD_LVOFF:
            value->rValue = model->BSIM4v7lvoff;
            return(OK);
        case BSIM4v7_MOD_LTVOFF:
            value->rValue = model->BSIM4v7ltvoff;
            return(OK);
        case BSIM4v7_MOD_LTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7ltnfactor;
            return(OK);
        case BSIM4v7_MOD_LTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7lteta0;
            return(OK);
        case BSIM4v7_MOD_LTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7ltvoffcv;
            return(OK);
        case BSIM4v7_MOD_LMINV:
            value->rValue = model->BSIM4v7lminv;
            return(OK);
        case BSIM4v7_MOD_LMINVCV:
            value->rValue = model->BSIM4v7lminvcv;
            return(OK);
        case BSIM4v7_MOD_LFPROUT:
            value->rValue = model->BSIM4v7lfprout;
            return(OK);
        case BSIM4v7_MOD_LPDITS:
            value->rValue = model->BSIM4v7lpdits;
            return(OK);
        case BSIM4v7_MOD_LPDITSD:
            value->rValue = model->BSIM4v7lpditsd;
            return(OK);
        case BSIM4v7_MOD_LDELTA:
            value->rValue = model->BSIM4v7ldelta;
            return(OK);
        case BSIM4v7_MOD_LRDSW:
            value->rValue = model->BSIM4v7lrdsw; 
            return(OK);             
        case BSIM4v7_MOD_LRDW:
            value->rValue = model->BSIM4v7lrdw;
            return(OK);
        case BSIM4v7_MOD_LRSW:
            value->rValue = model->BSIM4v7lrsw;
            return(OK);
        case BSIM4v7_MOD_LPRWB:
            value->rValue = model->BSIM4v7lprwb; 
            return(OK);             
        case BSIM4v7_MOD_LPRWG:
            value->rValue = model->BSIM4v7lprwg; 
            return(OK);             
        case BSIM4v7_MOD_LPRT:
            value->rValue = model->BSIM4v7lprt; 
            return(OK);              
        case BSIM4v7_MOD_LETA0:
            value->rValue = model->BSIM4v7leta0; 
            return(OK);               
        case BSIM4v7_MOD_LETAB:
            value->rValue = model->BSIM4v7letab; 
            return(OK);               
        case BSIM4v7_MOD_LPCLM:
            value->rValue = model->BSIM4v7lpclm; 
            return(OK);               
        case BSIM4v7_MOD_LPDIBL1:
            value->rValue = model->BSIM4v7lpdibl1; 
            return(OK);               
        case BSIM4v7_MOD_LPDIBL2:
            value->rValue = model->BSIM4v7lpdibl2; 
            return(OK);               
        case BSIM4v7_MOD_LPDIBLB:
            value->rValue = model->BSIM4v7lpdiblb; 
            return(OK);               
        case BSIM4v7_MOD_LPSCBE1:
            value->rValue = model->BSIM4v7lpscbe1; 
            return(OK);               
        case BSIM4v7_MOD_LPSCBE2:
            value->rValue = model->BSIM4v7lpscbe2; 
            return(OK);               
        case BSIM4v7_MOD_LPVAG:
            value->rValue = model->BSIM4v7lpvag; 
            return(OK);               
        case BSIM4v7_MOD_LWR:
            value->rValue = model->BSIM4v7lwr;
            return(OK);
        case BSIM4v7_MOD_LDWG:
            value->rValue = model->BSIM4v7ldwg;
            return(OK);
        case BSIM4v7_MOD_LDWB:
            value->rValue = model->BSIM4v7ldwb;
            return(OK);
        case BSIM4v7_MOD_LB0:
            value->rValue = model->BSIM4v7lb0;
            return(OK);
        case BSIM4v7_MOD_LB1:
            value->rValue = model->BSIM4v7lb1;
            return(OK);
        case BSIM4v7_MOD_LALPHA0:
            value->rValue = model->BSIM4v7lalpha0;
            return(OK);
        case BSIM4v7_MOD_LALPHA1:
            value->rValue = model->BSIM4v7lalpha1;
            return(OK);
        case BSIM4v7_MOD_LBETA0:
            value->rValue = model->BSIM4v7lbeta0;
            return(OK);
        case BSIM4v7_MOD_LAGIDL:
            value->rValue = model->BSIM4v7lagidl;
            return(OK);
        case BSIM4v7_MOD_LBGIDL:
            value->rValue = model->BSIM4v7lbgidl;
            return(OK);
        case BSIM4v7_MOD_LCGIDL:
            value->rValue = model->BSIM4v7lcgidl;
            return(OK);
	case BSIM4v7_MOD_LEGIDL:
            value->rValue = model->BSIM4v7legidl;
            return(OK);
        case BSIM4v7_MOD_LFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lfgidl;
            return(OK);
 	case BSIM4v7_MOD_LKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lkgidl;
            return(OK);
 	case BSIM4v7_MOD_LRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lrgidl;
            return(OK);
        case BSIM4v7_MOD_LAGISL:
            value->rValue = model->BSIM4v7lagisl;
            return(OK);
        case BSIM4v7_MOD_LBGISL:
            value->rValue = model->BSIM4v7lbgisl;
            return(OK);
        case BSIM4v7_MOD_LCGISL:
            value->rValue = model->BSIM4v7lcgisl;
            return(OK);
        case BSIM4v7_MOD_LEGISL:
            value->rValue = model->BSIM4v7legisl;
            return(OK);
        case BSIM4v7_MOD_LFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lfgisl;
            return(OK);
 	case BSIM4v7_MOD_LKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lkgisl;
            return(OK);
 	case BSIM4v7_MOD_LRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7lrgisl;
            return(OK);
        case BSIM4v7_MOD_LAIGC:
            value->rValue = model->BSIM4v7laigc;
            return(OK);
        case BSIM4v7_MOD_LBIGC:
            value->rValue = model->BSIM4v7lbigc;
            return(OK);
        case BSIM4v7_MOD_LCIGC:
            value->rValue = model->BSIM4v7lcigc;
            return(OK);
        case BSIM4v7_MOD_LAIGSD:
            value->rValue = model->BSIM4v7laigsd;
            return(OK);
        case BSIM4v7_MOD_LBIGSD:
            value->rValue = model->BSIM4v7lbigsd;
            return(OK);
        case BSIM4v7_MOD_LCIGSD:
            value->rValue = model->BSIM4v7lcigsd;
            return(OK);
        case BSIM4v7_MOD_LAIGS:
            value->rValue = model->BSIM4v7laigs;
            return(OK);
        case BSIM4v7_MOD_LBIGS:
            value->rValue = model->BSIM4v7lbigs;
            return(OK);
        case BSIM4v7_MOD_LCIGS:
            value->rValue = model->BSIM4v7lcigs;
            return(OK);
        case BSIM4v7_MOD_LAIGD:
            value->rValue = model->BSIM4v7laigd;
            return(OK);
        case BSIM4v7_MOD_LBIGD:
            value->rValue = model->BSIM4v7lbigd;
            return(OK);
        case BSIM4v7_MOD_LCIGD:
            value->rValue = model->BSIM4v7lcigd;
            return(OK);
        case BSIM4v7_MOD_LAIGBACC:
            value->rValue = model->BSIM4v7laigbacc;
            return(OK);
        case BSIM4v7_MOD_LBIGBACC:
            value->rValue = model->BSIM4v7lbigbacc;
            return(OK);
        case BSIM4v7_MOD_LCIGBACC:
            value->rValue = model->BSIM4v7lcigbacc;
            return(OK);
        case BSIM4v7_MOD_LAIGBINV:
            value->rValue = model->BSIM4v7laigbinv;
            return(OK);
        case BSIM4v7_MOD_LBIGBINV:
            value->rValue = model->BSIM4v7lbigbinv;
            return(OK);
        case BSIM4v7_MOD_LCIGBINV:
            value->rValue = model->BSIM4v7lcigbinv;
            return(OK);
        case BSIM4v7_MOD_LNIGC:
            value->rValue = model->BSIM4v7lnigc;
            return(OK);
        case BSIM4v7_MOD_LNIGBACC:
            value->rValue = model->BSIM4v7lnigbacc;
            return(OK);
        case BSIM4v7_MOD_LNIGBINV:
            value->rValue = model->BSIM4v7lnigbinv;
            return(OK);
        case BSIM4v7_MOD_LNTOX:
            value->rValue = model->BSIM4v7lntox;
            return(OK);
        case BSIM4v7_MOD_LEIGBINV:
            value->rValue = model->BSIM4v7leigbinv;
            return(OK);
        case BSIM4v7_MOD_LPIGCD:
            value->rValue = model->BSIM4v7lpigcd;
            return(OK);
        case BSIM4v7_MOD_LPOXEDGE:
            value->rValue = model->BSIM4v7lpoxedge;
            return(OK);
        case BSIM4v7_MOD_LPHIN:
            value->rValue = model->BSIM4v7lphin;
            return(OK);
        case BSIM4v7_MOD_LXRCRG1:
            value->rValue = model->BSIM4v7lxrcrg1;
            return(OK);
        case BSIM4v7_MOD_LXRCRG2:
            value->rValue = model->BSIM4v7lxrcrg2;
            return(OK);
        case BSIM4v7_MOD_LEU:
            value->rValue = model->BSIM4v7leu;
            return(OK);
		case BSIM4v7_MOD_LUCS:
            value->rValue = model->BSIM4v7lucs;
            return(OK);
        case BSIM4v7_MOD_LVFB:
            value->rValue = model->BSIM4v7lvfb;
            return(OK);

        case BSIM4v7_MOD_LCGSL:
            value->rValue = model->BSIM4v7lcgsl;
            return(OK);
        case BSIM4v7_MOD_LCGDL:
            value->rValue = model->BSIM4v7lcgdl;
            return(OK);
        case BSIM4v7_MOD_LCKAPPAS:
            value->rValue = model->BSIM4v7lckappas;
            return(OK);
        case BSIM4v7_MOD_LCKAPPAD:
            value->rValue = model->BSIM4v7lckappad;
            return(OK);
        case BSIM4v7_MOD_LCF:
            value->rValue = model->BSIM4v7lcf;
            return(OK);
        case BSIM4v7_MOD_LCLC:
            value->rValue = model->BSIM4v7lclc;
            return(OK);
        case BSIM4v7_MOD_LCLE:
            value->rValue = model->BSIM4v7lcle;
            return(OK);
        case BSIM4v7_MOD_LVFBCV:
            value->rValue = model->BSIM4v7lvfbcv;
            return(OK);
        case BSIM4v7_MOD_LACDE:
            value->rValue = model->BSIM4v7lacde;
            return(OK);
        case BSIM4v7_MOD_LMOIN:
            value->rValue = model->BSIM4v7lmoin;
            return(OK);
        case BSIM4v7_MOD_LNOFF:
            value->rValue = model->BSIM4v7lnoff;
            return(OK);
        case BSIM4v7_MOD_LVOFFCV:
            value->rValue = model->BSIM4v7lvoffcv;
            return(OK);
        case BSIM4v7_MOD_LVFBSDOFF:
            value->rValue = model->BSIM4v7lvfbsdoff;
            return(OK);
        case BSIM4v7_MOD_LTVFBSDOFF:
            value->rValue = model->BSIM4v7ltvfbsdoff;
            return(OK);

        case BSIM4v7_MOD_LLAMBDA:
            value->rValue = model->BSIM4v7llambda;
            return(OK);
        case BSIM4v7_MOD_LVTL:
            value->rValue = model->BSIM4v7lvtl;
            return(OK);
        case BSIM4v7_MOD_LXN:
            value->rValue = model->BSIM4v7lxn;
            return(OK);

	/* Width dependence */
        case  BSIM4v7_MOD_WCDSC :
          value->rValue = model->BSIM4v7wcdsc;
            return(OK);
        case  BSIM4v7_MOD_WCDSCB :
          value->rValue = model->BSIM4v7wcdscb;
            return(OK);
        case  BSIM4v7_MOD_WCDSCD :
          value->rValue = model->BSIM4v7wcdscd;
            return(OK);
        case  BSIM4v7_MOD_WCIT :
          value->rValue = model->BSIM4v7wcit;
            return(OK);
        case  BSIM4v7_MOD_WNFACTOR :
          value->rValue = model->BSIM4v7wnfactor;
            return(OK);
        case BSIM4v7_MOD_WXJ:
            value->rValue = model->BSIM4v7wxj;
            return(OK);
        case BSIM4v7_MOD_WVSAT:
            value->rValue = model->BSIM4v7wvsat;
            return(OK);
        case BSIM4v7_MOD_WAT:
            value->rValue = model->BSIM4v7wat;
            return(OK);
        case BSIM4v7_MOD_WA0:
            value->rValue = model->BSIM4v7wa0;
            return(OK);
        case BSIM4v7_MOD_WAGS:
            value->rValue = model->BSIM4v7wags;
            return(OK);
        case BSIM4v7_MOD_WA1:
            value->rValue = model->BSIM4v7wa1;
            return(OK);
        case BSIM4v7_MOD_WA2:
            value->rValue = model->BSIM4v7wa2;
            return(OK);
        case BSIM4v7_MOD_WKETA:
            value->rValue = model->BSIM4v7wketa;
            return(OK);   
        case BSIM4v7_MOD_WNSUB:
            value->rValue = model->BSIM4v7wnsub;
            return(OK);
        case BSIM4v7_MOD_WNDEP:
            value->rValue = model->BSIM4v7wndep;
            return(OK);
        case BSIM4v7_MOD_WNSD:
            value->rValue = model->BSIM4v7wnsd;
            return(OK);
        case BSIM4v7_MOD_WNGATE:
            value->rValue = model->BSIM4v7wngate;
            return(OK);
        case BSIM4v7_MOD_WGAMMA1:
            value->rValue = model->BSIM4v7wgamma1;
            return(OK);
        case BSIM4v7_MOD_WGAMMA2:
            value->rValue = model->BSIM4v7wgamma2;
            return(OK);
        case BSIM4v7_MOD_WVBX:
            value->rValue = model->BSIM4v7wvbx;
            return(OK);
        case BSIM4v7_MOD_WVBM:
            value->rValue = model->BSIM4v7wvbm;
            return(OK);
        case BSIM4v7_MOD_WXT:
            value->rValue = model->BSIM4v7wxt;
            return(OK);
        case  BSIM4v7_MOD_WK1:
          value->rValue = model->BSIM4v7wk1;
            return(OK);
        case  BSIM4v7_MOD_WKT1:
          value->rValue = model->BSIM4v7wkt1;
            return(OK);
        case  BSIM4v7_MOD_WKT1L:
          value->rValue = model->BSIM4v7wkt1l;
            return(OK);
        case  BSIM4v7_MOD_WKT2 :
          value->rValue = model->BSIM4v7wkt2;
            return(OK);
        case  BSIM4v7_MOD_WK2 :
          value->rValue = model->BSIM4v7wk2;
            return(OK);
        case  BSIM4v7_MOD_WK3:
          value->rValue = model->BSIM4v7wk3;
            return(OK);
        case  BSIM4v7_MOD_WK3B:
          value->rValue = model->BSIM4v7wk3b;
            return(OK);
        case  BSIM4v7_MOD_WW0:
          value->rValue = model->BSIM4v7ww0;
            return(OK);
        case  BSIM4v7_MOD_WLPE0:
          value->rValue = model->BSIM4v7wlpe0;
            return(OK);
        case  BSIM4v7_MOD_WDVTP0:
          value->rValue = model->BSIM4v7wdvtp0;
            return(OK);
        case  BSIM4v7_MOD_WDVTP1:
          value->rValue = model->BSIM4v7wdvtp1;
            return(OK);
        case  BSIM4v7_MOD_WDVTP2:
          value->rValue = model->BSIM4v7wdvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4v7_MOD_WDVTP3:
          value->rValue = model->BSIM4v7wdvtp3;
            return(OK);
        case  BSIM4v7_MOD_WDVTP4:
          value->rValue = model->BSIM4v7wdvtp4;
            return(OK);
        case  BSIM4v7_MOD_WDVTP5:
          value->rValue = model->BSIM4v7wdvtp5;
            return(OK);
        case  BSIM4v7_MOD_WLPEB:
          value->rValue = model->BSIM4v7wlpeb;
            return(OK);
        case  BSIM4v7_MOD_WDVT0:                
          value->rValue = model->BSIM4v7wdvt0;
            return(OK);
        case  BSIM4v7_MOD_WDVT1 :             
          value->rValue = model->BSIM4v7wdvt1;
            return(OK);
        case  BSIM4v7_MOD_WDVT2 :             
          value->rValue = model->BSIM4v7wdvt2;
            return(OK);
        case  BSIM4v7_MOD_WDVT0W :                
          value->rValue = model->BSIM4v7wdvt0w;
            return(OK);
        case  BSIM4v7_MOD_WDVT1W :             
          value->rValue = model->BSIM4v7wdvt1w;
            return(OK);
        case  BSIM4v7_MOD_WDVT2W :             
          value->rValue = model->BSIM4v7wdvt2w;
            return(OK);
        case  BSIM4v7_MOD_WDROUT :           
          value->rValue = model->BSIM4v7wdrout;
            return(OK);
        case  BSIM4v7_MOD_WDSUB :           
          value->rValue = model->BSIM4v7wdsub;
            return(OK);
        case BSIM4v7_MOD_WVTH0:
            value->rValue = model->BSIM4v7wvth0; 
            return(OK);
        case BSIM4v7_MOD_WUA:
            value->rValue = model->BSIM4v7wua; 
            return(OK);
        case BSIM4v7_MOD_WUA1:
            value->rValue = model->BSIM4v7wua1; 
            return(OK);
        case BSIM4v7_MOD_WUB:
            value->rValue = model->BSIM4v7wub;  
            return(OK);
        case BSIM4v7_MOD_WUB1:
            value->rValue = model->BSIM4v7wub1;  
            return(OK);
        case BSIM4v7_MOD_WUC:
            value->rValue = model->BSIM4v7wuc; 
            return(OK);
        case BSIM4v7_MOD_WUC1:
            value->rValue = model->BSIM4v7wuc1; 
            return(OK);
        case BSIM4v7_MOD_WUD:
            value->rValue = model->BSIM4v7wud; 
            return(OK);
        case BSIM4v7_MOD_WUD1:
            value->rValue = model->BSIM4v7wud1; 
            return(OK);
        case BSIM4v7_MOD_WUP:
            value->rValue = model->BSIM4v7wup; 
            return(OK);
        case BSIM4v7_MOD_WLP:
            value->rValue = model->BSIM4v7wlp; 
            return(OK);
        case BSIM4v7_MOD_WU0:
            value->rValue = model->BSIM4v7wu0;
            return(OK);
        case BSIM4v7_MOD_WUTE:
            value->rValue = model->BSIM4v7wute;
            return(OK);
        case BSIM4v7_MOD_WUCSTE:
            value->rValue = model->BSIM4v7wucste;
            return(OK);
        case BSIM4v7_MOD_WVOFF:
            value->rValue = model->BSIM4v7wvoff;
            return(OK);
        case BSIM4v7_MOD_WTVOFF:
            value->rValue = model->BSIM4v7wtvoff;
            return(OK);
        case BSIM4v7_MOD_WTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7wtnfactor;
            return(OK);
        case BSIM4v7_MOD_WTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7wteta0;
            return(OK);
        case BSIM4v7_MOD_WTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7wtvoffcv;
            return(OK);
        case BSIM4v7_MOD_WMINV:
            value->rValue = model->BSIM4v7wminv;
            return(OK);
        case BSIM4v7_MOD_WMINVCV:
            value->rValue = model->BSIM4v7wminvcv;
            return(OK);
        case BSIM4v7_MOD_WFPROUT:
            value->rValue = model->BSIM4v7wfprout;
            return(OK);
        case BSIM4v7_MOD_WPDITS:
            value->rValue = model->BSIM4v7wpdits;
            return(OK);
        case BSIM4v7_MOD_WPDITSD:
            value->rValue = model->BSIM4v7wpditsd;
            return(OK);
        case BSIM4v7_MOD_WDELTA:
            value->rValue = model->BSIM4v7wdelta;
            return(OK);
        case BSIM4v7_MOD_WRDSW:
            value->rValue = model->BSIM4v7wrdsw; 
            return(OK);             
        case BSIM4v7_MOD_WRDW:
            value->rValue = model->BSIM4v7wrdw;
            return(OK);
        case BSIM4v7_MOD_WRSW:
            value->rValue = model->BSIM4v7wrsw;
            return(OK);
        case BSIM4v7_MOD_WPRWB:
            value->rValue = model->BSIM4v7wprwb; 
            return(OK);             
        case BSIM4v7_MOD_WPRWG:
            value->rValue = model->BSIM4v7wprwg; 
            return(OK);             
        case BSIM4v7_MOD_WPRT:
            value->rValue = model->BSIM4v7wprt; 
            return(OK);              
        case BSIM4v7_MOD_WETA0:
            value->rValue = model->BSIM4v7weta0; 
            return(OK);               
        case BSIM4v7_MOD_WETAB:
            value->rValue = model->BSIM4v7wetab; 
            return(OK);               
        case BSIM4v7_MOD_WPCLM:
            value->rValue = model->BSIM4v7wpclm; 
            return(OK);               
        case BSIM4v7_MOD_WPDIBL1:
            value->rValue = model->BSIM4v7wpdibl1; 
            return(OK);               
        case BSIM4v7_MOD_WPDIBL2:
            value->rValue = model->BSIM4v7wpdibl2; 
            return(OK);               
        case BSIM4v7_MOD_WPDIBLB:
            value->rValue = model->BSIM4v7wpdiblb; 
            return(OK);               
        case BSIM4v7_MOD_WPSCBE1:
            value->rValue = model->BSIM4v7wpscbe1; 
            return(OK);               
        case BSIM4v7_MOD_WPSCBE2:
            value->rValue = model->BSIM4v7wpscbe2; 
            return(OK);               
        case BSIM4v7_MOD_WPVAG:
            value->rValue = model->BSIM4v7wpvag; 
            return(OK);               
        case BSIM4v7_MOD_WWR:
            value->rValue = model->BSIM4v7wwr;
            return(OK);
        case BSIM4v7_MOD_WDWG:
            value->rValue = model->BSIM4v7wdwg;
            return(OK);
        case BSIM4v7_MOD_WDWB:
            value->rValue = model->BSIM4v7wdwb;
            return(OK);
        case BSIM4v7_MOD_WB0:
            value->rValue = model->BSIM4v7wb0;
            return(OK);
        case BSIM4v7_MOD_WB1:
            value->rValue = model->BSIM4v7wb1;
            return(OK);
        case BSIM4v7_MOD_WALPHA0:
            value->rValue = model->BSIM4v7walpha0;
            return(OK);
        case BSIM4v7_MOD_WALPHA1:
            value->rValue = model->BSIM4v7walpha1;
            return(OK);
        case BSIM4v7_MOD_WBETA0:
            value->rValue = model->BSIM4v7wbeta0;
            return(OK);
        case BSIM4v7_MOD_WAGIDL:
            value->rValue = model->BSIM4v7wagidl;
            return(OK);
        case BSIM4v7_MOD_WBGIDL:
            value->rValue = model->BSIM4v7wbgidl;
            return(OK);
        case BSIM4v7_MOD_WCGIDL:
            value->rValue = model->BSIM4v7wcgidl;
            return(OK);
        case BSIM4v7_MOD_WEGIDL:
            value->rValue = model->BSIM4v7wegidl;
            return(OK);
        case BSIM4v7_MOD_WFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wfgidl;
            return(OK);
        case BSIM4v7_MOD_WKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wkgidl;
            return(OK);
        case BSIM4v7_MOD_WRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wrgidl;
            return(OK);
        case BSIM4v7_MOD_WAGISL:
            value->rValue = model->BSIM4v7wagisl;
            return(OK);
        case BSIM4v7_MOD_WBGISL:
            value->rValue = model->BSIM4v7wbgisl;
            return(OK);
        case BSIM4v7_MOD_WCGISL:
            value->rValue = model->BSIM4v7wcgisl;
            return(OK);
        case BSIM4v7_MOD_WEGISL:
            value->rValue = model->BSIM4v7wegisl;
            return(OK);
        case BSIM4v7_MOD_WFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wfgisl;
            return(OK);
        case BSIM4v7_MOD_WKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wkgisl;
            return(OK);
        case BSIM4v7_MOD_WRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7wrgisl;
            return(OK);
        case BSIM4v7_MOD_WAIGC:
            value->rValue = model->BSIM4v7waigc;
            return(OK);
        case BSIM4v7_MOD_WBIGC:
            value->rValue = model->BSIM4v7wbigc;
            return(OK);
        case BSIM4v7_MOD_WCIGC:
            value->rValue = model->BSIM4v7wcigc;
            return(OK);
        case BSIM4v7_MOD_WAIGSD:
            value->rValue = model->BSIM4v7waigsd;
            return(OK);
        case BSIM4v7_MOD_WBIGSD:
            value->rValue = model->BSIM4v7wbigsd;
            return(OK);
        case BSIM4v7_MOD_WCIGSD:
            value->rValue = model->BSIM4v7wcigsd;
            return(OK);
        case BSIM4v7_MOD_WAIGS:
            value->rValue = model->BSIM4v7waigs;
            return(OK);
        case BSIM4v7_MOD_WBIGS:
            value->rValue = model->BSIM4v7wbigs;
            return(OK);
        case BSIM4v7_MOD_WCIGS:
            value->rValue = model->BSIM4v7wcigs;
            return(OK);
        case BSIM4v7_MOD_WAIGD:
            value->rValue = model->BSIM4v7waigd;
            return(OK);
        case BSIM4v7_MOD_WBIGD:
            value->rValue = model->BSIM4v7wbigd;
            return(OK);
        case BSIM4v7_MOD_WCIGD:
            value->rValue = model->BSIM4v7wcigd;
            return(OK);
        case BSIM4v7_MOD_WAIGBACC:
            value->rValue = model->BSIM4v7waigbacc;
            return(OK);
        case BSIM4v7_MOD_WBIGBACC:
            value->rValue = model->BSIM4v7wbigbacc;
            return(OK);
        case BSIM4v7_MOD_WCIGBACC:
            value->rValue = model->BSIM4v7wcigbacc;
            return(OK);
        case BSIM4v7_MOD_WAIGBINV:
            value->rValue = model->BSIM4v7waigbinv;
            return(OK);
        case BSIM4v7_MOD_WBIGBINV:
            value->rValue = model->BSIM4v7wbigbinv;
            return(OK);
        case BSIM4v7_MOD_WCIGBINV:
            value->rValue = model->BSIM4v7wcigbinv;
            return(OK);
        case BSIM4v7_MOD_WNIGC:
            value->rValue = model->BSIM4v7wnigc;
            return(OK);
        case BSIM4v7_MOD_WNIGBACC:
            value->rValue = model->BSIM4v7wnigbacc;
            return(OK);
        case BSIM4v7_MOD_WNIGBINV:
            value->rValue = model->BSIM4v7wnigbinv;
            return(OK);
        case BSIM4v7_MOD_WNTOX:
            value->rValue = model->BSIM4v7wntox;
            return(OK);
        case BSIM4v7_MOD_WEIGBINV:
            value->rValue = model->BSIM4v7weigbinv;
            return(OK);
        case BSIM4v7_MOD_WPIGCD:
            value->rValue = model->BSIM4v7wpigcd;
            return(OK);
        case BSIM4v7_MOD_WPOXEDGE:
            value->rValue = model->BSIM4v7wpoxedge;
            return(OK);
        case BSIM4v7_MOD_WPHIN:
            value->rValue = model->BSIM4v7wphin;
            return(OK);
        case BSIM4v7_MOD_WXRCRG1:
            value->rValue = model->BSIM4v7wxrcrg1;
            return(OK);
        case BSIM4v7_MOD_WXRCRG2:
            value->rValue = model->BSIM4v7wxrcrg2;
            return(OK);
        case BSIM4v7_MOD_WEU:
            value->rValue = model->BSIM4v7weu;
            return(OK);
        case BSIM4v7_MOD_WUCS:
            value->rValue = model->BSIM4v7wucs;
            return(OK);
        case BSIM4v7_MOD_WVFB:
            value->rValue = model->BSIM4v7wvfb;
            return(OK);

        case BSIM4v7_MOD_WCGSL:
            value->rValue = model->BSIM4v7wcgsl;
            return(OK);
        case BSIM4v7_MOD_WCGDL:
            value->rValue = model->BSIM4v7wcgdl;
            return(OK);
        case BSIM4v7_MOD_WCKAPPAS:
            value->rValue = model->BSIM4v7wckappas;
            return(OK);
        case BSIM4v7_MOD_WCKAPPAD:
            value->rValue = model->BSIM4v7wckappad;
            return(OK);
        case BSIM4v7_MOD_WCF:
            value->rValue = model->BSIM4v7wcf;
            return(OK);
        case BSIM4v7_MOD_WCLC:
            value->rValue = model->BSIM4v7wclc;
            return(OK);
        case BSIM4v7_MOD_WCLE:
            value->rValue = model->BSIM4v7wcle;
            return(OK);
        case BSIM4v7_MOD_WVFBCV:
            value->rValue = model->BSIM4v7wvfbcv;
            return(OK);
        case BSIM4v7_MOD_WACDE:
            value->rValue = model->BSIM4v7wacde;
            return(OK);
        case BSIM4v7_MOD_WMOIN:
            value->rValue = model->BSIM4v7wmoin;
            return(OK);
        case BSIM4v7_MOD_WNOFF:
            value->rValue = model->BSIM4v7wnoff;
            return(OK);
        case BSIM4v7_MOD_WVOFFCV:
            value->rValue = model->BSIM4v7wvoffcv;
            return(OK);
        case BSIM4v7_MOD_WVFBSDOFF:
            value->rValue = model->BSIM4v7wvfbsdoff;
            return(OK);
        case BSIM4v7_MOD_WTVFBSDOFF:
            value->rValue = model->BSIM4v7wtvfbsdoff;
            return(OK);

        case BSIM4v7_MOD_WLAMBDA:
            value->rValue = model->BSIM4v7wlambda;
            return(OK);
        case BSIM4v7_MOD_WVTL:
            value->rValue = model->BSIM4v7wvtl;
            return(OK);
        case BSIM4v7_MOD_WXN:
            value->rValue = model->BSIM4v7wxn;
            return(OK);

	/* Cross-term dependence */
        case  BSIM4v7_MOD_PCDSC :
          value->rValue = model->BSIM4v7pcdsc;
            return(OK);
        case  BSIM4v7_MOD_PCDSCB :
          value->rValue = model->BSIM4v7pcdscb;
            return(OK);
        case  BSIM4v7_MOD_PCDSCD :
          value->rValue = model->BSIM4v7pcdscd;
            return(OK);
         case  BSIM4v7_MOD_PCIT :
          value->rValue = model->BSIM4v7pcit;
            return(OK);
        case  BSIM4v7_MOD_PNFACTOR :
          value->rValue = model->BSIM4v7pnfactor;
            return(OK);
        case BSIM4v7_MOD_PXJ:
            value->rValue = model->BSIM4v7pxj;
            return(OK);
        case BSIM4v7_MOD_PVSAT:
            value->rValue = model->BSIM4v7pvsat;
            return(OK);
        case BSIM4v7_MOD_PAT:
            value->rValue = model->BSIM4v7pat;
            return(OK);
        case BSIM4v7_MOD_PA0:
            value->rValue = model->BSIM4v7pa0;
            return(OK);
        case BSIM4v7_MOD_PAGS:
            value->rValue = model->BSIM4v7pags;
            return(OK);
        case BSIM4v7_MOD_PA1:
            value->rValue = model->BSIM4v7pa1;
            return(OK);
        case BSIM4v7_MOD_PA2:
            value->rValue = model->BSIM4v7pa2;
            return(OK);
        case BSIM4v7_MOD_PKETA:
            value->rValue = model->BSIM4v7pketa;
            return(OK);   
        case BSIM4v7_MOD_PNSUB:
            value->rValue = model->BSIM4v7pnsub;
            return(OK);
        case BSIM4v7_MOD_PNDEP:
            value->rValue = model->BSIM4v7pndep;
            return(OK);
        case BSIM4v7_MOD_PNSD:
            value->rValue = model->BSIM4v7pnsd;
            return(OK);
        case BSIM4v7_MOD_PNGATE:
            value->rValue = model->BSIM4v7pngate;
            return(OK);
        case BSIM4v7_MOD_PGAMMA1:
            value->rValue = model->BSIM4v7pgamma1;
            return(OK);
        case BSIM4v7_MOD_PGAMMA2:
            value->rValue = model->BSIM4v7pgamma2;
            return(OK);
        case BSIM4v7_MOD_PVBX:
            value->rValue = model->BSIM4v7pvbx;
            return(OK);
        case BSIM4v7_MOD_PVBM:
            value->rValue = model->BSIM4v7pvbm;
            return(OK);
        case BSIM4v7_MOD_PXT:
            value->rValue = model->BSIM4v7pxt;
            return(OK);
        case  BSIM4v7_MOD_PK1:
          value->rValue = model->BSIM4v7pk1;
            return(OK);
        case  BSIM4v7_MOD_PKT1:
          value->rValue = model->BSIM4v7pkt1;
            return(OK);
        case  BSIM4v7_MOD_PKT1L:
          value->rValue = model->BSIM4v7pkt1l;
            return(OK);
        case  BSIM4v7_MOD_PKT2 :
          value->rValue = model->BSIM4v7pkt2;
            return(OK);
        case  BSIM4v7_MOD_PK2 :
          value->rValue = model->BSIM4v7pk2;
            return(OK);
        case  BSIM4v7_MOD_PK3:
          value->rValue = model->BSIM4v7pk3;
            return(OK);
        case  BSIM4v7_MOD_PK3B:
          value->rValue = model->BSIM4v7pk3b;
            return(OK);
        case  BSIM4v7_MOD_PW0:
          value->rValue = model->BSIM4v7pw0;
            return(OK);
        case  BSIM4v7_MOD_PLPE0:
          value->rValue = model->BSIM4v7plpe0;
            return(OK);
        case  BSIM4v7_MOD_PLPEB:
          value->rValue = model->BSIM4v7plpeb;
            return(OK);
        case  BSIM4v7_MOD_PDVTP0:
          value->rValue = model->BSIM4v7pdvtp0;
            return(OK);
        case  BSIM4v7_MOD_PDVTP1:
          value->rValue = model->BSIM4v7pdvtp1;
            return(OK);
        case  BSIM4v7_MOD_PDVTP2:
          value->rValue = model->BSIM4v7pdvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4v7_MOD_PDVTP3:
          value->rValue = model->BSIM4v7pdvtp3;
            return(OK);
        case  BSIM4v7_MOD_PDVTP4:
          value->rValue = model->BSIM4v7pdvtp4;
            return(OK);
        case  BSIM4v7_MOD_PDVTP5:
          value->rValue = model->BSIM4v7pdvtp5;
            return(OK);
        case  BSIM4v7_MOD_PDVT0 :                
          value->rValue = model->BSIM4v7pdvt0;
            return(OK);
        case  BSIM4v7_MOD_PDVT1 :             
          value->rValue = model->BSIM4v7pdvt1;
            return(OK);
        case  BSIM4v7_MOD_PDVT2 :             
          value->rValue = model->BSIM4v7pdvt2;
            return(OK);
        case  BSIM4v7_MOD_PDVT0W :                
          value->rValue = model->BSIM4v7pdvt0w;
            return(OK);
        case  BSIM4v7_MOD_PDVT1W :             
          value->rValue = model->BSIM4v7pdvt1w;
            return(OK);
        case  BSIM4v7_MOD_PDVT2W :             
          value->rValue = model->BSIM4v7pdvt2w;
            return(OK);
        case  BSIM4v7_MOD_PDROUT :           
          value->rValue = model->BSIM4v7pdrout;
            return(OK);
        case  BSIM4v7_MOD_PDSUB :           
          value->rValue = model->BSIM4v7pdsub;
            return(OK);
        case BSIM4v7_MOD_PVTH0:
            value->rValue = model->BSIM4v7pvth0; 
            return(OK);
        case BSIM4v7_MOD_PUA:
            value->rValue = model->BSIM4v7pua; 
            return(OK);
        case BSIM4v7_MOD_PUA1:
            value->rValue = model->BSIM4v7pua1; 
            return(OK);
        case BSIM4v7_MOD_PUB:
            value->rValue = model->BSIM4v7pub;  
            return(OK);
        case BSIM4v7_MOD_PUB1:
            value->rValue = model->BSIM4v7pub1;  
            return(OK);
        case BSIM4v7_MOD_PUC:
            value->rValue = model->BSIM4v7puc; 
            return(OK);
        case BSIM4v7_MOD_PUC1:
            value->rValue = model->BSIM4v7puc1; 
            return(OK);
        case BSIM4v7_MOD_PUD:
            value->rValue = model->BSIM4v7pud; 
            return(OK);
        case BSIM4v7_MOD_PUD1:
            value->rValue = model->BSIM4v7pud1; 
            return(OK);
        case BSIM4v7_MOD_PUP:
            value->rValue = model->BSIM4v7pup; 
            return(OK);
        case BSIM4v7_MOD_PLP:
            value->rValue = model->BSIM4v7plp; 
            return(OK);
        case BSIM4v7_MOD_PU0:
            value->rValue = model->BSIM4v7pu0;
            return(OK);
        case BSIM4v7_MOD_PUTE:
            value->rValue = model->BSIM4v7pute;
            return(OK);
        case BSIM4v7_MOD_PUCSTE:
            value->rValue = model->BSIM4v7pucste;
            return(OK);
        case BSIM4v7_MOD_PVOFF:
            value->rValue = model->BSIM4v7pvoff;
            return(OK);
        case BSIM4v7_MOD_PTVOFF:
            value->rValue = model->BSIM4v7ptvoff;
            return(OK);
        case BSIM4v7_MOD_PTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7ptnfactor;
            return(OK);
        case BSIM4v7_MOD_PTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7pteta0;
            return(OK);
        case BSIM4v7_MOD_PTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4v7ptvoffcv;
            return(OK);
        case BSIM4v7_MOD_PMINV:
            value->rValue = model->BSIM4v7pminv;
            return(OK);
        case BSIM4v7_MOD_PMINVCV:
            value->rValue = model->BSIM4v7pminvcv;
            return(OK);
        case BSIM4v7_MOD_PFPROUT:
            value->rValue = model->BSIM4v7pfprout;
            return(OK);
        case BSIM4v7_MOD_PPDITS:
            value->rValue = model->BSIM4v7ppdits;
            return(OK);
        case BSIM4v7_MOD_PPDITSD:
            value->rValue = model->BSIM4v7ppditsd;
            return(OK);
        case BSIM4v7_MOD_PDELTA:
            value->rValue = model->BSIM4v7pdelta;
            return(OK);
        case BSIM4v7_MOD_PRDSW:
            value->rValue = model->BSIM4v7prdsw; 
            return(OK);             
        case BSIM4v7_MOD_PRDW:
            value->rValue = model->BSIM4v7prdw;
            return(OK);
        case BSIM4v7_MOD_PRSW:
            value->rValue = model->BSIM4v7prsw;
            return(OK);
        case BSIM4v7_MOD_PPRWB:
            value->rValue = model->BSIM4v7pprwb; 
            return(OK);             
        case BSIM4v7_MOD_PPRWG:
            value->rValue = model->BSIM4v7pprwg; 
            return(OK);             
        case BSIM4v7_MOD_PPRT:
            value->rValue = model->BSIM4v7pprt; 
            return(OK);              
        case BSIM4v7_MOD_PETA0:
            value->rValue = model->BSIM4v7peta0; 
            return(OK);               
        case BSIM4v7_MOD_PETAB:
            value->rValue = model->BSIM4v7petab; 
            return(OK);               
        case BSIM4v7_MOD_PPCLM:
            value->rValue = model->BSIM4v7ppclm; 
            return(OK);               
        case BSIM4v7_MOD_PPDIBL1:
            value->rValue = model->BSIM4v7ppdibl1; 
            return(OK);               
        case BSIM4v7_MOD_PPDIBL2:
            value->rValue = model->BSIM4v7ppdibl2; 
            return(OK);               
        case BSIM4v7_MOD_PPDIBLB:
            value->rValue = model->BSIM4v7ppdiblb; 
            return(OK);               
        case BSIM4v7_MOD_PPSCBE1:
            value->rValue = model->BSIM4v7ppscbe1; 
            return(OK);               
        case BSIM4v7_MOD_PPSCBE2:
            value->rValue = model->BSIM4v7ppscbe2; 
            return(OK);               
        case BSIM4v7_MOD_PPVAG:
            value->rValue = model->BSIM4v7ppvag; 
            return(OK);               
        case BSIM4v7_MOD_PWR:
            value->rValue = model->BSIM4v7pwr;
            return(OK);
        case BSIM4v7_MOD_PDWG:
            value->rValue = model->BSIM4v7pdwg;
            return(OK);
        case BSIM4v7_MOD_PDWB:
            value->rValue = model->BSIM4v7pdwb;
            return(OK);
        case BSIM4v7_MOD_PB0:
            value->rValue = model->BSIM4v7pb0;
            return(OK);
        case BSIM4v7_MOD_PB1:
            value->rValue = model->BSIM4v7pb1;
            return(OK);
        case BSIM4v7_MOD_PALPHA0:
            value->rValue = model->BSIM4v7palpha0;
            return(OK);
        case BSIM4v7_MOD_PALPHA1:
            value->rValue = model->BSIM4v7palpha1;
            return(OK);
        case BSIM4v7_MOD_PBETA0:
            value->rValue = model->BSIM4v7pbeta0;
            return(OK);
        case BSIM4v7_MOD_PAGIDL:
            value->rValue = model->BSIM4v7pagidl;
            return(OK);
        case BSIM4v7_MOD_PBGIDL:
            value->rValue = model->BSIM4v7pbgidl;
            return(OK);
        case BSIM4v7_MOD_PCGIDL:
            value->rValue = model->BSIM4v7pcgidl;
            return(OK);
        case BSIM4v7_MOD_PEGIDL:
            value->rValue = model->BSIM4v7pegidl;
            return(OK);
        case BSIM4v7_MOD_PFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7pfgidl;
            return(OK);
        case BSIM4v7_MOD_PKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7pkgidl;
            return(OK);
        case BSIM4v7_MOD_PRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7prgidl;
            return(OK);
        case BSIM4v7_MOD_PAGISL:
            value->rValue = model->BSIM4v7pagisl;
            return(OK);
        case BSIM4v7_MOD_PBGISL:
            value->rValue = model->BSIM4v7pbgisl;
            return(OK);
        case BSIM4v7_MOD_PCGISL:
            value->rValue = model->BSIM4v7pcgisl;
            return(OK);
        case BSIM4v7_MOD_PEGISL:
            value->rValue = model->BSIM4v7pegisl;
            return(OK);
        case BSIM4v7_MOD_PFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7pfgisl;
            return(OK);
        case BSIM4v7_MOD_PKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7pkgisl;
            return(OK);
        case BSIM4v7_MOD_PRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4v7prgisl;
            return(OK);
        case BSIM4v7_MOD_PAIGC:
            value->rValue = model->BSIM4v7paigc;
            return(OK);
        case BSIM4v7_MOD_PBIGC:
            value->rValue = model->BSIM4v7pbigc;
            return(OK);
        case BSIM4v7_MOD_PCIGC:
            value->rValue = model->BSIM4v7pcigc;
            return(OK);
        case BSIM4v7_MOD_PAIGSD:
            value->rValue = model->BSIM4v7paigsd;
            return(OK);
        case BSIM4v7_MOD_PBIGSD:
            value->rValue = model->BSIM4v7pbigsd;
            return(OK);
        case BSIM4v7_MOD_PCIGSD:
            value->rValue = model->BSIM4v7pcigsd;
            return(OK);
        case BSIM4v7_MOD_PAIGS:
            value->rValue = model->BSIM4v7paigs;
            return(OK);
        case BSIM4v7_MOD_PBIGS:
            value->rValue = model->BSIM4v7pbigs;
            return(OK);
        case BSIM4v7_MOD_PCIGS:
            value->rValue = model->BSIM4v7pcigs;
            return(OK);
        case BSIM4v7_MOD_PAIGD:
            value->rValue = model->BSIM4v7paigd;
            return(OK);
        case BSIM4v7_MOD_PBIGD:
            value->rValue = model->BSIM4v7pbigd;
            return(OK);
        case BSIM4v7_MOD_PCIGD:
            value->rValue = model->BSIM4v7pcigd;
            return(OK);
        case BSIM4v7_MOD_PAIGBACC:
            value->rValue = model->BSIM4v7paigbacc;
            return(OK);
        case BSIM4v7_MOD_PBIGBACC:
            value->rValue = model->BSIM4v7pbigbacc;
            return(OK);
        case BSIM4v7_MOD_PCIGBACC:
            value->rValue = model->BSIM4v7pcigbacc;
            return(OK);
        case BSIM4v7_MOD_PAIGBINV:
            value->rValue = model->BSIM4v7paigbinv;
            return(OK);
        case BSIM4v7_MOD_PBIGBINV:
            value->rValue = model->BSIM4v7pbigbinv;
            return(OK);
        case BSIM4v7_MOD_PCIGBINV:
            value->rValue = model->BSIM4v7pcigbinv;
            return(OK);
        case BSIM4v7_MOD_PNIGC:
            value->rValue = model->BSIM4v7pnigc;
            return(OK);
        case BSIM4v7_MOD_PNIGBACC:
            value->rValue = model->BSIM4v7pnigbacc;
            return(OK);
        case BSIM4v7_MOD_PNIGBINV:
            value->rValue = model->BSIM4v7pnigbinv;
            return(OK);
        case BSIM4v7_MOD_PNTOX:
            value->rValue = model->BSIM4v7pntox;
            return(OK);
        case BSIM4v7_MOD_PEIGBINV:
            value->rValue = model->BSIM4v7peigbinv;
            return(OK);
        case BSIM4v7_MOD_PPIGCD:
            value->rValue = model->BSIM4v7ppigcd;
            return(OK);
        case BSIM4v7_MOD_PPOXEDGE:
            value->rValue = model->BSIM4v7ppoxedge;
            return(OK);
        case BSIM4v7_MOD_PPHIN:
            value->rValue = model->BSIM4v7pphin;
            return(OK);
        case BSIM4v7_MOD_PXRCRG1:
            value->rValue = model->BSIM4v7pxrcrg1;
            return(OK);
        case BSIM4v7_MOD_PXRCRG2:
            value->rValue = model->BSIM4v7pxrcrg2;
            return(OK);
        case BSIM4v7_MOD_PEU:
            value->rValue = model->BSIM4v7peu;
            return(OK);
        case BSIM4v7_MOD_PUCS:
            value->rValue = model->BSIM4v7pucs;
            return(OK);
        case BSIM4v7_MOD_PVFB:
            value->rValue = model->BSIM4v7pvfb;
            return(OK);

        case BSIM4v7_MOD_PCGSL:
            value->rValue = model->BSIM4v7pcgsl;
            return(OK);
        case BSIM4v7_MOD_PCGDL:
            value->rValue = model->BSIM4v7pcgdl;
            return(OK);
        case BSIM4v7_MOD_PCKAPPAS:
            value->rValue = model->BSIM4v7pckappas;
            return(OK);
        case BSIM4v7_MOD_PCKAPPAD:
            value->rValue = model->BSIM4v7pckappad;
            return(OK);
        case BSIM4v7_MOD_PCF:
            value->rValue = model->BSIM4v7pcf;
            return(OK);
        case BSIM4v7_MOD_PCLC:
            value->rValue = model->BSIM4v7pclc;
            return(OK);
        case BSIM4v7_MOD_PCLE:
            value->rValue = model->BSIM4v7pcle;
            return(OK);
        case BSIM4v7_MOD_PVFBCV:
            value->rValue = model->BSIM4v7pvfbcv;
            return(OK);
        case BSIM4v7_MOD_PACDE:
            value->rValue = model->BSIM4v7pacde;
            return(OK);
        case BSIM4v7_MOD_PMOIN:
            value->rValue = model->BSIM4v7pmoin;
            return(OK);
        case BSIM4v7_MOD_PNOFF:
            value->rValue = model->BSIM4v7pnoff;
            return(OK);
        case BSIM4v7_MOD_PVOFFCV:
            value->rValue = model->BSIM4v7pvoffcv;
            return(OK);
        case BSIM4v7_MOD_PVFBSDOFF:
            value->rValue = model->BSIM4v7pvfbsdoff;
            return(OK);
        case BSIM4v7_MOD_PTVFBSDOFF:
            value->rValue = model->BSIM4v7ptvfbsdoff;
            return(OK);

        case BSIM4v7_MOD_PLAMBDA:
            value->rValue = model->BSIM4v7plambda;
            return(OK);
        case BSIM4v7_MOD_PVTL:
            value->rValue = model->BSIM4v7pvtl;
            return(OK);
        case BSIM4v7_MOD_PXN:
            value->rValue = model->BSIM4v7pxn;
            return(OK);

        case  BSIM4v7_MOD_TNOM :
          value->rValue = model->BSIM4v7tnom;
            return(OK);
        case BSIM4v7_MOD_CGSO:
            value->rValue = model->BSIM4v7cgso; 
            return(OK);
        case BSIM4v7_MOD_CGDO:
            value->rValue = model->BSIM4v7cgdo; 
            return(OK);
        case BSIM4v7_MOD_CGBO:
            value->rValue = model->BSIM4v7cgbo; 
            return(OK);
        case BSIM4v7_MOD_XPART:
            value->rValue = model->BSIM4v7xpart; 
            return(OK);
        case BSIM4v7_MOD_RSH:
            value->rValue = model->BSIM4v7sheetResistance; 
            return(OK);
        case BSIM4v7_MOD_JSS:
            value->rValue = model->BSIM4v7SjctSatCurDensity; 
            return(OK);
        case BSIM4v7_MOD_JSWS:
            value->rValue = model->BSIM4v7SjctSidewallSatCurDensity; 
            return(OK);
        case BSIM4v7_MOD_JSWGS:
            value->rValue = model->BSIM4v7SjctGateSidewallSatCurDensity;
            return(OK);
        case BSIM4v7_MOD_PBS:
            value->rValue = model->BSIM4v7SbulkJctPotential; 
            return(OK);
        case BSIM4v7_MOD_MJS:
            value->rValue = model->BSIM4v7SbulkJctBotGradingCoeff; 
            return(OK);
        case BSIM4v7_MOD_PBSWS:
            value->rValue = model->BSIM4v7SsidewallJctPotential; 
            return(OK);
        case BSIM4v7_MOD_MJSWS:
            value->rValue = model->BSIM4v7SbulkJctSideGradingCoeff; 
            return(OK);
        case BSIM4v7_MOD_CJS:
            value->rValue = model->BSIM4v7SunitAreaJctCap; 
            return(OK);
        case BSIM4v7_MOD_CJSWS:
            value->rValue = model->BSIM4v7SunitLengthSidewallJctCap; 
            return(OK);
        case BSIM4v7_MOD_PBSWGS:
            value->rValue = model->BSIM4v7SGatesidewallJctPotential; 
            return(OK);
        case BSIM4v7_MOD_MJSWGS:
            value->rValue = model->BSIM4v7SbulkJctGateSideGradingCoeff; 
            return(OK);
        case BSIM4v7_MOD_CJSWGS:
            value->rValue = model->BSIM4v7SunitLengthGateSidewallJctCap; 
            return(OK);
        case BSIM4v7_MOD_NJS:
            value->rValue = model->BSIM4v7SjctEmissionCoeff; 
            return(OK);
        case BSIM4v7_MOD_XTIS:
            value->rValue = model->BSIM4v7SjctTempExponent; 
            return(OK);
        case BSIM4v7_MOD_JSD:
            value->rValue = model->BSIM4v7DjctSatCurDensity;
            return(OK);
        case BSIM4v7_MOD_JSWD:
            value->rValue = model->BSIM4v7DjctSidewallSatCurDensity;
            return(OK);
        case BSIM4v7_MOD_JSWGD:
            value->rValue = model->BSIM4v7DjctGateSidewallSatCurDensity;
            return(OK);
        case BSIM4v7_MOD_PBD:
            value->rValue = model->BSIM4v7DbulkJctPotential;
            return(OK);
        case BSIM4v7_MOD_MJD:
            value->rValue = model->BSIM4v7DbulkJctBotGradingCoeff;
            return(OK);
        case BSIM4v7_MOD_PBSWD:
            value->rValue = model->BSIM4v7DsidewallJctPotential;
            return(OK);
        case BSIM4v7_MOD_MJSWD:
            value->rValue = model->BSIM4v7DbulkJctSideGradingCoeff;
            return(OK);
        case BSIM4v7_MOD_CJD:
            value->rValue = model->BSIM4v7DunitAreaJctCap;
            return(OK);
        case BSIM4v7_MOD_CJSWD:
            value->rValue = model->BSIM4v7DunitLengthSidewallJctCap;
            return(OK);
        case BSIM4v7_MOD_PBSWGD:
            value->rValue = model->BSIM4v7DGatesidewallJctPotential;
            return(OK);
        case BSIM4v7_MOD_MJSWGD:
            value->rValue = model->BSIM4v7DbulkJctGateSideGradingCoeff;
            return(OK);
        case BSIM4v7_MOD_CJSWGD:
            value->rValue = model->BSIM4v7DunitLengthGateSidewallJctCap;
            return(OK);
        case BSIM4v7_MOD_NJD:
            value->rValue = model->BSIM4v7DjctEmissionCoeff;
            return(OK);
        case BSIM4v7_MOD_XTID:
            value->rValue = model->BSIM4v7DjctTempExponent;
            return(OK);
        case BSIM4v7_MOD_LINTNOI:
            value->rValue = model->BSIM4v7lintnoi; 
            return(OK);
        case BSIM4v7_MOD_LINT:
            value->rValue = model->BSIM4v7Lint; 
            return(OK);
        case BSIM4v7_MOD_LL:
            value->rValue = model->BSIM4v7Ll;
            return(OK);
        case BSIM4v7_MOD_LLC:
            value->rValue = model->BSIM4v7Llc;
            return(OK);
        case BSIM4v7_MOD_LLN:
            value->rValue = model->BSIM4v7Lln;
            return(OK);
        case BSIM4v7_MOD_LW:
            value->rValue = model->BSIM4v7Lw;
            return(OK);
        case BSIM4v7_MOD_LWC:
            value->rValue = model->BSIM4v7Lwc;
            return(OK);
        case BSIM4v7_MOD_LWN:
            value->rValue = model->BSIM4v7Lwn;
            return(OK);
        case BSIM4v7_MOD_LWL:
            value->rValue = model->BSIM4v7Lwl;
            return(OK);
        case BSIM4v7_MOD_LWLC:
            value->rValue = model->BSIM4v7Lwlc;
            return(OK);
        case BSIM4v7_MOD_LMIN:
            value->rValue = model->BSIM4v7Lmin;
            return(OK);
        case BSIM4v7_MOD_LMAX:
            value->rValue = model->BSIM4v7Lmax;
            return(OK);
        case BSIM4v7_MOD_WINT:
            value->rValue = model->BSIM4v7Wint;
            return(OK);
        case BSIM4v7_MOD_WL:
            value->rValue = model->BSIM4v7Wl;
            return(OK);
        case BSIM4v7_MOD_WLC:
            value->rValue = model->BSIM4v7Wlc;
            return(OK);
        case BSIM4v7_MOD_WLN:
            value->rValue = model->BSIM4v7Wln;
            return(OK);
        case BSIM4v7_MOD_WW:
            value->rValue = model->BSIM4v7Ww;
            return(OK);
        case BSIM4v7_MOD_WWC:
            value->rValue = model->BSIM4v7Wwc;
            return(OK);
        case BSIM4v7_MOD_WWN:
            value->rValue = model->BSIM4v7Wwn;
            return(OK);
        case BSIM4v7_MOD_WWL:
            value->rValue = model->BSIM4v7Wwl;
            return(OK);
        case BSIM4v7_MOD_WWLC:
            value->rValue = model->BSIM4v7Wwlc;
            return(OK);
        case BSIM4v7_MOD_WMIN:
            value->rValue = model->BSIM4v7Wmin;
            return(OK);
        case BSIM4v7_MOD_WMAX:
            value->rValue = model->BSIM4v7Wmax;
            return(OK);

        /* stress effect */
        case BSIM4v7_MOD_SAREF:
            value->rValue = model->BSIM4v7saref;
            return(OK);
        case BSIM4v7_MOD_SBREF:
            value->rValue = model->BSIM4v7sbref;
            return(OK);
        case BSIM4v7_MOD_WLOD:
            value->rValue = model->BSIM4v7wlod;
            return(OK);
        case BSIM4v7_MOD_KU0:
            value->rValue = model->BSIM4v7ku0;
            return(OK);
        case BSIM4v7_MOD_KVSAT:
            value->rValue = model->BSIM4v7kvsat;
            return(OK);
        case BSIM4v7_MOD_KVTH0:
            value->rValue = model->BSIM4v7kvth0;
            return(OK);
        case BSIM4v7_MOD_TKU0:
            value->rValue = model->BSIM4v7tku0;
            return(OK);
        case BSIM4v7_MOD_LLODKU0:
            value->rValue = model->BSIM4v7llodku0;
            return(OK);
        case BSIM4v7_MOD_WLODKU0:
            value->rValue = model->BSIM4v7wlodku0;
            return(OK);
        case BSIM4v7_MOD_LLODVTH:
            value->rValue = model->BSIM4v7llodvth;
            return(OK);
        case BSIM4v7_MOD_WLODVTH:
            value->rValue = model->BSIM4v7wlodvth;
            return(OK);
        case BSIM4v7_MOD_LKU0:
            value->rValue = model->BSIM4v7lku0;
            return(OK);
        case BSIM4v7_MOD_WKU0:
            value->rValue = model->BSIM4v7wku0;
            return(OK);
        case BSIM4v7_MOD_PKU0:
            value->rValue = model->BSIM4v7pku0;
            return(OK);
        case BSIM4v7_MOD_LKVTH0:
            value->rValue = model->BSIM4v7lkvth0;
            return(OK);
        case BSIM4v7_MOD_WKVTH0:
            value->rValue = model->BSIM4v7wkvth0;
            return(OK);
        case BSIM4v7_MOD_PKVTH0:
            value->rValue = model->BSIM4v7pkvth0;
            return(OK);
        case BSIM4v7_MOD_STK2:
            value->rValue = model->BSIM4v7stk2;
            return(OK);
        case BSIM4v7_MOD_LODK2:
            value->rValue = model->BSIM4v7lodk2;
            return(OK);
        case BSIM4v7_MOD_STETA0:
            value->rValue = model->BSIM4v7steta0;
            return(OK);
        case BSIM4v7_MOD_LODETA0:
            value->rValue = model->BSIM4v7lodeta0;
            return(OK);

        /* Well Proximity Effect  */
        case BSIM4v7_MOD_WEB:
            value->rValue = model->BSIM4v7web;
            return(OK);
        case BSIM4v7_MOD_WEC:
            value->rValue = model->BSIM4v7wec;
            return(OK);
        case BSIM4v7_MOD_KVTH0WE:
            value->rValue = model->BSIM4v7kvth0we;
            return(OK);
        case BSIM4v7_MOD_K2WE:
            value->rValue = model->BSIM4v7k2we;
            return(OK);
        case BSIM4v7_MOD_KU0WE:
            value->rValue = model->BSIM4v7ku0we;
            return(OK);
        case BSIM4v7_MOD_SCREF:
            value->rValue = model->BSIM4v7scref;
            return(OK);
        case BSIM4v7_MOD_WPEMOD:
            value->rValue = model->BSIM4v7wpemod;
            return(OK);
        case BSIM4v7_MOD_LKVTH0WE:
            value->rValue = model->BSIM4v7lkvth0we;
            return(OK);
        case BSIM4v7_MOD_LK2WE:
            value->rValue = model->BSIM4v7lk2we;
            return(OK);
        case BSIM4v7_MOD_LKU0WE:
            value->rValue = model->BSIM4v7lku0we;
            return(OK);
        case BSIM4v7_MOD_WKVTH0WE:
            value->rValue = model->BSIM4v7wkvth0we;
            return(OK);
        case BSIM4v7_MOD_WK2WE:
            value->rValue = model->BSIM4v7wk2we;
            return(OK);
        case BSIM4v7_MOD_WKU0WE:
            value->rValue = model->BSIM4v7wku0we;
            return(OK);
        case BSIM4v7_MOD_PKVTH0WE:
            value->rValue = model->BSIM4v7pkvth0we;
            return(OK);
        case BSIM4v7_MOD_PK2WE:
            value->rValue = model->BSIM4v7pk2we;
            return(OK);
        case BSIM4v7_MOD_PKU0WE:
            value->rValue = model->BSIM4v7pku0we;
            return(OK);

        case BSIM4v7_MOD_NOIA:
            value->rValue = model->BSIM4v7oxideTrapDensityA;
            return(OK);
        case BSIM4v7_MOD_NOIB:
            value->rValue = model->BSIM4v7oxideTrapDensityB;
            return(OK);
        case BSIM4v7_MOD_NOIC:
            value->rValue = model->BSIM4v7oxideTrapDensityC;
            return(OK);
        case BSIM4v7_MOD_EM:
            value->rValue = model->BSIM4v7em;
            return(OK);
        case BSIM4v7_MOD_EF:
            value->rValue = model->BSIM4v7ef;
            return(OK);
        case BSIM4v7_MOD_AF:
            value->rValue = model->BSIM4v7af;
            return(OK);
        case BSIM4v7_MOD_KF:
            value->rValue = model->BSIM4v7kf;
            return(OK);

        case BSIM4v7_MOD_VGS_MAX:
            value->rValue = model->BSIM4v7vgsMax;
            return(OK);
        case BSIM4v7_MOD_VGD_MAX:
            value->rValue = model->BSIM4v7vgdMax;
            return(OK);
        case BSIM4v7_MOD_VGB_MAX:
            value->rValue = model->BSIM4v7vgbMax;
            return(OK);
        case BSIM4v7_MOD_VDS_MAX:
            value->rValue = model->BSIM4v7vdsMax;
            return(OK);
        case BSIM4v7_MOD_VBS_MAX:
            value->rValue = model->BSIM4v7vbsMax;
            return(OK);
        case BSIM4v7_MOD_VBD_MAX:
            value->rValue = model->BSIM4v7vbdMax;
            return(OK);

        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}



