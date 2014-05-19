/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4mask.c of BSIM4.8.0.
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
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4mAsk(
CKTcircuit *ckt,
GENmodel *inst,
int which,
IFvalue *value)
{
    BSIM4model *model = (BSIM4model *)inst;

    NG_IGNORE(ckt);

    switch(which) 
    {   case BSIM4_MOD_MOBMOD :
            value->iValue = model->BSIM4mobMod; 
            return(OK);
        case BSIM4_MOD_PARAMCHK :
            value->iValue = model->BSIM4paramChk; 
            return(OK);
        case BSIM4_MOD_BINUNIT :
            value->iValue = model->BSIM4binUnit; 
            return(OK);
        case BSIM4_MOD_CVCHARGEMOD :
            value->iValue = model->BSIM4cvchargeMod; 
            return(OK);
        case BSIM4_MOD_CAPMOD :
            value->iValue = model->BSIM4capMod; 
            return(OK);
        case BSIM4_MOD_DIOMOD :
            value->iValue = model->BSIM4dioMod;
            return(OK);
        case BSIM4_MOD_TRNQSMOD :
            value->iValue = model->BSIM4trnqsMod;
            return(OK);
        case BSIM4_MOD_ACNQSMOD :
            value->iValue = model->BSIM4acnqsMod;
            return(OK);
        case BSIM4_MOD_FNOIMOD :
            value->iValue = model->BSIM4fnoiMod; 
            return(OK);
        case BSIM4_MOD_TNOIMOD :
            value->iValue = model->BSIM4tnoiMod;
            return(OK);
        case BSIM4_MOD_RDSMOD :
            value->iValue = model->BSIM4rdsMod;
            return(OK);
        case BSIM4_MOD_RBODYMOD :
            value->iValue = model->BSIM4rbodyMod;
            return(OK);
        case BSIM4_MOD_RGATEMOD :
            value->iValue = model->BSIM4rgateMod;
            return(OK);
        case BSIM4_MOD_PERMOD :
            value->iValue = model->BSIM4perMod;
            return(OK);
        case BSIM4_MOD_GEOMOD :
            value->iValue = model->BSIM4geoMod;
            return(OK);
        case BSIM4_MOD_RGEOMOD :
            value->iValue = model->BSIM4rgeoMod;
            return(OK);
        case BSIM4_MOD_MTRLMOD :
            value->iValue = model->BSIM4mtrlMod;
            return(OK);
	case BSIM4_MOD_GIDLMOD :		/* v4.7 New GIDL/GISL*/
            value->iValue = model->BSIM4gidlMod;
            return(OK);
        case BSIM4_MOD_MTRLCOMPATMOD :
	    value->iValue = model->BSIM4mtrlCompatMod;
	    return(OK);
        case BSIM4_MOD_IGCMOD :
            value->iValue = model->BSIM4igcMod;
            return(OK);
        case BSIM4_MOD_IGBMOD :
            value->iValue = model->BSIM4igbMod;
            return(OK);
        case  BSIM4_MOD_TEMPMOD :
            value->iValue = model->BSIM4tempMod;
            return(OK);

        case  BSIM4_MOD_VERSION :
          value->sValue = model->BSIM4version;
            return(OK);
        case  BSIM4_MOD_TOXREF :
          value->rValue = model->BSIM4toxref;
          return(OK);
        case  BSIM4_MOD_EOT :
          value->rValue = model->BSIM4eot;
            return(OK);
        case  BSIM4_MOD_VDDEOT :
          value->rValue = model->BSIM4vddeot;
            return(OK);
		case  BSIM4_MOD_TEMPEOT :
          value->rValue = model->BSIM4tempeot;
            return(OK);
		case  BSIM4_MOD_LEFFEOT :
          value->rValue = model->BSIM4leffeot;
            return(OK);
		case  BSIM4_MOD_WEFFEOT :
          value->rValue = model->BSIM4weffeot;
            return(OK);
        case  BSIM4_MOD_ADOS :
          value->rValue = model->BSIM4ados;
            return(OK);
        case  BSIM4_MOD_BDOS :
          value->rValue = model->BSIM4bdos;
            return(OK);
        case  BSIM4_MOD_TOXE :
          value->rValue = model->BSIM4toxe;
            return(OK);
        case  BSIM4_MOD_TOXP :
          value->rValue = model->BSIM4toxp;
            return(OK);
        case  BSIM4_MOD_TOXM :
          value->rValue = model->BSIM4toxm;
            return(OK);
        case  BSIM4_MOD_DTOX :
          value->rValue = model->BSIM4dtox;
            return(OK);
        case  BSIM4_MOD_EPSROX :
          value->rValue = model->BSIM4epsrox;
            return(OK);
        case  BSIM4_MOD_CDSC :
          value->rValue = model->BSIM4cdsc;
            return(OK);
        case  BSIM4_MOD_CDSCB :
          value->rValue = model->BSIM4cdscb;
            return(OK);

        case  BSIM4_MOD_CDSCD :
          value->rValue = model->BSIM4cdscd;
            return(OK);

        case  BSIM4_MOD_CIT :
          value->rValue = model->BSIM4cit;
            return(OK);
        case  BSIM4_MOD_NFACTOR :
          value->rValue = model->BSIM4nfactor;
            return(OK);
        case BSIM4_MOD_XJ:
            value->rValue = model->BSIM4xj;
            return(OK);
        case BSIM4_MOD_VSAT:
            value->rValue = model->BSIM4vsat;
            return(OK);
        case BSIM4_MOD_VTL:
            value->rValue = model->BSIM4vtl;
            return(OK);
        case BSIM4_MOD_XN:
            value->rValue = model->BSIM4xn;
            return(OK);
        case BSIM4_MOD_LC:
            value->rValue = model->BSIM4lc;
            return(OK);
        case BSIM4_MOD_LAMBDA:
            value->rValue = model->BSIM4lambda;
            return(OK);
        case BSIM4_MOD_AT:
            value->rValue = model->BSIM4at;
            return(OK);
        case BSIM4_MOD_A0:
            value->rValue = model->BSIM4a0;
            return(OK);

        case BSIM4_MOD_AGS:
            value->rValue = model->BSIM4ags;
            return(OK);

        case BSIM4_MOD_A1:
            value->rValue = model->BSIM4a1;
            return(OK);
        case BSIM4_MOD_A2:
            value->rValue = model->BSIM4a2;
            return(OK);
        case BSIM4_MOD_KETA:
            value->rValue = model->BSIM4keta;
            return(OK);   
        case BSIM4_MOD_NSUB:
            value->rValue = model->BSIM4nsub;
            return(OK);
        case BSIM4_MOD_PHIG:
	    value->rValue = model->BSIM4phig;
	    return(OK);
        case BSIM4_MOD_EPSRGATE:
	    value->rValue = model->BSIM4epsrgate;
	    return(OK);
        case BSIM4_MOD_EASUB:
            value->rValue = model->BSIM4easub;
            return(OK);
        case BSIM4_MOD_EPSRSUB:
            value->rValue = model->BSIM4epsrsub;
            return(OK);
        case BSIM4_MOD_NI0SUB:
            value->rValue = model->BSIM4ni0sub;
            return(OK);
        case BSIM4_MOD_BG0SUB:
            value->rValue = model->BSIM4bg0sub;
            return(OK);
        case BSIM4_MOD_TBGASUB:
            value->rValue = model->BSIM4tbgasub;
            return(OK);
        case BSIM4_MOD_TBGBSUB:
            value->rValue = model->BSIM4tbgbsub;
            return(OK);
        case BSIM4_MOD_NDEP:
            value->rValue = model->BSIM4ndep;
            return(OK);
        case BSIM4_MOD_NSD:
            value->rValue = model->BSIM4nsd;
            return(OK);
        case BSIM4_MOD_NGATE:
            value->rValue = model->BSIM4ngate;
            return(OK);
        case BSIM4_MOD_GAMMA1:
            value->rValue = model->BSIM4gamma1;
            return(OK);
        case BSIM4_MOD_GAMMA2:
            value->rValue = model->BSIM4gamma2;
            return(OK);
        case BSIM4_MOD_VBX:
            value->rValue = model->BSIM4vbx;
            return(OK);
        case BSIM4_MOD_VBM:
            value->rValue = model->BSIM4vbm;
            return(OK);
        case BSIM4_MOD_XT:
            value->rValue = model->BSIM4xt;
            return(OK);
        case  BSIM4_MOD_K1:
          value->rValue = model->BSIM4k1;
            return(OK);
        case  BSIM4_MOD_KT1:
          value->rValue = model->BSIM4kt1;
            return(OK);
        case  BSIM4_MOD_KT1L:
          value->rValue = model->BSIM4kt1l;
            return(OK);
        case  BSIM4_MOD_KT2 :
          value->rValue = model->BSIM4kt2;
            return(OK);
        case  BSIM4_MOD_K2 :
          value->rValue = model->BSIM4k2;
            return(OK);
        case  BSIM4_MOD_K3:
          value->rValue = model->BSIM4k3;
            return(OK);
        case  BSIM4_MOD_K3B:
          value->rValue = model->BSIM4k3b;
            return(OK);
        case  BSIM4_MOD_W0:
          value->rValue = model->BSIM4w0;
            return(OK);
        case  BSIM4_MOD_LPE0:
          value->rValue = model->BSIM4lpe0;
            return(OK);
        case  BSIM4_MOD_LPEB:
          value->rValue = model->BSIM4lpeb;
            return(OK);
        case  BSIM4_MOD_DVTP0:
          value->rValue = model->BSIM4dvtp0;
            return(OK);
        case  BSIM4_MOD_DVTP1:
          value->rValue = model->BSIM4dvtp1;
            return(OK);
        case  BSIM4_MOD_DVTP2:
          value->rValue = model->BSIM4dvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4_MOD_DVTP3:
          value->rValue = model->BSIM4dvtp3;
            return(OK);
        case  BSIM4_MOD_DVTP4:
          value->rValue = model->BSIM4dvtp4;
            return(OK);
        case  BSIM4_MOD_DVTP5:
          value->rValue = model->BSIM4dvtp5;
            return(OK);
        case  BSIM4_MOD_DVT0 :                
          value->rValue = model->BSIM4dvt0;
            return(OK);
        case  BSIM4_MOD_DVT1 :             
          value->rValue = model->BSIM4dvt1;
            return(OK);
        case  BSIM4_MOD_DVT2 :             
          value->rValue = model->BSIM4dvt2;
            return(OK);
        case  BSIM4_MOD_DVT0W :                
          value->rValue = model->BSIM4dvt0w;
            return(OK);
        case  BSIM4_MOD_DVT1W :             
          value->rValue = model->BSIM4dvt1w;
            return(OK);
        case  BSIM4_MOD_DVT2W :             
          value->rValue = model->BSIM4dvt2w;
            return(OK);
        case  BSIM4_MOD_DROUT :           
          value->rValue = model->BSIM4drout;
            return(OK);
        case  BSIM4_MOD_DSUB :           
          value->rValue = model->BSIM4dsub;
            return(OK);
        case BSIM4_MOD_VTH0:
            value->rValue = model->BSIM4vth0; 
            return(OK);
        case BSIM4_MOD_EU:
            value->rValue = model->BSIM4eu;
            return(OK);
		 case BSIM4_MOD_UCS:
            value->rValue = model->BSIM4ucs;
            return(OK);
        case BSIM4_MOD_UA:
            value->rValue = model->BSIM4ua; 
            return(OK);
        case BSIM4_MOD_UA1:
            value->rValue = model->BSIM4ua1; 
            return(OK);
        case BSIM4_MOD_UB:
            value->rValue = model->BSIM4ub;  
            return(OK);
        case BSIM4_MOD_UB1:
            value->rValue = model->BSIM4ub1;  
            return(OK);
        case BSIM4_MOD_UC:
            value->rValue = model->BSIM4uc; 
            return(OK);
        case BSIM4_MOD_UC1:
            value->rValue = model->BSIM4uc1; 
            return(OK);
        case BSIM4_MOD_UD:
            value->rValue = model->BSIM4ud; 
            return(OK);
        case BSIM4_MOD_UD1:
            value->rValue = model->BSIM4ud1; 
            return(OK);
        case BSIM4_MOD_UP:
            value->rValue = model->BSIM4up; 
            return(OK);
        case BSIM4_MOD_LP:
            value->rValue = model->BSIM4lp; 
            return(OK);
        case BSIM4_MOD_U0:
            value->rValue = model->BSIM4u0;
            return(OK);
        case BSIM4_MOD_UTE:
            value->rValue = model->BSIM4ute;
            return(OK);
		 case BSIM4_MOD_UCSTE:
            value->rValue = model->BSIM4ucste;
            return(OK);
        case BSIM4_MOD_VOFF:
            value->rValue = model->BSIM4voff;
            return(OK);
        case BSIM4_MOD_TVOFF:
            value->rValue = model->BSIM4tvoff;
            return(OK);
        case BSIM4_MOD_TNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4tnfactor;
            return(OK);
        case BSIM4_MOD_TETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4teta0;
            return(OK);
        case BSIM4_MOD_TVOFFCV:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4tvoffcv;
            return(OK);
        case BSIM4_MOD_VFBSDOFF:
            value->rValue = model->BSIM4vfbsdoff;
            return(OK);
        case BSIM4_MOD_TVFBSDOFF:
            value->rValue = model->BSIM4tvfbsdoff;
            return(OK);
        case BSIM4_MOD_VOFFL:
            value->rValue = model->BSIM4voffl;
            return(OK);
        case BSIM4_MOD_VOFFCVL:
            value->rValue = model->BSIM4voffcvl;
            return(OK);
        case BSIM4_MOD_MINV:
            value->rValue = model->BSIM4minv;
            return(OK);
        case BSIM4_MOD_MINVCV:
            value->rValue = model->BSIM4minvcv;
            return(OK);
        case BSIM4_MOD_FPROUT:
            value->rValue = model->BSIM4fprout;
            return(OK);
        case BSIM4_MOD_PDITS:
            value->rValue = model->BSIM4pdits;
            return(OK);
        case BSIM4_MOD_PDITSD:
            value->rValue = model->BSIM4pditsd;
            return(OK);
        case BSIM4_MOD_PDITSL:
            value->rValue = model->BSIM4pditsl;
            return(OK);
        case BSIM4_MOD_DELTA:
            value->rValue = model->BSIM4delta;
            return(OK);
        case BSIM4_MOD_RDSW:
            value->rValue = model->BSIM4rdsw; 
            return(OK);
        case BSIM4_MOD_RDSWMIN:
            value->rValue = model->BSIM4rdswmin;
            return(OK);
        case BSIM4_MOD_RDWMIN:
            value->rValue = model->BSIM4rdwmin;
            return(OK);
        case BSIM4_MOD_RSWMIN:
            value->rValue = model->BSIM4rswmin;
            return(OK);
        case BSIM4_MOD_RDW:
            value->rValue = model->BSIM4rdw;
            return(OK);
        case BSIM4_MOD_RSW:
            value->rValue = model->BSIM4rsw;
            return(OK);
        case BSIM4_MOD_PRWG:
            value->rValue = model->BSIM4prwg; 
            return(OK);             
        case BSIM4_MOD_PRWB:
            value->rValue = model->BSIM4prwb; 
            return(OK);             
        case BSIM4_MOD_PRT:
            value->rValue = model->BSIM4prt; 
            return(OK);              
        case BSIM4_MOD_ETA0:
            value->rValue = model->BSIM4eta0; 
            return(OK);               
        case BSIM4_MOD_ETAB:
            value->rValue = model->BSIM4etab; 
            return(OK);               
        case BSIM4_MOD_PCLM:
            value->rValue = model->BSIM4pclm; 
            return(OK);               
        case BSIM4_MOD_PDIBL1:
            value->rValue = model->BSIM4pdibl1; 
            return(OK);               
        case BSIM4_MOD_PDIBL2:
            value->rValue = model->BSIM4pdibl2; 
            return(OK);               
        case BSIM4_MOD_PDIBLB:
            value->rValue = model->BSIM4pdiblb; 
            return(OK);               
        case BSIM4_MOD_PSCBE1:
            value->rValue = model->BSIM4pscbe1; 
            return(OK);               
        case BSIM4_MOD_PSCBE2:
            value->rValue = model->BSIM4pscbe2; 
            return(OK);               
        case BSIM4_MOD_PVAG:
            value->rValue = model->BSIM4pvag; 
            return(OK);               
        case BSIM4_MOD_WR:
            value->rValue = model->BSIM4wr;
            return(OK);
        case BSIM4_MOD_DWG:
            value->rValue = model->BSIM4dwg;
            return(OK);
        case BSIM4_MOD_DWB:
            value->rValue = model->BSIM4dwb;
            return(OK);
        case BSIM4_MOD_B0:
            value->rValue = model->BSIM4b0;
            return(OK);
        case BSIM4_MOD_B1:
            value->rValue = model->BSIM4b1;
            return(OK);
        case BSIM4_MOD_ALPHA0:
            value->rValue = model->BSIM4alpha0;
            return(OK);
        case BSIM4_MOD_ALPHA1:
            value->rValue = model->BSIM4alpha1;
            return(OK);
        case BSIM4_MOD_BETA0:
            value->rValue = model->BSIM4beta0;
            return(OK);
        case BSIM4_MOD_AGIDL:
            value->rValue = model->BSIM4agidl;
            return(OK);
        case BSIM4_MOD_BGIDL:
            value->rValue = model->BSIM4bgidl;
            return(OK);
        case BSIM4_MOD_CGIDL:
            value->rValue = model->BSIM4cgidl;
            return(OK);
        case BSIM4_MOD_EGIDL:
            value->rValue = model->BSIM4egidl;
            return(OK);
 	case BSIM4_MOD_FGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4fgidl;
            return(OK);
 	case BSIM4_MOD_KGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4kgidl;
            return(OK);
 	case BSIM4_MOD_RGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4rgidl;
            return(OK);
        case BSIM4_MOD_AGISL:
            value->rValue = model->BSIM4agisl;
            return(OK);
        case BSIM4_MOD_BGISL:
            value->rValue = model->BSIM4bgisl;
            return(OK);
        case BSIM4_MOD_CGISL:
            value->rValue = model->BSIM4cgisl;
            return(OK);
        case BSIM4_MOD_EGISL:
            value->rValue = model->BSIM4egisl;
            return(OK);
 	case BSIM4_MOD_FGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4fgisl;
            return(OK);
 	case BSIM4_MOD_KGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4kgisl;
            return(OK);
 	case BSIM4_MOD_RGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4rgisl;
            return(OK);
        case BSIM4_MOD_AIGC:
            value->rValue = model->BSIM4aigc;
            return(OK);
        case BSIM4_MOD_BIGC:
            value->rValue = model->BSIM4bigc;
            return(OK);
        case BSIM4_MOD_CIGC:
            value->rValue = model->BSIM4cigc;
            return(OK);
        case BSIM4_MOD_AIGSD:
            value->rValue = model->BSIM4aigsd;
            return(OK);
        case BSIM4_MOD_BIGSD:
            value->rValue = model->BSIM4bigsd;
            return(OK);
        case BSIM4_MOD_CIGSD:
            value->rValue = model->BSIM4cigsd;
            return(OK);
        case BSIM4_MOD_AIGS:
            value->rValue = model->BSIM4aigs;
            return(OK);
        case BSIM4_MOD_BIGS:
            value->rValue = model->BSIM4bigs;
            return(OK);
        case BSIM4_MOD_CIGS:
            value->rValue = model->BSIM4cigs;
            return(OK);
        case BSIM4_MOD_AIGD:
            value->rValue = model->BSIM4aigd;
            return(OK);
        case BSIM4_MOD_BIGD:
            value->rValue = model->BSIM4bigd;
            return(OK);
        case BSIM4_MOD_CIGD:
            value->rValue = model->BSIM4cigd;
            return(OK);
        case BSIM4_MOD_AIGBACC:
            value->rValue = model->BSIM4aigbacc;
            return(OK);
        case BSIM4_MOD_BIGBACC:
            value->rValue = model->BSIM4bigbacc;
            return(OK);
        case BSIM4_MOD_CIGBACC:
            value->rValue = model->BSIM4cigbacc;
            return(OK);
        case BSIM4_MOD_AIGBINV:
            value->rValue = model->BSIM4aigbinv;
            return(OK);
        case BSIM4_MOD_BIGBINV:
            value->rValue = model->BSIM4bigbinv;
            return(OK);
        case BSIM4_MOD_CIGBINV:
            value->rValue = model->BSIM4cigbinv;
            return(OK);
        case BSIM4_MOD_NIGC:
            value->rValue = model->BSIM4nigc;
            return(OK);
        case BSIM4_MOD_NIGBACC:
            value->rValue = model->BSIM4nigbacc;
            return(OK);
        case BSIM4_MOD_NIGBINV:
            value->rValue = model->BSIM4nigbinv;
            return(OK);
        case BSIM4_MOD_NTOX:
            value->rValue = model->BSIM4ntox;
            return(OK);
        case BSIM4_MOD_EIGBINV:
            value->rValue = model->BSIM4eigbinv;
            return(OK);
        case BSIM4_MOD_PIGCD:
            value->rValue = model->BSIM4pigcd;
            return(OK);
        case BSIM4_MOD_POXEDGE:
            value->rValue = model->BSIM4poxedge;
            return(OK);
        case BSIM4_MOD_PHIN:
            value->rValue = model->BSIM4phin;
            return(OK);
        case BSIM4_MOD_XRCRG1:
            value->rValue = model->BSIM4xrcrg1;
            return(OK);
        case BSIM4_MOD_XRCRG2:
            value->rValue = model->BSIM4xrcrg2;
            return(OK);
        case BSIM4_MOD_TNOIA:
            value->rValue = model->BSIM4tnoia;
            return(OK);
        case BSIM4_MOD_TNOIB:
            value->rValue = model->BSIM4tnoib;
            return(OK);
        case BSIM4_MOD_TNOIC:
            value->rValue = model->BSIM4tnoic;
            return(OK);
        case BSIM4_MOD_RNOIA:
            value->rValue = model->BSIM4rnoia;
            return(OK);
        case BSIM4_MOD_RNOIB:
            value->rValue = model->BSIM4rnoib;
            return(OK);
        case BSIM4_MOD_RNOIC:
            value->rValue = model->BSIM4rnoic;
            return(OK);
        case BSIM4_MOD_NTNOI:
            value->rValue = model->BSIM4ntnoi;
            return(OK);
        case BSIM4_MOD_IJTHDFWD:
            value->rValue = model->BSIM4ijthdfwd;
            return(OK);
        case BSIM4_MOD_IJTHSFWD:
            value->rValue = model->BSIM4ijthsfwd;
            return(OK);
        case BSIM4_MOD_IJTHDREV:
            value->rValue = model->BSIM4ijthdrev;
            return(OK);
        case BSIM4_MOD_IJTHSREV:
            value->rValue = model->BSIM4ijthsrev;
            return(OK);
        case BSIM4_MOD_XJBVD:
            value->rValue = model->BSIM4xjbvd;
            return(OK);
        case BSIM4_MOD_XJBVS:
            value->rValue = model->BSIM4xjbvs;
            return(OK);
        case BSIM4_MOD_BVD:
            value->rValue = model->BSIM4bvd;
            return(OK);
        case BSIM4_MOD_BVS:
            value->rValue = model->BSIM4bvs;
            return(OK);
        case BSIM4_MOD_VFB:
            value->rValue = model->BSIM4vfb;
            return(OK);

        case BSIM4_MOD_JTSS:
            value->rValue = model->BSIM4jtss;
            return(OK);
        case BSIM4_MOD_JTSD:
            value->rValue = model->BSIM4jtsd;
            return(OK);
        case BSIM4_MOD_JTSSWS:
            value->rValue = model->BSIM4jtssws;
            return(OK);
        case BSIM4_MOD_JTSSWD:
            value->rValue = model->BSIM4jtsswd;
            return(OK);
        case BSIM4_MOD_JTSSWGS:
            value->rValue = model->BSIM4jtsswgs;
            return(OK);
        case BSIM4_MOD_JTSSWGD:
            value->rValue = model->BSIM4jtsswgd;
            return(OK);
		case BSIM4_MOD_JTWEFF:
		    value->rValue = model->BSIM4jtweff;
			return(OK);
        case BSIM4_MOD_NJTS:
            value->rValue = model->BSIM4njts;
            return(OK);
        case BSIM4_MOD_NJTSSW:
            value->rValue = model->BSIM4njtssw;
            return(OK);
        case BSIM4_MOD_NJTSSWG:
            value->rValue = model->BSIM4njtsswg;
            return(OK);
        case BSIM4_MOD_NJTSD:
            value->rValue = model->BSIM4njtsd;
            return(OK);
        case BSIM4_MOD_NJTSSWD:
            value->rValue = model->BSIM4njtsswd;
            return(OK);
        case BSIM4_MOD_NJTSSWGD:
            value->rValue = model->BSIM4njtsswgd;
            return(OK);
        case BSIM4_MOD_XTSS:
            value->rValue = model->BSIM4xtss;
            return(OK);
        case BSIM4_MOD_XTSD:
            value->rValue = model->BSIM4xtsd;
            return(OK);
        case BSIM4_MOD_XTSSWS:
            value->rValue = model->BSIM4xtssws;
            return(OK);
        case BSIM4_MOD_XTSSWD:
            value->rValue = model->BSIM4xtsswd;
            return(OK);
        case BSIM4_MOD_XTSSWGS:
            value->rValue = model->BSIM4xtsswgs;
            return(OK);
        case BSIM4_MOD_XTSSWGD:
            value->rValue = model->BSIM4xtsswgd;
            return(OK);
        case BSIM4_MOD_TNJTS:
            value->rValue = model->BSIM4tnjts;
            return(OK);
        case BSIM4_MOD_TNJTSSW:
            value->rValue = model->BSIM4tnjtssw;
            return(OK);
        case BSIM4_MOD_TNJTSSWG:
            value->rValue = model->BSIM4tnjtsswg;
            return(OK);
        case BSIM4_MOD_TNJTSD:
            value->rValue = model->BSIM4tnjtsd;
            return(OK);
        case BSIM4_MOD_TNJTSSWD:
            value->rValue = model->BSIM4tnjtsswd;
            return(OK);
        case BSIM4_MOD_TNJTSSWGD:
            value->rValue = model->BSIM4tnjtsswgd;
            return(OK);
        case BSIM4_MOD_VTSS:
            value->rValue = model->BSIM4vtss;
            return(OK);
        case BSIM4_MOD_VTSD:
            value->rValue = model->BSIM4vtsd;
            return(OK);
        case BSIM4_MOD_VTSSWS:
            value->rValue = model->BSIM4vtssws;
            return(OK);
        case BSIM4_MOD_VTSSWD:
            value->rValue = model->BSIM4vtsswd;
            return(OK);
        case BSIM4_MOD_VTSSWGS:
            value->rValue = model->BSIM4vtsswgs;
            return(OK);
        case BSIM4_MOD_VTSSWGD:
            value->rValue = model->BSIM4vtsswgd;
            return(OK);

        case BSIM4_MOD_GBMIN:
            value->rValue = model->BSIM4gbmin;
            return(OK);
        case BSIM4_MOD_RBDB:
            value->rValue = model->BSIM4rbdb;
            return(OK);
        case BSIM4_MOD_RBPB:
            value->rValue = model->BSIM4rbpb;
            return(OK);
        case BSIM4_MOD_RBSB:
            value->rValue = model->BSIM4rbsb;
            return(OK);
        case BSIM4_MOD_RBPS:
            value->rValue = model->BSIM4rbps;
            return(OK);
        case BSIM4_MOD_RBPD:
            value->rValue = model->BSIM4rbpd;
            return(OK);

        case BSIM4_MOD_RBPS0:
            value->rValue = model->BSIM4rbps0;
            return(OK);
        case BSIM4_MOD_RBPSL:
            value->rValue = model->BSIM4rbpsl;
            return(OK);
        case BSIM4_MOD_RBPSW:
            value->rValue = model->BSIM4rbpsw;
            return(OK);
        case BSIM4_MOD_RBPSNF:
            value->rValue = model->BSIM4rbpsnf;
            return(OK);
        case BSIM4_MOD_RBPD0:
            value->rValue = model->BSIM4rbpd0;
            return(OK);
        case BSIM4_MOD_RBPDL:
            value->rValue = model->BSIM4rbpdl;
            return(OK);
        case BSIM4_MOD_RBPDW:
            value->rValue = model->BSIM4rbpdw;
            return(OK);
        case BSIM4_MOD_RBPDNF:
            value->rValue = model->BSIM4rbpdnf;
            return(OK);
        case BSIM4_MOD_RBPBX0:
            value->rValue = model->BSIM4rbpbx0;
            return(OK);
        case BSIM4_MOD_RBPBXL:
            value->rValue = model->BSIM4rbpbxl;
            return(OK);
        case BSIM4_MOD_RBPBXW:
            value->rValue = model->BSIM4rbpbxw;
            return(OK);
        case BSIM4_MOD_RBPBXNF:
            value->rValue = model->BSIM4rbpbxnf;
            return(OK);
        case BSIM4_MOD_RBPBY0:
            value->rValue = model->BSIM4rbpby0;
            return(OK);
        case BSIM4_MOD_RBPBYL:
            value->rValue = model->BSIM4rbpbyl;
            return(OK);
        case BSIM4_MOD_RBPBYW:
            value->rValue = model->BSIM4rbpbyw;
            return(OK);
        case BSIM4_MOD_RBPBYNF:
            value->rValue = model->BSIM4rbpbynf;
            return(OK);

        case BSIM4_MOD_RBSBX0:
            value->rValue = model->BSIM4rbsbx0;
            return(OK);
        case BSIM4_MOD_RBSBY0:
            value->rValue = model->BSIM4rbsby0;
            return(OK);
        case BSIM4_MOD_RBDBX0:
            value->rValue = model->BSIM4rbdbx0;
            return(OK);
        case BSIM4_MOD_RBDBY0:
            value->rValue = model->BSIM4rbdby0;
            return(OK);
        case BSIM4_MOD_RBSDBXL:
            value->rValue = model->BSIM4rbsdbxl;
            return(OK);
        case BSIM4_MOD_RBSDBXW:
            value->rValue = model->BSIM4rbsdbxw;
            return(OK);
        case BSIM4_MOD_RBSDBXNF:
            value->rValue = model->BSIM4rbsdbxnf;
            return(OK);
        case BSIM4_MOD_RBSDBYL:
            value->rValue = model->BSIM4rbsdbyl;
            return(OK);
        case BSIM4_MOD_RBSDBYW:
            value->rValue = model->BSIM4rbsdbyw;
            return(OK);
        case BSIM4_MOD_RBSDBYNF:
            value->rValue = model->BSIM4rbsdbynf;
            return(OK);


        case BSIM4_MOD_CGSL:
            value->rValue = model->BSIM4cgsl;
            return(OK);
        case BSIM4_MOD_CGDL:
            value->rValue = model->BSIM4cgdl;
            return(OK);
        case BSIM4_MOD_CKAPPAS:
            value->rValue = model->BSIM4ckappas;
            return(OK);
        case BSIM4_MOD_CKAPPAD:
            value->rValue = model->BSIM4ckappad;
            return(OK);
        case BSIM4_MOD_CF:
            value->rValue = model->BSIM4cf;
            return(OK);
        case BSIM4_MOD_CLC:
            value->rValue = model->BSIM4clc;
            return(OK);
        case BSIM4_MOD_CLE:
            value->rValue = model->BSIM4cle;
            return(OK);
        case BSIM4_MOD_DWC:
            value->rValue = model->BSIM4dwc;
            return(OK);
        case BSIM4_MOD_DLC:
            value->rValue = model->BSIM4dlc;
            return(OK);
        case BSIM4_MOD_XW:
            value->rValue = model->BSIM4xw;
            return(OK);
        case BSIM4_MOD_XL:
            value->rValue = model->BSIM4xl;
            return(OK);
        case BSIM4_MOD_DLCIG:
            value->rValue = model->BSIM4dlcig;
            return(OK);
        case BSIM4_MOD_DLCIGD:
            value->rValue = model->BSIM4dlcigd;
            return(OK);
        case BSIM4_MOD_DWJ:
            value->rValue = model->BSIM4dwj;
            return(OK);
        case BSIM4_MOD_VFBCV:
            value->rValue = model->BSIM4vfbcv; 
            return(OK);
        case BSIM4_MOD_ACDE:
            value->rValue = model->BSIM4acde;
            return(OK);
        case BSIM4_MOD_MOIN:
            value->rValue = model->BSIM4moin;
            return(OK);
        case BSIM4_MOD_NOFF:
            value->rValue = model->BSIM4noff;
            return(OK);
        case BSIM4_MOD_VOFFCV:
            value->rValue = model->BSIM4voffcv;
            return(OK);
        case BSIM4_MOD_DMCG:
            value->rValue = model->BSIM4dmcg;
            return(OK);
        case BSIM4_MOD_DMCI:
            value->rValue = model->BSIM4dmci;
            return(OK);
        case BSIM4_MOD_DMDG:
            value->rValue = model->BSIM4dmdg;
            return(OK);
        case BSIM4_MOD_DMCGT:
            value->rValue = model->BSIM4dmcgt;
            return(OK);
        case BSIM4_MOD_XGW:
            value->rValue = model->BSIM4xgw;
            return(OK);
        case BSIM4_MOD_XGL:
            value->rValue = model->BSIM4xgl;
            return(OK);
        case BSIM4_MOD_RSHG:
            value->rValue = model->BSIM4rshg;
            return(OK);
        case BSIM4_MOD_NGCON:
            value->rValue = model->BSIM4ngcon; 
            return(OK);
        case BSIM4_MOD_TCJ:
            value->rValue = model->BSIM4tcj;
            return(OK);
        case BSIM4_MOD_TPB:
            value->rValue = model->BSIM4tpb;
            return(OK);
        case BSIM4_MOD_TCJSW:
            value->rValue = model->BSIM4tcjsw;
            return(OK);
        case BSIM4_MOD_TPBSW:
            value->rValue = model->BSIM4tpbsw;
            return(OK);
        case BSIM4_MOD_TCJSWG:
            value->rValue = model->BSIM4tcjswg;
            return(OK);
        case BSIM4_MOD_TPBSWG:
            value->rValue = model->BSIM4tpbswg;
            return(OK);

	/* Length dependence */
        case  BSIM4_MOD_LCDSC :
          value->rValue = model->BSIM4lcdsc;
            return(OK);
        case  BSIM4_MOD_LCDSCB :
          value->rValue = model->BSIM4lcdscb;
            return(OK);
        case  BSIM4_MOD_LCDSCD :
          value->rValue = model->BSIM4lcdscd;
            return(OK);
        case  BSIM4_MOD_LCIT :
          value->rValue = model->BSIM4lcit;
            return(OK);
        case  BSIM4_MOD_LNFACTOR :
          value->rValue = model->BSIM4lnfactor;
            return(OK);
        case BSIM4_MOD_LXJ:
            value->rValue = model->BSIM4lxj;
            return(OK);
        case BSIM4_MOD_LVSAT:
            value->rValue = model->BSIM4lvsat;
            return(OK);
        case BSIM4_MOD_LAT:
            value->rValue = model->BSIM4lat;
            return(OK);
        case BSIM4_MOD_LA0:
            value->rValue = model->BSIM4la0;
            return(OK);
        case BSIM4_MOD_LAGS:
            value->rValue = model->BSIM4lags;
            return(OK);
        case BSIM4_MOD_LA1:
            value->rValue = model->BSIM4la1;
            return(OK);
        case BSIM4_MOD_LA2:
            value->rValue = model->BSIM4la2;
            return(OK);
        case BSIM4_MOD_LKETA:
            value->rValue = model->BSIM4lketa;
            return(OK);   
        case BSIM4_MOD_LNSUB:
            value->rValue = model->BSIM4lnsub;
            return(OK);
        case BSIM4_MOD_LNDEP:
            value->rValue = model->BSIM4lndep;
            return(OK);
        case BSIM4_MOD_LNSD:
            value->rValue = model->BSIM4lnsd;
            return(OK);
        case BSIM4_MOD_LNGATE:
            value->rValue = model->BSIM4lngate;
            return(OK);
        case BSIM4_MOD_LGAMMA1:
            value->rValue = model->BSIM4lgamma1;
            return(OK);
        case BSIM4_MOD_LGAMMA2:
            value->rValue = model->BSIM4lgamma2;
            return(OK);
        case BSIM4_MOD_LVBX:
            value->rValue = model->BSIM4lvbx;
            return(OK);
        case BSIM4_MOD_LVBM:
            value->rValue = model->BSIM4lvbm;
            return(OK);
        case BSIM4_MOD_LXT:
            value->rValue = model->BSIM4lxt;
            return(OK);
        case  BSIM4_MOD_LK1:
          value->rValue = model->BSIM4lk1;
            return(OK);
        case  BSIM4_MOD_LKT1:
          value->rValue = model->BSIM4lkt1;
            return(OK);
        case  BSIM4_MOD_LKT1L:
          value->rValue = model->BSIM4lkt1l;
            return(OK);
        case  BSIM4_MOD_LKT2 :
          value->rValue = model->BSIM4lkt2;
            return(OK);
        case  BSIM4_MOD_LK2 :
          value->rValue = model->BSIM4lk2;
            return(OK);
        case  BSIM4_MOD_LK3:
          value->rValue = model->BSIM4lk3;
            return(OK);
        case  BSIM4_MOD_LK3B:
          value->rValue = model->BSIM4lk3b;
            return(OK);
        case  BSIM4_MOD_LW0:
          value->rValue = model->BSIM4lw0;
            return(OK);
        case  BSIM4_MOD_LLPE0:
          value->rValue = model->BSIM4llpe0;
            return(OK);
        case  BSIM4_MOD_LLPEB:
          value->rValue = model->BSIM4llpeb;
            return(OK);
        case  BSIM4_MOD_LDVTP0:
          value->rValue = model->BSIM4ldvtp0;
            return(OK);
        case  BSIM4_MOD_LDVTP1:
          value->rValue = model->BSIM4ldvtp1;
            return(OK);
	case  BSIM4_MOD_LDVTP2:
          value->rValue = model->BSIM4ldvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4_MOD_LDVTP3:
          value->rValue = model->BSIM4ldvtp3;
            return(OK);
        case  BSIM4_MOD_LDVTP4:
          value->rValue = model->BSIM4ldvtp4;
            return(OK);
        case  BSIM4_MOD_LDVTP5:
          value->rValue = model->BSIM4ldvtp5;
            return(OK);
        case  BSIM4_MOD_LDVT0:                
          value->rValue = model->BSIM4ldvt0;
            return(OK);
        case  BSIM4_MOD_LDVT1 :             
          value->rValue = model->BSIM4ldvt1;
            return(OK);
        case  BSIM4_MOD_LDVT2 :             
          value->rValue = model->BSIM4ldvt2;
            return(OK);
        case  BSIM4_MOD_LDVT0W :                
          value->rValue = model->BSIM4ldvt0w;
            return(OK);
        case  BSIM4_MOD_LDVT1W :             
          value->rValue = model->BSIM4ldvt1w;
            return(OK);
        case  BSIM4_MOD_LDVT2W :             
          value->rValue = model->BSIM4ldvt2w;
            return(OK);
        case  BSIM4_MOD_LDROUT :           
          value->rValue = model->BSIM4ldrout;
            return(OK);
        case  BSIM4_MOD_LDSUB :           
          value->rValue = model->BSIM4ldsub;
            return(OK);
        case BSIM4_MOD_LVTH0:
            value->rValue = model->BSIM4lvth0; 
            return(OK);
        case BSIM4_MOD_LUA:
            value->rValue = model->BSIM4lua; 
            return(OK);
        case BSIM4_MOD_LUA1:
            value->rValue = model->BSIM4lua1; 
            return(OK);
        case BSIM4_MOD_LUB:
            value->rValue = model->BSIM4lub;  
            return(OK);
        case BSIM4_MOD_LUB1:
            value->rValue = model->BSIM4lub1;  
            return(OK);
        case BSIM4_MOD_LUC:
            value->rValue = model->BSIM4luc; 
            return(OK);
        case BSIM4_MOD_LUC1:
            value->rValue = model->BSIM4luc1; 
            return(OK);
        case BSIM4_MOD_LUD:
            value->rValue = model->BSIM4lud; 
            return(OK);
        case BSIM4_MOD_LUD1:
            value->rValue = model->BSIM4lud1; 
            return(OK);
        case BSIM4_MOD_LUP:
            value->rValue = model->BSIM4lup; 
            return(OK);
        case BSIM4_MOD_LLP:
            value->rValue = model->BSIM4llp; 
            return(OK);
        case BSIM4_MOD_LU0:
            value->rValue = model->BSIM4lu0;
            return(OK);
        case BSIM4_MOD_LUTE:
            value->rValue = model->BSIM4lute;
            return(OK);
		case BSIM4_MOD_LUCSTE:
            value->rValue = model->BSIM4lucste;
            return(OK);
        case BSIM4_MOD_LVOFF:
            value->rValue = model->BSIM4lvoff;
            return(OK);
        case BSIM4_MOD_LTVOFF:
            value->rValue = model->BSIM4ltvoff;
            return(OK);
        case BSIM4_MOD_LTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4ltnfactor;
            return(OK);
        case BSIM4_MOD_LTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4lteta0;
            return(OK);
        case BSIM4_MOD_LTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4ltvoffcv;
            return(OK);
        case BSIM4_MOD_LMINV:
            value->rValue = model->BSIM4lminv;
            return(OK);
        case BSIM4_MOD_LMINVCV:
            value->rValue = model->BSIM4lminvcv;
            return(OK);
        case BSIM4_MOD_LFPROUT:
            value->rValue = model->BSIM4lfprout;
            return(OK);
        case BSIM4_MOD_LPDITS:
            value->rValue = model->BSIM4lpdits;
            return(OK);
        case BSIM4_MOD_LPDITSD:
            value->rValue = model->BSIM4lpditsd;
            return(OK);
        case BSIM4_MOD_LDELTA:
            value->rValue = model->BSIM4ldelta;
            return(OK);
        case BSIM4_MOD_LRDSW:
            value->rValue = model->BSIM4lrdsw; 
            return(OK);             
        case BSIM4_MOD_LRDW:
            value->rValue = model->BSIM4lrdw;
            return(OK);
        case BSIM4_MOD_LRSW:
            value->rValue = model->BSIM4lrsw;
            return(OK);
        case BSIM4_MOD_LPRWB:
            value->rValue = model->BSIM4lprwb; 
            return(OK);             
        case BSIM4_MOD_LPRWG:
            value->rValue = model->BSIM4lprwg; 
            return(OK);             
        case BSIM4_MOD_LPRT:
            value->rValue = model->BSIM4lprt; 
            return(OK);              
        case BSIM4_MOD_LETA0:
            value->rValue = model->BSIM4leta0; 
            return(OK);               
        case BSIM4_MOD_LETAB:
            value->rValue = model->BSIM4letab; 
            return(OK);               
        case BSIM4_MOD_LPCLM:
            value->rValue = model->BSIM4lpclm; 
            return(OK);               
        case BSIM4_MOD_LPDIBL1:
            value->rValue = model->BSIM4lpdibl1; 
            return(OK);               
        case BSIM4_MOD_LPDIBL2:
            value->rValue = model->BSIM4lpdibl2; 
            return(OK);               
        case BSIM4_MOD_LPDIBLB:
            value->rValue = model->BSIM4lpdiblb; 
            return(OK);               
        case BSIM4_MOD_LPSCBE1:
            value->rValue = model->BSIM4lpscbe1; 
            return(OK);               
        case BSIM4_MOD_LPSCBE2:
            value->rValue = model->BSIM4lpscbe2; 
            return(OK);               
        case BSIM4_MOD_LPVAG:
            value->rValue = model->BSIM4lpvag; 
            return(OK);               
        case BSIM4_MOD_LWR:
            value->rValue = model->BSIM4lwr;
            return(OK);
        case BSIM4_MOD_LDWG:
            value->rValue = model->BSIM4ldwg;
            return(OK);
        case BSIM4_MOD_LDWB:
            value->rValue = model->BSIM4ldwb;
            return(OK);
        case BSIM4_MOD_LB0:
            value->rValue = model->BSIM4lb0;
            return(OK);
        case BSIM4_MOD_LB1:
            value->rValue = model->BSIM4lb1;
            return(OK);
        case BSIM4_MOD_LALPHA0:
            value->rValue = model->BSIM4lalpha0;
            return(OK);
        case BSIM4_MOD_LALPHA1:
            value->rValue = model->BSIM4lalpha1;
            return(OK);
        case BSIM4_MOD_LBETA0:
            value->rValue = model->BSIM4lbeta0;
            return(OK);
        case BSIM4_MOD_LAGIDL:
            value->rValue = model->BSIM4lagidl;
            return(OK);
        case BSIM4_MOD_LBGIDL:
            value->rValue = model->BSIM4lbgidl;
            return(OK);
        case BSIM4_MOD_LCGIDL:
            value->rValue = model->BSIM4lcgidl;
            return(OK);
	case BSIM4_MOD_LEGIDL:
            value->rValue = model->BSIM4legidl;
            return(OK);
        case BSIM4_MOD_LFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lfgidl;
            return(OK);
 	case BSIM4_MOD_LKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lkgidl;
            return(OK);
 	case BSIM4_MOD_LRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lrgidl;
            return(OK);
        case BSIM4_MOD_LAGISL:
            value->rValue = model->BSIM4lagisl;
            return(OK);
        case BSIM4_MOD_LBGISL:
            value->rValue = model->BSIM4lbgisl;
            return(OK);
        case BSIM4_MOD_LCGISL:
            value->rValue = model->BSIM4lcgisl;
            return(OK);
        case BSIM4_MOD_LEGISL:
            value->rValue = model->BSIM4legisl;
            return(OK);
        case BSIM4_MOD_LFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lfgisl;
            return(OK);
 	case BSIM4_MOD_LKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lkgisl;
            return(OK);
 	case BSIM4_MOD_LRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4lrgisl;
            return(OK);
        case BSIM4_MOD_LAIGC:
            value->rValue = model->BSIM4laigc;
            return(OK);
        case BSIM4_MOD_LBIGC:
            value->rValue = model->BSIM4lbigc;
            return(OK);
        case BSIM4_MOD_LCIGC:
            value->rValue = model->BSIM4lcigc;
            return(OK);
        case BSIM4_MOD_LAIGSD:
            value->rValue = model->BSIM4laigsd;
            return(OK);
        case BSIM4_MOD_LBIGSD:
            value->rValue = model->BSIM4lbigsd;
            return(OK);
        case BSIM4_MOD_LCIGSD:
            value->rValue = model->BSIM4lcigsd;
            return(OK);
        case BSIM4_MOD_LAIGS:
            value->rValue = model->BSIM4laigs;
            return(OK);
        case BSIM4_MOD_LBIGS:
            value->rValue = model->BSIM4lbigs;
            return(OK);
        case BSIM4_MOD_LCIGS:
            value->rValue = model->BSIM4lcigs;
            return(OK);
        case BSIM4_MOD_LAIGD:
            value->rValue = model->BSIM4laigd;
            return(OK);
        case BSIM4_MOD_LBIGD:
            value->rValue = model->BSIM4lbigd;
            return(OK);
        case BSIM4_MOD_LCIGD:
            value->rValue = model->BSIM4lcigd;
            return(OK);
        case BSIM4_MOD_LAIGBACC:
            value->rValue = model->BSIM4laigbacc;
            return(OK);
        case BSIM4_MOD_LBIGBACC:
            value->rValue = model->BSIM4lbigbacc;
            return(OK);
        case BSIM4_MOD_LCIGBACC:
            value->rValue = model->BSIM4lcigbacc;
            return(OK);
        case BSIM4_MOD_LAIGBINV:
            value->rValue = model->BSIM4laigbinv;
            return(OK);
        case BSIM4_MOD_LBIGBINV:
            value->rValue = model->BSIM4lbigbinv;
            return(OK);
        case BSIM4_MOD_LCIGBINV:
            value->rValue = model->BSIM4lcigbinv;
            return(OK);
        case BSIM4_MOD_LNIGC:
            value->rValue = model->BSIM4lnigc;
            return(OK);
        case BSIM4_MOD_LNIGBACC:
            value->rValue = model->BSIM4lnigbacc;
            return(OK);
        case BSIM4_MOD_LNIGBINV:
            value->rValue = model->BSIM4lnigbinv;
            return(OK);
        case BSIM4_MOD_LNTOX:
            value->rValue = model->BSIM4lntox;
            return(OK);
        case BSIM4_MOD_LEIGBINV:
            value->rValue = model->BSIM4leigbinv;
            return(OK);
        case BSIM4_MOD_LPIGCD:
            value->rValue = model->BSIM4lpigcd;
            return(OK);
        case BSIM4_MOD_LPOXEDGE:
            value->rValue = model->BSIM4lpoxedge;
            return(OK);
        case BSIM4_MOD_LPHIN:
            value->rValue = model->BSIM4lphin;
            return(OK);
        case BSIM4_MOD_LXRCRG1:
            value->rValue = model->BSIM4lxrcrg1;
            return(OK);
        case BSIM4_MOD_LXRCRG2:
            value->rValue = model->BSIM4lxrcrg2;
            return(OK);
        case BSIM4_MOD_LEU:
            value->rValue = model->BSIM4leu;
            return(OK);
		case BSIM4_MOD_LUCS:
            value->rValue = model->BSIM4lucs;
            return(OK);
        case BSIM4_MOD_LVFB:
            value->rValue = model->BSIM4lvfb;
            return(OK);

        case BSIM4_MOD_LCGSL:
            value->rValue = model->BSIM4lcgsl;
            return(OK);
        case BSIM4_MOD_LCGDL:
            value->rValue = model->BSIM4lcgdl;
            return(OK);
        case BSIM4_MOD_LCKAPPAS:
            value->rValue = model->BSIM4lckappas;
            return(OK);
        case BSIM4_MOD_LCKAPPAD:
            value->rValue = model->BSIM4lckappad;
            return(OK);
        case BSIM4_MOD_LCF:
            value->rValue = model->BSIM4lcf;
            return(OK);
        case BSIM4_MOD_LCLC:
            value->rValue = model->BSIM4lclc;
            return(OK);
        case BSIM4_MOD_LCLE:
            value->rValue = model->BSIM4lcle;
            return(OK);
        case BSIM4_MOD_LVFBCV:
            value->rValue = model->BSIM4lvfbcv;
            return(OK);
        case BSIM4_MOD_LACDE:
            value->rValue = model->BSIM4lacde;
            return(OK);
        case BSIM4_MOD_LMOIN:
            value->rValue = model->BSIM4lmoin;
            return(OK);
        case BSIM4_MOD_LNOFF:
            value->rValue = model->BSIM4lnoff;
            return(OK);
        case BSIM4_MOD_LVOFFCV:
            value->rValue = model->BSIM4lvoffcv;
            return(OK);
        case BSIM4_MOD_LVFBSDOFF:
            value->rValue = model->BSIM4lvfbsdoff;
            return(OK);
        case BSIM4_MOD_LTVFBSDOFF:
            value->rValue = model->BSIM4ltvfbsdoff;
            return(OK);

        case BSIM4_MOD_LLAMBDA:
            value->rValue = model->BSIM4llambda;
            return(OK);
        case BSIM4_MOD_LVTL:
            value->rValue = model->BSIM4lvtl;
            return(OK);
        case BSIM4_MOD_LXN:
            value->rValue = model->BSIM4lxn;
            return(OK);

	/* Width dependence */
        case  BSIM4_MOD_WCDSC :
          value->rValue = model->BSIM4wcdsc;
            return(OK);
        case  BSIM4_MOD_WCDSCB :
          value->rValue = model->BSIM4wcdscb;
            return(OK);
        case  BSIM4_MOD_WCDSCD :
          value->rValue = model->BSIM4wcdscd;
            return(OK);
        case  BSIM4_MOD_WCIT :
          value->rValue = model->BSIM4wcit;
            return(OK);
        case  BSIM4_MOD_WNFACTOR :
          value->rValue = model->BSIM4wnfactor;
            return(OK);
        case BSIM4_MOD_WXJ:
            value->rValue = model->BSIM4wxj;
            return(OK);
        case BSIM4_MOD_WVSAT:
            value->rValue = model->BSIM4wvsat;
            return(OK);
        case BSIM4_MOD_WAT:
            value->rValue = model->BSIM4wat;
            return(OK);
        case BSIM4_MOD_WA0:
            value->rValue = model->BSIM4wa0;
            return(OK);
        case BSIM4_MOD_WAGS:
            value->rValue = model->BSIM4wags;
            return(OK);
        case BSIM4_MOD_WA1:
            value->rValue = model->BSIM4wa1;
            return(OK);
        case BSIM4_MOD_WA2:
            value->rValue = model->BSIM4wa2;
            return(OK);
        case BSIM4_MOD_WKETA:
            value->rValue = model->BSIM4wketa;
            return(OK);   
        case BSIM4_MOD_WNSUB:
            value->rValue = model->BSIM4wnsub;
            return(OK);
        case BSIM4_MOD_WNDEP:
            value->rValue = model->BSIM4wndep;
            return(OK);
        case BSIM4_MOD_WNSD:
            value->rValue = model->BSIM4wnsd;
            return(OK);
        case BSIM4_MOD_WNGATE:
            value->rValue = model->BSIM4wngate;
            return(OK);
        case BSIM4_MOD_WGAMMA1:
            value->rValue = model->BSIM4wgamma1;
            return(OK);
        case BSIM4_MOD_WGAMMA2:
            value->rValue = model->BSIM4wgamma2;
            return(OK);
        case BSIM4_MOD_WVBX:
            value->rValue = model->BSIM4wvbx;
            return(OK);
        case BSIM4_MOD_WVBM:
            value->rValue = model->BSIM4wvbm;
            return(OK);
        case BSIM4_MOD_WXT:
            value->rValue = model->BSIM4wxt;
            return(OK);
        case  BSIM4_MOD_WK1:
          value->rValue = model->BSIM4wk1;
            return(OK);
        case  BSIM4_MOD_WKT1:
          value->rValue = model->BSIM4wkt1;
            return(OK);
        case  BSIM4_MOD_WKT1L:
          value->rValue = model->BSIM4wkt1l;
            return(OK);
        case  BSIM4_MOD_WKT2 :
          value->rValue = model->BSIM4wkt2;
            return(OK);
        case  BSIM4_MOD_WK2 :
          value->rValue = model->BSIM4wk2;
            return(OK);
        case  BSIM4_MOD_WK3:
          value->rValue = model->BSIM4wk3;
            return(OK);
        case  BSIM4_MOD_WK3B:
          value->rValue = model->BSIM4wk3b;
            return(OK);
        case  BSIM4_MOD_WW0:
          value->rValue = model->BSIM4ww0;
            return(OK);
        case  BSIM4_MOD_WLPE0:
          value->rValue = model->BSIM4wlpe0;
            return(OK);
        case  BSIM4_MOD_WDVTP0:
          value->rValue = model->BSIM4wdvtp0;
            return(OK);
        case  BSIM4_MOD_WDVTP1:
          value->rValue = model->BSIM4wdvtp1;
            return(OK);
        case  BSIM4_MOD_WDVTP2:
          value->rValue = model->BSIM4wdvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4_MOD_WDVTP3:
          value->rValue = model->BSIM4wdvtp3;
            return(OK);
        case  BSIM4_MOD_WDVTP4:
          value->rValue = model->BSIM4wdvtp4;
            return(OK);
        case  BSIM4_MOD_WDVTP5:
          value->rValue = model->BSIM4wdvtp5;
            return(OK);
        case  BSIM4_MOD_WLPEB:
          value->rValue = model->BSIM4wlpeb;
            return(OK);
        case  BSIM4_MOD_WDVT0:                
          value->rValue = model->BSIM4wdvt0;
            return(OK);
        case  BSIM4_MOD_WDVT1 :             
          value->rValue = model->BSIM4wdvt1;
            return(OK);
        case  BSIM4_MOD_WDVT2 :             
          value->rValue = model->BSIM4wdvt2;
            return(OK);
        case  BSIM4_MOD_WDVT0W :                
          value->rValue = model->BSIM4wdvt0w;
            return(OK);
        case  BSIM4_MOD_WDVT1W :             
          value->rValue = model->BSIM4wdvt1w;
            return(OK);
        case  BSIM4_MOD_WDVT2W :             
          value->rValue = model->BSIM4wdvt2w;
            return(OK);
        case  BSIM4_MOD_WDROUT :           
          value->rValue = model->BSIM4wdrout;
            return(OK);
        case  BSIM4_MOD_WDSUB :           
          value->rValue = model->BSIM4wdsub;
            return(OK);
        case BSIM4_MOD_WVTH0:
            value->rValue = model->BSIM4wvth0; 
            return(OK);
        case BSIM4_MOD_WUA:
            value->rValue = model->BSIM4wua; 
            return(OK);
        case BSIM4_MOD_WUA1:
            value->rValue = model->BSIM4wua1; 
            return(OK);
        case BSIM4_MOD_WUB:
            value->rValue = model->BSIM4wub;  
            return(OK);
        case BSIM4_MOD_WUB1:
            value->rValue = model->BSIM4wub1;  
            return(OK);
        case BSIM4_MOD_WUC:
            value->rValue = model->BSIM4wuc; 
            return(OK);
        case BSIM4_MOD_WUC1:
            value->rValue = model->BSIM4wuc1; 
            return(OK);
        case BSIM4_MOD_WUD:
            value->rValue = model->BSIM4wud; 
            return(OK);
        case BSIM4_MOD_WUD1:
            value->rValue = model->BSIM4wud1; 
            return(OK);
        case BSIM4_MOD_WUP:
            value->rValue = model->BSIM4wup; 
            return(OK);
        case BSIM4_MOD_WLP:
            value->rValue = model->BSIM4wlp; 
            return(OK);
        case BSIM4_MOD_WU0:
            value->rValue = model->BSIM4wu0;
            return(OK);
        case BSIM4_MOD_WUTE:
            value->rValue = model->BSIM4wute;
            return(OK);
        case BSIM4_MOD_WUCSTE:
            value->rValue = model->BSIM4wucste;
            return(OK);
        case BSIM4_MOD_WVOFF:
            value->rValue = model->BSIM4wvoff;
            return(OK);
        case BSIM4_MOD_WTVOFF:
            value->rValue = model->BSIM4wtvoff;
            return(OK);
        case BSIM4_MOD_WTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4wtnfactor;
            return(OK);
        case BSIM4_MOD_WTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4wteta0;
            return(OK);
        case BSIM4_MOD_WTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4wtvoffcv;
            return(OK);
        case BSIM4_MOD_WMINV:
            value->rValue = model->BSIM4wminv;
            return(OK);
        case BSIM4_MOD_WMINVCV:
            value->rValue = model->BSIM4wminvcv;
            return(OK);
        case BSIM4_MOD_WFPROUT:
            value->rValue = model->BSIM4wfprout;
            return(OK);
        case BSIM4_MOD_WPDITS:
            value->rValue = model->BSIM4wpdits;
            return(OK);
        case BSIM4_MOD_WPDITSD:
            value->rValue = model->BSIM4wpditsd;
            return(OK);
        case BSIM4_MOD_WDELTA:
            value->rValue = model->BSIM4wdelta;
            return(OK);
        case BSIM4_MOD_WRDSW:
            value->rValue = model->BSIM4wrdsw; 
            return(OK);             
        case BSIM4_MOD_WRDW:
            value->rValue = model->BSIM4wrdw;
            return(OK);
        case BSIM4_MOD_WRSW:
            value->rValue = model->BSIM4wrsw;
            return(OK);
        case BSIM4_MOD_WPRWB:
            value->rValue = model->BSIM4wprwb; 
            return(OK);             
        case BSIM4_MOD_WPRWG:
            value->rValue = model->BSIM4wprwg; 
            return(OK);             
        case BSIM4_MOD_WPRT:
            value->rValue = model->BSIM4wprt; 
            return(OK);              
        case BSIM4_MOD_WETA0:
            value->rValue = model->BSIM4weta0; 
            return(OK);               
        case BSIM4_MOD_WETAB:
            value->rValue = model->BSIM4wetab; 
            return(OK);               
        case BSIM4_MOD_WPCLM:
            value->rValue = model->BSIM4wpclm; 
            return(OK);               
        case BSIM4_MOD_WPDIBL1:
            value->rValue = model->BSIM4wpdibl1; 
            return(OK);               
        case BSIM4_MOD_WPDIBL2:
            value->rValue = model->BSIM4wpdibl2; 
            return(OK);               
        case BSIM4_MOD_WPDIBLB:
            value->rValue = model->BSIM4wpdiblb; 
            return(OK);               
        case BSIM4_MOD_WPSCBE1:
            value->rValue = model->BSIM4wpscbe1; 
            return(OK);               
        case BSIM4_MOD_WPSCBE2:
            value->rValue = model->BSIM4wpscbe2; 
            return(OK);               
        case BSIM4_MOD_WPVAG:
            value->rValue = model->BSIM4wpvag; 
            return(OK);               
        case BSIM4_MOD_WWR:
            value->rValue = model->BSIM4wwr;
            return(OK);
        case BSIM4_MOD_WDWG:
            value->rValue = model->BSIM4wdwg;
            return(OK);
        case BSIM4_MOD_WDWB:
            value->rValue = model->BSIM4wdwb;
            return(OK);
        case BSIM4_MOD_WB0:
            value->rValue = model->BSIM4wb0;
            return(OK);
        case BSIM4_MOD_WB1:
            value->rValue = model->BSIM4wb1;
            return(OK);
        case BSIM4_MOD_WALPHA0:
            value->rValue = model->BSIM4walpha0;
            return(OK);
        case BSIM4_MOD_WALPHA1:
            value->rValue = model->BSIM4walpha1;
            return(OK);
        case BSIM4_MOD_WBETA0:
            value->rValue = model->BSIM4wbeta0;
            return(OK);
        case BSIM4_MOD_WAGIDL:
            value->rValue = model->BSIM4wagidl;
            return(OK);
        case BSIM4_MOD_WBGIDL:
            value->rValue = model->BSIM4wbgidl;
            return(OK);
        case BSIM4_MOD_WCGIDL:
            value->rValue = model->BSIM4wcgidl;
            return(OK);
        case BSIM4_MOD_WEGIDL:
            value->rValue = model->BSIM4wegidl;
            return(OK);
        case BSIM4_MOD_WFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wfgidl;
            return(OK);
        case BSIM4_MOD_WKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wkgidl;
            return(OK);
        case BSIM4_MOD_WRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wrgidl;
            return(OK);
        case BSIM4_MOD_WAGISL:
            value->rValue = model->BSIM4wagisl;
            return(OK);
        case BSIM4_MOD_WBGISL:
            value->rValue = model->BSIM4wbgisl;
            return(OK);
        case BSIM4_MOD_WCGISL:
            value->rValue = model->BSIM4wcgisl;
            return(OK);
        case BSIM4_MOD_WEGISL:
            value->rValue = model->BSIM4wegisl;
            return(OK);
        case BSIM4_MOD_WFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wfgisl;
            return(OK);
        case BSIM4_MOD_WKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wkgisl;
            return(OK);
        case BSIM4_MOD_WRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4wrgisl;
            return(OK);
        case BSIM4_MOD_WAIGC:
            value->rValue = model->BSIM4waigc;
            return(OK);
        case BSIM4_MOD_WBIGC:
            value->rValue = model->BSIM4wbigc;
            return(OK);
        case BSIM4_MOD_WCIGC:
            value->rValue = model->BSIM4wcigc;
            return(OK);
        case BSIM4_MOD_WAIGSD:
            value->rValue = model->BSIM4waigsd;
            return(OK);
        case BSIM4_MOD_WBIGSD:
            value->rValue = model->BSIM4wbigsd;
            return(OK);
        case BSIM4_MOD_WCIGSD:
            value->rValue = model->BSIM4wcigsd;
            return(OK);
        case BSIM4_MOD_WAIGS:
            value->rValue = model->BSIM4waigs;
            return(OK);
        case BSIM4_MOD_WBIGS:
            value->rValue = model->BSIM4wbigs;
            return(OK);
        case BSIM4_MOD_WCIGS:
            value->rValue = model->BSIM4wcigs;
            return(OK);
        case BSIM4_MOD_WAIGD:
            value->rValue = model->BSIM4waigd;
            return(OK);
        case BSIM4_MOD_WBIGD:
            value->rValue = model->BSIM4wbigd;
            return(OK);
        case BSIM4_MOD_WCIGD:
            value->rValue = model->BSIM4wcigd;
            return(OK);
        case BSIM4_MOD_WAIGBACC:
            value->rValue = model->BSIM4waigbacc;
            return(OK);
        case BSIM4_MOD_WBIGBACC:
            value->rValue = model->BSIM4wbigbacc;
            return(OK);
        case BSIM4_MOD_WCIGBACC:
            value->rValue = model->BSIM4wcigbacc;
            return(OK);
        case BSIM4_MOD_WAIGBINV:
            value->rValue = model->BSIM4waigbinv;
            return(OK);
        case BSIM4_MOD_WBIGBINV:
            value->rValue = model->BSIM4wbigbinv;
            return(OK);
        case BSIM4_MOD_WCIGBINV:
            value->rValue = model->BSIM4wcigbinv;
            return(OK);
        case BSIM4_MOD_WNIGC:
            value->rValue = model->BSIM4wnigc;
            return(OK);
        case BSIM4_MOD_WNIGBACC:
            value->rValue = model->BSIM4wnigbacc;
            return(OK);
        case BSIM4_MOD_WNIGBINV:
            value->rValue = model->BSIM4wnigbinv;
            return(OK);
        case BSIM4_MOD_WNTOX:
            value->rValue = model->BSIM4wntox;
            return(OK);
        case BSIM4_MOD_WEIGBINV:
            value->rValue = model->BSIM4weigbinv;
            return(OK);
        case BSIM4_MOD_WPIGCD:
            value->rValue = model->BSIM4wpigcd;
            return(OK);
        case BSIM4_MOD_WPOXEDGE:
            value->rValue = model->BSIM4wpoxedge;
            return(OK);
        case BSIM4_MOD_WPHIN:
            value->rValue = model->BSIM4wphin;
            return(OK);
        case BSIM4_MOD_WXRCRG1:
            value->rValue = model->BSIM4wxrcrg1;
            return(OK);
        case BSIM4_MOD_WXRCRG2:
            value->rValue = model->BSIM4wxrcrg2;
            return(OK);
        case BSIM4_MOD_WEU:
            value->rValue = model->BSIM4weu;
            return(OK);
        case BSIM4_MOD_WUCS:
            value->rValue = model->BSIM4wucs;
            return(OK);
        case BSIM4_MOD_WVFB:
            value->rValue = model->BSIM4wvfb;
            return(OK);

        case BSIM4_MOD_WCGSL:
            value->rValue = model->BSIM4wcgsl;
            return(OK);
        case BSIM4_MOD_WCGDL:
            value->rValue = model->BSIM4wcgdl;
            return(OK);
        case BSIM4_MOD_WCKAPPAS:
            value->rValue = model->BSIM4wckappas;
            return(OK);
        case BSIM4_MOD_WCKAPPAD:
            value->rValue = model->BSIM4wckappad;
            return(OK);
        case BSIM4_MOD_WCF:
            value->rValue = model->BSIM4wcf;
            return(OK);
        case BSIM4_MOD_WCLC:
            value->rValue = model->BSIM4wclc;
            return(OK);
        case BSIM4_MOD_WCLE:
            value->rValue = model->BSIM4wcle;
            return(OK);
        case BSIM4_MOD_WVFBCV:
            value->rValue = model->BSIM4wvfbcv;
            return(OK);
        case BSIM4_MOD_WACDE:
            value->rValue = model->BSIM4wacde;
            return(OK);
        case BSIM4_MOD_WMOIN:
            value->rValue = model->BSIM4wmoin;
            return(OK);
        case BSIM4_MOD_WNOFF:
            value->rValue = model->BSIM4wnoff;
            return(OK);
        case BSIM4_MOD_WVOFFCV:
            value->rValue = model->BSIM4wvoffcv;
            return(OK);
        case BSIM4_MOD_WVFBSDOFF:
            value->rValue = model->BSIM4wvfbsdoff;
            return(OK);
        case BSIM4_MOD_WTVFBSDOFF:
            value->rValue = model->BSIM4wtvfbsdoff;
            return(OK);

        case BSIM4_MOD_WLAMBDA:
            value->rValue = model->BSIM4wlambda;
            return(OK);
        case BSIM4_MOD_WVTL:
            value->rValue = model->BSIM4wvtl;
            return(OK);
        case BSIM4_MOD_WXN:
            value->rValue = model->BSIM4wxn;
            return(OK);

	/* Cross-term dependence */
        case  BSIM4_MOD_PCDSC :
          value->rValue = model->BSIM4pcdsc;
            return(OK);
        case  BSIM4_MOD_PCDSCB :
          value->rValue = model->BSIM4pcdscb;
            return(OK);
        case  BSIM4_MOD_PCDSCD :
          value->rValue = model->BSIM4pcdscd;
            return(OK);
         case  BSIM4_MOD_PCIT :
          value->rValue = model->BSIM4pcit;
            return(OK);
        case  BSIM4_MOD_PNFACTOR :
          value->rValue = model->BSIM4pnfactor;
            return(OK);
        case BSIM4_MOD_PXJ:
            value->rValue = model->BSIM4pxj;
            return(OK);
        case BSIM4_MOD_PVSAT:
            value->rValue = model->BSIM4pvsat;
            return(OK);
        case BSIM4_MOD_PAT:
            value->rValue = model->BSIM4pat;
            return(OK);
        case BSIM4_MOD_PA0:
            value->rValue = model->BSIM4pa0;
            return(OK);
        case BSIM4_MOD_PAGS:
            value->rValue = model->BSIM4pags;
            return(OK);
        case BSIM4_MOD_PA1:
            value->rValue = model->BSIM4pa1;
            return(OK);
        case BSIM4_MOD_PA2:
            value->rValue = model->BSIM4pa2;
            return(OK);
        case BSIM4_MOD_PKETA:
            value->rValue = model->BSIM4pketa;
            return(OK);   
        case BSIM4_MOD_PNSUB:
            value->rValue = model->BSIM4pnsub;
            return(OK);
        case BSIM4_MOD_PNDEP:
            value->rValue = model->BSIM4pndep;
            return(OK);
        case BSIM4_MOD_PNSD:
            value->rValue = model->BSIM4pnsd;
            return(OK);
        case BSIM4_MOD_PNGATE:
            value->rValue = model->BSIM4pngate;
            return(OK);
        case BSIM4_MOD_PGAMMA1:
            value->rValue = model->BSIM4pgamma1;
            return(OK);
        case BSIM4_MOD_PGAMMA2:
            value->rValue = model->BSIM4pgamma2;
            return(OK);
        case BSIM4_MOD_PVBX:
            value->rValue = model->BSIM4pvbx;
            return(OK);
        case BSIM4_MOD_PVBM:
            value->rValue = model->BSIM4pvbm;
            return(OK);
        case BSIM4_MOD_PXT:
            value->rValue = model->BSIM4pxt;
            return(OK);
        case  BSIM4_MOD_PK1:
          value->rValue = model->BSIM4pk1;
            return(OK);
        case  BSIM4_MOD_PKT1:
          value->rValue = model->BSIM4pkt1;
            return(OK);
        case  BSIM4_MOD_PKT1L:
          value->rValue = model->BSIM4pkt1l;
            return(OK);
        case  BSIM4_MOD_PKT2 :
          value->rValue = model->BSIM4pkt2;
            return(OK);
        case  BSIM4_MOD_PK2 :
          value->rValue = model->BSIM4pk2;
            return(OK);
        case  BSIM4_MOD_PK3:
          value->rValue = model->BSIM4pk3;
            return(OK);
        case  BSIM4_MOD_PK3B:
          value->rValue = model->BSIM4pk3b;
            return(OK);
        case  BSIM4_MOD_PW0:
          value->rValue = model->BSIM4pw0;
            return(OK);
        case  BSIM4_MOD_PLPE0:
          value->rValue = model->BSIM4plpe0;
            return(OK);
        case  BSIM4_MOD_PLPEB:
          value->rValue = model->BSIM4plpeb;
            return(OK);
        case  BSIM4_MOD_PDVTP0:
          value->rValue = model->BSIM4pdvtp0;
            return(OK);
        case  BSIM4_MOD_PDVTP1:
          value->rValue = model->BSIM4pdvtp1;
            return(OK);
        case  BSIM4_MOD_PDVTP2:
          value->rValue = model->BSIM4pdvtp2;  /* New DIBL/Rout */
            return(OK);
        case  BSIM4_MOD_PDVTP3:
          value->rValue = model->BSIM4pdvtp3;
            return(OK);
        case  BSIM4_MOD_PDVTP4:
          value->rValue = model->BSIM4pdvtp4;
            return(OK);
        case  BSIM4_MOD_PDVTP5:
          value->rValue = model->BSIM4pdvtp5;
            return(OK);
        case  BSIM4_MOD_PDVT0 :                
          value->rValue = model->BSIM4pdvt0;
            return(OK);
        case  BSIM4_MOD_PDVT1 :             
          value->rValue = model->BSIM4pdvt1;
            return(OK);
        case  BSIM4_MOD_PDVT2 :             
          value->rValue = model->BSIM4pdvt2;
            return(OK);
        case  BSIM4_MOD_PDVT0W :                
          value->rValue = model->BSIM4pdvt0w;
            return(OK);
        case  BSIM4_MOD_PDVT1W :             
          value->rValue = model->BSIM4pdvt1w;
            return(OK);
        case  BSIM4_MOD_PDVT2W :             
          value->rValue = model->BSIM4pdvt2w;
            return(OK);
        case  BSIM4_MOD_PDROUT :           
          value->rValue = model->BSIM4pdrout;
            return(OK);
        case  BSIM4_MOD_PDSUB :           
          value->rValue = model->BSIM4pdsub;
            return(OK);
        case BSIM4_MOD_PVTH0:
            value->rValue = model->BSIM4pvth0; 
            return(OK);
        case BSIM4_MOD_PUA:
            value->rValue = model->BSIM4pua; 
            return(OK);
        case BSIM4_MOD_PUA1:
            value->rValue = model->BSIM4pua1; 
            return(OK);
        case BSIM4_MOD_PUB:
            value->rValue = model->BSIM4pub;  
            return(OK);
        case BSIM4_MOD_PUB1:
            value->rValue = model->BSIM4pub1;  
            return(OK);
        case BSIM4_MOD_PUC:
            value->rValue = model->BSIM4puc; 
            return(OK);
        case BSIM4_MOD_PUC1:
            value->rValue = model->BSIM4puc1; 
            return(OK);
        case BSIM4_MOD_PUD:
            value->rValue = model->BSIM4pud; 
            return(OK);
        case BSIM4_MOD_PUD1:
            value->rValue = model->BSIM4pud1; 
            return(OK);
        case BSIM4_MOD_PUP:
            value->rValue = model->BSIM4pup; 
            return(OK);
        case BSIM4_MOD_PLP:
            value->rValue = model->BSIM4plp; 
            return(OK);
        case BSIM4_MOD_PU0:
            value->rValue = model->BSIM4pu0;
            return(OK);
        case BSIM4_MOD_PUTE:
            value->rValue = model->BSIM4pute;
            return(OK);
        case BSIM4_MOD_PUCSTE:
            value->rValue = model->BSIM4pucste;
            return(OK);
        case BSIM4_MOD_PVOFF:
            value->rValue = model->BSIM4pvoff;
            return(OK);
        case BSIM4_MOD_PTVOFF:
            value->rValue = model->BSIM4ptvoff;
            return(OK);
        case BSIM4_MOD_PTNFACTOR:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4ptnfactor;
            return(OK);
        case BSIM4_MOD_PTETA0:		/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4pteta0;
            return(OK);
        case BSIM4_MOD_PTVOFFCV:	/* v4.7 temp dep of leakage current  */
            value->rValue = model->BSIM4ptvoffcv;
            return(OK);
        case BSIM4_MOD_PMINV:
            value->rValue = model->BSIM4pminv;
            return(OK);
        case BSIM4_MOD_PMINVCV:
            value->rValue = model->BSIM4pminvcv;
            return(OK);
        case BSIM4_MOD_PFPROUT:
            value->rValue = model->BSIM4pfprout;
            return(OK);
        case BSIM4_MOD_PPDITS:
            value->rValue = model->BSIM4ppdits;
            return(OK);
        case BSIM4_MOD_PPDITSD:
            value->rValue = model->BSIM4ppditsd;
            return(OK);
        case BSIM4_MOD_PDELTA:
            value->rValue = model->BSIM4pdelta;
            return(OK);
        case BSIM4_MOD_PRDSW:
            value->rValue = model->BSIM4prdsw; 
            return(OK);             
        case BSIM4_MOD_PRDW:
            value->rValue = model->BSIM4prdw;
            return(OK);
        case BSIM4_MOD_PRSW:
            value->rValue = model->BSIM4prsw;
            return(OK);
        case BSIM4_MOD_PPRWB:
            value->rValue = model->BSIM4pprwb; 
            return(OK);             
        case BSIM4_MOD_PPRWG:
            value->rValue = model->BSIM4pprwg; 
            return(OK);             
        case BSIM4_MOD_PPRT:
            value->rValue = model->BSIM4pprt; 
            return(OK);              
        case BSIM4_MOD_PETA0:
            value->rValue = model->BSIM4peta0; 
            return(OK);               
        case BSIM4_MOD_PETAB:
            value->rValue = model->BSIM4petab; 
            return(OK);               
        case BSIM4_MOD_PPCLM:
            value->rValue = model->BSIM4ppclm; 
            return(OK);               
        case BSIM4_MOD_PPDIBL1:
            value->rValue = model->BSIM4ppdibl1; 
            return(OK);               
        case BSIM4_MOD_PPDIBL2:
            value->rValue = model->BSIM4ppdibl2; 
            return(OK);               
        case BSIM4_MOD_PPDIBLB:
            value->rValue = model->BSIM4ppdiblb; 
            return(OK);               
        case BSIM4_MOD_PPSCBE1:
            value->rValue = model->BSIM4ppscbe1; 
            return(OK);               
        case BSIM4_MOD_PPSCBE2:
            value->rValue = model->BSIM4ppscbe2; 
            return(OK);               
        case BSIM4_MOD_PPVAG:
            value->rValue = model->BSIM4ppvag; 
            return(OK);               
        case BSIM4_MOD_PWR:
            value->rValue = model->BSIM4pwr;
            return(OK);
        case BSIM4_MOD_PDWG:
            value->rValue = model->BSIM4pdwg;
            return(OK);
        case BSIM4_MOD_PDWB:
            value->rValue = model->BSIM4pdwb;
            return(OK);
        case BSIM4_MOD_PB0:
            value->rValue = model->BSIM4pb0;
            return(OK);
        case BSIM4_MOD_PB1:
            value->rValue = model->BSIM4pb1;
            return(OK);
        case BSIM4_MOD_PALPHA0:
            value->rValue = model->BSIM4palpha0;
            return(OK);
        case BSIM4_MOD_PALPHA1:
            value->rValue = model->BSIM4palpha1;
            return(OK);
        case BSIM4_MOD_PBETA0:
            value->rValue = model->BSIM4pbeta0;
            return(OK);
        case BSIM4_MOD_PAGIDL:
            value->rValue = model->BSIM4pagidl;
            return(OK);
        case BSIM4_MOD_PBGIDL:
            value->rValue = model->BSIM4pbgidl;
            return(OK);
        case BSIM4_MOD_PCGIDL:
            value->rValue = model->BSIM4pcgidl;
            return(OK);
        case BSIM4_MOD_PEGIDL:
            value->rValue = model->BSIM4pegidl;
            return(OK);
        case BSIM4_MOD_PFGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4pfgidl;
            return(OK);
        case BSIM4_MOD_PKGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4pkgidl;
            return(OK);
        case BSIM4_MOD_PRGIDL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4prgidl;
            return(OK);
        case BSIM4_MOD_PAGISL:
            value->rValue = model->BSIM4pagisl;
            return(OK);
        case BSIM4_MOD_PBGISL:
            value->rValue = model->BSIM4pbgisl;
            return(OK);
        case BSIM4_MOD_PCGISL:
            value->rValue = model->BSIM4pcgisl;
            return(OK);
        case BSIM4_MOD_PEGISL:
            value->rValue = model->BSIM4pegisl;
            return(OK);
        case BSIM4_MOD_PFGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4pfgisl;
            return(OK);
        case BSIM4_MOD_PKGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4pkgisl;
            return(OK);
        case BSIM4_MOD_PRGISL:			/* v4.7 New GIDL/GISL*/
            value->rValue = model->BSIM4prgisl;
            return(OK);
        case BSIM4_MOD_PAIGC:
            value->rValue = model->BSIM4paigc;
            return(OK);
        case BSIM4_MOD_PBIGC:
            value->rValue = model->BSIM4pbigc;
            return(OK);
        case BSIM4_MOD_PCIGC:
            value->rValue = model->BSIM4pcigc;
            return(OK);
        case BSIM4_MOD_PAIGSD:
            value->rValue = model->BSIM4paigsd;
            return(OK);
        case BSIM4_MOD_PBIGSD:
            value->rValue = model->BSIM4pbigsd;
            return(OK);
        case BSIM4_MOD_PCIGSD:
            value->rValue = model->BSIM4pcigsd;
            return(OK);
        case BSIM4_MOD_PAIGS:
            value->rValue = model->BSIM4paigs;
            return(OK);
        case BSIM4_MOD_PBIGS:
            value->rValue = model->BSIM4pbigs;
            return(OK);
        case BSIM4_MOD_PCIGS:
            value->rValue = model->BSIM4pcigs;
            return(OK);
        case BSIM4_MOD_PAIGD:
            value->rValue = model->BSIM4paigd;
            return(OK);
        case BSIM4_MOD_PBIGD:
            value->rValue = model->BSIM4pbigd;
            return(OK);
        case BSIM4_MOD_PCIGD:
            value->rValue = model->BSIM4pcigd;
            return(OK);
        case BSIM4_MOD_PAIGBACC:
            value->rValue = model->BSIM4paigbacc;
            return(OK);
        case BSIM4_MOD_PBIGBACC:
            value->rValue = model->BSIM4pbigbacc;
            return(OK);
        case BSIM4_MOD_PCIGBACC:
            value->rValue = model->BSIM4pcigbacc;
            return(OK);
        case BSIM4_MOD_PAIGBINV:
            value->rValue = model->BSIM4paigbinv;
            return(OK);
        case BSIM4_MOD_PBIGBINV:
            value->rValue = model->BSIM4pbigbinv;
            return(OK);
        case BSIM4_MOD_PCIGBINV:
            value->rValue = model->BSIM4pcigbinv;
            return(OK);
        case BSIM4_MOD_PNIGC:
            value->rValue = model->BSIM4pnigc;
            return(OK);
        case BSIM4_MOD_PNIGBACC:
            value->rValue = model->BSIM4pnigbacc;
            return(OK);
        case BSIM4_MOD_PNIGBINV:
            value->rValue = model->BSIM4pnigbinv;
            return(OK);
        case BSIM4_MOD_PNTOX:
            value->rValue = model->BSIM4pntox;
            return(OK);
        case BSIM4_MOD_PEIGBINV:
            value->rValue = model->BSIM4peigbinv;
            return(OK);
        case BSIM4_MOD_PPIGCD:
            value->rValue = model->BSIM4ppigcd;
            return(OK);
        case BSIM4_MOD_PPOXEDGE:
            value->rValue = model->BSIM4ppoxedge;
            return(OK);
        case BSIM4_MOD_PPHIN:
            value->rValue = model->BSIM4pphin;
            return(OK);
        case BSIM4_MOD_PXRCRG1:
            value->rValue = model->BSIM4pxrcrg1;
            return(OK);
        case BSIM4_MOD_PXRCRG2:
            value->rValue = model->BSIM4pxrcrg2;
            return(OK);
        case BSIM4_MOD_PEU:
            value->rValue = model->BSIM4peu;
            return(OK);
        case BSIM4_MOD_PUCS:
            value->rValue = model->BSIM4pucs;
            return(OK);
        case BSIM4_MOD_PVFB:
            value->rValue = model->BSIM4pvfb;
            return(OK);

        case BSIM4_MOD_PCGSL:
            value->rValue = model->BSIM4pcgsl;
            return(OK);
        case BSIM4_MOD_PCGDL:
            value->rValue = model->BSIM4pcgdl;
            return(OK);
        case BSIM4_MOD_PCKAPPAS:
            value->rValue = model->BSIM4pckappas;
            return(OK);
        case BSIM4_MOD_PCKAPPAD:
            value->rValue = model->BSIM4pckappad;
            return(OK);
        case BSIM4_MOD_PCF:
            value->rValue = model->BSIM4pcf;
            return(OK);
        case BSIM4_MOD_PCLC:
            value->rValue = model->BSIM4pclc;
            return(OK);
        case BSIM4_MOD_PCLE:
            value->rValue = model->BSIM4pcle;
            return(OK);
        case BSIM4_MOD_PVFBCV:
            value->rValue = model->BSIM4pvfbcv;
            return(OK);
        case BSIM4_MOD_PACDE:
            value->rValue = model->BSIM4pacde;
            return(OK);
        case BSIM4_MOD_PMOIN:
            value->rValue = model->BSIM4pmoin;
            return(OK);
        case BSIM4_MOD_PNOFF:
            value->rValue = model->BSIM4pnoff;
            return(OK);
        case BSIM4_MOD_PVOFFCV:
            value->rValue = model->BSIM4pvoffcv;
            return(OK);
        case BSIM4_MOD_PVFBSDOFF:
            value->rValue = model->BSIM4pvfbsdoff;
            return(OK);
        case BSIM4_MOD_PTVFBSDOFF:
            value->rValue = model->BSIM4ptvfbsdoff;
            return(OK);

        case BSIM4_MOD_PLAMBDA:
            value->rValue = model->BSIM4plambda;
            return(OK);
        case BSIM4_MOD_PVTL:
            value->rValue = model->BSIM4pvtl;
            return(OK);
        case BSIM4_MOD_PXN:
            value->rValue = model->BSIM4pxn;
            return(OK);

        case  BSIM4_MOD_TNOM :
          value->rValue = model->BSIM4tnom;
            return(OK);
        case BSIM4_MOD_CGSO:
            value->rValue = model->BSIM4cgso; 
            return(OK);
        case BSIM4_MOD_CGDO:
            value->rValue = model->BSIM4cgdo; 
            return(OK);
        case BSIM4_MOD_CGBO:
            value->rValue = model->BSIM4cgbo; 
            return(OK);
        case BSIM4_MOD_XPART:
            value->rValue = model->BSIM4xpart; 
            return(OK);
        case BSIM4_MOD_RSH:
            value->rValue = model->BSIM4sheetResistance; 
            return(OK);
        case BSIM4_MOD_JSS:
            value->rValue = model->BSIM4SjctSatCurDensity; 
            return(OK);
        case BSIM4_MOD_JSWS:
            value->rValue = model->BSIM4SjctSidewallSatCurDensity; 
            return(OK);
        case BSIM4_MOD_JSWGS:
            value->rValue = model->BSIM4SjctGateSidewallSatCurDensity;
            return(OK);
        case BSIM4_MOD_PBS:
            value->rValue = model->BSIM4SbulkJctPotential; 
            return(OK);
        case BSIM4_MOD_MJS:
            value->rValue = model->BSIM4SbulkJctBotGradingCoeff; 
            return(OK);
        case BSIM4_MOD_PBSWS:
            value->rValue = model->BSIM4SsidewallJctPotential; 
            return(OK);
        case BSIM4_MOD_MJSWS:
            value->rValue = model->BSIM4SbulkJctSideGradingCoeff; 
            return(OK);
        case BSIM4_MOD_CJS:
            value->rValue = model->BSIM4SunitAreaJctCap; 
            return(OK);
        case BSIM4_MOD_CJSWS:
            value->rValue = model->BSIM4SunitLengthSidewallJctCap; 
            return(OK);
        case BSIM4_MOD_PBSWGS:
            value->rValue = model->BSIM4SGatesidewallJctPotential; 
            return(OK);
        case BSIM4_MOD_MJSWGS:
            value->rValue = model->BSIM4SbulkJctGateSideGradingCoeff; 
            return(OK);
        case BSIM4_MOD_CJSWGS:
            value->rValue = model->BSIM4SunitLengthGateSidewallJctCap; 
            return(OK);
        case BSIM4_MOD_NJS:
            value->rValue = model->BSIM4SjctEmissionCoeff; 
            return(OK);
        case BSIM4_MOD_XTIS:
            value->rValue = model->BSIM4SjctTempExponent; 
            return(OK);
        case BSIM4_MOD_JSD:
            value->rValue = model->BSIM4DjctSatCurDensity;
            return(OK);
        case BSIM4_MOD_JSWD:
            value->rValue = model->BSIM4DjctSidewallSatCurDensity;
            return(OK);
        case BSIM4_MOD_JSWGD:
            value->rValue = model->BSIM4DjctGateSidewallSatCurDensity;
            return(OK);
        case BSIM4_MOD_PBD:
            value->rValue = model->BSIM4DbulkJctPotential;
            return(OK);
        case BSIM4_MOD_MJD:
            value->rValue = model->BSIM4DbulkJctBotGradingCoeff;
            return(OK);
        case BSIM4_MOD_PBSWD:
            value->rValue = model->BSIM4DsidewallJctPotential;
            return(OK);
        case BSIM4_MOD_MJSWD:
            value->rValue = model->BSIM4DbulkJctSideGradingCoeff;
            return(OK);
        case BSIM4_MOD_CJD:
            value->rValue = model->BSIM4DunitAreaJctCap;
            return(OK);
        case BSIM4_MOD_CJSWD:
            value->rValue = model->BSIM4DunitLengthSidewallJctCap;
            return(OK);
        case BSIM4_MOD_PBSWGD:
            value->rValue = model->BSIM4DGatesidewallJctPotential;
            return(OK);
        case BSIM4_MOD_MJSWGD:
            value->rValue = model->BSIM4DbulkJctGateSideGradingCoeff;
            return(OK);
        case BSIM4_MOD_CJSWGD:
            value->rValue = model->BSIM4DunitLengthGateSidewallJctCap;
            return(OK);
        case BSIM4_MOD_NJD:
            value->rValue = model->BSIM4DjctEmissionCoeff;
            return(OK);
        case BSIM4_MOD_XTID:
            value->rValue = model->BSIM4DjctTempExponent;
            return(OK);
        case BSIM4_MOD_LINTNOI:
            value->rValue = model->BSIM4lintnoi; 
            return(OK);
        case BSIM4_MOD_LINT:
            value->rValue = model->BSIM4Lint; 
            return(OK);
        case BSIM4_MOD_LL:
            value->rValue = model->BSIM4Ll;
            return(OK);
        case BSIM4_MOD_LLC:
            value->rValue = model->BSIM4Llc;
            return(OK);
        case BSIM4_MOD_LLN:
            value->rValue = model->BSIM4Lln;
            return(OK);
        case BSIM4_MOD_LW:
            value->rValue = model->BSIM4Lw;
            return(OK);
        case BSIM4_MOD_LWC:
            value->rValue = model->BSIM4Lwc;
            return(OK);
        case BSIM4_MOD_LWN:
            value->rValue = model->BSIM4Lwn;
            return(OK);
        case BSIM4_MOD_LWL:
            value->rValue = model->BSIM4Lwl;
            return(OK);
        case BSIM4_MOD_LWLC:
            value->rValue = model->BSIM4Lwlc;
            return(OK);
        case BSIM4_MOD_LMIN:
            value->rValue = model->BSIM4Lmin;
            return(OK);
        case BSIM4_MOD_LMAX:
            value->rValue = model->BSIM4Lmax;
            return(OK);
        case BSIM4_MOD_WINT:
            value->rValue = model->BSIM4Wint;
            return(OK);
        case BSIM4_MOD_WL:
            value->rValue = model->BSIM4Wl;
            return(OK);
        case BSIM4_MOD_WLC:
            value->rValue = model->BSIM4Wlc;
            return(OK);
        case BSIM4_MOD_WLN:
            value->rValue = model->BSIM4Wln;
            return(OK);
        case BSIM4_MOD_WW:
            value->rValue = model->BSIM4Ww;
            return(OK);
        case BSIM4_MOD_WWC:
            value->rValue = model->BSIM4Wwc;
            return(OK);
        case BSIM4_MOD_WWN:
            value->rValue = model->BSIM4Wwn;
            return(OK);
        case BSIM4_MOD_WWL:
            value->rValue = model->BSIM4Wwl;
            return(OK);
        case BSIM4_MOD_WWLC:
            value->rValue = model->BSIM4Wwlc;
            return(OK);
        case BSIM4_MOD_WMIN:
            value->rValue = model->BSIM4Wmin;
            return(OK);
        case BSIM4_MOD_WMAX:
            value->rValue = model->BSIM4Wmax;
            return(OK);

        /* stress effect */
        case BSIM4_MOD_SAREF:
            value->rValue = model->BSIM4saref;
            return(OK);
        case BSIM4_MOD_SBREF:
            value->rValue = model->BSIM4sbref;
            return(OK);
        case BSIM4_MOD_WLOD:
            value->rValue = model->BSIM4wlod;
            return(OK);
        case BSIM4_MOD_KU0:
            value->rValue = model->BSIM4ku0;
            return(OK);
        case BSIM4_MOD_KVSAT:
            value->rValue = model->BSIM4kvsat;
            return(OK);
        case BSIM4_MOD_KVTH0:
            value->rValue = model->BSIM4kvth0;
            return(OK);
        case BSIM4_MOD_TKU0:
            value->rValue = model->BSIM4tku0;
            return(OK);
        case BSIM4_MOD_LLODKU0:
            value->rValue = model->BSIM4llodku0;
            return(OK);
        case BSIM4_MOD_WLODKU0:
            value->rValue = model->BSIM4wlodku0;
            return(OK);
        case BSIM4_MOD_LLODVTH:
            value->rValue = model->BSIM4llodvth;
            return(OK);
        case BSIM4_MOD_WLODVTH:
            value->rValue = model->BSIM4wlodvth;
            return(OK);
        case BSIM4_MOD_LKU0:
            value->rValue = model->BSIM4lku0;
            return(OK);
        case BSIM4_MOD_WKU0:
            value->rValue = model->BSIM4wku0;
            return(OK);
        case BSIM4_MOD_PKU0:
            value->rValue = model->BSIM4pku0;
            return(OK);
        case BSIM4_MOD_LKVTH0:
            value->rValue = model->BSIM4lkvth0;
            return(OK);
        case BSIM4_MOD_WKVTH0:
            value->rValue = model->BSIM4wkvth0;
            return(OK);
        case BSIM4_MOD_PKVTH0:
            value->rValue = model->BSIM4pkvth0;
            return(OK);
        case BSIM4_MOD_STK2:
            value->rValue = model->BSIM4stk2;
            return(OK);
        case BSIM4_MOD_LODK2:
            value->rValue = model->BSIM4lodk2;
            return(OK);
        case BSIM4_MOD_STETA0:
            value->rValue = model->BSIM4steta0;
            return(OK);
        case BSIM4_MOD_LODETA0:
            value->rValue = model->BSIM4lodeta0;
            return(OK);

        /* Well Proximity Effect  */
        case BSIM4_MOD_WEB:
            value->rValue = model->BSIM4web;
            return(OK);
        case BSIM4_MOD_WEC:
            value->rValue = model->BSIM4wec;
            return(OK);
        case BSIM4_MOD_KVTH0WE:
            value->rValue = model->BSIM4kvth0we;
            return(OK);
        case BSIM4_MOD_K2WE:
            value->rValue = model->BSIM4k2we;
            return(OK);
        case BSIM4_MOD_KU0WE:
            value->rValue = model->BSIM4ku0we;
            return(OK);
        case BSIM4_MOD_SCREF:
            value->rValue = model->BSIM4scref;
            return(OK);
        case BSIM4_MOD_WPEMOD:
            value->rValue = model->BSIM4wpemod;
            return(OK);
        case BSIM4_MOD_LKVTH0WE:
            value->rValue = model->BSIM4lkvth0we;
            return(OK);
        case BSIM4_MOD_LK2WE:
            value->rValue = model->BSIM4lk2we;
            return(OK);
        case BSIM4_MOD_LKU0WE:
            value->rValue = model->BSIM4lku0we;
            return(OK);
        case BSIM4_MOD_WKVTH0WE:
            value->rValue = model->BSIM4wkvth0we;
            return(OK);
        case BSIM4_MOD_WK2WE:
            value->rValue = model->BSIM4wk2we;
            return(OK);
        case BSIM4_MOD_WKU0WE:
            value->rValue = model->BSIM4wku0we;
            return(OK);
        case BSIM4_MOD_PKVTH0WE:
            value->rValue = model->BSIM4pkvth0we;
            return(OK);
        case BSIM4_MOD_PK2WE:
            value->rValue = model->BSIM4pk2we;
            return(OK);
        case BSIM4_MOD_PKU0WE:
            value->rValue = model->BSIM4pku0we;
            return(OK);

        case BSIM4_MOD_NOIA:
            value->rValue = model->BSIM4oxideTrapDensityA;
            return(OK);
        case BSIM4_MOD_NOIB:
            value->rValue = model->BSIM4oxideTrapDensityB;
            return(OK);
        case BSIM4_MOD_NOIC:
            value->rValue = model->BSIM4oxideTrapDensityC;
            return(OK);
        case BSIM4_MOD_EM:
            value->rValue = model->BSIM4em;
            return(OK);
        case BSIM4_MOD_EF:
            value->rValue = model->BSIM4ef;
            return(OK);
        case BSIM4_MOD_AF:
            value->rValue = model->BSIM4af;
            return(OK);
        case BSIM4_MOD_KF:
            value->rValue = model->BSIM4kf;
            return(OK);

        case BSIM4_MOD_VGS_MAX:
            value->rValue = model->BSIM4vgsMax;
            return(OK);
        case BSIM4_MOD_VGD_MAX:
            value->rValue = model->BSIM4vgdMax;
            return(OK);
        case BSIM4_MOD_VGB_MAX:
            value->rValue = model->BSIM4vgbMax;
            return(OK);
        case BSIM4_MOD_VDS_MAX:
            value->rValue = model->BSIM4vdsMax;
            return(OK);
        case BSIM4_MOD_VBS_MAX:
            value->rValue = model->BSIM4vbsMax;
            return(OK);
        case BSIM4_MOD_VBD_MAX:
            value->rValue = model->BSIM4vbdMax;
            return(OK);

        default:
            return(E_BADPARM);
    }
    /* NOTREACHED */
}



