/**** BSIM4v7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4mpar.c of BSIM4v7.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Authors: 2006- Mohan Dunga, Ali Niknejad, Chenming Hu
 * Authors: 2007- Mohan Dunga, Wenwei Yang, Ali Niknejad, Chenming Hu
 * Authors: 2008- Wenwei Yang, Ali Niknejad, Chenming Hu 
 * Project Director: Prof. Chenming Hu.
 * Modified by Xuemei Xi, 04/06/2001.
 * Modified by Xuemei Xi, 10/05/2001.
 * Modified by Xuemei Xi, 11/15/2002.
 * Modified by Xuemei Xi, 05/09/2003.
 * Modified by Xuemei Xi, 03/04/2004.
 * Modified by Xuemei Xi, Mohan Dunga, 07/29/2005.
 * Modified by Mohan Dunga, 12/13/2006
 * Modified by Mohan Dunga, Wenwei Yang, 05/18/2007.
 * Modified by Wenwei Yang, 07/31/2008.
 * Modified by Tanvir Morshed, Darsen Lu 03/27/2011
 **********/

#include "ngspice/ngspice.h"
#include "bsim4v7def.h"
#include "ngspice/ifsim.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/const.h"

int
BSIM4v7mParam(
int param,
IFvalue *value,
GENmodel *inMod)
{
    BSIM4v7model *mod = (BSIM4v7model*)inMod;
    switch(param)
    {   case  BSIM4v7_MOD_MOBMOD :
            mod->BSIM4v7mobMod = value->iValue;
            mod->BSIM4v7mobModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BINUNIT :
            mod->BSIM4v7binUnit = value->iValue;
            mod->BSIM4v7binUnitGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PARAMCHK :
            mod->BSIM4v7paramChk = value->iValue;
            mod->BSIM4v7paramChkGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CVCHARGEMOD :
            mod->BSIM4v7cvchargeMod = value->iValue;
            mod->BSIM4v7cvchargeModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CAPMOD :
            mod->BSIM4v7capMod = value->iValue;
            mod->BSIM4v7capModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DIOMOD :
            mod->BSIM4v7dioMod = value->iValue;
            mod->BSIM4v7dioModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RDSMOD :
            mod->BSIM4v7rdsMod = value->iValue;
            mod->BSIM4v7rdsModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TRNQSMOD :
            mod->BSIM4v7trnqsMod = value->iValue;
            mod->BSIM4v7trnqsModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_ACNQSMOD :
            mod->BSIM4v7acnqsMod = value->iValue;
            mod->BSIM4v7acnqsModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBODYMOD :
            mod->BSIM4v7rbodyMod = value->iValue;
            mod->BSIM4v7rbodyModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RGATEMOD :
            mod->BSIM4v7rgateMod = value->iValue;
            mod->BSIM4v7rgateModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PERMOD :
            mod->BSIM4v7perMod = value->iValue;
            mod->BSIM4v7perModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_GEOMOD :
            mod->BSIM4v7geoMod = value->iValue;
            mod->BSIM4v7geoModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RGEOMOD :
            mod->BSIM4v7rgeoMod = value->iValue;
            mod->BSIM4v7rgeoModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_FNOIMOD :
            mod->BSIM4v7fnoiMod = value->iValue;
            mod->BSIM4v7fnoiModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNOIMOD :
            mod->BSIM4v7tnoiMod = value->iValue;
            mod->BSIM4v7tnoiModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MTRLMOD :
            mod->BSIM4v7mtrlMod = value->iValue;
            mod->BSIM4v7mtrlModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MTRLCOMPATMOD :
            mod->BSIM4v7mtrlCompatMod = value->iValue;
            mod->BSIM4v7mtrlCompatModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_GIDLMOD :        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7gidlMod = value->iValue;
            mod->BSIM4v7gidlModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_IGCMOD :
            mod->BSIM4v7igcMod = value->iValue;
            mod->BSIM4v7igcModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_IGBMOD :
            mod->BSIM4v7igbMod = value->iValue;
            mod->BSIM4v7igbModGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TEMPMOD :
            mod->BSIM4v7tempMod = value->iValue;
            mod->BSIM4v7tempModGiven = TRUE;
            break;

        case  BSIM4v7_MOD_VERSION :
            mod->BSIM4v7version = value->sValue;
            mod->BSIM4v7versionGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TOXREF :
            mod->BSIM4v7toxref = value->rValue;
            mod->BSIM4v7toxrefGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EOT :
            mod->BSIM4v7eot = value->rValue;
            mod->BSIM4v7eotGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VDDEOT :
            mod->BSIM4v7vddeot = value->rValue;
            mod->BSIM4v7vddeotGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TEMPEOT :
            mod->BSIM4v7tempeot = value->rValue;
            mod->BSIM4v7tempeotGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LEFFEOT :
            mod->BSIM4v7leffeot = value->rValue;
            mod->BSIM4v7leffeotGiven = TRUE;
            break;
         case  BSIM4v7_MOD_WEFFEOT :
            mod->BSIM4v7weffeot = value->rValue;
            mod->BSIM4v7weffeotGiven = TRUE;
            break;
         case  BSIM4v7_MOD_ADOS :
            mod->BSIM4v7ados = value->rValue;
            mod->BSIM4v7adosGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BDOS :
            mod->BSIM4v7bdos = value->rValue;
            mod->BSIM4v7bdosGiven = TRUE;
            break;
       case  BSIM4v7_MOD_TOXE :
            mod->BSIM4v7toxe = value->rValue;
            mod->BSIM4v7toxeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TOXP :
            mod->BSIM4v7toxp = value->rValue;
            mod->BSIM4v7toxpGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TOXM :
            mod->BSIM4v7toxm = value->rValue;
            mod->BSIM4v7toxmGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DTOX :
            mod->BSIM4v7dtox = value->rValue;
            mod->BSIM4v7dtoxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EPSROX :
            mod->BSIM4v7epsrox = value->rValue;
            mod->BSIM4v7epsroxGiven = TRUE;
            break;

        case  BSIM4v7_MOD_CDSC :
            mod->BSIM4v7cdsc = value->rValue;
            mod->BSIM4v7cdscGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CDSCB :
            mod->BSIM4v7cdscb = value->rValue;
            mod->BSIM4v7cdscbGiven = TRUE;
            break;

        case  BSIM4v7_MOD_CDSCD :
            mod->BSIM4v7cdscd = value->rValue;
            mod->BSIM4v7cdscdGiven = TRUE;
            break;

        case  BSIM4v7_MOD_CIT :
            mod->BSIM4v7cit = value->rValue;
            mod->BSIM4v7citGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NFACTOR :
            mod->BSIM4v7nfactor = value->rValue;
            mod->BSIM4v7nfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_XJ:
            mod->BSIM4v7xj = value->rValue;
            mod->BSIM4v7xjGiven = TRUE;
            break;
        case BSIM4v7_MOD_VSAT:
            mod->BSIM4v7vsat = value->rValue;
            mod->BSIM4v7vsatGiven = TRUE;
            break;
        case BSIM4v7_MOD_A0:
            mod->BSIM4v7a0 = value->rValue;
            mod->BSIM4v7a0Given = TRUE;
            break;
        
        case BSIM4v7_MOD_AGS:
            mod->BSIM4v7ags= value->rValue;
            mod->BSIM4v7agsGiven = TRUE;
            break;
        
        case BSIM4v7_MOD_A1:
            mod->BSIM4v7a1 = value->rValue;
            mod->BSIM4v7a1Given = TRUE;
            break;
        case BSIM4v7_MOD_A2:
            mod->BSIM4v7a2 = value->rValue;
            mod->BSIM4v7a2Given = TRUE;
            break;
        case BSIM4v7_MOD_AT:
            mod->BSIM4v7at = value->rValue;
            mod->BSIM4v7atGiven = TRUE;
            break;
        case BSIM4v7_MOD_KETA:
            mod->BSIM4v7keta = value->rValue;
            mod->BSIM4v7ketaGiven = TRUE;
            break;    
        case BSIM4v7_MOD_NSUB:
            mod->BSIM4v7nsub = value->rValue;
            mod->BSIM4v7nsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_PHIG:
            mod->BSIM4v7phig = value->rValue;
            mod->BSIM4v7phigGiven = TRUE;
            break;
        case BSIM4v7_MOD_EPSRGATE:
            mod->BSIM4v7epsrgate = value->rValue;
            mod->BSIM4v7epsrgateGiven = TRUE;
            break;
        case BSIM4v7_MOD_EASUB:
            mod->BSIM4v7easub = value->rValue;
            mod->BSIM4v7easubGiven = TRUE;
            break;
        case BSIM4v7_MOD_EPSRSUB:
            mod->BSIM4v7epsrsub = value->rValue;
            mod->BSIM4v7epsrsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_NI0SUB:
            mod->BSIM4v7ni0sub = value->rValue;
            mod->BSIM4v7ni0subGiven = TRUE;
            break;
        case BSIM4v7_MOD_BG0SUB:
            mod->BSIM4v7bg0sub = value->rValue;
            mod->BSIM4v7bg0subGiven = TRUE;
            break;
        case BSIM4v7_MOD_TBGASUB:
            mod->BSIM4v7tbgasub = value->rValue;
            mod->BSIM4v7tbgasubGiven = TRUE;
            break;
        case BSIM4v7_MOD_TBGBSUB:
            mod->BSIM4v7tbgbsub = value->rValue;
            mod->BSIM4v7tbgbsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_NDEP:
            mod->BSIM4v7ndep = value->rValue;
            mod->BSIM4v7ndepGiven = TRUE;
            if (mod->BSIM4v7ndep > 1.0e20)
                mod->BSIM4v7ndep *= 1.0e-6;
            break;
        case BSIM4v7_MOD_NSD:
            mod->BSIM4v7nsd = value->rValue;
            mod->BSIM4v7nsdGiven = TRUE;
            if (mod->BSIM4v7nsd > 1.0e23)
                mod->BSIM4v7nsd *= 1.0e-6;
            break;
        case BSIM4v7_MOD_NGATE:
            mod->BSIM4v7ngate = value->rValue;
            mod->BSIM4v7ngateGiven = TRUE;
            if (mod->BSIM4v7ngate > 1.0e23)
                mod->BSIM4v7ngate *= 1.0e-6;
            break;
        case BSIM4v7_MOD_GAMMA1:
            mod->BSIM4v7gamma1 = value->rValue;
            mod->BSIM4v7gamma1Given = TRUE;
            break;
        case BSIM4v7_MOD_GAMMA2:
            mod->BSIM4v7gamma2 = value->rValue;
            mod->BSIM4v7gamma2Given = TRUE;
            break;
        case BSIM4v7_MOD_VBX:
            mod->BSIM4v7vbx = value->rValue;
            mod->BSIM4v7vbxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VBM:
            mod->BSIM4v7vbm = value->rValue;
            mod->BSIM4v7vbmGiven = TRUE;
            break;
        case BSIM4v7_MOD_XT:
            mod->BSIM4v7xt = value->rValue;
            mod->BSIM4v7xtGiven = TRUE;
            break;
        case  BSIM4v7_MOD_K1:
            mod->BSIM4v7k1 = value->rValue;
            mod->BSIM4v7k1Given = TRUE;
            break;
        case  BSIM4v7_MOD_KT1:
            mod->BSIM4v7kt1 = value->rValue;
            mod->BSIM4v7kt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_KT1L:
            mod->BSIM4v7kt1l = value->rValue;
            mod->BSIM4v7kt1lGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KT2:
            mod->BSIM4v7kt2 = value->rValue;
            mod->BSIM4v7kt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_K2:
            mod->BSIM4v7k2 = value->rValue;
            mod->BSIM4v7k2Given = TRUE;
            break;
        case  BSIM4v7_MOD_K3:
            mod->BSIM4v7k3 = value->rValue;
            mod->BSIM4v7k3Given = TRUE;
            break;
        case  BSIM4v7_MOD_K3B:
            mod->BSIM4v7k3b = value->rValue;
            mod->BSIM4v7k3bGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LPE0:
            mod->BSIM4v7lpe0 = value->rValue;
            mod->BSIM4v7lpe0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LPEB:
            mod->BSIM4v7lpeb = value->rValue;
            mod->BSIM4v7lpebGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP0:
            mod->BSIM4v7dvtp0 = value->rValue;
            mod->BSIM4v7dvtp0Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP1:
            mod->BSIM4v7dvtp1 = value->rValue;
            mod->BSIM4v7dvtp1Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP2:     /* New DIBL/Rout */
            mod->BSIM4v7dvtp2 = value->rValue;
            mod->BSIM4v7dvtp2Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP3:
            mod->BSIM4v7dvtp3 = value->rValue;
            mod->BSIM4v7dvtp3Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP4:
            mod->BSIM4v7dvtp4 = value->rValue;
            mod->BSIM4v7dvtp4Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVTP5:
            mod->BSIM4v7dvtp5 = value->rValue;
            mod->BSIM4v7dvtp5Given = TRUE;
            break;        
        case  BSIM4v7_MOD_W0:
            mod->BSIM4v7w0 = value->rValue;
            mod->BSIM4v7w0Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVT0:               
            mod->BSIM4v7dvt0 = value->rValue;
            mod->BSIM4v7dvt0Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVT1:             
            mod->BSIM4v7dvt1 = value->rValue;
            mod->BSIM4v7dvt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVT2:             
            mod->BSIM4v7dvt2 = value->rValue;
            mod->BSIM4v7dvt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_DVT0W:               
            mod->BSIM4v7dvt0w = value->rValue;
            mod->BSIM4v7dvt0wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DVT1W:             
            mod->BSIM4v7dvt1w = value->rValue;
            mod->BSIM4v7dvt1wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DVT2W:             
            mod->BSIM4v7dvt2w = value->rValue;
            mod->BSIM4v7dvt2wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DROUT:             
            mod->BSIM4v7drout = value->rValue;
            mod->BSIM4v7droutGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DSUB:             
            mod->BSIM4v7dsub = value->rValue;
            mod->BSIM4v7dsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_VTH0:
            mod->BSIM4v7vth0 = value->rValue;
            mod->BSIM4v7vth0Given = TRUE;
            break;
        case BSIM4v7_MOD_EU:
            mod->BSIM4v7eu = value->rValue;
            mod->BSIM4v7euGiven = TRUE;
            break;
                case BSIM4v7_MOD_UCS:
            mod->BSIM4v7ucs = value->rValue;
            mod->BSIM4v7ucsGiven = TRUE;
            break;
        case BSIM4v7_MOD_UA:
            mod->BSIM4v7ua = value->rValue;
            mod->BSIM4v7uaGiven = TRUE;
            break;
        case BSIM4v7_MOD_UA1:
            mod->BSIM4v7ua1 = value->rValue;
            mod->BSIM4v7ua1Given = TRUE;
            break;
        case BSIM4v7_MOD_UB:
            mod->BSIM4v7ub = value->rValue;
            mod->BSIM4v7ubGiven = TRUE;
            break;
        case BSIM4v7_MOD_UB1:
            mod->BSIM4v7ub1 = value->rValue;
            mod->BSIM4v7ub1Given = TRUE;
            break;
        case BSIM4v7_MOD_UC:
            mod->BSIM4v7uc = value->rValue;
            mod->BSIM4v7ucGiven = TRUE;
            break;
        case BSIM4v7_MOD_UC1:
            mod->BSIM4v7uc1 = value->rValue;
            mod->BSIM4v7uc1Given = TRUE;
            break;
        case  BSIM4v7_MOD_U0 :
            mod->BSIM4v7u0 = value->rValue;
            mod->BSIM4v7u0Given = TRUE;
            break;
        case  BSIM4v7_MOD_UTE :
            mod->BSIM4v7ute = value->rValue;
            mod->BSIM4v7uteGiven = TRUE;
            break;
        case  BSIM4v7_MOD_UCSTE :
            mod->BSIM4v7ucste = value->rValue;
            mod->BSIM4v7ucsteGiven = TRUE;
            break;
        case BSIM4v7_MOD_UD:
            mod->BSIM4v7ud = value->rValue;
            mod->BSIM4v7udGiven = TRUE;
            break;
        case BSIM4v7_MOD_UD1:
            mod->BSIM4v7ud1 = value->rValue;
            mod->BSIM4v7ud1Given = TRUE;
            break;
        case BSIM4v7_MOD_UP:
            mod->BSIM4v7up = value->rValue;
            mod->BSIM4v7upGiven = TRUE;
            break;
        case BSIM4v7_MOD_LP:
            mod->BSIM4v7lp = value->rValue;
            mod->BSIM4v7lpGiven = TRUE;
            break;
        case BSIM4v7_MOD_LUD:
            mod->BSIM4v7lud = value->rValue;
            mod->BSIM4v7ludGiven = TRUE;
            break;
        case BSIM4v7_MOD_LUD1:
            mod->BSIM4v7lud1 = value->rValue;
            mod->BSIM4v7lud1Given = TRUE;
            break;
        case BSIM4v7_MOD_LUP:
            mod->BSIM4v7lup = value->rValue;
            mod->BSIM4v7lupGiven = TRUE;
            break;
        case BSIM4v7_MOD_LLP:
            mod->BSIM4v7llp = value->rValue;
            mod->BSIM4v7llpGiven = TRUE;
            break;
        case BSIM4v7_MOD_WUD:
            mod->BSIM4v7wud = value->rValue;
            mod->BSIM4v7wudGiven = TRUE;
            break;
        case BSIM4v7_MOD_WUD1:
            mod->BSIM4v7wud1 = value->rValue;
            mod->BSIM4v7wud1Given = TRUE;
            break;
        case BSIM4v7_MOD_WUP:
            mod->BSIM4v7wup = value->rValue;
            mod->BSIM4v7wupGiven = TRUE;
            break;
        case BSIM4v7_MOD_WLP:
            mod->BSIM4v7wlp = value->rValue;
            mod->BSIM4v7wlpGiven = TRUE;
            break;
        case BSIM4v7_MOD_PUD:
            mod->BSIM4v7pud = value->rValue;
            mod->BSIM4v7pudGiven = TRUE;
            break;
        case BSIM4v7_MOD_PUD1:
            mod->BSIM4v7pud1 = value->rValue;
            mod->BSIM4v7pud1Given = TRUE;
            break;
        case BSIM4v7_MOD_PUP:
            mod->BSIM4v7pup = value->rValue;
            mod->BSIM4v7pupGiven = TRUE;
            break;
        case BSIM4v7_MOD_PLP:
            mod->BSIM4v7plp = value->rValue;
            mod->BSIM4v7plpGiven = TRUE;
            break;


        case BSIM4v7_MOD_VOFF:
            mod->BSIM4v7voff = value->rValue;
            mod->BSIM4v7voffGiven = TRUE;
            break;
        case BSIM4v7_MOD_TVOFF:
            mod->BSIM4v7tvoff = value->rValue;
            mod->BSIM4v7tvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_TNFACTOR:           /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7tnfactor = value->rValue;
            mod->BSIM4v7tnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_TETA0:                /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7teta0 = value->rValue;
            mod->BSIM4v7teta0Given = TRUE;
            break;
        case BSIM4v7_MOD_TVOFFCV:                /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7tvoffcv = value->rValue;
            mod->BSIM4v7tvoffcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_VOFFL:
            mod->BSIM4v7voffl = value->rValue;
            mod->BSIM4v7vofflGiven = TRUE;
            break;
        case BSIM4v7_MOD_VOFFCVL:
            mod->BSIM4v7voffcvl = value->rValue;
            mod->BSIM4v7voffcvlGiven = TRUE;
            break;
        case BSIM4v7_MOD_MINV:
            mod->BSIM4v7minv = value->rValue;
            mod->BSIM4v7minvGiven = TRUE;
            break;
        case BSIM4v7_MOD_MINVCV:
            mod->BSIM4v7minvcv = value->rValue;
            mod->BSIM4v7minvcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_FPROUT:
            mod->BSIM4v7fprout = value->rValue;
            mod->BSIM4v7fproutGiven = TRUE;
            break;
        case BSIM4v7_MOD_PDITS:
            mod->BSIM4v7pdits = value->rValue;
            mod->BSIM4v7pditsGiven = TRUE;
            break;
        case BSIM4v7_MOD_PDITSD:
            mod->BSIM4v7pditsd = value->rValue;
            mod->BSIM4v7pditsdGiven = TRUE;
            break;
        case BSIM4v7_MOD_PDITSL:
            mod->BSIM4v7pditsl = value->rValue;
            mod->BSIM4v7pditslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DELTA :
            mod->BSIM4v7delta = value->rValue;
            mod->BSIM4v7deltaGiven = TRUE;
            break;
        case BSIM4v7_MOD_RDSW:
            mod->BSIM4v7rdsw = value->rValue;
            mod->BSIM4v7rdswGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_RDSWMIN:
            mod->BSIM4v7rdswmin = value->rValue;
            mod->BSIM4v7rdswminGiven = TRUE;
            break;
        case BSIM4v7_MOD_RDWMIN:
            mod->BSIM4v7rdwmin = value->rValue;
            mod->BSIM4v7rdwminGiven = TRUE;
            break;
        case BSIM4v7_MOD_RSWMIN:
            mod->BSIM4v7rswmin = value->rValue;
            mod->BSIM4v7rswminGiven = TRUE;
            break;
        case BSIM4v7_MOD_RDW:
            mod->BSIM4v7rdw = value->rValue;
            mod->BSIM4v7rdwGiven = TRUE;
            break;
        case BSIM4v7_MOD_RSW:
            mod->BSIM4v7rsw = value->rValue;
            mod->BSIM4v7rswGiven = TRUE;
            break;
        case BSIM4v7_MOD_PRWG:
            mod->BSIM4v7prwg = value->rValue;
            mod->BSIM4v7prwgGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PRWB:
            mod->BSIM4v7prwb = value->rValue;
            mod->BSIM4v7prwbGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PRT:
            mod->BSIM4v7prt = value->rValue;
            mod->BSIM4v7prtGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_ETA0:
            mod->BSIM4v7eta0 = value->rValue;
            mod->BSIM4v7eta0Given = TRUE;
            break;                 
        case BSIM4v7_MOD_ETAB:
            mod->BSIM4v7etab = value->rValue;
            mod->BSIM4v7etabGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PCLM:
            mod->BSIM4v7pclm = value->rValue;
            mod->BSIM4v7pclmGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PDIBL1:
            mod->BSIM4v7pdibl1 = value->rValue;
            mod->BSIM4v7pdibl1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PDIBL2:
            mod->BSIM4v7pdibl2 = value->rValue;
            mod->BSIM4v7pdibl2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PDIBLB:
            mod->BSIM4v7pdiblb = value->rValue;
            mod->BSIM4v7pdiblbGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PSCBE1:
            mod->BSIM4v7pscbe1 = value->rValue;
            mod->BSIM4v7pscbe1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PSCBE2:
            mod->BSIM4v7pscbe2 = value->rValue;
            mod->BSIM4v7pscbe2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PVAG:
            mod->BSIM4v7pvag = value->rValue;
            mod->BSIM4v7pvagGiven = TRUE;
            break;                 
        case  BSIM4v7_MOD_WR :
            mod->BSIM4v7wr = value->rValue;
            mod->BSIM4v7wrGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DWG :
            mod->BSIM4v7dwg = value->rValue;
            mod->BSIM4v7dwgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DWB :
            mod->BSIM4v7dwb = value->rValue;
            mod->BSIM4v7dwbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_B0 :
            mod->BSIM4v7b0 = value->rValue;
            mod->BSIM4v7b0Given = TRUE;
            break;
        case  BSIM4v7_MOD_B1 :
            mod->BSIM4v7b1 = value->rValue;
            mod->BSIM4v7b1Given = TRUE;
            break;
        case  BSIM4v7_MOD_ALPHA0 :
            mod->BSIM4v7alpha0 = value->rValue;
            mod->BSIM4v7alpha0Given = TRUE;
            break;
        case  BSIM4v7_MOD_ALPHA1 :
            mod->BSIM4v7alpha1 = value->rValue;
            mod->BSIM4v7alpha1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PHIN :
            mod->BSIM4v7phin = value->rValue;
            mod->BSIM4v7phinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AGIDL :
            mod->BSIM4v7agidl = value->rValue;
            mod->BSIM4v7agidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BGIDL :
            mod->BSIM4v7bgidl = value->rValue;
            mod->BSIM4v7bgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGIDL :
            mod->BSIM4v7cgidl = value->rValue;
            mod->BSIM4v7cgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EGIDL :
            mod->BSIM4v7egidl = value->rValue;
            mod->BSIM4v7egidlGiven = TRUE;
            break;
          case  BSIM4v7_MOD_FGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7fgidl = value->rValue;
            mod->BSIM4v7fgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7kgidl = value->rValue;
            mod->BSIM4v7kgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7rgidl = value->rValue;
            mod->BSIM4v7rgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AGISL :
            mod->BSIM4v7agisl = value->rValue;
            mod->BSIM4v7agislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BGISL :
            mod->BSIM4v7bgisl = value->rValue;
            mod->BSIM4v7bgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGISL :
            mod->BSIM4v7cgisl = value->rValue;
            mod->BSIM4v7cgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EGISL :
            mod->BSIM4v7egisl = value->rValue;
            mod->BSIM4v7egislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_FGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7fgisl = value->rValue;
            mod->BSIM4v7fgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7kgisl = value->rValue;
            mod->BSIM4v7kgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7rgisl = value->rValue;
            mod->BSIM4v7rgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGC :
            mod->BSIM4v7aigc = value->rValue;
            mod->BSIM4v7aigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGC :
            mod->BSIM4v7bigc = value->rValue;
            mod->BSIM4v7bigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGC :
            mod->BSIM4v7cigc = value->rValue;
            mod->BSIM4v7cigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGSD :
            mod->BSIM4v7aigsd = value->rValue;
            mod->BSIM4v7aigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGSD :
            mod->BSIM4v7bigsd = value->rValue;
            mod->BSIM4v7bigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGSD :
            mod->BSIM4v7cigsd = value->rValue;
            mod->BSIM4v7cigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGS :
            mod->BSIM4v7aigs = value->rValue;
            mod->BSIM4v7aigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGS :
            mod->BSIM4v7bigs = value->rValue;
            mod->BSIM4v7bigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGS :
            mod->BSIM4v7cigs = value->rValue;
            mod->BSIM4v7cigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGD :
            mod->BSIM4v7aigd = value->rValue;
            mod->BSIM4v7aigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGD :
            mod->BSIM4v7bigd = value->rValue;
            mod->BSIM4v7bigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGD :
            mod->BSIM4v7cigd = value->rValue;
            mod->BSIM4v7cigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGBACC :
            mod->BSIM4v7aigbacc = value->rValue;
            mod->BSIM4v7aigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGBACC :
            mod->BSIM4v7bigbacc = value->rValue;
            mod->BSIM4v7bigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGBACC :
            mod->BSIM4v7cigbacc = value->rValue;
            mod->BSIM4v7cigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AIGBINV :
            mod->BSIM4v7aigbinv = value->rValue;
            mod->BSIM4v7aigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BIGBINV :
            mod->BSIM4v7bigbinv = value->rValue;
            mod->BSIM4v7bigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CIGBINV :
            mod->BSIM4v7cigbinv = value->rValue;
            mod->BSIM4v7cigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NIGC :
            mod->BSIM4v7nigc = value->rValue;
            mod->BSIM4v7nigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NIGBINV :
            mod->BSIM4v7nigbinv = value->rValue;
            mod->BSIM4v7nigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NIGBACC :
            mod->BSIM4v7nigbacc = value->rValue;
            mod->BSIM4v7nigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NTOX :
            mod->BSIM4v7ntox = value->rValue;
            mod->BSIM4v7ntoxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EIGBINV :
            mod->BSIM4v7eigbinv = value->rValue;
            mod->BSIM4v7eigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PIGCD :
            mod->BSIM4v7pigcd = value->rValue;
            mod->BSIM4v7pigcdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_POXEDGE :
            mod->BSIM4v7poxedge = value->rValue;
            mod->BSIM4v7poxedgeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XRCRG1 :
            mod->BSIM4v7xrcrg1 = value->rValue;
            mod->BSIM4v7xrcrg1Given = TRUE;
            break;
        case  BSIM4v7_MOD_XRCRG2 :
            mod->BSIM4v7xrcrg2 = value->rValue;
            mod->BSIM4v7xrcrg2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LAMBDA :
            mod->BSIM4v7lambda = value->rValue;
            mod->BSIM4v7lambdaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTL :
            mod->BSIM4v7vtl = value->rValue;
            mod->BSIM4v7vtlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XN:
            mod->BSIM4v7xn = value->rValue;
            mod->BSIM4v7xnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LC:
            mod->BSIM4v7lc = value->rValue;
            mod->BSIM4v7lcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNOIA :
            mod->BSIM4v7tnoia = value->rValue;
            mod->BSIM4v7tnoiaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNOIB :
            mod->BSIM4v7tnoib = value->rValue;
            mod->BSIM4v7tnoibGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNOIC :
            mod->BSIM4v7tnoic = value->rValue;
            mod->BSIM4v7tnoicGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RNOIA :
            mod->BSIM4v7rnoia = value->rValue;
            mod->BSIM4v7rnoiaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RNOIB :
            mod->BSIM4v7rnoib = value->rValue;
            mod->BSIM4v7rnoibGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RNOIC :
            mod->BSIM4v7rnoic = value->rValue;
            mod->BSIM4v7rnoicGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NTNOI :
            mod->BSIM4v7ntnoi = value->rValue;
            mod->BSIM4v7ntnoiGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VFBSDOFF:
            mod->BSIM4v7vfbsdoff = value->rValue;
            mod->BSIM4v7vfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TVFBSDOFF:
            mod->BSIM4v7tvfbsdoff = value->rValue;
            mod->BSIM4v7tvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LINTNOI:
            mod->BSIM4v7lintnoi = value->rValue;
            mod->BSIM4v7lintnoiGiven = TRUE;
            break;

        /* stress effect */
        case  BSIM4v7_MOD_SAREF :
            mod->BSIM4v7saref = value->rValue;
            mod->BSIM4v7sarefGiven = TRUE;
            break;
        case  BSIM4v7_MOD_SBREF :
            mod->BSIM4v7sbref = value->rValue;
            mod->BSIM4v7sbrefGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WLOD :
            mod->BSIM4v7wlod = value->rValue;
            mod->BSIM4v7wlodGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KU0 :
            mod->BSIM4v7ku0 = value->rValue;
            mod->BSIM4v7ku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_KVSAT :
            mod->BSIM4v7kvsat = value->rValue;
            mod->BSIM4v7kvsatGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KVTH0 :
            mod->BSIM4v7kvth0 = value->rValue;
            mod->BSIM4v7kvth0Given = TRUE;
            break;
        case  BSIM4v7_MOD_TKU0 :
            mod->BSIM4v7tku0 = value->rValue;
            mod->BSIM4v7tku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LLODKU0 :
            mod->BSIM4v7llodku0 = value->rValue;
            mod->BSIM4v7llodku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WLODKU0 :
            mod->BSIM4v7wlodku0 = value->rValue;
            mod->BSIM4v7wlodku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LLODVTH :
            mod->BSIM4v7llodvth = value->rValue;
            mod->BSIM4v7llodvthGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WLODVTH :
            mod->BSIM4v7wlodvth = value->rValue;
            mod->BSIM4v7wlodvthGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKU0 :
            mod->BSIM4v7lku0 = value->rValue;
            mod->BSIM4v7lku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WKU0 :
            mod->BSIM4v7wku0 = value->rValue;
            mod->BSIM4v7wku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PKU0 :
            mod->BSIM4v7pku0 = value->rValue;
            mod->BSIM4v7pku0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LKVTH0 :
            mod->BSIM4v7lkvth0 = value->rValue;
            mod->BSIM4v7lkvth0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WKVTH0 :
            mod->BSIM4v7wkvth0 = value->rValue;
            mod->BSIM4v7wkvth0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PKVTH0 :
            mod->BSIM4v7pkvth0 = value->rValue;
            mod->BSIM4v7pkvth0Given = TRUE;
            break;
        case  BSIM4v7_MOD_STK2 :
            mod->BSIM4v7stk2 = value->rValue;
            mod->BSIM4v7stk2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LODK2 :
            mod->BSIM4v7lodk2 = value->rValue;
            mod->BSIM4v7lodk2Given = TRUE;
            break;
        case  BSIM4v7_MOD_STETA0 :
            mod->BSIM4v7steta0 = value->rValue;
            mod->BSIM4v7steta0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LODETA0 :
            mod->BSIM4v7lodeta0 = value->rValue;
            mod->BSIM4v7lodeta0Given = TRUE;
            break;

        case  BSIM4v7_MOD_WEB :
            mod->BSIM4v7web = value->rValue;
            mod->BSIM4v7webGiven = TRUE;
            break;
        case BSIM4v7_MOD_WEC :
            mod->BSIM4v7wec = value->rValue;
            mod->BSIM4v7wecGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KVTH0WE :
            mod->BSIM4v7kvth0we = value->rValue;
            mod->BSIM4v7kvth0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_K2WE :
            mod->BSIM4v7k2we = value->rValue;
            mod->BSIM4v7k2weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KU0WE :
            mod->BSIM4v7ku0we = value->rValue;
            mod->BSIM4v7ku0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_SCREF :
            mod->BSIM4v7scref = value->rValue;
            mod->BSIM4v7screfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WPEMOD :
            mod->BSIM4v7wpemod = value->rValue;
            mod->BSIM4v7wpemodGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKVTH0WE :
            mod->BSIM4v7lkvth0we = value->rValue;
            mod->BSIM4v7lkvth0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LK2WE :
            mod->BSIM4v7lk2we = value->rValue;
            mod->BSIM4v7lk2weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKU0WE :
            mod->BSIM4v7lku0we = value->rValue;
            mod->BSIM4v7lku0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WKVTH0WE :
            mod->BSIM4v7wkvth0we = value->rValue;
            mod->BSIM4v7wkvth0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WK2WE :
            mod->BSIM4v7wk2we = value->rValue;
            mod->BSIM4v7wk2weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WKU0WE :
            mod->BSIM4v7wku0we = value->rValue;
            mod->BSIM4v7wku0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PKVTH0WE :
            mod->BSIM4v7pkvth0we = value->rValue;
            mod->BSIM4v7pkvth0weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PK2WE :
            mod->BSIM4v7pk2we = value->rValue;
            mod->BSIM4v7pk2weGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PKU0WE :
            mod->BSIM4v7pku0we = value->rValue;
            mod->BSIM4v7pku0weGiven = TRUE;
            break;

        case  BSIM4v7_MOD_BETA0 :
            mod->BSIM4v7beta0 = value->rValue;
            mod->BSIM4v7beta0Given = TRUE;
            break;
        case  BSIM4v7_MOD_IJTHDFWD :
            mod->BSIM4v7ijthdfwd = value->rValue;
            mod->BSIM4v7ijthdfwdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_IJTHSFWD :
            mod->BSIM4v7ijthsfwd = value->rValue;
            mod->BSIM4v7ijthsfwdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_IJTHDREV :
            mod->BSIM4v7ijthdrev = value->rValue;
            mod->BSIM4v7ijthdrevGiven = TRUE;
            break;
        case  BSIM4v7_MOD_IJTHSREV :
            mod->BSIM4v7ijthsrev = value->rValue;
            mod->BSIM4v7ijthsrevGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XJBVD :
            mod->BSIM4v7xjbvd = value->rValue;
            mod->BSIM4v7xjbvdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XJBVS :
            mod->BSIM4v7xjbvs = value->rValue;
            mod->BSIM4v7xjbvsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BVD :
            mod->BSIM4v7bvd = value->rValue;
            mod->BSIM4v7bvdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_BVS :
            mod->BSIM4v7bvs = value->rValue;
            mod->BSIM4v7bvsGiven = TRUE;
            break;
        
        /* reverse diode */
        case  BSIM4v7_MOD_JTSS :
            mod->BSIM4v7jtss = value->rValue;
            mod->BSIM4v7jtssGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JTSD :
            mod->BSIM4v7jtsd = value->rValue;
            mod->BSIM4v7jtsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JTSSWS :
            mod->BSIM4v7jtssws = value->rValue;
            mod->BSIM4v7jtsswsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JTSSWD :
            mod->BSIM4v7jtsswd = value->rValue;
            mod->BSIM4v7jtsswdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JTSSWGS :
            mod->BSIM4v7jtsswgs = value->rValue;
            mod->BSIM4v7jtsswgsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JTSSWGD :
            mod->BSIM4v7jtsswgd = value->rValue;
            mod->BSIM4v7jtsswgdGiven = TRUE;
            break;
        case BSIM4v7_MOD_JTWEFF :
            mod->BSIM4v7jtweff = value->rValue;
            mod->BSIM4v7jtweffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTS :
            mod->BSIM4v7njts = value->rValue;
            mod->BSIM4v7njtsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTSSW :
            mod->BSIM4v7njtssw = value->rValue;
            mod->BSIM4v7njtsswGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTSSWG :
            mod->BSIM4v7njtsswg = value->rValue;
            mod->BSIM4v7njtsswgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTSD :
            mod->BSIM4v7njtsd = value->rValue;
            mod->BSIM4v7njtsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTSSWD :
            mod->BSIM4v7njtsswd = value->rValue;
            mod->BSIM4v7njtsswdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJTSSWGD :
            mod->BSIM4v7njtsswgd = value->rValue;
            mod->BSIM4v7njtsswgdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSS :
            mod->BSIM4v7xtss = value->rValue;
            mod->BSIM4v7xtssGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSD :
            mod->BSIM4v7xtsd = value->rValue;
            mod->BSIM4v7xtsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSSWS :
            mod->BSIM4v7xtssws = value->rValue;
            mod->BSIM4v7xtsswsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSSWD :
            mod->BSIM4v7xtsswd = value->rValue;
            mod->BSIM4v7xtsswdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSSWGS :
            mod->BSIM4v7xtsswgs = value->rValue;
            mod->BSIM4v7xtsswgsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTSSWGD :
            mod->BSIM4v7xtsswgd = value->rValue;
            mod->BSIM4v7xtsswgdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTS :
            mod->BSIM4v7tnjts = value->rValue;
            mod->BSIM4v7tnjtsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTSSW :
            mod->BSIM4v7tnjtssw = value->rValue;
            mod->BSIM4v7tnjtsswGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTSSWG :
            mod->BSIM4v7tnjtsswg = value->rValue;
            mod->BSIM4v7tnjtsswgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTSD :
            mod->BSIM4v7tnjtsd = value->rValue;
            mod->BSIM4v7tnjtsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTSSWD :
            mod->BSIM4v7tnjtsswd = value->rValue;
            mod->BSIM4v7tnjtsswdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TNJTSSWGD :
            mod->BSIM4v7tnjtsswgd = value->rValue;
            mod->BSIM4v7tnjtsswgdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSS :
            mod->BSIM4v7vtss = value->rValue;
            mod->BSIM4v7vtssGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSD :
            mod->BSIM4v7vtsd = value->rValue;
            mod->BSIM4v7vtsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSSWS :
            mod->BSIM4v7vtssws = value->rValue;
            mod->BSIM4v7vtsswsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSSWD :
            mod->BSIM4v7vtsswd = value->rValue;
            mod->BSIM4v7vtsswdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSSWGS :
            mod->BSIM4v7vtsswgs = value->rValue;
            mod->BSIM4v7vtsswgsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VTSSWGD :
            mod->BSIM4v7vtsswgd = value->rValue;
            mod->BSIM4v7vtsswgdGiven = TRUE;
            break;

        case  BSIM4v7_MOD_VFB :
            mod->BSIM4v7vfb = value->rValue;
            mod->BSIM4v7vfbGiven = TRUE;
            break;

        case  BSIM4v7_MOD_GBMIN :
            mod->BSIM4v7gbmin = value->rValue;
            mod->BSIM4v7gbminGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBDB :
            mod->BSIM4v7rbdb = value->rValue;
            mod->BSIM4v7rbdbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPB :
            mod->BSIM4v7rbpb = value->rValue;
            mod->BSIM4v7rbpbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBSB :
            mod->BSIM4v7rbsb = value->rValue;
            mod->BSIM4v7rbsbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPS :
            mod->BSIM4v7rbps = value->rValue;
            mod->BSIM4v7rbpsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPD :
            mod->BSIM4v7rbpd = value->rValue;
            mod->BSIM4v7rbpdGiven = TRUE;
            break;

        case  BSIM4v7_MOD_RBPS0 :
            mod->BSIM4v7rbps0 = value->rValue;
            mod->BSIM4v7rbps0Given = TRUE;
            break;
        case  BSIM4v7_MOD_RBPSL :
            mod->BSIM4v7rbpsl = value->rValue;
            mod->BSIM4v7rbpslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPSW :
            mod->BSIM4v7rbpsw = value->rValue;
            mod->BSIM4v7rbpswGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPSNF :
            mod->BSIM4v7rbpsnf = value->rValue;
            mod->BSIM4v7rbpsnfGiven = TRUE;
            break;

        case  BSIM4v7_MOD_RBPD0 :
            mod->BSIM4v7rbpd0 = value->rValue;
            mod->BSIM4v7rbpd0Given = TRUE;
            break;
        case  BSIM4v7_MOD_RBPDL :
            mod->BSIM4v7rbpdl = value->rValue;
            mod->BSIM4v7rbpdlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPDW :
            mod->BSIM4v7rbpdw = value->rValue;
            mod->BSIM4v7rbpdwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPDNF :
            mod->BSIM4v7rbpdnf = value->rValue;
            mod->BSIM4v7rbpdnfGiven = TRUE;
            break;

        case  BSIM4v7_MOD_RBPBX0 :
            mod->BSIM4v7rbpbx0 = value->rValue;
            mod->BSIM4v7rbpbx0Given = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBXL :
            mod->BSIM4v7rbpbxl = value->rValue;
            mod->BSIM4v7rbpbxlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBXW :
            mod->BSIM4v7rbpbxw = value->rValue;
            mod->BSIM4v7rbpbxwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBXNF :
            mod->BSIM4v7rbpbxnf = value->rValue;
            mod->BSIM4v7rbpbxnfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBY0 :
            mod->BSIM4v7rbpby0 = value->rValue;
            mod->BSIM4v7rbpby0Given = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBYL :
            mod->BSIM4v7rbpbyl = value->rValue;
            mod->BSIM4v7rbpbylGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBYW :
            mod->BSIM4v7rbpbyw = value->rValue;
            mod->BSIM4v7rbpbywGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RBPBYNF :
            mod->BSIM4v7rbpbynf = value->rValue;
            mod->BSIM4v7rbpbynfGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSBX0 :
            mod->BSIM4v7rbsbx0 = value->rValue;
            mod->BSIM4v7rbsbx0Given = TRUE;
            break;
       case  BSIM4v7_MOD_RBSBY0 :
            mod->BSIM4v7rbsby0 = value->rValue;
            mod->BSIM4v7rbsby0Given = TRUE;
            break;
       case  BSIM4v7_MOD_RBDBX0 :
            mod->BSIM4v7rbdbx0 = value->rValue;
            mod->BSIM4v7rbdbx0Given = TRUE;
            break;
       case  BSIM4v7_MOD_RBDBY0 :
            mod->BSIM4v7rbdby0 = value->rValue;
            mod->BSIM4v7rbdby0Given = TRUE;
            break;


       case  BSIM4v7_MOD_RBSDBXL :
            mod->BSIM4v7rbsdbxl = value->rValue;
            mod->BSIM4v7rbsdbxlGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSDBXW :
            mod->BSIM4v7rbsdbxw = value->rValue;
            mod->BSIM4v7rbsdbxwGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSDBXNF :
            mod->BSIM4v7rbsdbxnf = value->rValue;
            mod->BSIM4v7rbsdbxnfGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSDBYL :
            mod->BSIM4v7rbsdbyl = value->rValue;
            mod->BSIM4v7rbsdbylGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSDBYW :
            mod->BSIM4v7rbsdbyw = value->rValue;
            mod->BSIM4v7rbsdbywGiven = TRUE;
            break;
       case  BSIM4v7_MOD_RBSDBYNF :
            mod->BSIM4v7rbsdbynf = value->rValue;
            mod->BSIM4v7rbsdbynfGiven = TRUE;
            break;
 
        case  BSIM4v7_MOD_CGSL :
            mod->BSIM4v7cgsl = value->rValue;
            mod->BSIM4v7cgslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGDL :
            mod->BSIM4v7cgdl = value->rValue;
            mod->BSIM4v7cgdlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CKAPPAS :
            mod->BSIM4v7ckappas = value->rValue;
            mod->BSIM4v7ckappasGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CKAPPAD :
            mod->BSIM4v7ckappad = value->rValue;
            mod->BSIM4v7ckappadGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CF :
            mod->BSIM4v7cf = value->rValue;
            mod->BSIM4v7cfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CLC :
            mod->BSIM4v7clc = value->rValue;
            mod->BSIM4v7clcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CLE :
            mod->BSIM4v7cle = value->rValue;
            mod->BSIM4v7cleGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DWC :
            mod->BSIM4v7dwc = value->rValue;
            mod->BSIM4v7dwcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DLC :
            mod->BSIM4v7dlc = value->rValue;
            mod->BSIM4v7dlcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XW :
            mod->BSIM4v7xw = value->rValue;
            mod->BSIM4v7xwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XL :
            mod->BSIM4v7xl = value->rValue;
            mod->BSIM4v7xlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DLCIG :
            mod->BSIM4v7dlcig = value->rValue;
            mod->BSIM4v7dlcigGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DLCIGD :
            mod->BSIM4v7dlcigd = value->rValue;
            mod->BSIM4v7dlcigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DWJ :
            mod->BSIM4v7dwj = value->rValue;
            mod->BSIM4v7dwjGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VFBCV :
            mod->BSIM4v7vfbcv = value->rValue;
            mod->BSIM4v7vfbcvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_ACDE :
            mod->BSIM4v7acde = value->rValue;
            mod->BSIM4v7acdeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MOIN :
            mod->BSIM4v7moin = value->rValue;
            mod->BSIM4v7moinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NOFF :
            mod->BSIM4v7noff = value->rValue;
            mod->BSIM4v7noffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_VOFFCV :
            mod->BSIM4v7voffcv = value->rValue;
            mod->BSIM4v7voffcvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DMCG :
            mod->BSIM4v7dmcg = value->rValue;
            mod->BSIM4v7dmcgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DMCI :
            mod->BSIM4v7dmci = value->rValue;
            mod->BSIM4v7dmciGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DMDG :
            mod->BSIM4v7dmdg = value->rValue;
            mod->BSIM4v7dmdgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_DMCGT :
            mod->BSIM4v7dmcgt = value->rValue;
            mod->BSIM4v7dmcgtGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XGW :
            mod->BSIM4v7xgw = value->rValue;
            mod->BSIM4v7xgwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XGL :
            mod->BSIM4v7xgl = value->rValue;
            mod->BSIM4v7xglGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RSHG :
            mod->BSIM4v7rshg = value->rValue;
            mod->BSIM4v7rshgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NGCON :
            mod->BSIM4v7ngcon = value->rValue;
            mod->BSIM4v7ngconGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TCJ :
            mod->BSIM4v7tcj = value->rValue;
            mod->BSIM4v7tcjGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TPB :
            mod->BSIM4v7tpb = value->rValue;
            mod->BSIM4v7tpbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TCJSW :
            mod->BSIM4v7tcjsw = value->rValue;
            mod->BSIM4v7tcjswGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TPBSW :
            mod->BSIM4v7tpbsw = value->rValue;
            mod->BSIM4v7tpbswGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TCJSWG :
            mod->BSIM4v7tcjswg = value->rValue;
            mod->BSIM4v7tcjswgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_TPBSWG :
            mod->BSIM4v7tpbswg = value->rValue;
            mod->BSIM4v7tpbswgGiven = TRUE;
            break;

        /* Length dependence */
        case  BSIM4v7_MOD_LCDSC :
            mod->BSIM4v7lcdsc = value->rValue;
            mod->BSIM4v7lcdscGiven = TRUE;
            break;


        case  BSIM4v7_MOD_LCDSCB :
            mod->BSIM4v7lcdscb = value->rValue;
            mod->BSIM4v7lcdscbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCDSCD :
            mod->BSIM4v7lcdscd = value->rValue;
            mod->BSIM4v7lcdscdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIT :
            mod->BSIM4v7lcit = value->rValue;
            mod->BSIM4v7lcitGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNFACTOR :
            mod->BSIM4v7lnfactor = value->rValue;
            mod->BSIM4v7lnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_LXJ:
            mod->BSIM4v7lxj = value->rValue;
            mod->BSIM4v7lxjGiven = TRUE;
            break;
        case BSIM4v7_MOD_LVSAT:
            mod->BSIM4v7lvsat = value->rValue;
            mod->BSIM4v7lvsatGiven = TRUE;
            break;
        
        
        case BSIM4v7_MOD_LA0:
            mod->BSIM4v7la0 = value->rValue;
            mod->BSIM4v7la0Given = TRUE;
            break;
        case BSIM4v7_MOD_LAGS:
            mod->BSIM4v7lags = value->rValue;
            mod->BSIM4v7lagsGiven = TRUE;
            break;
        case BSIM4v7_MOD_LA1:
            mod->BSIM4v7la1 = value->rValue;
            mod->BSIM4v7la1Given = TRUE;
            break;
        case BSIM4v7_MOD_LA2:
            mod->BSIM4v7la2 = value->rValue;
            mod->BSIM4v7la2Given = TRUE;
            break;
        case BSIM4v7_MOD_LAT:
            mod->BSIM4v7lat = value->rValue;
            mod->BSIM4v7latGiven = TRUE;
            break;
        case BSIM4v7_MOD_LKETA:
            mod->BSIM4v7lketa = value->rValue;
            mod->BSIM4v7lketaGiven = TRUE;
            break;    
        case BSIM4v7_MOD_LNSUB:
            mod->BSIM4v7lnsub = value->rValue;
            mod->BSIM4v7lnsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_LNDEP:
            mod->BSIM4v7lndep = value->rValue;
            mod->BSIM4v7lndepGiven = TRUE;
            if (mod->BSIM4v7lndep > 1.0e20)
                mod->BSIM4v7lndep *= 1.0e-6;
            break;
        case BSIM4v7_MOD_LNSD:
            mod->BSIM4v7lnsd = value->rValue;
            mod->BSIM4v7lnsdGiven = TRUE;
            if (mod->BSIM4v7lnsd > 1.0e23)
                mod->BSIM4v7lnsd *= 1.0e-6;
            break;
        case BSIM4v7_MOD_LNGATE:
            mod->BSIM4v7lngate = value->rValue;
            mod->BSIM4v7lngateGiven = TRUE;
            if (mod->BSIM4v7lngate > 1.0e23)
                mod->BSIM4v7lngate *= 1.0e-6;
            break;
        case BSIM4v7_MOD_LGAMMA1:
            mod->BSIM4v7lgamma1 = value->rValue;
            mod->BSIM4v7lgamma1Given = TRUE;
            break;
        case BSIM4v7_MOD_LGAMMA2:
            mod->BSIM4v7lgamma2 = value->rValue;
            mod->BSIM4v7lgamma2Given = TRUE;
            break;
        case BSIM4v7_MOD_LVBX:
            mod->BSIM4v7lvbx = value->rValue;
            mod->BSIM4v7lvbxGiven = TRUE;
            break;
        case BSIM4v7_MOD_LVBM:
            mod->BSIM4v7lvbm = value->rValue;
            mod->BSIM4v7lvbmGiven = TRUE;
            break;
        case BSIM4v7_MOD_LXT:
            mod->BSIM4v7lxt = value->rValue;
            mod->BSIM4v7lxtGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LK1:
            mod->BSIM4v7lk1 = value->rValue;
            mod->BSIM4v7lk1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LKT1:
            mod->BSIM4v7lkt1 = value->rValue;
            mod->BSIM4v7lkt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LKT1L:
            mod->BSIM4v7lkt1l = value->rValue;
            mod->BSIM4v7lkt1lGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKT2:
            mod->BSIM4v7lkt2 = value->rValue;
            mod->BSIM4v7lkt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LK2:
            mod->BSIM4v7lk2 = value->rValue;
            mod->BSIM4v7lk2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LK3:
            mod->BSIM4v7lk3 = value->rValue;
            mod->BSIM4v7lk3Given = TRUE;
            break;
        case  BSIM4v7_MOD_LK3B:
            mod->BSIM4v7lk3b = value->rValue;
            mod->BSIM4v7lk3bGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LLPE0:
            mod->BSIM4v7llpe0 = value->rValue;
            mod->BSIM4v7llpe0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LLPEB:
            mod->BSIM4v7llpeb = value->rValue;
            mod->BSIM4v7llpebGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP0:
            mod->BSIM4v7ldvtp0 = value->rValue;
            mod->BSIM4v7ldvtp0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP1:
            mod->BSIM4v7ldvtp1 = value->rValue;
            mod->BSIM4v7ldvtp1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP2:     /* New DIBL/Rout */
            mod->BSIM4v7ldvtp2 = value->rValue;
            mod->BSIM4v7ldvtp2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP3:
            mod->BSIM4v7ldvtp3 = value->rValue;
            mod->BSIM4v7ldvtp3Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP4:
            mod->BSIM4v7ldvtp4 = value->rValue;
            mod->BSIM4v7ldvtp4Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVTP5:
            mod->BSIM4v7ldvtp5 = value->rValue;
            mod->BSIM4v7ldvtp5Given = TRUE;
            break;
        case  BSIM4v7_MOD_LW0:
            mod->BSIM4v7lw0 = value->rValue;
            mod->BSIM4v7lw0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT0:               
            mod->BSIM4v7ldvt0 = value->rValue;
            mod->BSIM4v7ldvt0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT1:             
            mod->BSIM4v7ldvt1 = value->rValue;
            mod->BSIM4v7ldvt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT2:             
            mod->BSIM4v7ldvt2 = value->rValue;
            mod->BSIM4v7ldvt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT0W:               
            mod->BSIM4v7ldvt0w = value->rValue;
            mod->BSIM4v7ldvt0wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT1W:             
            mod->BSIM4v7ldvt1w = value->rValue;
            mod->BSIM4v7ldvt1wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDVT2W:             
            mod->BSIM4v7ldvt2w = value->rValue;
            mod->BSIM4v7ldvt2wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDROUT:             
            mod->BSIM4v7ldrout = value->rValue;
            mod->BSIM4v7ldroutGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDSUB:             
            mod->BSIM4v7ldsub = value->rValue;
            mod->BSIM4v7ldsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_LVTH0:
            mod->BSIM4v7lvth0 = value->rValue;
            mod->BSIM4v7lvth0Given = TRUE;
            break;
        case BSIM4v7_MOD_LUA:
            mod->BSIM4v7lua = value->rValue;
            mod->BSIM4v7luaGiven = TRUE;
            break;
        case BSIM4v7_MOD_LUA1:
            mod->BSIM4v7lua1 = value->rValue;
            mod->BSIM4v7lua1Given = TRUE;
            break;
        case BSIM4v7_MOD_LUB:
            mod->BSIM4v7lub = value->rValue;
            mod->BSIM4v7lubGiven = TRUE;
            break;
        case BSIM4v7_MOD_LUB1:
            mod->BSIM4v7lub1 = value->rValue;
            mod->BSIM4v7lub1Given = TRUE;
            break;
        case BSIM4v7_MOD_LUC:
            mod->BSIM4v7luc = value->rValue;
            mod->BSIM4v7lucGiven = TRUE;
            break;
        case BSIM4v7_MOD_LUC1:
            mod->BSIM4v7luc1 = value->rValue;
            mod->BSIM4v7luc1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LU0 :
            mod->BSIM4v7lu0 = value->rValue;
            mod->BSIM4v7lu0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LUTE :
            mod->BSIM4v7lute = value->rValue;
            mod->BSIM4v7luteGiven = TRUE;
            break;
                case  BSIM4v7_MOD_LUCSTE :
            mod->BSIM4v7lucste = value->rValue;
            mod->BSIM4v7lucsteGiven = TRUE;
            break;
        case BSIM4v7_MOD_LVOFF:
            mod->BSIM4v7lvoff = value->rValue;
            mod->BSIM4v7lvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_LTVOFF:
            mod->BSIM4v7ltvoff = value->rValue;
            mod->BSIM4v7ltvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_LTNFACTOR:           /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7ltnfactor = value->rValue;
            mod->BSIM4v7ltnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_LTETA0:                /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7lteta0 = value->rValue;
            mod->BSIM4v7lteta0Given = TRUE;
            break;
        case BSIM4v7_MOD_LTVOFFCV:        /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7ltvoffcv = value->rValue;
            mod->BSIM4v7ltvoffcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_LMINV:
            mod->BSIM4v7lminv = value->rValue;
            mod->BSIM4v7lminvGiven = TRUE;
            break;
        case BSIM4v7_MOD_LMINVCV:
            mod->BSIM4v7lminvcv = value->rValue;
            mod->BSIM4v7lminvcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_LFPROUT:
            mod->BSIM4v7lfprout = value->rValue;
            mod->BSIM4v7lfproutGiven = TRUE;
            break;
        case BSIM4v7_MOD_LPDITS:
            mod->BSIM4v7lpdits = value->rValue;
            mod->BSIM4v7lpditsGiven = TRUE;
            break;
        case BSIM4v7_MOD_LPDITSD:
            mod->BSIM4v7lpditsd = value->rValue;
            mod->BSIM4v7lpditsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDELTA :
            mod->BSIM4v7ldelta = value->rValue;
            mod->BSIM4v7ldeltaGiven = TRUE;
            break;
        case BSIM4v7_MOD_LRDSW:
            mod->BSIM4v7lrdsw = value->rValue;
            mod->BSIM4v7lrdswGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_LRDW:
            mod->BSIM4v7lrdw = value->rValue;
            mod->BSIM4v7lrdwGiven = TRUE;
            break;
        case BSIM4v7_MOD_LRSW:
            mod->BSIM4v7lrsw = value->rValue;
            mod->BSIM4v7lrswGiven = TRUE;
            break;
        case BSIM4v7_MOD_LPRWB:
            mod->BSIM4v7lprwb = value->rValue;
            mod->BSIM4v7lprwbGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_LPRWG:
            mod->BSIM4v7lprwg = value->rValue;
            mod->BSIM4v7lprwgGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_LPRT:
            mod->BSIM4v7lprt = value->rValue;
            mod->BSIM4v7lprtGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_LETA0:
            mod->BSIM4v7leta0 = value->rValue;
            mod->BSIM4v7leta0Given = TRUE;
            break;                 
        case BSIM4v7_MOD_LETAB:
            mod->BSIM4v7letab = value->rValue;
            mod->BSIM4v7letabGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_LPCLM:
            mod->BSIM4v7lpclm = value->rValue;
            mod->BSIM4v7lpclmGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_LPDIBL1:
            mod->BSIM4v7lpdibl1 = value->rValue;
            mod->BSIM4v7lpdibl1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_LPDIBL2:
            mod->BSIM4v7lpdibl2 = value->rValue;
            mod->BSIM4v7lpdibl2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_LPDIBLB:
            mod->BSIM4v7lpdiblb = value->rValue;
            mod->BSIM4v7lpdiblbGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_LPSCBE1:
            mod->BSIM4v7lpscbe1 = value->rValue;
            mod->BSIM4v7lpscbe1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_LPSCBE2:
            mod->BSIM4v7lpscbe2 = value->rValue;
            mod->BSIM4v7lpscbe2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_LPVAG:
            mod->BSIM4v7lpvag = value->rValue;
            mod->BSIM4v7lpvagGiven = TRUE;
            break;                 
        case  BSIM4v7_MOD_LWR :
            mod->BSIM4v7lwr = value->rValue;
            mod->BSIM4v7lwrGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDWG :
            mod->BSIM4v7ldwg = value->rValue;
            mod->BSIM4v7ldwgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LDWB :
            mod->BSIM4v7ldwb = value->rValue;
            mod->BSIM4v7ldwbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LB0 :
            mod->BSIM4v7lb0 = value->rValue;
            mod->BSIM4v7lb0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LB1 :
            mod->BSIM4v7lb1 = value->rValue;
            mod->BSIM4v7lb1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LALPHA0 :
            mod->BSIM4v7lalpha0 = value->rValue;
            mod->BSIM4v7lalpha0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LALPHA1 :
            mod->BSIM4v7lalpha1 = value->rValue;
            mod->BSIM4v7lalpha1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LBETA0 :
            mod->BSIM4v7lbeta0 = value->rValue;
            mod->BSIM4v7lbeta0Given = TRUE;
            break;
        case  BSIM4v7_MOD_LPHIN :
            mod->BSIM4v7lphin = value->rValue;
            mod->BSIM4v7lphinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAGIDL :
            mod->BSIM4v7lagidl = value->rValue;
            mod->BSIM4v7lagidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBGIDL :
            mod->BSIM4v7lbgidl = value->rValue;
            mod->BSIM4v7lbgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCGIDL :
            mod->BSIM4v7lcgidl = value->rValue;
            mod->BSIM4v7lcgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LEGIDL :
            mod->BSIM4v7legidl = value->rValue;
            mod->BSIM4v7legidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LFGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lfgidl = value->rValue;
            mod->BSIM4v7lfgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lkgidl = value->rValue;
            mod->BSIM4v7lkgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LRGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lrgidl = value->rValue;
            mod->BSIM4v7lrgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAGISL :
            mod->BSIM4v7lagisl = value->rValue;
            mod->BSIM4v7lagislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBGISL :
            mod->BSIM4v7lbgisl = value->rValue;
            mod->BSIM4v7lbgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCGISL :
            mod->BSIM4v7lcgisl = value->rValue;
            mod->BSIM4v7lcgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LEGISL :
            mod->BSIM4v7legisl = value->rValue;
            mod->BSIM4v7legislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LFGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lfgisl = value->rValue;
            mod->BSIM4v7lfgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LKGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lkgisl = value->rValue;
            mod->BSIM4v7lkgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LRGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7lrgisl = value->rValue;
            mod->BSIM4v7lrgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGC :
            mod->BSIM4v7laigc = value->rValue;
            mod->BSIM4v7laigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGC :
            mod->BSIM4v7lbigc = value->rValue;
            mod->BSIM4v7lbigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGC :
            mod->BSIM4v7lcigc = value->rValue;
            mod->BSIM4v7lcigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGSD :
            mod->BSIM4v7laigsd = value->rValue;
            mod->BSIM4v7laigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGSD :
            mod->BSIM4v7lbigsd = value->rValue;
            mod->BSIM4v7lbigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGSD :
            mod->BSIM4v7lcigsd = value->rValue;
            mod->BSIM4v7lcigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGS :
            mod->BSIM4v7laigs = value->rValue;
            mod->BSIM4v7laigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGS :
            mod->BSIM4v7lbigs = value->rValue;
            mod->BSIM4v7lbigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGS :
            mod->BSIM4v7lcigs = value->rValue;
            mod->BSIM4v7lcigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGD :
            mod->BSIM4v7laigd = value->rValue;
            mod->BSIM4v7laigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGD :
            mod->BSIM4v7lbigd = value->rValue;
            mod->BSIM4v7lbigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGD :
            mod->BSIM4v7lcigd = value->rValue;
            mod->BSIM4v7lcigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGBACC :
            mod->BSIM4v7laigbacc = value->rValue;
            mod->BSIM4v7laigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGBACC :
            mod->BSIM4v7lbigbacc = value->rValue;
            mod->BSIM4v7lbigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGBACC :
            mod->BSIM4v7lcigbacc = value->rValue;
            mod->BSIM4v7lcigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LAIGBINV :
            mod->BSIM4v7laigbinv = value->rValue;
            mod->BSIM4v7laigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LBIGBINV :
            mod->BSIM4v7lbigbinv = value->rValue;
            mod->BSIM4v7lbigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCIGBINV :
            mod->BSIM4v7lcigbinv = value->rValue;
            mod->BSIM4v7lcigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNIGC :
            mod->BSIM4v7lnigc = value->rValue;
            mod->BSIM4v7lnigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNIGBINV :
            mod->BSIM4v7lnigbinv = value->rValue;
            mod->BSIM4v7lnigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNIGBACC :
            mod->BSIM4v7lnigbacc = value->rValue;
            mod->BSIM4v7lnigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNTOX :
            mod->BSIM4v7lntox = value->rValue;
            mod->BSIM4v7lntoxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LEIGBINV :
            mod->BSIM4v7leigbinv = value->rValue;
            mod->BSIM4v7leigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LPIGCD :
            mod->BSIM4v7lpigcd = value->rValue;
            mod->BSIM4v7lpigcdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LPOXEDGE :
            mod->BSIM4v7lpoxedge = value->rValue;
            mod->BSIM4v7lpoxedgeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LXRCRG1 :
            mod->BSIM4v7lxrcrg1 = value->rValue;
            mod->BSIM4v7lxrcrg1Given = TRUE;
            break;
        case  BSIM4v7_MOD_LXRCRG2 :
            mod->BSIM4v7lxrcrg2 = value->rValue;
            mod->BSIM4v7lxrcrg2Given = TRUE;
            break;
        case  BSIM4v7_MOD_LLAMBDA :
            mod->BSIM4v7llambda = value->rValue;
            mod->BSIM4v7llambdaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LVTL :
            mod->BSIM4v7lvtl = value->rValue;
            mod->BSIM4v7lvtlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LXN:
            mod->BSIM4v7lxn = value->rValue;
            mod->BSIM4v7lxnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LVFBSDOFF:
            mod->BSIM4v7lvfbsdoff = value->rValue;
            mod->BSIM4v7lvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LTVFBSDOFF:
            mod->BSIM4v7ltvfbsdoff = value->rValue;
            mod->BSIM4v7ltvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LEU :
            mod->BSIM4v7leu = value->rValue;
            mod->BSIM4v7leuGiven = TRUE;
            break;
                case  BSIM4v7_MOD_LUCS :
            mod->BSIM4v7lucs = value->rValue;
            mod->BSIM4v7lucsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LVFB :
            mod->BSIM4v7lvfb = value->rValue;
            mod->BSIM4v7lvfbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCGSL :
            mod->BSIM4v7lcgsl = value->rValue;
            mod->BSIM4v7lcgslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCGDL :
            mod->BSIM4v7lcgdl = value->rValue;
            mod->BSIM4v7lcgdlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCKAPPAS :
            mod->BSIM4v7lckappas = value->rValue;
            mod->BSIM4v7lckappasGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCKAPPAD :
            mod->BSIM4v7lckappad = value->rValue;
            mod->BSIM4v7lckappadGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCF :
            mod->BSIM4v7lcf = value->rValue;
            mod->BSIM4v7lcfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCLC :
            mod->BSIM4v7lclc = value->rValue;
            mod->BSIM4v7lclcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LCLE :
            mod->BSIM4v7lcle = value->rValue;
            mod->BSIM4v7lcleGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LVFBCV :
            mod->BSIM4v7lvfbcv = value->rValue;
            mod->BSIM4v7lvfbcvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LACDE :
            mod->BSIM4v7lacde = value->rValue;
            mod->BSIM4v7lacdeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LMOIN :
            mod->BSIM4v7lmoin = value->rValue;
            mod->BSIM4v7lmoinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LNOFF :
            mod->BSIM4v7lnoff = value->rValue;
            mod->BSIM4v7lnoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LVOFFCV :
            mod->BSIM4v7lvoffcv = value->rValue;
            mod->BSIM4v7lvoffcvGiven = TRUE;
            break;

        /* Width dependence */
        case  BSIM4v7_MOD_WCDSC :
            mod->BSIM4v7wcdsc = value->rValue;
            mod->BSIM4v7wcdscGiven = TRUE;
            break;
       
       
         case  BSIM4v7_MOD_WCDSCB :
            mod->BSIM4v7wcdscb = value->rValue;
            mod->BSIM4v7wcdscbGiven = TRUE;
            break;
         case  BSIM4v7_MOD_WCDSCD :
            mod->BSIM4v7wcdscd = value->rValue;
            mod->BSIM4v7wcdscdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIT :
            mod->BSIM4v7wcit = value->rValue;
            mod->BSIM4v7wcitGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNFACTOR :
            mod->BSIM4v7wnfactor = value->rValue;
            mod->BSIM4v7wnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_WXJ:
            mod->BSIM4v7wxj = value->rValue;
            mod->BSIM4v7wxjGiven = TRUE;
            break;
        case BSIM4v7_MOD_WVSAT:
            mod->BSIM4v7wvsat = value->rValue;
            mod->BSIM4v7wvsatGiven = TRUE;
            break;


        case BSIM4v7_MOD_WA0:
            mod->BSIM4v7wa0 = value->rValue;
            mod->BSIM4v7wa0Given = TRUE;
            break;
        case BSIM4v7_MOD_WAGS:
            mod->BSIM4v7wags = value->rValue;
            mod->BSIM4v7wagsGiven = TRUE;
            break;
        case BSIM4v7_MOD_WA1:
            mod->BSIM4v7wa1 = value->rValue;
            mod->BSIM4v7wa1Given = TRUE;
            break;
        case BSIM4v7_MOD_WA2:
            mod->BSIM4v7wa2 = value->rValue;
            mod->BSIM4v7wa2Given = TRUE;
            break;
        case BSIM4v7_MOD_WAT:
            mod->BSIM4v7wat = value->rValue;
            mod->BSIM4v7watGiven = TRUE;
            break;
        case BSIM4v7_MOD_WKETA:
            mod->BSIM4v7wketa = value->rValue;
            mod->BSIM4v7wketaGiven = TRUE;
            break;    
        case BSIM4v7_MOD_WNSUB:
            mod->BSIM4v7wnsub = value->rValue;
            mod->BSIM4v7wnsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_WNDEP:
            mod->BSIM4v7wndep = value->rValue;
            mod->BSIM4v7wndepGiven = TRUE;
            if (mod->BSIM4v7wndep > 1.0e20)
                mod->BSIM4v7wndep *= 1.0e-6;
            break;
        case BSIM4v7_MOD_WNSD:
            mod->BSIM4v7wnsd = value->rValue;
            mod->BSIM4v7wnsdGiven = TRUE;
            if (mod->BSIM4v7wnsd > 1.0e23)
                mod->BSIM4v7wnsd *= 1.0e-6;
            break;
        case BSIM4v7_MOD_WNGATE:
            mod->BSIM4v7wngate = value->rValue;
            mod->BSIM4v7wngateGiven = TRUE;
            if (mod->BSIM4v7wngate > 1.0e23)
                mod->BSIM4v7wngate *= 1.0e-6;
            break;
        case BSIM4v7_MOD_WGAMMA1:
            mod->BSIM4v7wgamma1 = value->rValue;
            mod->BSIM4v7wgamma1Given = TRUE;
            break;
        case BSIM4v7_MOD_WGAMMA2:
            mod->BSIM4v7wgamma2 = value->rValue;
            mod->BSIM4v7wgamma2Given = TRUE;
            break;
        case BSIM4v7_MOD_WVBX:
            mod->BSIM4v7wvbx = value->rValue;
            mod->BSIM4v7wvbxGiven = TRUE;
            break;
        case BSIM4v7_MOD_WVBM:
            mod->BSIM4v7wvbm = value->rValue;
            mod->BSIM4v7wvbmGiven = TRUE;
            break;
        case BSIM4v7_MOD_WXT:
            mod->BSIM4v7wxt = value->rValue;
            mod->BSIM4v7wxtGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WK1:
            mod->BSIM4v7wk1 = value->rValue;
            mod->BSIM4v7wk1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WKT1:
            mod->BSIM4v7wkt1 = value->rValue;
            mod->BSIM4v7wkt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WKT1L:
            mod->BSIM4v7wkt1l = value->rValue;
            mod->BSIM4v7wkt1lGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WKT2:
            mod->BSIM4v7wkt2 = value->rValue;
            mod->BSIM4v7wkt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_WK2:
            mod->BSIM4v7wk2 = value->rValue;
            mod->BSIM4v7wk2Given = TRUE;
            break;
        case  BSIM4v7_MOD_WK3:
            mod->BSIM4v7wk3 = value->rValue;
            mod->BSIM4v7wk3Given = TRUE;
            break;
        case  BSIM4v7_MOD_WK3B:
            mod->BSIM4v7wk3b = value->rValue;
            mod->BSIM4v7wk3bGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WLPE0:
            mod->BSIM4v7wlpe0 = value->rValue;
            mod->BSIM4v7wlpe0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WLPEB:
            mod->BSIM4v7wlpeb = value->rValue;
            mod->BSIM4v7wlpebGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP0:
            mod->BSIM4v7wdvtp0 = value->rValue;
            mod->BSIM4v7wdvtp0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP1:
            mod->BSIM4v7wdvtp1 = value->rValue;
            mod->BSIM4v7wdvtp1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP2:     /* New DIBL/Rout */
            mod->BSIM4v7wdvtp2 = value->rValue;
            mod->BSIM4v7wdvtp2Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP3:
            mod->BSIM4v7wdvtp3 = value->rValue;
            mod->BSIM4v7wdvtp3Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP4:
            mod->BSIM4v7wdvtp4 = value->rValue;
            mod->BSIM4v7wdvtp4Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVTP5:
            mod->BSIM4v7wdvtp5 = value->rValue;
            mod->BSIM4v7wdvtp5Given = TRUE;
            break;        
        case  BSIM4v7_MOD_WW0:
            mod->BSIM4v7ww0 = value->rValue;
            mod->BSIM4v7ww0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT0:               
            mod->BSIM4v7wdvt0 = value->rValue;
            mod->BSIM4v7wdvt0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT1:             
            mod->BSIM4v7wdvt1 = value->rValue;
            mod->BSIM4v7wdvt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT2:             
            mod->BSIM4v7wdvt2 = value->rValue;
            mod->BSIM4v7wdvt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT0W:               
            mod->BSIM4v7wdvt0w = value->rValue;
            mod->BSIM4v7wdvt0wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT1W:             
            mod->BSIM4v7wdvt1w = value->rValue;
            mod->BSIM4v7wdvt1wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDVT2W:             
            mod->BSIM4v7wdvt2w = value->rValue;
            mod->BSIM4v7wdvt2wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDROUT:             
            mod->BSIM4v7wdrout = value->rValue;
            mod->BSIM4v7wdroutGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDSUB:             
            mod->BSIM4v7wdsub = value->rValue;
            mod->BSIM4v7wdsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_WVTH0:
            mod->BSIM4v7wvth0 = value->rValue;
            mod->BSIM4v7wvth0Given = TRUE;
            break;
        case BSIM4v7_MOD_WUA:
            mod->BSIM4v7wua = value->rValue;
            mod->BSIM4v7wuaGiven = TRUE;
            break;
        case BSIM4v7_MOD_WUA1:
            mod->BSIM4v7wua1 = value->rValue;
            mod->BSIM4v7wua1Given = TRUE;
            break;
        case BSIM4v7_MOD_WUB:
            mod->BSIM4v7wub = value->rValue;
            mod->BSIM4v7wubGiven = TRUE;
            break;
        case BSIM4v7_MOD_WUB1:
            mod->BSIM4v7wub1 = value->rValue;
            mod->BSIM4v7wub1Given = TRUE;
            break;
        case BSIM4v7_MOD_WUC:
            mod->BSIM4v7wuc = value->rValue;
            mod->BSIM4v7wucGiven = TRUE;
            break;
        case BSIM4v7_MOD_WUC1:
            mod->BSIM4v7wuc1 = value->rValue;
            mod->BSIM4v7wuc1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WU0 :
            mod->BSIM4v7wu0 = value->rValue;
            mod->BSIM4v7wu0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WUTE :
            mod->BSIM4v7wute = value->rValue;
            mod->BSIM4v7wuteGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WUCSTE :
            mod->BSIM4v7wucste = value->rValue;
            mod->BSIM4v7wucsteGiven = TRUE;
            break;
        case BSIM4v7_MOD_WVOFF:
            mod->BSIM4v7wvoff = value->rValue;
            mod->BSIM4v7wvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_WTVOFF:
            mod->BSIM4v7wtvoff = value->rValue;
            mod->BSIM4v7wtvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_WTNFACTOR:           /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7wtnfactor = value->rValue;
            mod->BSIM4v7wtnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_WTETA0:                /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7wteta0 = value->rValue;
            mod->BSIM4v7wteta0Given = TRUE;
            break;
        case BSIM4v7_MOD_WTVOFFCV:        /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7wtvoffcv = value->rValue;
            mod->BSIM4v7wtvoffcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_WMINV:
            mod->BSIM4v7wminv = value->rValue;
            mod->BSIM4v7wminvGiven = TRUE;
            break;
        case BSIM4v7_MOD_WMINVCV:
            mod->BSIM4v7wminvcv = value->rValue;
            mod->BSIM4v7wminvcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_WFPROUT:
            mod->BSIM4v7wfprout = value->rValue;
            mod->BSIM4v7wfproutGiven = TRUE;
            break;
        case BSIM4v7_MOD_WPDITS:
            mod->BSIM4v7wpdits = value->rValue;
            mod->BSIM4v7wpditsGiven = TRUE;
            break;
        case BSIM4v7_MOD_WPDITSD:
            mod->BSIM4v7wpditsd = value->rValue;
            mod->BSIM4v7wpditsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDELTA :
            mod->BSIM4v7wdelta = value->rValue;
            mod->BSIM4v7wdeltaGiven = TRUE;
            break;
        case BSIM4v7_MOD_WRDSW:
            mod->BSIM4v7wrdsw = value->rValue;
            mod->BSIM4v7wrdswGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_WRDW:
            mod->BSIM4v7wrdw = value->rValue;
            mod->BSIM4v7wrdwGiven = TRUE;
            break;
        case BSIM4v7_MOD_WRSW:
            mod->BSIM4v7wrsw = value->rValue;
            mod->BSIM4v7wrswGiven = TRUE;
            break;
        case BSIM4v7_MOD_WPRWB:
            mod->BSIM4v7wprwb = value->rValue;
            mod->BSIM4v7wprwbGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_WPRWG:
            mod->BSIM4v7wprwg = value->rValue;
            mod->BSIM4v7wprwgGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_WPRT:
            mod->BSIM4v7wprt = value->rValue;
            mod->BSIM4v7wprtGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_WETA0:
            mod->BSIM4v7weta0 = value->rValue;
            mod->BSIM4v7weta0Given = TRUE;
            break;                 
        case BSIM4v7_MOD_WETAB:
            mod->BSIM4v7wetab = value->rValue;
            mod->BSIM4v7wetabGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_WPCLM:
            mod->BSIM4v7wpclm = value->rValue;
            mod->BSIM4v7wpclmGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_WPDIBL1:
            mod->BSIM4v7wpdibl1 = value->rValue;
            mod->BSIM4v7wpdibl1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_WPDIBL2:
            mod->BSIM4v7wpdibl2 = value->rValue;
            mod->BSIM4v7wpdibl2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_WPDIBLB:
            mod->BSIM4v7wpdiblb = value->rValue;
            mod->BSIM4v7wpdiblbGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_WPSCBE1:
            mod->BSIM4v7wpscbe1 = value->rValue;
            mod->BSIM4v7wpscbe1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_WPSCBE2:
            mod->BSIM4v7wpscbe2 = value->rValue;
            mod->BSIM4v7wpscbe2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_WPVAG:
            mod->BSIM4v7wpvag = value->rValue;
            mod->BSIM4v7wpvagGiven = TRUE;
            break;                 
        case  BSIM4v7_MOD_WWR :
            mod->BSIM4v7wwr = value->rValue;
            mod->BSIM4v7wwrGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDWG :
            mod->BSIM4v7wdwg = value->rValue;
            mod->BSIM4v7wdwgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WDWB :
            mod->BSIM4v7wdwb = value->rValue;
            mod->BSIM4v7wdwbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WB0 :
            mod->BSIM4v7wb0 = value->rValue;
            mod->BSIM4v7wb0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WB1 :
            mod->BSIM4v7wb1 = value->rValue;
            mod->BSIM4v7wb1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WALPHA0 :
            mod->BSIM4v7walpha0 = value->rValue;
            mod->BSIM4v7walpha0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WALPHA1 :
            mod->BSIM4v7walpha1 = value->rValue;
            mod->BSIM4v7walpha1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WBETA0 :
            mod->BSIM4v7wbeta0 = value->rValue;
            mod->BSIM4v7wbeta0Given = TRUE;
            break;
        case  BSIM4v7_MOD_WPHIN :
            mod->BSIM4v7wphin = value->rValue;
            mod->BSIM4v7wphinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAGIDL :
            mod->BSIM4v7wagidl = value->rValue;
            mod->BSIM4v7wagidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBGIDL :
            mod->BSIM4v7wbgidl = value->rValue;
            mod->BSIM4v7wbgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCGIDL :
            mod->BSIM4v7wcgidl = value->rValue;
            mod->BSIM4v7wcgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WEGIDL :
            mod->BSIM4v7wegidl = value->rValue;
            mod->BSIM4v7wegidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WFGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wfgidl = value->rValue;
            mod->BSIM4v7wfgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WKGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wkgidl = value->rValue;
            mod->BSIM4v7wkgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WRGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wrgidl = value->rValue;
            mod->BSIM4v7wrgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAGISL :
            mod->BSIM4v7wagisl = value->rValue;
            mod->BSIM4v7wagislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBGISL :
            mod->BSIM4v7wbgisl = value->rValue;
            mod->BSIM4v7wbgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCGISL :
            mod->BSIM4v7wcgisl = value->rValue;
            mod->BSIM4v7wcgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WEGISL :
            mod->BSIM4v7wegisl = value->rValue;
            mod->BSIM4v7wegislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WFGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wfgisl = value->rValue;
            mod->BSIM4v7wfgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WKGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wkgisl = value->rValue;
            mod->BSIM4v7wkgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WRGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7wrgisl = value->rValue;
            mod->BSIM4v7wrgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGC :
            mod->BSIM4v7waigc = value->rValue;
            mod->BSIM4v7waigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGC :
            mod->BSIM4v7wbigc = value->rValue;
            mod->BSIM4v7wbigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGC :
            mod->BSIM4v7wcigc = value->rValue;
            mod->BSIM4v7wcigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGSD :
            mod->BSIM4v7waigsd = value->rValue;
            mod->BSIM4v7waigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGSD :
            mod->BSIM4v7wbigsd = value->rValue;
            mod->BSIM4v7wbigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGSD :
            mod->BSIM4v7wcigsd = value->rValue;
            mod->BSIM4v7wcigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGS :
            mod->BSIM4v7waigs = value->rValue;
            mod->BSIM4v7waigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGS :
            mod->BSIM4v7wbigs = value->rValue;
            mod->BSIM4v7wbigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGS :
            mod->BSIM4v7wcigs = value->rValue;
            mod->BSIM4v7wcigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGD :
            mod->BSIM4v7waigd = value->rValue;
            mod->BSIM4v7waigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGD :
            mod->BSIM4v7wbigd = value->rValue;
            mod->BSIM4v7wbigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGD :
            mod->BSIM4v7wcigd = value->rValue;
            mod->BSIM4v7wcigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGBACC :
            mod->BSIM4v7waigbacc = value->rValue;
            mod->BSIM4v7waigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGBACC :
            mod->BSIM4v7wbigbacc = value->rValue;
            mod->BSIM4v7wbigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGBACC :
            mod->BSIM4v7wcigbacc = value->rValue;
            mod->BSIM4v7wcigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WAIGBINV :
            mod->BSIM4v7waigbinv = value->rValue;
            mod->BSIM4v7waigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WBIGBINV :
            mod->BSIM4v7wbigbinv = value->rValue;
            mod->BSIM4v7wbigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCIGBINV :
            mod->BSIM4v7wcigbinv = value->rValue;
            mod->BSIM4v7wcigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNIGC :
            mod->BSIM4v7wnigc = value->rValue;
            mod->BSIM4v7wnigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNIGBINV :
            mod->BSIM4v7wnigbinv = value->rValue;
            mod->BSIM4v7wnigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNIGBACC :
            mod->BSIM4v7wnigbacc = value->rValue;
            mod->BSIM4v7wnigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNTOX :
            mod->BSIM4v7wntox = value->rValue;
            mod->BSIM4v7wntoxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WEIGBINV :
            mod->BSIM4v7weigbinv = value->rValue;
            mod->BSIM4v7weigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WPIGCD :
            mod->BSIM4v7wpigcd = value->rValue;
            mod->BSIM4v7wpigcdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WPOXEDGE :
            mod->BSIM4v7wpoxedge = value->rValue;
            mod->BSIM4v7wpoxedgeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WXRCRG1 :
            mod->BSIM4v7wxrcrg1 = value->rValue;
            mod->BSIM4v7wxrcrg1Given = TRUE;
            break;
        case  BSIM4v7_MOD_WXRCRG2 :
            mod->BSIM4v7wxrcrg2 = value->rValue;
            mod->BSIM4v7wxrcrg2Given = TRUE;
            break;
        case  BSIM4v7_MOD_WLAMBDA :
            mod->BSIM4v7wlambda = value->rValue;
            mod->BSIM4v7wlambdaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WVTL :
            mod->BSIM4v7wvtl = value->rValue;
            mod->BSIM4v7wvtlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WXN:
            mod->BSIM4v7wxn = value->rValue;
            mod->BSIM4v7wxnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WVFBSDOFF:
            mod->BSIM4v7wvfbsdoff = value->rValue;
            mod->BSIM4v7wvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WTVFBSDOFF:
            mod->BSIM4v7wtvfbsdoff = value->rValue;
            mod->BSIM4v7wtvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WEU :
            mod->BSIM4v7weu = value->rValue;
            mod->BSIM4v7weuGiven = TRUE;
            break;
                 case  BSIM4v7_MOD_WUCS :
            mod->BSIM4v7wucs = value->rValue;
            mod->BSIM4v7wucsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WVFB :
            mod->BSIM4v7wvfb = value->rValue;
            mod->BSIM4v7wvfbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCGSL :
            mod->BSIM4v7wcgsl = value->rValue;
            mod->BSIM4v7wcgslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCGDL :
            mod->BSIM4v7wcgdl = value->rValue;
            mod->BSIM4v7wcgdlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCKAPPAS :
            mod->BSIM4v7wckappas = value->rValue;
            mod->BSIM4v7wckappasGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCKAPPAD :
            mod->BSIM4v7wckappad = value->rValue;
            mod->BSIM4v7wckappadGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCF :
            mod->BSIM4v7wcf = value->rValue;
            mod->BSIM4v7wcfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCLC :
            mod->BSIM4v7wclc = value->rValue;
            mod->BSIM4v7wclcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WCLE :
            mod->BSIM4v7wcle = value->rValue;
            mod->BSIM4v7wcleGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WVFBCV :
            mod->BSIM4v7wvfbcv = value->rValue;
            mod->BSIM4v7wvfbcvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WACDE :
            mod->BSIM4v7wacde = value->rValue;
            mod->BSIM4v7wacdeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WMOIN :
            mod->BSIM4v7wmoin = value->rValue;
            mod->BSIM4v7wmoinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WNOFF :
            mod->BSIM4v7wnoff = value->rValue;
            mod->BSIM4v7wnoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WVOFFCV :
            mod->BSIM4v7wvoffcv = value->rValue;
            mod->BSIM4v7wvoffcvGiven = TRUE;
            break;

        /* Cross-term dependence */
        case  BSIM4v7_MOD_PCDSC :
            mod->BSIM4v7pcdsc = value->rValue;
            mod->BSIM4v7pcdscGiven = TRUE;
            break;


        case  BSIM4v7_MOD_PCDSCB :
            mod->BSIM4v7pcdscb = value->rValue;
            mod->BSIM4v7pcdscbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCDSCD :
            mod->BSIM4v7pcdscd = value->rValue;
            mod->BSIM4v7pcdscdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIT :
            mod->BSIM4v7pcit = value->rValue;
            mod->BSIM4v7pcitGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNFACTOR :
            mod->BSIM4v7pnfactor = value->rValue;
            mod->BSIM4v7pnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_PXJ:
            mod->BSIM4v7pxj = value->rValue;
            mod->BSIM4v7pxjGiven = TRUE;
            break;
        case BSIM4v7_MOD_PVSAT:
            mod->BSIM4v7pvsat = value->rValue;
            mod->BSIM4v7pvsatGiven = TRUE;
            break;


        case BSIM4v7_MOD_PA0:
            mod->BSIM4v7pa0 = value->rValue;
            mod->BSIM4v7pa0Given = TRUE;
            break;
        case BSIM4v7_MOD_PAGS:
            mod->BSIM4v7pags = value->rValue;
            mod->BSIM4v7pagsGiven = TRUE;
            break;
        case BSIM4v7_MOD_PA1:
            mod->BSIM4v7pa1 = value->rValue;
            mod->BSIM4v7pa1Given = TRUE;
            break;
        case BSIM4v7_MOD_PA2:
            mod->BSIM4v7pa2 = value->rValue;
            mod->BSIM4v7pa2Given = TRUE;
            break;
        case BSIM4v7_MOD_PAT:
            mod->BSIM4v7pat = value->rValue;
            mod->BSIM4v7patGiven = TRUE;
            break;
        case BSIM4v7_MOD_PKETA:
            mod->BSIM4v7pketa = value->rValue;
            mod->BSIM4v7pketaGiven = TRUE;
            break;    
        case BSIM4v7_MOD_PNSUB:
            mod->BSIM4v7pnsub = value->rValue;
            mod->BSIM4v7pnsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_PNDEP:
            mod->BSIM4v7pndep = value->rValue;
            mod->BSIM4v7pndepGiven = TRUE;
            if (mod->BSIM4v7pndep > 1.0e20)
                mod->BSIM4v7pndep *= 1.0e-6;
            break;
        case BSIM4v7_MOD_PNSD:
            mod->BSIM4v7pnsd = value->rValue;
            mod->BSIM4v7pnsdGiven = TRUE;
            if (mod->BSIM4v7pnsd > 1.0e23)
                mod->BSIM4v7pnsd *= 1.0e-6;
            break;
        case BSIM4v7_MOD_PNGATE:
            mod->BSIM4v7pngate = value->rValue;
            mod->BSIM4v7pngateGiven = TRUE;
            if (mod->BSIM4v7pngate > 1.0e23)
                mod->BSIM4v7pngate *= 1.0e-6;
            break;
        case BSIM4v7_MOD_PGAMMA1:
            mod->BSIM4v7pgamma1 = value->rValue;
            mod->BSIM4v7pgamma1Given = TRUE;
            break;
        case BSIM4v7_MOD_PGAMMA2:
            mod->BSIM4v7pgamma2 = value->rValue;
            mod->BSIM4v7pgamma2Given = TRUE;
            break;
        case BSIM4v7_MOD_PVBX:
            mod->BSIM4v7pvbx = value->rValue;
            mod->BSIM4v7pvbxGiven = TRUE;
            break;
        case BSIM4v7_MOD_PVBM:
            mod->BSIM4v7pvbm = value->rValue;
            mod->BSIM4v7pvbmGiven = TRUE;
            break;
        case BSIM4v7_MOD_PXT:
            mod->BSIM4v7pxt = value->rValue;
            mod->BSIM4v7pxtGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PK1:
            mod->BSIM4v7pk1 = value->rValue;
            mod->BSIM4v7pk1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PKT1:
            mod->BSIM4v7pkt1 = value->rValue;
            mod->BSIM4v7pkt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PKT1L:
            mod->BSIM4v7pkt1l = value->rValue;
            mod->BSIM4v7pkt1lGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PKT2:
            mod->BSIM4v7pkt2 = value->rValue;
            mod->BSIM4v7pkt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_PK2:
            mod->BSIM4v7pk2 = value->rValue;
            mod->BSIM4v7pk2Given = TRUE;
            break;
        case  BSIM4v7_MOD_PK3:
            mod->BSIM4v7pk3 = value->rValue;
            mod->BSIM4v7pk3Given = TRUE;
            break;
        case  BSIM4v7_MOD_PK3B:
            mod->BSIM4v7pk3b = value->rValue;
            mod->BSIM4v7pk3bGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PLPE0:
            mod->BSIM4v7plpe0 = value->rValue;
            mod->BSIM4v7plpe0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PLPEB:
            mod->BSIM4v7plpeb = value->rValue;
            mod->BSIM4v7plpebGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP0:
            mod->BSIM4v7pdvtp0 = value->rValue;
            mod->BSIM4v7pdvtp0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP1:
            mod->BSIM4v7pdvtp1 = value->rValue;
            mod->BSIM4v7pdvtp1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP2:     /* New DIBL/Rout */
            mod->BSIM4v7pdvtp2 = value->rValue;
            mod->BSIM4v7pdvtp2Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP3:
            mod->BSIM4v7pdvtp3 = value->rValue;
            mod->BSIM4v7pdvtp3Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP4:
            mod->BSIM4v7pdvtp4 = value->rValue;
            mod->BSIM4v7pdvtp4Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVTP5:
            mod->BSIM4v7pdvtp5 = value->rValue;
            mod->BSIM4v7pdvtp5Given = TRUE;
            break;
        case  BSIM4v7_MOD_PW0:
            mod->BSIM4v7pw0 = value->rValue;
            mod->BSIM4v7pw0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT0:               
            mod->BSIM4v7pdvt0 = value->rValue;
            mod->BSIM4v7pdvt0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT1:             
            mod->BSIM4v7pdvt1 = value->rValue;
            mod->BSIM4v7pdvt1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT2:             
            mod->BSIM4v7pdvt2 = value->rValue;
            mod->BSIM4v7pdvt2Given = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT0W:               
            mod->BSIM4v7pdvt0w = value->rValue;
            mod->BSIM4v7pdvt0wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT1W:             
            mod->BSIM4v7pdvt1w = value->rValue;
            mod->BSIM4v7pdvt1wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDVT2W:             
            mod->BSIM4v7pdvt2w = value->rValue;
            mod->BSIM4v7pdvt2wGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDROUT:             
            mod->BSIM4v7pdrout = value->rValue;
            mod->BSIM4v7pdroutGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDSUB:             
            mod->BSIM4v7pdsub = value->rValue;
            mod->BSIM4v7pdsubGiven = TRUE;
            break;
        case BSIM4v7_MOD_PVTH0:
            mod->BSIM4v7pvth0 = value->rValue;
            mod->BSIM4v7pvth0Given = TRUE;
            break;
        case BSIM4v7_MOD_PUA:
            mod->BSIM4v7pua = value->rValue;
            mod->BSIM4v7puaGiven = TRUE;
            break;
        case BSIM4v7_MOD_PUA1:
            mod->BSIM4v7pua1 = value->rValue;
            mod->BSIM4v7pua1Given = TRUE;
            break;
        case BSIM4v7_MOD_PUB:
            mod->BSIM4v7pub = value->rValue;
            mod->BSIM4v7pubGiven = TRUE;
            break;
        case BSIM4v7_MOD_PUB1:
            mod->BSIM4v7pub1 = value->rValue;
            mod->BSIM4v7pub1Given = TRUE;
            break;
        case BSIM4v7_MOD_PUC:
            mod->BSIM4v7puc = value->rValue;
            mod->BSIM4v7pucGiven = TRUE;
            break;
        case BSIM4v7_MOD_PUC1:
            mod->BSIM4v7puc1 = value->rValue;
            mod->BSIM4v7puc1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PU0 :
            mod->BSIM4v7pu0 = value->rValue;
            mod->BSIM4v7pu0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PUTE :
            mod->BSIM4v7pute = value->rValue;
            mod->BSIM4v7puteGiven = TRUE;
            break;
                 case  BSIM4v7_MOD_PUCSTE :
            mod->BSIM4v7pucste = value->rValue;
            mod->BSIM4v7pucsteGiven = TRUE;
            break;
        case BSIM4v7_MOD_PVOFF:
            mod->BSIM4v7pvoff = value->rValue;
            mod->BSIM4v7pvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_PTVOFF:
            mod->BSIM4v7ptvoff = value->rValue;
            mod->BSIM4v7ptvoffGiven = TRUE;
            break;
        case BSIM4v7_MOD_PTNFACTOR:           /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7ptnfactor = value->rValue;
            mod->BSIM4v7ptnfactorGiven = TRUE;
            break;
        case BSIM4v7_MOD_PTETA0:                /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7pteta0 = value->rValue;
            mod->BSIM4v7pteta0Given = TRUE;
            break;
        case BSIM4v7_MOD_PTVOFFCV:        /* v4.7 temp dep of leakage current  */
            mod->BSIM4v7ptvoffcv = value->rValue;
            mod->BSIM4v7ptvoffcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_PMINV:
            mod->BSIM4v7pminv = value->rValue;
            mod->BSIM4v7pminvGiven = TRUE;
            break;
        case BSIM4v7_MOD_PMINVCV:
            mod->BSIM4v7pminvcv = value->rValue;
            mod->BSIM4v7pminvcvGiven = TRUE;
            break;
        case BSIM4v7_MOD_PFPROUT:
            mod->BSIM4v7pfprout = value->rValue;
            mod->BSIM4v7pfproutGiven = TRUE;
            break;
        case BSIM4v7_MOD_PPDITS:
            mod->BSIM4v7ppdits = value->rValue;
            mod->BSIM4v7ppditsGiven = TRUE;
            break;
        case BSIM4v7_MOD_PPDITSD:
            mod->BSIM4v7ppditsd = value->rValue;
            mod->BSIM4v7ppditsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDELTA :
            mod->BSIM4v7pdelta = value->rValue;
            mod->BSIM4v7pdeltaGiven = TRUE;
            break;
        case BSIM4v7_MOD_PRDSW:
            mod->BSIM4v7prdsw = value->rValue;
            mod->BSIM4v7prdswGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PRDW:
            mod->BSIM4v7prdw = value->rValue;
            mod->BSIM4v7prdwGiven = TRUE;
            break;
        case BSIM4v7_MOD_PRSW:
            mod->BSIM4v7prsw = value->rValue;
            mod->BSIM4v7prswGiven = TRUE;
            break;
        case BSIM4v7_MOD_PPRWB:
            mod->BSIM4v7pprwb = value->rValue;
            mod->BSIM4v7pprwbGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PPRWG:
            mod->BSIM4v7pprwg = value->rValue;
            mod->BSIM4v7pprwgGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PPRT:
            mod->BSIM4v7pprt = value->rValue;
            mod->BSIM4v7pprtGiven = TRUE;
            break;                     
        case BSIM4v7_MOD_PETA0:
            mod->BSIM4v7peta0 = value->rValue;
            mod->BSIM4v7peta0Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PETAB:
            mod->BSIM4v7petab = value->rValue;
            mod->BSIM4v7petabGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PPCLM:
            mod->BSIM4v7ppclm = value->rValue;
            mod->BSIM4v7ppclmGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PPDIBL1:
            mod->BSIM4v7ppdibl1 = value->rValue;
            mod->BSIM4v7ppdibl1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PPDIBL2:
            mod->BSIM4v7ppdibl2 = value->rValue;
            mod->BSIM4v7ppdibl2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PPDIBLB:
            mod->BSIM4v7ppdiblb = value->rValue;
            mod->BSIM4v7ppdiblbGiven = TRUE;
            break;                 
        case BSIM4v7_MOD_PPSCBE1:
            mod->BSIM4v7ppscbe1 = value->rValue;
            mod->BSIM4v7ppscbe1Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PPSCBE2:
            mod->BSIM4v7ppscbe2 = value->rValue;
            mod->BSIM4v7ppscbe2Given = TRUE;
            break;                 
        case BSIM4v7_MOD_PPVAG:
            mod->BSIM4v7ppvag = value->rValue;
            mod->BSIM4v7ppvagGiven = TRUE;
            break;                 
        case  BSIM4v7_MOD_PWR :
            mod->BSIM4v7pwr = value->rValue;
            mod->BSIM4v7pwrGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDWG :
            mod->BSIM4v7pdwg = value->rValue;
            mod->BSIM4v7pdwgGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PDWB :
            mod->BSIM4v7pdwb = value->rValue;
            mod->BSIM4v7pdwbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PB0 :
            mod->BSIM4v7pb0 = value->rValue;
            mod->BSIM4v7pb0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PB1 :
            mod->BSIM4v7pb1 = value->rValue;
            mod->BSIM4v7pb1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PALPHA0 :
            mod->BSIM4v7palpha0 = value->rValue;
            mod->BSIM4v7palpha0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PALPHA1 :
            mod->BSIM4v7palpha1 = value->rValue;
            mod->BSIM4v7palpha1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PBETA0 :
            mod->BSIM4v7pbeta0 = value->rValue;
            mod->BSIM4v7pbeta0Given = TRUE;
            break;
        case  BSIM4v7_MOD_PPHIN :
            mod->BSIM4v7pphin = value->rValue;
            mod->BSIM4v7pphinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAGIDL :
            mod->BSIM4v7pagidl = value->rValue;
            mod->BSIM4v7pagidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBGIDL :
            mod->BSIM4v7pbgidl = value->rValue;
            mod->BSIM4v7pbgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCGIDL :
            mod->BSIM4v7pcgidl = value->rValue;
            mod->BSIM4v7pcgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PEGIDL :
            mod->BSIM4v7pegidl = value->rValue;
            mod->BSIM4v7pegidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PFGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7pfgidl = value->rValue;
            mod->BSIM4v7pfgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PKGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7pkgidl = value->rValue;
            mod->BSIM4v7pkgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PRGIDL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7prgidl = value->rValue;
            mod->BSIM4v7prgidlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAGISL :
            mod->BSIM4v7pagisl = value->rValue;
            mod->BSIM4v7pagislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBGISL :
            mod->BSIM4v7pbgisl = value->rValue;
            mod->BSIM4v7pbgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCGISL :
            mod->BSIM4v7pcgisl = value->rValue;
            mod->BSIM4v7pcgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PEGISL :
            mod->BSIM4v7pegisl = value->rValue;
            mod->BSIM4v7pegislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PFGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7pfgisl = value->rValue;
            mod->BSIM4v7pfgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PKGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7pkgisl = value->rValue;
            mod->BSIM4v7pkgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PRGISL :                        /* v4.7 New GIDL/GISL */
            mod->BSIM4v7prgisl = value->rValue;
            mod->BSIM4v7prgislGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGC :
            mod->BSIM4v7paigc = value->rValue;
            mod->BSIM4v7paigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGC :
            mod->BSIM4v7pbigc = value->rValue;
            mod->BSIM4v7pbigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGC :
            mod->BSIM4v7pcigc = value->rValue;
            mod->BSIM4v7pcigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGSD :
            mod->BSIM4v7paigsd = value->rValue;
            mod->BSIM4v7paigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGSD :
            mod->BSIM4v7pbigsd = value->rValue;
            mod->BSIM4v7pbigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGSD :
            mod->BSIM4v7pcigsd = value->rValue;
            mod->BSIM4v7pcigsdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGS :
            mod->BSIM4v7paigs = value->rValue;
            mod->BSIM4v7paigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGS :
            mod->BSIM4v7pbigs = value->rValue;
            mod->BSIM4v7pbigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGS :
            mod->BSIM4v7pcigs = value->rValue;
            mod->BSIM4v7pcigsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGD :
            mod->BSIM4v7paigd = value->rValue;
            mod->BSIM4v7paigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGD :
            mod->BSIM4v7pbigd = value->rValue;
            mod->BSIM4v7pbigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGD :
            mod->BSIM4v7pcigd = value->rValue;
            mod->BSIM4v7pcigdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGBACC :
            mod->BSIM4v7paigbacc = value->rValue;
            mod->BSIM4v7paigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGBACC :
            mod->BSIM4v7pbigbacc = value->rValue;
            mod->BSIM4v7pbigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGBACC :
            mod->BSIM4v7pcigbacc = value->rValue;
            mod->BSIM4v7pcigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PAIGBINV :
            mod->BSIM4v7paigbinv = value->rValue;
            mod->BSIM4v7paigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBIGBINV :
            mod->BSIM4v7pbigbinv = value->rValue;
            mod->BSIM4v7pbigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCIGBINV :
            mod->BSIM4v7pcigbinv = value->rValue;
            mod->BSIM4v7pcigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNIGC :
            mod->BSIM4v7pnigc = value->rValue;
            mod->BSIM4v7pnigcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNIGBINV :
            mod->BSIM4v7pnigbinv = value->rValue;
            mod->BSIM4v7pnigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNIGBACC :
            mod->BSIM4v7pnigbacc = value->rValue;
            mod->BSIM4v7pnigbaccGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNTOX :
            mod->BSIM4v7pntox = value->rValue;
            mod->BSIM4v7pntoxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PEIGBINV :
            mod->BSIM4v7peigbinv = value->rValue;
            mod->BSIM4v7peigbinvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PPIGCD :
            mod->BSIM4v7ppigcd = value->rValue;
            mod->BSIM4v7ppigcdGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PPOXEDGE :
            mod->BSIM4v7ppoxedge = value->rValue;
            mod->BSIM4v7ppoxedgeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PXRCRG1 :
            mod->BSIM4v7pxrcrg1 = value->rValue;
            mod->BSIM4v7pxrcrg1Given = TRUE;
            break;
        case  BSIM4v7_MOD_PXRCRG2 :
            mod->BSIM4v7pxrcrg2 = value->rValue;
            mod->BSIM4v7pxrcrg2Given = TRUE;
            break;
        case  BSIM4v7_MOD_PLAMBDA :
            mod->BSIM4v7plambda = value->rValue;
            mod->BSIM4v7plambdaGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PVTL :
            mod->BSIM4v7pvtl = value->rValue;
            mod->BSIM4v7pvtlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PXN:
            mod->BSIM4v7pxn = value->rValue;
            mod->BSIM4v7pxnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PVFBSDOFF:
            mod->BSIM4v7pvfbsdoff = value->rValue;
            mod->BSIM4v7pvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PTVFBSDOFF:
            mod->BSIM4v7ptvfbsdoff = value->rValue;
            mod->BSIM4v7ptvfbsdoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PEU :
            mod->BSIM4v7peu = value->rValue;
            mod->BSIM4v7peuGiven = TRUE;
            break;
                case  BSIM4v7_MOD_PUCS :
            mod->BSIM4v7pucs = value->rValue;
            mod->BSIM4v7pucsGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PVFB :
            mod->BSIM4v7pvfb = value->rValue;
            mod->BSIM4v7pvfbGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCGSL :
            mod->BSIM4v7pcgsl = value->rValue;
            mod->BSIM4v7pcgslGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCGDL :
            mod->BSIM4v7pcgdl = value->rValue;
            mod->BSIM4v7pcgdlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCKAPPAS :
            mod->BSIM4v7pckappas = value->rValue;
            mod->BSIM4v7pckappasGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCKAPPAD :
            mod->BSIM4v7pckappad = value->rValue;
            mod->BSIM4v7pckappadGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCF :
            mod->BSIM4v7pcf = value->rValue;
            mod->BSIM4v7pcfGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCLC :
            mod->BSIM4v7pclc = value->rValue;
            mod->BSIM4v7pclcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PCLE :
            mod->BSIM4v7pcle = value->rValue;
            mod->BSIM4v7pcleGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PVFBCV :
            mod->BSIM4v7pvfbcv = value->rValue;
            mod->BSIM4v7pvfbcvGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PACDE :
            mod->BSIM4v7pacde = value->rValue;
            mod->BSIM4v7pacdeGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PMOIN :
            mod->BSIM4v7pmoin = value->rValue;
            mod->BSIM4v7pmoinGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PNOFF :
            mod->BSIM4v7pnoff = value->rValue;
            mod->BSIM4v7pnoffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PVOFFCV :
            mod->BSIM4v7pvoffcv = value->rValue;
            mod->BSIM4v7pvoffcvGiven = TRUE;
            break;

        case  BSIM4v7_MOD_TNOM :
            mod->BSIM4v7tnom = value->rValue + CONSTCtoK;
            mod->BSIM4v7tnomGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGSO :
            mod->BSIM4v7cgso = value->rValue;
            mod->BSIM4v7cgsoGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGDO :
            mod->BSIM4v7cgdo = value->rValue;
            mod->BSIM4v7cgdoGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CGBO :
            mod->BSIM4v7cgbo = value->rValue;
            mod->BSIM4v7cgboGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XPART :
            mod->BSIM4v7xpart = value->rValue;
            mod->BSIM4v7xpartGiven = TRUE;
            break;
        case  BSIM4v7_MOD_RSH :
            mod->BSIM4v7sheetResistance = value->rValue;
            mod->BSIM4v7sheetResistanceGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSS :
            mod->BSIM4v7SjctSatCurDensity = value->rValue;
            mod->BSIM4v7SjctSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSWS :
            mod->BSIM4v7SjctSidewallSatCurDensity = value->rValue;
            mod->BSIM4v7SjctSidewallSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSWGS :
            mod->BSIM4v7SjctGateSidewallSatCurDensity = value->rValue;
            mod->BSIM4v7SjctGateSidewallSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBS :
            mod->BSIM4v7SbulkJctPotential = value->rValue;
            mod->BSIM4v7SbulkJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJS :
            mod->BSIM4v7SbulkJctBotGradingCoeff = value->rValue;
            mod->BSIM4v7SbulkJctBotGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBSWS :
            mod->BSIM4v7SsidewallJctPotential = value->rValue;
            mod->BSIM4v7SsidewallJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJSWS :
            mod->BSIM4v7SbulkJctSideGradingCoeff = value->rValue;
            mod->BSIM4v7SbulkJctSideGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJS :
            mod->BSIM4v7SunitAreaJctCap = value->rValue;
            mod->BSIM4v7SunitAreaJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJSWS :
            mod->BSIM4v7SunitLengthSidewallJctCap = value->rValue;
            mod->BSIM4v7SunitLengthSidewallJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJS :
            mod->BSIM4v7SjctEmissionCoeff = value->rValue;
            mod->BSIM4v7SjctEmissionCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBSWGS :
            mod->BSIM4v7SGatesidewallJctPotential = value->rValue;
            mod->BSIM4v7SGatesidewallJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJSWGS :
            mod->BSIM4v7SbulkJctGateSideGradingCoeff = value->rValue;
            mod->BSIM4v7SbulkJctGateSideGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJSWGS :
            mod->BSIM4v7SunitLengthGateSidewallJctCap = value->rValue;
            mod->BSIM4v7SunitLengthGateSidewallJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTIS :
            mod->BSIM4v7SjctTempExponent = value->rValue;
            mod->BSIM4v7SjctTempExponentGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSD :
            mod->BSIM4v7DjctSatCurDensity = value->rValue;
            mod->BSIM4v7DjctSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSWD :
            mod->BSIM4v7DjctSidewallSatCurDensity = value->rValue;
            mod->BSIM4v7DjctSidewallSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_JSWGD :
            mod->BSIM4v7DjctGateSidewallSatCurDensity = value->rValue;
            mod->BSIM4v7DjctGateSidewallSatCurDensityGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBD :
            mod->BSIM4v7DbulkJctPotential = value->rValue;
            mod->BSIM4v7DbulkJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJD :
            mod->BSIM4v7DbulkJctBotGradingCoeff = value->rValue;
            mod->BSIM4v7DbulkJctBotGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBSWD :
            mod->BSIM4v7DsidewallJctPotential = value->rValue;
            mod->BSIM4v7DsidewallJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJSWD :
            mod->BSIM4v7DbulkJctSideGradingCoeff = value->rValue;
            mod->BSIM4v7DbulkJctSideGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJD :
            mod->BSIM4v7DunitAreaJctCap = value->rValue;
            mod->BSIM4v7DunitAreaJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJSWD :
            mod->BSIM4v7DunitLengthSidewallJctCap = value->rValue;
            mod->BSIM4v7DunitLengthSidewallJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NJD :
            mod->BSIM4v7DjctEmissionCoeff = value->rValue;
            mod->BSIM4v7DjctEmissionCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_PBSWGD :
            mod->BSIM4v7DGatesidewallJctPotential = value->rValue;
            mod->BSIM4v7DGatesidewallJctPotentialGiven = TRUE;
            break;
        case  BSIM4v7_MOD_MJSWGD :
            mod->BSIM4v7DbulkJctGateSideGradingCoeff = value->rValue;
            mod->BSIM4v7DbulkJctGateSideGradingCoeffGiven = TRUE;
            break;
        case  BSIM4v7_MOD_CJSWGD :
            mod->BSIM4v7DunitLengthGateSidewallJctCap = value->rValue;
            mod->BSIM4v7DunitLengthGateSidewallJctCapGiven = TRUE;
            break;
        case  BSIM4v7_MOD_XTID :
            mod->BSIM4v7DjctTempExponent = value->rValue;
            mod->BSIM4v7DjctTempExponentGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LINT :
            mod->BSIM4v7Lint = value->rValue;
            mod->BSIM4v7LintGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LL :
            mod->BSIM4v7Ll = value->rValue;
            mod->BSIM4v7LlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LLC :
            mod->BSIM4v7Llc = value->rValue;
            mod->BSIM4v7LlcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LLN :
            mod->BSIM4v7Lln = value->rValue;
            mod->BSIM4v7LlnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LW :
            mod->BSIM4v7Lw = value->rValue;
            mod->BSIM4v7LwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LWC :
            mod->BSIM4v7Lwc = value->rValue;
            mod->BSIM4v7LwcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LWN :
            mod->BSIM4v7Lwn = value->rValue;
            mod->BSIM4v7LwnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LWL :
            mod->BSIM4v7Lwl = value->rValue;
            mod->BSIM4v7LwlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LWLC :
            mod->BSIM4v7Lwlc = value->rValue;
            mod->BSIM4v7LwlcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LMIN :
            mod->BSIM4v7Lmin = value->rValue;
            mod->BSIM4v7LminGiven = TRUE;
            break;
        case  BSIM4v7_MOD_LMAX :
            mod->BSIM4v7Lmax = value->rValue;
            mod->BSIM4v7LmaxGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WINT :
            mod->BSIM4v7Wint = value->rValue;
            mod->BSIM4v7WintGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WL :
            mod->BSIM4v7Wl = value->rValue;
            mod->BSIM4v7WlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WLC :
            mod->BSIM4v7Wlc = value->rValue;
            mod->BSIM4v7WlcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WLN :
            mod->BSIM4v7Wln = value->rValue;
            mod->BSIM4v7WlnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WW :
            mod->BSIM4v7Ww = value->rValue;
            mod->BSIM4v7WwGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WWC :
            mod->BSIM4v7Wwc = value->rValue;
            mod->BSIM4v7WwcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WWN :
            mod->BSIM4v7Wwn = value->rValue;
            mod->BSIM4v7WwnGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WWL :
            mod->BSIM4v7Wwl = value->rValue;
            mod->BSIM4v7WwlGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WWLC :
            mod->BSIM4v7Wwlc = value->rValue;
            mod->BSIM4v7WwlcGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WMIN :
            mod->BSIM4v7Wmin = value->rValue;
            mod->BSIM4v7WminGiven = TRUE;
            break;
        case  BSIM4v7_MOD_WMAX :
            mod->BSIM4v7Wmax = value->rValue;
            mod->BSIM4v7WmaxGiven = TRUE;
            break;

        case  BSIM4v7_MOD_NOIA :
            mod->BSIM4v7oxideTrapDensityA = value->rValue;
            mod->BSIM4v7oxideTrapDensityAGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NOIB :
            mod->BSIM4v7oxideTrapDensityB = value->rValue;
            mod->BSIM4v7oxideTrapDensityBGiven = TRUE;
            break;
        case  BSIM4v7_MOD_NOIC :
            mod->BSIM4v7oxideTrapDensityC = value->rValue;
            mod->BSIM4v7oxideTrapDensityCGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EM :
            mod->BSIM4v7em = value->rValue;
            mod->BSIM4v7emGiven = TRUE;
            break;
        case  BSIM4v7_MOD_EF :
            mod->BSIM4v7ef = value->rValue;
            mod->BSIM4v7efGiven = TRUE;
            break;
        case  BSIM4v7_MOD_AF :
            mod->BSIM4v7af = value->rValue;
            mod->BSIM4v7afGiven = TRUE;
            break;
        case  BSIM4v7_MOD_KF :
            mod->BSIM4v7kf = value->rValue;
            mod->BSIM4v7kfGiven = TRUE;
            break;

        case BSIM4v7_MOD_VGS_MAX:
            mod->BSIM4v7vgsMax = value->rValue;
            mod->BSIM4v7vgsMaxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VGD_MAX:
            mod->BSIM4v7vgdMax = value->rValue;
            mod->BSIM4v7vgdMaxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VGB_MAX:
            mod->BSIM4v7vgbMax = value->rValue;
            mod->BSIM4v7vgbMaxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VDS_MAX:
            mod->BSIM4v7vdsMax = value->rValue;
            mod->BSIM4v7vdsMaxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VBS_MAX:
            mod->BSIM4v7vbsMax = value->rValue;
            mod->BSIM4v7vbsMaxGiven = TRUE;
            break;
        case BSIM4v7_MOD_VBD_MAX:
            mod->BSIM4v7vbdMax = value->rValue;
            mod->BSIM4v7vbdMaxGiven = TRUE;
            break;

        case  BSIM4v7_MOD_NMOS  :
            if(value->iValue) {
                mod->BSIM4v7type = 1;
                mod->BSIM4v7typeGiven = TRUE;
            }
            break;
        case  BSIM4v7_MOD_PMOS  :
            if(value->iValue) {
                mod->BSIM4v7type = - 1;
                mod->BSIM4v7typeGiven = TRUE;
            }
            break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}


