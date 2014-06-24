/**** BSIM4.7.0 Released by Darsen Lu 04/08/2011 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4set.c of BSIM4.7.0.
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
#include "ngspice/jobdefs.h"
#include "ngspice/ftedefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "bsim4v7def.h"
#include "ngspice/const.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

#ifdef USE_OMP
#include "ngspice/cpextern.h"
#endif

#define MAX_EXP 5.834617425e14
#define MIN_EXP 1.713908431e-15
#define EXP_THRESHOLD 34.0
#define EPS0 8.85418e-12
#define EPSSI 1.03594e-10
#define PI 3.141592654
#define Charge_q 1.60219e-19

int
BSIM4v7setup(
SMPmatrix *matrix,
GENmodel *inModel,
CKTcircuit *ckt,
int *states)
{
BSIM4v7model *model = (BSIM4v7model*)inModel;
BSIM4v7instance *here;
int error;
CKTnode *tmp;
int    noiseAnalGiven = 0, createNode;  /* Criteria for new node creation */
double Rtot, DMCGeff, DMCIeff, DMDGeff;
JOB   *job;

#ifdef USE_OMP
unsigned int idx, InstCount;
BSIM4v7instance **InstArray;
#endif

    /* Search for a noise analysis request */
    for (job = ((TSKtask *)ft_curckt->ci_curTask)->jobs;job;job = job->JOBnextJob) {
        if(strcmp(job->JOBname,"Noise Analysis")==0) {
            noiseAnalGiven = 1;
            break;
        }
    }

    /*  loop through all the BSIM4v7 device models */
    for( ; model != NULL; model = model->BSIM4v7nextModel )
    {   /* process defaults of model parameters */
        if (!model->BSIM4v7typeGiven)
            model->BSIM4v7type = NMOS;     

        if (!model->BSIM4v7mobModGiven) 
            model->BSIM4v7mobMod = 0;
        else if ((model->BSIM4v7mobMod != 0) && (model->BSIM4v7mobMod != 1)
                 && (model->BSIM4v7mobMod != 2)&& (model->BSIM4v7mobMod != 3))
        {   model->BSIM4v7mobMod = 0;
            printf("Warning: mobMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7binUnitGiven) 
            model->BSIM4v7binUnit = 1;
        if (!model->BSIM4v7paramChkGiven) 
            model->BSIM4v7paramChk = 1;

        if (!model->BSIM4v7dioModGiven)
            model->BSIM4v7dioMod = 1;
        else if ((model->BSIM4v7dioMod != 0) && (model->BSIM4v7dioMod != 1)
            && (model->BSIM4v7dioMod != 2))
        {   model->BSIM4v7dioMod = 1;
            printf("Warning: dioMod has been set to its default value: 1.\n");
        }

        if (!model->BSIM4v7cvchargeModGiven) 
            model->BSIM4v7cvchargeMod = 0;
        if (!model->BSIM4v7capModGiven) 
            model->BSIM4v7capMod = 2;
        else if ((model->BSIM4v7capMod != 0) && (model->BSIM4v7capMod != 1)
            && (model->BSIM4v7capMod != 2))
        {   model->BSIM4v7capMod = 2;
            printf("Warning: capMod has been set to its default value: 2.\n");
        }

        if (!model->BSIM4v7rdsModGiven)
            model->BSIM4v7rdsMod = 0;
        else if ((model->BSIM4v7rdsMod != 0) && (model->BSIM4v7rdsMod != 1))
        {   model->BSIM4v7rdsMod = 0;
            printf("Warning: rdsMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7rbodyModGiven)
            model->BSIM4v7rbodyMod = 0;
        else if ((model->BSIM4v7rbodyMod != 0) && (model->BSIM4v7rbodyMod != 1) && (model->BSIM4v7rbodyMod != 2))
        {   model->BSIM4v7rbodyMod = 0;
            printf("Warning: rbodyMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7rgateModGiven)
            model->BSIM4v7rgateMod = 0;
        else if ((model->BSIM4v7rgateMod != 0) && (model->BSIM4v7rgateMod != 1)
            && (model->BSIM4v7rgateMod != 2) && (model->BSIM4v7rgateMod != 3))
        {   model->BSIM4v7rgateMod = 0;
            printf("Warning: rgateMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7perModGiven)
            model->BSIM4v7perMod = 1;
        else if ((model->BSIM4v7perMod != 0) && (model->BSIM4v7perMod != 1))
        {   model->BSIM4v7perMod = 1;
            printf("Warning: perMod has been set to its default value: 1.\n");
        }

        if (!model->BSIM4v7geoModGiven)
            model->BSIM4v7geoMod = 0;

        if (!model->BSIM4v7rgeoModGiven)
            model->BSIM4v7rgeoMod = 0;
        else if ((model->BSIM4v7rgeoMod != 0) && (model->BSIM4v7rgeoMod != 1))
        {   model->BSIM4v7rgeoMod = 1;
            printf("Warning: rgeoMod has been set to its default value: 1.\n");
        }
        if (!model->BSIM4v7fnoiModGiven) 
            model->BSIM4v7fnoiMod = 1;
        else if ((model->BSIM4v7fnoiMod != 0) && (model->BSIM4v7fnoiMod != 1))
        {   model->BSIM4v7fnoiMod = 1;
            printf("Warning: fnoiMod has been set to its default value: 1.\n");
        }
        if (!model->BSIM4v7tnoiModGiven)
            model->BSIM4v7tnoiMod = 0; /* WDLiu: tnoiMod=1 needs to set internal S/D nodes */
        else if ((model->BSIM4v7tnoiMod != 0) && (model->BSIM4v7tnoiMod != 1) && (model->BSIM4v7tnoiMod != 2))  /* v4.7 */
        {   model->BSIM4v7tnoiMod = 0;
            printf("Warning: tnoiMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7trnqsModGiven)
            model->BSIM4v7trnqsMod = 0; 
        else if ((model->BSIM4v7trnqsMod != 0) && (model->BSIM4v7trnqsMod != 1))
        {   model->BSIM4v7trnqsMod = 0;
            printf("Warning: trnqsMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7acnqsModGiven)
            model->BSIM4v7acnqsMod = 0;
        else if ((model->BSIM4v7acnqsMod != 0) && (model->BSIM4v7acnqsMod != 1))
        {   model->BSIM4v7acnqsMod = 0;
            printf("Warning: acnqsMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7mtrlModGiven)   
            model->BSIM4v7mtrlMod = 0;
        else if((model->BSIM4v7mtrlMod != 0) && (model->BSIM4v7mtrlMod != 1))
        {
            model->BSIM4v7mtrlMod = 0;
            printf("Warning: mtrlMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7mtrlCompatModGiven)   
            model->BSIM4v7mtrlCompatMod = 0;
        else if((model->BSIM4v7mtrlCompatMod != 0) && (model->BSIM4v7mtrlCompatMod != 1))
        {
            model->BSIM4v7mtrlCompatMod = 0;
            printf("Warning: mtrlCompatMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7igcModGiven)
            model->BSIM4v7igcMod = 0;
        else if ((model->BSIM4v7igcMod != 0) && (model->BSIM4v7igcMod != 1)
                  && (model->BSIM4v7igcMod != 2))
        {   model->BSIM4v7igcMod = 0;
            printf("Warning: igcMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7igbModGiven)
            model->BSIM4v7igbMod = 0;
        else if ((model->BSIM4v7igbMod != 0) && (model->BSIM4v7igbMod != 1))
        {   model->BSIM4v7igbMod = 0;
            printf("Warning: igbMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7tempModGiven)
            model->BSIM4v7tempMod = 0;
        else if ((model->BSIM4v7tempMod != 0) && (model->BSIM4v7tempMod != 1) 
                  && (model->BSIM4v7tempMod != 2) && (model->BSIM4v7tempMod != 3))
        {   model->BSIM4v7tempMod = 0;
            printf("Warning: tempMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4v7versionGiven) 
            model->BSIM4v7version = copy("4.7.0");
        if (!model->BSIM4v7toxrefGiven)
            model->BSIM4v7toxref = 30.0e-10;
        if (!model->BSIM4v7eotGiven)
            model->BSIM4v7eot = 15.0e-10;
        if (!model->BSIM4v7vddeotGiven)
            model->BSIM4v7vddeot = (model->BSIM4v7type == NMOS) ? 1.5 : -1.5;
        if (!model->BSIM4v7tempeotGiven)
            model->BSIM4v7tempeot = 300.15;
        if (!model->BSIM4v7leffeotGiven)
            model->BSIM4v7leffeot = 1;
        if (!model->BSIM4v7weffeotGiven)
            model->BSIM4v7weffeot = 10;        
        if (!model->BSIM4v7adosGiven)
            model->BSIM4v7ados = 1.0;
        if (!model->BSIM4v7bdosGiven)
            model->BSIM4v7bdos = 1.0;
        if (!model->BSIM4v7toxeGiven)
            model->BSIM4v7toxe = 30.0e-10;
        if (!model->BSIM4v7toxpGiven)
            model->BSIM4v7toxp = model->BSIM4v7toxe;
        if (!model->BSIM4v7toxmGiven)
            model->BSIM4v7toxm = model->BSIM4v7toxe;
        if (!model->BSIM4v7dtoxGiven)
            model->BSIM4v7dtox = 0.0;
        if (!model->BSIM4v7epsroxGiven)
            model->BSIM4v7epsrox = 3.9;

        if (!model->BSIM4v7cdscGiven)
            model->BSIM4v7cdsc = 2.4e-4;   /* unit Q/V/m^2  */
        if (!model->BSIM4v7cdscbGiven)
            model->BSIM4v7cdscb = 0.0;   /* unit Q/V/m^2  */    
            if (!model->BSIM4v7cdscdGiven)
            model->BSIM4v7cdscd = 0.0;   /* unit Q/V/m^2  */
        if (!model->BSIM4v7citGiven)
            model->BSIM4v7cit = 0.0;   /* unit Q/V/m^2  */
        if (!model->BSIM4v7nfactorGiven)
            model->BSIM4v7nfactor = 1.0;
        if (!model->BSIM4v7xjGiven)
            model->BSIM4v7xj = .15e-6;
        if (!model->BSIM4v7vsatGiven)
            model->BSIM4v7vsat = 8.0e4;    /* unit m/s */ 
        if (!model->BSIM4v7atGiven)
            model->BSIM4v7at = 3.3e4;    /* unit m/s */ 
        if (!model->BSIM4v7a0Given)
            model->BSIM4v7a0 = 1.0;  
        if (!model->BSIM4v7agsGiven)
            model->BSIM4v7ags = 0.0;
        if (!model->BSIM4v7a1Given)
            model->BSIM4v7a1 = 0.0;
        if (!model->BSIM4v7a2Given)
            model->BSIM4v7a2 = 1.0;
        if (!model->BSIM4v7ketaGiven)
            model->BSIM4v7keta = -0.047;    /* unit  / V */
        if (!model->BSIM4v7nsubGiven)
            model->BSIM4v7nsub = 6.0e16;   /* unit 1/cm3 */
        if (!model->BSIM4v7phigGiven)
            model->BSIM4v7phig = 4.05;  
        if (!model->BSIM4v7epsrgateGiven)
            model->BSIM4v7epsrgate = 11.7;  
        if (!model->BSIM4v7easubGiven)
            model->BSIM4v7easub = 4.05;  
        if (!model->BSIM4v7epsrsubGiven)
            model->BSIM4v7epsrsub = 11.7; 
        if (!model->BSIM4v7ni0subGiven)
            model->BSIM4v7ni0sub = 1.45e10;   /* unit 1/cm3 */
        if (!model->BSIM4v7bg0subGiven)
            model->BSIM4v7bg0sub =  1.16;     /* unit eV */
        if (!model->BSIM4v7tbgasubGiven)
            model->BSIM4v7tbgasub = 7.02e-4;  
        if (!model->BSIM4v7tbgbsubGiven)
            model->BSIM4v7tbgbsub = 1108.0;  
        if (!model->BSIM4v7ndepGiven)
            model->BSIM4v7ndep = 1.7e17;   /* unit 1/cm3 */
        if (!model->BSIM4v7nsdGiven)
            model->BSIM4v7nsd = 1.0e20;   /* unit 1/cm3 */
        if (!model->BSIM4v7phinGiven)
            model->BSIM4v7phin = 0.0; /* unit V */
        if (!model->BSIM4v7ngateGiven)
            model->BSIM4v7ngate = 0;   /* unit 1/cm3 */
        if (!model->BSIM4v7vbmGiven)
            model->BSIM4v7vbm = -3.0;
        if (!model->BSIM4v7xtGiven)
            model->BSIM4v7xt = 1.55e-7;
        if (!model->BSIM4v7kt1Given)
            model->BSIM4v7kt1 = -0.11;      /* unit V */
        if (!model->BSIM4v7kt1lGiven)
            model->BSIM4v7kt1l = 0.0;      /* unit V*m */
        if (!model->BSIM4v7kt2Given)
            model->BSIM4v7kt2 = 0.022;      /* No unit */
        if (!model->BSIM4v7k3Given)
            model->BSIM4v7k3 = 80.0;      
        if (!model->BSIM4v7k3bGiven)
            model->BSIM4v7k3b = 0.0;      
        if (!model->BSIM4v7w0Given)
            model->BSIM4v7w0 = 2.5e-6;    
        if (!model->BSIM4v7lpe0Given)
            model->BSIM4v7lpe0 = 1.74e-7;     
        if (!model->BSIM4v7lpebGiven)
            model->BSIM4v7lpeb = 0.0;
        if (!model->BSIM4v7dvtp0Given)
            model->BSIM4v7dvtp0 = 0.0;
        if (!model->BSIM4v7dvtp1Given)
            model->BSIM4v7dvtp1 = 0.0;
        if (!model->BSIM4v7dvtp2Given)        /* New DIBL/Rout */
            model->BSIM4v7dvtp2 = 0.0;
        if (!model->BSIM4v7dvtp3Given)
            model->BSIM4v7dvtp3 = 0.0;
        if (!model->BSIM4v7dvtp4Given)
            model->BSIM4v7dvtp4 = 0.0;
        if (!model->BSIM4v7dvtp5Given)
            model->BSIM4v7dvtp5 = 0.0;
        if (!model->BSIM4v7dvt0Given)
            model->BSIM4v7dvt0 = 2.2;    
        if (!model->BSIM4v7dvt1Given)
            model->BSIM4v7dvt1 = 0.53;      
        if (!model->BSIM4v7dvt2Given)
            model->BSIM4v7dvt2 = -0.032;   /* unit 1 / V */     

        if (!model->BSIM4v7dvt0wGiven)
            model->BSIM4v7dvt0w = 0.0;    
        if (!model->BSIM4v7dvt1wGiven)
            model->BSIM4v7dvt1w = 5.3e6;    
        if (!model->BSIM4v7dvt2wGiven)
            model->BSIM4v7dvt2w = -0.032;   

        if (!model->BSIM4v7droutGiven)
            model->BSIM4v7drout = 0.56;     
        if (!model->BSIM4v7dsubGiven)
            model->BSIM4v7dsub = model->BSIM4v7drout;     
        if (!model->BSIM4v7vth0Given)
            model->BSIM4v7vth0 = (model->BSIM4v7type == NMOS) ? 0.7 : -0.7;
        if (!model->BSIM4v7vfbGiven)
            model->BSIM4v7vfb = -1.0;
        if (!model->BSIM4v7euGiven)
            model->BSIM4v7eu = (model->BSIM4v7type == NMOS) ? 1.67 : 1.0;
        if (!model->BSIM4v7ucsGiven)
            model->BSIM4v7ucs = (model->BSIM4v7type == NMOS) ? 1.67 : 1.0;
        if (!model->BSIM4v7uaGiven)
            model->BSIM4v7ua = ((model->BSIM4v7mobMod == 2)) ? 1.0e-15 : 1.0e-9; /* unit m/V */
        if (!model->BSIM4v7ua1Given)
            model->BSIM4v7ua1 = 1.0e-9;      /* unit m/V */
        if (!model->BSIM4v7ubGiven)
            model->BSIM4v7ub = 1.0e-19;     /* unit (m/V)**2 */
        if (!model->BSIM4v7ub1Given)
            model->BSIM4v7ub1 = -1.0e-18;     /* unit (m/V)**2 */
        if (!model->BSIM4v7ucGiven)
            model->BSIM4v7uc = (model->BSIM4v7mobMod == 1) ? -0.0465 : -0.0465e-9;   
        if (!model->BSIM4v7uc1Given)
            model->BSIM4v7uc1 = (model->BSIM4v7mobMod == 1) ? -0.056 : -0.056e-9;   
        if (!model->BSIM4v7udGiven)
            model->BSIM4v7ud = 0.0;     /* unit m**(-2) */
        if (!model->BSIM4v7ud1Given)
            model->BSIM4v7ud1 = 0.0;     
        if (!model->BSIM4v7upGiven)
            model->BSIM4v7up = 0.0;     
        if (!model->BSIM4v7lpGiven)
            model->BSIM4v7lp = 1.0e-8;     
        if (!model->BSIM4v7u0Given)
            model->BSIM4v7u0 = (model->BSIM4v7type == NMOS) ? 0.067 : 0.025;
        if (!model->BSIM4v7uteGiven)
            model->BSIM4v7ute = -1.5; 
        if (!model->BSIM4v7ucsteGiven)
            model->BSIM4v7ucste = -4.775e-3;
        if (!model->BSIM4v7voffGiven)
            model->BSIM4v7voff = -0.08;
        if (!model->BSIM4v7vofflGiven)
            model->BSIM4v7voffl = 0.0;
        if (!model->BSIM4v7voffcvlGiven)
            model->BSIM4v7voffcvl = 0.0;
        if (!model->BSIM4v7minvGiven)
            model->BSIM4v7minv = 0.0;
        if (!model->BSIM4v7minvcvGiven)
            model->BSIM4v7minvcv = 0.0;
        if (!model->BSIM4v7fproutGiven)
            model->BSIM4v7fprout = 0.0;
        if (!model->BSIM4v7pditsGiven)
            model->BSIM4v7pdits = 0.0;
        if (!model->BSIM4v7pditsdGiven)
            model->BSIM4v7pditsd = 0.0;
        if (!model->BSIM4v7pditslGiven)
            model->BSIM4v7pditsl = 0.0;
        if (!model->BSIM4v7deltaGiven)  
           model->BSIM4v7delta = 0.01;
        if (!model->BSIM4v7rdswminGiven)
            model->BSIM4v7rdswmin = 0.0;
        if (!model->BSIM4v7rdwminGiven)
            model->BSIM4v7rdwmin = 0.0;
        if (!model->BSIM4v7rswminGiven)
            model->BSIM4v7rswmin = 0.0;
        if (!model->BSIM4v7rdswGiven)
            model->BSIM4v7rdsw = 200.0; /* in ohm*um */     
        if (!model->BSIM4v7rdwGiven)
            model->BSIM4v7rdw = 100.0;
        if (!model->BSIM4v7rswGiven)
            model->BSIM4v7rsw = 100.0;
        if (!model->BSIM4v7prwgGiven)
            model->BSIM4v7prwg = 1.0; /* in 1/V */
        if (!model->BSIM4v7prwbGiven)
            model->BSIM4v7prwb = 0.0;      
        if (!model->BSIM4v7prtGiven)
        if (!model->BSIM4v7prtGiven)
            model->BSIM4v7prt = 0.0;      
        if (!model->BSIM4v7eta0Given)
            model->BSIM4v7eta0 = 0.08;      /* no unit  */ 
        if (!model->BSIM4v7etabGiven)
            model->BSIM4v7etab = -0.07;      /* unit  1/V */ 
        if (!model->BSIM4v7pclmGiven)
            model->BSIM4v7pclm = 1.3;      /* no unit  */ 
        if (!model->BSIM4v7pdibl1Given)
            model->BSIM4v7pdibl1 = 0.39;    /* no unit  */
        if (!model->BSIM4v7pdibl2Given)
            model->BSIM4v7pdibl2 = 0.0086;    /* no unit  */ 
        if (!model->BSIM4v7pdiblbGiven)
            model->BSIM4v7pdiblb = 0.0;    /* 1/V  */ 
        if (!model->BSIM4v7pscbe1Given)
            model->BSIM4v7pscbe1 = 4.24e8;     
        if (!model->BSIM4v7pscbe2Given)
            model->BSIM4v7pscbe2 = 1.0e-5;    
        if (!model->BSIM4v7pvagGiven)
            model->BSIM4v7pvag = 0.0;     
        if (!model->BSIM4v7wrGiven)  
            model->BSIM4v7wr = 1.0;
        if (!model->BSIM4v7dwgGiven)  
            model->BSIM4v7dwg = 0.0;
        if (!model->BSIM4v7dwbGiven)  
            model->BSIM4v7dwb = 0.0;
        if (!model->BSIM4v7b0Given)
            model->BSIM4v7b0 = 0.0;
        if (!model->BSIM4v7b1Given)  
            model->BSIM4v7b1 = 0.0;
        if (!model->BSIM4v7alpha0Given)  
            model->BSIM4v7alpha0 = 0.0;
        if (!model->BSIM4v7alpha1Given)
            model->BSIM4v7alpha1 = 0.0;
        if (!model->BSIM4v7beta0Given)  
            model->BSIM4v7beta0 = 0.0;
        if (!model->BSIM4v7gidlModGiven)
            model->BSIM4v7gidlMod = 0;         /* v4.7 New GIDL/GISL */
        if (!model->BSIM4v7agidlGiven)
            model->BSIM4v7agidl = 0.0;
        if (!model->BSIM4v7bgidlGiven)
            model->BSIM4v7bgidl = 2.3e9; /* V/m */
        if (!model->BSIM4v7cgidlGiven)
            model->BSIM4v7cgidl = 0.5; /* V^3 */
        if (!model->BSIM4v7egidlGiven)
            model->BSIM4v7egidl = 0.8; /* V */
        if (!model->BSIM4v7rgidlGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7rgidl = 1.0;
        if (!model->BSIM4v7kgidlGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7kgidl = 0.0;
        if (!model->BSIM4v7fgidlGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7fgidl = 0.0;
        if (!model->BSIM4v7agislGiven)
        {
            if (model->BSIM4v7agidlGiven)
                model->BSIM4v7agisl = model->BSIM4v7agidl;
            else
                model->BSIM4v7agisl = 0.0;
        }
        if (!model->BSIM4v7bgislGiven)
        {
            if (model->BSIM4v7bgidlGiven)
                model->BSIM4v7bgisl = model->BSIM4v7bgidl;
            else
                model->BSIM4v7bgisl = 2.3e9; /* V/m */
        }
        if (!model->BSIM4v7cgislGiven)
        {
            if (model->BSIM4v7cgidlGiven)
                model->BSIM4v7cgisl = model->BSIM4v7cgidl;
            else
                model->BSIM4v7cgisl = 0.5; /* V^3 */
        }
        if (!model->BSIM4v7egislGiven)
        {
            if (model->BSIM4v7egidlGiven)
                model->BSIM4v7egisl = model->BSIM4v7egidl;
            else
                model->BSIM4v7egisl = 0.8; /* V */
        }
        if (!model->BSIM4v7rgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7rgisl = model->BSIM4v7rgidl;
        if (!model->BSIM4v7kgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7kgisl = model->BSIM4v7kgidl;
        if (!model->BSIM4v7fgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4v7fgisl = model->BSIM4v7fgidl;
        if (!model->BSIM4v7aigcGiven)
            model->BSIM4v7aigc = (model->BSIM4v7type == NMOS) ? 1.36e-2 : 9.80e-3;
        if (!model->BSIM4v7bigcGiven)
            model->BSIM4v7bigc = (model->BSIM4v7type == NMOS) ? 1.71e-3 : 7.59e-4;
        if (!model->BSIM4v7cigcGiven)
            model->BSIM4v7cigc = (model->BSIM4v7type == NMOS) ? 0.075 : 0.03;
        if (model->BSIM4v7aigsdGiven)
        {
            model->BSIM4v7aigs = model->BSIM4v7aigd = model->BSIM4v7aigsd;
        }
        else
        {
            model->BSIM4v7aigsd = (model->BSIM4v7type == NMOS) ? 1.36e-2 : 9.80e-3;
            if (!model->BSIM4v7aigsGiven)
                model->BSIM4v7aigs = (model->BSIM4v7type == NMOS) ? 1.36e-2 : 9.80e-3;
            if (!model->BSIM4v7aigdGiven)
                model->BSIM4v7aigd = (model->BSIM4v7type == NMOS) ? 1.36e-2 : 9.80e-3;
        }
        if (model->BSIM4v7bigsdGiven)
        {
            model->BSIM4v7bigs = model->BSIM4v7bigd = model->BSIM4v7bigsd;
        }
        else
        {
            model->BSIM4v7bigsd = (model->BSIM4v7type == NMOS) ? 1.71e-3 : 7.59e-4;
            if (!model->BSIM4v7bigsGiven)
                model->BSIM4v7bigs = (model->BSIM4v7type == NMOS) ? 1.71e-3 : 7.59e-4;
            if (!model->BSIM4v7bigdGiven)
                model->BSIM4v7bigd = (model->BSIM4v7type == NMOS) ? 1.71e-3 : 7.59e-4;
        }
        if (model->BSIM4v7cigsdGiven)
        {
            model->BSIM4v7cigs = model->BSIM4v7cigd = model->BSIM4v7cigsd;
        }
        else
        {
             model->BSIM4v7cigsd = (model->BSIM4v7type == NMOS) ? 0.075 : 0.03;
           if (!model->BSIM4v7cigsGiven)
                model->BSIM4v7cigs = (model->BSIM4v7type == NMOS) ? 0.075 : 0.03;
            if (!model->BSIM4v7cigdGiven)
                model->BSIM4v7cigd = (model->BSIM4v7type == NMOS) ? 0.075 : 0.03;
        }
        if (!model->BSIM4v7aigbaccGiven)
            model->BSIM4v7aigbacc = 1.36e-2;
        if (!model->BSIM4v7bigbaccGiven)
            model->BSIM4v7bigbacc = 1.71e-3;
        if (!model->BSIM4v7cigbaccGiven)
            model->BSIM4v7cigbacc = 0.075;
        if (!model->BSIM4v7aigbinvGiven)
            model->BSIM4v7aigbinv = 1.11e-2;
        if (!model->BSIM4v7bigbinvGiven)
            model->BSIM4v7bigbinv = 9.49e-4;
        if (!model->BSIM4v7cigbinvGiven)
            model->BSIM4v7cigbinv = 0.006;
        if (!model->BSIM4v7nigcGiven)
            model->BSIM4v7nigc = 1.0;
        if (!model->BSIM4v7nigbinvGiven)
            model->BSIM4v7nigbinv = 3.0;
        if (!model->BSIM4v7nigbaccGiven)
            model->BSIM4v7nigbacc = 1.0;
        if (!model->BSIM4v7ntoxGiven)
            model->BSIM4v7ntox = 1.0;
        if (!model->BSIM4v7eigbinvGiven)
            model->BSIM4v7eigbinv = 1.1;
        if (!model->BSIM4v7pigcdGiven)
            model->BSIM4v7pigcd = 1.0;
        if (!model->BSIM4v7poxedgeGiven)
            model->BSIM4v7poxedge = 1.0;
        if (!model->BSIM4v7xrcrg1Given)
            model->BSIM4v7xrcrg1 = 12.0;
        if (!model->BSIM4v7xrcrg2Given)
            model->BSIM4v7xrcrg2 = 1.0;
        if (!model->BSIM4v7ijthsfwdGiven)
            model->BSIM4v7ijthsfwd = 0.1; /* unit A */
        if (!model->BSIM4v7ijthdfwdGiven)
            model->BSIM4v7ijthdfwd = model->BSIM4v7ijthsfwd;
        if (!model->BSIM4v7ijthsrevGiven)
            model->BSIM4v7ijthsrev = 0.1; /* unit A */
        if (!model->BSIM4v7ijthdrevGiven)
            model->BSIM4v7ijthdrev = model->BSIM4v7ijthsrev;
        if (!model->BSIM4v7tnoiaGiven)
            model->BSIM4v7tnoia = 1.5;
        if (!model->BSIM4v7tnoibGiven)
            model->BSIM4v7tnoib = 3.5;
        if (!model->BSIM4v7tnoicGiven)
            model->BSIM4v7tnoic = 0.0;
        if (!model->BSIM4v7rnoiaGiven)
            model->BSIM4v7rnoia = 0.577;
        if (!model->BSIM4v7rnoibGiven)
            model->BSIM4v7rnoib = 0.5164;
        if (!model->BSIM4v7rnoicGiven)
            model->BSIM4v7rnoic = 0.395;
        if (!model->BSIM4v7ntnoiGiven)
            model->BSIM4v7ntnoi = 1.0;
        if (!model->BSIM4v7lambdaGiven)
            model->BSIM4v7lambda = 0.0;
        if (!model->BSIM4v7vtlGiven)
            model->BSIM4v7vtl = 2.0e5;    /* unit m/s */ 
        if (!model->BSIM4v7xnGiven)
            model->BSIM4v7xn = 3.0;   
        if (!model->BSIM4v7lcGiven)
            model->BSIM4v7lc = 5.0e-9;   
        if (!model->BSIM4v7vfbsdoffGiven)  
            model->BSIM4v7vfbsdoff = 0.0;  /* unit v */  
        if (!model->BSIM4v7tvfbsdoffGiven)
            model->BSIM4v7tvfbsdoff = 0.0;  
        if (!model->BSIM4v7tvoffGiven)
            model->BSIM4v7tvoff = 0.0;  
        if (!model->BSIM4v7tnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4v7tnfactor = 0.0;  
        if (!model->BSIM4v7teta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7teta0 = 0.0;  
        if (!model->BSIM4v7tvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7tvoffcv = 0.0;  
 
        if (!model->BSIM4v7lintnoiGiven)
            model->BSIM4v7lintnoi = 0.0;  /* unit m */  

        if (!model->BSIM4v7xjbvsGiven)
            model->BSIM4v7xjbvs = 1.0; /* no unit */
        if (!model->BSIM4v7xjbvdGiven)
            model->BSIM4v7xjbvd = model->BSIM4v7xjbvs;
        if (!model->BSIM4v7bvsGiven)
            model->BSIM4v7bvs = 10.0; /* V */
        if (!model->BSIM4v7bvdGiven)
            model->BSIM4v7bvd = model->BSIM4v7bvs;

        if (!model->BSIM4v7gbminGiven)
            model->BSIM4v7gbmin = 1.0e-12; /* in mho */
        if (!model->BSIM4v7rbdbGiven)
            model->BSIM4v7rbdb = 50.0; /* in ohm */
        if (!model->BSIM4v7rbpbGiven)
            model->BSIM4v7rbpb = 50.0;
        if (!model->BSIM4v7rbsbGiven)
            model->BSIM4v7rbsb = 50.0;
        if (!model->BSIM4v7rbpsGiven)
            model->BSIM4v7rbps = 50.0;
        if (!model->BSIM4v7rbpdGiven)
            model->BSIM4v7rbpd = 50.0;

        if (!model->BSIM4v7rbps0Given)
            model->BSIM4v7rbps0 = 50.0;
        if (!model->BSIM4v7rbpslGiven)
            model->BSIM4v7rbpsl = 0.0;
        if (!model->BSIM4v7rbpswGiven)
            model->BSIM4v7rbpsw = 0.0;
        if (!model->BSIM4v7rbpsnfGiven)
            model->BSIM4v7rbpsnf = 0.0;

        if (!model->BSIM4v7rbpd0Given)
            model->BSIM4v7rbpd0 = 50.0;
        if (!model->BSIM4v7rbpdlGiven)
            model->BSIM4v7rbpdl = 0.0;
        if (!model->BSIM4v7rbpdwGiven)
            model->BSIM4v7rbpdw = 0.0;
        if (!model->BSIM4v7rbpdnfGiven)
            model->BSIM4v7rbpdnf = 0.0;

        if (!model->BSIM4v7rbpbx0Given)
            model->BSIM4v7rbpbx0 = 100.0;
        if (!model->BSIM4v7rbpbxlGiven)
            model->BSIM4v7rbpbxl = 0.0;
        if (!model->BSIM4v7rbpbxwGiven)
            model->BSIM4v7rbpbxw = 0.0;
        if (!model->BSIM4v7rbpbxnfGiven)
            model->BSIM4v7rbpbxnf = 0.0;
        if (!model->BSIM4v7rbpby0Given)
            model->BSIM4v7rbpby0 = 100.0;
        if (!model->BSIM4v7rbpbylGiven)
            model->BSIM4v7rbpbyl = 0.0;
        if (!model->BSIM4v7rbpbywGiven)
            model->BSIM4v7rbpbyw = 0.0;
        if (!model->BSIM4v7rbpbynfGiven)
            model->BSIM4v7rbpbynf = 0.0;


        if (!model->BSIM4v7rbsbx0Given)
            model->BSIM4v7rbsbx0 = 100.0;
        if (!model->BSIM4v7rbsby0Given)
            model->BSIM4v7rbsby0 = 100.0;
        if (!model->BSIM4v7rbdbx0Given)
            model->BSIM4v7rbdbx0 = 100.0;
        if (!model->BSIM4v7rbdby0Given)
            model->BSIM4v7rbdby0 = 100.0;


        if (!model->BSIM4v7rbsdbxlGiven)
            model->BSIM4v7rbsdbxl = 0.0;
        if (!model->BSIM4v7rbsdbxwGiven)
            model->BSIM4v7rbsdbxw = 0.0;
        if (!model->BSIM4v7rbsdbxnfGiven)
            model->BSIM4v7rbsdbxnf = 0.0;
        if (!model->BSIM4v7rbsdbylGiven)
            model->BSIM4v7rbsdbyl = 0.0;
        if (!model->BSIM4v7rbsdbywGiven)
            model->BSIM4v7rbsdbyw = 0.0;
        if (!model->BSIM4v7rbsdbynfGiven)
            model->BSIM4v7rbsdbynf = 0.0;

        if (!model->BSIM4v7cgslGiven)  
            model->BSIM4v7cgsl = 0.0;
        if (!model->BSIM4v7cgdlGiven)  
            model->BSIM4v7cgdl = 0.0;
        if (!model->BSIM4v7ckappasGiven)  
            model->BSIM4v7ckappas = 0.6;
        if (!model->BSIM4v7ckappadGiven)
            model->BSIM4v7ckappad = model->BSIM4v7ckappas;
        if (!model->BSIM4v7clcGiven)  
            model->BSIM4v7clc = 0.1e-6;
        if (!model->BSIM4v7cleGiven)  
            model->BSIM4v7cle = 0.6;
        if (!model->BSIM4v7vfbcvGiven)  
            model->BSIM4v7vfbcv = -1.0;
        if (!model->BSIM4v7acdeGiven)
            model->BSIM4v7acde = 1.0;
        if (!model->BSIM4v7moinGiven)
            model->BSIM4v7moin = 15.0;
        if (!model->BSIM4v7noffGiven)
            model->BSIM4v7noff = 1.0;
        if (!model->BSIM4v7voffcvGiven)
            model->BSIM4v7voffcv = 0.0;
        if (!model->BSIM4v7dmcgGiven)
            model->BSIM4v7dmcg = 0.0;
        if (!model->BSIM4v7dmciGiven)
            model->BSIM4v7dmci = model->BSIM4v7dmcg;
        if (!model->BSIM4v7dmdgGiven)
            model->BSIM4v7dmdg = 0.0;
        if (!model->BSIM4v7dmcgtGiven)
            model->BSIM4v7dmcgt = 0.0;
        if (!model->BSIM4v7xgwGiven)
            model->BSIM4v7xgw = 0.0;
        if (!model->BSIM4v7xglGiven)
            model->BSIM4v7xgl = 0.0;
        if (!model->BSIM4v7rshgGiven)
            model->BSIM4v7rshg = 0.1;
        if (!model->BSIM4v7ngconGiven)
            model->BSIM4v7ngcon = 1.0;
        if (!model->BSIM4v7tcjGiven)
            model->BSIM4v7tcj = 0.0;
        if (!model->BSIM4v7tpbGiven)
            model->BSIM4v7tpb = 0.0;
        if (!model->BSIM4v7tcjswGiven)
            model->BSIM4v7tcjsw = 0.0;
        if (!model->BSIM4v7tpbswGiven)
            model->BSIM4v7tpbsw = 0.0;
        if (!model->BSIM4v7tcjswgGiven)
            model->BSIM4v7tcjswg = 0.0;
        if (!model->BSIM4v7tpbswgGiven)
            model->BSIM4v7tpbswg = 0.0;

        /* Length dependence */
        if (!model->BSIM4v7lcdscGiven)
            model->BSIM4v7lcdsc = 0.0;
        if (!model->BSIM4v7lcdscbGiven)
            model->BSIM4v7lcdscb = 0.0;
            if (!model->BSIM4v7lcdscdGiven) 
            model->BSIM4v7lcdscd = 0.0;
        if (!model->BSIM4v7lcitGiven)
            model->BSIM4v7lcit = 0.0;
        if (!model->BSIM4v7lnfactorGiven)
            model->BSIM4v7lnfactor = 0.0;
        if (!model->BSIM4v7lxjGiven)
            model->BSIM4v7lxj = 0.0;
        if (!model->BSIM4v7lvsatGiven)
            model->BSIM4v7lvsat = 0.0;
        if (!model->BSIM4v7latGiven)
            model->BSIM4v7lat = 0.0;
        if (!model->BSIM4v7la0Given)
            model->BSIM4v7la0 = 0.0; 
        if (!model->BSIM4v7lagsGiven)
            model->BSIM4v7lags = 0.0;
        if (!model->BSIM4v7la1Given)
            model->BSIM4v7la1 = 0.0;
        if (!model->BSIM4v7la2Given)
            model->BSIM4v7la2 = 0.0;
        if (!model->BSIM4v7lketaGiven)
            model->BSIM4v7lketa = 0.0;
        if (!model->BSIM4v7lnsubGiven)
            model->BSIM4v7lnsub = 0.0;
        if (!model->BSIM4v7lndepGiven)
            model->BSIM4v7lndep = 0.0;
        if (!model->BSIM4v7lnsdGiven)
            model->BSIM4v7lnsd = 0.0;
        if (!model->BSIM4v7lphinGiven)
            model->BSIM4v7lphin = 0.0;
        if (!model->BSIM4v7lngateGiven)
            model->BSIM4v7lngate = 0.0;
        if (!model->BSIM4v7lvbmGiven)
            model->BSIM4v7lvbm = 0.0;
        if (!model->BSIM4v7lxtGiven)
            model->BSIM4v7lxt = 0.0;
        if (!model->BSIM4v7lkt1Given)
            model->BSIM4v7lkt1 = 0.0; 
        if (!model->BSIM4v7lkt1lGiven)
            model->BSIM4v7lkt1l = 0.0;
        if (!model->BSIM4v7lkt2Given)
            model->BSIM4v7lkt2 = 0.0;
        if (!model->BSIM4v7lk3Given)
            model->BSIM4v7lk3 = 0.0;      
        if (!model->BSIM4v7lk3bGiven)
            model->BSIM4v7lk3b = 0.0;      
        if (!model->BSIM4v7lw0Given)
            model->BSIM4v7lw0 = 0.0;    
        if (!model->BSIM4v7llpe0Given)
            model->BSIM4v7llpe0 = 0.0;
        if (!model->BSIM4v7llpebGiven)
            model->BSIM4v7llpeb = 0.0; 
        if (!model->BSIM4v7ldvtp0Given)
            model->BSIM4v7ldvtp0 = 0.0;
        if (!model->BSIM4v7ldvtp1Given)
            model->BSIM4v7ldvtp1 = 0.0;
        if (!model->BSIM4v7ldvtp2Given)        /* New DIBL/Rout */
            model->BSIM4v7ldvtp2 = 0.0;
        if (!model->BSIM4v7ldvtp3Given)
            model->BSIM4v7ldvtp3 = 0.0;
        if (!model->BSIM4v7ldvtp4Given)
            model->BSIM4v7ldvtp4 = 0.0;
        if (!model->BSIM4v7ldvtp5Given)
            model->BSIM4v7ldvtp5 = 0.0;
        if (!model->BSIM4v7ldvt0Given)
            model->BSIM4v7ldvt0 = 0.0;    
        if (!model->BSIM4v7ldvt1Given)
            model->BSIM4v7ldvt1 = 0.0;      
        if (!model->BSIM4v7ldvt2Given)
            model->BSIM4v7ldvt2 = 0.0;
        if (!model->BSIM4v7ldvt0wGiven)
            model->BSIM4v7ldvt0w = 0.0;    
        if (!model->BSIM4v7ldvt1wGiven)
            model->BSIM4v7ldvt1w = 0.0;      
        if (!model->BSIM4v7ldvt2wGiven)
            model->BSIM4v7ldvt2w = 0.0;
        if (!model->BSIM4v7ldroutGiven)
            model->BSIM4v7ldrout = 0.0;     
        if (!model->BSIM4v7ldsubGiven)
            model->BSIM4v7ldsub = 0.0;
        if (!model->BSIM4v7lvth0Given)
           model->BSIM4v7lvth0 = 0.0;
        if (!model->BSIM4v7luaGiven)
            model->BSIM4v7lua = 0.0;
        if (!model->BSIM4v7lua1Given)
            model->BSIM4v7lua1 = 0.0;
        if (!model->BSIM4v7lubGiven)
            model->BSIM4v7lub = 0.0;
        if (!model->BSIM4v7lub1Given)
            model->BSIM4v7lub1 = 0.0;
        if (!model->BSIM4v7lucGiven)
            model->BSIM4v7luc = 0.0;
        if (!model->BSIM4v7luc1Given)
            model->BSIM4v7luc1 = 0.0;
        if (!model->BSIM4v7ludGiven)
            model->BSIM4v7lud = 0.0;
        if (!model->BSIM4v7lud1Given)
            model->BSIM4v7lud1 = 0.0;
        if (!model->BSIM4v7lupGiven)
            model->BSIM4v7lup = 0.0;
        if (!model->BSIM4v7llpGiven)
            model->BSIM4v7llp = 0.0;
        if (!model->BSIM4v7lu0Given)
            model->BSIM4v7lu0 = 0.0;
        if (!model->BSIM4v7luteGiven)
            model->BSIM4v7lute = 0.0;  
          if (!model->BSIM4v7lucsteGiven)
            model->BSIM4v7lucste = 0.0;                 
        if (!model->BSIM4v7lvoffGiven)
            model->BSIM4v7lvoff = 0.0;
        if (!model->BSIM4v7lminvGiven)
            model->BSIM4v7lminv = 0.0;
        if (!model->BSIM4v7lminvcvGiven)
            model->BSIM4v7lminvcv = 0.0;
        if (!model->BSIM4v7lfproutGiven)
            model->BSIM4v7lfprout = 0.0;
        if (!model->BSIM4v7lpditsGiven)
            model->BSIM4v7lpdits = 0.0;
        if (!model->BSIM4v7lpditsdGiven)
            model->BSIM4v7lpditsd = 0.0;
        if (!model->BSIM4v7ldeltaGiven)  
            model->BSIM4v7ldelta = 0.0;
        if (!model->BSIM4v7lrdswGiven)
            model->BSIM4v7lrdsw = 0.0;
        if (!model->BSIM4v7lrdwGiven)
            model->BSIM4v7lrdw = 0.0;
        if (!model->BSIM4v7lrswGiven)
            model->BSIM4v7lrsw = 0.0;
        if (!model->BSIM4v7lprwbGiven)
            model->BSIM4v7lprwb = 0.0;
        if (!model->BSIM4v7lprwgGiven)
            model->BSIM4v7lprwg = 0.0;
        if (!model->BSIM4v7lprtGiven)
            model->BSIM4v7lprt = 0.0;
        if (!model->BSIM4v7leta0Given)
            model->BSIM4v7leta0 = 0.0;
        if (!model->BSIM4v7letabGiven)
            model->BSIM4v7letab = -0.0;
        if (!model->BSIM4v7lpclmGiven)
            model->BSIM4v7lpclm = 0.0; 
        if (!model->BSIM4v7lpdibl1Given)
            model->BSIM4v7lpdibl1 = 0.0;
        if (!model->BSIM4v7lpdibl2Given)
            model->BSIM4v7lpdibl2 = 0.0;
        if (!model->BSIM4v7lpdiblbGiven)
            model->BSIM4v7lpdiblb = 0.0;
        if (!model->BSIM4v7lpscbe1Given)
            model->BSIM4v7lpscbe1 = 0.0;
        if (!model->BSIM4v7lpscbe2Given)
            model->BSIM4v7lpscbe2 = 0.0;
        if (!model->BSIM4v7lpvagGiven)
            model->BSIM4v7lpvag = 0.0;     
        if (!model->BSIM4v7lwrGiven)  
            model->BSIM4v7lwr = 0.0;
        if (!model->BSIM4v7ldwgGiven)  
            model->BSIM4v7ldwg = 0.0;
        if (!model->BSIM4v7ldwbGiven)  
            model->BSIM4v7ldwb = 0.0;
        if (!model->BSIM4v7lb0Given)
            model->BSIM4v7lb0 = 0.0;
        if (!model->BSIM4v7lb1Given)  
            model->BSIM4v7lb1 = 0.0;
        if (!model->BSIM4v7lalpha0Given)  
            model->BSIM4v7lalpha0 = 0.0;
        if (!model->BSIM4v7lalpha1Given)
            model->BSIM4v7lalpha1 = 0.0;
        if (!model->BSIM4v7lbeta0Given)  
            model->BSIM4v7lbeta0 = 0.0;
        if (!model->BSIM4v7lagidlGiven)
            model->BSIM4v7lagidl = 0.0;
        if (!model->BSIM4v7lbgidlGiven)
            model->BSIM4v7lbgidl = 0.0;
        if (!model->BSIM4v7lcgidlGiven)
            model->BSIM4v7lcgidl = 0.0;
        if (!model->BSIM4v7legidlGiven)
            model->BSIM4v7legidl = 0.0;
        if (!model->BSIM4v7lrgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7lrgidl = 0.0;
        if (!model->BSIM4v7lkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7lkgidl = 0.0;
        if (!model->BSIM4v7lfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7lfgidl = 0.0;

        if (!model->BSIM4v7lagislGiven)
        {
            if (model->BSIM4v7lagidlGiven)
                model->BSIM4v7lagisl = model->BSIM4v7lagidl;
            else
                model->BSIM4v7lagisl = 0.0;
        }
        if (!model->BSIM4v7lbgislGiven)
        {
            if (model->BSIM4v7lbgidlGiven)
                model->BSIM4v7lbgisl = model->BSIM4v7lbgidl;
            else
                model->BSIM4v7lbgisl = 0.0;
        }
        if (!model->BSIM4v7lcgislGiven)
        {
            if (model->BSIM4v7lcgidlGiven)
                model->BSIM4v7lcgisl = model->BSIM4v7lcgidl;
            else
                model->BSIM4v7lcgisl = 0.0;
        }
        if (!model->BSIM4v7legislGiven)
        {
            if (model->BSIM4v7legidlGiven)
                model->BSIM4v7legisl = model->BSIM4v7legidl;
            else
                model->BSIM4v7legisl = 0.0; 
        }
        if (!model->BSIM4v7lrgislGiven)         /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7lrgidlGiven)
                model->BSIM4v7lrgisl = model->BSIM4v7lrgidl;
        }
        if (!model->BSIM4v7lkgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7lkgidlGiven)
                model->BSIM4v7lkgisl = model->BSIM4v7lkgidl;
        }
        if (!model->BSIM4v7lfgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7lfgidlGiven)
                model->BSIM4v7lfgisl = model->BSIM4v7lfgidl;
        }
        if (!model->BSIM4v7laigcGiven)
            model->BSIM4v7laigc = 0.0;
        if (!model->BSIM4v7lbigcGiven)
            model->BSIM4v7lbigc = 0.0;
        if (!model->BSIM4v7lcigcGiven)
            model->BSIM4v7lcigc = 0.0;
        if (!model->BSIM4v7aigsdGiven && (model->BSIM4v7aigsGiven || model->BSIM4v7aigdGiven))
        {
            if (!model->BSIM4v7laigsGiven)
                model->BSIM4v7laigs = 0.0;
            if (!model->BSIM4v7laigdGiven)
                model->BSIM4v7laigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7laigsdGiven)
               model->BSIM4v7laigsd = 0.0;
           model->BSIM4v7laigs = model->BSIM4v7laigd = model->BSIM4v7laigsd;
        }
        if (!model->BSIM4v7bigsdGiven && (model->BSIM4v7bigsGiven || model->BSIM4v7bigdGiven))
        {
            if (!model->BSIM4v7lbigsGiven)
                model->BSIM4v7lbigs = 0.0;
            if (!model->BSIM4v7lbigdGiven)
                model->BSIM4v7lbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7lbigsdGiven)
               model->BSIM4v7lbigsd = 0.0;
           model->BSIM4v7lbigs = model->BSIM4v7lbigd = model->BSIM4v7lbigsd;
        }
        if (!model->BSIM4v7cigsdGiven && (model->BSIM4v7cigsGiven || model->BSIM4v7cigdGiven))
        {
            if (!model->BSIM4v7lcigsGiven)
                model->BSIM4v7lcigs = 0.0;
            if (!model->BSIM4v7lcigdGiven)
                model->BSIM4v7lcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7lcigsdGiven)
               model->BSIM4v7lcigsd = 0.0;
           model->BSIM4v7lcigs = model->BSIM4v7lcigd = model->BSIM4v7lcigsd;
        }
        if (!model->BSIM4v7laigbaccGiven)
            model->BSIM4v7laigbacc = 0.0;
        if (!model->BSIM4v7lbigbaccGiven)
            model->BSIM4v7lbigbacc = 0.0;
        if (!model->BSIM4v7lcigbaccGiven)
            model->BSIM4v7lcigbacc = 0.0;
        if (!model->BSIM4v7laigbinvGiven)
            model->BSIM4v7laigbinv = 0.0;
        if (!model->BSIM4v7lbigbinvGiven)
            model->BSIM4v7lbigbinv = 0.0;
        if (!model->BSIM4v7lcigbinvGiven)
            model->BSIM4v7lcigbinv = 0.0;
        if (!model->BSIM4v7lnigcGiven)
            model->BSIM4v7lnigc = 0.0;
        if (!model->BSIM4v7lnigbinvGiven)
            model->BSIM4v7lnigbinv = 0.0;
        if (!model->BSIM4v7lnigbaccGiven)
            model->BSIM4v7lnigbacc = 0.0;
        if (!model->BSIM4v7lntoxGiven)
            model->BSIM4v7lntox = 0.0;
        if (!model->BSIM4v7leigbinvGiven)
            model->BSIM4v7leigbinv = 0.0;
        if (!model->BSIM4v7lpigcdGiven)
            model->BSIM4v7lpigcd = 0.0;
        if (!model->BSIM4v7lpoxedgeGiven)
            model->BSIM4v7lpoxedge = 0.0;
        if (!model->BSIM4v7lxrcrg1Given)
            model->BSIM4v7lxrcrg1 = 0.0;
        if (!model->BSIM4v7lxrcrg2Given)
            model->BSIM4v7lxrcrg2 = 0.0;
        if (!model->BSIM4v7leuGiven)
            model->BSIM4v7leu = 0.0;
                if (!model->BSIM4v7lucsGiven)
            model->BSIM4v7lucs = 0.0;
        if (!model->BSIM4v7lvfbGiven)
            model->BSIM4v7lvfb = 0.0;
        if (!model->BSIM4v7llambdaGiven)
            model->BSIM4v7llambda = 0.0;
        if (!model->BSIM4v7lvtlGiven)
            model->BSIM4v7lvtl = 0.0;  
        if (!model->BSIM4v7lxnGiven)
            model->BSIM4v7lxn = 0.0;  
        if (!model->BSIM4v7lvfbsdoffGiven)
            model->BSIM4v7lvfbsdoff = 0.0;   
        if (!model->BSIM4v7ltvfbsdoffGiven)
            model->BSIM4v7ltvfbsdoff = 0.0;  
        if (!model->BSIM4v7ltvoffGiven)
            model->BSIM4v7ltvoff = 0.0;  
        if (!model->BSIM4v7ltnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4v7ltnfactor = 0.0;  
        if (!model->BSIM4v7lteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7lteta0 = 0.0;  
        if (!model->BSIM4v7ltvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7ltvoffcv = 0.0;  


        if (!model->BSIM4v7lcgslGiven)  
            model->BSIM4v7lcgsl = 0.0;
        if (!model->BSIM4v7lcgdlGiven)  
            model->BSIM4v7lcgdl = 0.0;
        if (!model->BSIM4v7lckappasGiven)  
            model->BSIM4v7lckappas = 0.0;
        if (!model->BSIM4v7lckappadGiven)
            model->BSIM4v7lckappad = 0.0;
        if (!model->BSIM4v7lclcGiven)  
            model->BSIM4v7lclc = 0.0;
        if (!model->BSIM4v7lcleGiven)  
            model->BSIM4v7lcle = 0.0;
        if (!model->BSIM4v7lcfGiven)  
            model->BSIM4v7lcf = 0.0;
        if (!model->BSIM4v7lvfbcvGiven)  
            model->BSIM4v7lvfbcv = 0.0;
        if (!model->BSIM4v7lacdeGiven)
            model->BSIM4v7lacde = 0.0;
        if (!model->BSIM4v7lmoinGiven)
            model->BSIM4v7lmoin = 0.0;
        if (!model->BSIM4v7lnoffGiven)
            model->BSIM4v7lnoff = 0.0;
        if (!model->BSIM4v7lvoffcvGiven)
            model->BSIM4v7lvoffcv = 0.0;

        /* Width dependence */
        if (!model->BSIM4v7wcdscGiven)
            model->BSIM4v7wcdsc = 0.0;
        if (!model->BSIM4v7wcdscbGiven)
            model->BSIM4v7wcdscb = 0.0;  
            if (!model->BSIM4v7wcdscdGiven)
            model->BSIM4v7wcdscd = 0.0;
        if (!model->BSIM4v7wcitGiven)
            model->BSIM4v7wcit = 0.0;
        if (!model->BSIM4v7wnfactorGiven)
            model->BSIM4v7wnfactor = 0.0;
        if (!model->BSIM4v7wxjGiven)
            model->BSIM4v7wxj = 0.0;
        if (!model->BSIM4v7wvsatGiven)
            model->BSIM4v7wvsat = 0.0;
        if (!model->BSIM4v7watGiven)
            model->BSIM4v7wat = 0.0;
        if (!model->BSIM4v7wa0Given)
            model->BSIM4v7wa0 = 0.0; 
        if (!model->BSIM4v7wagsGiven)
            model->BSIM4v7wags = 0.0;
        if (!model->BSIM4v7wa1Given)
            model->BSIM4v7wa1 = 0.0;
        if (!model->BSIM4v7wa2Given)
            model->BSIM4v7wa2 = 0.0;
        if (!model->BSIM4v7wketaGiven)
            model->BSIM4v7wketa = 0.0;
        if (!model->BSIM4v7wnsubGiven)
            model->BSIM4v7wnsub = 0.0;
        if (!model->BSIM4v7wndepGiven)
            model->BSIM4v7wndep = 0.0;
        if (!model->BSIM4v7wnsdGiven)
            model->BSIM4v7wnsd = 0.0;
        if (!model->BSIM4v7wphinGiven)
            model->BSIM4v7wphin = 0.0;
        if (!model->BSIM4v7wngateGiven)
            model->BSIM4v7wngate = 0.0;
        if (!model->BSIM4v7wvbmGiven)
            model->BSIM4v7wvbm = 0.0;
        if (!model->BSIM4v7wxtGiven)
            model->BSIM4v7wxt = 0.0;
        if (!model->BSIM4v7wkt1Given)
            model->BSIM4v7wkt1 = 0.0; 
        if (!model->BSIM4v7wkt1lGiven)
            model->BSIM4v7wkt1l = 0.0;
        if (!model->BSIM4v7wkt2Given)
            model->BSIM4v7wkt2 = 0.0;
        if (!model->BSIM4v7wk3Given)
            model->BSIM4v7wk3 = 0.0;      
        if (!model->BSIM4v7wk3bGiven)
            model->BSIM4v7wk3b = 0.0;      
        if (!model->BSIM4v7ww0Given)
            model->BSIM4v7ww0 = 0.0;    
        if (!model->BSIM4v7wlpe0Given)
            model->BSIM4v7wlpe0 = 0.0;
        if (!model->BSIM4v7wlpebGiven)
            model->BSIM4v7wlpeb = 0.0; 
        if (!model->BSIM4v7wdvtp0Given)
            model->BSIM4v7wdvtp0 = 0.0;
        if (!model->BSIM4v7wdvtp1Given)
            model->BSIM4v7wdvtp1 = 0.0;
        if (!model->BSIM4v7wdvtp2Given)        /* New DIBL/Rout */
            model->BSIM4v7wdvtp2 = 0.0;
        if (!model->BSIM4v7wdvtp3Given)
            model->BSIM4v7wdvtp3 = 0.0;
        if (!model->BSIM4v7wdvtp4Given)
            model->BSIM4v7wdvtp4 = 0.0;
        if (!model->BSIM4v7wdvtp5Given)
            model->BSIM4v7wdvtp5 = 0.0;
        if (!model->BSIM4v7wdvt0Given)
            model->BSIM4v7wdvt0 = 0.0;    
        if (!model->BSIM4v7wdvt1Given)
            model->BSIM4v7wdvt1 = 0.0;      
        if (!model->BSIM4v7wdvt2Given)
            model->BSIM4v7wdvt2 = 0.0;
        if (!model->BSIM4v7wdvt0wGiven)
            model->BSIM4v7wdvt0w = 0.0;    
        if (!model->BSIM4v7wdvt1wGiven)
            model->BSIM4v7wdvt1w = 0.0;      
        if (!model->BSIM4v7wdvt2wGiven)
            model->BSIM4v7wdvt2w = 0.0;
        if (!model->BSIM4v7wdroutGiven)
            model->BSIM4v7wdrout = 0.0;     
        if (!model->BSIM4v7wdsubGiven)
            model->BSIM4v7wdsub = 0.0;
        if (!model->BSIM4v7wvth0Given)
           model->BSIM4v7wvth0 = 0.0;
        if (!model->BSIM4v7wuaGiven)
            model->BSIM4v7wua = 0.0;
        if (!model->BSIM4v7wua1Given)
            model->BSIM4v7wua1 = 0.0;
        if (!model->BSIM4v7wubGiven)
            model->BSIM4v7wub = 0.0;
        if (!model->BSIM4v7wub1Given)
            model->BSIM4v7wub1 = 0.0;
        if (!model->BSIM4v7wucGiven)
            model->BSIM4v7wuc = 0.0;
        if (!model->BSIM4v7wuc1Given)
            model->BSIM4v7wuc1 = 0.0;
        if (!model->BSIM4v7wudGiven)
            model->BSIM4v7wud = 0.0;
        if (!model->BSIM4v7wud1Given)
            model->BSIM4v7wud1 = 0.0;
        if (!model->BSIM4v7wupGiven)
            model->BSIM4v7wup = 0.0;
        if (!model->BSIM4v7wlpGiven)
            model->BSIM4v7wlp = 0.0;
        if (!model->BSIM4v7wu0Given)
            model->BSIM4v7wu0 = 0.0;
        if (!model->BSIM4v7wuteGiven)
                model->BSIM4v7wute = 0.0; 
        if (!model->BSIM4v7wucsteGiven)
                model->BSIM4v7wucste = 0.0;                 
        if (!model->BSIM4v7wvoffGiven)
                model->BSIM4v7wvoff = 0.0;
        if (!model->BSIM4v7wminvGiven)
            model->BSIM4v7wminv = 0.0;
        if (!model->BSIM4v7wminvcvGiven)
            model->BSIM4v7wminvcv = 0.0;
        if (!model->BSIM4v7wfproutGiven)
            model->BSIM4v7wfprout = 0.0;
        if (!model->BSIM4v7wpditsGiven)
            model->BSIM4v7wpdits = 0.0;
        if (!model->BSIM4v7wpditsdGiven)
            model->BSIM4v7wpditsd = 0.0;
        if (!model->BSIM4v7wdeltaGiven)  
            model->BSIM4v7wdelta = 0.0;
        if (!model->BSIM4v7wrdswGiven)
            model->BSIM4v7wrdsw = 0.0;
        if (!model->BSIM4v7wrdwGiven)
            model->BSIM4v7wrdw = 0.0;
        if (!model->BSIM4v7wrswGiven)
            model->BSIM4v7wrsw = 0.0;
        if (!model->BSIM4v7wprwbGiven)
            model->BSIM4v7wprwb = 0.0;
        if (!model->BSIM4v7wprwgGiven)
            model->BSIM4v7wprwg = 0.0;
        if (!model->BSIM4v7wprtGiven)
            model->BSIM4v7wprt = 0.0;
        if (!model->BSIM4v7weta0Given)
            model->BSIM4v7weta0 = 0.0;
        if (!model->BSIM4v7wetabGiven)
            model->BSIM4v7wetab = 0.0;
        if (!model->BSIM4v7wpclmGiven)
            model->BSIM4v7wpclm = 0.0; 
        if (!model->BSIM4v7wpdibl1Given)
            model->BSIM4v7wpdibl1 = 0.0;
        if (!model->BSIM4v7wpdibl2Given)
            model->BSIM4v7wpdibl2 = 0.0;
        if (!model->BSIM4v7wpdiblbGiven)
            model->BSIM4v7wpdiblb = 0.0;
        if (!model->BSIM4v7wpscbe1Given)
            model->BSIM4v7wpscbe1 = 0.0;
        if (!model->BSIM4v7wpscbe2Given)
            model->BSIM4v7wpscbe2 = 0.0;
        if (!model->BSIM4v7wpvagGiven)
            model->BSIM4v7wpvag = 0.0;     
        if (!model->BSIM4v7wwrGiven)  
            model->BSIM4v7wwr = 0.0;
        if (!model->BSIM4v7wdwgGiven)  
            model->BSIM4v7wdwg = 0.0;
        if (!model->BSIM4v7wdwbGiven)  
            model->BSIM4v7wdwb = 0.0;
        if (!model->BSIM4v7wb0Given)
            model->BSIM4v7wb0 = 0.0;
        if (!model->BSIM4v7wb1Given)  
            model->BSIM4v7wb1 = 0.0;
        if (!model->BSIM4v7walpha0Given)  
            model->BSIM4v7walpha0 = 0.0;
        if (!model->BSIM4v7walpha1Given)
            model->BSIM4v7walpha1 = 0.0;
        if (!model->BSIM4v7wbeta0Given)  
            model->BSIM4v7wbeta0 = 0.0;
        if (!model->BSIM4v7wagidlGiven)
            model->BSIM4v7wagidl = 0.0;
        if (!model->BSIM4v7wbgidlGiven)
            model->BSIM4v7wbgidl = 0.0;
        if (!model->BSIM4v7wcgidlGiven)
            model->BSIM4v7wcgidl = 0.0;
        if (!model->BSIM4v7wegidlGiven)
            model->BSIM4v7wegidl = 0.0;
        if (!model->BSIM4v7wrgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7wrgidl = 0.0;
        if (!model->BSIM4v7wkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7wkgidl = 0.0;
        if (!model->BSIM4v7wfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7wfgidl = 0.0;

        if (!model->BSIM4v7wagislGiven)
        {
            if (model->BSIM4v7wagidlGiven)
                model->BSIM4v7wagisl = model->BSIM4v7wagidl;
            else
                model->BSIM4v7wagisl = 0.0;
        }
        if (!model->BSIM4v7wbgislGiven)
        {
            if (model->BSIM4v7wbgidlGiven)
                model->BSIM4v7wbgisl = model->BSIM4v7wbgidl;
            else
                model->BSIM4v7wbgisl = 0.0;
        }
        if (!model->BSIM4v7wcgislGiven)
        {
            if (model->BSIM4v7wcgidlGiven)
                model->BSIM4v7wcgisl = model->BSIM4v7wcgidl;
            else
                model->BSIM4v7wcgisl = 0.0;
        }
        if (!model->BSIM4v7wegislGiven)
        {
            if (model->BSIM4v7wegidlGiven)
                model->BSIM4v7wegisl = model->BSIM4v7wegidl;
            else
                model->BSIM4v7wegisl = 0.0; 
        }
        if (!model->BSIM4v7wrgislGiven)         /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7wrgidlGiven)
                model->BSIM4v7wrgisl = model->BSIM4v7wrgidl;
        }
        if (!model->BSIM4v7wkgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7wkgidlGiven)
                model->BSIM4v7wkgisl = model->BSIM4v7wkgidl;
        }
        if (!model->BSIM4v7wfgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7wfgidlGiven)
                model->BSIM4v7wfgisl = model->BSIM4v7wfgidl;
        }
        if (!model->BSIM4v7waigcGiven)
            model->BSIM4v7waigc = 0.0;
        if (!model->BSIM4v7wbigcGiven)
            model->BSIM4v7wbigc = 0.0;
        if (!model->BSIM4v7wcigcGiven)
            model->BSIM4v7wcigc = 0.0;
        if (!model->BSIM4v7aigsdGiven && (model->BSIM4v7aigsGiven || model->BSIM4v7aigdGiven))
        {
            if (!model->BSIM4v7waigsGiven)
                model->BSIM4v7waigs = 0.0;
            if (!model->BSIM4v7waigdGiven)
                model->BSIM4v7waigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7waigsdGiven)
               model->BSIM4v7waigsd = 0.0;
           model->BSIM4v7waigs = model->BSIM4v7waigd = model->BSIM4v7waigsd;
        }
        if (!model->BSIM4v7bigsdGiven && (model->BSIM4v7bigsGiven || model->BSIM4v7bigdGiven))
        {
            if (!model->BSIM4v7wbigsGiven)
                model->BSIM4v7wbigs = 0.0;
            if (!model->BSIM4v7wbigdGiven)
                model->BSIM4v7wbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7wbigsdGiven)
               model->BSIM4v7wbigsd = 0.0;
           model->BSIM4v7wbigs = model->BSIM4v7wbigd = model->BSIM4v7wbigsd;
        }
        if (!model->BSIM4v7cigsdGiven && (model->BSIM4v7cigsGiven || model->BSIM4v7cigdGiven))
        {
            if (!model->BSIM4v7wcigsGiven)
                model->BSIM4v7wcigs = 0.0;
            if (!model->BSIM4v7wcigdGiven)
                model->BSIM4v7wcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7wcigsdGiven)
               model->BSIM4v7wcigsd = 0.0;
           model->BSIM4v7wcigs = model->BSIM4v7wcigd = model->BSIM4v7wcigsd;
        }
        if (!model->BSIM4v7waigbaccGiven)
            model->BSIM4v7waigbacc = 0.0;
        if (!model->BSIM4v7wbigbaccGiven)
            model->BSIM4v7wbigbacc = 0.0;
        if (!model->BSIM4v7wcigbaccGiven)
            model->BSIM4v7wcigbacc = 0.0;
        if (!model->BSIM4v7waigbinvGiven)
            model->BSIM4v7waigbinv = 0.0;
        if (!model->BSIM4v7wbigbinvGiven)
            model->BSIM4v7wbigbinv = 0.0;
        if (!model->BSIM4v7wcigbinvGiven)
            model->BSIM4v7wcigbinv = 0.0;
        if (!model->BSIM4v7wnigcGiven)
            model->BSIM4v7wnigc = 0.0;
        if (!model->BSIM4v7wnigbinvGiven)
            model->BSIM4v7wnigbinv = 0.0;
        if (!model->BSIM4v7wnigbaccGiven)
            model->BSIM4v7wnigbacc = 0.0;
        if (!model->BSIM4v7wntoxGiven)
            model->BSIM4v7wntox = 0.0;
        if (!model->BSIM4v7weigbinvGiven)
            model->BSIM4v7weigbinv = 0.0;
        if (!model->BSIM4v7wpigcdGiven)
            model->BSIM4v7wpigcd = 0.0;
        if (!model->BSIM4v7wpoxedgeGiven)
            model->BSIM4v7wpoxedge = 0.0;
        if (!model->BSIM4v7wxrcrg1Given)
            model->BSIM4v7wxrcrg1 = 0.0;
        if (!model->BSIM4v7wxrcrg2Given)
            model->BSIM4v7wxrcrg2 = 0.0;
        if (!model->BSIM4v7weuGiven)
            model->BSIM4v7weu = 0.0;
            if (!model->BSIM4v7wucsGiven)
            model->BSIM4v7wucs = 0.0;
        if (!model->BSIM4v7wvfbGiven)
            model->BSIM4v7wvfb = 0.0;
        if (!model->BSIM4v7wlambdaGiven)
            model->BSIM4v7wlambda = 0.0;
        if (!model->BSIM4v7wvtlGiven)
            model->BSIM4v7wvtl = 0.0;  
        if (!model->BSIM4v7wxnGiven)
            model->BSIM4v7wxn = 0.0;  
        if (!model->BSIM4v7wvfbsdoffGiven)
            model->BSIM4v7wvfbsdoff = 0.0;   
        if (!model->BSIM4v7wtvfbsdoffGiven)
            model->BSIM4v7wtvfbsdoff = 0.0;  
        if (!model->BSIM4v7wtvoffGiven)
            model->BSIM4v7wtvoff = 0.0;  
        if (!model->BSIM4v7wtnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4v7wtnfactor = 0.0;  
        if (!model->BSIM4v7wteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7wteta0 = 0.0;  
        if (!model->BSIM4v7wtvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7wtvoffcv = 0.0;  

        if (!model->BSIM4v7wcgslGiven)  
            model->BSIM4v7wcgsl = 0.0;
        if (!model->BSIM4v7wcgdlGiven)  
            model->BSIM4v7wcgdl = 0.0;
        if (!model->BSIM4v7wckappasGiven)  
            model->BSIM4v7wckappas = 0.0;
        if (!model->BSIM4v7wckappadGiven)
            model->BSIM4v7wckappad = 0.0;
        if (!model->BSIM4v7wcfGiven)  
            model->BSIM4v7wcf = 0.0;
        if (!model->BSIM4v7wclcGiven)  
            model->BSIM4v7wclc = 0.0;
        if (!model->BSIM4v7wcleGiven)  
            model->BSIM4v7wcle = 0.0;
        if (!model->BSIM4v7wvfbcvGiven)  
            model->BSIM4v7wvfbcv = 0.0;
        if (!model->BSIM4v7wacdeGiven)
            model->BSIM4v7wacde = 0.0;
        if (!model->BSIM4v7wmoinGiven)
            model->BSIM4v7wmoin = 0.0;
        if (!model->BSIM4v7wnoffGiven)
            model->BSIM4v7wnoff = 0.0;
        if (!model->BSIM4v7wvoffcvGiven)
            model->BSIM4v7wvoffcv = 0.0;

        /* Cross-term dependence */
        if (!model->BSIM4v7pcdscGiven)
            model->BSIM4v7pcdsc = 0.0;
        if (!model->BSIM4v7pcdscbGiven)
            model->BSIM4v7pcdscb = 0.0;   
            if (!model->BSIM4v7pcdscdGiven)
            model->BSIM4v7pcdscd = 0.0;
        if (!model->BSIM4v7pcitGiven)
            model->BSIM4v7pcit = 0.0;
        if (!model->BSIM4v7pnfactorGiven)
            model->BSIM4v7pnfactor = 0.0;
        if (!model->BSIM4v7pxjGiven)
            model->BSIM4v7pxj = 0.0;
        if (!model->BSIM4v7pvsatGiven)
            model->BSIM4v7pvsat = 0.0;
        if (!model->BSIM4v7patGiven)
            model->BSIM4v7pat = 0.0;
        if (!model->BSIM4v7pa0Given)
            model->BSIM4v7pa0 = 0.0; 
            
        if (!model->BSIM4v7pagsGiven)
            model->BSIM4v7pags = 0.0;
        if (!model->BSIM4v7pa1Given)
            model->BSIM4v7pa1 = 0.0;
        if (!model->BSIM4v7pa2Given)
            model->BSIM4v7pa2 = 0.0;
        if (!model->BSIM4v7pketaGiven)
            model->BSIM4v7pketa = 0.0;
        if (!model->BSIM4v7pnsubGiven)
            model->BSIM4v7pnsub = 0.0;
        if (!model->BSIM4v7pndepGiven)
            model->BSIM4v7pndep = 0.0;
        if (!model->BSIM4v7pnsdGiven)
            model->BSIM4v7pnsd = 0.0;
        if (!model->BSIM4v7pphinGiven)
            model->BSIM4v7pphin = 0.0;
        if (!model->BSIM4v7pngateGiven)
            model->BSIM4v7pngate = 0.0;
        if (!model->BSIM4v7pvbmGiven)
            model->BSIM4v7pvbm = 0.0;
        if (!model->BSIM4v7pxtGiven)
            model->BSIM4v7pxt = 0.0;
        if (!model->BSIM4v7pkt1Given)
            model->BSIM4v7pkt1 = 0.0; 
        if (!model->BSIM4v7pkt1lGiven)
            model->BSIM4v7pkt1l = 0.0;
        if (!model->BSIM4v7pkt2Given)
            model->BSIM4v7pkt2 = 0.0;
        if (!model->BSIM4v7pk3Given)
            model->BSIM4v7pk3 = 0.0;      
        if (!model->BSIM4v7pk3bGiven)
            model->BSIM4v7pk3b = 0.0;      
        if (!model->BSIM4v7pw0Given)
            model->BSIM4v7pw0 = 0.0;    
        if (!model->BSIM4v7plpe0Given)
            model->BSIM4v7plpe0 = 0.0;
        if (!model->BSIM4v7plpebGiven)
            model->BSIM4v7plpeb = 0.0;
        if (!model->BSIM4v7pdvtp0Given)
            model->BSIM4v7pdvtp0 = 0.0;
        if (!model->BSIM4v7pdvtp1Given)
            model->BSIM4v7pdvtp1 = 0.0;
        if (!model->BSIM4v7pdvtp2Given)        /* New DIBL/Rout */
            model->BSIM4v7pdvtp2 = 0.0;
        if (!model->BSIM4v7pdvtp3Given)
            model->BSIM4v7pdvtp3 = 0.0;
        if (!model->BSIM4v7pdvtp4Given)
            model->BSIM4v7pdvtp4 = 0.0;
        if (!model->BSIM4v7pdvtp5Given)
            model->BSIM4v7pdvtp5 = 0.0;
        if (!model->BSIM4v7pdvt0Given)
            model->BSIM4v7pdvt0 = 0.0;    
        if (!model->BSIM4v7pdvt1Given)
            model->BSIM4v7pdvt1 = 0.0;      
        if (!model->BSIM4v7pdvt2Given)
            model->BSIM4v7pdvt2 = 0.0;
        if (!model->BSIM4v7pdvt0wGiven)
            model->BSIM4v7pdvt0w = 0.0;    
        if (!model->BSIM4v7pdvt1wGiven)
            model->BSIM4v7pdvt1w = 0.0;      
        if (!model->BSIM4v7pdvt2wGiven)
            model->BSIM4v7pdvt2w = 0.0;
        if (!model->BSIM4v7pdroutGiven)
            model->BSIM4v7pdrout = 0.0;     
        if (!model->BSIM4v7pdsubGiven)
            model->BSIM4v7pdsub = 0.0;
        if (!model->BSIM4v7pvth0Given)
           model->BSIM4v7pvth0 = 0.0;
        if (!model->BSIM4v7puaGiven)
            model->BSIM4v7pua = 0.0;
        if (!model->BSIM4v7pua1Given)
            model->BSIM4v7pua1 = 0.0;
        if (!model->BSIM4v7pubGiven)
            model->BSIM4v7pub = 0.0;
        if (!model->BSIM4v7pub1Given)
            model->BSIM4v7pub1 = 0.0;
        if (!model->BSIM4v7pucGiven)
            model->BSIM4v7puc = 0.0;
        if (!model->BSIM4v7puc1Given)
            model->BSIM4v7puc1 = 0.0;
        if (!model->BSIM4v7pudGiven)
            model->BSIM4v7pud = 0.0;
        if (!model->BSIM4v7pud1Given)
            model->BSIM4v7pud1 = 0.0;
        if (!model->BSIM4v7pupGiven)
            model->BSIM4v7pup = 0.0;
        if (!model->BSIM4v7plpGiven)
            model->BSIM4v7plp = 0.0;
        if (!model->BSIM4v7pu0Given)
            model->BSIM4v7pu0 = 0.0;
        if (!model->BSIM4v7puteGiven)
            model->BSIM4v7pute = 0.0;  
     if (!model->BSIM4v7pucsteGiven)
            model->BSIM4v7pucste = 0.0;                 
        if (!model->BSIM4v7pvoffGiven)
            model->BSIM4v7pvoff = 0.0;
        if (!model->BSIM4v7pminvGiven)
            model->BSIM4v7pminv = 0.0;
        if (!model->BSIM4v7pminvcvGiven)
            model->BSIM4v7pminvcv = 0.0;
        if (!model->BSIM4v7pfproutGiven)
            model->BSIM4v7pfprout = 0.0;
        if (!model->BSIM4v7ppditsGiven)
            model->BSIM4v7ppdits = 0.0;
        if (!model->BSIM4v7ppditsdGiven)
            model->BSIM4v7ppditsd = 0.0;
        if (!model->BSIM4v7pdeltaGiven)  
            model->BSIM4v7pdelta = 0.0;
        if (!model->BSIM4v7prdswGiven)
            model->BSIM4v7prdsw = 0.0;
        if (!model->BSIM4v7prdwGiven)
            model->BSIM4v7prdw = 0.0;
        if (!model->BSIM4v7prswGiven)
            model->BSIM4v7prsw = 0.0;
        if (!model->BSIM4v7pprwbGiven)
            model->BSIM4v7pprwb = 0.0;
        if (!model->BSIM4v7pprwgGiven)
            model->BSIM4v7pprwg = 0.0;
        if (!model->BSIM4v7pprtGiven)
            model->BSIM4v7pprt = 0.0;
        if (!model->BSIM4v7peta0Given)
            model->BSIM4v7peta0 = 0.0;
        if (!model->BSIM4v7petabGiven)
            model->BSIM4v7petab = 0.0;
        if (!model->BSIM4v7ppclmGiven)
            model->BSIM4v7ppclm = 0.0; 
        if (!model->BSIM4v7ppdibl1Given)
            model->BSIM4v7ppdibl1 = 0.0;
        if (!model->BSIM4v7ppdibl2Given)
            model->BSIM4v7ppdibl2 = 0.0;
        if (!model->BSIM4v7ppdiblbGiven)
            model->BSIM4v7ppdiblb = 0.0;
        if (!model->BSIM4v7ppscbe1Given)
            model->BSIM4v7ppscbe1 = 0.0;
        if (!model->BSIM4v7ppscbe2Given)
            model->BSIM4v7ppscbe2 = 0.0;
        if (!model->BSIM4v7ppvagGiven)
            model->BSIM4v7ppvag = 0.0;     
        if (!model->BSIM4v7pwrGiven)  
            model->BSIM4v7pwr = 0.0;
        if (!model->BSIM4v7pdwgGiven)  
            model->BSIM4v7pdwg = 0.0;
        if (!model->BSIM4v7pdwbGiven)  
            model->BSIM4v7pdwb = 0.0;
        if (!model->BSIM4v7pb0Given)
            model->BSIM4v7pb0 = 0.0;
        if (!model->BSIM4v7pb1Given)  
            model->BSIM4v7pb1 = 0.0;
        if (!model->BSIM4v7palpha0Given)  
            model->BSIM4v7palpha0 = 0.0;
        if (!model->BSIM4v7palpha1Given)
            model->BSIM4v7palpha1 = 0.0;
        if (!model->BSIM4v7pbeta0Given)  
            model->BSIM4v7pbeta0 = 0.0;
        if (!model->BSIM4v7pagidlGiven)
            model->BSIM4v7pagidl = 0.0;
        if (!model->BSIM4v7pbgidlGiven)
            model->BSIM4v7pbgidl = 0.0;
        if (!model->BSIM4v7pcgidlGiven)
            model->BSIM4v7pcgidl = 0.0;
        if (!model->BSIM4v7pegidlGiven)
            model->BSIM4v7pegidl = 0.0;
        if (!model->BSIM4v7prgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7prgidl = 0.0;
        if (!model->BSIM4v7pkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7pkgidl = 0.0;
        if (!model->BSIM4v7pfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4v7pfgidl = 0.0;

        if (!model->BSIM4v7pagislGiven)
        {
            if (model->BSIM4v7pagidlGiven)
                model->BSIM4v7pagisl = model->BSIM4v7pagidl;
            else
                model->BSIM4v7pagisl = 0.0;
        }
        if (!model->BSIM4v7pbgislGiven)
        {
            if (model->BSIM4v7pbgidlGiven)
                model->BSIM4v7pbgisl = model->BSIM4v7pbgidl;
            else
                model->BSIM4v7pbgisl = 0.0;
        }
        if (!model->BSIM4v7pcgislGiven)
        {
            if (model->BSIM4v7pcgidlGiven)
                model->BSIM4v7pcgisl = model->BSIM4v7pcgidl;
            else
                model->BSIM4v7pcgisl = 0.0;
        }
        if (!model->BSIM4v7pegislGiven)
        {
            if (model->BSIM4v7pegidlGiven)
                model->BSIM4v7pegisl = model->BSIM4v7pegidl;
            else
                model->BSIM4v7pegisl = 0.0; 
        }
        if (!model->BSIM4v7prgislGiven)         /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7prgidlGiven)
                model->BSIM4v7prgisl = model->BSIM4v7prgidl;
        }
        if (!model->BSIM4v7pkgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7pkgidlGiven)
                model->BSIM4v7pkgisl = model->BSIM4v7pkgidl;
        }
        if (!model->BSIM4v7pfgislGiven)        /* v4.7 New GIDL/GISL */
        {
            if (model->BSIM4v7pfgidlGiven)
                model->BSIM4v7pfgisl = model->BSIM4v7pfgidl;
        }
        if (!model->BSIM4v7paigcGiven)
            model->BSIM4v7paigc = 0.0;
        if (!model->BSIM4v7pbigcGiven)
            model->BSIM4v7pbigc = 0.0;
        if (!model->BSIM4v7pcigcGiven)
            model->BSIM4v7pcigc = 0.0;
        if (!model->BSIM4v7aigsdGiven && (model->BSIM4v7aigsGiven || model->BSIM4v7aigdGiven))
        {
            if (!model->BSIM4v7paigsGiven)
                model->BSIM4v7paigs = 0.0;
            if (!model->BSIM4v7paigdGiven)
                model->BSIM4v7paigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7paigsdGiven)
               model->BSIM4v7paigsd = 0.0;
           model->BSIM4v7paigs = model->BSIM4v7paigd = model->BSIM4v7paigsd;
        }
        if (!model->BSIM4v7bigsdGiven && (model->BSIM4v7bigsGiven || model->BSIM4v7bigdGiven))
        {
            if (!model->BSIM4v7pbigsGiven)
                model->BSIM4v7pbigs = 0.0;
            if (!model->BSIM4v7pbigdGiven)
                model->BSIM4v7pbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7pbigsdGiven)
               model->BSIM4v7pbigsd = 0.0;
           model->BSIM4v7pbigs = model->BSIM4v7pbigd = model->BSIM4v7pbigsd;
        }
        if (!model->BSIM4v7cigsdGiven && (model->BSIM4v7cigsGiven || model->BSIM4v7cigdGiven))
        {
            if (!model->BSIM4v7pcigsGiven)
                model->BSIM4v7pcigs = 0.0;
            if (!model->BSIM4v7pcigdGiven)
                model->BSIM4v7pcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4v7pcigsdGiven)
               model->BSIM4v7pcigsd = 0.0;
           model->BSIM4v7pcigs = model->BSIM4v7pcigd = model->BSIM4v7pcigsd;
        }
        if (!model->BSIM4v7paigbaccGiven)
            model->BSIM4v7paigbacc = 0.0;
        if (!model->BSIM4v7pbigbaccGiven)
            model->BSIM4v7pbigbacc = 0.0;
        if (!model->BSIM4v7pcigbaccGiven)
            model->BSIM4v7pcigbacc = 0.0;
        if (!model->BSIM4v7paigbinvGiven)
            model->BSIM4v7paigbinv = 0.0;
        if (!model->BSIM4v7pbigbinvGiven)
            model->BSIM4v7pbigbinv = 0.0;
        if (!model->BSIM4v7pcigbinvGiven)
            model->BSIM4v7pcigbinv = 0.0;
        if (!model->BSIM4v7pnigcGiven)
            model->BSIM4v7pnigc = 0.0;
        if (!model->BSIM4v7pnigbinvGiven)
            model->BSIM4v7pnigbinv = 0.0;
        if (!model->BSIM4v7pnigbaccGiven)
            model->BSIM4v7pnigbacc = 0.0;
        if (!model->BSIM4v7pntoxGiven)
            model->BSIM4v7pntox = 0.0;
        if (!model->BSIM4v7peigbinvGiven)
            model->BSIM4v7peigbinv = 0.0;
        if (!model->BSIM4v7ppigcdGiven)
            model->BSIM4v7ppigcd = 0.0;
        if (!model->BSIM4v7ppoxedgeGiven)
            model->BSIM4v7ppoxedge = 0.0;
        if (!model->BSIM4v7pxrcrg1Given)
            model->BSIM4v7pxrcrg1 = 0.0;
        if (!model->BSIM4v7pxrcrg2Given)
            model->BSIM4v7pxrcrg2 = 0.0;
        if (!model->BSIM4v7peuGiven)
            model->BSIM4v7peu = 0.0;
                if (!model->BSIM4v7pucsGiven)
            model->BSIM4v7pucs = 0.0;
        if (!model->BSIM4v7pvfbGiven)
            model->BSIM4v7pvfb = 0.0;
        if (!model->BSIM4v7plambdaGiven)
            model->BSIM4v7plambda = 0.0;
        if (!model->BSIM4v7pvtlGiven)
            model->BSIM4v7pvtl = 0.0;  
        if (!model->BSIM4v7pxnGiven)
            model->BSIM4v7pxn = 0.0;  
        if (!model->BSIM4v7pvfbsdoffGiven)
            model->BSIM4v7pvfbsdoff = 0.0;   
        if (!model->BSIM4v7ptvfbsdoffGiven)
            model->BSIM4v7ptvfbsdoff = 0.0;  
        if (!model->BSIM4v7ptvoffGiven)
            model->BSIM4v7ptvoff = 0.0;  
        if (!model->BSIM4v7ptnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4v7ptnfactor = 0.0;  
        if (!model->BSIM4v7pteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7pteta0 = 0.0;  
        if (!model->BSIM4v7ptvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4v7ptvoffcv = 0.0;  

        if (!model->BSIM4v7pcgslGiven)  
            model->BSIM4v7pcgsl = 0.0;
        if (!model->BSIM4v7pcgdlGiven)  
            model->BSIM4v7pcgdl = 0.0;
        if (!model->BSIM4v7pckappasGiven)  
            model->BSIM4v7pckappas = 0.0;
        if (!model->BSIM4v7pckappadGiven)
            model->BSIM4v7pckappad = 0.0;
        if (!model->BSIM4v7pcfGiven)  
            model->BSIM4v7pcf = 0.0;
        if (!model->BSIM4v7pclcGiven)  
            model->BSIM4v7pclc = 0.0;
        if (!model->BSIM4v7pcleGiven)  
            model->BSIM4v7pcle = 0.0;
        if (!model->BSIM4v7pvfbcvGiven)  
            model->BSIM4v7pvfbcv = 0.0;
        if (!model->BSIM4v7pacdeGiven)
            model->BSIM4v7pacde = 0.0;
        if (!model->BSIM4v7pmoinGiven)
            model->BSIM4v7pmoin = 0.0;
        if (!model->BSIM4v7pnoffGiven)
            model->BSIM4v7pnoff = 0.0;
        if (!model->BSIM4v7pvoffcvGiven)
            model->BSIM4v7pvoffcv = 0.0;

        if (!model->BSIM4v7gamma1Given)
            model->BSIM4v7gamma1 = 0.0;
        if (!model->BSIM4v7lgamma1Given)
            model->BSIM4v7lgamma1 = 0.0;
        if (!model->BSIM4v7wgamma1Given)
            model->BSIM4v7wgamma1 = 0.0;
        if (!model->BSIM4v7pgamma1Given)
            model->BSIM4v7pgamma1 = 0.0;
        if (!model->BSIM4v7gamma2Given)
            model->BSIM4v7gamma2 = 0.0;
        if (!model->BSIM4v7lgamma2Given)
            model->BSIM4v7lgamma2 = 0.0;
        if (!model->BSIM4v7wgamma2Given)
            model->BSIM4v7wgamma2 = 0.0;
        if (!model->BSIM4v7pgamma2Given)
            model->BSIM4v7pgamma2 = 0.0;
        if (!model->BSIM4v7vbxGiven)
            model->BSIM4v7vbx = 0.0;
        if (!model->BSIM4v7lvbxGiven)
            model->BSIM4v7lvbx = 0.0;
        if (!model->BSIM4v7wvbxGiven)
            model->BSIM4v7wvbx = 0.0;
        if (!model->BSIM4v7pvbxGiven)
            model->BSIM4v7pvbx = 0.0;

        /* unit degree celcius */
        if (!model->BSIM4v7tnomGiven)  
            model->BSIM4v7tnom = ckt->CKTnomTemp; 
        if (!model->BSIM4v7LintGiven)  
           model->BSIM4v7Lint = 0.0;
        if (!model->BSIM4v7LlGiven)  
           model->BSIM4v7Ll = 0.0;
        if (!model->BSIM4v7LlcGiven)
           model->BSIM4v7Llc = model->BSIM4v7Ll;
        if (!model->BSIM4v7LlnGiven)  
           model->BSIM4v7Lln = 1.0;
        if (!model->BSIM4v7LwGiven)  
           model->BSIM4v7Lw = 0.0;
        if (!model->BSIM4v7LwcGiven)
           model->BSIM4v7Lwc = model->BSIM4v7Lw;
        if (!model->BSIM4v7LwnGiven)  
           model->BSIM4v7Lwn = 1.0;
        if (!model->BSIM4v7LwlGiven)  
           model->BSIM4v7Lwl = 0.0;
        if (!model->BSIM4v7LwlcGiven)
           model->BSIM4v7Lwlc = model->BSIM4v7Lwl;
        if (!model->BSIM4v7LminGiven)  
           model->BSIM4v7Lmin = 0.0;
        if (!model->BSIM4v7LmaxGiven)  
           model->BSIM4v7Lmax = 1.0;
        if (!model->BSIM4v7WintGiven)  
           model->BSIM4v7Wint = 0.0;
        if (!model->BSIM4v7WlGiven)  
           model->BSIM4v7Wl = 0.0;
        if (!model->BSIM4v7WlcGiven)
           model->BSIM4v7Wlc = model->BSIM4v7Wl;
        if (!model->BSIM4v7WlnGiven)  
           model->BSIM4v7Wln = 1.0;
        if (!model->BSIM4v7WwGiven)  
           model->BSIM4v7Ww = 0.0;
        if (!model->BSIM4v7WwcGiven)
           model->BSIM4v7Wwc = model->BSIM4v7Ww;
        if (!model->BSIM4v7WwnGiven)  
           model->BSIM4v7Wwn = 1.0;
        if (!model->BSIM4v7WwlGiven)  
           model->BSIM4v7Wwl = 0.0;
        if (!model->BSIM4v7WwlcGiven)
           model->BSIM4v7Wwlc = model->BSIM4v7Wwl;
        if (!model->BSIM4v7WminGiven)  
           model->BSIM4v7Wmin = 0.0;
        if (!model->BSIM4v7WmaxGiven)  
           model->BSIM4v7Wmax = 1.0;
        if (!model->BSIM4v7dwcGiven)  
           model->BSIM4v7dwc = model->BSIM4v7Wint;
        if (!model->BSIM4v7dlcGiven)  
           model->BSIM4v7dlc = model->BSIM4v7Lint;
        if (!model->BSIM4v7xlGiven)  
           model->BSIM4v7xl = 0.0;
        if (!model->BSIM4v7xwGiven)  
           model->BSIM4v7xw = 0.0;
        if (!model->BSIM4v7dlcigGiven)
           model->BSIM4v7dlcig = model->BSIM4v7Lint;
        if (!model->BSIM4v7dlcigdGiven)
        {
           if (model->BSIM4v7dlcigGiven) 
               model->BSIM4v7dlcigd = model->BSIM4v7dlcig;
           else             
               model->BSIM4v7dlcigd = model->BSIM4v7Lint;
        }
        if (!model->BSIM4v7dwjGiven)
           model->BSIM4v7dwj = model->BSIM4v7dwc;
        if (!model->BSIM4v7cfGiven)
           model->BSIM4v7cf = 2.0 * model->BSIM4v7epsrox * EPS0 / PI
                          * log(1.0 + 0.4e-6 / model->BSIM4v7toxe);

        if (!model->BSIM4v7xpartGiven)
            model->BSIM4v7xpart = 0.0;
        if (!model->BSIM4v7sheetResistanceGiven)
            model->BSIM4v7sheetResistance = 0.0;

        if (!model->BSIM4v7SunitAreaJctCapGiven)
            model->BSIM4v7SunitAreaJctCap = 5.0E-4;
        if (!model->BSIM4v7DunitAreaJctCapGiven)
            model->BSIM4v7DunitAreaJctCap = model->BSIM4v7SunitAreaJctCap;
        if (!model->BSIM4v7SunitLengthSidewallJctCapGiven)
            model->BSIM4v7SunitLengthSidewallJctCap = 5.0E-10;
        if (!model->BSIM4v7DunitLengthSidewallJctCapGiven)
            model->BSIM4v7DunitLengthSidewallJctCap = model->BSIM4v7SunitLengthSidewallJctCap;
        if (!model->BSIM4v7SunitLengthGateSidewallJctCapGiven)
            model->BSIM4v7SunitLengthGateSidewallJctCap = model->BSIM4v7SunitLengthSidewallJctCap ;
        if (!model->BSIM4v7DunitLengthGateSidewallJctCapGiven)
            model->BSIM4v7DunitLengthGateSidewallJctCap = model->BSIM4v7SunitLengthGateSidewallJctCap;
        if (!model->BSIM4v7SjctSatCurDensityGiven)
            model->BSIM4v7SjctSatCurDensity = 1.0E-4;
        if (!model->BSIM4v7DjctSatCurDensityGiven)
            model->BSIM4v7DjctSatCurDensity = model->BSIM4v7SjctSatCurDensity;
        if (!model->BSIM4v7SjctSidewallSatCurDensityGiven)
            model->BSIM4v7SjctSidewallSatCurDensity = 0.0;
        if (!model->BSIM4v7DjctSidewallSatCurDensityGiven)
            model->BSIM4v7DjctSidewallSatCurDensity = model->BSIM4v7SjctSidewallSatCurDensity;
        if (!model->BSIM4v7SjctGateSidewallSatCurDensityGiven)
            model->BSIM4v7SjctGateSidewallSatCurDensity = 0.0;
        if (!model->BSIM4v7DjctGateSidewallSatCurDensityGiven)
            model->BSIM4v7DjctGateSidewallSatCurDensity = model->BSIM4v7SjctGateSidewallSatCurDensity;
        if (!model->BSIM4v7SbulkJctPotentialGiven)
            model->BSIM4v7SbulkJctPotential = 1.0;
        if (!model->BSIM4v7DbulkJctPotentialGiven)
            model->BSIM4v7DbulkJctPotential = model->BSIM4v7SbulkJctPotential;
        if (!model->BSIM4v7SsidewallJctPotentialGiven)
            model->BSIM4v7SsidewallJctPotential = 1.0;
        if (!model->BSIM4v7DsidewallJctPotentialGiven)
            model->BSIM4v7DsidewallJctPotential = model->BSIM4v7SsidewallJctPotential;
        if (!model->BSIM4v7SGatesidewallJctPotentialGiven)
            model->BSIM4v7SGatesidewallJctPotential = model->BSIM4v7SsidewallJctPotential;
        if (!model->BSIM4v7DGatesidewallJctPotentialGiven)
            model->BSIM4v7DGatesidewallJctPotential = model->BSIM4v7SGatesidewallJctPotential;
        if (!model->BSIM4v7SbulkJctBotGradingCoeffGiven)
            model->BSIM4v7SbulkJctBotGradingCoeff = 0.5;
        if (!model->BSIM4v7DbulkJctBotGradingCoeffGiven)
            model->BSIM4v7DbulkJctBotGradingCoeff = model->BSIM4v7SbulkJctBotGradingCoeff;
        if (!model->BSIM4v7SbulkJctSideGradingCoeffGiven)
            model->BSIM4v7SbulkJctSideGradingCoeff = 0.33;
        if (!model->BSIM4v7DbulkJctSideGradingCoeffGiven)
            model->BSIM4v7DbulkJctSideGradingCoeff = model->BSIM4v7SbulkJctSideGradingCoeff;
        if (!model->BSIM4v7SbulkJctGateSideGradingCoeffGiven)
            model->BSIM4v7SbulkJctGateSideGradingCoeff = model->BSIM4v7SbulkJctSideGradingCoeff;
        if (!model->BSIM4v7DbulkJctGateSideGradingCoeffGiven)
            model->BSIM4v7DbulkJctGateSideGradingCoeff = model->BSIM4v7SbulkJctGateSideGradingCoeff;
        if (!model->BSIM4v7SjctEmissionCoeffGiven)
            model->BSIM4v7SjctEmissionCoeff = 1.0;
        if (!model->BSIM4v7DjctEmissionCoeffGiven)
            model->BSIM4v7DjctEmissionCoeff = model->BSIM4v7SjctEmissionCoeff;
        if (!model->BSIM4v7SjctTempExponentGiven)
            model->BSIM4v7SjctTempExponent = 3.0;
        if (!model->BSIM4v7DjctTempExponentGiven)
            model->BSIM4v7DjctTempExponent = model->BSIM4v7SjctTempExponent;

        if (!model->BSIM4v7jtssGiven)
            model->BSIM4v7jtss = 0.0;
        if (!model->BSIM4v7jtsdGiven)
            model->BSIM4v7jtsd = model->BSIM4v7jtss;
        if (!model->BSIM4v7jtsswsGiven)
            model->BSIM4v7jtssws = 0.0;
        if (!model->BSIM4v7jtsswdGiven)
            model->BSIM4v7jtsswd = model->BSIM4v7jtssws;
        if (!model->BSIM4v7jtsswgsGiven)
            model->BSIM4v7jtsswgs = 0.0;
        if (!model->BSIM4v7jtsswgdGiven)
            model->BSIM4v7jtsswgd = model->BSIM4v7jtsswgs;
                if (!model->BSIM4v7jtweffGiven)
                    model->BSIM4v7jtweff = 0.0;
        if (!model->BSIM4v7njtsGiven)
            model->BSIM4v7njts = 20.0;
        if (!model->BSIM4v7njtsswGiven)
            model->BSIM4v7njtssw = 20.0;
        if (!model->BSIM4v7njtsswgGiven)
            model->BSIM4v7njtsswg = 20.0;
        if (!model->BSIM4v7njtsdGiven)
        {
            if (model->BSIM4v7njtsGiven)
                model->BSIM4v7njtsd =  model->BSIM4v7njts;
            else
              model->BSIM4v7njtsd = 20.0;
        }
        if (!model->BSIM4v7njtsswdGiven)
        {
            if (model->BSIM4v7njtsswGiven)
                model->BSIM4v7njtsswd =  model->BSIM4v7njtssw;
            else
              model->BSIM4v7njtsswd = 20.0;
        }
        if (!model->BSIM4v7njtsswgdGiven)
        {
            if (model->BSIM4v7njtsswgGiven)
                model->BSIM4v7njtsswgd =  model->BSIM4v7njtsswg;
            else
              model->BSIM4v7njtsswgd = 20.0;
        }
        if (!model->BSIM4v7xtssGiven)
            model->BSIM4v7xtss = 0.02;
        if (!model->BSIM4v7xtsdGiven)
            model->BSIM4v7xtsd = model->BSIM4v7xtss;
        if (!model->BSIM4v7xtsswsGiven)
            model->BSIM4v7xtssws = 0.02;
        if (!model->BSIM4v7xtsswdGiven)
            model->BSIM4v7xtsswd = model->BSIM4v7xtssws;
        if (!model->BSIM4v7xtsswgsGiven)
            model->BSIM4v7xtsswgs = 0.02;
        if (!model->BSIM4v7xtsswgdGiven)
            model->BSIM4v7xtsswgd = model->BSIM4v7xtsswgs;
        if (!model->BSIM4v7tnjtsGiven)
            model->BSIM4v7tnjts = 0.0;
        if (!model->BSIM4v7tnjtsswGiven)
            model->BSIM4v7tnjtssw = 0.0;
        if (!model->BSIM4v7tnjtsswgGiven)
            model->BSIM4v7tnjtsswg = 0.0;
        if (!model->BSIM4v7tnjtsdGiven)
        {
            if (model->BSIM4v7tnjtsGiven)
                model->BSIM4v7tnjtsd =  model->BSIM4v7tnjts;
            else
              model->BSIM4v7tnjtsd = 0.0;
        }
        if (!model->BSIM4v7tnjtsswdGiven)
        {
            if (model->BSIM4v7tnjtsswGiven)
                model->BSIM4v7tnjtsswd =  model->BSIM4v7tnjtssw;
            else
              model->BSIM4v7tnjtsswd = 0.0;
        }
        if (!model->BSIM4v7tnjtsswgdGiven)
        {
            if (model->BSIM4v7tnjtsswgGiven)
                model->BSIM4v7tnjtsswgd =  model->BSIM4v7tnjtsswg;
            else
              model->BSIM4v7tnjtsswgd = 0.0;
        }
        if (!model->BSIM4v7vtssGiven)
            model->BSIM4v7vtss = 10.0;
        if (!model->BSIM4v7vtsdGiven)
            model->BSIM4v7vtsd = model->BSIM4v7vtss;
        if (!model->BSIM4v7vtsswsGiven)
            model->BSIM4v7vtssws = 10.0;
        if (!model->BSIM4v7vtsswdGiven)
            model->BSIM4v7vtsswd = model->BSIM4v7vtssws;
        if (!model->BSIM4v7vtsswgsGiven)
            model->BSIM4v7vtsswgs = 10.0;
        if (!model->BSIM4v7vtsswgdGiven)
            model->BSIM4v7vtsswgd = model->BSIM4v7vtsswgs;

        if (!model->BSIM4v7oxideTrapDensityAGiven)
        {   if (model->BSIM4v7type == NMOS)
                model->BSIM4v7oxideTrapDensityA = 6.25e41;
            else
                model->BSIM4v7oxideTrapDensityA= 6.188e40;
        }
        if (!model->BSIM4v7oxideTrapDensityBGiven)
        {   if (model->BSIM4v7type == NMOS)
                model->BSIM4v7oxideTrapDensityB = 3.125e26;
            else
                model->BSIM4v7oxideTrapDensityB = 1.5e25;
        }
        if (!model->BSIM4v7oxideTrapDensityCGiven)
            model->BSIM4v7oxideTrapDensityC = 8.75e9;
        if (!model->BSIM4v7emGiven)
            model->BSIM4v7em = 4.1e7; /* V/m */
        if (!model->BSIM4v7efGiven)
            model->BSIM4v7ef = 1.0;
        if (!model->BSIM4v7afGiven)
            model->BSIM4v7af = 1.0;
        if (!model->BSIM4v7kfGiven)
            model->BSIM4v7kf = 0.0;

        if (!model->BSIM4v7vgsMaxGiven)
            model->BSIM4v7vgsMax = 1e99;
        if (!model->BSIM4v7vgdMaxGiven)
            model->BSIM4v7vgdMax = 1e99;
        if (!model->BSIM4v7vgbMaxGiven)
            model->BSIM4v7vgbMax = 1e99;
        if (!model->BSIM4v7vdsMaxGiven)
            model->BSIM4v7vdsMax = 1e99;
        if (!model->BSIM4v7vbsMaxGiven)
            model->BSIM4v7vbsMax = 1e99;
        if (!model->BSIM4v7vbdMaxGiven)
            model->BSIM4v7vbdMax = 1e99;

        /* stress effect */
        if (!model->BSIM4v7sarefGiven)
            model->BSIM4v7saref = 1e-6; /* m */
        if (!model->BSIM4v7sbrefGiven)
            model->BSIM4v7sbref = 1e-6;  /* m */
        if (!model->BSIM4v7wlodGiven)
            model->BSIM4v7wlod = 0;  /* m */
        if (!model->BSIM4v7ku0Given)
            model->BSIM4v7ku0 = 0; /* 1/m */
        if (!model->BSIM4v7kvsatGiven)
            model->BSIM4v7kvsat = 0;
        if (!model->BSIM4v7kvth0Given) /* m */
            model->BSIM4v7kvth0 = 0;
        if (!model->BSIM4v7tku0Given)
            model->BSIM4v7tku0 = 0;
        if (!model->BSIM4v7llodku0Given)
            model->BSIM4v7llodku0 = 0;
        if (!model->BSIM4v7wlodku0Given)
            model->BSIM4v7wlodku0 = 0;
        if (!model->BSIM4v7llodvthGiven)
            model->BSIM4v7llodvth = 0;
        if (!model->BSIM4v7wlodvthGiven)
            model->BSIM4v7wlodvth = 0;
        if (!model->BSIM4v7lku0Given)
            model->BSIM4v7lku0 = 0;
        if (!model->BSIM4v7wku0Given)
            model->BSIM4v7wku0 = 0;
        if (!model->BSIM4v7pku0Given)
            model->BSIM4v7pku0 = 0;
        if (!model->BSIM4v7lkvth0Given)
            model->BSIM4v7lkvth0 = 0;
        if (!model->BSIM4v7wkvth0Given)
            model->BSIM4v7wkvth0 = 0;
        if (!model->BSIM4v7pkvth0Given)
            model->BSIM4v7pkvth0 = 0;
        if (!model->BSIM4v7stk2Given)
            model->BSIM4v7stk2 = 0;
        if (!model->BSIM4v7lodk2Given)
            model->BSIM4v7lodk2 = 1.0;
        if (!model->BSIM4v7steta0Given)
            model->BSIM4v7steta0 = 0;
        if (!model->BSIM4v7lodeta0Given)
            model->BSIM4v7lodeta0 = 1.0;

        /* Well Proximity Effect  */
        if (!model->BSIM4v7webGiven)
            model->BSIM4v7web = 0.0; 
        if (!model->BSIM4v7wecGiven)
            model->BSIM4v7wec = 0.0;
        if (!model->BSIM4v7kvth0weGiven)
            model->BSIM4v7kvth0we = 0.0; 
        if (!model->BSIM4v7k2weGiven)
            model->BSIM4v7k2we = 0.0; 
        if (!model->BSIM4v7ku0weGiven)
            model->BSIM4v7ku0we = 0.0; 
        if (!model->BSIM4v7screfGiven)
            model->BSIM4v7scref = 1.0E-6; /* m */
        if (!model->BSIM4v7wpemodGiven)
            model->BSIM4v7wpemod = 0; 
        else if ((model->BSIM4v7wpemod != 0) && (model->BSIM4v7wpemod != 1))
        {   model->BSIM4v7wpemod = 0;
            printf("Warning: wpemod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4v7lkvth0weGiven)
            model->BSIM4v7lkvth0we = 0; 
        if (!model->BSIM4v7lk2weGiven)
            model->BSIM4v7lk2we = 0;
        if (!model->BSIM4v7lku0weGiven)
            model->BSIM4v7lku0we = 0;
        if (!model->BSIM4v7wkvth0weGiven)
            model->BSIM4v7wkvth0we = 0; 
        if (!model->BSIM4v7wk2weGiven)
            model->BSIM4v7wk2we = 0;
        if (!model->BSIM4v7wku0weGiven)
            model->BSIM4v7wku0we = 0;
        if (!model->BSIM4v7pkvth0weGiven)
            model->BSIM4v7pkvth0we = 0; 
        if (!model->BSIM4v7pk2weGiven)
            model->BSIM4v7pk2we = 0;
        if (!model->BSIM4v7pku0weGiven)
            model->BSIM4v7pku0we = 0;

        DMCGeff = model->BSIM4v7dmcg - model->BSIM4v7dmcgt;
        DMCIeff = model->BSIM4v7dmci;
        DMDGeff = model->BSIM4v7dmdg - model->BSIM4v7dmcgt;

        /*
         * End processing models and begin to loop
         * through all the instances of the model
         */

        for (here = model->BSIM4v7instances; here != NULL ;
             here=here->BSIM4v7nextInstance) 
        {
            /* allocate a chunk of the state vector */
            here->BSIM4v7states = *states;
            *states += BSIM4v7numStates;

            /* perform the parameter defaulting */
            if (!here->BSIM4v7lGiven)
                here->BSIM4v7l = 5.0e-6;
            if (!here->BSIM4v7wGiven)
                here->BSIM4v7w = 5.0e-6;
            if (!here->BSIM4v7mGiven)
                here->BSIM4v7m = 1.0;
            if (!here->BSIM4v7nfGiven)
                here->BSIM4v7nf = 1.0;
            if (!here->BSIM4v7minGiven)
                here->BSIM4v7min = 0; /* integer */
            if (!here->BSIM4v7icVDSGiven)
                here->BSIM4v7icVDS = 0.0;
            if (!here->BSIM4v7icVGSGiven)
                here->BSIM4v7icVGS = 0.0;
            if (!here->BSIM4v7icVBSGiven)
                here->BSIM4v7icVBS = 0.0;
            if (!here->BSIM4v7drainAreaGiven)
                here->BSIM4v7drainArea = 0.0;
            if (!here->BSIM4v7drainPerimeterGiven)
                here->BSIM4v7drainPerimeter = 0.0;
            if (!here->BSIM4v7drainSquaresGiven)
                here->BSIM4v7drainSquares = 1.0;
            if (!here->BSIM4v7sourceAreaGiven)
                here->BSIM4v7sourceArea = 0.0;
            if (!here->BSIM4v7sourcePerimeterGiven)
                here->BSIM4v7sourcePerimeter = 0.0;
            if (!here->BSIM4v7sourceSquaresGiven)
                here->BSIM4v7sourceSquares = 1.0;

            if (!here->BSIM4v7rbdbGiven)
                here->BSIM4v7rbdb = model->BSIM4v7rbdb; /* in ohm */
            if (!here->BSIM4v7rbsbGiven)
                here->BSIM4v7rbsb = model->BSIM4v7rbsb;
            if (!here->BSIM4v7rbpbGiven)
                here->BSIM4v7rbpb = model->BSIM4v7rbpb;
            if (!here->BSIM4v7rbpsGiven)
                here->BSIM4v7rbps = model->BSIM4v7rbps;
            if (!here->BSIM4v7rbpdGiven)
                here->BSIM4v7rbpd = model->BSIM4v7rbpd;
            if (!here->BSIM4v7delvtoGiven)
                here->BSIM4v7delvto = 0.0;
            if (!here->BSIM4v7xgwGiven)
                here->BSIM4v7xgw = model->BSIM4v7xgw;
            if (!here->BSIM4v7ngconGiven)
                here->BSIM4v7ngcon = model->BSIM4v7ngcon;

                    
            /* Process instance model selectors, some
             * may override their global counterparts
             */
            if (!here->BSIM4v7rbodyModGiven)
                here->BSIM4v7rbodyMod = model->BSIM4v7rbodyMod;
            else if ((here->BSIM4v7rbodyMod != 0) && (here->BSIM4v7rbodyMod != 1) && (here->BSIM4v7rbodyMod != 2))
            {   here->BSIM4v7rbodyMod = model->BSIM4v7rbodyMod;
                printf("Warning: rbodyMod has been set to its global value %d.\n",
                model->BSIM4v7rbodyMod);
            }

            if (!here->BSIM4v7rgateModGiven)
                here->BSIM4v7rgateMod = model->BSIM4v7rgateMod;
            else if ((here->BSIM4v7rgateMod != 0) && (here->BSIM4v7rgateMod != 1)
                && (here->BSIM4v7rgateMod != 2) && (here->BSIM4v7rgateMod != 3))
            {   here->BSIM4v7rgateMod = model->BSIM4v7rgateMod;
                printf("Warning: rgateMod has been set to its global value %d.\n",
                model->BSIM4v7rgateMod);
            }

            if (!here->BSIM4v7geoModGiven)
                here->BSIM4v7geoMod = model->BSIM4v7geoMod;

            if (!here->BSIM4v7rgeoModGiven)
                here->BSIM4v7rgeoMod = model->BSIM4v7rgeoMod;
            else if ((here->BSIM4v7rgeoMod != 0) && (here->BSIM4v7rgeoMod != 1))
            {   here->BSIM4v7rgeoMod = model->BSIM4v7rgeoMod;
                printf("Warning: rgeoMod has been set to its global value %d.\n",
                model->BSIM4v7rgeoMod);
            }

            if (!here->BSIM4v7trnqsModGiven)
                here->BSIM4v7trnqsMod = model->BSIM4v7trnqsMod;
            else if ((here->BSIM4v7trnqsMod != 0) && (here->BSIM4v7trnqsMod != 1))
            {   here->BSIM4v7trnqsMod = model->BSIM4v7trnqsMod;
                printf("Warning: trnqsMod has been set to its global value %d.\n",
                model->BSIM4v7trnqsMod);
            }

            if (!here->BSIM4v7acnqsModGiven)
                here->BSIM4v7acnqsMod = model->BSIM4v7acnqsMod;
            else if ((here->BSIM4v7acnqsMod != 0) && (here->BSIM4v7acnqsMod != 1))
            {   here->BSIM4v7acnqsMod = model->BSIM4v7acnqsMod;
                printf("Warning: acnqsMod has been set to its global value %d.\n",
                model->BSIM4v7acnqsMod);
            }

            /* stress effect */
            if (!here->BSIM4v7saGiven)
                here->BSIM4v7sa = 0.0;
            if (!here->BSIM4v7sbGiven)
                here->BSIM4v7sb = 0.0;
            if (!here->BSIM4v7sdGiven)
                here->BSIM4v7sd = 2 * model->BSIM4v7dmcg;
            /* Well Proximity Effect  */
            if (!here->BSIM4v7scaGiven)
                here->BSIM4v7sca = 0.0;
            if (!here->BSIM4v7scbGiven)
                here->BSIM4v7scb = 0.0;
            if (!here->BSIM4v7sccGiven)
                here->BSIM4v7scc = 0.0;
            if (!here->BSIM4v7scGiven)
                here->BSIM4v7sc = 0.0; /* m */

            /* process drain series resistance */
            createNode = 0;
            if ( (model->BSIM4v7rdsMod != 0)
                            || (model->BSIM4v7tnoiMod == 1 && noiseAnalGiven))
            {
               createNode = 1;
            } else if (model->BSIM4v7sheetResistance > 0)
            {
                     if (here->BSIM4v7drainSquaresGiven
                                       && here->BSIM4v7drainSquares > 0)
                     {
                          createNode = 1;
                     } else if (!here->BSIM4v7drainSquaresGiven
                                       && (here->BSIM4v7rgeoMod != 0))
                     {
                          BSIM4v7RdseffGeo(here->BSIM4v7nf*here->BSIM4v7m, here->BSIM4v7geoMod,
                                  here->BSIM4v7rgeoMod, here->BSIM4v7min,
                                  here->BSIM4v7w, model->BSIM4v7sheetResistance,
                                  DMCGeff, DMCIeff, DMDGeff, 0, &Rtot);
                          if(Rtot > 0)
                             createNode = 1;
                     }
            }
            if ( createNode != 0 )
            {   if ( here->BSIM4v7dNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"drain");
                if(error) return(error);
                here->BSIM4v7dNodePrime = tmp->number;
                if (ckt->CKTcopyNodesets) {
                  CKTnode *tmpNode;
                  IFuid tmpName;    
                  if (CKTinst2Node(ckt,here,1,&tmpNode,&tmpName)==OK) {
                     if (tmpNode->nsGiven) {
                       tmp->nodeset=tmpNode->nodeset; 
                       tmp->nsGiven=tmpNode->nsGiven; 
                     }
                  }
                }
            }
            }
            else
            {   here->BSIM4v7dNodePrime = here->BSIM4v7dNode;
            }
            
            /* process source series resistance */
            createNode = 0;
            if ( (model->BSIM4v7rdsMod != 0)
                            || (model->BSIM4v7tnoiMod == 1 && noiseAnalGiven))
            {
               createNode = 1;
            } else if (model->BSIM4v7sheetResistance > 0)
            {
                     if (here->BSIM4v7sourceSquaresGiven
                                        && here->BSIM4v7sourceSquares > 0)
                     {
                          createNode = 1;
                     } else if (!here->BSIM4v7sourceSquaresGiven
                                        && (here->BSIM4v7rgeoMod != 0))
                     {
                          BSIM4v7RdseffGeo(here->BSIM4v7nf*here->BSIM4v7m, here->BSIM4v7geoMod,
                                  here->BSIM4v7rgeoMod, here->BSIM4v7min,
                                  here->BSIM4v7w, model->BSIM4v7sheetResistance,
                                  DMCGeff, DMCIeff, DMDGeff, 1, &Rtot);
                          if(Rtot > 0)
                             createNode = 1;
                     }
            }
            if ( createNode != 0 )
            {   if ( here->BSIM4v7sNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"source");
                if(error) return(error);
                here->BSIM4v7sNodePrime = tmp->number;
                if (ckt->CKTcopyNodesets) {
                  CKTnode *tmpNode;
                  IFuid tmpName;
                  if (CKTinst2Node(ckt,here,3,&tmpNode,&tmpName)==OK) {
                     if (tmpNode->nsGiven) {
                       tmp->nodeset=tmpNode->nodeset; 
                       tmp->nsGiven=tmpNode->nsGiven; 
                     }
                  }
                }
            }
            }
            else
                here->BSIM4v7sNodePrime = here->BSIM4v7sNode;

            if ( here->BSIM4v7rgateMod > 0 )
            {   if ( here->BSIM4v7gNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"gate");
                if(error) return(error);
                   here->BSIM4v7gNodePrime = tmp->number;
                if (ckt->CKTcopyNodesets) {
                  CKTnode *tmpNode;
                  IFuid tmpName;
                  if (CKTinst2Node(ckt,here,2,&tmpNode,&tmpName)==OK) {
                    if (tmpNode->nsGiven) {
                      tmp->nodeset=tmpNode->nodeset; 
                      tmp->nsGiven=tmpNode->nsGiven; 
                    }
                  }
                }
            }
            }
            else
                here->BSIM4v7gNodePrime = here->BSIM4v7gNodeExt;

            if ( here->BSIM4v7rgateMod == 3 )
            {   if ( here->BSIM4v7gNodeMid == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"midgate");
                if(error) return(error);
                   here->BSIM4v7gNodeMid = tmp->number;
            }
            }
            else
                here->BSIM4v7gNodeMid = here->BSIM4v7gNodeExt;
            

            /* internal body nodes for body resistance model */
            if ((here->BSIM4v7rbodyMod ==1) || (here->BSIM4v7rbodyMod ==2))
            {   if (here->BSIM4v7dbNode == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"dbody");
                    if(error) return(error);
                    here->BSIM4v7dbNode = tmp->number;
                }
                if (here->BSIM4v7bNodePrime == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"body");
                    if(error) return(error);
                    here->BSIM4v7bNodePrime = tmp->number;
                    if (ckt->CKTcopyNodesets) {
                      CKTnode *tmpNode;
                      IFuid tmpName;
                      if (CKTinst2Node(ckt,here,4,&tmpNode,&tmpName)==OK) {
                        if (tmpNode->nsGiven) {
                          tmp->nodeset=tmpNode->nodeset; 
                          tmp->nsGiven=tmpNode->nsGiven; 
                        }
                      }
                    }
                }
                if (here->BSIM4v7sbNode == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"sbody");
                    if(error) return(error);
                    here->BSIM4v7sbNode = tmp->number;
                }
            }
            else
                here->BSIM4v7dbNode = here->BSIM4v7bNodePrime = here->BSIM4v7sbNode
                                  = here->BSIM4v7bNode;

            /* NQS node */
            if ( here->BSIM4v7trnqsMod )
            {   if ( here->BSIM4v7qNode == 0 ) 
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4v7name,"charge");
                if(error) return(error);
                here->BSIM4v7qNode = tmp->number;
            }
            }
            else 
                here->BSIM4v7qNode = 0;
      

/* set Sparse Matrix Pointers 
 * macro to make elements with built-in out-of-memory test */
#define TSTALLOC(ptr,first,second) \
do { if((here->ptr = SMPmakeElt(matrix,here->first,here->second))==(double *)NULL){\
    return(E_NOMEM);\
} } while(0)

            TSTALLOC(BSIM4v7DPbpPtr, BSIM4v7dNodePrime, BSIM4v7bNodePrime);
            TSTALLOC(BSIM4v7GPbpPtr, BSIM4v7gNodePrime, BSIM4v7bNodePrime);
            TSTALLOC(BSIM4v7SPbpPtr, BSIM4v7sNodePrime, BSIM4v7bNodePrime);

            TSTALLOC(BSIM4v7BPdpPtr, BSIM4v7bNodePrime, BSIM4v7dNodePrime);
            TSTALLOC(BSIM4v7BPgpPtr, BSIM4v7bNodePrime, BSIM4v7gNodePrime);
            TSTALLOC(BSIM4v7BPspPtr, BSIM4v7bNodePrime, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7BPbpPtr, BSIM4v7bNodePrime, BSIM4v7bNodePrime);

            TSTALLOC(BSIM4v7DdPtr, BSIM4v7dNode, BSIM4v7dNode);
            TSTALLOC(BSIM4v7GPgpPtr, BSIM4v7gNodePrime, BSIM4v7gNodePrime);
            TSTALLOC(BSIM4v7SsPtr, BSIM4v7sNode, BSIM4v7sNode);
            TSTALLOC(BSIM4v7DPdpPtr, BSIM4v7dNodePrime, BSIM4v7dNodePrime);
            TSTALLOC(BSIM4v7SPspPtr, BSIM4v7sNodePrime, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7DdpPtr, BSIM4v7dNode, BSIM4v7dNodePrime);
            TSTALLOC(BSIM4v7GPdpPtr, BSIM4v7gNodePrime, BSIM4v7dNodePrime);
            TSTALLOC(BSIM4v7GPspPtr, BSIM4v7gNodePrime, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7SspPtr, BSIM4v7sNode, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7DPspPtr, BSIM4v7dNodePrime, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7DPdPtr, BSIM4v7dNodePrime, BSIM4v7dNode);
            TSTALLOC(BSIM4v7DPgpPtr, BSIM4v7dNodePrime, BSIM4v7gNodePrime);
            TSTALLOC(BSIM4v7SPgpPtr, BSIM4v7sNodePrime, BSIM4v7gNodePrime);
            TSTALLOC(BSIM4v7SPsPtr, BSIM4v7sNodePrime, BSIM4v7sNode);
            TSTALLOC(BSIM4v7SPdpPtr, BSIM4v7sNodePrime, BSIM4v7dNodePrime);

            TSTALLOC(BSIM4v7QqPtr, BSIM4v7qNode, BSIM4v7qNode);
            TSTALLOC(BSIM4v7QbpPtr, BSIM4v7qNode, BSIM4v7bNodePrime) ;
            TSTALLOC(BSIM4v7QdpPtr, BSIM4v7qNode, BSIM4v7dNodePrime);
            TSTALLOC(BSIM4v7QspPtr, BSIM4v7qNode, BSIM4v7sNodePrime);
            TSTALLOC(BSIM4v7QgpPtr, BSIM4v7qNode, BSIM4v7gNodePrime);
            TSTALLOC(BSIM4v7DPqPtr, BSIM4v7dNodePrime, BSIM4v7qNode);
            TSTALLOC(BSIM4v7SPqPtr, BSIM4v7sNodePrime, BSIM4v7qNode);
            TSTALLOC(BSIM4v7GPqPtr, BSIM4v7gNodePrime, BSIM4v7qNode);

            if (here->BSIM4v7rgateMod != 0)
            {   TSTALLOC(BSIM4v7GEgePtr, BSIM4v7gNodeExt, BSIM4v7gNodeExt);
                TSTALLOC(BSIM4v7GEgpPtr, BSIM4v7gNodeExt, BSIM4v7gNodePrime);
                TSTALLOC(BSIM4v7GPgePtr, BSIM4v7gNodePrime, BSIM4v7gNodeExt);
                TSTALLOC(BSIM4v7GEdpPtr, BSIM4v7gNodeExt, BSIM4v7dNodePrime);
                TSTALLOC(BSIM4v7GEspPtr, BSIM4v7gNodeExt, BSIM4v7sNodePrime);
                TSTALLOC(BSIM4v7GEbpPtr, BSIM4v7gNodeExt, BSIM4v7bNodePrime);

                TSTALLOC(BSIM4v7GMdpPtr, BSIM4v7gNodeMid, BSIM4v7dNodePrime);
                TSTALLOC(BSIM4v7GMgpPtr, BSIM4v7gNodeMid, BSIM4v7gNodePrime);
                TSTALLOC(BSIM4v7GMgmPtr, BSIM4v7gNodeMid, BSIM4v7gNodeMid);
                TSTALLOC(BSIM4v7GMgePtr, BSIM4v7gNodeMid, BSIM4v7gNodeExt);
                TSTALLOC(BSIM4v7GMspPtr, BSIM4v7gNodeMid, BSIM4v7sNodePrime);
                TSTALLOC(BSIM4v7GMbpPtr, BSIM4v7gNodeMid, BSIM4v7bNodePrime);
                TSTALLOC(BSIM4v7DPgmPtr, BSIM4v7dNodePrime, BSIM4v7gNodeMid);
                TSTALLOC(BSIM4v7GPgmPtr, BSIM4v7gNodePrime, BSIM4v7gNodeMid);
                TSTALLOC(BSIM4v7GEgmPtr, BSIM4v7gNodeExt, BSIM4v7gNodeMid);
                TSTALLOC(BSIM4v7SPgmPtr, BSIM4v7sNodePrime, BSIM4v7gNodeMid);
                TSTALLOC(BSIM4v7BPgmPtr, BSIM4v7bNodePrime, BSIM4v7gNodeMid);
            }        

            if ((here->BSIM4v7rbodyMod ==1) || (here->BSIM4v7rbodyMod ==2))
            {   TSTALLOC(BSIM4v7DPdbPtr, BSIM4v7dNodePrime, BSIM4v7dbNode);
                TSTALLOC(BSIM4v7SPsbPtr, BSIM4v7sNodePrime, BSIM4v7sbNode);

                TSTALLOC(BSIM4v7DBdpPtr, BSIM4v7dbNode, BSIM4v7dNodePrime);
                TSTALLOC(BSIM4v7DBdbPtr, BSIM4v7dbNode, BSIM4v7dbNode);
                TSTALLOC(BSIM4v7DBbpPtr, BSIM4v7dbNode, BSIM4v7bNodePrime);
                TSTALLOC(BSIM4v7DBbPtr, BSIM4v7dbNode, BSIM4v7bNode);

                TSTALLOC(BSIM4v7BPdbPtr, BSIM4v7bNodePrime, BSIM4v7dbNode);
                TSTALLOC(BSIM4v7BPbPtr, BSIM4v7bNodePrime, BSIM4v7bNode);
                TSTALLOC(BSIM4v7BPsbPtr, BSIM4v7bNodePrime, BSIM4v7sbNode);

                TSTALLOC(BSIM4v7SBspPtr, BSIM4v7sbNode, BSIM4v7sNodePrime);
                TSTALLOC(BSIM4v7SBbpPtr, BSIM4v7sbNode, BSIM4v7bNodePrime);
                TSTALLOC(BSIM4v7SBbPtr, BSIM4v7sbNode, BSIM4v7bNode);
                TSTALLOC(BSIM4v7SBsbPtr, BSIM4v7sbNode, BSIM4v7sbNode);

                TSTALLOC(BSIM4v7BdbPtr, BSIM4v7bNode, BSIM4v7dbNode);
                TSTALLOC(BSIM4v7BbpPtr, BSIM4v7bNode, BSIM4v7bNodePrime);
                TSTALLOC(BSIM4v7BsbPtr, BSIM4v7bNode, BSIM4v7sbNode);
                TSTALLOC(BSIM4v7BbPtr, BSIM4v7bNode, BSIM4v7bNode);
            }

            if (model->BSIM4v7rdsMod)
            {   TSTALLOC(BSIM4v7DgpPtr, BSIM4v7dNode, BSIM4v7gNodePrime);
                TSTALLOC(BSIM4v7DspPtr, BSIM4v7dNode, BSIM4v7sNodePrime);
                TSTALLOC(BSIM4v7DbpPtr, BSIM4v7dNode, BSIM4v7bNodePrime);
                TSTALLOC(BSIM4v7SdpPtr, BSIM4v7sNode, BSIM4v7dNodePrime);
                TSTALLOC(BSIM4v7SgpPtr, BSIM4v7sNode, BSIM4v7gNodePrime);
                TSTALLOC(BSIM4v7SbpPtr, BSIM4v7sNode, BSIM4v7bNodePrime);
            }
        }
    }

#ifdef USE_OMP
    InstCount = 0;
    model = (BSIM4v7model*)inModel;
    /* loop through all the BSIM4v7 device models 
       to count the number of instances */
    
    for( ; model != NULL; model = model->BSIM4v7nextModel )
    {
        /* loop through all the instances of the model */
        for (here = model->BSIM4v7instances; here != NULL ;
             here=here->BSIM4v7nextInstance) 
        { 
            InstCount++;
        }
    }
    InstArray = TMALLOC(BSIM4v7instance*, InstCount);
    model = (BSIM4v7model*)inModel;
    idx = 0;
    for( ; model != NULL; model = model->BSIM4v7nextModel )
    {
        /* loop through all the instances of the model */
        for (here = model->BSIM4v7instances; here != NULL ;
             here=here->BSIM4v7nextInstance) 
        { 
            InstArray[idx] = here;
            idx++;
        }
        /* set the array pointer and instance count into each model */
        model->BSIM4v7InstCount = InstCount;
        model->BSIM4v7InstanceArray = InstArray;		
    }
#endif

    return(OK);
}  

int
BSIM4v7unsetup(
GENmodel *inModel,
CKTcircuit *ckt)
{
#ifndef HAS_BATCHSIM
    BSIM4v7model *model;
    BSIM4v7instance *here;

    for (model = (BSIM4v7model *)inModel; model != NULL;
            model = model->BSIM4v7nextModel)
    {
        for (here = model->BSIM4v7instances; here != NULL;
                here=here->BSIM4v7nextInstance)
        {
            if (here->BSIM4v7dNodePrime
                    && here->BSIM4v7dNodePrime != here->BSIM4v7dNode)
            {
                CKTdltNNum(ckt, here->BSIM4v7dNodePrime);
                here->BSIM4v7dNodePrime = 0;
            }
            if (here->BSIM4v7sNodePrime
                    && here->BSIM4v7sNodePrime != here->BSIM4v7sNode)
            {
                CKTdltNNum(ckt, here->BSIM4v7sNodePrime);
                here->BSIM4v7sNodePrime = 0;
            }
        }
    }
#endif
    return OK;
}
