/**** BSIM4.8.0 Released by Navid Paydavosi 11/01/2013 ****/

/**********
 * Copyright 2006 Regents of the University of California. All rights reserved.
 * File: b4set.c of BSIM4.8.0.
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
 * Modified by Pankaj Kumar Thakur, 07//2012
 * Modified by Navid Paydavosi, 06/12/2013
**********/

#include "ngspice/ngspice.h"
#include "ngspice/jobdefs.h"
#include "ngspice/ftedefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
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
BSIM4setup(
SMPmatrix *matrix,
GENmodel *inModel,
CKTcircuit *ckt,
int *states)
{
BSIM4model *model = (BSIM4model*)inModel;
BSIM4instance *here;
int error;
CKTnode *tmp;
int    noiseAnalGiven = 0, createNode;  /* Criteria for new node creation */
double Rtot, DMCGeff, DMCIeff, DMDGeff;
JOB   *job;

#ifdef USE_OMP
unsigned int idx, InstCount;
BSIM4instance **InstArray;
#endif

    /* Search for a noise analysis request */
    for (job = ((TSKtask *)ft_curckt->ci_curTask)->jobs;job;job = job->JOBnextJob) {
        if(strcmp(job->JOBname,"Noise Analysis")==0) {
            noiseAnalGiven = 1;
            break;
        }
    }

    /*  loop through all the BSIM4 device models */
    for( ; model != NULL; model = model->BSIM4nextModel )
    {   /* process defaults of model parameters */
        if (!model->BSIM4typeGiven)
            model->BSIM4type = NMOS;     

        if (!model->BSIM4mobModGiven) 
            model->BSIM4mobMod = 0;
        else if ((model->BSIM4mobMod != 0) && (model->BSIM4mobMod != 1)&& (model->BSIM4mobMod != 2)&& (model->BSIM4mobMod != 3)
                 && (model->BSIM4mobMod != 4) && (model->BSIM4mobMod != 5) && (model->BSIM4mobMod != 6)) /* Synopsys 08/30/2013 modify */
        {   model->BSIM4mobMod = 0;
            printf("Warning: mobMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4binUnitGiven) 
            model->BSIM4binUnit = 1;
        if (!model->BSIM4paramChkGiven) 
            model->BSIM4paramChk = 1;

        if (!model->BSIM4dioModGiven)
            model->BSIM4dioMod = 1;
        else if ((model->BSIM4dioMod != 0) && (model->BSIM4dioMod != 1)
            && (model->BSIM4dioMod != 2))
        {   model->BSIM4dioMod = 1;
            printf("Warning: dioMod has been set to its default value: 1.\n");
        }

        if (!model->BSIM4cvchargeModGiven) 
            model->BSIM4cvchargeMod = 0;
        if (!model->BSIM4capModGiven) 
            model->BSIM4capMod = 2;
        else if ((model->BSIM4capMod != 0) && (model->BSIM4capMod != 1)
            && (model->BSIM4capMod != 2))
        {   model->BSIM4capMod = 2;
            printf("Warning: capMod has been set to its default value: 2.\n");
        }

        if (!model->BSIM4rdsModGiven)
            model->BSIM4rdsMod = 0;
        else if ((model->BSIM4rdsMod != 0) && (model->BSIM4rdsMod != 1))
        {   model->BSIM4rdsMod = 0;
            printf("Warning: rdsMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4rbodyModGiven)
            model->BSIM4rbodyMod = 0;
        else if ((model->BSIM4rbodyMod != 0) && (model->BSIM4rbodyMod != 1) && (model->BSIM4rbodyMod != 2))
        {   model->BSIM4rbodyMod = 0;
            printf("Warning: rbodyMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4rgateModGiven)
            model->BSIM4rgateMod = 0;
        else if ((model->BSIM4rgateMod != 0) && (model->BSIM4rgateMod != 1)
            && (model->BSIM4rgateMod != 2) && (model->BSIM4rgateMod != 3))
        {   model->BSIM4rgateMod = 0;
            printf("Warning: rgateMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4perModGiven)
            model->BSIM4perMod = 1;
        else if ((model->BSIM4perMod != 0) && (model->BSIM4perMod != 1))
        {   model->BSIM4perMod = 1;
            printf("Warning: perMod has been set to its default value: 1.\n");
        }

        if (!model->BSIM4geoModGiven)
            model->BSIM4geoMod = 0;

        if (!model->BSIM4rgeoModGiven)
            model->BSIM4rgeoMod = 0;
        else if ((model->BSIM4rgeoMod != 0) && (model->BSIM4rgeoMod != 1))
        {   model->BSIM4rgeoMod = 1;
            printf("Warning: rgeoMod has been set to its default value: 1.\n");
        }
        if (!model->BSIM4fnoiModGiven) 
            model->BSIM4fnoiMod = 1;
        else if ((model->BSIM4fnoiMod != 0) && (model->BSIM4fnoiMod != 1))
        {   model->BSIM4fnoiMod = 1;
            printf("Warning: fnoiMod has been set to its default value: 1.\n");
        }
        if (!model->BSIM4tnoiModGiven)
            model->BSIM4tnoiMod = 0; /* WDLiu: tnoiMod=1 needs to set internal S/D nodes */
        else if ((model->BSIM4tnoiMod != 0) && (model->BSIM4tnoiMod != 1) && (model->BSIM4tnoiMod != 2))  /* v4.7 */
        {   model->BSIM4tnoiMod = 0;
            printf("Warning: tnoiMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4trnqsModGiven)
            model->BSIM4trnqsMod = 0; 
        else if ((model->BSIM4trnqsMod != 0) && (model->BSIM4trnqsMod != 1))
        {   model->BSIM4trnqsMod = 0;
            printf("Warning: trnqsMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4acnqsModGiven)
            model->BSIM4acnqsMod = 0;
        else if ((model->BSIM4acnqsMod != 0) && (model->BSIM4acnqsMod != 1))
        {   model->BSIM4acnqsMod = 0;
            printf("Warning: acnqsMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4mtrlModGiven)   
            model->BSIM4mtrlMod = 0;
        else if((model->BSIM4mtrlMod != 0) && (model->BSIM4mtrlMod != 1))
        {
            model->BSIM4mtrlMod = 0;
            printf("Warning: mtrlMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4mtrlCompatModGiven)   
            model->BSIM4mtrlCompatMod = 0;
        else if((model->BSIM4mtrlCompatMod != 0) && (model->BSIM4mtrlCompatMod != 1))
        {
            model->BSIM4mtrlCompatMod = 0;
            printf("Warning: mtrlCompatMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4igcModGiven)
            model->BSIM4igcMod = 0;
        else if ((model->BSIM4igcMod != 0) && (model->BSIM4igcMod != 1)
                  && (model->BSIM4igcMod != 2))
        {   model->BSIM4igcMod = 0;
            printf("Warning: igcMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4igbModGiven)
            model->BSIM4igbMod = 0;
        else if ((model->BSIM4igbMod != 0) && (model->BSIM4igbMod != 1))
        {   model->BSIM4igbMod = 0;
            printf("Warning: igbMod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4tempModGiven)
            model->BSIM4tempMod = 0;
        else if ((model->BSIM4tempMod != 0) && (model->BSIM4tempMod != 1) 
                  && (model->BSIM4tempMod != 2) && (model->BSIM4tempMod != 3))
        {   model->BSIM4tempMod = 0;
            printf("Warning: tempMod has been set to its default value: 0.\n");
        }

        if (!model->BSIM4versionGiven) 
            model->BSIM4version = copy("4.8.0");
        if (!model->BSIM4toxrefGiven)
            model->BSIM4toxref = 30.0e-10;
        if (!model->BSIM4eotGiven)
            model->BSIM4eot = 15.0e-10;
        if (!model->BSIM4vddeotGiven)
            model->BSIM4vddeot = (model->BSIM4type == NMOS) ? 1.5 : -1.5;
        if (!model->BSIM4tempeotGiven)
            model->BSIM4tempeot = 300.15;
        if (!model->BSIM4leffeotGiven)
            model->BSIM4leffeot = 1;
        if (!model->BSIM4weffeotGiven)
            model->BSIM4weffeot = 10;        
        if (!model->BSIM4adosGiven)
            model->BSIM4ados = 1.0;
        if (!model->BSIM4bdosGiven)
            model->BSIM4bdos = 1.0;
        if (!model->BSIM4toxeGiven)
            model->BSIM4toxe = 30.0e-10;
        if (!model->BSIM4toxpGiven)
            model->BSIM4toxp = model->BSIM4toxe;
        if (!model->BSIM4toxmGiven)
            model->BSIM4toxm = model->BSIM4toxe;
        if (!model->BSIM4dtoxGiven)
            model->BSIM4dtox = 0.0;
        if (!model->BSIM4epsroxGiven)
            model->BSIM4epsrox = 3.9;

        if (!model->BSIM4cdscGiven)
            model->BSIM4cdsc = 2.4e-4;   /* unit Q/V/m^2  */
        if (!model->BSIM4cdscbGiven)
            model->BSIM4cdscb = 0.0;   /* unit Q/V/m^2  */    
            if (!model->BSIM4cdscdGiven)
            model->BSIM4cdscd = 0.0;   /* unit Q/V/m^2  */
        if (!model->BSIM4citGiven)
            model->BSIM4cit = 0.0;   /* unit Q/V/m^2  */
        if (!model->BSIM4nfactorGiven)
            model->BSIM4nfactor = 1.0;
        if (!model->BSIM4xjGiven)
            model->BSIM4xj = .15e-6;
        if (!model->BSIM4vsatGiven)
            model->BSIM4vsat = 8.0e4;    /* unit m/s */ 
        if (!model->BSIM4atGiven)
            model->BSIM4at = 3.3e4;    /* unit m/s */ 
        if (!model->BSIM4a0Given)
            model->BSIM4a0 = 1.0;  
        if (!model->BSIM4agsGiven)
            model->BSIM4ags = 0.0;
        if (!model->BSIM4a1Given)
            model->BSIM4a1 = 0.0;
        if (!model->BSIM4a2Given)
            model->BSIM4a2 = 1.0;
        if (!model->BSIM4ketaGiven)
            model->BSIM4keta = -0.047;    /* unit  / V */
        if (!model->BSIM4nsubGiven)
            model->BSIM4nsub = 6.0e16;   /* unit 1/cm3 */
        if (!model->BSIM4phigGiven)
            model->BSIM4phig = 4.05;  
        if (!model->BSIM4epsrgateGiven)
            model->BSIM4epsrgate = 11.7;  
        if (!model->BSIM4easubGiven)
            model->BSIM4easub = 4.05;  
        if (!model->BSIM4epsrsubGiven)
            model->BSIM4epsrsub = 11.7; 
        if (!model->BSIM4ni0subGiven)
            model->BSIM4ni0sub = 1.45e10;   /* unit 1/cm3 */
        if (!model->BSIM4bg0subGiven)
            model->BSIM4bg0sub =  1.16;     /* unit eV */
        if (!model->BSIM4tbgasubGiven)
            model->BSIM4tbgasub = 7.02e-4;  
        if (!model->BSIM4tbgbsubGiven)
            model->BSIM4tbgbsub = 1108.0;  
        if (!model->BSIM4ndepGiven)
            model->BSIM4ndep = 1.7e17;   /* unit 1/cm3 */
        if (!model->BSIM4nsdGiven)
            model->BSIM4nsd = 1.0e20;   /* unit 1/cm3 */
        if (!model->BSIM4phinGiven)
            model->BSIM4phin = 0.0; /* unit V */
        if (!model->BSIM4ngateGiven)
            model->BSIM4ngate = 0;   /* unit 1/cm3 */
        if (!model->BSIM4vbmGiven)
            model->BSIM4vbm = -3.0;
        if (!model->BSIM4xtGiven)
            model->BSIM4xt = 1.55e-7;
        if (!model->BSIM4kt1Given)
            model->BSIM4kt1 = -0.11;      /* unit V */
        if (!model->BSIM4kt1lGiven)
            model->BSIM4kt1l = 0.0;      /* unit V*m */
        if (!model->BSIM4kt2Given)
            model->BSIM4kt2 = 0.022;      /* No unit */
        if (!model->BSIM4k3Given)
            model->BSIM4k3 = 80.0;      
        if (!model->BSIM4k3bGiven)
            model->BSIM4k3b = 0.0;      
        if (!model->BSIM4w0Given)
            model->BSIM4w0 = 2.5e-6;    
        if (!model->BSIM4lpe0Given)
            model->BSIM4lpe0 = 1.74e-7;     
        if (!model->BSIM4lpebGiven)
            model->BSIM4lpeb = 0.0;
        if (!model->BSIM4dvtp0Given)
            model->BSIM4dvtp0 = 0.0;
        if (!model->BSIM4dvtp1Given)
            model->BSIM4dvtp1 = 0.0;
        if (!model->BSIM4dvtp2Given)        /* New DIBL/Rout */
            model->BSIM4dvtp2 = 0.0;
        if (!model->BSIM4dvtp3Given)
            model->BSIM4dvtp3 = 0.0;
        if (!model->BSIM4dvtp4Given)
            model->BSIM4dvtp4 = 0.0;
        if (!model->BSIM4dvtp5Given)
            model->BSIM4dvtp5 = 0.0;
        if (!model->BSIM4dvt0Given)
            model->BSIM4dvt0 = 2.2;    
        if (!model->BSIM4dvt1Given)
            model->BSIM4dvt1 = 0.53;      
        if (!model->BSIM4dvt2Given)
            model->BSIM4dvt2 = -0.032;   /* unit 1 / V */     

        if (!model->BSIM4dvt0wGiven)
            model->BSIM4dvt0w = 0.0;    
        if (!model->BSIM4dvt1wGiven)
            model->BSIM4dvt1w = 5.3e6;    
        if (!model->BSIM4dvt2wGiven)
            model->BSIM4dvt2w = -0.032;   

        if (!model->BSIM4droutGiven)
            model->BSIM4drout = 0.56;     
        if (!model->BSIM4dsubGiven)
            model->BSIM4dsub = model->BSIM4drout;     
        if (!model->BSIM4vth0Given)
            model->BSIM4vth0 = (model->BSIM4type == NMOS) ? 0.7 : -0.7;
        if (!model->BSIM4vfbGiven)
            model->BSIM4vfb = -1.0;
        if (!model->BSIM4euGiven)
            model->BSIM4eu = (model->BSIM4type == NMOS) ? 1.67 : 1.0;
        if (!model->BSIM4ucsGiven)
            model->BSIM4ucs = (model->BSIM4type == NMOS) ? 1.67 : 1.0;
        if (!model->BSIM4uaGiven)
            model->BSIM4ua = ((model->BSIM4mobMod == 2)) ? 1.0e-15 : 1.0e-9; /* unit m/V */
        if (!model->BSIM4ua1Given)
            model->BSIM4ua1 = 1.0e-9;      /* unit m/V */
        if (!model->BSIM4ubGiven)
            model->BSIM4ub = 1.0e-19;     /* unit (m/V)**2 */
        if (!model->BSIM4ub1Given)
            model->BSIM4ub1 = -1.0e-18;     /* unit (m/V)**2 */
        if (!model->BSIM4ucGiven)
            model->BSIM4uc = (model->BSIM4mobMod == 1) ? -0.0465 : -0.0465e-9;   
        if (!model->BSIM4uc1Given)
            model->BSIM4uc1 = (model->BSIM4mobMod == 1) ? -0.056 : -0.056e-9;   
        if (!model->BSIM4udGiven)
            model->BSIM4ud = 0.0;     /* unit m**(-2) */
        if (!model->BSIM4ud1Given)
            model->BSIM4ud1 = 0.0;     
        if (!model->BSIM4upGiven)
            model->BSIM4up = 0.0;     
        if (!model->BSIM4lpGiven)
            model->BSIM4lp = 1.0e-8;     
        if (!model->BSIM4u0Given)
            model->BSIM4u0 = (model->BSIM4type == NMOS) ? 0.067 : 0.025;
        if (!model->BSIM4uteGiven)
            model->BSIM4ute = -1.5; 
        if (!model->BSIM4ucsteGiven)
            model->BSIM4ucste = -4.775e-3;
        if (!model->BSIM4voffGiven)
            model->BSIM4voff = -0.08;
        if (!model->BSIM4vofflGiven)
            model->BSIM4voffl = 0.0;
        if (!model->BSIM4voffcvlGiven)
            model->BSIM4voffcvl = 0.0;
        if (!model->BSIM4minvGiven)
            model->BSIM4minv = 0.0;
        if (!model->BSIM4minvcvGiven)
            model->BSIM4minvcv = 0.0;
        if (!model->BSIM4fproutGiven)
            model->BSIM4fprout = 0.0;
        if (!model->BSIM4pditsGiven)
            model->BSIM4pdits = 0.0;
        if (!model->BSIM4pditsdGiven)
            model->BSIM4pditsd = 0.0;
        if (!model->BSIM4pditslGiven)
            model->BSIM4pditsl = 0.0;
        if (!model->BSIM4deltaGiven)  
           model->BSIM4delta = 0.01;
        if (!model->BSIM4rdswminGiven)
            model->BSIM4rdswmin = 0.0;
        if (!model->BSIM4rdwminGiven)
            model->BSIM4rdwmin = 0.0;
        if (!model->BSIM4rswminGiven)
            model->BSIM4rswmin = 0.0;
        if (!model->BSIM4rdswGiven)
            model->BSIM4rdsw = 200.0; /* in ohm*um */     
        if (!model->BSIM4rdwGiven)
            model->BSIM4rdw = 100.0;
        if (!model->BSIM4rswGiven)
            model->BSIM4rsw = 100.0;
        if (!model->BSIM4prwgGiven)
            model->BSIM4prwg = 1.0; /* in 1/V */
        if (!model->BSIM4prwbGiven)
            model->BSIM4prwb = 0.0;      
        if (!model->BSIM4prtGiven)
            model->BSIM4prt = 0.0;      
        if (!model->BSIM4eta0Given)
            model->BSIM4eta0 = 0.08;      /* no unit  */ 
        if (!model->BSIM4etabGiven)
            model->BSIM4etab = -0.07;      /* unit  1/V */ 
        if (!model->BSIM4pclmGiven)
            model->BSIM4pclm = 1.3;      /* no unit  */ 
        if (!model->BSIM4pdibl1Given)
            model->BSIM4pdibl1 = 0.39;    /* no unit  */
        if (!model->BSIM4pdibl2Given)
            model->BSIM4pdibl2 = 0.0086;    /* no unit  */ 
        if (!model->BSIM4pdiblbGiven)
            model->BSIM4pdiblb = 0.0;    /* 1/V  */ 
        if (!model->BSIM4pscbe1Given)
            model->BSIM4pscbe1 = 4.24e8;     
        if (!model->BSIM4pscbe2Given)
            model->BSIM4pscbe2 = 1.0e-5;    
        if (!model->BSIM4pvagGiven)
            model->BSIM4pvag = 0.0;     
        if (!model->BSIM4wrGiven)  
            model->BSIM4wr = 1.0;
        if (!model->BSIM4dwgGiven)  
            model->BSIM4dwg = 0.0;
        if (!model->BSIM4dwbGiven)  
            model->BSIM4dwb = 0.0;
        if (!model->BSIM4b0Given)
            model->BSIM4b0 = 0.0;
        if (!model->BSIM4b1Given)  
            model->BSIM4b1 = 0.0;
        if (!model->BSIM4alpha0Given)  
            model->BSIM4alpha0 = 0.0;
        if (!model->BSIM4alpha1Given)
            model->BSIM4alpha1 = 0.0;
        if (!model->BSIM4beta0Given)  
            model->BSIM4beta0 = 0.0;
        if (!model->BSIM4gidlModGiven)
            model->BSIM4gidlMod = 0;         /* v4.7 New GIDL/GISL */
        if (!model->BSIM4agidlGiven)
            model->BSIM4agidl = 0.0;
        if (!model->BSIM4bgidlGiven)
            model->BSIM4bgidl = 2.3e9; /* V/m */
        if (!model->BSIM4cgidlGiven)
            model->BSIM4cgidl = 0.5; /* V^3 */
        if (!model->BSIM4egidlGiven)
            model->BSIM4egidl = 0.8; /* V */
        if (!model->BSIM4rgidlGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4rgidl = 1.0;
        if (!model->BSIM4kgidlGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4kgidl = 0.0;
        if (!model->BSIM4fgidlGiven)          /* v4.7 New GIDL/GISL */
        /*model->BSIM4fgidl = 0.0;*/
        /* Default value of fgdil set to 1 in BSIM4.8.0*/
            model->BSIM4fgidl = 1.0;

        /*if (!model->BSIM4agislGiven)
        {
            if (model->BSIM4agidlGiven)
                model->BSIM4agisl = model->BSIM4agidl;
            else
                model->BSIM4agisl = 0.0;
        }*/

        /*Default value of agidl being 0, agisl set as follows */

        /*if (!model->BSIM4bgislGiven)
        {
            if (model->BSIM4bgidlGiven)
                model->BSIM4bgisl = model->BSIM4bgidl;
            else
                model->BSIM4bgisl = 2.3e9;
        }*/

        /*Default value of bgidl being 2.3e9, bgisl set as follows */

        /*if (!model->BSIM4cgislGiven)
        {
            if (model->BSIM4cgidlGiven)
                model->BSIM4cgisl = model->BSIM4cgidl;
            else
                model->BSIM4cgisl = 0.5;
        }*/

        /*Default value of cgidl being 0.5, cgisl set as follows */

        /*if (!model->BSIM4egislGiven)
        {
            if (model->BSIM4egidlGiven)
                model->BSIM4egisl = model->BSIM4egidl;
            else
                model->BSIM4egisl = 0.8;
        }*/

        /*Default value of agisl, bgisl, cgisl, egisl, rgisl, kgisl, and fgisl are set as follows */
        if (!model->BSIM4agislGiven)
           model->BSIM4agisl = model->BSIM4agidl;
        if (!model->BSIM4bgislGiven)
           model->BSIM4bgisl = model->BSIM4bgidl;
        if (!model->BSIM4cgislGiven)
           model->BSIM4cgisl = model->BSIM4cgidl;
        if (!model->BSIM4egislGiven)
           model->BSIM4egisl = model->BSIM4egidl;
        if (!model->BSIM4rgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4rgisl = model->BSIM4rgidl;
        if (!model->BSIM4kgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4kgisl = model->BSIM4kgidl;
        if (!model->BSIM4fgislGiven)          /* v4.7 New GIDL/GISL */
            model->BSIM4fgisl = model->BSIM4fgidl;
        if (!model->BSIM4aigcGiven)
            model->BSIM4aigc = (model->BSIM4type == NMOS) ? 1.36e-2 : 9.80e-3;
        if (!model->BSIM4bigcGiven)
            model->BSIM4bigc = (model->BSIM4type == NMOS) ? 1.71e-3 : 7.59e-4;
        if (!model->BSIM4cigcGiven)
            model->BSIM4cigc = (model->BSIM4type == NMOS) ? 0.075 : 0.03;
        if (model->BSIM4aigsdGiven)
        {
            model->BSIM4aigs = model->BSIM4aigd = model->BSIM4aigsd;
        }
        else
        {
            model->BSIM4aigsd = (model->BSIM4type == NMOS) ? 1.36e-2 : 9.80e-3;
            if (!model->BSIM4aigsGiven)
                model->BSIM4aigs = (model->BSIM4type == NMOS) ? 1.36e-2 : 9.80e-3;
            if (!model->BSIM4aigdGiven)
                model->BSIM4aigd = (model->BSIM4type == NMOS) ? 1.36e-2 : 9.80e-3;
        }
        if (model->BSIM4bigsdGiven)
        {
            model->BSIM4bigs = model->BSIM4bigd = model->BSIM4bigsd;
        }
        else
        {
            model->BSIM4bigsd = (model->BSIM4type == NMOS) ? 1.71e-3 : 7.59e-4;
            if (!model->BSIM4bigsGiven)
                model->BSIM4bigs = (model->BSIM4type == NMOS) ? 1.71e-3 : 7.59e-4;
            if (!model->BSIM4bigdGiven)
                model->BSIM4bigd = (model->BSIM4type == NMOS) ? 1.71e-3 : 7.59e-4;
        }
        if (model->BSIM4cigsdGiven)
        {
            model->BSIM4cigs = model->BSIM4cigd = model->BSIM4cigsd;
        }
        else
        {
             model->BSIM4cigsd = (model->BSIM4type == NMOS) ? 0.075 : 0.03;
           if (!model->BSIM4cigsGiven)
                model->BSIM4cigs = (model->BSIM4type == NMOS) ? 0.075 : 0.03;
            if (!model->BSIM4cigdGiven)
                model->BSIM4cigd = (model->BSIM4type == NMOS) ? 0.075 : 0.03;
        }
        if (!model->BSIM4aigbaccGiven)
            model->BSIM4aigbacc = 1.36e-2;
        if (!model->BSIM4bigbaccGiven)
            model->BSIM4bigbacc = 1.71e-3;
        if (!model->BSIM4cigbaccGiven)
            model->BSIM4cigbacc = 0.075;
        if (!model->BSIM4aigbinvGiven)
            model->BSIM4aigbinv = 1.11e-2;
        if (!model->BSIM4bigbinvGiven)
            model->BSIM4bigbinv = 9.49e-4;
        if (!model->BSIM4cigbinvGiven)
            model->BSIM4cigbinv = 0.006;
        if (!model->BSIM4nigcGiven)
            model->BSIM4nigc = 1.0;
        if (!model->BSIM4nigbinvGiven)
            model->BSIM4nigbinv = 3.0;
        if (!model->BSIM4nigbaccGiven)
            model->BSIM4nigbacc = 1.0;
        if (!model->BSIM4ntoxGiven)
            model->BSIM4ntox = 1.0;
        if (!model->BSIM4eigbinvGiven)
            model->BSIM4eigbinv = 1.1;
        if (!model->BSIM4pigcdGiven)
            model->BSIM4pigcd = 1.0;
        if (!model->BSIM4poxedgeGiven)
            model->BSIM4poxedge = 1.0;
        if (!model->BSIM4xrcrg1Given)
            model->BSIM4xrcrg1 = 12.0;
        if (!model->BSIM4xrcrg2Given)
            model->BSIM4xrcrg2 = 1.0;
        if (!model->BSIM4ijthsfwdGiven)
            model->BSIM4ijthsfwd = 0.1; /* unit A */
        if (!model->BSIM4ijthdfwdGiven)
            model->BSIM4ijthdfwd = model->BSIM4ijthsfwd;
        if (!model->BSIM4ijthsrevGiven)
            model->BSIM4ijthsrev = 0.1; /* unit A */
        if (!model->BSIM4ijthdrevGiven)
            model->BSIM4ijthdrev = model->BSIM4ijthsrev;
        if (!model->BSIM4tnoiaGiven)
            model->BSIM4tnoia = 1.5;
        if (!model->BSIM4tnoibGiven)
            model->BSIM4tnoib = 3.5;
        if (!model->BSIM4tnoicGiven)
            model->BSIM4tnoic = 0.0;
        if (!model->BSIM4rnoiaGiven)
            model->BSIM4rnoia = 0.577;
        if (!model->BSIM4rnoibGiven)
            model->BSIM4rnoib = 0.5164;
        if (!model->BSIM4rnoicGiven)
            model->BSIM4rnoic = 0.395;
        if (!model->BSIM4ntnoiGiven)
            model->BSIM4ntnoi = 1.0;
        if (!model->BSIM4lambdaGiven)
            model->BSIM4lambda = 0.0;
        if (!model->BSIM4vtlGiven)
            model->BSIM4vtl = 2.0e5;    /* unit m/s */ 
        if (!model->BSIM4xnGiven)
            model->BSIM4xn = 3.0;   
        if (!model->BSIM4lcGiven)
            model->BSIM4lc = 5.0e-9;   
        if (!model->BSIM4vfbsdoffGiven)  
            model->BSIM4vfbsdoff = 0.0;  /* unit v */  
        if (!model->BSIM4tvfbsdoffGiven)
            model->BSIM4tvfbsdoff = 0.0;  
        if (!model->BSIM4tvoffGiven)
            model->BSIM4tvoff = 0.0;  
        if (!model->BSIM4tnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4tnfactor = 0.0;  
        if (!model->BSIM4teta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4teta0 = 0.0;  
        if (!model->BSIM4tvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4tvoffcv = 0.0;  
 
        if (!model->BSIM4lintnoiGiven)
            model->BSIM4lintnoi = 0.0;  /* unit m */  

        if (!model->BSIM4xjbvsGiven)
            model->BSIM4xjbvs = 1.0; /* no unit */
        if (!model->BSIM4xjbvdGiven)
            model->BSIM4xjbvd = model->BSIM4xjbvs;
        if (!model->BSIM4bvsGiven)
            model->BSIM4bvs = 10.0; /* V */
        if (!model->BSIM4bvdGiven)
            model->BSIM4bvd = model->BSIM4bvs;

        if (!model->BSIM4gbminGiven)
            model->BSIM4gbmin = 1.0e-12; /* in mho */
        if (!model->BSIM4rbdbGiven)
            model->BSIM4rbdb = 50.0; /* in ohm */
        if (!model->BSIM4rbpbGiven)
            model->BSIM4rbpb = 50.0;
        if (!model->BSIM4rbsbGiven)
            model->BSIM4rbsb = 50.0;
        if (!model->BSIM4rbpsGiven)
            model->BSIM4rbps = 50.0;
        if (!model->BSIM4rbpdGiven)
            model->BSIM4rbpd = 50.0;

        if (!model->BSIM4rbps0Given)
            model->BSIM4rbps0 = 50.0;
        if (!model->BSIM4rbpslGiven)
            model->BSIM4rbpsl = 0.0;
        if (!model->BSIM4rbpswGiven)
            model->BSIM4rbpsw = 0.0;
        if (!model->BSIM4rbpsnfGiven)
            model->BSIM4rbpsnf = 0.0;

        if (!model->BSIM4rbpd0Given)
            model->BSIM4rbpd0 = 50.0;
        if (!model->BSIM4rbpdlGiven)
            model->BSIM4rbpdl = 0.0;
        if (!model->BSIM4rbpdwGiven)
            model->BSIM4rbpdw = 0.0;
        if (!model->BSIM4rbpdnfGiven)
            model->BSIM4rbpdnf = 0.0;

        if (!model->BSIM4rbpbx0Given)
            model->BSIM4rbpbx0 = 100.0;
        if (!model->BSIM4rbpbxlGiven)
            model->BSIM4rbpbxl = 0.0;
        if (!model->BSIM4rbpbxwGiven)
            model->BSIM4rbpbxw = 0.0;
        if (!model->BSIM4rbpbxnfGiven)
            model->BSIM4rbpbxnf = 0.0;
        if (!model->BSIM4rbpby0Given)
            model->BSIM4rbpby0 = 100.0;
        if (!model->BSIM4rbpbylGiven)
            model->BSIM4rbpbyl = 0.0;
        if (!model->BSIM4rbpbywGiven)
            model->BSIM4rbpbyw = 0.0;
        if (!model->BSIM4rbpbynfGiven)
            model->BSIM4rbpbynf = 0.0;


        if (!model->BSIM4rbsbx0Given)
            model->BSIM4rbsbx0 = 100.0;
        if (!model->BSIM4rbsby0Given)
            model->BSIM4rbsby0 = 100.0;
        if (!model->BSIM4rbdbx0Given)
            model->BSIM4rbdbx0 = 100.0;
        if (!model->BSIM4rbdby0Given)
            model->BSIM4rbdby0 = 100.0;


        if (!model->BSIM4rbsdbxlGiven)
            model->BSIM4rbsdbxl = 0.0;
        if (!model->BSIM4rbsdbxwGiven)
            model->BSIM4rbsdbxw = 0.0;
        if (!model->BSIM4rbsdbxnfGiven)
            model->BSIM4rbsdbxnf = 0.0;
        if (!model->BSIM4rbsdbylGiven)
            model->BSIM4rbsdbyl = 0.0;
        if (!model->BSIM4rbsdbywGiven)
            model->BSIM4rbsdbyw = 0.0;
        if (!model->BSIM4rbsdbynfGiven)
            model->BSIM4rbsdbynf = 0.0;

        if (!model->BSIM4cgslGiven)  
            model->BSIM4cgsl = 0.0;
        if (!model->BSIM4cgdlGiven)  
            model->BSIM4cgdl = 0.0;
        if (!model->BSIM4ckappasGiven)  
            model->BSIM4ckappas = 0.6;
        if (!model->BSIM4ckappadGiven)
            model->BSIM4ckappad = model->BSIM4ckappas;
        if (!model->BSIM4clcGiven)  
            model->BSIM4clc = 0.1e-6;
        if (!model->BSIM4cleGiven)  
            model->BSIM4cle = 0.6;
        if (!model->BSIM4vfbcvGiven)  
            model->BSIM4vfbcv = -1.0;
        if (!model->BSIM4acdeGiven)
            model->BSIM4acde = 1.0;
        if (!model->BSIM4moinGiven)
            model->BSIM4moin = 15.0;
        if (!model->BSIM4noffGiven)
            model->BSIM4noff = 1.0;
        if (!model->BSIM4voffcvGiven)
            model->BSIM4voffcv = 0.0;
        if (!model->BSIM4dmcgGiven)
            model->BSIM4dmcg = 0.0;
        if (!model->BSIM4dmciGiven)
            model->BSIM4dmci = model->BSIM4dmcg;
        if (!model->BSIM4dmdgGiven)
            model->BSIM4dmdg = 0.0;
        if (!model->BSIM4dmcgtGiven)
            model->BSIM4dmcgt = 0.0;
        if (!model->BSIM4xgwGiven)
            model->BSIM4xgw = 0.0;
        if (!model->BSIM4xglGiven)
            model->BSIM4xgl = 0.0;
        if (!model->BSIM4rshgGiven)
            model->BSIM4rshg = 0.1;
        if (!model->BSIM4ngconGiven)
            model->BSIM4ngcon = 1.0;
        if (!model->BSIM4tcjGiven)
            model->BSIM4tcj = 0.0;
        if (!model->BSIM4tpbGiven)
            model->BSIM4tpb = 0.0;
        if (!model->BSIM4tcjswGiven)
            model->BSIM4tcjsw = 0.0;
        if (!model->BSIM4tpbswGiven)
            model->BSIM4tpbsw = 0.0;
        if (!model->BSIM4tcjswgGiven)
            model->BSIM4tcjswg = 0.0;
        if (!model->BSIM4tpbswgGiven)
            model->BSIM4tpbswg = 0.0;

        /* Length dependence */
        if (!model->BSIM4lcdscGiven)
            model->BSIM4lcdsc = 0.0;
        if (!model->BSIM4lcdscbGiven)
            model->BSIM4lcdscb = 0.0;
            if (!model->BSIM4lcdscdGiven) 
            model->BSIM4lcdscd = 0.0;
        if (!model->BSIM4lcitGiven)
            model->BSIM4lcit = 0.0;
        if (!model->BSIM4lnfactorGiven)
            model->BSIM4lnfactor = 0.0;
        if (!model->BSIM4lxjGiven)
            model->BSIM4lxj = 0.0;
        if (!model->BSIM4lvsatGiven)
            model->BSIM4lvsat = 0.0;
        if (!model->BSIM4latGiven)
            model->BSIM4lat = 0.0;
        if (!model->BSIM4la0Given)
            model->BSIM4la0 = 0.0; 
        if (!model->BSIM4lagsGiven)
            model->BSIM4lags = 0.0;
        if (!model->BSIM4la1Given)
            model->BSIM4la1 = 0.0;
        if (!model->BSIM4la2Given)
            model->BSIM4la2 = 0.0;
        if (!model->BSIM4lketaGiven)
            model->BSIM4lketa = 0.0;
        if (!model->BSIM4lnsubGiven)
            model->BSIM4lnsub = 0.0;
        if (!model->BSIM4lndepGiven)
            model->BSIM4lndep = 0.0;
        if (!model->BSIM4lnsdGiven)
            model->BSIM4lnsd = 0.0;
        if (!model->BSIM4lphinGiven)
            model->BSIM4lphin = 0.0;
        if (!model->BSIM4lngateGiven)
            model->BSIM4lngate = 0.0;
        if (!model->BSIM4lvbmGiven)
            model->BSIM4lvbm = 0.0;
        if (!model->BSIM4lxtGiven)
            model->BSIM4lxt = 0.0;
        if (!model->BSIM4lk1Given)
            model->BSIM4lk1 = 0.0;
        if (!model->BSIM4lkt1Given)
            model->BSIM4lkt1 = 0.0; 
        if (!model->BSIM4lkt1lGiven)
            model->BSIM4lkt1l = 0.0;
        if (!model->BSIM4lkt2Given)
            model->BSIM4lkt2 = 0.0;
        if (!model->BSIM4lk2Given)
            model->BSIM4lk2 = 0.0;
        if (!model->BSIM4lk3Given)
            model->BSIM4lk3 = 0.0;      
        if (!model->BSIM4lk3bGiven)
            model->BSIM4lk3b = 0.0;      
        if (!model->BSIM4lw0Given)
            model->BSIM4lw0 = 0.0;    
        if (!model->BSIM4llpe0Given)
            model->BSIM4llpe0 = 0.0;
        if (!model->BSIM4llpebGiven)
            model->BSIM4llpeb = 0.0; 
        if (!model->BSIM4ldvtp0Given)
            model->BSIM4ldvtp0 = 0.0;
        if (!model->BSIM4ldvtp1Given)
            model->BSIM4ldvtp1 = 0.0;
        if (!model->BSIM4ldvtp2Given)        /* New DIBL/Rout */
            model->BSIM4ldvtp2 = 0.0;
        if (!model->BSIM4ldvtp3Given)
            model->BSIM4ldvtp3 = 0.0;
        if (!model->BSIM4ldvtp4Given)
            model->BSIM4ldvtp4 = 0.0;
        if (!model->BSIM4ldvtp5Given)
            model->BSIM4ldvtp5 = 0.0;
        if (!model->BSIM4ldvt0Given)
            model->BSIM4ldvt0 = 0.0;    
        if (!model->BSIM4ldvt1Given)
            model->BSIM4ldvt1 = 0.0;      
        if (!model->BSIM4ldvt2Given)
            model->BSIM4ldvt2 = 0.0;
        if (!model->BSIM4ldvt0wGiven)
            model->BSIM4ldvt0w = 0.0;    
        if (!model->BSIM4ldvt1wGiven)
            model->BSIM4ldvt1w = 0.0;      
        if (!model->BSIM4ldvt2wGiven)
            model->BSIM4ldvt2w = 0.0;
        if (!model->BSIM4ldroutGiven)
            model->BSIM4ldrout = 0.0;     
        if (!model->BSIM4ldsubGiven)
            model->BSIM4ldsub = 0.0;
        if (!model->BSIM4lvth0Given)
           model->BSIM4lvth0 = 0.0;
        if (!model->BSIM4luaGiven)
            model->BSIM4lua = 0.0;
        if (!model->BSIM4lua1Given)
            model->BSIM4lua1 = 0.0;
        if (!model->BSIM4lubGiven)
            model->BSIM4lub = 0.0;
        if (!model->BSIM4lub1Given)
            model->BSIM4lub1 = 0.0;
        if (!model->BSIM4lucGiven)
            model->BSIM4luc = 0.0;
        if (!model->BSIM4luc1Given)
            model->BSIM4luc1 = 0.0;
        if (!model->BSIM4ludGiven)
            model->BSIM4lud = 0.0;
        if (!model->BSIM4lud1Given)
            model->BSIM4lud1 = 0.0;
        if (!model->BSIM4lupGiven)
            model->BSIM4lup = 0.0;
        if (!model->BSIM4llpGiven)
            model->BSIM4llp = 0.0;
        if (!model->BSIM4lu0Given)
            model->BSIM4lu0 = 0.0;
        if (!model->BSIM4luteGiven)
            model->BSIM4lute = 0.0;  
          if (!model->BSIM4lucsteGiven)
            model->BSIM4lucste = 0.0;                 
        if (!model->BSIM4lvoffGiven)
            model->BSIM4lvoff = 0.0;
        if (!model->BSIM4lminvGiven)
            model->BSIM4lminv = 0.0;
        if (!model->BSIM4lminvcvGiven)
            model->BSIM4lminvcv = 0.0;
        if (!model->BSIM4lfproutGiven)
            model->BSIM4lfprout = 0.0;
        if (!model->BSIM4lpditsGiven)
            model->BSIM4lpdits = 0.0;
        if (!model->BSIM4lpditsdGiven)
            model->BSIM4lpditsd = 0.0;
        if (!model->BSIM4ldeltaGiven)  
            model->BSIM4ldelta = 0.0;
        if (!model->BSIM4lrdswGiven)
            model->BSIM4lrdsw = 0.0;
        if (!model->BSIM4lrdwGiven)
            model->BSIM4lrdw = 0.0;
        if (!model->BSIM4lrswGiven)
            model->BSIM4lrsw = 0.0;
        if (!model->BSIM4lprwbGiven)
            model->BSIM4lprwb = 0.0;
        if (!model->BSIM4lprwgGiven)
            model->BSIM4lprwg = 0.0;
        if (!model->BSIM4lprtGiven)
            model->BSIM4lprt = 0.0;
        if (!model->BSIM4leta0Given)
            model->BSIM4leta0 = 0.0;
        if (!model->BSIM4letabGiven)
            model->BSIM4letab = -0.0;
        if (!model->BSIM4lpclmGiven)
            model->BSIM4lpclm = 0.0; 
        if (!model->BSIM4lpdibl1Given)
            model->BSIM4lpdibl1 = 0.0;
        if (!model->BSIM4lpdibl2Given)
            model->BSIM4lpdibl2 = 0.0;
        if (!model->BSIM4lpdiblbGiven)
            model->BSIM4lpdiblb = 0.0;
        if (!model->BSIM4lpscbe1Given)
            model->BSIM4lpscbe1 = 0.0;
        if (!model->BSIM4lpscbe2Given)
            model->BSIM4lpscbe2 = 0.0;
        if (!model->BSIM4lpvagGiven)
            model->BSIM4lpvag = 0.0;     
        if (!model->BSIM4lwrGiven)  
            model->BSIM4lwr = 0.0;
        if (!model->BSIM4ldwgGiven)  
            model->BSIM4ldwg = 0.0;
        if (!model->BSIM4ldwbGiven)  
            model->BSIM4ldwb = 0.0;
        if (!model->BSIM4lb0Given)
            model->BSIM4lb0 = 0.0;
        if (!model->BSIM4lb1Given)  
            model->BSIM4lb1 = 0.0;
        if (!model->BSIM4lalpha0Given)  
            model->BSIM4lalpha0 = 0.0;
        if (!model->BSIM4lalpha1Given)
            model->BSIM4lalpha1 = 0.0;
        if (!model->BSIM4lbeta0Given)  
            model->BSIM4lbeta0 = 0.0;
        if (!model->BSIM4lagidlGiven)
            model->BSIM4lagidl = 0.0;
        if (!model->BSIM4lbgidlGiven)
            model->BSIM4lbgidl = 0.0;
        if (!model->BSIM4lcgidlGiven)
            model->BSIM4lcgidl = 0.0;
        if (!model->BSIM4legidlGiven)
            model->BSIM4legidl = 0.0;
        if (!model->BSIM4lrgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4lrgidl = 0.0;
        if (!model->BSIM4lkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4lkgidl = 0.0;
        if (!model->BSIM4lfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4lfgidl = 0.0;
        /*if (!model->BSIM4lagislGiven)
        {
            if (model->BSIM4lagidlGiven)
                model->BSIM4lagisl = model->BSIM4lagidl;
            else
                model->BSIM4lagisl = 0.0;
        }
        if (!model->BSIM4lbgislGiven)
        {
            if (model->BSIM4lbgidlGiven)
                model->BSIM4lbgisl = model->BSIM4lbgidl;
            else
                model->BSIM4lbgisl = 0.0;
        }
        if (!model->BSIM4lcgislGiven)
        {
            if (model->BSIM4lcgidlGiven)
                model->BSIM4lcgisl = model->BSIM4lcgidl;
            else
                model->BSIM4lcgisl = 0.0;
        }
        if (!model->BSIM4legislGiven)
        {
            if (model->BSIM4legidlGiven)
                model->BSIM4legisl = model->BSIM4legidl;
            else
                model->BSIM4legisl = 0.0; 
        }*/
        /*if (!model->BSIM4lrgislGiven)
        {
            if (model->BSIM4lrgidlGiven)
                model->BSIM4lrgisl = model->BSIM4lrgidl;
        }
        if (!model->BSIM4lkgislGiven)
        {
            if (model->BSIM4lkgidlGiven)
                model->BSIM4lkgisl = model->BSIM4lkgidl;
        }
        if (!model->BSIM4lfgislGiven)
        {
            if (model->BSIM4lfgidlGiven)
                model->BSIM4lfgisl = model->BSIM4lfgidl;
        }*/

        /*Default value of lagisl, lbgisl, lcgisl, legisl, lrgisl, lkgisl, and lfgisl are set as follows */
        if (!model->BSIM4lagislGiven)
           model->BSIM4lagisl = model->BSIM4lagidl;
        if (!model->BSIM4lbgislGiven)
           model->BSIM4lbgisl = model->BSIM4lbgidl;
        if (!model->BSIM4lcgislGiven)
           model->BSIM4lcgisl = model->BSIM4lcgidl;
        if (!model->BSIM4legislGiven)
           model->BSIM4legisl = model->BSIM4legidl;
        if (!model->BSIM4lrgislGiven) 	/* v4.7 New GIDL/GISL */
           model->BSIM4lrgisl = model->BSIM4lrgidl;
        if (!model->BSIM4lkgislGiven)	/* v4.7 New GIDL/GISL */
            model->BSIM4lkgisl = model->BSIM4lkgidl;
        if (!model->BSIM4lfgislGiven)	/* v4.7 New GIDL/GISL */
            model->BSIM4lfgisl = model->BSIM4lfgidl;

        if (!model->BSIM4laigcGiven)
            model->BSIM4laigc = 0.0;
        if (!model->BSIM4lbigcGiven)
            model->BSIM4lbigc = 0.0;
        if (!model->BSIM4lcigcGiven)
            model->BSIM4lcigc = 0.0;
        if (!model->BSIM4aigsdGiven && (model->BSIM4aigsGiven || model->BSIM4aigdGiven))
        {
            if (!model->BSIM4laigsGiven)
                model->BSIM4laigs = 0.0;
            if (!model->BSIM4laigdGiven)
                model->BSIM4laigd = 0.0;
        }
        else
        {
           if (!model->BSIM4laigsdGiven)
               model->BSIM4laigsd = 0.0;
           model->BSIM4laigs = model->BSIM4laigd = model->BSIM4laigsd;
        }
        if (!model->BSIM4bigsdGiven && (model->BSIM4bigsGiven || model->BSIM4bigdGiven))
        {
            if (!model->BSIM4lbigsGiven)
                model->BSIM4lbigs = 0.0;
            if (!model->BSIM4lbigdGiven)
                model->BSIM4lbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4lbigsdGiven)
               model->BSIM4lbigsd = 0.0;
           model->BSIM4lbigs = model->BSIM4lbigd = model->BSIM4lbigsd;
        }
        if (!model->BSIM4cigsdGiven && (model->BSIM4cigsGiven || model->BSIM4cigdGiven))
        {
            if (!model->BSIM4lcigsGiven)
                model->BSIM4lcigs = 0.0;
            if (!model->BSIM4lcigdGiven)
                model->BSIM4lcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4lcigsdGiven)
               model->BSIM4lcigsd = 0.0;
           model->BSIM4lcigs = model->BSIM4lcigd = model->BSIM4lcigsd;
        }
        if (!model->BSIM4laigbaccGiven)
            model->BSIM4laigbacc = 0.0;
        if (!model->BSIM4lbigbaccGiven)
            model->BSIM4lbigbacc = 0.0;
        if (!model->BSIM4lcigbaccGiven)
            model->BSIM4lcigbacc = 0.0;
        if (!model->BSIM4laigbinvGiven)
            model->BSIM4laigbinv = 0.0;
        if (!model->BSIM4lbigbinvGiven)
            model->BSIM4lbigbinv = 0.0;
        if (!model->BSIM4lcigbinvGiven)
            model->BSIM4lcigbinv = 0.0;
        if (!model->BSIM4lnigcGiven)
            model->BSIM4lnigc = 0.0;
        if (!model->BSIM4lnigbinvGiven)
            model->BSIM4lnigbinv = 0.0;
        if (!model->BSIM4lnigbaccGiven)
            model->BSIM4lnigbacc = 0.0;
        if (!model->BSIM4lntoxGiven)
            model->BSIM4lntox = 0.0;
        if (!model->BSIM4leigbinvGiven)
            model->BSIM4leigbinv = 0.0;
        if (!model->BSIM4lpigcdGiven)
            model->BSIM4lpigcd = 0.0;
        if (!model->BSIM4lpoxedgeGiven)
            model->BSIM4lpoxedge = 0.0;
        if (!model->BSIM4lxrcrg1Given)
            model->BSIM4lxrcrg1 = 0.0;
        if (!model->BSIM4lxrcrg2Given)
            model->BSIM4lxrcrg2 = 0.0;
        if (!model->BSIM4leuGiven)
            model->BSIM4leu = 0.0;
                if (!model->BSIM4lucsGiven)
            model->BSIM4lucs = 0.0;
        if (!model->BSIM4lvfbGiven)
            model->BSIM4lvfb = 0.0;
        if (!model->BSIM4llambdaGiven)
            model->BSIM4llambda = 0.0;
        if (!model->BSIM4lvtlGiven)
            model->BSIM4lvtl = 0.0;  
        if (!model->BSIM4lxnGiven)
            model->BSIM4lxn = 0.0;  
        if (!model->BSIM4lvfbsdoffGiven)
            model->BSIM4lvfbsdoff = 0.0;   
        if (!model->BSIM4ltvfbsdoffGiven)
            model->BSIM4ltvfbsdoff = 0.0;  
        if (!model->BSIM4ltvoffGiven)
            model->BSIM4ltvoff = 0.0;  
        if (!model->BSIM4ltnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4ltnfactor = 0.0;  
        if (!model->BSIM4lteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4lteta0 = 0.0;  
        if (!model->BSIM4ltvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4ltvoffcv = 0.0;  


        if (!model->BSIM4lcgslGiven)  
            model->BSIM4lcgsl = 0.0;
        if (!model->BSIM4lcgdlGiven)  
            model->BSIM4lcgdl = 0.0;
        if (!model->BSIM4lckappasGiven)  
            model->BSIM4lckappas = 0.0;
        if (!model->BSIM4lckappadGiven)
            model->BSIM4lckappad = 0.0;
        if (!model->BSIM4lclcGiven)  
            model->BSIM4lclc = 0.0;
        if (!model->BSIM4lcleGiven)  
            model->BSIM4lcle = 0.0;
        if (!model->BSIM4lcfGiven)  
            model->BSIM4lcf = 0.0;
        if (!model->BSIM4lvfbcvGiven)  
            model->BSIM4lvfbcv = 0.0;
        if (!model->BSIM4lacdeGiven)
            model->BSIM4lacde = 0.0;
        if (!model->BSIM4lmoinGiven)
            model->BSIM4lmoin = 0.0;
        if (!model->BSIM4lnoffGiven)
            model->BSIM4lnoff = 0.0;
        if (!model->BSIM4lvoffcvGiven)
            model->BSIM4lvoffcv = 0.0;

        /* Width dependence */
        if (!model->BSIM4wcdscGiven)
            model->BSIM4wcdsc = 0.0;
        if (!model->BSIM4wcdscbGiven)
            model->BSIM4wcdscb = 0.0;  
            if (!model->BSIM4wcdscdGiven)
            model->BSIM4wcdscd = 0.0;
        if (!model->BSIM4wcitGiven)
            model->BSIM4wcit = 0.0;
        if (!model->BSIM4wnfactorGiven)
            model->BSIM4wnfactor = 0.0;
        if (!model->BSIM4wxjGiven)
            model->BSIM4wxj = 0.0;
        if (!model->BSIM4wvsatGiven)
            model->BSIM4wvsat = 0.0;
        if (!model->BSIM4watGiven)
            model->BSIM4wat = 0.0;
        if (!model->BSIM4wa0Given)
            model->BSIM4wa0 = 0.0; 
        if (!model->BSIM4wagsGiven)
            model->BSIM4wags = 0.0;
        if (!model->BSIM4wa1Given)
            model->BSIM4wa1 = 0.0;
        if (!model->BSIM4wa2Given)
            model->BSIM4wa2 = 0.0;
        if (!model->BSIM4wketaGiven)
            model->BSIM4wketa = 0.0;
        if (!model->BSIM4wnsubGiven)
            model->BSIM4wnsub = 0.0;
        if (!model->BSIM4wndepGiven)
            model->BSIM4wndep = 0.0;
        if (!model->BSIM4wnsdGiven)
            model->BSIM4wnsd = 0.0;
        if (!model->BSIM4wphinGiven)
            model->BSIM4wphin = 0.0;
        if (!model->BSIM4wngateGiven)
            model->BSIM4wngate = 0.0;
        if (!model->BSIM4wvbmGiven)
            model->BSIM4wvbm = 0.0;
        if (!model->BSIM4wxtGiven)
            model->BSIM4wxt = 0.0;
        if (!model->BSIM4wk1Given)
            model->BSIM4wk1 = 0.0;
        if (!model->BSIM4wkt1Given)
            model->BSIM4wkt1 = 0.0; 
        if (!model->BSIM4wkt1lGiven)
            model->BSIM4wkt1l = 0.0;
        if (!model->BSIM4wkt2Given)
            model->BSIM4wkt2 = 0.0;
        if (!model->BSIM4wk2Given)
            model->BSIM4wk2 = 0.0;
        if (!model->BSIM4wk3Given)
            model->BSIM4wk3 = 0.0;      
        if (!model->BSIM4wk3bGiven)
            model->BSIM4wk3b = 0.0;      
        if (!model->BSIM4ww0Given)
            model->BSIM4ww0 = 0.0;    
        if (!model->BSIM4wlpe0Given)
            model->BSIM4wlpe0 = 0.0;
        if (!model->BSIM4wlpebGiven)
            model->BSIM4wlpeb = 0.0; 
        if (!model->BSIM4wdvtp0Given)
            model->BSIM4wdvtp0 = 0.0;
        if (!model->BSIM4wdvtp1Given)
            model->BSIM4wdvtp1 = 0.0;
        if (!model->BSIM4wdvtp2Given)        /* New DIBL/Rout */
            model->BSIM4wdvtp2 = 0.0;
        if (!model->BSIM4wdvtp3Given)
            model->BSIM4wdvtp3 = 0.0;
        if (!model->BSIM4wdvtp4Given)
            model->BSIM4wdvtp4 = 0.0;
        if (!model->BSIM4wdvtp5Given)
            model->BSIM4wdvtp5 = 0.0;
        if (!model->BSIM4wdvt0Given)
            model->BSIM4wdvt0 = 0.0;    
        if (!model->BSIM4wdvt1Given)
            model->BSIM4wdvt1 = 0.0;      
        if (!model->BSIM4wdvt2Given)
            model->BSIM4wdvt2 = 0.0;
        if (!model->BSIM4wdvt0wGiven)
            model->BSIM4wdvt0w = 0.0;    
        if (!model->BSIM4wdvt1wGiven)
            model->BSIM4wdvt1w = 0.0;      
        if (!model->BSIM4wdvt2wGiven)
            model->BSIM4wdvt2w = 0.0;
        if (!model->BSIM4wdroutGiven)
            model->BSIM4wdrout = 0.0;     
        if (!model->BSIM4wdsubGiven)
            model->BSIM4wdsub = 0.0;
        if (!model->BSIM4wvth0Given)
           model->BSIM4wvth0 = 0.0;
        if (!model->BSIM4wuaGiven)
            model->BSIM4wua = 0.0;
        if (!model->BSIM4wua1Given)
            model->BSIM4wua1 = 0.0;
        if (!model->BSIM4wubGiven)
            model->BSIM4wub = 0.0;
        if (!model->BSIM4wub1Given)
            model->BSIM4wub1 = 0.0;
        if (!model->BSIM4wucGiven)
            model->BSIM4wuc = 0.0;
        if (!model->BSIM4wuc1Given)
            model->BSIM4wuc1 = 0.0;
        if (!model->BSIM4wudGiven)
            model->BSIM4wud = 0.0;
        if (!model->BSIM4wud1Given)
            model->BSIM4wud1 = 0.0;
        if (!model->BSIM4wupGiven)
            model->BSIM4wup = 0.0;
        if (!model->BSIM4wlpGiven)
            model->BSIM4wlp = 0.0;
        if (!model->BSIM4wu0Given)
            model->BSIM4wu0 = 0.0;
        if (!model->BSIM4wuteGiven)
                model->BSIM4wute = 0.0; 
        if (!model->BSIM4wucsteGiven)
                model->BSIM4wucste = 0.0;                 
        if (!model->BSIM4wvoffGiven)
                model->BSIM4wvoff = 0.0;
        if (!model->BSIM4wminvGiven)
            model->BSIM4wminv = 0.0;
        if (!model->BSIM4wminvcvGiven)
            model->BSIM4wminvcv = 0.0;
        if (!model->BSIM4wfproutGiven)
            model->BSIM4wfprout = 0.0;
        if (!model->BSIM4wpditsGiven)
            model->BSIM4wpdits = 0.0;
        if (!model->BSIM4wpditsdGiven)
            model->BSIM4wpditsd = 0.0;
        if (!model->BSIM4wdeltaGiven)  
            model->BSIM4wdelta = 0.0;
        if (!model->BSIM4wrdswGiven)
            model->BSIM4wrdsw = 0.0;
        if (!model->BSIM4wrdwGiven)
            model->BSIM4wrdw = 0.0;
        if (!model->BSIM4wrswGiven)
            model->BSIM4wrsw = 0.0;
        if (!model->BSIM4wprwbGiven)
            model->BSIM4wprwb = 0.0;
        if (!model->BSIM4wprwgGiven)
            model->BSIM4wprwg = 0.0;
        if (!model->BSIM4wprtGiven)
            model->BSIM4wprt = 0.0;
        if (!model->BSIM4weta0Given)
            model->BSIM4weta0 = 0.0;
        if (!model->BSIM4wetabGiven)
            model->BSIM4wetab = 0.0;
        if (!model->BSIM4wpclmGiven)
            model->BSIM4wpclm = 0.0; 
        if (!model->BSIM4wpdibl1Given)
            model->BSIM4wpdibl1 = 0.0;
        if (!model->BSIM4wpdibl2Given)
            model->BSIM4wpdibl2 = 0.0;
        if (!model->BSIM4wpdiblbGiven)
            model->BSIM4wpdiblb = 0.0;
        if (!model->BSIM4wpscbe1Given)
            model->BSIM4wpscbe1 = 0.0;
        if (!model->BSIM4wpscbe2Given)
            model->BSIM4wpscbe2 = 0.0;
        if (!model->BSIM4wpvagGiven)
            model->BSIM4wpvag = 0.0;     
        if (!model->BSIM4wwrGiven)  
            model->BSIM4wwr = 0.0;
        if (!model->BSIM4wdwgGiven)  
            model->BSIM4wdwg = 0.0;
        if (!model->BSIM4wdwbGiven)  
            model->BSIM4wdwb = 0.0;
        if (!model->BSIM4wb0Given)
            model->BSIM4wb0 = 0.0;
        if (!model->BSIM4wb1Given)  
            model->BSIM4wb1 = 0.0;
        if (!model->BSIM4walpha0Given)  
            model->BSIM4walpha0 = 0.0;
        if (!model->BSIM4walpha1Given)
            model->BSIM4walpha1 = 0.0;
        if (!model->BSIM4wbeta0Given)  
            model->BSIM4wbeta0 = 0.0;
        if (!model->BSIM4wagidlGiven)
            model->BSIM4wagidl = 0.0;
        if (!model->BSIM4wbgidlGiven)
            model->BSIM4wbgidl = 0.0;
        if (!model->BSIM4wcgidlGiven)
            model->BSIM4wcgidl = 0.0;
        if (!model->BSIM4wegidlGiven)
            model->BSIM4wegidl = 0.0;
        if (!model->BSIM4wrgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4wrgidl = 0.0;
        if (!model->BSIM4wkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4wkgidl = 0.0;
        if (!model->BSIM4wfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4wfgidl = 0.0;

        /*if (!model->BSIM4wagislGiven)
        {
            if (model->BSIM4wagidlGiven)
                model->BSIM4wagisl = model->BSIM4wagidl;
            else
                model->BSIM4wagisl = 0.0;
        }
        if (!model->BSIM4wbgislGiven)
        {
            if (model->BSIM4wbgidlGiven)
                model->BSIM4wbgisl = model->BSIM4wbgidl;
            else
                model->BSIM4wbgisl = 0.0;
        }
        if (!model->BSIM4wcgislGiven)
        {
            if (model->BSIM4wcgidlGiven)
                model->BSIM4wcgisl = model->BSIM4wcgidl;
            else
                model->BSIM4wcgisl = 0.0;
        }
        if (!model->BSIM4wegislGiven)
        {
            if (model->BSIM4wegidlGiven)
                model->BSIM4wegisl = model->BSIM4wegidl;
            else
                model->BSIM4wegisl = 0.0; 
        }*/
        /*if (!model->BSIM4wrgislGiven)
        {
            if (model->BSIM4wrgidlGiven)
                model->BSIM4wrgisl = model->BSIM4wrgidl;
        }
        if (!model->BSIM4wkgislGiven)
        {
            if (model->BSIM4wkgidlGiven)
                model->BSIM4wkgisl = model->BSIM4wkgidl;
        }
        if (!model->BSIM4wfgislGiven)
        {
            if (model->BSIM4wfgidlGiven)
                model->BSIM4wfgisl = model->BSIM4wfgidl;
        }*/

        /*Default value of wagisl, wbgisl, wcgisl, wegisl, wrgisl, wkgisl, and wfgisl are set as follows */
        if (!model->BSIM4wagislGiven)
           model->BSIM4wagisl = model->BSIM4wagidl;
        if (!model->BSIM4wbgislGiven)
           model->BSIM4wbgisl = model->BSIM4wbgidl;
        if (!model->BSIM4wcgislGiven)
           model->BSIM4wcgisl = model->BSIM4wcgidl;
        if (!model->BSIM4wegislGiven)
           model->BSIM4wegisl = model->BSIM4wegidl;
        if (!model->BSIM4wrgislGiven) 	/* v4.7 New GIDL/GISL */
            model->BSIM4wrgisl = model->BSIM4wrgidl;
        if (!model->BSIM4wkgislGiven)	/* v4.7 New GIDL/GISL */
            model->BSIM4wkgisl = model->BSIM4wkgidl;
        if (!model->BSIM4wfgislGiven)	/* v4.7 New GIDL/GISL */
            model->BSIM4wfgisl = model->BSIM4wfgidl;

        if (!model->BSIM4waigcGiven)
            model->BSIM4waigc = 0.0;
        if (!model->BSIM4wbigcGiven)
            model->BSIM4wbigc = 0.0;
        if (!model->BSIM4wcigcGiven)
            model->BSIM4wcigc = 0.0;
        if (!model->BSIM4aigsdGiven && (model->BSIM4aigsGiven || model->BSIM4aigdGiven))
        {
            if (!model->BSIM4waigsGiven)
                model->BSIM4waigs = 0.0;
            if (!model->BSIM4waigdGiven)
                model->BSIM4waigd = 0.0;
        }
        else
        {
           if (!model->BSIM4waigsdGiven)
               model->BSIM4waigsd = 0.0;
           model->BSIM4waigs = model->BSIM4waigd = model->BSIM4waigsd;
        }
        if (!model->BSIM4bigsdGiven && (model->BSIM4bigsGiven || model->BSIM4bigdGiven))
        {
            if (!model->BSIM4wbigsGiven)
                model->BSIM4wbigs = 0.0;
            if (!model->BSIM4wbigdGiven)
                model->BSIM4wbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4wbigsdGiven)
               model->BSIM4wbigsd = 0.0;
           model->BSIM4wbigs = model->BSIM4wbigd = model->BSIM4wbigsd;
        }
        if (!model->BSIM4cigsdGiven && (model->BSIM4cigsGiven || model->BSIM4cigdGiven))
        {
            if (!model->BSIM4wcigsGiven)
                model->BSIM4wcigs = 0.0;
            if (!model->BSIM4wcigdGiven)
                model->BSIM4wcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4wcigsdGiven)
               model->BSIM4wcigsd = 0.0;
           model->BSIM4wcigs = model->BSIM4wcigd = model->BSIM4wcigsd;
        }
        if (!model->BSIM4waigbaccGiven)
            model->BSIM4waigbacc = 0.0;
        if (!model->BSIM4wbigbaccGiven)
            model->BSIM4wbigbacc = 0.0;
        if (!model->BSIM4wcigbaccGiven)
            model->BSIM4wcigbacc = 0.0;
        if (!model->BSIM4waigbinvGiven)
            model->BSIM4waigbinv = 0.0;
        if (!model->BSIM4wbigbinvGiven)
            model->BSIM4wbigbinv = 0.0;
        if (!model->BSIM4wcigbinvGiven)
            model->BSIM4wcigbinv = 0.0;
        if (!model->BSIM4wnigcGiven)
            model->BSIM4wnigc = 0.0;
        if (!model->BSIM4wnigbinvGiven)
            model->BSIM4wnigbinv = 0.0;
        if (!model->BSIM4wnigbaccGiven)
            model->BSIM4wnigbacc = 0.0;
        if (!model->BSIM4wntoxGiven)
            model->BSIM4wntox = 0.0;
        if (!model->BSIM4weigbinvGiven)
            model->BSIM4weigbinv = 0.0;
        if (!model->BSIM4wpigcdGiven)
            model->BSIM4wpigcd = 0.0;
        if (!model->BSIM4wpoxedgeGiven)
            model->BSIM4wpoxedge = 0.0;
        if (!model->BSIM4wxrcrg1Given)
            model->BSIM4wxrcrg1 = 0.0;
        if (!model->BSIM4wxrcrg2Given)
            model->BSIM4wxrcrg2 = 0.0;
        if (!model->BSIM4weuGiven)
            model->BSIM4weu = 0.0;
            if (!model->BSIM4wucsGiven)
            model->BSIM4wucs = 0.0;
        if (!model->BSIM4wvfbGiven)
            model->BSIM4wvfb = 0.0;
        if (!model->BSIM4wlambdaGiven)
            model->BSIM4wlambda = 0.0;
        if (!model->BSIM4wvtlGiven)
            model->BSIM4wvtl = 0.0;  
        if (!model->BSIM4wxnGiven)
            model->BSIM4wxn = 0.0;  
        if (!model->BSIM4wvfbsdoffGiven)
            model->BSIM4wvfbsdoff = 0.0;   
        if (!model->BSIM4wtvfbsdoffGiven)
            model->BSIM4wtvfbsdoff = 0.0;  
        if (!model->BSIM4wtvoffGiven)
            model->BSIM4wtvoff = 0.0;  
        if (!model->BSIM4wtnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4wtnfactor = 0.0;  
        if (!model->BSIM4wteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4wteta0 = 0.0;  
        if (!model->BSIM4wtvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4wtvoffcv = 0.0;  

        if (!model->BSIM4wcgslGiven)  
            model->BSIM4wcgsl = 0.0;
        if (!model->BSIM4wcgdlGiven)  
            model->BSIM4wcgdl = 0.0;
        if (!model->BSIM4wckappasGiven)  
            model->BSIM4wckappas = 0.0;
        if (!model->BSIM4wckappadGiven)
            model->BSIM4wckappad = 0.0;
        if (!model->BSIM4wcfGiven)  
            model->BSIM4wcf = 0.0;
        if (!model->BSIM4wclcGiven)  
            model->BSIM4wclc = 0.0;
        if (!model->BSIM4wcleGiven)  
            model->BSIM4wcle = 0.0;
        if (!model->BSIM4wvfbcvGiven)  
            model->BSIM4wvfbcv = 0.0;
        if (!model->BSIM4wacdeGiven)
            model->BSIM4wacde = 0.0;
        if (!model->BSIM4wmoinGiven)
            model->BSIM4wmoin = 0.0;
        if (!model->BSIM4wnoffGiven)
            model->BSIM4wnoff = 0.0;
        if (!model->BSIM4wvoffcvGiven)
            model->BSIM4wvoffcv = 0.0;

        /* Cross-term dependence */
        if (!model->BSIM4pcdscGiven)
            model->BSIM4pcdsc = 0.0;
        if (!model->BSIM4pcdscbGiven)
            model->BSIM4pcdscb = 0.0;   
            if (!model->BSIM4pcdscdGiven)
            model->BSIM4pcdscd = 0.0;
        if (!model->BSIM4pcitGiven)
            model->BSIM4pcit = 0.0;
        if (!model->BSIM4pnfactorGiven)
            model->BSIM4pnfactor = 0.0;
        if (!model->BSIM4pxjGiven)
            model->BSIM4pxj = 0.0;
        if (!model->BSIM4pvsatGiven)
            model->BSIM4pvsat = 0.0;
        if (!model->BSIM4patGiven)
            model->BSIM4pat = 0.0;
        if (!model->BSIM4pa0Given)
            model->BSIM4pa0 = 0.0; 
            
        if (!model->BSIM4pagsGiven)
            model->BSIM4pags = 0.0;
        if (!model->BSIM4pa1Given)
            model->BSIM4pa1 = 0.0;
        if (!model->BSIM4pa2Given)
            model->BSIM4pa2 = 0.0;
        if (!model->BSIM4pketaGiven)
            model->BSIM4pketa = 0.0;
        if (!model->BSIM4pnsubGiven)
            model->BSIM4pnsub = 0.0;
        if (!model->BSIM4pndepGiven)
            model->BSIM4pndep = 0.0;
        if (!model->BSIM4pnsdGiven)
            model->BSIM4pnsd = 0.0;
        if (!model->BSIM4pphinGiven)
            model->BSIM4pphin = 0.0;
        if (!model->BSIM4pngateGiven)
            model->BSIM4pngate = 0.0;
        if (!model->BSIM4pvbmGiven)
            model->BSIM4pvbm = 0.0;
        if (!model->BSIM4pxtGiven)
            model->BSIM4pxt = 0.0;
        if (!model->BSIM4pk1Given)
            model->BSIM4pk1 = 0.0;
        if (!model->BSIM4pkt1Given)
            model->BSIM4pkt1 = 0.0; 
        if (!model->BSIM4pkt1lGiven)
            model->BSIM4pkt1l = 0.0;
        if (!model->BSIM4pk2Given)
            model->BSIM4pk2 = 0.0;
        if (!model->BSIM4pkt2Given)
            model->BSIM4pkt2 = 0.0;
        if (!model->BSIM4pk3Given)
            model->BSIM4pk3 = 0.0;      
        if (!model->BSIM4pk3bGiven)
            model->BSIM4pk3b = 0.0;      
        if (!model->BSIM4pw0Given)
            model->BSIM4pw0 = 0.0;    
        if (!model->BSIM4plpe0Given)
            model->BSIM4plpe0 = 0.0;
        if (!model->BSIM4plpebGiven)
            model->BSIM4plpeb = 0.0;
        if (!model->BSIM4pdvtp0Given)
            model->BSIM4pdvtp0 = 0.0;
        if (!model->BSIM4pdvtp1Given)
            model->BSIM4pdvtp1 = 0.0;
        if (!model->BSIM4pdvtp2Given)        /* New DIBL/Rout */
            model->BSIM4pdvtp2 = 0.0;
        if (!model->BSIM4pdvtp3Given)
            model->BSIM4pdvtp3 = 0.0;
        if (!model->BSIM4pdvtp4Given)
            model->BSIM4pdvtp4 = 0.0;
        if (!model->BSIM4pdvtp5Given)
            model->BSIM4pdvtp5 = 0.0;
        if (!model->BSIM4pdvt0Given)
            model->BSIM4pdvt0 = 0.0;    
        if (!model->BSIM4pdvt1Given)
            model->BSIM4pdvt1 = 0.0;      
        if (!model->BSIM4pdvt2Given)
            model->BSIM4pdvt2 = 0.0;
        if (!model->BSIM4pdvt0wGiven)
            model->BSIM4pdvt0w = 0.0;    
        if (!model->BSIM4pdvt1wGiven)
            model->BSIM4pdvt1w = 0.0;      
        if (!model->BSIM4pdvt2wGiven)
            model->BSIM4pdvt2w = 0.0;
        if (!model->BSIM4pdroutGiven)
            model->BSIM4pdrout = 0.0;     
        if (!model->BSIM4pdsubGiven)
            model->BSIM4pdsub = 0.0;
        if (!model->BSIM4pvth0Given)
           model->BSIM4pvth0 = 0.0;
        if (!model->BSIM4puaGiven)
            model->BSIM4pua = 0.0;
        if (!model->BSIM4pua1Given)
            model->BSIM4pua1 = 0.0;
        if (!model->BSIM4pubGiven)
            model->BSIM4pub = 0.0;
        if (!model->BSIM4pub1Given)
            model->BSIM4pub1 = 0.0;
        if (!model->BSIM4pucGiven)
            model->BSIM4puc = 0.0;
        if (!model->BSIM4puc1Given)
            model->BSIM4puc1 = 0.0;
        if (!model->BSIM4pudGiven)
            model->BSIM4pud = 0.0;
        if (!model->BSIM4pud1Given)
            model->BSIM4pud1 = 0.0;
        if (!model->BSIM4pupGiven)
            model->BSIM4pup = 0.0;
        if (!model->BSIM4plpGiven)
            model->BSIM4plp = 0.0;
        if (!model->BSIM4pu0Given)
            model->BSIM4pu0 = 0.0;
        if (!model->BSIM4puteGiven)
            model->BSIM4pute = 0.0;  
     if (!model->BSIM4pucsteGiven)
            model->BSIM4pucste = 0.0;                 
        if (!model->BSIM4pvoffGiven)
            model->BSIM4pvoff = 0.0;
        if (!model->BSIM4pminvGiven)
            model->BSIM4pminv = 0.0;
        if (!model->BSIM4pminvcvGiven)
            model->BSIM4pminvcv = 0.0;
        if (!model->BSIM4pfproutGiven)
            model->BSIM4pfprout = 0.0;
        if (!model->BSIM4ppditsGiven)
            model->BSIM4ppdits = 0.0;
        if (!model->BSIM4ppditsdGiven)
            model->BSIM4ppditsd = 0.0;
        if (!model->BSIM4pdeltaGiven)  
            model->BSIM4pdelta = 0.0;
        if (!model->BSIM4prdswGiven)
            model->BSIM4prdsw = 0.0;
        if (!model->BSIM4prdwGiven)
            model->BSIM4prdw = 0.0;
        if (!model->BSIM4prswGiven)
            model->BSIM4prsw = 0.0;
        if (!model->BSIM4pprwbGiven)
            model->BSIM4pprwb = 0.0;
        if (!model->BSIM4pprwgGiven)
            model->BSIM4pprwg = 0.0;
        if (!model->BSIM4pprtGiven)
            model->BSIM4pprt = 0.0;
        if (!model->BSIM4peta0Given)
            model->BSIM4peta0 = 0.0;
        if (!model->BSIM4petabGiven)
            model->BSIM4petab = 0.0;
        if (!model->BSIM4ppclmGiven)
            model->BSIM4ppclm = 0.0; 
        if (!model->BSIM4ppdibl1Given)
            model->BSIM4ppdibl1 = 0.0;
        if (!model->BSIM4ppdibl2Given)
            model->BSIM4ppdibl2 = 0.0;
        if (!model->BSIM4ppdiblbGiven)
            model->BSIM4ppdiblb = 0.0;
        if (!model->BSIM4ppscbe1Given)
            model->BSIM4ppscbe1 = 0.0;
        if (!model->BSIM4ppscbe2Given)
            model->BSIM4ppscbe2 = 0.0;
        if (!model->BSIM4ppvagGiven)
            model->BSIM4ppvag = 0.0;     
        if (!model->BSIM4pwrGiven)  
            model->BSIM4pwr = 0.0;
        if (!model->BSIM4pdwgGiven)  
            model->BSIM4pdwg = 0.0;
        if (!model->BSIM4pdwbGiven)  
            model->BSIM4pdwb = 0.0;
        if (!model->BSIM4pb0Given)
            model->BSIM4pb0 = 0.0;
        if (!model->BSIM4pb1Given)  
            model->BSIM4pb1 = 0.0;
        if (!model->BSIM4palpha0Given)  
            model->BSIM4palpha0 = 0.0;
        if (!model->BSIM4palpha1Given)
            model->BSIM4palpha1 = 0.0;
        if (!model->BSIM4pbeta0Given)  
            model->BSIM4pbeta0 = 0.0;
        if (!model->BSIM4pagidlGiven)
            model->BSIM4pagidl = 0.0;
        if (!model->BSIM4pbgidlGiven)
            model->BSIM4pbgidl = 0.0;
        if (!model->BSIM4pcgidlGiven)
            model->BSIM4pcgidl = 0.0;
        if (!model->BSIM4pegidlGiven)
            model->BSIM4pegidl = 0.0;
        if (!model->BSIM4prgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4prgidl = 0.0;
        if (!model->BSIM4pkgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4pkgidl = 0.0;
        if (!model->BSIM4pfgidlGiven)        /* v4.7 New GIDL/GISL */
            model->BSIM4pfgidl = 0.0;

        /*if (!model->BSIM4pagislGiven)
        {
            if (model->BSIM4pagidlGiven)
                model->BSIM4pagisl = model->BSIM4pagidl;
            else
                model->BSIM4pagisl = 0.0;
        }
        if (!model->BSIM4pbgislGiven)
        {
            if (model->BSIM4pbgidlGiven)
                model->BSIM4pbgisl = model->BSIM4pbgidl;
            else
                model->BSIM4pbgisl = 0.0;
        }
        if (!model->BSIM4pcgislGiven)
        {
            if (model->BSIM4pcgidlGiven)
                model->BSIM4pcgisl = model->BSIM4pcgidl;
            else
                model->BSIM4pcgisl = 0.0;
        }
        if (!model->BSIM4pegislGiven)
        {
            if (model->BSIM4pegidlGiven)
                model->BSIM4pegisl = model->BSIM4pegidl;
            else
                model->BSIM4pegisl = 0.0; 
        }*/

        /*if (!model->BSIM4prgislGiven)
        {
            if (model->BSIM4prgidlGiven)
                model->BSIM4prgisl = model->BSIM4prgidl;
        }
        if (!model->BSIM4pkgislGiven)
        {
            if (model->BSIM4pkgidlGiven)
                model->BSIM4pkgisl = model->BSIM4pkgidl;
        }
        if (!model->BSIM4pfgislGiven)
        {
            if (model->BSIM4pfgidlGiven)
                model->BSIM4pfgisl = model->BSIM4pfgidl;
        }*/

        /*Default value of pagisl, pbgisl, pcgisl, pegisl, prgisl, pkgisl, and pfgisl are set as follows */
        if (!model->BSIM4pagislGiven)
            model->BSIM4pagisl = model->BSIM4pagidl;
        if (!model->BSIM4pbgislGiven)
            model->BSIM4pbgisl = model->BSIM4pbgidl;
        if (!model->BSIM4pcgislGiven)
            model->BSIM4pcgisl = model->BSIM4pcgidl;
        if (!model->BSIM4pegislGiven)
            model->BSIM4pegisl = model->BSIM4pegidl;
        if (!model->BSIM4prgislGiven)   /* v4.7 New GIDL/GISL */
            model->BSIM4prgisl = model->BSIM4prgidl;
        if (!model->BSIM4pkgislGiven)   /* v4.7 New GIDL/GISL */
            model->BSIM4pkgisl = model->BSIM4pkgidl;
        if (!model->BSIM4pfgislGiven)   /* v4.7 New GIDL/GISL */
            model->BSIM4pfgisl = model->BSIM4pfgidl;
        if (!model->BSIM4paigcGiven)
            model->BSIM4paigc = 0.0;
        if (!model->BSIM4pbigcGiven)
            model->BSIM4pbigc = 0.0;
        if (!model->BSIM4pcigcGiven)
            model->BSIM4pcigc = 0.0;
        if (!model->BSIM4aigsdGiven && (model->BSIM4aigsGiven || model->BSIM4aigdGiven))
        {
            if (!model->BSIM4paigsGiven)
                model->BSIM4paigs = 0.0;
            if (!model->BSIM4paigdGiven)
                model->BSIM4paigd = 0.0;
        }
        else
        {
           if (!model->BSIM4paigsdGiven)
               model->BSIM4paigsd = 0.0;
           model->BSIM4paigs = model->BSIM4paigd = model->BSIM4paigsd;
        }
        if (!model->BSIM4bigsdGiven && (model->BSIM4bigsGiven || model->BSIM4bigdGiven))
        {
            if (!model->BSIM4pbigsGiven)
                model->BSIM4pbigs = 0.0;
            if (!model->BSIM4pbigdGiven)
                model->BSIM4pbigd = 0.0;
        }
        else
        {
           if (!model->BSIM4pbigsdGiven)
               model->BSIM4pbigsd = 0.0;
           model->BSIM4pbigs = model->BSIM4pbigd = model->BSIM4pbigsd;
        }
        if (!model->BSIM4cigsdGiven && (model->BSIM4cigsGiven || model->BSIM4cigdGiven))
        {
            if (!model->BSIM4pcigsGiven)
                model->BSIM4pcigs = 0.0;
            if (!model->BSIM4pcigdGiven)
                model->BSIM4pcigd = 0.0;
        }
        else
        {
           if (!model->BSIM4pcigsdGiven)
               model->BSIM4pcigsd = 0.0;
           model->BSIM4pcigs = model->BSIM4pcigd = model->BSIM4pcigsd;
        }
        if (!model->BSIM4paigbaccGiven)
            model->BSIM4paigbacc = 0.0;
        if (!model->BSIM4pbigbaccGiven)
            model->BSIM4pbigbacc = 0.0;
        if (!model->BSIM4pcigbaccGiven)
            model->BSIM4pcigbacc = 0.0;
        if (!model->BSIM4paigbinvGiven)
            model->BSIM4paigbinv = 0.0;
        if (!model->BSIM4pbigbinvGiven)
            model->BSIM4pbigbinv = 0.0;
        if (!model->BSIM4pcigbinvGiven)
            model->BSIM4pcigbinv = 0.0;
        if (!model->BSIM4pnigcGiven)
            model->BSIM4pnigc = 0.0;
        if (!model->BSIM4pnigbinvGiven)
            model->BSIM4pnigbinv = 0.0;
        if (!model->BSIM4pnigbaccGiven)
            model->BSIM4pnigbacc = 0.0;
        if (!model->BSIM4pntoxGiven)
            model->BSIM4pntox = 0.0;
        if (!model->BSIM4peigbinvGiven)
            model->BSIM4peigbinv = 0.0;
        if (!model->BSIM4ppigcdGiven)
            model->BSIM4ppigcd = 0.0;
        if (!model->BSIM4ppoxedgeGiven)
            model->BSIM4ppoxedge = 0.0;
        if (!model->BSIM4pxrcrg1Given)
            model->BSIM4pxrcrg1 = 0.0;
        if (!model->BSIM4pxrcrg2Given)
            model->BSIM4pxrcrg2 = 0.0;
        if (!model->BSIM4peuGiven)
            model->BSIM4peu = 0.0;
                if (!model->BSIM4pucsGiven)
            model->BSIM4pucs = 0.0;
        if (!model->BSIM4pvfbGiven)
            model->BSIM4pvfb = 0.0;
        if (!model->BSIM4plambdaGiven)
            model->BSIM4plambda = 0.0;
        if (!model->BSIM4pvtlGiven)
            model->BSIM4pvtl = 0.0;  
        if (!model->BSIM4pxnGiven)
            model->BSIM4pxn = 0.0;  
        if (!model->BSIM4pvfbsdoffGiven)
            model->BSIM4pvfbsdoff = 0.0;   
        if (!model->BSIM4ptvfbsdoffGiven)
            model->BSIM4ptvfbsdoff = 0.0;  
        if (!model->BSIM4ptvoffGiven)
            model->BSIM4ptvoff = 0.0;  
        if (!model->BSIM4ptnfactorGiven)         /* v4.7 temp dep of leakage current  */
            model->BSIM4ptnfactor = 0.0;  
        if (!model->BSIM4pteta0Given)                /* v4.7 temp dep of leakage current  */
            model->BSIM4pteta0 = 0.0;  
        if (!model->BSIM4ptvoffcvGiven)                /* v4.7 temp dep of leakage current  */
            model->BSIM4ptvoffcv = 0.0;  

        if (!model->BSIM4pcgslGiven)  
            model->BSIM4pcgsl = 0.0;
        if (!model->BSIM4pcgdlGiven)  
            model->BSIM4pcgdl = 0.0;
        if (!model->BSIM4pckappasGiven)  
            model->BSIM4pckappas = 0.0;
        if (!model->BSIM4pckappadGiven)
            model->BSIM4pckappad = 0.0;
        if (!model->BSIM4pcfGiven)  
            model->BSIM4pcf = 0.0;
        if (!model->BSIM4pclcGiven)  
            model->BSIM4pclc = 0.0;
        if (!model->BSIM4pcleGiven)  
            model->BSIM4pcle = 0.0;
        if (!model->BSIM4pvfbcvGiven)  
            model->BSIM4pvfbcv = 0.0;
        if (!model->BSIM4pacdeGiven)
            model->BSIM4pacde = 0.0;
        if (!model->BSIM4pmoinGiven)
            model->BSIM4pmoin = 0.0;
        if (!model->BSIM4pnoffGiven)
            model->BSIM4pnoff = 0.0;
        if (!model->BSIM4pvoffcvGiven)
            model->BSIM4pvoffcv = 0.0;

        if (!model->BSIM4gamma1Given)
            model->BSIM4gamma1 = 0.0;
        if (!model->BSIM4lgamma1Given)
            model->BSIM4lgamma1 = 0.0;
        if (!model->BSIM4wgamma1Given)
            model->BSIM4wgamma1 = 0.0;
        if (!model->BSIM4pgamma1Given)
            model->BSIM4pgamma1 = 0.0;
        if (!model->BSIM4gamma2Given)
            model->BSIM4gamma2 = 0.0;
        if (!model->BSIM4lgamma2Given)
            model->BSIM4lgamma2 = 0.0;
        if (!model->BSIM4wgamma2Given)
            model->BSIM4wgamma2 = 0.0;
        if (!model->BSIM4pgamma2Given)
            model->BSIM4pgamma2 = 0.0;
        if (!model->BSIM4vbxGiven)
            model->BSIM4vbx = 0.0;
        if (!model->BSIM4lvbxGiven)
            model->BSIM4lvbx = 0.0;
        if (!model->BSIM4wvbxGiven)
            model->BSIM4wvbx = 0.0;
        if (!model->BSIM4pvbxGiven)
            model->BSIM4pvbx = 0.0;

        /* unit degree celcius */
        if (!model->BSIM4tnomGiven)  
            model->BSIM4tnom = ckt->CKTnomTemp; 
        if (!model->BSIM4LintGiven)  
           model->BSIM4Lint = 0.0;
        if (!model->BSIM4LlGiven)  
           model->BSIM4Ll = 0.0;
        if (!model->BSIM4LlcGiven)
           model->BSIM4Llc = model->BSIM4Ll;
        if (!model->BSIM4LlnGiven)  
           model->BSIM4Lln = 1.0;
        if (!model->BSIM4LwGiven)  
           model->BSIM4Lw = 0.0;
        if (!model->BSIM4LwcGiven)
           model->BSIM4Lwc = model->BSIM4Lw;
        if (!model->BSIM4LwnGiven)  
           model->BSIM4Lwn = 1.0;
        if (!model->BSIM4LwlGiven)  
           model->BSIM4Lwl = 0.0;
        if (!model->BSIM4LwlcGiven)
           model->BSIM4Lwlc = model->BSIM4Lwl;
        if (!model->BSIM4LminGiven)  
           model->BSIM4Lmin = 0.0;
        if (!model->BSIM4LmaxGiven)  
           model->BSIM4Lmax = 1.0;
        if (!model->BSIM4WintGiven)  
           model->BSIM4Wint = 0.0;
        if (!model->BSIM4WlGiven)  
           model->BSIM4Wl = 0.0;
        if (!model->BSIM4WlcGiven)
           model->BSIM4Wlc = model->BSIM4Wl;
        if (!model->BSIM4WlnGiven)  
           model->BSIM4Wln = 1.0;
        if (!model->BSIM4WwGiven)  
           model->BSIM4Ww = 0.0;
        if (!model->BSIM4WwcGiven)
           model->BSIM4Wwc = model->BSIM4Ww;
        if (!model->BSIM4WwnGiven)  
           model->BSIM4Wwn = 1.0;
        if (!model->BSIM4WwlGiven)  
           model->BSIM4Wwl = 0.0;
        if (!model->BSIM4WwlcGiven)
           model->BSIM4Wwlc = model->BSIM4Wwl;
        if (!model->BSIM4WminGiven)  
           model->BSIM4Wmin = 0.0;
        if (!model->BSIM4WmaxGiven)  
           model->BSIM4Wmax = 1.0;
        if (!model->BSIM4dwcGiven)  
           model->BSIM4dwc = model->BSIM4Wint;
        if (!model->BSIM4dlcGiven)  
           model->BSIM4dlc = model->BSIM4Lint;
        if (!model->BSIM4xlGiven)  
           model->BSIM4xl = 0.0;
        if (!model->BSIM4xwGiven)  
           model->BSIM4xw = 0.0;
        if (!model->BSIM4dlcigGiven)
           model->BSIM4dlcig = model->BSIM4Lint;
        if (!model->BSIM4dlcigdGiven)
        {
           if (model->BSIM4dlcigGiven) 
               model->BSIM4dlcigd = model->BSIM4dlcig;
           else             
               model->BSIM4dlcigd = model->BSIM4Lint;
        }
        if (!model->BSIM4dwjGiven)
           model->BSIM4dwj = model->BSIM4dwc;
        if (!model->BSIM4cfGiven)
           model->BSIM4cf = 2.0 * model->BSIM4epsrox * EPS0 / PI
                          * log(1.0 + 0.4e-6 / model->BSIM4toxe);

        if (!model->BSIM4xpartGiven)
            model->BSIM4xpart = 0.0;
        if (!model->BSIM4sheetResistanceGiven)
            model->BSIM4sheetResistance = 0.0;

        if (!model->BSIM4SunitAreaJctCapGiven)
            model->BSIM4SunitAreaJctCap = 5.0E-4;
        if (!model->BSIM4DunitAreaJctCapGiven)
            model->BSIM4DunitAreaJctCap = model->BSIM4SunitAreaJctCap;
        if (!model->BSIM4SunitLengthSidewallJctCapGiven)
            model->BSIM4SunitLengthSidewallJctCap = 5.0E-10;
        if (!model->BSIM4DunitLengthSidewallJctCapGiven)
            model->BSIM4DunitLengthSidewallJctCap = model->BSIM4SunitLengthSidewallJctCap;
        if (!model->BSIM4SunitLengthGateSidewallJctCapGiven)
            model->BSIM4SunitLengthGateSidewallJctCap = model->BSIM4SunitLengthSidewallJctCap ;
        if (!model->BSIM4DunitLengthGateSidewallJctCapGiven)
            model->BSIM4DunitLengthGateSidewallJctCap = model->BSIM4SunitLengthGateSidewallJctCap;
        if (!model->BSIM4SjctSatCurDensityGiven)
            model->BSIM4SjctSatCurDensity = 1.0E-4;
        if (!model->BSIM4DjctSatCurDensityGiven)
            model->BSIM4DjctSatCurDensity = model->BSIM4SjctSatCurDensity;
        if (!model->BSIM4SjctSidewallSatCurDensityGiven)
            model->BSIM4SjctSidewallSatCurDensity = 0.0;
        if (!model->BSIM4DjctSidewallSatCurDensityGiven)
            model->BSIM4DjctSidewallSatCurDensity = model->BSIM4SjctSidewallSatCurDensity;
        if (!model->BSIM4SjctGateSidewallSatCurDensityGiven)
            model->BSIM4SjctGateSidewallSatCurDensity = 0.0;
        if (!model->BSIM4DjctGateSidewallSatCurDensityGiven)
            model->BSIM4DjctGateSidewallSatCurDensity = model->BSIM4SjctGateSidewallSatCurDensity;
        if (!model->BSIM4SbulkJctPotentialGiven)
            model->BSIM4SbulkJctPotential = 1.0;
        if (!model->BSIM4DbulkJctPotentialGiven)
            model->BSIM4DbulkJctPotential = model->BSIM4SbulkJctPotential;
        if (!model->BSIM4SsidewallJctPotentialGiven)
            model->BSIM4SsidewallJctPotential = 1.0;
        if (!model->BSIM4DsidewallJctPotentialGiven)
            model->BSIM4DsidewallJctPotential = model->BSIM4SsidewallJctPotential;
        if (!model->BSIM4SGatesidewallJctPotentialGiven)
            model->BSIM4SGatesidewallJctPotential = model->BSIM4SsidewallJctPotential;
        if (!model->BSIM4DGatesidewallJctPotentialGiven)
            model->BSIM4DGatesidewallJctPotential = model->BSIM4SGatesidewallJctPotential;
        if (!model->BSIM4SbulkJctBotGradingCoeffGiven)
            model->BSIM4SbulkJctBotGradingCoeff = 0.5;
        if (!model->BSIM4DbulkJctBotGradingCoeffGiven)
            model->BSIM4DbulkJctBotGradingCoeff = model->BSIM4SbulkJctBotGradingCoeff;
        if (!model->BSIM4SbulkJctSideGradingCoeffGiven)
            model->BSIM4SbulkJctSideGradingCoeff = 0.33;
        if (!model->BSIM4DbulkJctSideGradingCoeffGiven)
            model->BSIM4DbulkJctSideGradingCoeff = model->BSIM4SbulkJctSideGradingCoeff;
        if (!model->BSIM4SbulkJctGateSideGradingCoeffGiven)
            model->BSIM4SbulkJctGateSideGradingCoeff = model->BSIM4SbulkJctSideGradingCoeff;
        if (!model->BSIM4DbulkJctGateSideGradingCoeffGiven)
            model->BSIM4DbulkJctGateSideGradingCoeff = model->BSIM4SbulkJctGateSideGradingCoeff;
        if (!model->BSIM4SjctEmissionCoeffGiven)
            model->BSIM4SjctEmissionCoeff = 1.0;
        if (!model->BSIM4DjctEmissionCoeffGiven)
            model->BSIM4DjctEmissionCoeff = model->BSIM4SjctEmissionCoeff;
        if (!model->BSIM4SjctTempExponentGiven)
            model->BSIM4SjctTempExponent = 3.0;
        if (!model->BSIM4DjctTempExponentGiven)
            model->BSIM4DjctTempExponent = model->BSIM4SjctTempExponent;

        if (!model->BSIM4jtssGiven)
            model->BSIM4jtss = 0.0;
        if (!model->BSIM4jtsdGiven)
            model->BSIM4jtsd = model->BSIM4jtss;
        if (!model->BSIM4jtsswsGiven)
            model->BSIM4jtssws = 0.0;
        if (!model->BSIM4jtsswdGiven)
            model->BSIM4jtsswd = model->BSIM4jtssws;
        if (!model->BSIM4jtsswgsGiven)
            model->BSIM4jtsswgs = 0.0;
        if (!model->BSIM4jtsswgdGiven)
            model->BSIM4jtsswgd = model->BSIM4jtsswgs;
                if (!model->BSIM4jtweffGiven)
                    model->BSIM4jtweff = 0.0;
        if (!model->BSIM4njtsGiven)
            model->BSIM4njts = 20.0;
        if (!model->BSIM4njtsswGiven)
            model->BSIM4njtssw = 20.0;
        if (!model->BSIM4njtsswgGiven)
            model->BSIM4njtsswg = 20.0;
        if (!model->BSIM4njtsdGiven)
        {
            if (model->BSIM4njtsGiven)
                model->BSIM4njtsd =  model->BSIM4njts;
            else
              model->BSIM4njtsd = 20.0;
        }
        if (!model->BSIM4njtsswdGiven)
        {
            if (model->BSIM4njtsswGiven)
                model->BSIM4njtsswd =  model->BSIM4njtssw;
            else
              model->BSIM4njtsswd = 20.0;
        }
        if (!model->BSIM4njtsswgdGiven)
        {
            if (model->BSIM4njtsswgGiven)
                model->BSIM4njtsswgd =  model->BSIM4njtsswg;
            else
              model->BSIM4njtsswgd = 20.0;
        }
        if (!model->BSIM4xtssGiven)
            model->BSIM4xtss = 0.02;
        if (!model->BSIM4xtsdGiven)
            model->BSIM4xtsd = model->BSIM4xtss;
        if (!model->BSIM4xtsswsGiven)
            model->BSIM4xtssws = 0.02;
        if (!model->BSIM4xtsswdGiven)
            model->BSIM4xtsswd = model->BSIM4xtssws;
        if (!model->BSIM4xtsswgsGiven)
            model->BSIM4xtsswgs = 0.02;
        if (!model->BSIM4xtsswgdGiven)
            model->BSIM4xtsswgd = model->BSIM4xtsswgs;
        if (!model->BSIM4tnjtsGiven)
            model->BSIM4tnjts = 0.0;
        if (!model->BSIM4tnjtsswGiven)
            model->BSIM4tnjtssw = 0.0;
        if (!model->BSIM4tnjtsswgGiven)
            model->BSIM4tnjtsswg = 0.0;
        if (!model->BSIM4tnjtsdGiven)
        {
            if (model->BSIM4tnjtsGiven)
                model->BSIM4tnjtsd =  model->BSIM4tnjts;
            else
              model->BSIM4tnjtsd = 0.0;
        }
        if (!model->BSIM4tnjtsswdGiven)
        {
            if (model->BSIM4tnjtsswGiven)
                model->BSIM4tnjtsswd =  model->BSIM4tnjtssw;
            else
              model->BSIM4tnjtsswd = 0.0;
        }
        if (!model->BSIM4tnjtsswgdGiven)
        {
            if (model->BSIM4tnjtsswgGiven)
                model->BSIM4tnjtsswgd =  model->BSIM4tnjtsswg;
            else
              model->BSIM4tnjtsswgd = 0.0;
        }
        if (!model->BSIM4vtssGiven)
            model->BSIM4vtss = 10.0;
        if (!model->BSIM4vtsdGiven)
            model->BSIM4vtsd = model->BSIM4vtss;
        if (!model->BSIM4vtsswsGiven)
            model->BSIM4vtssws = 10.0;
        if (!model->BSIM4vtsswdGiven)
            model->BSIM4vtsswd = model->BSIM4vtssws;
        if (!model->BSIM4vtsswgsGiven)
            model->BSIM4vtsswgs = 10.0;
        if (!model->BSIM4vtsswgdGiven)
            model->BSIM4vtsswgd = model->BSIM4vtsswgs;

        if (!model->BSIM4oxideTrapDensityAGiven)
        {   if (model->BSIM4type == NMOS)
                model->BSIM4oxideTrapDensityA = 6.25e41;
            else
                model->BSIM4oxideTrapDensityA= 6.188e40;
        }
        if (!model->BSIM4oxideTrapDensityBGiven)
        {   if (model->BSIM4type == NMOS)
                model->BSIM4oxideTrapDensityB = 3.125e26;
            else
                model->BSIM4oxideTrapDensityB = 1.5e25;
        }
        if (!model->BSIM4oxideTrapDensityCGiven)
            model->BSIM4oxideTrapDensityC = 8.75e9;
        if (!model->BSIM4emGiven)
            model->BSIM4em = 4.1e7; /* V/m */
        if (!model->BSIM4efGiven)
            model->BSIM4ef = 1.0;
        if (!model->BSIM4afGiven)
            model->BSIM4af = 1.0;
        if (!model->BSIM4kfGiven)
            model->BSIM4kf = 0.0;

        if (!model->BSIM4vgsMaxGiven)
            model->BSIM4vgsMax = 1e99;
        if (!model->BSIM4vgdMaxGiven)
            model->BSIM4vgdMax = 1e99;
        if (!model->BSIM4vgbMaxGiven)
            model->BSIM4vgbMax = 1e99;
        if (!model->BSIM4vdsMaxGiven)
            model->BSIM4vdsMax = 1e99;
        if (!model->BSIM4vbsMaxGiven)
            model->BSIM4vbsMax = 1e99;
        if (!model->BSIM4vbdMaxGiven)
            model->BSIM4vbdMax = 1e99;

        /* stress effect */
        if (!model->BSIM4sarefGiven)
            model->BSIM4saref = 1e-6; /* m */
        if (!model->BSIM4sbrefGiven)
            model->BSIM4sbref = 1e-6;  /* m */
        if (!model->BSIM4wlodGiven)
            model->BSIM4wlod = 0;  /* m */
        if (!model->BSIM4ku0Given)
            model->BSIM4ku0 = 0; /* 1/m */
        if (!model->BSIM4kvsatGiven)
            model->BSIM4kvsat = 0;
        if (!model->BSIM4kvth0Given) /* m */
            model->BSIM4kvth0 = 0;
        if (!model->BSIM4tku0Given)
            model->BSIM4tku0 = 0;
        if (!model->BSIM4llodku0Given)
            model->BSIM4llodku0 = 0;
        if (!model->BSIM4wlodku0Given)
            model->BSIM4wlodku0 = 0;
        if (!model->BSIM4llodvthGiven)
            model->BSIM4llodvth = 0;
        if (!model->BSIM4wlodvthGiven)
            model->BSIM4wlodvth = 0;
        if (!model->BSIM4lku0Given)
            model->BSIM4lku0 = 0;
        if (!model->BSIM4wku0Given)
            model->BSIM4wku0 = 0;
        if (!model->BSIM4pku0Given)
            model->BSIM4pku0 = 0;
        if (!model->BSIM4lkvth0Given)
            model->BSIM4lkvth0 = 0;
        if (!model->BSIM4wkvth0Given)
            model->BSIM4wkvth0 = 0;
        if (!model->BSIM4pkvth0Given)
            model->BSIM4pkvth0 = 0;
        if (!model->BSIM4stk2Given)
            model->BSIM4stk2 = 0;
        if (!model->BSIM4lodk2Given)
            model->BSIM4lodk2 = 1.0;
        if (!model->BSIM4steta0Given)
            model->BSIM4steta0 = 0;
        if (!model->BSIM4lodeta0Given)
            model->BSIM4lodeta0 = 1.0;

        /* Well Proximity Effect  */
        if (!model->BSIM4webGiven)
            model->BSIM4web = 0.0; 
        if (!model->BSIM4wecGiven)
            model->BSIM4wec = 0.0;
        if (!model->BSIM4kvth0weGiven)
            model->BSIM4kvth0we = 0.0; 
        if (!model->BSIM4k2weGiven)
            model->BSIM4k2we = 0.0; 
        if (!model->BSIM4ku0weGiven)
            model->BSIM4ku0we = 0.0; 
        if (!model->BSIM4screfGiven)
            model->BSIM4scref = 1.0E-6; /* m */
        if (!model->BSIM4wpemodGiven)
            model->BSIM4wpemod = 0; 
        else if ((model->BSIM4wpemod != 0) && (model->BSIM4wpemod != 1))
        {   model->BSIM4wpemod = 0;
            printf("Warning: wpemod has been set to its default value: 0.\n");
        }
        if (!model->BSIM4lkvth0weGiven)
            model->BSIM4lkvth0we = 0; 
        if (!model->BSIM4lk2weGiven)
            model->BSIM4lk2we = 0;
        if (!model->BSIM4lku0weGiven)
            model->BSIM4lku0we = 0;
        if (!model->BSIM4wkvth0weGiven)
            model->BSIM4wkvth0we = 0; 
        if (!model->BSIM4wk2weGiven)
            model->BSIM4wk2we = 0;
        if (!model->BSIM4wku0weGiven)
            model->BSIM4wku0we = 0;
        if (!model->BSIM4pkvth0weGiven)
            model->BSIM4pkvth0we = 0; 
        if (!model->BSIM4pk2weGiven)
            model->BSIM4pk2we = 0;
        if (!model->BSIM4pku0weGiven)
            model->BSIM4pku0we = 0;

        DMCGeff = model->BSIM4dmcg - model->BSIM4dmcgt;
        DMCIeff = model->BSIM4dmci;
        DMDGeff = model->BSIM4dmdg - model->BSIM4dmcgt;

        /*
         * End processing models and begin to loop
         * through all the instances of the model
         */

        for (here = model->BSIM4instances; here != NULL ;
             here=here->BSIM4nextInstance) 
        {
            /* allocate a chunk of the state vector */
            here->BSIM4states = *states;
            *states += BSIM4numStates;

            /* perform the parameter defaulting */
            if (!here->BSIM4lGiven)
                here->BSIM4l = 5.0e-6;
            if (!here->BSIM4wGiven)
                here->BSIM4w = 5.0e-6;
            if (!here->BSIM4mGiven)
                here->BSIM4m = 1.0;
            if (!here->BSIM4nfGiven)
                here->BSIM4nf = 1.0;
            if (!here->BSIM4minGiven)
                here->BSIM4min = 0; /* integer */
            if (!here->BSIM4icVDSGiven)
                here->BSIM4icVDS = 0.0;
            if (!here->BSIM4icVGSGiven)
                here->BSIM4icVGS = 0.0;
            if (!here->BSIM4icVBSGiven)
                here->BSIM4icVBS = 0.0;
            if (!here->BSIM4drainAreaGiven)
                here->BSIM4drainArea = 0.0;
            if (!here->BSIM4drainPerimeterGiven)
                here->BSIM4drainPerimeter = 0.0;
            if (!here->BSIM4drainSquaresGiven)
                here->BSIM4drainSquares = 1.0;
            if (!here->BSIM4sourceAreaGiven)
                here->BSIM4sourceArea = 0.0;
            if (!here->BSIM4sourcePerimeterGiven)
                here->BSIM4sourcePerimeter = 0.0;
            if (!here->BSIM4sourceSquaresGiven)
                here->BSIM4sourceSquares = 1.0;

            if (!here->BSIM4rbdbGiven)
                here->BSIM4rbdb = model->BSIM4rbdb; /* in ohm */
            if (!here->BSIM4rbsbGiven)
                here->BSIM4rbsb = model->BSIM4rbsb;
            if (!here->BSIM4rbpbGiven)
                here->BSIM4rbpb = model->BSIM4rbpb;
            if (!here->BSIM4rbpsGiven)
                here->BSIM4rbps = model->BSIM4rbps;
            if (!here->BSIM4rbpdGiven)
                here->BSIM4rbpd = model->BSIM4rbpd;
            if (!here->BSIM4delvtoGiven)
                here->BSIM4delvto = 0.0;
            if (!here->BSIM4xgwGiven)
                here->BSIM4xgw = model->BSIM4xgw;
            if (!here->BSIM4ngconGiven)
                here->BSIM4ngcon = model->BSIM4ngcon;

                    
            /* Process instance model selectors, some
             * may override their global counterparts
             */
            if (!here->BSIM4rbodyModGiven)
                here->BSIM4rbodyMod = model->BSIM4rbodyMod;
            else if ((here->BSIM4rbodyMod != 0) && (here->BSIM4rbodyMod != 1) && (here->BSIM4rbodyMod != 2))
            {   here->BSIM4rbodyMod = model->BSIM4rbodyMod;
                printf("Warning: rbodyMod has been set to its global value %d.\n",
                model->BSIM4rbodyMod);
            }

            if (!here->BSIM4rgateModGiven)
                here->BSIM4rgateMod = model->BSIM4rgateMod;
            else if ((here->BSIM4rgateMod != 0) && (here->BSIM4rgateMod != 1)
                && (here->BSIM4rgateMod != 2) && (here->BSIM4rgateMod != 3))
            {   here->BSIM4rgateMod = model->BSIM4rgateMod;
                printf("Warning: rgateMod has been set to its global value %d.\n",
                model->BSIM4rgateMod);
            }

            if (!here->BSIM4geoModGiven)
                here->BSIM4geoMod = model->BSIM4geoMod;

            if (!here->BSIM4rgeoModGiven)
                here->BSIM4rgeoMod = model->BSIM4rgeoMod;
            else if ((here->BSIM4rgeoMod != 0) && (here->BSIM4rgeoMod != 1))
            {   here->BSIM4rgeoMod = model->BSIM4rgeoMod;
                printf("Warning: rgeoMod has been set to its global value %d.\n",
                model->BSIM4rgeoMod);
            }

            if (!here->BSIM4trnqsModGiven)
                here->BSIM4trnqsMod = model->BSIM4trnqsMod;
            else if ((here->BSIM4trnqsMod != 0) && (here->BSIM4trnqsMod != 1))
            {   here->BSIM4trnqsMod = model->BSIM4trnqsMod;
                printf("Warning: trnqsMod has been set to its global value %d.\n",
                model->BSIM4trnqsMod);
            }

            if (!here->BSIM4acnqsModGiven)
                here->BSIM4acnqsMod = model->BSIM4acnqsMod;
            else if ((here->BSIM4acnqsMod != 0) && (here->BSIM4acnqsMod != 1))
            {   here->BSIM4acnqsMod = model->BSIM4acnqsMod;
                printf("Warning: acnqsMod has been set to its global value %d.\n",
                model->BSIM4acnqsMod);
            }

            /* stress effect */
            if (!here->BSIM4saGiven)
                here->BSIM4sa = 0.0;
            if (!here->BSIM4sbGiven)
                here->BSIM4sb = 0.0;
            if (!here->BSIM4sdGiven)
                here->BSIM4sd = 2 * model->BSIM4dmcg;
            /* Well Proximity Effect  */
            if (!here->BSIM4scaGiven)
                here->BSIM4sca = 0.0;
            if (!here->BSIM4scbGiven)
                here->BSIM4scb = 0.0;
            if (!here->BSIM4sccGiven)
                here->BSIM4scc = 0.0;
            if (!here->BSIM4scGiven)
                here->BSIM4sc = 0.0; /* m */

            /* process drain series resistance */
            createNode = 0;
            if ( (model->BSIM4rdsMod != 0)
                            || (model->BSIM4tnoiMod == 1 && noiseAnalGiven))
            {
               createNode = 1;
            } else if (model->BSIM4sheetResistance > 0)
            {
                     if (here->BSIM4drainSquaresGiven
                                       && here->BSIM4drainSquares > 0)
                     {
                          createNode = 1;
                     } else if (!here->BSIM4drainSquaresGiven
                                       && (here->BSIM4rgeoMod != 0))
                     {
                          BSIM4RdseffGeo(here->BSIM4nf*here->BSIM4m, here->BSIM4geoMod,
                                  here->BSIM4rgeoMod, here->BSIM4min,
                                  here->BSIM4w, model->BSIM4sheetResistance,
                                  DMCGeff, DMCIeff, DMDGeff, 0, &Rtot);
                          if(Rtot > 0)
                             createNode = 1;
                     }
            }
            if ( createNode != 0 )
            {   if ( here->BSIM4dNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"drain");
                if(error) return(error);
                here->BSIM4dNodePrime = tmp->number;
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
            {   here->BSIM4dNodePrime = here->BSIM4dNode;
            }
            
            /* process source series resistance */
            createNode = 0;
            if ( (model->BSIM4rdsMod != 0)
                            || (model->BSIM4tnoiMod == 1 && noiseAnalGiven))
            {
               createNode = 1;
            } else if (model->BSIM4sheetResistance > 0)
            {
                     if (here->BSIM4sourceSquaresGiven
                                        && here->BSIM4sourceSquares > 0)
                     {
                          createNode = 1;
                     } else if (!here->BSIM4sourceSquaresGiven
                                        && (here->BSIM4rgeoMod != 0))
                     {
                          BSIM4RdseffGeo(here->BSIM4nf*here->BSIM4m, here->BSIM4geoMod,
                                  here->BSIM4rgeoMod, here->BSIM4min,
                                  here->BSIM4w, model->BSIM4sheetResistance,
                                  DMCGeff, DMCIeff, DMDGeff, 1, &Rtot);
                          if(Rtot > 0)
                             createNode = 1;
                     }
            }
            if ( createNode != 0 )
            {   if ( here->BSIM4sNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"source");
                if(error) return(error);
                here->BSIM4sNodePrime = tmp->number;
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
                here->BSIM4sNodePrime = here->BSIM4sNode;

            if ( here->BSIM4rgateMod > 0 )
            {   if ( here->BSIM4gNodePrime == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"gate");
                if(error) return(error);
                   here->BSIM4gNodePrime = tmp->number;
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
                here->BSIM4gNodePrime = here->BSIM4gNodeExt;

            if ( here->BSIM4rgateMod == 3 )
            {   if ( here->BSIM4gNodeMid == 0 )
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"midgate");
                if(error) return(error);
                   here->BSIM4gNodeMid = tmp->number;
            }
            }
            else
                here->BSIM4gNodeMid = here->BSIM4gNodeExt;
            

            /* internal body nodes for body resistance model */
            if ((here->BSIM4rbodyMod ==1) || (here->BSIM4rbodyMod ==2))
            {   if (here->BSIM4dbNode == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"dbody");
                    if(error) return(error);
                    here->BSIM4dbNode = tmp->number;
                }
                if (here->BSIM4bNodePrime == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"body");
                    if(error) return(error);
                    here->BSIM4bNodePrime = tmp->number;
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
                if (here->BSIM4sbNode == 0)
                {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"sbody");
                    if(error) return(error);
                    here->BSIM4sbNode = tmp->number;
                }
            }
            else
                here->BSIM4dbNode = here->BSIM4bNodePrime = here->BSIM4sbNode
                                  = here->BSIM4bNode;

            /* NQS node */
            if ( here->BSIM4trnqsMod )
            {   if ( here->BSIM4qNode == 0 ) 
            {   error = CKTmkVolt(ckt,&tmp,here->BSIM4name,"charge");
                if(error) return(error);
                here->BSIM4qNode = tmp->number;
            }
            }
            else 
                here->BSIM4qNode = 0;
      

/* set Sparse Matrix Pointers 
 * macro to make elements with built-in out-of-memory test */
#define TSTALLOC(ptr,first,second) \
do { if((here->ptr = SMPmakeElt(matrix,here->first,here->second))==(double *)NULL){\
    return(E_NOMEM);\
} } while(0)

            TSTALLOC(BSIM4DPbpPtr, BSIM4dNodePrime, BSIM4bNodePrime);
            TSTALLOC(BSIM4GPbpPtr, BSIM4gNodePrime, BSIM4bNodePrime);
            TSTALLOC(BSIM4SPbpPtr, BSIM4sNodePrime, BSIM4bNodePrime);

            TSTALLOC(BSIM4BPdpPtr, BSIM4bNodePrime, BSIM4dNodePrime);
            TSTALLOC(BSIM4BPgpPtr, BSIM4bNodePrime, BSIM4gNodePrime);
            TSTALLOC(BSIM4BPspPtr, BSIM4bNodePrime, BSIM4sNodePrime);
            TSTALLOC(BSIM4BPbpPtr, BSIM4bNodePrime, BSIM4bNodePrime);

            TSTALLOC(BSIM4DdPtr, BSIM4dNode, BSIM4dNode);
            TSTALLOC(BSIM4GPgpPtr, BSIM4gNodePrime, BSIM4gNodePrime);
            TSTALLOC(BSIM4SsPtr, BSIM4sNode, BSIM4sNode);
            TSTALLOC(BSIM4DPdpPtr, BSIM4dNodePrime, BSIM4dNodePrime);
            TSTALLOC(BSIM4SPspPtr, BSIM4sNodePrime, BSIM4sNodePrime);
            TSTALLOC(BSIM4DdpPtr, BSIM4dNode, BSIM4dNodePrime);
            TSTALLOC(BSIM4GPdpPtr, BSIM4gNodePrime, BSIM4dNodePrime);
            TSTALLOC(BSIM4GPspPtr, BSIM4gNodePrime, BSIM4sNodePrime);
            TSTALLOC(BSIM4SspPtr, BSIM4sNode, BSIM4sNodePrime);
            TSTALLOC(BSIM4DPspPtr, BSIM4dNodePrime, BSIM4sNodePrime);
            TSTALLOC(BSIM4DPdPtr, BSIM4dNodePrime, BSIM4dNode);
            TSTALLOC(BSIM4DPgpPtr, BSIM4dNodePrime, BSIM4gNodePrime);
            TSTALLOC(BSIM4SPgpPtr, BSIM4sNodePrime, BSIM4gNodePrime);
            TSTALLOC(BSIM4SPsPtr, BSIM4sNodePrime, BSIM4sNode);
            TSTALLOC(BSIM4SPdpPtr, BSIM4sNodePrime, BSIM4dNodePrime);

            TSTALLOC(BSIM4QqPtr, BSIM4qNode, BSIM4qNode);
            TSTALLOC(BSIM4QbpPtr, BSIM4qNode, BSIM4bNodePrime) ;
            TSTALLOC(BSIM4QdpPtr, BSIM4qNode, BSIM4dNodePrime);
            TSTALLOC(BSIM4QspPtr, BSIM4qNode, BSIM4sNodePrime);
            TSTALLOC(BSIM4QgpPtr, BSIM4qNode, BSIM4gNodePrime);
            TSTALLOC(BSIM4DPqPtr, BSIM4dNodePrime, BSIM4qNode);
            TSTALLOC(BSIM4SPqPtr, BSIM4sNodePrime, BSIM4qNode);
            TSTALLOC(BSIM4GPqPtr, BSIM4gNodePrime, BSIM4qNode);

            if (here->BSIM4rgateMod != 0)
            {   TSTALLOC(BSIM4GEgePtr, BSIM4gNodeExt, BSIM4gNodeExt);
                TSTALLOC(BSIM4GEgpPtr, BSIM4gNodeExt, BSIM4gNodePrime);
                TSTALLOC(BSIM4GPgePtr, BSIM4gNodePrime, BSIM4gNodeExt);
                TSTALLOC(BSIM4GEdpPtr, BSIM4gNodeExt, BSIM4dNodePrime);
                TSTALLOC(BSIM4GEspPtr, BSIM4gNodeExt, BSIM4sNodePrime);
                TSTALLOC(BSIM4GEbpPtr, BSIM4gNodeExt, BSIM4bNodePrime);
                TSTALLOC(BSIM4GMdpPtr, BSIM4gNodeMid, BSIM4dNodePrime);
                TSTALLOC(BSIM4GMgpPtr, BSIM4gNodeMid, BSIM4gNodePrime);
                TSTALLOC(BSIM4GMgmPtr, BSIM4gNodeMid, BSIM4gNodeMid);
                TSTALLOC(BSIM4GMgePtr, BSIM4gNodeMid, BSIM4gNodeExt);
                TSTALLOC(BSIM4GMspPtr, BSIM4gNodeMid, BSIM4sNodePrime);
                TSTALLOC(BSIM4GMbpPtr, BSIM4gNodeMid, BSIM4bNodePrime);
                TSTALLOC(BSIM4DPgmPtr, BSIM4dNodePrime, BSIM4gNodeMid);
                TSTALLOC(BSIM4GPgmPtr, BSIM4gNodePrime, BSIM4gNodeMid);
                TSTALLOC(BSIM4GEgmPtr, BSIM4gNodeExt, BSIM4gNodeMid);
                TSTALLOC(BSIM4SPgmPtr, BSIM4sNodePrime, BSIM4gNodeMid);
                TSTALLOC(BSIM4BPgmPtr, BSIM4bNodePrime, BSIM4gNodeMid);
            }        

            if ((here->BSIM4rbodyMod ==1) || (here->BSIM4rbodyMod ==2))
            {   TSTALLOC(BSIM4DPdbPtr, BSIM4dNodePrime, BSIM4dbNode);
                TSTALLOC(BSIM4SPsbPtr, BSIM4sNodePrime, BSIM4sbNode);

                TSTALLOC(BSIM4DBdpPtr, BSIM4dbNode, BSIM4dNodePrime);
                TSTALLOC(BSIM4DBdbPtr, BSIM4dbNode, BSIM4dbNode);
                TSTALLOC(BSIM4DBbpPtr, BSIM4dbNode, BSIM4bNodePrime);
                TSTALLOC(BSIM4DBbPtr, BSIM4dbNode, BSIM4bNode);

                TSTALLOC(BSIM4BPdbPtr, BSIM4bNodePrime, BSIM4dbNode);
                TSTALLOC(BSIM4BPbPtr, BSIM4bNodePrime, BSIM4bNode);
                TSTALLOC(BSIM4BPsbPtr, BSIM4bNodePrime, BSIM4sbNode);

                TSTALLOC(BSIM4SBspPtr, BSIM4sbNode, BSIM4sNodePrime);
                TSTALLOC(BSIM4SBbpPtr, BSIM4sbNode, BSIM4bNodePrime);
                TSTALLOC(BSIM4SBbPtr, BSIM4sbNode, BSIM4bNode);
                TSTALLOC(BSIM4SBsbPtr, BSIM4sbNode, BSIM4sbNode);

                TSTALLOC(BSIM4BdbPtr, BSIM4bNode, BSIM4dbNode);
                TSTALLOC(BSIM4BbpPtr, BSIM4bNode, BSIM4bNodePrime);
                TSTALLOC(BSIM4BsbPtr, BSIM4bNode, BSIM4sbNode);
                TSTALLOC(BSIM4BbPtr, BSIM4bNode, BSIM4bNode);
            }

            if (model->BSIM4rdsMod)
            {   TSTALLOC(BSIM4DgpPtr, BSIM4dNode, BSIM4gNodePrime);
                TSTALLOC(BSIM4DspPtr, BSIM4dNode, BSIM4sNodePrime);
                TSTALLOC(BSIM4DbpPtr, BSIM4dNode, BSIM4bNodePrime);
                TSTALLOC(BSIM4SdpPtr, BSIM4sNode, BSIM4dNodePrime);
                TSTALLOC(BSIM4SgpPtr, BSIM4sNode, BSIM4gNodePrime);
                TSTALLOC(BSIM4SbpPtr, BSIM4sNode, BSIM4bNodePrime);
            }
        }
    }

#ifdef USE_OMP
    InstCount = 0;
    model = (BSIM4model*)inModel;
    /* loop through all the BSIM4 device models 
       to count the number of instances */
    
    for( ; model != NULL; model = model->BSIM4nextModel )
    {
        /* loop through all the instances of the model */
        for (here = model->BSIM4instances; here != NULL ;
             here=here->BSIM4nextInstance) 
        { 
            InstCount++;
        }
    }
    InstArray = TMALLOC(BSIM4instance*, InstCount);
    model = (BSIM4model*)inModel;
    idx = 0;
    for( ; model != NULL; model = model->BSIM4nextModel )
    {
        /* loop through all the instances of the model */
        for (here = model->BSIM4instances; here != NULL ;
             here=here->BSIM4nextInstance) 
        { 
            InstArray[idx] = here;
            idx++;
        }
        /* set the array pointer and instance count into each model */
        model->BSIM4InstCount = InstCount;
        model->BSIM4InstanceArray = InstArray;		
    }
#endif

    return(OK);
}  

int
BSIM4unsetup(
GENmodel *inModel,
CKTcircuit *ckt)
{
#ifndef HAS_BATCHSIM
    BSIM4model *model;
    BSIM4instance *here;

    for (model = (BSIM4model *)inModel; model != NULL;
            model = model->BSIM4nextModel)
    {
        for (here = model->BSIM4instances; here != NULL;
                here=here->BSIM4nextInstance)
        {
            if (here->BSIM4dNodePrime
                    && here->BSIM4dNodePrime != here->BSIM4dNode)
            {
                CKTdltNNum(ckt, here->BSIM4dNodePrime);
                here->BSIM4dNodePrime = 0;
            }
            if (here->BSIM4sNodePrime
                    && here->BSIM4sNodePrime != here->BSIM4sNode)
            {
                CKTdltNNum(ckt, here->BSIM4sNodePrime);
                here->BSIM4sNodePrime = 0;
            }
        }
    }
#endif
    return OK;
}
