/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2012 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 1  SUBVERSION : 2  REVISION : 4 )
 Model Parameter VERSION : 1.23
 FILE : hsmhvnoi.c

 DATE : 2013.04.30

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

#include "ngspice/ngspice.h"
#include "hsmhvdef.h"
#include "ngspice/cktdefs.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/noisedef.h"
#include "ngspice/suffix.h"
#include "ngspice/const.h"  /* jwan */
#include "hsmhvevalenv.h"
/* #include "hsmhvmacro.h" */

/*
 * HSMHVnoise (mode, operation, firstModel, ckt, data, OnDens)
 *    This routine names and evaluates all of the noise sources
 *    associated with MOSFET's.  It starts with the model *firstModel and
 *    traverses all of its insts.  It then proceeds to any other models
 *    on the linked list.  The total output noise density generated by
 *    all of the MOSFET's is summed with the variable "OnDens".
 */

int HSMHVnoise (
     int mode, int operation,
     GENmodel *inModel,
     CKTcircuit *ckt,
     register Ndata *data,
     double *OnDens)
{
  register HSMHVmodel *model = (HSMHVmodel *)inModel;
  register HSMHVinstance *here;
  char name[N_MXVLNTH];
  double tempOnoise=0.0 ;
  double tempInoise=0.0 ;
  double noizDens[HSMHVNSRCS] ;
  double lnNdens[HSMHVNSRCS] ;
  register int i;
  double G =0.0 ;
  double TTEMP = 0.0 ;

  /* define the names of the noise sources */
  static char * HSMHVnNames[HSMHVNSRCS] = {
    /* Note that we have to keep the order
       consistent with the index definitions 
       in hsmhvdefs.h */
    ".rd",              /* noise due to rd */
    ".rs",              /* noise due to rs */
    ".id",              /* noise due to id */
    ".1ovf",            /* flicker (1/f) noise */
    ".ign",             /* induced gate noise component at the drain node */
    ""                  /* total transistor noise */
  };
  
  for ( ;model != NULL; model = model->HSMHVnextModel ) {
    for ( here = model->HSMHVinstances; here != NULL;
	  here = here->HSMHVnextInstance ) {
      switch (operation) {
      case N_OPEN:
	/* see if we have to to produce a summary report */
	/* if so, name all the noise generators */
	  
	if (((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0) {
	  switch (mode) {
	  case N_DENS:
	    for ( i = 0; i < HSMHVNSRCS; i++ ) { 
	      (void) sprintf(name, "onoise.%s%s", 
			     here->HSMHVname, HSMHVnNames[i]);
	      data->namelist = TREALLOC(IFuid, data->namelist, data->numPlots + 1);
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, NULL);
	    }
	    break;
	  case INT_NOIZ:
	    for ( i = 0; i < HSMHVNSRCS; i++ ) {
	      (void) sprintf(name, "onoise_total.%s%s", 
			     here->HSMHVname, HSMHVnNames[i]);
	      data->namelist = TREALLOC(IFuid, data->namelist, data->numPlots + 1);
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, NULL);
	      
	      (void) sprintf(name, "inoise_total.%s%s", 
			     here->HSMHVname, HSMHVnNames[i]);
	      data->namelist = TREALLOC(IFuid, data->namelist, data->numPlots + 1);
	      if (!data->namelist)
		return(E_NOMEM);
	      (*(SPfrontEnd->IFnewUid)) 
		(ckt, &(data->namelist[data->numPlots++]),
		 (IFuid) NULL, name, UID_OTHER, NULL);
	    }
	    break;
	  }
	}
	break;
      case N_CALC:
	switch (mode) {
	case N_DENS:

         /* temperature */
         TTEMP = ckt->CKTtemp;
         if ( here->HSMHV_dtemp_Given ) { TTEMP = TTEMP + here->HSMHV_dtemp ; }
         TTEMP = TTEMP + *(ckt->CKTstate0 + here->HSMHVdeltemp) ;

         /* rs/rd thermal noise */
	  if ( model->HSMHV_corsrd == 1 || model->HSMHV_corsrd == 3 ) {
	    NevalSrc(&noizDens[HSMHVRDNOIZ], (double*) NULL,
		     ckt, N_GAIN,
		     here->HSMHVdNodePrime, here->HSMHVdNode,
		     (double) 0.0);
	    noizDens[HSMHVRDNOIZ] *= 4 * C_KB * TTEMP * here->HSMHVdrainConductance ;
            lnNdens[HSMHVRDNOIZ] = log( MAX(noizDens[HSMHVRDNOIZ],N_MINLOG) );
	    
	    NevalSrc(&noizDens[HSMHVRSNOIZ], (double*) NULL,
		     ckt, N_GAIN,
		     here->HSMHVsNodePrime, here->HSMHVsNode,
		     (double) 0.0);
	    noizDens[HSMHVRSNOIZ] *= 4 * C_KB * TTEMP * here->HSMHVsourceConductance ;
            lnNdens[HSMHVRSNOIZ] = log( MAX(noizDens[HSMHVRSNOIZ],N_MINLOG) );
	  } else {
	    noizDens[HSMHVRDNOIZ] = 0e0 ;
	    lnNdens[HSMHVRDNOIZ] = N_MINLOG ;
	    noizDens[HSMHVRSNOIZ] = 0e0 ;
	    lnNdens[HSMHVRSNOIZ] = N_MINLOG ;
	  }

	  /* channel thermal noise */
	  NevalSrc(&noizDens[HSMHVIDNOIZ], (double*) NULL,
		   ckt, N_GAIN,
		   here->HSMHVdNodePrime, here->HSMHVsNodePrime,
		   (double) 0.0);
	  switch( model->HSMHV_noise ) {
	  case 1:
	    /* HiSIMHV model */
	    G = here->HSMHV_noithrml ;
	    noizDens[HSMHVIDNOIZ] *= 4 * C_KB * TTEMP * G ;
	    lnNdens[HSMHVIDNOIZ] = log( MAX(noizDens[HSMHVIDNOIZ],N_MINLOG) );
	    break;
	  }

	  /* flicker noise */
	  NevalSrc(&noizDens[HSMHVFLNOIZ], (double*) NULL,
		   ckt, N_GAIN,
		   here->HSMHVdNodePrime, here->HSMHVsNodePrime, 
		   (double) 0.0);
	  switch ( model->HSMHV_noise ) {
	  case 1:
	    /* HiSIM model */
	    noizDens[HSMHVFLNOIZ] *= here->HSMHV_noiflick / pow(data->freq, model->HSMHV_falph) ; 
	    lnNdens[HSMHVFLNOIZ] = log(MAX(noizDens[HSMHVFLNOIZ], N_MINLOG));
	    break;
	  }	  

	  /* induced gate noise */
	  NevalSrc(&noizDens[HSMHVIGNOIZ], (double*) NULL,
		   ckt, N_GAIN, 
		   here->HSMHVdNodePrime, here->HSMHVsNodePrime, 
		   (double) 0.0);
	  switch ( model->HSMHV_noise ) {
	  case 1:
	    /* HiSIM model */
	    noizDens[HSMHVIGNOIZ] *= here->HSMHV_noiigate * here->HSMHV_noicross * here->HSMHV_noicross * data->freq * data->freq;
	    lnNdens[HSMHVIGNOIZ] = log(MAX(noizDens[HSMHVIGNOIZ], N_MINLOG));
	    break;
	  }

	  /* total */
	  noizDens[HSMHVTOTNOIZ] = noizDens[HSMHVRDNOIZ] + noizDens[HSMHVRSNOIZ]
	    + noizDens[HSMHVIDNOIZ] + noizDens[HSMHVFLNOIZ] + noizDens[HSMHVIGNOIZ];
	  lnNdens[HSMHVTOTNOIZ] = log(MAX(noizDens[HSMHVTOTNOIZ], N_MINLOG));
	  
	  *OnDens += noizDens[HSMHVTOTNOIZ];
	  
	  if ( data->delFreq == 0.0 ) {
	    /* if we haven't done any previous 
	       integration, we need to initialize our
	       "history" variables.
	    */
	    
	    for ( i = 0; i < HSMHVNSRCS; i++ ) 
	      here->HSMHVnVar[LNLSTDENS][i] = lnNdens[i];
	    
	    /* clear out our integration variables
	       if it's the first pass
	    */
	    if (data->freq == ((NOISEAN*) ckt->CKTcurJob)->NstartFreq) {
	      for (i = 0; i < HSMHVNSRCS; i++) {
		here->HSMHVnVar[OUTNOIZ][i] = 0.0;
		here->HSMHVnVar[INNOIZ][i] = 0.0;
	      }
	    }
	  }
	  else {
	    /* data->delFreq != 0.0,
	       we have to integrate.
	    */
	    for ( i = 0; i < HSMHVNSRCS; i++ ) {
	      if ( i != HSMHVTOTNOIZ ) {
		tempOnoise = 
		  Nintegrate(noizDens[i], lnNdens[i],
			     here->HSMHVnVar[LNLSTDENS][i], data);
		tempInoise = 
		  Nintegrate(noizDens[i] * data->GainSqInv, 
			     lnNdens[i] + data->lnGainInv,
			     here->HSMHVnVar[LNLSTDENS][i] + data->lnGainInv,
			     data);
		here->HSMHVnVar[LNLSTDENS][i] = lnNdens[i];
		data->outNoiz += tempOnoise;
		data->inNoise += tempInoise;
		if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
		  here->HSMHVnVar[OUTNOIZ][i] += tempOnoise;
		  here->HSMHVnVar[OUTNOIZ][HSMHVTOTNOIZ] += tempOnoise;
		  here->HSMHVnVar[INNOIZ][i] += tempInoise;
		  here->HSMHVnVar[INNOIZ][HSMHVTOTNOIZ] += tempInoise;
		}
	      }
	    }
	  }
	  if ( data->prtSummary ) {
	    for (i = 0; i < HSMHVNSRCS; i++) {
	      /* print a summary report */
	      data->outpVector[data->outNumber++] = noizDens[i];
	    }
	  }
	  break;
	case INT_NOIZ:
	  /* already calculated, just output */
	  if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
	    for ( i = 0; i < HSMHVNSRCS; i++ ) {
	      data->outpVector[data->outNumber++] = here->HSMHVnVar[OUTNOIZ][i];
	      data->outpVector[data->outNumber++] = here->HSMHVnVar[INNOIZ][i];
	    }
	  }
	  break;
	}
	break;
      case N_CLOSE:
	/* do nothing, the main calling routine will close */
	return (OK);
	break;   /* the plots */
      }       /* switch (operation) */
    }    /* for here */
  }    /* for model */
  
  return(OK);
}



