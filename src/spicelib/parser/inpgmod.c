/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles, 1991 David A. Gates
Modified: 2001 Paolo Nenzi (Cider Integration)
**********/

#include "ngspice/ngspice.h"
#include "ngspice/inpdefs.h"
#include "ngspice/ifsim.h"
#include "ngspice/cpstd.h"
#include "ngspice/fteext.h"
#include "inpxx.h"
#include <errno.h>

#ifdef CIDER
/* begin Cider Integration */
#include "ngspice/numcards.h"
#include "ngspice/carddefs.h"
#include "ngspice/numgen.h"
#include "ngspice/suffix.h"

extern IFcardInfo *INPcardTab[];
extern int INPnumCards;
#define E_MISSING        -1
#define E_AMBIGUOUS        -2

static int INPparseNumMod( CKTcircuit* ckt, INPmodel *model, INPtables *tab, char **errMessage );
static int INPfindCard( char *name, IFcardInfo *table[], int numCards );
static int INPfindParm( char *name, IFparm *table, int numParms );

/* end Cider Integration */
#endif /* CIDER */

extern INPmodel *modtab;

/*
  code moved from INPgetMod
 */
static int
create_model( CKTcircuit* ckt, INPmodel* modtmp, INPtables* tab )
{
  IFvalue* val;
  char*    err = NULL;
  char*    line;
  char*    parm;
  char*    temp;
  int      error = 0;
  int      j;
  char *endptr; double dval;

  /* not already defined, so create & give parameters */
  error = ft_sim->newModel (ckt, modtmp->INPmodType, &(modtmp->INPmodfast), modtmp->INPmodName);

  if (error) return error;

  /* parameter isolation, identification, binding */

#ifdef CIDER
  /* begin cider integration */
  /* Handle Numerical Models Differently */
  if ( modtmp->INPmodType == INPtypelook("NUMD") ||
       modtmp->INPmodType == INPtypelook("NBJT") ||
       modtmp->INPmodType == INPtypelook("NUMD2") ||
       modtmp->INPmodType == INPtypelook("NBJT2") ||
       modtmp->INPmodType == INPtypelook("NUMOS") ) {
    error = INPparseNumMod( ckt, modtmp, tab, &err );
    if (error) return error;
  } else {
    /* It's an analytical model */
#endif /* CIDER */

    line = modtmp->INPmodLine->line;

#ifdef TRACE
    /* SDB debug statement */
    printf("In INPgetMod, inserting new model into table.  line = %s . . . \n", line);
#endif

    INPgetTok(&line, &parm, 1);        /* throw away '.model' */
    tfree(parm);
    INPgetTok(&line, &parm, 1);        /* throw away 'modname' */
    tfree(parm);
    while (*line != 0) {
      INPgetTok(&line, &parm, 1);
      if (!*parm)
        continue;

      for (j = 0; j < *(ft_sim->devices[modtmp->INPmodType]->numModelParms); j++) {

        if (strcmp(parm, "txl") == 0) {
          if (strcmp("cpl", ft_sim->devices[modtmp->INPmodType]->modelParms[j].keyword) == 0) {
            strcpy(parm, "cpl");
          }
        }

        if (strcmp(parm, ft_sim->devices[modtmp->INPmodType]->modelParms[j].keyword) == 0) {

          val = INPgetValue(ckt, &line, ft_sim->devices[modtmp->INPmodType]->modelParms[j].dataType, tab);

          error = ft_sim->setModelParm (ckt, modtmp->INPmodfast,
                                             ft_sim->devices[modtmp->INPmodType]->modelParms[j].id,
                                             val, NULL);
          if (error)
            return error;
          break;
        }
      } /* end for(j = 0 . . .*/

      if (strcmp(parm, "level") == 0) {
        /* just grab the level number and throw away */
        /* since we already have that info from pass1 */
        val = INPgetValue(ckt, &line, IF_REAL, tab);
      } else if (j >= *(ft_sim->devices[modtmp->INPmodType]->numModelParms)) {

        /* want only the parameter names in output - not the values */
        errno = 0;    /* To distinguish success/failure after call */
        dval = strtod(parm, &endptr);
        /* Check for various possible errors */
        if ((errno == ERANGE && dval == HUGE_VAL) || errno != 0) {
            perror("strtod");
            controlled_exit(EXIT_FAILURE);
        }
        if (endptr == parm) { /* it was no number - it is really a string */
          temp = tprintf("unrecognized parameter (%s) - ignored", parm);
          err = INPerrCat(err, temp);
        }
      }
      FREE(parm);
    }
#ifdef CIDER
  } /* analytical vs. numerical model parsing */
#endif
  modtmp->INPmodLine->error = err;

  return 0;
}

static bool
parse_line( char* line, char* tokens[], int num_tokens, double values[], bool found[] )
{
  char*     token = NULL;
  int       get_index = -1;
  int       i;
  bool      flag = TRUE;
  int       error;

  for ( i = 0; i < num_tokens; i++ )
    found[i] = FALSE;

  while( *line != '\0' ) {

    if ( get_index != -1 ) {
      values[get_index] = INPevaluate( &line, &error, 1 );
      found[get_index]  = TRUE;
      get_index         = -1;
      continue;
    }

    INPgetNetTok( &line, &token, 1 );

    for ( i = 0; i < num_tokens; i++ )
      if ( strcmp( tokens[i], token ) == 0 ) get_index = i;
    txfree(token);
  }

  for ( i = 0; i < num_tokens; i++ )
    if ( found[i] == FALSE ) {
      flag = FALSE;
      break;
    }

  return flag;
}

static bool
is_equal( double result, double expectedResult )
{
  if (fabs(result - expectedResult) < 1e-15) return TRUE;
  else                                       return FALSE;
}

static bool
in_range( double value, double min, double max )
{
  if ( (is_equal( value, min ) == TRUE) || /* the standard binning rule is: min <= value < max */
       (min < value && value < max) ) return TRUE;
  else                                return FALSE;
}

char*
INPgetModBin( CKTcircuit* ckt, char* name, INPmodel** model, INPtables* tab, char* line )
{
  INPmodel*    modtmp;
  double       l, w, lmin, lmax, wmin, wmax;
  double       parse_values[4];
  bool         parse_found[4];
  static char* instance_tokens[] = { "l", "w" };
  static char* model_tokens[]    = { "lmin", "lmax", "wmin", "wmax" };
  int          error;
  double       scale;
  char *err = NULL;

  if (!cp_getvar("scale", CP_REAL, &scale))
      scale = 1;

  *model = NULL;

  if ( parse_line( line, instance_tokens, 2, parse_values, parse_found ) != TRUE )
    return NULL;

  l = parse_values[0]*scale;
  w = parse_values[1]*scale;

  for ( modtmp = modtab; modtmp != NULL; modtmp = modtmp->INPnextModel ) {
    if ( /* This is the list of binable models */
           modtmp->INPmodType != INPtypelook ("BSIM3")
        && modtmp->INPmodType != INPtypelook ("BSIM3v32")
        && modtmp->INPmodType != INPtypelook ("BSIM3v0")
        && modtmp->INPmodType != INPtypelook ("BSIM3v1")
        && modtmp->INPmodType != INPtypelook ("BSIM4")
        && modtmp->INPmodType != INPtypelook ("BSIM4v5")
        && modtmp->INPmodType != INPtypelook ("BSIM4v6")
        && modtmp->INPmodType != INPtypelook ("BSIM4v7")
        && modtmp->INPmodType != INPtypelook ("HiSIM2")
        && modtmp->INPmodType != INPtypelook ("HiSIMHV")
       ) continue; /* We left the loop if the model is not in the list */

    if (modtmp->INPmodType < 0) {  /* First check for illegal model type */
      /* illegal device type, so can't handle */
      *model = NULL;
      err = tprintf("Unknown device type for model %s \n", name);
      return (err);
    } /* end of checking for illegal model */

    if ( parse_line( modtmp->INPmodLine->line, model_tokens, 4, parse_values, parse_found ) != TRUE )
      continue;

    lmin = parse_values[0]; lmax = parse_values[1];
    wmin = parse_values[2]; wmax = parse_values[3];

    if ( strncmp( modtmp->INPmodName, name, strlen( name ) ) == 0 &&
      in_range( l, lmin, lmax ) && in_range( w, wmin, wmax ) ) {
        if ( !modtmp->INPmodfast ) {
            error = create_model( ckt, modtmp, tab );
            if ( error ) return NULL;
        }
        *model = modtmp;
        return NULL;
    }
  }
  return NULL;
}

char *INPgetMod(CKTcircuit *ckt, char *name, INPmodel ** model, INPtables * tab)
{
  INPmodel *modtmp;
  char *err = NULL;
  int error;

#ifdef TRACE
  /* SDB debug statement */
  printf("In INPgetMod, examining model %s . . . \n", name);
#endif

  for (modtmp = modtab; modtmp != NULL; modtmp = modtmp->INPnextModel) {

#ifdef TRACE
    /* SDB debug statement */
    printf("In INPgetMod, comparing %s against stored model %s . . . \n", name, modtmp->INPmodName);
#endif

    if (strcmp(modtmp->INPmodName, name) == 0) {
      /* found the model in question - now instantiate if necessary */
      /* and return an appropriate pointer to it */

      if (modtmp->INPmodType < 0) {    /* First check for illegal model type */
        /* illegal device type, so can't handle */
        *model = NULL;
        err = tprintf("Unknown device type for model %s \n", name);

#ifdef TRACE
        /* SDB debug statement */
        printf("In INPgetMod, illegal device type for model %s . . . \n", name);
#endif

        return (err);
      }  /* end of checking for illegal model */

      if (! modtmp->INPmodfast) {   /* Check if model is already defined */
        error = create_model( ckt, modtmp, tab );
        if ( error ) return INPerror(error);
      }
      *model = modtmp;
      return (NULL);
    }
  }
  /* didn't find model - ERROR  - return model */
  *model = NULL;
  err = tprintf("Unable to find definition of model %s - default assumed \n", name);

#ifdef TRACE
  /* SDB debug statement */
  printf("In INPgetMod, didn't find model for %s, using default . . . \n", name);
#endif

  return (err);
}

#ifdef CIDER
/*
 * Parse a numerical model by running through the list of original
 * input cards which make up the model
 * Given:
 * 1. First card looks like: .model modname modtype <level=val>
 * 2. Other cards look like: +<whitespace>? where ? tells us what
 * to do with the next card:
 *    '#$*' = comment card
 *    '+'   = continue previous card
 *    other = new card
 */
static int
INPparseNumMod( CKTcircuit* ckt, INPmodel *model, INPtables *tab, char **errMessage )
{
    card *txtCard;        /* Text description of a card */
    GENcard *tmpCard;        /* Processed description of a card */
    IFcardInfo *info;        /* Info about the type of card located */
    char *line;
    char *cardName = NULL;        /* name of a card */
    char *parm;                /* name of a parameter */
    int cardType;        /* type/index for the current card */
    int cardNum = 0;        /* number of this card in the overall line */
    int lastType = E_MISSING;        /* type of previous card */
    char *err = NULL, *tmp;    /* Strings for error messages */
    IFvalue *value;
    int error, idx, invert;

    /* Chase down to the top of the list of actual cards */
    txtCard = model->INPmodLine->actualLine;

    /* Skip the first card if it exists since there's nothing interesting */
    /* txtCard will be empty if the numerical model is empty */
    if (txtCard) txtCard = txtCard->nextcard;

    /* Now parse each remaining card */
    while (txtCard) {
        line = txtCard->line;
        cardType = E_MISSING;
        cardNum++;

        /* Skip the initial '+' and any whitespace. */
        line++;
        while (*line == ' ' || *line == '\t')
            line++;

        switch (*line) {
        case '*':
        case '$':
        case '#':
        case '\0':
        case '\n':
            /* comment or empty cards */
            lastType = E_MISSING;
            break;
        case '+':
            /* continuation card */
            if (lastType >= 0) {
                cardType = lastType;
                while (*line == '+') line++;        /* Skip leading '+'s */
            } else {
                tmp = tprintf("Error on card %d : illegal continuation \'+\' - ignored",
                                cardNum);
                err = INPerrCat(err,tmp);
                lastType = E_MISSING;
                break;
            }
            /* FALL THRU when continuing a card */
        default:
            if (cardType == E_MISSING) {
                /* new command card */
                if (cardName) FREE(cardName);        /* get rid of old card name */
                INPgetTok(&line,&cardName,1);        /* get new card name */
                if (*cardName) {                 /* Found a name? */
                    cardType = INPfindCard(cardName,INPcardTab,INPnumCards);
                    if (cardType >= 0) {
                        /* Add card structure to model */
                        info = INPcardTab[cardType];
                        error = info->newCard (&tmpCard, model->INPmodfast );
                        if (error) return(error);
                    /* Handle parameter-less cards */
                    } else if (cinprefix( cardName, "title", 3 ) ) {
                        /* Do nothing */
                    } else if (cinprefix( cardName, "comment", 3 ) ) {
                        /* Do nothing */
                    } else if (cinprefix( cardName, "end", 3 ) ) {
                        /* Terminate parsing */
                        txtCard = ((card *) 0);
                        cardType = E_MISSING;
                    } else {
                        /* Error */
                        tmp = tprintf("Error on card %d : unrecognized name (%s) - ignored",
                                        cardNum, cardName);
                        err = INPerrCat(err,tmp);
                    }
                }
            }
            if (cardType >= 0) { /* parse the rest of this line */
                while (*line) {
                    /* Strip leading carat from booleans */
                    if (*line == '^') {
                      invert = TRUE;
                      line++;        /* Skip the '^' */
                    } else {
                      invert = FALSE;
                    }
                    INPgetTok(&line,&parm,1);
                    if (!*parm)
                      break;
                    idx = INPfindParm(parm, info->cardParms, info->numParms);
                    if (idx == E_MISSING) {
                        /* parm not found */
                        tmp = tprintf("Error on card %d : unrecognized parameter (%s) - ignored",
                                        cardNum, parm);
                        err = INPerrCat(err, tmp);
                    } else if (idx == E_AMBIGUOUS) {
                        /* parm ambiguous */
                        tmp = tprintf("Error on card %d : ambiguous parameter (%s) - ignored",
                                        cardNum, parm);
                        err = INPerrCat(err, tmp);
                    } else {
                        value = INPgetValue( ckt, &line,
                            info->cardParms[idx].dataType, tab );
                        if (invert) { /* invert if it's a boolean entry */
                          if ((info->cardParms[idx].dataType & IF_VARTYPES)
                              == IF_FLAG) {
                            value->iValue = 0;
                          } else {
                            tmp = tprintf("Error on card %d : non-boolean parameter (%s) - \'^\' ignored",
                                            cardNum, parm);
                            err = INPerrCat(err, tmp);
                          }
                        }
                        error = info->setCardParm (
                            info->cardParms[idx].id, value, tmpCard );
                        if (error) return(error);
                    }
                    FREE(parm);
                }
            }
            lastType = cardType;
            break;
        }
        if (txtCard) txtCard = txtCard->nextcard;
    }
    *errMessage = err;
    return( 0 );
}

/*
 * Locate the best match to a card name in an IFcardInfo table
 */
static int
INPfindCard( char *name, IFcardInfo *table[], int numCards )
{
    int test;
    int match, bestMatch;
    int best;
    int length;

    length = (int) strlen(name);

    /* compare all the names in the card table to this name */
    best = E_MISSING;
    bestMatch = 0;
    for ( test = 0; test < numCards; test++ ) {
        match = cimatch( name, table[test]->name );
        if ((match == bestMatch ) && (match > 0)){
            best = E_AMBIGUOUS;
        } else if ((match > bestMatch) && (match == length)) {
            best = test;
            bestMatch = match;
        }
    }
    return(best);
}

/*
 * Locate the best match to a parameter name in an IFparm table
 */
static int
INPfindParm( char *name, IFparm *table, int numParms )
{
    int test, best;
    int match, bestMatch;
    int id, bestId;
    int length;

    length = (int) strlen(name);

    /* compare all the names in the parameter table to this name */
    best = E_MISSING;
    bestId = -1;
    bestMatch = 0;
    for ( test = 0; test < numParms; test++ ) {
        match = cimatch( name, table[test].keyword );
        if ( (match == length) && (match == (int) strlen(table[test].keyword)) ) {
            /* exact match */
            best = test;
            /* all done */
            break;
        }
        id = table[test].id;
        if ((match == bestMatch) && (match > 0) && (id != bestId)) {
            best = E_AMBIGUOUS;
        } else if ((match > bestMatch) && (match == length)) {
            bestMatch = match;
            bestId = id;
            best = test;
        }
    }
    return(best);
}

#endif /* CIDER */
