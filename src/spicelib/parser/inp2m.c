/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1988 Thomas L. Quarles
Modified: 2001 Paolo Nenzi (Cider Integration)
**********/

#include "ngspice/ngspice.h"
#include <stdio.h>
#include "ngspice/ifsim.h"
#include "ngspice/inpdefs.h"
#include "ngspice/inpmacs.h"
#include "ngspice/fteext.h"
#include "inpxx.h"

void
INP2M (CKTcircuit *ckt, INPtables * tab, card * current)
{

    /* Mname <node> <node> <node> <node> <model> [L=<val>]
     *       [W=<val>] [AD=<val>] [AS=<val>] [PD=<val>]
     *       [PS=<val>] [NRD=<val>] [NRS=<val>] [OFF]
     *       [IC=<val>,<val>,<val>]
     */

    int type;                  /* the type the model says it is */
    char *line;                /* the part of the current line left to parse */
    char *name;                /* the resistor's name */
    char *nname1;              /* the first node's name */
    char *nname2;              /* the second node's name */
    char *nname3;              /* the third node's name */
    char *nname4;              /* the fourth node's name */
    char *nname5;              /* the fifth node's name */
    char *nname6;              /* the sixt node's name */
    char *nname7;              /* the seventh node's name */
    char *save;                /* saj - used to save the posn of the start of
                                the parameters if the model is a mosfet*/
    CKTnode *node1;            /* the first node's node pointer */
    CKTnode *node2;            /* the second node's node pointer */
    CKTnode *node3;            /* the third node's node pointer */
    CKTnode *node4;            /* the fourth node's node pointer */
    CKTnode *node5;            /* the fifth node's node pointer */
    CKTnode *node6;            /* the sixth node's node pointer */
    CKTnode *node7;            /* the seventh node's node pointer */
    int error;                 /* error code temporary */
    int nodeflag;              /* flag indicating 4 or 5 (or 6 or 7) nodes */
    GENinstance *fast;         /* pointer to the actual instance */
    int waslead;               /* flag to indicate that funny unlabeled number was found */
    double leadval;            /* actual value of unlabeled number */
    char *model;               /* the name of the model */
    INPmodel *thismodel;       /* pointer to model description for user's model */
    GENmodel *mdfast;          /* pointer to the actual model */
    IFuid uid;                 /* uid for default model */
    char* err_msg;

#ifdef TRACE
    printf("INP2M: Parsing '%s'\n",current->line);
#endif

    nodeflag = 0;                  /*  initially specify a 4 terminal device  */
    line = current->line;
    INPgetTok (&line, &name, 1);
    INPinsert (&name, tab);
    INPgetNetTok (&line, &nname1, 1);
    INPtermInsert (ckt, &nname1, tab, &node1);
    INPgetNetTok (&line, &nname2, 1);
    INPtermInsert (ckt, &nname2, tab, &node2);
    INPgetNetTok (&line, &nname3, 1);
    INPtermInsert (ckt, &nname3, tab, &node3);
    INPgetNetTok (&line, &nname4, 1);
    INPtermInsert (ckt, &nname4, tab, &node4);

    node5 = NULL;
    node6 = NULL;
    node7 = NULL;

    /*  See if 5th token after device specification is a model name  */

    INPgetNetTok (&line, &nname5, 1);      /*  get 5th token  */
    save = line; /* saj - save the posn for later if
                  the default mosfet model is used */
    thismodel = NULL;
#ifdef TRACE
    printf("INP2M: checking for 4 node device\n");
#endif
    err_msg = INPgetMod (ckt, nname5, &thismodel, tab);
    tfree(err_msg);
    /* check if using model binning -- pass in line since need 'l' and 'w' */
    if ( thismodel == NULL ) {
        INPgetModBin( ckt, nname5, &thismodel, tab, line );
    }

    if (thismodel == NULL) {
        /*  5th token is not a model in the table  */
        nodeflag = 1;          /*  now specify a 5 node device  */
        INPgetNetTok (&line, &nname6, 1);   /*  get next token  */
        thismodel = NULL;
#ifdef TRACE
        printf("INP2M: checking for 5 node device\n");
#endif
        err_msg = INPgetMod (ckt, nname6, &thismodel, tab);
        tfree(err_msg);
        if (thismodel == NULL) {
            /*  6th token is not a model in the table  */
            nodeflag = 2;              /*  now specify a 6 node device  */
            INPgetNetTok (&line, &nname7, 1);       /*  get next token  */
            thismodel = NULL;
#ifdef TRACE
            printf("INP2M: checking for 6 node device\n");
#endif
            err_msg = INPgetMod (ckt, nname7, &thismodel, tab);
            tfree(err_msg);
            if (thismodel == NULL) {
                /*  7th token is not a model in the table  */
                nodeflag = 3;              /*  now specify a 7 node device  */
                INPgetTok (&line, &model, 1);           /*  get model name  */
#ifdef TRACE
                printf("INP2M: checking for 7 node device\n");
#endif
                err_msg = INPgetMod (ckt, model, &thismodel, tab);      /*  get pointer to the model  */
                tfree(err_msg);
                if (thismodel != NULL) {
                    if ((thismodel->INPmodType != INPtypelook ("B4SOI"))   &&
                            (thismodel->INPmodType != INPtypelook ("B3SOIPD")) &&
                            (thismodel->INPmodType != INPtypelook ("B3SOIFD")) &&
                            (thismodel->INPmodType != INPtypelook ("B3SOIDD"))
                       ) {
                        /*  if model is not variable node B3SOIPD/FD/DD model, error!  */
                        LITERR ("only level 55-58: B3SOI(PD|FD|DD) and B4SOI can have 7 nodes");
                        return;
                    } else {
                        /*  if looking at B3SOIPD/FD/DD or B4SOI model, allocate the 7th node  */
                        INPtermInsert (ckt, &nname5, tab, &node5);
                        INPtermInsert (ckt, &nname6, tab, &node6);
                        INPtermInsert (ckt, &nname7, tab, &node7);
                    }
                }
                /* saj - unbreak the default model creation*/
                else {
#ifdef TRACE
                    printf("INP2M: couldn't workout number of nodes, assuming 4\n");
#endif
                    nodeflag = 0;                  /*  now reset to a 4 node device  */
                    model = nname5;/* mosfet*/
                    line = save;   /* reset the posn to what it sould be */
                }
                /*saj*/
            } else {
                /*  7th token is a model - only have 6 terminal device  */
                if ((thismodel->INPmodType != INPtypelook ("B4SOI"))   &&
                        (thismodel->INPmodType != INPtypelook ("B3SOIPD")) &&
                        (thismodel->INPmodType != INPtypelook ("B3SOIFD")) &&
                        (thismodel->INPmodType != INPtypelook ("B3SOIDD")) &&
                        (thismodel->INPmodType != INPtypelook ("HiSIMHV")) &&
                        (thismodel->INPmodType != INPtypelook ("SOI3"))
                   ) {
                    /*  if model is not variable node B3SOIPD/FD/DD or STAG model, error!  */
                    LITERR ("only level 55-58,61,62: B3SOI(PD|FD|DD), B4SOI, STAG (SOI3) and HiSIMHV can have 6 nodes");
                    return;
                } else {
                    /*  if looking at B3SOIPD/FD/DD, B4SOI, STAG (SOI3) or HiSIMHV model, allocate the 6th node  */
                    INPtermInsert (ckt, &nname5, tab, &node5);
                    INPtermInsert (ckt, &nname6, tab, &node6);
                    model = nname7;
                }
            }
        } else {
            /*  6th token is a model - only have 5 terminal device  */
            if ((thismodel->INPmodType != INPtypelook ("B4SOI"))   &&
                    (thismodel->INPmodType != INPtypelook ("B3SOIPD")) &&
                    (thismodel->INPmodType != INPtypelook ("B3SOIFD")) &&
                    (thismodel->INPmodType != INPtypelook ("B3SOIDD")) &&
                    (thismodel->INPmodType != INPtypelook ("HiSIMHV")) &&
                    (thismodel->INPmodType != INPtypelook ("SOI3"))
               ) {
                /*  if model is not variable node B3SOIPD/FD/DD  model, error!  */
                LITERR ("only level 55-58,61,62: B3SOI(PD|FD|DD), B4SOI, STAG (SOI3) and HiSIMHV can have 5 nodes");
                return;
            } else {
                /*  if looking at B3SOIPD/FD/DD, B4SOI, STAG (SOI3) or HiSIMHV model, allocate the 5th node  */
                INPtermInsert (ckt, &nname5, tab, &node5);
                model = nname6;      /*  make model point to the correct token  */
            }
        }
    }

    else {
        /*  5th token is a model - only have 4 terminal device  */
        model = nname5;        /*  make model point to the correct token  */
    }


    INPinsert (&model, tab);
    thismodel = NULL;
#ifdef TRACE
    printf("INP2M: Looking up model\n");
#endif

    err_msg = INPgetMod (ckt, model, &thismodel, tab);
    if ( thismodel == NULL ) {
        INPgetModBin( ckt, model, &thismodel, tab, save );
        if ( thismodel == NULL )
            current->error = err_msg;
        else
            tfree(err_msg);
    }
    else
        tfree(err_msg);

    if (thismodel != NULL) {
        if (thismodel->INPmodType != INPtypelook ("Mos1")
                && thismodel->INPmodType != INPtypelook ("Mos2")
                && thismodel->INPmodType != INPtypelook ("Mos3")
                && thismodel->INPmodType != INPtypelook ("Mos5")
                && thismodel->INPmodType != INPtypelook ("Mos6")
                && thismodel->INPmodType != INPtypelook ("Mos8")
                && thismodel->INPmodType != INPtypelook ("Mos9")
                && thismodel->INPmodType != INPtypelook ("BSIM1")
                && thismodel->INPmodType != INPtypelook ("BSIM2")
                && thismodel->INPmodType != INPtypelook ("BSIM3")
                && thismodel->INPmodType != INPtypelook ("BSIM3v32")
                && thismodel->INPmodType != INPtypelook ("B4SOI")
                && thismodel->INPmodType != INPtypelook ("B3SOIPD")
                && thismodel->INPmodType != INPtypelook ("B3SOIFD")
                && thismodel->INPmodType != INPtypelook ("B3SOIDD")
                && thismodel->INPmodType != INPtypelook ("BSIM4")
                && thismodel->INPmodType != INPtypelook ("BSIM4v5")
                && thismodel->INPmodType != INPtypelook ("BSIM4v6")
                && thismodel->INPmodType != INPtypelook ("BSIM4v7")
                && thismodel->INPmodType != INPtypelook ("BSIM3v0")
                && thismodel->INPmodType != INPtypelook ("BSIM3v1")
                && thismodel->INPmodType != INPtypelook ("SOI3")
#ifdef CIDER
                && thismodel->INPmodType != INPtypelook ("NUMOS")
#endif
#ifdef ADMS
                && thismodel->INPmodType != INPtypelook ("ekv")
                && thismodel->INPmodType != INPtypelook ("psp102")
#endif
                && thismodel->INPmodType != INPtypelook ("HiSIM2")
                && thismodel->INPmodType != INPtypelook ("HiSIMHV")
           ) {
            LITERR ("incorrect model type");
            return;
        }
        type = thismodel->INPmodType;
        mdfast = (thismodel->INPmodfast);
    } else {
        type = INPtypelook ("Mos1");
        if (type < 0) {
            LITERR ("Device type MOS1 not supported by this binary\n");
            return;
        }
        if (!tab->defMmod) {
            /* create default M model */
            IFnewUid (ckt, &uid, NULL, "M", UID_MODEL, NULL);
            IFC (newModel, (ckt, type, &(tab->defMmod), uid));
        }
        mdfast = tab->defMmod;
    }
    IFC (newInstance, (ckt, mdfast, &fast, name));
    IFC (bindNode, (ckt, fast, 1, node1));
    IFC (bindNode, (ckt, fast, 2, node2));
    IFC (bindNode, (ckt, fast, 3, node3));
    IFC (bindNode, (ckt, fast, 4, node4));
    /* use type - not thismodel->INPmodType as it might not exist! */
    if ((type == INPtypelook ("B4SOI"))   ||
            (type == INPtypelook ("B3SOIPD")) ||
            (type == INPtypelook ("B3SOIFD")) ||
            (type == INPtypelook ("B3SOIDD")) ||
            (type == INPtypelook ("HiSIMHV")) ||
            (type == INPtypelook ("SOI3"))) {
        switch (nodeflag) {
        case 0:
            fast->GENnode5 = -1;
            fast->GENnode6 = -1;
            fast->GENnode7 = -1;
            break;
        case 1:
            IFC (bindNode, (ckt, fast, 5, node5));
            fast->GENnode6 = -1;
            fast->GENnode7 = -1;
            break;
        case 2:
            IFC (bindNode, (ckt, fast, 5, node5));
            IFC (bindNode, (ckt, fast, 6, node6));
            fast->GENnode7 = -1;
            break;
        case 3:
            IFC (bindNode, (ckt, fast, 5, node5));
            IFC (bindNode, (ckt, fast, 6, node6));
            IFC (bindNode, (ckt, fast, 7, node7));
            break;
        default:
            break;
        }
    }

    PARSECALL ((&line, ckt, type, fast, &leadval, &waslead, tab));
    if (waslead) {
        LITERR (" error:  no unlabeled parameter permitted on mosfet\n");
    }
}
