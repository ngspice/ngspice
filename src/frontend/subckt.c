/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Wayne A. Christopher, U. C. Berkeley CAD Group
Modified: 2000 AlansFixes
**********/

/*------------------------------------------------------------------------------
 * encapsulated string assembly in translate() and finishLine()
 *   this string facility (bxx_buffer) mainly abstracts away buffer allocation.
 *   this fixes a buffer overflow in finishLine, caused by lengthy descriptions
 *   of the kind:
 *     B1  1 2  I=v(1)+v(2)+v(3)+...
 * Larice, 22nd Aug 2009
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 * Added changes supplied by by H.Tanaka with some tidy up of comments, debug
 * statements, and variables. This fixes a problem with nested .subsck elements
 * that accessed .model lines. Code not ideal, but it seems to work okay.
 * Also took opportunity to tidy a few other items (unused variables etc.), plus
 * fix a few spelling errors in the comments, and a memory leak.
 * SJB 25th March 2005
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 * re-written by SDB during 4.2003 to enable SPICE2 POLY statements to be processed
 * properly.  This is particularly important for dependent sources, whose argument
 * list changes when POLY is used.
 * Major changes include:
 * -- Added lots of comments which (hopefully) elucidate the steps taken
 *    by the program during its processing.
 * -- Re-wrote translate, which does the processing of each card.
 * Please direct comments/questions/complaints to Stuart Brorson:
 * mailto:sdb@cloud9.net
 *-----------------------------------------------------------------------------*/

/*
 * Expand subcircuits. This is very spice-dependent. Bug fixes by Norbert
 * Jeske on 10/5/85.
 */

/*======================================================================*
 * Expand all subcircuits in the deck. This handles imbedded .subckt
 * definitions. The variables substart, subend, and subinvoke can be used
 * to redefine the controls used. The syntax is invariant though.
 * NOTE: the deck must be passed without the title line.
 * What we do is as follows: first make one pass through the circuit
 * and collect all of the subcircuits. Then, whenever a line that starts
 * with 'x' is found, copy the subcircuit associated with that name and
 * splice it in. A few of the problems: the nodes in the spliced-in
 * stuff must be unique, so when we copy it, append "subcktname:" to
 * each node. If we are in a nested subcircuit, use foo:bar:...:node.
 * Then we have to systematically change all references to the renamed
 * nodes. On top of that, we have to know how many args BJT's have,
 * so we have to keep track of model names.
 *======================================================================*/
/*#define TRACE*/
#include "ngspice/ngspice.h"
#include "ngspice/cpdefs.h"
#include "ngspice/ftedefs.h"
#include "ngspice/fteinp.h"

#include <stdarg.h>

#ifdef XSPICE
/* gtri - add - wbk - 11/9/90 - include MIF function prototypes */
#include "ngspice/mifproto.h"
/* gtri - end - wbk - 11/9/90 */
#endif

#include "subckt.h"
#include "variable.h"

#include "numparam/numpaif.h"

extern void line_free_x(struct line * deck, bool recurse);

#define line_free(line, flag)                   \
    do {                                        \
        line_free_x(line, flag);                \
        line = NULL;                            \
    } while(0)



struct subs;
static struct line *doit(struct line *deck, wordlist *modnames);
static int translate(struct line *deck, char *formal, char *actual, char *scname,
                     char *subname, struct subs *subs, wordlist const *modnames);
struct bxx_buffer;
static void finishLine(struct bxx_buffer *dst, char *src, char *scname);
static int settrans(char *formal, char *actual, char *subname);
static char *gettrans(const char *name, const char *name_end);
static int numnodes(char *name, struct subs *subs, wordlist const *modnames);
static int  numdevs(char *s);
static wordlist *modtranslate(struct line *deck, char *subname, wordlist  **const modnames);
static void devmodtranslate(struct line *deck, char *subname, wordlist * const orig_modnames);
static int inp_numnodes(char c);

/*---------------------------------------------------------------------
 * table is used in settrans and gettrans -- it holds the netnames used
 * in the .subckt definition (t_old), and in the subcircuit invocation
 * (t_new)
 *--------------------------------------------------------------------*/
static struct tab {
    char *t_old;
    char *t_new;
} table[512];   /* That had better be enough. */


/*---------------------------------------------------------------------
 *  subs is the linked list which holds the .subckt definitions
 *  found during processing.
 *--------------------------------------------------------------------*/
struct subs {
    char *su_name;          /* The .subckt name. */
    char *su_args;          /* The .subckt arguments, space separated. */
    int su_numargs;
    struct line *su_def;    /* Pointer to the .subckt definition. */
    struct subs *su_next;
};


/* orig_modnames is the list of original model names, modnames is the
 * list of translated names (i.e. after subckt expansion)
 */

static bool nobjthack = FALSE;
/* flag indicating use of the experimental numparams library */
static bool use_numparams = FALSE;

static char start[32], sbend[32], invoke[32], model[32];

static char node[128][128];
static int numgnode;

/*-------------------------------------------------------------------
  inp_subcktexpand is the top level function which translates
  .subckts into mainlined code.   Note that there are several things
  we need to do:  1. Find all .subckt definitions & stick them
  into a list.  2. Find all subcircuit invocations (refdes X)
  and replace them with the .subckt definition stored earlier.
  3. Do parameter substitution.

  The algorithm is as follows:
  1.  Define some aliases for .subckt, .ends, etc.
  2.  First numparam pass: substitute paramterized tokens by
  intermediate values 1000000001 etc.
  3.  Make a list node[] of global nodes
  4.  Clean up parens around netnames
  5.  Call doit, which does the actual translation.
  6.  Second numparam pass: Do final substitution
  7.  Check the results & return.
  inp_subcktexpand takes as argument a pointer to deck, and
  it returns a pointer to the same deck after the new subcircuits
  are spliced in.
  -------------------------------------------------------------------*/
struct line *
inp_subcktexpand(struct line *deck) {
    struct line *ll, *c;
    char *s;
    int ok = 0;
    char *t;
    int i;
    wordlist *modnames = NULL;

    if (!cp_getvar("substart", CP_STRING, start))
        (void) strcpy(start, ".subckt");
    if (!cp_getvar("subend", CP_STRING, sbend))
        (void) strcpy(sbend, ".ends");
    if (!cp_getvar("subinvoke", CP_STRING, invoke))
        (void) strcpy(invoke, "x");
    if (!cp_getvar("modelcard", CP_STRING, model))
        (void) strcpy(model, ".model");
    if (!cp_getvar("modelline", CP_STRING, model))
        (void) strcpy(model, ".model");
    nobjthack = cp_getvar("nobjthack", CP_BOOL, NULL);

    use_numparams = cp_getvar("numparams", CP_BOOL, NULL);

    use_numparams = TRUE;

    /*  deck has .control sections already removed, but not comments */
    if (use_numparams) {

#ifdef TRACE
        fprintf(stderr, "Numparams is processing this deck:\n");
        for (c = deck; c; c = c->li_next)
            fprintf(stderr, "%3d:%s\n", c->li_linenum, c->li_line);
#endif

        ok = nupa_signal(NUPADECKCOPY, NULL);
        /* get the subckt/model names from the deck */
        for (c = deck; c; c = c->li_next) {  /* first Numparam pass */
            if (ciprefix(".subckt", c->li_line))
                nupa_scan(c->li_line, c->li_linenum, TRUE);
            if (ciprefix(".model", c->li_line))
                nupa_scan(c->li_line, c->li_linenum, FALSE);
        }
        for (c = deck; c; c = c->li_next)  /* first Numparam pass */
            c->li_line = nupa_copy(c->li_line, c->li_linenum);
        /* now copy instances */

#ifdef TRACE
        fprintf(stderr, "Numparams transformed deck:\n");
        for (c = deck; c; c = c->li_next)
            fprintf(stderr, "%3d:%s\n", c->li_linenum, c->li_line);
#endif

    }

    /* Get all the model names so we can deal with BJTs, etc.
     *  Stick all the model names into the doubly-linked wordlist modnames.
     */
    {
        int nest = 0;
        for (c = deck; c; c = c->li_next) {

            if (ciprefix(".subckt", c->li_line))
                nest++;
            else if (ciprefix(".ends", c->li_line))
                nest--;
            else if (nest > 0)
                continue;

            if (ciprefix(model, c->li_line)) {
                s = c->li_line;
                txfree(gettok(&s)); /* discard the model keyword */
                modnames = wl_cons(gettok(&s), modnames);
            } /* model name finding routine */
        }
    }

#ifdef TRACE
    {
        wordlist *w;
        printf("Models found:\n");
        for (w = modnames; w; w = w->wl_next)
            printf("%s\n", w->wl_word);
    }
#endif

    /* Added by H.Tanaka to find global nodes */
    for (i = 0; i < 128; i++)
        strcpy(node[i], ""); /* Clear global node holder */

    for (c = deck; c; c = c->li_next)
        if (ciprefix(".global", c->li_line)) {
            s = c->li_line;
            txfree(gettok(&s));
            numgnode = 0;
            while (*s) {
                i = 0;
                t = s;
                for (/*s*/; *s && !isspace(*s); s++)
                    i++;
                strncpy(node[numgnode], t, (size_t) i);
                if (i>0 && t[i-1] != '\0')
                    node[numgnode][i] = '\0';
                while (isspace(*s))
                    s++;
                numgnode++;
            } /* node[] holds name of global node */
#ifdef TRACE
            printf("***Global node option has been found.***\n");
            for (i = 0; i<numgnode; i++)
                printf("***Global node no.%d is %s.***\n", i, node[i]);
            printf("\n");
#endif
            c->li_line[0] = '*'; /* comment it out */
        }

    /* Let's do a few cleanup things... Get rid of ( ) around node lists... */
    for (c = deck; c; c = c->li_next) {    /* iterate on lines in deck */

        char *s = c->li_line;

        if (*s == '*')           /* skip comment */
            continue;

        if (ciprefix(start, s)) {   /* if we find .subckt . . . */
#ifdef TRACE
            /* SDB debug statement */
            printf("In inp_subcktexpand, found a .subckt: %s\n", s);
#endif
            while (*s && *s != '(') /* search opening paren */
                s++;

            if (*s == '(') {
                int level = 0;
                do {
                    /* strip outer parens '(' ')', just the first pair */
                    if (*s == '('  &&  level++ == 0) {
                        *s = ' ';
                    }
                    if (*s == ')'  &&  --level == 0) {
                        *s = ' ';
                        break;
                    }
                } while(*s++);
            }
        } else if  (*s == '.') {
            continue;   /* skip .commands */
        } else {        /* any other line . . . */
            while (*s && !isspace(*s)) /* skip first token */
                s++;
            while (*s && isspace(*s)) /* skip whitespace */
                s++;

            if (*s == '(') {
                int level = 0;
                do {
                    /* strip outer parens '(' ')', just the first pair, why ? */
                    if (*s == '('  &&  level++ == 0) {
                        *s = ' ';
                    }
                    if (*s == ')'  &&  --level == 0) {
                        *s = ' ';
                        break;
                    }
                } while(*s++);
            } /* if (*s == '(' . . . */
        } /* any other line */
    }   /*  for (c = deck . . . */

#ifdef TRACE
    /* SDB debug statement */
    printf("In inp_subcktexpand, about to call doit.\n");
#endif

    /* doit does the actual splicing in of the .subckt . . .  */
    ll = doit(deck, modnames);

    wl_free(modnames);

    /* Count numbers of line in deck after expansion */
    if (ll) {
        dynMaxckt = 0; /* number of lines in deck after expansion */
        for (c = ll; c; c = c->li_next)
            dynMaxckt++;
    }

    /* Now check to see if there are still subckt instances undefined... */
    for (c = ll; c; c = c->li_next)
        if (ciprefix(invoke, c->li_line)) {
            fprintf(cp_err, "Error: unknown subckt: %s\n", c->li_line);
            if (use_numparams)
                ok = ok && nupa_signal(NUPAEVALDONE, NULL);
            return NULL;
        }

    if (use_numparams) {
        /* the NUMPARAM final line translation pass */
        ok = ok && nupa_signal(NUPASUBDONE, NULL);
        for (c = ll; c; c = c->li_next)
            /* 'param' .meas statements can have dependencies on measurement values */
            /* need to skip evaluating here and evaluate after other .meas statements */
            if (ciprefix(".meas", c->li_line)) {
                if (!strstr(c->li_line, "param"))
                    nupa_eval(c->li_line, c->li_linenum, c->li_linenum_orig);
            } else {
                nupa_eval(c->li_line, c->li_linenum, c->li_linenum_orig);
            }

#ifdef TRACE
        fprintf(stderr, "Numparams converted deck:\n");
        for (c = ll; c; c = c->li_next)
            fprintf(stderr, "%3d:%s\n", c->li_linenum, c->li_line);
#endif

        /*nupa_list_params(stdout);*/
        nupa_copy_inst_dico();
        ok = ok && nupa_signal(NUPAEVALDONE, NULL);
    }

    return (ll);  /* return the spliced deck.  */
}


#define MAXNEST 21
/*-------------------------------------------------------------------*/
/*  doit does the actual substitution of .subckts.                   */
/*  It takes two passes:  the first extracts .subckts                */
/*  and sticks pointer to them into the linked list sss.  It does    */
/*  the extraction recursively.  Then, it look for subcircuit        */
/*  invocations and substitutes the stored .subckt into              */
/*  the main circuit file.                                           */
/*  It takes as argument a pointer to the deck, and returns a        */
/*  pointer to the deck after the subcircuit has been spliced in.    */
/*-------------------------------------------------------------------*/
static struct line *
doit(struct line *deck, wordlist *modnames) {
    struct subs *sss = NULL;   /*  *sss temporarily hold decks to substitute  */
    int numpasses = MAXNEST;
    bool gotone;
    int error;

    /* Save all the old stuff... */
    struct subs *subs = NULL;
    wordlist *orig_modnames;
    wordlist *xmodnames = modnames;

#ifdef TRACE
    /* SDB debug statement */
    {
        struct line *c;
        printf("In doit, about to start first pass through deck.\n");
        for (c = deck; c; c = c->li_next)
            printf("   %s\n", c->li_line);
    }
#endif

    {
        /* First pass: xtract all the .subckts and stick pointers to them into sss.  */

        struct line *last = deck;
        struct line *lc   = NULL;

        while (last) {

            struct line *c, *lcc;

            if (ciprefix(sbend, last->li_line)) {         /* if line == .ends  */
                fprintf(cp_err, "Error: misplaced %s line: %s\n", sbend,
                        last->li_line);
                return (NULL);
            }

            if (ciprefix(start, last->li_line)) {  /* if line == .subckt  */
                if (last->li_next == NULL) {            /* first check that next line is non null */
                    fprintf(cp_err, "Error: no %s line.\n", sbend);
                    return (NULL);
                }

                /* Here we loop through the deck looking for .subckt and .ends cards.
                 * At the end of this section, last will point to the location of the
                 * .subckt card, and c will point to the location of the .ends card.
                 * and lcc->li_next === c, thus lcc will be the last body card
                 */
                {
                    int nest = 1;
                    lcc = last;
                    c = lcc->li_next;

                    while (c) {

                        if (ciprefix(sbend, c->li_line)) /* found a .ends */
                            nest--;
                        else if (ciprefix(start, c->li_line))  /* found a .subckt */
                            nest++;

                        if (!nest)
                            break;

                        lcc = c;
                        c = lcc->li_next;
                    }
                }

                /* Check to see if we have looped through remainder of deck without finding .ends */
                if (!c) {
                    fprintf(cp_err, "Error: no %s line.\n", sbend);
                    return (NULL);
                }

                /* last is the opening .subckt card */
                /* c    is the terminating .ends card */
                /* lcc  is one card before, which is the last body card */

                sss = alloc(struct subs);

                if (use_numparams == FALSE)
                    lcc->li_next = NULL;    /* shouldn't we free some memory here????? */

                /* cut the whole .subckt ... .ends sequence from the deck chain */
                if (lc)
                    lc->li_next = c->li_next;
                else
                    deck = c->li_next;

                /*  Now put the .subckt definition found into sss  */
                sss->su_def = last->li_next;

                {
                    char *s = last->li_line;
                    txfree(gettok(&s));
                    sss->su_name = gettok(&s);
                    sss->su_args = copy(s);
                    /* count the number of args in the .subckt line */
                    sss->su_numargs = 0;
                    for (;;) {
                        while (isspace(*s))
                            s++;
                        if (*s == '\0')
                            break;
                        while (*s && !isspace(*s))
                            s++;
                        sss->su_numargs ++;
                    }
                }

                sss->su_next = subs;
                subs = sss;            /* Now that sss is built, assign it to subs */

                line_free_x(last, FALSE);
                last = c->li_next;

                /*gp */
                c->li_next = NULL;  /* Numparam needs line c */
                c->li_line[0] = '*'; /* comment it out */
            } else { /*  line is neither .ends nor .subckt.  */
                /* make lc point to this card, and advance last to next card. */
                lc = last;
                last = last->li_next;
            }
        } /* for (last = deck . . . .  */
    }


    /* At this point, sss holds the .subckt definition found, subs holds
     * all .subckt defs found, including this one
     */

    if (!subs)            /* we have found no subckts.  Just return.  */
        return (deck);

    /* Otherwise, expand sub-subcircuits recursively. */
    for (sss = subs; sss; sss = sss->su_next)  /* iterate through the list of subcircuits */
        if ((sss->su_def = doit(sss->su_def, modnames)) == NULL)
            return (NULL);

#ifdef TRACE
    /* SDB debug statement */
    {
        struct line *c;
        printf("In doit, about to start second pass through deck.\n");
        for (c = deck; c; c = c->li_next)
            printf("   %s\n", c->li_line);
    }
#endif

    error = 0;
    /* Second pass: do the replacements. */
    do {                    /*  while (!error && numpasses-- && gotone)  */
        struct line *c = deck;
        struct line *lc = NULL;
        gotone = FALSE;
        while (c) {
            if (ciprefix(invoke, c->li_line)) {  /* found reference to .subckt (i.e. component with refdes X)  */

                char *tofree, *tofree2, *s, *t;
                char *scname, *subname;
                struct line *lcc;

                gotone = TRUE;
                t = tofree = s = copy(c->li_line);       /*  s & t hold copy of component line  */

                /*  make scname point to first non-whitepace chars after refdes invocation
                 * e.g. if invocation is Xreference, *scname = reference
                 */
                tofree2 = scname = gettok(&s);
                /*scname += strlen(invoke);   */
                while ((*scname == ' ') || (*scname == '\t') || (*scname == ':'))
                    scname++;

                /*  Now set s to point to last non-space chars in line (i.e.
                 *   the name of the model invoked
                 */
                while (*s)
                    s++;
                s--;
                while ((*s == ' ') || (*s == '\t'))
                    *s-- = '\0';
                while ((*s != ' ') && (*s != '\t'))
                    s--;
                s++;

                /* iterate through .subckt list and look for .subckt name invoked */
                for (sss = subs; sss; sss = sss->su_next)
                    if (eq(sss->su_name, s))
                        break;


                /* At this point, sss points to the .subckt invoked,
                 * and scname points to the netnames
                 * involved.
                 */


                /* If no .subckt is found, don't complain -- this might be an
                 * instance of a subckt that is defined above at higher level.
                 */
                if (!sss) {
                    lc = c;
                    c = c->li_next;
                    tfree(tofree);
                    tfree(tofree2);
                    continue;
                }

                /* Now we have to replace this line with the
                 * macro definition.
                 */
                subname = copy(sss->su_name);

                /*  make lcc point to a copy of the .subckt definition  */
                lcc = inp_deckcopy(sss->su_def);

                /* Change the names of .models found in .subckts . . .  */
                /* this translates the model name in the .model line */
                orig_modnames = modtranslate(lcc, scname, &modnames);
                if (orig_modnames)
                    devmodtranslate(lcc, scname, orig_modnames); /* This translates the model name on all components in the deck */
                wl_free(orig_modnames);

                {
                    char *s = sss->su_args;
                    txfree(gettok(&t));  /* Throw out the subcircuit refdes */

                    /* now invoke translate, which handles the remainder of the
                     * translation.
                     */
                    if (!translate(lcc, s, t, scname, subname, subs, modnames))
                        error = 1;
                    tfree(subname);
                }

                /* Now splice the decks together. */
                {
                    struct line *savenext = c->li_next;
                    if (use_numparams == FALSE) {
                        /* old style: c will drop a dangling pointer: memory leak  */
                        if (lc)
                            lc->li_next = lcc;
                        else
                            deck = lcc;
                    } else {
                        /* ifdef NUMPARAMS, keep the invoke line as a comment  */
                        c->li_next = lcc;
                        c->li_line[0] = '*'; /* comment it out */
                    }
                    while (lcc->li_next)
                        lcc = lcc->li_next;
                    lcc->li_next = c->li_next;
                    lcc->li_next = savenext;
                }
                c = lcc->li_next;
                lc = lcc;
                tfree(tofree);
                tfree(tofree2);
            }     /* if (ciprefix(invoke, c->li_line)) . . . */
            else {
                lc = c;
                c = c->li_next;
            }
        }
    } while (!error && numpasses-- && gotone);


    if (!numpasses) {
        fprintf(cp_err, "Error: infinite subckt recursion\n");
        error = 1;
    }

#ifdef TRACE
    /* Added by H.Tanaka to display converted deck */
    {
        struct line *c = deck;
        printf("Converted deck\n");
        for (; c; c = c->li_next)
            printf("%s\n", c->li_line);
    }
    {
        wordlist *w = modnames;
        printf("Models:\n");
        for (; w; w = w->wl_next)
            printf("%s\n", w->wl_word);
    }
#endif

    wl_delete_slice(modnames, xmodnames);

    if (error)
        return NULL;    /* error message already reported; should free() */


    /*
      struct subs {
      char *su_name;
      char *su_args;
      int su_numargs;
      struct line *su_def;
      struct subs *su_next;
      };
    */
    while (subs) {
        struct subs *sss2 = subs;
        subs = subs->su_next;
        tfree(sss2->su_name);
        tfree(sss2->su_args);
        line_free(sss2->su_def, TRUE);
        tfree(sss2);
    }

    return (deck);
}


/*-------------------------------------------------------------------*/
/* Copy a deck, including the actual lines.                          */
/*-------------------------------------------------------------------*/
struct line *
inp_deckcopy(struct line *deck) {
    struct line *d = NULL, *nd = NULL;

    while (deck) {
        if (nd) {
            d->li_next = alloc(struct line);
            d = d->li_next;
        } else {
            nd = d = alloc(struct line);
        }
        d->li_linenum = deck->li_linenum;
        d->li_line = copy(deck->li_line);
        if (deck->li_error)
            d->li_error = copy(deck->li_error);
        d->li_actual = inp_deckcopy(deck->li_actual);
        deck = deck->li_next;
    }
    return (nd);
}


/*-------------------------------------------------------------------
 * struct bxx_buffer,
 *   a string assembly facility.
 *
 * usage:
 *
 *   struct bxx_buffer thing;
 *   bxx_init(&thing);
 *   ...
 *   while (...) {
 *     bxx_rewind(&thing);
 *     ...
 *     bxx_putc(&thing, ...)
 *     bxx_printf(&thing, ...)
 *     bxx_put_cstring(&thing, ...)
 *     bxx_put_substring(&thing, ...)
 *     ...
 *     strcpy(bxx_buffer(&thing)
 *   }
 *   ..
 *   bxx_free(&thing)
 *
 * main aspect:
 *   reallocates/extends its buffer itself.
 *
 * note:
 *   during asssembly the internal buffer is
 *   not necessarily '\0' terminated.
 *   but will be when bxx_buffer() is invoked
 */

struct bxx_buffer {
    char *dst;
    char *limit;
    char *buffer;
};

/* must be a power of 2 */
static const int bxx_chunksize = 1024;

static void
bxx_init(struct bxx_buffer *t)
{
    /* assert(0 == (bxx_chunksize & (bxx_chunksize - 1))); */

    t->buffer = TMALLOC(char, bxx_chunksize);

    t->dst   = t->buffer;
    t->limit = t->buffer + bxx_chunksize;
}


static void
bxx_free(struct bxx_buffer *t)
{
    tfree(t->buffer);
}


static void
bxx_rewind(struct bxx_buffer *t)
{
    t->dst = t->buffer;
}


static void
bxx_extend(struct bxx_buffer *t, int howmuch)
{
    int pos = (int)(t->dst   - t->buffer);
    int len = (int)(t->limit - t->buffer);

    /* round up */
    howmuch +=  (bxx_chunksize - 1);
    howmuch &= ~(bxx_chunksize - 1);

    len += howmuch;

    t->buffer = TREALLOC(char, t->buffer, len);

    t->dst   = t->buffer + pos;
    t->limit = t->buffer + len;
}


static void
bxx_printf(struct bxx_buffer *t, const char *fmt, ...)
{
    va_list ap;

    for (;;) {
        int ret;
        int size = (int)(t->limit - t->dst);
        va_start(ap, fmt);
        ret  = vsnprintf(t->dst, (size_t) size, fmt, ap);
        va_end(ap);
        if (ret == -1) {
            bxx_extend(t, bxx_chunksize);
        } else if (ret >= size) {
            bxx_extend(t, ret - size + 1);
        } else {
            t->dst += ret;
            break;
        }
    }

    va_end(ap);
}


static inline char
bxx_putc(struct bxx_buffer *t, char c)
{
    if (t->dst >= t->limit)
        bxx_extend(t, 1);
    return *(t->dst)++ = c;
}


static void
bxx_put_cstring(struct bxx_buffer *t, const char *cstring)
{
    while (*cstring)
        bxx_putc(t, *cstring++);
}


static void
bxx_put_substring(struct bxx_buffer *t, const char *str, const char *end)
{
    while (str < end)
        bxx_putc(t, *str++);
}


static char *
bxx_buffer(struct bxx_buffer *t)
{
    if ((t->dst == t->buffer) || (t->dst[-1] != '\0'))
        bxx_putc(t, '\0');
    return t->buffer;
}


/*------------------------------------------------------------------------------------------*
 * Translate all of the device names and node names in the .subckt deck. They are
 * pre-pended with subname:, unless they are in the formal list, in which case
 * they are replaced with the corresponding entry in the actual list.
 * The one special case is node 0 -- this is always ground and we don't
 * touch it.
 *
 * Variable name meanings:
 * *deck = pointer to subcircuit definition (lcc) (struct line)
 * formal = copy of the .subckt definition line (e.g. ".subckt subcircuitname 1 2 3") (string)
 * actual = copy of the .subcircuit invocation line (e.g. "Xexample 4 5 6 subcircuitname") (string)
 * scname = refdes (- first letter) used at invocation (e.g. "example") (string)
 * subname = copy of the subcircuit name
 *-------------------------------------------------------------------------------------------*/
static int
translate(struct line *deck, char *formal, char *actual, char *scname, char *subname, struct subs *subs, wordlist const *modnames)
{
    struct line *c;
    struct bxx_buffer buffer;
    char *next_name, dev_type, *name, *s, *t, ch, *nametofree, *paren_ptr, *new_str;
    int nnodes, i, dim;
    int rtn = 0;

    bxx_init(&buffer);

    /* settrans builds the table holding the translated netnames.  */
    i = settrans(formal, actual, subname);
    if (i < 0) {
        fprintf(stderr,
                "Too few parameters for subcircuit type \"%s\" (instance: x%s)\n",
                subname, scname);
        goto quit;
    } else if (i > 0) {
        fprintf(stderr,
                "Too many parameters for subcircuit type \"%s\" (instance: x%s)\n",
                subname, scname);
        goto quit;
    }

    /* now iterate through the .subckt deck and translate the cards. */
    for (c = deck; c; c = c->li_next) {
        dev_type = *(c->li_line);

#ifdef TRACE
        /* SDB debug statement */
        printf("\nIn translate, examining line (dev_type: %c, subname: %s, instance: %s) %s \n", dev_type, subname, scname, c->li_line);
#endif

        if (ciprefix(".ic", c->li_line) || ciprefix(".nodeset", c->li_line)) {
            paren_ptr = s = c->li_line;
            while ((paren_ptr = strstr(paren_ptr, "("))  != NULL) {
                *paren_ptr = '\0';
                paren_ptr++;
                name = paren_ptr;

                if ((paren_ptr = strstr(paren_ptr, ")")) == NULL) {
                    *(name-1) = '(';
                    fprintf(cp_err, "Error: missing closing ')' for .ic|.nodeset statement %s\n", c->li_line);
                    goto quit;
                }
                *paren_ptr = '\0';
                t = gettrans(name, NULL);

                if (t) {
                    new_str = tprintf("%s(%s)%s", s, t, paren_ptr+1);
                } else {
                    new_str = tprintf("%s(%s.%s)%s", s, scname, name, paren_ptr+1);
                }

                paren_ptr = new_str + strlen(s) + 1;

                tfree(s);
                s = new_str;
            }
            c->li_line = s;
            continue;
        }

        /* Rename the device. */
        switch (dev_type) {
        case '\0':
        case '*':
        case '$':
        case '.':
            /* Just a pointer to the line into s and then break */
            bxx_rewind(&buffer);
            s = c->li_line;
            break;

#ifdef XSPICE
            /*===================  case A  ====================*/
            /* gtri - add - wbk - 10/23/90 - process A devices specially */
            /* since they have a more involved and variable length node syntax */

        case 'a':
        case 'A':

            /* translate the instance name according to normal rules */
            s = c->li_line;
            name = MIFgettok(&s);

            bxx_rewind(&buffer);
            bxx_printf(&buffer, "a.%s.%s ", scname, name);


            /* Now translate the nodes, looking ahead one token to recognize */
            /* when we reach the model name which should not be translated   */
            /* here.                                                         */

            next_name = MIFgettok(&s);

            for (;;) {
                /* rotate the tokens and get the the next one */
                if (name)
                    tfree(name);
                name = next_name;
                next_name = MIFgettok(&s);

                /* if next token is NULL, name holds the model name, so exit */
                if (next_name == NULL)
                    break;

                /* Process the token in name.  If it is special, then don't */
                /* translate it.                                            */
                switch (*name) {
                case '[':
                case ']':
                case '~':
                    bxx_printf(&buffer, "%s ", name);
                    break;

                case '%':
                    bxx_printf(&buffer, "%%");
                    /* don't translate the port type identifier */
                    if (name)
                        tfree(name);
                    name = next_name;
                    next_name = MIFgettok(&s);
                    bxx_printf(&buffer, "%s ", name);
                    break;

                default:

                    /* must be a node name at this point, so translate it */
                    t = gettrans(name, NULL);
                    if (t) {
                        bxx_printf(&buffer, "%s ", t);
                    } else {
                        bxx_printf(&buffer, "%s.%s ", scname, name);
                    }
                    break;

                } /* switch */

            } /* while */

            /* copy in the last token, which is the model name */
            if (name) {
                bxx_printf(&buffer, "%s ", name);
                tfree(name);
            }
            /* Set s to null string for compatibility with code */
            /* after switch statement                           */
            s = "";
            break; /* case 'a' */

            /* gtri - end - wbk - 10/23/90 */
#endif

            /*================   case E, F, G, H  ================*/
            /* This section handles controlled sources and allows for SPICE2 POLY attributes.
             * This is a new section, added by SDB to handle POLYs in sources.  Significant
             * changes were made in here.
             * 4.21.2003 -- SDB.  mailto:sdb@cloud9.net
             */
        case 'E':
        case 'e':
        case 'F':
        case 'f':
        case 'G':
        case 'g':
        case 'H':
        case 'h':

            s = c->li_line;       /* s now holds the SPICE line */
            t = name = gettok(&s);    /* name points to the refdes  */
            if (!name)
                continue;
            if (!*name) {
                tfree(name);
                continue;
            }

            /* Here's where we translate the refdes to e.g. F:subcircuitname:57
             * and stick the translated name into buffer.
             */
            ch = *name;           /* ch identifies the type of component */
            bxx_rewind(&buffer);
            bxx_printf(&buffer, "%c.%s.%s ", ch, scname, name);
            tfree(t);

            /* Next iterate over all nodes (netnames) found and translate them. */
            nnodes = numnodes(c->li_line, subs, modnames);

            while (nnodes-- > 0) {
                name = gettok_node(&s);
                if (name == NULL) {
                    fprintf(cp_err, "Error: too few nodes: %s\n",
                            c->li_line);
                    goto quit;
                }

                /* call gettrans and see if netname was used in the invocation */
                t = gettrans(name, NULL);

                if (t) {   /* the netname was used during the invocation; print it into the buffer */
                    bxx_printf(&buffer, "%s ", t);
                } else {
                    /* net netname was not used during the invocation; place a
                     * translated name into the buffer.*/
                    bxx_printf(&buffer, "%s.%s ", scname, name);
                }
                tfree(name);
            }  /* while (nnodes-- . . . . */


            /*  Next we handle the POLY (if any) */
            /* get next token */
            t = s;
            next_name = gettok_noparens(&t);
            if ((strcmp(next_name, "POLY") == 0) ||
                (strcmp(next_name, "poly") == 0)) {         /* found POLY . . . . */

#ifdef TRACE
                /* SDB debug statement */
                printf("In translate, looking at e, f, g, h found poly\n");
#endif

                /* move pointer ahead of (  */
                if (get_l_paren(&s) == 1) {
                    fprintf(cp_err, "Error: no left paren after POLY %s\n",
                            c->li_line);
                    tfree(next_name);
                    goto quit;
                }

                nametofree = gettok_noparens(&s);
                dim = atoi(nametofree);  /* convert returned string to int */
                tfree(nametofree);

                /* move pointer ahead of ) */
                if (get_r_paren(&s) == 1) {
                    fprintf(cp_err, "Error: no right paren after POLY %s\n",
                            c->li_line);
                    tfree(next_name);
                    goto quit;
                }

                /* Write POLY(dim) into buffer */
                bxx_printf(&buffer, "POLY( %d ) ", dim);

            } /* if ((strcmp(next_name, "POLY") == 0) . . .  */
            else
                dim = 1;    /* only one controlling source . . . */
            tfree(next_name);

            /* Now translate the controlling source/nodes */
            nnodes = dim * numdevs(c->li_line);
            while (nnodes-- > 0) {
                nametofree = name = gettok_node(&s);   /* name points to the returned token  */
                if (name == NULL) {
                    fprintf(cp_err, "Error: too few devs: %s\n", c->li_line);
                    goto quit;
                }

                if ((dev_type == 'f') ||
                    (dev_type == 'F') ||
                    (dev_type == 'h') ||
                    (dev_type == 'H')) {

                    /* Handle voltage source name */

#ifdef TRACE
                    /* SDB debug statement */
                    printf("In translate, found type f or h\n");
#endif

                    ch = *name;         /*  ch is the first char of the token.  */

                    bxx_printf(&buffer, "%c.%s.%s ", ch, scname, name);
                    /* From Vsense and Urefdes creates V.Urefdes.sense */
                } else {                            /* Handle netname */

#ifdef TRACE
                    /* SDB debug statement */
                    printf("In translate, found type e or g\n");
#endif

                    /* call gettrans and see if netname was used in the invocation */
                    t = gettrans(name, NULL);

                    if (t) {   /* the netname was used during the invocation; print it into the buffer */
                        bxx_printf(&buffer, "%s ", t);
                    } else {
                        /* net netname was not used during the invocation; place a
                         * translated name into the buffer.
                         */
                        bxx_printf(&buffer, "%s.%s ", scname, name);
                        /* From netname and Urefdes creates Urefdes:netname */
                    }
                }
                tfree(nametofree);
            }      /* while (nnodes--. . . . */

            /* Now write out remainder of line (polynomial coeffs) */
            finishLine(&buffer, s, scname);
            s = "";
            break;


            /*=================   Default case  ===================*/
        default:            /* this section handles ordinary components */
            s = c->li_line;
            nametofree = name = gettok_node(&s);  /* changed to gettok_node to handle netlists with ( , ) */
            if (!name)
                continue;
            if (!*name) {
                tfree(name);
                continue;
            }

            /* Here's where we translate the refdes to e.g. R:subcircuitname:57
             * and stick the translated name into buffer.
             */
            ch = *name;

            bxx_rewind(&buffer);

            if (ch != 'x')
                bxx_printf(&buffer, "%c.%s.%s ", ch, scname, name);
            else
                bxx_printf(&buffer, "%s.%s ", scname, name);

            tfree(nametofree);

            /* Next iterate over all nodes (netnames) found and translate them. */
            nnodes = numnodes(c->li_line, subs, modnames);
            while (nnodes-- > 0) {
                name = gettok_node(&s);
                if (name == NULL) {
                    fprintf(cp_err, "Error: too few nodes: %s\n", c->li_line);
                    goto quit;
                }

                /* call gettrans and see if netname was used in the invocation */
                t = gettrans(name, NULL);

                if (t) {   /* the netname was used during the invocation; print it into the buffer */
                    bxx_printf(&buffer, "%s ", t);
                } else {
                    /* net netname was not used during the invocation; place a
                     * translated name into the buffer.
                     */
                    bxx_printf(&buffer, "%s.%s ", scname, name);
                }
                tfree(name);
            }  /* while (nnodes-- . . . . */

            /* Now translate any devices (i.e. controlling sources).
             * This may be superfluous because we handle dependent
             * source devices above . . . .
             */
            nnodes = numdevs(c->li_line);
            while (nnodes-- > 0) {
                t = name = gettok_node(&s);
                if (name == NULL) {
                    fprintf(cp_err, "Error: too few devs: %s\n", c->li_line);
                    goto quit;
                }
                ch = *name;

                if (ch != 'x')
                    bxx_printf(&buffer, "%c.%s.%s ", ch, scname, name);
                else
                    bxx_printf(&buffer, "%s ", scname);

                tfree(t);
            } /* while (nnodes--. . . . */

            /* Now we finish off the line.  For most components (R, C, etc),
             * this involves adding the component value to the buffer.
             * We also scan through the line for v(something) and
             * i(something)...
             */
            finishLine(&buffer, s, scname);
            s = "";

        } /* switch (c->li_line . . . . */

        bxx_printf(&buffer, "%s", s);
        tfree(c->li_line);
        c->li_line = copy(bxx_buffer(&buffer));

#ifdef TRACE
        /* SDB debug statement */
        printf("In translate, translated line = %s \n", c->li_line);
#endif

    }  /* for (c = deck . . . . */
    rtn = 1;
quit:
    for (i = 0; ; i++) {
        if (!table[i].t_old && !table[i].t_new)
            break;
        FREE(table[i].t_old);
        FREE(table[i].t_new);
    }

    bxx_free(&buffer);
    return rtn;
}


/*-------------------------------------------------------------------*
 * finishLine now doesn't handle current or voltage sources.
 * Therefore, it just writes out the final netnames, if required.
 * Changes made by SDB on 4.29.2003.
 *-------------------------------------------------------------------*/
static void
finishLine(struct bxx_buffer *t, char *src, char *scname)
{
    char *buf, *buf_end, which;
    char *s;
    int lastwasalpha;

    lastwasalpha = 0;
    while (*src) {
        /* Find the next instance of "<non-alpha>[vi]<opt spaces>(" in
         * this string.
         */
        if (((*src != 'v') && (*src != 'V') &&
             (*src != 'i') && (*src != 'I')) ||
            lastwasalpha) {
            lastwasalpha = isalpha(*src);
            bxx_putc(t, *src++);
            continue;
        }
        for (s = src + 1; *s && isspace(*s); s++)
            ;
        if (!*s || (*s != '(')) {
            lastwasalpha = isalpha(*src);
            bxx_putc(t, *src++);
            continue;
        }
        lastwasalpha = 0;
        bxx_putc(t, which = *src);
        src = s;
        bxx_putc(t, *src++);
        while (isspace(*src))
            src++;
        for (buf = src; *src && !isspace(*src) && *src != ',' && *src != ')'; )
            src++;
        buf_end = src;

        if ((which == 'v') || (which == 'V'))
            s = gettrans(buf, buf_end);
        else
            s = NULL;

        if (s) {
            bxx_put_cstring(t, s);
        } else {  /* just a normal netname . . . . */
            /*
              i(vname) -> i(v.subckt.vname)
              i(ename) -> i(e.subckt.ename)
              i(hname) -> i(h.subckt.hname)
              i(bname) -> i(b.subckt.hname)
            */
            if ((which == 'i' || which == 'I') &&
                (buf[0] == 'v' || buf[0] == 'V' || buf[0] == 'e' || buf[0] == 'h'
                 || buf[0] == 'b' || buf[0] == 'B')) {
                bxx_putc(t, buf[0]);
                bxx_putc(t, '.');
                /*i = 1; */
            } /* else {
                 i = 0;
                 } */
            bxx_put_cstring(t, scname);
            bxx_putc(t, '.');
            bxx_put_substring(t, buf, buf_end);
        }

        /* translate the reference node, as in the "2" in "v(4,2)" */

        if ((which == 'v') || (which == 'V')) {
            while (*src && (isspace(*src) || *src == ',')) {
                src++;
            }
            if (*src && *src != ')') {
                for (buf = src; *src && !isspace(*src) && (*src != ')'); )
                    src++;
                s = gettrans(buf, buf_end = src);
                bxx_putc(t, ',');
                if (s) {
                    bxx_put_cstring(t, s);
                } else {
                    bxx_put_cstring(t, scname);
                    bxx_putc(t, '.');
                    bxx_put_substring(t, buf, buf_end);
                }
            }
        }
    }
}


/*------------------------------------------------------------------------------*
 * settrans builds the table which holds the old and new netnames.
 * it also compares the number of nets present in the .subckt definition against
 * the number of nets present in the subcircuit invocation.  It returns 0 if they
 * match, otherwise, it returns an error.
 *
 * Variable definitions:
 * formal = copy of the .subckt definition line (e.g. ".subckt subcircuitname 1 2 3") (string)
 * actual = copy of the .subcircuit invocation line (e.g. "Xexample 4 5 6 subcircuitname") (string)
 * subname = copy of the subcircuit name
 *------------------------------------------------------------------------------*/
static int
settrans(char *formal, char *actual, char *subname)
{
    int i;

    bzero(table, sizeof(*table));

    for (i = 0; ; i++) {
        table[i].t_old = gettok(&formal);
        table[i].t_new = gettok(&actual);

        if (table[i].t_new == NULL) {
            return -1;          /* Too few actual / too many formal */
        } else if (table[i].t_old == NULL) {
            if (eq(table[i].t_new, subname))
                break;
            else
                return 1;       /* Too many actual / too few formal */
        }
    }
    return 0;
}


/* compare a substring, with a '\0' terminated string
 *   the substring itself is required to be free of a '\0'
 */

static int
eq_substr(const char *str, const char *end, const char *cstring)
{
    while (str < end)
        if (*str++ != *cstring++)
            return 0;
    return (*cstring == '\0');
}


/*------------------------------------------------------------------------------*
 * gettrans returns the name of the top level net if it is in the list,
 * otherwise it returns NULL.
 *------------------------------------------------------------------------------*/
static char *
gettrans(const char *name, const char *name_end)
{
    int i;

    if (!name_end)
        name_end = strchr(name, '\0');

#ifdef XSPICE
    /* gtri - wbk - 2/27/91 - don't translate the reserved word 'null' */
    if (eq_substr(name, name_end, "null"))
        return ("null");
    /* gtri - end */
#endif

    /* Added by H.Tanaka to translate global nodes */
    for (i = 0; i<numgnode; i++)
        if (eq_substr(name, name_end, node[i]))
            return (node[i]);

    if (eq_substr(name, name_end, "0"))
        return ("0");

    for (i = 0; table[i].t_old; i++)
        if (eq_substr(name, name_end, table[i].t_old))
            return (table[i].t_new);

    return (NULL);
}


/*
  check if current token matches model bin name -- <token>.[0-9]+
*/
static bool
model_bin_match(char *token, char *model_name)
{
    /* find last dot in model_name */
    char *dot_char = strrchr(model_name, '.');
    bool  flag = FALSE;
    /* check if token equals the substring before last dot in model_name */
    if (dot_char) {
        char *mtoken = copy_substring(model_name, dot_char);
        if (cieq(mtoken, token)) {
            flag = TRUE;
            dot_char++;
            /* check if model_name has binning info (trailing digit(s)) */
            while (*dot_char != '\0') {
                if (!isdigit(*dot_char)) {
                    flag = FALSE;
                    break;
                }
                dot_char++;
            }
        }
        tfree(mtoken);
    }
    return flag;
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
static int
numnodes(char *name, struct subs *subs, wordlist const *modnames)
{
    /* gtri - comment - wbk - 10/23/90 - Do not modify this routine for */
    /* 'A' type devices since the callers will not know how to find the */
    /* nodes even if they know how many there are.  Modify the callers  */
    /* instead.                                                         */
    /* gtri - end - wbk - 10/23/90 */
    char c;
    struct subs *sss;
    char *s, *t, buf[4 * BSIZE_SP];
    const wordlist *wl;
    int n, i, gotit;

    while (*name && isspace(*name))
        name++;

    c = *name;
    if (isupper(c))
        c = (char) tolower(c);

    (void) strncpy(buf, name, sizeof(buf));
    s = buf;
    if (c == 'x') {     /* Handle this ourselves. */
        while (*s)
            s++;
        s--;
        while ((*s == ' ') || (*s == '\t'))
            *s-- = '\0';
        while ((*s != ' ') && (*s != '\t'))
            s--;
        s++;
        for (sss = subs; sss; sss = sss->su_next)
            if (eq(sss->su_name, s))
                return (sss->su_numargs);
        /*
         * number of nodes not known so far.
         * lets count the nodes ourselves,
         * assuming `buf' looks like this:
         *    xname n1 n2 ... nn subname
         */
        {
            int nodes = -2;
            for (s = buf; *s; ) {
                nodes++;
                while (*s && !isspace(*s))
                    s++;
                while (isspace(*s))
                    s++;
            }
            return (nodes);
        }
    }

    n = inp_numnodes(c);

    /* Added this code for variable number of nodes on BSIM3SOI/CPL devices  */
    /* The consequence of this code is that the value returned by the    */
    /* inp_numnodes(c) call must be regarded as "maximum number of nodes */
    /* for a given device type.                                          */
    /* Paolo Nenzi Jan-2001                                              */

    /* I hope that works, this code is very very untested */

    if ((c == 'm') || (c == 'p')) {              /* IF this is a mos or cpl */
        i = 0;
        s = buf;
        gotit = 0;
        txfree(gettok(&s));          /* Skip component name */
        while ((i < n) && (*s) && !gotit) {
            t = gettok_node(&s);       /* get nodenames . . .  */
            for (wl = modnames; wl; wl = wl->wl_next) {
                /* also need to check if binnable device mos model */
                if (eq(t, wl->wl_word) || model_bin_match(t, wl->wl_word))
                    gotit = 1;
            }
            i++;
            tfree(t);
        } /* while . . . . */

        /* Note: node checks must be done on #_of_node-1 because the */
        /* "while" cycle increments the counter even when a model is */
        /* recognized. This code may be better!                      */

        if (i < 5) {
            fprintf(cp_err, "Error: too few nodes for MOS or CPL: %s\n", name);
            return (0);
        }
        return (i-1); /* compensate the unnecessary increment in the while cycle */
    } /* if (c == 'm' . . .  */

    if (nobjthack || (c != 'q'))
        return (n);

    for (s = buf, i = 0; *s && (i < 4); i++)
        txfree(gettok(&s));

    if (i == 3)
        return (3);
    else if (i < 4) {
        fprintf(cp_err, "Error: too few nodes for BJT: %s\n", name);
        return (0);
    }

    /* Now, is this a model? */
    t = gettok(&s);
    wl = wl_find(t, modnames);
    tfree(t);
    return wl ? 3 : 4;
}


/*-------------------------------------------------------------------*
 *  This function returns the number of controlling voltage sources
 *  (for F, H) or controlling nodes (for G, E)  attached to a dependent
 *  source.
 *-------------------------------------------------------------------*/
static int
numdevs(char *s)
{

    while (*s && isspace(*s))
        s++;
    switch (*s) {
    case 'K':
    case 'k':
        return (2);

        /* two nodes per voltage controlled source */
    case 'G':
    case 'g':
    case 'E':
    case 'e':
        return (2);

        /* one source per current controlled source */
    case 'F':
    case 'f':
    case 'H':
    case 'h':
        /* 2 lines here added to fix w bug, NCF 1/31/95 */
    case 'W':
    case 'w':
        return (1);

    default:
        return (0);
    }
}


/*----------------------------------------------------------------------*
 *  modtranslate --  translates .model lines found in subckt definitions.
 *  Calling arguments are:
 *  *c = pointer to the .subckt definition (linked list)
 *  *subname = pointer to the subcircuit name used at the subcircuit invocation (string)
 *  modtranslate returns the list of model names which have been translated
 *----------------------------------------------------------------------*/
static wordlist *
modtranslate(struct line *c, char *subname, wordlist ** const modnames)
{
    wordlist *orig_modnames = NULL;

    for (; c; c = c->li_next)
        if (ciprefix(".model", c->li_line)) {
            char *model_name, *new_model_name;
            char *t = c->li_line;

#ifdef TRACE
            printf("modtranslate(), translating:\n"
                   "  \"%s\" -->\n", t);
#endif

            /* swallow ".model" */
            txfree(gettok(&t));

            model_name = gettok(&t);

            new_model_name = tprintf("%s:%s", subname, model_name);

            /* remember the translation */
            orig_modnames = wl_cons(model_name, orig_modnames);
            *modnames = wl_cons(new_model_name, *modnames);

            /* perform the actual translation of this .model line */
            t = tprintf(".model %s %s", new_model_name, t);
            tfree(c->li_line);
            c->li_line = t;

#ifdef TRACE
            printf("  \"%s\"\n", t);
            printf("  mapped modelname \"%s\" --> \"%s\"\n",
                   model_name, new_model_name);
#endif

        }

    return orig_modnames;
}


/*-------------------------------------------------------------------*
 *  Devmodtranslate scans through the deck, and translates the
 *  name of the model in a line held in a .subckt.  For example:
 *  before:   .subckt U1 . . . .
 *            Q1 c b e 2N3904
 *  after:    Q1 c b e U1:2N3904
 *-------------------------------------------------------------------*/
static void
devmodtranslate(struct line *deck, char *subname, wordlist * const orig_modnames)
{
    struct line *s;
    int found;

    for (s = deck; s; s = s->li_next) {

        char *buffer, *t, c, *name, *next_name;
        wordlist *wlsub;

        t = s->li_line;

#ifdef TRACE
        /* SDB debug stuff */
        printf("In devmodtranslate, examining line %s.\n", t);
#endif

        while (*t && isspace(*t))
            t++;
        c = *t;                           /* set c to first char in line. . . . */
        if (isupper(c))
            c = (char) tolower(c);

        buffer = TMALLOC(char, strlen(t) + strlen(subname) + 4);

        switch (c) {

#ifdef XSPICE

        case 'a':
            /*  Code for codemodels (dev prefix "A") added by SDB on 6.10.2004.
             *  The algorithm is simple.  We don't know how many nodes or sources are attached,
             *  but the name of the model is always last.  Therefore, just iterate through all
             *  tokens until the last one is reached.  Then translate it.
             */

#ifdef TRACE
            /* SDB debug statement */
            printf("In devmodtranslate, found codemodel, line= %s\n", t);
#endif

            /* first do refdes. */
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);

            /* now do remainder of line. */
            next_name = gettok(&t);
            for (;;) {
                name = next_name;
                next_name = gettok(&t);

                if (next_name == NULL) {
                    /* if next_name is NULL, we are at the line end.
                     * name holds the model name.  Therefore, break */
                    break;

                } else {
                    /* next_name holds something.  Write name into the buffer and continue. */
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                    tfree(name);
                }
            }  /* while  */


            /*
             *  Note that we compare against orig_modnames,
             *    which is the list of untranslated names of models.
             */
            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            tfree(name);

#ifdef TRACE
            /* SDB debug statement */
            printf("In devmodtranslate, translated codemodel line= %s\n", buffer);
#endif

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

#endif /* XSPICE */

        case 'r':
        case 'c':
        case 'l':
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get first netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get second netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);

            if (*t) {    /* if there is a model, process it. . . . */
                name = gettok(&t);
                wlsub = wl_find(name, orig_modnames);

                if (!wlsub)
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                else
                    (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);
                tfree(name);
            }

            if (*t) {
                name = gettok(&t);
                wlsub = wl_find(name, orig_modnames);

                if (!wlsub)
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                else
                    (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);
                tfree(name);
            }

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

        case 'd':
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get first attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get second attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);

            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            tfree(name);
            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

        case 'u': /* urc transmissionline */
            /* 3 terminal devices */
        case 'w': /* current controlled switch */
        case 'j': /* jfet */
        case 'z': /* hfet, mesa */
            name = gettok(&t);
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok(&t);
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);

            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            tfree(name);
            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

            /* 4 terminal devices */
        case 'o':    /* ltra */
        case 's':    /* vc switch */
        case 'y':    /* txl */
            /*  Changed gettok() to gettok_node() on 12.2.2003 by SDB
                to enable parsing lines like "S1 10 11 (80,51) SLATCH1"
                which occur in real Analog Devices SPICE models.
            */
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get first attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get second attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get third attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get fourth attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);

            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            tfree(name);
            break;

             /* 4-7 terminal mos devices */
        case 'm':
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get first attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get second attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get third attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get fourth attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok(&t);

            found = 0;
            while (!found) {
                /* Now, is this a subcircuit model? */
                for (wlsub = orig_modnames; wlsub; wlsub = wlsub->wl_next) {
                    /* FIXME, probably too unspecific */
                    int i = (int) strlen(wlsub->wl_word);
                    int j = 0; /* Now, have we a binned model? */
                    char* dot_char;
                    if ((dot_char = strstr(wlsub->wl_word, ".")) != NULL) {
                        dot_char++;
                        j++;
                        while (*dot_char != '\0') {
                            if (!isdigit(*dot_char)) {
                                break;
                            }
                            dot_char++;
                            j++;
                        }
                    }
                    if (strncmp(name, wlsub->wl_word, (size_t) (i - j)) == 0) {
                        found = 1;
                        break;
                    }
                }
                if (!found) { /* name was not a model - was a netname */
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                    tfree(name);
                    name = gettok(&t);
                    if (name == NULL) {
                        name = copy(""); /* allow 'tfree' */
                        break;
                    }
                }
            }  /* while  */

            if (!found)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            tfree(name);
            break;

            /* 3-5 terminal bjt devices */
        case 'q':
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get first attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get second attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* get third attached netname */
            (void) sprintf(buffer + strlen(buffer), "%s ", name);
            tfree(name);
            name = gettok_node(&t);  /* this can be either a model name or a node name. */

            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                if (*t) { /* There is another token - perhaps a model */
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                    tfree(name);
                    name = gettok(&t);
                    wlsub = wl_find(name, orig_modnames);
                }

#ifdef ADMS
            if (!wlsub)
                if (*t) { /* There is another token - perhaps a model */
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                    tfree(name);
                    name = gettok(&t);
                    wlsub = wl_find(name, orig_modnames);
                }
#endif

            if (!wlsub) /* Fallback w/o subckt name before */
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            tfree(name);

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

            /* 4-18 terminal devices */
        case 'p': /* cpl */
            name = gettok(&t);  /* get refdes */
            (void) sprintf(buffer, "%s ", name);
            tfree(name);

            /* now do remainder of line. */
            next_name = gettok(&t);
            for (;;) {
                name = next_name;
                next_name = gettok(&t);
                if (!next_name || strstr(next_name, "len")) {
                    /* if next_name is NULL or len or length, we are at the line end.
                     * name holds the model name.  Therefore, break */
                    break;
                } else {
                    /* next_name holds something.  Write name into the buffer and continue. */
                    (void) sprintf(buffer + strlen(buffer), "%s ", name);
                    tfree(name);
                }
            }  /* while  */

            wlsub = wl_find(name, orig_modnames);

            if (!wlsub)
                (void) sprintf(buffer + strlen(buffer), "%s ", name);
            else
                (void) sprintf(buffer + strlen(buffer), "%s:%s ", subname, name);

            tfree(name);

            (void) strcat(buffer, t);
            tfree(s->li_line);
            s->li_line = buffer;
            break;

        default:
            tfree(buffer);
            break;
        }
    }
}


/*----------------------------------------------------------------------*
 * inp_numnodes returns the number of nodes (netnames) attached to the
 * component.
 * This is a spice-dependent thing.  It should probably go somewhere
 * else, but...  Note that we pretend that dependent sources and mutual
 * inductors have more nodes than they really do...
 *----------------------------------------------------------------------*/
static int
inp_numnodes(char c)
{
    if (isupper(c))
        c = (char) tolower(c);
    switch (c) {
    case ' ':
    case '\t':
    case '.':
    case 'x':
    case '*':
    case '$':
        return (0);

    case 'b':
        return (2);
    case 'c':
        return (2);
    case 'd':
        return (2);
    case 'e':
        return (2); /* changed from 4 to 2 by SDB on 4.22.2003 to enable POLY */
    case 'f':
        return (2);
    case 'g':
        return (2); /* changed from 4 to 2 by SDB on 4.22.2003 to enable POLY */
    case 'h':
        return (2);
    case 'i':
        return (2);
    case 'j':
        return (3);
    case 'k':
        return (0);
    case 'l':
        return (2);
    case 'm':
        return (7); /* This means that 7 is the maximun number of nodes */
    case 'o':
        return (4);
    case 'p':
        return (18);/* 16 lines + 2 gnd is the maximum number of nodes for CPL */
    case 'q':
        return (5);
    case 'r':
        return (2);
    case 's':
        return (4);
    case 't':
        return (4);
    case 'u':
        return (3);
    case 'v':
        return (2);
    case 'w':
        return (2); /* change 3 to 2 here to fix w bug, NCF 1/31/95 */
    case 'y':
        return (4);
    case 'z':
        return (3);

    default:
        fprintf(cp_err, "Warning: unknown device type: %c\n", c);
        return (2);
    }
}
