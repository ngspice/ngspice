/*       spicenum.c                Copyright (C)  2002    Georg Post
 *
 *  This file is part of Numparam, see:  readme.txt
 *  Free software under the terms of the GNU Lesser General Public License
 */

/* number parameter add-on for Spice.
   to link with mystring.o, xpressn.o (math formula interpreter),
   and with Spice frontend src/lib/fte.a .
   Interface function nupa_signal to tell us about automaton states.
Buglist (some are 'features'):
  blank lines get category '*'
  inserts conditional blanks before or after  braces
  between .control and .endc, flags all lines as 'category C', dont touch.
  there are reserved magic numbers (1e9 + n) as placeholders
  control lines must not contain {} .
  ignores the '.option numparam' line planned to trigger the actions
  operation of .include certainly doesnt work
  there are frozen maxima for source and expanded circuit size.
Todo:
  add support for nested .if .elsif .else .endif controls.
*/

#include "ngspice/ngspice.h"

#include "general.h"
#include "numparam.h"

#include "ngspice/fteext.h"


extern bool ft_batchmode;

void dump_symbols(tdico *dico_p);

char *nupa_inst_name;

/* number of parameter substitutions, available only after the substitution */
extern long dynsubst; /* spicenum.c:144 */

/* number of lines in input deck */
extern int dynmaxline; /* inpcom.c:1529 */

/* Uncomment this line to allow debug tracing */
/* #define TRACE_NUMPARAMS */

/*  the nupa_signal arguments sent from Spice:

    sig=1: Start of the subckt expansion.
    sig=2: Stop of the subckt expansion.
    sig=3: Stop of the evaluation phase.
    sig=0: Start of a deck copy operation

    After sig=1 until sig=2, nupa_copy does no transformations.
    At sig=2, we prepare for nupa_eval loop.
    After sig=3, we assume the initial state (clean).

    In Clean state, a lot of deckcopy operations come in and we
    overwrite any line pointers, or we start a new set after each sig=0 ?
    Anyway, we neutralize all & and .param lines  (category[] array!)
    and we substitute all {} &() and &id placeholders by dummy identifiers.
    those look like numparm__________XXXXXXXX (8 hexadecimal digits)

*/
/**********  string handling ***********/

static long placeholder = 0;


static void
stripsomespace(SPICE_DSTRINGPTR dstr_p, unsigned char incontrol)
{
    /* if s starts with one of some markers, strip leading space */
    int i, ls;
    char *sstr;                 /* string contained in s */
    SPICE_DSTRING markers;

    spice_dstring_init(&markers);
    scopys(&markers, "*.&+#$");

    if (!incontrol)
        sadd(&markers, "xX");

    sstr = spice_dstring_value(dstr_p);
    ls = spice_dstring_length(dstr_p);

    i = 0;
    while ((i < ls) && (sstr[i] <= ' '))
        i++;

    if ((i > 0) && (i < ls) && (cpos(sstr[i], spice_dstring_value(&markers)) >= 0))
        pscopy(dstr_p, sstr, i, ls);
}


static int
stripbraces(SPICE_DSTRINGPTR dstr_p)
/* puts the funny placeholders. returns the number of {...} substitutions */
{
    int n, i, nest, ls, j;
    char *s;                    /* value of dynamic string */
    char *t_p;                  /* value of t dynamic string */
    SPICE_DSTRING tstr;         /* temporary dynamic string */

    n = 0;
    spice_dstring_init(&tstr);
    s = spice_dstring_value(dstr_p);
    ls = spice_dstring_length(dstr_p);
    i = 0;

    while (i < ls) {

        if (s[i] == '{') {

            /* something to strip */
            j = i + 1;
            nest = 1;
            n++;

            while ((nest > 0) && (j < ls)) {
                if (s[j] == '{')
                    nest++;
                else if (s[j] == '}')
                    nest--;
                j++;
            }

            pscopy(&tstr, s, 0, i);
            placeholder++;

            t_p = spice_dstring_value(&tstr);

            if (t_p[i - 1] > ' ')
                cadd(&tstr, ' ');

            cadd(&tstr, ' ');
            {
                char buf[25+1];
                sprintf(buf, "numparm__________%08lx", placeholder);
                sadd(&tstr, buf);
            }
            cadd(&tstr, ' ');

            if (s[j] >= ' ')
                cadd(&tstr, ' ');

            i = spice_dstring_length(&tstr);
            pscopy(dstr_p, s, j, ls);
            sadd(&tstr, s);
            scopyd(dstr_p, &tstr);
            s = spice_dstring_value(dstr_p);
            ls = spice_dstring_length(dstr_p);

        } else {

            i++;

        }
    }

    dynsubst = placeholder;
    spice_dstring_free(&tstr);

    return n;
}


static int
findsubname(tdico *dico, SPICE_DSTRINGPTR dstr_p)
/* truncate the parameterized subckt call to regular old Spice */
/* scan a string from the end, skipping non-idents and {expressions} */
/* then truncate s after the last subckt(?) identifier */
{
    SPICE_DSTRING name;         /* extract a name */
    char *s;                    /* current dstring */
    int h, j, k, nest, ls;
    int found;

    h = 0;

    ls = spice_dstring_length(dstr_p);
    s = spice_dstring_value(dstr_p);
    k = ls - 1;                 /* now a C - string */
    found = 0;
    spice_dstring_init(&name);

    while ((k >= 0) && (!found)) {

        /* skip space, then non-space */
        while ((k >= 0) && (s[k] <= ' '))
            k--;

        h = k + 1;              /* at h: space */
        while ((k >= 0) && (s[k] > ' ')) {

            if (s[k] == '}') {
                nest = 1;
                k--;

                while ((nest > 0) && (k >= 0)) {
                    if (s[k] == '{')
                        nest--;
                    else if (s[k] == '}')
                        nest++;

                    k--;
                }
                h = k + 1;      /* h points to '{' */

            } else {
                k--;
            }
        }

        found = (k >= 0) && alfanum(s[k + 1]); /* suppose an identifier */
        if (found) {
            /* check for known subckt name */
            spice_dstring_reinit(&name);
            j = k + 1;
            while (alfanum(s[j])) {
                cadd(&name, upcase(s[j]));
                j++;
            }
            found = (getidtype(dico, spice_dstring_value(&name)) == 'U');
        }
    }

    if (found && (h < ls))
        pscopy(dstr_p, s, 0, h);

    return h;
}


static void
modernizeex(SPICE_DSTRINGPTR dstr_p)
/* old style expressions &(..) and &id --> new style with braces. */
{
    int i, state, ls;
    char c, d;
    char *s;                    /* current string */
    SPICE_DSTRING t;            /* temporary dyna string */

    i = 0;
    state = 0;
    ls = spice_dstring_length(dstr_p);
    s = spice_dstring_value(dstr_p);

    /* check if string might need modernizing */
    if (!memchr(s, Intro, (size_t) ls))
        return;

    spice_dstring_init(&t);

    while (i < ls) {
        c = s[i];
        d = s[i + 1];
        if ((!state) && (c == Intro) && (i > 0)) {
            if (d == '(') {
                state = 1;
                i++;
                c = '{';
            } else if (alfa(d)) {
                cadd(&t, '{');
                i++;
                while (alfanum(s[i])) {
                    cadd(&t, s[i]);
                    i++;
                }
                c = '}';
                i--;
            }
        } else if (state) {
            if (c == '(')
                state++;
            else if (c == ')')
                state--;

            if (!state)         /* replace--) by terminator */
                c = '}';
        }

        cadd(&t, c);
        i++;
    }

    scopyd(dstr_p, &t);
    spice_dstring_free(&t);
}


static char
transform(tdico *dico, SPICE_DSTRINGPTR dstr_p, unsigned char nostripping,
          SPICE_DSTRINGPTR u_p)
/*         line s is categorized and crippled down to basic Spice
 *         returns in u control word following dot, if any
 *
 * any + line is copied as-is.
 * any & or .param line is commented-out.
 * any .subckt line has params section stripped off
 * any X line loses its arguments after sub-circuit name
 * any &id or &() or {} inside line gets a 10-digit substitute.
 *
 * strip  the new syntax off the codeline s, and
 * return the line category as follows:
 *   '*'  comment line
 *   '+'  continuation line
 *   ' '  other untouched netlist or command line
 *   'P'  parameter line, commented-out; (name,linenr)-> symbol table.
 *   'S'  subckt entry line, stripped;   (name,linenr)-> symbol table.
 *   'U'  subckt exit line
 *   'X'  subckt call line, stripped
 *   'C'  control entry line
 *   'E'  control exit line
 *   '.'  any other dot line
 *   'B'  netlist (or .model ?) line that had Braces killed
 */
{
    int k, a, n;
    char *s;                    /* dstring value of dstr_p */
    char *t;                    /* dstring value of tstr */
    char category;
    SPICE_DSTRING tstr;         /* temporary string */

    spice_dstring_init(&tstr);
    spice_dstring_reinit(u_p);
    stripsomespace(dstr_p, nostripping);
    modernizeex(dstr_p);        /* required for stripbraces count */

    s = spice_dstring_value(dstr_p);

    if (s[0] == '.') {
        /* check PS parameter format */
        scopy_up(&tstr, spice_dstring_value(dstr_p));
        k = 1;

        t = spice_dstring_value(&tstr);
        while (t[k] > ' ') {
            cadd(u_p, t[k]);
            k++;
        }

        if (ci_prefix(".PARAM", t) == 1) {
            /* comment it out */
            /* s[0] = '*'; */
            category = 'P';
        } else if (ci_prefix(".SUBCKT", t) == 1) {
            /* split off any "params" tail */
            a = spos_("PARAMS:", t);
            if (a >= 0)
                pscopy(dstr_p, s, 0, a);
            category = 'S';
        } else if (ci_prefix(".CONTROL", t) == 1) {
            category = 'C';
        } else if (ci_prefix(".ENDC", t) == 1) {
            category = 'E';
        } else if (ci_prefix(".ENDS", t) == 1) {
            category = 'U';
        } else {
            category = '.';
            n = stripbraces(dstr_p);
            if (n > 0)
                category = 'B'; /* priority category ! */
        }
    } else if (s[0] == Intro) {
        /* private style preprocessor line */
        s[0] = '*';
        category = 'P';
    } else if (upcase(s[0]) == 'X') {
        /* strip actual parameters */
        findsubname(dico, dstr_p); /* i= index following last identifier in s */
        category = 'X';
    } else if (s[0] == '+') {   /* continuation line */
        category = '+';
    } else if (cpos(s[0], "*$#") < 0) {
        /* not a comment line! */
        n = stripbraces(dstr_p);
        if (n > 0)
            category = 'B';     /* line that uses braces */
        else
            category = ' ';     /* ordinary code line */
    } else {
        category = '*';
    }

    spice_dstring_free(&tstr);
    return category;
}


/************ core of numparam **************/

/* some day, all these nasty globals will go into the tdico structure
   and everything will get hidden behind some "handle" ...
   For the time being we will rename this variable to end in S so we know
   they are statics within this file for easier reading of the code.
*/

static int linecountS = 0;      /* global: number of lines received via nupa_copy */
static int evalcountS = 0;      /* number of lines through nupa_eval() */
static int nblogS = 0;          /* serial number of (debug) logfile */
static unsigned char inexpansionS = 0;  /* flag subckt expansion phase */
static unsigned char incontrolS = 0;    /* flag control code sections */
static unsigned char dologfileS = 0;    /* for debugging */
static unsigned char firstsignalS = 1;
static FILE *logfileS = NULL;
static tdico *dicoS = NULL;


/*  already part of dico : */
/*
  Open ouput to a log file.
  takes no action if logging is disabled.
  Open the log if not already open.
*/
static void
putlogfile(char c, int num, char *t)
{
    if (!dologfileS)
        return;

    if (!logfileS) {
        char *fname = tprintf("logfile.%d", ++nblogS);
        logfileS = fopen(fname, "w");
        tfree(fname);
    }

    if (logfileS)
        fprintf(logfileS, "%c%d: %s\n", c, num, t);
}


static void
nupa_init(char *srcfile)
{
    int i;

    /* init the symbol table and so on, before the first  nupa_copy. */
    evalcountS = 0;
    linecountS = 0;
    incontrolS = 0;
    placeholder = 0;
    dicoS = (tdico *) new(sizeof(tdico));
    initdico(dicoS);

    dicoS->dynrefptr = TMALLOC(char*, dynmaxline + 1);
    dicoS->dyncategory = TMALLOC(char, dynmaxline + 1);

    for (i = 0; i <= dynmaxline; i++) {
        dicoS->dynrefptr[i] = NULL;
        dicoS->dyncategory[i] = '?';
    }

    if (srcfile != NULL)
        scopys(&dicoS->srcfile, srcfile);
}


/* free dicoS (called from com_remcirc()) */
void
nupa_del_dicoS(void)
{
    int i;

    if(!dicoS)
        return;

    for (i = dynmaxline; i >= 0; i--)
        dispose(dicoS->dynrefptr[i]);

    dispose(dicoS->dynrefptr);
    dispose(dicoS->dyncategory);
    dispose(dicoS->inst_name);
    dispose(dicoS->local_symbols);
    nghash_free(dicoS->global_symbols, del_attrib, NULL);
    dispose(dicoS);
    dicoS = NULL;
}


static void
nupa_done(void)
{
    /* int i; not needed so far, see below */
    SPICE_DSTRING rep;          /* dynamic report */
    int dictsize, nerrors;

    spice_dstring_init(&rep);

    if (logfileS != NULL) {
        fclose(logfileS);
        logfileS = NULL;
    }

    nerrors = dicoS->errcount;
    dictsize = donedico(dicoS);

    /* We cannot remove dicoS here because numparam is used by
       the .measure statements, which are invoked only after the
       simulation has finished. */

    if (nerrors) {
        /* debug: ask if spice run really wanted */
        sadd(&rep, " Copies=");
        nadd(&rep, linecountS);
        sadd(&rep, " Evals=");
        nadd(&rep, evalcountS);
        sadd(&rep, " Placeholders=");
        nadd(&rep, placeholder);
        sadd(&rep, " Symbols=");
        nadd(&rep, dictsize);
        sadd(&rep, " Errors=");
        nadd(&rep, nerrors);
        cadd(&rep, '\n');
        printf("%s", spice_dstring_value(&rep));
        if (ft_batchmode)
            controlled_exit(EXIT_FAILURE);
        for (;;) {
            int c;
            printf("Numparam expansion errors: Run Spice anyway? y/n ?\n");
            c = yes_or_no();
            if (c == 'n' || c == EOF)
                controlled_exit(EXIT_FAILURE);
            if (c == 'y')
                break;
        }
    }

    linecountS = 0;
    evalcountS = 0;
    placeholder = 0;
    /* release symbol table data */
}


/* SJB - Scan the line for subcircuits */
void
nupa_scan(char *s, int linenum, int is_subckt)
{
    if (is_subckt)
        defsubckt(dicoS, s, linenum, 'U');
    else
        defsubckt(dicoS, s, linenum, 'O');
}


/* -----------------------------------------------------------------
 * Dump the contents of a symbol table.
 * ----------------------------------------------------------------- */
static void
dump_symbol_table(tdico *dico_p, NGHASHPTR htable_p, FILE *cp_out)
{
    char *name;                 /* current symbol */
    entry *entry_p;             /* current entry */
    NGHASHITER iter;            /* hash iterator - thread safe */

    NGHASH_FIRST(&iter);
    for (entry_p = (entry *) nghash_enumerateRE(htable_p, &iter);
         entry_p;
         entry_p = (entry *) nghash_enumerateRE(htable_p, &iter))
    {
        if (entry_p->tp == 'R') {
            spice_dstring_reinit(& dico_p->lookup_buf);
            scopy_lower(& dico_p->lookup_buf, entry_p->symbol);
            name = spice_dstring_value(& dico_p->lookup_buf);
            fprintf(cp_out, "       ---> %s = %g\n", name, entry_p->vl);
            spice_dstring_free(& dico_p->lookup_buf);
        }
    }
}


/* -----------------------------------------------------------------
 * Dump the contents of the symbol table.
 * ----------------------------------------------------------------- */
void
nupa_list_params(FILE *cp_out)
{
    int depth;                  /* nested subcircit depth */
    tdico *dico_p;              /* local copy for speed */
    NGHASHPTR htable_p;         /* current hash table */

    dico_p = dicoS;
    if (dico_p == NULL) {
        fprintf(cp_err, "\nWarning: No symbol table available for 'listing param'\n");
        return;
    }

    fprintf(cp_out, "\n\n");

    /* -----------------------------------------------------------------
     * Print out the locally defined symbols from highest to lowest priority.
     * If there are no parameters, the hash table will not be allocated as
     * we use lazy allocation to save memory.
     * ----------------------------------------------------------------- */
    for (depth = dico_p->stack_depth; depth > 0; depth--) {
        htable_p = dico_p->local_symbols[depth];
        if (htable_p) {
            fprintf(cp_out, " local symbol definitions for:%s\n", dico_p->inst_name[depth]);
            dump_symbol_table(dico_p, htable_p, cp_out);
        }
    }

    /* -----------------------------------------------------------------
     * Finally dump the global symbols.
     * ----------------------------------------------------------------- */
    fprintf(cp_out, " global symbol definitions:\n");
    dump_symbol_table(dico_p, dico_p->global_symbols, cp_out);
}


/* -----------------------------------------------------------------
 * Lookup a parameter value in the symbol tables.   This involves
 * multiple lookups in various hash tables in order to get the scope
 * correct.  Each subcircuit instance will have its own local hash
 * table if it has parameters.   We can return whenever we get a hit.
 * Otherwise, we have to exhaust all of the tables including the global
 * table.
 * ----------------------------------------------------------------- */
double
nupa_get_param(char *param_name, int *found)
{
    int depth;                  /* nested subcircit depth */
    char *up_name;              /* current parameter upper case */
    entry *entry_p;             /* current entry */
    tdico *dico_p;              /* local copy for speed */
    NGHASHPTR htable_p;         /* current hash table */
    double result = 0;          /* parameter value */

    dico_p = dicoS;
    spice_dstring_reinit(& dico_p->lookup_buf);
    scopy_up(& dico_p->lookup_buf, param_name);
    up_name = spice_dstring_value(& dico_p->lookup_buf);

    *found = 0;
    for (depth = dico_p->stack_depth; depth > 0; depth--) {
        htable_p = dico_p->local_symbols[depth];
        if (htable_p) {
            entry_p = (entry *) nghash_find(htable_p, up_name);
            if (entry_p) {
                result = entry_p->vl;
                *found = 1;
                break;
            }
        }
    }

    if (!(*found)) {
        /* No luck.  Try the global table. */
        entry_p = (entry *) nghash_find(dico_p->global_symbols, up_name);
        if (entry_p) {
            result = entry_p->vl;
            *found = 1;
        }
    }

    spice_dstring_free(& dico_p->lookup_buf);
    return result;
}


void
nupa_add_param(char *param_name, double value)
{
    char *up_name;              /* current parameter upper case */
    entry *entry_p;             /* current entry */
    tdico *dico_p;              /* local copy for speed */
    NGHASHPTR htable_p;         /* hash table of interest */

    dico_p = dicoS;
    /* -----------------------------------------------------------------
     * We use a dynamic string here because most of the time we will
     * be using short names and no memory allocation will occur.
     * ----------------------------------------------------------------- */
    spice_dstring_reinit(& dico_p->lookup_buf);
    scopy_up(& dico_p->lookup_buf, param_name);
    up_name = spice_dstring_value(& dico_p->lookup_buf);

    if (dico_p->stack_depth > 0) {
        /* can't be lazy anymore */
        if (!(dico_p->local_symbols[dico_p->stack_depth]))
            dico_p->local_symbols[dico_p->stack_depth] = nghash_init(NGHASH_MIN_SIZE);
        htable_p = dico_p->local_symbols[dico_p->stack_depth];
    } else {
        /* global symbol */
        htable_p = dico_p->global_symbols;
    }

    entry_p = attrib(dico_p, htable_p, up_name, 'N');
    if (entry_p) {
        entry_p->vl = value;
        entry_p->tp = 'R';
        entry_p->ivl = 0;
        entry_p->sbbase = NULL;
    }

    spice_dstring_free(& dico_p->lookup_buf);
}


void
nupa_add_inst_param(char *param_name, double value)
{
    char *up_name;              /* current parameter upper case */
    entry *entry_p;             /* current entry */
    tdico *dico_p;              /* local copy for speed */

    dico_p = dicoS;
    spice_dstring_reinit(& dico_p->lookup_buf);
    scopy_up(& dico_p->lookup_buf, param_name);
    up_name = spice_dstring_value(& dico_p->lookup_buf);

    if (!(dico_p->inst_symbols))
        dico_p->inst_symbols = nghash_init(NGHASH_MIN_SIZE);

    entry_p = attrib(dico_p, dico_p->inst_symbols, up_name, 'N');
    if (entry_p) {
        entry_p->vl = value;
        entry_p->tp = 'R';
        entry_p->ivl = 0;
        entry_p->sbbase = NULL;
    }

    spice_dstring_free(& dico_p->lookup_buf);
}


/* -----------------------------------------------------------------
 * This function copies any definitions in the inst_symbols hash
 * table which are qualified symbols and makes them available at
 * the global level.  Afterwards, the inst_symbols table is freed.
 * ----------------------------------------------------------------- */
void
nupa_copy_inst_dico(void)
{
    entry *entry_p;             /* current entry */
    tdico *dico_p;              /* local copy for speed */
    NGHASHITER iter;            /* hash iterator - thread safe */

    dico_p = dicoS;
    if (dico_p->inst_symbols) {
        /* We we perform this operation we should be in global scope */
        if (dico_p->stack_depth > 0)
            fprintf(stderr, "stack depth should be zero.\n");

        NGHASH_FIRST(&iter);
        for (entry_p = (entry *) nghash_enumerateRE(dico_p->inst_symbols, &iter);
             entry_p;
             entry_p = (entry *) nghash_enumerateRE(dico_p->inst_symbols, &iter))
        {
            nupa_add_param(entry_p->symbol, entry_p->vl);
            dico_free_entry(entry_p);
        }

        nghash_free(dico_p->inst_symbols, NULL, NULL);
        dico_p->inst_symbols = NULL;
    }
}


char *
nupa_copy(char *s, int linenum)
/* returns a copy (not quite) of s in freshly allocated memory.
   linenum, for info only, is the source line number.
   origin pointer s is kept, memory is freed later in nupa_done.
   must abort all Spice if malloc() fails.
   :{ called for the first time sequentially for all spice deck lines.
   :{ then called again for all X invocation lines, top-down for
   subckts defined at the outer level, but bottom-up for local
   subcircuit expansion, but has no effect in that phase.
   we steal a copy of the source line pointer.
   - comment-out a .param or & line
   - substitute placeholders for all {..} --> 10-digit numeric values.
*/
{
    char *t;
    int ls;
    char c, d;
    SPICE_DSTRING u;
    SPICE_DSTRING keywd;

    spice_dstring_init(&u);
    spice_dstring_init(&keywd);
    ls = length(s);

    while ((ls > 0) && (s[ls - 1] <= ' '))
        ls--;

    pscopy(&u, s, 0, ls);       /* strip trailing space, CrLf and so on */
    dicoS->srcline = linenum;

    if ((!inexpansionS) && (linenum >= 0) && (linenum <= dynmaxline)) {
        linecountS++;
        dicoS->dynrefptr[linenum] = s;
        c = transform(dicoS, &u, incontrolS, &keywd);
        if (c == 'C')
            incontrolS = 1;
        else if (c == 'E')
            incontrolS = 0;

        if (incontrolS)
            c = 'C';            /* force it */

        d = dicoS->dyncategory[linenum]; /* warning if already some strategic line! */

        if ((d == 'P') || (d == 'S') || (d == 'X'))
            fprintf(stderr,
                    " Numparam warning: overwriting P,S or X line (linenum == %d).\n",
                    linenum);
        dicoS->dyncategory[linenum] = c;
    } /* keep a local copy and mangle the string */

    ls = spice_dstring_length(&u);
    t = strdup(spice_dstring_value(&u));

    if (t == NULL) {
        fputs("Fatal: String malloc crash in nupa_copy()\n", stderr);
        controlled_exit(EXIT_FAILURE);
    } else {
        if (!inexpansionS)
            putlogfile(dicoS->dyncategory[linenum], linenum, t);
    }

    spice_dstring_free(&u);
    return t;
}


int
nupa_eval(char *s, int linenum, int orig_linenum)
/* s points to a partially transformed line.
   compute variables if linenum points to a & or .param line.
   if ( the original is an X line,  compute actual params.;
   } else {  substitute any &(expr) with the current values.
   All the X lines are preserved (commented out) in the expanded circuit.
*/
{
    int idef;                   /* subckt definition line */
    char c, keep, *ptr;
    SPICE_DSTRING subname;      /* dynamic string for subcircuit name */
    bool err = 1;

    spice_dstring_init(&subname);
    dicoS->srcline = linenum;
    dicoS->oldline = orig_linenum;

    c = dicoS->dyncategory[linenum];

#ifdef TRACE_NUMPARAMS
    fprintf(stderr, "** SJB - in nupa_eval()\n");
    fprintf(stderr, "** SJB - processing line %3d: %s\n", linenum, s);
    fprintf(stderr, "** SJB - category '%c'\n", c);
#endif

    if (c == 'P') {                     /* evaluate parameters */
        // err = nupa_substitute(dico, dico->dynrefptr[linenum], s, 0);
        nupa_assignment(dicoS, dicoS->dynrefptr[linenum], 'N');
    } else if (c == 'B') {              /* substitute braces line */
        err = nupa_substitute(dicoS, dicoS->dynrefptr[linenum], s, 0);
    } else if (c == 'X') {
        /* compute args of subcircuit, if required */
        ptr = s;
        while (!isspace(*ptr))
            ptr++;
        keep = *ptr;
        *ptr = '\0';
        nupa_inst_name = strdup(s);
        *nupa_inst_name = 'x';
        *ptr = keep;

        strtoupper(nupa_inst_name);

        idef = findsubckt(dicoS, s, &subname);
        if (idef > 0)
            nupa_subcktcall(dicoS, dicoS->dynrefptr[idef], dicoS->dynrefptr[linenum], 0);
        else
            putlogfile('?', linenum, "  illegal subckt call.");
    } else if (c == 'U') {              /*  release local symbols = parameters */
        nupa_subcktexit(dicoS);
    }

    putlogfile('e', linenum, s);
    evalcountS++;

#ifdef TRACE_NUMPARAMS
    fprintf(stderr, "** SJB - leaving nupa_eval(): %s   %d\n", s, err);
    printf("** SJB -                  --> %s\n", s);
    printf("** SJB - leaving nupa_eval()\n\n");
#endif

    if (err)
        return 0;
    else
        return 1;
}


int
nupa_signal(int sig, char *info)
/* warning: deckcopy may come inside a recursion ! substart no! */
/* info is context-dependent string data */
{
    putlogfile('!', sig, " Nupa Signal");

    if (sig == NUPADECKCOPY) {
        if (firstsignalS) {
            nupa_init(info);
            firstsignalS = 0;
        }
    } else if (sig == NUPASUBSTART) {
        inexpansionS = 1;
    } else if (sig == NUPASUBDONE) {
        inexpansionS = 0;
        nupa_inst_name = NULL;
    } else if (sig == NUPAEVALDONE) {
        nupa_done();
        firstsignalS = 1;
    }

    return 1;
}


#ifdef USING_NUPATEST
/* This is use only by the nupatest program */
tdico *
nupa_fetchinstance(void)
{
    return dico;
}
#endif


void dump_symbols(tdico *dico_p)
{
    NG_IGNORE(dico_p);

    fprintf(stderr, "Symbol table\n");
    nupa_list_params(stderr);
}
