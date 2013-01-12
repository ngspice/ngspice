/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Wayne A. Christopher
**********/

/*
   For dealing with spice input decks and command scripts

   Central function is inp_readall()
 */

#include "ngspice/ngspice.h"

#include "ngspice/cpdefs.h"
#include "ngspice/ftedefs.h"
#include "ngspice/fteext.h"
#include "ngspice/dvec.h"
#include "ngspice/fteinp.h"
#include "ngspice/compatmode.h"

#include "inpcom.h"
#include "variable.h"
#include "../misc/util.h" /* ngdirname() */
#include "ngspice/stringutil.h"
#include "ngspice/wordlist.h"

#ifdef XSPICE
/* gtri - add - 12/12/90 - wbk - include new stuff */
#include "ngspice/ipctiein.h"
#include "ngspice/enh.h"
/* gtri - end - 12/12/90 */
#endif

/* SJB - Uncomment this line for debug tracing */
/*#define TRACE*/

/* globals -- wanted to avoid complicating inp_readall interface */
static char *library_file[1000];
static char *library_name[1000][1000];
struct line *library_ll_ptr[1000][1000];
struct line *libraries[1000];
int         num_libraries;
int         num_lib_names[1000];
static      char *global;
static char *subckt_w_params[1000];
static int  num_subckt_w_params;
static char *func_names[1000];
static char *func_params[1000][1000];
static char *func_macro[5000];
static int  num_functions;
static int  num_parameters[1000];
static COMPATMODE_T inp_compat_mode;

/* Collect information for dynamic allocation of numparam arrays */
/* number of lines in input deck */
int dynmaxline;  /* inpcom.c 1529 */
/* number of lines in deck after expansion */
int dynMaxckt = 0; /* subckt.c 307 */
/* number of parameter substitutions */
long dynsubst; /* spicenum.c 221 */


static char *readline(FILE *fd);
static int  get_number_terminals(char *c);
static void inp_stripcomments_deck(struct line *deck);
static void inp_stripcomments_line(char *s);
static void inp_fix_for_numparam(struct line *deck);
static void inp_remove_excess_ws(struct line *deck);
static void inp_determine_libraries(struct line *deck, char *lib_name);
static void inp_init_lib_data(void);
static void inp_grab_func(struct line *deck);
static void inp_fix_inst_calls_for_numparam(struct line *deck);
static void inp_expand_macros_in_func(void);
static void inp_expand_macros_in_deck(struct line *deck);
static void inp_fix_param_values(struct line *deck);
static void inp_reorder_params(struct line *deck, struct line *list_head, struct line *end);
static int  inp_split_multi_param_lines(struct line *deck, int line_number);
static void inp_sort_params(struct line *start_card, struct line *end_card, struct line *card_bf_start, struct line *s_c, struct line *e_c);
static char *inp_remove_ws(char *s);
static void inp_compat(struct line *deck);
static void inp_bsource_compat(struct line *deck);

static bool chk_for_line_continuation(char *line);
static void comment_out_unused_subckt_models(struct line *start_card, int no_of_lines);
static void inp_fix_macro_param_func_paren_io(struct line *begin_card);
static char *inp_fix_ternary_operator_str(char *line, bool all);
static void inp_fix_ternary_operator(struct line *start_card);
static void inp_fix_gnd_name(struct line *deck);
static void inp_chk_for_multi_in_vcvs(struct line *deck, int *line_number);
static void inp_add_control_section(struct line *deck, int *line_number);
static char *get_quoted_token(char *string, char **token);
static void replace_token(char *string, char *token, int where, int total);
static void inp_add_series_resistor(struct line *deck);

static char *skip_back_non_ws(char *d) { while (*d && !isspace(*d)) d--; return d; }
static char *skip_non_ws(char *d)      { while (*d && !isspace(*d)) d++; return d; }
static char *skip_back_ws(char *d)     { while (isspace(*d))        d--; return d; }
static char *skip_ws(char *d)          { while (isspace(*d))        d++; return d; }

#ifndef XSPICE
static void inp_poly_err(struct line *deck);
#endif


static int
find_lib_name(int i, char *s) {
    int j;
    for (j = 0; j < num_lib_names[i]; j++)
        if (strcmp(library_name[i][j], s) == 0)
            break;
    return -1;
}


static void
new_lib_name(int i, char *y, struct line *c) {
    library_ll_ptr[i][num_lib_names[i]] = c;
    library_name[i][num_lib_names[i]++] = strdup(y);
}


static int
expand_libs(int line_number)
{
    struct line *tmp_ptr = NULL, *prev;
    bool found_lib_name = FALSE;
    int i;
    for (i = 0; i < num_libraries; i++) {
        struct line *working = libraries[i];
        while (working) {
            char *buffer = working->li_line;

            if (found_lib_name && ciprefix(".endl", buffer)) {
                struct line *tmp_ptr2;
                /* Make the .endl a comment */
                *buffer = '*';
                found_lib_name = FALSE;

                /* set pointer and continue to avoid deleting below */
                tmp_ptr2         = working->li_next;
                working->li_next = tmp_ptr;
                working          = tmp_ptr2;

                /* end          = working;
                 * working      = working->li_next;
                 * end->li_next = NULL; */

                continue;
            }

            if (ciprefix(".lib", buffer)) {

                char keep_char;
                int j;
                char *s, *t;

                if (found_lib_name == TRUE) {
                    fprintf(stderr, "ERROR: .lib is missing .endl!\n");
                    controlled_exit(EXIT_FAILURE);
                }

                s = skip_non_ws(buffer);       /* skip over .lib                        */
                while (isspace(*s) || isquote(*s))
                    s++;    /* advance past space chars              */
                for (t = s; *t && !isspace(*t) && !isquote(*t); t++)
                    ;       /* skip to end of word                   */
                keep_char = *t;
                *t = '\0';
                /* see if library we want to copy */

                j = find_lib_name(i, s);

                found_lib_name = (j >= 0);

                if (found_lib_name) {
                    struct line *start_lib = working;
                    int line_number_lib;

                    /* make the .lib a comment */
                    *buffer = '*';

                    tmp_ptr = library_ll_ptr[i][j]->li_next;
                    library_ll_ptr[i][j]->li_next = working;

                    /* renumber lines */
                    line_number_lib = 1;
                    for (start_lib = working; !ciprefix(".endl", start_lib->li_line); start_lib = start_lib->li_next) {
                        start_lib->li_linenum = line_number++;
                        start_lib->li_linenum_orig = line_number_lib++;
                    }
                    start_lib->li_linenum = line_number++;  // renumber endl line
                    start_lib->li_linenum_orig = line_number_lib++;
                }
                *t = keep_char;
            }
            prev = working;
            working = working->li_next;

            if (found_lib_name == FALSE) {
                tfree(prev->li_line);
                tfree(prev);
            }
        } /* end while */
    } /* end for */

    return line_number;
}


/*-------------------------------------------------------------------------
 Read the entire input file and return  a pointer to the first line of
 the linked list of 'card' records in data.  The pointer is stored in
 *data.
 Called from fcn inp_spsource() in inp.c to load circuit or command files.
 Called from fcn com_alter_mod() in device.c to load model files.
 Called from here to load .library or .include files.

 Procedure:
 read in all lines & put them in the struct cc
   read next line
   process .TITLE line
     store contents in string new_title
   process .lib lines
     read file and library name, open file using fcn inp_pathopen()
     read file contents and put into struct libraries[], one entry per .lib line
   process .inc lines
     read file and library name, open file using fcn inp_pathopen()
     read file contents and add lines to cc
   make line entry lower case
   allow for shell end of line continuation (\\)
     add '+' to beginning of next line
   add line entry to list cc
 add '.global gnd'
 add libraries
   find library section
   add lines
 add .end card
 strip end-of-line comments
 make continuation lines a single line
*** end of processing for command files ***
 start preparation of input deck for numparam
   ...
 debug printout to debug-out.txt
 *-------------------------------------------------------------------------*/
void
inp_readall(FILE *fp, struct line **data, int call_depth, char *dir_name, bool comfile)
/* fp: in, pointer to file to be read,
   data: out, linked list of cards
   call_depth: in, nested call to fcn
   dir_name: in, name of directory of file to be read
   comfile: in, TRUE if command file (e.g. spinit, .spiceinit
*/
{
    struct line *end = NULL, *cc = NULL, *prev, *working, *newcard, *global_card;
    char *buffer = NULL, c;
    /* segfault fix */
#ifdef XSPICE
    char big_buff[5000];
    int line_count = 0;
    Ipc_Status_t    ipc_status;
    char            ipc_buffer[1025];  /* Had better be big enough */
    int             ipc_len;
#endif
    char *new_title = NULL;
    int line_number = 1; /* sjb - renamed to avoid confusion with struct line */
    int line_number_orig = 1, line_number_inc = 1;
    unsigned int no_braces = 0; /* number of '{' */

    size_t max_line_length; /* max. line length in input deck */

    FILE *fdo;
    struct line *tmp_ptr1 = NULL;

    bool found_end = FALSE, shell_eol_continuation = FALSE;

    if (call_depth == 0) {
        num_subckt_w_params = 0;
        num_libraries       = 0;
        num_functions       = 0;
        global              = NULL;
        found_end           = FALSE;
        inp_compat_mode = ngspice_compat_mode();
    }

    /* First read in all lines & put them in the struct cc */
    for (;;) {

    /*   gtri - modify - 12/12/90 - wbk - read from mailbox if ipc enabled   */
#ifdef XSPICE

        /* If IPC is not enabled, do equivalent of what SPICE did before */
        if (! g_ipc.enabled) {
            if (call_depth == 0 && line_count == 0) {
                line_count++;
                if (fgets(big_buff, 5000, fp))
                    buffer = copy(big_buff);
            } else {
                buffer = readline(fp);
                if (!buffer)
                    break;
            }
        } else {
            /* else, get the line from the ipc channel. */
            /* We assume that newlines are not sent by the client */
            /* so we add them here */
            ipc_status = ipc_get_line(ipc_buffer, &ipc_len, IPC_WAIT);
            if (ipc_status == IPC_STATUS_END_OF_DECK) {
                buffer = NULL;
                break;
            } else if (ipc_status == IPC_STATUS_OK) {
                buffer = TMALLOC(char, strlen(ipc_buffer) + 3);
                strcpy(buffer, ipc_buffer);
                strcat(buffer, "\n");
            } else {            /* No good way to report this so just die */
                controlled_exit(EXIT_FAILURE);
            }
        }

        /* gtri - end - 12/12/90 */
#else
        buffer = readline(fp);
        if(!buffer)
            break;
#endif

#ifdef TRACE
        /* SDB debug statement */
        printf("in inp_readall, just read   %s", buffer);
#endif

        if (!buffer)
            continue;

        /* OK -- now we have loaded the next line into 'buffer'.  Process it. */
        /* If input line is blank, ignore it & continue looping.  */
        if ((strcmp(buffer, "\n") == 0) || (strcmp(buffer, "\r\n") == 0))
            if (call_depth != 0 || (call_depth == 0 && cc != NULL)) {
                line_number_orig++;
                tfree(buffer);  /* was allocated by readline() */
                continue;
            }

        if (*buffer == '@') {
            tfree(buffer);      /* was allocated by readline() */
            break;
        }

        /* now handle .title statement */
        if (ciprefix(".title", buffer)) {
            char *s;
            s = skip_non_ws(buffer);               /* skip over .title */
            s = skip_ws(s);            /* advance past space chars */

            /* only the last title line remains valid */
            if (new_title != NULL)
                tfree(new_title);
            new_title = copy(s);
            if ((s = strchr(new_title, '\n')) != NULL)
                *s = ' ';
            *buffer = '*';      /* change .TITLE line to comment line */
        }

        /* now handle .lib statements */
        if (ciprefix(".lib", buffer)) {

            char *y = NULL;     /* filename */
            char *z = NULL;     /* libname */
            char *s, *t;

            inp_stripcomments_line(buffer);
            s = skip_non_ws(buffer);               /* skip over .lib           */

            s = strdup(s);

            t = get_quoted_token(s, &y);

            if (!y) {
                fprintf(cp_err, "Error: .lib filename missing\n");
                tfree(buffer);  /* was allocated by readline()           */
                controlled_exit(EXIT_FAILURE);
            }

            t = get_quoted_token(t, &z);

            if (z && (inp_compat_mode == COMPATMODE_ALL ||
                      inp_compat_mode == COMPATMODE_HS  ||
                      inp_compat_mode == COMPATMODE_NATIVE)) { /* .lib <file name> <lib name> */

                char *copyy = NULL;
                int i;

                if (*y == '~') {
                    copyy = cp_tildexpand(y); /* allocates memory, but can also return NULL */
                    if (copyy)
                        y = copyy; /* reuse y, but remember, buffer still points to allocated memory */
                }

                for (i = 0; i < num_libraries; i++)
                    if (cieq(library_file[i], y))
                        break;

                if (i >= num_libraries) {

                    bool dir_name_flag = FALSE;
                    FILE *newfp = inp_pathopen(y, "r");

                    if (!newfp) {
                        char big_buff2[5000];

                        if (dir_name)
                            sprintf(big_buff2, "%s/%s", dir_name, y);
                        else
                            sprintf(big_buff2, "./%s", y);

                        newfp = inp_pathopen(big_buff2, "r");
                        if (!newfp) {
                            fprintf(cp_err, "Error: Could not find library file %s\n", y);
                            if (copyy)
                                tfree(copyy); /* allocated by the cp_tildexpand() above */
                            tfree(s);
                            tfree(buffer);
                            controlled_exit(EXIT_FAILURE);
                        }

                        dir_name_flag = TRUE;
                    }

                    library_file[num_libraries++] = strdup(y);

                    if (dir_name_flag == FALSE) {
                        char *y_dir_name = ngdirname(y);
                        inp_readall(newfp, &libraries[num_libraries-1], call_depth+1, y_dir_name, FALSE);
                        tfree(y_dir_name);
                    } else {
                        inp_readall(newfp, &libraries[num_libraries-1], call_depth+1, dir_name, FALSE);
                    }

                    fclose(newfp);
                }

                if (copyy)
                    tfree(copyy);   /* allocated by the cp_tildexpand() above */
                tfree(s);

                /* Make the .lib a comment */
                *buffer = '*';
            } else if (inp_compat_mode == COMPATMODE_PS) {
                /* .lib <file name> (no lib name given ) */
                fprintf(cp_err, "Warning: library name missing in line\n  %s", buffer);
                fprintf(cp_err, "  File included as:   .inc %s\n", s);
                memcpy(buffer, ".inc", 4);
            }

        }   /*  end of .lib handling  */

        /* now handle .include statements */
        if (ciprefix(".include", buffer) || ciprefix(".inc", buffer)) {

            char *copyy = NULL;
            char *y = NULL;
            char *s, *t;

            inp_stripcomments_line(buffer);

            s = skip_non_ws(buffer);               /* advance past non-space chars */

            t = get_quoted_token(s, &y);

            if (!y) {
                fprintf(cp_err, "Error: .include filename missing\n");
                tfree(buffer);  /* was allocated by readline() */
                controlled_exit(EXIT_FAILURE);
            }

            if (*y == '~') {
                copyy = cp_tildexpand(y); /* allocates memory, but can also return NULL */
                if (copyy)
                    y = copyy; /* reuse y, but remember, buffer still points to allocated memory */
            }

            {
                bool dir_name_flag = FALSE;
                FILE *newfp = inp_pathopen(y, "r");

                if (!newfp) {
                    char big_buff2[5000];

                    /* open file specified by  .include statement */
                    if (dir_name)
                        sprintf(big_buff2, "%s/%s", dir_name, y);
                    else
                        sprintf(big_buff2, "./%s", y);

                    newfp = inp_pathopen(big_buff2, "r");
                    if (!newfp) {
                        perror(y);
                        if (copyy)
                            tfree(copyy);       /* allocated by the cp_tildexpand() above */
                        fprintf(cp_err, "Error: .include statement failed.\n");
                        tfree(buffer);          /* allocated by readline() above */
                        controlled_exit(EXIT_FAILURE);
                    }

                    dir_name_flag = TRUE;
                }

                if (dir_name_flag == FALSE) {
                    char *y_dir_name = ngdirname(y);
                    inp_readall(newfp, &newcard, call_depth+1, y_dir_name, FALSE);  /* read stuff in include file into netlist */
                    tfree(y_dir_name);
                } else {
                    inp_readall(newfp, &newcard, call_depth+1, dir_name, FALSE);  /* read stuff in include file into netlist */
                }

                (void) fclose(newfp);
            }

            if (copyy)
                tfree(copyy);           /* allocated by the cp_tildexpand() above */

            /* Make the .include a comment */
            *buffer = '*';

            /* now check if this is the first pass (i.e. end points to null) */
            if (end) {                             /* end already exists */
                end->li_next = alloc(struct line);  /* create next card */
                end = end->li_next;                 /* make end point to next card */
            } else {
                end = cc = alloc(struct line);   /* create the deck & end.  cc will
                                                    point to beginning of deck, end to
                                                    the end */
            }

            /* now fill out rest of struct end. */
            end->li_next = NULL;
            end->li_error = NULL;
            end->li_actual = NULL;
            end->li_line = copy(buffer);
            end->li_linenum = end->li_linenum_orig = line_number++;
            if (newcard) {
                end->li_next = newcard;
                /* Renumber the lines */
                line_number_inc = 1;
                for (end = newcard; end && end->li_next; end = end->li_next) {
                    end->li_linenum = line_number++;
                    end->li_linenum_orig = line_number_inc++;
                }
                end->li_linenum = line_number++;        /* SJB - renumber the last line */
                end->li_linenum_orig = line_number_inc++;       /* SJB - renumber the last line */
            }

            /* Fix the buffer up a bit. */
            (void) strncpy(buffer + 1, "end of:", 7);
        }   /*  end of .include handling  */

        /* loop through 'buffer' until end is reached.  Then test for
         premature end.  If premature end is reached, spew
         error and zap the line. */
        if (!ciprefix("write", buffer)) {    // exclude 'write' command so filename case preserved
            char *s;
            for (s = buffer; *s && (*s != '\n'); s++)
                *s = (char) tolower(*s);
            if (!*s) {
                // fprintf(cp_err, "Warning: premature EOF\n");
            }
            *s = '\0';      /* Zap the newline. */

            if ((s-1) >= buffer && *(s-1) == '\r') /* Zop the carriage return under windows */
                *(s-1) = '\0';
        }

        /* find the true .end command out of .endc, .ends, .endl, .end (comments may follow) */
        if (ciprefix(".end", buffer))
            if ((*(buffer+4) == '\0') || (isspace(*(buffer+4)))) {
                found_end = TRUE;
                *buffer   = '*';
            }

        if (shell_eol_continuation) {
            char *new_buffer = TMALLOC(char, strlen(buffer) + 2);
            sprintf(new_buffer, "+%s", buffer);

            tfree(buffer);
            buffer = new_buffer;
        }

        shell_eol_continuation = chk_for_line_continuation(buffer);

        /* now check if this is the first pass (i.e. end points to null) */
        if (end) {                             /* end already exists */
            end->li_next = alloc(struct line); /* create next card */
            end = end->li_next;                /* point to next card */
        } else {                           /* End doesn't exist.  Create it. */
            end = cc = alloc(struct line); /* note that cc points to beginning
                                              of deck, end to the end */
        }

        /* now put buffer into li */
        end->li_next = NULL;
        end->li_error = NULL;
        end->li_actual = NULL;
        end->li_line = copy(buffer);
        end->li_linenum = line_number++;
        end->li_linenum_orig = line_number_orig++;
        tfree(buffer);
    }  /* end while ((buffer = readline(fp)) != NULL) */

    if (!end) { /* No stuff here */
        *data = NULL;
        return;
    }

    if (call_depth == 0 && found_end == TRUE) {
        struct line *prev;
        if (global == NULL) {
            global = TMALLOC(char, strlen(".global gnd") + 1);
            sprintf(global, ".global gnd");
        }
        global_card = alloc(struct line);
        global_card->li_error   = NULL;
        global_card->li_actual  = NULL;
        global_card->li_line    = global;
        global_card->li_linenum = 1;

        prev                 = cc->li_next;
        cc->li_next          = global_card;
        global_card->li_next = prev;

        inp_init_lib_data();
        inp_determine_libraries(cc, NULL);
    }

    /*
       add libraries
    */
    if (call_depth == 0) {
        line_number = expand_libs(line_number);
    }

    /*
      add a terminal ".end" card
    */

    if (call_depth == 0) {
        if (found_end == TRUE) {
            end->li_next = alloc(struct line);    /* create next card */
            end = end->li_next;                   /* point to next card */

            /* now put buffer into li */
            end->li_next = NULL;
            end->li_error = NULL;
            end->li_actual = NULL;
            end->li_line = copy(".end");
            end->li_linenum = end->li_linenum_orig = line_number++;
            end->li_linenum_orig = line_number_orig++;
        }
    }

    /* Replace first line with the new title, if available */
    if (new_title != NULL) {
        if (cc->li_line) tfree(cc->li_line);
        cc->li_line = new_title;
    }

    /* Now clean up li: remove comments & stitch together continuation lines. */
    working = cc->li_next;      /* cc points to head of deck.  Start with the
                                   next card. */

    /* sjb - strip or convert end-of-line comments.
        This must be cone before stitching continuation lines.
        If the line only contains an end-of-line comment then it is converted
        into a normal comment with a '*' at the start.  This will then get
        stripped in the following code. */
    inp_stripcomments_deck(working);

    prev = NULL;
    while (working) {
        char *s;
        for (s = working->li_line; (c = *s) != '\0' && c <= ' '; s++)
            ;

#ifdef TRACE
        /* SDB debug statement */
        printf("In inp_readall, processing linked list element line = %d, s = %s . . . \n", working->li_linenum, s);
#endif

        switch (c) {
        case '#':
        case '$':
        case '*':
        case '\0':
            /* this used to be commented out.  Why? */
            /*  prev = NULL; */
            working = working->li_next;  /* for these chars, go to next card */
            break;

        case '+':   /* handle continuation */
            if (!prev) {
                working->li_error = copy("Illegal continuation line: ignored.");
                working = working->li_next;
                break;
            }

            /* create buffer and write last and current line into it. */
            buffer = TMALLOC(char, strlen(prev->li_line) + strlen(s) + 2);
            (void) sprintf(buffer, "%s %s", prev->li_line, s + 1);

            s = prev->li_line;
            prev->li_line = buffer;
            prev->li_next = working->li_next;
            working->li_next = NULL;
            if (prev->li_actual) {
                for (end = prev->li_actual; end->li_next; end = end->li_next)
                    ;
                end->li_next = working;
                tfree(s);
            } else {
                newcard = alloc(struct line);
                newcard->li_linenum = prev->li_linenum;
                newcard->li_line = s;
                newcard->li_next = working;
                newcard->li_error = NULL;
                newcard->li_actual = NULL;
                prev->li_actual = newcard;
            }
            working = prev->li_next;
            break;

        default:  /* regular one-line card */
            prev = working;
            working = working->li_next;
            break;
        }
    }

    /* The following processing of an input file is not required for command files
       like spinit or .spiceinit, so return command files here. */
    if (comfile) {
        /* save the return value (via **data) */
        *data = cc;
        return;
    }

    working = cc->li_next;

    inp_fix_for_numparam(working);
    inp_remove_excess_ws(working);

    if (call_depth == 0) {

        comment_out_unused_subckt_models(working, line_number);

        line_number = inp_split_multi_param_lines(working, line_number);

        inp_fix_macro_param_func_paren_io(working);
        inp_fix_ternary_operator(working);
        inp_grab_func(working);

        inp_expand_macros_in_func();
        inp_expand_macros_in_deck(working);
        inp_fix_param_values(working);

        /* get end card as last card in list; end card pntr does not appear to always
           be correct at this point */
        for (newcard = working; newcard != NULL; newcard = newcard->li_next)
            end = newcard;

        inp_reorder_params(working, cc, end);
        inp_fix_inst_calls_for_numparam(working);
        inp_fix_gnd_name(working);
        inp_chk_for_multi_in_vcvs(working, &line_number);

        if (cp_getvar("addcontrol", CP_BOOL, NULL))
            inp_add_control_section(working, &line_number);
#ifndef XSPICE
        inp_poly_err(working);
#endif
        if (inp_compat_mode != COMPATMODE_SPICE3) {
            /* Do all the compatibility stuff here */
            working = cc->li_next;
            /* E, G, L, R, C compatibility transformations */
            inp_compat(working);
            working = cc->li_next;
            /* B source numparam compatibility transformation */
            inp_bsource_compat(working);
        }

        inp_add_series_resistor(working);
    }

    /* save the return value (via **data) */
    *data = cc;

    /* get max. line length and number of lines in input deck,
       and renumber the lines,
       count the number of '{' per line as an upper estimate of the number
       of parameter substitutions in a line*/
    dynmaxline = 0;
    max_line_length = 0;
    for (tmp_ptr1 = cc; tmp_ptr1 != NULL; tmp_ptr1 = tmp_ptr1->li_next) {
        char *s;
        unsigned int braces_per_line = 0;
        /* count number of lines */
        dynmaxline++;
        /* renumber the lines of the processed input deck */
        tmp_ptr1->li_linenum = dynmaxline;
        if (max_line_length < strlen(tmp_ptr1->li_line))
            max_line_length = strlen(tmp_ptr1->li_line);
        /* count '{' */
        for (s = tmp_ptr1->li_line; *s; s++)
            if (*s == '{')
                braces_per_line++;
        if (no_braces <  braces_per_line)
            no_braces = braces_per_line;
    }

    if (ft_ngdebug) {
        /*debug: print into file*/
        fdo = fopen("debug-out.txt", "w");
        for (tmp_ptr1 = cc; tmp_ptr1 != NULL; tmp_ptr1 = tmp_ptr1->li_next)
            fprintf(fdo, "%d  %d  %s\n", tmp_ptr1->li_linenum_orig, tmp_ptr1->li_linenum, tmp_ptr1->li_line);

        (void) fclose(fdo);
        fprintf(stdout, "max line length %d, max subst. per line %d, number of lines %d\n",
                (int) max_line_length, no_braces, dynmaxline);
    }
}


/*-------------------------------------------------------------------------*
  Look up the variable sourcepath and try everything in the list in order
  if the file isn't in . and it isn't an abs path name.
  For MS Windows: First try the path of the source file.
 *-------------------------------------------------------------------------*/
FILE *
inp_pathopen(char *name, char *mode)
{
    FILE *fp;
    char buf[BSIZE_SP];
    struct variable *v;

#if defined(HAS_WINDOWS)
    char buf2[BSIZE_SP];
    /* search in the path where the source (input) file has been found,
       but only if "name" is just a file name */
    if (!strchr(name, DIR_TERM) && !strchr(name, DIR_TERM_LINUX) && cp_getvar("sourcefile", CP_STRING, buf2)) {
        /* If pathname is found, get path.
        (char *dirname(const char *name) might have been used here) */
        if (substring(DIR_PATHSEP, buf2) || substring(DIR_PATHSEP_LINUX, buf2)) {
            int i, j = 0;
            for (i = 0; i < BSIZE_SP-1; i++) {
                if (buf2[i] == '\0')
                    break;
                if ((buf2[i] == DIR_TERM) || (buf2[i] == DIR_TERM_LINUX))
                    j = i;
            }
            buf2[j+1] = '\0'; /* include path separator */
        }
        /* add file name */
        strcat(buf2, name);
        /* try to open file */
        if ((fp = fopen(buf2, mode)) != NULL)
            return (fp);
    }
    /* If this is an abs pathname, or there is no sourcepath var, just
     * do an fopen.
     */
    if (strchr(name, DIR_TERM) || strchr(name, DIR_TERM_LINUX) ||
        !cp_getvar("sourcepath", CP_LIST, &v))
        return (fopen(name, mode));
#else


    /* If this is an abs pathname, or there is no sourcepath var, just
     * do an fopen.
     */
    if (strchr(name, DIR_TERM) || !cp_getvar("sourcepath", CP_LIST, &v))
        return (fopen(name, mode));
#endif

    while (v) {
        switch (v->va_type) {
        case CP_STRING:
            cp_wstrip(v->va_string);
            (void) sprintf(buf, "%s%s%s", v->va_string, DIR_PATHSEP, name);
            break;
        case CP_NUM:
            (void) sprintf(buf, "%d%s%s", v->va_num, DIR_PATHSEP, name);
            break;
        case CP_REAL:           /* This is foolish */
            (void) sprintf(buf, "%g%s%s", v->va_real, DIR_PATHSEP, name);
            break;
        default: {
            fprintf(stderr, "ERROR: enumeration value `CP_BOOL' or `CP_LIST' not handled in inp_pathopen\nAborting...\n");
            controlled_exit(EXIT_FAILURE);
        }
        }
        if ((fp = fopen(buf, mode)) != NULL)
            return (fp);
        v = v->va_next;
    }
    return (NULL);
}

/*-------------------------------------------------------------------------*
 *  This routine reads a line (of arbitrary length), up to a '\n' or 'EOF' *
 *  and returns a pointer to the resulting null terminated string.         *
 *  The '\n' if found, is included in the returned string.                 *
 *  From: jason@ucbopal.BERKELEY.EDU (Jason Venner)                        *
 *  Newsgroups: net.sources                                                *
 *-------------------------------------------------------------------------*/
#define STRGROW 256

static char *
readline(FILE *fd)
{
    int c;
    int memlen;
    char *strptr;
    int strlen;

    strptr = NULL;
    strlen = 0;
    memlen = STRGROW;
    strptr = TMALLOC(char, memlen);
    memlen -= 1;                /* Save constant -1's in while loop */
    while ((c = getc(fd)) != EOF) {
        if (strlen == 0 && (c == '\t' || c == ' ')) /* Leading spaces away */
            continue;
        strptr[strlen] = (char) c;
        strlen++;
        if (strlen >= memlen) {
            memlen += STRGROW;
            if ((strptr = TREALLOC(char, strptr, memlen + 1)) == NULL)
                return (NULL);
        }
        if (c == '\n')
            break;
    }
    if (!strlen) {
        tfree(strptr);
        return (NULL);
    }
    // strptr[strlen] = '\0';
    /* Trim the string */
    strptr = TREALLOC(char, strptr, strlen + 1);
    strptr[strlen] = '\0';
    return (strptr);
}


/* replace "gnd" by " 0 "
   Delimiters of gnd may be ' ' or ',' or '(' or ')' */
static void
inp_fix_gnd_name(struct line *deck)
{
    struct line *c = deck;
    char *gnd;

    while (c != NULL) {
        gnd = c->li_line;
        // if there is a comment or no gnd, go to next line
        if ((*gnd == '*') || (strstr(gnd, "gnd") == NULL)) {
            c = c->li_next;
            continue;
        }
        // replace "?gnd?" by "? 0 ?", ? being a ' '  ','  '('  ')'.
        while ((gnd = strstr(gnd, "gnd")) != NULL) {
            if ((isspace(*(gnd-1)) || *(gnd-1) == '(' || *(gnd-1) == ',') &&
                (isspace(*(gnd+3)) || *(gnd+3) == ')' || *(gnd+3) == ',')) {
                memcpy(gnd, " 0 ", 3);
            }
            gnd += 3;
        }
        // now remove the extra white spaces around 0
        c->li_line = inp_remove_ws(c->li_line);
        c = c->li_next;
    }
}

static struct line*
create_new_card(char *card_str, int *line_number) {
    char        *str = strdup(card_str);
    struct line *newcard = alloc(struct line);

    newcard->li_line    = str;
    newcard->li_linenum = *line_number;
    newcard->li_error   = NULL;
    newcard->li_actual  = NULL;

    *line_number = *line_number + 1;

    return newcard;
}


static void
inp_chk_for_multi_in_vcvs(struct line *deck, int *line_number)
{
    struct line *c, *a_card, *model_card, *next_card;
    char *line, *bool_ptr, *str_ptr1, *str_ptr2, keep, *comma_ptr, *xy_values1[5], *xy_values2[5];
    char *node_str, *ctrl_node_str, *xy_str1, *model_name, *fcn_name;
    char big_buf[1000];
    int  xy_count1 = 0, xy_count2 = 0, skip_control = 0;

    for (c = deck; c != NULL; c = c->li_next) {
        line = c->li_line;

        /* there is no e source inside .control ... .endc */
        if (ciprefix(".control", line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }

        if (*line == 'e') {
            if ((bool_ptr = strstr(line, "nand(")) != NULL ||
                (bool_ptr = strstr(line, "and(")) != NULL ||
                (bool_ptr = strstr(line, "nor(")) != NULL ||
                (bool_ptr = strstr(line, "or(")) != NULL)
            {
                str_ptr1 = skip_non_ws(line);
                model_name = copy_substring(line, str_ptr1);

                str_ptr1 = skip_ws(str_ptr1);
                str_ptr2 = skip_back_ws(bool_ptr - 1) + 1;
                keep = *str_ptr2;
                *str_ptr2 = '\0';
                node_str  = strdup(str_ptr1);
                *str_ptr2 = keep;

                str_ptr1 = bool_ptr + 1;
                while (*str_ptr1 != '(')
                    str_ptr1++;
                fcn_name = copy_substring(bool_ptr, str_ptr1);
                str_ptr1  = strchr(str_ptr1, ')');
                comma_ptr = str_ptr2 = strchr(line, ',');
                if ((str_ptr1 == NULL)|| (str_ptr1 == NULL)) {
                    fprintf(stderr, "ERROR: mal formed line: %s\n", line);
                    controlled_exit(EXIT_FAILURE);
                }
                str_ptr2 = skip_back_ws(str_ptr2 - 1);
                str_ptr1 = skip_ws(str_ptr1 + 1);
                if (*str_ptr2 == '}') {
                    while (*str_ptr2 != '{')
                        str_ptr2--;
                    xy_str1 = str_ptr2;
                    str_ptr2 = skip_back_ws(str_ptr2 - 1) + 1;
                } else {
                    str_ptr2 = skip_back_non_ws(str_ptr2);
                    xy_str1 = str_ptr2 + 1;
                    str_ptr2 = skip_back_ws(str_ptr2) + 1;
                }
                keep = *str_ptr2;
                *str_ptr2 = '\0';
                ctrl_node_str = strdup(str_ptr1);
                *str_ptr2 = keep;

                str_ptr1 = skip_ws(comma_ptr + 1);
                if (*str_ptr1 == '{') {
                    while (*str_ptr1 != '}')
                        str_ptr1++;
                    str_ptr1++;
                } else {
                    str_ptr1 = skip_non_ws(str_ptr1);
                }
                keep = *str_ptr1;
                *str_ptr1 = '\0';
                xy_count1 = get_comma_separated_values(xy_values1, xy_str1);
                *str_ptr1 = keep;

                str_ptr1 = skip_ws(str_ptr1);
                xy_count2 = get_comma_separated_values(xy_values2, str_ptr1);

                // place restrictions on only having 2 point values; this can change later
                if (xy_count1 != 2 && xy_count2 != 2)
                    fprintf(stderr, "ERROR: only expecting 2 pair values for multi-input vcvs!\n");

                sprintf(big_buf, "%s %%vd[ %s ] %%vd( %s ) %s",
                        model_name, ctrl_node_str, node_str, model_name);
                a_card = create_new_card(big_buf, line_number);
                *a_card->li_line = 'a';

                sprintf(big_buf, ".model %s multi_input_pwl ( x = [%s %s] y = [%s %s] model = \"%s\" )",
                         model_name, xy_values1[0], xy_values2[0],
                         xy_values1[1], xy_values2[1], fcn_name);
                model_card = create_new_card(big_buf, line_number);

                tfree(model_name);
                tfree(node_str);
                tfree(fcn_name);
                tfree(ctrl_node_str);
                tfree(xy_values1[0]);
                tfree(xy_values1[1]);
                tfree(xy_values2[0]);
                tfree(xy_values2[1]);

                *c->li_line = '*';
                next_card   = c->li_next;
                c->li_next  = a_card;
                a_card->li_next     = model_card;
                model_card->li_next = next_card;
            }
        }
    }
}


static void
inp_add_control_section(struct line *deck, int *line_number)
{
    struct line *c, *newcard, *prev_card = NULL;
    bool        found_control = FALSE, found_run = FALSE;
    bool        found_end = FALSE;
    char        *op_line  = NULL, rawfile[1000], *line;

    for (c = deck; c != NULL; c = c->li_next) {
        if (*c->li_line == '*')
            continue;
        if (ciprefix(".op ", c->li_line)) {
            *c->li_line = '*';
            op_line = c->li_line + 1;
        }
        if (ciprefix(".end", c->li_line))
            found_end = TRUE;
        if (found_control && ciprefix("run", c->li_line))
            found_run = TRUE;

        if (ciprefix(".control", c->li_line))
            found_control = TRUE;
        if (ciprefix(".endc", c->li_line)) {
            found_control = FALSE;

            if (!found_run) {
                newcard            = create_new_card("run", line_number);
                prev_card->li_next = newcard;
                newcard->li_next   = c;
                prev_card          = newcard;
                found_run          = TRUE;
            }
            if (cp_getvar("rawfile", CP_STRING, rawfile)) {
                line = TMALLOC(char, strlen("write") + strlen(rawfile) + 2);
                sprintf(line, "write %s", rawfile);
                newcard            = create_new_card(line, line_number);
                prev_card->li_next = newcard;
                newcard->li_next   = c;
                prev_card          = newcard;
                tfree(line);
            }
        }
        prev_card = c;
    }
    // check if need to add control section
    if (!found_run && found_end) {
        prev_card        = deck->li_next;
        newcard          = create_new_card(".endc", line_number);
        deck->li_next    = newcard;
        newcard->li_next = prev_card;

        if (cp_getvar("rawfile", CP_STRING, rawfile)) {
            line = TMALLOC(char, strlen("write") + strlen(rawfile) + 2);
            sprintf(line, "write %s", rawfile);
            prev_card        = deck->li_next;
            newcard          = create_new_card(line, line_number);
            deck->li_next    = newcard;
            newcard->li_next = prev_card;
            tfree(line);
        }
        if (op_line != NULL) {
            prev_card        = deck->li_next;
            newcard          = create_new_card(op_line, line_number);
            deck->li_next    = newcard;
            newcard->li_next = prev_card;
        }

        prev_card        = deck->li_next;
        newcard          = create_new_card("run", line_number);
        deck->li_next    = newcard;
        newcard->li_next = prev_card;

        prev_card        = deck->li_next;
        newcard          = create_new_card(".control", line_number);
        deck->li_next    = newcard;
        newcard->li_next = prev_card;
    }
}


// look for shell-style end-of-line continuation '\\'
static bool
chk_for_line_continuation(char *line)
{
    char *ptr = line + strlen(line) - 1;

    if (*line != '*' && *line != '$') {
        while (ptr >= line && *ptr && isspace(*ptr))
            ptr--;

        if ((ptr-1) >= line && *ptr == '\\' && *(ptr-1) && *(ptr-1) == '\\') {
            *ptr = ' ';
            *(ptr-1) = ' ';
            return TRUE;
        }
    }

    return FALSE;
}


//
// change .macro --> .subckt
//        .eom   --> .ends
//        .subckt name 1 2 3 params: w=9u l=180n --> .subckt name 1 2 3 w=9u l=180n
//        .subckt name (1 2 3) --> .subckt name 1 2 3
//        x1 (1 2 3)      --> x1 1 2 3
//        .param func1(x,y) = {x*y} --> .func func1(x,y) {x*y}

static void
inp_fix_macro_param_func_paren_io(struct line *begin_card)
{
    struct line *card;
    char        *str_ptr, *new_str;
    bool        is_func = FALSE;

    for (card = begin_card; card != NULL; card = card->li_next) {

        if (*card->li_line == '*')
            continue;

        if (ciprefix(".macro", card->li_line) || ciprefix(".eom", card->li_line)) {
            str_ptr = skip_non_ws(card->li_line);

            if (ciprefix(".macro", card->li_line)) {
                new_str = TMALLOC(char, strlen(".subckt") + strlen(str_ptr) + 1);
                sprintf(new_str, ".subckt%s", str_ptr);
            } else {
                new_str = TMALLOC(char, strlen(".ends") + strlen(str_ptr) + 1);
                sprintf(new_str, ".ends%s", str_ptr);
            }

            tfree(card->li_line);
            card->li_line = new_str;
        }

        if (ciprefix(".subckt", card->li_line) || ciprefix("x", card->li_line)) {
            /* remove () */
            str_ptr = skip_non_ws(card->li_line);  // skip over .subckt, instance name
            str_ptr = skip_ws(str_ptr);
            if (ciprefix(".subckt", card->li_line)) {
                str_ptr = skip_non_ws(str_ptr);  // skip over subckt name
                str_ptr = skip_ws(str_ptr);
            }
            if (*str_ptr == '(') {
                *str_ptr = ' ';
                while (*str_ptr && *str_ptr != '\0') {
                    if (*str_ptr == ')') {
                        *str_ptr = ' ';
                        break;
                    }
                    str_ptr++;
                }
                card->li_line = inp_remove_ws(card->li_line); /* remove the extra white spaces just introduced */
            }
        }
        is_func = FALSE;
        if (ciprefix(".param", card->li_line)) {
            str_ptr = skip_non_ws(card->li_line);  // skip over .param
            str_ptr = skip_ws(str_ptr);
            while (!isspace(*str_ptr) && *str_ptr != '=') {
                if (*str_ptr == '(')
                    is_func = TRUE;
                str_ptr++;
            }

            if (is_func) {
                str_ptr = strchr(card->li_line, '=');
                if (str_ptr)
                    *str_ptr = ' ';
                str_ptr = card->li_line + 1;
                *str_ptr = 'f';
                *(str_ptr+1) = 'u';
                *(str_ptr+2) = 'n';
                *(str_ptr+3) = 'c';
                *(str_ptr+4) = ' ';
            }
        }
    }
}


static char *
get_instance_subckt(char *line)
{
    char *end_ptr, *inst_name_ptr;
    char *equal_ptr = strchr(line, '=');

    // see if instance has parameters
    if (equal_ptr) {
        end_ptr = skip_back_ws(equal_ptr - 1);
        end_ptr = skip_back_non_ws(end_ptr) + 1;
    } else {
        end_ptr = line + strlen(line);
    }

    end_ptr = skip_back_ws(end_ptr - 1) + 1;

    inst_name_ptr = skip_back_non_ws(end_ptr - 1) + 1;
    return copy_substring(inst_name_ptr, end_ptr);
}


static char*
get_subckt_model_name(char *line)
{
    char *name, *end_ptr;

    name = skip_non_ws(line);   // eat .subckt|.model
    name = skip_ws(name);

    end_ptr = skip_non_ws(name);
    return copy_substring(name, end_ptr);
}


static char*
get_model_name(char *line, int num_terminals)
{
    char *beg_ptr, *end_ptr;
    int  i = 0;

    beg_ptr = skip_non_ws(line); /* eat device name */
    beg_ptr = skip_ws(beg_ptr);

    for (i = 0; i < num_terminals; i++) { /* skip the terminals */
        beg_ptr = skip_non_ws(beg_ptr);
        beg_ptr = skip_ws(beg_ptr);
    }
    if (*line == 'r')           /* special dealing for r models */
        if ((*beg_ptr == '+') || (*beg_ptr == '-') || isdigit(*beg_ptr)) { /* looking for a value before model */
            beg_ptr = skip_non_ws(beg_ptr); /* skip the value */
            beg_ptr = skip_ws(beg_ptr);
        }

    end_ptr = skip_non_ws(beg_ptr);
    return copy_substring(beg_ptr, end_ptr);
}


static char*
get_model_type(char *line)
{
    char *model_type, *beg_ptr;
    if (!(ciprefix(".model", line)))
        return NULL;
    beg_ptr = skip_non_ws(line); /* eat .model */
    beg_ptr = skip_ws(beg_ptr);
    beg_ptr = skip_non_ws(beg_ptr); /* eat model name */
    beg_ptr = skip_ws(beg_ptr);
    model_type = gettok(&beg_ptr);
    return model_type;
}


static char *
get_adevice_model_name(char *line)
{
    char *ptr_end, *ptr_beg;

    ptr_end = skip_back_ws(line + strlen(line) - 1) + 1;
    ptr_beg = skip_back_non_ws(ptr_end - 1) + 1;
    return copy_substring(ptr_beg, ptr_end);
}


static void
get_subckts_for_subckt(struct line *start_card, char *subckt_name,
                        char *used_subckt_names[], int *num_used_subckt_names,
                        char *used_model_names[], int *num_used_model_names,
                        bool has_models)
{
    struct line *card;
    char *line = NULL, *curr_subckt_name, *inst_subckt_name, *model_name, *new_names[100];
    bool found_subckt = FALSE, have_subckt = FALSE, found_model = FALSE;
    int  i, num_terminals = 0, tmp_cnt = 0;

    for (card = start_card; card != NULL; card = card->li_next) {
        line = card->li_line;

        if (*line == '*')
            continue;

        if ((ciprefix(".ends", line) || ciprefix(".eom", line)) && found_subckt)
            break;

        if (ciprefix(".subckt", line) || ciprefix(".macro", line)) {
            curr_subckt_name = get_subckt_model_name(line);

            if (strcmp(curr_subckt_name, subckt_name) == 0)
                found_subckt = TRUE;

            tfree(curr_subckt_name);
        }
        if (found_subckt) {
            if (*line == 'x') {
                inst_subckt_name = get_instance_subckt(line);
                have_subckt = FALSE;
                for (i = 0; i < *num_used_subckt_names; i++)
                    if (strcmp(used_subckt_names[i], inst_subckt_name) == 0)
                        have_subckt = TRUE;
                if (!have_subckt) {
                    new_names[tmp_cnt++] = used_subckt_names[*num_used_subckt_names] = inst_subckt_name;
                    *num_used_subckt_names = *num_used_subckt_names + 1;
                } else {
                    tfree(inst_subckt_name);
                }
            } else if (*line == 'a') {
                model_name = get_adevice_model_name(line);
                found_model = FALSE;
                for (i = 0; i < *num_used_model_names; i++)
                    if (strcmp(used_model_names[i], model_name) == 0)
                        found_model = TRUE;
                if (!found_model) {
                    used_model_names[*num_used_model_names] = model_name;
                    *num_used_model_names = *num_used_model_names + 1;
                } else {
                    tfree(model_name);
                }
            } else if (has_models) {
                num_terminals = get_number_terminals(line);

                if (num_terminals != 0) {
                    char *tmp_name, *tmp_name1;
                    tmp_name1 = tmp_name = model_name = get_model_name(line, num_terminals);

                    if (isalpha(*model_name) ||
                            /* first character is digit, second is alpha, third is digit,
                            e.g. 1N4002 */
                            ((strlen(model_name) > 2) && isdigit(*tmp_name) &&
                             isalpha(*(++tmp_name)) && isdigit(*(++tmp_name))) ||
                            /* first character is is digit, second is alpha, third is alpha, fourth is digit
                            e.g. 2SK456 */
                            ((strlen(model_name) > 3) && isdigit(*tmp_name1) && isalpha(*(++tmp_name1)) &&
                             isalpha(*(++tmp_name1)) && isdigit(*(++tmp_name1)))) {
                        found_model = FALSE;
                        for (i = 0; i < *num_used_model_names; i++)
                            if (strcmp(used_model_names[i], model_name) == 0) found_model = TRUE;
                        if (!found_model) {
                            used_model_names[*num_used_model_names] = model_name;
                            *num_used_model_names = *num_used_model_names + 1;
                        } else {
                            tfree(model_name);
                        }
                    } else {
                        tfree(model_name);
                    }
                }
            }
        }
    }
    // now make recursive call on instances just found above
    for (i = 0; i < tmp_cnt; i++)
        get_subckts_for_subckt(start_card, new_names[i], used_subckt_names, num_used_subckt_names,
                                used_model_names, num_used_model_names, has_models);
}


/*
  check if current token matches model bin name -- <token>.[0-9]+
 */
static bool
model_bin_match(char *token, char *model_name)
{
    char *dot_char;
    bool  flag = FALSE;

    if (strncmp(model_name, token, strlen(token)) == 0)
        if ((dot_char = strchr(model_name, '.')) != NULL) {
            flag = TRUE;
            dot_char++;
            while (*dot_char != '\0') {
                if (!isdigit(*dot_char)) {
                    flag = FALSE;
                    break;
                }
                dot_char++;
            }
        }

    return flag;
}


/*
  iterate through the deck and comment out unused subckts, models
  (don't want to waste time processing everything)
  also comment out .param lines with no parameters defined
 */
static void
comment_out_unused_subckt_models(struct line *start_card, int no_of_lines)
{
    struct line *card;
    char **used_subckt_names, **used_model_names, *line = NULL, *subckt_name, *model_name;
    int  num_used_subckt_names = 0, num_used_model_names = 0, i = 0, num_terminals = 0, tmp_cnt = 0;
    bool processing_subckt = FALSE, found_subckt = FALSE, remove_subckt = FALSE, found_model = FALSE, has_models = FALSE;
    int skip_control = 0, nested_subckt = 0;

    /* generate arrays of *char for subckt or model names. Start
    with 1000, but increase, if number of lines in deck is larger */
    if (no_of_lines < 1000)
        no_of_lines = 1000;
    used_subckt_names = TMALLOC(char*, no_of_lines);
    used_model_names = TMALLOC(char*, no_of_lines);

    for (card = start_card; card != NULL; card = card->li_next) {
        if (ciprefix(".model", card->li_line))
            has_models = TRUE;
        if (ciprefix(".cmodel", card->li_line))
            has_models = TRUE;
        if (ciprefix(".param", card->li_line) && !strchr(card->li_line, '='))
            *card->li_line = '*';
    }

    for (card = start_card; card != NULL; card = card->li_next) {
        line = card->li_line;

        if (*line == '*')
            continue;

        /* there is no .subckt, .model or .param inside .control ... .endc */
        if (ciprefix(".control", line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }

        if (ciprefix(".subckt", line) || ciprefix(".macro", line))
            processing_subckt = TRUE;
        if (ciprefix(".ends", line) || ciprefix(".eom", line))
            processing_subckt = FALSE;
        if (!processing_subckt) {
            if (*line == 'x') {
                subckt_name = get_instance_subckt(line);
                found_subckt = FALSE;
                for (i = 0; i < num_used_subckt_names; i++)
                    if (strcmp(used_subckt_names[i], subckt_name) == 0) found_subckt = TRUE;
                if (!found_subckt) {
                    used_subckt_names[num_used_subckt_names++] = subckt_name;
                    tmp_cnt++;
                } else {
                    tfree(subckt_name);
                }
            } else if (*line == 'a') {
                model_name = get_adevice_model_name(line);
                found_model = FALSE;
                for (i = 0; i < num_used_model_names; i++)
                    if (strcmp(used_model_names[i], model_name) == 0)
                        found_model = TRUE;
                if (!found_model)
                    used_model_names[num_used_model_names++] = model_name;
                else
                    tfree(model_name);
            } else if (has_models) {
                /* This is a preliminary version, until we have found a reliable
                method to detect the model name out of the input line (Many
                options have to be taken into account.). */
                num_terminals = get_number_terminals(line);
                if (num_terminals != 0) {
                    bool model_ok = FALSE;
                    char *tmp_name, *tmp_name1;
                    tmp_name = tmp_name1 = model_name = get_model_name(line, num_terminals);
                    /* first character of model name is character from alphabet */
                    if (isalpha(*model_name))
                        model_ok = TRUE;
                    /* first character is digit, second is alpha, third is digit,
                    e.g. 1N4002 */
                    else if ((strlen(model_name) > 2) && isdigit(*tmp_name) &&
                             isalpha(*(++tmp_name)) && isdigit(*(++tmp_name)))
                        model_ok = TRUE;
                    /* first character is is digit, second is alpha, third is alpha, fourth is digit
                    e.g. 2SK456 */
                    else if ((strlen(model_name) > 3) && isdigit(*tmp_name1) &&
                             isalpha(*(++tmp_name1)) && isalpha(*(++tmp_name1)) &&
                             isdigit(*(++tmp_name1)))
                        model_ok = TRUE;
                    /* Check if model has already been recognized, if not, add its name to
                    list used_model_names[i] */
                    if (model_ok) {
                        found_model = FALSE;
                        for (i = 0; i < num_used_model_names; i++)
                            if (strcmp(used_model_names[i], model_name) == 0)
                                found_model = TRUE;
                        if (!found_model)
                            used_model_names[num_used_model_names++] = model_name;
                        else
                            tfree(model_name);
                    } else {
                        tfree(model_name);
                    }
                }
            } /* if (has_models)  */
        } /* if (!processing_subckt) */
    } /* for loop through all cards */
    for (i = 0; i < tmp_cnt; i++)
        get_subckts_for_subckt
            (start_card, used_subckt_names[i],
             used_subckt_names, &num_used_subckt_names,
             used_model_names, &num_used_model_names, has_models);

    /* comment out any unused subckts, currently only at top level */
    for (card = start_card; card != NULL; card = card->li_next) {
        line = card->li_line;

        if (*line == '*')
            continue;

        if (ciprefix(".subckt", line) || ciprefix(".macro", line)) {
            nested_subckt++;
            subckt_name = get_subckt_model_name(line);
            if (nested_subckt == 1) {
                /* check if unused, only at top level */
                remove_subckt = TRUE;
                for (i = 0; i < num_used_subckt_names; i++)
                    if (strcmp(used_subckt_names[i], subckt_name) == 0)
                        remove_subckt = FALSE;
            }
            tfree(subckt_name);
        }
        if (ciprefix(".ends", line) || ciprefix(".eom", line)) {
            nested_subckt--;
            if (remove_subckt)
                *line = '*';
            if (nested_subckt == 0)
                remove_subckt = FALSE;
        }
        if (remove_subckt)
            *line = '*';
        else if (has_models &&
                 (ciprefix(".model", line) || ciprefix(".cmodel", line)))
        {
            char *model_type = get_model_type(line);
            model_name = get_subckt_model_name(line);
            /* keep R, L, C models because in addition to no. of terminals the value may be given,
               as in RE1 1 2 800 newres dtemp=5, so model name may be token no. 4 or 5,
               and, if 5, will not be detected by get_subckt_model_name()*/
            if (cieq(model_type, "c") ||
                cieq(model_type, "l") ||
                cieq(model_type, "r"))
            {
                found_model = TRUE;
            } else {
                found_model = FALSE;
                for (i = 0; i < num_used_model_names; i++)
                    if (strcmp(used_model_names[i], model_name) == 0 ||
                        model_bin_match(used_model_names[i], model_name))
                        found_model = TRUE;
            }
            if (model_type)
                tfree(model_type);
            if (!found_model)
                *line = '*';
            tfree(model_name);
        }
    }
    for (i = 0; i < num_used_subckt_names; i++)
        tfree(used_subckt_names[i]);
    for (i = 0; i < num_used_model_names;  i++)
        tfree(used_model_names[i]);
    tfree(used_subckt_names);
    tfree(used_model_names);
}


/* replace ternary operator ? : by fcn ternary_fcn() in .param, .func, and .meas lines,
   if all is FALSE, for all lines if all is TRUE */
static char*
inp_fix_ternary_operator_str(char *line, bool all)
{
    char *conditional, *if_str, *else_str, *question, *colon, keep, *str_ptr, *str_ptr2, *new_str;
    char *paren_ptr = NULL, *end_str = NULL, *beg_str = NULL;
    int  count = 0;

    if (!strchr(line, '?') && !strchr(line, ':'))
        return line;

    if (all || ciprefix(".param", line) ||
        ciprefix(".func", line) || ciprefix(".meas", line))
    {
        if (ciprefix(".param", line) || ciprefix(".meas", line))
            str_ptr = strchr(line, '=');
        else
            str_ptr = strchr(line, ')');

        if ((str_ptr == NULL) && all == FALSE) {
            fprintf(stderr, "ERROR: mal formed .param, .func or .meas line:\n  %s\n", line);
            controlled_exit(EXIT_FAILURE);
        }

        if ((str_ptr == NULL) && all == TRUE) {
            fprintf(stderr, "ERROR: mal formed expression in line:\n  %s\n", line);
            fprintf(stderr, "  We need parentheses around 'if' clause and nested ternary functions\n");
            fprintf(stderr, "  like:  Rtern4 1 0 '(ut > 0.7) ? 2k : ((ut < 0.3) ? 500 : 1k)'\n");
            controlled_exit(EXIT_FAILURE);
        }

        str_ptr = skip_ws(str_ptr + 1);
        if (*str_ptr == '{')
            str_ptr = skip_ws(str_ptr + 1);
    } else {
        return line;
    }

    // get conditional
    question = strchr(str_ptr, '?');
    str_ptr2 = skip_back_ws(question - 1);
    if (*str_ptr2 == ')') {
        count = 1;
        str_ptr = str_ptr2;
        while ((count != 0) && (str_ptr != line)) {
            str_ptr--;
            if (*str_ptr == '(')
                count--;
            if (*str_ptr == ')')
                count++;
        }
    }
    str_ptr2++;
    keep = *str_ptr2;
    *str_ptr2 = '\0';
    conditional = strdup(str_ptr);
    *str_ptr2 = keep;

    // get beginning (whatever is left before the conditional)
    keep = *str_ptr;
    *str_ptr = '\0';
    beg_str = strdup(line);
    *str_ptr = keep;

    // get if
    str_ptr = skip_ws(question + 1);
    if (*str_ptr == '(') {
        // find closing paren
        count    = 1;
        str_ptr2 = str_ptr/* + 1*/;
        while (count != 0 && *str_ptr2 != '\0') {
            str_ptr2++;
            if (*str_ptr2 == '(')
                count++;
            if (*str_ptr2 == ')')
                count--;
        }
        if (count != 0) {
            fprintf(stderr, "ERROR: problem parsing 'if' of ternary string %s!\n", line);
            controlled_exit(EXIT_FAILURE);
        }
        colon = str_ptr2 + 1;
        while (*colon != ':' && *colon != '\0')
            colon++;
        if (*colon != ':') {
            fprintf(stderr, "ERROR: problem parsing ternary string (finding ':') %s!\n", line);
            controlled_exit(EXIT_FAILURE);
        }
        str_ptr2 = skip_back_ws(colon - 1);
    } else if ((colon = strchr(str_ptr, ':')) != NULL) {
        str_ptr2 = skip_back_ws(colon - 1);
    } else {
        fprintf(stderr, "ERROR: problem parsing ternary string (missing ':') %s!\n", line);
        controlled_exit(EXIT_FAILURE);
    }
    str_ptr2++;
    keep = *str_ptr2;
    *str_ptr2 = '\0';
    if_str = inp_fix_ternary_operator_str(strdup(str_ptr), all);
    *str_ptr2 = keep;

    // get else
    str_ptr = skip_ws(colon + 1);
    paren_ptr = strchr(question, '(');
    if (paren_ptr) {
        // find end paren ')'
        bool found_paren = FALSE;
        count = 0;
        str_ptr2 = str_ptr;
        while (*str_ptr2 != '\0') {
            if (*str_ptr2 == '(') {
                count++;
                found_paren = TRUE;
            }
            if (*str_ptr2 == ')')
                count--;
            str_ptr2++;
            if (found_paren && count == 0)
                break;
        }
        if (found_paren && count != 0) {
            fprintf(stderr, "ERROR: problem parsing ternary line %s!\n", line);
            controlled_exit(EXIT_FAILURE);
        }
        else_str = inp_fix_ternary_operator_str(copy_substring(str_ptr, str_ptr2), all);
        if (*str_ptr2 != '}')
            end_str = inp_fix_ternary_operator_str(strdup(str_ptr2+1), all);
        else
            end_str = strdup(str_ptr2);
    } else {
        if ((str_ptr2 = strchr(str_ptr, '}')) != NULL) {
            else_str = inp_fix_ternary_operator_str(copy_substring(str_ptr, str_ptr2), all);
            end_str = strdup(str_ptr2);
        } else {
            else_str = inp_fix_ternary_operator_str(strdup(str_ptr), all);
        }
    }

    if (end_str != NULL) {
        if (beg_str != NULL) {
            new_str = TMALLOC(char, strlen(beg_str) + strlen("ternary_fcn") + strlen(conditional) + strlen(if_str) + strlen(else_str) + strlen(end_str) + 5);
            sprintf(new_str, "%sternary_fcn(%s,%s,%s)%s", beg_str, conditional, if_str, else_str, end_str);
        } else {
            new_str = TMALLOC(char, strlen("ternary_fcn") + strlen(conditional) + strlen(if_str) + strlen(else_str) + strlen(end_str) + 5);
            sprintf(new_str, "ternary_fcn(%s,%s,%s)%s", conditional, if_str, else_str, end_str);
        }
    } else {
        if (beg_str != NULL) {
            new_str = TMALLOC(char, strlen(beg_str) + strlen("ternary_fcn") + strlen(conditional) + strlen(if_str) + strlen(else_str) + 5);
            sprintf(new_str, "%sternary_fcn(%s,%s,%s)", beg_str, conditional, if_str, else_str);
        } else {
            new_str = TMALLOC(char, strlen("ternary_fcn") + strlen(conditional) + strlen(if_str) + strlen(else_str) + 5);
            sprintf(new_str, "ternary_fcn(%s,%s,%s)", conditional, if_str, else_str);
        }
    }

    tfree(line);
    tfree(conditional);
    tfree(if_str);
    tfree(else_str);
    if (beg_str != NULL)
        tfree(beg_str);
    if (end_str != NULL)
        tfree(end_str);

    return new_str;
}


static void
inp_fix_ternary_operator(struct line *start_card)
{
    struct line *card;
    char        *line;
    bool found_control = FALSE;

    for (card = start_card; card != NULL; card = card->li_next) {
        line = card->li_line;

        /* exclude replacement of ternary function between .control and .endc */
        if (ciprefix(".control", line))
            found_control = TRUE;
        if (ciprefix(".endc", line))
            found_control = FALSE;
        if (found_control)
            continue;

        /* ternary operator for B source done elsewhere */
        if (*line == 'B' || *line == 'b')
            continue;
        if (*line == '*')
            continue;
        /* .param, .func, and .meas lines handled here (2nd argument FALSE) */
        if (strchr(line, '?') && strchr(line, ':'))
            card->li_line = inp_fix_ternary_operator_str(line, FALSE);
    }
}


/*-------------------------------------------------------------------------*
   removes  " " quotes, returns lower case letters,
   replaces non-printable characterss with '_'                                                                       *
 *-------------------------------------------------------------------------*/
void
inp_casefix(char *string)
{
#ifdef HAVE_CTYPE_H
    if (string)
        while (*string) {
#ifdef HAS_ASCII
            /* ((*string) & 0177): mask off all but the first seven bits, 0177: octal */
            *string = (char) strip(*string);
#endif
            if (*string == '"') {
                *string++ = ' ';
                while (*string && *string != '"')
                    string++;
                if (*string == '\0')
                    continue; /* needed if string is "something ! */
                if (*string == '"')
                    *string = ' ';
            }
            if (!isspace(*string) && !isprint(*string))
                *string = '_';
            if (isupper(*string))
                *string = (char) tolower(*string);
            string++;
        }
#endif
}


/* Strip all end-of-line comments from a deck */
static void
inp_stripcomments_deck(struct line *deck)
{
    struct line *c = deck;
    while (c != NULL) {
        inp_stripcomments_line(c->li_line);
        c = c->li_next;
    }
}


/*
 * SJB 21 April 2005
 * Added support for end-of-line comments that begin with any of the following:
 *   ';'
 *   '$ '
 *   '//' (like in c++ and as per the numparam code)
 *   '--' (as per the numparam code)
 * Any following text to the end of the line is ignored.
 * Note requirement for $ to be followed by a space. This is to avoid conflict
 * with use in front of a variable.
 * Comments on a contunuation line (i.e. line begining with '+') are allowed
 * and are removed before lines are stitched.
 * Lines that contain only an end-of-line comment with or without leading white
 * space are also allowed.

   If there is only white space before the end-of-line comment the
   the whole line is converted to a normal comment line (i.e. one that
   begins with a '*').
   BUG: comment characters in side of string literals are not ignored. */
static void
inp_stripcomments_line(char *s)
{
    char c = ' ';               /* anything other than a comment character */
    char *d = s;
    if (*s == '\0')
        return;                 /* empty line */
    if (*s == '*')
        return;                 /* line is already a comment */
    /* look for comments */
    while ((c = *d) != '\0') {
        d++;
        if (*d == ';') {
            break;
        } else if ((c == '$') && (*d == ' ')) {
            d--;                /* move d back to first comment character */
            break;
        } else if ((*d == c) && ((c == '/') || (c == '-'))) {
            d--;                /* move d back to first comment character */
            break;
        }
    }
    /* d now points to the first comment character of the null at the string end */

    /* check for special case of comment at start of line */
    if (d == s) {
        *s = '*';               /* turn into normal comment */
        return;
    }

    if (d>s) {
        d--;
        /* d now points to character just before comment */

        /* eat white space at end of line */
        while (d >= s) {
            if ((*d != ' ') && (*d != '\t'))
                break;
            d--;
        }
        d++;
        /* d now points to the first white space character before the
           end-of-line or end-of-line comment, or it points to the first
           end-of-line comment character, or to the begining of the line */
    }

    /* Check for special case of comment at start of line
       with or without preceeding white space */
    if (d <= s) {
        *s = '*'; /* turn the whole line into normal comment */
        return;
    }

    *d = '\0'; /* terminate line in new location */
}


static void
inp_change_quotes(char *s)
{
    bool first_quote = FALSE;

    while (*s) {
        if (*s == '\'') {
            if (first_quote == FALSE) {
                *s = '{';
                first_quote = TRUE;
            } else {
                *s = '}';
                first_quote = FALSE;
            }
        }
        s++;
    }
}


static char*
inp_fix_subckt(char *s)
{
    struct line *head = NULL, *newcard = NULL, *start_card = NULL, *end_card = NULL, *prev_card = NULL, *c = NULL;
    char *equal, *beg, *buffer, *ptr1, *ptr2, *str, *new_str = NULL;
    char keep;
    int  num_params = 0, i = 0, bracedepth = 0;
    /* find first '=' */
    equal = strchr(s, '=');
    if (!strstr(s, "params:") && equal != NULL) {
        /* get subckt name (ptr1 will point to name) */
        ptr1 = skip_non_ws(s);
        ptr1 = skip_ws(ptr1);
        for (ptr2 = ptr1; *ptr2 && !isspace(*ptr2) && !isquote(*ptr2); ptr2++)
            ;

        subckt_w_params[num_subckt_w_params++] = copy_substring(ptr1, ptr2);

        /* go to beginning of first parameter word  */
        /* s    will contain only subckt definition */
        /* beg  will point to start of param list   */
        beg = skip_back_ws(equal - 1);
        beg = skip_back_non_ws(beg);
        *beg = '\0';
        beg++;

        head = alloc(struct line);
        /* create list of parameters that need to get sorted */
        while (*beg && (ptr1 = strchr(beg, '=')) != NULL) {
#ifndef NOBRACE
            /* alternative patch to cope with spaces:
               get expression between braces {...} */
            ptr2 = skip_ws(ptr1 + 1);
            ptr1 = skip_back_ws(ptr1 - 1);
            ptr1 = skip_back_non_ws(ptr1) + 1;
            /* ptr1 points to beginning of parameter */

            /* if parameter is an expression and starts with '{', find closing '}'
               Braces maybe nested (will they ever be ?). */
            if (*ptr2 == '{') {
                bracedepth = 1;
                while (bracedepth > 0) {
                    ptr2++;
                    if (*ptr2 == '{')
                        bracedepth++;
                    else if (*ptr2 == '}')
                        bracedepth--;
                    else if (*ptr2 == '\0') {
                        fprintf(stderr, "Error: Missing } in line %s\n", s);
                        controlled_exit(EXIT_FAILURE);
                    }
                }
                ptr2++;/* ptr2 points past end of parameter {...} */
            }
            else
            /* take only the next token (separated by space) as the parameter */
                ptr2 = skip_non_ws(ptr2); /* ptr2 points past end of parameter       */

            keep  = *ptr2;
            if (keep == '\0') {
                /* End of string - don't go over end. */
                beg = ptr2;
            } else {
                *ptr2 = '\0';
                beg   = ptr2+1;
            }

            newcard = alloc(struct line);
            str     = strdup(ptr1);

            newcard->li_line = str;
            newcard->li_next = NULL;

            if (start_card == NULL)
                head->li_next = start_card = newcard;
            else
                prev_card->li_next = newcard;

            prev_card = end_card = newcard;
            num_params++;
        }
#else
            /* patch provided by Ivan Riis Nielsen */
            bool done = FALSE;
            int buf_len = 32, buf_idx = 0;
            char *buf = TMALLOC(char, buf_len), *p1 = beg, *p2 = beg;

            while (!done) {
                while (*p2 && !isspace(*p2)) {

                    if (buf_idx >= buf_len) {
                        buf_len *= 2;
                        buf = TREALLOC(char, buf, buf_len);
                    }
                    buf[buf_idx++] = *(p2++);
                }
                p1 = skip_ws(p2);
                if (*p1 == '\0' || !(strchr("+-*/<>=(!,{", p2[-1]) || strchr("+-*/<>=()!,}", *p1))) {
                    if (buf_idx >= buf_len) {
                        buf_len *= 2;
                        buf = TREALLOC(char, buf, buf_len);
                    }
                    buf[buf_idx++] = '\0';
                    beg = p1;

                    newcard = alloc(struct line);
                    str     = buf;

                    newcard->li_line = str;
                    newcard->li_next = NULL;

                    if (start_card == NULL)
                        head->li_next = start_card = newcard;
                    else
                        prev_card->li_next = newcard;

                    prev_card = end_card = newcard;
                    num_params++;

                    done = TRUE;
                }
                else
                    p2 = p1;
            }
        }
#endif
        /* now sort parameters in order of dependencies */
        inp_sort_params(start_card, end_card, head, start_card, end_card);

        /* create new ordered parameter string for subckt call */
        c = head->li_next;
        tfree(head);
        for (i = 0; i < num_params && c != NULL; i++) {
            if (new_str == NULL) {
                new_str = strdup(c->li_line);
            } else {
                str     = new_str;
                new_str = TMALLOC(char, strlen(str) + strlen(c->li_line) + 2);
                sprintf(new_str, "%s %s", str, c->li_line);
                tfree(str);
            }
            tfree(c->li_line);
            head = c;
            c = c->li_next;
            tfree(head);
        }

        /* create buffer and insert params: */
        buffer = TMALLOC(char, strlen(s) + 9 + strlen(new_str) + 1);
        sprintf(buffer, "%s params: %s", s, new_str);

        tfree(s);
        tfree(new_str);

        s = buffer;
    }
    return s;
}


static char*
inp_remove_ws(char *s)
{
    char *big_buff;
    int  big_buff_index = 0;
    char *buffer, *curr;
    bool is_expression = FALSE;

    big_buff = TMALLOC(char, strlen(s) + 2);
    curr = s;

    while (*curr != '\0') {
        if (*curr == '{')
            is_expression = TRUE;
        if (*curr == '}')
            is_expression = FALSE;

        big_buff[big_buff_index++] = *curr;
        if (*curr == '=' || (is_expression && (is_arith_char(*curr) || *curr == ','))) {
            curr = skip_ws(curr + 1);

            if (*curr == '{')
                is_expression = TRUE;
            if (*curr == '}')
                is_expression = FALSE;

            big_buff[big_buff_index++] = *curr;
        }
        if (*curr != '\0')
            curr++;
        if (isspace(*curr)) {
            curr = skip_ws(curr);
            if (is_expression) {
                if (*curr != '=' && !is_arith_char(*curr) && *curr != ',')
                    big_buff[big_buff_index++] = ' ';
            } else {
                if (*curr != '=')
                    big_buff[big_buff_index++] = ' ';
            }
        }
    }

    big_buff[big_buff_index] = '\0';

    buffer = copy(big_buff);

    tfree(s);
    tfree(big_buff);

    return buffer;
}


/*
  change quotes from '' to {}
  .subckt name 1 2 3 params: l=1 w=2 --> .subckt name 1 2 3 l=1 w=2
   x1 1 2 3 params: l=1 w=2 --> x1 1 2 3 l=1 w=2
   modify .subckt lines by calling inp_fix_subckt()
   No changes to lines in .control section !
*/
static void
inp_fix_for_numparam(struct line *deck)
{
    bool found_control = FALSE;
    struct line *c = deck;
    char *str_ptr;

    while (c != NULL) {
        if (ciprefix("*lib", c->li_line) || ciprefix("*inc", c->li_line)) {
            c = c->li_next;
            continue;
        }

        /* exclude lines between .control and .endc from getting quotes changed */
        if (ciprefix(".control", c->li_line))
            found_control = TRUE;
        if (ciprefix(".endc", c->li_line))
            found_control = FALSE;
        if (found_control) {
            c = c->li_next;
            continue;
        }

        inp_change_quotes(c->li_line);

        if ((inp_compat_mode == COMPATMODE_ALL) || (inp_compat_mode == COMPATMODE_PS))
            if (ciprefix(".subckt", c->li_line) || ciprefix("x", c->li_line)) {
               /* remove params: */
                str_ptr = strstr(c->li_line, "params:");
                if (str_ptr)
                    memcpy(str_ptr, "       ", 7);
            }

        if (ciprefix(".subckt", c->li_line))
            c->li_line = inp_fix_subckt(c->li_line);

        c = c->li_next;
    }
}


static void
inp_remove_excess_ws(struct line *deck)
{
    struct line *c = deck;
    bool found_control = FALSE;
    while (c != NULL) {
        if (*c->li_line == '*') {
            c = c->li_next;
            continue;
        }

        /* exclude echo lines between .control and .endc from removing white spaces */
        if (ciprefix(".control", c->li_line))
            found_control = TRUE;
        if (ciprefix(".endc", c->li_line))
            found_control = FALSE;
        if ((found_control) && (ciprefix("echo", c->li_line))) {
            c = c->li_next;
            continue;
        }
        c->li_line = inp_remove_ws(c->li_line); /* freed in fcn */
        c = c->li_next;
    }
}


static void
inp_determine_libraries(struct line *deck, char *lib_name)
{
    struct line *c = deck;
    char *z, keep_char1, keep_char2;
    int i;
    bool read_line = FALSE;

    if (lib_name == NULL)
        read_line = TRUE;

    while (c != NULL) {
        char *line = c->li_line;

        if (ciprefix(".endl", line) && lib_name != NULL)
            read_line = FALSE;

        if (ciprefix("*lib", line) || ciprefix(".lib", line)) {
            char *s, *t, *y;
            s = skip_non_ws(line);
            while (isspace(*s) || isquote(*s))
                s++;
            for (t = s; *t && !isspace(*t) && !isquote(*t); t++)
                ;
            y = t;
            while (isspace(*y) || isquote(*y))
                y++;

            /* .lib <lib name> */
            if (!*y) {
                keep_char1 = *t;
                *t = '\0';

                if (lib_name != NULL && strcmp(lib_name, s) == 0)
                    read_line = TRUE;
                *t = keep_char1;
            }
            /* .lib <file name> <lib name> */
            else if (read_line == TRUE) {

                char *copys = NULL;

                for (z = y; *z && !isspace(*z) && !isquote(*z); z++)
                    ;
                keep_char1 = *t;
                keep_char2 = *z;
                *t = '\0';
                *z = '\0';

                if (*s == '~') {
                    copys = cp_tildexpand(s);
                    if (copys)
                        s = copys;
                }
                for (i = 0; i < num_libraries; i++)
                    if (cieq(library_file[i], s)) {
                        if (find_lib_name(i, y) < 0) {
                            new_lib_name(i, y, c);
                            /* see if other libraries referenced */
                            inp_determine_libraries(libraries[i], y);
                        }
                    }
                *line = '*';  /* comment out .lib line */
                *t = keep_char1;
                *z = keep_char2;
                /* FIXME, copys not freed ?! */
            }
        }
        c = c->li_next;
    }
}


static void
inp_init_lib_data(void)
{
    int i;

    for (i = 0; i < num_libraries; i++)
        num_lib_names[i] = 0;
}


static char*
inp_get_subckt_name(char *s)
{
    char *subckt_name, *end_ptr = strchr(s, '=');

    if (end_ptr) {
        end_ptr = skip_back_ws(end_ptr - 1);
        end_ptr = skip_back_non_ws(end_ptr) + 1;
    } else {
        end_ptr = s + strlen(s);
    }

    end_ptr = skip_back_ws(end_ptr - 1) + 1;
    subckt_name = skip_back_non_ws(end_ptr - 1) + 1;

    return copy_substring(subckt_name, end_ptr);
}


static int
inp_get_params(char *line, char *param_names[], char *param_values[])
{
    char *equal_ptr = strchr(line, '=');
    char *end, *name, *value;
    int  num_params = 0;
    char tmp_str[1000];
    char keep;
    bool is_expression = FALSE;

    while ((equal_ptr = strchr(line, '=')) != NULL) {

        // check for equality '=='
        if (*(equal_ptr+1) == '=') {
            line = equal_ptr+2;
            continue;
        }
        // check for '!=', '<=', '>='
        if (*(equal_ptr-1) == '!' || *(equal_ptr-1) == '<' || *(equal_ptr-1) == '>') {
            line = equal_ptr+1;
            continue;
        }

        is_expression = FALSE;

        /* get parameter name */
        name = skip_back_ws(equal_ptr - 1);
        end  = name + 1;
        name = skip_back_non_ws(name) + 1;

        param_names[num_params++] = copy_substring(name, end);

        /* get parameter value */
        value = skip_ws(equal_ptr + 1);

        if (*value == '{')
            is_expression = TRUE;
        end = value;
        if (is_expression)
            while (*end && *end != '}')
                end++;
        else
            end = skip_non_ws(end);

        if (is_expression)
            end++;
        keep = *end;
        *end = '\0';

        if (*value != '{' &&
            !(isdigit(*value) || (*value == '.' && isdigit(*(value+1))))) {
            sprintf(tmp_str, "{%s}", value);
            value = tmp_str;
        }
        param_values[num_params-1] = strdup(value);
        *end = keep;

        line = end;
    }

    return num_params;
}


static char*
inp_fix_inst_line(char *inst_line,
                   int num_subckt_params, char *subckt_param_names[], char *subckt_param_values[],
                   int num_inst_params, char *inst_param_names[], char *inst_param_values[])
{
    char *end, *inst_name, *inst_name_end;
    char *curr_line = inst_line, *new_line = NULL;
    int i, j;

    inst_name_end = skip_non_ws(inst_line);
    inst_name = copy_substring(inst_line, inst_name_end);

    end = strchr(inst_line, '=');
    if (end) {
        end = skip_back_ws(end - 1);
        end = skip_back_non_ws(end);
        *end = '\0';
    }

    for (i = 0; i < num_subckt_params; i++)
        for (j = 0; j < num_inst_params; j++)
            if (strcmp(subckt_param_names[i], inst_param_names[j]) == 0) {
                tfree(subckt_param_values[i]);
                subckt_param_values[i] = strdup(inst_param_values[j]);
            }

    for (i = 0; i < num_subckt_params; i++) {
        new_line = TMALLOC(char, strlen(curr_line) + strlen(subckt_param_values[i]) + 2);
        sprintf(new_line, "%s %s", curr_line, subckt_param_values[i]);

        tfree(curr_line);
        tfree(subckt_param_names[i]);
        tfree(subckt_param_values[i]);

        curr_line = new_line;
    }

    for (i = 0; i < num_inst_params; i++) {
        tfree(inst_param_names[i]);
        tfree(inst_param_values[i]);
    }

    tfree(inst_name);

    return curr_line;
}


static bool
found_mult_param(int num_params, char *param_names[])
{
    bool found_mult = FALSE;
    int i;

    for (i = 0; i < num_params; i++)
        if (strcmp(param_names[i], "m") == 0)
            found_mult = TRUE;

    return found_mult;
}


static int
inp_fix_subckt_multiplier(struct line *subckt_card,
                          int num_subckt_params, char *subckt_param_names[], char *subckt_param_values[])
{
    struct line *card;
    char *new_str;

    subckt_param_names[num_subckt_params]  = strdup("m");
    subckt_param_values[num_subckt_params] = strdup("1");
    num_subckt_params = num_subckt_params + 1;

    if (!strstr(subckt_card->li_line, "params:")) {
        new_str = TMALLOC(char, strlen(subckt_card->li_line) + 13);
        sprintf(new_str, "%s params: m=1", subckt_card->li_line);
        subckt_w_params[num_subckt_w_params++] = get_subckt_model_name(subckt_card->li_line);
    } else {
        new_str = TMALLOC(char, strlen(subckt_card->li_line) + 5);
        sprintf(new_str, "%s m=1", subckt_card->li_line);
    }

    tfree(subckt_card->li_line);
    subckt_card->li_line = new_str;

    for (card = subckt_card->li_next;
         card != NULL && !ciprefix(".ends", card->li_line);
         card = card->li_next) {
        /* no 'm' for B line or comment line */
        if ((*(card->li_line) == '*') || (*(card->li_line) == 'b'))
            continue;
        /* no 'm' for model cards */
        if (ciprefix(".model", card->li_line))
            continue;
        new_str = TMALLOC(char, strlen(card->li_line) + 7);
        sprintf(new_str, "%s m={m}", card->li_line);

        tfree(card->li_line);
        card->li_line = new_str;
    }

    return num_subckt_params;
}


static void
inp_fix_inst_calls_for_numparam(struct line *deck)
{
    struct line *c = deck;
    struct line *d, *p = NULL;
    char *inst_line;
    char *subckt_line;
    char *subckt_name;
    char *subckt_param_names[1000];
    char *subckt_param_values[1000];
    char *inst_param_names[1000];
    char *inst_param_values[1000];
    char name_w_space[1000];
    int  num_subckt_params = 0;
    int  num_inst_params   = 0;
    int  i, j, k;
    bool flag = FALSE;
    bool found_subckt = FALSE;
    bool found_param_match = FALSE;

    // first iterate through instances and find occurences where 'm' multiplier needs to be
    // added to the subcircuit -- subsequent instances will then need this parameter as well
    for (c = deck; c != NULL; c = c->li_next) {
        inst_line = c->li_line;

        if (*inst_line == '*')
            continue;

        if (ciprefix("x", inst_line)) {
            num_inst_params = inp_get_params(inst_line, inst_param_names, inst_param_values);
            subckt_name     = inp_get_subckt_name(inst_line);

            if (found_mult_param(num_inst_params, inst_param_names)) {
                flag = FALSE;
                // iterate through the deck to find the subckt (last one defined wins)
                d = deck;
                while (d != NULL) {
                    subckt_line = d->li_line;
                    if (ciprefix(".subckt", subckt_line)) {
                        subckt_line = skip_non_ws(subckt_line);
                        subckt_line = skip_ws(subckt_line);

                        sprintf(name_w_space, "%s ", subckt_name);
                        if (strncmp(subckt_line, name_w_space, strlen(name_w_space)) == 0) {
                            p = d;
                            flag = TRUE;
                        }
                    }
                    d = d->li_next;
                }
                if (flag) {
                    num_subckt_params = inp_get_params(p->li_line, subckt_param_names, subckt_param_values);

                    if (num_subckt_params == 0 || !found_mult_param(num_subckt_params, subckt_param_names))
                        inp_fix_subckt_multiplier(p, num_subckt_params, subckt_param_names, subckt_param_values);
                }
            }
            tfree(subckt_name);
            if (flag)
                for (i = 0; i < num_subckt_params; i++) {
                    tfree(subckt_param_names[i]);
                    tfree(subckt_param_values[i]);
                }
            for (i = 0; i < num_inst_params; i++) {
                tfree(inst_param_names[i]);
                tfree(inst_param_values[i]);
            }
        }
    }

    c = deck;
    while (c != NULL) {
        inst_line = c->li_line;

        if (*inst_line == '*') {
            c = c->li_next;
            continue;
        }
        if (ciprefix("x", inst_line)) {
            subckt_name = inp_get_subckt_name(inst_line);
            for (i = 0; i < num_subckt_w_params; i++) {
                if (strcmp(subckt_w_params[i], subckt_name) == 0) {
                    sprintf(name_w_space, "%s ", subckt_name);

                    /* find .subckt line */
                    found_subckt = FALSE;

                    d = deck;
                    while (d != NULL) {
                        subckt_line = d->li_line;
                        if (ciprefix(".subckt", subckt_line)) {
                            subckt_line = skip_non_ws(subckt_line);
                            subckt_line = skip_ws(subckt_line);

                            if (strncmp(subckt_line, name_w_space, strlen(name_w_space)) == 0) {
                                num_subckt_params = inp_get_params(subckt_line, subckt_param_names, subckt_param_values);
                                num_inst_params   = inp_get_params(inst_line, inst_param_names, inst_param_values);

                                // make sure that if have inst params that one matches subckt
                                found_param_match = FALSE;
                                if (num_inst_params == 0) {
                                    found_param_match = TRUE;
                                } else {
                                    for (j = 0; j < num_inst_params; j++) {
                                        for (k = 0; k < num_subckt_params; k++)
                                            if (strcmp(subckt_param_names[k], inst_param_names[j]) == 0) {
                                                found_param_match = TRUE;
                                                break;
                                            }
                                        if (found_param_match)
                                            break;
                                    }
                                }

                                if (!found_param_match) {
                                    // comment out .subckt and continue
                                    while (d != NULL && !ciprefix(".ends", d->li_line)) {
                                        *(d->li_line) = '*';
                                        d = d->li_next;
                                    }
                                    *(d->li_line) = '*';
                                    d = d->li_next;
                                    continue;
                                }

                                c->li_line = inp_fix_inst_line(inst_line, num_subckt_params, subckt_param_names, subckt_param_values, num_inst_params, inst_param_names, inst_param_values);
                                found_subckt = TRUE;
                                for (i = 0; i < num_subckt_params; i++) {
                                    tfree(subckt_param_names[i]);
                                    tfree(subckt_param_values[i]);
                                }
                                for (i = 0; i < num_inst_params; i++) {
                                    tfree(inst_param_names[i]);
                                    tfree(inst_param_values[i]);
                                }
                            }
                        }
                        if (found_subckt)
                            break;
                        d = d->li_next;
                    }
                    break;
                }
            }
            tfree(subckt_name);
        }
        c = c->li_next;
    }
}


static void
inp_get_func_from_line(char *line)
{
    char *ptr, *end;
    char keep;
    char temp_buf[5000];
    int  num_params = 0;
    int  str_len = 0;
    int  i = 0;

    /* get function name */
    line = skip_non_ws(line);
    line = skip_ws(line);
    end = line;
    while (!isspace(*end) && *end != '(')
        end++;
    keep = *end;
    *end = '\0';

    /* see if already encountered this function */
    for (i = 0; i < num_functions; i++)
        if (strcmp(func_names[i], line) == 0)
            break;

    func_names[num_functions++] = strdup(line);
    *end = keep;

    num_params = 0;

    /* get function parameters */
    while (*end != '(')
        end++;
    while (*end != ')') {
        end = skip_ws(end + 1);
        ptr = end;
        while (!isspace(*end) && *end != ',' && *end != ')')
            end++;
        if (end > ptr)
            func_params[num_functions-1][num_params++] = copy_substring(ptr, end);
    }
    num_parameters[num_functions-1] = num_params;

    /* get function macro */
    str_len = 0;
    while (*end != '{')
        end++;
    end++;
    while (*end  &&  *end != '}') {
        if (!isspace(*end))
            temp_buf[str_len++] = *end;
        end++;
    }
    temp_buf[str_len++] = '\0';

    func_macro[num_functions-1] = strdup(temp_buf);
}


//
// only grab global functions; skip subckt functions
//
static void
inp_grab_func(struct line *deck)
{
    struct line *c        = deck;
    bool        is_subckt = FALSE;

    while (c != NULL) {
        if (*c->li_line == '*') {
            c = c->li_next;
            continue;
        }
        if (ciprefix(".subckt", c->li_line))
            is_subckt = TRUE;
        if (ciprefix(".ends", c->li_line))
            is_subckt = FALSE;

        if (!is_subckt && ciprefix(".func", c->li_line)) {
            inp_get_func_from_line(c->li_line);
            *c->li_line = '*';
        }
        c = c->li_next;
    }
}


static void
inp_grab_subckt_func(struct line *subckt)
{
    struct line *c = subckt;

    while (!ciprefix(".ends", c->li_line)) {
        if (ciprefix(".func", c->li_line)) {
            inp_get_func_from_line(c->li_line);
            *c->li_line = '*';
        }
        c = c->li_next;
    }
}


static char*
inp_do_macro_param_replace(int fcn_number, char *params[])
{
    char *param_ptr, *curr_ptr, *new_str, *curr_str = NULL, *search_ptr;
    char keep, before, after;
    int  i;

    if (num_parameters[fcn_number] == 0)
        return strdup(func_macro[fcn_number]);

    for (i = 0; i < num_parameters[fcn_number]; i++) {
        if (curr_str == NULL) {
            search_ptr = curr_ptr = func_macro[fcn_number];
        } else {
            search_ptr = curr_ptr = curr_str;
            curr_str = NULL;
        }
        while ((param_ptr = strstr(search_ptr, func_params[fcn_number][i])) != NULL) {

            /* make sure actually have the parameter name */
            if (param_ptr == search_ptr) /* no valid 'before' */
                before = '\0';
            else
                before = *(param_ptr-1);
            after  = *(param_ptr+strlen(func_params[fcn_number][i]));
            if (!(is_arith_char(before) || isspace(before) ||
                  before == ',' || before == '=' || (param_ptr-1) < curr_ptr) ||
                !(is_arith_char(after) || isspace(after) ||
                  after == ',' || after  == '=' || after  == '\0'))
            {
                search_ptr = param_ptr + 1;
                continue;
            }

            keep       = *param_ptr;
            *param_ptr = '\0';

            {
                size_t curr_str_len = curr_str ? strlen(curr_str) : 0;
                size_t len = strlen(curr_ptr) + strlen(params[i]) + 1;
                if (str_has_arith_char(params[i])) {
                    curr_str = TREALLOC(char, curr_str, curr_str_len + len + 2);
                    sprintf(curr_str + curr_str_len, "%s(%s)", curr_ptr, params[i]);
                } else {
                    curr_str = TREALLOC(char, curr_str, curr_str_len + len);
                    sprintf(curr_str + curr_str_len, "%s%s", curr_ptr, params[i]);
                }
            }

            *param_ptr = keep;
            search_ptr = curr_ptr = param_ptr + strlen(func_params[fcn_number][i]);
        }
        if (param_ptr == NULL) {
            if (curr_str == NULL) {
                curr_str = curr_ptr;
            } else {
                new_str = TMALLOC(char, strlen(curr_str) + strlen(curr_ptr) + 1);
                sprintf(new_str, "%s%s", curr_str, curr_ptr);
                tfree(curr_str);
                curr_str = new_str;
            }
        }
    }
    return curr_str;
}


static char*
inp_expand_macro_in_str(char *str)
{
    int  i;
    char *c;
    char *open_paren_ptr, *close_paren_ptr, *fcn_name, *params[1000];
    char *curr_ptr, *macro_str, *curr_str = NULL;
    int  num_parens, num_params;
    char *orig_ptr = str, *search_ptr = str, *orig_str = strdup(str);
    char keep;

    // printf("%s: enter(\"%s\")\n", __FUNCTION__, str);
    while ((open_paren_ptr = strchr(search_ptr, '(')) != NULL) {

        fcn_name = open_paren_ptr;
        while (--fcn_name >= search_ptr)
            if (!isalnum(*fcn_name) && *fcn_name != '_')
                break;
        fcn_name++;

        search_ptr = open_paren_ptr + 1;
        if (open_paren_ptr == fcn_name)
            continue;

        *open_paren_ptr = '\0';

        for (i = 0; i < num_functions; i++)
            if (strcmp(func_names[i], fcn_name) == 0)
                break;

        *open_paren_ptr = '(';

        if (i >= num_functions)
            continue;

        /* find the closing paren */
        num_parens = 1;
        for (c = open_paren_ptr + 1; *c; c++) {
            if (*c == '(')
                num_parens++;
            if (*c == ')' && --num_parens == 0)
                break;
        }

        if (num_parens) {
            fprintf(stderr, "ERROR: did not find closing parenthesis for function call in str: %s\n", orig_str);
            controlled_exit(EXIT_FAILURE);
        }

        close_paren_ptr = c;

        /*
         * if (ciprefix("v(", curr_ptr)) {
         *     // look for any commas and change to ' '
         *     char *str_ptr = curr_ptr;
         *     while (*str_ptr != '\0' && *str_ptr != ')') {
         *         if (*str_ptr == ',' || *str_ptr == '(')
         *             *str_ptr = ' '; str_ptr++; }
         *     if (*str_ptr == ')')
         *         *str_ptr = ' ';
         * }
         */

        /* get the parameters */
        curr_ptr = open_paren_ptr + 1;

        for (num_params = 0; curr_ptr < close_paren_ptr; curr_ptr++) {
            char *beg_parameter;
            int num_parens;
            if (isspace(*curr_ptr))
                continue;
            beg_parameter = curr_ptr;
            num_parens = 0;
            for (; curr_ptr < close_paren_ptr; curr_ptr++) {
                if (*curr_ptr == '(')
                    num_parens++;
                if (*curr_ptr == ')')
                    num_parens--;
                if (*curr_ptr == ',' && num_parens == 0)
                    break;
            }
            params[num_params++] =
                inp_expand_macro_in_str(copy_substring(beg_parameter, curr_ptr));
        }

        if (num_parameters[i] != num_params) {
            fprintf(stderr, "ERROR: parameter mismatch for function call in str: %s\n", orig_str);
            controlled_exit(EXIT_FAILURE);
        }

        macro_str = inp_do_macro_param_replace(i, params);
        keep  = *fcn_name;
        *fcn_name = '\0';
        {
            size_t curr_str_len = curr_str ? strlen(curr_str) : 0;
            size_t len = strlen(str) + strlen(macro_str) + 3;
            curr_str = TREALLOC(char, curr_str, curr_str_len + len);
            sprintf(curr_str + curr_str_len, "%s(%s)", str, macro_str);
        }
        *fcn_name = keep;

        search_ptr = str = close_paren_ptr + 1;
    }

    if (curr_str == NULL) {
        curr_str = orig_ptr;
    } else {
        if (str != NULL) {
            size_t curr_str_len = strlen(curr_str);
            size_t len = strlen(str) + 1;
            curr_str = TREALLOC(char, curr_str, curr_str_len + len);
            sprintf(curr_str + curr_str_len, "%s", str);
        }
        tfree(orig_ptr);
    }

    tfree(orig_str);
    // printf("%s: --> \"%s\"\n", __FUNCTION__, curr_str);

    return curr_str;
}


static void
inp_expand_macros_in_func(void)
{
    int  i;

    for (i = 0; i < num_functions; i++)
        func_macro[i] = inp_expand_macro_in_str(func_macro[i]);
}


static void
inp_expand_macros_in_deck(struct line *deck)
{
    struct line *c = deck;
    int         prev_num_functions = 0, i, j;

    while (c != NULL) {
        if (*c->li_line == '*') {
            c = c->li_next;
            continue;
        }
        if (ciprefix(".subckt", c->li_line)) {
            prev_num_functions = num_functions;
            inp_grab_subckt_func(c);
            if (prev_num_functions != num_functions)
                inp_expand_macros_in_func();
        }
        if (ciprefix(".ends", c->li_line)) {
            if (prev_num_functions != num_functions) {
                for (i = prev_num_functions; i < num_functions; i++) {
                    tfree(func_names[i]);
                    tfree(func_macro[i]);
                    for (j = 0; j < num_parameters[i]; j++)
                        tfree(func_params[i][j]);
                    num_functions = prev_num_functions;
                }
            }
        }

        if (*c->li_line != '*')
            c->li_line = inp_expand_macro_in_str(c->li_line);
        c = c->li_next;
    }
}


/* Put {} around tokens for handling in numparam.
   Searches for the next '=' in the line to become active.
   Several exceptions (eg. no 'set' or 'b' lines, no .cmodel lines,
   no lines between .control and .endc, no .option lines).
   Special handling of vectors with [] and complex values with < >

   h_vogt 20 April 2008
 * For xspice and num_pram compatibility .cmodel added
 * .cmodel will be replaced by .model in inp_fix_param_values()
 * and then the entire line is skipped (will not be changed by this function).
 * Usage of numparam requires {} around the parameters in the .cmodel line.
 * May be obsolete?
 */
static void
inp_fix_param_values(struct line *deck)
{
    struct line *c = deck;
    char *line, *beg_of_str, *end_of_str, *old_str, *equal_ptr, *new_str;
    char *vec_str, *tmp_str, *natok, *buffer, *newvec, *whereisgt;
    bool control_section = FALSE;
    wordlist *nwl;
    int parens;

    while (c != NULL) {
        line = c->li_line;

        if (*line == '*' || (ciprefix(".param", line) && strchr(line, '{'))) {
            c = c->li_next;
            continue;
        }

        if (ciprefix(".control", line)) {
            control_section = TRUE;
            c = c->li_next;
            continue;
        }
        if (ciprefix(".endc", line)) {
            control_section = FALSE;
            c = c->li_next;
            continue;
        }
        if (control_section || ciprefix(".option", line)) {
            c = c->li_next;    /* no handling of params in "option" lines */
            continue;
        }
        if (ciprefix("set", line)) {
            c = c->li_next;    /* no handling of params in "set" lines */
            continue;
        }
        if (*line == 'b') {
            c = c->li_next;    /* no handling of params in B source lines */
            continue;
        }

        /* for xspice .cmodel: replace .cmodel with .model and skip entire line) */
        if (ciprefix(".cmodel", line)) {
            *(++line) = 'm';
            *(++line) = 'o';
            *(++line) = 'd';
            *(++line) = 'e';
            *(++line) = 'l';
            *(++line) = ' ';
            c = c->li_next;
            continue;
        }

        /* exclude CIDER models */
        if (ciprefix(".model", line) && (strstr(line, "numos") || strstr(line, "numd") ||
                                         strstr(line, "nbjt") || strstr(line, "nbjt2") ||
                                         strstr(line, "numd2"))) {
            c = c->li_next;
            continue;
        }
        /* exclude CIDER devices with ic.file parameter */
        if (strstr(line, "ic.file")) {
            c = c->li_next;
            continue;
        }

        while ((equal_ptr = strchr(line, '=')) != NULL) {

            // special case: .MEASURE {DC|AC|TRAN} result FIND out_variable WHEN out_variable2=out_variable3
            // no braces around out_variable3. out_variable3 may be v(...) or i(...)
            if (ciprefix(".meas", line))
                if (((*(equal_ptr+1) == 'v') || (*(equal_ptr+1) == 'i')) &&
                    (*(equal_ptr+2) == '('))
                {
                    // find closing ')' and skip token v(...) or i(...)
                    while (*equal_ptr != ')' && *(equal_ptr+1) != '\0')
                        equal_ptr++;
                    line = equal_ptr + 1;
                    continue;
                }

            // skip over equality '=='
            if (*(equal_ptr+1) == '=') {
                line += 2;
                continue;
            }
            // check for '!=', '<=', '>='
            if (*(equal_ptr-1) == '!' || *(equal_ptr-1) == '<' || *(equal_ptr-1) == '>')
            {
                line += 1;
                continue;
            }

            beg_of_str = skip_ws(equal_ptr + 1);
            /* all cases where no {} have to be put around selected token */
            if (isdigit(*beg_of_str) ||
                *beg_of_str == '{' ||
                *beg_of_str == '.' ||
                *beg_of_str == '"' ||
                (*beg_of_str == '-' && isdigit(*(beg_of_str+1))) ||
                (*beg_of_str == '-' && *(beg_of_str+1) == '.' && isdigit(*(beg_of_str+2))) ||
                ciprefix("true", beg_of_str) ||
                ciprefix("false", beg_of_str))
            {
                line = equal_ptr + 1;
            } else if (*beg_of_str == '[') {
                /* A vector following the '=' token: code to put curly brackets around all params
                inside a pair of square brackets */
                end_of_str = beg_of_str;
                while (*end_of_str != ']')
                    end_of_str++;
                /* string xx yyy from vector [xx yyy] */
                tmp_str = vec_str = copy_substring(beg_of_str + 1, end_of_str);

                /* work on vector elements inside [] */
                nwl = NULL;
                for (;;) {
                    natok = gettok(&vec_str);
                    if (!natok)
                        break;

                    buffer = TMALLOC(char, strlen(natok) + 4);
                    if (isdigit(*natok) || *natok == '{' || *natok == '.' ||
                        *natok == '"' || (*natok == '-' && isdigit(*(natok+1))) ||
                        ciprefix("true", natok) || ciprefix("false", natok) ||
                        eq(natok, "<") || eq(natok, ">"))
                    {
                        (void) sprintf(buffer, "%s", natok);
                        /* A complex value found inside a vector [< x1 y1> <x2 y2>] */
                        /* < xx and yy > have been dealt with before */
                        /* <xx */
                    } else if (*natok == '<') {
                        if (isdigit(*(natok + 1)) ||
                            (*(natok + 1) == '-' && isdigit(*(natok + 2))))
                        {
                            (void) sprintf(buffer, "%s", natok);
                        } else {
                            *natok = '{';
                            (void) sprintf(buffer, "<%s}", natok);
                        }
                        /* yy> */
                    } else if (strchr(natok, '>')) {
                        if (isdigit(*natok) || (*natok == '-' && isdigit(*(natok+1)))) {
                            (void) sprintf(buffer, "%s", natok);
                        } else {
                            whereisgt = strchr(natok, '>');
                            *whereisgt = '}';
                            (void) sprintf(buffer, "{%s>", natok);
                        }
                        /* all other tokens */
                    } else {
                        (void) sprintf(buffer, "{%s}", natok);
                    }
                    tfree(natok);
                    nwl = wl_cons(copy(buffer), nwl);
                    tfree(buffer);
                }
                tfree(tmp_str);
                nwl = wl_reverse(nwl);
                /* new vector elements */
                newvec = wl_flatten(nwl);
                wl_free(nwl);
                /* insert new vector into actual line */
                *equal_ptr  = '\0';
                new_str = TMALLOC(char, strlen(c->li_line) + strlen(newvec) + strlen(end_of_str + 1) + 5);
                sprintf(new_str, "%s=[%s] %s", c->li_line, newvec, end_of_str+1);
                tfree(newvec);

                old_str    = c->li_line;
                c->li_line = new_str;
                line = new_str + strlen(old_str) + 1;
                tfree(old_str);
            } else if (*beg_of_str == '<') {
                /* A complex value following the '=' token: code to put curly brackets around all params
                inside a pair < > */
                end_of_str = beg_of_str;
                while (*end_of_str != '>')
                    end_of_str++;
                /* string xx yyy from vector [xx yyy] */
                vec_str = copy_substring(beg_of_str + 1, end_of_str);

                /* work on tokens inside <> */
                nwl = NULL;
                for (;;) {
                    natok = gettok(&vec_str);
                    if (!natok)
                        break;

                    buffer = TMALLOC(char, strlen(natok) + 4);
                    if (isdigit(*natok) || *natok == '{' || *natok == '.' ||
                        *natok == '"' || (*natok == '-' && isdigit(*(natok+1))) ||
                        ciprefix("true", natok) || ciprefix("false", natok))
                    {
                        (void) sprintf(buffer, "%s", natok);
                    } else {
                        (void) sprintf(buffer, "{%s}", natok);
                    }
                    tfree(natok);
                    nwl = wl_cons(copy(buffer), nwl);
                    tfree(buffer);
                }
                nwl = wl_reverse(nwl);
                /* new elements of complex variable */
                newvec = wl_flatten(nwl);
                wl_free(nwl);
                /* insert new complex value into actual line */
                *equal_ptr  = '\0';
                new_str = TMALLOC(char, strlen(c->li_line) + strlen(newvec) + strlen(end_of_str + 1) + 5);
                sprintf(new_str, "%s=<%s> %s", c->li_line, newvec, end_of_str+1);
                tfree(newvec);

                old_str    = c->li_line;
                c->li_line = new_str;
                line = new_str + strlen(old_str) + 1;
                tfree(old_str);
            } else {
                /* put {} around token to be accepted as numparam */
                end_of_str = beg_of_str;
                parens = 0;
                while (*end_of_str != '\0' &&
                       (!isspace(*end_of_str) || (parens > 0)))
                {
                    if (*end_of_str == '(')
                        parens++;
                    if (*end_of_str == ')')
                        parens--;
                    end_of_str++;
                }

                *equal_ptr  = '\0';

                if (*end_of_str == '\0') {
                    new_str = TMALLOC(char, strlen(c->li_line) + strlen(beg_of_str) + 4);
                    sprintf(new_str, "%s={%s}", c->li_line, beg_of_str);

                } else {
                    *end_of_str = '\0';
                    new_str = TMALLOC(char, strlen(c->li_line) + strlen(beg_of_str) + strlen(end_of_str + 1) + 5);
                    sprintf(new_str, "%s={%s} %s", c->li_line, beg_of_str, end_of_str+1);
                }
                old_str    = c->li_line;
                c->li_line = new_str;

                line = new_str + strlen(old_str) + 1;
                tfree(old_str);
            }
        }
        c = c->li_next;
    }
}


static char*
get_param_name(char *line)
{
    char *name = NULL, *equal_ptr, *beg;
    char keep;

    if ((equal_ptr = strchr(line, '=')) != NULL) {
        equal_ptr = skip_back_ws(equal_ptr - 1) + 1;

        beg = equal_ptr - 1;
        while (!isspace(*beg) && beg != line)
            beg--;
        if (beg != line)
            beg++;
        keep       = *equal_ptr;
        *equal_ptr = '\0';
        name       = strdup(beg);
        *equal_ptr = keep;
    } else {
        fprintf(stderr, "ERROR: could not find '=' on parameter line '%s'!\n", line);
        controlled_exit(EXIT_FAILURE);
    }
    return name;
}


static char*
get_param_str(char *line)
{
    char *equal_ptr = strchr(line, '=');

    if (equal_ptr)
        return skip_ws(equal_ptr + 1);
    else
        return line;
}


static int
//inp_get_param_level(int param_num, char *depends_on[12000][100], char *param_names[12000], char *param_strs[12000], int total_params, int *level)
inp_get_param_level(int param_num, char ***depends_on, char **param_names, char **param_strs, int total_params, int *level)
{
    int index1 = 0, comp_level = 0, temp_level = 0;
    int index2 = 0;

    if (*(level+param_num) != -1)
        return *(level+param_num);

    while (depends_on[param_num][index1] != NULL) {
        index2 = 0;
        while (index2 <= total_params &&
               param_names[index2] != depends_on[param_num][index1])
            index2++;

        if (index2 > total_params) {
            fprintf(stderr, "ERROR: unable to find dependency parameter for %s!\n", param_names[param_num]);
            controlled_exit(EXIT_FAILURE);
        }
        temp_level = inp_get_param_level(index2, depends_on, param_names, param_strs, total_params, level);
        temp_level++;

        if (temp_level > comp_level)
            comp_level = temp_level;
        index1++;
    }
    *(level+param_num) = comp_level;
    return comp_level;
}


static int
get_number_terminals(char *c)
{
    int i, j, k;
    char *name[12];
    char nam_buf[33];
    bool area_found = FALSE;

    switch (*c) {
    case 'r':
    case 'c':
    case 'l':
    case 'k':
    case 'f':
    case 'h':
    case 'b':
    case 'v':
    case 'i':
    case 'd':
        return 2;
        break;
    case 'u':
    case 'j':
    case 'w':
    case 'z':
        return 3;
        break;
    case 't':
    case 'o':
    case 'g':
    case 'e':
    case 's':
    case 'y':
        return 4;
        break;
    case 'm': /* recognition of 4, 5, 6, or 7 nodes for SOI devices needed */
        i = 0;
        /* find the first token with "off" or "=" in the line*/
        while ((i < 20) && (*c != '\0')) {
            char *inst = gettok_instance(&c);
            strncpy(nam_buf, inst, 32);
            txfree(inst);
            if (strstr(nam_buf, "off") || strchr(nam_buf, '='))
                break;
            i++;
        }
        return i-2;
        break;
    case 'p': /* recognition of up to 100 cpl nodes */
        i = j = 0;
        /* find the last token in the line*/
        while ((i < 100) && (*c != '\0')) {
            char *tmp_inst = gettok_instance(&c);
            strncpy(nam_buf, tmp_inst, 32);
            tfree(tmp_inst);
            if (strchr(nam_buf, '='))
                j++;
            i++;
        }
        if (i == 100)
            return 0;
        return i - j - 2;
        break;
    case 'q': /* recognition of 3/4 terminal bjt's needed */
        /* QXXXXXXX NC NB NE <NS> MNAME <AREA> <OFF> <IC=VBE, VCE> <TEMP=T> */
        /* 12 tokens maximum */
        i = j = 0;
        while ((i < 12) && (*c != '\0')) {
            char *comma;
            name[i] = gettok_instance(&c);
            if (strstr(name[i], "off") || strchr(name[i], '='))
                j++;
            /* If we have IC=VBE, VCE instead of IC=VBE,VCE we need to inc j */
            if ((comma = strchr(name[i], ',')) != NULL && (*(++comma) == '\0'))
                j++;
            /* If we have IC=VBE , VCE ("," is a token) we need to inc j */
            if (eq(name[i], ","))
                j++;
            i++;
        }
        i--;
        area_found = FALSE;
        for (k = i; k > i-j-1; k--) {
            bool only_digits = TRUE;
            char *nametmp = name[k];
            /* MNAME has to contain at least one alpha character. AREA may be assumed
            if we have a token with only digits, and where the previous token does not
            end with a ',' */
            while (*nametmp) {
                if (isalpha(*nametmp) || (*nametmp == ','))
                    only_digits = FALSE;
                nametmp++;
            }
            if (only_digits && (strchr(name[k-1], ',') == NULL))
                area_found = TRUE;
        }
        for (k = i; k >= 0; k--)
            tfree(name[k]);
        if (area_found) {
            return i-j-2;
        } else {
            return i-j-1;
        }
        break;
    default:
        return 0;
        break;
    }
}


/* sort parameters based on parameter dependencies */
static void
inp_sort_params(struct line *start_card, struct line *end_card, struct line *card_bf_start, struct line *s_c, struct line *e_c)
{
    char *param_name = NULL, *param_str = NULL, *param_ptr = NULL;
    int  i, j, num_params = 0, ind = 0, max_level = 0, num_terminals = 0, ioff = 1;
    bool in_control = FALSE;

    bool found_in_list = FALSE;

    struct line *ptr;
    char *curr_line;
    char *str_ptr, *beg, *end, *new_str;
    int  skipped = 0;
    int arr_size = 12000;

    int *level;
    int *param_skip;
    char **param_names;
    char **param_strs;
    char ***depends_on;
    struct line **ptr_array;
    struct line **ptr_array_ordered;

    NG_IGNORE(end_card);

    if (start_card == NULL)
        return;

    /* determine the number of lines with .param */
    ptr = start_card;
    while (ptr != NULL) {
        if (strchr(ptr->li_line, '='))
            num_params++;
        ptr = ptr->li_next;
    }

    arr_size = num_params;
    num_params = 0; /* This is just to keep the code in row 2907ff. */

    // dynamic memory allocation
    level = TMALLOC(int, arr_size);
    param_skip = TMALLOC(int, arr_size);
    param_names = TMALLOC(char *, arr_size);
    param_strs = TMALLOC(char *, arr_size);

    /* array[row][column] -> depends_on[array_size][100] */
    /* rows */
    depends_on = TMALLOC(char **, arr_size);
    /* columns */
    for (i = 0; i < arr_size; i++)
        depends_on[i] = TMALLOC(char *, 100);

    ptr_array = TMALLOC(struct line *, arr_size);
    ptr_array_ordered = TMALLOC(struct line *, arr_size);

    ptr = start_card;
    while (ptr != NULL) {
        // ignore .param lines without '='
        if (strchr(ptr->li_line, '=')) {
            depends_on[num_params][0] = NULL;
            level[num_params]         = -1;
            param_names[num_params]   = get_param_name(ptr->li_line); /* strdup in fcn */
            param_strs[num_params]    = strdup(get_param_str(ptr->li_line));

            ptr_array[num_params++]   = ptr;
        }
        ptr = ptr->li_next;
    }
    // look for duplicately defined parameters and mark earlier one to skip
    // param list is ordered as defined in netlist
    for (i = 0; i < num_params; i++)
        param_skip[i] = 0;

    for (i = 0; i < num_params; i++)
        for (j = num_params-1; j >= 0; j--)
            if (i != j && i < j && strcmp(param_names[i], param_names[j]) == 0) {
                // skip earlier one in list
                param_skip[i] = 1;
                skipped++;
            }

    for (i = 0; i < num_params; i++) {
        if (param_skip[i] == 1)
            continue;

        param_name = param_names[i];
        for (j = 0; j < num_params; j++) {
            if (j == i)
                continue;

            param_str = param_strs[j];

            while ((param_ptr = strstr(param_str, param_name)) != NULL) {
                ioff = (strchr(param_ptr, '}') != NULL ? 1 : 0);  /* want prevent wrong memory access below */
                /* looking for curly braces or other string limiter */
                if ((!isalnum(*(param_ptr-ioff)) && *(param_ptr-ioff) != '_' &&
                     !isalnum(*(param_ptr+strlen(param_name))) &&
                     *(param_ptr+strlen(param_name)) != '_') ||
                    strcmp(param_ptr, param_name) == 0)
                { /* this are cases without curly braces */
                    ind = 0;
                    found_in_list = FALSE;
                    while (depends_on[j][ind] != NULL) {
                        if (strcmp(param_name, depends_on[j][ind]) == 0) {
                            found_in_list = TRUE;
                            break;
                        }
                        ind++;
                    }
                    if (!found_in_list) {
                        depends_on[j][ind++] = param_name;
                        depends_on[j][ind]   = NULL;
                    }
                    break;
                }
                param_str = param_ptr + strlen(param_name);
            }
        }
    }

    for (i = 0; i < num_params; i++) {
        level[i] = inp_get_param_level(i, depends_on, param_names, param_strs, num_params, level);
        if (level[i] > max_level)
            max_level = level[i];
    }

    /* look for unquoted parameters and quote them */
    ptr = s_c;
    in_control = FALSE;
    while (ptr != NULL && ptr != e_c) {
        curr_line = ptr->li_line;

        if (ciprefix(".control", curr_line)) {
            in_control = TRUE;
            ptr = ptr->li_next;
            continue;
        }
        if (ciprefix(".endc", curr_line)) {
            in_control = FALSE;
            ptr = ptr->li_next;
            continue;
        }
        if (in_control || curr_line[0] == '.' || curr_line[0] == '*') {
            ptr = ptr->li_next;
            continue;
        }

        num_terminals = get_number_terminals(curr_line);

        if (num_terminals <= 0) {
            ptr = ptr->li_next;
            continue;
        }

        for (i = 0; i < num_params; i++) {
            str_ptr = curr_line;

            for (j = 0; j < num_terminals+1; j++) {
                str_ptr = skip_non_ws(str_ptr);
                str_ptr = skip_ws(str_ptr);
            }

            while ((str_ptr = strstr(str_ptr, param_names[i])) != NULL) {
                /* make sure actually have the parameter name */
                char before = *(str_ptr-1);
                char after  = *(str_ptr+strlen(param_names[i]));
                if (!(is_arith_char(before) || isspace(before) || (str_ptr-1) < curr_line) ||
                        !(is_arith_char(after)  || isspace(after)  || after == '\0')) {
                    str_ptr = str_ptr + 1;
                    continue;
                }
                beg = str_ptr - 1;
                end = str_ptr + strlen(param_names[i]);
                if ((isspace(*beg) || *beg == '=') &&
                    (isspace(*end) || *end == '\0' || *end == ')'))
                { /* ')' added, any risks? h_vogt */
                    *str_ptr = '\0';
                    if (*end != '\0') {
                        new_str = TMALLOC(char, strlen(curr_line) + strlen(param_names[i]) + strlen(end) + 3);
                        sprintf(new_str, "%s{%s}%s", curr_line, param_names[i], end);
                    } else {
                        new_str = TMALLOC(char, strlen(curr_line) + strlen(param_names[i]) + 3);
                        sprintf(new_str, "%s{%s}", curr_line, param_names[i]);
                    }
                    str_ptr = new_str + strlen(curr_line) + strlen(param_names[i]);

                    tfree(ptr->li_line);
                    curr_line = ptr->li_line = new_str;
                }
                str_ptr++;
            }
        }
        ptr = ptr->li_next;
    }

    ind = 0;
    for (i = 0; i <= max_level; i++)
        for (j = num_params-1; j >= 0; j--)
            if (level[j] == i)
                if (param_skip[j] == 0)
                    ptr_array_ordered[ind++] = ptr_array[j];

    num_params -= skipped;
    if (ind != num_params) {
        fprintf(stderr, "ERROR: found wrong number of parameters during levelization ( %d instead of %d parameter s)!\n", ind, num_params);
        controlled_exit(EXIT_FAILURE);
    }

    /* fix next ptrs */
    ptr                                      = card_bf_start->li_next;
    card_bf_start->li_next                   = ptr_array_ordered[0];
    ptr_array_ordered[num_params-1]->li_next = ptr;
    for (i = 0; i < num_params-1; i++)
        ptr_array_ordered[i]->li_next = ptr_array_ordered[i+1];

    // clean up memory
    for (i = 0; i < num_params; i++) {
        tfree(param_names[i]);
        tfree(param_strs[i]);
    }

    tfree(level);
    tfree(param_skip);
    tfree(param_names);
    tfree(param_strs);

    for (i = 0; i< arr_size; i++)
        tfree(depends_on[i]);
    tfree(depends_on);

    tfree(ptr_array);
    tfree(ptr_array_ordered);

}


static void
inp_add_params_to_subckt(struct line *subckt_card)
{
    struct line *card        = subckt_card->li_next;
    char        *curr_line   = card->li_line;
    char        *subckt_line = subckt_card->li_line;
    char        *new_line, *param_ptr, *subckt_name, *end_ptr;

    while (card != NULL && ciprefix(".param", curr_line)) {
        param_ptr = strchr(curr_line, ' ');
        param_ptr = skip_ws(param_ptr);

        if (!strstr(subckt_line, "params:")) {
            new_line = TMALLOC(char, strlen(subckt_line) + strlen("params: ") + strlen(param_ptr) + 2);
            sprintf(new_line, "%s params: %s", subckt_line, param_ptr);

            subckt_name = skip_non_ws(subckt_card->li_line);
            subckt_name = skip_ws(subckt_name);
            end_ptr = skip_non_ws(subckt_name);
            subckt_w_params[num_subckt_w_params++] = copy_substring(subckt_name, end_ptr);
        } else {
            new_line = TMALLOC(char, strlen(subckt_line) + strlen(param_ptr) + 2);
            sprintf(new_line, "%s %s", subckt_line, param_ptr);
        }

        tfree(subckt_line);
        subckt_card->li_line = subckt_line = new_line;

        *curr_line = '*';

        card      = card->li_next;
        curr_line = card->li_line;
    }
}


/*
 * process a sequence of decks
 *   starting from a         `.suckt' deck
 *   upto the corresponding  `.ends'  deck
 * return a pointer to the terminating `.ends' deck
 *
 * recursivly descend
 *   when another `.subckt' is found
 *
 * parameters are removed from the main list
 *   and collected into a local list `first_param_card'
 * then processed and reinserted into the main list
 *
 */

static struct line *
inp_reorder_params_subckt(struct line *subckt_card)
{
    struct line *first_param_card = NULL;
    struct line *last_param_card = NULL;

    struct line *prev_card = subckt_card;
    struct line *c         = subckt_card->li_next;

    /* move .param lines to beginning of deck */
    while (c != NULL) {

        char *curr_line = c->li_line;

        if (*curr_line == '*') {
            c = c->li_next;
            continue;
        }

        if (ciprefix(".subckt", curr_line)) {
            prev_card = inp_reorder_params_subckt(c);
            c         = prev_card->li_next;
            continue;
        }

        if (ciprefix(".ends", curr_line)) {
            if (first_param_card) {
                inp_sort_params(first_param_card, last_param_card, subckt_card, subckt_card, c);
                inp_add_params_to_subckt(subckt_card);
            }
            return c;
        }

        if (ciprefix(".param", curr_line)) {
            if (first_param_card)
                last_param_card->li_next = c;
            else
                first_param_card = c;

            last_param_card    = c;
            prev_card->li_next = c->li_next;
            c                  = c->li_next;

            last_param_card->li_next = NULL;
            continue;
        }

        prev_card = c;
        c         = c->li_next;
    }

    /* the terminating `.ends' deck wasn't found */
    controlled_exit(EXIT_FAILURE);
    return NULL;
}


static void
inp_reorder_params(struct line *deck, struct line *list_head, struct line *end)
{
    struct line *first_param_card = NULL;
    struct line *last_param_card = NULL;

    struct line *prev_card = list_head;
    struct line *c = deck;

    /* move .param lines to beginning of deck */
    while (c != NULL) {

        char *curr_line = c->li_line;

        if (*curr_line == '*') {
            c = c->li_next;
            continue;
        }

        if (ciprefix(".subckt", curr_line)) {
            prev_card = inp_reorder_params_subckt(c);
            c         = prev_card->li_next;
            continue;
        }

        /* check for an unexpected extra `.ends' deck */
        if (ciprefix(".ends", curr_line)) {
            fprintf(stderr, "Error: Unexpected extra .ends in line:\n  %s.\n", curr_line);
            controlled_exit(EXIT_FAILURE);
        }

        if (ciprefix(".param", curr_line)) {
            if (first_param_card)
                last_param_card->li_next = c;
            else
                first_param_card = c;

            last_param_card    = c;
            prev_card->li_next = c->li_next;
            c                  = c->li_next;

            last_param_card->li_next = NULL;
            continue;
        }

        prev_card = c;
        c = c->li_next;
    }

    inp_sort_params(first_param_card, last_param_card, list_head, deck, end);
}


// iterate through deck and find lines with multiply defined parameters
//
// split line up into multiple lines and place those new lines immediately
// afetr the current multi-param line in the deck
static int
inp_split_multi_param_lines(struct line *deck, int line_num)
{
    struct line *card = deck, *param_end = NULL, *param_beg = NULL, *prev = NULL, *tmp_ptr;
    char *curr_line, *equal_ptr, *beg_param, *end_param, *new_line;
    char *array[5000];
    int  counter = 0, i;
    bool get_expression = FALSE, get_paren_expression = FALSE;
    char keep;

    while (card != NULL) {
        curr_line = card->li_line;

        if (*curr_line == '*') {
            card = card->li_next;
            continue;
        }

        if (ciprefix(".param", curr_line)) {
            counter = 0;
            while ((equal_ptr = strchr(curr_line, '=')) != NULL) {
                // check for equality '=='
                if (*(equal_ptr+1) == '=') {
                    curr_line = equal_ptr+2;
                    continue;
                }
                // check for '!=', '<=', '>='
                if (*(equal_ptr-1) == '!' || *(equal_ptr-1) == '<' || *(equal_ptr-1) == '>') {
                    curr_line = equal_ptr+1;
                    continue;
                }
                counter++;
                curr_line = equal_ptr + 1;
            }
            if (counter <= 1) {
                card = card->li_next;
                continue;
            }

            // need to split multi param line
            curr_line = card->li_line;
            counter   = 0;
            while (curr_line < card->li_line+strlen(card->li_line) && (equal_ptr = strchr(curr_line, '=')) != NULL) {
                // check for equality '=='
                if (*(equal_ptr+1) == '=') {
                    curr_line = equal_ptr+2;
                    continue;
                }
                // check for '!=', '<=', '>='
                if (*(equal_ptr-1) == '!' || *(equal_ptr-1) == '<' || *(equal_ptr-1) == '>') {
                    curr_line = equal_ptr+1;
                    continue;
                }

                beg_param = skip_back_ws(equal_ptr - 1);
                beg_param = skip_back_non_ws(beg_param);
                end_param = skip_ws(equal_ptr + 1);
                while (*end_param != '\0' && (!isspace(*end_param) || get_expression || get_paren_expression)) {
                    if (*end_param == '{')
                        get_expression = TRUE;
                    if (*end_param == '(')
                        get_paren_expression = TRUE;
                    if (*end_param == '}')
                        get_expression = FALSE;
                    if (*end_param == ')')
                        get_paren_expression = FALSE;
                    end_param++;
                }
                beg_param++;
                keep       = *end_param;
                *end_param = '\0';
                new_line   = TMALLOC(char, strlen(".param ") + strlen(beg_param) + 1);
                sprintf(new_line, ".param %s", beg_param);
                array[counter++] = new_line;
                *end_param = keep;
                curr_line = end_param;
            }
            tmp_ptr = card->li_next;
            for (i = 0; i < counter; i++) {
                if (param_end) {
                    param_end->li_next = alloc(struct line);
                    param_end          = param_end->li_next;
                } else {
                    param_end = param_beg = alloc(struct line);
                }
                param_end->li_next    = NULL;
                param_end->li_error   = NULL;
                param_end->li_actual  = NULL;
                param_end->li_line    = array[i];
                param_end->li_linenum = line_num++;
            }
            // comment out current multi-param line
            *(card->li_line)   = '*';
            // insert new param lines immediately after current line
            tmp_ptr            = card->li_next;
            card->li_next      = param_beg;
            param_end->li_next = tmp_ptr;
            // point 'card' pointer to last in scalar list
            card               = param_end;

            param_beg = param_end = NULL;
        } // if (ciprefix(".param", curr_line))
        prev = card;
        card = card->li_next;
    } // while (card != NULL)
    if (param_end) {
        prev->li_next = param_beg;
        prev          = param_end;
    }
    return line_num;
}


/* ps compatibility:
   ECOMP 3 0 TABLE {V(1,2)} = (-1 0V) (1, 10V)
   -->
   ECOMP 3 0 int3 int0 1
   BECOMP int3 int0 V = pwl(V(1,2), -2, 0, -1, 0 , 1, 10V, 2, 10V)

   GD16 16 1 TABLE {V(16,1)} ((-100V,-1pV)(0,0)(1m,1u)(2m,1m))
   -->
   GD16 16 1 int_16 int_1 1
   BGD16 int_16 int_1 V = pwl (v(16,1) , -100V , -1pV , 0 , 0 , 1m , 1u , 2m , 1m)
*/

/* hs compatibility:
   Exxx n1 n2 VCVS n3 n4 gain --> Exxx n1 n2 n3 n4 gain
   Gxxx n1 n2 VCCS n3 n4 tr --> Gxxx n1 n2 n3 n4 tr

   Two step approach to keep the original names for reuse,
   i.e. for current measurements like i(Exxx):
   Exxx n1 n2 VOL = {equation}
   -->
   Exxx n1 n2 int1 0 1
   BExxx int1 0 V = {equation}

   Gxxx n1 n2 CUR = {equation}
   -->
   Gxxx n1 n2 int1 0 1
   BGxxx int1 0 V = {equation}

   Do the following transformations only if {equation} contains
   simulation output like v(node), v(node1, node2), i(branch).
   Otherwise let do numparam the substitutions (R=const is handled
   in inp2r.c).

   Rxxx n1 n2 R = {equation} or Rxxx n1 n2 {equation}
   -->
   BRxxx n1 n2 I = V(n1,n2)/{equation}

   Unfortunately the capability for ac noise calculation of
   resistance may be lost.

   Cxxx n1 n2 C = {equation} or Cxxx n1 n2 {equation}
      -->
      Exxx  n-aux 0  n1 n2  1
      Cxxx  n-aux 0         1
      Bxxx  n2 n1  I = i(Exxx) * equation

   Lxxx n1 n2 L = {equation} or Lxxx n1 n2 {equation}
      -->
      Fxxx n-aux 0  Bxxx -1
      Lxxx n-aux 0      1
      Bxxx n1 n2 V = v(n-aux) * 1e-16

     */
static void
inp_compat(struct line *deck)
{
    char *str_ptr, *cut_line, *title_tok, *node1, *node2;
    char *out_ptr, *exp_ptr, *beg_ptr, *end_ptr, *copy_ptr, *del_ptr;
    char *xline;
    size_t xlen, i, pai = 0, paui = 0, ii;
    char *ckt_array[100];
    struct line *new_line, *tmp_ptr;

    struct line  *param_end = NULL, *param_beg = NULL;
    struct line *card;
    int skip_control = 0;

    char *equation, *tc1_ptr = NULL, *tc2_ptr = NULL;
    double tc1 = 0.0, tc2 = 0.0;

    for (card = deck; card; card = card->li_next) {

        char *curr_line = card->li_line;

        /* exclude any command inside .control ... .endc */
        if (ciprefix(".control", curr_line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", curr_line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }

        if (*curr_line == '*')
            continue;

        if (*curr_line == 'e') {
            /*    Exxx n1 n2 VCVS n3 n4 gain --> Exxx n1 n2 n3 n4 gain
                  remove vcvs */
            replace_token(curr_line, "vcvs", 4, 7);

            /* Exxx n1 n2 value={equation}
               -->
               Exxx n1 n2   vol={equation} */
            if ((str_ptr = strstr(curr_line, "value=")) != NULL) {
                *str_ptr = ' ';
                *(str_ptr + 1) = ' ';
                *(str_ptr + 2) = 'v';
                *(str_ptr + 3) = 'o';
                *(str_ptr + 4) = 'l';
            }
            /* Exxx n1 n2 TABLE {expression} = (x0, y0) (x1, y1) (x2, y2)
               -->
               Exxx n1 n2 int1 0 1
               BExxx int1 0 V = pwl (expression, x0-(x2-x0)/2, y0, x0, y0, x1, y1, x2, y2, x2+(x2-x0)/2, y2)
             */
            if ((str_ptr = strstr(curr_line, "table")) != NULL) {
                char *expression, *firstno, *ffirstno, *secondno, *midline, *lastno, *lastlastno;
                double fnumber, lnumber, delta;
                int nerror;
                cut_line = curr_line;
                /* title and nodes */
                title_tok = gettok(&cut_line);
                node1 =  gettok(&cut_line);
                node2 =  gettok(&cut_line);
                // Exxx  n1 n2 int1 0 1
                xlen = 2*strlen(title_tok) + strlen(node1) + strlen(node2)
                       + 20 - 4*2 + 1;
                ckt_array[0] = TMALLOC(char, xlen);
                sprintf(ckt_array[0], "%s %s %s %s_int1 0 1",
                        title_tok, node1, node2, title_tok);
                // get the expression
                str_ptr = gettok(&cut_line); /* ignore 'table' */
                if (cieq(str_ptr, "table")) {
                    tfree(str_ptr);
                } else {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                str_ptr =  gettok_char(&cut_line, '{', FALSE, FALSE);
                expression = gettok_char(&cut_line, '}', TRUE, TRUE); /* expression */
                if ((!expression) || (!str_ptr)) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                else
                    tfree(str_ptr);
                /* remove '{' and '}' from expression */
                if ((str_ptr = strchr(expression, '{')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(expression, '}')) != NULL)
                    *str_ptr = ' ';
                /* cut_line may now have a '=', if yes, it will have '{' and '}'
                   (braces around token after '=') */
                if ((str_ptr = strchr(cut_line, '=')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(cut_line, '{')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(cut_line, '}')) != NULL)
                    *str_ptr = ' ';
                /* get first two numbers to establish extrapolation */
                str_ptr = cut_line;
                ffirstno = gettok_node(&cut_line);
                if (!ffirstno) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                firstno = copy(ffirstno);
                fnumber = INPevaluate(&ffirstno, &nerror, TRUE);
                secondno = gettok_node(&cut_line);
                midline = cut_line;
                cut_line = strrchr(str_ptr, '(');
                /* replace '(' with ',' and ')' with ' ' */
                for (; *str_ptr; str_ptr++)
                    if (*str_ptr == '(')
                        *str_ptr = ',';
                    else if (*str_ptr == ')')
                        *str_ptr = ' ';
                /* scan for last two numbers */
                lastno = gettok_node(&cut_line);
                lnumber = INPevaluate(&lastno, &nerror, FALSE);
                /* check for max-min and take half the difference for delta */
                delta = (lnumber-fnumber)/2.;
                lastlastno = gettok_node(&cut_line);
                if (!secondno || (*midline == 0) || (delta <= 0.) || !lastlastno) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                xlen = 2*strlen(title_tok) + strlen(expression) + 14 + strlen(firstno) +
                    2*strlen(secondno) + strlen(midline) + 14 +
                    strlen(lastlastno) + 50;
                ckt_array[1] = TMALLOC(char, xlen);
                sprintf(ckt_array[1], "b%s %s_int1 0 v = pwl(%s, %e, %s, %s, %s, %s, %e, %s)",
                        title_tok, title_tok, expression, fnumber-delta, secondno, firstno, secondno,
                        midline, lnumber + delta, lastlastno);

                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (i = 0; i < 2; i++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[i];
                    param_end->li_linenum = 0;
                }
                // comment out current variable e line
                *(card->li_line)   = '*';
                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
                tfree(firstno);
                tfree(lastlastno);
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
            }
            /* Exxx n1 n2 VOL = {equation}
               -->
               Exxx n1 n2 int1 0 1
               BExxx int1 0 V = {equation}
            */
            if ((str_ptr = strstr(curr_line, "vol")) != NULL) {
                cut_line = curr_line;
                /* title and nodes */
                title_tok = gettok(&cut_line);
                node1 =  gettok(&cut_line);
                node2 =  gettok(&cut_line);
                /* Find equation, starts with '{', till end of line */
                str_ptr = strchr(cut_line, '{');
                if (str_ptr == NULL) {
                    fprintf(stderr, "ERROR: mal formed E line: %s\n", curr_line);
                    controlled_exit(EXIT_FAILURE);
                }

                // Exxx  n1 n2 int1 0 1
                xlen = 2*strlen(title_tok) + strlen(node1) + strlen(node2)
                       + 20 - 4*2 + 1;
                ckt_array[0] = TMALLOC(char, xlen);
                sprintf(ckt_array[0], "%s %s %s %s_int1 0 1",
                        title_tok, node1, node2, title_tok);
                // BExxx int1 0 V = {equation}
                xlen = 2*strlen(title_tok) + strlen(str_ptr)
                       + 20 - 3*2 + 1;
                ckt_array[1] = TMALLOC(char, xlen);
                sprintf(ckt_array[1], "b%s %s_int1 0 v = %s",
                        title_tok, title_tok, str_ptr);

                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (i = 0; i < 2; i++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[i];
                    param_end->li_linenum = 0;
                }
                // comment out current variable e line
                *(card->li_line)   = '*';
                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
            }
        } else if (*curr_line == 'g') {
            /* Gxxx n1 n2 VCCS n3 n4 tr --> Gxxx n1 n2 n3 n4 tr
               remove vccs */
            replace_token(curr_line, "vccs", 4, 7);

            /* Gxxx n1 n2 value={equation}
               -->
               Gxxx n1 n2   cur={equation} */
            if ((str_ptr = strstr(curr_line, "value=")) != NULL) {
                *str_ptr = ' ';
                *(str_ptr + 1) = ' ';
                *(str_ptr + 2) = 'c';
                *(str_ptr + 3) = 'u';
                *(str_ptr + 4) = 'r';
            }

            /* Gxxx n1 n2 TABLE {expression} = (x0, y0) (x1, y1) (x2, y2)
               -->
               Gxxx n1 n2 int1 0 1
               BGxxx int1 0 V = pwl (expression, x0-(x2-x0)/2, y0, x0, y0, x1, y1, x2, y2, x2+(x2-x0)/2, y2)
             */
            if ((str_ptr = strstr(curr_line, "table")) != NULL) {
                char *expression, *firstno, *ffirstno, *secondno, *midline, *lastno, *lastlastno;
                double fnumber, lnumber, delta;
                int nerror;
                cut_line = curr_line;
                /* title and nodes */
                title_tok = gettok(&cut_line);
                node1 =  gettok(&cut_line);
                node2 =  gettok(&cut_line);
                // Gxxx  n1 n2 int1 0 1
                xlen = 2*strlen(title_tok) + strlen(node1) + strlen(node2)
                       + 20 - 4*2 + 1;
                ckt_array[0] = TMALLOC(char, xlen);
                sprintf(ckt_array[0], "%s %s %s %s_int1 0 1",
                        title_tok, node1, node2, title_tok);
                // get the expression
                str_ptr = gettok(&cut_line); /* ignore 'table' */
                if (cieq(str_ptr, "table")) {
                    tfree(str_ptr);
                } else {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                str_ptr =  gettok_char(&cut_line, '{', FALSE, FALSE);
                expression = gettok_char(&cut_line, '}', TRUE, TRUE); /* expression */
                if ((!expression) || (!str_ptr)) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                else
                    tfree(str_ptr);
                /* remove '{' and '}' from expression */
                if ((str_ptr = strchr(expression, '{')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(expression, '}')) != NULL)
                    *str_ptr = ' ';
                /* cut_line may now have a '=', if yes, it will have '{' and '}'
                   (braces around token after '=') */
                if ((str_ptr = strchr(cut_line, '=')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(cut_line, '{')) != NULL)
                    *str_ptr = ' ';
                if ((str_ptr = strchr(cut_line, '}')) != NULL)
                    *str_ptr = ' ';
                /* get first two numbers to establish extrapolation */
                str_ptr = cut_line;
                ffirstno = gettok_node(&cut_line);
                if (!ffirstno) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                firstno = copy(ffirstno);
                fnumber = INPevaluate(&ffirstno, &nerror, TRUE);
                secondno = gettok_node(&cut_line);
                midline = cut_line;
                cut_line = strrchr(str_ptr, '(');
                /* replace '(' with ',' and ')' with ' ' */
                for (; *str_ptr; str_ptr++)
                    if (*str_ptr == '(')
                        *str_ptr = ',';
                    else if (*str_ptr == ')')
                        *str_ptr = ' ';
                /* scan for last two numbers */
                lastno = gettok_node(&cut_line);
                lnumber = INPevaluate(&lastno, &nerror, FALSE);
                /* check for max-min and take half the difference for delta */
                delta = (lnumber-fnumber)/2.;
                lastlastno = gettok_node(&cut_line);
                if (!secondno || (*midline == 0) || (delta <= 0.) || !lastlastno) {
                    fprintf(stderr, "Error: bad sytax in line %d\n  %s\n",
                        card->li_linenum_orig, card->li_line);
                    controlled_exit(EXIT_BAD);
                }
                /* BGxxx int1 0 V = pwl (expression, x0-(x2-x0)/2, y0, x0, y0, x1, y1, x2, y2, x2+(x2-x0)/2, y2) */
                xlen = 2*strlen(title_tok) + strlen(expression) + 14 + strlen(firstno) +
                    2*strlen(secondno) + strlen(midline) + 14 +
                    strlen(lastlastno) + 50;
                ckt_array[1] = TMALLOC(char, xlen);
                sprintf(ckt_array[1], "b%s %s_int1 0 v = pwl(%s, %e, %s, %s, %s, %s, %e, %s)",
                        title_tok, title_tok, expression, fnumber-delta, secondno, firstno, secondno,
                        midline, lnumber + delta, lastlastno);

                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (i = 0; i < 2; i++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[i];
                    param_end->li_linenum = 0;
                }
                // comment out current variable e line
                *(card->li_line)   = '*';
                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
                tfree(firstno);
                tfree(lastlastno);
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
            }
            /*
               Gxxx n1 n2 CUR = {equation}
               -->
               Gxxx n1 n2 int1 0 1
               BGxxx int1 0 V = {equation}
            */
            if ((str_ptr = strstr(curr_line, "cur")) != NULL) {
                cut_line = curr_line;
                /* title and nodes */
                title_tok = gettok(&cut_line);
                node1 =  gettok(&cut_line);
                node2 =  gettok(&cut_line);
                /* Find equation, starts with '{', till end of line */
                str_ptr = strchr(cut_line, '{');
                if (str_ptr == NULL) {
                    fprintf(stderr, "ERROR: mal formed G line: %s\n", curr_line);
                    controlled_exit(EXIT_FAILURE);
                }
                // Gxxx  n1 n2 int1 0 1
                xlen = 2*strlen(title_tok) + strlen(node1) + strlen(node2)
                       + 20 - 4*2 + 1;
                ckt_array[0] = TMALLOC(char, xlen);
                sprintf(ckt_array[0], "%s %s %s %s_int1 0 1",
                        title_tok, node1, node2, title_tok);
                // BGxxx int1 0 V = {equation}
                xlen = 2*strlen(title_tok) + strlen(str_ptr)
                       + 20 - 3*2 + 1;
                ckt_array[1] = TMALLOC(char, xlen);
                sprintf(ckt_array[1], "b%s %s_int1 0 v = %s",
                        title_tok, title_tok, str_ptr);

                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (i = 0; i < 2; i++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[i];
                    param_end->li_linenum = 0;
                }
                // comment out current variable g line
                *(card->li_line)   = '*';
                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
            }
        }

        /* F element compatibility */
        else if (*curr_line == 'f') {
            /* Fxxx n1 n2 CCCS vnam gain --> Fxxx n1 n2 vnam gain
               remove cccs */
            replace_token(curr_line, "cccs", 4, 6);
        }
        /* H element compatibility */
        else if (*curr_line == 'h') {
            /* Hxxx n1 n2 CCVS vnam transres --> Hxxx n1 n2 vnam transres
               remove cccs */
            replace_token(curr_line, "ccvs", 4, 6);
        }

        /* Rxxx n1 n2 R = {equation} or Rxxx n1 n2 {equation}
           -->
           BRxxx pos neg I = V(pos, neg)/{equation}
        */
        else if (*curr_line == 'r') {
            cut_line = curr_line;
            /* make BRxxx pos neg I = V(pos, neg)/{equation}*/
            title_tok = gettok(&cut_line);
            node1 =  gettok(&cut_line);
            node2 =  gettok(&cut_line);
            /* check only after skipping Rname and nodes, either may contain time (e.g. Rtime)*/
            if ((!strstr(cut_line, "v(")) &&  (!strstr(cut_line, "i(")) &&
                (!strstr(cut_line, "temper")) &&  (!strstr(cut_line, "hertz")) &&
                (!strstr(cut_line, "time"))) {
                /* no handling in B-Source, so we have to prepare ternary fcn
                   for numparam */
                if (strchr(curr_line, '?') && strchr(curr_line, ':'))
                    card->li_line = inp_fix_ternary_operator_str(curr_line, TRUE);
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
                continue;
            }

            /* Find equation, starts with '{', till end of line */
            str_ptr = strchr(cut_line, '{');
            if (str_ptr == NULL) {
                fprintf(stderr, "ERROR: mal formed R line: %s\n", curr_line);
                controlled_exit(EXIT_FAILURE);
            }
            equation = gettok_char(&str_ptr, '}', TRUE, TRUE);
            str_ptr = strstr(cut_line, "tc1");
            if (str_ptr) {
                /* We need to have 'tc1=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc1_ptr = strchr(str_ptr, '=');
                    if (tc1_ptr)
                        tc1 = atof(tc1_ptr+1);
                }
            }
            str_ptr = strstr(cut_line, "tc2");
            if (str_ptr) {
                /* We need to have 'tc2=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc2_ptr = strchr(str_ptr, '=');
                    if (tc2_ptr)
                        tc2 = atof(tc2_ptr+1);
                }
            }
            if ((tc1_ptr == NULL) && (tc2_ptr == NULL)) {
                xlen = strlen(title_tok) + strlen(node1) + strlen(node2) +
                       strlen(node1) + strlen(node2) + strlen(equation)  +
                       28 - 6*2 + 1;
                xline = TMALLOC(char, xlen);
                sprintf(xline, "b%s %s %s i = v(%s, %s)/(%s)", title_tok, node1, node2,
                        node1, node2, equation);
            } else if (tc2_ptr == NULL) {
                xlen = strlen(title_tok) + strlen(node1) + strlen(node2) +
                       strlen(node1) + strlen(node2) + strlen(equation)  +
                       28 - 6*2 + 1 + 21 + 13;
                xline = TMALLOC(char, xlen);
                sprintf(xline, "b%s %s %s i = v(%s, %s)/(%s) tc1=%15.8e reciproctc=1", title_tok, node1, node2,
                        node1, node2, equation, tc1);
            } else {
                xlen = strlen(title_tok) + strlen(node1) + strlen(node2) +
                       strlen(node1) + strlen(node2) + strlen(equation)  +
                       28 - 6*2 + 1 + 21 + 21 + 13;
                xline = TMALLOC(char, xlen);
                sprintf(xline, "b%s %s %s i = v(%s, %s)/(%s) tc1=%15.8e tc2=%15.8e reciproctc=1", title_tok, node1, node2,
                        node1, node2, equation, tc1, tc2);
            }
            tc1_ptr = NULL;
            tc2_ptr = NULL;
            new_line = alloc(struct line);
            new_line->li_next    = NULL;
            new_line->li_error   = NULL;
            new_line->li_actual  = NULL;
            new_line->li_line    = xline;
            new_line->li_linenum = 0;
            // comment out current old R line
            *(card->li_line)   = '*';
            // insert new B source line immediately after current line
            tmp_ptr           = card->li_next;
            card->li_next     = new_line;
            new_line->li_next = tmp_ptr;
            // point 'card' pointer to the new line
            card              = new_line;
            tfree(title_tok);
            tfree(node1);
            tfree(node2);
            tfree(equation);
        }
        /* Cxxx n1 n2 C = {equation} or Cxxx n1 n2 {equation}
           -->
           Exxx  n-aux 0  n1 n2  1
           Cxxx  n-aux 0         1
           Bxxx  n2 n1  I = i(Exxx) * equation
         */
        else if (*curr_line == 'c') {
            cut_line = curr_line;
            title_tok = gettok(&cut_line);
            node1 =  gettok(&cut_line);
            node2 =  gettok(&cut_line);
            /* check only after skipping Cname and nodes, either may contain time (e.g. Ctime)*/
            if ((!strstr(cut_line, "v(")) &&  (!strstr(cut_line, "i(")) &&
                (!strstr(cut_line, "temper")) &&  (!strstr(cut_line, "hertz")) &&
                (!strstr(cut_line, "time")))
            {
                /* no handling in B-Source, so we have to prepare ternary fcn
                   for numparam */
                if (strchr(curr_line, '?') && strchr(curr_line, ':'))
                    card->li_line = inp_fix_ternary_operator_str(curr_line, TRUE);
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
                continue;
            }

            /* Find equation, starts with '{', till end of line */
            str_ptr = strchr(cut_line, '{');
            if (str_ptr == NULL) {
                fprintf(stderr, "ERROR: mal formed C line: %s\n", curr_line);
                controlled_exit(EXIT_FAILURE);
            }
            equation = gettok_char(&str_ptr, '}', TRUE, TRUE);
            str_ptr = strstr(cut_line, "tc1");
            if (str_ptr) {
                /* We need to have 'tc1=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc1_ptr = strchr(str_ptr, '=');
                    if (tc1_ptr)
                        tc1 = atof(tc1_ptr+1);
                }
            }
            str_ptr = strstr(cut_line, "tc2");
            if (str_ptr) {
                /* We need to have 'tc2=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc2_ptr = strchr(str_ptr, '=');
                    if (tc2_ptr)
                        tc2 = atof(tc2_ptr+1);
                }
            }
            // Exxx  n-aux 0  n1 n2  1
            xlen = 2*strlen(title_tok) + strlen(node1) + strlen(node2)
                   + 21 - 4*2 + 1;
            ckt_array[0] = TMALLOC(char, xlen);
            sprintf(ckt_array[0], "e%s %s_int2 0 %s %s 1",
                    title_tok, title_tok, node1, node2);
            // Cxxx  n-aux 0  1
            xlen = 2*strlen(title_tok)
                   + 15 - 2*2 + 1;
            ckt_array[1] = TMALLOC(char, xlen);
            sprintf(ckt_array[1], "c%s %s_int2 0 1", title_tok, title_tok);
            // Bxxx  n2 n1  I = i(Exxx) * equation
            if ((tc1_ptr == NULL) && (tc2_ptr == NULL)) {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 27 - 2*5 + 1;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s i = i(e%s) * (%s)",
                        title_tok, node2, node1, title_tok, equation);
            } else if (tc2_ptr == NULL) {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 27 - 2*5 + 1 + 21 + 13;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s i = i(e%s) * (%s) tc1=%15.8e reciproctc=1",
                        title_tok, node2, node1, title_tok, equation, tc1);
            } else {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 27 - 2*5 + 1 + 21 + 21 + 13;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s i = i(e%s) * (%s) tc1=%15.8e tc2=%15.8e reciproctc=1",
                        title_tok, node2, node1, title_tok, equation, tc1, tc2);
            }
            tc1_ptr = NULL;
            tc2_ptr = NULL;
            // insert new B source line immediately after current line
            tmp_ptr = card->li_next;
            for (i = 0; i < 3; i++) {
                if (param_end) {
                    param_end->li_next = alloc(struct line);
                    param_end          = param_end->li_next;
                } else {
                    param_end = param_beg = alloc(struct line);
                }
                param_end->li_next    = NULL;
                param_end->li_error   = NULL;
                param_end->li_actual  = NULL;
                param_end->li_line    = ckt_array[i];
                param_end->li_linenum = 0;
            }
            // comment out current variable capacitor line
            *(card->li_line)   = '*';
            // insert new param lines immediately after current line
            tmp_ptr            = card->li_next;
            card->li_next      = param_beg;
            param_end->li_next = tmp_ptr;
            // point 'card' pointer to last in scalar list
            card               = param_end;

            param_beg = param_end = NULL;
            tfree(title_tok);
            tfree(node1);
            tfree(node2);
            tfree(equation);
        }

        /* Lxxx n1 n2 L = {equation} or Lxxx n1 n2 {equation}
           -->
           Fxxx n-aux 0  Bxxx -1
           Lxxx n-aux 0      1
           Bxxx n1 n2 V = v(n-aux) * equation
         */
        else if (*curr_line == 'l') {
            cut_line = curr_line;
            /* title and nodes */
            title_tok = gettok(&cut_line);
            node1 =  gettok(&cut_line);
            node2 =  gettok(&cut_line);
            if ((!strstr(cut_line, "v(")) &&  (!strstr(cut_line, "i(")) &&
                (!strstr(cut_line, "temper")) &&  (!strstr(cut_line, "hertz")) &&
                (!strstr(cut_line, "time")))
            {
                /* no handling in B-Source, so we have to prepare ternary fcn
                   for numparam */
                if (strchr(curr_line, '?') && strchr(curr_line, ':'))
                    card->li_line = inp_fix_ternary_operator_str(curr_line, TRUE);
                tfree(title_tok);
                tfree(node1);
                tfree(node2);
                continue;
            }

            /* Find equation, starts with '{', till end of line */
            str_ptr = strchr(cut_line, '{');
            if (str_ptr == NULL) {
                fprintf(stderr, "ERROR: mal formed L line: %s\n", curr_line);
                controlled_exit(EXIT_FAILURE);
            }
            equation = gettok_char(&str_ptr, '}', TRUE, TRUE);
            str_ptr = strstr(cut_line, "tc1");
            if (str_ptr) {
                /* We need to have 'tc1=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc1_ptr = strchr(str_ptr, '=');
                    if (tc1_ptr)
                        tc1 = atof(tc1_ptr+1);
                }
            }
            str_ptr = strstr(cut_line, "tc2");
            if (str_ptr) {
                /* We need to have 'tc2=something */
                if (str_ptr[3] && (isspace(str_ptr[3]) || (str_ptr[3] == '='))) {
                    tc2_ptr = strchr(str_ptr, '=');
                    if (tc2_ptr)
                        tc2 = atof(tc2_ptr+1);
                }
            }
            // Fxxx  n-aux 0  Bxxx  1
            xlen = 3*strlen(title_tok)
                   + 20 - 3*2 + 1;
            ckt_array[0] = TMALLOC(char, xlen);
            sprintf(ckt_array[0], "f%s %s_int2 0 b%s -1",
                    title_tok, title_tok, title_tok);
            // Lxxx  n-aux 0  1
            xlen = 2*strlen(title_tok)
                   + 15 - 2*2 + 1;
            ckt_array[1] = TMALLOC(char, xlen);
            sprintf(ckt_array[1], "l%s %s_int2 0 1", title_tok, title_tok);
            // Bxxx  n1 n2  V = v(n-aux) * equation
            if ((tc1_ptr == NULL) && (tc2_ptr == NULL)) {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 31 - 2*5 + 1;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s v = v(%s_int2) * (%s)",
                        title_tok, node1, node2, title_tok, equation);
            } else if (tc2_ptr == NULL) {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 31 - 2*5 + 1 + 21 + 13;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s v = v(%s_int2) * (%s) tc1=%15.8e reciproctc=0",
                        title_tok, node2, node1, title_tok, equation, tc1);
            } else {
                xlen = 2*strlen(title_tok) + strlen(node2) + strlen(node1)
                       + strlen(equation) + 31 - 2*5 + 1 + 21 + 21 + 13;
                ckt_array[2] = TMALLOC(char, xlen);
                sprintf(ckt_array[2], "b%s %s %s v = v(%s_int2) * (%s) tc1=%15.8e tc2=%15.8e reciproctc=0",
                        title_tok, node2, node1, title_tok, equation, tc1, tc2);
            }
            tc1_ptr = NULL;
            tc2_ptr = NULL;
            // insert new B source line immediately after current line
            tmp_ptr = card->li_next;
            for (i = 0; i < 3; i++) {
                if (param_end) {
                    param_end->li_next = alloc(struct line);
                    param_end          = param_end->li_next;
                } else {
                    param_end = param_beg = alloc(struct line);
                }
                param_end->li_next    = NULL;
                param_end->li_error   = NULL;
                param_end->li_actual  = NULL;
                param_end->li_line    = ckt_array[i];
                param_end->li_linenum = 0;
            }
            // comment out current variable inductor line
            *(card->li_line)   = '*';
            // insert new param lines immediately after current line
            tmp_ptr            = card->li_next;
            card->li_next      = param_beg;
            param_end->li_next = tmp_ptr;
            // point 'card' pointer to last in scalar list
            card               = param_end;

            param_beg = param_end = NULL;
            tfree(title_tok);
            tfree(node1);
            tfree(node2);
            tfree(equation);
        }
        /* .probe -> .save
           .print, .plot, .save, .four,
           An ouput vector may be replaced by the following:
           myoutput=par('expression')
           .meas
           A vector out_variable may be replaced by
           par('expression')
        */
        else if (*curr_line == '.') {
            // replace .probe by .save
            if ((str_ptr = strstr(curr_line, ".probe")) != NULL)
                memcpy(str_ptr, ".save ", 6);

            /* Various formats for measure statement:
             * .MEASURE {DC|AC|TRAN} result WHEN out_variable=val
             * + <TD=td> <FROM=val> <TO=val>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
             *
             * .MEASURE {DC|AC|TRAN} result WHEN out_variable=out_variable2
             * + <TD=td> <FROM=val> <TO=val>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
             *
             * .MEASURE {DC|AC|TRAN} result FIND out_variable WHEN out_variable2=val
             * + <TD=td> <FROM=val> <TO=val>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
             *
             * .MEASURE {DC|AC|TRAN} result FIND out_variable WHEN out_variable2=out_variable3
             * + <TD=td>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
             *
             * .MEASURE {DC|AC|TRAN} result FIND out_variable AT=val
             * + <FROM=val> <TO=val>
             *
             * .MEASURE {DC|AC|TRAN} result {AVG|MIN|MAX|MIN_AT|MAX_AT|PP|RMS} out_variable
             * + <TD=td> <FROM=val> <TO=val>
             *
             * .MEASURE {DC|AC|TRAN} result INTEG<RAL> out_variable
             * + <TD=td> <FROM=val> <TO=val>
             *
             * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable AT=val
             *
             * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable WHEN out_variable2=val
             * + <TD=td>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
             *
             * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable WHEN out_variable2=out_variable3
             * + <TD=td>
             * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>

             The user may set any out_variable to par(' expr ').
             We have to replace this by v(pa_xx) and generate a B source line.

             * ----------------------------------------------------------------- */
            if (ciprefix(".meas", curr_line)) {
                if (strstr(curr_line, "par(") == NULL)
                    continue;
                cut_line = curr_line;
                // search for 'par('
                while ((str_ptr = strstr(cut_line, "par(")) != NULL) {
                    if (pai > 99) {
                        fprintf(stderr, "ERROR: More than 99 function calls to par()\n");
                        fprintf(stderr, "  Limited to 99 per input file\n");
                        controlled_exit(EXIT_FAILURE);
                    }

                    // we have ' par({ ... })', the right delimeter is a ' ' or '='
                    if (ciprefix(" par({", (str_ptr-1))) {
                        // find expression
                        beg_ptr = end_ptr = str_ptr + 5;
                        while ((*end_ptr != ' ') && (*end_ptr != '=') && (*end_ptr != '\0'))
                            end_ptr++;
                        exp_ptr = copy_substring(beg_ptr, end_ptr-2);
                        cut_line = str_ptr;
                        // generate node
                        out_ptr = TMALLOC(char, 6);
                        sprintf(out_ptr, "pa_%02d", (int)pai);
                        // Bout_ptr  out_ptr 0  V = v(expr_ptr)
                        xlen = 2*strlen(out_ptr) + strlen(exp_ptr)+ 15 - 2*3 + 1;
                        ckt_array[pai] = TMALLOC(char, xlen);
                        sprintf(ckt_array[pai], "b%s %s 0 v = %s",
                                out_ptr, out_ptr, exp_ptr);
                        ckt_array[++pai] = NULL;
                        // length of the replacement V(out_ptr)
                        xlen = strlen(out_ptr) + 4;
                        del_ptr = copy_ptr = TMALLOC(char, xlen);
                        sprintf(copy_ptr, "v(%s)", out_ptr);
                        // length of the replacement part in original line
                        xlen = strlen(exp_ptr) + 7;
                        // copy the replacement without trailing '\0'
                        for (ii = 0; ii < xlen; ii++)
                            if (*copy_ptr)
                                *cut_line++ = *copy_ptr++;
                            else
                                *cut_line++ = ' ';

                        tfree(del_ptr);
                        tfree(exp_ptr);
                        tfree(out_ptr);
                    }
                    // or we have '={par({ ... })}', the right delimeter is a ' '
                    else if (ciprefix("={par({", (str_ptr-2))) {
                        // find expression
                        beg_ptr = end_ptr = str_ptr + 5;
                        while ((*end_ptr != ' ') && (*end_ptr != '\0'))
                            end_ptr++;
                        exp_ptr = copy_substring(beg_ptr, end_ptr-3);
                        // generate node
                        out_ptr = TMALLOC(char, 6);
                        sprintf(out_ptr, "pa_%02d", (int)pai);
                        // Bout_ptr  out_ptr 0  V = v(expr_ptr)
                        xlen = 2*strlen(out_ptr) + strlen(exp_ptr)+ 15 - 2*3 + 1;
                        ckt_array[pai] = TMALLOC(char, xlen);
                        sprintf(ckt_array[pai], "b%s %s 0 v = %s",
                                out_ptr, out_ptr, exp_ptr);
                        ckt_array[++pai] = NULL;
                        // length of the replacement V(out_ptr)
                        xlen = strlen(out_ptr) + 4;
                        del_ptr = copy_ptr = TMALLOC(char, xlen);
                        sprintf(copy_ptr, "v(%s)", out_ptr);
                        // length of the replacement part in original line
                        xlen = strlen(exp_ptr) + 9;
                        // skip '='
                        cut_line++;
                        // copy the replacement without trailing '\0'
                        for (ii = 0; ii < xlen; ii++)
                            if (*copy_ptr)
                                *cut_line++ = *copy_ptr++;
                            else *cut_line++ = ' ';

                        tfree(del_ptr);
                        tfree(exp_ptr);
                        tfree(out_ptr);
                    } else {
                        // nothing to replace
                        cut_line = str_ptr + 1;
                        continue;
                    }

                } // while 'par'
                // no replacement done, go to next line
                if (pai == paui)
                    continue;
                // remove white spaces
                card->li_line = inp_remove_ws(curr_line);
                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (ii = paui; ii < pai; ii++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[ii];
                    param_end->li_linenum = 0;
                }

                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
                paui = pai;
            } else if ((ciprefix(".save", curr_line)) ||
                       (ciprefix(".four", curr_line)) ||
                       (ciprefix(".print", curr_line)) ||
                       (ciprefix(".plot", curr_line)))
            {
                if (strstr(curr_line, "par(") == NULL)
                    continue;
                cut_line = curr_line;
                // search for 'par('
                while ((str_ptr = strstr(cut_line, "par(")) != NULL) {
                    if (pai > 99) {
                        fprintf(stderr, "ERROR: More than 99 function calls to par()\n");
                        fprintf(stderr, "  Limited to 99 per input file\n");
                        controlled_exit(EXIT_FAILURE);
                    }

                    // we have ' par({ ... })'
                    if (ciprefix(" par({", (str_ptr-1))) {

                        // find expression
                        beg_ptr = end_ptr = str_ptr + 5;
                        while ((*end_ptr != ' ') && (*end_ptr != '\0'))
                            end_ptr++;
                        exp_ptr = copy_substring(beg_ptr, end_ptr-2);
                        cut_line = str_ptr;
                        // generate node
                        out_ptr = TMALLOC(char, 6);
                        sprintf(out_ptr, "pa_%02d", (int)pai);
                        // Bout_ptr  out_ptr 0  V = v(expr_ptr)
                        xlen = 2*strlen(out_ptr) + strlen(exp_ptr)+ 15 - 2*3 + 1;
                        ckt_array[pai] = TMALLOC(char, xlen);
                        sprintf(ckt_array[pai], "b%s %s 0 v = %s",
                                out_ptr, out_ptr, exp_ptr);
                        ckt_array[++pai] = NULL;
                        // length of the replacement V(out_ptr)
                        xlen = strlen(out_ptr) + 1;
                        del_ptr = copy_ptr = TMALLOC(char, xlen);
                        sprintf(copy_ptr, "%s", out_ptr);
                        // length of the replacement part in original line
                        xlen = strlen(exp_ptr) + 7;
                        // copy the replacement without trailing '\0'
                        for (ii = 0; ii < xlen; ii++)
                            if (*copy_ptr)
                                *cut_line++ = *copy_ptr++;
                            else
                                *cut_line++ = ' ';

                        tfree(del_ptr);
                        tfree(exp_ptr);
                        tfree(out_ptr);
                    }
                    // or we have '={par({ ... })}'
                    else if (ciprefix("={par({", (str_ptr-2))) {

                        // find myoutput
                        beg_ptr = end_ptr = str_ptr - 2;
                        while (*beg_ptr != ' ')
                            beg_ptr--;
                        out_ptr = copy_substring(beg_ptr + 1, end_ptr);
                        cut_line = beg_ptr + 1;
                        // find expression
                        beg_ptr = end_ptr = str_ptr + 5;
                        while ((*end_ptr != ' ') && (*end_ptr != '\0'))
                            end_ptr++;
                        exp_ptr = copy_substring(beg_ptr, end_ptr-3);
                        // Bout_ptr  out_ptr 0  V = v(expr_ptr)
                        xlen = 2*strlen(out_ptr) + strlen(exp_ptr)+ 15 - 2*3 + 1;
                        ckt_array[pai] = TMALLOC(char, xlen);
                        sprintf(ckt_array[pai], "b%s %s 0 v = %s",
                                out_ptr, out_ptr, exp_ptr);
                        ckt_array[++pai] = NULL;
                        // length of the replacement V(out_ptr)
                        xlen = strlen(out_ptr) + 1;
                        del_ptr = copy_ptr = TMALLOC(char, xlen);
                        sprintf(copy_ptr, "%s", out_ptr);
                        // length of the replacement part in original line
                        xlen = strlen(out_ptr) + strlen(exp_ptr) + 10;
                        // copy the replacement without trailing '\0'
                        for (ii = 0; ii < xlen; ii++)
                            if (*copy_ptr)
                                *cut_line++ = *copy_ptr++;
                            else *cut_line++ = ' ';

                        tfree(del_ptr);
                        tfree(exp_ptr);
                        tfree(out_ptr);
                    }
                    // nothing to replace
                    else
                        cut_line = str_ptr + 1;
                } // while 'par('
                // no replacement done, go to next line
                if (pai == paui)
                    continue;
                // remove white spaces
                card->li_line = inp_remove_ws(curr_line);
                // insert new B source line immediately after current line
                tmp_ptr = card->li_next;
                for (ii = paui; ii < pai; ii++) {
                    if (param_end) {
                        param_end->li_next = alloc(struct line);
                        param_end          = param_end->li_next;
                    } else {
                        param_end = param_beg = alloc(struct line);
                    }
                    param_end->li_next    = NULL;
                    param_end->li_error   = NULL;
                    param_end->li_actual  = NULL;
                    param_end->li_line    = ckt_array[ii];
                    param_end->li_linenum = 0;
                }
                // comment out current variable capacitor line
                // *(ckt_array[0])   = '*';
                // insert new param lines immediately after current line
                tmp_ptr            = card->li_next;
                card->li_next      = param_beg;
                param_end->li_next = tmp_ptr;
                // point 'card' pointer to last in scalar list
                card               = param_end;

                param_beg = param_end = NULL;
                paui = pai;
                // continue;
            } // if .print etc.
        } // if ('.')
    }
}

/* replace a token (length 4 char) in string by spaces, if it is found
    at the correct position and the total number of tokens is o.k. */
static void
replace_token(char *string, char *token, int wherereplace, int total)
{
    int count = 0, i;
    char *actstring = string;

    /* token to be replaced not in string */
    if (strstr(string, token) == NULL)
        return;

    /* get total number of tokens */
    while (*actstring) {
        txfree(gettok(&actstring));
        count++;
    }
    /* If total number of tokens correct */
    if (count == total) {
        actstring = string;
        for (i = 1; i < wherereplace; i++)
            txfree(gettok(&actstring));
        /* If token to be replaced at right position */
        if (ciprefix(token, actstring)) {
            *actstring = ' ';
            *(actstring + 1) = ' ';
            *(actstring + 2) = ' ';
            *(actstring + 3) = ' ';
        }
    }
}


/* lines for B sources: no parsing in numparam code, just replacement of parameters.
   Parsing done in B source parser.
   To achive this, do the following:
   Remove all '{' and '}' --> no parsing of equations in numparam
   Place '{' and '}' directly around all potential parameters,
   thus skip function names like exp (search for exp(to detect fcn name),
   functions containing nodes like v(node), v(node1, node2), i(branch)
   and other keywords. --> Only parameter replacement in numparam

*/

static void
inp_bsource_compat(struct line *deck)
{
    char *equal_ptr, *str_ptr, *tmp_char, *new_str, *final_str;
    char actchar;
    struct line *card, *new_line, *tmp_ptr;
    wordlist *wl = NULL, *wlist = NULL;
    char buf[512];
    size_t i, xlen, ustate = 0;
    int skip_control = 0;
    int error1;

    for (card = deck; card; card = card->li_next) {

        char *curr_line = card->li_line;

        /* exclude any command inside .control ... .endc */
        if (ciprefix(".control", curr_line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", curr_line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }

        if (*curr_line == 'b') {
            /* remove white spaces of everything inside {}*/
            card->li_line = inp_remove_ws(card->li_line);
            curr_line = card->li_line;
            /* store starting point for later parsing, beginning of {expression} */
            equal_ptr = strchr(curr_line, '=');
            /* check for errors */
            if (equal_ptr == NULL) {
                fprintf(stderr, "ERROR: mal formed B line: %s\n", curr_line);
                controlled_exit(EXIT_FAILURE);
            }
            /* find the m={m} token and remove it */
            if ((str_ptr = strstr(curr_line, "m={m}")) != NULL)
                memcpy(str_ptr, "     ", 5);
            /* scan the line and remove all '{' and '}' */
            str_ptr = curr_line;
            while (*str_ptr) {
                if ((*str_ptr == '{') || (*str_ptr == '}'))
                    *str_ptr = ' ';
                str_ptr++;
            }
            /* scan the expression */
            str_ptr = equal_ptr + 1;
            while (*str_ptr != '\0') {
                str_ptr = skip_ws(str_ptr);
                if (*str_ptr == '\0')
                    break;
                actchar = *str_ptr;
                wl_append_word(&wlist, &wl, NULL);
                if ((actchar == ',') || (actchar == '(') || (actchar == ')') ||
                    (actchar == '*') || (actchar == '/') || (actchar == '^') ||
                    (actchar == '+') || (actchar == '?') || (actchar == ':'))
                {
                    if ((actchar == '*') && (*(str_ptr+1) == '*')) {
                        actchar = '^';
                        str_ptr++;
                    }
                    buf[0] = actchar;
                    buf[1] = '\0';
                    wl->wl_word = copy(buf);
                    str_ptr++;
                    if (actchar == ')')
                        ustate = 0;
                    else
                        ustate = 1; /* we have an operator */
                } else if ((actchar == '>') || (actchar == '<') ||
                           (actchar == '!') || (actchar == '='))
                {
                    /* >=, <=, !=, ==, <>, ... */
                    char *beg = str_ptr++;
                    if ((*str_ptr == '=') || (*str_ptr == '<') || (*str_ptr == '>'))
                        str_ptr++;
                    wl->wl_word = copy_substring(beg, str_ptr);
                    ustate = 1; /* we have an operator */
                } else if ((actchar == '|') || (actchar == '&')) {
                    char *beg = str_ptr++;
                    if ((*str_ptr == '|') || (*str_ptr == '&'))
                        str_ptr++;
                    wl->wl_word = copy_substring(beg, str_ptr);
                    ustate = 1; /* we have an operator */
                } else if ((actchar == '-') && (ustate == 0)) {
                    buf[0] = actchar;
                    buf[1] = '\0';
                    wl->wl_word = copy(buf);
                    str_ptr++;
                    ustate = 1; /* we have an operator */
                } else if ((actchar == '-') && (ustate == 1)) {
                    wl->wl_word = copy("");
                    str_ptr++;
                    ustate = 2; /* place a '-' in front of token */
                } else if (isalpha(actchar)) {
                    /* unary -, change sign */
                    if (ustate == 2) {
                        i = 1;
                        buf[0] = '-';
                    } else {
                        i = 0;
                    }

                    if (((actchar == 'v') || (actchar == 'i')) && (*(str_ptr+1) == '(')) {
                        while (*str_ptr != ')') {
                            buf[i] = *str_ptr;
                            i++;
                            str_ptr++;
                        }
                        buf[i] = *str_ptr;
                        buf[i+1] = '\0';
                        wl->wl_word = copy(buf);
                        str_ptr++;
                    } else {
                        while (isalnum(*str_ptr) ||
                               (*str_ptr == '!') || (*str_ptr == '#') ||
                               (*str_ptr == '$') || (*str_ptr == '%') ||
                               (*str_ptr == '_') || (*str_ptr == '[') ||
                               (*str_ptr == ']'))
                        {
                            buf[i] = *str_ptr;
                            i++;
                            str_ptr++;
                        }
                        buf[i] = '\0';
                        /* no parens {} around time, hertz, temper, the constants
                           pi and e which are defined in inpptree.c, around pwl and temp. coeffs */
                        if ((*str_ptr == '(') ||
                            cieq(buf, "hertz") || cieq(buf, "temper") ||
                            cieq(buf, "time") || cieq(buf, "pi") ||
                            cieq(buf, "e") || cieq(buf, "pwl"))
                        {
                            /* special handling of pwl lines:
                               Put braces around tokens and around expressions, use ','
                               as separator like:
                               pwl(i(Vin), {x0-1},{y0},
                                  {x0},{y0},{x1},{y1}, {x2},{y2},{x3},{y3},
                                  {x3+1},{y3})
                            */
                            /*
                             * if (cieq(buf, "pwl")) {
                             *     // go past i(Vin)
                             *     i = 3;
                             *     while (*str_ptr != ')') {
                             *         buf[i] = *str_ptr;
                             *         i++;
                             *         str_ptr++;
                             *     }
                             *     buf[i] = *str_ptr;
                             *     i++;
                             *     str_ptr++;
                             *     // find first ','
                             *     while (*str_ptr != ',') {
                             *         buf[i] = *str_ptr;
                             *         i++;
                             *         str_ptr++;
                             *     }
                             *     buf[i] = *str_ptr;
                             *     i++;
                             *     buf[i] = '{';
                             *     i++;
                             *     str_ptr++;
                             *     while (*str_ptr != ')') {
                             *         if (*str_ptr == ',') {
                             *             buf[i] = '}';
                             *             i++;
                             *             buf[i] = ',';
                             *             i++;
                             *             buf[i] = '{';
                             *             i++;
                             *             str_ptr++;
                             *         }
                             *         else {
                             *             buf[i] = *str_ptr;
                             *             i++;
                             *             str_ptr++;
                             *         }
                             *     }
                             *     buf[i] = '}';
                             *     i++;
                             *     buf[i] = *str_ptr;
                             *     i++;
                             *     buf[i] = '\0';
                             *     str_ptr++;
                             * }
                             */
                            wl->wl_word = copy(buf);

                        } else if (cieq(buf, "tc1") || cieq(buf, "tc2") ||
                                   cieq(buf, "reciproctc"))
                        {

                            str_ptr = skip_ws(str_ptr);
                            /* no {} around tc1 = or tc2 = , these are temp coeffs. */
                            if (str_ptr[0] == '='  &&  str_ptr[1] != '=') {
                                buf[i++] = '=';
                                buf[i] = '\0';
                                str_ptr++;
                                wl->wl_word = copy(buf);
                            } else {
                                xlen = strlen(buf);
                                tmp_char = TMALLOC(char, xlen + 3);
                                sprintf(tmp_char, "{%s}", buf);
                                wl->wl_word = tmp_char;
                            }

                        } else {
                            /* {} around all other tokens */
                            xlen = strlen(buf);
                            tmp_char = TMALLOC(char, xlen + 3);
                            sprintf(tmp_char, "{%s}", buf);
                            wl->wl_word = tmp_char;
                        }
                    }
                    ustate = 0; /* we have a number */
                } else if (isdigit(actchar)) {
                    /* allow 100p, 5MEG etc. */
                    double dvalue = INPevaluate(&str_ptr, &error1, 0);
                    char   cvalue[19];
                    /* unary -, change sign */
                    if (ustate == 2)
                        dvalue *= -1;
                    sprintf(cvalue, "%18.10e", dvalue);
                    wl->wl_word = copy(cvalue);
                    ustate = 0; /* we have a number */
                    /* skip the `unit', FIXME INPevaluate() should do this */
                    while (isalpha(*str_ptr))
                        str_ptr++;
                } else { /* strange char */
                    printf("Preparing B line for numparam\nWhat is this?\n%s\n", str_ptr);
                    buf[0] = *str_ptr;
                    buf[1] = '\0';
                    wl->wl_word = copy(buf);
                    str_ptr++;
                }
            }

            new_str = wl_flatten(wlist);
            wl_free(wlist);
            wlist = NULL;
            wl = NULL;

            tmp_char = copy(curr_line);
            equal_ptr = strchr(tmp_char, '=');
            if (str_ptr == NULL) {
                fprintf(stderr, "ERROR: mal formed B line:\n  %s\n", curr_line);
                controlled_exit(EXIT_FAILURE);
            }
            /* cut the tmp_char after the equal sign */
            *(equal_ptr + 1) = '\0';
            xlen = strlen(tmp_char) + strlen(new_str) + 2;
            final_str = TMALLOC(char, xlen);
            sprintf(final_str, "%s %s", tmp_char, new_str);

            new_line = alloc(struct line);
            new_line->li_next    = NULL;
            new_line->li_error   = NULL;
            new_line->li_actual  = NULL;
            new_line->li_line    = final_str;
            /* Copy old line numbers into new B source line */
            new_line->li_linenum = card->li_linenum;
            new_line->li_linenum_orig = card->li_linenum_orig;
            // comment out current line (old B source line)
            *(card->li_line)   = '*';
            // insert new B source line immediately after current line
            tmp_ptr           = card->li_next;
            card->li_next     = new_line;
            new_line->li_next = tmp_ptr;
            // point 'card' pointer to the new line
            card              = new_line;

            tfree(new_str);
            tfree(tmp_char);
        } /* end of if 'b' */
    } /* end of for loop */
}


/*
 * destructively fetch a token from the input string
 *   token is either quoted, or a plain nonwhitespace sequence
 * function will return the place from where to continue
 */

static char *
get_quoted_token(char *string, char **token)
{
    char *s = skip_ws(string);

    if (!*s)            /* nothing found */
        return string;

    if (isquote(*s)) {

        char *t = ++s;

        while (*t && !isquote(*t))
            t++;

        if (!*t) {        /* teriminator quote not found */
            *token = NULL;
            return string;
        }

        *t++ = '\0';

        *token = s;
        return t;

    } else {

        char *t = skip_non_ws(s);

        if (t == s) {     /* nothing found */
            *token = NULL;
            return string;
        }

        if (*t)
            *t++ = '\0';

        *token = s;
        return t;
    }
}


/* Option RSERIES=rval
 * Lxxx n1 n2 Lval
 * -->
 * Lxxx n1 n2_intern__ Lval
 * RLxxx_n2_intern__ n2_intern__ n2 rval
*/

static void
inp_add_series_resistor(struct line *deck)
{
    size_t skip_control = 0, xlen, i;
    bool has_rseries = FALSE;
    struct line *card;
    char *tmp_p, *title_tok, *node1, *node2, *rval = NULL;
    char *ckt_array[10];
    struct line  *tmp_ptr, *param_end = NULL, *param_beg = NULL;

    for (card = deck; card; card = card->li_next) {
        char *curr_line = card->li_line;
        if (*curr_line == '*')
            continue;
        if (strstr(curr_line, "option") && strstr(curr_line, "rseries"))
            has_rseries = TRUE;
        else
            continue;
        tmp_p = strstr(curr_line, "rseries");
        tmp_p += 7;
        /* default to "1e-3" if no value given */
        if (ciprefix("=", tmp_p)) {
            tmp_p = strchr(tmp_p, '=') + 1;
            rval = gettok(&tmp_p);
        }
        else
            rval = copy("1e-3");
    }

    if (!has_rseries || !rval)
        return;

    fprintf(stdout,
            "\nOption rseries given: \n"
            "resistor %s Ohms added in series to each inductor L\n\n", rval);

    for (card = deck; card; card = card->li_next) {
        char *cut_line;
        char *curr_line = cut_line = card->li_line;

        /* exclude any command inside .control ... .endc */
        if (ciprefix(".control", curr_line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", curr_line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }

        if (ciprefix("l", curr_line)) {
            title_tok = gettok(&cut_line);
            node1 =  gettok(&cut_line);
            node2 =  gettok(&cut_line);
            /* new L line */
            xlen = strlen(curr_line) + 10;
            ckt_array[0] = TMALLOC(char, xlen);
            sprintf(ckt_array[0], "%s %s %s_intern__ %s",
                    title_tok, node1, node2, cut_line);
            /* new R line */
            xlen = strlen(curr_line) + 19;
            ckt_array[1] = TMALLOC(char, xlen);
            sprintf(ckt_array[1], "R%s_intern__ %s_intern__ %s %s",
                    title_tok, node2, node2, rval);
            /* assemble new L and R lines */
            tmp_ptr = card->li_next;
            for (i = 0; i < 2; i++) {
                if (param_end) {
                    param_end->li_next = alloc(struct line);
                    param_end          = param_end->li_next;
                } else {
                    param_end = param_beg = alloc(struct line);
                }
                param_end->li_next    = NULL;
                param_end->li_error   = NULL;
                param_end->li_actual  = NULL;
                param_end->li_line    = ckt_array[i];
                param_end->li_linenum = 0;
            }
            // comment out current L line
            *(card->li_line)   = '*';
            // insert new new L and R lines immediately after current line
            tmp_ptr            = card->li_next;
            card->li_next      = param_beg;
            param_end->li_next = tmp_ptr;
            // point 'card' pointer to last in scalar list
            card               = param_end;
            param_beg = param_end = NULL;
            tfree(title_tok);
            tfree(node1);
            tfree(node2);
        }
    }

    tfree(rval);
}

/* If XSPICE option is not selected, run this function to alert and exit
   if the 'poly' option is found in e, g, f, or h controlled sources. */

#ifndef XSPICE
static void
inp_poly_err(struct line *deck)
{
    struct line *card;

    size_t skip_control = 0;
    for (card = deck; card; card = card->li_next) {
        char *curr_line = card->li_line;
        if (*curr_line == '*')
            continue;
        /* exclude any command inside .control ... .endc */
        if (ciprefix(".control", curr_line)) {
            skip_control ++;
            continue;
        } else if (ciprefix(".endc", curr_line)) {
            skip_control --;
            continue;
        } else if (skip_control > 0) {
            continue;
        }
        /* get the fourth token in a controlled source line and exit,
        if it is 'poly' */
        if ((ciprefix("e", curr_line)) || (ciprefix("g", curr_line)) ||
            (ciprefix("f", curr_line)) || (ciprefix("h", curr_line))) {
            txfree(gettok(&curr_line));
            txfree(gettok(&curr_line));
            txfree(gettok(&curr_line));
            if (ciprefix("poly", curr_line)) {
                fprintf(stderr,
                    "\nError: XSPICE is required to run the 'poly' option in line %d\n",
                    card->li_linenum_orig);
                fprintf(stderr, "  %s\n", card->li_line);
                fprintf(stderr, "\nSee manual chapt. 31 for installation instructions\n");
                controlled_exit(EXIT_BAD);
            }
        }
    }
}
#endif
