/*       mystring.c                Copyright (C)  2002    Georg Post
 *
 *  This file is part of Numparam, see:  readme.txt
 *  Free software under the terms of the GNU Lesser General Public License
 */

#include "ngspice/ngspice.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include "ngspice/config.h"
#include "ngspice/memory.h"
#ifdef HAS_WINGUI
#include "ngspice/wstdio.h"
#endif

#include "general.h"


/*
 * fetch a human answer to a y/n question from stdin
 * insist on a single non white-space char on a '\n' terminated line
 * return this char or '\n' or EOF
 * return '\0' if the answer doesn't fit this pattern
 */

int
yes_or_no(void)
{
    int first;

    do {
        first = getchar();
        if (first == '\n' || first == EOF)
            return first;
    } while (isspace(first));

    for (;;) {
        int c = getchar();
        if (c == EOF)
            return c;
        if (c == '\n')
            return tolower(first);
        if (!isspace(c))
            first = '\0';
    }
}


bool
ci_prefix(const char *p, const char *s)
{
    while (*p) {
        if ((isupper(*p) ? tolower(*p) : *p) !=
            (isupper(*s) ? tolower(*s) : *s))
            return (0);
        p++;
        s++;
    }

    return (1);
}


/*******  Strings ************
 *  are 0-terminated char arrays with a 2-byte trailer: max length.
 *  the string mini-library is "overflow-safe" under these conditions:
 *    use Str(n,s) macro: define and initialize a string s of maxlen n<255
 *    use sini() to initialize empty strings;  sfix() for non-empty ones.
 *    the Sini() macro does automatic sizing, for automatic char arrays
 *    to allocate a string on the heap, use newstring(n).
 *    use maxlen() and length() to retrieve string max and actual length
 *    use: cadd, cins, sadd, sins, scopy, pscopy to manipulate them
 *    never put '\x0' characters inside strings !
 *
 *    the 'killer idea' is the following:
 *    on string overflow and/or on heap allocation failure, a program
 *    MUST die.   Now we only die on a heap failure as with dynamic
 *    string we cannot have a string overflow.
 */

void
sfix(SPICE_DSTRINGPTR dstr_p, int len)
/* suppose s is allocated and filled with non-zero stuff */
{
    /* This function will now eliminate the max field.   The length of
     * the string is going to be i-1 and a null is going to be added
     * at the ith position to be compatible with old codel.  Also no
     * null characters will be present in the string leading up to the
     * NULL so this will make it a valid string. */
    int j;
    char *s;

    spice_dstring_setlength(dstr_p, len);
    s = spice_dstring_value(dstr_p);
    for (j = 0; j < len; j++)   /* eliminate null characters ! */
        if (s[j] == 0)
            s[j] = 1;
}


int
length(const char *s)
{
    return (int) strlen(s);
}


/* -----------------------------------------------------------------
 * Function: add string t to dynamic string dstr_p.
 * ----------------------------------------------------------------- */
bool
sadd(SPICE_DSTRINGPTR dstr_p, const char *t)
{
    spice_dstring_append(dstr_p, t, -1);
    return 1;
}


/* -----------------------------------------------------------------
 * Function: add character c to dynamic string dstr_p.
 * ----------------------------------------------------------------- */
bool
cadd(SPICE_DSTRINGPTR dstr_p, char c)
{
    char tmp_str[2];
    tmp_str[0] = c;
    tmp_str[1] = '\0';
    spice_dstring_append(dstr_p, tmp_str, -1);
    return 1;
}


/* -----------------------------------------------------------------
 * Function: insert character c at front of dynamic string dstr_p
 * ----------------------------------------------------------------- */
bool
cins(SPICE_DSTRINGPTR dstr_p, char c)
{
    int i;
    int ls;
    char *s_p;

    ls = spice_dstring_length(dstr_p);
    spice_dstring_setlength(dstr_p, ls+2); /* make sure we have space for char + EOS */
    s_p = spice_dstring_value(dstr_p);
    for (i = ls + 1; i >= 0; i--)
        s_p[i + 1] = s_p[i];
    s_p[0] = c;
    return 1;
}


/* -----------------------------------------------------------------
 * Function: insert string t at front of dynamic string dstr_p
 * ----------------------------------------------------------------- */
bool
sins(SPICE_DSTRINGPTR dstr_p, const char *t)
{
    int i;
    int ls;
    int lt;
    char *s_p;

    ls = spice_dstring_length(dstr_p);
    lt = length(t);
    spice_dstring_setlength(dstr_p, ls+lt+1); /* make sure we have space for string + EOS */
    s_p = spice_dstring_value(dstr_p);
    for (i = ls + 1; i >= 0; i--)
        s_p[i + lt] = s_p[i];

    for (i = 0; i < lt; i++)
        s_p[i] = t[i];
    return 1;
}


int
cpos(char c, char *s)
/* return position of c in s, or 0 if not found.
 * BUG, Pascal inherited: first char is at 1, not 0 !
 * No longer!  Now position is C-based to make life easier.
 */
{
    int i = 0;
    while ((s[i] != c) && (s[i] != '\0'))
        i++;

    if (s[i] == c)
        return i;
    else
        return -1;
}


char
upcase(char c)
{
    if ((c >= 'a') && (c <= 'z'))
        return (char) (c + 'A' - 'a');
    else
        return c;
}


/* -----------------------------------------------------------------
 * Create copy of the dynamic string.  Dynamic strings are always NULL
 * terminated.
 * ----------------------------------------------------------------- */
bool
scopyd(SPICE_DSTRINGPTR s, SPICE_DSTRINGPTR t)  /* returns success flag */
{
    spice_dstring_reinit(s);
    spice_dstring_append(s, spice_dstring_value(t), -1);
    return 1; /* Dstrings expand to any length */
}


/* -----------------------------------------------------------------
 * Create copy of the string in the dynamic string.  Dynamic strings
 * are always NULLterminated.
 * ----------------------------------------------------------------- */
bool
scopys(SPICE_DSTRINGPTR s, const char *t)     /* returns success flag */
{
    spice_dstring_reinit(s);
    spice_dstring_append(s, t, -1);
    return 1; /* Dstrings expand to any length */
}


/* -----------------------------------------------------------------
 * Create an upper case copy of a string and store it in a dynamic string.
 * Dynamic strings are always NULL * terminated.
 * ----------------------------------------------------------------- */
bool
scopy_up(SPICE_DSTRINGPTR dstr_p, const char *str)    /* returns success flag */
{
    char up[2];                 /* short string */
    const char *ptr;            /* position in string */

    spice_dstring_reinit(dstr_p);
    up[1] = '\0';
    for (ptr = str; ptr && *ptr; ptr++) {
        up[0] = upcase(*ptr);
        spice_dstring_append(dstr_p, up, 1);
    }
    return 1; /* Dstrings expand to any length */
}


/* -----------------------------------------------------------------
 * Create a lower case copy of a string and store it in a dynamic string.
 * Dynamic strings are always NULL * terminated.
 * ----------------------------------------------------------------- */
bool
scopy_lower(SPICE_DSTRINGPTR dstr_p, const char *str) /* returns success flag */
{
    char low[2];                /* short string */
    const char *ptr;            /* position in string */

    spice_dstring_reinit(dstr_p);
    low[1] = '\0';
    for (ptr = str; ptr && *ptr; ptr++) {
        low[0] = lowcase(*ptr);
        spice_dstring_append(dstr_p, low, 1);
    }
    return 1; /* Dstrings expand to any length */
}


bool
ccopy(SPICE_DSTRINGPTR dstr_p, char c)  /* returns success flag */
{
    char *s_p;                  /* current string */

    sfix(dstr_p, 1);
    s_p = spice_dstring_value(dstr_p);
    s_p[0] = c;
    return 1;
}


char *
pscopy(SPICE_DSTRINGPTR dstr_p, const char *t, int start, int leng)
/* partial string copy, with C-based start - Because we now have a 0 based
 * start and string may copy outselves, we may need to restore the first
 * character of the original dstring because resetting string will wipe
 * out first character. */
{
    int i;                      /* counter */
    int stop;                   /* stop value */
    char *s_p;                  /* value of dynamic string */

    stop = length(t);

    if (start < stop) {         /* nothing! */

        if ((start + leng - 1) > stop) {
            // leng = stop - start + 1;
            leng = stop - start;
        }

        _spice_dstring_setlength(dstr_p, leng);
        s_p = spice_dstring_value(dstr_p);

        for (i = 0; i < leng; i++)
            s_p[i] = t[start + i];

        s_p[leng] = '\0';

    } else {

        s_p = spice_dstring_reinit(dstr_p);

    }

    return s_p;
}


char *
pscopy_up(SPICE_DSTRINGPTR dstr_p, const char *t, int start, int leng)
/* partial string copy to upper case, with C convention for start. */
{
    int i;                      /* counter */
    int stop;                   /* stop value */
    char *s_p;                  /* value of dynamic string */

    stop = length(t);

    if (start < stop) {         /* nothing! */

        if ((start + leng - 1) > stop) {
            // leng = stop - start + 1;
            leng = stop - start;
        }

        _spice_dstring_setlength(dstr_p, leng);
        s_p = spice_dstring_value(dstr_p);

        for (i = 0; i < leng; i++)
            s_p[i] = upcase(t[start + i]);

        s_p[leng] = '\0';

    } else {

        s_p = spice_dstring_reinit(dstr_p);

    }

    return s_p;
}


bool
nadd(SPICE_DSTRINGPTR dstr_p, long n)
/* append a decimal integer to a string */
{
    int d[25];
    int j, k;
    char sg;                    /* the sign */
    char load_str[2];           /* used to load dstring */
    k = 0;

    if (n < 0) {
        n = -n;
        sg = '-';
    } else {
        sg = '+';
    }

    while (n > 0) {
        d[k] = (int)(n % 10);
        k++;
        n = n / 10;
    }

    if (k == 0) {
        cadd(dstr_p, '0');
    } else {
        load_str[1] = '\0';
        if (sg == '-') {
            load_str[0] = sg;
            spice_dstring_append(dstr_p, load_str, 1);
        }
        for (j = k - 1; j >= 0; j--) {
            load_str[0] = (char) ('0' + d[j]);
            spice_dstring_append(dstr_p, load_str, 1);
        }
    }

    return 1;
}


bool
naddll(SPICE_DSTRINGPTR dstr_p, long long n)
/* append a decimal integer (but a long long) to a string */
{
    int d[25];
    int j, k;
    char sg;                    /* the sign */
    char load_str[2];           /* used to load dstring */
    k = 0;

    if (n < 0) {
        n = -n;
        sg = '-';
    } else {
        sg = '+';
    }

    while (n > 0) {
        d[k] = (int) (n % 10);
        k++;
        n = n / 10;
    }

    if (k == 0) {
        cadd(dstr_p, '0');
    } else {
        load_str[1] = '\0';
        if (sg == '-') {
            load_str[0] = sg;
            spice_dstring_append(dstr_p, load_str, 1);
        }
        for (j = k - 1; j >= 0; j--) {
            load_str[0] = (char) ('0' + d[j]);
            spice_dstring_append(dstr_p, load_str, 1);
        }
    }

    return 1;
}


void
stri(long n, SPICE_DSTRINGPTR dstr_p)
/* convert integer to string */
{
    spice_dstring_reinit(dstr_p);
    nadd(dstr_p, n);
}


bool
steq(const char *a, const char *b)  /* string a==b test */
{
    return strcmp(a, b) == 0;
}


bool
stne(const char *s, const char *t)
{
    return strcmp(s, t) != 0;
}


char
lowcase(char c)
{
    if ((c >= 'A') && (c <= 'Z'))
        return (char) (c - 'A' + 'a');
    else
        return c;
}


bool
alfa(char c)
{
    return
        ((c >= 'a') && (c <= 'z')) ||
        ((c >= 'A') && (c <= 'Z')) ||
        c == '_' || c == '[' || c == ']';
}


bool
num(char c)
{
    return (c >= '0') && (c <= '9');
}


bool
alfanum(char c)
{
    return alfa(c) || ((c >= '0') && (c <= '9'));
}


int
freadstr(FILE * f, SPICE_DSTRINGPTR dstr_p)
/* read a line from a file.
   was BUG: long lines truncated without warning, ctrl chars are dumped.
   Bug no more as we can only run out of memory.  Removed max argument.
*/
{
    char c;
    char str_load[2];
    int len = 0;

    str_load[0] = '\0';
    str_load[1] = '\0';
    spice_dstring_reinit(dstr_p);

    do
    {
        c = (char) fgetc(f);    /*  tab is the only control char accepted */
        if (((c >= ' ') || (c < 0) || (c == '\t'))) {
            str_load[0] = c;
            spice_dstring_append(dstr_p, str_load, 1);
        }
    }
    while (!feof(f) && (c != '\n'));

    return len;
}


char *
stupcase(char *s)
{
    int i = 0;

    while (s[i] != '\0') {
        s[i] = upcase(s[i]);
        i++;
    }

    return s;
}


/*****  pointer tricks: app won't use naked malloc(), free() ****/

void
dispose(void *p)
{
    if (p != NULL)
        free(p);
}


void *
new(size_t sz)
{
    void *p = tmalloc(sz);

    if (p == NULL) {            /* fatal error */
        printf(" new() failure. Program halted.\n");
        controlled_exit(EXIT_FAILURE);
    }

    return p;
}


/***** elementary math *******/

double
absf(double x)
{
    if (x < 0.0)
        return -x;
    else
        return x;
}


long
absi(long i)
{
    if (i >= 0)
        return (i);
    else
        return (-i);
}


int
spos_(char *sub, const char *s)
{
    const char *ptr = strstr(s, sub);

    if (ptr)
        return (int) (ptr - s);
    else
        return -1;
}


#ifndef HAVE_LIBM

long
np_round(double x)
/* using <math.h>, it would be simpler: floor(x+0.5), see below */
{
    double u;
    long z;
    int n;
//  Str (40, s);
    SPICE_DSTRING s;

    spice_dstring_init(&s);
    u = 2e9;

    if (x > u)
        x = u;
    else if (x < -u)
        x = -u;

    n = sprintf(s, "%-12.0f", x);
    s[n] = 0;
    sscanf(s, "%ld", &z);

    return z;
}


long
np_trunc(double x)
{
    long n = np_round(x);

    if ((n > x) && (x >= 0.0))
        n--;
    else if ((n < x) && (x < 0.0))
        n++;

    return n;
}


#else /* use floor() and ceil() */


long
np_round(double r)
{
    return (long) floor(r + 0.5);
}


long
np_trunc(double r)
{
    if (r >= 0.0)
        return (long) floor(r);
    else
        return (long) ceil(r);
}

#endif
