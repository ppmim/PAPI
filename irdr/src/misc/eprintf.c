/* eprintf.c -- error handling tasks */

/*
 * eprintf functions are based on similar functions presented in 'The
 * Practice of Programming' by Brian W. Kernighan and Rob Pike.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include "eprintf.h"

/* eprintf: print error message and exit */
void eprintf(char *fmt, ...)
{
    va_list args;

    fflush(stdout);

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    if (fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':')
        fprintf(stderr, " %s", strerror(errno));

    fprintf(stderr, "\n");

    exit(1);
}

/* emalloc: malloc and report if error */
void *emalloc(size_t n)
{
    void *p = malloc(n);

    if (p == NULL)
        eprintf("emalloc: malloc of %u bytes failed:", n);

    return p;
}

/* erealloc: realloc and report if error */
void *erealloc(void *vp, size_t n)
{
    void *p = realloc(vp, n);

    if (p == NULL)
        eprintf("erealloc: realloc of %u bytes failed:", n);

    return p;
}

/* estrdup: duplicate a string, report if error */
char *estrdup(char *s)
{
    char *t = (char *) malloc(strlen(s)+1);

    if (t == NULL)
        eprintf("estrdup(\"%.20s\") failed:", s);

    strcpy(t, s);

    return t;
}

/* ecalloc: calloc and report if error */
void *ecalloc(size_t nmemb, size_t size)
{
    void *p = calloc(nmemb, size);

    if (p == NULL)
        eprintf("ecalloc: calloc of %u bytes failed:", size * nmemb);

    return p;
}
