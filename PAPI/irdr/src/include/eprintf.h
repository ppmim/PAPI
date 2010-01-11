/* eprintf.h -- header file for eprintf.c */

/*
 * eprintf functions are based on similar functions presented in 'The
 * Practice of Programming' by Brian W. Kernighan and Rob Pike.
 */

extern void eprintf(char *fmt, ...);
extern char *estrdup(char *s);
extern void *emalloc(size_t n);
extern void *ecalloc(size_t nmemb, size_t size);
extern void *erealloc(void *vp, size_t n);
