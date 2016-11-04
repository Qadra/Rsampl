#ifndef DEFSH
#define DEFSH

#include <R_ext/RS.h> /* For Calloc() etc. */
#include <R_ext/Random.h> /* For unif_rand() */

#ifndef XALLOC
#define XALLOC(n, t) Calloc(n, t)
#endif

#ifndef XFREE
#define XFREE(p) Free(p)
#endif

#ifndef XMEMZERO
#define XMEMZERO(p, n) Memzero(p, n)
#endif

#ifndef XREALLOC
#define XREALLOC(p, n, t) Realloc(p, n, t)
#endif

#ifndef XRANDFUN
#define XRANDFUN() unif_rand()
#endif

#endif
