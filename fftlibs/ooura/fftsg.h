#ifndef __OOURA_H__
#define __OOURA_H__

#if defined(OOFLOAT_SINGLE)
    #define OOFLOAT float
#elif defined(OOFLOAT_DOUBLE)
    #define OOFLOAT double
#elif defined(OOFLOAT_LONGDOUBLE)
    #define OOFLOAT long double
#endif

void cdft(int, int, OOFLOAT *, int *, OOFLOAT *);
void rdft(int, int, OOFLOAT *, int *, OOFLOAT *);
void ddct(int, int, OOFLOAT *, int *, OOFLOAT *);
void ddst(int, int, OOFLOAT *, int *, OOFLOAT *);
void dfct(int, OOFLOAT *, OOFLOAT *, int *, OOFLOAT *);
void dfst(int, OOFLOAT *, OOFLOAT *, int *, OOFLOAT *);

#endif // __OOURA_H__
