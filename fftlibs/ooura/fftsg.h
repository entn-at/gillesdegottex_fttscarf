#ifndef __OOURA_H__
#define __OOURA_H__

#ifdef OOFLOAT_FLOAT
    #define OOFLOAT float
#else
    #define OOFLOAT double
#endif

void cdft(int, int, OOFLOAT *, int *, OOFLOAT *);
void rdft(int, int, OOFLOAT *, int *, OOFLOAT *);
void ddct(int, int, OOFLOAT *, int *, OOFLOAT *);
void ddst(int, int, OOFLOAT *, int *, OOFLOAT *);
void dfct(int, OOFLOAT *, OOFLOAT *, int *, OOFLOAT *);
void dfst(int, OOFLOAT *, OOFLOAT *, int *, OOFLOAT *);

#endif // __OOURA_H__
