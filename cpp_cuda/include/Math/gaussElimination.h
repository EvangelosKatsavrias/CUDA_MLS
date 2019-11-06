#ifndef GAUSSELIMINATIONHEADER
#define GAUSSELIMINATIONHEADER

template<class T>
void gaussElimination_backSubst(T* A, T* b, T* x, int n);
template<class T>
void gaussElimination_frontSubst(T* A, T* b, T* x, int n);


#endif
