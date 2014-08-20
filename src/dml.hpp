#ifndef _DML_H_
#define _DML_H_

#include <array>

typedef struct {
	double A,B,C,D;
	double x;
} dml_s;

dml_s dml_init(double* parms, double* X, double* coef);
double dml_comp(dml_s dml, double*coef);
dml_s dml_init(double* parms, const std::array<int, 24>& X, int count, double* coef);

#endif
