#ifndef _DML_H_
#define _DML_H_

#include <array>

typedef struct {
	double A,B,C,D;
	double x;
} dml_s;

dml_s dml_init(float* parms, float* X, float* coef);
double dml_comp(dml_s dml, float*coef);
dml_s dml_init(float* parms, const std::array<int, 24>& X, int count, float* coef);

#endif
