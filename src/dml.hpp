#ifndef _DML_H_
#define _DML_H_

typedef struct {
	double A,B,C,D;
	double x;
} dml_s;

dml_s dml_init(float* parms, float* X, float* coef);
double dml_comp(dml_s dml, float*coef);

#endif
