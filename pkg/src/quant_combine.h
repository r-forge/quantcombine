#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <R_ext/RS.h>
#include <R.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>


void disco_chisq(int *freqs,int *n_genes,int *labels,int *n_labels,double *res);

void derivatives1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);
void derivatives2(double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);

void newton1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels,long double tol);
void newton2(double *f1,double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels,long double tol);

