#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <R_ext/RS.h>
#include <R.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

#define BETA1_SIZE 4
#define BETA2_SIZE 3
#define TOL 0.000000000001


void disco_chisq(int *freqs,int *n_genes,int *labels,int *n_labels,double *res);

void derivatives1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);
void derivatives2(double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);

void newton1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);
void newton2(double *f1,double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels);

void get_quantile_freqs(int *scores,int *n_genes,int *n_grp1,int *n_grp2,int *labels,int *n_labels,int *freqs);

void get_quantile_scores(double *exprs,int *n_genes,int *n_grp1,int *n_grp2,int *labels,int *n_labels,double *quantile_thresholds,int *scores);
