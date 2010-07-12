#include "quant_combine.h"

void derivatives1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels) {
  
  double nr[*n_labels][2];

  double ar[2] = {0,0};
  double ar1[2] = {0,0};
  double ar2[2] = {0,0};
  double ar3[2] = {0,0};
  double ar4[2] = {0,0};

  int i=0,j=0;
  
  for(i=0;i<*n_labels;i++) {
	  	  
	  nr[i][0] = exp( 
		  ( *(labels+i) * *(beta1) ) 
		  + ( pow(*(labels+i),2) * *(beta1+2) ) 
		  );

	  nr[i][1] = exp( 
		  ( *(labels+i) * *(beta1+1) ) 
		  + ( pow(*(labels+i),2) * *(beta1+3) ) 
		  );

	  ar[0] += nr[i][0];
	  ar[1] += nr[i][1];

	  ar1[0] += labels[i] * nr[i][0];
	  ar1[1] += labels[i] * nr[i][1];

	  ar2[0] += pow(labels[i],2) * nr[i][0];
	  ar2[1] += pow(labels[i],2) * nr[i][1];

	  ar3[0] += pow(labels[i],3) * nr[i][0];
	  ar3[1] += pow(labels[i],3) * nr[i][1];

	  ar4[0] += pow(labels[i],4) * nr[i][0];
	  ar4[1] += pow(labels[i],4) * nr[i][1];
  }

  for(i=0;i<BETA1_SIZE;i++) {

	  *(f1+i) = 0;
	  
	  for(j=0;j<BETA1_SIZE;j++) {
		  *(d1 + BETA1_SIZE*i + j) = 0;
	  }
  }


  for(i=0;i<2;i++) {

	  *(f1+i) = *(colsums_freqs_x_labels+i) - *(colsums_freqs+i) * (ar1[i]/ar[i]);
	  
	  *(f1+i+2) = *(colsums_freqs_x_labels_x_labels+i) - *(colsums_freqs+i) * (ar2[i]/ar[i]);

	  *(d1+BETA1_SIZE*i+i) = -*(colsums_freqs+i) * ((-pow(ar1[i],2) + (ar[i]*ar2[i]))/pow(ar[i],2) );
		  
	  *(d1 + BETA1_SIZE*i + 2+i) = -*(colsums_freqs+i) * ( ( (-ar2[i]*ar1[i]) + (ar[i]*ar3[i]))/pow(ar[i],2));
	  
	  *(d1 + BETA1_SIZE*(2+i) + i) = *(d1 + BETA1_SIZE*i + 2+i);
	  
	  *(d1 + BETA1_SIZE*(2+i) + 2+i) = -*(colsums_freqs+i) * ( ( -pow(ar2[i],2) + (ar[i]*ar4[i]) )/pow(ar[i],2));

  }
  
  //f1 and d1 are used
}


void derivatives2(double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels) {

  
  double anr[*n_labels][2];

  double aar[2] = {0,0};
  double aar1[2] = {0,0};
  double aar2[2] = {0,0};
  double aar3[2] = {0,0};
  double aar4[2] = {0,0};

  int i=0,j=0;
  
  for(i=0;i<*n_labels;i++) {
	  	  
	  anr[i][0] = exp( ( *(labels+i) * *(beta2) ) + ( pow(*(labels+i),2) * *(beta2+1) ) );
	  anr[i][1] = exp( ( *(labels+i) * *(beta2) ) + ( pow(*(labels+i),2) * *(beta2+2) ) );

	  aar[0] += anr[i][0];
	  aar[1] += anr[i][1];

	  aar1[0] += labels[i] * anr[i][0];
	  aar1[1] += labels[i] * anr[i][1];

	  aar2[0] += pow(labels[i],2) * anr[i][0];
	  aar2[1] += pow(labels[i],2) * anr[i][1];

	  aar3[0] += pow(labels[i],3) * anr[i][0];
	  aar3[1] += pow(labels[i],3) * anr[i][1];

	  aar4[0] += pow(labels[i],4) * anr[i][0];
	  aar4[1] += pow(labels[i],4) * anr[i][1];
  }

  for(i=0;i<BETA2_SIZE;i++) {

	  *(f2+i) = 0;
	  
	  for(j=0;j<BETA2_SIZE;j++) {
		  *(d2 + BETA2_SIZE*i + j) = 0;
	  }
  }


  double acc1 = 0;
  double acc2 = 0;
  double xbb = 0;
  

  for(i=0;i<2;i++) {

	  acc1 = acc1 + *(colsums_freqs+i) * (aar1[i]/aar[i]);

	  acc2 = acc2 + *(colsums_freqs+i) * ( (-pow(aar1[i],2) + (aar[i] * aar2[i]))) / pow(aar[i],2);

	  *(f2+i+1) = *(colsums_freqs_x_labels_x_labels+i) - *(colsums_freqs+i) * (aar2[i]/aar[i]);

	  *(d2 + 4 + 4*i) = -*(colsums_freqs+i) * ( ( -pow(aar2[i],2) + (aar[i] * aar4[i]) ) / pow(aar[i],2));

	  *(d2 + i + 1) = -*(colsums_freqs+i) * ( ( (-aar2[i]*aar1[i]) + (aar[i]*aar3[i]) ) / pow(aar[i],2) ); 

	  *(d2 + 3*i + 3) = *(d2 + i + 1);
	  
	  xbb = xbb + *(colsums_freqs_x_labels+i);
  }
  
  *(f2) = xbb - acc1;
  *(d2) = -acc2;


  //f2 and d2 are used
}

