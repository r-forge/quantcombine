#include "quant_combine.h"

void newton1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels) {

	int i=0,j=0;

	double abs_f1 = fabs(*(f1));

	//for dgesv
	int N = BETA1_SIZE;
	int NRHS = BETA1_SIZE;
	int LDA = BETA1_SIZE;
	int IPIV[BETA1_SIZE] = {0,0,0,0};
	int LDB = BETA1_SIZE;
	int INFO = 0;
	double B[] = { 	//BETA1_SIZE x BETA1_SIZE
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1
	};	
	
	//for dgemm
	char transa = 'N';
	char transb = 'N';
	int m = BETA1_SIZE;
	int n = 1;
	int k = BETA1_SIZE;
	double alpha = 1.0;
	//a is B
	int lda = BETA1_SIZE;
	//b is *(f1)
	int ldb = BETA1_SIZE;
	double beta = 1.0;
	double c[] = {0,0,0,0};	//BETA1_SIZE
	int ldc = BETA1_SIZE;
	

	while( abs_f1 > TOL ) {

		derivatives1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);

		abs_f1 = fabs(*(f1));


		for(i=0;i<BETA1_SIZE;i++) {

			for(j=0;j<BETA1_SIZE;j++) {

				*(d1 + BETA1_SIZE*i + j) *= -1;	 //multiply d1 by -1
				
				//reset B for dgesv
				if(i==j) {
					*(B + BETA1_SIZE*i + j) = 1;
				}
				else {
					*(B + BETA1_SIZE*i + j) = 0;
				}
			}
			IPIV[i] = 0; //reset of dgesv
			*(c+i) = 0; //reset for dgemm
		}

		INFO = 0; //reset for dgesv

		F77_CALL(dgesv)(&N,&NRHS,d1,&LDA,IPIV,B,&LDB,&INFO);
		
		if(INFO != 0) {
			error("Error occured in dgesv LAPACK routine in newton1 C function, which is part of the discoChiSq R function; return value is %d\n",INFO);
		}
		
		F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,B,&lda,f1,&ldb,&beta,c,&ldc);

		for(i=0;i<BETA1_SIZE;i++) {
			*(beta1+i) += *(c+i);
		}
	}
}



void newton2(double *f1,double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels) {

	int i=0,j=0;
	
	double abs_f1 = fabs(*(f1));
	double abs_f2 = fabs(*(f2));

	//for dgesv
	int N = BETA2_SIZE;
	int NRHS = BETA2_SIZE;
	int LDA = BETA2_SIZE;
	int IPIV[BETA2_SIZE] = {0,0,0};
	int LDB = BETA2_SIZE;
	int INFO = 0;
	double B[] = { //BETA2_SIZE x BETA2_SIZE
		1,0,0,
		0,1,0,
		0,0,1
	};

	//for dgemm
	char transa = 'N';
	char transb = 'N';
	int m = BETA2_SIZE;
	int n = 1;
	int k = BETA2_SIZE;
	double alpha = 1.0;
	//a is B
	int lda = BETA2_SIZE;
	//b is *(f2)
	int ldb = BETA2_SIZE;
	double beta = 1.0;
	double c[] = {0,0,0}; //BETA2_SIZE
	int ldc = BETA2_SIZE;

	if( abs_f1 > TOL ) {
		
		while( abs_f2 > TOL ) {
			
			derivatives2(f2,d2,labels,n_labels,beta2,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);

			abs_f2 = fabs(*(f2));

			for(i=0;i<BETA2_SIZE;i++) {
				
				for(j=0;j<BETA2_SIZE;j++) {
					
					*(d2 + BETA2_SIZE*i + j) *= -1;
					
					//reset B for dgesv
					if(i==j) {
						*(B + BETA2_SIZE*i + j) = 1;
					}
					else {
						*(B + BETA2_SIZE*i + j) = 0;
					}
					IPIV[i] = 0; //reset of dgesv
					*(c+i) = 0; //reset for dgemm
				}
			}
			
			INFO = 0; //reset for dgesv

			F77_CALL(dgesv)(&N,&NRHS,d2,&LDA,IPIV,B,&LDB,&INFO);
		
			if(INFO != 0) {
				error("Error occured in dgesv LAPACK routine in newton2 C function, which is part of discoChiSq R function; return value is %d\n",INFO);
			}
						
			F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,B,&lda,f2,&ldb,&beta,c,&ldc);
			
			for(i=0;i<BETA2_SIZE;i++) {
				*(beta2+i) += *(c+i);
			}
		}
	}
}

