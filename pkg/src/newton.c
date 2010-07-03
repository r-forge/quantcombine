#include "quant_combine.h"
#include <R_ext/Lapack.h>

void newton1(double *f1,double *d1,int *labels,int *n_labels,double *beta1,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels,long double tol) {

	double delta[4];

	//printf("initially in newton1, *(f1) is %2.13f\n",*(f1));

	double abs_f1 = fabs(*(f1));

	while( abs_f1 > tol ) {

		//double abs_f1 = 0;
		//abs_f1 = abs(*(f1));
		//printf("absf1 is %2.13f\n",abs_f1);
		//printf("fabs(*(f1)), which is %2.13f, is still > tol, which is %3.14f, so about to run derivatives1\n",abs_f1,tol);

		derivatives1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);
		
		//printf("finished running derivatives1\n");
		//printf("f1 now looks like %2.13f,%2.13f,%2.13f,%2.13f\n",*(f1),*(f1+1),*(f1+2),*(f1+3));

		int i=0,j=0;

		//printf("d1 initially looks like (after derivatives 1)\n");
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) printf("%2.13f\t",*(d1+4*i+j));
			//putchar('\n');
		}


		for(i=0;i<4;i++) {
			for(j=0;j<4;j++) {
				if (*(d1+4*i+j) == 0) {
					//*(d1 + 4*i + j) = 0;
					//printf("d1 i=%d,j=%d is 0\n",i,j);
				}
				else {
					*(d1 + 4*i + j) *= -1;
				}
			}
		}

		//printf("after multiplying by -1, d1 looks like\n");
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) printf("%5.8f\t",*(d1+4*i+j));
			//putchar('\n');
		}
		
		int N = 4;
		int NRHS = 4 ;
		int LDA = 4;
		int IPIV[4] = {0,0,0,0};
		int LDB = 4;
		int INFO;

		double B[] = {
			1,0,0,0,
			0,1,0,0,
			0,0,1,0,
			0,0,0,1
		};

		//printf("B initially looks like\n");
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) printf("%5.8f\t",B[i*4+j]);
			//putchar('\n');
		}

		//printf("populated data for dgesv, now will run\n");

		F77_CALL(dgesv)(&N,&NRHS,d1,&LDA,IPIV,B,&LDB,&INFO);
		
		//printf("finished running dgesv, INFO is %d, and values in IPIV are %d,%d,%d,%d\n",INFO,IPIV[0],IPIV[1],IPIV[2],IPIV[3]);
		
		//printf("after dgesv, B looks like\n");
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) printf("%5.8f\t",B[i*4+j]);
			//putchar('\n');
		}
		//printf("after dgesv, d1 looks like\n");
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) printf("%3.8f\t",*(d1+i*4+j));
			//putchar('\n');
		}


		if(INFO != 0) {
			error("Error occured in dgesv LAPACK routine in newton1 C function, return value is %d\n",INFO);
		}
		/*
		delta[0] = IPIV[0] * *(f1);
		delta[1] = IPIV[1] * *(f1+1);
		delta[2] = IPIV[2] * *(f1+2);
		delta[3] = IPIV[3] * *(f1+3);
		*/

		int f[4];
		char transa = 'N';
		char transb = 'N';
		int m = 4;
		int n = 1;
		int k = 4;
		double alpha = 1.0;
		//a is B
		int lda = 4;
		//b is *(f1)
		int ldb = 4;
		double beta = 1.0;
		//double  c[4][1] ;
		double c[] = {
			0,
			0,
			0,
			0
		};
		int ldc = 4;

		abs_f1 = fabs(*(f1));
		//printf("now fabs(*(f1)) is %f\n",abs_f1);

		//printf("before dgemm, f1 is %f,%f,%f,%f\n",*(f1),*(f1+1),*(f1+2),*(f1+3));
		
		
		
		F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,B,&lda,f1,&ldb,&beta,c,&ldc);

		//printf("after dgemm, f1 is %f,%f,%f,%f\n",*(f1),*(f1+1),*(f1+2),*(f1+3));

		for(i=0;i<4;i++) {
			delta[i] = *(c + i); //id[3] * *(f1+3);
		}

		//printf("delta contains %f,%f,%f,%f\n",delta[0],delta[1],delta[2],delta[3]);

		*(beta1) += delta[0];
		*(beta1+1) += delta[1];
		*(beta1+2) += delta[2];
		*(beta1+3) += delta[3];

		//printf("beta1 contains %f,%f,%f,%f\n",*(beta1),*(beta1+1),*(beta1+2),*(beta1+3));


	}
	//printf("now fabs(*(f1)), which is %f, is no longer > tol, which is %3.14f, so will return out of newton1 to disco_chisq\n",abs_f1,tol);
}



void newton2(double *f1,double *f2,double *d2,int *labels,int *n_labels,double *beta2,int *colsums_freqs,int *colsums_freqs_x_labels,int *colsums_freqs_x_labels_x_labels,long double tol) {
	
	double delta[3];
	
	//printf("initially in newton2, *(f1) is %2.13f\n",*(f1));
	//printf("initially in newton2, *(f2) is %2.13f\n",*(f2));

	double abs_f1 = fabs(*(f1));
	double abs_f2 = fabs(*(f2));

	if( abs_f1 > tol ) {
		
		//printf("fabs(*(f1)), which is %2.13f, is > tol, which is %3.14f, so about to loop while f2>tol\n",abs_f1,tol);

		while( abs_f2 > tol ) {
			
			//printf("fabs(*(f2)), which is %2.13f, is still > tol, which is %3.14f, so about to run derivatives2\n",abs_f2,tol);

			derivatives2(f2,d2,labels,n_labels,beta2,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);
			//printf("finished running derivatives2\n");
			//printf("f2 now looks like %2.13f,%2.13f,%2.13f\n",*(f2),*(f2+1),*(f2+2));

			//printf("before multiplying by 1, d2 looks like:\n");
			
			int i=0,j=0;
			for(i=0;i<3;i++) {
				for(j=0;j<3;j++) {
					//printf("d2[i,j] is %f\n",*(d2 + 3*i + j));
				}
			}

			for(i=0;i<3;i++) {
				for(j=0;j<3;j++) {
					if(*(d2+3*i+j) == 0) {
						//*(d2+3*i+j) = 0;
					}
					else {
						*(d2 + 3*i + j) *= -1;
					}
				}
			}
			
			//printf("after multiplying by -1, d2 looks like\n");
			for(i=0;i<3;i++) {
				//for(j=0;j<3;j++) printf("%5.8f\t",*(d2+3*i+j));
				//putchar('\n');
			}

			int N = 3;
			int NRHS = 3;
			int LDA = 3;
			int IPIV[3] = {0,0,0};
			int LDB = 3;
			int INFO;
			
			double B[] = {
				1,0,0,
				0,1,0,
				0,0,1
			};

			//printf("B initially looks like\n");
			for(i=0;i<3;i++) {
				//for(j=0;j<3;j++) printf("%5.8f\t",B[i*3+j]);
				//putchar('\n');
			}
			
			F77_CALL(dgesv)(&N,&NRHS,d2,&LDA,IPIV,B,&LDB,&INFO);
		
			//printf("finished running dgesv, INFO is %d, and values in IPIV are %d,%d,%d\n",INFO,IPIV[0],IPIV[1],IPIV[2]);
			
			//printf("after dgesv, B looks like\n");
			for(i=0;i<3;i++) {
				//for(j=0;j<3;j++) printf("%5.8f\t",B[i*3+j]);
				//putchar('\n');
			}
			
			if(INFO != 0) {
				error("Error occured in dgesv LAPACK routine in newton2 C function, return value is %d\n",INFO);
			}
			
			int f[3];
			char transa = 'N';
			char transb = 'N';
			int m = 3;
			int n = 1;
			int k = 3;
			double alpha = 1.0;
			//a is B
			int lda = 3;
			//b is *(f1)
			int ldb = 3;
			double beta = 1.0;
			double c[] = {
				0,
				0,
				0
			};
			int ldc = 3;
			
			//printf("before dgemm, f2 is %f,%f,%f\n",*(f2),*(f2+1),*(f2+2));
			
			
			
			F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,B,&lda,f2,&ldb,&beta,c,&ldc);
			
			for(i=0;i<3;i++) {
				//printf("c[i] is %f\n",*(c+i));
				delta[i] = *(c + i);
			}
			
			//printf("delta contains %f,%f,%f\n",delta[0],delta[1],delta[2]);
			
			*(beta2) += delta[0];
			*(beta2+1) += delta[1];
			*(beta2+2) += delta[2];
			
			//printf("beta2 contains %f,%f,%f\n",*(beta2),*(beta2+1),*(beta2+2));

			abs_f2 = fabs(*(f2));
		}
		//printf("now fabs(*(f2)), which is %2.13f, is no longer > tol, which is %2.14f, so will return out of newton2 to disco_chisq\n",abs_f2,tol);
		
	}
	else {
		//printf("fabs(*(f1), which is %2.13f is not > tol,which is %2.14f\n",abs_f1,tol);
	}
}

