#include "quant_combine.h"

void disco_chisq(int *freqs,int *n_genes,int *labels,int *n_labels,double *res) {

	int g = 0,i = 0, j = 0;

	int freq_length = 2 * *n_labels;

	double beta1[BETA1_SIZE] = {0,0,0,0};
	double beta2[BETA2_SIZE] = {0,0,0};
	double f1[] = {0,0,0,0};
	double f2[] = {0,0,0};
	double d1[] = {
		0,0,0,0,
		0,0,0,0,
		0,0,0,0,
		0,0,0,0
	};
	double d2[] = {
		0,0,0,
		0,0,0,
		0,0,0
	};


	for(g=0;g<*n_genes;g++) {

		int gene_offset = g * freq_length;
		int res_offset = 9 * g;
		
		int colsums_freqs[2] = {0,0};
		int colsums_freqs_x_labels[2] = {0,0};
		int colsums_freqs_x_labels_x_labels[2] = {0,0};
		
		for(i=0;i<*n_labels;i++) {

			colsums_freqs[0] += *(freqs+i+gene_offset);
			colsums_freqs[1] += *(freqs+i+*n_labels+gene_offset);
			
			colsums_freqs_x_labels[0] += *(freqs+i+gene_offset) * *(labels+i);
			colsums_freqs_x_labels[1] += *(freqs+i+*n_labels+gene_offset) * *(labels+i);
			
			colsums_freqs_x_labels_x_labels[0] += *(freqs+i+gene_offset) * *(labels+i) * *(labels+i);
			colsums_freqs_x_labels_x_labels[1] += *(freqs+i+*n_labels+gene_offset) * *(labels+i) * *(labels+i);
		}
		
		for(i=0;i<BETA1_SIZE;i++) {

			beta1[i] = 0;
			*(f1+i) = 0;
	
			for(j=0;j<BETA1_SIZE;j++) {
				*(d1 + BETA1_SIZE*i + j) = 0;
			}
		}

		for(i=0;i<BETA2_SIZE;i++) {

			beta2[i] = 0;
			*(f2+i) = 0;
			
			for(j=0;j<BETA2_SIZE;j++) {
				*(d2 + BETA2_SIZE*i + j) = 0;
			}		
		}



		derivatives1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);

		double f1b[] = {f1[0],f1[1],f1[2]};
		
		newton1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);
		
		f2[0] = f1b[0];
		
		newton2(f1b,f2,d2,labels,n_labels,beta2,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);
		
		double log1[2] = {0,0};
		double larp1[2] = {0,0};	
		
		for(i=0;i<2;i++) {
			
			double alpha[2];
			double bes[2];
			
			if (i==0) {
				alpha[0] = beta1[0];
				alpha[1] = beta1[1];
				
				bes[0] = beta1[2];
				bes[1] = beta1[3];
			}
			else if (i==1) {
				alpha[0] = beta2[0];
				alpha[1] = beta2[0];
				
				bes[0] = beta2[1];
				bes[1] = beta2[2];
			}
			
			double br[2] = {0,0};
			
			for(j=0;j<*n_labels;j++) {
				
				br[0] += exp( (labels[j] * alpha[0]) + (pow(labels[j],2) * bes[0] ) );
				br[1] += exp( (labels[j] * alpha[1]) + (pow(labels[j],2) * bes[1] ) );
			}
			
			larp1[i] = ((double)colsums_freqs[0] * log(br[0])) + ((double)colsums_freqs[1] * log(br[1]));
			
			log1[i] = ( (colsums_freqs_x_labels[0] * alpha[0]) + (colsums_freqs_x_labels[1] * alpha[1]) ) 
				+ ( (colsums_freqs_x_labels_x_labels[0] * bes[0]) + (colsums_freqs_x_labels_x_labels[1] * bes[1]) )
				- larp1[i];
		}
		
		double chisq = -2 * (log1[1] - log1[0]);
		
		double pvalue = 1 - pchisq(chisq,1,1,0);
		
		
		//return beta1,beta2,chisq,pvalue
		*(res + res_offset) = pvalue;
		*(res + 1 + res_offset) = chisq;
		*(res + 2 + res_offset) = beta1[0];
		*(res + 3 + res_offset) = beta1[1];
		*(res + 4 + res_offset) = beta1[2];
		*(res + 5 + res_offset) = beta1[3];
		*(res + 6 + res_offset) = beta2[0];
		*(res + 7 + res_offset) = beta2[1];
		*(res + 8 + res_offset) = beta2[2];
	}
}
