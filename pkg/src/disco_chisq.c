#include "quant_combine.h"

void disco_chisq(int *freqs,int *n_genes,int *labels,int *n_labels,double *res) {

	int freq_length = 2 * *n_labels;
	int g = 0;

	for(g=0;g<*n_genes;g++) {
		
		int gene_offset = g * freq_length;
		int res_offset = 9 * g;

		//printf("*n_labels is %d\n",*n_labels);
		int i = 0, j = 0;
		
		for(i=0;i<freq_length;i++) {
			//printf("freq is %d\n",*(freqs+i+gene_offset ));
		}
		
		int colsums_freqs[2] = {0,0};
		int colsums_freqs_x_labels[2] = {0,0};
		int colsums_freqs_x_labels_x_labels[2] = {0,0};
		
		//printf("just allocated colsums_freqs");
		for(i=0;i<*n_labels;i++) {
			//printf("i is %d\n",i);
			
			//printf("*(freqs+i) is %d, which will be added to colsums_freqs[0]\n",*(freqs+i+gene_offset));
			colsums_freqs[0] += *(freqs+i+gene_offset);
			//printf("colsums_freqs[0] is %d\n",colsums_freqs[0]);
			
			//printf("*(freqs+i+*n_labels) is %d, which will be added to colsums_freqs[1]\n",*(freqs+i+*n_labels+gene_offset));
			colsums_freqs[1] += *(freqs+i+*n_labels+gene_offset);
			//printf("colsums_freqs[1] is %d\n",colsums_freqs[1]);
			
			colsums_freqs_x_labels[0] += *(freqs+i+gene_offset) * *(labels+i);
			colsums_freqs_x_labels[1] += *(freqs+i+*n_labels+gene_offset) * *(labels+i);
			
			colsums_freqs_x_labels_x_labels[0] += *(freqs+i+gene_offset) * *(labels+i) * *(labels+i);
			colsums_freqs_x_labels_x_labels[1] += *(freqs+i+*n_labels+gene_offset) * *(labels+i) * *(labels+i);
			
		}
		//printf("finished creating colsums_freqs");
		
		
		double beta1[4] = {0,0,0,0};
		double beta2[3] = {0,0,0};
		double f1[] = {0,0,0,0};
		
		//double d1[4][4];
		double f2[] = {0,0,0};
		//double f2b[] = {f2[0],f2[1],f2[2]};
		//double d2[3][3];
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
		
		long double tol = 0.000000000001;
		
		derivatives1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels);
		//printf("finished derivatives first time\n");
		
		double f1b[] = {f1[0],f1[1],f1[2]};
		//printf("f1b is %2.13f,%2.13f,%2.13f,%2.13f\n",f1b[0],f1b[1],f1b[2],f1b[3]);
		
		double d1b[] = {
			d1[0],d1[1],d1[2],d1[3],
			d1[4],d1[5],d1[6],d1[7],
			d1[8],d1[9],d1[10],d1[11],
			d1[12],d1[13],d1[14],d1[15]
		};
		
		for(i=0;i<4;i++) {
			//for(j=0;j<4;j++) {printf("d1b[i,j] is %f\t",*(d1b+4*i+j)); }
			//putchar('\n');
		}
		
		
		for(i=0;i<4;i++) {
			for(j=0;j<4;j++) {
				
			}
		}
		
		//printf("f1b is %f,%f,%f,%f\n",f1b[0],f1b[1],f1b[2],f1b[3]);
		
		newton1(f1,d1,labels,n_labels,beta1,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels,tol);
		//printf("finished newton1\n");
		
		f2[0] = f1b[0];
		
		newton2(f1b,f2,d2,labels,n_labels,beta2,colsums_freqs,colsums_freqs_x_labels,colsums_freqs_x_labels_x_labels,tol);
		//printf("finished newton2 time\n");
		
		//printf("beta1 is %f,%f,%f,%f\n",beta1[0],beta1[1],beta1[2],beta1[3]);
		//printf("beta2 is %f,%f,%f\n",beta2[0],beta2[1],beta2[2]);
		
		
		double log1[2] = {0,0};
		//double nr[n_labels][2];
		double larp1[2] = {0,0};	
		
		for(i=0;i<2;i++) {
			//printf("i is %d\n",i);
			
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
			
			//printf("alpha is %f,%f\n",alpha[0],alpha[1]);
			//printf("bes is %f,%f\n",bes[0],bes[1]);
			
			double br[2] = {0,0};
			
			for(j=0;j<*n_labels;j++) {
				//printf("j is %d\n",j);
				
				br[0] += exp( (labels[j] * alpha[0]) + (pow(labels[j],2) * bes[0] ) );
				br[1] += exp( (labels[j] * alpha[1]) + (pow(labels[j],2) * bes[1] ) );
			}
			
			//printf("br[0] is %f, br[1] is %f\n",br[0],br[1]);
			
			//printf("colsums_freqs is %d and %d\n",colsums_freqs[0],colsums_freqs[1]);
			//printf("log(br) is %f and %f\n",log(br[0]),log(br[1]));
			
			double p1 = (double)colsums_freqs[0] * log(br[0]);
			double p2 = (double)colsums_freqs[1] * log(br[1]);
			//printf("colsums_freqs[0] * log(br[0] is %f\n",p1);				  
			//printf("colsums_freqs[1] * log(br[1] is %f\n",p2);				  

			larp1[i] = ((double)colsums_freqs[0] * log(br[0])) + ((double)colsums_freqs[1] * log(br[1]));
			
			//printf("larp1[%d] is:%f\n",i,larp1[i]);
			
			log1[i] = ( (colsums_freqs_x_labels[0] * alpha[0]) + (colsums_freqs_x_labels[1] * alpha[1]) ) 
				+ ( (colsums_freqs_x_labels_x_labels[0] * bes[0]) + (colsums_freqs_x_labels_x_labels[1] * bes[1]) )
				- larp1[i];
			//printf("log1[%d] is:%f\n",i,log1[i]);
		}
		
		//printf("larp1[0] is %f and larp1[1] is %f\n",larp1[0],larp1[1]);
		//printf("log1[0] is %f and log1[1] is %f\n",log1[0],log1[1]);
		
		
		double chisq = -2 * (log1[1] - log1[0]);
		
		double pvalue = 1 - pchisq(chisq,1,1,0);
		
		//printf("about to store results\n");
		
		//return beta1,beta1,chisq,pvalue
		*(res + res_offset) = chisq;
		*(res + 1 + res_offset) = pvalue;
		*(res + 2 + res_offset) = beta1[0];
		*(res + 3 + res_offset) = beta1[1];
		*(res + 4 + res_offset) = beta1[2];
		*(res + 5 + res_offset) = beta1[3];
		*(res + 6 + res_offset) = beta2[0];
		*(res + 7 + res_offset) = beta2[1];
		*(res + 8 + res_offset) = beta2[2];
		
		for(i=0;i<9;i++) {
			//printf("*(res + %d) is %f\n",i,*(res+i+res_offset));
		}
	}
}
