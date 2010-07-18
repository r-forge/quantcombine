#include "quant_combine.h"

void quantile(int *n_quantiles,int n_samples,double *gene_exprs_sorted,double *quantile_thresholds,double *quantiles);


void get_quantile_scores(double *exprs,int *n_genes,int *n_grp1,int *n_grp2,int *labels,int *n_labels,double *quantile_thresholds,int *scores) {

	int g = 0,i = 0,j = 0;
	int n_samples = *n_grp1 + *n_grp2;
	int n_quantiles = *n_labels - 1;
	
	double grp1_exprs[*n_grp1];
	double grp2_exprs[*n_grp2];

	double gene_exprs[n_samples];
	double gene_exprs_sorted[n_samples];

	double quantiles[n_quantiles];


	for(g=0;g<*n_genes;g++) {

		int gene_offset = g * n_samples;
		int NA_count = 0;


		for(i=0;i<*n_grp1;i++) {
			int NA_count_grp1 = 0;

			if(isnan(*(exprs + gene_offset + i))) {

				grp1_exprs[i] = NA_VALUE;
				gene_exprs[i] = NA_VALUE;	
				NA_count++;
				NA_count_grp1++;
			}
			else {
				grp1_exprs[i] = *(exprs + gene_offset + i);
				gene_exprs[i] = grp1_exprs[i];	
				gene_exprs_sorted[(i - NA_count_grp1)] = grp1_exprs[i];
			}
		}

		for(i=0;i<*n_grp2;i++) {
			int NA_count_grp2 = 0;

			if(isnan(*(exprs + gene_offset + *n_grp1 + i))) {

				grp2_exprs[i] = NA_VALUE;
				gene_exprs[(i+*n_grp1)] = NA_VALUE;
				NA_count++;
				NA_count_grp2++;
			}
			else {
				grp2_exprs[i] = *(exprs + gene_offset + *n_grp1 + i);
				gene_exprs[(i+*n_grp1)] = grp2_exprs[i];
				gene_exprs_sorted[(i + *n_grp1 - NA_count_grp2)] = grp2_exprs[i];
			}
		}

		quantile(&n_quantiles,(n_samples-NA_count),gene_exprs_sorted,quantile_thresholds,quantiles);


		for(i=0;i<*n_grp1;i++) {
			
			if(grp1_exprs[i] == NA_VALUE) {
				*(scores+i+gene_offset) = NA_VALUE;
			}
			else {

				for(j=0;j<n_quantiles;j++) {
					
					if( (j==0) && (grp1_exprs[i] < quantiles[j]) ) {
						*(scores+i+gene_offset) = *(labels);
					}
					else if ( (grp1_exprs[i] >= quantiles[(j-1)]  ) && (grp1_exprs[i] < quantiles[j]) ) {
						*(scores+i+gene_offset) = *(labels+j);
					}
					else if( (j==(n_quantiles-1)) && (grp1_exprs[i] >= quantiles[j]) ) {
						*(scores+i+gene_offset) = *(labels+n_quantiles);
					}
				}
			}
		}
		for(i=0;i<*n_grp2;i++) {
			
			if(grp2_exprs[i] == NA_VALUE) {
				*(scores+i+*n_grp1+gene_offset) = NA_VALUE;
			}
			else {
				for(j=0;j<n_quantiles;j++) {
					
					if( (j==0) && (grp2_exprs[i] < quantiles[j]) ) {
						*(scores+i+*n_grp1+gene_offset) = *(labels);
					}
					else if ( (grp2_exprs[i] >= quantiles[(j-1)]  ) && (grp2_exprs[i] < quantiles[j]) ) {
						*(scores+i+*n_grp1+gene_offset) = *(labels+j);
					}
					else if( (j==(n_quantiles-1)) && (grp2_exprs[i] >= quantiles[j]) ) {
						*(scores+i+*n_grp1+gene_offset) = *(labels+n_quantiles);
					}
				}
			}
		}
	}
}

void insertion_sort(double *a,int length) {

	int i;
	for(i=0;i<length;i++) {

		int j;
		double v=*(a+i);

		for(j= (i-1);j>=0;j--) {
			if(*(a+j) <= v) break;
			*(a+j+1) = *(a+j);
		}
		*(a+j+1) = v;
	}
}

void quantile(int *n_quantiles,int n_samples,double *gene_exprs_sorted,double *quantile_thresholds,double *quantiles) {

	int i=0,j=0;

	double index[*n_quantiles];
	
	double lo[*n_quantiles];
	
	double hi[*n_quantiles];
	
	double h[*n_quantiles];

	insertion_sort(gene_exprs_sorted,n_samples);

	for(i=0;i<*n_quantiles;i++) {
		
		index[i] = 1 + (n_samples - 1) * *(quantile_thresholds+i);
		
		lo[i] = floor(index[i]);
		
		hi[i] = ceil(index[i]);



		for(j=0;j<n_samples;j++) {
			
			if((j+1) == (int)lo[i]) {

				*(quantiles+i) = *(gene_exprs_sorted+j);
			}
		}

		h[i] = index[i] - lo[i];

		if(h[i] != 0) {

			*(quantiles+i) = (1 - h[i]) * *(quantiles+i) + h[i] * *(gene_exprs_sorted + ((int)hi[i]-1) );
		}
	}
}




