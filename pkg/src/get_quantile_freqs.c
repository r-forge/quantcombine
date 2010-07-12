#include "quant_combine.h"

void get_quantile_freqs(int *scores,int *n_genes,int *n_grp1,int *n_grp2,int *labels,int *n_labels,int *freqs) {

	int g = 0,i = 0,j = 0;
	int n_samples = *n_grp1 + *n_grp2;
	int freqs_length = 2 * *n_labels;

	int grp1_scores[*n_grp1];
	int grp2_scores[*n_grp2];

	int grp1_freqs[*n_labels];
	int grp2_freqs[*n_labels];

	for(g=0;g<*n_genes;g++) {

		int gene_offset = g * n_samples;
		int freqs_offset = g * freqs_length;
		
		for(i=0;i<*n_grp1;i++) {

			grp1_scores[i] = *(scores + gene_offset + i);
		}

		for(i=0;i<*n_grp2;i++) {

			grp2_scores[i] = *(scores + gene_offset + *n_grp1 + i);
		}

		for(i=0;i<*n_labels;i++) {
			
			grp1_freqs[i] = 0;
			grp2_freqs[i] = 0;

			for(j=0;j<*n_grp1;j++) {

				if(*(labels+i) == grp1_scores[j]) {
					grp1_freqs[i] += 1;
				}
			}

			for(j=0;j<*n_grp2;j++) {
				
				if(*(labels+i) == grp2_scores[j]) {
					grp2_freqs[i] += 1;
				}
			}
			
			*(freqs + freqs_offset + i) = grp1_freqs[i];
			
			*(freqs + freqs_offset + i + *n_labels) = grp2_freqs[i];
		}
	}
}
