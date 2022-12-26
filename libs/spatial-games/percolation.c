#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spatial_games.h"
#include "percolation.h"

int uf_find(int *labels, int x){
	// finds x equivalence class' label
	int z;
	int y = x;

	// encontra a raiz
	while(labels[y] != y){
		y = labels[y];
	}

	// colapsa a arvore no caminho
	// pra todos do caminho apontarem
	// para a raiz
	while(labels[x] != x){
		z = labels[x];
		labels[x] = y;
		x = z;
	}

	return y;
}

void uf_union(int *labels, int x, int y){
	// aponta a maior label para a menor
	int a, b;

	a = uf_find(labels, x);
	b = uf_find(labels, y);

	// pego o maior
	if(a > b){
		labels[a] = labels[b];
	}
	else{
		labels[b] = labels[a];
	}
}

/*********** USER FUNCTIONS ************/

void initHK(HKstuff *hk, int L, char *histname, char *percolname){
	int L2 = L*L;
	hk->L = L;
	hk->L2 = L2;
	hk->clusters = (int *) malloc(hk->L2 * sizeof(int));
	hk->histname = histname;
	hk->percolname = percolname;
	hk->biggest_labels[0] = -1;
	hk->biggest_labels[1] = -1;
}

void freeHK(HKstuff *hk){
	free(hk->clusters);
}

void hk(Lattice grid, HKstuff *hks, int opt){
	int *labels;
	int i;
	int L = grid.L;
	int L2 = L * L;
	int n2, n3;	// direita e abaixo

	labels = (int *) calloc(L2, sizeof(int));

	// pigeon hole
	for(i = 0; i < L2; i++){
		labels[i] = i;
	}

	if(opt == 1 || opt == 0){
		for(i = 0; i < L2; i++){
			n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
			n3 = ((i/L + 1 + L2) % L)*L + (i % L);
			if(grid.players[i] == grid.players[n2] && grid.players[i] == opt) uf_union(labels, i, n2);
			if(grid.players[i] == grid.players[n3] && grid.players[i] == opt) uf_union(labels, i, n3);
		}
		// escrevo a rede etiquetada
		for(i = 0; i < L2; i++){
			if(grid.players[i] == opt){
				// pego a raiz de cada coordenada
				hks->clusters[i] = uf_find(labels, i);
			}
			else{
				hks->clusters[i] = -1;
			}
		}
	}
	else if(opt == 2){
		for(i = 0; i < L2; i++){
			n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
			n3 = ((i/L + 1 + L2) % L)*L + (i % L);
			if(grid.players[i] == grid.players[n2]) uf_union(labels, i, n2);
			if(grid.players[i] == grid.players[n3]) uf_union(labels, i, n3);
		}

		for(i = 0; i < L2; i++){
			hks->clusters[i] = uf_find(labels, i);
		}
	}

	free(labels);

}

void bond_hk(Lattice grid, HKstuff *hks){
	int *labels;
	int i;
	int L = grid.L;
	int L2 = L*L;
	int n2, n3;

	labels = (int *) calloc(L2, sizeof(int));

	for(i = 0; i < L2; i++) labels[i] = i;
	
	for(i = 0; i < L2; i++){
		n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
		n3 = ((i/L + 1 + L2) % L)*L + (i % L);
		if(grid.links[i][1] == 1) uf_union(labels, i, n2);
		if(grid.links[i][2] == 1) uf_union(labels, i, n3);
	}

	for(i = 0; i < L2; i++) hks->clusters[i] = uf_find(labels, i);

	free(labels);
}

int does_it_percolate(int *clusters, int L, int label){
	int i, j, ok1, ok2, res;

	for(i = 0; i < L; i++){
		ok1 = 0;
		for(j = i; j < L*L; j += L){
			if(clusters[j] == label){
				ok1 = 1;
				break;
			}
		}
		if(ok1 == 0) break;
	}

	for(i = 0; i < L*L; i += L){
		ok2 = 0;
		for(j = i; j < i + L; j++){
			if(clusters[j] == label){
				ok2 = 1;
				break;
			}
		}
		if(ok2 == 0) break;
	}

	res = ok1 + ok2;
	return res;
}

void percolation_measures(HKstuff *hks){
	int i;
	int L = hks->L;
	int L2 = L*L;
	int ok = 0;
	int percol = 0;
	int biggest, biggest_label;
	int second_biggest, second_biggest_label;
	int cluster_count;
	int *size, *hist;
	FILE *arq;

	size = (int *) calloc(L2, sizeof(int));
	hist = (int *) calloc(L2 + 1, sizeof(int));

	// sizes
	for(i = 0; i < L2; i++){
		if(hks->clusters[i] != -1){
			size[hks->clusters[i]] += 1;
		}
	}

	biggest = 0;
	biggest_label = 0;
	second_biggest = 0;
	second_biggest_label = 0;
	cluster_count = 0;
	arq = fopen(hks->percolname, "w");
	for(i = 0; i < L2; i++){
		if(size[i] > 0){
			hist[size[i]] += 1;	// distrib
			if(size[i] >= biggest){
				second_biggest = biggest;
				second_biggest_label = biggest_label;
				biggest = size[i];
				biggest_label = i;
			} else if(size[i] >= second_biggest){
				second_biggest = size[i];
				second_biggest_label = i;
			}
			if(size[i] >= L){
				percol = does_it_percolate(hks->clusters, L, i);
			}
			cluster_count++;
		}

		if(percol != 0 && ok == 0){
			fprintf(arq, "Percolation status: %d\n", percol);
			ok = 1;
		}
	}
	if(ok == 0){
		fprintf(arq, "Percolation status: %d\n", 0);
	}
	fprintf(arq, "Biggest cluster label: %d\n", biggest_label);
	fprintf(arq, "Biggest cluster size: %d\n", biggest);
	fprintf(arq, "Second biggest cluster label: %d\n", second_biggest_label);
	fprintf(arq, "Second biggest cluster size: %d\n", second_biggest);
	fprintf(arq, "Cluster count = %d\n", cluster_count);
	fclose(arq);
	
	arq = fopen(hks->histname, "w");
	for(i = 0; i < L2 + 1; i++){
		if(hist[i] != 0)
			fprintf(arq, "%d %d\n", i, hist[i]);
	}
	fclose(arq);

	free(size);
	free(hist);

	hks->biggest_labels[0] = biggest_label;
	hks->biggest_labels[1] = second_biggest_label;

}

