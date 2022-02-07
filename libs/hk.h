#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////// Function headers //////////
int uf_find(int x, int *labels);
void uf_union(int x, int y, int *labels);
void hk(int *hk_grid, int *grid, int L, int opt);
void bond_hk(int *hk_grid, int *grid, int *links, int L);
void percolation_measures(int *clusters, int L, char *histname, char *percolname, int *bl);
int does_it_percolate(int *clusters, int L, int label);
//////////////////////////////////////

//////////// Declarations ////////////
int uf_find(int x, int *labels){
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

void uf_union(int x, int y, int *labels){
	// aponta a maior label para a menor
	int a, b;

	a = uf_find(x, labels);
	b = uf_find(y, labels);

	// pego o maior
	if(a > b){
		labels[a] = labels[b];
	}
	else{
		labels[b] = labels[a];
	}
}

void hk(int *hk_grid, int *grid, int L, int opt){
	int *labels;
	int i;
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
			if(grid[i] == grid[n2] && grid[i] == opt) uf_union(i, n2, labels);
			if(grid[i] == grid[n3] && grid[i] == opt) uf_union(i, n3, labels);
		}

		// escrevo a rede etiquetada
		for(i = 0; i < L2; i++){
			if(grid[i] == opt){
				// pego a raiz de cada coordenada
				hk_grid[i] = uf_find(i, labels);
			}
			else{
				hk_grid[i] = -1;
			}
		}
	}
	else if(opt == 2){
		for(i = 0; i < L2; i++){
			n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
			n3 = ((i/L + 1 + L2) % L)*L + (i % L);
			if(grid[i] == grid[n2]) uf_union(i, n2, labels);
			if(grid[i] == grid[n3]) uf_union(i, n3, labels);
		}

		for(i = 0; i < L2; i++){
			hk_grid[i] = uf_find(i, labels);
		}
	}

	free(labels);

}

void percolation_measures(int *clusters, int L, char *histname, char *percolname, int *bl){
	int i, L2 = L*L;
	int ok = 0;
	int percol = 0;
	int biggest, biggest_label, second_biggest, second_biggest_label, cluster_count;
	int *size, *hist;
	FILE *arq;

	size = (int *) calloc(L2, sizeof(int));
	hist = (int *) calloc(L2 + 1, sizeof(int));

	// sizes
	for(i = 0; i < L2; i++){
		if(clusters[i] != -1){
			size[clusters[i]] += 1;
		}
	}

	biggest = 0;
	biggest_label = 0;
	second_biggest = 0;
	second_biggest_label = 0;
	cluster_count = 0;
	arq = fopen(percolname, "w");
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
				percol = does_it_percolate(clusters, L, i);
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
	
	arq = fopen(histname, "w");
	for(i = 0; i < L2 + 1; i++){
		if(hist[i] != 0)
			fprintf(arq, "%d %d\n", i, hist[i]);
	}
	fclose(arq);

	free(size);
	free(hist);

	bl[0] = biggest_label;
	bl[1] = second_biggest_label;

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

void bond_hk(int *hk_grid, int *grid, int *links, int L){
	int *labels;
	int i;
	int L2 = L*L;
	int n2, n3;

	labels = (int *) calloc(L2, sizeof(int));

	for(i = 0; i < L2; i++) labels[i] = i;
	
	for(i = 0; i < L2; i++){
		n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
		n3 = ((i/L + 1 + L2) % L)*L + (i % L);
		if(links[4*i + 1] == 1) uf_union(i, n2, labels);
		if(links[4*i + 2] == 1) uf_union(i, n3, labels);
	}

	for(i = 0; i < L2; i++) hk_grid[i] = uf_find(i, labels);

	free(labels);
}

