#ifndef _PERCOLATION_H_
#define _PERCOLATION_H_

typedef struct{
	unsigned int L;
	unsigned int L2;
	int *clusters;
	char *histname;
	char *percolname;
	int biggest_labels[2];
} HKstuff;

void initHK(HKstuff *hk, int L, char *histname, char *percolname);

void freeHK(HKstuff *hk);

void hk(Lattice grid, HKstuff *hks, int opt);

void bond_hk(Lattice grid, HKstuff *hks);

int does_it_percolate(int *clusters, int L, int label);

void percolation_measures(HKstuff *hks);

#endif
