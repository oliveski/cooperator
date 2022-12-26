/*
 * gcc -Wall simulacao.c spatial_games.c -o exe -lm
 * -lgsl -lgslcblas
 * ./exe L r semente tau M crit
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <lat2eps.h>
#include <gsl/gsl_rng.h>
#include "spatial_games.h"
#include "percolation.h"
#include "lattice.h"

#define SNAP 1

extern gsl_rng *r;

void printParams(Params p){
	printf("L = %d\n", p.L);
	printf("L2 = %d\n", p.L2);
	printf("r = %lf\n", p.r);
	printf("seed = %d\n", p.seed);
	printf("tau = %lf\n", p.tau);
	printf("M = %d\n", p.M);
	printf("crit = %lf\n", p.crit);
	printf("| %lf\t%lf |\n", p.payoffs[0][0], p.payoffs[0][1]);
	printf("| %lf\t%lf |\n", p.payoffs[1][0], p.payoffs[1][1]);
}

int main(int argc, char *argv[]){

	Lattice net;
	Memory mem;
	Params p;

	getParams(&p, argc, argv);
	/* printParams(p); */

	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, p.seed);

	initLattice(&net, p.L);
	initMemory(&mem, net, p);

	///////// debug /////////
	/* printNet(net); */
	/* printMemory(mem, p); */
	////////////////////////

	int t;	// time
	char series_filename[1000];
	sprintf(series_filename, "r%.4lf_L%d_S%d_M%d_tau%.4f.dat", p.r, p.L, p.seed, p.M, p.tau);
	FILE *series = fopen(series_filename, "w");

	// Simulation loop //
	for(int mcs = 1; mcs <= MCS; mcs++){
		// index of history
		t = mcs % p.M;
		diluteLinks(&net, mem, p, t);
		for(int i = 0; i < p.L2; i++)
			playGame(&net, &mem, p, t);
		// if we should measure, we measure
		if(mcs >= TRAN && mcs % WSIZE == 0){
			// Test for extinction
			if(net.coop_count == p.L2 || net.coop_count == 0){
				dumpExtinct(net, series, mcs);
				fclose(series);
				freeLattice(&net);
				freeMemory(&mem, p);
				gsl_rng_free(r);
				return 0;
			}
			fprintf(series, "%d %d %d\n", mcs-TRAN, net.coop_count, net.link_count);
		}
	}
	fclose(series);

	// Percolation measures //
	HKstuff hks;
	char histname[100], percolname[100];
	sprintf(histname, "clstr_hstgrm_r%.4lf_L%d_S%d_M%d_tau%.4f.dat",
			p.r, p.L, p.seed, p.M, p.tau);
	sprintf(percolname, "prcltn_1_r%.4lf_L%d_S%d_M%d_tau%.4f.dat",
			p.r, p.L, p.seed, p.M, p.tau);
	// state 1
	initHK(&hks, p.L, histname, percolname);
	hk(net, &hks, 1);
	percolation_measures(&hks);
	freeHK(&hks);

	// bond
	sprintf(histname, "bond_clstr_hstgrm_r%.4lf_L%d_S%d_M%d_tau%.4f.dat",
			p.r, p.L, p.seed, p.M, p.tau);
	sprintf(percolname, "bond_prcltn_r%.4lf_L%d_S%d_M%d_tau%.4f.dat",
			p.r, p.L, p.seed, p.M, p.tau);
	initHK(&hks, p.L, histname, percolname);
	bond_hk(net, &hks);
	percolation_measures(&hks);
	freeHK(&hks);

	// snapshots
	int *lpp = linksPerPlayer(net);
	char snpsht[100];
	sprintf(snpsht, "snpsht_r%.4lf_L%d_S%d_M%d_tau%.4f.eps",
			p.r, p.L, p.seed, p.M, p.tau);

	openlat(p.L, 4);
	printlat(net, hks, lpp, snpsht, 4);


	///////// debug ////////
	/* printNet(net); */
	/* printMemory(mem, p); */
	////////////////////////

	freeLattice(&net);
	freeMemory(&mem, p);
	gsl_rng_free(r);

	return 0;
}
