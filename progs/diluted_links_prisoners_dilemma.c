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

#define k 0.1

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
}

int main(int argc, char *argv[]){

	Lattice net;
	Memory mem;
	Params p;

	getParams(&p, argc, argv);
	printParams(p);

	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, p.seed);

	initLattice(&net, p.L);
	printf("Lattice initiated\n");
	initMemory(&mem, p);
	printf("Memory initiated\n");

	int count = 0;
	for(int i = 0; i < p.L*p.L; i++){
		printf("%d ", net.players[i]);
		count++;
		if(count % p.L == 0) printf("\n");
	}
	printf("\n");
	printf("\n");
	count = 0;
	for(int i = 0; i < p.L*p.L; i++) count += net.players[i];
	printf("Total = %d\n", count);

	printf("\n");
	printf("\n");

	count = 0;
	for(int i = 0; i < p.L*p.L; i++){
		printf("%d ", net.neigh[i][0]);
		count++;
		if(count % p.L == 0) printf("\n");
	}

	printf("\n");
	printf("\n");

	count = 0;
	for(int i = 0; i < p.L*p.L; i++){
		printf("%d ", net.links[i][0] + net.links[i][1] + net.links[i][2] + net.links[i][3]);
		count++;
		if(count % p.L == 0) printf("\n");
	}

	freeLattice(&net);
	freeMemory(&mem, p);

	return 0;
}
