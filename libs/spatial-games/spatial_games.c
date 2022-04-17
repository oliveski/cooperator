#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "spatial_games.h"

gsl_rng *r;

void getParams(Params *p, int argc, char *argv[]){
	p->L = atoi(argv[1]);
	p->L2 = p->L * p->L;
	p->r = atof(argv[2]);
	p->seed = atoi(argv[3]);
	p->tau = atof(argv[4]);
	p->M = atoi(argv[5]);
	p->crit = atof(argv[6]);
}

double powerlawRng(Params p){
	double x = gsl_rng_uniform(r);
	double rn;

	rn = x*pow(p.M + 1, 1 - p.tau) - x + 1;
	rn = pow(rn, 1.0/(1 - p.tau));
	
	return rn;
}

int int_powerlawRng(Params p){
	double rn = powerlawRng(p);
	rn = (int) rn;

	return rn;
}

void initLattice(Lattice *net, int L){

	int L2 = L*L;
	net->L = L;
	net->L2 = L2;
	net->players = (int *) malloc(net->L2 * sizeof(int));
	net->neigh = malloc(4 * net->L2 * sizeof(int));
	net->links = malloc(4 * net->L2 * sizeof(int));

	int n1, n2, n3, n4;
	for(int i = 0; i < L2; i++){
		net->players[i] = gsl_rng_uniform(r) * 2;
		n1 = ((i/L - 1 + L2) % L)*L + (i % L);
		n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
		n3 = ((i/L + 1 + L2) % L)*L + (i % L);
		n4 = ((i/L) % L)*L + ((i - 1 + L) % L);

		net->neigh[i][0] = n1;
		net->neigh[i][1] = n2;
		net->neigh[i][2] = n3;
		net->neigh[i][3] = n4;

		net->links[i][0] = 1;
		net->links[i][1] = 1;
		net->links[i][2] = 1;
		net->links[i][3] = 1;
	}

}

void freeLattice(Lattice *net, int L){
	int L2 = L*L;
	free(net->players);
	free(net->neigh);
	free(net->links);
}
