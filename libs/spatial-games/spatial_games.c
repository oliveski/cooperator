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
	p->payoffs[1][1] = 1;		// R
	p->payoffs[1][0] = - p->r;	// S
	p->payoffs[0][1] = 1 + p->r;	// T
	p->payoffs[0][0] = 0;		// P
	printf("R = %lf\n", p->payoffs[1][1]);
	printf("S = %lf\n", p->payoffs[1][0]);
	printf("T = %lf\n", p->payoffs[0][1]);
	printf("P = %lf\n", p->payoffs[0][0]);
}

double W(double Px, double Py){
	double diff_Py_Px = Py - Px;
	double expr;

	expr = 1 + exp(-(diff_Py_Px)/k);
	return (1.0 / expr);
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

	net->coop_count = 0;
	net->link_count = 0;
	int n1, n2, n3, n4;
	for(int i = 0; i < L2; i++){
		net->players[i] = gsl_rng_uniform(r) * 2;
		net->coop_count += net->players[i];
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
		net->link_count += 4;
	}

}

void updateCounts(Lattice *net){
	net->coop_count = 0;
	net->link_count = 0;
	
	for(int i = 0; i < net->L2; i++){
		net->coop_count += net->players[i];
		net->link_count += net->links[i][0];
		net->link_count += net->links[i][1];
		net->link_count += net->links[i][2];
		net->link_count += net->links[i][3];
	}
}

void freeLattice(Lattice *net){
	free(net->players);
	free(net->neigh);
	free(net->links);
}

void initMemory(Memory *mem, Lattice net, Params p){
	mem->hist = (int **) malloc(p.L2 * sizeof(int *));
	mem->n = (int *) malloc(p.L2 * sizeof(int));

	for(int i = 0; i < p.L2; i++){
		mem->hist[i] = (int *) malloc(p.M * sizeof(int));
	}

	for(int i = 0; i < p.L2; i++){
		mem->n[i] = int_powerlawRng(p);
#ifdef DELTA
		mem->n[i] = p.M;
#endif
		for(int j = 1; j < p.M; j++){
			// at t=0 my history is equal to my strategy
			// and t>0 is initiated as if I had cooperated
			mem->hist[i][0] = net.players[i];
			mem->hist[i][j] = 1;
		}
	}
}

void freeMemory(Memory *mem, Params p){
	free(mem->n);
	for(int i = 0; i < p.L2; i++) free(mem->hist[i]);
	free(mem->hist);
}

void linkAvaliation(Lattice *net, Memory mem, Params p, int t, int player){
	double frac;	// fraction of cooperation in the history
	int neigh;	// neighbours coordinates

	// if the fraction of cooperation in the neighbours history
	// is less or equal to crit we cut ties
	for(int i = 0; i < 4; i++){
		if(net->links[player][i]){
			frac = 0;
			neigh = net->neigh[player][i];
			for(int j = 0; j < mem.n[player]; j++)
				frac += (double) mem.hist[neigh][(t - j + p.M) % p.M] / mem.n[player];
			if(frac <= p.crit){
				net->links[player][i] = 0;	// i -> j
				net->links[neigh][(i + 2)%4] = 0; // j -> i
				net->link_count -= 2;
			}
		}
	}
}

double payoffCalculation(Lattice net, Params p, int player){
	int player_state = net.players[player];
	int link;
	int neigh_state;
	double payoff;
	double pay;

	payoff = 0;
	/* printf("Payoff para %d\n", player); */
	for(int i = 0; i < 4; i++){
		/* printf("viz %d: ", i); */
		link = net.links[player][i];
		/* printf("link=%d ", link); */
		neigh_state = net.players[net.neigh[player][i]];
		/* printf("nstat=%d ", neigh_state); */
		pay = p.payoffs[player_state][neigh_state];
		/* printf("pay=%lf\n", pay); */
		payoff += link * pay;
		/* printf("\n"); */
	}
	/* printf("payoff de %d: %lf\n", player, payoff); */

	return payoff;
}

void strategyUpdate(Lattice *net, double Px, double Py, int player, int neigh){
	double prob = W(Px, Py);
	double random = gsl_rng_uniform(r);

	if(random < prob){
		net->players[player] = net->players[net->neigh[player][neigh]];
		/* printf("[%d] trocou\n", player); */
	}
}

void printNet(Lattice net){
	for(int i = 0; i < net.L; i++){
		for(int j = 0; j < net.L; j++){
			printf("%d ", net.players[i*net.L + j]);
		}
		printf("\n");
	}
	printf("coop_count = %d\n", net.coop_count);
	printf("link_count = %d\n", net.link_count);
}

void printMemory(Memory mem, Params p){
	for(int i = 0; i < p.L2; i++){
		printf("[%d]:", i);
		for(int j = 0; j < p.M; j++){
			printf(" %d", mem.hist[i][j]);
		}
		printf("\n");
	}
}

void playGame(Lattice *net, Memory *mem, Params p, int t){
	/* printf("MCS %d\n", t); */
	int player = (int) gsl_rng_uniform_int(r, p.L2);
	int init_state = net->players[player];
	/* printf("Sorteado %d (%d)\n", player, init_state); */
	int neigh = (int) gsl_rng_uniform_int(r, 4);
	/* printf("vizinho %d\n", neigh); */

	double payoff_i = payoffCalculation(*net, p, player);
	double payoff_j = payoffCalculation(*net, p, net->neigh[player][neigh]);
	/* printf("Pi = %.4lf; Pj = %.4lf\n", payoff_i, payoff_j); */

	strategyUpdate(net, payoff_i, payoff_j, player, neigh);

	mem->hist[player][t] = net->players[player];
	net->coop_count += net->players[player] - init_state;
	
	/* printf("#### %d ####\n", t); */
	/* printNet(*net); */
	/* printf("###########\n"); */
}

void diluteLinks(Lattice *net, Memory mem, Params p, int t){
	for(int i = 0; i < p.L2; i++){
		for(int j = 0; j < 4; j++){
			net->links[i][j] = 1;
		}
	}

	net->link_count = 4*p.L2;

	for(int i = 0; i < p.L2; i++)
		linkAvaliation(net, mem, p, t, i);

}

int *linksPerPlayer(Lattice net){
	int L2 = net.L2;
	int *lpp = (int *) malloc(L2 * sizeof(int));
	for(int i = 0; i < L2; i++) lpp[i] = 0;

	for(int i = 0; i < L2; i++){
		for(int j = 0; j < 4; j++){
			lpp[i] += net.links[i][j];
		}
	}
	return lpp;
}

void dumpExtinct(Lattice net, FILE *series, int mcs){
	while(mcs < MCS){
		fprintf(series, "%d %d %d\n", mcs-TRAN, net.coop_count, net.link_count);
		mcs += 50;
	}
}
