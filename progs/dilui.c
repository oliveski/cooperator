/************ compilar **************/
// gcc -Wall dilui.c -o exe -lm -llat2eps//
/************ rodar ***********/
// ./exe L r semente tau M crit //
/*****************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <lat2eps.h>
#include "pointers.h"
#include "mc.h"
#include "hk.h"
#include "jlat.h"

#define k 0.1
#define TRAN 100000
#define MCS 500000
#define WSIZE 50
#define c 1

/************* SWITCHES ***********/
#define INICIAL 1
#define SNAP	1
/**********************************/

/********* VARIAVEIS GLOBAIS ******/
int L;			// Tamanho da rede
int L2;
int M;			// Memória máxima
double tau;		// expoente da distribuição
double crit;		// frac mínima para cooperar
/**********************************/

/*********** DECLARAÇÕES *************/

void inicializa(int *jogadores, int *vizinhos, int *links,
		int *n, int *hist, int *count);

double calcula_payoff(int *jogadores, int *vizinhos, int *links,
		double payoffs[2][2], int jogador);

double teste_calcula_payoff(int *jogadores, int *vizinhos, int *links,
		double payoffs[2][2], int jogador);

double W(double Px, double Py);

void troca_estrategia(int *jogadores, int *vizinhos,
		double Px, double Py, int jogador, int vizinho);

void avalia_links(int *vizinhos, int *links, int *n, int *hist, int *count_links, int t, int jogador);

void conta_coop(int *jogadores, int *count);

void conta_links(int *links, int *count_links);

void atualiza_coop(int *count, int jogador_antigo, int jogador_novo, 
		int atcount[2][2]);

void jogo_tran(int *jogadores, int *vizinhos, int *links, int *hist,
		int t, double payoffs[2][2]);

void jogo_mcs(int *jogadores, int *vizinhos, int *links, int *hist,
		int t, double payoffs[2][2], int *count, int atcount[2][2]);

void escreve_sistema(int *jogadores, int *vizinhos, int *links);

void diluir(int *vizinhos, int *links, int *n, int *hist, int *count_links, int t);

double power_law(double x);

void escreve_memorias(int *hist);

double rng();

int int_rng();

void dist_n(int *n);

int transientcheck(int *jogadores, int *vizinhos, int *links, int *hist,
		int *n, FILE *series, int count, int count_links);

int measurecheck(int *jogadores, int *vizinhos, int *links, int *hist,
		int *n, FILE *series, int count, int count_links, int measure);

void gnuplot_view(int *jogadores, FILE *gnu_p);

int *lpjogador(int *links);

#ifdef SNAP
void old_openlat(int opt);

void old_printlat(int *lat_pt, char *snpsht, int opt);
#endif

/**************************************/

int main(int argc, char *argv[]){

	int count;		// conta cooperadores
	int count_links;	// conta links
	int *jogadores;		// matriz de jogadores
	int *vizinhos;		// topologia
	int *links;		// quais vizinhos estão ligados
	int *n;			// n de cada jogador
	int *hist;		// historico dos jogadores
	double payoffs[2][2];	// matriz de payoffs
	int atcount[2][2];	// atualiza count
#ifdef SNAP
	int *bl;		// biggest clusters
	int *lpj;		// links per site
	int *clusters0;
	int *clusters1;
	int *clusters2;
	int *bonds;
#endif
	int L_arg;
	double tau_arg;
	double crit_arg;
	int M_arg;
	int semente;
	int mcs, measure;
	int i, t;
	double R, T, S, P, r;	// payoffs normalizados pra R
	FILE *series;		// serie temporal
	FILE *gnu_p;		// pipe pro gnuplot
	char nome_arq[1000];

	L_arg = atoi(argv[1]);
	L = L_arg;		// recebe o tamanho da rede
	L2 = L*L;
	
	r = atof(argv[2]);	// recebe r

	semente = atoi(argv[3]);
	srand(semente);		// inicia o gerador
	start_randomic(semente); // inicia FRANDOM

	tau_arg = atof(argv[4]);
	tau = tau_arg;
	M_arg = atoi(argv[5]);
	M = M_arg;
	crit_arg = atof(argv[6]);
	crit = crit_arg;	// recebe o crit

	R = 1;
	T = 1 + r;
	S = -r;
	P = 0;
	/* payoffs[0] = R; */
	/* payoffs[1] = S; */
	/* payoffs[2] = T; */
	/* payoffs[3] = P; */

	payoffs[1][1]= R;
	payoffs[1][0]= S;
	payoffs[0][1]= T;
	payoffs[0][0]= P;

	atcount[0][0] = 0;
	atcount[0][1] = 1;
	atcount[1][0] = -1;
	atcount[1][1] = 0;

	jogadores = create_int_pointer(L2);
	vizinhos = create_int_pointer(4*L2);	// topologia de vizinhos
	links = create_int_pointer(4*L2);	// Jogadores conectados
	n = create_int_pointer(L2);		// extensão da memória
	hist = create_int_pointer(M*L2);	// histórico de jogadas

#ifdef SNAP
	// hk_grids for hk() and bond_hk()
	clusters0 = create_int_pointer(L2);
	clusters1 = create_int_pointer(L2);
	clusters2 = create_int_pointer(L2);
	bonds = create_int_pointer(L2);
	bl = create_int_pointer(2);	// biggest labels
#endif

	inicializa(jogadores, vizinhos, links, n, hist, &count);

#ifdef VIEW
	gnu_p = popen("gnuplot --persist", "w");
	fprintf(gnu_p, "set terminal wxt enhanced background rgb 'black'\n");
	gnuplot_view(jogadores, gnu_p);
#endif

	//Transiente
	for(mcs = 1; mcs <= TRAN; mcs++){
		t = mcs % M;
		diluir(vizinhos, links, n, hist, &count_links, t);
		for(i = 0; i < L2; i++){
			jogo_tran(jogadores, vizinhos, links, hist, t, payoffs);
		}
	}

	// config de teste
	/* for(i = 0; i < L2; i++){ */
	/* 	jogadores[i] = 0; */
	/* } */
	
	/* jogadores[0] = 1; */
	/* jogadores[2] = 1; */
	/* jogadores[3] = 1; */
	/* jogadores[3 + L] = 1; */
	/* jogadores[5] = 1; */
	/* jogadores[6] = 1; */
	/* jogadores[6 + L] = 1; */
	/* jogadores[96] = 1; */
	/* jogadores[93] = 1; */
	/* jogadores[44] = 1; */
	/* jogadores[34] = 1; */
	/* jogadores[54] = 1; */
	/* jogadores[45] = 1; */
	/* jogadores[43] = 1; */

	conta_coop(jogadores, &count);
	printf("contagem após transiente = %d\n", count);

	conta_links(links, &count_links);
	printf("contagem dos links após o transiente = %d\n", count_links);

	sprintf(nome_arq, "r%.4lf_L%d_S%d_M%d_tau%.4f.dat", r, L, semente, M, tau);
	series = fopen(nome_arq, "w");

	//Testa se houve extinção
	if(transientcheck(jogadores, vizinhos, links, hist,
				n, series, count, count_links)){
		return 0;
	}
	
#ifdef VIEW
	gnuplot_view(jogadores, gnu_p);
#endif

	//Medidas
	for(measure = 0; measure < MCS/WSIZE; measure++){
		for(mcs = t; mcs < WSIZE + t; mcs++){
			t = mcs % M;
			diluir(vizinhos, links, n, hist, &count_links, t);
			for(i = 0; i < L2; i++){
				jogo_mcs(jogadores, vizinhos, links, hist, t,
						payoffs, &count, atcount);
			}
		#ifdef VIEW
			gnuplot_view(jogadores, gnu_p);
		#endif
		}
		if(measurecheck(jogadores, vizinhos, links, hist,
					n, series, count, count_links, measure)){
		#ifdef SNAP
			openlat(L, 4);
			char snpsht[1000];
			sprintf(snpsht, "snapshot_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			lpj = lpjogador(links);
			printlat(jogadores, snpsht, L, 4, bl, lpj);

			openlat(L, 2);
			sprintf(snpsht, "links_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			printlat(links, snpsht, L, 2, bl, lpj);

			hk(clusters0, jogadores, L, 0);
			hk(clusters1, jogadores, L, 1);
			hk(clusters2, jogadores, L, 2);
			char histname[1000];
			char percolname[1000];

			openlat(L, 3);
			sprintf(snpsht, "clusters0_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			sprintf(histname, "hist0s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			sprintf(percolname, "per0s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			percolation_measures(clusters0, L, histname, percolname, bl);
			printlat(clusters0, snpsht, L, 3, bl, lpj);

			openlat(L, 3);
			sprintf(snpsht, "clusters1_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			sprintf(histname, "hist1s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			sprintf(percolname, "per1s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			percolation_measures(clusters1, L, histname, percolname, bl);
			printlat(clusters1, snpsht, L, 3, bl, lpj);

			openlat(L, 3);
			sprintf(snpsht, "clusters2_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			sprintf(histname, "hist2s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			sprintf(percolname, "per2s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			percolation_measures(clusters2, L, histname, percolname, bl);
			printlat(clusters2, snpsht, L, 3, bl, lpj);

			openlat(L, 3);
			sprintf(snpsht, "bonds_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
			sprintf(histname, "histbonds_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			sprintf(percolname, "perbonds_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
			bond_hk(bonds, jogadores, links, L);
			percolation_measures(bonds, L, histname, percolname, bl);
			printlat(bonds, snpsht, L, 3, bl, lpj);

			free(jogadores);
			free(links);
			free(clusters0);
			free(clusters1);
			free(clusters2);
			free(bl);
			free(lpj);
		#endif
			return 0;
		}
		fprintf(series, "%d %d %d\n", measure*WSIZE, count, count_links);
	}

	fclose(series);

#ifdef VIEW
	fclose(gnu_p);
#endif

#ifdef SNAP
	openlat(L, 4);
	char snpsht[1000];
	sprintf(snpsht, "snapshot_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	lpj = lpjogador(links);
	printlat(jogadores, snpsht, L, 4, bl, lpj);

	openlat(L, 2);
	sprintf(snpsht, "links_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	printlat(links, snpsht, L, 2, bl, lpj);

	hk(clusters0, jogadores, L, 0);
	hk(clusters1, jogadores, L, 1);
	hk(clusters2, jogadores, L, 2);
	char histname[1000];
	char percolname[1000];

	openlat(L, 3);
	sprintf(snpsht, "clusters0_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	sprintf(histname, "hist0s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	sprintf(percolname, "per0s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	percolation_measures(clusters0, L, histname, percolname, bl);
	printlat(clusters0, snpsht, L, 3, bl, lpj);

	openlat(L, 3);
	sprintf(snpsht, "clusters1_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	sprintf(histname, "hist1s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	sprintf(percolname, "per1s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	percolation_measures(clusters1, L, histname, percolname, bl);
	printlat(clusters1, snpsht, L, 3, bl, lpj);

	openlat(L, 3);
	sprintf(snpsht, "clusters2_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	sprintf(histname, "hist2s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	sprintf(percolname, "per2s_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	percolation_measures(clusters2, L, histname, percolname, bl);
	printlat(clusters2, snpsht, L, 3, bl, lpj);

	openlat(L, 3);
	sprintf(snpsht, "bonds_r%.4lf_L%d_S%d_M%d.eps", r, L, semente, M);
	sprintf(histname, "histbonds_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	sprintf(percolname, "perbonds_r%.4lf_L%d_S%d_M%d.dsf", r, L, semente, M);
	bond_hk(bonds, jogadores, links, L);
	percolation_measures(bonds, L, histname, percolname, bl);
	printlat(bonds, snpsht, L, 3, bl, lpj);

	free(clusters0);
	free(clusters1);
	free(clusters2);
	free(bl);
	free(lpj);

#endif

#ifdef VIEW
	fclose(gnu_p);
#endif

	/*** LIBERA POINTERS ***/
	free(jogadores);
	free(vizinhos);
	free(links);
	free(n);
	free(hist);

	return 0;

}

/***************** FUNÇÕES **********************/

void inicializa(int *jogadores, int *vizinhos, int *links,
		int *n, int *hist, int *count){

	int i, j, h, l;
	int n1, n2, n3, n4;

	*count = 0;
	j = 0;
	l = 0;
	for(i = 0; i < L2; i++){
		jogadores[i] = FRANDOM * 2;
		n[i] = int_rng();
		n1 = ((i/L - 1 + L2) % L)*L + (i % L);
		n2 = ((i/L) % L)*L + ((i + 1 + L) % L);
		n3 = ((i/L + 1 + L2) % L)*L + (i % L);
		n4 = ((i/L) % L)*L + ((i - 1 + L) % L);
		vizinhos[j + i + 0] = n1;
		vizinhos[j + i + 1] = n2;
		vizinhos[j + i + 2] = n3;
		vizinhos[j + i + 3] = n4;
		links[j + i + 0] = 1;
		links[j + i + 1] = 1;
		links[j + i + 2] = 1;
		links[j + i + 3] = 1;

		// como inicializar o estado da memoria
#ifdef UM
		hist[l + i + 0] = jogadores[i];
		for(h = 1; h < M; h++){
			hist[l + i + h] = 1;
		}
#endif
#ifdef INICIAL
		for(h = 0; h < M; h++){
			hist[l + i + h] = jogadores[i];
		}
#endif
 
		j = j + 3;
		l = l + M - 1;

		*count = *count + jogadores[i];
	}

}

/* double calcula_payoff(int *jogadores, int *vizinhos, int *links, */
/* 		double *payoffs[4], int jogador){ */
/* 	int i; */
/* 	int Nc = 0; */
/* 	int Nd = 0; */
/* 	int Z = 0; */
/* 	int link; */
/* 	int neigh_state; */
/* 	double payoff; */

/* 	payoff = 0; */

/* 	for(i = 0; i < 4; i++){ */
/* 		link = links[4*jogador + i]; */
/* 		neigh_state = jogadores[vizinhos[4*jogador + i]]; */
/* 		Nc = Nc + link*neigh_state; */
/* 		Z += link;	// number of neighs in the game */
/* 	} */

/* 	Nd = Z - Nc; */

/* 	if(jogadores[jogador] == 1){	// C */
/* 		payoff = Nc*payoffs[0] + Nd*payoffs[1]; */
/* 	} */
/* 	else{				// D */
/* 		payoff = Nc*payoffs[2] + Nd*payoffs[3]; */
/* 	} */

/* 	return payoff; */
/* } */

double teste_calcula_payoff(int *jogadores, int *vizinhos, int *links,
		double payoffs[2][2], int jogador){
	int i;
	int link;
	int player_state;
	int neigh_state;
	double payoff;
	double pay;

	payoff = 0;

	for(i = 0; i < 4; i++){
		link = links[4*jogador + i];
		player_state = jogadores[jogador];
		neigh_state = jogadores[vizinhos[4*jogador + i]];
		pay = payoffs[player_state][neigh_state];
		payoff += link * pay;
	}

	return payoff;
}

double W(double Px, double Py){

	double diff_Py_Px;
	double expr;

	diff_Py_Px = Py - Px;
	expr = 1 + exp(-(diff_Py_Px)/k);

	return(1.0/(expr));

}

void troca_estrategia(int *jogadores, int *vizinhos, 
		double Px, double Py, int jogador, int vizinho){

	double random, prob;

	prob = W(Px, Py);
	random = FRANDOM;

	if(random < prob){
		jogadores[jogador] = jogadores[vizinhos[4*jogador + vizinho]];
	}

}

void avalia_links(int *vizinhos, int *links, int *n, int *hist, int *count_links, int t, int jogador){

	int i, j;
	double frac;
	int hist_viz;
	int T = t % M;

	for(i = 0; i < 4; i++){
		if(links[4*jogador + i]){
			frac = 0;
			hist_viz = M*vizinhos[4*jogador + i];
			for(j = 0; j < n[jogador]; j++){
				frac += (double) hist[hist_viz + (T - j + M)%M] / n[jogador];
			}
			if(frac <= crit){
				links[4*jogador + i] = 0; // j -> i
				links[4*vizinhos[4*jogador + i] + ((i + 2) % 4)] = 0; // i -> j
				*count_links -= 2;
				/* printf("%d e %d desligaram\n", 4*jogador + i, 4*vizinhos[4*jogador + i] + ((i + 2) % 4)); */
			}
		}
	}

}

void conta_coop(int *jogadores, int *count){

	int i;

	*count = 0;
	for(i = 0; i < L2; i++){
		*count = *count + jogadores[i];
	}

}

void conta_links(int *links, int *count_links){

	int i;

	*count_links = 0;
	for(i = 0; i < 4*L2; i++){
		*count_links = *count_links + links[i];
	}

}

void atualiza_coop(int *count, int jogador_antigo, int jogador_novo, 
		int atcount[2][2]){

	*count = *count + atcount[jogador_antigo][jogador_novo];	

}

void jogo_tran(int *jogadores, int *vizinhos, int *links, int *hist,
		int t, double payoffs[2][2]){

	int jogador, vizinho;
	double payoff_i, payoff_j;

	jogador = FRANDOM * L2;
	/* printf("MCS %d: jogador %d: %d -> ", t, jogador, jogadores[jogador]); */
	vizinho = FRANDOM * 4;

	payoff_i = teste_calcula_payoff(jogadores, vizinhos, links,
			payoffs, jogador);

	payoff_j = teste_calcula_payoff(jogadores, vizinhos, links,
			payoffs, vizinhos[4*jogador + vizinho]);

	troca_estrategia(jogadores, vizinhos, payoff_i, payoff_j,
			jogador, vizinho);
	/* printf("%d\n", jogadores[jogador]); */
	hist[M*jogador + t] = jogadores[jogador];
	/* hist[4*jogador + (mem[jogador] % M)] = jogadores[jogador]; */
	/* mem[jogador]++; */

}

void jogo_mcs(int *jogadores, int *vizinhos, int *links, int *hist,
		int t, double payoffs[2][2], int *count, int atcount[2][2]){

	int jogador, vizinho;
	int jogador_antigo;
	double payoff_i, payoff_j;

	// escolhe o sitio focal
	jogador = FRANDOM * L2;
	jogador_antigo = jogadores[jogador];
	// escolhe o vizinho
	vizinho = FRANDOM * 4;

	payoff_i = teste_calcula_payoff(jogadores, vizinhos, links,
			payoffs, jogador);

	payoff_j = teste_calcula_payoff(jogadores, vizinhos, links,
			payoffs, vizinhos[4*jogador + vizinho]);

	troca_estrategia(jogadores, vizinhos, payoff_i, payoff_j,
			jogador, vizinho);
	
	hist[M*jogador + t] = jogadores[jogador];
	/* hist[M*jogador + ((T - i + M) % M)] = jogadores[jogador]; */
	/* hist[4*jogador + (mem[jogador] % M)] = jogadores[jogador]; */
	atualiza_coop(count, jogador_antigo, jogadores[jogador], atcount);

}

void escreve_sistema(int *jogadores, int *vizinhos, int *links){
	int i, j;

	for(i = 0; i < L2; i++){
		printf("[%d] = %d; ", i, jogadores[i]);
		for(j = 0; j < 4; j++){
			printf("links[%d] = %d; ", 4*i + j, links[4*i + j]);
		}
		printf("\n");
	}
}

void escreve_memorias(int *hist){
	int i, j;

	for(i = 0; i < L2; i++){
		printf("[%d]\t", i);
		for(j = 0; j < M; j++){
			printf("%d ", hist[M*i + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void diluir(int *vizinhos, int *links, int *n, int *hist, int *count_links, int t){

	int i;

	for(i = 0; i < 4*L2; i++){
		links[i] = 1;
	}

	*count_links = 4*L2;

	for(i = 0; i < L2; i++){
		avalia_links(vizinhos, links, n, hist, count_links, t, i);
	}

}

double power_law(double x){
	return c*(pow(x, -tau));
}


double rng(){
	// Inversa da cumulativa
	double x;
	double rn;

	x = FRANDOM;
	/* x = 1.0*rand()/RAND_MAX; */
	rn = x*pow(M + 1, 1 - tau) - x + 1;
	rn = pow(rn, 1.0/(1 - tau));

	return rn;
}

int int_rng(){
	double rn;

	rn = rng();
	rn = (int) rn;

	return rn;
}

void dist_n(int *n){
	FILE *arq;
	int i;

	arq = fopen("distrib.dsf", "w");
	for(i = 0; i < L2; i++){
		fprintf(arq, "%d\n", n[i]);
	}
	fclose(arq);
}

int transientcheck(int *jogadores, int *vizinhos, int *links, int *hist,
		int *n, FILE *series, int count, int count_links){

	int measures;

	if(count == 0 || count == L2){
		printf("TERMINATED! count = %d\n", count);
		for(measures = 0; measures < MCS/WSIZE; measures++){
			fprintf(series, "%d %d %d\n", measures*WSIZE, count, count_links);
		}
		fclose(series);
		free(jogadores);
		free(vizinhos);
		free(links);
		free(n);
		free(hist);
		return 1;
	}
	else{
		return 0;
	}
}

int measurecheck(int *jogadores, int *vizinhos, int *links, int *hist,
		int *n, FILE *series, int count, int count_links, int measure){

	int i;

	if(count == 0 || count == L2){
		printf("TERMINATED! count = %d; measure = %d\n", count, measure*WSIZE);
		for(i = measure; i < MCS/WSIZE; i++){
			fprintf(series, "%d %d %d\n", i*WSIZE, count, count_links);
		}
		fclose(series);
#ifndef SNAP
		free(jogadores);
		free(links);
#endif
		free(vizinhos);
		free(n);
		free(hist);
		return 1;
	}

	return 0;

}

void gnuplot_view(int *jogadores, FILE *gnu_p){

	int i;
	int j;

	fprintf(gnu_p, "set grid\n");
	fprintf(gnu_p, "set size square\n");
	fprintf(gnu_p, "set xrange[-1:%d]\n", L + 1);
	fprintf(gnu_p, "set yrange[-1:%d]\n", L + 1);
	fprintf(gnu_p, "unset key\n");
	fprintf(gnu_p, "plot \"-\" u 1:($3 == 1? $2 : 1/0) pt 22 lt 22, \"-\" u 1:($3 == 0? $2 : 1/0) pt 22 lt 23\n");
	for(i = 0; i < L; i++){
		for(j = 0; j < L; j++){
			fprintf(gnu_p, "%d %d %d\n", j, i, jogadores[i*L + j]);
		}
	}

	fprintf(gnu_p, "e\n");

	for(i = 0; i < L; i++){
		for(j = 0; j < L; j++){
			fprintf(gnu_p, "%d %d %d\n", j, i, jogadores[i*L + j]);
		}
	}

	fprintf(gnu_p, "e\n");
}

int *lpjogador(int *links){
	// conta quantos links
	// cada jogador tem
	
	int i, j;
	int *lpj;

	lpj = create_int_pointer(L2);

	for(i = 0; i < L2; i++){
		for(j = 0; j < 4; j++){
			lpj[i] += links[4*i + j];
		}
	}

	return lpj;
}

/* #ifdef SNAP */
/* void old_openlat(int opt){ */
/* 	// opens lattice for snapshots */
/* 	// opt = 1 for C/D */
/* 	// opt = 2 for #links/site */

/* 	lat2eps_init(L, L); */
/* 	if(opt == 1){ */
/* 		lat2eps_set_color(0,0xFFFFFF); // white */
/* 		/1* lat2eps_set_color(1,0x000000); // black *1/ */
/* 		/1* lat2eps_set_color(1,0xff7cc5); // some pink *1/ */
/* 		/1* lat2eps_set_color(1,0xeeacf2); // perfume *1/ */
/* 		/1* lat2eps_set_color(1,0xacf2d2); // teal? *1/ */
/* 		lat2eps_set_color(1,0x54e4a2); // teal?? */
/* 		/1* lat2eps_set_color(1,0xc97424); // ochre *1/ */
/* 		/1* lat2eps_set_color(1,0x81160c); // dark burgundy *1/ */
/* 	} */
/* 	if(opt == 2){ */
/* 		lat2eps_set_color(0,0xFFFFFF);	// white */
/* 		/1* lat2eps_set_color(1,0xe5e590);	// weird yellow *1/ */
/* 		/1* lat2eps_set_color(2,0xd2d23d); *1/ */
/* 		/1* lat2eps_set_color(3,0x898920); *1/ */
/* 		/1* lat2eps_set_color(4,0x36360d); *1/ */
/* 		lat2eps_set_color(1,0xc2d6ff);	// blue */
/* 		lat2eps_set_color(2,0x8fa3ec); */
/* 		lat2eps_set_color(3,0x5c70b9); */
/* 		lat2eps_set_color(4,0x293d86); */

/* 	} */
/* } */

/* void old_printlat(int *lat_pt, char *snpsht, int opt){ */
/* 	// prints the lattice to eps file */
/* 	// opt = 1 for C/D */
/* 	// opt = 2 for #links/site */
/* 	int i; */
/* 	if(opt == 1){ */
/* 		for(i = 0; i < L2; i++){ */
/* 			lat2eps_set_site(i % L, L - (i / L) - 1, lat_pt[i]); */
/* 		} */
/* 		lat2eps_gen_eps(snpsht, 0, 0, L, L, 0, 3); */
/* 		lat2eps_release(); */
/* 	} */
/* 	if(opt == 2){ */
/* 		int lpj; */
/* 		for(i = 0; i < L2; i++){ */
/* 			lpj = lpjogador(lat_pt, i); */
/* 			lat2eps_set_site(i % L, L - (i / L) - 1, lpj); */
/* 		} */
/* 		lat2eps_gen_eps(snpsht, 0, 0, L, L, 0, 3); */
/* 		lat2eps_release(); */
/* 	} */
/* } */
/* #endif */
