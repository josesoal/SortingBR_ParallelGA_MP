/*
 ============================================================================
 Author :
 Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 Group of Theory of Computation
 Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h> //malloc, free
#include <math.h> //log2
#include "operadores.h"
#include "structs_ga.h"
#include "calc_fitness.h"
#include "ordenar_pob.h"

# include "mpi.h" //processamento paralelo

void poblacion_inicial(poblacion *pob, permutacion *perm, float porcentaje_reemp, float porcentaje_selec, parameters *params){
	int i,j,n,prob;
	
    //Calcular tamanho de la poblacion
	n = perm->tamanho;
	if ( params->size_population > -1 ) {
		pob->tamanhoPoblacion = params->size_population;
	}
	else {
		pob->tamanhoPoblacion =  (int) (n * (log(n) / log(2))); //tamanho de la poblacion = (nlogn), log base 2
	}
	pob->basePoblacion = 0;
	pob->topePoblacion = pob->tamanhoPoblacion;

	pob->limiteSeleccion = pob->tamanhoPoblacion * porcentaje_selec;//limite de pob q seran "seleccionados" para crossover
	pob->baseReemplazo = pob->tamanhoPoblacion * porcentaje_reemp;//base de pob a partir de donde se hara el reemplazo

	pob->baseDescendencia = pob->tamanhoPoblacion;
	pob->topeDescendencia = pob->tamanhoPoblacion + pob->limiteSeleccion;
	pob->tamanhoTotal = pob->topeDescendencia;// tamanho total = poblacion + descendencia
    
    /*alocar memoria*/
	pob->cromosomas = malloc(pob->tamanhoTotal*sizeof(cromosoma)); //alocar memoria para poblacion total
    if (pob->cromosomas == NULL) 
        printf("pob->cromosomas NULL en poblacion_inicial() de operadores.c\n");
    for(i=0;i<pob->tamanhoTotal;i++){
        pob->cromosomas[i].genes = malloc(n * sizeof(int));//alocar memoria para cada cromosoma
        if (pob->cromosomas[i].genes == NULL) 
            printf("pob->cromosomas[i].genes NULL en poblacion_inicial() de operadores.c\n");
    }

	/*crear poblacion a partir de permutacion*/
	for(i=pob->basePoblacion; i < pob->topePoblacion; i++){
		for(j = 0; j < perm->tamanho; j++){
			prob = rand() % 2;//valores entre 0 y 1
			if (prob == 1)//si prob 1, entonces gen positivo
				pob->cromosomas[i].genes[j] = perm->pi[j];
			else//si prob 0, entonces gen negativo
				pob->cromosomas[i].genes[j] = perm->pi[j] * -1;
		}
		pob->cromosomas[i].tamanho = perm->tamanho;
	}
	/*agregar tamanho a desc*/
	for(i=pob->baseDescendencia; i<pob->topeDescendencia; i++){
		pob->cromosomas[i].tamanho = perm->tamanho;
	}
}

void fitness(poblacion *pob, int nro_generacion, int *num_eval_fit_function){
	/*calcular fitness para cada elemento de la poblacion
	 * utilizando el algoritmo lineal para distancia de reversion*/
	int i;
	if (nro_generacion > 1){ //fitness sobre descendencia
		for(i=pob->baseDescendencia; i < pob->topeDescendencia; i++){
			pob->cromosomas[i].fitness = calcular_fitness(pob->cromosomas[i], num_eval_fit_function);
		}
	}
	else{ //nro_generacion == 1, fitness sobre poblacion
		for(i=pob->basePoblacion; i < pob->topePoblacion; i++)
			pob->cromosomas[i].fitness = calcular_fitness(pob->cromosomas[i], num_eval_fit_function);
	}
}

void seleccion(poblacion *pob,  int nro_generacion, int *mejor_solucion){ //cromosoma *mejor_solucion){
	/*ordenar la poblacion*/
	//quickSort(pob,pob->basePoblacion,pob->topePoblacion-1);
	
    countingSort(pob,pob->cromosomas[0].tamanho+2);

	/*actualizar mejor solucion*/
	//if (pob->cromosomas[0].fitness < mejor_solucion->fitness)
		//*mejor_solucion = pob->cromosomas[0];
    if (pob->cromosomas[0].fitness < *mejor_solucion)
        *mejor_solucion = pob->cromosomas[0].fitness;

}

void crossover3(poblacion *pob, permutacion *perm, float prob_crossover){
	int i,j,k,pos_padre1, pos_padre2, punto_crossover;
	float prob;

	k = pob->baseDescendencia;//indice para colocar en descendencia
	/*croosver sobre poblacion actual, hasta el limite de seleccion*/
	for(i=pob->basePoblacion; i<pob->limiteSeleccion; i+=2){
		//si falta ultima posicion para colocar descendencia, terminar
		if (k+1 == pob->tamanhoTotal)
			break;

		//escoger posiciones de los padres
		pos_padre1 = rand() % pob->limiteSeleccion;//valores entre 0 y limiteSeleccion-1
		pos_padre2 = rand() % pob->limiteSeleccion;//valores entre 0 y limiteSeleccion-1
		//buscar pos2 que no sea igual a pos1, por lo menos 5 veces
		for(j=0; j<5; j++){
			if (pos_padre1 != pos_padre2) break;
			pos_padre2 = rand() % pob->limiteSeleccion;//valores entre 0 y limiteSeleccion-1
		}

		prob = ((rand() % 100) + 1)/100.0; //de 0.01, 0.02 ... a 1.0
		//printf("----->%f\n",prob);
		if (prob <= prob_crossover){//crossover
			punto_crossover = rand() % perm->tamanho;//valores entre 0 tamanho de permutacion - 1
			//aplicar single_point_crossover a la primera parte antes del punto de crossover
			for(j=0; j<punto_crossover; j++){
				pob->cromosomas[k].genes[j] = pob->cromosomas[pos_padre1].genes[j];
				pob->cromosomas[k+1].genes[j] = pob->cromosomas[pos_padre2].genes[j];
			}
			//aplicar single_point_crossover a la segunda parte desp. del punto de crossover
			for(j=punto_crossover; j<perm->tamanho; j++){
				pob->cromosomas[k+1].genes[j] = pob->cromosomas[pos_padre1].genes[j];
				pob->cromosomas[k].genes[j] = pob->cromosomas[pos_padre2].genes[j];
			}
			k = k + 2;
		}
	}
	pob->topeDescendencia = k;//actualizar tope descendencia
}

void mutacion(poblacion *pob, permutacion *perm, float prob_mutacion){
	int i,j;
	float prob;
	/*mutacion sobre descendencia*/
	for(i=pob->baseDescendencia; i<pob->topeDescendencia; i++){
		for(j=0; j<perm->tamanho; j++){
			prob = ((rand() % 100) + 1)/100.0; //de 0.01, 0.02 ... a 1.0
			if(prob <= prob_mutacion){//mutar gen (invertir signo)
				pob->cromosomas[i].genes[j] = pob->cromosomas[i].genes[j] * -1;
			}

		}
	}

}

void reemplazo(poblacion *pob){
	int i,j,pos_reemp;
	/*reemplaza la descendencia en poblacion actual*/
	for(i=pob->baseDescendencia; i<pob->topeDescendencia; i++){
		//buscar pos de reemplazo
		int tamanhoReemp = pob->topePoblacion - pob->baseReemplazo;
		pos_reemp = (rand() % tamanhoReemp) + pob->baseReemplazo;//(valores entre 0 y tamanhoReemp-1)+base
		//realizar reemplazo
		if (pob->cromosomas[i].fitness < pob->cromosomas[pos_reemp].fitness){
			pob->cromosomas[pos_reemp].fitness = pob->cromosomas[i].fitness;
            
            for(j=0; j<pob->cromosomas[0].tamanho; j++)
                pob->cromosomas[pos_reemp].genes[j] = pob->cromosomas[i].genes[j];
            
        }
	}
}
