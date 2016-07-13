/*
 ============================================================================
 Author :
 Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 Group of Theory of Computation
 Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "calc_fitness.h"
#include "structs_ga.h"

#include "invdist.h"
#include "structs.h"

#include <mpi.h>

int calcular_fitness_mpi(msgCromosoma crom){
	int i,buffer;
	struct genome_struct *genome_list;//lista de genomas
	int NumGenomes = 2;//numero de genomas
	int NumGenes;//numero de genes

	//leer numero de genes
	NumGenes = crom.tamanho;
	NumGenes++;//aumentar en uno el numero de genes, para colocar el "1" al inicio

	//alocar lista de genomas
	genome_list = ( struct genome_struct * ) malloc ( NumGenomes * sizeof ( struct genome_struct ) );

	//leer primer genoma
	genome_list[0].genome_num = 1;
	genome_list[0].genes = ( int * ) malloc ( NumGenes * sizeof ( int ) );//alocar espacio para genes

	genome_list[0].genes[0] = 1;//el primer gen es 1
	for(i=1; i<NumGenes; i++){
		buffer = crom.genes[i-1];
		if (buffer >= 1)//si es positivo incrementar en uno
			genome_list[0].genes[i] = buffer + 1;
		else//si es negativo decrementar uno
			genome_list[0].genes[i] = buffer - 1;
	}

	//generar segundo genoma(identidad)
	genome_list[1].genome_num = 2;
	genome_list[1].genes = ( int * ) malloc ( NumGenes * sizeof ( int ) );//alocar espacio para genes

	for(i=0; i<NumGenes; i++){
		genome_list[1].genes[i] = i+1;
	}

	//calcular la distancia de reversion y mostrar resultado
	int score = invdist_circular_nomem ( genome_list, genome_list + 1, NumGenes);
	//printf ( "score = %d\n", score );

	//liberar memoria
	free(genome_list[0].genes);
	free(genome_list[1].genes);
	free(genome_list);

	return score;
}

int calcular_fitness(cromosoma crom, int *num_eval_fit_function){
	int i,buffer;
	struct genome_struct *genome_list;//lista de genomas
	int NumGenomes = 2;//numero de genomas
	int NumGenes;//numero de genes
    
	//leer numero de genes
	NumGenes = crom.tamanho;
	NumGenes++;//aumentar en uno el numero de genes, para colocar el "1" al inicio
    
	//alocar lista de genomas
	genome_list = ( struct genome_struct * ) malloc ( NumGenomes * sizeof ( struct genome_struct ) );
    
	//leer primer genoma
	genome_list[0].genome_num = 1;
	genome_list[0].genes = ( int * ) malloc ( NumGenes * sizeof ( int ) );//alocar espacio para genes
    
	genome_list[0].genes[0] = 1;//el primer gen es 1
	for(i=1; i<NumGenes; i++){
		buffer = crom.genes[i-1];
		if (buffer >= 1)//si es positivo incrementar en uno
			genome_list[0].genes[i] = buffer + 1;
		else//si es negativo decrementar uno
			genome_list[0].genes[i] = buffer - 1;
	}
    
	//generar segundo genoma(identidad)
	genome_list[1].genome_num = 2;
	genome_list[1].genes = ( int * ) malloc ( NumGenes * sizeof ( int ) );//alocar espacio para genes
    
	for(i=0; i<NumGenes; i++){
		genome_list[1].genes[i] = i+1;
	}
    
	//calcular la distancia de reversion y mostrar resultado
	int score = invdist_circular_nomem ( genome_list, genome_list + 1, NumGenes);
	//printf ( "score = %d\n", score );
    
    /* increment the number of evaluations of fitness function */
    *num_eval_fit_function = *num_eval_fit_function + 1;
    
	//liberar memoria
	free(genome_list[0].genes);
	free(genome_list[1].genes);
	free(genome_list);
    
	return score;
}


