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
#include "ordenar_pob.h"
#include "structs_ga.h"

# include "mpi.h" //processamento paralelo

/*void swap2(poblacion *pob,int i, int j) {
    cromosoma cromosomaT;
    cromosomaT = pob->cromosomas[i];
    pob->cromosomas[i] = pob->cromosomas[j];
    pob->cromosomas[j] = cromosomaT;
}

int partition(poblacion *pob, int left, int right) {
    int i, j;

    i = left;
    for (j = left + 1; j <= right; ++j) {
        if (pob->cromosomas[j].fitness < pob->cromosomas[left].fitness) {
            ++i;
            swap2(pob,i,j);
        }
    }
    swap2(pob,left,i);

    return i;
}

void quickSort(poblacion *pob, int left, int right) {
    int r;

    if (right > left) {
        r = partition(pob, left, right);
        quickSort(pob, left, r - 1);
        quickSort(pob, r + 1, right);
    }
}*/

void countingSort(poblacion *pob, int k){
	int i,j;
	int c[k];
	cromosoma *a;

	/*alocar memoria para poblacion temporal "a" */
	a = malloc(pob->tamanhoPoblacion*sizeof(cromosoma));
	if (a == NULL) printf("a NULL en countingSort() de ordenar_pob.c \n");
    
    for(i=0;i<pob->tamanhoPoblacion;i++){
        a[i].genes = malloc(pob->cromosomas[0].tamanho * sizeof(int));//alocar memoria para cada cromosoma
        if (a[i].genes == NULL) 
            printf("a[i].genes NULL en countingSort() de ordenar_pob.c \n");
    }
    
    /*hacer una copia poblacion*/
	for(i=0; i<pob->tamanhoPoblacion; i++){
		a[i].fitness = pob->cromosomas[i].fitness;
        for(j=0; j<pob->cromosomas[0].tamanho; j++)
            a[i].genes[j] = pob->cromosomas[i].genes[j]; 
	}
    
	/*counting sort*/
	for(i=0; i<k; i++){
		c[i] = 0;
	}

	for(i=0; i<pob->tamanhoPoblacion; i++){
		c[a[i].fitness] = c[a[i].fitness] + 1;
	}

	for(i=1; i<k; i++){
		c[i] = c[i] + c[i-1];
	}

	for(i=pob->tamanhoPoblacion-1; i>=0; i--){
		//copiar en poblacion
        pob->cromosomas[c[a[i].fitness]-1].fitness = a[i].fitness;
        for(j=0; j<pob->cromosomas[0].tamanho; j++)
            pob->cromosomas[c[a[i].fitness]-1].genes[j] = a[i].genes[j];
        
        
		c[a[i].fitness] = c[a[i].fitness] - 1;
	}
    
    
    /*Liberar memoria*/
    for(i=0;i<pob->tamanhoPoblacion;i++){
        free(a[i].genes);
    }
	free(a);
}


