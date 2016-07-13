/*
 ============================================================================
 Author :
 Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 Group of Theory of Computation
 Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef STRUCTS_GA_H_
#define STRUCTS_GA_H_

#define MAX_N 1000//200 //tamanho maximo de la permutacion inicial

/*define types of print modes */
#define SHOW_FINAL_RESULT   0 //show final results (score and number of evaluations)
#define BEST_BY_GENERATION  1 //show best solution by generation
#define EVALS_BY_GENERATION 2 //show number of evaluations fit funct by generation

typedef struct{
    unsigned int seed;
    int size_population;
    int number_generations;
    int num_eval_fit_function;
    int number_reversals;
    int print_mode; //types of print modes
}parameters;

typedef struct{
	int izq;
	int der;
}gen;

typedef struct{
    int pi[MAX_N]; //permutacio
    int tamanho; //tamanho da permutacion
}permutacion;

typedef struct{
	int *genes;//int genes[MAX_N];
	int tamanho;
	int fitness;
}cromosoma;

typedef struct{
	cromosoma *cromosomas;
	int tamanhoPoblacion;//tamanho de la poblacion, con indice base en 0
	int limiteSeleccion;//Top 50% de la pob o cualquier otro porcentaje
	int baseReemplazo;//A partir de donde se hara el reemplazo

	int basePoblacion;//base de poblacion actual
	int topePoblacion;//tope de poblacion actual

	int baseDescendencia;//base de la descendencia
	int topeDescendencia;//tope de la descendencia

	int tamanhoTotal;//tamanho total = 2 * tamanhoPoblacion (sin usar!!)
}poblacion;

/*****Para MPI****/
typedef struct{
    int index;
    int result;
}msgFitness;

typedef struct{
    int work_number;//MPI: worknumber es indice de cromosoma en la poblacion
    int tamanho;
    int genes[MAX_N];
}msgCromosoma;

#endif /* STRUCTS_GA_H_ */
