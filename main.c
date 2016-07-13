/*
 ============================================================================
 Project    : Parallel Genetic Algorithm for the SUPBR Problem
                (version with multiple populations)
 File       : main.c
 Author     : Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
                Group of Theory of Computation
                Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>//time
#include <mpi.h>
#include <math.h>
#include <string.h>

#include "operadores.h"
#include "structs_ga.h"
#include "calc_fitness.h"


#define WORKTAG 1
#define DIETAG 2
#define EVALTAG 3
#define GOTAG 4

#define DEBUG 0
#define TAM_BUFFER 1000

void readCommandLine(int argc, char *argv[], parameters *params, int numprocs);
void leer_permutacion(permutacion *perm);
void leer_permutacion_de_archivo(char* nombre_archivo, permutacion *perm);
void mostrar_poblacion(poblacion *pob, permutacion *perm, int nro_generacion);
void master(parameters *params);
void slave(parameters *params);

int main(int argc,char **argv)
{
    int myid, numprocs, proc;
    double wtime;
    parameters params;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    params.seed = time(NULL);
    params.size_population = -1;
    params.number_generations = -1;
    params.num_eval_fit_function = -1;
    params.number_reversals = -1;
    params.print_mode = SHOW_FINAL_RESULT;

    readCommandLine(argc, argv, &params, numprocs);

    if(myid == 0){
        wtime = MPI_Wtime ( );//get initial time
        master(&params);
        wtime = MPI_Wtime ( ) - wtime;//get total time
        //printf ( "\tTotal time = %.16g seconds.\n\n", wtime );
        for (proc = 1; proc<numprocs; proc++){
            MPI_Send(0, 0, MPI_BYTE, proc, DIETAG, MPI_COMM_WORLD);
        }
    }
    else{
        slave(&params);
       
    }
    MPI_Finalize();     
    return 0;
}

void master(parameters *params)
{
    int fit_evals[2];//[0] for fitness, and [1] for evaluations of fitness
    int num_generaciones, num_evaluations, num_used_generations;
    permutacion perm;
    int solution, recv_fitness, recv_evals, i, k, myid, numprocs, proc;
    MPI_Status status;    
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    
    //Read Permutation
    leer_permutacion_de_archivo("i", &perm);
    
    //Send Permutation to each slave
    for (proc = 1; proc<numprocs; proc++){
        MPI_Send(&perm, 
                 sizeof(permutacion), 
                 MPI_BYTE, 
                 proc,
                 WORKTAG, 
                 MPI_COMM_WORLD);
    }
    
    //Receive fitness from each slave, for each generation
    solution = perm.tamanho; //assume the solution is the length of permutation
    num_generaciones = perm.tamanho;
    if (params->number_generations > -1) 
        num_generaciones = params->number_generations;

    for (k = 1; k <= num_generaciones; k++) {
        //send empty message to all slaves to start working
        for (i = 0; i < numprocs - 1; i++){
            MPI_Send(0, 0, MPI_BYTE, i % ( numprocs - 1 ) + 1, GOTAG, MPI_COMM_WORLD);
        }
        //receive fitness from all slaves
        num_evaluations = 0;
        for (i = 0; i < numprocs - 1; i++){
            MPI_Recv(&fit_evals, 
                    2, 
                    MPI_INT,  
                    MPI_ANY_SOURCE, //from any slave ID 
                    EVALTAG,
                    MPI_COMM_WORLD, 
                    &status);
            recv_fitness = fit_evals[0];
            recv_evals = fit_evals[1];

            num_evaluations = num_evaluations + recv_evals;
            if (recv_fitness < solution) solution = recv_fitness;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (params->print_mode == BEST_BY_GENERATION){
            printf("%d ",solution); //print best solution by generation
        }
        else if (params->print_mode == EVALS_BY_GENERATION){
            printf("%d ", num_evaluations); //print number of evaluations
        }
        //other stop conditions
        if (params->num_eval_fit_function > -1 &&
            num_evaluations >= params->num_eval_fit_function){
            break;
        }
        if (solution <= params->number_reversals){
            break;
        } 
    }

    if (k > num_generaciones)
        num_used_generations = num_generaciones;
    else
        num_used_generations = k;

    if (params->print_mode == SHOW_FINAL_RESULT){
        printf("%d\n",solution);//show result (number of reversals)
        printf("%d\n",num_evaluations);//show number of evalutions of fit_function
        printf("%d\n",num_used_generations);//show number of used generations
    }
}

void slave(parameters *params)
{
    int fit_evals[2];//[0] for fitness, and [1] for evaluations of fitness
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Status status; 
    srand ( params->seed * (myid+1) );//initialize random generator
    
    /*NOTE: The parameters are the same proposed in the paper:
     [Parallelization and virtualization of genetic algorithms for 
     sorting permutations by reversals] */
    
    /*Parametros*/
	float prob_crossover = 0.9;//probabilidad de crossover: prob_crossover
	float prob_mutacion = 0.02;//probabilidad de mutacion: prob_mutacion
	float porcentaje_selec = 0.6;//porcentaje de poblacion que seran seleccionados para crossover (60%)
	float porcentaje_reemp = 1-0.6;//porcentaje de poblacion a partir de donde se realizara el reemplazo (60%)
    //de la descendencia
	/*Variables*/
	int i, nro_generacion;
	permutacion perm;//permutacion
	poblacion pob;//poblacion
	int mejor_solucion;//mejor solucion encontrada
	int num_generaciones;
    int num_eval_fit_function = 0;
    int wait_for_dietag = 1;
    
	/*Receive permutation from master*/
    MPI_Recv(&perm, 
             sizeof(permutacion), 
             MPI_BYTE,  
             0, //from master 
             WORKTAG,
             MPI_COMM_WORLD, 
             &status);
    
    /*Iniciar ga*/
    MPI_Recv(0, 0, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (status.MPI_TAG == DIETAG) return; //no allocated memory
    
    poblacion_inicial(&pob,&perm,porcentaje_reemp, porcentaje_selec, params);//crear poblacion inicial
    fitness(&pob,1,&num_eval_fit_function);//calcular fitness pob inicial
    mejor_solucion = pob.cromosomas[0].fitness;
    num_generaciones = perm.tamanho;//num de generaciones = a tamanho de permutacion
    if (params->number_generations > -1) num_generaciones = params->number_generations;
   
    fit_evals[0] = mejor_solucion;
    fit_evals[1] = num_eval_fit_function;
    MPI_Send(&fit_evals, 
             2, 
             MPI_INT, 
             0, //Master ID  
             EVALTAG, 
             MPI_COMM_WORLD);//send a solution to the master
    MPI_Barrier(MPI_COMM_WORLD);

    for(nro_generacion=2; nro_generacion<=num_generaciones; nro_generacion++){
     
        MPI_Recv(0, 0, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            wait_for_dietag = 0;//is dead dont need to wait
            break;
        }

        seleccion(&pob,nro_generacion,&mejor_solucion);//seleccion:ordenar poblacion
        crossover3(&pob, &perm, prob_crossover);//crossover:generar descendencia
        mutacion(&pob, &perm, prob_mutacion);//mutacion:actuar sobre toda la descendencia (de tam_desc hasta tamanho de descencia)
        fitness(&pob,nro_generacion,&num_eval_fit_function);//calcular fitness de descendencia
        reemplazo(&pob);//reemplazar descendencia en nueva generacion

        fit_evals[0] = mejor_solucion;
        fit_evals[1] = num_eval_fit_function;
        MPI_Send(&fit_evals, 
                2, 
                MPI_INT, 
                0, //Master ID  
                EVALTAG, 
                MPI_COMM_WORLD);//send a solution to the master
        MPI_Barrier(MPI_COMM_WORLD);        
    }

    if ( wait_for_dietag == 1 ) 
        MPI_Recv(0, 0, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    /*Liberar memoria*/
    for(i=0;i<pob.tamanhoTotal;i++){
        free(pob.cromosomas[i].genes);
    }
    free(pob.cromosomas);
}

void readCommandLine(int argc, char *argv[], parameters *params, int numprocs) 
{
    int i, g_option, t_option;
    char option;

    g_option = 0;
    t_option = 0;
    /* read parameters from command line */
    i = 1;
    while ( i < argc) {
        if ( argv[ i ][ 0 ] == '-' && (i + 1) < argc ) { 
            option = argv[ i ][ 1 ]; 
            switch ( option ) {
                case 's':
                    params->seed = atoi( argv[ i + 1 ] );
                    break;
                case 'g':
                    g_option = 1;
                    params->number_generations = atoi( argv[ i + 1 ] );
                    break;
                case 't':
                    t_option = 1;
                    params->number_generations = atoi( argv[ i + 1 ] );
                    //divide the number of generations among the processes
                    //and, use the ceil function over this division
                    if (params->number_generations % (numprocs-1) == 0)
                        params->number_generations = 
                            params->number_generations / (numprocs-1);
                    else
                        params->number_generations = 
                            (params->number_generations / (numprocs-1)) + 1;
                    break;
                case 'e':
                    params->num_eval_fit_function = atoi( argv[ i + 1 ] );
                    break;
                case 'r':
                    params->number_reversals = atoi( argv[ i + 1 ] );
                    break;
                case 'm':
                    if ( strcmp( argv[ i + 1], "final_result" ) == 0 ) {
                        params->print_mode = SHOW_FINAL_RESULT;
                    }
                    else if ( strcmp( argv[ i + 1], "best_by_gen" ) == 0 ) {
                        params->print_mode = BEST_BY_GENERATION;
                    }
                    else if ( strcmp( argv[ i + 1], "eval_by_gen" ) == 0 ) {
                        params->print_mode = EVALS_BY_GENERATION;
                    }
                    else{
                        fprintf( stderr, " stderr: incorrect mode (-m).\n" );
                        exit( EXIT_FAILURE );
                    }
                    break;
                case 'p':
                    params->size_population = atoi( argv[ i + 1 ] );
                    break;
                default:
                    fprintf( stderr, " stderr: incorrect option: %c.\n", option );
                    exit( EXIT_FAILURE );
            }
            i = i + 2;  
        }
        else{
            fprintf( stderr, " stderr: incorrect options or parameters.\n" );
            exit( EXIT_FAILURE );
        }
    }//end-while

    /* verify that option "-g" and "-t" are not used at the same time */
    if ( g_option == 1 && t_option == 1 ) {
        fprintf( stderr, " stderr: options \"-g\" and \"-t\" can not be used at the same time.\n" );
        exit( EXIT_FAILURE );
    }

}

void leer_permutacion(permutacion *perm){
	int i,n;
	scanf("%d",&n);//leer tamanho de permutacion
	perm->tamanho = n;
	for(i=0; i<n; i++){
		scanf("%d",&perm->pi[i]);
	}
}

void leer_permutacion_de_archivo(char* nombre_archivo, permutacion *perm){
    int i,n;
    FILE *filePointer;
    if ((filePointer = fopen(nombre_archivo,"r")) == NULL){
        printf("Input file could not be opened\n");
    }
    else{
        fscanf(filePointer,"%d",&n);//leer tamanho de permutacion
        perm->tamanho = n;
        for(i=0; i<n; i++){
            fscanf(filePointer,"%d",&perm->pi[i]);
        }
    }
    fclose(filePointer);
}

void mostrar_poblacion(poblacion *pob, permutacion *perm, int nro_generacion){
	
    int i,j;
    printf("*** Generacion %d ***\n", nro_generacion);
    //imprimir poblacion actual
    for(i=pob->basePoblacion; i<pob->topePoblacion; i++){
        printf("{");
        for(j=0; j<perm->tamanho-1; j++){
            printf("%d, ", pob->cromosomas[i].genes[j]);
        }
        printf("%d", pob->cromosomas[i].genes[perm->tamanho-1]);
        printf("} fit: %d\n",pob->cromosomas[i].fitness);
    }
    //imprimir descendencia
    printf("-------------------\n");
    for(i=pob->baseDescendencia; i<pob->topeDescendencia; i++){
        printf("{");
        for(j=0; j<perm->tamanho-1; j++){
            printf("%d, ", pob->cromosomas[i].genes[j]);
        }
        printf("%d", pob->cromosomas[i].genes[perm->tamanho-1]);
        printf("} fit: %d\n",pob->cromosomas[i].fitness);
    }
	
}


