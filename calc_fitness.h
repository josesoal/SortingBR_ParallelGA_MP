/*
 ============================================================================
 Author :
 Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 Group of Theory of Computation
 Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef CALC_FITNESS_H_
#define CALC_FITNESS_H_

#include "structs_ga.h"

int calcular_fitness_mpi(msgCromosoma crom);
int calcular_fitness(cromosoma crom, int *num_eval_fit_function);


#endif /* CALC_FITNESS_H_ */
