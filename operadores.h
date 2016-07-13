/*
 ============================================================================
 Author :
 Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 Group of Theory of Computation
 Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef OPERADORES_H_
#define OPERADORES_H_

#include "structs_ga.h"

void poblacion_inicial(poblacion *pob, permutacion *perm, float porcentaje_reemp, float porcentaje_selec, parameters *params);
void fitness(poblacion *pob, int nro_generacion, int *num_eval_fit_function);
void seleccion(poblacion *pob, int nro_generacion, int *mejor_solucion);
void crossover3(poblacion *pob, permutacion *perm, float prob_crossover);
void mutacion(poblacion *pob, permutacion *perm, float prob_mutacion);
void reemplazo(poblacion *pob);
#endif /* OPERADORES_H_ */
