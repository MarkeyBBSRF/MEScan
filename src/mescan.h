#ifndef MESCAN_H
#define MESCAN_H

#define _GNU_SOURCE
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

typedef struct {
  int nDepth;
  int nRow;
  int nCol;
  int *** Matrix;
  int *** bMatrix;
} mut_categs;

typedef struct {
  int nDepth;
  int nRow;
  int nCol;
  double *** Matrix;
} mut_rate_categs;


/* mutation data matrix */
typedef struct {
  int nRow;
  int nCol;
  int ** Matrix;
  int ** bMatrix;
} mutation_data;

/* mutation rate matirx */
typedef struct {
  int nRow;
  int nCol;
  double ** Matrix;
} mutation_rate;


/* geneset index */
typedef struct {
  size_t len;
  int *idx;
} geneset_index;


geneset_index *init_index(int *indx, size_t length);
mutation_data vect_to_mutation_data(const int *vect, int *size);
mutation_rate vect_to_mutation_rate(const double *vect, int *size);
mutation_data * init_mutation_data(int nrow, int ncol);
mutation_rate * init_mutation_rate(int nrow, int ncol);
void free_geneset_idx(geneset_index * g);
void update_geneset_idx(geneset_index * dest, geneset_index * src);

extern int* preAddPatCoverage;
geneset_index *get_geneset_index_from_array(int * array, int size, int num);
int compare_int_array(int *a, int*b, int length);
mutation_rate * get_complement_mut_rate(mutation_rate matrix);
double *get_column_prod(mutation_rate *matrix, geneset_index *G);
double get_tg_score(mutation_rate dmatrix, mutation_rate * comp, mutation_data imatrix, geneset_index *G, double tuning);
void r_calculate_TG(
  int *rMutMat, int *sizeMutMat,
  double *rRateMat, int *sizeRateMat,
  int *rCombMat, int *sizeCombMat,
  double *rTuning, double *rOut_TG
);
int sample_next_gene(int min_num, int max_num);
geneset_index *get_next_geneset_idx(int* curr_set, int next_gene_idx, int size, int replace);

double flipcoin();
int rand_int();
unsigned long get_seed(unsigned long a, unsigned long b, unsigned long c);
void write_file_int_vect(int * vec, int size, FILE *file);

#endif
