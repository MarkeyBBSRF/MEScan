#include "mescan.h"

int * preAddPatCoverage;

// Return TG scores for the gene sets in rCombMat
// rMutMat: mutation matrix passed by R
// rRateMat: mutation rate matrix passed by R
// rCombMat: gene sets combinations matrix passed by R
// rOut_TG: vector stores the calculated TG scores
void r_calculate_TG(
  int *rMutMat, int *sizeMutMat,
  double *rRateMat, int *sizeRateMat,
  int *rCombMat, int *sizeCombMat,
  double *rTuning, double *rOut_TG
) {
  /* pass input to c */
  mutation_data MutMat = vect_to_mutation_data(rMutMat, sizeMutMat);
  mutation_rate RateMat = vect_to_mutation_rate(rRateMat, sizeRateMat);
  mutation_data combMat = vect_to_mutation_data(rCombMat, sizeCombMat);

  if (MutMat.nRow != RateMat.nRow || MutMat.nCol != RateMat.nCol) {
    perror("Dimensions differ in MutMat and RateMat");
    exit(EXIT_FAILURE);
  }
  int Gsize = combMat.nCol;
  int nG = combMat.nRow;
  int G[Gsize];
  if (Gsize < 2) {
    perror("The gene set must have at least 2 genes");
    exit(EXIT_FAILURE);
  }

  /* 1-Mutation_Rate */
  mutation_rate *compRateMat = get_complement_mut_rate(RateMat);
  double tg_tmp;
  int i;
  geneset_index * Gi;
  for (i = 0; i < nG; i++) {
    memcpy(G, combMat.Matrix[i], sizeof(G));
    Gi = init_index(G, Gsize);
    // print_int_vect(G,Gi->len);
    tg_tmp = get_tg_score(RateMat, compRateMat, MutMat, Gi, *rTuning);
    rOut_TG[i] = tg_tmp;
  }
  free(compRateMat->Matrix[0]);
  free(compRateMat);
}

// Return non-zero index from the binary array
// For example:
// For the binary array [0,1,0,0,1,1] return index: [1,4,5]
geneset_index * get_geneset_index_from_array(int * array, int size, int num) {
  geneset_index *tmpI;
  int i;
  int count = 0;
  int *tmp = malloc(size * sizeof(int));
  for ( i = 0; i < size; i++) {
    if (array[i] == num) {
      tmp[count] = i;
      count++;
    }
  }
  if (count == 0) {
    free(tmp);
    return NULL;
  } else {
    tmpI = init_index(tmp, count);
    return tmpI;
  }

}

// Return 1 minus mutation rate
mutation_rate * get_complement_mut_rate(mutation_rate matrix) {
  mutation_rate *newMatrix = init_mutation_rate(matrix.nRow, matrix.nCol);
  int i;
  int j;
  for (i = 0; i < matrix.nRow; i++) {
    for (j = 0; j < matrix.nCol; j++) {
      newMatrix->Matrix[i][j] = 1 - matrix.Matrix[i][j];
    }
  }
  return newMatrix;
}

// Calculate product by column
double *get_column_prod(mutation_rate *matrix, geneset_index *G) {
  double * prod = malloc(matrix->nCol * sizeof(double));
  int i;
  int j;
  for (j = 0; j < matrix->nCol; j++) {
    prod[j] = matrix->Matrix[G->idx[0]][j];

  }

  for (j = 0; j < matrix->nCol; j++) {
    for (i = 1; i < G->len; i++) {
      prod[j] = prod[j] * matrix->Matrix[G->idx[i]][j];
    }

  }
  return prod;
}

// calculate tg score
double get_tg_score(mutation_rate dmatrix,
                    mutation_rate * compRateMat,
                    mutation_data imatrix, geneset_index *G, double tuning)
{
  int i;
  int j;
  int * colSum = (int *)malloc(sizeof(int) * imatrix.nCol);
  for (i = 0; i < imatrix.nCol; ++i) {
    colSum[i] = 0;
    for (j = 0; j < G->len; ++j) {
      colSum[i] = colSum[i] + imatrix.bMatrix[G->idx[j]][i];
    }
  }

  // count ME size for each gene
  int * nME = (int *)malloc(sizeof(int) * G->len);
  double totalnME = 0;

  for (j = 0; j < G->len; ++j) {
    nME[j] = 0;
    for (i = 0; i < imatrix.nCol; ++i) {
      if (colSum[i] == 1 && imatrix.bMatrix[G->idx[j]][i] > 0)
        nME[j]++;
    }

    if (nME[j] == 0) {
      free(nME);
      return 0;
    }
    totalnME = totalnME + 1.0 / nME[j];
  }
  int mCount = 0;
  double *compRateMatProd = get_column_prod(compRateMat, G);
  double *weight = (double *)malloc(sizeof(double) * G->len);
  for (j = 0; j < G->len; ++j) {
    weight[j] = (1.0 / nME[j]) / totalnME;
  }
  double uSum = 0;
  // Tvec
  double T = 0;
  // ETvect
  double eT = 0;
  // VarTvec
  double varT = 0;
  // T_G
  double T_G = 0;
  double thetaij = 0;
  double tmp1;
  double Tnom = 0;
  double Tdnom = 0;
  double varT1 = 0;
  double varT2;
  double dtij;
  double Sij;

  for (j = 0; j < dmatrix.nCol; j++) {
    T = 0;
    eT = 0;
    uSum = 0;
    mCount = 0;
    Tnom = 0;
    Tdnom = 0;
    varT1 = 0;
    varT2 = 0;

    for (i = 0; i < G->len; i++) {
      mCount = imatrix.Matrix[G->idx[i]][j];
      uSum = uSum + mCount; // track # of mutations for patient(ij)
      thetaij = compRateMatProd[j] / compRateMat->Matrix[G->idx[i]][j] * dmatrix.Matrix[G->idx[i]][j];
      tmp1 = thetaij * (1 - thetaij) * (1 - thetaij) * (1 - thetaij);
      dtij = sqrt(tmp1) + tuning;
      Tnom = Tnom + weight[i] * (mCount - thetaij) * mCount;
      Tdnom = Tdnom + mCount * dtij;
      Sij = weight[i] * thetaij * (1 - thetaij) / dtij; //TempMat in R
      eT = eT + Sij; // ETvect in R
      varT1 = varT1 + weight[i] * weight[i] * tmp1 / dtij / dtij;
      varT2 = varT2 + Sij * Sij;
    }

    if (Tdnom == 0 || uSum != 1) {
      T = 0;
    } else {
      T = Tnom / Tdnom; // Tvect in R
    }
    varT = varT + varT1 + varT2 - eT * eT;
    T_G = T_G + (T - eT);

  }

  T_G = T_G / sqrt(varT);
  free(compRateMatProd);
  free(colSum);
  free(nME);
  return T_G;
}


int sample_next_gene(int min_num, int max_num)
{
  int result = 0, low_num = 0, hi_num = 0;
  if (min_num < max_num)
  {
    low_num = min_num;
    hi_num = max_num + 1;
  } else {
    low_num = max_num + 1;
    hi_num = min_num;
  }

  result = (rand() % (hi_num - low_num)) + low_num;
  return result;
}


// Compare two pre-sorted int array
// Return 1 if they are the same, 0 otherwise
int compare_int_array(int *a, int *b, int length) {
  int i;
  int same = 1;
  for (i = 0; i < length; i++) {
    if (a[i] != b[i]) {
      same = 0;
      break;
    }
  }
  return same;

}

int* get_sample_coverage_exclude(mutation_data matrix, geneset_index * geneIndex, int numPatient, int exclude) {
  if (geneIndex->len == 0) {
    perror("geneset_index has ZERO gene index called in --> get_sample_coverage");
    exit(EXIT_FAILURE);
  }
  int *ptCoverage = malloc(numPatient * sizeof(int));
  int i, j;
  for (j = 0; j < numPatient; j++) {
    for (i = 0; i < geneIndex->len ; i++) {
      if ( i != exclude) {
        if (matrix.Matrix[geneIndex->idx[i]][j] > 0) {
          ptCoverage[j] = 1;
          break;
        } else {
          ptCoverage[j] = 0;
        }
      }

    }
  }
  return ptCoverage;
}

int* get_sample_coverage(mutation_data matrix, geneset_index * geneIndex, int numPatient) {
  if (geneIndex->len == 0) {
    perror("geneset_index is 0 in --> get_sample_coverage");
    exit(EXIT_FAILURE);
  }
  int *ptCoverage = malloc(numPatient * sizeof(int));
  int i, j;
  for (j = 0; j < numPatient; j++) {
    for (i = 0; i < geneIndex->len ; i++) {
      if (matrix.Matrix[geneIndex->idx[i]][j] > 0) {
        ptCoverage[j] = 1;
        break;
      } else {
        ptCoverage[j] = 0;
      }
    }
  }
  return ptCoverage;
}


geneset_index *get_next_geneset_idx(int* curr_set, int next_gene_idx, int size, int replace) {
  /* turn on/off candidate gene */
  curr_set[next_gene_idx] = 1 - curr_set[next_gene_idx];
  if (replace != -1) {
    curr_set[replace] = 1 - curr_set[replace];
  }
  geneset_index *next_geneset_idx = get_geneset_index_from_array(curr_set, size, 1);
  if (next_geneset_idx == NULL) {
    return NULL;
  } else {
    return next_geneset_idx;
  }
}



geneset_index *init_index(int *indx, size_t length) {
  geneset_index *nIdx;
  nIdx = (geneset_index *)malloc(sizeof(geneset_index));
  nIdx->len = length;
  nIdx->idx = indx;
  return nIdx;
}

/*******************Data Reformat*******************/

/* convert vect input from R to 2d array */
mutation_data vect_to_mutation_data(const int *vect, int *size) {
  mutation_data data;
  int i, j;
  data.nRow = size[0];
  data.nCol = size[1];
  data.Matrix = malloc(data.nRow * sizeof(int*));
  data.bMatrix = malloc(data.nRow * sizeof(int*));
  int* values = calloc(data.nRow * data.nCol, sizeof(int));
  int* binaries = calloc(data.nRow * data.nCol, sizeof(int));
  for (i = 0; i < data.nRow; ++i)
  {
    data.Matrix[i] = values + i * data.nCol;
    data.bMatrix[i] = binaries + i * data.nCol;
  }
  if (data.Matrix == NULL)
  {
    printf("Error in allocating 2D memory\n");
    exit( EXIT_FAILURE);
  }
  int ndata = 0;
  for (i = 0; i < data.nRow; i++) {
    for (j = 0; j < data.nCol; j++) {
      data.Matrix[i][j] =  vect[j * data.nRow + i];
      if (data.Matrix[i][j] > 0) {
        data.bMatrix[i][j] = 1;
      } else {
        data.bMatrix[i][j] = 0;
      }
      ndata++;

    }
  }

  ndata -= data.nCol * data.nRow;
  if (ndata != 0)
  {
    printf("%s data values in input\n", ndata > 0 ? "Too few\n" : "Too many\n");
    exit(EXIT_FAILURE);
  }
  else {
    return data;
  }

}

mutation_rate vect_to_mutation_rate(const double *vect, int *size) {
  mutation_rate data;
  int i;
  int j;
  data.nRow = size[0];
  data.nCol = size[1];

  data.Matrix = malloc(data.nRow * sizeof(double*));
  double* values = calloc(data.nRow * data.nCol, sizeof(double));
  for (i = 0; i < data.nRow; ++i)
  {
    data.Matrix[i] = values + i * data.nCol;
  }
  if (data.Matrix == NULL)
  {
    printf("Error in allocating 2D memory\n");
    exit( EXIT_FAILURE);
  }
  int ndata = 0;
  for (i = 0; i < data.nRow; i++) {
    for (j = 0; j < data.nCol; j++) {
      data.Matrix[i][j] =  vect[j * data.nRow + i];
      ndata++;
    }
  }

  ndata -= data.nCol * data.nRow;

  if (ndata != 0)
  {
    printf("%s data values in input\n", ndata > 0 ? "Too few\n" : "Too many\n");
    exit(EXIT_FAILURE);
  }
  else {
    return data;
  }
}

mutation_rate * init_mutation_rate(int nrow, int ncol) {
  mutation_rate *newM;
  int i;
  newM = (mutation_rate*)malloc(sizeof(mutation_rate));
  newM->nRow = nrow;
  newM->nCol = ncol;
  newM->Matrix = malloc(nrow * sizeof(double*));
  double* values = calloc(nrow * ncol, sizeof(double));
  for (i = 0; i < nrow; ++i)
  {
    newM->Matrix[i] = values + i * ncol;
  }
  return newM;
}

mutation_data *init_mutation_data(int nrow, int ncol) {
  mutation_data *newM;
  int i;
  newM = (mutation_data*)malloc(sizeof(mutation_data));
  newM->nRow = nrow;
  newM->nCol = ncol;
  newM->Matrix = malloc(nrow * sizeof(int*));
  newM->bMatrix = malloc(nrow * sizeof(int*));
  int* values = calloc(nrow * ncol, sizeof(int));
  for (i = 0; i < nrow; ++i)
  {
    newM->Matrix[i] = values + i * ncol;
    newM->bMatrix[i] = values + i * ncol;
  }
  return newM;
}

void update_geneset_idx(geneset_index * dest, geneset_index * src) {
  // if(dest->idx != NULL){
  // 	free(dest -> idx);
  // }
  dest->idx = malloc(src->len * sizeof(int));
  dest->len = src->len;
  for (int i = 0; i < src->len; i++) {
    dest->idx[i] = src->idx[i];
  }
}

void free_geneset_idx(geneset_index * g) {
  free(g->idx);
  free(g);
}
