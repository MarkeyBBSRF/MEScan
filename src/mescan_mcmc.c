#include "mescan.h"


void mcmc(int *r_mutation_data, int *r_mudata_data_dim,
          double *r_mutation_rate, int *r_mutation_rate_dim,
          int *r_num_iters, int *r_burning, double *r_mcmc_tuning, double *r_tg_lambda,
          double *r_tg_scores, char **r_geneset_file, int *r_geneset_max_tg,
          int *r_iter_max_tg, int *r_geneset_size, int *r_seed,
          int *r_debug, int *genes) {

  unsigned long seed;
  if (r_seed[0] == -1) {
    seed = get_seed(clock(), time(NULL), getpid());
    srand(seed);
  } else {
    seed = r_seed[0];
    srand(r_seed[0]);
  }

  clock_t begin = clock(), end;

  const short int SIZE = r_geneset_size[0];

  /* convert to matrix */
  mutation_data mut_data = vect_to_mutation_data(r_mutation_data, r_mudata_data_dim);
  mutation_rate mut_rate_data = vect_to_mutation_rate(r_mutation_rate, r_mutation_rate_dim);

  if (mut_data.nRow != mut_rate_data.nRow || mut_data.nCol != mut_rate_data.nCol) {
    perror("Dimensions differ in mut_data and mut_rate_data");
    exit(EXIT_FAILURE);
  }

  /* initilized mcmc params */
  const int num_genes = mut_data.nRow;
  double tg_lambda = r_tg_lambda[0];
  double mcmc_tuning = r_mcmc_tuning[0];
  int num_inters = r_num_iters[0];
  int burning = r_burning[0];
  int trans_counter = 0;
  int stay_counter = 0;
  int step = num_inters / 68;
  int steps = 0;

  /* 1-mut_rate_data */
  mutation_rate *complement_rate = get_complement_mut_rate(mut_rate_data);

  /* dump gene set index to file */
  FILE * geneset_filePtr;
  geneset_filePtr = fopen(*r_geneset_file, "a");
  if (geneset_filePtr == NULL) {
    perror(*r_geneset_file);
    exit (EXIT_FAILURE);
  }
  FILE * state_filePtr;
  char* stat = "stat.txt";
  state_filePtr = fopen(stat, "w");
  if (state_filePtr == NULL) {
    perror(stat);
    exit (EXIT_FAILURE);
  }
  fprintf(state_filePtr, "%lu\n", seed);
  fprintf(state_filePtr, "%d\n", SIZE);

  /****************/
  /*mcmc procedure*/
  /****************/
  /* initialized the gene set to start */

  int i, j, next_gene, max_tg_iter = -1, to_replace;
  int *curr_geneset = malloc(num_genes * sizeof(int));
  int *next_geneset = malloc(num_genes * sizeof(int));
  int *track = malloc(num_genes * sizeof(int));
  double tg_score, tg_curr, tg_next, max_tgscore, rate, randrate;

  int iter = 0;
  memset(curr_geneset, 0, num_genes * sizeof(int));
  memset(next_geneset, 0, num_genes * sizeof(int));
  memset(track, 0, num_genes * sizeof(int));

  geneset_index* curr_geneset_idx = (geneset_index *)malloc(sizeof(geneset_index));
  geneset_index* next_geneset_idx  = (geneset_index *)malloc(sizeof(geneset_index));

  if (genes[0] != -1) {
    curr_geneset_idx = init_index(genes, SIZE);
    tg_next = get_tg_score(mut_rate_data, complement_rate, mut_data, curr_geneset_idx, tg_lambda);
  } else {
    do {
      memset(curr_geneset, 0, num_genes * sizeof(int));
      memset(next_geneset, 0, num_genes * sizeof(int));
      do {
        next_gene = sample_next_gene(0, num_genes - 1);
        next_geneset_idx = get_next_geneset_idx(next_geneset, next_gene , num_genes, -1);
      } while (next_geneset_idx == NULL || next_geneset_idx->len < SIZE);
      tg_next = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);

    }
    while (tg_next < 0);
    printf("TC_init = %.3f\n", tg_next);
    for (i = 0; i < next_geneset_idx->len; i++) {
      track[next_geneset_idx->idx[i]] = 1;
    }
    tg_next = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);
    max_tgscore = tg_next;
    // update curr_geneset
    memcpy(curr_geneset, next_geneset, num_genes * sizeof(int));
    update_geneset_idx(curr_geneset_idx, next_geneset_idx);

    if (iter > (burning - 1)) {
      if (r_debug[0] == 1) {
        fprintf(geneset_filePtr, "%.3f\t", tg_next);
      }
      fprintf(geneset_filePtr, "%.3f\t", tg_next);

      write_file_int_vect(curr_geneset_idx->idx, curr_geneset_idx->len, geneset_filePtr);
    }
    free_geneset_idx(next_geneset_idx);
  }

  tg_curr = tg_next;

  /* For each iteration, sample a gene G
  	1. if G in current_set, resample G;
  	2. if G not in current_set, randomly sample a gene G' in current_set, and
  	replace G' with G.
  */
  for (iter = 1; iter < num_inters; iter++) {
    memcpy(next_geneset, curr_geneset, num_genes * sizeof(int));
    // find a gene to replace
    do {
      next_gene = sample_next_gene(0, num_genes - 1);
    } while (track[next_gene] == 1);

    to_replace = sample_next_gene(0, curr_geneset_idx->len - 1);
    to_replace = curr_geneset_idx->idx[to_replace];
    next_geneset_idx = get_next_geneset_idx(next_geneset, next_gene, num_genes, to_replace);
    tg_score = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);

    rate = tg_score / tg_curr;
    rate = pow(rate, mcmc_tuning);
    randrate = flipcoin();


    if (rate < 1 && rate < randrate) {
      /* stay */
      stay_counter++;
      if (r_debug[0] == 1) {
        fprintf(geneset_filePtr, "%s\t", "stay");
      }
      track[next_gene] = 0;

    } else {
      /* trans */
      memcpy(curr_geneset, next_geneset, num_genes * sizeof(int));
      update_geneset_idx(curr_geneset_idx, next_geneset_idx);
      tg_curr = tg_score;
      trans_counter++;
      if (r_debug[0] == 1) {
        fprintf(geneset_filePtr, "%s\t", "tran");
      }
      track[next_gene] = 1;
      track[to_replace] = 0;
    }

    if (iter >= burning) {
      if (r_debug[0] == 1) {
        fprintf(geneset_filePtr, "%.3f\t%.3f\t%.3f\t", rate, randrate, tg_score);
      }
      fprintf(geneset_filePtr, "%.3f\t", tg_curr);
      write_file_int_vect(curr_geneset_idx->idx, curr_geneset_idx->len, geneset_filePtr);
      if (tg_curr > max_tgscore) {
        max_tgscore = tg_curr;
        max_tg_iter = iter;
        memcpy(r_geneset_max_tg, curr_geneset, num_genes * sizeof(int));
      }
    }

    if (iter % step == 0) {
      steps = (int) 67. * iter / num_inters;
      for (j = 0; j < steps; j++) printf("+");
      for (j = 0; j < 67 - steps; j++) printf(" ");
      printf(" [ %d%% ]\r", (int) (100. * iter / num_inters));
      fflush(stdout);
    }

    free_geneset_idx(next_geneset_idx);
  }
  end = clock();

  write_file_int_vect(curr_geneset_idx->idx, curr_geneset_idx->len, state_filePtr);

  for (j = 0; j < 67; j++) printf("+");
  printf(" [ 100%% ]\n");
  printf("MCMC parameters:\n  - geneset_size = %d\n  - lambda = %.3f\n  - iteration = %d\n  - burning = %d\n  - mcmc_tuning = %.3f\n\n",
         SIZE, tg_lambda, r_num_iters[0], r_burning[0], r_mcmc_tuning[0]);
  printf("  - Max tg_score = %.3f at iteration: %d\n", max_tgscore, max_tg_iter);
  memcpy(r_iter_max_tg, &max_tg_iter, sizeof(int));

  printf("  - runtime: %.3f mins\n", (double)(end - begin) / CLOCKS_PER_SEC / 60);
  printf("  - Stay = %d\n", stay_counter);
  printf("  - Trans = %d\n", trans_counter);


  free(next_geneset);
  free(curr_geneset);
  fclose(geneset_filePtr);
  fclose(state_filePtr);

}


void select_tuning(int *r_mutation_data, int *r_mudata_data_dim,
                   double *r_mutation_rate, int *r_mutation_rate_dim,
                   int *r_num_iters, double *r_mcmc_tuning, double *r_tg_lambda,
                   int *r_seed, int *min_gs, int *max_gs) {
  mutation_data mut_data = vect_to_mutation_data(r_mutation_data, r_mudata_data_dim);
  mutation_rate mut_rate_data = vect_to_mutation_rate(r_mutation_rate, r_mutation_rate_dim);
  mutation_rate *complement_rate = get_complement_mut_rate(mut_rate_data);
  const int num_genes = mut_data.nRow;
  double tg_lambda = r_tg_lambda[0];


  int num_inters = r_num_iters[0];
  int s = 0;
  double ratio = 100;
  for (s = min_gs[0]; s < (max_gs[0] + 1); s++) {
    int mcmc_tuning = r_mcmc_tuning[0];
    ratio = 100;
    while (ratio > 35.0 || ratio < 20.0) {
      //printf("gs = %d, tuning = %.3f, ratio = %.3f\n",s, mcmc_tuning, ratio);
      if (ratio > 35.0)
        mcmc_tuning = mcmc_tuning + 1;
      if (ratio < 20.0)
        mcmc_tuning = mcmc_tuning - 1;;

      int trans_counter = 0;
      int stay_counter = 0;
      unsigned long seed;
      if (r_seed[0] == -1) {
        seed = get_seed(clock(), time(NULL), getpid());
        srand(seed);
      } else {
        srand(r_seed[0]);
      }
      int SIZE = s;
      int i, next_gene, to_replace;
      int *curr_geneset = malloc(num_genes * sizeof(int));
      int *next_geneset = malloc(num_genes * sizeof(int));
      int *track = malloc(num_genes * sizeof(int));
      double *tg_arr = malloc(num_inters * sizeof(double));
      double tg_score, rate, randrate;

      int iter = 0;
      memset(curr_geneset, 0, num_genes * sizeof(int));
      memset(next_geneset, 0, num_genes * sizeof(int));
      memset(track, 0, num_genes * sizeof(int));

      geneset_index* curr_geneset_idx = (geneset_index *)malloc(sizeof(geneset_index));
      geneset_index* next_geneset_idx  = (geneset_index *)malloc(sizeof(geneset_index));

      /* Initialized MCMC change to the SIZE */
      do {
        memset(curr_geneset, 0, num_genes * sizeof(int));
        memset(next_geneset, 0, num_genes * sizeof(int));
        do {
          next_gene = sample_next_gene(0, num_genes - 1);
          next_geneset_idx = get_next_geneset_idx(next_geneset, next_gene , num_genes, -1);
        } while (next_geneset_idx == NULL || next_geneset_idx->len < SIZE);
        tg_arr[iter] = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);
        // printf("TC_init = %.3f\n",tg_arr[iter]);
      }
      while (tg_arr[iter] < 0);
      for (i = 0; i < next_geneset_idx->len; i++) {
        track[next_geneset_idx->idx[i]] = 1;
      }


      tg_arr[iter] = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);

      // update curr_geneset
      memcpy(curr_geneset, next_geneset, num_genes * sizeof(int));
      update_geneset_idx(curr_geneset_idx, next_geneset_idx);

      free_geneset_idx(next_geneset_idx);

      for (iter = 1; iter < num_inters; iter++) {
        memcpy(next_geneset, curr_geneset, num_genes * sizeof(int));
        // find a gene to replace
        do {
          next_gene = sample_next_gene(0, num_genes - 1);
          //printf("nextgene = %d\n",next_gene);
        } while (track[next_gene] == 1);

        to_replace = sample_next_gene(0, curr_geneset_idx->len - 1);
        to_replace = curr_geneset_idx->idx[to_replace];
        next_geneset_idx = get_next_geneset_idx(next_geneset, next_gene, num_genes, to_replace);
        tg_score = get_tg_score(mut_rate_data, complement_rate, mut_data, next_geneset_idx, tg_lambda);
        //printf("tg=%.3f\n",tg_score);
        // rate = mcmc_tuning*(tg_score/tg_arr[iter - 1]);
        // rate = pow(rate,4);
        // adjust the tunning in the burning part,check the tran ratio
        // for every 10000 iteration.

        rate = tg_score / tg_arr[iter - 1];
        rate = pow(rate, mcmc_tuning);
        randrate = flipcoin();


        if (rate < 1 && rate < randrate) {
          /* stay */
          // MIY - current is curr_geneset
          // memcpy(MIY,curr_geneset,num_genes*sizeof(int));
          tg_arr[iter] = tg_arr[iter - 1];
          stay_counter++;

          track[next_gene] = 0;

        } else {
          /* trans */
          // MIY - current is next_geneset
          memcpy(curr_geneset, next_geneset, num_genes * sizeof(int));
          update_geneset_idx(curr_geneset_idx, next_geneset_idx);
          tg_arr[iter] = tg_score;
          trans_counter++;

          track[next_gene] = 1;
          track[to_replace] = 0;
        }
        //printf("tran = %d, stay = %d",trans_counter,stay_counter);
        ratio = trans_counter * 1.0 / (trans_counter + stay_counter) * 100.0;

        free_geneset_idx(next_geneset_idx);


      }
      free(next_geneset);
      free(curr_geneset);

    }
    printf("%d %d %.3f\n", s, mcmc_tuning, ratio);
  }

}

