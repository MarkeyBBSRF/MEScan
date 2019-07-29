#include "mescan.h"

unsigned long get_seed(unsigned long a, unsigned long b, unsigned long c)
{
  a = a - b;  a = a - c;  a = a ^ (c >> 13);
  b = b - c;  b = b - a;  b = b ^ (a << 8);
  c = c - a;  c = c - b;  c = c ^ (b >> 13);
  a = a - b;  a = a - c;  a = a ^ (c >> 12);
  b = b - c;  b = b - a;  b = b ^ (a << 16);
  c = c - a;  c = c - b;  c = c ^ (b >> 5);
  a = a - b;  a = a - c;  a = a ^ (c >> 3);
  b = b - c;  b = b - a;  b = b ^ (a << 10);
  c = c - a;  c = c - b;  c = c ^ (b >> 15);
  return c;
}

double flipcoin()
{
  return rand() / (RAND_MAX + 1.);
}

int rand_int()
{
  return rand() % 100 + 1;
  // return (int)(100.0*rand())/(RAND_MAX+1.0) + 1;
}


void write_file_int_vect(int * vec, int size, FILE *file) {
  int i;
  for (i = 0; i < size; i++) {
    fprintf(file, "%d ", vec[i]);
  }
  fputs("\n", file);
}

