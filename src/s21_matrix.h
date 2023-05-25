#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SUCCESS 1
#define FAILURE 0
#define S21_EPSILON 1e-8
#define FRACTIONAL_BUFFER_SIZE 11

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

void zeroing_matrix_result(matrix_t *result);
int matrix_cell_comparison(double A, double B);
int eq_fractional_part(char *A, char *B);
double count_determinant(matrix_t *A);
void minor(int r, int c, matrix_t *A, matrix_t *B);
double count_third_order_matrix_determinant(matrix_t *A);
void inverse(matrix_t *A);
int eq_size_matrix(matrix_t *A, matrix_t *B);

#endif  // S21_MATRIX_H