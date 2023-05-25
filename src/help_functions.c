#include "s21_matrix.h"

void zeroing_matrix_result(matrix_t *result) {
  result->matrix = NULL;
  result->rows = 0;
  result->columns = 0;
}

double count_determinant(matrix_t *A) {
  double res = 0.0;
  double tmp = 1.0;
  double calc = 0.0;
  matrix_t B;
  s21_create_matrix(A->rows - 1, A->rows - 1, &B);

  for (int row = 0; row < A->rows; row++) {
    if (A->matrix[row][0]) {
      minor(row, 0, A, &B);
      s21_determinant(&B, &tmp);
      calc = row % 2 == 0 ? 1.0 : -1.0;
      res += calc * A->matrix[row][0] * tmp;
    }
  }

  s21_remove_matrix(&B);

  return res;
}

double count_third_order_matrix_determinant(matrix_t *A) {
  double res = 0.0;

  res += A->matrix[0][0] * A->matrix[1][1] * A->matrix[2][2];
  res += A->matrix[0][1] * A->matrix[1][2] * A->matrix[2][0];
  res += A->matrix[0][2] * A->matrix[1][0] * A->matrix[2][1];
  res -= A->matrix[0][2] * A->matrix[1][1] * A->matrix[2][0];
  res -= A->matrix[1][0] * A->matrix[0][1] * A->matrix[2][2];
  res -= A->matrix[0][0] * A->matrix[2][1] * A->matrix[1][2];

  return res;
}

int eq_fractional_part(char *A, char *B) {
  int res = 1;

  for (int i = 0; i < FRACTIONAL_BUFFER_SIZE; i++) {
    A[i] != B[i] ? res = 0, i = FRACTIONAL_BUFFER_SIZE : res;
  }

  return res;
}

int matrix_cell_comparison(double A, double B) {
  int res = 1;
  double A_plus_accuracy = A + S21_EPSILON;
  double B_plus_accuracy = B + S21_EPSILON;
  int whole_part_A = (int)A_plus_accuracy;
  int whole_part_B = (int)B_plus_accuracy;
  double fractional_part_A = A - (double)whole_part_A + S21_EPSILON;
  double fractional_part_B = B - (double)whole_part_B + S21_EPSILON;
  char fractional_A_to_string[FRACTIONAL_BUFFER_SIZE] = {0};
  char fractional_B_to_string[FRACTIONAL_BUFFER_SIZE] = {0};

  whole_part_A != whole_part_B ? res = 0 : res;
  if (res) {
    sprintf(fractional_A_to_string, "%.7f", fractional_part_A);
    sprintf(fractional_B_to_string, "%.7f", fractional_part_B);
    res = eq_fractional_part(fractional_A_to_string, fractional_B_to_string);
  }

  return res;
}

void minor(int r, int c, matrix_t *A, matrix_t *B) {
  int count_row = 0;
  int count_column = 0;

  for (int row = 0; row < A->rows; row++) {
    for (int column = 0; column < A->columns; column++) {
      if (row != r && column != c) {
        B->matrix[count_row][count_column] = A->matrix[row][column];
        count_column++;
      }
    }
    if (row != r) {
      count_row++;
      count_column = 0;
    }
  }
}

int eq_size_matrix(matrix_t *A, matrix_t *B) {
  int res = 0;

  if (A->rows < 1 || A->columns < 1) res = 1;
  if (!res && (B->rows < 1 || B->columns < 1)) res = 1;
  if (!res && (A->rows != B->rows)) res = 2;
  if (!res && (A->columns != B->columns)) res = 2;

  return res;
}