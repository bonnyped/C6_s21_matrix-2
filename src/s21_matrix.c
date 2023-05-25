#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  if (rows > 0 && columns > 0) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
    }
  } else {
    res = 1;
  }

  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
  }
  zeroing_matrix_result(A);
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;

  if (A->rows < 1 || A->columns < 1 || B->rows < 1 || B->columns < 1)
    res = FAILURE;
  if (res && A->rows == B->rows && A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (!matrix_cell_comparison(A->matrix[i][j], B->matrix[i][j]))
          res = FAILURE;
      }
    }
  } else
    res = FAILURE;

  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  res = eq_size_matrix(A, B);
  if (!res) {
    if (!(res = s21_create_matrix(A->rows, A->columns, result))) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  }

  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  res = eq_size_matrix(A, B);
  if (!res) {
    if (!(res = s21_create_matrix(A->rows, A->columns, result))) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    }
  }

  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  if (!(res = s21_create_matrix(A->rows, A->columns, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }

  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  if (A->rows < 1 || B->columns < 1)
    res = 1;
  else if (!res && A->columns != B->rows)
    res = 2;
  else if (!res && !(res = s21_create_matrix(A->rows, B->columns, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }

  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = 0;
  zeroing_matrix_result(result);

  if (!(res = s21_create_matrix(A->columns, A->rows, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }

  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = 0;
  double resultat = 0.0;
  matrix_t minor_matrix = {0};
  zeroing_matrix_result(result);

  if (A->rows == 1 && A->columns == 1) res = 2;
  if (A->rows == A->columns && !res) {
    if (!(res = s21_create_matrix(A->rows, A->columns, result))) {
      s21_create_matrix(A->rows - 1, A->columns - 1, &minor_matrix);
      for (int row = 0; row < A->rows; row++) {
        for (int column = 0; column < A->columns; column++) {
          minor(row, column, A, &minor_matrix);
          s21_determinant(&minor_matrix, &resultat);
          if (resultat != 0) {
            result->matrix[row][column] =
                ((row + column) % 2 == 0 ? 1.0 : -1.0) * resultat;
          }
        }
      }
      s21_remove_matrix(&minor_matrix);
    }
  } else
    res = 2;

  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int res = 0;
  *result = 0.0;

  if (A->rows != A->columns) res = 2;
  if (A->rows < 1 || A->columns < 1) res = 1;
  if (!res) {
    if (A->rows == 1) {
      *result = **A->matrix;
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
    } else if (A->rows == 3) {
      *result = count_third_order_matrix_determinant(A);
    } else if (A->rows > 3) {
      *result = count_determinant(A);
    }
  }

  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = 0;
  double determinant = 0.0;
  zeroing_matrix_result(result);

  if (A->rows == 1 && A->columns == 1) {
    s21_create_matrix(A->rows, A->columns, result);
    s21_determinant(A, &determinant);
    if (fabs(determinant) > S21_EPSILON) {
      result->matrix[0][0] = 1 / determinant;
    } else
      res = 2;
  } else {
    if (!res && (A->rows < 1 || A->columns < 1)) res = 1;
    if (!res && !(res = s21_determinant(A, &determinant))) {
      if (A->rows == A->columns && fabs(determinant) > S21_EPSILON) {
        matrix_t tmp = {0};
        s21_calc_complements(A, result);
        s21_transpose(result, &tmp);
        s21_remove_matrix(result);
        determinant = 1 / determinant;
        s21_mult_number(&tmp, determinant, result);
        s21_remove_matrix(&tmp);
      } else
        res = 2;
    }
  }

  return res;
}
