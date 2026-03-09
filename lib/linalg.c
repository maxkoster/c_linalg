#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "linalg.h"

double random_double(double min, double max){
    double range = max - min;
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

mat_type* create_matrix(unsigned int n_rows, unsigned int n_cols){
    if (n_rows == 0) {
        printf("Cannot create a matrix with 0 rows");
        return NULL;
    }
    if (n_cols == 0) {
        printf("Cannot create a matrix with 0 columns");
        return NULL;
    }
    mat_type *matrix = calloc(1,sizeof(*matrix));
    matrix->rows = n_rows;
    matrix->cols = n_cols;
    matrix->is_square = (n_cols == n_rows) ? 1 : 0;
    matrix->data = calloc(matrix->rows,sizeof(*matrix->data));
    for (int i = 0; i < matrix->rows; i++)
    {
        matrix->data[i] = calloc(matrix->cols, sizeof(**matrix->data));
    }
    
    return matrix;
}

mat_type* create_random_matrix(unsigned int num_rows, unsigned int num_cols, double min, double max){
    mat_type *random_matrix = create_matrix(num_rows, num_cols);

    srand(time(NULL));
    // Set values in matrix
    for (int i = 0; i < random_matrix->rows; i++){
        for (int j = 0; j < random_matrix->cols; j++){
            random_matrix->data[i][j] = random_double(min, max);
        }
    }

    return random_matrix;
}

mat_type* create_square_matrix(unsigned int size){
    mat_type *square_matrix = create_matrix(size,size);
    return square_matrix;
}

mat_type* create_eye_matrix(unsigned int size){
    mat_type *eye_matrix = create_square_matrix(size);

    for (int i = 0; i < eye_matrix->rows; i++){
        eye_matrix->data[i][i] = 1;
    }
    return eye_matrix;
}

mat_type* create_matrix_from_file(FILE *file){
    unsigned int num_rows = 0, num_cols = 0;
    fscanf(file, "%d", &num_rows);
    fscanf(file, "%d", &num_cols);
    mat_type *matrix = create_matrix(num_rows, num_cols);
    for(int i = 0; i < matrix->rows; i++) {
        for(int j = 0; j < num_cols; j++) {
            fscanf(file, "%lf\t", &matrix->data[i][j]);
    }
  }
  return matrix;
}

int check_matrix_dimension_equality(mat_type *matrix_a, mat_type *matrix_b){
    return (matrix_a->rows == matrix_b->rows) && 
        (matrix_a->cols == matrix_b->cols);
}

int check_matrix_equality(mat_type *matrix_a, mat_type *matrix_b, double tolerance){
    
    if (!check_matrix_dimension_equality(matrix_a, matrix_b)){
        printf("Matrix A of shape (%i,%i) and Matrix B of shape (%i,%i) have different dimensions\n",
            matrix_a->rows, matrix_a->cols, matrix_b->rows, matrix_b->cols);
        return 0;
    }

    for (int i = 0; i < matrix_a->rows; i++){
        for (int j = 0; j < matrix_a->cols; j++){
            if (fabs(matrix_a->data[i][j] - matrix_b->data[i][j]) > tolerance){
                return 0;
            }
        }
    }

    return 1;
}

mat_type* get_matrix_column(mat_type *matrix, unsigned int column_index){
    mat_type *column = create_matrix(matrix->rows, 1);

    if (column_index >= matrix->cols) {
        printf("Selected column index %i is out of range. Value should be between 0 and %i", column_index, matrix->cols);
        return column;
    }

    for (int i = 0; i < matrix->rows; i++) {
        column->data[i][0] = matrix->data[i][column_index];
    }

    return column;
}

mat_type* get_matrix_row(mat_type *matrix, unsigned int row_index){
    mat_type *row = create_matrix(1,matrix->cols);

    if (row_index >= matrix->rows) {
        printf("Selected row index %i is out of range. Value should be between 0 and %i", row_index, matrix->rows);
        return row;
    }

    memcpy(row->data[0], matrix->data[row_index], row->cols * sizeof(matrix->data[row_index]));
    return row;
}

void mult_mat_row_scalar(mat_type *matrix, unsigned int row_index, double scalar){
    for (int i = 0; i < matrix->cols; i++){
        matrix->data[row_index][i] *=  scalar;
    }
}

void mult_mat_col_scalar(mat_type *matrix, unsigned int column_index, double scalar){
    for (int i = 0; i < matrix->rows; i++){
        matrix->data[i][column_index] *=  scalar;
    }
}

void swap_rows(mat_type *matrix, unsigned int row_index_1, unsigned int row_index_2){
    double temp_row[matrix->rows];
    memcpy(temp_row, matrix->data[row_index_1], matrix->cols * sizeof(matrix->data[0]));
    memcpy(matrix->data[row_index_1], matrix->data[row_index_2], matrix->cols * sizeof(matrix->data[0]));
    memcpy(matrix->data[row_index_2], temp_row, matrix->cols * sizeof(matrix->data[0]));
}

// Add values of row at index row_index_1 to row at row_index_2 with a multiplier
void mat_add_rows(mat_type *matrix, unsigned int row_index_1, unsigned int row_index_2, double multiplier){
    for (int i = 0; i < matrix->cols; i++){
        matrix->data[row_index_2][i] += multiplier * matrix->data[row_index_1][i];
    }
}

// Find row with pivot point, for stability takes largest number from column.
int find_pivot_row(mat_type *matrix, unsigned int column_idx, unsigned int row_idx, double precision){
    int pivot_idx = -1;
    double max_pivot = precision;
    for (int i = row_idx; i < matrix->rows; i++){
        if ( fabs(matrix->data[i][column_idx]) > max_pivot ){
            pivot_idx = i;
            max_pivot = matrix->data[i][column_idx];
        }
    }
    return pivot_idx;
}

void row_echelon_form(mat_type *matrix, double precision){
    int pivot_idx;
    int i = 0, j = 0;
    while (i < matrix->rows && j < matrix->cols){
        pivot_idx = find_pivot_row(matrix, j, i, precision);
        if (pivot_idx < 0){
            j++;
            continue;
        }
        mult_mat_row_scalar(matrix, pivot_idx, 1/(matrix->data[pivot_idx][j]));
        if (pivot_idx != i){
            swap_rows(matrix, pivot_idx, i);                        
        }

        // subtract row at pivot from other rows
        for (int row_index = i + 1; row_index < matrix->rows; row_index++){
            if (matrix->data[row_index][j] > precision || matrix->data[row_index][j] < -precision){
                mat_add_rows(matrix, i, row_index, - matrix->data[row_index][j]);
            }
        }
        i++;
        j++;
    }
}

void reduced_row_echelon_form(mat_type *matrix, double precision){
    int pivot_idx;
    int i = 0, j = 0;
    while (i < matrix->rows && j < matrix->cols){
        pivot_idx = find_pivot_row(matrix, j, i, precision);
        if (pivot_idx < 0){
            j++;
            continue;
        }
        mult_mat_row_scalar(matrix, pivot_idx, 1/(matrix->data[pivot_idx][j]));
        if (pivot_idx != i){
            swap_rows(matrix, pivot_idx, i);                        
        }

        // subtract row at pivot from other rows
        for (int row_index = 0; row_index < matrix->rows; row_index++){
            if (row_index == i){
                continue;
            }
            mat_add_rows(matrix, i, row_index, - matrix->data[row_index][j]);
        }
        i++;
        j++;
    }
}

mat_type* transpose_mat(mat_type *matrix){
    mat_type *out_matrix = create_matrix(matrix->cols, matrix->rows);
    return out_matrix;
}

int check_dims_matmul(mat_type *matrix_1, mat_type *matrix_2){
    if (matrix_1->cols != matrix_2->rows){
        return 0;
    }
    return 1;
}

mat_type* mat_mul(mat_type *matrix_1, mat_type *matrix_2){
    if (!check_dims_matmul(matrix_1, matrix_2)){
        printf("Dimensions do not match nr. of cols of matrix 1: %i, does not match nr. of rows of matrix 2: %i", matrix_1->cols, matrix_2->rows);
        return NULL;
    }
    mat_type *matrix_r = create_matrix(matrix_1->rows, matrix_2->cols);
    
    for (int col = 0; col < matrix_2->cols; col++){
        for (int row = 0; row < matrix_1->rows; row++){
            for (int i = 0; i < matrix_2->rows; i++){
                matrix_r->data[row][col] += matrix_1->data[row][i] * matrix_2->data[i][col];
            }
        }
    }
    return matrix_r;
}

double mat_determinant(mat_type *matrix, double precision){
    // Create upper diagonal matrix
    // Swapping rows = det * -1
    if (matrix->rows != matrix->cols){
        printf("Matrix is not a square matrix, dimensions are %i x %i\n", matrix->rows, matrix->cols);
        exit(1);
    }
    int cofactor = 1;
    double determinant = 1;
    for (int i = 0; i < matrix->rows; i++){
        int pivot_idx = find_pivot_row(matrix, i, i, precision);

        if (pivot_idx < 0){
            return 0;
        }
        
        if (pivot_idx != i){
            cofactor *= -1;
            swap_rows(matrix, pivot_idx, i);
        }

        determinant *= matrix->data[i][i];

        // subtract row at pivot from other rows
        for (int row_index = i + 1; row_index < matrix->rows; row_index++){
            double coefficient = matrix->data[row_index][i]/matrix->data[i][i];
            mat_add_rows(matrix, i, row_index, - coefficient);
        }
    }
    return cofactor * determinant;
}

void print_matrix(mat_type *matrix){
    for (int i = 0; i < matrix->rows; i++){     
        printf("|");
        for (int j = 0; j < matrix->cols; j++){
            if (j == matrix->cols-1){
                printf("%f|", matrix->data[i][j]);
            }
            else{
                printf("%f  ", matrix->data[i][j]);
            }
        }
        printf("\n");
    }
}


void destroy_matrix(mat_type *matrix){
    for (int i = 0; i < matrix->rows; i++){
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}

void print_mem_adresses(mat_type *matrix){
    for (int i = 0; i< matrix->rows; i++){
        printf("%p\n", &matrix->data[i]);
        for (int j = 0; j < matrix->cols; j++){
            printf("%i %p ", i, &matrix->data[i][j]);
        }
        printf("\n");
    }
}