#include <stdio.h>
#include "lib/linalg.h"

int main(){
    
    FILE *matrix_file = fopen("matrix.data", "r");
    mat_type *test_matrix = create_matrix_from_file(matrix_file);
    printf("%f\n", mat_determinant(test_matrix, 0001));
    destroy_matrix(test_matrix);
    // mat_type *eye_matrix_4 = create_eye_matrix(4);
    // mat_type *test_column = get_matrix_column(test_matrix, 1);
    // mat_type *test_row = get_matrix_row(test_matrix,2);
    // mat_type *test_matmul = mat_mul(test_matrix, eye_matrix_4);
    // mat_type *test_matmul_2 = mat_mul(test_matrix, test_column);
    
    // printf("original matrix\n");
    // print_matrix(test_matrix);
    // printf("\n");

    // printf("matrix after row addition\n");
    // mat_add_rows(test_matrix, 0,1,2.3);
    // print_matrix(test_matrix);
    // printf("\n");

    // print_matrix(test_matmul);

    // printf("\n");
    // print_matrix(test_matmul_2);
    // printf("\n");

    // print_matrix(test_column);
    // printf("\n");
    // print_matrix(test_row);
    // printf("The matrices are equal: %i\n", check_matrix_equality(test_matrix,test_matrix, 0.0));
    // printf("The matrices are equal: %i\n", check_matrix_equality(test_matrix,test_column, 0.0));
    
    // destroy_matrix(eye_matrix_4);
    // destroy_matrix(test_matmul);
    // destroy_matrix(test_matmul_2);
    // destroy_matrix(test_column);
    // destroy_matrix(test_row);
    
    return 0;
}