typedef struct mat_type_s{
    unsigned int rows;
    unsigned int cols;
    double **data;
    int is_square;
} mat_type;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

mat_type* create_matrix(unsigned int num_rows, unsigned int num_cols);

mat_type* create_random_matrix(unsigned int num_rows, unsigned int num_cols, double min, double max);

mat_type* create_square_matrix(unsigned int size);

mat_type* create_eye_matrix(unsigned int size);

mat_type* create_matrix_from_file(FILE *file);

int check_matrix_dimension_equality(mat_type *matrix_a, mat_type *matrix_b);

int check_matrix_equality(mat_type *matrix_a, mat_type *matrix_b, double tolerance);

mat_type* get_matrix_column(mat_type *matrix, unsigned int column_index);

mat_type* get_matrix_row(mat_type *matrix, unsigned int row_index);

void mult_mat_row_scalar(mat_type *matrix, unsigned int row_index, double scalar);

void mult_mat_col_scalar(mat_type *matrix, unsigned int col_index, double scalar);

void swap_rows(mat_type *matrix, unsigned int row_index_1, unsigned int row_index_2);

void mat_add_rows(mat_type *matrix, unsigned int row_index_1, unsigned int row_index_2, double multiplier);

int find_pivot_row(mat_type *matrix, unsigned int column_idx, unsigned int row_idx, double precision);

void row_echelon_form(mat_type *matrix, double precision);

void reduced_row_echelon_form(mat_type *matrix, double precision);

int check_dims_matmul(mat_type *matrix_1, mat_type *matrix_2);

mat_type* mat_mul(mat_type *matrix_1, mat_type *matrix_2);

double mat_determinant(mat_type *matrix, double precision);

void print_matrix(mat_type *matrix);

void print_mem_adresses(mat_type *matrix);

void destroy_matrix(mat_type *matrix);