#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "s21_matrix_oop.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

unsigned int seed = 0;
using namespace std;

TEST(TestGroupName, eq_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = rand_val;
            k += 0.000001;
        }
    }
    ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(TestGroupName, not_eq_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val;
            k += 0.000001;
        }
    }
    ASSERT_TRUE(A.EqMatrix(B) == FAILURE);
}

TEST(TestGroupName, not_rows_eq_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows + 10, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = i * 10.54;
            A(i, j) = rand_val;
        }
    }
    for (int i = 0; i < rows + 10; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = i * 10.54;
            B(i, j) = rand_val;
        }
    }
    ASSERT_TRUE(A.EqMatrix(B) == FAILURE);
}

TEST(TestGroupName, sum_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            Check_m(i, j) = 3 * rand_val + 0.015;
            k += 0.000001;
        }
    }
    A.SumMatrix(B);
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, sub_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            Check_m(i, j) = -rand_val - 0.015;
            k += 0.000001;
        }
    }
    A.SubMatrix(B);  
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, mult_number) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix Check_m(rows, cols);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            Check_m(i, j) = rand_val * 0.345;
            k += 0.000001;
        }
    }
    A.MulNumber(0.345);   
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, mult_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix B(cols, rows + 10);
    S21Matrix Check_m(rows, rows + 10);


    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            k += 0.000001;
        }
    }
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows + 10; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            B(i, j) = rand_val + k + 0.32;
            k += 0.000001;
        }
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows + 10; j++) {
            Check_m(i, j) = 0;
            for (int s = 0; s < cols; s++) {
                Check_m(i, j) += A(i, s) * B(s, j);
            }
        }
    }

    A.MulMatrix(B);
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, transpose_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix A_trans(cols, rows);
    S21Matrix Check_m(cols, rows);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            Check_m(j, i) = rand_val;
            k += 0.000001;
        }
    }
    A_trans = A.Transpose();   
    ASSERT_TRUE(A_trans.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, calc_complements_mat) {
    S21Matrix A(3, 3);
    S21Matrix A_calc(3, 3);
    S21Matrix check(3, 3);

    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 3;
    A(1, 0) = 0;
    A(1, 1) = 4;
    A(1, 2) = 2;
    A(2, 0) = 5;
    A(2, 1) = 2;
    A(2, 2) = 1;

    check(0, 0) = 0;
    check(0, 1) = 10;
    check(0, 2) = -20;
    check(1, 0) = 4;
    check(1, 1) = -14;
    check(1, 2) = 8;
    check(2, 0) = -8;
    check(2, 1) = -2;
    check(2, 2) = 4;

    A_calc = A.CalcComplements();   
    ASSERT_TRUE(A_calc.EqMatrix(check) == SUCCESS);
}

TEST(TestGroupName, determinant_mat) {
    S21Matrix A(3, 3);
    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 3;
    A(1, 0) = 0;
    A(1, 1) = 4;
    A(1, 2) = 2;
    A(2, 0) = 5;
    A(2, 1) = 2;
    A(2, 2) = 1;

    double result = A.Determinant();   
    ASSERT_TRUE(result == -40);
}

TEST(TestGroupName, inverse_mat) {
    S21Matrix A(3, 3);
    S21Matrix A_inverse(3, 3);
    S21Matrix check(3, 3);

    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 3;
    A(1, 0) = 0;
    A(1, 1) = 4;
    A(1, 2) = 2;
    A(2, 0) = 5;
    A(2, 1) = 2;
    A(2, 2) = 1;

    check(0, 0) = 0;
    check(0, 1) = -0.1;
    check(0, 2) = 0.2;
    check(1, 0) = -0.25;
    check(1, 1) = 0.35;
    check(1, 2) = 0.05;
    check(2, 0) = 0.5;
    check(2, 1) = -0.2;
    check(2, 2) = -0.1;

    A_inverse = A.InverseMatrix();   
    ASSERT_TRUE(A_inverse.EqMatrix(check) == SUCCESS);
}

TEST(TestGroupName, operator_plus) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            k += 0.000001;
        }
    }

    Check_m = A + B;
    A.SumMatrix(B);    
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, operator_minus) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            k += 0.000001;
        }
    }

    Check_m = A - B;
    A.SubMatrix(B);    
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, operator_mult_matrix) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix B(cols, rows + 10);
    S21Matrix Check_m(rows, rows + 10);


    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            k += 0.000001;
        }
    }
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows + 10; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            B(i, j) = rand_val + k + 0.32;
            k += 0.000001;
        }
    }

    Check_m = A * B;
    A.MulMatrix(B);    
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, operator_mult_number_1) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix Check_m(rows, cols);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            k += 0.000001;
        }
    }
    Check_m = A * 0.345;
    A.MulNumber(0.345);
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, operator_mult_number_2) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix Check_m(rows, cols);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            k += 0.000001;
        }
    }
    Check_m = 0.345 * A;
    A.MulNumber(0.345);    
    ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(TestGroupName, operator_eq_1) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = rand_val;
            k += 0.000001;
        }
    }   
    ASSERT_TRUE(A.EqMatrix(B) == (A == B));
}

TEST(TestGroupName, operator_eq_2) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val;
            k += 0.000001;
        }
    }   
    ASSERT_TRUE(A.EqMatrix(B) == (A == B));
}

TEST(TestGroupName, operator_plus_eq) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            Check_m(i, j) = 3 * rand_val + 0.015;
            k += 0.000001;
        }
    }
    S21Matrix Copy_A(A);
    Copy_A += B;
    A.SumMatrix(B);   
    ASSERT_TRUE((Copy_A == A) == SUCCESS);
}

TEST(TestGroupName, operator_minus_eq) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;

    S21Matrix A(rows, cols);
    S21Matrix B(rows, cols);
    S21Matrix Check_m(rows, cols);
    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            B(i, j) = 2 * rand_val + 0.015;
            Check_m(i, j) = 3 * rand_val + 0.015;
            k += 0.000001;
        }
    }
    S21Matrix Copy_A(A);
    Copy_A -= B;
    A.SubMatrix(B);   
    ASSERT_TRUE((Copy_A == A) == SUCCESS);
}

TEST(TestGroupName, operator_mult_number_eq) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix Check_m(rows, cols);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            Check_m(i, j) = rand_val * 0.345;
            k += 0.000001;
        }
    }
    S21Matrix Copy_A(A);
    Copy_A *= 0.345;
    A.MulNumber(0.345);   
    ASSERT_TRUE((A == Copy_A) == SUCCESS);
}

TEST(TestGroupName, operator_mult_matrix_eq) {
    const int rows = rand_r(&seed) % 100 + 1;
    const int cols = rand_r(&seed) % 100 + 1;
    S21Matrix A(rows, cols);
    S21Matrix B(cols, rows + 10);
    S21Matrix Check_m(rows, rows + 10);

    double k = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            A(i, j) = rand_val;
            k += 0.000001;
        }
    }
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows + 10; j++) {
            double rand_val = rand_r(&seed) % 2001 - 1000 + k;
            B(i, j) = rand_val + k + 0.32;
            k += 0.000001;
        }
    }

    S21Matrix Copy_A(A);
    Copy_A *= B;
    A.MulMatrix(B);    
    ASSERT_TRUE((A == Copy_A) == SUCCESS);
}


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
