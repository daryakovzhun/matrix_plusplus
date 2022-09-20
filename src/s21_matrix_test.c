#include <check.h>

#include "s21_matrix_oop.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

unsigned int seed = 0;

START_TEST(eq_matrix) {

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
    ck_assert_int_eq(A.EqMatrix(B), SUCCESS);
}
END_TEST

START_TEST(not_eq_matrix) {

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
    ck_assert_int_eq(A.EqMatrix(B), FAILURE);
}
END_TEST

START_TEST(not_rows_eq_matrix) {

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
    ck_assert_int_eq(A.EqMatrix(B), FAILURE);
}
END_TEST

START_TEST(sum_matrix) {

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
    ck_assert_int_eq(A.EqMatrix(Check_m), SUCCESS);
}
END_TEST

START_TEST(sub_matrix) {

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
    ck_assert_int_eq(A.EqMatrix(Check_m), SUCCESS);
}
END_TEST

START_TEST(mult_number) {
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

    ck_assert_int_eq(A.EqMatrix(Check_m), SUCCESS);
}
END_TEST

START_TEST(mult_matrix) {
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
    ck_assert_int_eq(A.EqMatrix(Check_m), SUCCESS);
}
END_TEST

START_TEST(transpose_matrix) {
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

    ck_assert_int_eq(A_trans.EqMatrix(Check_m), SUCCESS);
}
END_TEST

START_TEST(calc_complements_mat) {
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

    ck_assert_int_eq(A_calc.EqMatrix(check), SUCCESS);
}
END_TEST

START_TEST(determinant_mat) {
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
    ck_assert_int_eq(result, -40);
}
END_TEST

START_TEST(inverse_mat) {
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

    ck_assert_int_eq(A_inverse.EqMatrix(check), SUCCESS);
}
END_TEST


int main(void) {
    Suite *s1 = suite_create("Test_matrix");
    TCase *tc1 = tcase_create("Test_matrix");
    SRunner *sr = srunner_create(s1);
    int nf;

    suite_add_tcase(s1, tc1);
    tcase_add_loop_test(tc1, eq_matrix, 0, 100);
    tcase_add_loop_test(tc1, not_eq_matrix, 0, 100);
    tcase_add_loop_test(tc1, not_rows_eq_matrix, 0, 100);
    tcase_add_loop_test(tc1, sum_matrix, 0, 100);
    tcase_add_loop_test(tc1, sub_matrix, 0, 100);
    tcase_add_loop_test(tc1, mult_number, 0, 100);
    tcase_add_loop_test(tc1, mult_matrix, 0, 100);
    tcase_add_loop_test(tc1, transpose_matrix, 0, 100);
    tcase_add_loop_test(tc1, calc_complements_mat, 0, 100);
    tcase_add_loop_test(tc1, determinant_mat, 0, 100);
    tcase_add_loop_test(tc1, inverse_mat, 0, 100);


    srunner_set_fork_status(sr, CK_NOFORK);
    srunner_run_all(sr, CK_ENV);
    nf = srunner_ntests_failed(sr);
    srunner_free(sr);

    return nf == 0 ? 0 : 1;
}
