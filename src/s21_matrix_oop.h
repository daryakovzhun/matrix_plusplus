#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <exception>

#define SUCCESS 1
#define FAILURE 0

class S21Matrix {
    friend void calc_minor(double** matrix, int rows, int cols, int i_remove, int j_remove, S21Matrix* minor);

    public:
        S21Matrix();
        S21Matrix(int rows, int cols);
        S21Matrix(const S21Matrix& other);
        S21Matrix(S21Matrix&& other);

        ~S21Matrix();

        bool EqMatrix(const S21Matrix& other);
        void SumMatrix(const S21Matrix& other);
        void SubMatrix(const S21Matrix& other);
        void MulNumber(const double num);
        void MulMatrix(const S21Matrix& other);
        S21Matrix Transpose();
        S21Matrix CalcComplements();
        double Determinant();
        S21Matrix InverseMatrix();

    private:
        int rows_, cols_;         
        double **matrix_;
};

#endif  //  SRC_S21_MATRIX_H_
