#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {
    rows_ = 10;
    cols_ = 10;
    matrix_ = new double*[rows_];
 
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[cols_];
    }
}

S21Matrix::S21Matrix(int rows, int cols) {
    rows_ = rows;
    cols_ = cols;
    matrix_ = new double*[rows_];
 
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[cols_];
    }
}

S21Matrix::S21Matrix(const S21Matrix& other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double*[rows_];
 
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[cols_];
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i][j] = other.matrix_[i][j];
        }
    }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
}
        
S21Matrix::~S21Matrix() {
    for (int i = 0; i < rows_; i++) {
        delete[] matrix_[i];
    }
    delete[] matrix_;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {

}

void S21Matrix::SumMatrix(const S21Matrix& other) {

}

void S21Matrix::SubMatrix(const S21Matrix& other) {

}

void S21Matrix::MulNumber(const double num) {

}

void S21Matrix::MulMatrix(const S21Matrix& other) {

}

S21Matrix S21Matrix::Transpose() {

}

S21Matrix S21Matrix::CalcComplements() {

}

double S21Matrix::Determinant() {

}

S21Matrix S21Matrix::InverseMatrix() {

}




int main () {
    // S21Matrix A;
    return 0;
}