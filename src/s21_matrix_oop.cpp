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
    bool error = SUCCESS;
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        error = FAILURE;
    }

    for (int i = 0; error != FAILURE && i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            if (fabs(matrix_[i][j] - other.matrix_[i][j]) >= 1e-8) {
                error = FAILURE;
                break;
            }
        }
    }

    return error;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::exception();
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
        }
    }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::exception();
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
        }
    }
}

void S21Matrix::MulNumber(const double num) {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] * num;
        }
    }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
    if (cols_ != other.rows_) {
        throw std::exception();
    }

    S21Matrix result(rows_, other.cols_);

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < other.cols_; j++) {
            result.matrix_[i][j] = 0;
            for (int k = 0; k < cols_; k++) {
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }

    rows_ = result.rows_;
    cols_ = result.cols_;
    matrix_ = new double*[result.rows_];
    for (int i = 0; i < result.rows_; i++) {
        matrix_[i] = new double[result.cols_];
        for (int j = 0; j < result.cols_; j++) {
            matrix_[i][j] = result.matrix_[i][j];
        }
    }

    //  проверить на течку
}

S21Matrix S21Matrix::Transpose() {
    S21Matrix matrix_trans(cols_, rows_);

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_trans.matrix_[j][i] = matrix_[i][j];
        }
    }

    return matrix_trans;
}

void calc_minor(double** matrix, int rows, int cols, int i_remove, int j_remove, S21Matrix* minor) {
    int k = 0, s;
    for (int i = 0; i < rows; i++) {
        if (i != i_remove) {
            s = 0;
            for (int j = 0; j < cols; j++) {
                if (j != j_remove) {
                    minor->matrix_[k][s] = matrix[i][j];
                    s++;
                }
            }
            k++;
        }
    }
}


S21Matrix S21Matrix::CalcComplements() {
    // if (rows_ != cols_) {
    //     throw std::exception();
    // }

    // S21Matrix result(rows_, cols_);
}

double S21Matrix::Determinant() {
    if (rows_ != cols_) {
        throw std::exception();
    }

    double det = 0;

    if (rows_ == 1) {
        det = matrix_[0][0];
    } else  if (rows_ == 2) {
        det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else {
        S21Matrix minor(rows_ - 1, cols_ - 1);
        int k = 1;
        for (int i = 0; i < rows_; i++) {
            calc_minor(matrix_, rows_, cols_, i, 0, &minor);
            det += k * matrix_[i][0] * minor.Determinant();
            k *= -1;
        }
        minor.~S21Matrix();  //  ???
    }

    return det;

}

S21Matrix S21Matrix::InverseMatrix() {

}

int main() {
  // S21Matrix A;
  return 0;
}