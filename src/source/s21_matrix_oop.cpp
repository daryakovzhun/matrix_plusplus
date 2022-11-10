#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : S21Matrix(10, 10) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  matrix_ = new double *[rows_];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double *[rows_];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    memcpy(matrix_[i], other.matrix_[i], cols_ * sizeof(double));
  }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}

int S21Matrix::getRows() { return rows_; }

int S21Matrix::getCols() { return cols_; }

void S21Matrix::setRows(const int rows) {
  S21Matrix buffer(*this);
  this->~S21Matrix();
  matrix_ = new double *[rows];
  rows_ = rows;
  cols_ = buffer.cols_;

  for (int i = 0; i < rows; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) {
      if (i < buffer.rows_) {
        matrix_[i][j] = buffer.matrix_[i][j];
      } else {
        matrix_[i][j] = 0;
      }
    }
  }
}

void S21Matrix::setCols(const int cols) {
  S21Matrix buffer(*this);
  this->~S21Matrix();
  matrix_ = new double *[buffer.rows_];
  rows_ = buffer.rows_;
  cols_ = cols;

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) {
      if (j < buffer.cols_) {
        matrix_[i][j] = buffer.matrix_[i][j];
      } else {
        matrix_[i][j] = 0;
      }
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
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

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
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

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
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
  *this = result;
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

void calc_minor(double **matrix, int rows, int cols, int i_remove, int j_remove,
                S21Matrix *minor) {
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
  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  S21Matrix result(rows_, cols_);
  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      calc_minor(matrix_, rows_, cols_, i, j, &minor);
      result.matrix_[i][j] = minor.Determinant() * pow(-1, i + j);
    }
  }

  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  double det = 0;

  if (rows_ == 1) {
    det = matrix_[0][0];
  } else if (rows_ == 2) {
    det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    S21Matrix minor(rows_ - 1, cols_ - 1);
    int k = 1;
    for (int i = 0; i < rows_; i++) {
      calc_minor(matrix_, rows_, cols_, i, 0, &minor);
      det += k * matrix_[i][0] * minor.Determinant();
      k *= -1;
    }
  }

  return det;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (fabs(det) <= 1e-7) {
    throw std::out_of_range("Incorrect input, determinat = 0");
  }

  S21Matrix result = this->CalcComplements().Transpose();
  result.MulNumber(1 / det);

  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix &o) {
  S21Matrix result(*this);
  result.SumMatrix(o);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &o) {
  S21Matrix result(*this);
  result.SubMatrix(o);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &o) {
  S21Matrix result(*this);
  result.MulMatrix(o);
  return result;
}

S21Matrix operator*(const S21Matrix &lhs, const double num) {
  S21Matrix result(lhs);
  result.MulNumber(num);
  return result;
}

S21Matrix operator*(const double num, const S21Matrix &rhs) {
  S21Matrix result(rhs);
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix &o) { return this->EqMatrix(o); }

S21Matrix S21Matrix::operator=(const S21Matrix &rhs) {
  if (this == &rhs) {
    return *this;
  }
  rows_ = rhs.rows_;
  cols_ = rhs.cols_;

  this->~S21Matrix();

  matrix_ = new double *[rows_];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    memcpy(matrix_[i], rhs.matrix_[i], cols_ * sizeof(double));
  }

  return *this;
}

S21Matrix S21Matrix::operator+=(const S21Matrix &rhs) {
  this->SumMatrix(rhs);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &rhs) {
  this->SubMatrix(rhs);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &rhs) {
  this->MulMatrix(rhs);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  return matrix_[i][j];
}
