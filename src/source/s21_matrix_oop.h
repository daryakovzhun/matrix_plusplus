#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <string.h>

#include <exception>
#include <stdexcept>

#define SUCCESS 1
#define FAILURE 0

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);

  ~S21Matrix();

  int getRows();
  int getCols();
  void setRows(const int rows);
  void setCols(const int cols);

  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix operator+(const S21Matrix &o);
  S21Matrix operator-(const S21Matrix &o);
  S21Matrix operator*(const S21Matrix &o);
  bool operator==(const S21Matrix &o);
  S21Matrix operator=(const S21Matrix &rhs);
  S21Matrix operator+=(const S21Matrix &rhs);
  S21Matrix operator-=(const S21Matrix &rhs);
  S21Matrix operator*=(const double num);
  S21Matrix operator*=(const S21Matrix &rhs);
  double &operator()(int x, int y);

  friend void calc_minor(double **matrix, int rows, int cols, int i_remove,
                         int j_remove, S21Matrix *minor);

  friend S21Matrix operator*(const S21Matrix &lhs, const double num);
  friend S21Matrix operator*(const double num, const S21Matrix &rhs);

 private:
  int rows_, cols_;
  double **matrix_;
};

#endif  //  SRC_S21_MATRIX_H_
