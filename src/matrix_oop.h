#ifndef SRC_MATRIX_OOP_H
#define SRC_MATRIX_OOP_H

#include <cmath>
#include <iostream>

class Matrix {
 public:
  // constructors
  Matrix();
  Matrix(int rows, int cols);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other) noexcept;
  Matrix& operator=(const Matrix& other);
  Matrix& operator=(Matrix&& other) noexcept;

  // destructor
  ~Matrix();

  // matrix operations
  bool EqMatrix(const Matrix& other) const;
  void SumMatrix(const Matrix& other);
  void SubMatrix(const Matrix& other);
  void MulNumber(const double number);
  void MulMatrix(const Matrix& other);
  Matrix Transpose() const;
  Matrix CalcComplements() const;
  Matrix CreateMinor(const int i, const int j) const;
  double Determinant() const;
  Matrix InverseMatrix() const;

  // overloads
  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix operator*(const Matrix& other) const;
  Matrix operator*(double number) const;
  friend Matrix operator*(const double number, const Matrix& other);
  bool operator==(Matrix const& other) const;
  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  Matrix& operator*=(const Matrix& other);
  Matrix& operator*=(const double number) noexcept;
  const double& operator()(int row_index, int col_index) const;
  double& operator()(int row_index, int col_index);

  // mutators & accessors
  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows_number);
  void SetCols(int cols_number);

  // other functions
  friend std::ostream& operator<<(std::ostream& out, const Matrix& p);
  void Delete();

 private:
  // data members
  int rows_, cols_;
  double** matrix_;
};

#endif  // SRC_MATRIX_OOP_H