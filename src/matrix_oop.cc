#include "matrix_oop.h"

Matrix::Matrix() : rows_(1), cols_(1), matrix_(nullptr) {
  matrix_ = new double*[rows_]();
  matrix_[0] = new double[cols_]();
}

Matrix::Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  if (rows <= 0 || cols <= 0) {
    throw std::exception();
  }
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_]();
  }
}

Matrix::Matrix(const Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

Matrix::Matrix(Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
}

Matrix& Matrix::operator=(const Matrix& other) {
  if (this != &other) {
    Delete();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {
  if (this != &other) {
    Delete();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
  }
  return *this;
}

Matrix::~Matrix() {
  if (matrix_) {
    Delete();
  }
}

bool Matrix::EqMatrix(const Matrix& other) const {
  bool result = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    result = false;
  } else {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) >= 1e-07) {
          result = false;
        }
      }
    }
  }
  return result;
}

void Matrix::SumMatrix(const Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::exception();
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void Matrix::SubMatrix(const Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::exception();
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void Matrix::MulNumber(const double number) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= number;
    }
  }
}

void Matrix::MulMatrix(const Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::exception();
  }
  Matrix matrix_tmp(rows_, other.cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      double result = 0;
      for (int k = 0; k < cols_; ++k) {
        result += matrix_[i][k] * other.matrix_[k][j];
      }
      matrix_tmp.matrix_[i][j] = result;
    }
  }
  *this = std::move(matrix_tmp);
}

Matrix Matrix::Transpose() const {
  Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

Matrix Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::exception();
  }
  Matrix res(rows_, cols_);
  int correction_r = (rows_ > 1) ? 1 : 0;
  int correction_c = (cols_ > 1) ? 1 : 0;
  Matrix minor(rows_ - correction_r, cols_ - correction_c);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (rows_ > 1)
        minor = this->CreateMinor(i, j);
      else
        minor.matrix_[i][j] = 1;
      double det = minor.Determinant();
      res.matrix_[i][j] = pow(-1, (i + j)) * det;
    }
  }
  return res;
}

Matrix Matrix::CreateMinor(const int i, const int j) const {
  Matrix res(cols_ - 1, rows_ - 1);
  for (int r = 0; r < rows_; r++) {
    for (int c = 0; c < cols_; c++) {
      if (r != i && c != j) {
        int correction_r = (r > i) ? -1 : 0;
        int correction_c = (c > j) ? -1 : 0;
        res.matrix_[r + correction_r][c + correction_c] = matrix_[r][c];
      }
    }
  }
  return res;
}

double Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::exception();
  }
  Matrix tmp(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  double result = 1;
  for (int i = 0; i < tmp.rows_; i++) {
    if (tmp.matrix_[i][i] == 0) {
      for (int j = i + 1; j < tmp.rows_; j++) {
        if (tmp.matrix_[j][i] != 0) {
          std::swap(tmp.matrix_[i], tmp.matrix_[j]);
          result *= -1;
        }
      }
    }
    if (tmp.matrix_[i][i] == 0) {
      result = 0;
      break;
    }
    result *= tmp.matrix_[i][i];
    double divider = tmp.matrix_[i][i];
    for (int j = i; j < tmp.cols_; j++) {
      tmp.matrix_[i][j] = tmp.matrix_[i][j] / divider;
    }
    for (int j = i + 1; j < tmp.rows_; j++) {
      double multiplier = tmp.matrix_[j][i] / tmp.matrix_[i][i];
      for (int k = i; k < tmp.cols_; k++) {
        tmp.matrix_[j][k] -= tmp.matrix_[i][k] * multiplier;
      }
    }
  }
  return result;
}

Matrix Matrix::InverseMatrix() const {
  double determinant = this->Determinant();
  if (determinant == 0) {
    throw std::exception();
  }
  Matrix result(rows_, cols_);
  Matrix matrix_complements = std::move(this->CalcComplements());
  Matrix matrix_complements_transpose =
      std::move(matrix_complements.Transpose());
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[i][j] =
          (1 / determinant) * matrix_complements_transpose.matrix_[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
  Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
  Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
  Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

Matrix Matrix::operator*(double number) const {
  Matrix result(*this);
  result.MulNumber(number);
  return result;
}

Matrix operator*(const double number, const Matrix& other) {
  Matrix result(other);
  result.MulNumber(number);
  return result;
}

bool Matrix::operator==(Matrix const& other) const {
  return EqMatrix(other);
}

Matrix& Matrix::operator+=(const Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

Matrix& Matrix::operator*=(const Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

Matrix& Matrix::operator*=(const double number) noexcept {
  this->MulNumber(number);
  return *this;
}

const double& Matrix::operator()(int row_index, int col_index) const {
  if (row_index < 0 || row_index >= rows_ || col_index < 0 ||
      col_index >= cols_) {
    throw std::exception();
  }
  return matrix_[row_index][col_index];
}

double& Matrix::operator()(int row_index, int col_index) {
  if (row_index < 0 || row_index >= rows_ || col_index < 0 ||
      col_index >= cols_) {
    throw std::exception();
  }
  return matrix_[row_index][col_index];
}

int Matrix::GetRows() const { return rows_; }

int Matrix::GetCols() const { return cols_; }

void Matrix::SetRows(int rows_number) {
  if (rows_number < 1) {
    throw std::exception();
  }
  Matrix result(rows_number, cols_);
  for (int i = 0; i < rows_number; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (rows_ > i && cols_ > j) result.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(result);
}

void Matrix::SetCols(int cols_number) {
  if (cols_number < 1) {
    throw std::exception();
  }
  Matrix result(rows_, cols_number);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_number; ++j) {
      if (rows_ > i && cols_ > j) result.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(result);
}

std::ostream& operator<<(std::ostream& out, const Matrix& p) {
  for (int i = 0; i < p.rows_; ++i) {
    for (int j = 0; j < p.cols_; ++j) {
      out << p.matrix_[i][j] << "\t";
    }
    out << '\n';
  }
  return out;
}

void Matrix::Delete() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}