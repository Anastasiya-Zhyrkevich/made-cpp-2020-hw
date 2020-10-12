#include <math.h> 
#include "matrix.h"

using namespace task;

Matrix::Matrix() {
  rows = 1;
  cols = 1;
  data = new double* [rows];

  data[0] = new double[1];
  data[0][0] = 1.;
}

Matrix::Matrix(size_t rows, size_t cols){
  this->rows = rows;
  this->cols = cols;

  this->data = new double* [rows];

  for (size_t i = 0; i < rows; i++) {
    data[i] = new double[cols];
  }

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (i == j) {
        data[i][j] = 1.0;
      }
      else {
        data[i][j] = 0.0;
      }
    }
  }
}

Matrix::Matrix(const Matrix& copy){
  this->rows = copy.rows;
  this->cols = copy.cols;

  this->data = new double* [rows];

  for (size_t i = 0; i < rows; i++) {
    data[i] = new double[cols];
    for (size_t j = 0; j < cols; j++) {
      data[i][j] = copy.data[i][j];
    }
  }

}

Matrix& Matrix::operator=(const Matrix& copy){
  if (this == &copy) {
    return *this;
  }

  clear_data();

  this->rows = copy.rows;
  this->cols = copy.cols;

  this->data = new double* [rows];

  for (size_t i = 0; i < rows; i++) {
    data[i] = new double[cols];
    for (size_t j = 0; j < cols; j++) {
      data[i][j] = copy.data[i][j];
    }
  }
  return *this;
}

double* Matrix::operator[](size_t row){
  if (row >= rows) {
    throw OutOfBoundsException();
  }
  return data[row];
}

const double * Matrix::operator[](size_t row) const{
  if (row >= rows) {
    throw OutOfBoundsException();
  }
  return data[row];
}

Matrix::~Matrix() {
  clear_data();
}

double& Matrix::get(size_t row, size_t col) {
  if (row >= rows || col >= cols) {
    throw OutOfBoundsException();
  }
  return data[row][col];
}

const double& Matrix::get(size_t row, size_t col) const{
  if (row >= rows || col >= cols) {
    throw OutOfBoundsException();
  }
  return data[row][col];
}

void Matrix::set(size_t row, size_t col, const double& value) {
  if (row >= rows || col >= cols) {
    throw OutOfBoundsException();
  }
  data[row][col] = value;
}

void Matrix::clear_data() {
  for (size_t i = 0; i < this->rows; i++) {
    delete[] data[i];
  }
  delete[] data;
}

void Matrix::resize(size_t new_rows, size_t new_cols) {
  double ** new_data = new double*[new_rows];
  for (size_t i = 0; i < new_rows; i++) {
    new_data[i] = new double[new_cols];
    for (size_t j = 0; j < new_cols; j++) {
      new_data[i][j] = 0.0;
    }
  }

  for (size_t i = 0; i < std::min(new_rows, rows); i++) {
    for(size_t j = 0; j < std::min(new_cols, cols); j++) {
      new_data[i][j] = data[i][j];
    }
  }
  clear_data();

  rows = new_rows;
  cols = new_cols;
  data = new_data;
}

std::vector<double> Matrix::getRow(size_t row) {
  if (row >= rows) {
    throw OutOfBoundsException();
  }
  std::vector<double> rowV(cols);
  for (size_t j = 0; j < cols; j++) {
    rowV[j] = data[row][j];
  }
  return rowV;
}

std::vector<double> Matrix::getColumn(size_t column) {
  if (column >= cols) {
    throw OutOfBoundsException();
  }
  std::vector<double> colV(rows);
  for (size_t j = 0; j < rows; j++) {
    colV[j] = data[j][column];
  }
  return colV;
}

Matrix& Matrix::operator+=(const Matrix& a){
  if (rows != a.rows || cols != a.cols) {
    throw SizeMismatchException();
  }
  for (size_t i =0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      data[i][j] += a.data[i][j];
    }
  }
  return *(this);
}
Matrix& Matrix::operator-=(const Matrix& a){
  if (rows != a.rows || cols != a.cols) {
    throw SizeMismatchException();
  }
  for (size_t i =0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      data[i][j] -= a.data[i][j];
    }
  }
  return *(this);
}
Matrix& Matrix::operator*=(const Matrix& a){
  if (cols != a.rows) {
    throw SizeMismatchException();
  }

  size_t r = rows;
  size_t c = cols;
  Matrix left(*(this));
  this->resize(rows, a.cols);

  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < a.cols; j++) {
      double value = 0.0;
      for (size_t k = 0; k < c; k++) {
        value += left.data[i][k] * a.data[k][j];
      }
      data[i][j] = value;   
    }
  }
  return *(this);
}

Matrix& Matrix::operator*=(const double& number){
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      data[i][j] *= number;
    }
  }
  return *(this);
}

Matrix Matrix::operator+(const Matrix& a) const{
  Matrix result(*this);
  result += a;
  return result;
}
Matrix Matrix::operator-(const Matrix& a) const{
  Matrix result(*this);
  result -= a;
  return result;
}
Matrix Matrix::operator*(const Matrix& a) const{
  Matrix result(*this);
  result *= a;
  return result;
}
Matrix Matrix::operator*(const double& a) const{
  Matrix result(*this);
  result *= a;
  return result;
}

Matrix Matrix::operator-() const {
  Matrix result(*this);
  result *= (-1);
  return result;
}
Matrix Matrix::operator+() const {
  Matrix result(*this);
  return result;
}

double Matrix::det() const{
  if (rows != cols) {
    throw SizeMismatchException();
  }

  if (rows == 1) {
    return data[0][0];
  }

  if (rows == 2) {
    return data[0][0] * data[1][1] - data[0][1] * data[1][0];
  }

  Matrix tmp(rows - 1, cols - 1);

  double answer = 0.0;
  double coeff = 1.0;
  for (size_t j = 0; j < cols; j++) {
    // fill tmp 
    for (size_t i = 1; i < rows; i++) {
      for (size_t jj = 0; jj < cols; jj++) {
        if (jj == j) {
          continue;
        }
        size_t j_coeff = jj;
        if (j_coeff > j) {
          j_coeff -= 1;
        }
        tmp.data[i - 1][j_coeff] = data[i][jj];
      }
    }

    answer += coeff * data[0][j] * tmp.det();
    coeff *= (-1);
  }

  return answer;
}
void Matrix::transpose(){
  double ** new_data = new double*[cols];
  for (size_t i = 0; i < cols; i++) {
    new_data[i] = new double[rows];
    for (size_t j = 0; j < rows; j++) {
      new_data[i][j] = data[j][i];
    }
  }

  clear_data();

  size_t new_rows = cols;
  size_t new_cols = rows;
  rows = new_rows;
  cols = new_cols;
  data = new_data;
}

Matrix Matrix::transposed() const{
  Matrix result(*this);
  result.transpose();
  return result;
}

double Matrix::trace() const{
  if (rows != cols) {
    throw SizeMismatchException();
  }
  double result = 0.0;
  for (size_t i = 0; i < std::min(rows, cols); i++) {
    result += data[i][i];
  }
}

bool Matrix::operator==(const Matrix& a) const {
  if (rows != a.rows || cols != a.cols) {
    return false;
  }

  for(size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (fabs(data[i][j] - a.data[i][j]) > EPS) {
        return false;
      }
    }
  }
  return true;
}

bool Matrix::operator!=(const Matrix& a) const {
  return (!(*this == a));
}

Matrix operator*(const double& a, const Matrix& b) {
  Matrix result(b);
  result *= a;
  return result;
}

size_t Matrix::getRows() const{
  return rows;
}
size_t Matrix::getCols() const {
  return cols;
}

std::ostream& operator<<(std::ostream& output, const Matrix& matrix){
  for(size_t i = 0; i < matrix.getRows(); i++) {
    for (size_t j = 0; j < matrix.getCols(); j++) {
      output << matrix.get(i, j) << ' ';
    }
  }
  output << std::endl;

  return output;
}

std::istream& operator>>(std::istream& input, Matrix& matrix) {
  int r,c; 
  input >> r >> c;

  matrix.resize(r, c);
  for(size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < c; j++) {
      double value;
      input >> value;
      matrix.set(i, j, value);
    }
  }
  return input;
}


// Your code goes here...
