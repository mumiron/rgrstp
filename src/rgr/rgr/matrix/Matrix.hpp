#pragma once

#include <iostream>
#include <vector>

template <class T> class Matrix {

private:
  unsigned column_{0};
  unsigned row_{0};
  unsigned matrixSize_{0};
  T *matrix_;

  struct MatrixIterator {
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type *;
    using reference = value_type &;

    explicit MatrixIterator(pointer ptr = nullptr) : ptr_(ptr) {}

    reference operator*() const { return *ptr_; }

    pointer operator->() { return ptr_; }

    MatrixIterator &operator++() {
      ptr_++;
      return *this;
    }

    MatrixIterator &operator--() {
      ptr_--;
      return *this;
    }

    MatrixIterator operator++(int) { // NOLINT
      MatrixIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    MatrixIterator operator--(int) { // NOLINT
      MatrixIterator tmp = *this;
      --(*this);
      return tmp;
    }

    MatrixIterator operator+(size_t n) const {
      return MatrixIterator(ptr_ + n);
    }

    MatrixIterator &operator+=(size_t n) {
      ptr_ += n;
      return *this;
    }

    MatrixIterator operator-(size_t n) const {
      return MatrixIterator(ptr_ - n);
    }

    MatrixIterator &operator-=(size_t n) const {
      ptr_ -= n;
      return *this;
    }

    difference_type operator-(const MatrixIterator &iterator) {
      return ptr_ - iterator.ptr_;
    }

    bool operator==(const MatrixIterator &other) const {
      return ptr_ == other.ptr_;
    }

    bool operator!=(const MatrixIterator &other) const {
      return ptr_ != other.ptr_;
    }

    bool operator>(const MatrixIterator &other) const {
      return ptr_ > other.ptr_;
    }

    bool operator<(const MatrixIterator &other) const {
      return ptr_ < other.ptr_;
    }

    bool operator>=(const MatrixIterator &other) const {
      return ptr_ >= other.ptr_;
    }

    bool operator<=(const MatrixIterator &other) const {
      return ptr_ <= other.ptr_;
    }

    reference operator[](size_t n) const { return *(ptr_ + n); }

  private:
    pointer ptr_;
  };

public:
  MatrixIterator begin() const { return MatrixIterator(matrix_); }

  MatrixIterator end() const { return MatrixIterator(matrix_ + matrixSize_); }

  unsigned getColumn() const { return column_; }

  unsigned getRow() const { return row_; }

  unsigned int getMatrixSize() const { return matrixSize_; }

  T *getMatrix() const { return matrix_; }

  Matrix();

  Matrix(unsigned column, unsigned row);

  Matrix(unsigned column, unsigned row, const std::vector<T> &inputValues);

  Matrix(const Matrix<T> &copyMatrix);

  Matrix<T> &operator=(const Matrix<T> &copyMatrix);

  Matrix(Matrix<T> &&moveMatrix) noexcept;

  Matrix<T> &operator=(Matrix<T> &&moveMatrix) noexcept;

  ~Matrix() { delete[] matrix_; }

  const T &operator()(unsigned rows, unsigned cols) const {
    if (rows >= row_ || cols >= column_) {
      throw std::out_of_range("&operator() : Index out of bounds");
    }
    return matrix_[column_ * rows + cols];
  }

  T &operator()(unsigned rows, unsigned cols) {
    if (rows >= row_ || cols >= column_) {
      throw std::out_of_range("&operator() : Index out of bounds");
    }
    return matrix_[column_ * rows + cols];
  }

  Matrix<T> &operator+=(const Matrix<T> &bMatrix);

  Matrix<T> &operator-=(const Matrix<T> &bMatrix);

  Matrix<T> &operator*=(const T &scalar);

  Matrix<T> transpose();
};

template <class T> Matrix<T>::Matrix() : matrix_(nullptr) {}

template <class T>
Matrix<T>::Matrix(unsigned column, unsigned row)
    : column_(column), row_(row) { // пустая
  if (row_ == 0 || column_ == 0) {
    throw std::invalid_argument("Matrix constructor has 0 size");
  }
  matrixSize_ = column_ * row;
  matrix_ = new T[matrixSize_]();
}

template <class T>
Matrix<T>::Matrix(unsigned column, unsigned row, // с польз числ
                  const std::vector<T> &inputValues)
    : column_(column), row_(row) {
  if (row_ == 0 || column_ == 0) {
    throw std::invalid_argument("Matrix constructor has 0 size");
  }
  matrixSize_ = column_ * row_;
  matrix_ = new T[matrixSize_];

  std::copy(inputValues.begin(), inputValues.end(), matrix_);
}

template <class T> // копирование
Matrix<T>::Matrix(const Matrix<T> &copyMatrix)
    : column_(copyMatrix.column_), row_(copyMatrix.row_),
      matrixSize_(copyMatrix.matrixSize_), matrix_(new T[matrixSize_]) {

  std::copy(copyMatrix.begin(), copyMatrix.end(), matrix_);
}

template <class T> // копир через присваивание
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &copyMatrix) {
  if (&copyMatrix == *this) {
    return *this;
  }
  Matrix tmp(copyMatrix);
  std::swap(tmp);
  return *this;
}

template <class T> // перемещение (-1 мат)
Matrix<T>::Matrix(Matrix<T> &&moveMatrix) noexcept
    : column_(std::move(moveMatrix.column_)), row_(std::move(moveMatrix.row_)),
      matrixSize_(std::move(moveMatrix.matrixSize_)),
      matrix_(std::move(moveMatrix.matrix_)) {
  moveMatrix.column_ = 0;
  moveMatrix.row_ = 0;
  moveMatrix.matrix_ = nullptr;
}

template <class T> // перемещ через присв
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&moveMatrix) noexcept {
  Matrix tmp(std::move(moveMatrix));
  std::swap(tmp);
  return *this;
}

template <class T> // перегруж оператор вывода
std::ostream &operator<<(std::ostream &os, const Matrix<T> &outputMatrix) {
  for (size_t i = 0; i < outputMatrix.getRow(); ++i) {
    for (size_t j = 0; j < outputMatrix.getColumn(); ++j) {
      os << outputMatrix.getMatrix()[outputMatrix.getColumn() * i + j] << " ";
    }
    os << "\n";
  }
  return os;
}

template <class T> // оператор сравнения (равенство)
bool operator==(const Matrix<T> &aMatrix, const Matrix<T> &bMatrix) {
  if (aMatrix.getRow() != bMatrix.getRow() &&
      aMatrix.getColumn() != bMatrix.getColumn()) {
    throw std::invalid_argument("operator== : Different dimensions");
  }
  for (size_t i = 0; i < aMatrix.getMatrixSize(); ++i) {
    if (aMatrix.getMatrix()[i] != bMatrix.getMatrix()[i]) {
      return false;
    }
  }
  return true;
}

template <class T> // неравнество
bool operator!=(const Matrix<T> &aMatrix, const Matrix<T> &bMatrix) {
  if (aMatrix.getRow() != bMatrix.getRow() &&
      aMatrix.getColumn() != bMatrix.getColumn()) {
    throw std::invalid_argument("bool operator!= : Different dimensions");
  }
  for (size_t i = 0; i < aMatrix.getMatrixSize(); ++i) {
    if (aMatrix.getMatrix()[i] != bMatrix.getMatrix()[i]) {
      return true;
    }
  }
  return false;
}

template <class T> // сложение двух матриц в новую
Matrix<T> operator+(const Matrix<T> &aMatrix, const Matrix<T> bMatrix) {
  if (aMatrix.getRow() != bMatrix.getRow() &&
      aMatrix.getColumn() != bMatrix.getColumn()) {
    throw std::invalid_argument("operator+ : Different dimensions");
  }
  Matrix<T> matrixProduct(aMatrix.getRow(), aMatrix.getColumn()); // резмат
  for (size_t i = 0; i < aMatrix.getMatrixSize(); ++i) {
    matrixProduct.getMatrix()[i] =
        aMatrix.getMatrix()[i] + bMatrix.getMatrix()[i];
  }
  return matrixProduct;
}

template <class T> // сложение двух в первую матрицу
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &bMatrix) {
  if (row_ != bMatrix.row_ && column_ != bMatrix.column_) {
    throw std::invalid_argument("operator+= : Different dimensions");
  }
  for (size_t i = 0; i < matrixSize_; ++i) {
    matrix_[i] += bMatrix.matrix_[i];
  }
  return *this;
}

template <class T> // вычитание двух в новую мат
Matrix<T> operator-(const Matrix<T> &aMatrix, const Matrix<T> bMatrix) {
  if (aMatrix.getRow() != bMatrix.getRow() &&
      aMatrix.getColumn() != bMatrix.getColumn()) {
    throw std::invalid_argument("operator- : Different dimensions");
  }
  Matrix<T> matrixProduct(aMatrix.getRow(), aMatrix.getColumn());
  for (size_t i = 0; i < aMatrix.getMatrixSize(); ++i) {
    matrixProduct.getMatrix()[i] =
        aMatrix.getMatrix()[i] - bMatrix.getMatrix()[i];
  }
  return matrixProduct;
}

template <class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &bMatrix) { // вычитание в одну
  if (row_ != bMatrix.row_ && column_ != bMatrix.column_) {
    throw std::invalid_argument("operator+= : Different dimensions");
  }
  for (size_t i = 0; i < matrixSize_; ++i) {
    matrix_[i] -= bMatrix.matrix_[i];
  }
  return *this;
}

template <class T> // умножение двух матриц
Matrix<T> operator*(const Matrix<T> &aMatrix, const Matrix<T> bMatrix) {
  if (aMatrix.getColumn() != bMatrix.getRow()) {
    throw std::invalid_argument("operator* : The number of A-matrix column has "
                                "to be equal to the number of"
                                "B-matrix row");
  }
  Matrix<T> matrixProduct(bMatrix.getColumn(), aMatrix.getRow());

  for (size_t i = 0; i < aMatrix.getRow(); ++i) {
    for (size_t j = 0; j < bMatrix.getColumn(); ++j) {
      for (size_t k = 0; k < aMatrix.getRow(); ++k) {
        matrixProduct(i, j) += aMatrix(i, k) * bMatrix(k, j);
      }
    }
  }
  return matrixProduct;
}

template <class T>
Matrix<T> operator*(const Matrix<T> &aMatrix,
                    T scalar) { // умнож на число в новую
  unsigned size = aMatrix.getRow() * aMatrix.getColumn();
  Matrix<T> matrixProduct(aMatrix.getColumn(), aMatrix.getRow());
  for (size_t i = 0; i < size; ++i) {
    matrixProduct.getMatrix()[i] = aMatrix.getMatrix()[i] * scalar;
  }
  return matrixProduct;
}

template <class T> // умнож на число в стар
Matrix<T> &Matrix<T>::operator*=(const T &scalar) {
  unsigned size = getRow() * getColumn();
  for (size_t i = 0; i < size; ++i) {
    getMatrix()[i] *= scalar;
  }
  return *this;
}

template <class T> // транспонирование
Matrix<T> Matrix<T>::transpose() {
  Matrix<T> matrixProduct(column_, row_);
  for (size_t i = 0; i < row_; ++i) {
    for (size_t j = 0; j < column_; ++j) {
      matrixProduct(i, j) = matrix_[column_ * j + i];
    }
  }
  return matrixProduct;
}
