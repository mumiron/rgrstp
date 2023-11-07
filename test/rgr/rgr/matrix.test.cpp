#include <gtest/gtest.h>

#include <rgr/matrix/Matrix.hpp>

TEST(Matrix, Copy) {
  std::vector<int> aVector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix<int> aMatrix(3, 3, aVector);
  const auto &bMatrix(aMatrix);
  EXPECT_EQ(aMatrix, bMatrix);
}

TEST(Matrix, CopyAssigment) {
  std::vector<int> aVector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix<int> aMatrix(3, 3, aVector);
  const auto &bMatrix = aMatrix;
  EXPECT_EQ(aMatrix, bMatrix);
}

TEST(Matrix, Move) {
  Matrix<int> aMatrix(2, 2, {1, 2, 3, 4});
  const unsigned aMatrixCols = aMatrix.getColumn();
  const unsigned aMatrixRows = aMatrix.getRow();
  auto *aMatrixData = aMatrix.getMatrix();

  auto bMatrix(std::move(aMatrix));

  EXPECT_EQ(bMatrix.getColumn(), aMatrixCols);
  EXPECT_EQ(bMatrix.getRow(), aMatrixRows);
  EXPECT_EQ(bMatrix.getMatrix(), aMatrixData);
}

TEST(Matrix, MoveAssignment) {
  std::vector<int> bVector = {1, 2, 3, 4};
  Matrix<int> bMatrix(4, 1, bVector);
  Matrix<int> aMatrix = Matrix<int>(4, 1, bVector);

  EXPECT_EQ(aMatrix, bMatrix);
}

TEST(Matrix, Equality) {                                    //сравнение на равенство
  std::vector<int> aVector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix<int> aMatrix(3, 3, aVector);
  Matrix<int> bMatrix(3, 3, aVector);

  EXPECT_EQ(aMatrix, bMatrix);

  std::vector<int> bVector = {1, 1, 1, 1, 1, 1, 1, 1, 1};

  const Matrix<int> &cMatrix = aMatrix;
  Matrix<int> dMatrix(3, 3, bVector);

  EXPECT_NE(cMatrix, dMatrix);
}

TEST(Matrix, Subscript) {                                 //доступ к элементу
  Matrix<int> aMatrix(2, 2);
  aMatrix(1, 1) = 5;
  EXPECT_EQ(aMatrix(1, 1), 5);
}

TEST(Matrix, Addition) {                           //сложение
  std::vector<int> aVector = {1, 2, 3};
  std::vector<int> bVector = {4, 5, 6};
  std::vector<int> expectedVector = {5, 7, 9, 0};

  Matrix<int> aMatrix(2, 2, aVector);
  Matrix<int> bMatrix(2, 2, bVector);
  Matrix<int> expectedMatrix(2, 2, expectedVector);

  EXPECT_EQ(aMatrix + bMatrix, expectedMatrix);
  EXPECT_EQ(aMatrix += bMatrix, expectedMatrix);
}

TEST(Matrix, Subtraction) {             //вычит
  std::vector<int> aVector = {1, 2, 3};
  std::vector<int> bVector = {4, 5, 6};

  Matrix<int> aMatrix(2, 2, aVector);
  Matrix<int> bMatrix(2, 2, bVector);
  Matrix<int> expectedMatrix(2, 2, {3, 3, 3, 0});

  EXPECT_EQ(bMatrix - aMatrix, expectedMatrix);
  EXPECT_EQ(bMatrix -= aMatrix, expectedMatrix);
}

TEST(Matrix, Multiplication) {       //умнож
  std::vector<int> aVector = {1, 2, 3, 4};
  std::vector<int> bVector = {5, 6, 7, 8};
  std::vector<int> expectedVector = {19, 22, 43, 50};

  Matrix<int> aMatrix(2, 2, aVector);
  Matrix<int> bMatrix(2, 2, bVector);
  Matrix<int> expectedMatrix(2, 2, expectedVector);

  EXPECT_EQ(aMatrix * bMatrix, expectedMatrix);

  Matrix<int> expMatrix(2, 2, {3, 6, 9, 12});

  EXPECT_EQ(aMatrix * 3, expMatrix);
  EXPECT_EQ(aMatrix *= 3, expMatrix);
}

TEST(Matrix, Transpose) {
  std::vector<int> aVector = {1, 2, 3, 4};
  std::vector<int> transposeVector = {1, 3, 2, 4};
  Matrix<int> aMatrix(2, 2, aVector);
  Matrix<int> transposeMatrix(2, 2, transposeVector);

  EXPECT_EQ(aMatrix.transpose(), transposeMatrix);
}

TEST(Matrix, Iterator) {              //стл алг
  std::vector<int> aVector = {1, 2, 3, 3, 4, 1, 2, 3, 4};
  std::vector<int> bVector = {1, 2, 3};

  Matrix<int> aMatrix(3, 3, aVector);
  Matrix<int> bMatrix(3, 1, bVector);

  auto find_end_result = std::find_end(aMatrix.begin(), aMatrix.end(),
                                       bMatrix.begin(), bMatrix.end());
  auto result = std::distance(aMatrix.begin(), find_end_result);
  EXPECT_EQ(result, 5);

  auto is_even = [](int i) { return i % 2 == 0; };
  auto find_if_result = std::find_if(bMatrix.begin(), bMatrix.end(), is_even);
  EXPECT_EQ(*find_if_result, 2);

  Matrix<int> filledMatrix(3, 1, {1, 1, 1});
  std::fill(bMatrix.begin(), bMatrix.end(), 1);
  EXPECT_EQ(bMatrix, filledMatrix);

  Matrix<int> sortedMatrix(3, 3, {1, 1, 2, 2, 3, 3, 3, 4, 4});
  std::sort(aMatrix.begin(), aMatrix.end());
  EXPECT_EQ(aMatrix, sortedMatrix);
}

TEST(Matrix, Exception) {               //на исключение
  Matrix<int> aMatrix(2, 2);
  EXPECT_THROW(aMatrix(4, 6) = 5, std::out_of_range);

  EXPECT_THROW(Matrix<int> bMatrix(0, 0), std::invalid_argument);
  EXPECT_THROW(Matrix<int> cMatrix(0, 0, {1, 2, 3, 4}), std::invalid_argument);

  Matrix<int> dMatrix(4, 4);
  EXPECT_THROW(aMatrix + dMatrix, std::invalid_argument);
  EXPECT_THROW(aMatrix += dMatrix, std::invalid_argument);
  EXPECT_THROW(aMatrix - dMatrix, std::invalid_argument);
  EXPECT_THROW(aMatrix -= dMatrix, std::invalid_argument);
  EXPECT_THROW(aMatrix * dMatrix, std::invalid_argument);
}
