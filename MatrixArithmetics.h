/*
  ==============================================================================

    MatrixArithmetics.h
    Created: 25 Nov 2022 5:02:28pm
    Author:  Selim Sheta
    
    MatrixArithmetics is a header only "micro-library" containing useful basic functions 
    for manipulating vectors and matrices. Functions starting with "vector" can only be applied to
    a vector, and functions starting with "matrix" can only be applied to a matrix. 

    The library is designed for maximum flexibility, prioritising no compile errors over
    "mathematical safety". For example, you can do a dot product between two vectors of different sizes. 
    The result will obviously be wrong, but you won't get errors or warnings.

    Strongly recommended: "using namespace arrmath" at the top of your file.
    Notes: 
    - This library is not particularly fast, to be avoided for time-critical operations.
    - This library hasn't been fully tested as of 14/12/2022
    - To be tested : matrix functions, coord system functions
    - To be added : more functions, support for pass-by-reference for all functions

  ==============================================================================
*/

#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>

namespace arrmath {

    static const bool DISABLE_WARNINGS = true;

    template<typename T>
    using matrix = std::vector<std::vector<T>>;

    template<typename T>
    using vector = std::vector<T>;

    // index pair for matrices
    using mtxIndex = std::pair<size_t, size_t>;

    enum class CoordSystem { cart = 0, pol = 1, sph = 2 };

    // create vector of zeros
    template<typename T>
    vector<T> vectorOfZeros(size_t size);

    // get the maximum value in a vector
    template<typename T>
    T vectorMax(vector<T> vec);

    // get the minimum value in a vector
    template<typename T>
    T vectorMin(vector<T> vec);

    // get the minimum and maximum values in a vector with one operation. 
    // The result is a pair in the form {min, max}
    template<typename T>
    std::pair<T, T> vectorMinMax(vector<T> vec);

    // get the difference between the largest and smallest value in a vector
    template<typename T>
    T vectorSpan(vector<T> vec);

    // get Index of largest element in a vector
    template<typename T>
    size_t vectorMaxIndex(vector<T> vec);

    // get Index of smallest element in a vector
    template<typename T>
    size_t vectorMinIndex(vector<T> vec);

    // Apply zero-averaging to a vector
    // This scales the positive and/or negative elements in a vector
    // So that the average of all values is zero.
    template<typename T>
    vector<T> vectorZeroAvg(vector<T> vec);

    // Element-wise mutliplication of two vectors
    template<typename T>
    vector<T> vectorElemMult(vector<T> vec1, vector<T> vec2);

    // Element-wise addition of two vectors
    template<typename T>
    vector<T> vectorElemAdd(vector<T> vec1, vector<T> vec2);

    // Multiply vector by scalar
    template<typename T>
    vector<T> vectorScale(vector<T> vec, T scalar);

    // Transpose a vector
    // Note - This will result in a "column vector":
    // [1, 2, 3] (1 row, 3 columns) => [[1], [2], [3]] (3 rows, 1 column).
    template<typename T>
    matrix<T> vectorTranspose(vector<T> vec);

    // Compute the dot product of two vectors
    template<typename T>
    T vectorDotProduct(vector<T> vec1, vector<T> vec2);

    // Compute the cross product of two vectors
    template<typename T>
    vector<T> vectorCrossProduct(vector<T> vec1, vector<T> vec2);

    // Compute the magnitude (Euclidean norm) of a vector
    template<typename T>
    T vectorMagnitude(vector<T> vec);

    // Divide all the values in the vector by its magnitude
    template<typename T>
    vector<T> vectorNormalise(vector<T> vec);

    // Compute the mean of all values in a vector
    template<typename T>
    T vectorAverage(vector<T> vec);

    // Compute the average of two points
    template<typename T>
    vector<T> vectorAverage(vector<T> vec1, vector<T> vec2);

    // Compute the average of multiple points stored in a matrix.
    // All points must have the same number of coordinates. The points must be
    // rows in the matrix.
    template<typename T>
    vector<T> vectorAverage(matrix<T> mtx);

    // Compute the distance between two vectors. Use the coordinateSystem argument to specify
    // if the input vectors are cartesian (cart), polar (pol), or spherical (sph). 
    // By default it's assumed that both vectors are cartesian. 
    template<typename T>
    T vectorGetDistance(vector<T> vec1, vector<T> vec2, CoordSystem coordinateSystem = CoordSystem::cart);

    // Apply arbitrary function to all values in a vector
    template<typename T>
    vector<T> vectorFunc(vector<T> vec, std::function<T(T)> func);

    // Applies quick sort to the vector and returns the sorted version in ascending order. 
    // EG: [0, 100, 50, 25, 75] => [0, 25, 50, 75, 100]
    template<typename T>
    vector<size_t> vectorQuickSort(vector<T> vec);

    // Sorts vector & outputs the indices in order.
    // EG: [0, 100, 50, 25, 75] => [0, 3, 2, 4, 1]
    template<typename T>
    vector<size_t> vectorSortIndices(vector<T> vec);

    // create matrix of zeros
    template<typename T>
    matrix<T> matrixOfZeros(size_t rows, size_t cols);

    // get the maximum value in a matrix
    template<typename T>
    T matrixMax(matrix<T> mtx);

    // get the minimum value in a matrix
    template<typename T>
    T matrixMin(matrix<T> mtx);

    // get the minimum and maximum values in a matrix in one operation.
    template<typename T>
    std::pair<T,T> matrixMinMax(matrix<T> mtx);

    // Element-wise mutliplication of two matrices
    template<typename T>
    matrix<T> matrixElemMult(matrix<T> mtx1, matrix<T> mtx2);

    // Element-wise addition of two matrices
    template<typename T>
    matrix<T> matrixElemAdd(matrix<T> mtx1, matrix<T> mtx2);

    // Multiply matrix by scalar
    template<typename T>
    matrix<T> matrixScale(matrix<T> mtx, T scalar);

    // Multiply two matrices using regular matrix multiplication rules
    template<typename T>
    matrix<T> matrixMultiply(matrix<T> mtx1, matrix<T> mtx2);

    // Transpose a matrix
    template<typename T>
    matrix<T> matrixTranspose(matrix<T> mtx);

    // Remap all values to fit in the desired range, by default [0,1]
    template<typename T>
    matrix<T> matrixRemap(matrix<T> mtx, T min = 0.0, T max = 1.0);

    // Compute the average of all values in a matrix
    template<typename T>
    T matrixAverage(matrix<T> mtx);

    // apply arbitrary function to all values in a matrix
    template<typename T>
    matrix<T> matrixFunc(matrix<T> mtx, std::function<T(T)> func);

    // Substract the smallest value in a matrix from all its element
    template<typename T>
    matrix<T> matrixShiftToZero(matrix<T> mtx);

    // Apply quick sort to a matrix.
    template<typename T>
    matrix<T> matrixQuickSort(matrix<T> mtx);

    // Sorts matrix & outputs the indices in order.
    template<typename T>
    vector<mtxIndex> matrixSortIndices(matrix<T> mtx);

    // Convert polar coordinates to a cartesian vector
    template<typename T>
    vector<T> polToCart(T theta, T radius);

    // Convert polar vector to cartesian vector
    template<typename T>
    vector<T> polToCart(vector<T> vec);

    // Convert polar coordinates stored in separate vectors to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> polToCart(vector<T> theta, vector<T> radius, bool transpose = false);

    // Convert polar vectors stored in a matrix to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> polToCart(matrix<T> mtx, bool transpose = false);

    // convert spherical coordinates to a cartesian vector
    template<typename T>
    vector<T> sphToCart(T theta, T phi, T radius);

    // convert spherical vector to cartesian vector
    template<typename T>
    vector<T> sphToCart(vector<T> vec);

    // Convert spherical coordinates stored in separate vectors to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> sphToCart(vector<T> theta, vector<T> phi, vector<T> radius, bool transpose = false);

    // Convert spherical vectors stored in a matrix to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> sphToCart(matrix<T> mtx, bool transpose = false);

    //////////// QUICKSORT UTILITIES ///////////

    namespace utils {

        // value, index i, index j
        template<typename T>
        using valIdx = std::pair<T, std::pair<size_t, size_t>>;

        template<typename T>
        int partition(vector<valIdx<T>>& arr, int start, int end)
        {
            T pivot = arr[start].first;

            int count = 0;
            for (int i = start + 1; i <= end; i++) {
                if (arr[i].first <= pivot) count++;
            }

            // Giving pivot element its correct position
            int pivotIndex = start + count;
            std::swap(arr[pivotIndex], arr[start]);

            // Sorting left and right parts of the pivot element
            int i = start;
            int j = end;

            while (i < pivotIndex && j > pivotIndex) {

                while (arr[i].first <= pivot) {
                    i++;
                }

                while (arr[j].first > pivot) {
                    j--;
                }

                if (i < pivotIndex && j > pivotIndex) {
                    std::swap(arr[i++], arr[j--]);
                }
            }

            return pivotIndex;
        }

        template<typename T>
        void quickSort(vector<valIdx<T>>& arr, int start, int end)
        {
            // base case
            if (start >= end)
                return;

            // partitioning the array
            int p = partition(arr, start, end);

            // Sorting the left part
            quickSort(arr, start, p - 1);

            // Sorting the right part
            quickSort(arr, p + 1, end);
        }
    }
   

    //=================//
    //  VECTOR MATHS   //
    //=================//

    // create vector of zeros
    template<typename T>
    vector<T> vectorOfZeros(size_t size) {
        return vector<T>(size, static_cast<T>(0.0));
    }

    // get the maximum value in a vector
    template<typename T>
    T vectorMax(vector<T> vec) {
        T max = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            max = (vec[i] > max) ? vec[i] : max;
        }
        return max;
    }

    // get the minimum value in a vector
    template<typename T>
    T vectorMin(vector<T> vec) {
        T min = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            min = (vec[i] < min) ? vec[i] : min;
        }
        return min;
    }

    // get the minimum and maximum values in a vector with one operation. 
    // The result is a pair in the form {min, max}
    template<typename T>
    std::pair<T,T> vectorMinMax(vector<T> vec) {
        T min = vec[0];
        T max = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            min = (vec[i] < min) ? vec[i] : min;
            max = (vec[i] > max) ? vec[i] : max;
        }
        return std::pair<T,T>{ min, max };
    }

    // get the difference between the largest and smallest value in a vector
    template<typename T>
    T vectorSpan(vector<T> vec) {
        T max = vec[0];
        T min = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i] > max) {
                max = vec[i];
            }
            else if (vec[i] < min) {
                min = vec[i];
            }
        }
        return abs(max - min);
    }

    // get Index of largest element in a vector
    template<typename T>
    size_t vectorMaxIndex(vector<T> vec) {
        size_t index = 0;
        T max = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i] > max) {
                max = vec[i];
                index = i;
            }
        }
        return index;
    }

    // get Index of smallest element in a vector
    template<typename T>
    size_t vectorMinIndex(vector<T> vec) {
        size_t index = 0;
        T min = vec[0];
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i] < min) {
                min = vec[i];
                index = i;
            }
        }
        return index;
    }

    // Apply zero-averaging to a vector
    // This scales the positive and/or negative elements in a vector
    // So that the average of all values is zero.
    template<typename T>
    vector<T> vectorZeroAvg(vector<T> vec) {
        T vposSum = 0.0;
        T vnegSum = 0.0;
        for (size_t i = 0; i < vec.size(); i++) {
            if (vec[i] > 0) {
                vposSum += vec[i];
            }
            else {
                vnegSum -= vec[i];
            }
        }
        if (vposSum == 0 || vnegSum == 0) {
            return vectorOfZeros<T>(vec.size());
        }
        else if (vposSum > vnegSum) {
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] > 0) {
                    vec[i] *= (vnegSum / vposSum);
                }
            }
            return vec;
        }
        else if (vposSum < vnegSum) {
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] < 0) {
                    vec[i] *= (vposSum / vnegSum);
                }
            }
            return vec;
        }
        else {
            return vec;
        }
    }

    // Element-wise mutliplication of two vectors
    template<typename T>
    vector<T> vectorElemMult(vector<T> vec1, vector<T> vec2) {
        if (vec1.size() != vec2.size() && !DISABLE_WARNINGS) {
            std::cout << "elementMult | Warning: dimensions don't match." << std::endl;
        }
        size_t oplength = std::min(vec1.size(), vec2.size());
        for (size_t i = 0; i < oplength; i++) {
            vec1[i] *= vec2[i];
        }
        return vec1;
    }

    // Element-wise addition of two vectors
    template<typename T>
    vector<T> vectorElemAdd(vector<T> vec1, vector<T> vec2) {
        if (vec1.size() != vec2.size() && !DISABLE_WARNINGS) {
            std::cout << "elementAdd | Warning: dimensions don't match." << std::endl;
        }
        size_t oplength = std::min(vec1.size(), vec2.size());
        for (size_t i = 0; i < oplength; i++) {
            vec1[i] += vec2[i];
        }
        return vec1;
    }

    // Multiply vector by scalar
    template<typename T>
    vector<T> vectorScale(vector<T> vec, T scalar) {
        for (size_t i = 0; i < vec.size(); i++) {
            vec[i] *= scalar;
        }
        return vec;
    }

    // Transpose a vector
    // Note - This will result in a "column vector":
    // [1, 2, 3] (1 row, 3 columns) => [[1], [2], [3]] (3 rows, 1 column).
    template<typename T>
    matrix<T> vectorTranspose(vector<T> vec) {
        matrix<T> result = matrixOfZeros<T>(vec.size(), 1);
        for (size_t i = 0; i < vec.size(); i++) {
            result[i][0] = vec[i];
        }
        return result;
    }

    // Compute the dot product of two vectors
    template<typename T>
    T vectorDotProduct(vector<T> vec1, vector<T> vec2) {
        T result = static_cast<T>(0.0);
        if (vec1.size() != vec2.size() && !DISABLE_WARNINGS) {
            std::cout << "vectorDotProduct | Warning : dimensions don't match." << std::endl;
        }
        size_t oplength = std::min(vec1.size(), vec2.size());
        for (size_t i = 0; i < oplength; i++) {
            result += vec1[i] * vec2[i];
        }
        return result;
    }

    // Compute the cross product of two vectors
    template<typename T>
    vector<T> vectorCrossProduct(vector<T> vec1, vector<T> vec2) {
        if ((vec1.size() != vec2.size() || vec1.size() != 3) && !DISABLE_WARNINGS) {
            std::cout << "vectorCrossProduct | Warning : dimensions don't match. Both vectors must have 3 values." << std::endl;
        }
        while (vec1.size() < 3) {
            vec1.push_back(static_cast<T>(0.0));
        }
        while (vec2.size() < 3) {
            vec2.push_back(static_cast<T>(0.0));
        }
        
        return vector<T>{vec1[1]*vec2[2] - vec1[2] * vec2[1], vec1[2] * vec2[0] - vec1[0] * vec2[2], vec1[0] * vec2[1] - vec1[1] * vec2[0]};
    }

    // Compute the magnitude (Euclidean norm) of a vector
    template<typename T>
    T vectorMagnitude(vector<T> vec) {
        T result = static_cast<T>(0.0);
        for (size_t i = 0; i < vec.size(); i++) {
            result += vec[i] * vec[i];
        }
        return sqrt(result);
    }

    // Divide all the values in the vector by its magnitude
    template<typename T>
    vector<T> vectorNormalise(vector<T> vec) {
        T mag = vectorMagnitude<T>(vec);
        if (mag != static_cast<T>(0.0)) return vectorScale<T>(vec, static_cast<T>(1.0) / mag);
        else return vec;
    }

    // Compute the mean of all values in a vector
    template<typename T>
    T vectorAverage(vector<T> vec) {
        T result = static_cast<T>(0.0);
        size_t numOfValues = 0;
        for (size_t i = 0; i < vec.size(); i++) {
            result += vec[i];
            numOfValues++;
        }
        if (numOfValues == 0) {
            return static_cast<T>(0.0);
        }
        return result / static_cast<T>(numOfValues);
    }

    // Compute the average of two points
    template<typename T>
    vector<T> vectorAverage(vector<T> vec1, vector<T> vec2) {
        for (size_t i = 0; i < vec1.size(); i++) {
            vec1[i] += vec2[i];
        }
        return vectorScale<T>(vec1, static_cast<T>(0.5));
    }

    // Compute the average of multiple points stored in a matrix.
    // All points must have the same number of coordinates. The points must be
    // rows in the matrix.
    template<typename T>
    vector<T> vectorAverage(matrix<T> mtx) {
        size_t numOfVectors = mtx.size();
        size_t lengthOfVectors = mtx[0].size();
        vector<T> result = vectorOfZeros<T>(lengthOfVectors);
        for (size_t i = 0; i < numOfVectors; i++) {
            for (size_t j = 0; j < lengthOfVectors; j++) {
                result[j] += mtx[i][j];
            }
        }
        return vectorScale<T>(result, static_cast<T>(1.0) / static_cast<T>(numOfVectors));
    }

    // Compute the distance between two vectors. Use the coordinateSystem argument to specify
    // if the input vectors are cartesian (cart), polar (pol), or spherical (sph). 
    // By default it's assumed that both vectors are cartesian. 
    template<typename T>
    T vectorGetDistance(vector<T> vec1, vector<T> vec2, CoordSystem coordinateSystem) {
        size_t length = std::min(vec1.size(), vec2.size());
        if (coordinateSystem == CoordSystem::pol) {
            vec1 = polToCart<T>(vec1);
            vec2 = polToCart<T>(vec2);
            length = 2;
        }
        else if (coordinateSystem == CoordSystem::sph) {
            vec1 = sphToCart<T>(vec1);
            vec2 = sphToCart<T>(vec2);
            length = 3;
        }
        T sum = static_cast<T>(0.0);
        for (size_t i = 0; i < length; i++) {
            sum += pow(vec1[i] - vec2[i], 2);
        }
        return sqrt(sum);
    }

    // apply arbitrary function to all values in a vector
    template<typename T>
    vector<T> vectorFunc(vector<T> vec, std::function<T(T)> func) {
        for (size_t i = 0; i < vec.size(); i++) {
            vec[i] = func(vec[i]);
        }
        return vec;
    }

    // Applies quick sort to the vector and returns the sorted version in ascending order. 
    // EG: [0, 100, 50, 25, 75] => [0, 25, 50, 75, 100]
    template<typename T>
    vector<T> vectorQuickSort(vector<T> vec) {
        size_t length = vec.size();
        vector<utils::valIdx<T>> vec2;
        for (size_t i{}; i < length; i++) {
            vec2.push_back({ vec[i], {i, 0} });
        }
        utils::quickSort(vec2, 0, (int)length - 1);
        vector<T> result;
        for (size_t i{}; i < length; i++) {
            result.push_back(vec2[i].first);
        }
        return result;
    }

    // Sorts the vector, outputs the indices in order.
    // EG: [0, 100, 50, 25, 75] => [0, 3, 2, 4, 1]
    template<typename T>
    vector<size_t> vectorSortIndices(vector<T> vec){
        size_t length = vec.size();
        vector<utils::valIdx<T>> vec2;
        for (size_t i{}; i < length; i++){
            vec2.push_back({ vec[i], {i, 0} }); 
        }
        utils::quickSort(vec2, 0, (int)length - 1);
        vector<size_t> result;
        for (size_t i{}; i<length; i++){
            result.push_back(vec2[i].second.first);
        }
        return result;
    }

    //=================//
    //  MATRIX MATHS   //
    //=================//
 
    // create matrix of zeros
    template<typename T>
    matrix<T> matrixOfZeros(size_t rows, size_t cols)
    {
        matrix<T> newMatrix;
        for (size_t i = 0; i < rows; i++) {
            newMatrix.push_back(vector<T>(cols, static_cast<T>(0.0)));
        }
        return newMatrix;
    }

    // get the maximum value in a matrix
    template<typename T>
    T matrixMax(matrix<T> mtx) {
        T max = mtx[0][0];
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[i].size(); j++) {
                max = (mtx[i][j] > max) ? mtx[i][j] : max;
            }
        }
        return max;
    }

    // get the minimum value in a matrix
    template<typename T>
    T matrixMin(matrix<T> mtx) {
        T min = mtx[0][0];
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[i].size(); j++) {
                min = (mtx[i][j] < min) ? mtx[i][j] : min;
            }
        }
        return min;
    }

    // get the minimum and maximum values in a matrix with one operation.
    template<typename T>
    std::pair<T,T> matrixMinMax(matrix<T> mtx) {
        T min = mtx[0][0];
        T max = mtx[0][0];
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[i].size(); j++) {
                min = (mtx[i][j] < min) ? mtx[i][j] : min;
                max = (mtx[i][j] > max) ? mtx[i][j] : max;
            }
        }
        return std::pair<T,T>{ min, max };
    }

    // Element-wise mutliplication of two matrices
    template<typename T>
    matrix<T> matrixElemMult(matrix<T> mtx1, matrix<T> mtx2) {
        size_t nRows1 = mtx1.size();
        size_t nCols1 = mtx1[0].size();
        size_t nRows2 = mtx2.size();
        size_t nCols2 = mtx2[0].size();
        if (nRows1 != nRows2 || nCols1 != nCols2 && !DISABLE_WARNINGS) {
            std::cout << "elementMult | Warning : dimensions don't match." << std::endl;
        }
        size_t opRows = std::min(nRows1, nRows2);
        size_t opCols = std::min(nCols1, nCols2);
        for (size_t i = 0; i < opRows; i++) {
            for (size_t j = 0; j < opCols; j++) {
                mtx1[i][j] *= mtx2[i][j];
            }
        }
        return mtx1;
    }

    // Element-wise addition of two matrices
    template<typename T>
    matrix<T> matrixElemAdd(matrix<T> mtx1, matrix<T> mtx2) {
        size_t nRows1 = mtx1.size();
        size_t nCols1 = mtx1[0].size();
        size_t nRows2 = mtx2.size();
        size_t nCols2 = mtx2[0].size();
        if (nRows1 != nRows2 || nCols1 != nCols2 && !DISABLE_WARNINGS) {
            std::cout << "elementAdd | Warning : dimensions don't match." << std::endl;
        }
        size_t opRows = std::min(nRows1, nRows2);
        size_t opCols = std::min(nCols1, nCols2);
        for (size_t i = 0; i < opRows; i++) {
            for (size_t j = 0; j < opCols; j++) {
                mtx1[i][j] += mtx2[i][j];
            }
        }
        return mtx1;
    }

    // Multiply matrix by scalar
    template<typename T>
    matrix<T> matrixScale(matrix<T> mtx, T scalar) {
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[i].size(); j++) {
                mtx[i][j] *= scalar;
            }
        }
        return mtx;
    }

    // Multiply two matrices using regular matrix multiplication rules
    template<typename T>
    matrix<T> matrixMultiply(matrix<T> mtx1, matrix<T> mtx2) {
        size_t nRows1 = mtx1.size();
        size_t nCols1 = mtx1[0].size();
        size_t nRows2 = mtx2.size();
        size_t nCols2 = mtx2[0].size();
        if (nCols1 != nRows2 && !DISABLE_WARNINGS) {
            std::cout << "matrixMultiply | Warning : dimensions not fit for matrix multiplication. Returning first argument." << std::endl;
            return mtx1;
        }
        matrix<T> result = matrixOfZeros<T>(nRows1, nCols2);

        matrix<T> mtxT = matrixTranspose<T>(mtx2);

        for (size_t i = 0; i < nRows1; i++) {
            for (size_t j = 0; j < nCols2; j++) {
                result[i][j] = vectorDotProduct<T>(mtx1[i], mtxT[j]);
            }
        }

        return result;
    }

    // Transpose a matrix
    template<typename T>
    matrix<T> matrixTranspose(matrix<T> mtx) {
        size_t nCols = mtx.size();
        size_t nRows = mtx[0].size();
        matrix<T> result = matrixOfZeros<T>(nRows, nCols);
        for (size_t i = 0; i < nRows; i++) {
            for (size_t j = 0; j < nCols; j++) {
                result[i][j] = mtx[j][i];
            }
        }
        return result;
    }

    // Linear mapping of all values to fit in the desired range, by default [0,1]
    template<typename T>
    matrix<T> matrixRemap(matrix<T> mtx, T min, T max){
        std::pair<T,T> minMax = matrixMinMax<T>(mtx);
        if ((minMax.second - minMax.first) != static_cast<T>(0.0)){
            T m = (max - min) / (minMax.second - minMax.first);
            T b = min - minMax.first * m;
            return matrixFunc<T>(mtx, [m,b](T x) { return m*x+b; });
        }
        return matrixFunc<T>(mtx, [](T x) { return static_cast<T>(0.0); });
    }

    // Compute the average of all values in a matrix
    template<typename T>
    T matrixAverage(matrix<T> mtx) {
        T result = static_cast<T>(0.0);
        size_t numOfValues = 0;
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[j].size(); j++) {
                result += mtx[i][j];
                numOfValues++;
            }
        }
        if (numOfValues == 0) {
            return static_cast<T>(0.0);
        }
        return result / static_cast<T>(numOfValues);
    }

    // apply arbitrary function to all values in a matrix
    template<typename T>
    matrix<T> matrixFunc(matrix<T> mtx, std::function<T(T)> func) {
        for (size_t i = 0; i < mtx.size(); i++) {
            for (size_t j = 0; j < mtx[i].size(); j++) {
                mtx[i][j] = func(mtx[i][j]);
            }
        }
        return mtx;
    }

    // Substract the smallest value in a matrix from all its element
    template<typename T>
    matrix<T> matrixShiftToZero(matrix<T> mtx) {
        T correction = matrixMin<T>(mtx);
        if (correction != static_cast<T>(0.0)){
            for (size_t i = 0; i < mtx.size(); i++) {
                for (size_t j = 0; j < mtx[i].size(); j++) {
                    mtx[i][j] -= correction;
                }
            }
        }
        return mtx;
    }

    // Apply quick sort to a matrix.
    template<typename T>
    matrix<T> matrixQuickSort(matrix<T> mtx) {
        size_t numRows = mtx.size();
        size_t numCols = mtx[0].size();
        size_t length = numRows * numCols;
        vector<utils::valIdx<T>> vec;
        for (size_t i{}; i < mtx.size(); i++) {
            for (size_t j{}; j < mtx[i].size(); j++) {
                vec.push_back({ mtx[i][j], {i, j} });
            }
        }
        utils::quickSort(vec, 0, (int)length - 1);
        matrix<T> result;
        size_t idx = 0;
        for (size_t i{}; i < length; i++) {
            for (size_t j{}; j < mtx[i].size(); j++) {
                result[i][j] = vec[idx++].first;
            }
        }
        return result;
    }

    // Sort a matrix, output the indices in order.
    template<typename T>
    vector<mtxIndex> matrixSortIndices(matrix<T> mtx){
        size_t length = mtx.size() * mtx[0].size();
        vector<utils::valIdx<T>> vec;
        for (size_t i{}; i < mtx.size(); i++){
            for (size_t j{}; j < mtx[i].size(); j++) {
                vec.push_back({ mtx[i][j], {i, j} });
            }
        }
        utils::quickSort(vec, 0, (int)length - 1);
        vector<mtxIndex> result;
        for (size_t i{}; i<length; i++){
            result.push_back({ vec[i].second.first, vec[i].second.second });
        }
        return result;
    }

    //==========================================//
    //  CONVERSION BETWEEN COORDINATE SYSTEMS   //
    //==========================================//

    // Convert polar coordinates to a cartesian vector
    template<typename T>
    vector<T> polToCart(T theta, T radius) {
        vector<T> result = vectorOfZeros<T>(2);
        result[0] = radius * cos(theta);
        result[1] = radius * sin(theta);
        return result;
    }

    // Convert polar vector to cartesian vector
    template<typename T>
    vector<T> polToCart(vector<T> vec) {
        vector<T> result = vectorOfZeros<T>(2);
        result[0] = vec[1] * cos(vec[0]);
        result[1] = vec[1] * sin(vec[0]);
        return result;
    }

    // Convert polar coordinates stored in separate vectors to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> polToCart(vector<T> theta, vector<T> radius, bool transpose) {
        matrix<T> result = matrixOfZeros<T>(theta.size(), 2);
        for (size_t i = 0; i < theta.size(); i++) {
            result[i][0] = radius[i] * cos(theta[i]);
            result[i][1] = radius[i] * sin(theta[i]);
        }
        if (transpose) {
            result = matrixTranspose<T>(result);
        }
        return result;
    }

    // Convert polar vectors stored in a matrix, to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> polToCart(matrix<T> mtx, bool transpose) {
        matrix<T> result = matrixOfZeros<T>(mtx.size(), 2);
        for (size_t i = 0; i < mtx.size(); i++) {
            result[i][0] = mtx[1][i] * cos(mtx[0][i]);
            result[i][1] = mtx[1][i] * sin(mtx[0][i]);
        }
        if (transpose) {
            result = matrixTranspose<T>(result);
        }
        return result;
    }

    // convert spherical coordinates to a cartesian vector
    template<typename T>
    vector<T> sphToCart(T theta, T phi, T radius) {
        vector<T> result = vectorOfZeros<T>(3);
        T rCosElev = radius * cos(phi);
        result[0] = rCosElev * cos(theta);
        result[1] = rCosElev * sin(theta);
        result[2] = radius * sin(phi);
        return result;
    }

    // convert spherical vector to cartesian vector
    template<typename T>
    vector<T> sphToCart(vector<T> vec) {
        vector<T> result = vectorOfZeros<T>(3);
        T rCosElev = vec[2] * cos(vec[1]);
        result[0] = rCosElev * cos(vec[0]);
        result[1] = rCosElev * sin(vec[0]);
        result[2] = vec[2] * sin(vec[1]);
        return result;
    }

    // Convert spherical coordinates stored in separate vectors, to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> sphToCart(vector<T> theta, vector<T> phi, vector<T> radius, bool transpose) {
        matrix<T> result = matrixOfZeros<T>(theta.size(), 3);
        for (size_t i = 0; i < theta.size(); i++) {
            T rCosElev = radius[i] * cos(phi[i]);
            result[i][0] = rCosElev * cos(theta[i]);
            result[i][1] = rCosElev * sin(theta[i]);
            result[i][2] = radius[i] * sin(phi[i]);
        }
        if (transpose) {
            result = matrixTranspose<T>(result);
        }
        return result;
    }

    // Convert spherical vectors stored in a matrix to cartesian vectors stored
    // in a matrix.
    template<typename T>
    matrix<T> sphToCart(matrix<T> mtx, bool transpose) {
        matrix<T> result = matrixOfZeros<T>(mtx.size(), 3);
        for (size_t i = 0; i < mtx.size(); i++) {
            T rCosElev = mtx[i][2] * cos(mtx[i][1]);
            result[i][0] = rCosElev * cos(mtx[i][0]);
            result[i][1] = rCosElev * sin(mtx[i][0]);
            result[i][2] = mtx[i][2] * sin(mtx[i][1]);
        }
        if (transpose) {
            result = matrixTranspose<T>(result);
        }
        return result;
    }

    // Wow, you reached the end of the file. Here's a secret function...
    /*namespace scrt {
        void matrixFun() {
           const char d=32;matrix<char>m=arrmath::matrixOfZeros<char>(9,14);for(size_t i{};i<m.size();i++){for(size_t j{};j<m[i].size();
            j++){m[i][j]=d*((int)(((i==1||i==3)&&(j==3||j==10))||((i==2||i==6)&&((j>1&&j<5)||(j>8&&j<12)))||((i==5)&&(j==1||j==2||
            j==11||j==12))||((i==7)&&(j>3&&j<10)))+1);std::cout<<m[i][j];}std::cout<<std::endl;}
        }
    }*/
}
