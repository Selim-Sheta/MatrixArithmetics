# MatrixArithmetics
A header-only "micro library" containing useful basic functions for vector and matrix manipulation.

Functions starting with "vector" can only be applied to a vector, and functions starting with "matrix" can only be applied to a matrix. 

The library is designed for maximum flexibility, prioritising no compile errors or warning.

Notes: 
- recommended: "using namespace arrmath" at the top of your file.
- This library is not particularly fast, to be avoided for time-critical operations.
- This library hasn't been fully tested as of 14/12/2022
- to be tested : matrix functions, coord system functions

## Functions

### Vector Maths:
- vectorOfZeros : create a vector of zeros
- vectorMax : find largest value
- vectorMin : find smallest value
- vectorMinMax : get smallest and larget value as an std::pair
- vectorSpan : find difference between largest and smallest value
- vectorMaxIndex : find index of largest value
- vectorMinIndex : find index of smallest value
- vectorZeroAvg : apply zero-averaging
- vectorElemMult : element-wise multiplication of two vectors
- vectorElemAdd : element-wise addition of two vectors
- vectorScale : scale by a constant
- vectorTranspose : convert to column vector
- vectorDotProduct : compute dot product
- vectorMagnitude : compute Euclidean norm (magnitude)
- vectorNormalise : divide all values in the vector by its magnitude
- vectorAverage : compute average of two or more points in space
- vectorGetDistance : compute distance between two points in space
- vectorFunc : apply arbitrary lambda function to all elements in a vector
- vectorSortIndices : sort vector and return indices arranged in order

### Matrix Maths:
- matrixOfZeros : create a matrix of zeros
- matrixMax : find largest value
- matrixMin : find smallest value
- matrixMinMax : get smallest and larget value as an std::pair
- matrixElemMult : element-wise multiplication of two matrices
- matrixElemAdd : element-wise addition of two vectors
- matrixScale : scale by a constant
- matrixMultiply : multiply two matrices following standard matrix-multiplication rules
- matrixTranspose : tranpose matrix
- matrixFunc : apply arbitrary lambda function to all elements in a matrix
- matrixShiftToZero : substract the smallest value in a matrix from all its element
- matrixSortIndices : sort vector and return indices arranged in order

### Coordinate System Conversion:
- polToCart : convert polar coordinates to cartesian
- sphToCart : convert spherical coordinates to cartesian

To be added soon:
- vectorQuickSort
- vectorCrossProduct
- vectorAvg
- vectorShiftToZero
- vectorAngleBetween
- matrixQuickSort
- matrixMaxIndex
- matrixMinIndex
- matrixSpan
- matrixAvg
- matrixZeroAvg
- matrixRotate
- matrixInvert
- matrixFlip
- matrixIdentity
- matrixDeterm
- cartToPol
- cartToSph
- polToSph
- sphToPol
