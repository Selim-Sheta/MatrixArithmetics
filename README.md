# MatrixArithmetics
A header-only "micro library" containing useful basic functions for vector and matrix manipulation.

Functions starting with "vector" can only be applied to a vector, and functions starting with "matrix" can only be applied to a matrix. 

The library is designed for maximum flexibility, prioritising no compile errors over maths riggor. For example, you can do a dot product between two vectors of different sizes. The result will obviously be wrong, but you won't get errors or warnings.

Notes: 
- Strongly recommended: "using namespace arrmath" at the top of your file.
- This library is not particularly fast, to be avoided for time-critical operations.
- This library hasn't been fully tested as of 26/11/2022

## Functions

### Vector Maths:
- vectorOfZeros : create a vector of zeros
- vectorMax : find largest value
- vectorMin : find smallest value
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
- vectorAverage : compute average of two or more points in space
- vectorGetDistance : compute distance between two points in space

### Matrix Maths:
- matrixOfZeros : create a matrix of zeros
- matrixMax : find largest value
- matrixMin : find smallest value
- matrixElemMult : element-wise multiplication of two matrices
- matrixElemAdd : element-wise addition of two vectors
- matrixScale : scale by a constant
- matrixMultiply : multiply two matrices following standard matrix-multiplication rules
- matrixTranspose : tranpose matrix

### Coordinate System Conversion:
- polToCart : convert polar coordinates to cartesian
- sphToCart : convert spherical coordinates to cartesian

To be added:
- vectorQuickSort
- vectorCrossProduct
- matrixSpan
- matrixMaxIndex
- matrixMinIndex
- matrixRotate
- cartToPol
- cartToSph
