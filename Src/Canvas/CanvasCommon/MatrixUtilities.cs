using System;

namespace CanvasCommon
{
    public static class MatrixUtilities
    {
        /// <summary>
        /// Allocates/creates a 2D matrix initialized to all 0.0, assume rows and cols > 0
        /// </summary>
        public static double[][] MatrixCreate(int rows, int cols)
        {
            
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

        /// <summary>
        /// Allocates/creates a 3D matrix initialized to all 0.0, assume rows and cols > 0
        /// </summary>
        public static double[][][] MatrixCreate(int rows, int cols1, int cols2)
        {

            double[][][] result = new double[rows][][];
            for (int i = 0; i < rows; ++i) { 
                result[i] = new double[cols1][];
                for (int j = 0; j < cols1; ++j)
                    result[i][j] = new double[cols2];
            }
            return result;
        }

        /// <summary>
        /// Allocates/creates a duplicate of a matrix assumes matrix is not null
        /// </summary>
        public static double[][] MatrixDuplicate(double[][] matrix)
        {          
            double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; ++i) // copy the values
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }

        /// <summary>
        /// LUP decomposition with pivoting
        /// </summary>
        /// <returns>
        /// result is L (with 1s on diagonal) and U; perm holds row permutations; 
        /// toggle is +1 or -1 (even or odd)
        /// </returns>
        public static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
        {
            // rerturns: 
            int rows = matrix.Length;
            int cols = matrix[0].Length; // assume all rows have the same number of columns so just use row [0].
            if (rows != cols)
                throw new Exception("Attempt to MatrixDecompose a non-square mattrix");

            int n = rows; 
            double[][] result = MatrixDuplicate(matrix); 

            perm = new int[n]; // set up row permutation result
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1; // toggle tracks row swaps. +1 -> even, -1 -> odd. used by MatrixDeterminant

            for (int j = 0; j < n - 1; ++j) // each column
            {
                double colMax = Math.Abs(result[j][j]); // find largest value in col j
                int pRow = j;
                for (int i = j + 1; i < n; ++i)
                {
                    if (result[i][j] > colMax)
                    {
                        colMax = result[i][j];
                        pRow = i;
                    }
                }

                if (pRow != j) // if largest value not on pivot, swap rows
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[pRow]; // and swap perm info
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }

                if (Math.Abs(result[j][j]) < 1.0E-20) // if diagonal after swap is zero . . .
                    return null; // consider a throw

                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i][k] -= result[i][j] * result[j][k];
                    }
                }
            } 

            return result;
        }

        public static double MatrixDeterminant(double[][] matrix)
        {
            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
            if (lum == null)
                throw new Exception("Unable to compute MatrixDeterminant");
            double result = toggle;
            for (int i = 0; i < lum.Length; ++i)
                result *= lum[i][i];
            return result;
        }

        public static double[][] MatrixInverse(double[][] matrix)
        {
            int n = matrix.Length;
            double[][] result = MatrixDuplicate(matrix);

            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b); // 

                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }

        /// <summary>
        /// before calling this helper, permute b using the perm array from MatrixDecompose that generated lu Matrix
        /// </summary>
        public static double[] HelperSolve(double[][] luMatrix, double[] b) // helper
        {
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }

            return x;
        }

        /// <summary>
        ///  result of multiplying an n x m matrix by a m x 1 column vector (yielding an n x 1 column vector)
        /// </summary>
        public static double[] MatrixVectorProduct(double[][] matrix, double[] vector)
        {           
            int mRows = matrix.Length;
            int mCols = matrix[0].Length;
            int vRows = vector.Length;
            if (mCols != vRows)
                throw new Exception("Non-conformable matrix and vector in MatrixVectorProduct");
            double[] result = new double[mRows]; // an n x m matrix times a m x 1 column vector is a n x 1 column vector
            for (int i = 0; i < mRows; ++i)
                for (int j = 0; j < mCols; ++j)
                    result[i] += matrix[i][j] * vector[j];
            return result;
        }
    }
}