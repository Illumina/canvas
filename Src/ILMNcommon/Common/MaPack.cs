using System;
using System.Text;

// -------------------------------------------------------------------------
// Lutz Roeder's .NET Mapack, adapted from Mapack for COM and Jama routines.
// Copyright (C) 2001-2003 Lutz Roeder. All rights reserved.
// http://www.aisto.com/roeder/dotnet
// roeder@aisto.com
// -------------------------------------------------------------------------

namespace Illumina.Common.LinearAlgebra
{
    /// <summary>Matrix provides the fundamental operations of numerical linear algebra.</summary>
    public interface IMatrix
    {
        /// <summary>Returns the number of columns.</summary>
        int Rows { get; }

        /// <summary>Returns the number of columns.</summary>
        int Columns { get; }

        /// <summary>Access the value at the given location.</summary>
        double this[int i, int j] { get; set; }

        /// <summary>Inverse of the matrix if matrix is square, pseudoinverse otherwise.</summary>
        IMatrix Inverse { get; }

        /// <summary>Determinant if matrix is square.</summary>
        double Determinant { get; }

        /// <summary>Returns the One Norm for the matrix.</summary>
        /// <value>The maximum column sum.</value>
        double Norm1 { get; }

        double Norm { get; }

        /// <summary>Returns the Infinity Norm for the matrix.</summary>
        /// <value>The maximum row sum.</value>
        double InfinityNorm { get; }

        /// <summary>Returns the Frobenius Norm for the matrix.</summary>
        /// <value>The square root of sum of squares of all elements.</value>
        double FrobeniusNorm { get; }

        /// <summary>
        ///     Return <see langword="true" /> if the matrix is a square matrix.
        /// </summary>
        bool IsSquare { get; }

        /// <summary>
        ///     Returns <see langword="true" /> if the matrix is symmetric.
        /// </summary>
        bool IsSymmetric { get; }

        /// <summary>Returns the trace of the matrix.</summary>
        /// <returns>Sum of the diagonal elements.</returns>
        double Trace { get; }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="startRow">Start row index.</param>
        /// <param name="endRow">End row index;</param>
        /// <param name="startColumn">Start column index;</param>
        /// <param name="endColumn">End column index;</param>
        IMatrix Submatrix(int startRow, int endRow, int startColumn, int endColumn);

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="r">Array of row indices;</param>
        /// <param name="c">Array of row indices;</param>
        IMatrix Submatrix(int[] r, int[] c);

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="startRow">Starttial row index.</param>
        /// <param name="endRow">End row index.</param>
        /// <param name="c">Array of row indices.</param>
        IMatrix Submatrix(int startRow, int endRow, int[] c);

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="r">Array of row indices.</param>
        /// <param name="startColumn">Start column index.</param>
        /// <param name="endColumn">End column index.</param>
        IMatrix Submatrix(int[] r, int startColumn, int endColumn);

        /// <summary>Creates a copy of the matrix.</summary>
        IMatrix Clone();

        /// <summary>Returns the transposed matrix.</summary>
        IMatrix Transpose();

        /// <summary>Matrix addition.</summary>
        IMatrix Addition(IMatrix b);

        IMatrix Sum();

        /// <summary>Matrix-matrix multiplication.</summary>
        IMatrix Multiply(IMatrix b);

        IMatrix DotMultiply(IMatrix b);

        /// <summary>Matrix-scalar multiplication.</summary>
        IMatrix Multiply(double s);

        /// <summary>Matrix subtraction.</summary>
        IMatrix Subtraction(IMatrix b);

        /// <summary>Returns the LHS solution vetor if the matrix is square or the least squares solution otherwise.</summary>
        IMatrix Solve(IMatrix rhs);

        IMatrix Divide(IMatrix rhs);

        /// <summary>Returns the cholesky decomposition for this matrix.</summary>
        ICholeskyDecomposition GetCholeskyDecomposition();

        /// <summary>Returns the LU decomposition for this matrix.</summary>
        ILuDecomposition GetLuDecomposition();

        /// <summary>Returns the singular value decomposition for this matrix.</summary>
        ISingularValueDecomposition GetSingularValueDecomposition();

        /// <summary>Returns the QR decomposition for this matrix.</summary>
        IQrDecomposition GetQrDecomposition();

        /// <summary>Returns the eigenvalue decomposition for this matrix.</summary>
        IEigenvalueDecomposition GetEigenvalueDecomposition();
    }

    /// <summary>
    ///     Cholesky Decomposition of a symmetric, positive definite matrix.
    /// </summary>
    /// <remarks>
    ///     For a symmetric, positive definite matrix <c>A</c>, the Cholesky decomposition is a
    ///     lower triangular matrix <c>L</c> so that <c>A = L * L'</c>.
    ///     If the matrix is not symmetric or positive definite, the constructor returns a partial
    ///     decomposition and sets two internal variables that can be queried using the
    ///     <see cref="IsSymmetric" /> and <see cref="IsPositiveDefinite" /> properties.
    /// </remarks>
    public interface ICholeskyDecomposition
    {
        /// <summary>
        ///     Returns <see langword="true" /> if the matrix is positive definite.
        /// </summary>
        Boolean IsPositiveDefinite { get; }

        /// <summary>
        ///     Returns <see langword="true" /> if the matrix is symmetric.
        /// </summary>
        Boolean IsSymmetric { get; }

        /// <summary>
        ///     Returns the left triangular factor <c>L</c> so that <c>A = L * L'</c>.
        /// </summary>
        IMatrix LeftTriangularFactor { get; }

        /// <summary>
        ///     Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// <param name="rhs">
        ///     Right hand side matrix with as many rows as <c>A</c> and any number of columns.
        /// </param>
        /// <returns>
        ///     Matrix <c>X</c> so that <c>L * L' * X = B</c>.
        /// </returns>
        /// <exception cref="T:System.ArgumentException">Matrix dimensions do not match.</exception>
        /// <exception cref="T:System.InvalidOperationException">Matrix is not symmetrix and positive definite.</exception>
        IMatrix Solve(IMatrix rhs);
    }

    /// <summary>
    ///     LU decomposition of a rectangular matrix.
    /// </summary>
    /// <remarks>
    ///     For an m-by-n matrix <c>A</c> with m >= n, the LU decomposition is an m-by-n
    ///     unit lower triangular matrix <c>L</c>, an n-by-n upper triangular matrix <c>U</c>,
    ///     and a permutation vector <c>piv</c> of length m so that <c>A(piv)=L*U</c>.
    ///     If m &lt; n, then <c>L</c> is m-by-m and <c>U</c> is m-by-n.
    ///     The LU decompostion with pivoting always exists, even if the matrix is
    ///     singular, so the constructor will never fail.  The primary use of the
    ///     LU decomposition is in the solution of square systems of simultaneous
    ///     linear equations. This will fail if <see cref="IsNonSingular" /> returns <see langword="false" />.
    /// </remarks>
    public interface ILuDecomposition
    {
        /// <summary>Returns if the matrix is non-singular.</summary>
        bool IsNonSingular { get; }

        /// <summary>Returns the determinant of the matrix.</summary>
        double Determinant { get; }

        /// <summary>
        ///     Returns the lower triangular factor <c>L</c> with <c>A=LU</c>.
        /// </summary>
        IMatrix LowerTriangularFactor { get; }

        /// <summary>
        ///     Returns the lower triangular factor <c>L</c> with <c>A=LU</c>.
        /// </summary>
        IMatrix UpperTriangularFactor { get; }

        /// <summary>Returns the pivot permuation vector.</summary>
        double[] PivotPermutationVector { get; }

        /// <summary>
        ///     Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// <param name="b">
        ///     Right hand side matrix with as many rows as <c>A</c> and any number of columns.
        /// </param>
        /// <returns>
        ///     Matrix <c>X</c> so that <c>L * U * X = B</c>.
        /// </returns>
        IMatrix Solve(IMatrix b);
    }

    /// <summary>
    ///     QR decomposition for a rectangular matrix.
    /// </summary>
    /// <remarks>
    ///     For an m-by-n matrix <c>A</c> with <c>m &gt;= n</c>, the QR decomposition is an m-by-n
    ///     orthogonal matrix <c>Q</c> and an n-by-n upper triangular
    ///     matrix <c>R</c> so that <c>A = Q * R</c>.
    ///     The QR decompostion always exists, even if the matrix does not have
    ///     full rank, so the constructor will never fail.  The primary use of the
    ///     QR decomposition is in the least squares solution of nonsquare systems
    ///     of simultaneous linear equations.
    ///     This will fail if <see cref="IsFullRank" /> returns <see langword="false" />.
    /// </remarks>
    public interface IQrDecomposition
    {
        /// <summary>
        ///     Shows if the matrix <c>A</c> is of full rank.
        /// </summary>
        /// <value>
        ///     The value is <see langword="true" /> if <c>R</c>, and hence <c>A</c>, has full rank.
        /// </value>
        bool IsFullRank { get; }

        /// <summary>
        ///     Returns the upper triangular factor <c>R</c>.
        /// </summary>
        IMatrix UpperTriangularFactor { get; }

        /// <summary>
        ///     Returns the orthogonal factor <c>Q</c>.
        /// </summary>
        IMatrix OrthogonalFactor { get; }

        /// <summary>
        ///     Least squares solution of <c>A * X = B</c>
        /// </summary>
        /// <param name="rhs">
        ///     Right-hand-side matrix with as many rows as <c>A</c> and any number of columns.
        /// </param>
        /// <returns>
        ///     A matrix that minimized the two norm of <c>Q * R * X - B</c>.
        /// </returns>
        /// <exception cref="T:System.ArgumentException">Matrix row dimensions must be the same.</exception>
        /// <exception cref="T:System.InvalidOperationException">Matrix is rank deficient.</exception>
        IMatrix Solve(IMatrix rhs);
    }

    /// <summary>
    ///     Singular Value Decomposition for a rectangular matrix.
    /// </summary>
    /// <remarks>
    ///     For an m-by-n matrix <c>A</c> with <c>m >= n</c>, the singular value decomposition is
    ///     an m-by-n orthogonal matrix <c>U</c>, an n-by-n diagonal matrix <c>S</c>, and
    ///     an n-by-n orthogonal matrix <c>V</c> so that <c>A = U * S * V'</c>.
    ///     The singular values, <c>sigma[k] = S[k,k]</c>, are ordered so that
    ///     <c>sigma[0] >= sigma[1] >= ... >= sigma[n-1]</c>.
    ///     The singular value decompostion always exists, so the constructor will
    ///     never fail. The matrix condition number and the effective numerical
    ///     rank can be computed from this decomposition.
    /// </remarks>
    public interface ISingularValueDecomposition
    {
        /// <summary>
        ///     Returns the condition number <c>max(S) / min(S)</c>.
        /// </summary>
        double Condition { get; }

        /// <summary>Returns the Two norm.</summary>
        double Norm2 { get; }

        /// <summary>Returns the effective numerical matrix rank.</summary>
        /// <value>Number of non-negligible singular values.</value>
        int Rank { get; }

        /// <summary>Return the one-dimensional array of singular values.</summary>
        double[] Diagonal { get; }
    }

    /// <summary>
    ///     Determines the eigenvalues and eigenvectors of a real square matrix.
    /// </summary>
    /// <remarks>
    ///     If <c>A</c> is symmetric, then <c>A = V * D * V'</c> and <c>A = V * V'</c>
    ///     where the eigenvalue matrix <c>D</c> is diagonal and the eigenvector matrix <c>V</c> is orthogonal.
    ///     If <c>A</c> is not symmetric, the eigenvalue matrix <c>D</c> is block diagonal
    ///     with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    ///     <c>lambda+i*mu</c>, in 2-by-2 blocks, <c>[lambda, mu; -mu, lambda]</c>.
    ///     The columns of <c>V</c> represent the eigenvectors in the sense that <c>A * V = V * D</c>.
    ///     The matrix V may be badly conditioned, or even singular, so the validity of the equation
    ///     <c>A=V*D*inverse(V)</c> depends upon the condition of <c>V</c>.
    /// </remarks>
    public interface IEigenvalueDecomposition
    {
        /// <summary>Returns the real parts of the eigenvalues.</summary>
        double[] RealEigenvalues { get; }

        /// <summary>Returns the imaginary parts of the eigenvalues.</summary>
        double[] ImaginaryEigenvalues { get; }

        /// <summary>Returns the eigenvector matrix.</summary>
        IMatrix EigenvectorMatrix { get; }

        /// <summary>Returns the block diagonal eigenvalue matrix.</summary>
        IMatrix DiagonalMatrix { get; }
    }

    /// <summary>Matrix provides the fundamental operations of numerical linear algebra.</summary>
    public class Matrix : IMatrix
    {
        //public static byte[] K = new byte[] { 0x7F, 0x0A, 0x49, 0x8D, 0xD1, 0xD8, 0x19, 0xAB };
        private readonly int _columns;
        private readonly double[][] _data;
        private readonly int _rows;

        /// <summary>Constructs an empty matrix of the given size.</summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        public Matrix(int rows, int columns)
        {
            _rows = rows;
            _columns = columns;
            _data = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                _data[i] = new double[columns];
            }
        }

        /// <summary>Constructs a matrix of the given size and assigns a given value to all diagonal elements.</summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="value">Value to assign to the diagnoal elements.</param>
        public Matrix(int rows, int columns, double value)
        {
            _rows = rows;
            _columns = columns;
            _data = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                _data[i] = new double[columns];
            }
            for (int i = 0; i < rows; i++)
            {
                _data[i][i] = value;
            }
        }

        public Matrix(float[] data, int rowCount, int columnCount)
        {
            if (rowCount*columnCount != data.Length) throw new ArgumentException();
            _rows = rowCount;
            _columns = columnCount;
            _data = new double[_rows][];
            int overallIndex = 0;
            for (int y = 0; y < _rows; y++)
            {
                _data[y] = new double[_columns];
                for (int x = 0; x < _columns; x++)
                {
                    _data[y][x] = data[overallIndex++];
                }
            }
        }

        /// <summary>Constructs a matrix from the given array.</summary>
        /// <param name="data">The array the matrix gets constructed from.</param>
        public Matrix(double[][] data)
        {
            _rows = data.Length;
            _columns = data[0].Length;

            for (int i = 0; i < _rows; i++)
            {
                if (data[i].Length != _columns)
                {
                    throw new ArgumentException();
                }
            }

            this._data = data;
        }

        private double[][] Array
        {
            get { return _data; }
        }

        /// <summary>Returns the number of columns.</summary>
        public int Rows
        {
            get { return _rows; }
        }

        /// <summary>Returns the number of columns.</summary>
        public int Columns
        {
            get { return _columns; }
        }

        /// <summary>
        ///     Return <see langword="true" /> if the matrix is a square matrix.
        /// </summary>
        public bool IsSquare
        {
            get { return (_rows == _columns); }
        }

        /// <summary>
        ///     Returns <see langword="true" /> if the matrix is symmetric.
        /// </summary>
        public bool IsSymmetric
        {
            get
            {
                if (IsSquare)
                {
                    for (int i = 0; i < _rows; i++)
                    {
                        for (int j = 0; j <= i; j++)
                        {
                            if (_data[i][j] != _data[j][i])
                            {
                                return false;
                            }
                        }
                    }

                    return true;
                }

                return false;
            }
        }

        /// <summary>Access the value at the given location.</summary>
        public double this[int i, int j]
        {
            set { _data[i][j] = value; }
            get { return _data[i][j]; }
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="i0">Starttial row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        public IMatrix Submatrix(int i0, int i1, int j0, int j1)
        {
            Matrix xx = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
            double[][] x = xx.Array;
            for (int i = i0; i <= i1; i++)
                for (int j = j0; j <= j1; j++)
                    x[i - i0][j - j0] = _data[i][j];
            return xx;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="c">Array of row indices</param>
        public IMatrix Submatrix(int[] r, int[] c)
        {
            Matrix xx = new Matrix(r.Length, c.Length);
            double[][] x = xx.Array;
            for (int i = 0; i < r.Length; i++)
                for (int j = 0; j < c.Length; j++)
                    x[i][j] = _data[r[i]][c[j]];
            return xx;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="i0">Starttial row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="c">Array of row indices</param>
        public IMatrix Submatrix(int i0, int i1, int[] c)
        {
            Matrix xx = new Matrix(i1 - i0 + 1, c.Length);
            double[][] x = xx.Array;
            for (int i = i0; i <= i1; i++)
                for (int j = 0; j < c.Length; j++)
                    x[i - i0][j] = _data[i][c[j]];
            return xx;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        public IMatrix Submatrix(int[] r, int j0, int j1)
        {
            Matrix xx = new Matrix(r.Length, j1 - j0 + 1);
            double[][] x = xx.Array;
            for (int i = 0; i < r.Length; i++)
                for (int j = j0; j <= j1; j++)
                    x[i][j - j0] = _data[r[i]][j];
            return xx;
        }

        /// <summary>Creates a copy of the matrix.</summary>
        public IMatrix Clone()
        {
            Matrix xx = new Matrix(_rows, _columns);
            double[][] x = xx.Array;
            for (int i = 0; i < _rows; i++)
                for (int j = 0; j < _columns; j++)
                    x[i][j] = _data[i][j];
            return xx;
        }

        /// <summary>Returns the transposed matrix.</summary>
        public IMatrix Transpose()
        {
            Matrix xx = new Matrix(_columns, _rows);
            double[][] x = xx.Array;
            for (int i = 0; i < _rows; i++)
                for (int j = 0; j < _columns; j++)
                    x[j][i] = _data[i][j];
            return xx;
        }


        public double Norm
        {
            get
            {
                ISingularValueDecomposition svd = GetSingularValueDecomposition();
                double[] vals = svd.Diagonal;
                double max = vals[0];
                for (int i = 1; i < vals.Length; i++)
                    if (vals[i] > max)
                        max = vals[i];
                return max;
            }
        }

        /// <summary>Returns the One Norm for the matrix.</summary>
        /// <value>The maximum column sum.</value>
        public double Norm1
        {
            get
            {
                double f = 0;
                for (int j = 0; j < _columns; j++)
                {
                    double s = 0;
                    for (int i = 0; i < _rows; i++)
                        s += Math.Abs(_data[i][j]);
                    f = Math.Max(f, s);
                }
                return f;
            }
        }

        /// <summary>Returns the Infinity Norm for the matrix.</summary>
        /// <value>The maximum row sum.</value>
        public double InfinityNorm
        {
            get
            {
                double f = 0;
                for (int i = 0; i < _rows; i++)
                {
                    double s = 0;
                    for (int j = 0; j < _columns; j++)
                        s += Math.Abs(_data[i][j]);
                    f = Math.Max(f, s);
                }
                return f;
            }
        }

        /// <summary>Returns the Frobenius Norm for the matrix.</summary>
        /// <value>The square root of sum of squares of all elements.</value>
        public double FrobeniusNorm
        {
            get
            {
                double f = 0;
                for (int i = 0; i < _rows; i++)
                {
                    for (int j = 0; j < _columns; j++)
                    {
                        f = MathHelper.Hypotenuse(f, _data[i][j]);
                    }
                }

                return f;
            }
        }

        /// <summary>Matrix addition.</summary>
        public IMatrix Addition(IMatrix b)
        {
            if ((_rows != b.Rows) || (_columns != b.Columns))
            {
                throw new ArgumentException("Matrix dimension do not match.");
            }
            Matrix xx = new Matrix(_rows, _columns);
            double[][] x = xx.Array;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _columns; j++)
                {
                    x[i][j] = _data[i][j] + b[i, j];
                }
            }
            return xx;
        }

        /// <summary>Matrix subtraction.</summary>
        public IMatrix Subtraction(IMatrix b)
        {
            if ((_rows != b.Rows) || (_columns != b.Columns))
            {
                throw new ArgumentException("Matrix dimension do not match.");
            }
            Matrix xx = new Matrix(_rows, _columns);
            double[][] x = xx.Array;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _columns; j++)
                {
                    x[i][j] = _data[i][j] - b[i, j];
                }
            }
            return xx;
        }

        /// <summary>Matrix-scalar multiplication.</summary>
        public IMatrix Multiply(double s)
        {
            Matrix xx = new Matrix(_rows, _columns);
            double[][] x = xx.Array;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _columns; j++)
                {
                    x[i][j] = _data[i][j]*s;
                }
            }

            return xx;
        }


        public IMatrix Sum()
        {
            if (_rows == 1)
            {
                Matrix m1 = new Matrix(1, 1);
                double sum = 0;
                for (int i = 0; i < _columns; i++)
                    sum += this[0, i];
                m1[0, 0] = sum;
                return m1;
            }

            Matrix m = new Matrix(1, _columns);

            for (int c = 0; c < _columns; c++)
            {
                double sum = 0;
                for (int r = 0; r < _rows; r++)
                    sum += this[r, c];
                m[0, c] = sum;
            }
            return m;
        }

        public IMatrix DotMultiply(IMatrix b)
        {
            if (b.Rows != Rows)
            {
                throw new ArgumentException("Matrix dimensions are not valid.");
            }
            if (b.Columns != Columns)
            {
                throw new ArgumentException("Matrix dimensions are not valid.");
            }

            Matrix xx = new Matrix(_rows, _columns);

            for (int i = 0; i < _rows; i++)
                for (int j = 0; j < _columns; j++)
                    xx[i, j] = this[i, j]*b[i, j];

            return xx;
        }

        /// <summary>Matrix-matrix multiplication.</summary>
        public IMatrix Multiply(IMatrix b)
        {
            if (b.Rows != _columns)
            {
                throw new ArgumentException("Matrix dimensions are not valid.");
            }

            int lColumns = b.Columns;
            Matrix xx = new Matrix(_rows, lColumns);
            double[][] x = xx.Array;

            int size = _columns;
            double[] column = new double[size];
            for (int j = 0; j < lColumns; j++)
            {
                for (int k = 0; k < size; k++)
                {
                    column[k] = b[k, j];
                }
                for (int i = 0; i < _rows; i++)
                {
                    double[] row = _data[i];
                    double s = 0;
                    for (int k = 0; k < size; k++)
                    {
                        s += row[k]*column[k];
                    }
                    x[i][j] = s;
                }
            }

            return xx;
        }

        /// <summary>Returns the LHS solution vetor if the matrix is square or the least squares solution otherwise.</summary>
        public IMatrix Solve(IMatrix rhs)
        {
            return (_rows == _columns) ? GetLuDecomposition().Solve(rhs) : GetQrDecomposition().Solve(rhs);
        }

        /// <summary> returns this / rhs </summary>
        public IMatrix Divide(IMatrix rhs)
        {
            IMatrix bt = rhs.Transpose();
            IMatrix at = Transpose();
            return bt.Solve(at).Transpose();
        }

        /// <summary>Inverse of the matrix if matrix is square, pseudoinverse otherwise.</summary>
        public IMatrix Inverse
        {
            get { return Solve(Diagonal(_rows, _rows, 1.0)); }
        }

        /// <summary>Determinant if matrix is square.</summary>
        public double Determinant
        {
            get { return GetLuDecomposition().Determinant; }
        }

        /// <summary>Returns the trace of the matrix.</summary>
        /// <returns>Sum of the diagonal elements.</returns>
        public double Trace
        {
            get
            {
                double trace = 0;
                for (int i = 0; i < Math.Min(_rows, _columns); i++)
                {
                    trace += _data[i][i];
                }
                return trace;
            }
        }

        /// <summary>Returns the cholesky decomposition for this matrix.</summary>
        public ICholeskyDecomposition GetCholeskyDecomposition()
        {
            return new CholeskyDecomposition(this);
        }

        /// <summary>Returns the LU decomposition for this matrix.</summary>
        public ILuDecomposition GetLuDecomposition()
        {
            return new LuDecomposition(this);
        }

        /// <summary>Returns the singular value decomposition for this matrix.</summary>
        public ISingularValueDecomposition GetSingularValueDecomposition()
        {
            return new SingularValueDecomposition(this);
        }

        /// <summary>Returns the QR decomposition for this matrix.</summary>
        public IQrDecomposition GetQrDecomposition()
        {
            return new QrDecomposition(this);
        }

        /// <summary>Returns the eigenvalue decomposition for this matrix.</summary>
        public IEigenvalueDecomposition GetEigenvalueDecomposition()
        {
            return new EigenvalueDecomposition(this);
        }

        /// <summary>Unary minus.</summary>
        public IMatrix UnaryMinus()
        {
            int lRows = _rows;
            int lColumns = _columns;
            Matrix xx = new Matrix(lRows, lColumns);
            double[][] x = xx.Array;
            for (int i = 0; i < lRows; i++)
                for (int j = 0; j < lColumns; j++)
                    x[i][j] = -_data[i][j];
            return xx;
        }

        /// <summary>Returns a matrix filled with random values.</summary>
        public static IMatrix Random(int rows, int columns)
        {
            Matrix xx = new Matrix(rows, columns);
            double[][] x = xx.Array;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    x[i][j] = MathHelper.Random();
                }
            }
            return xx;
        }

        /// <summary>Returns a diagonal matrix of the given size.</summary>
        public static IMatrix Diagonal(int rows, int columns, double value)
        {
            Matrix xx = new Matrix(rows, columns);
            double[][] x = xx.Array;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    x[i][j] = ((i == j) ? value : 0.0);
                }
            }
            return xx;
        }

        /// <summary>Returns the matrix in a textual form.</summary>
        public override string ToString()
        {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _columns; j++)
                    builder.Append(_data[i][j] + " ");

                builder.Append(Environment.NewLine);
            }

            return builder.ToString();
        }

        private class CholeskyDecomposition : ICholeskyDecomposition
        {
            private readonly Matrix _l;
            private readonly bool _isPositiveDefinite;
            private readonly bool _isSymmetric;

            public CholeskyDecomposition(Matrix pA)
            {
                if (!pA.IsSquare)
                {
                    throw new ArgumentNullException("Matrix is not square.");
                }

                int dimension = pA.Rows;
                _l = new Matrix(dimension, dimension);

                double[][] a = pA.Array;
                double[][] l = _l.Array;

                _isPositiveDefinite = true;
                _isSymmetric = true;

                for (int j = 0; j < dimension; j++)
                {
                    double[] lrowj = l[j];
                    double d = 0.0;
                    for (int k = 0; k < j; k++)
                    {
                        double[] lrowk = l[k];
                        double s = 0.0;
                        for (int i = 0; i < k; i++)
                        {
                            s += lrowk[i]*lrowj[i];
                        }
                        lrowj[k] = s = (a[j][k] - s)/l[k][k];
                        d = d + s*s;
                        _isSymmetric = _isSymmetric & (a[k][j] == a[j][k]);
                    }

                    d = a[j][j] - d;
                    _isPositiveDefinite = _isPositiveDefinite & (d > 0.0);
                    l[j][j] = Math.Sqrt(Math.Max(d, 0.0));
                    for (int k = j + 1; k < dimension; k++)
                        l[j][k] = 0.0;
                }
            }

            public bool IsSymmetric
            {
                get { return _isSymmetric; }
            }

            public bool IsPositiveDefinite
            {
                get { return _isPositiveDefinite; }
            }

            public IMatrix LeftTriangularFactor
            {
                get { return _l; }
            }

            public IMatrix Solve(IMatrix rhs)
            {
                if (rhs.Rows != _l.Rows)
                {
                    throw new ArgumentException("Matrix dimensions do not match.");
                }
                if (!_isSymmetric)
                {
                    throw new InvalidOperationException("Matrix is not symmetric.");
                }
                if (!_isPositiveDefinite)
                {
                    throw new InvalidOperationException("Matrix is not positive definite.");
                }

                int dimension = _l.Rows;
                int count = rhs.Columns;

                IMatrix b = rhs.Clone();
                double[][] l = _l.Array;

                // Solve L*Y = B;
                for (int k = 0; k < _l.Rows; k++)
                {
                    for (int i = k + 1; i < dimension; i++)
                    {
                        for (int j = 0; j < count; j++)
                        {
                            b[i, j] -= b[k, j]*l[i][k];
                        }
                    }

                    for (int j = 0; j < count; j++)
                    {
                        b[k, j] /= l[k][k];
                    }
                }

                // Solve L'*X = Y;
                for (int k = dimension - 1; k >= 0; k--)
                {
                    for (int j = 0; j < count; j++)
                    {
                        b[k, j] /= l[k][k];
                    }

                    for (int i = 0; i < k; i++)
                    {
                        for (int j = 0; j < count; j++)
                        {
                            b[i, j] -= b[k, j]*l[k][i];
                        }
                    }
                }

                return b;
            }
        }

        private class EigenvalueDecomposition : IEigenvalueDecomposition
        {
            private readonly Matrix _h; // storage of nonsymmetric Hessenberg form.
            private readonly Matrix _v; // storage of eigenvectors.
            private readonly double[] _d; // storage of eigenvalues.
            private readonly double[] _e; // storage of eigenvalues.
            private readonly bool _isSymmetric;
            private readonly int _n; // matrix dimension
            private readonly double[] _ort; // storage for nonsymmetric algorithm.
            private double _cdivi;
            private double _cdivr;

            public EigenvalueDecomposition(Matrix a)
            {
                if (a.Rows != a.Columns) throw new ArgumentException("Matrix is not a square matrix.");

                _n = a.Columns;
                _v = new Matrix(_n, _n);
                _d = new double[_n];
                _e = new double[_n];

                // Check for symmetry.
                _isSymmetric = a.IsSymmetric;

                if (_isSymmetric)
                {
                    for (int i = 0; i < _n; i++)
                        for (int j = 0; j < _n; j++)
                            _v[i, j] = a[i, j];

                    // Tridiagonalize.
                    Tred2();
                    // Diagonalize.
                    Tql2();
                }
                else
                {
                    _h = new Matrix(_n, _n);
                    _ort = new double[_n];

                    for (int j = 0; j < _n; j++)
                        for (int i = 0; i < _n; i++)
                            _h[i, j] = a[i, j];

                    // Reduce to Hessenberg form.
                    Orthes();

                    // Reduce Hessenberg to real Schur form.
                    Hqr2();
                }
            }

            public double[] RealEigenvalues
            {
                get { return _d; }
            }

            public double[] ImaginaryEigenvalues
            {
                get { return _e; }
            }

            public IMatrix EigenvectorMatrix
            {
                get { return _v; }
            }

            public IMatrix DiagonalMatrix
            {
                get
                {
                    Matrix xx = new Matrix(_n, _n);
                    double[][] x = xx.Array;

                    for (int i = 0; i < _n; i++)
                    {
                        for (int j = 0; j < _n; j++)
                            x[i][j] = 0.0;

                        x[i][i] = _d[i];
                        if (_e[i] > 0)
                        {
                            x[i][i + 1] = _e[i];
                        }
                        else if (_e[i] < 0)
                        {
                            x[i][i - 1] = _e[i];
                        }
                    }

                    return xx;
                }
            }

            private void Tred2()
            {
                // Symmetric Householder reduction to tridiagonal form.
                // This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
                // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
                for (int j = 0; j < _n; j++)
                    _d[j] = _v[_n - 1, j];

                // Householder reduction to tridiagonal form.
                for (int i = _n - 1; i > 0; i--)
                {
                    // Scale to avoid under/overflow.
                    double scale = 0.0;
                    double h = 0.0;
                    for (int k = 0; k < i; k++)
                        scale = scale + Math.Abs(_d[k]);

                    if (scale == 0.0)
                    {
                        _e[i] = _d[i - 1];
                        for (int j = 0; j < i; j++)
                        {
                            _d[j] = _v[i - 1, j];
                            _v[i, j] = 0.0;
                            _v[j, i] = 0.0;
                        }
                    }
                    else
                    {
                        // Generate Householder vector.
                        for (int k = 0; k < i; k++)
                        {
                            _d[k] /= scale;
                            h += _d[k]*_d[k];
                        }

                        double f = _d[i - 1];
                        double g = Math.Sqrt(h);
                        if (f > 0) g = -g;

                        _e[i] = scale*g;
                        h = h - f*g;
                        _d[i - 1] = f - g;
                        for (int j = 0; j < i; j++)
                            _e[j] = 0.0;

                        // Apply similarity transformation to remaining columns.
                        for (int j = 0; j < i; j++)
                        {
                            f = _d[j];
                            _v[j, i] = f;
                            g = _e[j] + _v[j, j]*f;
                            for (int k = j + 1; k <= i - 1; k++)
                            {
                                g += _v[k, j]*_d[k];
                                _e[k] += _v[k, j]*f;
                            }
                            _e[j] = g;
                        }

                        f = 0.0;
                        for (int j = 0; j < i; j++)
                        {
                            _e[j] /= h;
                            f += _e[j]*_d[j];
                        }

                        double hh = f/(h + h);
                        for (int j = 0; j < i; j++)
                            _e[j] -= hh*_d[j];

                        for (int j = 0; j < i; j++)
                        {
                            f = _d[j];
                            g = _e[j];
                            for (int k = j; k <= i - 1; k++)
                                _v[k, j] -= (f*_e[k] + g*_d[k]);

                            _d[j] = _v[i - 1, j];
                            _v[i, j] = 0.0;
                        }
                    }
                    _d[i] = h;
                }

                // Accumulate transformations.
                for (int i = 0; i < _n - 1; i++)
                {
                    _v[_n - 1, i] = _v[i, i];
                    _v[i, i] = 1.0;
                    double h = _d[i + 1];
                    if (h != 0.0)
                    {
                        for (int k = 0; k <= i; k++)
                            _d[k] = _v[k, i + 1]/h;

                        for (int j = 0; j <= i; j++)
                        {
                            double g = 0.0;
                            for (int k = 0; k <= i; k++)
                                g += _v[k, i + 1]*_v[k, j];
                            for (int k = 0; k <= i; k++)
                                _v[k, j] -= g*_d[k];
                        }
                    }

                    for (int k = 0; k <= i; k++)
                        _v[k, i + 1] = 0.0;
                }

                for (int j = 0; j < _n; j++)
                {
                    _d[j] = _v[_n - 1, j];
                    _v[_n - 1, j] = 0.0;
                }

                _v[_n - 1, _n - 1] = 1.0;
                _e[0] = 0.0;
            }

            private void Tql2()
            {
                // Symmetric tridiagonal QL algorithm.
                // This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson, 
                // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
                for (int i = 1; i < _n; i++)
                    _e[i - 1] = _e[i];

                _e[_n - 1] = 0.0;

                double f = 0.0;
                double tst1 = 0.0;
                double eps = Math.Pow(2.0, -52.0);

                for (int l = 0; l < _n; l++)
                {
                    // Find small subdiagonal element.
                    tst1 = Math.Max(tst1, Math.Abs(_d[l]) + Math.Abs(_e[l]));
                    int m = l;
                    while (m < _n)
                    {
                        if (Math.Abs(_e[m]) <= eps*tst1)
                            break;
                        m++;
                    }

                    // If m == l, d[l] is an eigenvalue, otherwise, iterate.
                    if (m > l)
                    {
                        int iter = 0;
                        do
                        {
                            iter = iter + 1; // (Could check iteration count here.)

                            // Compute implicit shift
                            double g = _d[l];
                            double p = (_d[l + 1] - g)/(2.0*_e[l]);
                            double r = MathHelper.Hypotenuse(p, 1.0);
                            if (p < 0) r = -r;

                            _d[l] = _e[l]/(p + r);
                            _d[l + 1] = _e[l]*(p + r);
                            double dl1 = _d[l + 1];
                            double h = g - _d[l];
                            for (int i = l + 2; i < _n; i++)
                                _d[i] -= h;
                            f = f + h;

                            // Implicit QL transformation.
                            p = _d[m];
                            double c = 1.0;
                            double c2 = c;
                            double c3 = c;
                            double el1 = _e[l + 1];
                            double s = 0.0;
                            double s2 = 0.0;
                            for (int i = m - 1; i >= l; i--)
                            {
                                c3 = c2;
                                c2 = c;
                                s2 = s;
                                g = c*_e[i];
                                h = c*p;
                                r = MathHelper.Hypotenuse(p, _e[i]);
                                _e[i + 1] = s*r;
                                s = _e[i]/r;
                                c = p/r;
                                p = c*_d[i] - s*g;
                                _d[i + 1] = h + s*(c*g + s*_d[i]);

                                // Accumulate transformation.
                                for (int k = 0; k < _n; k++)
                                {
                                    h = _v[k, i + 1];
                                    _v[k, i + 1] = s*_v[k, i] + c*h;
                                    _v[k, i] = c*_v[k, i] - s*h;
                                }
                            }

                            p = -s*s2*c3*el1*_e[l]/dl1;
                            _e[l] = s*p;
                            _d[l] = c*p;

                            // Check for convergence.
                        } while (Math.Abs(_e[l]) > eps*tst1);
                    }
                    _d[l] = _d[l] + f;
                    _e[l] = 0.0;
                }

                // Sort eigenvalues and corresponding vectors.
                for (int i = 0; i < _n - 1; i++)
                {
                    int k = i;
                    double p = _d[i];
                    for (int j = i + 1; j < _n; j++)
                    {
                        if (_d[j] < p)
                        {
                            k = j;
                            p = _d[j];
                        }
                    }

                    if (k != i)
                    {
                        _d[k] = _d[i];
                        _d[i] = p;
                        for (int j = 0; j < _n; j++)
                        {
                            p = _v[j, i];
                            _v[j, i] = _v[j, k];
                            _v[j, k] = p;
                        }
                    }
                }
            }

            private void Orthes()
            {
                // Nonsymmetric reduction to Hessenberg form.
                // This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
                // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
                const int low = 0;
                int high = _n - 1;

                for (int m = low + 1; m <= high - 1; m++)
                {
                    // Scale column.

                    double scale = 0.0;
                    for (int i = m; i <= high; i++)
                        scale = scale + Math.Abs(_h[i, m - 1]);

                    if (scale != 0.0)
                    {
                        // Compute Householder transformation.
                        double h = 0.0;
                        for (int i = high; i >= m; i--)
                        {
                            _ort[i] = _h[i, m - 1]/scale;
                            h += _ort[i]*_ort[i];
                        }

                        double g = Math.Sqrt(h);
                        if (_ort[m] > 0) g = -g;

                        h = h - _ort[m]*g;
                        _ort[m] = _ort[m] - g;

                        // Apply Householder similarity transformation
                        // H = (I - u * u' / h) * H * (I - u * u') / h)
                        for (int j = m; j < _n; j++)
                        {
                            double f = 0.0;
                            for (int i = high; i >= m; i--)
                                f += _ort[i]*_h[i, j];

                            f = f/h;
                            for (int i = m; i <= high; i++)
                                _h[i, j] -= f*_ort[i];
                        }

                        for (int i = 0; i <= high; i++)
                        {
                            double f = 0.0;
                            for (int j = high; j >= m; j--)
                                f += _ort[j]*_h[i, j];

                            f = f/h;
                            for (int j = m; j <= high; j++)
                                _h[i, j] -= f*_ort[j];
                        }

                        _ort[m] = scale*_ort[m];
                        _h[m, m - 1] = scale*g;
                    }
                }

                // Accumulate transformations (Algol's ortran).
                for (int i = 0; i < _n; i++)
                    for (int j = 0; j < _n; j++)
                        _v[i, j] = (i == j ? 1.0 : 0.0);

                for (int m = high - 1; m >= low + 1; m--)
                {
                    if (_h[m, m - 1] != 0.0)
                    {
                        for (int i = m + 1; i <= high; i++)
                            _ort[i] = _h[i, m - 1];

                        for (int j = m; j <= high; j++)
                        {
                            double g = 0.0;
                            for (int i = m; i <= high; i++)
                                g += _ort[i]*_v[i, j];

                            // Double division avoids possible underflow.
                            g = (g/_ort[m])/_h[m, m - 1];
                            for (int i = m; i <= high; i++)
                                _v[i, j] += g*_ort[i];
                        }
                    }
                }
            }

            private void Cdiv(double xr, double xi, double yr, double yi)
            {
                // Complex scalar division.
                double r, dd;
                if (Math.Abs(yr) > Math.Abs(yi))
                {
                    r = yi/yr;
                    dd = yr + r*yi;
                    _cdivr = (xr + r*xi)/dd;
                    _cdivi = (xi - r*xr)/dd;
                }
                else
                {
                    r = yr/yi;
                    dd = yi + r*yr;
                    _cdivr = (r*xr + xi)/dd;
                    _cdivi = (r*xi - xr)/dd;
                }
            }

            private void Hqr2()
            {
                // Nonsymmetric reduction from Hessenberg to real Schur form.   
                // This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
                // Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
                int nn = _n;
                int nNormalized = nn - 1;
                const int low = 0;
                int high = nn - 1;
                double eps = Math.Pow(2.0, -52.0);
                double exshift = 0.0;
                double p = 0, q = 0, r = 0, s = 0, z = 0;
                double w, x, y;

                // Store roots isolated by balanc and compute matrix norm
                double norm = 0.0;
                for (int i = 0; i < nn; i++)
                {
                    if (i < low | i > high)
                    {
                        _d[i] = _h[i, i];
                        _e[i] = 0.0;
                    }

                    for (int j = Math.Max(i - 1, 0); j < nn; j++)
                        norm = norm + Math.Abs(_h[i, j]);
                }

                // Outer loop over eigenvalue index
                int iter = 0;
                while (nNormalized >= low)
                {
                    // Look for single small sub-diagonal element
                    int l = nNormalized;
                    while (l > low)
                    {
                        s = Math.Abs(_h[l - 1, l - 1]) + Math.Abs(_h[l, l]);
                        if (s == 0.0) s = norm;
                        if (Math.Abs(_h[l, l - 1]) < eps*s)
                            break;

                        l--;
                    }

                    // Check for convergence
                    if (l == nNormalized)
                    {
                        // One root found
                        _h[nNormalized, nNormalized] = _h[nNormalized, nNormalized] + exshift;
                        _d[nNormalized] = _h[nNormalized, nNormalized];
                        _e[nNormalized] = 0.0;
                        nNormalized--;
                        iter = 0;
                    }
                    else if (l == nNormalized - 1)
                    {
                        // Two roots found
                        w = _h[nNormalized, nNormalized - 1]*_h[nNormalized - 1, nNormalized];
                        p = (_h[nNormalized - 1, nNormalized - 1] - _h[nNormalized, nNormalized])/2.0;
                        q = p*p + w;
                        z = Math.Sqrt(Math.Abs(q));
                        _h[nNormalized, nNormalized] = _h[nNormalized, nNormalized] + exshift;
                        _h[nNormalized - 1, nNormalized - 1] = _h[nNormalized - 1, nNormalized - 1] + exshift;
                        x = _h[nNormalized, nNormalized];

                        if (q >= 0)
                        {
                            // Real pair
                            z = (p >= 0) ? (p + z) : (p - z);
                            _d[nNormalized - 1] = x + z;
                            _d[nNormalized] = _d[nNormalized - 1];
                            if (z != 0.0)
                                _d[nNormalized] = x - w/z;
                            _e[nNormalized - 1] = 0.0;
                            _e[nNormalized] = 0.0;
                            x = _h[nNormalized, nNormalized - 1];
                            s = Math.Abs(x) + Math.Abs(z);
                            p = x/s;
                            q = z/s;
                            r = Math.Sqrt(p*p + q*q);
                            p = p/r;
                            q = q/r;

                            // Row modification
                            for (int j = nNormalized - 1; j < nn; j++)
                            {
                                z = _h[nNormalized - 1, j];
                                _h[nNormalized - 1, j] = q*z + p*_h[nNormalized, j];
                                _h[nNormalized, j] = q*_h[nNormalized, j] - p*z;
                            }

                            // Column modification
                            for (int i = 0; i <= nNormalized; i++)
                            {
                                z = _h[i, nNormalized - 1];
                                _h[i, nNormalized - 1] = q*z + p*_h[i, nNormalized];
                                _h[i, nNormalized] = q*_h[i, nNormalized] - p*z;
                            }

                            // Accumulate transformations
                            for (int i = low; i <= high; i++)
                            {
                                z = _v[i, nNormalized - 1];
                                _v[i, nNormalized - 1] = q*z + p*_v[i, nNormalized];
                                _v[i, nNormalized] = q*_v[i, nNormalized] - p*z;
                            }
                        }
                        else
                        {
                            // Complex pair
                            _d[nNormalized - 1] = x + p;
                            _d[nNormalized] = x + p;
                            _e[nNormalized - 1] = z;
                            _e[nNormalized] = -z;
                        }

                        nNormalized = nNormalized - 2;
                        iter = 0;
                    }
                    else
                    {
                        // No convergence yet    

                        // Form shift
                        x = _h[nNormalized, nNormalized];
                        y = 0.0;
                        w = 0.0;
                        if (l < nNormalized)
                        {
                            y = _h[nNormalized - 1, nNormalized - 1];
                            w = _h[nNormalized, nNormalized - 1]*_h[nNormalized - 1, nNormalized];
                        }

                        // Wilkinson's original ad hoc shift
                        if (iter == 10)
                        {
                            exshift += x;
                            for (int i = low; i <= nNormalized; i++)
                                _h[i, i] -= x;

                            s = Math.Abs(_h[nNormalized, nNormalized - 1]) + Math.Abs(_h[nNormalized - 1, nNormalized - 2]);
                            x = y = 0.75*s;
                            w = -0.4375*s*s;
                        }

                        // MATLAB's new ad hoc shift
                        if (iter == 30)
                        {
                            s = (y - x)/2.0;
                            s = s*s + w;
                            if (s > 0)
                            {
                                s = Math.Sqrt(s);
                                if (y < x) s = -s;
                                s = x - w/((y - x)/2.0 + s);
                                for (int i = low; i <= nNormalized; i++)
                                    _h[i, i] -= s;
                                exshift += s;
                                x = y = w = 0.964;
                            }
                        }

                        iter = iter + 1;

                        // Look for two consecutive small sub-diagonal elements
                        int m = nNormalized - 2;
                        while (m >= l)
                        {
                            z = _h[m, m];
                            r = x - z;
                            s = y - z;
                            p = (r*s - w)/_h[m + 1, m] + _h[m, m + 1];
                            q = _h[m + 1, m + 1] - z - r - s;
                            r = _h[m + 2, m + 1];
                            s = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                            p = p/s;
                            q = q/s;
                            r = r/s;
                            if (m == l)
                                break;
                            if (Math.Abs(_h[m, m - 1])*(Math.Abs(q) + Math.Abs(r)) <
                                eps*(Math.Abs(p)*(Math.Abs(_h[m - 1, m - 1]) + Math.Abs(z) + Math.Abs(_h[m + 1, m + 1]))))
                                break;
                            m--;
                        }

                        for (int i = m + 2; i <= nNormalized; i++)
                        {
                            _h[i, i - 2] = 0.0;
                            if (i > m + 2)
                                _h[i, i - 3] = 0.0;
                        }

                        // Double QR step involving rows l:n and columns m:n
                        for (int k = m; k <= nNormalized - 1; k++)
                        {
                            bool notlast = (k != nNormalized - 1);
                            if (k != m)
                            {
                                p = _h[k, k - 1];
                                q = _h[k + 1, k - 1];
                                r = (notlast ? _h[k + 2, k - 1] : 0.0);
                                x = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                                if (x != 0.0)
                                {
                                    p = p/x;
                                    q = q/x;
                                    r = r/x;
                                }
                            }

                            if (x == 0.0) break;

                            s = Math.Sqrt(p*p + q*q + r*r);
                            if (p < 0) s = -s;

                            if (s != 0)
                            {
                                if (k != m)
                                    _h[k, k - 1] = -s*x;
                                else if (l != m)
                                    _h[k, k - 1] = -_h[k, k - 1];

                                p = p + s;
                                x = p/s;
                                y = q/s;
                                z = r/s;
                                q = q/p;
                                r = r/p;

                                // Row modification
                                for (int j = k; j < nn; j++)
                                {
                                    p = _h[k, j] + q*_h[k + 1, j];
                                    if (notlast)
                                    {
                                        p = p + r*_h[k + 2, j];
                                        _h[k + 2, j] = _h[k + 2, j] - p*z;
                                    }

                                    _h[k, j] = _h[k, j] - p*x;
                                    _h[k + 1, j] = _h[k + 1, j] - p*y;
                                }

                                // Column modification
                                for (int i = 0; i <= Math.Min(nNormalized, k + 3); i++)
                                {
                                    p = x*_h[i, k] + y*_h[i, k + 1];
                                    if (notlast)
                                    {
                                        p = p + z*_h[i, k + 2];
                                        _h[i, k + 2] = _h[i, k + 2] - p*r;
                                    }

                                    _h[i, k] = _h[i, k] - p;
                                    _h[i, k + 1] = _h[i, k + 1] - p*q;
                                }

                                // Accumulate transformations
                                for (int i = low; i <= high; i++)
                                {
                                    p = x*_v[i, k] + y*_v[i, k + 1];
                                    if (notlast)
                                    {
                                        p = p + z*_v[i, k + 2];
                                        _v[i, k + 2] = _v[i, k + 2] - p*r;
                                    }

                                    _v[i, k] = _v[i, k] - p;
                                    _v[i, k + 1] = _v[i, k + 1] - p*q;
                                }
                            }
                        }
                    }
                }

                // Backsubstitute to find vectors of upper triangular form
                if (norm == 0.0) return;

                for (nNormalized = nn - 1; nNormalized >= 0; nNormalized--)
                {
                    p = _d[nNormalized];
                    q = _e[nNormalized];

                    // Real vector
                    double t;
                    if (q == 0)
                    {
                        int l = nNormalized;
                        _h[nNormalized, nNormalized] = 1.0;
                        for (int i = nNormalized - 1; i >= 0; i--)
                        {
                            w = _h[i, i] - p;
                            r = 0.0;
                            for (int j = l; j <= nNormalized; j++)
                                r = r + _h[i, j]*_h[j, nNormalized];

                            if (_e[i] < 0.0)
                            {
                                z = w;
                                s = r;
                            }
                            else
                            {
                                l = i;
                                if (_e[i] == 0.0)
                                {
                                    _h[i, nNormalized] = (w != 0.0) ? (-r/w) : (-r/(eps*norm));
                                }
                                else
                                {
                                    // Solve real equations
                                    x = _h[i, i + 1];
                                    y = _h[i + 1, i];
                                    q = (_d[i] - p)*(_d[i] - p) + _e[i]*_e[i];
                                    t = (x*s - z*r)/q;
                                    _h[i, nNormalized] = t;
                                    _h[i + 1, nNormalized] = (Math.Abs(x) > Math.Abs(z)) ? ((-r - w*t)/x) : ((-s - y*t)/z);
                                }

                                // Overflow control
                                t = Math.Abs(_h[i, nNormalized]);
                                if ((eps*t)*t > 1)
                                    for (int j = i; j <= nNormalized; j++)
                                        _h[j, nNormalized] = _h[j, nNormalized]/t;
                            }
                        }
                    }
                    else if (q < 0)
                    {
                        // Complex vector
                        int l = nNormalized - 1;

                        // Last vector component imaginary so matrix is triangular
                        if (Math.Abs(_h[nNormalized, nNormalized - 1]) > Math.Abs(_h[nNormalized - 1, nNormalized]))
                        {
                            _h[nNormalized - 1, nNormalized - 1] = q/_h[nNormalized, nNormalized - 1];
                            _h[nNormalized - 1, nNormalized] = -(_h[nNormalized, nNormalized] - p)/_h[nNormalized, nNormalized - 1];
                        }
                        else
                        {
                            Cdiv(0.0, -_h[nNormalized - 1, nNormalized], _h[nNormalized - 1, nNormalized - 1] - p, q);
                            _h[nNormalized - 1, nNormalized - 1] = _cdivr;
                            _h[nNormalized - 1, nNormalized] = _cdivi;
                        }

                        _h[nNormalized, nNormalized - 1] = 0.0;
                        _h[nNormalized, nNormalized] = 1.0;
                        for (int i = nNormalized - 2; i >= 0; i--)
                        {
                            double ra = 0.0;
                            double sa = 0.0;
                            for (int j = l; j <= nNormalized; j++)
                            {
                                ra = ra + _h[i, j]*_h[j, nNormalized - 1];
                                sa = sa + _h[i, j]*_h[j, nNormalized];
                            }

                            w = _h[i, i] - p;

                            if (_e[i] < 0.0)
                            {
                                z = w;
                                r = ra;
                                s = sa;
                            }
                            else
                            {
                                l = i;
                                if (_e[i] == 0)
                                {
                                    Cdiv(-ra, -sa, w, q);
                                    _h[i, nNormalized - 1] = _cdivr;
                                    _h[i, nNormalized] = _cdivi;
                                }
                                else
                                {
                                    // Solve complex equations
                                    x = _h[i, i + 1];
                                    y = _h[i + 1, i];
                                    double vr = (_d[i] - p)*(_d[i] - p) + _e[i]*_e[i] - q*q;
                                    double vi = (_d[i] - p)*2.0*q;
                                    if (vr == 0.0 & vi == 0.0)
                                        vr = eps*norm*
                                             (Math.Abs(w) + Math.Abs(q) + Math.Abs(x) + Math.Abs(y) + Math.Abs(z));
                                    Cdiv(x*r - z*ra + q*sa, x*s - z*sa - q*ra, vr, vi);
                                    _h[i, nNormalized - 1] = _cdivr;
                                    _h[i, nNormalized] = _cdivi;
                                    if (Math.Abs(x) > (Math.Abs(z) + Math.Abs(q)))
                                    {
                                        _h[i + 1, nNormalized - 1] = (-ra - w*_h[i, nNormalized - 1] + q*_h[i, nNormalized])/x;
                                        _h[i + 1, nNormalized] = (-sa - w*_h[i, nNormalized] - q*_h[i, nNormalized - 1])/x;
                                    }
                                    else
                                    {
                                        Cdiv(-r - y*_h[i, nNormalized - 1], -s - y*_h[i, nNormalized], z, q);
                                        _h[i + 1, nNormalized - 1] = _cdivr;
                                        _h[i + 1, nNormalized] = _cdivi;
                                    }
                                }

                                // Overflow control
                                t = Math.Max(Math.Abs(_h[i, nNormalized - 1]), Math.Abs(_h[i, nNormalized]));
                                if ((eps*t)*t > 1)
                                    for (int j = i; j <= nNormalized; j++)
                                    {
                                        _h[j, nNormalized - 1] = _h[j, nNormalized - 1]/t;
                                        _h[j, nNormalized] = _h[j, nNormalized]/t;
                                    }
                            }
                        }
                    }
                }

                // Vectors of isolated roots
                for (int i = 0; i < nn; i++)
                    if (i < low | i > high)
                        for (int j = i; j < nn; j++)
                            _v[i, j] = _h[i, j];

                // Back transformation to get eigenvectors of original matrix
                for (int j = nn - 1; j >= low; j--)
                    for (int i = low; i <= high; i++)
                    {
                        z = 0.0;
                        for (int k = low; k <= Math.Min(j, high); k++)
                            z = z + _v[i, k]*_h[k, j];
                        _v[i, j] = z;
                    }
            }
        }

        private class LuDecomposition : ILuDecomposition
        {
            private readonly Matrix _lu;
            private readonly int _pivotSign;
            private readonly int[] _pivotVector;

            public LuDecomposition(Matrix a)
            {
                _lu = (Matrix) a.Clone();
                double[][] lu = _lu.Array;
                int rows = a.Rows;
                int columns = a.Columns;
                _pivotVector = new int[rows];
                for (int i = 0; i < rows; i++)
                    _pivotVector[i] = i;
                _pivotSign = 1;
                double[] luColj = new double[rows];

                // Outer loop.
                for (int j = 0; j < columns; j++)
                {
                    // Make a copy of the j-th column to localize references.
                    for (int i = 0; i < rows; i++)
                        luColj[i] = lu[i][j];

                    // Apply previous transformations.
                    for (int i = 0; i < rows; i++)
                    {
                        double[] luRowi = lu[i];

                        // Most of the time is spent in the following dot product.
                        int kmax = Math.Min(i, j);
                        double s = 0.0;
                        for (int k = 0; k < kmax; k++)
                            s += luRowi[k]*luColj[k];
                        luRowi[j] = luColj[i] -= s;
                    }

                    // Find pivot and exchange if necessary.
                    int p = j;
                    for (int i = j + 1; i < rows; i++)
                        if (Math.Abs(luColj[i]) > Math.Abs(luColj[p]))
                            p = i;

                    if (p != j)
                    {
                        for (int k = 0; k < columns; k++)
                        {
                            double t = lu[p][k];
                            lu[p][k] = lu[j][k];
                            lu[j][k] = t;
                        }

                        int v = _pivotVector[p];
                        _pivotVector[p] = _pivotVector[j];
                        _pivotVector[j] = v;

                        _pivotSign = -_pivotSign;
                    }

                    // Compute multipliers.

                    if (j < rows & lu[j][j] != 0.0)
                    {
                        for (int i = j + 1; i < rows; i++)
                        {
                            lu[i][j] /= lu[j][j];
                        }
                    }
                }
            }

            public bool IsNonSingular
            {
                get
                {
                    for (int j = 0; j < _lu.Columns; j++)
                        if (_lu[j, j] == 0)
                            return false;
                    return true;
                }
            }

            public double Determinant
            {
                get
                {
                    if (_lu.Rows != _lu.Columns) throw new ArgumentException("Matrix must be square.");
                    double determinant = _pivotSign;
                    for (int j = 0; j < _lu.Columns; j++)
                        determinant *= _lu[j, j];
                    return determinant;
                }
            }

            public IMatrix LowerTriangularFactor
            {
                get
                {
                    int rows = _lu.Rows;
                    int columns = _lu.Columns;
                    Matrix x = new Matrix(rows, columns);
                    for (int i = 0; i < rows; i++)
                        for (int j = 0; j < columns; j++)
                            if (i > j)
                                x[i, j] = _lu[i, j];
                            else if (i == j)
                                x[i, j] = 1.0;
                            else
                                x[i, j] = 0.0;
                    return x;
                }
            }

            public IMatrix UpperTriangularFactor
            {
                get
                {
                    int rows = _lu.Rows;
                    int columns = _lu.Columns;
                    Matrix x = new Matrix(rows, columns);
                    for (int i = 0; i < rows; i++)
                        for (int j = 0; j < columns; j++)
                            if (i <= j)
                                x[i, j] = _lu[i, j];
                            else
                                x[i, j] = 0.0;
                    return x;
                }
            }

            public double[] PivotPermutationVector
            {
                get
                {
                    int rows = _lu.Rows;
                    double[] p = new double[rows];
                    for (int i = 0; i < rows; i++)
                        p[i] = _pivotVector[i];
                    return p;
                }
            }

            public IMatrix Solve(IMatrix b)
            {
                if (b.Rows != _lu.Rows) throw new ArgumentException("Invalid matrix dimensions.");
                if (!IsNonSingular) throw new InvalidOperationException("Matrix is singular");

                // Copy right hand side with pivoting
                int count = b.Columns;
                IMatrix x = b.Submatrix(_pivotVector, 0, count - 1);

                //int rows = LU.Rows;
                int columns = _lu.Columns;
                double[][] lu = _lu.Array;

                // Solve L*Y = B(piv,:)
                for (int k = 0; k < columns; k++)
                {
                    for (int i = k + 1; i < columns; i++)
                    {
                        for (int j = 0; j < count; j++)
                        {
                            x[i, j] -= x[k, j]*lu[i][k];
                        }
                    }
                }

                // Solve U*X = Y;
                for (int k = columns - 1; k >= 0; k--)
                {
                    for (int j = 0; j < count; j++)
                    {
                        x[k, j] /= lu[k][k];
                    }

                    for (int i = 0; i < k; i++)
                    {
                        for (int j = 0; j < count; j++)
                        {
                            x[i, j] -= x[k, j]*lu[i][k];
                        }
                    }
                }

                return x;
            }
        }

        private static class MathHelper
        {
            private static readonly Random StaticRandom = new Random();

            public static double Random()
            {
                return StaticRandom.NextDouble();
            }

            public static double Hypotenuse(double a, double b)
            {
                if (Math.Abs(a) > Math.Abs(b))
                {
                    double r = b/a;
                    return Math.Abs(a)*Math.Sqrt(1 + r*r);
                }

                if (b != 0)
                {
                    double r = a/b;
                    return Math.Abs(b)*Math.Sqrt(1 + r*r);
                }

                return 0.0;
            }
        }

        private class QrDecomposition : IQrDecomposition
        {
            private readonly Matrix _qr;
            private readonly double[] _rdiag;

            public QrDecomposition(Matrix a)
            {
                _qr = (Matrix) a.Clone();
                double[][] qr = _qr.Array;
                int m = a.Rows;
                int n = a.Columns;
                _rdiag = new double[n];

                for (int k = 0; k < n; k++)
                {
                    // Compute 2-norm of k-th column without under/overflow.
                    double nrm = 0;
                    for (int i = k; i < m; i++)
                        nrm = MathHelper.Hypotenuse(nrm, qr[i][k]);

                    if (nrm != 0.0)
                    {
                        // Form k-th Householder vector.
                        if (qr[k][k] < 0)
                            nrm = -nrm;
                        for (int i = k; i < m; i++)
                            qr[i][k] /= nrm;
                        qr[k][k] += 1.0;

                        // Apply transformation to remaining columns.
                        for (int j = k + 1; j < n; j++)
                        {
                            double s = 0.0;
                            for (int i = k; i < m; i++)
                                s += qr[i][k]*qr[i][j];
                            s = -s/qr[k][k];
                            for (int i = k; i < m; i++)
                                qr[i][j] += s*qr[i][k];
                        }
                    }
                    _rdiag[k] = -nrm;
                }
            }

            public IMatrix Solve(IMatrix rhs)
            {
                if (rhs.Rows != _qr.Rows) throw new ArgumentException("Matrix row dimensions must agree.");
                if (!IsFullRank) throw new InvalidOperationException("Matrix is rank deficient.");

                // Copy right hand side
                int count = rhs.Columns;
                IMatrix x = rhs.Clone();
                int m = _qr.Rows;
                int n = _qr.Columns;
                double[][] qr = _qr.Array;

                // Compute Y = transpose(Q)*B
                for (int k = 0; k < n; k++)
                {
                    for (int j = 0; j < count; j++)
                    {
                        double s = 0.0;
                        for (int i = k; i < m; i++)
                            s += qr[i][k]*x[i, j];
                        s = -s/qr[k][k];
                        for (int i = k; i < m; i++)
                            x[i, j] += s*qr[i][k];
                    }
                }

                // Solve R*X = Y;
                for (int k = n - 1; k >= 0; k--)
                {
                    for (int j = 0; j < count; j++)
                        x[k, j] /= _rdiag[k];

                    for (int i = 0; i < k; i++)
                        for (int j = 0; j < count; j++)
                            x[i, j] -= x[k, j]*qr[i][k];
                }

                return x.Submatrix(0, n - 1, 0, count - 1);
            }

            public bool IsFullRank
            {
                get
                {
                    int columns = _qr.Columns;
                    for (int j = 0; j < columns; j++)
                        if (_rdiag[j] == 0)
                            return false;
                    return true;
                }
            }

            public IMatrix UpperTriangularFactor
            {
                get
                {
                    int n = _qr.Columns;
                    Matrix xx = new Matrix(n, n);
                    double[][] x = xx.Array;
                    double[][] qr = _qr.Array;
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                            if (i < j)
                                x[i][j] = qr[i][j];
                            else if (i == j)
                                x[i][j] = _rdiag[i];
                            else
                                x[i][j] = 0.0;

                    return xx;
                }
            }

            public IMatrix OrthogonalFactor
            {
                get
                {
                    Matrix xx = new Matrix(_qr.Rows, _qr.Columns);
                    double[][] x = xx.Array;
                    double[][] qr = _qr.Array;
                    for (int k = _qr.Columns - 1; k >= 0; k--)
                    {
                        for (int i = 0; i < _qr.Rows; i++)
                            x[i][k] = 0.0;

                        x[k][k] = 1.0;
                        for (int j = k; j < _qr.Columns; j++)
                        {
                            if (qr[k][k] != 0)
                            {
                                double s = 0.0;
                                for (int i = k; i < _qr.Rows; i++)
                                    s += qr[i][k]*x[i][j];
                                s = -s/qr[k][k];
                                for (int i = k; i < _qr.Rows; i++)
                                    x[i][j] += s*qr[i][k];
                            }
                        }
                    }
                    return xx;
                }
            }
        }

        private class SingularValueDecomposition : ISingularValueDecomposition
        {
            private readonly Matrix _u;
            private readonly Matrix _v;
            private readonly int _m;
            private readonly int _n;
            private readonly double[] _s; // singular values

            public SingularValueDecomposition(Matrix pA)
            {
                Matrix copy = (Matrix) pA.Clone();
                double[][] a = copy.Array;
                _m = pA.Rows;
                _n = pA.Columns;
                int nu = Math.Min(_m, _n);
                _s = new double[Math.Min(_m + 1, _n)];
                _u = new Matrix(_m, nu);
                _v = new Matrix(_n, _n);
                double[][] u = _u.Array;
                double[][] v = _v.Array;
                double[] e = new double[_n];
                double[] work = new double[_m];
                const bool wantu = true;
                const bool wantv = true;

                // Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
                int nct = Math.Min(_m - 1, _n);
                int nrt = Math.Max(0, Math.Min(_n - 2, _m));
                for (int k = 0; k < Math.Max(nct, nrt); k++)
                {
                    if (k < nct)
                    {
                        // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
                        // Compute 2-norm of k-th column without under/overflow.
                        _s[k] = 0;
                        for (int i = k; i < _m; i++)
                            _s[k] = MathHelper.Hypotenuse(_s[k], a[i][k]);

                        if (_s[k] != 0.0)
                        {
                            if (a[k][k] < 0.0)
                                _s[k] = -_s[k];

                            for (int i = k; i < _m; i++)
                                a[i][k] /= _s[k];

                            a[k][k] += 1.0;
                        }
                        _s[k] = -_s[k];
                    }

                    for (int j = k + 1; j < _n; j++)
                    {
                        if ((k < nct) & (_s[k] != 0.0))
                        {
                            // Apply the transformation.
                            double t = 0;
                            for (int i = k; i < _m; i++)
                                t += a[i][k]*a[i][j];
                            t = -t/a[k][k];
                            for (int i = k; i < _m; i++)
                                a[i][j] += t*a[i][k];
                        }

                        // Place the k-th row of A into e for the subsequent calculation of the row transformation.
                        e[j] = a[k][j];
                    }

                    if (wantu & (k < nct))
                    {
                        // Place the transformation in U for subsequent back
                        // multiplication.
                        for (int i = k; i < _m; i++)
                            u[i][k] = a[i][k];
                    }

                    if (k < nrt)
                    {
                        // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
                        // Compute 2-norm without under/overflow.
                        e[k] = 0;
                        for (int i = k + 1; i < _n; i++)
                            e[k] = MathHelper.Hypotenuse(e[k], e[i]);

                        if (e[k] != 0.0)
                        {
                            if (e[k + 1] < 0.0)
                                e[k] = -e[k];

                            for (int i = k + 1; i < _n; i++)
                                e[i] /= e[k];

                            e[k + 1] += 1.0;
                        }

                        e[k] = -e[k];
                        if ((k + 1 < _m) & (e[k] != 0.0))
                        {
                            // Apply the transformation.
                            for (int i = k + 1; i < _m; i++)
                                work[i] = 0.0;

                            for (int j = k + 1; j < _n; j++)
                                for (int i = k + 1; i < _m; i++)
                                    work[i] += e[j]*a[i][j];

                            for (int j = k + 1; j < _n; j++)
                            {
                                double t = -e[j]/e[k + 1];
                                for (int i = k + 1; i < _m; i++)
                                    a[i][j] += t*work[i];
                            }
                        }

                        if (wantv)
                        {
                            // Place the transformation in V for subsequent back multiplication.
                            for (int i = k + 1; i < _n; i++)
                                v[i][k] = e[i];
                        }
                    }
                }

                // Set up the final bidiagonal matrix or order p.
                int p = Math.Min(_n, _m + 1);
                if (nct < _n) _s[nct] = a[nct][nct];
                if (_m < p) _s[p - 1] = 0.0;
                if (nrt + 1 < p) e[nrt] = a[nrt][p - 1];
                e[p - 1] = 0.0;

                // If required, generate U.
                if (wantu)
                {
                    for (int j = nct; j < nu; j++)
                    {
                        for (int i = 0; i < _m; i++)
                            u[i][j] = 0.0;
                        u[j][j] = 1.0;
                    }

                    for (int k = nct - 1; k >= 0; k--)
                    {
                        if (_s[k] != 0.0)
                        {
                            for (int j = k + 1; j < nu; j++)
                            {
                                double t = 0;
                                for (int i = k; i < _m; i++)
                                    t += u[i][k]*u[i][j];

                                t = -t/u[k][k];
                                for (int i = k; i < _m; i++)
                                    u[i][j] += t*u[i][k];
                            }

                            for (int i = k; i < _m; i++)
                                u[i][k] = -u[i][k];

                            u[k][k] = 1.0 + u[k][k];
                            for (int i = 0; i < k - 1; i++)
                                u[i][k] = 0.0;
                        }
                        else
                        {
                            for (int i = 0; i < _m; i++)
                                u[i][k] = 0.0;
                            u[k][k] = 1.0;
                        }
                    }
                }

                // If required, generate V.
                if (wantv)
                {
                    for (int k = _n - 1; k >= 0; k--)
                    {
                        if ((k < nrt) & (e[k] != 0.0))
                        {
                            for (int j = k + 1; j < nu; j++)
                            {
                                double t = 0;
                                for (int i = k + 1; i < _n; i++)
                                    t += v[i][k]*v[i][j];

                                t = -t/v[k + 1][k];
                                for (int i = k + 1; i < _n; i++)
                                    v[i][j] += t*v[i][k];
                            }
                        }

                        for (int i = 0; i < _n; i++)
                            v[i][k] = 0.0;
                        v[k][k] = 1.0;
                    }
                }

                // Main iteration loop for the singular values.
                int pp = p - 1;
                int iter = 0;
                double eps = Math.Pow(2.0, -52.0);
                while (p > 0)
                {
                    int k, kase;

                    // Here is where a test for too many iterations would go.
                    // This section of the program inspects for
                    // negligible elements in the s and e arrays.  On
                    // completion the variables kase and k are set as follows.
                    // kase = 1     if s(p) and e[k-1] are negligible and k<p
                    // kase = 2     if s(k) is negligible and k<p
                    // kase = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
                    // kase = 4     if e(p-1) is negligible (convergence).
                    for (k = p - 2; k >= -1; k--)
                    {
                        if (k == -1)
                            break;

                        if (Math.Abs(e[k]) <= eps*(Math.Abs(_s[k]) + Math.Abs(_s[k + 1])))
                        {
                            e[k] = 0.0;
                            break;
                        }
                    }

                    if (k == p - 2)
                    {
                        kase = 4;
                    }
                    else
                    {
                        int ks;
                        for (ks = p - 1; ks >= k; ks--)
                        {
                            if (ks == k)
                                break;

                            double t = (ks != p ? Math.Abs(e[ks]) : 0.0) + (ks != k + 1 ? Math.Abs(e[ks - 1]) : 0.0);
                            if (Math.Abs(_s[ks]) <= eps*t)
                            {
                                _s[ks] = 0.0;
                                break;
                            }
                        }

                        if (ks == k)
                            kase = 3;
                        else if (ks == p - 1)
                            kase = 1;
                        else
                        {
                            kase = 2;
                            k = ks;
                        }
                    }

                    k++;

                    // Perform the task indicated by kase.
                    switch (kase)
                    {
                            // Deflate negligible s(p).
                        case 1:
                            {
                                double f = e[p - 2];
                                e[p - 2] = 0.0;
                                for (int j = p - 2; j >= k; j--)
                                {
                                    double t = MathHelper.Hypotenuse(_s[j], f);
                                    double cs = _s[j]/t;
                                    double sn = f/t;
                                    _s[j] = t;
                                    if (j != k)
                                    {
                                        f = -sn*e[j - 1];
                                        e[j - 1] = cs*e[j - 1];
                                    }

                                    if (wantv)
                                    {
                                        for (int i = 0; i < _n; i++)
                                        {
                                            t = cs*v[i][j] + sn*v[i][p - 1];
                                            v[i][p - 1] = -sn*v[i][j] + cs*v[i][p - 1];
                                            v[i][j] = t;
                                        }
                                    }
                                }
                            }
                            break;

                            // Split at negligible s(k).
                        case 2:
                            {
                                double f = e[k - 1];
                                e[k - 1] = 0.0;
                                for (int j = k; j < p; j++)
                                {
                                    double t = MathHelper.Hypotenuse(_s[j], f);
                                    double cs = _s[j]/t;
                                    double sn = f/t;
                                    _s[j] = t;
                                    f = -sn*e[j];
                                    e[j] = cs*e[j];
                                    if (wantu)
                                    {
                                        for (int i = 0; i < _m; i++)
                                        {
                                            t = cs*u[i][j] + sn*u[i][k - 1];
                                            u[i][k - 1] = -sn*u[i][j] + cs*u[i][k - 1];
                                            u[i][j] = t;
                                        }
                                    }
                                }
                            }
                            break;

                            // Perform one qr step.
                        case 3:
                            {
                                // Calculate the shift.
                                double scale =
                                    Math.Max(
                                        Math.Max(
                                            Math.Max(Math.Max(Math.Abs(_s[p - 1]), Math.Abs(_s[p - 2])),
                                                     Math.Abs(e[p - 2])), Math.Abs(_s[k])), Math.Abs(e[k]));
                                double sp = _s[p - 1]/scale;
                                double spm1 = _s[p - 2]/scale;
                                double epm1 = e[p - 2]/scale;
                                double sk = _s[k]/scale;
                                double ek = e[k]/scale;
                                double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
                                double c = (sp*epm1)*(sp*epm1);
                                double shift = 0.0;
                                if ((b != 0.0) | (c != 0.0))
                                {
                                    shift = Math.Sqrt(b*b + c);
                                    if (b < 0.0)
                                        shift = -shift;
                                    shift = c/(b + shift);
                                }

                                double f = (sk + sp)*(sk - sp) + shift;
                                double g = sk*ek;

                                // Chase zeros.
                                for (int j = k; j < p - 1; j++)
                                {
                                    double t = MathHelper.Hypotenuse(f, g);
                                    double cs = f/t;
                                    double sn = g/t;
                                    if (j != k)
                                        e[j - 1] = t;
                                    f = cs*_s[j] + sn*e[j];
                                    e[j] = cs*e[j] - sn*_s[j];
                                    g = sn*_s[j + 1];
                                    _s[j + 1] = cs*_s[j + 1];
                                    if (wantv)
                                    {
                                        for (int i = 0; i < _n; i++)
                                        {
                                            t = cs*v[i][j] + sn*v[i][j + 1];
                                            v[i][j + 1] = -sn*v[i][j] + cs*v[i][j + 1];
                                            v[i][j] = t;
                                        }
                                    }

                                    t = MathHelper.Hypotenuse(f, g);
                                    cs = f/t;
                                    sn = g/t;
                                    _s[j] = t;
                                    f = cs*e[j] + sn*_s[j + 1];
                                    _s[j + 1] = -sn*e[j] + cs*_s[j + 1];
                                    g = sn*e[j + 1];
                                    e[j + 1] = cs*e[j + 1];
                                    if (wantu && (j < _m - 1))
                                    {
                                        for (int i = 0; i < _m; i++)
                                        {
                                            t = cs*u[i][j] + sn*u[i][j + 1];
                                            u[i][j + 1] = -sn*u[i][j] + cs*u[i][j + 1];
                                            u[i][j] = t;
                                        }
                                    }
                                }

                                e[p - 2] = f;
                                iter = iter + 1;
                            }
                            break;

                            // Convergence.
                        case 4:
                            {
                                // Make the singular values positive.
                                if (_s[k] <= 0.0)
                                {
                                    _s[k] = (_s[k] < 0.0 ? -_s[k] : 0.0);
                                    if (wantv)
                                        for (int i = 0; i <= pp; i++)
                                            v[i][k] = -v[i][k];
                                }

                                // Order the singular values.
                                while (k < pp)
                                {
                                    if (_s[k] >= _s[k + 1])
                                        break;

                                    double t = _s[k];
                                    _s[k] = _s[k + 1];
                                    _s[k + 1] = t;
                                    if (wantv && (k < _n - 1))
                                        for (int i = 0; i < _n; i++)
                                        {
                                            t = v[i][k + 1];
                                            v[i][k + 1] = v[i][k];
                                            v[i][k] = t;
                                        }

                                    if (wantu && (k < _m - 1))
                                        for (int i = 0; i < _m; i++)
                                        {
                                            t = u[i][k + 1];
                                            u[i][k + 1] = u[i][k];
                                            u[i][k] = t;
                                        }

                                    k++;
                                }

                                iter = 0;
                                p--;
                            }
                            break;
                    }
                }
            }

            public double Condition
            {
                get { return _s[0]/_s[Math.Min(_m, _n) - 1]; }
            }

            public double Norm2
            {
                get { return _s[0]; }
            }

            public int Rank
            {
                get
                {
                    double eps = Math.Pow(2.0, -52.0);
                    double tol = Math.Max(_m, _n)*_s[0]*eps;
                    int r = 0;
                    for (int i = 0; i < _s.Length; i++)
                        if (_s[i] > tol)
                            r++;
                    return r;
                }
            }

            public double[] Diagonal
            {
                get { return _s; }
            }
        }
    }
}