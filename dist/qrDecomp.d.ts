import { Vector, Matrix } from 'minimatrix';
/** Class for QR decomposition and linear equation solving. */
export declare class QRSolver {
    /**
     * Decomposes the matrix using Householder reflections into
     * an orthogonal matrix Q and a right triangular matrix R.
     */
    static decompose<T extends Matrix>(A: T): {
        Q: T;
        R: T;
    };
    static solveLinear<T extends Matrix, U extends Vector>(Q: T, R: T, b: U): U;
}
