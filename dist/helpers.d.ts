import { Vector, Matrix } from 'minimatrix';
/**
 * Swaps rows in-place in the matrix. Zero is the first row.
 */
export declare const swapMatrixArrayRows: (A: number[], n: number, i: number, j: number) => number[];
/** Swaps two values in an array. */
export declare const swapArrayElements: (A: number[], i: number, j: number) => number[];
/** Swaps components in a vector. */
export declare const swapVectorComponents: (A: Vector, i: number, j: number) => Vector;
/** Finds the largest element in the row of the matrix. */
export declare const findLargestInRow: (M: Matrix, i: number) => number;
/** Finds the largest element in the column of the matrix. */
export declare const findLargestInCol: (M: Matrix, i: number) => number;
