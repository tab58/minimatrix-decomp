import { Vector } from 'minimatrix';
/**
 * Swaps rows in-place in the matrix. Zero is the first row.
 */
export declare const swapMatrixArrayRows: (A: number[], n: number, i: number, j: number) => number[];
/** Swaps two values in an array. */
export declare const swapArrayElements: (A: number[], i: number, j: number) => number[];
/** Swaps components in a vector. */
export declare const swapVectorComponents: (A: Vector, i: number, j: number) => Vector;
/** Finds the largest element in the row of the matrix. */
export declare const findLargestInRow: (m: number[], n: number, i: number) => number;
/** Finds the largest element in the column of the matrix. */
export declare const findLargestInCol: (m: number[], n: number, i: number) => number;
/** 2-norm function that avoids overflow and underflow. */
export declare const hypot: (a: number, b: number) => number;
/**
 * Get the values of the Householder reflection for the vector.
 * @param x The array of values for the vector to reflect.
 * @param n The number of values to take from the array.
 * @param offset The starting offset for the array.
 */
export declare const getHouseholderVectorValues: (x: number[], n: number, offset?: number) => number[];
