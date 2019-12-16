import { Vector, Matrix } from 'minimatrix';

/**
 * Swaps rows in-place in the matrix. Zero is the first row.
 */
export const swapMatrixArrayRows = (A: number[], n: number, i: number, j: number): number[] => {
  if (i > n - 1 || j > n - 1) {
    throw new Error(`swapRows(): row index out of bounds.`);
  }
  if (i !== j) {
    for (let k = 0; k < n; ++k) {
      const offset = k * n;
      const tmp = A[i + offset];
      A[i + offset] = A[j + offset];
      A[j + offset] = tmp;
    }
  }
  return A;
}

/** Swaps two values in an array. */
export const swapArrayElements = (A: number[], i: number, j: number): number[] => {
  if (i !== j) {
    const tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
  }
  return A;
};

/** Swaps components in a vector. */
export const swapVectorComponents = (A: Vector, i: number, j: number): Vector => {
  if (i !== j) {
    const tmp = A.getComponent(i);
    A.setComponent(i, A.getComponent(j));
    A.setComponent(j, tmp);
  }
  return A;
};

/** Finds the largest element in the row of the matrix. */
export const findLargestInRow = (m: number[], n: number, i: number): number => {
  const rowOffset = i;

  let colOffset = 0;
  let lrgElem = Math.abs(m[rowOffset]);
  let lrgCol = colOffset;
  for (colOffset = 1; colOffset < n; ++colOffset) {
    const val = Math.abs(m[rowOffset + colOffset * n]);
    if (val > lrgElem) {
      lrgCol = colOffset;
      lrgElem = val;
    }
  }
  return lrgCol;
};

/** Finds the largest element in the column of the matrix. */
export const findLargestInCol = (m: number[], n: number, i: number): number => {
  const rowOffset = i * n;

  let columnOffset = 0;
  let lrgElem = Math.abs(m[rowOffset]);
  let lrgRow = columnOffset;
  for (columnOffset = 1; columnOffset < n; ++columnOffset) {
    const val = Math.abs(m[rowOffset + columnOffset]);
    if (val > lrgElem) {
      lrgRow = columnOffset;
      lrgElem = val;
    }
  }
  return lrgRow;
};