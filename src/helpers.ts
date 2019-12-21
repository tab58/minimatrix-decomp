import { Vector2, Vector3, Vector4, Matrix2, Matrix3, Matrix4, LinAlgHelpers } from 'minimatrix';

export type Vector = Vector2 | Vector3 | Vector4;
export type Matrix = Matrix2 | Matrix3 | Matrix4;

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

/** 2-norm function that avoids overflow and underflow. */
export const hypot = (a: number, b: number): number => {
  if (a === 0 && b === 0) {
    return 0;
  }
  const x = Math.abs(a);
  const y = Math.abs(b);
  const t = Math.min(x, y);
  const u = Math.max(x, y);
  const w = t / u;
  return u * Math.sqrt(1 + w * w);
}

export class Householder {
  /**
   * Get the values of the unit Householder reflection vector for the matrix vector.
   * @param x The column of the matrix.
   * @param offset The starting component offset for the vector.
   */
  public static getVector (x: Vector): { v: Vector, beta: number } {
    const n = x.dimension;
    if (n === 0) { throw new Error(`getVector(): length of vector is zero.`); }
    if (n > 4 || n < 2) { throw new Error(`getVector(): length of vector is not Vector2, Vector3, or Vector4.`); }

    const v = x.clone();
    const x0 = v.getComponent(0);
    v.setComponent(0, 0);
    const sigma = v.lengthSq();
    v.setComponent(0, 1);
    
    let beta = 0;
    if (sigma === 0 && x0 >= 0) {
      beta = 0;
    } else if (sigma === 0 && x0 < 0) {
      beta = -2;
    } else {
      const mu = hypot(sigma, x0);
      const v0 = (x0 <= 0)
        ? x0 - mu
        : -sigma / (x0 + mu);
      v.setComponent(0, v0);
      beta = (2 * v0 * v0) / (sigma + (v0 * v0));
      v.divideScalar(v0);
    }
    return { v, beta };
  }

  /**
   * Pre-multiplies the Householder matrix P = I - beta * v * v^T to a matrix.
   * @param A The matrix to multiply.
   * @param x The Householder vector.
   */
  public static premultiply (A: Matrix, x: Vector, beta: number): Matrix {    
    const n = A.colDimension;
    const m = x.dimension;
    if (m > n) { throw new Error(`premultiply(): Householder vector dimension is larger than matrix column dimension.`); }
    if (n > 4 || n < 2) { throw new Error(`premultiply(): length of matrix column is not 2, 3, or 4.`); }

    // set the Householder vector at the end of the outer product vector.    
    const v = A.getColumn(0).setScalar(0);
    for (let i = m - 1; i >= 0; --i) { v.setComponent(n - m + i, x.getComponent(i)); }
    
    // PA = A - (beta * v) (v^T * A)
    const b = LinAlgHelpers.transformRowVector(A, v.clone());
    const a = v.multiplyScalar(beta);
    LinAlgHelpers.addMatrixOuterProduct(A, a, b, -1);
    return A;
  }

  /**
   * Post-multiplies the Householder matrix P = I - beta * v * v^T to a matrix.
   * @param m The matrix elements.
   * @param n The matrix column dimension.
   * @param col The column number.
   */
  public static postmultiplyMatrixColumnReflection (A: Matrix, x: Vector, beta: number): Matrix {  
    const n = A.rowDimension;
    const m = x.dimension;
    if (m > n) { throw new Error(`premultiply(): Householder vector dimension is larger than matrix row dimension.`); }
    if (n > 4 || n < 2) { throw new Error(`premultiply(): length of matrix row is not 2, 3, or 4.`); }

    // set the Householder vector at the end of the outer product vector.    
    const v = A.getRow(0).setScalar(0);
    for (let i = m - 1; i >= 0; --i) { v.setComponent(n - m + i, x.getComponent(i)); }
    
    // AP = A - (A * v) (beta * v)^T
    const a = LinAlgHelpers.transformVector(A, v.clone());
    const b = v.multiplyScalar(beta);
    LinAlgHelpers.addMatrixOuterProduct(A, a, b, -1);
    return A;
  }
}

/**
 * Calculates the parameters of a vector [r, 0]^T to rotate it to [a, b]^T.
 * @param a The first parameter of the vector.
 * @param b 
 * @param csr 
 */
function rotg (a: number, b: number, csr: number[] = []): number[] {
  // Based on Algorithm 4 from "Discontinuous Plane
  // Rotations and the Symmetric Eigenvalue Problem"
  // by Anderson, 2000.
  var c = 0;
  var s = 0;
  var r = 0;
  var t = 0;
  var u = 0;

  if (b === 0) {
    c = Math.sign(a);
    s = 0;
    r = Math.abs(a);
  } else if (a === 0) {
    c = 0;
    s = Math.sign(b);
    r = Math.abs(b);
  } else if (Math.abs(a) > Math.abs(b)) {
    t = b / a;
    u = Math.sign(a) * Math.sqrt(1 + t * t);
    c = 1 / u;
    s = t * c;
    r = a * u;
  } else {
    t = a / b;
    u = Math.sign(a) * Math.sqrt(1 + t * t);
    s = 1 / u;
    c = t * s;
    r = b * u;
  }
  // try to save some unnecessary object creation
  csr[0] = c;
  csr[1] = s;
  csr[2] = r;
  return csr;
}

/**
 * Applies a Givens rotation to a matrix.
 * @param M 
 * @param j 
 */
export function applyGivensRotation (M: Matrix, a: number, b: number, i: number, j: number): Matrix {
  const [c, s] = rotg(a, b);
  M.set(i, i, c);
  M.set(i, j, s);
  M.set(j, i , -s);
  M.set(j, j, c);
  return M;
}