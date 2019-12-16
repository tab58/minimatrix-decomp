import { Vector, Matrix, LinAlgHelpers } from 'minimatrix';

/** 2-norm function that avoids overflow and underflow. */
const hypot = (a: number, b: number): number => {
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

/**
 * Get the values of the Householder reflection for the vector.
 * @param x The array of values for the vector to reflect.
 * @param n The number of values to take from the array.
 * @param offset The starting offset for the array.
 */
const getHouseholderVectorValues = (x: number[], n: number, offset: number = 0): number[] => {
  if (n === 0) {
    throw new Error(`getHouseholderVector(): length of vector is zero.`);
  }
  const v = [];
  let alpha = 0;
  for (let i = 0; i < n; ++i) {
    const xi = x[offset + i];
    v.push(xi);
    alpha = hypot(alpha, xi);
  }
  const [ x0 ] = v;
  // construct v = u / |u|, u = x - (alpha)*e1, alpha = |x|
  v[0] -= Math.abs(alpha);
  const vLen = 1.0 / Math.sqrt(2 * alpha * (alpha - x0));
  for (let i = 0; i < n; ++i) {
    v[i] *= vLen;
  }
  return v;
}

/** Class for QR decomposition and linear equation solving. */
export class QRSolver {
  /** 
   * Decomposes the matrix using Householder reflections into
   * an orthogonal matrix Q and a right triangular matrix R.
   */
  public static decompose <T extends Matrix>(A: T): { Q: T, R: T } {
    const n = A.colDimension;
    const Q = A.clone().identity();
    const R = A.clone();
    const x = [];
    for (let k = 0; k < n; ++k) { x.push(0); }
    
    // cycle through columns
    const H = A.clone();
    for (let j = 0; j < n; ++j) {
      H.identity();
      const nrows = n - j;
      if (nrows > 1) {
        // get householder vector
        for (let i = j; i < n; ++i) { x[i] = R.get(i, j); }
        const vv = getHouseholderVectorValues(x, nrows, j);
        const v = LinAlgHelpers.vectorFromValues(vv, nrows, 0);

        // calculate the Householder matrix H
        const HH = LinAlgHelpers.getOuterProduct(v).multiplyScalar(-2);
        for (let i = 0; i < nrows; ++i) {
          for (let k = 0; k < nrows; ++k) {
            const t = (i === k ? 1 : 0) + HH.get(i, k); // HH = I - 2 * v * v^T
            H.set(j + i, j + k, t);
          }
        }

        // R' = H * R, Q' = Q * H^T
        R.premultiply(H);
        H.transpose();
        Q.multiply(H);
      }
    }
    return { Q, R };
  }

  public static solveLinear <T extends Matrix, U extends Vector>(Q: T, R: T, b: U): U {
    const n = b.dimension;
    // y = Q^T * b
    const Qt = Q.clone().transpose();
    const x = LinAlgHelpers.transformVector(Qt, b) as U;

    // R * x = y
    for (let i = n - 1; i >= 0; --i) {
      const y = x.getComponent(i);
      const rii = R.get(i, i);
      let sum = 0;
      for (let k = i + 1; k < n; ++k) { sum += R.get(i, k) * x.getComponent(k); }
      x.setComponent(i, (y - sum) / rii);
    }
    return x;
  }
}