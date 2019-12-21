import { findLargestInRow, swapArrayElements, swapMatrixArrayRows } from './helpers';

export class LUSolver {
  public static decompose <T extends Matrix>(AA: T): { P: number[], A: T } {
    // adapted Crout decomposition with ideas from:
    // https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/class06.pdf
  
    const A = AA.clone();
    const a = A.toArray([], 0); // indexed via a_ij = a[i + j * n]
    const n = A.colDimension;
    const P = [];
    const rowScalers = [];
    for (let i = 0; i < n; ++i) {
      P[i] = i;
      const col = findLargestInRow(a, n, i);
      const scaler = a[col * n + i]; // row scaling
      if (scaler === 0) {
        throw new Error('LUSolver.decompose(): matrix is singular and cannot be LU factorized.');
      }
      rowScalers[i] = scaler; // don't actually want to scale the matrix or else (PA != LU).
    }
  
    // iterate over columns to reduce the matrix
    // implicitly reduce the matrix
    for (let j = 0; j < n; ++j) {
      for (let i = 0; i < j; ++i) {
        // compute upper tri, computes values across the row
        let sum = 0;
        for (let k = 0; k < i; ++k) {
          sum += a[i + k * n] * a[k + j * n]; // avoid big + small roundoffs
        }
        a[i + j * n] -= sum;
      }

      let pivotLrgElem = 0;
      let pivotIndex = j;
      for (let i = j; i < n; ++i) {
        // compute lower tri and diagonal, computes values down the column
        let sum = 0;
        for (let k = 0; k < j; ++k) {
          sum += a[i + k * n] * a[k + j * n]; // avoid big + small roundoffs
        }
        a[i + j * n] -= sum;
  
        // find the pivot element
        const pivotTest = Math.abs(a[i + j * n]) / rowScalers[i];
        if (pivotTest > pivotLrgElem) {
          pivotLrgElem = pivotTest;
          pivotIndex = i;
        }
      }
      // swap the row on the pivot
      swapMatrixArrayRows(a, n, j, pivotIndex);
      swapArrayElements(P, j, pivotIndex);
      swapArrayElements(rowScalers, j, pivotIndex);
      // divide by the pivot
      if (j < n - 1) {
        const pivotScale = a[j + j * n];
        for (let i = j + 1; i < n; ++i) {
          a[i + j * n] /= pivotScale;
        }
      }
    }

    A.fromArray(a, 0);
    return {
      P,
      A
    };
  }

  public static solveLinear <T extends Matrix, U extends Vector>(A: T, P: number[], b: U): U {
    // Since PA = LU, then L(Ux) = Pb
    const a = A.toArray([], 0);
    const n = Math.floor(Math.sqrt(a.length));
    const x = b.clone().setScalar(0);
    // L * y  = P * b, solve for y.
    // Implicit 1's on the diagonal.
    for (let i = 0; i < n; ++i) {
      let sum = 0;
      for (let j = 0; j < i; ++j) {
        sum += a[i + j * n] * x.getComponent(j);
      }
      const xi = b.getComponent(P[i]) - sum;
      x.setComponent(i, xi);
    }

    // U * x = y  ==> i = n-1 -> 0, j = 0
    for (let i = n - 1; i >= 0; --i) {
      let sum = 0;
      for (let j = i + 1; j < n; ++j) {
        sum += a[i + j * n] * x.getComponent(j);
      }
      const scale = a[i + i * n];
      if (scale === 0) {
        console.warn(`LUSolver.solve(): x[${i}] is free.`);
      }
      x.setComponent(i, (x.getComponent(i) - sum) / scale);
    }
    return x;
  }
}