'use strict';

// default numerical tolerance
// const TOL = 1e-14;
const consoleWarning = console.warn;

const helpers = {
  swapValuesInArray: (A, i, j) => {
    if (i !== j) {
      const tmp = A[i];
      A[i] = A[j];
      A[j] = tmp;
    }
    return A;
  },
  swapValuesInVector: (A, i, j) => {
    if (i !== j) {
      const tmp = A.getComponent(i);
      A.setComponent(i, A.getComponent(j));
      A.setComponent(j, tmp);
    }
    return A;
  },
  swapRowsInMatrixArray: (A, i, j) => {
    if (i !== j) {
      const n = Math.sqrt(A.length); // assuming M is square because it's a minimatrix
      for (let k = 0; k < n; ++k) {
        const offset = k * n;
        const tmp = A[i + offset];
        A[i + offset] = A[j + offset];
        A[j + offset] = tmp;
      }
    }
    return A;
  },
  findLargestInRow: (M, i) => {
    // get the row:
    const m = M.elements;
    const n = Math.sqrt(m.length); // assuming M is square because it's a minimatrix
    const offset = i;

    let j = 0;
    let lrgElem = Math.abs(m[offset]);
    let lrgCol = j;
    for (j = 1; j < n; ++j) {
      const val = Math.abs(m[offset + j * n]);
      if (val > lrgElem) {
        lrgCol = j;
        lrgElem = val;
      }
    }
    return lrgCol;
  },
  findLargestInCol: (M, i) => {
    // get the column:
    const m = M.elements;
    const n = Math.sqrt(m.length); // assuming M is square because it's a minimatrix
    const offset = i * n;

    let j = 0;
    let lrgElem = Math.abs(m[offset]);
    let lrgRow = j;
    for (j = 1; j < n; ++j) {
      const val = Math.abs(m[offset + j]);
      if (val > lrgElem) {
        lrgRow = j;
        lrgElem = val;
      }
    }
    return lrgRow;
  }
};

module.exports = {
  luDecomposition: (AA) => {
    // adapted Crout decomposition with ideas from:
    // https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/class06.pdf

    const A = AA.clone();
    const a = A.elements; // indexed via a_ij = a[i + j * n]
    const n = Math.floor(Math.sqrt(a.length));
    const P = [];
    const rowScalers = [];
    for (let i = 0; i < n; ++i) {
      P[i] = i;
      const col = helpers.findLargestInRow(A, i);
      const scaler = a[col * n + i]; // row scaling
      if (scaler === 0) {
        throw new Error('luDecompAndSolve(): matrix is singular and cannot be LU factorized.');
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
        a[i + j * n] = a[i + j * n] - sum;
      }
      let pivotLrgElem = 0;
      let pivotIndex = j;
      for (let i = j; i < n; ++i) {
        // compute lower tri and diagonal, computes values down the column
        let sum = 0;
        for (let k = 0; k < j; ++k) {
          sum += a[i + k * n] * a[k + j * n]; // avoid big + small roundoffs
        }
        a[i + j * n] = a[i + j * n] - sum;

        // find the pivot element
        const pivotTest = rowScalers[i] * Math.abs(a[i + j * n]);
        if (pivotTest > pivotLrgElem) {
          pivotLrgElem = pivotTest;
          pivotIndex = i;
        }
      }
      helpers.swapRowsInMatrixArray(a, j, pivotIndex);
      helpers.swapValuesInArray(P, j, pivotIndex);
      helpers.swapValuesInArray(rowScalers, j, pivotIndex);
      if (j < n - 1) {
        const pivotScale = a[j + j * n];
        for (let i = j + 1; i < n; ++i) {
          a[i + j * n] /= pivotScale;
        }
      }
    }
    return {
      P,
      A
    };
  },
  luSolve: (A, P, b) => {
    // Since PA = LU, then L(U x) = Pb
    const a = A.elements;
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
        consoleWarning(`luSolve(): x[${i}] is free.`);
      }
      x.setComponent(i, (x.getComponent(i) - sum) / scale);
    }
    return x;
  }
};
