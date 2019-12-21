import { Matrix } from 'minimatrix';

function determinant2x2 (a11: number, a12: number, a21: number, a22: number) {
  return a11 * a22 - a12 * a21;
}

export class SVDSolver {
  // public static decompose (A: Matrix): { U: Matrix, S: Matrix, Vt: Matrix } {
  //   const n = A.colDimension;
  //   // calculate Wilkinson horizontal shift
  //   const mu = determinant2x2(A.get(n - 1, n - 1), A.get(n - 1, n), A.get(n, n - 1), A.get(n, n));
  // }
}