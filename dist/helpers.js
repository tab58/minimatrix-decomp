"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * Swaps rows in-place in the matrix. Zero is the first row.
 */
exports.swapMatrixArrayRows = (A, n, i, j) => {
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
};
/** Swaps two values in an array. */
exports.swapArrayElements = (A, i, j) => {
    if (i !== j) {
        const tmp = A[i];
        A[i] = A[j];
        A[j] = tmp;
    }
    return A;
};
/** Swaps components in a vector. */
exports.swapVectorComponents = (A, i, j) => {
    if (i !== j) {
        const tmp = A.getComponent(i);
        A.setComponent(i, A.getComponent(j));
        A.setComponent(j, tmp);
    }
    return A;
};
/** Finds the largest element in the row of the matrix. */
exports.findLargestInRow = (M, i) => {
    const m = M.toArray([], 0);
    const n = M.colDimension;
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
exports.findLargestInCol = (M, i) => {
    const m = M.toArray([], 0);
    const n = M.colDimension;
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
//# sourceMappingURL=helpers.js.map