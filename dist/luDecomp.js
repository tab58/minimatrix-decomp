"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("./helpers");
var LUSolver = /** @class */ (function () {
    function LUSolver() {
    }
    LUSolver.decompose = function (AA) {
        // adapted Crout decomposition with ideas from:
        // https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/class06.pdf
        var A = AA.clone();
        var a = A.toArray([], 0); // indexed via a_ij = a[i + j * n]
        var n = A.colDimension;
        var P = [];
        var rowScalers = [];
        for (var i = 0; i < n; ++i) {
            P[i] = i;
            var col = helpers_1.findLargestInRow(a, n, i);
            var scaler = a[col * n + i]; // row scaling
            if (scaler === 0) {
                throw new Error('LUSolver.decompose(): matrix is singular and cannot be LU factorized.');
            }
            rowScalers[i] = scaler; // don't actually want to scale the matrix or else (PA != LU).
        }
        // iterate over columns to reduce the matrix
        // implicitly reduce the matrix
        for (var j = 0; j < n; ++j) {
            for (var i = 0; i < j; ++i) {
                // compute upper tri, computes values across the row
                var sum = 0;
                for (var k = 0; k < i; ++k) {
                    sum += a[i + k * n] * a[k + j * n]; // avoid big + small roundoffs
                }
                a[i + j * n] -= sum;
            }
            var pivotLrgElem = 0;
            var pivotIndex = j;
            for (var i = j; i < n; ++i) {
                // compute lower tri and diagonal, computes values down the column
                var sum = 0;
                for (var k = 0; k < j; ++k) {
                    sum += a[i + k * n] * a[k + j * n]; // avoid big + small roundoffs
                }
                a[i + j * n] -= sum;
                // find the pivot element
                var pivotTest = Math.abs(a[i + j * n]) / rowScalers[i];
                if (pivotTest > pivotLrgElem) {
                    pivotLrgElem = pivotTest;
                    pivotIndex = i;
                }
            }
            // swap the row on the pivot
            helpers_1.swapMatrixArrayRows(a, n, j, pivotIndex);
            helpers_1.swapArrayElements(P, j, pivotIndex);
            helpers_1.swapArrayElements(rowScalers, j, pivotIndex);
            // divide by the pivot
            if (j < n - 1) {
                var pivotScale = a[j + j * n];
                for (var i = j + 1; i < n; ++i) {
                    a[i + j * n] /= pivotScale;
                }
            }
        }
        A.fromArray(a, 0);
        return {
            P: P,
            A: A
        };
    };
    LUSolver.solveLinear = function (A, P, b) {
        // Since PA = LU, then L(Ux) = Pb
        var a = A.toArray([], 0);
        var n = Math.floor(Math.sqrt(a.length));
        var x = b.clone().setScalar(0);
        // L * y  = P * b, solve for y.
        // Implicit 1's on the diagonal.
        for (var i = 0; i < n; ++i) {
            var sum = 0;
            for (var j = 0; j < i; ++j) {
                sum += a[i + j * n] * x.getComponent(j);
            }
            var xi = b.getComponent(P[i]) - sum;
            x.setComponent(i, xi);
        }
        // U * x = y  ==> i = n-1 -> 0, j = 0
        for (var i = n - 1; i >= 0; --i) {
            var sum = 0;
            for (var j = i + 1; j < n; ++j) {
                sum += a[i + j * n] * x.getComponent(j);
            }
            var scale = a[i + i * n];
            if (scale === 0) {
                console.warn("LUSolver.solve(): x[" + i + "] is free.");
            }
            x.setComponent(i, (x.getComponent(i) - sum) / scale);
        }
        return x;
    };
    return LUSolver;
}());
exports.LUSolver = LUSolver;
//# sourceMappingURL=luDecomp.js.map