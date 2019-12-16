"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var minimatrix_1 = require("minimatrix");
var helpers_1 = require("./helpers");
/** Class for QR decomposition and linear equation solving. */
var QRSolver = /** @class */ (function () {
    function QRSolver() {
    }
    /**
     * Decomposes the matrix using Householder reflections into
     * an orthogonal matrix Q and a right triangular matrix R.
     */
    QRSolver.decompose = function (A) {
        var n = A.colDimension;
        var Q = A.clone().identity();
        var R = A.clone();
        var x = [];
        for (var k = 0; k < n; ++k) {
            x.push(0);
        }
        // cycle through columns
        var H = A.clone();
        for (var j = 0; j < n; ++j) {
            H.identity();
            var nrows = n - j;
            if (nrows > 1) {
                // get householder vector
                for (var i = j; i < n; ++i) {
                    x[i] = R.get(i, j);
                }
                var vv = helpers_1.getHouseholderVectorValues(x, nrows, j);
                var v = minimatrix_1.LinAlgHelpers.vectorFromValues(vv, nrows, 0);
                // calculate the Householder matrix H
                var HH = minimatrix_1.LinAlgHelpers.getOuterProduct(v).multiplyScalar(-2);
                for (var i = 0; i < nrows; ++i) {
                    for (var k = 0; k < nrows; ++k) {
                        var t = (i === k ? 1 : 0) + HH.get(i, k); // HH = I - 2 * v * v^T
                        H.set(j + i, j + k, t);
                    }
                }
                // R' = H * R, Q' = Q * H^T
                R.premultiply(H);
                H.transpose();
                Q.multiply(H);
            }
        }
        return { Q: Q, R: R };
    };
    QRSolver.solveLinear = function (Q, R, b) {
        var n = b.dimension;
        // y = Q^T * b
        var Qt = Q.clone().transpose();
        var x = minimatrix_1.LinAlgHelpers.transformVector(Qt, b);
        // R * x = y
        for (var i = n - 1; i >= 0; --i) {
            var y = x.getComponent(i);
            var rii = R.get(i, i);
            var sum = 0;
            for (var k = i + 1; k < n; ++k) {
                sum += R.get(i, k) * x.getComponent(k);
            }
            x.setComponent(i, (y - sum) / rii);
        }
        return x;
    };
    return QRSolver;
}());
exports.QRSolver = QRSolver;
//# sourceMappingURL=qrDecomp.js.map