'use strict';
/* globals it describe */

const _Math = require('minimatrix');
const roots = require('./index.js');
const expect = require('chai').expect;

const TOL = 1e-14;

describe('LU Decomposition', () => {
  it('should solve correctly', () => {
    const AA = new _Math.Matrix3();
    AA.set(3, -0.1, -0.2, 0.1, 7, -0.3, 0.3, -0.2, 10);
    const b = new _Math.Vector3();
    b.x = 7.85;
    b.y = -19.3;
    b.z = 71.4;
    const luInfo = roots.luDecomposition(AA);
    const P = luInfo.P;
    const A = luInfo.A;
    const x = roots.luSolve(A, P, b);
    expect(_Math.abs(x.x - 3) < TOL).to.be.true;
    expect(_Math.abs(x.y - -2.5) < TOL).to.be.true;
    expect(_Math.abs(x.z - 7) < TOL).to.be.true;
  });
});
