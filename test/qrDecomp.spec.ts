import { expect } from 'chai';
import { Vector2, Vector3, Vector4, Matrix2, Matrix3, Matrix4 } from 'minimatrix';
import { QRSolver } from '../src/index';

const EPS = 1.5e-14;
describe('QR Decomposition', (): void => {
  it('should solve an Ax=b problem for a 2x2 matrix', (): void => {
    const M = new Matrix2();
    M.setElements(1, 2, 3, -5);
    const { Q, R } = QRSolver.decompose(M);
    const b = new Vector2(4, 1);
    const x = QRSolver.solveLinear(Q, R, b);
    expect(x.dimension).to.be.eql(2);
    expect(Math.abs(x.getComponent(0) - 2)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(1) - 1)).to.be.lessThan(EPS);
  });
  it('should solve an Ax=b problem for a 3x3 matrix', (): void => {
    const M = new Matrix3();
    M.setElements(3, -0.1, -0.2, 0.1, 7, -0.3, 0.3, -0.2, 10);
    const { Q, R } = QRSolver.decompose(M);
    const b = new Vector3(7.85, -19.3, 71.4);
    const x = QRSolver.solveLinear(Q, R, b);
    expect(x.dimension).to.be.eql(3);
    expect(Math.abs(x.getComponent(0) - 3)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(1) - -2.5)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(2) - 7)).to.be.lessThan(EPS);
  });
  it('should solve an Ax=b problem for a 4x4 matrix', (): void => {
    const M = new Matrix4();
    M.setElements(2, -1, 1, 1, 1, 2, -1, -1, -1, 2, 2, 2, 1, -1, 2, 1);
    const { Q, R } = QRSolver.decompose(M);
    const b = new Vector4(6, 3, 14, 8);
    const x = QRSolver.solveLinear(Q, R, b);
    expect(x.dimension).to.be.eql(4);
    expect(Math.abs(x.getComponent(0) - 2)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(1) - 3)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(2) - 4)).to.be.lessThan(EPS);
    expect(Math.abs(x.getComponent(3) - 1)).to.be.lessThan(EPS);
  });
});