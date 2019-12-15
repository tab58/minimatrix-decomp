import { Vector, Matrix } from 'minimatrix';
export declare class LUSolver {
    static decompose(AA: Matrix): {
        P: number[];
        A: Matrix;
    };
    static solve(A: Matrix, P: number[], b: Vector): Vector;
}
