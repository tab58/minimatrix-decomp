import { Vector, Matrix } from 'minimatrix';
export declare class LUSolver {
    static decompose<T extends Matrix>(AA: T): {
        P: number[];
        A: T;
    };
    static solveLinear<T extends Matrix, U extends Vector>(A: T, P: number[], b: U): U;
}
