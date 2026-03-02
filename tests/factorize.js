
import { factorize } from './fftpack.ts'
import { deepEqualTolerance } from '../tests/deepEqualTolerance.js'
function eq(a, b, tol = 2e-15) { return deepEqualTolerance(a, b, tol) }

eq(factorize(2), [2, 1, 2])
eq(factorize(3), [3, 1, 3])
eq(factorize(4), [4, 1, 4])
eq(factorize(5), [5, 1, 5])
eq(factorize(6), [6, 2, 2, 3])
eq(factorize(7), [7, 1, 7])
eq(factorize(8), [8, 2, 2, 4])
eq(factorize(9), [9, 2, 3, 3])
eq(factorize(10), [10, 2, 2, 5])
eq(factorize(11), [11, 1, 11])
eq(factorize(12), [12, 2, 3, 4])
eq(factorize(13), [13, 1, 13])
eq(factorize(14), [14, 2, 2, 7])
eq(factorize(36), [36, 3, 3, 3, 4])

