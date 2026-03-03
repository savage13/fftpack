
import assert from 'node:assert'

const TOL = 1e-15

export function deepEqualTolerance(a, b, tol = TOL) {
    try {
        __deepEqualTolerance(a, b, tol)
    } catch (error) {
        console.log(`Error: ${error.message}`)
        assert.deepStrictEqual(a, b)
    }
}
function __deepEqualTolerance(a, b, tol) {
    if(typeof a != typeof b) {
        let ta = typeof a
        let tb = typeof b
        throw new Error(`inconsisent types ${ta} ${tb}`)
    }
    if(Array.isArray(a) != Array.isArray(b)) {
        throw new Error('inconsisent types Arrays')
    }
    if(Array.isArray(a)) {
        if(a.length != b.length) {
            throw new Error("arrays inconsistent length")
        }
        for(let i = 0; i < a.length; i++) {
            __deepEqualTolerance(a[i],b[i],tol)
        }
    } else if(typeof a == 'string') {
        if(a != b) {
            throw new Error(`${a} != ${b}`)
        }
    } else if(typeof a == 'number') {
        if(Number.isNaN(a) && Number.isNaN(b)) {
        } else if(Number.isNaN(a) != Number.isNaN(b)) {
            throw new Error(`NaN a: ${a} b: ${b}`)
        } else {
            let diff = Math.abs(a - b)
            if(diff > tol) {
                throw new Error(`abs(${a} - ${b}) (${diff}) > ${tol}`)
            }
        }
    } else if(typeof a == 'object') {
        let key_a = Object.keys(a).sort()
        let key_b = Object.keys(a).sort()
        if(key_a.length != key_b.length) {
            throw new Error("objects with inconsistent key lengths")
        }
        __deepEqualTolerance(key_a, key_b,tol)
        for(const key of key_a) {
            __deepEqualTolerance(a[key], b[key], tol)
        }
    } else {
        if(a != b) {
            throw new Error(`assertion failed ${a} ${b}`)
        }
    }
}
