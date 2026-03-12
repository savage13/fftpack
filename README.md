fftpack in javascript
---------------------
Implementation of fftpack in javascript, see https://www.netlib.org/fftpack/

Tries to follow the fftpack C implementation


Example
-------

     import { fft } from './fftpack.js'

     let data = Array(100).fill(0)
     data[50] = 1.0

     // Fast Fourier Transform if input signal
     let Data = fft(data)


Testing
-------
Tests are based on scipy version 1.17.1 (Feb 2026) fftpack

Functions
---------

     function fft(x: number[]): number[]
         Compute the Fast Fourier Transform of a real input signal x
         Converts to a complex number and computes 
         Arguments:
         - x - input signal ( number[] ) length n (Real)
         Returns:
         - X - complex (re, im) pairs of transform (number[]) length 2 * n

     function ifft(x: number[]): number[]
         Compute the Inverse Fast Fourier Transform of a complex signal
         Arguments:
         - x - input signal ( number[] ) length 2*n (Real, Imag, Real, Imag ...)
         Returns:
         - x - complex signal ( number() ) length 2*n (Real, Imag, Real, Imag ...)

     function fftfreq(n: number, dt: number): number[]
         Sample frequencies for the Fourier Transform of length n and sample spacing dt
         Argument:
         - n - Input signal length
         - dt - Sample spacing or time step
         Returns:
         - f - fourier sample frequencies
               f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (dt*n)   if n is even
               f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (dt*n)   if n is odd

License
-------
BSD 2-Clause license

fftpack is in the public domain. This implementation should not be
at odds with the original fftpack (Fortran) or its C version

