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

     function fft(x)
         Compute the Fast Fourier Transform of input signal x
         Arguments:
         - x - input signal ( number[] ) length n
         Returns:
         - X - complex (re, im) pairs of transform (number[]) length 2 * n


License
-------
BSD 2-Clause license

fftpack is in the public domain. This implementation should not be 
at odds with the original fftpack or its C version

