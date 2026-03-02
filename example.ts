import { fft } from './fftpack.js'

let data = Array(100).fill(0)
data[50] = 1.0

// Fast Fourier Transform if input signal
let Data = fft(data)
console.log(Data)


