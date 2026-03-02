
import json
import numpy as np
from scipy import fftpack

def imp(n, k=0):
    v = np.zeros(n)
    v[k] = 1
    return v
def zflat(z):
    out = []
    for v in z:
        out.extend([v.real, v.imag])
    return out
tests = []

n = 8
for n in range(2,256+1):
    for i in np.linspace(0,n,5, dtype='int',endpoint=False):
        input = imp(n, i)
        tests.append({'input': list(input), 'output': zflat( fftpack.fft(input) )})
    for k in [1,2,-1,5]:
        input = np.ones(n) * k
        tests.append({'input': list(input), 'output': zflat( fftpack.fft(input) )})

    #input = np.linspace(0, 10, n,dtype=np.float64)
    #tests.append({'input': list(input), 'output': zflat( fftpack.fft(input) )})

    #input = np.linspace(0, 10, n,dtype=np.float64) * -1
    #tests.append({'input': list(input), 'output': zflat( fftpack.fft(input) )})

    for i in range(10):
        input = np.random.rand(n) * 2 - 1.0
        tests.append({'input': list(input), 'output': zflat( fftpack.fft(input) )})

json.dump(tests, open('cfft.json','w'), indent=2)
