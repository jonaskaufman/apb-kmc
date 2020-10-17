import sys
import numpy as np
from grid import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize

def fit_func(x, a):
    return a*x*x

data = np.loadtxt(sys.argv[1], skiprows=1)
wavelength = data[:,0]
inverse_wavelength = 1/wavelength
relaxation_time = data[:,1]
inverse_time = 1/relaxation_time

plt.plot(wavelength, relaxation_time, 'ko')
plt.xlabel('wavelength')
plt.ylabel('relaxation time')
plt.tight_layout()
plt.show()


opt, cov = scipy.optimize.curve_fit(fit_func, inverse_wavelength, inverse_time, p0=[1])
print(opt)
samples = np.linspace(0, max(inverse_wavelength), 10)
fitted = [fit_func(x, *opt) for x in samples]
plt.plot(samples, fitted, 'r-') 

plt.plot(inverse_wavelength, inverse_time, 'ko')
plt.xlabel('inverse wavelength')
plt.ylabel('inverse relaxation time')
plt.tight_layout()
plt.show()
