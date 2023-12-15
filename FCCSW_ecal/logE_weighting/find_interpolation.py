import matplotlib.pyplot as plt
import math as ma
from scipy.optimize import curve_fit
import numpy as np


plt.style.use('bmh')


 
min_w0 = [2.5, 2.75, 4.5, 5.75, 6.0, 6.0, 6.25, 6.25, 6.25, 6.25, 6.5, 6.25, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5] 
energies = [0.1, 0.5, 2.0, 7.0, 12.0, 17.0, 22.0, 25.0, 32.0, 35.0, 37.0, 42.0, 45.0, 47.0, 52.0, 55.0, 57.0, 62.0, 65.0, 67.0, 72.0, 75.0, 77.0, 82.0, 85.0, 87.0, 92.0, 95.0, 97.0]


# f_inv = lambda e, alpha, beta: alpha-beta/np.sqrt(e)
f_inv = lambda e, alpha, beta: alpha-beta/np.sqrt(e)
f_exp = lambda e, alpha, beta: alpha*(1-np.exp(-beta*np.sqrt(e)))
f_exp2 = lambda e, alpha, beta, n: alpha*(1-np.exp(-beta*(e**n)))


params_inv, _ = curve_fit(f_inv, energies, min_w0, p0=[6.5, 0.5], maxfev=2000)
params_exp, _ = curve_fit(f_exp, energies, min_w0, p0=[6.5, 0.5], maxfev=2000)
params_exp2, _ = curve_fit(f_exp2, energies, min_w0, p0=[6.5, 0.5, 0.5], maxfev=2000)

print(params_exp2)

print(params_inv)

plt.figure(figsize=(9,7))
plt.suptitle("w0_min value with the energy",
    fontsize=15,
    fontweight='bold',
    )
plt.title("w0_min -> w0 such as theta_res(E) is the smallest", style='italic', fontsize = 11, loc='left')
plt.plot(energies, min_w0, 'x', label='data')
# plt.plot(np.linspace(0, 90, 180), f_inv(np.linspace(0.1, 97, 180), *params_inv), label="invert.")
plt.plot(np.linspace(0, 97, 180), f_exp(np.linspace(0.1, 97, 180), *params_exp), label="exp(sqrt)")
plt.plot(np.linspace(0, 97, 180), f_exp2(np.linspace(0.1, 97, 180), *params_exp2), label="exp(sqrt)")
plt.xlabel('E part. (GeV)', loc='right')
plt.ylabel('min w')
plt.legend()
plt.show()
plt.savefig('plots/w0_min_energy_interpolation.png')
