import sys
sys.path.append(r'C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Python Code Lib')
import numpy as np
import pandas as pd
import re
import numpy as np
from matplotlib import pyplot as plt 
import your_dad
import scipy
from scipy.interpolate import griddata
import matplotlib.patches as patches
from scipy.optimize import curve_fit
import os

file_path = r'measured_refractive_index\ICP-PECVD Long Anneal.txt'
title = "ICP-PECVD"
output_file = r'temp.txt'
plotly = 0
output=0
temp=0
if title is None:
    title = os.path.basename(file_path)
data = np.loadtxt(file_path)
wave = data[:, 0] # Wavelength data
n = data[:, 1]     # Refractive index data
xi=300
xf=450
print(wave[xi],wave[xf])
fitted_wavelength, fitted_n, _ = your_dad.fit_Sellmeier(n, wave,index_begin=xi,index_end=xf,step=0.5)
plt.figure(dpi=150, figsize=(8, 6))
plt.plot(wave, n, 'bo', markersize=3, label=f'Measured Data_{title}')
plt.plot(fitted_wavelength, fitted_n, 'r-', linewidth=2, label=f'Fitted Sellmeier Curve_{title}')
if temp:
    file_path = r'Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt'
    data = np.loadtxt(file_path)
    wave = data[:, 0] # Wavelength data
    n = data[:, 1]     # Refractive index data
    plt.plot(wave, n, 'go', markersize=3, label='2%')

    


plt.xlabel('Wavelength (nm)', fontsize=12)
plt.ylabel('Refractive Index $n$', fontsize=12)
plt.title(title, fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

if plotly:
    your_dad.plotly_plot(fitted_wavelength, fitted_n, wave, n, x_label="Wavelength (nm)", y_label="Refractive Index $n$", title="Refractive Index Plot")

if output:
    with open(output_file, 'w') as f:
        for wavelength, refractive_index in zip(fitted_wavelength, fitted_n):
            f.write(f"{wavelength:.9f}\t{refractive_index:.9f}\n")