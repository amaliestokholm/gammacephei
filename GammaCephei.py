# coding: utf-8

# # Gamma Cephei
# Gamma Cephei is a spectroscopic binary with a planet of minimum mass of 1.7 sin(i) jupiter mass.
# Gamma Cephei A is a K1 subgiant star, while Gamma Cephei B is a red dwarf.

import numpy as np

# Stellar parameters found in the literature
pi = 70.91  # mas. From Hipparcos, new estimates

# Colors from SIMBAD
U = 5.190
B = 4.250
V = 3.220
R = 2.6
I = 2.1
J = 1.66
H = 1.19
K = 1.04

# Metallicities from the literature: here the lowest, median, and highest value since 1990.
FeHs = [-0.05, 0.04, 0.17]

BC = -0.24  # Estimate from Kaler 1989
M = 1  # in solar masses

# Solar values from Prsa et al 2016
teff_sun = 5772
logg_sun = 4.43

# Distance
d = 1 / (pi * 10 ** (-3))
print('The distance is %.3f pc' % d)

# Effective temperature
# We assume a metallicity similar to the Sun and use Alonso et al. 1999
def alonso1999(X, FeH, a):
    theta = (a[0] + a[1] * X
             + a[2] * X ** 2
             + a[3] * X * FeH
             + a[4] * FeH
             + a[5] * FeH ** 2)
    teff = 5040 / theta
    teffs.append(teff)
    return teff


def print_teff(color, teff, sigma):
    print('Effective temperature from %s is %.0f +/- %.0f K'
        % (color, teff, sigma))

for FeH in FeHs:
    teffs = []
    # Using the relation with the most stars in the sample
    print('If [Fe/H] = %.2f:' % FeH)
    if U and V is not None:
        #a_uv = [0.6388, 0.4065, -0.1117, -2.308e-3, -7.783e-2, 1.200e-2]
        a_uv = [0.8323, 9.374e-2, 1.184e-2, 2.351e-2, -0.1392, -1.944e-2]
        sigma_uv = 80
        teff_uv = alonso1999(U - V, FeH, a_uv)
        print_teff('U-V', teff_uv, sigma_uv)

    if B and V is not None:
        a_bv = [0.6177, 0.4354, -4.025e-3, 5.204e-2, -0.1127, -1.385e-2]
        sigma_bv = 96
        teff_bv = alonso1999(B - V, FeH, a_bv)
        print_teff('B-V', teff_bv, sigma_bv)

    if V and R is not None:
        a_vr = [0.4972, 0.8841, -0.1904, -1.197e-2, -1.025e-2, -5.500e-3]
        sigma_vr = 150
        teff_vr = alonso1999(V - R, FeH, a_vr)
        print_teff('V-R', teff_vr, sigma_vr)

    if V and K is not None:
        a_vk = [0.5558, 0.2105, 1.981e-3, -9.965e-3, 1.325e-2, -2.726e-3]
        sigma_vk = 25
        teff_vk = alonso1999(V - K, FeH, a_vk)
        print_teff('V-K', teff_vk, sigma_vk)

    if J and H is not None:
        a_jh = [0.5977, 1.015, -1.020e-1, -1.029e-2, 3.006e-2, 1.013e-2]
        sigma_jh = 170
        teff_jh = alonso1999(J - H, FeH, a_jh)
        print_teff('J-H', teff_jh, sigma_jh)

    if J and K is not None:
        a_jh = [0.5816, 0.9134, -0.1443, 0, 0, 0]
        sigma_jh = 125
        teff_jh = alonso1999(J - H, FeH, a_jh)
        print_teff('J-H', teff_jh, sigma_jh)

        
    teffs = np.asarray(teffs)
    mean = np.mean(teffs)
    median = np.median(teffs)
    std = np.std(teffs)
    print('The mean of all calculated effective temperatures is %.0f K' % mean)
    print('The median of all calculated effective temperatures is %.0f K' % median)
    print('The std of all calculated effective temperatures is %.0f K' % std)
    print('')

# Surface gravity
logg_ratio = np.log10(M) + 4 * np.log10(median / teff_sun) + 0.4 * V + 0.4 * BC + 2 * np.log10(pi) + 0.12
logg = logg_ratio * logg_sun
print('The surface gravity is %.1f' % logg)
