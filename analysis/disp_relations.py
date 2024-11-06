"""
Plot dispersion relations of schemes

Author: JWT
"""
import numpy as np
import os
import matplotlib.pyplot as plt

path = os.path.dirname(__file__)

# Parameters
lw = 1.5
CFL = 0.4  # Courant number
k_dx = np.linspace(0., np.pi, 100)  # Normalised wavenumber k * Δx

# Compute the numerical phase speed for CTCS scheme
numerical_phase_speed_ctcs_1 = (1 / CFL) * np.arctan2(CFL * np.sin(k_dx), np.sqrt(1 - (CFL * np.sin(k_dx)) ** 2))
numerical_phase_speed_ctcs_2 = -(1 / CFL) * np.arctan2(CFL * np.sin(k_dx), np.sqrt(1 - (CFL * np.sin(k_dx)) ** 2))

# Compute the numerical phase speed for FTFS scheme
numerical_phase_speed_ftfs = (1 / CFL) * np.arctan2(CFL * np.sin(k_dx), (1 - CFL + CFL * np.cos(k_dx)))

# Compute the numerical phase speed for BTFS scheme
numerical_phase_speed_btfs = (1 / CFL) * np.arctan2((2 * CFL * np.sin(k_dx / 2)**2), (1 + 2 * CFL * np.sin(k_dx / 2) * np.cos(k_dx / 2)))

# Exact phase speed (should be 1 for comparison)
exact_phase_speed = k_dx  # np.ones_like(k_dx)

# Plotting
plt.rc('text', usetex=True)
plt.figure(figsize=(10, 6))
plt.plot(k_dx, numerical_phase_speed_ctcs_1, 'k--', label="CTCS", lw=lw)
plt.plot(k_dx, numerical_phase_speed_ctcs_2, 'k--', lw=lw)

plt.plot(k_dx, numerical_phase_speed_ftfs, 'b:', label="FTFS", lw=lw)
plt.plot(k_dx, numerical_phase_speed_btfs, 'r-.', label="BTFS", lw=lw)
plt.plot(k_dx, numerical_phase_speed_ftfs + numerical_phase_speed_btfs, 'gx--', label="FTFS + BTFS", lw=lw, ms=3)

plt.plot(k_dx, exact_phase_speed, 'k-', label="True")
plt.xlim([0., np.pi])
plt.ylim([-1.1, 3.1])
plt.xlabel(r'$k \Delta x$', usetex=True)
plt.xticks(np.linspace(0., np.pi, 5), ['0.', r'$\pi / 4$', r'$\pi / 2$', r'$3 \pi / 4$', r'$\pi$'], usetex=True)
plt.ylabel(r'$\omega_{n} / c$', usetex=True)
plt.title("Dispersion relations, numerical vs exact normalised angular frequencies")
plt.legend()
plt.tick_params(axis='both', direction='in')
plt.savefig(f'{path}/disp_relations.eps', dpi=1000)
plt.show()
