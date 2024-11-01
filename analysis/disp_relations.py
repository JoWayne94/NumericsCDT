"""
Plot dispersion relations of schemes

Author: JWT
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
lw = 1.
CFL = 0.4  # Courant number
k_dx = np.linspace(0., np.pi, 100)  # Normalized wavenumber k * Δx

# Compute the numerical phase speed for CTCS scheme
numerical_phase_speed_ctcs_1 = (1 / CFL) * np.arctan2(CFL * np.sin(k_dx), np.sqrt(1 - (CFL * np.sin(k_dx)) ** 2))
numerical_phase_speed_ctcs_2 = -(1 / CFL) * np.arctan2(CFL * np.sin(k_dx), np.sqrt(1 - (CFL * np.sin(k_dx)) ** 2))

# Compute the numerical phase speed for FTFS scheme
numerical_phase_speed_ftfs_1 = (1 / CFL) * np.arctan2(CFL * np.sin(k_dx), (1 - CFL + CFL * np.cos(k_dx)))


# Exact phase speed (should be 1 for comparison)
exact_phase_speed = k_dx  # np.ones_like(k_dx)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(k_dx, numerical_phase_speed_ctcs_1, 'k-', label="CTCS", lw=lw)
plt.plot(k_dx, numerical_phase_speed_ctcs_2, 'k-', lw=lw)

plt.plot(k_dx, numerical_phase_speed_ftfs_1, 'b-.', label="FTFS", lw=lw)

plt.plot(k_dx, exact_phase_speed, 'k--', label="True")
plt.xlim([0., np.pi])
plt.ylim([-2., 2.])
plt.xlabel(r'$\mathbf{k \Delta x}$', usetex=True)
plt.xticks(np.linspace(0., np.pi, 5), ['0.', r'$\pi / 4$', r'$\pi / 2$', r'$3 \pi / 4$', r'$\pi$'], usetex=True)
plt.ylabel(r'$\mathbf{\omega_{n} / c}$', usetex=True)
plt.title("Dispersion relation, numerical vs exact phase speeds")
plt.legend()
plt.axhline(0, linestyle=':', color='black')
# plt.grid()
plt.tick_params(axis='both', direction='in')
plt.show()
