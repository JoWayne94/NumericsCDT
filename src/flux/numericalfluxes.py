"""
Numerical flux definitions

Author: JWT
"""
import numpy as np


def upwind(values: np.array, f: np.array, dx=None, dt=None):

    # fluxes = np.empty(len(f) - 1)
    #
    # for i in range(len(f) - 1):
    #
    #     velocity = f[i]
    #     # velocity = 0.5 * (v[i] + v[i + 1])
    #     # velocity = max(v[i], v[i + 1] * -1.)
    #
    #     if velocity >= 0.:
    #         fluxes[i] = velocity * values[i]
    #     elif velocity < 0.:
    #         fluxes[i] = velocity * values[i + 1]
    #     else:
    #         raise ValueError('Flow velocity undefined.')

    # fluxes_l = np.where(values[:-1] >= 0., f[:-1], 0.)
    # fluxes_r = np.where(values[:-1] < 0., f[1:], 0.)
    fluxes_l = np.where(0.5 * (values[:-1] + values[1:]) >= 0., f[:-1], 0.)
    fluxes_r = np.where(0.5 * (values[:-1] + values[1:]) < 0., f[1:], 0.)

    return fluxes_l + fluxes_r


def lax_friedrichs(values: np.array, f: np.array, dx=None, dt=None):

    # for i in range(n + 1):
    #
    #     fluxes[i] = 0.5 * (dx / dt) * (values[i] - values[i + 1]) + \
    #                 0.5 * (f[i] + f[i + 1])

    return 0.5 * (dx / dt) * (values[:-1] - values[1:]) + 0.5 * (f[:-1] + f[1:])


def richtmyer(values: np.array, f: np.array, dx=None, dt=None):

    # for i in range(n + 1):
    #
    #     # velocity = 0.5 * (values[i] + values[i + 1])
    #
    #     u_half = 0.5 * (values[i] + values[i + 1]) - \
    #              0.5 * (dt / dx) * (f[i + 1] - f[i])
    #
    #     fluxes[i] = 0.5 * u_half * u_half

    u_half = 0.5 * (values[:-1] + values[1:]) - 0.5 * (dt / dx) * (f[1:] - f[:-1])

    return 0.5 * u_half * u_half


def force(values: np.array, f: np.array, dx=None, dt=None):

    fluxes_lf = lax_friedrichs(values, f, dx, dt)
    fluxes_ri = richtmyer(values, f, dx, dt)

    return 0.5 * (fluxes_lf + fluxes_ri)
