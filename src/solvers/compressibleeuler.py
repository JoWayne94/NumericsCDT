"""
Inviscid, compressible Euler equations set

\partial U / \partial t + div(F) = 0

where u is

U = [rho, rho * u, rho * E]^T,

and the flux, F, is

F = [rho * u, rho * u^2 + p, rho * u * E + u * p]^T.

Author: JWT
"""
import numpy as np


"""
Macros
"""
gamma = 1.4 # 5. / 3.
gamma2 = gamma - 1
mu = gamma2 / (2. * gamma)
eToPres = lambda rho, e: (gamma - 1.) * rho * e
presToe = lambda rho, p: p / ((gamma - 1.) * rho)
eToE = lambda rho, u, p: presToe(rho, p) + u * u / 2.


def speed_of_sound(p, rho):
    """
    Ideal EOS, calculate the speeds of sound in the system
    """

    return np.sqrt(gamma * p / rho)

def prim_to_cons(rho, u, p):

    return rho, rho * u, rho * eToE(rho, u, p)

def cons_to_prim(rho, rhou, rhoE):

    u = rhou / rho
    e = (rhoE / rho) - u * u / 2.

    return rho, u, eToPres(rho, e)

def prim_to_flux(rho, u, p):

    return rho * u, rho * u * u + p, rho * u * eToE(rho, u, p) + u * p

def cons_to_flux(rho, rhou, rhoE):

    Vtmp1, Vtmp2, Vtmp3 = cons_to_prim(rho, rhou, rhoE)
    Ftmp1, Ftmp2, Ftmp3 = prim_to_flux(Vtmp1, Vtmp2, Vtmp3)

    return Ftmp1, Ftmp2, Ftmp3


def wave_estimates(u_1l, u_2l, u_3l, u_1r, u_2r, u_3r):

    rho_l, u_l, p_l = cons_to_prim(u_1l, u_2l, u_3l)
    rho_r, u_r, p_r = cons_to_prim(u_1r, u_2r, u_3r)

    """ Ideal Gas EoS """
    a_l = speed_of_sound(p_l, rho_l)
    a_r = speed_of_sound(p_r, rho_r)

    return (np.array([min(u_l[i], u_r[i]) - max(a_l[i], a_r[i]) for i in range(len(u_l))]),
            np.array([max(u_l[i], u_r[i]) + max(a_l[i], a_r[i]) for i in range(len(u_l))]))


def hll(u_1l, u_2l, u_3l, u_1r, u_2r, u_3r):

    s_l, s_r = wave_estimates(u_1l, u_2l, u_3l, u_1r, u_2r, u_3r)

    f_1l, f_2l, f_3l  = cons_to_flux(u_1l, u_2l, u_3l)
    f_1r, f_2r, f_3r = cons_to_flux(u_1r, u_2r, u_3r)

    f_1hll = np.empty_like(f_1l)
    f_2hll = np.empty_like(f_2l)
    f_3hll = np.empty_like(f_3l)

    for i in range(len(s_l)):

        if s_l[i] >= 0.:
            f_1hll[i] = f_1l[i]
            f_2hll[i] = f_2l[i]
            f_3hll[i] = f_3l[i]

        elif s_r[i] <= 0.:
            f_1hll[i] = f_1r[i]
            f_2hll[i] = f_2r[i]
            f_3hll[i] = f_3r[i]

        else:
            f_1hll[i] = (s_r[i] * f_1l[i] - s_l[i] * f_1r[i] + s_l[i] * s_r[i] * (u_1r[i] - u_1l[i])) / (s_r[i] - s_l[i])
            f_2hll[i] = (s_r[i] * f_2l[i] - s_l[i] * f_2r[i] + s_l[i] * s_r[i] * (u_2r[i] - u_2l[i])) / (s_r[i] - s_l[i])
            f_3hll[i] = (s_r[i] * f_3l[i] - s_l[i] * f_3r[i] + s_l[i] * s_r[i] * (u_3r[i] - u_3l[i])) / (s_r[i] - s_l[i])

    return f_1hll, f_2hll, f_3hll


def compressible_euler_1d(dt, velocity_vector, numerical_flux, u_1, u_2=None, u_3=None):
    """
    @:brief Inviscid, compressible Euler flux computations

    :param dt:              Time-step size
    :param u_1:             Grid density object
    :param u_2:             Grid density * velocity object in the x_1 direction
    :param u_3:             Grid density * energy object
    :param numerical_flux:  Numerical flux of choice
    :param velocity_vector: Constant in time driving velocity field
    :return: 0.5 * u * u
    """

    Ftmp1, Ftmp2, Ftmp3 = cons_to_flux(u_1.variables_data.values, u_2.variables_data.values, u_3.variables_data.values)

    # u_1.variables_data.x_fluxes = numerical_flux(u_1.variables_data.x_fluxes,
    #                                               u_1.variables_data.values, u_1.nx,
    #                                               Ftmp1, u_1.dx, dt)
    #
    # u_2.variables_data.x_fluxes = numerical_flux(u_2.variables_data.x_fluxes,
    #                                              u_2.variables_data.values, u_2.nx,
    #                                              Ftmp2, u_2.dx, dt)
    #
    # u_3.variables_data.x_fluxes = numerical_flux(u_3.variables_data.x_fluxes,
    #                                              u_3.variables_data.values, u_3.nx,
    #                                              Ftmp3, u_3.dx, dt)

    (u_1.variables_data.x_fluxes,
     u_2.variables_data.x_fluxes,
     u_3.variables_data.x_fluxes) = hll(u_1.variables_data.values[:-1], u_2.variables_data.values[:-1],
                                        u_3.variables_data.values[:-1], u_1.variables_data.values[1:],
                                        u_2.variables_data.values[1:], u_3.variables_data.values[1:])
