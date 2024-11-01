"""
CABARET implementations

Author: JWT
"""
from src.CABARET.mesh.mesh1D import Mesh1D


def cabaret_1d(msh: Mesh1D):

    for cell in range(1, msh.nx + 1):

        msh.cells_data[cell].predictor()
        msh.cells_data[cell].second_order_extrapolation()
        msh.cells_data[cell].nonlinear_limiter()

    msh.flux_reconstruction()

    for cell in range(1, msh.nx + 1):

        msh.cells_data[cell].corrector()
        msh.cells_data[cell].faces = msh.cells_data[cell].faces_new.copy()

