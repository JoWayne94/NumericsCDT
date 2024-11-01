"""
Solver options and case setup

Author: JWT
"""
import numpy as np

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[0.], [2. * np.pi]]
NO_OF_GRID_POINTS = [101]
BOUNDARY_CONDITIONS_TYPE = ["Zero Neumann", "Zero Neumann"]
FINAL_TIME = 1.

"""======SOLVER OPTIONS======"""
CFL = 1.
SPATIAL_SCHEME = "Finite difference"
TEMPORAL_SCHEME = "Forward Euler"
FLUX_SCHEME = "Upwind"
EQUATION_TYPE = "Inviscid Burgers"
