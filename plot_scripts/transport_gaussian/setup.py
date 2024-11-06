"""
Solver options and case setup

Author: JWT
"""
import numpy as np

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[-1.], [1.]]
NO_OF_GRID_POINTS = [65]
BOUNDARY_CONDITIONS_TYPE = ["Periodic", "Periodic"]
FINAL_TIME = 2.

"""======SOLVER OPTIONS======"""
CFL = 1.
SPATIAL_SCHEME = "Finite difference"
TEMPORAL_SCHEME = "Forward Euler"
FLUX_SCHEME = "Upwind"
EQUATION_TYPE = "Scalar transport"