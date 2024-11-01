"""
Solver options and case setup
"""
import numpy as np

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 2
DOMAIN_BOUNDARIES = [[0., 0.], [1., 0.], [1., 1.], [0, 1.]]
NO_OF_GRID_POINTS = [101, 101]
BOUNDARY_CONDITIONS_TYPE = ["Periodic", "Periodic", "Periodic", "Periodic"]
FINAL_TIME = 10.

"""======SOLVER OPTIONS======"""
CFL = 0.5
SPATIAL_SCHEME = "Finite difference"
TEMPORAL_SCHEME = "Forward Euler"
FLUX_SCHEME = "Upwind"
EQUATION_TYPE = "Scalar transport"
