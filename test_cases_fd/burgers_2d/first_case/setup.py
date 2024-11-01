"""
Solver options and case setup
"""
import numpy as np

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 2
DOMAIN_BOUNDARIES = [[-6., -1.5], [6., -1.5], [6., 1.5], [-6, 1.5]]
NO_OF_GRID_POINTS = [361, 91]
BOUNDARY_CONDITIONS_TYPE = ["Zero Neumann", "Zero Neumann", "Zero Neumann", "Zero Neumann"]
FINAL_TIME = 2.5

"""======SOLVER OPTIONS======"""
CFL = 0.5
SPATIAL_SCHEME = "Finite difference"
TEMPORAL_SCHEME = "Forward Euler"
FLUX_SCHEME = "Upwind"
EQUATION_TYPE = "Inviscid Burgers"
