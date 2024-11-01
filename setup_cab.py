"""
Solver options and case setup

Author: JWT
"""
import numpy as np

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[-1.], [1.]]
NO_OF_CELLS = [400]
BOUNDARY_CONDITIONS_TYPE = ["Periodic", "Periodic"]
FINAL_TIME = 2.

"""======SOLVER OPTIONS======"""
CFL = 0.9
EQUATION_TYPE = "Scalar transport"
NO_OF_VARIABLES = 1
