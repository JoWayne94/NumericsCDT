"""
Solver options and case setup

Author: JWT
"""
import numpy as np


"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[0.], [2. * np.pi]]
NO_OF_CELLS = [100]
BOUNDARY_CONDITIONS_TYPE = ["Periodic", "Periodic"]
FINAL_TIME = 2. * np.pi

"""======SOLVER OPTIONS======"""
CFL = 1.
EQUATION_TYPE = "Scalar transport"
NO_OF_VARIABLES = 1
