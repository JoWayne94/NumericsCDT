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
FINAL_TIME = 1.

"""======SOLVER OPTIONS======"""
CFL = 1.0
EQUATION_TYPE = "Inviscid Burgers"
NO_OF_VARIABLES = 1
