"""
Solver options and case setup

Author: JWT
"""


"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[0.], [2.]]
NO_OF_CELLS = [100]
BOUNDARY_CONDITIONS_TYPE = ["Periodic", "Periodic"]
FINAL_TIME = 2.

"""======SOLVER OPTIONS======"""
CFL = 0.45
EQUATION_TYPE = "Scalar transport"
NO_OF_VARIABLES = 1
