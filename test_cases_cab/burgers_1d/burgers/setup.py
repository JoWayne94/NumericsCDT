"""
Solver options and case setup

Author: JWT
"""

"""======SPACE-TIME DOMAIN======"""
NO_OF_DIMENSIONS = 1
DOMAIN_BOUNDARIES = [[-1.], [2.]]
NO_OF_CELLS = [100]
BOUNDARY_CONDITIONS_TYPE = ["Zero Neumann", "Zero Neumann"]
FINAL_TIME = 1.

"""======SOLVER OPTIONS======"""
CFL = 0.5
EQUATION_TYPE = "Inviscid Burgers"
NO_OF_VARIABLES = 1
