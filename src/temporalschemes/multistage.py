"""
Time stepping methods for the semi-discrete form

Note:

    1. for 1d, set y first derivatives as zeros

Author: JWT
"""


def forward_euler(grid, dt):

    grid.variables_data.__add__(dt * grid.variables_data.x_first_derivative)
    grid.variables_data.__add__(dt * grid.variables_data.y_first_derivative)
    grid.variables_data.enforce_boundary_conditions()
