"""
Grid subclass in two-dimensional space.

Author: JWT
"""
import numpy as np
from src.grid.grid import Grid


class Grid2D(Grid):
    """
    @:brief Contains data of a grid (variable) object
    """
    class GridPointsData:
        """
        @:brief Points data subclass
        """
        def __init__(self):
            self._x_coordinates = None
            self._y_coordinates = None

        @property
        def x_coordinates(self):
            return self._x_coordinates

        @x_coordinates.setter
        def x_coordinates(self, value):
            self._x_coordinates = value

        @property
        def y_coordinates(self):
            return self._y_coordinates

        @y_coordinates.setter
        def y_coordinates(self, value):
            self._y_coordinates = value

    class GridVariableData:
        """
        @:brief Variables data subclass
        """
        def __init__(self, no_of_points, bcs_type):
            self._values = np.empty((no_of_points[1] + 2, no_of_points[0] + 2))
            self._x_fluxes = np.empty((no_of_points[1], no_of_points[0] + 1))
            self._y_fluxes = np.empty((no_of_points[1] + 1, no_of_points[0]))
            self._x_first_derivative = np.empty((no_of_points[1], no_of_points[0]))
            self._y_first_derivative = np.empty((no_of_points[1], no_of_points[0]))
            self._bcs_type = bcs_type

        @property
        def values(self):
            return self._values

        @values.setter
        def values(self, value):
            self._values = value

        @property
        def x_fluxes(self):
            return self._x_fluxes

        @x_fluxes.setter
        def x_fluxes(self, value):
            self._x_fluxes = value

        @property
        def y_fluxes(self):
            return self._y_fluxes

        @y_fluxes.setter
        def y_fluxes(self, value):
            self._y_fluxes = value

        @property
        def x_first_derivative(self):
            return self._x_first_derivative

        @x_first_derivative.setter
        def x_first_derivative(self, value):
            self._x_first_derivative = value

        @property
        def y_first_derivative(self):
            return self._y_first_derivative

        @y_first_derivative.setter
        def y_first_derivative(self, value):
            self._y_first_derivative = value

        def __add__(self, other):
            self._values[1:-1, 1:-1] += other

        def enforce_boundary_conditions(self):

            # Left most domain boundary
            if self._bcs_type[0] == "Periodic":
                self._values[:, 0] = self._values[:, -3]
            elif self._bcs_type[0] == "Zero Neumann":
                self._values[:, 0] = self._values[:, 1]
            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[0]} for left boundary not implemented")

            # Bottom domain boundary
            if self._bcs_type[1] == "Periodic":
                self._values[-1, :] = self._values[2, :]
            elif self._bcs_type[1] == "Zero Neumann":
                self._values[-1, :] = self._values[-2, :]
            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[1]} for bottom boundary not implemented")

            # Top domain boundary
            if self._bcs_type[2] == "Periodic":
                self._values[0, :] = self._values[-3, :]
            elif self._bcs_type[2] == "Zero Neumann":
                self._values[0, :] = self._values[1, :]
            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[2]} for right boundary not implemented")

            # Right most domain boundary
            if self._bcs_type[3] == "Periodic":
                self._values[:, -1] = self._values[:, 2]
            elif self._bcs_type[3] == "Zero Neumann":
                self._values[:, -1] = self._values[:, -2]
            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[3]} for top boundary not implemented")


    def __init__(self, ngp, boundary_coords, bcs_type):
        """
        @:brief Main constructor

        :param ngp: Number of grid points
        :type ngp: list[nd]
        :param boundary_coords: Physical boundary coordinates
        :type boundary_coords: list[nd * 2][nd]
        :param bcs_type: Type of boundary conditions
        :type bcs_type: str
        """
        self._points_data = self.GridPointsData()
        self._variables_data = self.GridVariableData(ngp, bcs_type)
        self._ngp = ngp

        # Points data
        self._dx = (boundary_coords[1][0] - boundary_coords[0][0]) / (ngp[0] - 1)
        self._points_data.x_coordinates = np.linspace(boundary_coords[0][0] - self._dx,
                                                      boundary_coords[1][0] + self._dx, ngp[0] + 2)

        self._dy = (boundary_coords[3][1] - boundary_coords[0][1]) / (ngp[1] - 1)
        self._points_data.y_coordinates = np.linspace(boundary_coords[0][1] - self._dy,
                                                      boundary_coords[3][1] + self._dy, ngp[1] + 2)

        self._n = [self.nx, self.ny]

    @property
    def nx(self):
        return self._ngp[0]

    @property
    def ny(self):
        return self._ngp[1]

    @property
    def nz(self):
        return 0

    @property
    def dx(self):
        return self._dx

    @property
    def dy(self):
        return self._dy

    @property
    def dz(self):
        return 0.

    @property
    def points_data(self):
        return self._points_data

    @property
    def variables_data(self):
        return self._variables_data

    def reset_fluxes(self):
        self.variables_data.x_fluxes = np.zeros((self.ny, self.nx + 1))
        self.variables_data.y_fluxes = np.zeros((self.ny + 1, self.nx))
