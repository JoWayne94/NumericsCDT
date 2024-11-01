"""
Grid subclass in one-dimensional space.

Author: JWT
"""
import numpy as np
from src.grid.grid import Grid


class Grid1D(Grid):
    """
    @:brief Contains data of a grid (variable) object
    """
    class GridPointsData:
        """
        @:brief Points data subclass
        """
        def __init__(self):
            self._x_coordinates = None

        @property
        def x_coordinates(self):
            return self._x_coordinates

        @x_coordinates.setter
        def x_coordinates(self, value):
            self._x_coordinates = value

    class GridVariableData:
        """
        @:brief Variables data subclass
        """
        def __init__(self, no_of_points, bcs_type):
            self._values = np.empty(no_of_points[0] + 2)
            self._x_fluxes = np.empty(no_of_points[0] + 1)
            self._x_first_derivative = np.empty(no_of_points[0])
            self._bcs_type = bcs_type

        @property
        def values(self):
            return self._values

        @values.setter
        def values(self, value):
            # self._values[:len(value)] = value
            self._values = value

        @property
        def x_fluxes(self):
            return self._x_fluxes

        @x_fluxes.setter
        def x_fluxes(self, value):
            self._x_fluxes = value

        @property
        def x_first_derivative(self):
            return self._x_first_derivative

        @x_first_derivative.setter
        def x_first_derivative(self, value):
            self._x_first_derivative = value

        @property
        def y_first_derivative(self):
            return 0.

        def __add__(self, other):
            self._values[1:-1] += other

        def enforce_boundary_conditions(self):

            # Left boundary
            if self._bcs_type[0] == "Periodic":

                self._values[0] = self._values[-3]
                # self._fluxes[0] = self._values[-2]

            elif self._bcs_type[0] == "Zero Neumann":

                self._values[0] = self._values[1]
                # self._fluxes[0] = self._fluxes[1]

            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[0]} not implemented")

            # Right boundary
            if self._bcs_type[1] == "Periodic":

                self._values[-1] = self._values[2]
                # self._values[-1] = self._values[1]

            elif self._bcs_type[1] == "Zero Neumann":

                self._values[-1] = self._values[-2]
                # self._fluxes[-1] = self._fluxes[-2]

            else:
                raise NotImplementedError(f"Type of BC {self._bcs_type[1]} not implemented")


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
        self._points_data.coordinates = np.linspace(boundary_coords[0][0] - self._dx, boundary_coords[1][0] + self._dx, ngp[0] + 2)

    @property
    def nx(self):
        return self._ngp[0]

    @property
    def ny(self):
        return 0

    @property
    def nz(self):
        return 0

    @property
    def dx(self):
        return self._dx

    @property
    def dy(self):
        return 0.

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
        self.variables_data.x_fluxes = np.zeros(self.nx + 1)
