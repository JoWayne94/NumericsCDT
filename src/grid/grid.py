"""
Abstract base class for grid objects

Author: JWT
"""
from abc import ABC, abstractmethod


class Grid(ABC):
    """
    @:brief Base class to represent a grid (variable) object in physical space
    """
    class GridPointsData:
        """
        @:brief Points data subclass
        """
        @property
        @abstractmethod
        def x_coordinates(self):
            """
            @:brief Point coordinates in first principle direction
            """
            pass

        @x_coordinates.setter
        @abstractmethod
        def x_coordinates(self, value):
            return None

    class GridVariableData:
        """
        @:brief Variables data subclass
        """
        @property
        @abstractmethod
        def values(self):
            """
            @:brief Variable values
            """
            pass

        @values.setter
        @abstractmethod
        def values(self, value):
            return None

        @property
        @abstractmethod
        def x_fluxes(self):
            """
            @:brief Flux values in first principle direction
            """
            pass

        @x_fluxes.setter
        @abstractmethod
        def x_fluxes(self, value):
            return None

        @property
        @abstractmethod
        def x_first_derivative(self):
            """
            @:brief First-order derivative in space for the x_1 direction
            """
            pass

        @x_first_derivative.setter
        @abstractmethod
        def x_first_derivative(self, value):
            return None

        @abstractmethod
        def __add__(self, other):
            pass

        @abstractmethod
        def enforce_boundary_conditions(self):
            """
            @:brief Enforce boundary conditions
            :return: void
            """
            pass

    @property
    @abstractmethod
    def nx(self):
        """
        @:brief Number of grid points in the x_1 direction
        """
        pass

    @property
    @abstractmethod
    def ny(self):
        """
        @:brief Number of grid points in the x_2 direction
        """
        pass

    @property
    @abstractmethod
    def nz(self):
        """
        @:brief Number of grid points in the x_3 direction
        """
        pass

    @property
    @abstractmethod
    def dx(self):
        """
        @:brief Grid spacing in the x_1 direction
        """
        pass

    @property
    @abstractmethod
    def dy(self):
        """
        @:brief Grid spacing in the x_2 direction
        """
        pass

    @property
    @abstractmethod
    def dz(self):
        """
        @:brief Grid spacing in the x_3 direction
        """
        pass

    @property
    @abstractmethod
    def points_data(self):
        """
        @:brief Grid points data
        """
        pass

    @property
    @abstractmethod
    def variables_data(self):
        """
        @:brief Variables data
        """
        pass

    @abstractmethod
    def reset_fluxes(self):
        """
        @:brief Zero all flux values in numpy array
        """
        pass
