"""
Abstract base class for CABARET mesh objects

Notes:

    1. In the future for unstructured grids, geometry data could contain more information
    2. Investigate class methods

Author: JWT
"""
from abc import ABC, abstractmethod


class Mesh(ABC):
    """
    @:brief Base class to represent a mesh object
    """
    class GeometryData:
        """
        @:brief Geometry data subclass
        """
        @property
        @abstractmethod
        def bcs_type(self):
            """
            @:brief Boundary conditions type
            """
            pass


    class SpaceTimeCell:
        """
        @:brief A space-time cell object to represent variables data
        """
        @property
        @abstractmethod
        def coordinates(self):
            """
            @:brief Cell centre coordinates [x_1, x_2, x_3]
            """
            pass

        @coordinates.setter
        @abstractmethod
        def coordinates(self, value):
            return None

        @property
        @abstractmethod
        def dx(self):
            """
            @:brief Cell size in the x_1 direction
            """
            pass

        @dx.setter
        @abstractmethod
        def dx(self, value):
            return None

        @property
        @abstractmethod
        def dy(self):
            """
            @:brief Cell size in the x_2 direction
            """
            pass

        @dy.setter
        @abstractmethod
        def dy(self, value):
            return None

        @property
        @abstractmethod
        def dz(self):
            """
            @:brief Cell size in the x_3 direction
            """
            pass

        @dz.setter
        @abstractmethod
        def dz(self, value):
            return None

        @property
        @abstractmethod
        def dt(self):
            """
            @:brief Local time-step size
            """
            pass

        @dt.setter
        @abstractmethod
        def dt(self, value):
            return None

        @property
        @abstractmethod
        def centre(self):
            """
            @:brief Variable values at cell centre
            """
            pass

        @centre.setter
        @abstractmethod
        def centre(self, value):
            return None

        @property
        @abstractmethod
        def faces(self):
            """
            @:brief Variable values at cell faces for flux values
            """
            pass

        @faces.setter
        @abstractmethod
        def faces(self, value):
            return None

        @property
        @abstractmethod
        def faces_new(self):
            """
            @:brief Variable values at cell faces at time-level n + 1
            """
            pass

        @faces_new.setter
        @abstractmethod
        def faces_new(self, value):
            return None

        @property
        @abstractmethod
        def centre_half(self):
            """
            @:brief Variable values at cell centre and time-level n + 1/2
            """
            pass

        @centre_half.setter
        @abstractmethod
        def centre_half(self, value):
            return None

        @abstractmethod
        def predictor(self):
            """
            @:brief Predictor step to get cell centred values at time-level n + 1/2
            """
            pass

        @abstractmethod
        def second_order_extrapolation(self):
            """
            @:brief Extrapolate values from t = n and n + 1/2 to get fluxes at t = n + 1
            """
            pass

        @abstractmethod
        def corrector(self):
            """
            @:brief Corrector step to get cell centred values at t = n + 1
            """
            pass


    @property
    @abstractmethod
    def nx(self):
        """
        @:brief Number of cells in the x_1 direction
        """
        pass

    @property
    @abstractmethod
    def ny(self):
        """
        @:brief Number of cells in the x_2 direction
        """
        pass

    @property
    @abstractmethod
    def nz(self):
        """
        @:brief Number of cells in the x_3 direction
        """
        pass

    @property
    @abstractmethod
    def geometry_data(self):
        """
        @:brief Geometry data
        """
        pass

    @property
    @abstractmethod
    def cells_data(self):
        """
        @:brief A collection of space-time cells
        """
        pass

    @abstractmethod
    def enforce_boundary_conditions(self):
        """
        @:brief Enforce boundary conditions
        :return:
        """
        pass
