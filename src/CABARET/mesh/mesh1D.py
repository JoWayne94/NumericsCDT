"""
Mesh subclass in one-dimensional space.

Author: JWT
"""
import numpy as np
from src.CABARET.mesh.mesh import Mesh


class Mesh1D(Mesh):
    """
    @:brief Contains data of a mesh (variable) object
    """
    class GeometryData:
        """
        @:brief Geometry data subclass
        """
        def __init__(self, bcs_type):
            self._bcs_type = bcs_type

        @property
        def bcs_type(self):
            return self._bcs_type


    class SpaceTimeCell(object):
        """
        @:brief A space-time cell object to represent variables data
        """
        def __init__(self, no_of_variables, flux_def):
            self._coordinates = None
            self._dx = None
            self._dy = None
            self._dz = None
            self._dt = None

            self._centre = np.empty(no_of_variables)
            self._c_centre = np.empty(no_of_variables)
            self._centre_half = np.empty(no_of_variables)

            self._faces = np.empty((2, no_of_variables))
            self._c_faces = np.empty((2, no_of_variables))
            self._faces_new = np.empty((2, no_of_variables))

            self._flux_def = flux_def

        @property
        def coordinates(self):
            return self._coordinates

        @coordinates.setter
        def coordinates(self, value):
            self._coordinates = value

        @property
        def dx(self):
            return self._dx

        @dx.setter
        def dx(self, value):
            self._dx = value

        @property
        def dy(self):
            return 0.

        @dy.setter
        def dy(self, value):
            self._dy = value

        @property
        def dz(self):
            return 0.

        @dz.setter
        def dz(self, value):
            self._dz = value

        @property
        def dt(self):
            return self._dt

        @dt.setter
        def dt(self, value):
            self._dt = value

        @property
        def centre(self):
            return self._centre

        @centre.setter
        def centre(self, value: np.array):
            self._centre = value

        @property
        def c_centre(self):
            return self._c_centre

        @c_centre.setter
        def c_centre(self, value: np.array):
            self._c_centre = value

        @property
        def faces(self):
            return self._faces

        @faces.setter
        def faces(self, value: np.array):
            self._faces = value

        @property
        def c_faces(self):
            return self._c_faces

        @c_faces.setter
        def c_faces(self, value: np.array):
            self._c_faces = value

        @property
        def faces_new(self):
            return self._faces_new

        @faces_new.setter
        def faces_new(self, value: np.array):
            self._faces_new = value

        @property
        def centre_half(self):
            return self._centre_half

        @centre_half.setter
        def centre_half(self, value: np.array):
            self._centre_half = value

        def predictor(self):

            self.centre_half = self.centre - 0.5 * self.dt * (self._flux_def(self.c_faces[1, :], self.faces[1, :])
                                                            - self._flux_def(self.c_faces[0, :], self.faces[0, :])) / self.dx

        def second_order_extrapolation(self):

            self.faces_new[0, :] = 2. * self.centre_half - self.faces[1, :]
            self.faces_new[1, :] = 2. * self.centre_half - self.faces[0, :]

        def nonlinear_limiter(self):

            _M = np.maximum.reduce([self.faces[0, :], self.faces[1, :], self.centre_half])
            _m = np.minimum.reduce([self.faces[0, :], self.faces[1, :], self.centre_half])

            for face in range(self.faces.shape[0]):

                self.faces_new[face, :] = np.minimum.reduce([self.faces_new[face, :], _M])
                self.faces_new[face, :] = np.maximum.reduce([self.faces_new[face, :], _m])

        def corrector(self):

            self.centre = self.centre_half - 0.5 * self.dt * (self._flux_def(self.c_faces[1, :], self.faces_new[1, :])
                                                            - self._flux_def(self.c_faces[0, :], self.faces_new[0, :])) / self.dx


    def __init__(self, nc, nv, boundary_coords, bcs_type, flux_def):
        """
        @:brief Main constructor

        :param nc: Number of cells
        :type nc: list[nd]
        :param nv: Number of variables
        :type nv: int
        :param boundary_coords: Physical boundary coordinates
        :type boundary_coords: list[nd * 2][nd]
        :param bcs_type: Type of boundary conditions
        :type bcs_type: str
        :param flux_def: Interface flux definition
        :type flux_def: method(., .)
        """
        self._nc = nc
        self._geometry_data = self.GeometryData(bcs_type)
        self._cells_data = np.ndarray((nc[0] + 2,), dtype=self.SpaceTimeCell)
        dx = (boundary_coords[1][0] - boundary_coords[0][0]) / nc[0]

        self._x_coordinates = np.linspace(boundary_coords[0][0] - 0.5 * dx, boundary_coords[1][0] + 0.5 * dx, nc[0] + 2)

        # Wave speeds at cell faces at half time-steps: np.array((nv, nc[0] + 1,))
        self._c = None

        # Cells data
        for cell in range(nc[0] + 2):

            self._cells_data[cell] = self.SpaceTimeCell(nv, flux_def)
            self._cells_data[cell].coordinates = [self._x_coordinates[cell]]
            self._cells_data[cell].dx = dx


    @property
    def nx(self):
        return self._nc[0]

    @property
    def ny(self):
        return 0

    @property
    def nz(self):
        return 0

    @property
    def geometry_data(self):
        return self._geometry_data

    @property
    def cells_data(self):
        return self._cells_data

    @property
    def x_coordinates(self):
        return self._x_coordinates

    @property
    def cell_centre_values(self):
        return np.array([self.cells_data[cell].centre for cell in range(self.nx + 2)])

    @property
    def c_centre_values(self):
        return np.array([self.cells_data[cell].c_centre for cell in range(self.nx + 2)])

    def enforce_boundary_conditions(self):

        # Left boundary
        if self.geometry_data.bcs_type[0] == "Periodic":

            self.cells_data[0] = self.cells_data[-2]

        elif self.geometry_data.bcs_type[0] == "Zero Neumann":

            self.cells_data[0] = self.cells_data[1]

        else:
            raise NotImplementedError(f"Type of BC {self.geometry_data.bcs_type[0]} not implemented")

        # Right boundary
        if self.geometry_data.bcs_type[1] == "Periodic":

            self.cells_data[-1] = self.cells_data[1]

        elif self.geometry_data.bcs_type[1] == "Zero Neumann":

            self.cells_data[-1] = self.cells_data[-2]

        else:
            raise NotImplementedError(f"Type of BC {self.geometry_data.bcs_type[1]} not implemented")


    def flux_reconstruction(self):

        """ Wave speed estimation """
        # Wave speed at cell centres are pointed to ICs and cell centred values,
        # needs to be improved in the future to half time-steps or better estimates
        self._c = 0.5 * (self.c_centre_values[:-1] + self.c_centre_values[1:])

        """ Upwind reconstruction """
        for face in range(self.nx + 1):

            if self._c[face][0] >= 0.:
                self.cells_data[face + 1].faces_new[0] = self.cells_data[face].faces_new[1].copy()
            else:
                self.cells_data[face].faces_new[1] = self.cells_data[face + 1].faces_new[0].copy()

