"""
Inviscid (scalar) Burgers' 2D test case

Author: JWT
"""

# Import third-party libraries
import sys, os, imageio
import matplotlib.pyplot as plt

# Import user inputs
from setup import *

# Configure system path
if (sys.platform[:3] == 'win') or (sys.platform[:3] == 'Win'):
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..\..')))
else:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

# Import grid classes
from src.grid.grid1D import Grid1D
from src.grid.grid2D import Grid2D

# Import flux definitions & spatial and temporal schemes
from src.spatialschemes.finitedifference import *
from src.temporalschemes.multistage import *
from src.flux.numericalfluxes import *

# Import equation types
from src.solvers.scalartransport import *
from src.solvers.inviscidburgers import *


if __name__ == '__main__':
    """
    main()
    """

    if NO_OF_DIMENSIONS == 1:
        grid = Grid1D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_1d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_1d
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_1d
        else:
            raise NotImplementedError('Equation type is not implemented.')

    elif NO_OF_DIMENSIONS == 2:
        grid = Grid2D

        if SPATIAL_SCHEME == 'Finite difference':
            spatial_scheme = conservative_fd_2d
        else:
            raise NotImplementedError('Spatial discretisation not implemented.')

        if EQUATION_TYPE == 'Scalar transport':
            eqn_type = scalar_transport_2d
        elif EQUATION_TYPE == 'Inviscid Burgers':
            eqn_type = inviscid_burgers_2d
        else:
            raise NotImplementedError('Equation type is not implemented.')

    else:
        raise NotImplementedError('Number of spatial dimensions is not supported.')

    if FLUX_SCHEME == 'Upwind':
        numerical_flux = upwind
    else:
        raise NotImplementedError('Numerical flux of choice is not implemented.')


    """ Create variable object """
    var = grid(NO_OF_GRID_POINTS, DOMAIN_BOUNDARIES, BOUNDARY_CONDITIONS_TYPE)

    '''======PLEASE INPUT YOUR INITIAL CONDITIONS BELOW======'''
    # x coordinates
    x = var.points_data.x_coordinates
    # y coordinates
    y = var.points_data.y_coordinates
    X, Y = np.meshgrid(x, y)
    # Initial datum/conditions defined with a numpy array
    ic = np.where(X < 0., -1., 1.)
    '''======================================================'''

    """ Define driving velocity field  """
    var.variables_data.values = ic
    v = np.ones((var.ny + 2, var.nx + 2, 2))
    dt = CFL * min(var.dx, var.dy) / np.max(np.abs(ic))


    """ Plotting parameters and IC visualisations """
    # Directory to save frames
    frame_dir = 'frames'
    os.makedirs(frame_dir, exist_ok=True)

    # Initialize empty list to store frame filenames
    frame_files = []
    lw = 1.5
    x_left = DOMAIN_BOUNDARIES[0][0]
    x_right = DOMAIN_BOUNDARIES[1][0]
    y_bottom = DOMAIN_BOUNDARIES[0][1]
    y_top = DOMAIN_BOUNDARIES[2][1]
    z_min = np.min(ic)
    z_max = np.max(ic)

    # Customise your visuals here
    # plt.rcParams.update({
    #     "text.usetex": True,
    #     "font.family": "Helvetica"
    # })
    fig = plt.figure(figsize=(10, 5))

    # Subplot 1: imshow
    ax1 = fig.add_subplot(1, 2, 1)
    im = ax1.imshow(ic, extent=[x_left, x_right, y_bottom, y_top], origin='lower', cmap='viridis', aspect='auto',
                    vmin=z_min, vmax=z_max)
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$y$')
    cbar = fig.colorbar(im, ax=ax1, label=r'$u$')  # Individual color bar for imshow

    # Subplot 2: plot_surface with adjustable camera angle
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax2.plot_surface(X, Y, ic, cmap='RdBu_r', edgecolor='none')

    # Adjust the camera angle (elevation, azimuth)
    ax2.view_init(elev=30, azim=45)  # Change these values to adjust the view
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$y$')
    ax2.set_zlabel(r'$u$')
    ax2.set_xlim([x_right, x_left])
    ax2.set_ylim([y_top, y_bottom])
    ax2.set_xlim([x_right, x_left])
    ax2.set_zlim([z_min, z_max])
    plt.suptitle(r'Initial datum', fontsize=16, usetex=True)

    # Show the plot
    plt.tight_layout()  # rect=[0, 0, 1, 0.95]
    plt.pause(1)


    """======START TIME EVOLUTION======"""
    t = 0 # Initial time
    n = 0 # Time iteration count
    while abs(FINAL_TIME - t) > 1.e-8:

        var.reset_fluxes()

        if t + dt >= FINAL_TIME: dt = FINAL_TIME - t

        """ Equation type """
        eqn_type(dt, v, numerical_flux, var)

        """ Spatial discretisation """
        spatial_scheme(var, 1.)

        """ Time evolution """
        forward_euler(var, dt)

        # Replot
        plt.cla()
        fig = plt.figure(figsize=(10, 5))

        # Subplot 1: imshow
        ax1 = fig.add_subplot(1, 2, 1)
        im = ax1.imshow(ic, extent=[x_left, x_right, y_bottom, y_top], origin='lower', cmap='viridis', aspect='auto',
                        vmin=z_min, vmax=z_max)
        ax1.set_xlabel(r'$x$')
        ax1.set_ylabel(r'$y$')
        ax1.set_xlim([x_left, x_right])
        ax1.set_ylim([y_bottom, y_top])
        cbar = fig.colorbar(im, ax=ax1, label=r'$u$')  # Individual color bar for imshow

        # Subplot 2: plot_surface with adjustable camera angle
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        surf = ax2.plot_surface(X, Y, ic, cmap='RdBu_r', edgecolor='none')

        # Adjust the camera angle (elevation, azimuth)
        ax2.view_init(elev=30, azim=45)  # Change these values to adjust the view
        ax2.set_xlabel(r'$x$')
        ax2.set_ylabel(r'$y$')
        ax2.set_zlabel(r'$u$')
        ax2.set_xlim([x_right, x_left])
        ax2.set_ylim([y_top, y_bottom])
        ax2.set_zlim([z_min, z_max])
        plt.suptitle(r'$CFL = $' + str(CFL) + ', Time = ' + '%.2f' % (t + dt) + ' s', fontsize=16, usetex=True)

        # Show the plot
        plt.tight_layout()
        plt.pause(0.01)

        # Save current frame as an image file
        frame_filename = os.path.join(frame_dir, f'frame_{n:03d}.png')
        plt.savefig(frame_filename)
        frame_files.append(frame_filename)

        t += dt
        n += 1
        dt = CFL * min(var.dx, var.dy) / np.max(np.abs(var.variables_data.values))

    # Create the GIF
    with imageio.get_writer('burgers.gif', mode='I', duration=0.002) as writer:
        for filename in frame_files:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Optional: Clean up the frame files after GIF creation
    for filename in frame_files:
        os.remove(filename)

    """ Analytical solution """
#     Analytical solution exists so find it and plot here
