This is a repository for the MFC CDT Numerical Analysis assignment. It contains a conservative finite difference method and the CABARET; written in Python. Currently, it solves first-order, scalar, nonlinear partial differential equations in one spatial dimension.

[//]: # (### Setup)

[//]: # (Python 3.7 or higher is required. The following libraries should also be installed &#40;tested version numbers provided&#41;:)

[//]: # (  - NumPy 1.17.4, 1.22.3)

[//]: # (  - Matplotlib 3.3.1, 3.5.1)

[//]: # (  - SciPy 1.4.1, 1.7.3)

[//]: # (  - LaTeX &#40;used for post-processing&#41;)

[//]: # ()
[//]: # (For convenience, the Quail src directory can be added to PATH. The driver script &#40;`quail`&#41; is located in this directory.)

[//]: # (```sh)

[//]: # ($ export PATH=$PATH:/your/quail/directory/src)

[//]: # (```)

[//]: # (The above line can also be added to the appropriate file &#40;e.g., `~/.bashrc`, `~/.bash_profile`, `~/.profile`&#41; and sourced.)


[//]: # (### Using Quail )

[//]: # (A suite of example 1D and 2D cases for scalar equations, the compressible Euler equations, the compressible Navier-Stokes equations, and chemically reacting flow is available in the `examples` directory. Also showcased are the ADERDG scheme and various stabilization methods &#40;positivity-preserving limiter, WENO limiter, and artificial viscosity&#41;. For instance, to run the 2D isentropic vortex case, do the following:)

[//]: # (```sh)

[//]: # ($ cd examples/euler/2D/isentropic_vortex/)

[//]: # ($ quail isentropic_vortex.py)

[//]: # (```)

[//]: # (Depending on the configuration of your machine, the above command &#40;`quail isentropic_vortex.py`&#41; may not work. In that case, try replacing the first line of `src/quail` with `#!/usr/bin/env python3`. If that still doesn't work, run the following command instead:)

[//]: # (```sh)

[//]: # ($ python /your/quail/directory/src/quail isentropic_vortex.py)

[//]: # (```)

[//]: # (Note that this command doesn't require the Quail src directory to be added to PATH.)

[//]: # ()
[//]: # (Additional tools for performing dissipation and dispersion analysis and plotting basis functions are available in the `tools` directory. To perform said analysis, do the following:)

[//]: # (```sh)

[//]: # ($ cd tools/dissipation_dispersion_analysis/)

[//]: # ($ python plot_dissipation_dispersion_relations.py )

[//]: # (```)

[//]: # (Settings can be changed directly in `plot_dissipation_dispersion_relations.py`.)

[//]: # (To plot 1D basis functions, do the following:)

[//]: # (```sh)

[//]: # ($ cd tools/plot_basis_functions/)

[//]: # ($ python plot_segment_basis_fcn.py  )

[//]: # (```)

[//]: # (Settings can be changed directly in `plot_segment_basis_fcn.py`. Basis functions for triangles and quadrilaterals can also be plotted.)




### Contact details
Please contact Jo Wayne Tan (jwt617@ic.ac.uk) for any questions or issues.
