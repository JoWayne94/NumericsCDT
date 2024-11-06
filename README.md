This is a repository for the MFC CDT Numerical Analysis assignment. It contains a conservative finite difference method and the CABARET; written in Python. Currently, it solves first-order, scalar, nonlinear partial differential equations in one spatial dimension.

### Setup

Python 3.7 or higher is required. The libraries required to be installed are also listed in `requirements.txt` with their tested version numbers provided.

### Generating figures 

The scripts generating plots and figures are located in the `analysis` and `plot_scripts` directories. In order to obtain them, run the following shell script:

```sh
$ chmod +x plots.sh
$ ./plots.sh
```

The figures will pop up on your screen, and at the same time saved to your local machine; they are located in the individual subdirectories respectively (e.g., `/analysis`, `/plot_scripts/transport_sine`, `/burgers_tophat`).

### Others

If you want to run individual simulation, the `main_.py` and `setup_.py` scripts are a good starting point and they act as templates for you to customise your test case using different methods and setup specifications. Otherwise, you can also navigate to the `test_cases_` directories to look at how various test cases are set up, and if you want to try them, just run the script that is not `setup.py`. For example:  

```sh
$ cd test_cases_fd/transport_1d/sine/

$ python sine.py
```

The configuration of the simulation can be altered in the `setup.py` scripts that are located in the same directory as the respective python test scripts.

### Solver capabilities

Spatial dimensions:

  - Conservative finite difference - 1D, 2D (minor changes required after some development works)
  - CABARET - 1D

Equation types:

  - Scalar transport
  - Inviscid Burgers'

Time-stepping methods:

  - Forward Euler (only for finite difference)

Boundary conditions types:

  - Periodic
  - Zero Neumann

### Contact details
Please contact Jo Wayne Tan (jwt617@ic.ac.uk) for any questions or issues.
