#!/bin/sh
# This script is used to plot all the figures included in the report.

cd analysis
python convergence_test.py
python disp_relations.py
cd ..
cd plot_scripts/transport_gaussian
python gaussian.py
cd ..
cd transport_sine
python sine.py
cd ..
cd transport_tophat
python top_hat.py
cd ..
cd burgers_sine
python sine.py
cd ..
cd burgers_tophat
python tophat.py
