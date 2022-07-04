# Fluorescence recovery after photobleaching (FRAP) 

This is a repository of code developed in the Membrane Enzymology group of University of Groningen for the analysis of protein diffusion in bacterial cytoplasm by FRAP.
The repository contains Python scripts working within ImageJ (Fiji). The analysis requires FRAP acquisitions on single cells.     

* Create_selection: working with files ending with .lsm. It finds the first photo-bleached frame, measures the background value, and finds the selection line through the cell by fitting an ellipse to the thresholded image. Output is a .csv file containing selection information.
* Adjust_selection: visually adjusting the selection line aligned along the cell.
* Create_kymograph: creating kymograph based on the fluorescent intensity profile
* 1d_heat_diff_fit: simulating 1D diffusion simulation in Python is based on the heat diffusion equation of the Crank-Nicolson scheme (doi: 10.1017/S0305004100023197) and calculating diffusion coefficients. This code needs be run outside ImageJ (Fiji).
* plotting array: plotting the fluorescent intensity along the selection line from the experimental data, simulation data and residuals.   
