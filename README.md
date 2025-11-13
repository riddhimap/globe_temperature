# Calculating Globe Temperature

This repository contains the Python code for calculating Wet Bulb Globe Temperature (WBGT), modified from the Cython code under the pyWBGT package (Kong & Huber, 2022; Liljegren et al., 2008). From this code, individual components of WBGT such as Globe Temperature can be calculated. A notebook showing the same is included.

## Details of the files included

- **calculate_GT_example.ipynb**  
  *A Jupyter notebook demonstrating how to calculate Globe Temperature (GT) using the `_wbgt.py` script and CORDEX data. It shows data loading, pre-processing, and the use of xclim.*

- **_wbgt.py**  
  *A Python module implementing the calculation of WBGT, Globe Temperature (Tg), and Natural Wet Bulb Temperature (Tnwb) following Kong et al. (2022) and using physical equations. Used as a backend for the notebook and other scripts.*

- **environment.txt**
  *A text file listing the main Python packages required and their versions to reproduce the computational environment needed to run the two files above.*
  *You can install these using `pip install -r environment.txt`.*

## Usage

- See `calculate_GT_example.ipynb` for a step-by-step example.

## Data
These scripts are designed for use with EURO-CORDEX NetCDF data, but given that you have all the variables needed to calculate GT, you can use this notebook to calculate GT from other climate model/reanalysis data.

---

**Contact:**  
For questions, please contact Riddhima Puri or open an issue on this repository.
