https://doi.org/10.5281/zenodo.17800947
# Calculating Globe Temperature from Regional Climate Model Data

This repository contains the Python code for calculating Wet Bulb Globe Temperature (WBGT), modified from the Cython code under the pyWBGT package (Kong & Huber, 2022; Qinkong, 2022; Liljegren et al., 2008). From this code, individual components of WBGT such as Globe Temperature can be calculated. A notebook showing the same is included.

**These scripts support the analysis done in Puri et al. (2026).**

**Citation**: If you reuse the code, please cite 
 _Puri, R., Nam, C. C. W., & Bouwer, L. M. (2026). Rethinking infrastructure thermal exposure beyond air temperature: globe temperature from regional climate models. Environmental Research Letters, 21(12), 124006. https://doi.org/10.1088/1748-9326/ae7390_ 
 and 
 _Riddhima Puri. (2025). Globe Temperature from Regional Climate Models. Zenodo. https://doi.org/10.5281/ZENODO.17800947_


## Details of the files included

- **calculate_GT_example.ipynb**  
  *A Jupyter notebook demonstrating how to calculate Globe Temperature (GT) using the `_wbgt.py` script and EURO-CORDEX data. It shows data loading, pre-processing, and the use of xclim.*

- **_wbgt.py**  
  *A Python module implementing the calculation of Wet Bulb Globe Temperature (WBGT), Globe Temperature (Tg), and Natural Wet Bulb Temperature (Tnwb), integrated with the xclim open-source library (Bourgault et al., 2023). It is modified from the open-source pyWBGT package (Kong and Huber, 2022; Qinkong, 2022), based on the original Liljegren et al. (2008) code.
Used as a backend for the notebook and other scripts.*

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

---

**References**
- Bourgault, P., Huard, D., Smith, T. J., Logan, T., Aoun, A., Lavoie, J., ... & Whelan, C. (2023). xclim: Xarray-based climate data analytics. Journal of Open Source Software, 8(85), 5415. doi: 10.21105/joss.05415
- Kong, Q., & Huber, M. (2022). Explicit calculations of wet-bulb globe temperature compared with approximations and why it matters for labor productivity. Earth’s Future, 10 (3). doi: 10.1029/2021EF002334
- Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurement. Journal of Occupational and Environmental Hygiene, 5 (10), 645–655. doi: 10.1080/15459620802310770
- Puri, R., Nam, C. C. W., & Bouwer, L. M. (2026). Rethinking infrastructure thermal exposure beyond air temperature: globe temperature from regional climate models. Environmental Research Letters, 21(12), 124006. doi: 10.1088/1748-9326/ae7390
- Qinkong. (2022). Explicit calculations of wet bulb globe temperature compared with approximations and why it matters for labor productivity (v1.0.0). Zenodo. doi: 10.5281/zenodo.5980536
