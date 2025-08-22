# Neuronal heterogeneity of normalization strength in a circuit model

> Song, D., Ruff, D., Cohen, M., & Huang, C. (2024). Neuronal heterogeneity of normalization strength in a circuit model. bioRxiv.
> DOI: [https://doi.org/10.1101/2024.11.22.624903](https://doi.org/10.1101/2024.11.22.624903)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.6084/m9.figshare.29940062-blue)](https://doi.org/10.6084/m9.figshare.29940062)



This repository contains the source code and data to reproduce the results presented in our paper. We provide scripts for simulation, analysis, and figure generation.

For a comprehensive introduction and tutorial on the high-performance C implementations of exponential integrate-and-fire network simulations, please consult this [repository](https://github.com/deyingsong/EIF-network-simulation).

### Repository Structure

```text
.
├── README.md
├── data
├── sim_script
├── src
│   ├── simulation
│   │   ├── +spiking_simulation
│   │   ├── C_source_codes
│   │   └── +physiological_properties
│   ├── utils
│   │   ├── +utils_simulation
│   │   ├── +utils_analysis
│   │   └── +utils_plot
│   └── figure
│       ├── +MainFigure
│       └── +SuppFigure
├── plot
├── results
└── LICENSE
```


## Installation

### Prerequisites
- **MATLAB** R2021b (or compatible)  
  Required toolboxes:  
  - Statistics and Machine Learning Toolbox  
  - Image Processing Toolbox  
  - Curve Fitting Toolbox  

- **C Compiler** (for building MEX files):  
  - **Windows** → Visual Studio 2017/2019/2022, or MinGW-w64 GCC 6.3+  
  - **macOS** → Xcode 12 or 13 (Clang)  
  - **Linux** → GCC 6.3 – 9.x  

### Setup
1. Run the setup script:  
   ```matlab
   setup
   ```
2.  Download the data:

## How to Reproduce Our Results

### Quick Test
To verify your installation, run:
```matlab
demo
```

### Full replicaation
- Main figures:
    ```matlab
    plot/runMainFigures
    ```
- Supplementary figures:
    ```matlab
    plot/runSuppFigures
    ```
- Simulations
    ```matlab
    sim_script/sim_default
    ```

## **Licensing and Contact**

* **License**: MIT.
* **Contact**: please email deyingsong [at] cmu [dot] edu for questions

