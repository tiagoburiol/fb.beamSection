# 🌀 napy.torsion

A Saint-Venant Torsion Module for the NAPy Finite Element Framework

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

## Overview
This repository contains the **Saint-Venant Torsion Module**, developed for the **NAPy Finite Element Framework**, a Python-based environment for numerical simulation in solid and structural mechanics, which is an ongoing development project.


**NAPy**, short for **Numeric AirPlane**, was originally conceived as a computational platform for the numerical analysis of aerospace structures. Over time, it has evolved into a general-purpose finite element framework designed for research, education, and engineering applications.

## Features
- Implements the **classical Saint-Venant torsion theory**.  
- Computes **shear stress distribution**, **torsional stiffness**, and **warping functions**.  
- Fully written in **Python**, using an **object-oriented and modular design**.  
- Integrates with NAPy’s components for **mesh generation**, **boundary conditions**, and **post-processing**.  
- Supports **structured and unstructured meshes**.  
- Ready for extension to **nonlinear torsion**, **coupled problems**, and **aeroelastic analyses**.

## Implementation Notes
The finite element formulation adopted in this module follows the standard weak form of Saint-Venant’s torsion problem.  
The resulting linear system is assembled using NAPy’s internal FEM infrastructure, ensuring compatibility with other NAPy modules.

## Requirements
- Python ≥ 3.10  
- NumPy  
- SciPy  
- Matplotlib  
- (Optional) meshio or gmsh for mesh import/export

## Usage

See the [multicell profile example](https://github.com/tiagoburiol/napy.torsion/blob/main/notebooks/multicell_profile_example.ipynb).

![Distribution of the warping function and Shear stress distribution](https://github.com/tiagoburiol/napy.torsion/blob/main/imagens/fig1_readme.png?raw=true)


Or

```bash
# Clone the repository
git clone https://github.com/yourusername/napy-torsion.git
cd napy-torsion


# Example: open the multicell profile example using jupyter
jupyter-notebook notebooks/multicell_profile_example.ipynb
```

## Citation
If you use this module in your research, please cite it as:

> Pedro F. M. Pires, Pedro F. M., Buriol, Tiago M. e Santos, Tiago dos (2025). *Saint-Venant Torsion Module for the NAPy Finite Element Framework*. GitHub repository.  
> [https://github.com/yourusername/napy-torsion](https://github.com/yourusername/napy-torsion)

## License

This project is distributed under the MIT License — see the LICENSE
 file for details.

## Future Work

- Nonlinear torsion formulation
- Coupling with elasticity and thermal modules
- Validation with analytical and experimental results
- GUI integration for torsion analysis setup