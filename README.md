# 🌀 fb.beamSection

A Saint-Venant Torsion Module for the FlightBEND Finite Element Framework

[![PyPI version](https://badge.fury.io/py/FlightBEND.svg)](https://pypi.org/project/FlightBEND/1.0.0/) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17460543.svg)](https://doi.org/10.5281/zenodo.17460543) 

## Overview

This repository contains the **Saint-Venant Torsion Module**, developed for the **FlightBEND Finite Element Framework**, a Python-based environment for numerical simulation in solid and structural mechanics, which is an ongoing development project.

The name **FlightBEND** is a tribute to Bernoulli, Euler, Saint-Venant, and Santos Dumont. The project was originally conceived as a computational platform for the numerical analysis of aerospace structures. Over time, it has evolved into a general-purpose finite element framework designed for research, education, and engineering applications.

This project has been carried out at the Numerical Analysis Laboratory (LANum) of the [Federal University of Santa Maria - UFSM](https://www.ufsm.br/), involving graduate students from the [Mechanical Engineering Program - PGMec](https://www.ufsm.br/cursos/pos-graduacao/santa-maria/pgmec) and faculty members from the Departments of Mechanical Engineering and Mathematics.

## Features

- Implements the **classical Saint-Venant torsion theory**.  
- Computes **shear stress distribution**, **torsional stiffness**, and **warping functions**.  
- Fully written in **Python**, using an **object-oriented and modular design**.  
- Integrates with FlightBEND’s components for **mesh generation**, **boundary conditions**, and **post-processing**.  
- Supports **structured and unstructured meshes**.  
- Ready for extension to **nonlinear torsion**, **coupled problems**, and **aeroelastic analyses**.

## Implementation Notes

The finite element formulation adopted in this module follows the standard weak form of Saint-Venant’s torsion problem.  
The resulting linear system is assembled using FlightBEND’s internal FEM infrastructure, ensuring compatibility with other FlightBEND modules.

## Requirements

- Python ≥ 3.10  
- NumPy  
- SciPy  
- Matplotlib  

## Usage

See the [multicell profile example](https://github.com/tiagoburiol/fb.beamSection/blob/main/tutorials/multicell_profile_example.ipynb).

![Distribution of the warping function and Shear stress distribution](https://raw.githubusercontent.com/tiagoburiol/fb.beamSection/refs/heads/main/images/fig1_readme.png)

Or

```bash
# Clone the repository
git clone https://github.com/tiagoburiol/fb.beamSection.git
cd fb.beamSection


# Example: open the multicell profile example using jupyter
jupyter-notebook notebooks/multicell_profile_example.ipynb
```

Other examples of usage can be found in [fb.beamSection/examples](https://github.com/tiagoburiol/fb.beamSection/tree/main/examples). The figure below is from the [wing_torsion_box.ipynb](https://github.com/tiagoburiol/fb.beamSection/blob/main/examples/5_wing_torsion_box/wing_torsion_box.ipynb) example.

![Warping function](https://raw.githubusercontent.com/tiagoburiol/fb.beamSection/refs/heads/main/images/fig2_readme.png)

## Install on pip

You can also install the `flightbend` package direcly from PyPl on the comand prompt (Windows) or Linux terminal:

```powershell
pip install flightbend
```

*Note : Clone the repository to get the most up-to-date version. PyPl updates may be infrequent*. 

## Verification

The `.inp` files for the ABAQUS 3D exemples from the paper may be found at the [verification_files](https://github.com/tiagoburiol/fb.beamSection/blob/main/verification_files) folder in a compacted `.7z` format.

## Mesh generation

For all examples we used the mesh generation from the [GiD](https://www.gidsimulation.com/) pre-processing software. If you wish to create `.py` files like the ones from our examples you'll also need to install the [MATFem](https://www.gidsimulation.com/downloads/educational-finite-element-codes-matfem/) plugin and download our customized [.bas](https://github.com/tiagoburiol/fb.beamSection/blob/main/mesh_files) file. The plugin allows for element material assigment and the `.bas` file exports the mesh as a python dictionary saved in the variable `mesh_data` containing the keys `'materials'`, `'coordinates'` and `'elements'`.  Check our tutorial on how to use GiD for this application on the [notebooks folder](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks).

## Citation

If you use this module in your research, please cite it as:

> Pedro F. M. Pires, Pedro F. M., Buriol, Tiago M. e Santos, Tiago dos (2025). *Saint-Venant Torsion Module for the FlightBEND Finite Element Framework*. GitHub repository.  
> [https://github.com/tiagoburiol/fb.beamSection](https://github.com/tiagoburiol/fb.beamSection)

## License

This project is distributed under the MIT License — see the LICENSE
 file for details.

## Future Work

- Nonlinear torsion formulation
- Coupling with elasticity and thermal modules
- Validation with analytical and experimental results
- GUI integration for torsion analysis setup
