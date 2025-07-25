# PorousAirfoil_PotentialFlow

This repository implements a **Source Panel Method with Vortex correction (SPVP)** to simulate airflow around NACA airfoils and compute pressure distributions, lift coefficients (CL), and moment coefficients (CM). The tool also compares results against reference data generated using **XFOIL**.

---

## 🗂️ Repository Structure

- `AoA_Sweep.py`  
  This script compares the aerodynamic performance of a porous and a solid airfoil using the SPVP method across a range of angles of attack. It computes and plots the lift (CL), drag (CD), and lift-to-drag ratio (CL/CD) for both configurations.

- `COMPARAISON.py`  
    This script compares pressure coefficients, lift (CL), and moment (CM) between the SPVP panel method and XFOIL for a NACA airfoil.  
    It includes both local Cp comparisons at a fixed angle of attack and global CL/CM trends over a range of angles.


- `SPVP_Airfoil.py`  
  This script implements the Source and Vortex Panel Method (SPVP) to compute the pressure distribution, lift, moment, and drag coefficients of a NACA 4-digit airfoil. It supports both porous and solid configurations and provides detailed visualizations of aerodynamic quantities.

- `PLOT.py`  
    This module provides a suite of visualization tools for the SPVP method, including plots for airfoil geometry, normal vectors, pressure coefficients, streamlines, and pressure contours.

- `COMPUTATION/`  
  - `COMPUTE.py`: It includes functions to build influence matrices, solve the linear system, and compute lift, drag, and moment coefficients based on airfoil geometry and flow characteristics.
  - `Hydraulic_Resistance.py`: Models resistances for porous airfoil simulation.

- `GEOMETRY/`  
  - `GEOMETRY.py`: Constructs panel geometry for the airfoil.  
  - `Hydraulic_GEOMETRY.py`: Adapts geometry for porous surfaces.  
  - `NACA.py`: Generates 4-digit NACA airfoil coordinates.  
  - `PANEL_DIRECTION.py`: Changes the panels directions if needed.

- `X_FOIL/`  
  - `X_FOIL.py`: Reads XFOIL-generated data files.  
  - `X_FOIL_DATA/`: Directory containing `.dat` files exported from XFOIL.

---

## ⚙️ Features

- **SPVP Solver**  
  Models inviscid, incompressible flow using source and vortex panels.

- **Lift and Moment Computation**  
  Integrates pressure distribution to compute CL and CM.

- **XFOIL Integration**  
  Imports `.dat` files generated with XFOIL for validation and comparison.

- **Porous Surface Modeling (Experimental)**  
  Adds porous resistance through hydraulic modeling components.

---

## 🚀 Getting Started

### 1. Requirements

Install required dependencies:

```bash
pip install numpy matplotlib
