
# TubulePy
TubulePy allows calculation of microtubule protofilament skew angles (_&theta;_) and helical parameters for different protofilament numbers (_N_) and helical-start numbers (_S_), taking into account variations in their lattice parameters: the tubulin monomer subunit repeat (_a_), their rise along the left-handed helices (_r_), and the separation between protofilaments (_&delta;x_).
<div align="center">
<img width="293" height="442" alt="LAM" src="https://github.com/user-attachments/assets/64617499-8e6e-4a62-ab61-ca6dafa0ff47" />
</div>

Contact: denis.chretien@cnrs.fr

# Citing TubulePy
To cite TubulePy in your publications, use at least one of these references:

Chrétien, D., & Fuller, S. D. (2000). Microtubules switch occasionally into unfavorable configurations during elongation. Journal of Molecular Biology, 298(4), 663‑676. https://doi.org/10.1006/jmbi.2000.3696

Guyomar _et al._ (in preparation, to be updated).

# Installation
Download TubulePy.py and place it in a dedicated folder. Download the three PDB files (Alpha.pdb, Alpha_Beta.pdb and Alpha_Beta_kin.pdb) and place them in a subfolder named PDB_models.

•	It is recommended to create a TubulePy environment.

•	Python 3.9+ versions recommended (runs fine on 3.11/3.13).

•	Needs Tkinter available with the Python build (standard on most macOS/Windows installers; on some Linux distros: run 'install python3-tk').

•	External app: Latest version of [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) installed and its ChimeraX binary accessible (path is user-selectable in the app).

•	Packages for full feature set (projections and gallery): numpy, mrcfile, Pillow.

•	Optional packages: matplotlib (needed for the _N~&theta;_ plot and as a fallback saver if Pillow is missing).

•	One-line install (inside your desired env): run 'pip install numpy mrcfile Pillow matplotlib’ in the terminal before running TubulePy.

•	Once in your environment, type at the terminal: python TubulePy.py

# TubulePy interface
<img width="4241" height="1999" alt="Interface" src="https://github.com/user-attachments/assets/a23562b8-a14e-4941-a702-ed52af59dde4" />

TubulePy is composed of two tabs (Parameters and ChimeraX). The first tab allows calculation of protofilament skew angles and helical parameters of microtubule N_S configurations. The second tab allows computation of PDB-based microtubule models, MRC derived 3D models, as well as projected densities using ChimeraX.

Instructions are provided in the document TubulePy.pdf.

# History
TubulePy is based on the Lattice Accommodation Model (LAM) that allows prediction of microtubule protofilament skew angles for any N_S configuration [(Chrétien and Fuller, 2000)](https://doi.org/10.1006/jmbi.2000.3696). Microtubules with fractional helical-starts mimicking unique C- or D-type lattices are computed (Guyomar _et al._, in preparation).

TubulePy has been written using [Microsoft Visual Studio Code](https://github.com/microsoft/vscode), with the help of IA Agents to build the Graphical User Interface.
