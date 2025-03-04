# MLDockKit
This is a simple platform for computing Lipinsky's Rule of five using the rdkit package, predicting pIC50 of canonical SMILES that are potential targets against Oestrogen receptor alpha protein as ant-prostate cancer agaents using a preformatted RandomForest model, and docking of the canonical SMILE with the Oestrogen receptor alpha protein using Audodock Vina package. 
### Purpose of the Package
The purpose of the package is to provide a unified platform for computing prostate cancer drug likeness indicess and performing docking on the same compounds. 
### Features
Important chemoinformatics features of Oestrogen receptor alpha antagonists such as:
    - Lipinsky descriptors
    - Prediction of pIC50
    - Docking and visiualization 
### Getting Started
The package is found on pypi hence can be installed with pip
#### Pre-requisites
Installation of Vina requires [boost](https://www.boost.org/doc/libs/1_83_0/tools/build/doc/html/index.html#bbv2.installation) and [swig](https://www.swig.org/)  
Importantly: Install pymol as the first package in the new environment, this is due to environmental package conflict.
#### Installation
It is important to ensure that all the required dependencies are installed in your working environment. It would be much easier if you create a conda environment before installation of packages. The following packages are required, **pymol**, **rdkit**, **pandas**, **padelpy**, **joblib**, **meeko**, **Autodock Vina**, **java**, **scipy**, **MDAnalysis** and **scikit-learn**.
```bash
conda create -n MLDockKit python=3.11
conda activate MLDockKit
```
MLDockKit requires [meeko](https://github.com/forlilab/Meeko) for ligand preparation. The current version of Meeko is compatible only with Python 3.11 or lower. We will update MLDockKit to support higher Python versions once Meeko is updated.

Then, install pymol before installing other packages:
```bash
conda install -c conda-forge pymol-open-source openbabel vina numpy

conda install -c cyclus java-jre

pip install MLDockKit
```

### Run MLDockKit pipeline

```python
from MLDockKit import MLDockKit

MLDockKit(
    smiles="Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2", 
    output_file="docking_results.txt",
    presentation="cartoon", 
    label_residues=True, 
    show_iAA=False
)
```

### Params:

```bash
1. smiles (str): SMILES string for ligand.
2. output_file (str): File path for saving output.
3. presentation (str): How to display the receptor [[e.g., 'surface', 'sticks', 'spheres', 'cartoon', etc.] default: 'cartoon')].
4. label_residues (bool): Option to label residues (default: True).
5. show_iAA (bool): Option to show interacting amino acids only (default: True).
```
### Output
The pipeline's output is an MLDockKit_output.txt file which contains **Lipinsky descriptos**, **predicted pIC50 value** and the **docking score**. Docking image is rentered in pymol for further analysis by the user. Also, the ligand's and protein's **.sfd** and **.pdpqt** files are rentered in the user's working directory.

### Acknowledgment
Autodock Vina and pymol were greatily used in writing the codes for molecular docking and visualization. If you use these functions in your work, please cite the original publications for [vina](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203) and [pymol](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ab82608e9a44c17b60d7f908565fba628295dc72#page=44)

We extracted part of Angel Ruiz Moreno's Jupyter_Dock [Jupyter Dock](https://github.com/AngelRuizMoreno/Jupyter_Dock) to include it in our visualization function. 

### Contribution
We welcome any contributions. Should you notice a bug, please let us know through issues in the, [GitHub Issue Tracker](https://github.com/clabe-wekesa/MLDockKit/issues)

### Authors
Edwin mwakio, [Clabe Wekesa](https://www.ice.mpg.de/246268/group-members) and [Patrick Okoth](https://mmust.ac.ke/staffprofiles/index.php/dr-patrick-okoth)  
Department of Biological Sciences, [Masinde Muliro University of Science and Technology](https://www.mmust.ac.ke/)
