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
It is important to ensure that all the required dependencies are installed in your working environment. It would be much easier if you create a conda environment before installation of packages. The following packages are required, **pymol**, **rdkit**, **pandas**, **gemmi**, **padelpy**, **joblib**,**openbabel**, **Autodock Vina**, **java**, **scipy**, **MDAnalysis** and **scikit-learn**.
```bash
conda create -n MLDockKit
conda activate MLDockKit
```

Then, install pymol before installing other packages:
```bash
conda install -c conda-forge pymol-open-source vina mdanalysis 

conda install -c cyclus java-jre

pip install MLDockKit
```

MLDockKit requires the installation of [Meeko](https://github.com/forlilab/Meeko) for ligand preparation. Please note that the current release of Meeko is only compatible with Python 3.11 or lower. If you're using a Python version higher than 3.11, you must install the [development version from the source](https://meeko.readthedocs.io/en/release-doc/installation.html#from-source) for compatibility within the MLDockKit environment.

To install Meeko in the MLDockKit environment, follow these steps:
```bash
git clone https://github.com/forlilab/Meeko.git
cd Meeko
git checkout develop
pip install .
```

### Run MLDockKit pipeline

```python
>>>from MLDockKit import MLDockKit

>>>MLDockKit(smiles,output_file,presentation,show_iAA,color,num_poses,exhaustiveness)
```

### Params:

```bash
1. smiles (str): SMILES string for ligand.
2. output_file (str): File path for saving output [default: "MLDockKit_results.txt"].
3. presentation (str): How to display the receptor [[e.g., 'surface', 'sticks', 'spheres', 'cartoon', etc.] default: 'cartoon')].
4. show_iAA (bool): Option to show interacting amino acid residues only (default: True).
5. color (str): Color for the receptor and interacting residues (default: 'green')
6. num_poses (int): Number of docking poses to generate (default: 10).
7. exhaustiveness (int): Controls search thoroughness (default: 10).
```

### Example running with default settings

```python
from MLDockKit import MLDockKit

MLDockKit(smiles="Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```

### Output
The pipeline's output is an MLDockKit_output.txt file which contains **Lipinsky descriptos**, **predicted pIC50 value** and the **docking score**. Docking image is rentered in pymol for further analysis by the user. Also, the ligand's and protein's **.sfd** and **.pdpqt** files are rentered in the user's working directory.

### Acknowledgment
Autodock Vina and pymol were greatily used in writing the codes for molecular docking and visualization. If you use these functions in your work, please cite the original publications for [vina](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203) and [pymol](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ab82608e9a44c17b60d7f908565fba628295dc72#page=44)


### Contribution
We welcome any contributions. Should you notice a bug, please let us know through issues in the, [GitHub Issue Tracker](https://github.com/clabe-wekesa/MLDockKit/issues)

### Authors
Edwin mwakio, [Clabe Wekesa](https://www.ice.mpg.de/246268/group-members) and [Patrick Okoth](https://mmust.ac.ke/staffprofiles/index.php/dr-patrick-okoth)  
Department of Biological Sciences, [Masinde Muliro University of Science and Technology](https://www.mmust.ac.ke/)
 
