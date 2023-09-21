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
It is important to ensure that all the required dependencies are installed in your working environment. It would be much easier if you create a conda environment before installation of packages. The following packages are required, **pymol**, **rdkit**, **pandas**, **padelpy**, **joblib**, **meeko**, **Autodock Vina**, **java**, **scipy**, and **scikit-learn**.
```bash
conda create -n MLDockKit
conda activate MLDockKit
```
Then, install pymol before installing other packages:
```bash
conda install -c conda-forge pymol-open-source

conda install -c conda-forge openbabel
conda install -c cyclus java-jre

pip install -U numpy vina

pip install rdkit pandas padelpy joblib meeko scikit-learn scipy

pip install MLDockKit
```

### Usage
#### Computating Lipinsky descriptors
```python
>>>from MLDockKit import calculate_lipinski_descriptors
>>>calculate_lipinski_descriptors("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Predictioning pIC50
```python
>>>from MLDockKit import predict_pIC50
>>>predict_pIC50("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Docking with protein, pdb_id 5gs4
```python
>>>from MLDockKit import prot_lig_docking
>>>prot_lig_docking("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
```
#### Visualization of docking results
This opens pymol for visulization and analysis. If you need help on pymol analysis please have a look on [pymol documentation](https://fitzkee.chemistry.msstate.edu/sites/default/files/ch8990/pymol-tutorial.pdf)

```python
>>>from MLDockKit import vizualize_dock_results
>>>vizualize_dock_results()
```

### Acknowledgment
Autodock Vina and pymol were greatily used in writing the codes for molecular docking and visualization. If you use these functions in your work, please cite the original publications for [vina](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203) and [pymol](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=ab82608e9a44c17b60d7f908565fba628295dc72#page=44)

We extracted part of Angel Ruiz Moreno's Jupyter_Dock [Jupyter Dock](https://github.com/AngelRuizMoreno/Jupyter_Dock) to include it in our visualization function. 

### Contribution
We welcome any contributions. Should you notice a bug, please let us know through issues in the, [GitHub Issue Tracker](https://github.com/clabe-wekesa/MLDockKit/issues)

### Authors
Edwin mwakio, [Dr. Clabe Wekesa](https://www.ice.mpg.de/246268/group-members) and [Dr. Patrick Okoth](https://mmust.ac.ke/staffprofiles/index.php/dr-patrick-okoth),  
Department of Biological Sciences, [Masinde Muliro University of Science and Technology](https://www.mmust.ac.ke/)
 
