#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Descriptors
from padelpy import padeldescriptor
import joblib
import csv
from rdkit.Chem import AllChem, SDWriter, SDMolSupplier
from pymol import cmd
from vina import Vina
import os
import subprocess


# constants
docking_protein = "5gs4.pdbqt"
prediction_model = "padel_model.joblib"
current_directory = os.getcwd()
file_paths = ["ligand_clean.sdf", "ligand.pdbqt"]


def delete_files_with_extension(directory, extensions):
    """
    Delete files with specified extensions in a directory.
    """
    for file in os.listdir(directory):
        if any(file.endswith(ext) for ext in extensions):
            os.remove(os.path.join(directory, file))

# Delete files in the working directory
current_directory = os.getcwd()
#delete_files_with_extension(current_directory, [".sdf", ".pdbqt"])
delete_files_with_extension(current_directory, [".sdf",".pdbqt"])            


def getbox(selection="sele", extending=6.0, software="vina"):
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)

    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)

    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX = (maxX + minX) / 2
    CenterY = (maxY + minY) / 2
    CenterZ = (maxZ + minZ) / 2

    cmd.delete("all")

    return {"center_x": CenterX, "center_y": CenterY, "center_z": CenterZ}, {
        "size_x": SizeX,
        "size_y": SizeY,
        "size_z": SizeZ,
    }


def pdbqt_to_sdf(pdbqt_file=None, output=None):
    results = [m for m in pybel.readfile(filename=pdbqt_file, format="pdbqt")]
    out = pybel.Outputfile(filename=output, format="sdf", overwrite=True)
    for pose in results:
        pose.data.update({"Pose": pose.data["MODEL"]})
        pose.data.update({"Score": pose.data["REMARK"].split()[2]})
        del pose.data["MODEL"], pose.data["REMARK"], pose.data["TORSDO"]

        out.write(pose)
    out.close()


def calculate_lipinski_descriptors(smiles):
    """Lipinski descriptors: A set of molecular properties used to assess the drug-likeness or pharmacokinetic profile of a chemical compound

    Params
    ------
    smiles: string: An rdkit valid canonical SMILES or chemical structure a compound.

    Usage
    -----
    from MLDockKit import calculate_lipinski_descriptors

    calculate_lipinski_descriptors("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("You entered an invalid SMILES string")

    else:
        descriptors = {
            "Molecular Weight": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "Num H Donors": Descriptors.NumHDonors(mol),
            "Num H Acceptors": Descriptors.NumHAcceptors(mol),
            "Num Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Carbon Count": Descriptors.HeavyAtomCount(mol),
            "Oxygen Count": sum(
                1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8
            ),
        }

        aliases = {
            "Molecular Weight": "Molecular Weight",
            "LogP": "LogP",
            "Num H Donors": "Number Hydrogen Bond Donors",
            "Num H Acceptors": "Number of Hydrogen Bond Acceptors",
            "Num Rotatable Bonds": "Number of Rotatable Bonds",
            "Carbon Count": "Carbon Count",
            "Oxygen Count": "Oxygen Count",
        }

        formatted_descriptors = ""
        for key, value in descriptors.items():
            formatted_descriptors += f"{aliases[key]}: {value}\n"

        return formatted_descriptors


def predict_pIC50(smiles):
    """Prediction model is based on RandomForest regression constructed using a collection of all known cannonical SMILES that interact with Oestrogen Receptor alpha protein stored in the ChEMBL database.

    Params
    ------
    smiles: string: An rdkit valid canonical SMILES or chemical structure a compound.

    Usage
    -----
    from MLDockKit import predict_pIC50

    predict_pIC50("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    """
    # Get the directory of the currently executing script

    script_dir = os.path.dirname(__file__)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("You entered an invalid SMILES string")

    # Convert SMILES to molecule object
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Write the molecule to an SDF file
    sdf_file = os.path.join(script_dir, "molecule.smi")
    data = [[smiles + "\t" + "Compound_name"]]
    with open(sdf_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(data)

    # Process the fingerprints
    padeldescriptor(
        mol_dir=sdf_file,
        d_file=os.path.join(script_dir, "descriptors.csv"),
        detectaromaticity=True,
        standardizenitro=True,
        standardizetautomers=True,
        removesalt=True,
        fingerprints=True,
    )
    data = pd.read_csv(os.path.join(script_dir, "descriptors.csv"))
    X = data.drop(columns=["Name"])

    # Specify the path to the "padel_model.joblib" file
    prediction_model = os.path.join(script_dir, "padel_model.joblib")
    loaded_model = joblib.load(prediction_model)
    y_pred = loaded_model.predict(X)
    predicted_value = y_pred[0]
    predicted_value = format(predicted_value, ".2f")
    return f"Predicted pIC50: {predicted_value}"

def prot_lig_docking(smiles):
    """
    Docking procedure is performed by Autodock Vina on the Oestrogen Receptor alpha protein, pdb_id: 5gs4.
    
    Params
    ------
    smiles: string, an rdkit valid canonical SMILES or chemical structure a compound.
    
    Returns
    ------
    str: Docking score or an error message.
    """
    # Get the directory of the currently executing script
    script_dir = os.path.dirname(__file__)
    current_directory = os.getcwd()

    # Convert SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Error: Invalid SMILES string"

    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Save the ligand to an SDF file
    sdf_file = os.path.join(current_directory, "ligand_clean.sdf")
    writer = SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

    # Prepare ligand using Open Babel
    ligand_pdbqt = os.path.join(current_directory, "ligand.pdbqt")
    mol = next(pybel.readfile("sdf", sdf_file))
    mol.write("pdbqt", ligand_pdbqt, overwrite=True)

    # Load the docking protein
    docking_protein = os.path.join(script_dir, "5gs4.pdbqt")
    cmd.load(filename=docking_protein, format="pdb", object="prot")
    cmd.load(filename=sdf_file, format="sdf", object="lig")
    center, size = getbox(selection="lig", extending=5.0, software="vina")
    cmd.delete("all")

    # Initialize Vina
    v = Vina(sf_name="vina")
    v.set_receptor(docking_protein)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(
        center=[center["center_x"], center["center_y"], center["center_z"]],
        box_size=[size["size_x"], size["size_y"], size["size_z"]],
    )

    # Perform docking
    v.dock(exhaustiveness=10, n_poses=10)

    # Save docking results
    vina_out_file = os.path.join(current_directory, "5gs4_ligand_vina_out.pdbqt")
    sdf_file = os.path.join(current_directory, "5gs4_ligand_vina_out.sdf")
    v.write_poses(vina_out_file, n_poses=10, overwrite=True)
    pdbqt_to_sdf(
        pdbqt_file=vina_out_file,
        output=sdf_file,
    )

    # Process docking results
    results = Chem.SDMolSupplier(sdf_file)
    if not results or len(results) == 0:
        return "No docking results found"

    try:
        docking_score = results[0].GetProp("Score")  # Extract docking score
        return f"Docking score: {docking_score}"
    except Exception as e:
        return f"Error during docking: {str(e)}"


# Run PyMOL with the specified PDBQT files
def vizualize_dock_results():
    """Visualization and protein-ligand interaction in pymol. Users should only run this function after running prot_lig_docking function.

    Usage
    -----
    from MLDockKit import vizualize_dock_results

    vizualize_dock_results()
    """
    script_dir = os.path.dirname(__file__)
    current_directory = os.getcwd()
    docking_protein = os.path.join(script_dir, "5gs4.pdbqt")
    vina_out_file = os.path.join(current_directory, "5gs4_ligand_vina_out.pdbqt")
    return subprocess.run(["pymol", "-Q", docking_protein, vina_out_file])


def MLDockKit(smiles, output_file="MLDockKit_output.txt"):
    """
    Perform the entire molecular modeling pipeline:
    1. Calculate Lipinski descriptors
    2. Predict pIC50
    3. Perform protein-ligand docking
    4. Visualize docking results

    Params:
    smiles (str): An rdkit valid canonical SMILES or chemical structure of a compound.
    output_file (str): File path to save the output.

    Returns:
    str: Summary of the pipeline execution.
    """
    try:
        with open(output_file, "w") as f:
            f.write("." * 200 + "\n")
            # Calculate Lipinski descriptors
            lipinski_descriptors = calculate_lipinski_descriptors(smiles)
            f.write("Lipinski Descriptors"+ "\n")
            f.write(str(lipinski_descriptors))
            f.write("\n" + "." * 200 + "\n")
            print("\n" +'###Computation of Lipinsky descriptors complete'+"\n")
            
            # Predict pIC50
            pIC50_prediction = predict_pIC50(smiles)
            f.write(pIC50_prediction + "\n")
            f.write("\n" + "." * 200 + "\n")
            print('###Prediction of pIC50 complete'+"\n")

            # Perform protein-ligand docking
            docking_result = prot_lig_docking(smiles)
            #f.write("\n")
            f.write(docking_result + "\n")
            f.write("\n" + "." * 200 + "\n")
            print("\n" + '###Docking process complete'+"\n")
            print("##MLDockKit output is saved to " + output_file + "and image rentered in pymol"+"\n")
            #f.write("\n"+"\n")

            # Delete files in the script directory
            script_dir = os.path.dirname(__file__)
            delete_files_with_extension(script_dir, [".smi", ".csv"])
            # Delete files in the working directory
        

        # Iterate over files in the directory
        for file in os.listdir(current_directory):
            if file in file_paths:
                file_path = os.path.join(current_directory, file)  # Get the full path
                os.remove(file_path)    


        # Visualize docking results
        vizualize_dock_results()

    except Exception as e:
        return f"Error occurred: {str(e)}"



