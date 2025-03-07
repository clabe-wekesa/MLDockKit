#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import MDAnalysis as mda
import pymol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from padelpy import padeldescriptor
import joblib
import csv
from rdkit.Chem import AllChem, SDWriter, SDMolSupplier
from pymol import cmd,stored
from vina import Vina
import os
import subprocess


# constants
docking_protein = "7te7_prepared.pdbqt"
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


def prepare_ligand(input_sdf: str, output_pdbqt: str):
    # Read the input molecule
    mol = Chem.MolFromMolFile(input_sdf)
    if mol is None:
        raise ValueError(f"Invalid SDF file: {input_sdf}")
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    # Save as pdbqt
    with open(output_pdbqt, 'w') as pdbqt_file:
        pdbqt_file.write(Chem.MolToPDBBlock(mol))


def calculate_lipinski_descriptors(smiles):
    """
    Calculate Lipinski descriptors and assess Rule of 5 and Veber’s Rule violations.

    Parameters
    ----------
    smiles : str
        An RDKit valid canonical SMILES or chemical structure of a compound.

    Returns
    -------
    str
        A formatted string of Lipinski descriptors with Rule of 5, Veber’s Rule, 
        and interpretations for Carbon and Oxygen counts.
    """

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string. Please provide a valid SMILES notation.")

    # Calculate Lipinski descriptors
    molecular_weight = Descriptors.MolWt(mol)
    logP = Descriptors.MolLogP(mol)
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Evaluate Rule of 5 criteria
    mw_status = "Rule of 5 not violated (less than or equal to 500)" if molecular_weight <= 500 else "Rule of 5 violated (greater than 500)"
    logP_status = "Rule of 5 not violated (less than or equal to 5)" if logP <= 5 else "Rule of 5 violated (greater than 5)"
    h_donor_status = "Rule of 5 not violated (less than or equal to 5)" if num_h_donors <= 5 else "Rule of 5 violated (greater than 5)"
    h_acceptor_status = "Rule of 5 not violated (less than or equal to 10)" if num_h_acceptors <= 10 else "Rule of 5 violated (greater than 10)"

    # Format output
    formatted_descriptors = (
        f"Molecular Weight = {molecular_weight:.2f}: {mw_status}\n"
        f"LogP = {logP:.2f}: {logP_status}\n"
        f"Number of Hydrogen Bond Donors = {num_h_donors}: {h_donor_status}\n"
        f"Number of Hydrogen Bond Acceptors = {num_h_acceptors}: {h_acceptor_status}\n"
    )

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
    predicted_value = float(predicted_value)  # Convert back to float
    if predicted_value > 6:
        return f"Predicted pIC50 = {predicted_value}: Active (Strong Inhibition [greater than 6])"
    elif 5 <= predicted_value <= 6:  # Corrected range condition
        return f"Predicted pIC50 = {predicted_value}: Moderate Inhibition (from 5 to 6)"
    else:
        return f"Predicted pIC50 = {predicted_value}: Inactive (Weak or No Inhibition [less than 5])"   

def prot_lig_docking(smiles, num_poses=10, exhaustiveness=10):
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
    sdf_file = os.path.join(current_directory, "ligand_initial.sdf")
    writer = SDWriter(sdf_file)
    writer.write(mol)
    writer.close()

    # Prepare ligand using meeko
    ligand_pdbqt = os.path.join(current_directory, "ligand_prepared.pdbqt")
    subprocess.run(["mk_prepare_ligand.py", "-i", sdf_file, "-o", ligand_pdbqt], check=True)
    

    # Load the docking protein
    docking_protein = os.path.join(script_dir, "7te7_prepared.pdbqt")
    original_protein = os.path.join(script_dir, "7te7_original.pdb")
    original_structure = mda.Universe(original_protein)
    ligand_mda = original_structure.select_atoms("resname I0V")

    # Get the center of the ligand as the "pocket center"
    pocket_center = ligand_mda.center_of_geometry()
    ligand_box = ligand_mda.positions.max(axis=0) - ligand_mda.positions.min(axis=0) + 5

    ## convert ligand_box to list
    pocket_center = pocket_center.tolist()
    ligand_box = ligand_box.tolist() 

    # Initialize Vina
    v = Vina(sf_name="vina")
    v.set_receptor(docking_protein)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=pocket_center, box_size=ligand_box)

    # Perform docking
    v.dock(exhaustiveness=exhaustiveness, n_poses=num_poses)

    # Save docking results
    vina_out_file = os.path.join(current_directory, "ligand_docked.pdbqt")
    sdf_file = os.path.join(current_directory, "ligand_docked.sdf")
    v.write_poses(vina_out_file, n_poses=num_poses, overwrite=True)
    subprocess.run(["mk_export.py", vina_out_file, "-s", sdf_file], check=True)

    # Process docking results
    try:
        docking_score = v.score()  # Extract docking score
        docking_score_value = docking_score[0]  # Get the actual score
        
        # Classify Docking Score
        if docking_score_value <= -9.0:
            docking_interpretation = "Very Strong Binding (less than -9.0 kcal/mol)"
        elif -9.0 < docking_score_value <= -7.0:
            docking_interpretation = "Strong Binding (Between -9.0 to -7.0 kcal/mol)"
        elif -7.0 < docking_score_value <= -5.0:
            docking_interpretation = "Moderate Binding (Between -7.0 to -5.0 kcal/mol)"
        else:
            docking_interpretation = "Weak or No Significant Binding (more than -5.0 kcal/mol)"

        return f"Docking score = {docking_score_value:.3f} kcal/mol: {docking_interpretation}"

    except Exception as e:
        return f"Error during docking: {str(e)}"


def visualize_dock_results(presentation='cartoon', show_iAA=True, color='green'):
    """Visualizes a docking result in PyMOL, keeping only interacting residues.

    Parameters:
    - presentation (str): How to display the receptor (default: 'cartoon').
    - show_iAA (bool): Option to show interacting amino acids only (default: True).
    - color (str): Color for the receptor and interacting residues (default: 'green').
    """

    script_dir = os.path.dirname(__file__)
    current_directory = os.getcwd()
    receptor_file = os.path.join(script_dir, "7te7_H_no_HOHandI0V.pdb")
    ligand_file = os.path.join(current_directory, "ligand_docked.sdf")

    # Ensure required files are available
    if not os.path.exists(receptor_file) or not os.path.exists(ligand_file):
        raise ValueError("Both receptor_file and ligand_file must exist.")

    pymol.finish_launching()

    # Load receptor and ligand
    cmd.load(receptor_file, "receptor")
    cmd.load(ligand_file, "ligand")
    
    # Remove water molecules
    cmd.remove("resn HOH")

    # Show ligand
    cmd.show("sticks", "ligand")
    cmd.color("magenta", "ligand")
    cmd.set("stick_radius", 0.3, "ligand")

    if show_iAA:
        # Hide everything except interacting residues and ligand
        cmd.hide("everything")
        
        # Show only interacting residues within 5Å of the ligand
        cmd.select("interacting_residues", "byres receptor within 5 of ligand")
        cmd.show('sticks', "interacting_residues")
        
        # Apply the user-provided color to interacting residues
        cmd.color(color, "interacting_residues")
        
        # Ensure the ligand remains visible
        cmd.show("sticks", "ligand")
        cmd.color("magenta", "ligand")

        # Store unique interacting residues
        stored.residues = []
        cmd.iterate("interacting_residues and name CA", 
                    "stored.residues.append((resi, resn))")

        # Ensure residues were found
        if stored.residues:
            for resi, resn in set(stored.residues):
                cmd.label(f"resi {resi} and name CA", f'"{resn}-{resi}"')
        else:
            print("No interacting residues found for labeling.")

    else:
        cmd.dss("receptor")  # Assign secondary structure if missing
        cmd.show(presentation, "receptor")
        cmd.color(color, "receptor")  # Apply the color to the receptor

    # Zoom in on the ligand
    cmd.zoom("ligand")


def MLDockKit(smiles, 
              output_file="MLDockKit_output.txt", 
              presentation='cartoon', 
              show_iAA=True, 
              color='green', 
              num_poses=10, 
              exhaustiveness=10
              ):
    """
    Perform the entire molecular modeling pipeline:
    1. Calculate Lipinski descriptors
    2. Predict pIC50
    3. Perform protein-ligand docking
    4. Visualize docking results

    Params:
    smiles (str): SMILES string for ligand.
    output_file (str): File path for saving output.
    presentation (str): How to display the receptor [[e.g., 'surface', 'sticks', 'spheres', 'cartoon', etc.] (default: 'cartoon')].
    show_iAA (bool): Option to show interacting amino acids only (default: True).
    color (str): Color for the receptor and interacting residues (default: 'green').
    num_poses (int): Number of docking poses to generate (default: 10).
    exhaustiveness (int): Controls search thoroughness (default: 10).

    Returns:
    str: Summary of the pipeline execution.
    """
    try:
        with open(output_file, "w") as f:
            f.write("." * 200 + "\n")
            # Calculate Lipinski descriptors
            lipinski_descriptors = calculate_lipinski_descriptors(smiles)
            f.write("Lipinski Descriptors/Rule of 5"+ "\n")
            f.write(str(lipinski_descriptors))
            f.write("\n" + "." * 200 + "\n")
            print("\n" +'###Computation of Lipinsky descriptors complete'+"\n")
            
            # Predict pIC50
            pIC50_prediction = predict_pIC50(smiles)
            f.write(pIC50_prediction + "\n")
            f.write("\n" + "." * 200 + "\n")
            print('###Prediction of pIC50 complete'+"\n")

            # Perform protein-ligand docking
            docking_result = prot_lig_docking(smiles, num_poses=num_poses, exhaustiveness=exhaustiveness)
            f.write(docking_result + "\n")
            f.write("\n" + "." * 200 + "\n")
            print("\n" + '###Docking process complete'+"\n")
            print("##MLDockKit output is saved to " + output_file + "and image rendered in pymol"+"\n")

            # Delete files in the script directory
            script_dir = os.path.dirname(__file__)
            current_directory = os.getcwd()
            delete_files_with_extension(script_dir, [".smi", ".csv"])
            # Remove ligand files
            subprocess.run(["rm", "ligand_initial.sdf", "ligand_prepared.pdbqt"], check=True)
            # Copy receptor file (fixed missing comma in subprocess.run argument list)
            subprocess.run(["cp", os.path.join(script_dir, "7te7_H_no_HOHandI0V.pdb"), "7te7_receptor.pdb"], check=True)


        # Visualize docking results, passing the user-defined parameters
        visualize_dock_results(presentation=presentation,show_iAA=show_iAA,color=color)
        
    except Exception as e:
        return f"Error occurred: {str(e)}"


            