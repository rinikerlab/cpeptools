from cpeptools import generate_conformers_with_eccentricity
from rdkit import Chem
from rdkit.Chem import AllChem

def test_confgen():
    CsA_smiles = "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C"
    CsA = Chem.MolFromSmiles(CsA_smiles, sanitize = True)
    ref = CsA
    m = Chem.MolFromPDBFile('/home/shuwang/Documents/Modelling/CP/Datum/macrocycles_benchmark/csa_reference.pdb', removeHs = False, sanitize = True)
    CsA = AllChem.AssignBondOrdersFromTemplate(ref, m)
    assert generate_conformers_with_eccentricity(CsA, 1, 0)
