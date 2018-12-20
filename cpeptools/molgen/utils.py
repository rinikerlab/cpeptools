import copy
from rdkit import Chem

_aa321 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

_bondtypes = {1: Chem.BondType.SINGLE,
              1.5: Chem.BondType.AROMATIC,
              2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE,
              4: Chem.BondType.QUADRUPLE,
              5: Chem.BondType.QUINTUPLE,
              6: Chem.BondType.HEXTUPLE,
              7: Chem.BondType.ONEANDAHALF,}
def wrap_mol_derivative(func):
    def wrapper(*args, **kwargs):
        # print(args) #first arg is mol
        original_mol_derivative = copy.deepcopy(args[0]) #don't overwrite

        # print(type(original_mol_derivative))
        rdmol = func(*args, **kwargs)
        if not hasattr(original_mol_derivative, "load_rdmol"):
            return rdmol
        else:
            original_mol_derivative.load_rdmol(rdmol)
            return original_mol_derivative
    return wrapper
#####################################################
#####################################################
#####################################################
#####################################################

def _select(selection_string):
    """
    TODOs:
        - include chain name and chain id
    """
    logic_words = {"(", ")", "==", "!=", "and", "or"} #True, False
    rdkit_func_keyword = ".GetPDBResidueInfo()." #FIXME matches better than this?

    selection_string = selection_string.replace("is not", "!=")
    selection_string = selection_string.replace("is", "==")
    selection_string = selection_string.replace("(", " ( ")
    selection_string = selection_string.replace(")", " ) ")

    selection_string = selection_string.lower() #should be after is/ is not replacement

    selection_string = selection_string.replace("resname", "atm.GetPDBResidueInfo().GetResidueName().strip().lower() ")
    selection_string = selection_string.replace("resnum", "str(atm.GetPDBResidueInfo().GetResidueNumber()) ")
    # selection_string = selection_string.replace("index", "str(atm.GetIdx()) ") #headache as there is offset bettwen pdb index and rdkit mol atom index
    selection_string = selection_string.replace("name", "atm.GetPDBResidueInfo().GetName().strip().lower() ")

    selection_string = selection_string.split()
    for idx,val in enumerate(selection_string):
        if val not in logic_words and rdkit_func_keyword not in val:
            selection_string[idx] = "\"{}\"".format(selection_string[idx])
    return "{}".format(" ".join(selection_string))

def select_atoms(mol, selection):
    selection = _select(selection)
    return list(eval("filter(lambda atm : {} , mol.GetAtoms())".format(selection)))

@wrap_mol_derivative
def remove_atoms(mol, selection):
    """
    """
    to_remove = select_atoms(mol, selection)
    to_remove = [i.GetIdx() for i in to_remove]
    to_remove.sort()
    if len(to_remove) == 0:
        print("No atom fits selection")
        return mol

    # old_mol = mol
    mol = copy.deepcopy(mol)
    mol = Chem.RWMol(mol)
    for idx in reversed(to_remove):
        mol.RemoveAtom(idx)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise ValueError("Removing bridging atoms are not allowed")
        # return old_mol
    return mol.GetMol()



@wrap_mol_derivative
def reorder_atoms(mol):
    """change index of the atoms to ensure atoms are ordered by ascending residue number
    """
    order = [(i.GetPDBResidueInfo().GetName().strip(), i.GetPDBResidueInfo().GetResidueNumber()) for i in mol.GetAtoms()] # currently the atom name is not used in sorting
    order = [i[0] for i in sorted(enumerate(order), key  = lambda x : x[1][1])]
    return Chem.RenumberAtoms(mol, order)


@wrap_mol_derivative
def add_atoms(mol, atom, selection, bondtype = 1, aromatic = False):
    """
    - will add an atom to all sites that fits the @selection
    - thus all newly added atoms will have the same pdb record (except index)
    - only use it for adding an atom to the same residue as selection
    """
    bondtype = _bondtypes[bondtype]
    to_add = select_atoms(mol, selection)
    mol = copy.deepcopy(mol)
    mol = Chem.RWMol(mol)
    for atm1 in to_add:
        atm2 = create_atom(atom.GetAtomicNum(), atom.GetPDBResidueInfo().GetName())

        atm2.GetPDBResidueInfo().SetIsHeteroAtom(False)
        atm2.GetPDBResidueInfo().SetResidueName(atm1.GetPDBResidueInfo().GetResidueName())
        atm2.GetPDBResidueInfo().SetResidueNumber(atm1.GetPDBResidueInfo().GetResidueNumber())
        atm2.GetPDBResidueInfo().SetChainId(atm1.GetPDBResidueInfo().GetChainId()),
        mol.AddAtom(atm2)
        # print(atm1.GetIdx(), atm2.GetIdx(), bondtype)
        mol.AddBond(atm1.GetIdx(), mol.GetNumAtoms() - 1, bondtype)
        mol.UpdatePropertyCache()
        # mol = add_bond(atm1, atm2, bondtype, same_mol = True)
        # atom.GetPDBResidueInfo().SetName(" {: <3s}".format(name))
        #TODO reorder
    Chem.GetSSSR(mol)
    return reorder_atoms(mol)


def create_atom(number, name, aromatic = False):
    name = name.strip()
    atm = Chem.Atom(number)
    atm.SetIsAromatic(aromatic)
    atm.SetMonomerInfo(Chem.AtomPDBResidueInfo())
    atm.GetPDBResidueInfo().SetName(" {: <3s}".format(name)) #means padding with the space character
    atm.GetPDBResidueInfo().SetOccupancy(0)
    atm.GetPDBResidueInfo().SetTempFactor(0)
    return atm

# #FIXME
# def add_bond(atm1, atm2, bondtype = 1, same_mol = False): #TODO both idx needs to be in the same molecule
#     """
#     - !!!! should only be used after all the atom names are sorted/solved
#     - atm1 atm2 order matters
#     - usually user case of additional inter-residue linkage, e.g. disulphide bridges
#
#     - GetOwningMol() creates a new copy of molecule, this is desired
#     """
#
#     if same_mol:
#         mol1 = atm1.GetOwningMol()
#         mol2 = atm2.GetOwningMol()
#         if Chem.MolToSmiles(atm1.GetOwningMol()) != Chem.MolToSmiles(atm2.GetOwningMol()):
#             raise AssertionError("atm1 and atm2 cannot come from the same molecule as they have non-equivelent SMILES code")
#
#         mol = atm1.GetOwningMol()
#         idx1, idx2 = atm1.GetIdx(), atm2.GetIdx()
#     else:
#         mol1 = atm1.GetOwningMol()
#         mol2 = atm2.GetOwningMol()
#         # resnum ordering
#         max_resnum = max([atm.GetPDBResidueInfo().GetResidueNumber() for atm in mol1.GetAtoms()])
#         for atm in mol2.GetAtoms():
#             atm.GetPDBResidueInfo().SetResidueNumber(max_resnum + atm.GetPDBResidueInfo().GetResidueNumber())
#
#         mol = Chem.CombineMols(mol1, mol2)
#         idx1, idx2 = atm1.GetIdx(), atm2.GetIdx() + mol1.GetNumAtoms()
#
#     bondtype = _bondtypes[bondtype]
#     mol = Chem.RWMol(mol)
#     mol.AddBond(idx1, idx2, bondtype)
#     mol.UpdatePropertyCache()
#     Chem.GetSSSR(mol)
#     return mol.GetMol()

@wrap_mol_derivative
def add_bond(mol, selection, atom, bondtype = 1, same_mol = False): #TODO both idx needs to be in the same molecule
    """
    - !!!! should only be used after all the atom names are sorted/solved
    - atm1 atm2 order matters
    - usually user case of additional inter-residue linkage, e.g. disulphide bridges

    - GetOwningMol() creates a new copy of molecule, this is desired
    - addes to mol
    """
    to_add_bond = select_atoms(mol, selection)
    if len(to_add_bond) < 1:
        raise ValueError("No atom fits selection in mol")
    if same_mol:
        mol2 = atom.GetOwningMol()
        if Chem.MolToSmiles(mol) != Chem.MolToSmiles(mol2):
            raise AssertionError("atom cannot belong to mol as they have non-equivelent SMILES code")

        idx2 = atom.GetIdx()
    else:
        mol2 = atom.GetOwningMol()
        # resnum ordering
        max_resnum = max([atm.GetPDBResidueInfo().GetResidueNumber() for atm in mol.GetAtoms()])
        for atm in mol2.GetAtoms():
            atm.GetPDBResidueInfo().SetResidueNumber(max_resnum + atm.GetPDBResidueInfo().GetResidueNumber())

        idx2 =  atom.GetIdx() + mol.GetNumAtoms() #this needs to be before `mol = Chem.CombineMols(mol, mol2)`
        mol = Chem.CombineMols(mol, mol2)


    bondtype = _bondtypes[bondtype]
    mol = Chem.RWMol(mol)
    for atm in to_add_bond:
        idx1 = atm.GetIdx()
        mol.AddBond(idx1, idx2, bondtype)
        mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

def break_bond(atm1, atm2):
    raise NotImplementedError()


##############
######################################

def make_peptide_bond(mol, res = None, start = "N", end = "C", delete = "OXT"):
    """Performs one condesation rxn between a molecule and residue/itself
    default creates peptide bond

    Parameters
    ----------
        mol : rdkmol
            Main molecule on which reaction is performed
        res : rdkmol
            None or a single residue, when it is None, self-condensation is performed on @mol
        start : str
            atom name of one of the two atoms to which connection is established
        end : str
            atom name of the other atom to which connection is established
        delete : str
            default to hydroxy oxygen thats eliminated during condensation

    Returns
    -------
    mol : rdkmol
        modified molecule
    """
    startIdx, endIdx, deleteIdx  = -1, -1, -1
    for idx, atm in enumerate(mol.GetAtoms()):  # get the last occurence of end and delete atomname
        if atm.GetPDBResidueInfo().GetName().strip() == end:
            endIdx = idx
        elif atm.GetPDBResidueInfo().GetName().strip() == delete:
            deleteIdx = idx

    if res is not None: #residue addition
        lastResNum = -1
        for idx, atm in enumerate(mol.GetAtoms()):
            lastResNum = atm.GetPDBResidueInfo().GetResidueNumber()
        lastResNum += 1
        for idx, atm in enumerate(res.GetAtoms()): #get the last occurence of start atomname
            atm.GetPDBResidueInfo().SetResidueNumber(lastResNum)
            if atm.GetPDBResidueInfo().GetName().strip() == start:
                startIdx = idx
        startIdx += mol.GetNumAtoms()
        mol = Chem.CombineMols(mol, res)
    else: #cyclisation
        for idx, atm in enumerate(mol.GetAtoms()): #get the first occurence of start atomname
            if atm.GetPDBResidueInfo().GetName().strip() == start:
                startIdx = idx
                break

    mol = Chem.RWMol(mol)
    mol.AddBond(startIdx, endIdx, Chem.BondType.SINGLE)
    mol.RemoveAtom(deleteIdx)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

def remove_terminal_atom(mol, condition):
    """Removes atom(s) on the molecule based on expression

    Parameters
    ----------
        mol : rdkmol
            Main molecule on which deletion is performed
        condition : function
            returns boolean values for each atom in @mol, those with True are deleted

    Returns
    -------
    mol : rdkmol
        modified molecule
    """
    mol = Chem.RWMol(mol)
    todelete = [atm.GetIdx() for atm in filter(condition, copy.deepcopy(mol).GetAtoms())]
    for idx in reversed(todelete):
        mol.RemoveAtom(idx)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()


#more general than add_atom_at_site()
def add_terminal_atom(mol, condition, newAtomicNumber, newAtomInfo = lambda atm : Chem.AtomMonomerInfo(), aromatic = lambda atm: False, bondType = 1): #TODO make bondType also a function?
    """create new atom(s) based on expression

    Parameters:
    ----------
    mol : rdkmol
        Main moelcule one which addition is perforemd
    condition : function
        outputs how many atoms to add at each atom site of @mol
    newAtomicNumber : function
        outputs element types to add at each atom site of @mol
    newAtomInfo : function
        outputs the rdkit AtomInfo block information add at each atom site of @mol
    aromatic : function
        True/False value for atom aromaticity
    bondType : float

    Returns
    -------
    mol : rdkmol
        modified molecule

    Example:
    --------
    residue = Chem.MolFromSequence("KK", flavor = 0)
    residue = add_terminal_atom(mol = residue,
        condition = lambda atm: atm.GetPDBResidueInfo().GetName().strip() == "N",
        newAtomName = lambda atm : "CN",
        newAtomicNumber = lambda atm: 6)
    print(Chem.MolToPDBBlock(residue))
    """
    def mapper(condition, mol):
        for idx, repeats in enumerate(map(condition, mol.GetAtoms())):
            for j in range(repeats):
                yield mol.GetAtomWithIdx(idx)


    bondType = _bondtypes[bondType]
    mol = Chem.RWMol(mol)
    for atom in mapper(condition, copy.deepcopy(mol)):
        newatom = Chem.Atom(newAtomicNumber(atom))
        newatom.SetIsAromatic(aromatic(atom))
        newatom.SetMonomerInfo(newAtomInfo(atom))
        # newatom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName = " {: <3s}".format(newAtomName(atom)),
        # residueName = atom.GetPDBResidueInfo().GetResidueName(),
        # residueNumber = atom.GetPDBResidueInfo().GetResidueNumber(),
        # chainId = atom.GetPDBResidueInfo().GetChainId(),
        # ))

        mol.AddAtom(newatom)
        mol.AddBond(atom.GetIdx(), mol.GetNumAtoms() - 1, bondType)
    # mol.UpdatePropertyCache(strict = False)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

def add_atom_at_site(mol, site, name, atomicNumber = 6, aromatic = False, bondType = 1):
    """
    ! Assumes the molecule has PDBResidueInfo()
    """
    bondType = _bondtypes[bondType]
    tmpmol = copy.deepcopy(mol)
    mol = Chem.RWMol(mol)
    for atom in tmpmol.GetAtoms():
        # print(atom.GetPDBResidueInfo().GetName())
        if atom.GetPDBResidueInfo().GetName().strip() == site:
            atom.SetIsAromatic(aromatic)
            atom.SetAtomicNum(atomicNumber)
            atom.GetPDBResidueInfo().SetName(" {: <3s}".format(name))
            mol.AddAtom(atom)
            mol.AddBond(atom.GetIdx(), mol.GetNumAtoms() - 1, bondType)
            mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()
