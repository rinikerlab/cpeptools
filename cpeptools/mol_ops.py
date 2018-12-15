from rdkit import Chem
# from rdkit.Chem import AllChem

def get_carbonyl_O(mol):
    return [i[0] for i in get_atom_mapping(mol, "[C]=[O:1]")]
def get_amine_H(mol):
    return [i[0] for i in get_atom_mapping(mol, "[N]-[H:1]")]
def get_1_4_pairs(mol, smirks_1_4 = "[O:1]=[C:2]@;-[NX3:3]-[CX4H3:4]" ):
    return get_atom_mapping(mol, smirks_1_4)
def get_1_5_pairs(mol, smirks_1_5 = "[O:1]=[C]@;-[CX4H1,CX4H2]-[NX3H1]-[H:5]"):
    return get_atom_mapping(mol, smirks_1_5)

def get_atom_mapping(mol, smirks):
    qmol = Chem.MolFromSmarts(smirks)
    ind_map = {}
    for atom in qmol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num:
            ind_map[map_num - 1] = atom.GetIdx()
    map_list = [ind_map[x] for x in sorted(ind_map)]
    matches = list()
    for match in mol.GetSubstructMatches(qmol, uniquify = False) :
        mas = [match[x] for x in map_list]
        matches.append(tuple(mas))
    return matches


def get_rings(mol):
    for ring in mol.GetRingInfo().AtomRings():
        yield ring

def get_largest_ring(mol):
    out = []
    for r in get_rings(mol):
        if len(r) > len(out):
            out = r
    return list(out)
