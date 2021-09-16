import pandas as pd
import copy
from . import PandasTools
from rdkit import Chem
from rdkit.Chem import BRICS
from ..mol_ops import get_rings, get_largest_ring
from functools import reduce

def featurise_sidechain_torsion(mol): #XXX name right? here not really torsions indices
    for atm in mol.GetAtoms():
        atm.SetUnsignedProp("idx" , atm.GetIdx())

    sidechains = get_sidechains(mol)
    sidechains = (Chem.RemoveHs(i) for i in sidechains)

    out = []
    for noh_mol in sidechains:
        out.append((noh_mol, decide_paths(noh_mol)))
    
    for m, arr in out:
        m.__sssAtoms = reduce(lambda a,b : a + b, arr)
    return out


def linear_path_to_torsion_indices(path):
    return [path[i:i + 4] for i in range(len(path) -3)]

def build_df(sidechains_paths_arr):
    mol_column = []
    torsion_column = []
    for sidechain, paths in sidechains_paths_arr:
        for path in paths:
            torsions = linear_path_to_torsion_indices(path)
            for t in torsions:
                m = copy.deepcopy(sidechain)
                m.__sssAtoms = t

                real_indices = [atm.GetProp("idx") if atm.GetIsotope() == 0 else atm.GetIsotope() for atm in [sidechain.GetAtomWithIdx(i) for i in t]]
                torsion_column.append("({})".format(",".join([str(i) for i in real_indices]))) #FIXME
                mol_column.append(m)

    df = pd.DataFrame({"Torsion Indices" : torsion_column, 
    "Depiction" : mol_column})
    df.drop_duplicates(subset = "Torsion Indices", inplace = True)
    PandasTools.RenderImagesInAllDataFrames(images=True)
    return  df

def export_sidechain_torsions(mol, discard = []):
    torsion_set = set([]) 
    torsion_output = [] #I stil want to keep the order


    indices = get_largest_ring(mol)
    indices = indices + indices[:3]
    macrocycle_ring_torsion_indices = [tuple(indices[i:i + 4]) for i in range(len(indices) -3)]

    for real_indices in macrocycle_ring_torsion_indices:
        if real_indices in discard or real_indices in torsion_set:
            continue
        else:
            torsion_set.add(real_indices)
            torsion_output.append(real_indices)

    sidechains_paths_arr = featurise_sidechain_torsion(mol)
    for sidechain, paths in sidechains_paths_arr: #add ring1-ring2-side1-side2 torsion for each sidechain
        ring2, side1, side2, *_ = paths[0]  #only one path is enough per sidechain
        ring2 = int(sidechain.GetAtomWithIdx(ring2).GetIsotope())
        side1 = int(sidechain.GetAtomWithIdx(side1).GetProp("idx"))
        side2 = int(sidechain.GetAtomWithIdx(side2).GetProp("idx"))

        ring1 = [atm.GetIdx() for atm in mol.GetAtomWithIdx(ring2).GetNeighbors() if atm.GetIdx() != side1][0]

        torsion_set.add((ring1, ring2, side1, side2))
        torsion_output.append((ring1, ring2, side1, side2))


    for sidechain, paths in sidechains_paths_arr:
        for path in paths:
            torsions = linear_path_to_torsion_indices(path)
            for t in torsions:
                real_indices = tuple(int(atm.GetProp("idx")) if atm.GetIsotope() == 0 else atm.GetIsotope() for atm in [sidechain.GetAtomWithIdx(i) for i in t])

                if real_indices in discard or real_indices in torsion_set:
                    continue
                else:
                    torsion_set.add(real_indices)
                    torsion_output.append(real_indices)
    return  torsion_output



def decide_paths(noh_mol):
    anchor = [atm.GetIdx() for atm in noh_mol.GetAtoms() if atm.GetAtomicNum() == 0][0] #the break point

    out = [Chem.GetShortestPath(noh_mol, anchor, atm.GetIdx()) for atm in noh_mol.GetAtoms() if atm.GetIdx() != anchor]
    out = [",".join([str(j) for j in i]) for i in out] #each path as , concat string

    #############################################3
    #remove all subpaths
    ignore = set([])
    for i in range(len(out)):
        for j in range(len(out)):
            if i == j:
                continue
            if out[i].startswith(out[j]): #one is subpath of another
                ignore.add(j)
    out = [v.split(",") for i,v in enumerate(out) if i not in ignore]
    out = sorted(out, key = lambda x: len(x)) #each path back to list sorted by length

    #############################################3
    #remove short branching along the path
    ignore = set([])
    for i in range(0, len(out) - 1):
        for j in range(i + 1, len(out)):

            tmp =  out[i][:-1]
            #branching cannot occur at the end of the longer path
            if len(tmp) <= len(out[j]) - 2 and tmp == out[j][:len(tmp)]: 
                ignore.add(i)
    out = [[int(k) for k in v] for i,v in enumerate(out) if i not in ignore]  #idx converted back to int

    #############################################3
    #definitely keeping the longest path(s)
    #for shorter path, if last or second atom is in ring (<= 6 atoms ring), and if longer paths have element also in the same ring, remove it.
    ignore = set([])
    rings = [set(i) for i in noh_mol.GetRingInfo().AtomRings() if len(i) <= 6 ]
    for i in range(0, len(out) - 1):
        for j in range(i + 1, len(out)):
            if len(out[i]) == max(map(lambda x: len(x), out)): #always keep longest path(s)
                continue
            if noh_mol.GetAtomWithIdx(out[i][-1]).IsInRing():
                for r in rings:
                    if out[i][-1] in r:
                        for k in out[j]:
                            if k in r:
                                ignore.add(i)
                                break
            elif noh_mol.GetAtomWithIdx(out[i][-2]).IsInRing():
                for r in rings:
                    if out[i][-2] in r:
                        for k in out[j]:
                            if k in r:
                                ignore.add(i)
                                break
    out = [v for i,v in enumerate(out) if i not in ignore]   
    
    #############################################3
    #filtering on all paths that are of the max length
    ignore = set([])
    for i in range(len(out)):
        for j in range(len(out)):
            if len(out[i]) != max(map(lambda x: len(x), out)): #only looking for longest paths
                continue
            
            if out[i][:-1] == out[j][:-1]: #only the tip differ, then pick the heavier atom
                if noh_mol.GetAtomWithIdx(out[i][-1]).GetAtomicNum() < noh_mol.GetAtomWithIdx(out[j][-1]).GetAtomicNum():
                    ignore.add(i)
                    if j in ignore:
                        ignore.remove(j)
                else:
                    ignore.add(j)
                    if i in ignore:
                        ignore.remove(i)
                    
    # for i in range(0, len(out) - 1):
    #     for j in range(i + 1, len(out)):
    #         if len(out[i]) != max(map(lambda x: len(x), out)): #only looking for longest paths
    #             continue
    #         if noh_mol.GetAtomWithIdx(out[i][-1]).IsInRing(): 
    #             for r in rings:
    #                 if out[i][-1] in r:
    #                     if not np.any([k in r for k in out[j]]):
    #                         keep.add(i) #add i if it does not share the same ring along the path to any already kept path
    out = tuple([v for i,v in enumerate(out) if i not in ignore])
    return out

def get_sidechains(mol):
    #obtain a rough idea of what consists of the core, BRICS quite useful
    # breakage_atompair_list = [i[0] for i in list(BRICS.FindBRICSBonds(mol))] #important list!
    # breakage_bondid_list = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in breakage_atompair_list]

    #break up molecule at the core/side chain boundaries
    *non_core_rings, core_ring_indices = sorted([set(i) for i in get_rings(mol)], key = lambda x: len(x))
    for r in non_core_rings:
        if len(r.intersection(core_ring_indices)) >= 2: #XXX fused ring, should it be == 2?
            core_ring_indices |= r
    
    
    breakage_bondid_list = []
    for bond in mol.GetBonds():
        idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (idx1 in core_ring_indices and idx2 in core_ring_indices) or (idx1 not in core_ring_indices and idx2 not in core_ring_indices):
            continue
        
        idx1, idx2 = (idx1, idx2) if idx1 in core_ring_indices else (idx2, idx1)
        if not mol.GetAtomWithIdx(idx2).IsInRing():
            breakage_bondid_list.append(bond.GetIdx())
            
            
    fragmented_mol = Chem.FragmentOnBonds(mol, breakage_bondid_list, addDummies = True)
    #largest piece is the core
    tmp = sorted(Chem.GetMolFrags(fragmented_mol, asMols = True), key = lambda x: x.GetNumAtoms())
    *sidechains, core = tmp
    
    # sidechains = Chem.ReplaceCore(mol,core,labelByIndex = True)#XXX returns None on some cores
    # sidechains = Chem.GetMolFrags(sidechains, asMols = True)

    #remove sidechains that are too short
    sidechains = [frag for frag in sidechains if len([i for i in frag.GetAtoms() if i.GetAtomicNum() > 1]) > 2] #select only those that has more than just having one heavy atom

    return sidechains 









##############################################################################
# #FIXME
# indices = get_1_4_pairs(look_up[j], "[O:1]=[C:2]@;-[NX3:3]-[CX4H3:4]") + get_1_4_pairs(look_up[j], "[O:1]=[C:2]@;-[NX3:3]-[H:4]")
# out = md.compute_dihedrals(traj, indices)

# raise NotImplementedError()
# #TODO per-residue Ramanchandran plots
