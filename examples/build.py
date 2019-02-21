from cpeptools import molgen
from rdkit.Chem import AllChem
from rdkit import Chem


f1 = molgen.Fragment("N[C@@H](C)CCCC(=O)")
f1.show()

p1 = molgen.Peptide("AFG")
p1.show()

# print(Chem.MolToPDBBlock(p1.peptide))
type(p1)
p1 = molgen.remove_atoms(p1, "name is OXT")
p1.show()

# atm1 = molgen.select_atoms(p1, "name is N and resnum is 1")[0]
atm2 = molgen.select_atoms(f1, "name is C6")[0]

p1 = molgen.add_bond(p1, "name is N and resnum is 1", atm2)

p1.show()
print(Chem.MolToPDBBlock(p1)) #FIXME resnum seems to increase


atm2 = molgen.select_atoms(p1, "name is N0 and resnum is 4")[0]
p1= molgen.add_bond(p1, "name is C and resnum is  3", atm2, same_mol = True)
p1.show()


print(Chem.MolToPDBBlock(p1)) #FIXME resnum seems to increase



#### Build a D-hydroxproline residue
f1 = molgen.Fragment("C1[C@@H](CN[C@H]1C)O") #somehow defining both stereocenter does not work, only one of them work (below)
f1 = molgen.Fragment("C1C(CN[C@H]1C)O")
f1.show()
f1.set_resname("DHP")
f1.rename("C0", "CB")
f1.rename("C1", "CG")
f1.rename("C2", "CD2")
f1.rename("N3", "N")
f1.rename("C4", "CA")
f1.rename("C5", "C")
f1.rename("O6", "OD1")
f1 = molgen.add_atoms(f1, molgen.create_atom(1, "HD1"), "name is OD1")
f1.show()

print(Chem.MolToPDBBlock(f1))
AllChem.EmbedMultipleConfs(f1,1)
