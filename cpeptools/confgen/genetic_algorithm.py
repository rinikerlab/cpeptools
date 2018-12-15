"""
1. when charges too strong the charged partners overlap
    - play around with the weights between the bound matrix and charges
1. fine-tuning the angle per iteration
1. mutation:
    -
"""
class GA:
    def __init__(self, objectivefn, N):
        self.objectivefn = objectivefn
        self.N = N

    def createInitialPop(self, db, targetlensmi):
        # The population will always be in sorted order by decreasing
        # objectivefn() (i.e. most similar first). This will simplify
        # top-slicing and tournament selection.
        pop = []
        while len(pop) < self.N:
            smi = random.choice(db)
            mol = pybel.readstring("smi", smi)
            # Select molecules with a similar SMILES length but with a low value of the objectivefn
            if HasCommonValence(mol.OBMol) and abs(targetlensmi - len(smi)) < 10 and self.objectivefn(smi) < 0.2:
                mol.OBMol.SetTitle("")
                pop.append(mol.write("smi", opt={"i":True}).rstrip()) # random order, leave out stereo
        self.pop = sorted(pop, key=lambda x:self.objectivefn(x), reverse=True)

    def createChildren(self):
        # tournament selection of 2 smis
        children = []
        mrange = range(self.N)
        while len(children) < self.N:
            chosenA = sorted(random.sample(mrange, 3))[0]
            chosenB = chosenA
            while chosenB == chosenA:
                chosenB = sorted(random.sample(mrange, 3))[0]
            # unleash the mutants
            mutantsmi, mol = get_mutant(self.pop[chosenA], self.pop[chosenB])
            if not mol:
                continue
            children.append(mutantsmi)
        self.children = sorted(children, key=lambda x:self.objectivefn(x), reverse=True)

    def createNextGen(self):
        # top-slice the existing population and the children
        self.pop = sorted(self.pop[:int(self.N/2)] + self.children[:int(self.N/2)],
                          key=lambda x:self.objectivefn(x), reverse=True)

    def report(self):
        for p in self.pop[:10]:
            print(p, self.objectivefn(p))
        print()
