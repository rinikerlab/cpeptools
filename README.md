cPepTools
==============================
[//]: # (Badges)
<!-- [![Travis Build Status](https://travis-ci.org/hjuinj/cPepTools.png)](https://travis-ci.org/hjuinj/cPepTools) -->
[![Travis Build Status](https://travis-ci.org/hjuinj/cpeptools.png)](https://travis-ci.org/hjuinj/cpeptools)
[![codecov](https://codecov.io/gh/hjuinj/cPepTools/branch/master/graph/badge.svg)](https://codecov.io/gh/hjuinj/cPepTools/branch/master)

This package includes various tools that is being used on the ongoing computational modelling studies of macrocycles, with particular keen on cyclic peptidic species.

The functionalities are primarily split into three submodules:
- Molgen : to construct an insillico (cyclic)molecular entity that contains the atom names (for atom type assignment) and connectivity information including bond orders.
- Confgen : to construct 3D coordinates for atoms in a molecular entity
- Metrics : to analyse a given set of 3D structures 

The tools rely heavily on RDKit and mdtraj functionalities and enable interconversions between these two molecular representations.

### Installation

This package only supports Python 3.5 or above.

The cPeptools packages is built on various packages that is conda-installable. Therefore, either [anaconda or miniconda](https://docs.anaconda.com/anaconda/install/) is needed as prerequisites.

Once conda is available, execute the following commands to install all essential dependencies:
```
conda install -c anaconda numpy scipy scikit-learn
conda install -c omnia mdtraj
conda install -c rdkit rdkit
```
Please ensure your rdkit version is v.2020.03 or above.

After all dependencies are installed, navigate to a directory where you wish to clone this repository and execute:
```
git clone https://github.com/rinikerlab/cpeptools.git
cd cpeptools
python setup.py build install
```

### Usage
Published use involving `cPeptools` is documented here.

- [Improving Conformer Generation for Small Rings and Macrocycles Based on Distance Geometry and Experimental Torsional-Angle Preferences](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025) 

In the method section of this publication, we describe the use of elliptical 2D geometry to confine the conformational space of conformer generation. Figure 2 is repreduced below, and we would like to bias sampling of the rightmost-like cyclic molecules. This requires the useage of this package.
<div align="center">
<img src="https://pubs.acs.org/na101/home/literatum/publisher/achs/journals/content/jcisd8/0/jcisd8.ahead-of-print/acs.jcim.0c00025/20200331/images/medium/ci0c00025_0003.gif" alt="depict_conf"></img>
</div>

The following code snippet is the minimum code needed:
```python
import cpeptools
from rdkit import Chem
from rdkit.Chem import AllChem


mol = Chem.MolFromSmiles("C-C(-C)-C-[C@@H]1-N-C(=O)-[C@H](-C)-N(-C)-C(=O)-[C@H]2-C-C-C-N-2-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](C)-N(-C)-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](-C)-N(-C)-C(=O)-[C@H]2-C-C-C-N-2-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](-C)-N(-C)-C-1=O")
mol = Chem.AddHs(mol)

params = AllChem.ETKDGv3() #contains the modified macrocycle torsion workflow

bmat = cpeptools.bound_matrix_from_ellipse(mol, angle = 0, eccentricity = 0.99)
params.SetBoundsMat(bmat)

AllChem.EmbedMultipleConfs(mol, 3 , params) #conformers written to mol object
```
Check the documentation for [bound_matrix_from_ellipse](https://github.com/rinikerlab/cpeptools/blob/master/cpeptools/confgen/utils.py#L11) for details of additional parameters.

The usage of custom pairwise Coulombic interactions (CPCIs) is also described in the method section. The following code gets CPCIS specific for amide interactions:
```python
amide_cpci = cpeptools.get_amide_pairwise_coulomb_interaction(mol)
params.SetCPCI(amide_cpci)

AllChem.EmbedMultipleConfs(mol, 3 , params) #conformers written to mol object
```
Check the documentation for [get_amide_pairwise_coulomb_interaction](https://github.com/rinikerlab/cpeptools/blob/master/cpeptools/confgen/utils.py#L67) for details of additional parameters.


Another useage described in the method section is to map a 3D ring structure to a 2D represenation with calculated metrics, figure 3 of the paper is reproduced below:
<div align="center">
<img src="https://pubs.acs.org/na101/home/literatum/publisher/achs/journals/content/jcisd8/0/jcisd8.ahead-of-print/acs.jcim.0c00025/20200331/images/medium/ci0c00025_0016.gif" alt="depict_conf"></img>
</div>

The following code depicts a 3D conformer (left) as a 2D representation (right) with calculated metrics:
```python
from cpeptools import metrics
import mdtraj as md

smiles = "C-C(-C)-C-[C@@H]1-N-C(=O)-[C@H](-C)-N(-C)-C(=O)-[C@H]2-C-C-C-N-2-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](C)-N(-C)-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](-C)-N(-C)-C(=O)-[C@H]2-C-C-C-N-2-C(=O)-[C@H](-C-C(-C)-C)-N-C(=O)-[C@H](-C)-N(-C)-C-1=O"
metrics.show_ellipse_fitting("./examples/data/Decapeptide.pdb", smiles = smiles)
```
See doc for [show_ellipse_fitting](https://github.com/rinikerlab/cpeptools/blob/master/cpeptools/metrics/bb_eccentricity.py#L118) for options for additional parameters. To only calculate the values without drawing the figure see [calculate_eccentricity](https://github.com/rinikerlab/cpeptools/blob/master/cpeptools/metrics/bb_eccentricity.py#L49).

To cite this work:
```
@article{Wang2020,
  doi = {10.1021/acs.jcim.0c00025},
  url = {https://doi.org/10.1021/acs.jcim.0c00025},
  year = {2020},
  month = mar,
  publisher = {American Chemical Society ({ACS})},
  author = {Shuzhe Wang and Jagna Witek and Gregory A. Landrum and Sereina Riniker},
  title = {Improving Conformer Generation for Small Rings and Macrocycles Based on Distance Geometry and Experimental Torsional-Angle Preferences},
  journal = {Journal of Chemical Information and Modeling}
}
```

### Copyright

Copyright (c) 2018, shuzhe Wang


#### Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
