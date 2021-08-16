"""Molecule enumerator blog:
https://greglandrum.github.io/rdkit-blog/tutorial/substructure/2021/05/13/intro-to-the-molecule-enumerator.html
"""

from rdkit import Chem
from rdkit.Chem import rdMolEnumerator
from rdkit.Chem import rdTautomerQuery
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.drawOptions.minFontSize = 10
Draw.SetComicMode(IPythonConsole.drawOptions)
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)
import rdkit

print(rdkit.__version__)
import time

print(time.asctime())

qry = Chem.MolFromMolBlock(
    """
  Mrv2108 08012107372D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 11 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.4167 7.8734 0 0
M  V30 2 C -3.7503 7.1034 0 0
M  V30 3 C -3.7503 5.5633 0 0
M  V30 4 N -2.4167 4.7933 0 0
M  V30 5 C -1.083 5.5633 0 0
M  V30 6 C -1.083 7.1034 0 0
M  V30 7 N 0.3973 7.5279 0 0
M  V30 8 N 0.3104 5.0377 0 0
M  V30 9 C 1.2585 6.2511 0 0
M  V30 10 * 0.3539 6.2828 0 0
M  V30 11 C 1.5089 8.2833 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 6
M  V30 8 1 5 8
M  V30 9 1 8 9
M  V30 10 2 7 9
M  V30 11 1 10 11 ENDPTS=(2 8 7) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
"""
)

print(qry)
bndl = rdMolEnumerator.Enumerate(qry)
print(bndl)

for m in bndl:
    print(m)
