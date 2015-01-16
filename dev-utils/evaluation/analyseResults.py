# -*- coding: utf-8 -*-

import sys
import glob
import os

pathOfTrees = sys.argv[1]
print("Looking for trees into "+pathOfTrees)
os.chdir(pathOfTrees)
treeFiles = glob.glob("*ReconciledTree")
for currTree in treeFiles:
  print("."),