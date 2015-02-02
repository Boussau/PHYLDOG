# -*- coding: utf-8 -*-

import sys
import glob
import os
import script_tree
import re


# FUNCTIONS

def concatenateTrees(path,suffix,resultFile):
  filesToConcatenate = glob.glob(path+"*"+suffix)
  with open(resultFile, 'w') as outfile:
    for currTree in filesToConcatenate:
      print("."),
      with open(currTree) as infile:
	outfile.write(infile.read())
      sys.stdout.flush()
  print("\nRead "+ str(len(filesToConcatenate)) + " files")
  if len(filesToConcatenate) == 0:
    print("No trees in " + path + " with suffix ``" + suffix + "''. I give up this experiment.\n____________________\n\n")
    sys.exit(0)

def collectLk_beforeAndAfter(dest_file,run_path):
  print("Collecting likelihood in "+run_path)
  regexp_lk = re.compile("/([^\/.]+).opt total logLk: -([0-9\.]+) ;")
  regexp_initlk = re.compile("/([^\/.]+).opt ; Initial sequence and reconciliation likelihood: -([0-9\.]+)")
  lkfile = open(dest_file,"w")
  clientFiles = glob.glob(run_path+"Client_*.out")
  dict_families_total_llk = dict()
  dict_families_init_loglk = dict()
  for currClientFile in clientFiles:
    with open(currClientFile, 'r') as currClient:
      for clientLine in currClient:
	found = regexp_lk.search(clientLine)
	if found:
	  family = found.group(1)
	  lk = found.group(2)
	  dict_families_total_llk[family] = lk
	found = regexp_initlk.search(clientLine)
	if found:
	  family = found.group(1)
	  lk = found.group(2)
	  dict_families_init_loglk[family] = lk
	  # puting it in lk, in case no lk improvement is made
	  dict_families_total_llk[family] = lk
	  
	  
  for family in dict_families_total_llk.keys():
    lkfile.write(family + "," + dict_families_init_loglk[family] + "," + dict_families_total_llk[family] +"\n")



os.system("rm RESULTS -rf")

pathOfTrees = sys.argv[1]
print("________________\n\n\nLooking for trees into "+pathOfTrees)

#finding exp name
exp_name = pathOfTrees.split("/")[-2]
print("Experiment name: "+ exp_name)

collectLk_beforeAndAfter("./lk_before_and_after.csv",pathOfTrees.replace("result","run"))
collectLk_beforeAndAfter("./lk_ensembl.csv",pathOfTrees.replace("result","run").replace(exp_name,"exp_ensemblTrees"))

concatenateTrees(pathOfTrees,"ReconciledTree","reconciledTrees.txt")
concatenateTrees(pathOfTrees,"StartingTree","startingTrees.txt")


tree_file = open("reconciledTrees.txt","r").readlines()


# Now Eric’s scripts build DeCo instance

print "READING COMPARA FILE"
fichier_ensembl = open("./Compara.73.protein.nhx.emf","r").readlines()

genes = {}
l = 0
while l < len(fichier_ensembl):
	if fichier_ensembl[l][:3] == "SEQ":
		mots = fichier_ensembl[l].split()
		species = mots[1].split("_")[0].title()+"_"+mots[1].split("_")[1]
		genes[mots[2]] = [species,mots[3],int(mots[4]),int(mots[5]),mots[6]] 
	l = l + 1
	
name = "reconciledTrees"
adjacence_file = open("deco_adjacences_"+name,"w")
gene_file = open("deco_genes_"+name,"w")

restricted_genes = []
for line in tree_file:
  tree = script_tree.readTree(line)
  leaves = script_tree.getLeavesNames(tree)
  for l in leaves:
    restricted_genes.append(l)

print len(tree_file),len(restricted_genes)
#print restricted_genes

def compare(x,y):
  #print x,y
  if genes[x][0] != genes[y][0]:
    return cmp(genes[x][0],genes[y][0])
  elif genes[x][1] != genes[y][1]:
    return cmp(genes[x][1],genes[y][1])
  else:
    return cmp(genes[x][2],genes[y][2])

restricted_genes.sort(lambda x,y: compare(x,y))

for g in range(len(restricted_genes)):
  gene_file.write(genes[restricted_genes[g]][0]+" "+restricted_genes[g]+"\n")
  if (g + 1 < len(restricted_genes) and
      genes[restricted_genes[g]][0] == genes[restricted_genes[g+1]][0] and
      genes[restricted_genes[g]][1] == genes[restricted_genes[g+1]][1]):
    adjacence_file.write(restricted_genes[g]+" "+restricted_genes[g+1]+"\n")
    
    
    
adjacence_file.close()
gene_file.close()
    
file_config = open("deco_config_"+name,"w")
file_config.write("trees_file reconciledTrees.txt\n")
file_config.write("genes_file deco_genes_"+name+"\n")
file_config.write("species_file species_tree.tree\n")
file_config.write("adjacencies_file deco_adjacences_"+name+"\n")
file_config.write("exp_name deco_"+name+"\n")
file_config.write("directory RESULTS\n")
file_config.write("ReconcilDone false\n")
file_config.write("INPUT_FORMAT 0\n")
file_config.write("OUTPUT_FORMAT 1\n")
file_config.write("sep |\n")
file_config.write("Adj_percentage 0\n")
file_config.write("Gain 2\n")
file_config.write("Break 1\n\n")
file_config.close()

name="startingTrees"
file_config = open("deco_config_startingTrees","w")
file_config.write("trees_file startingTrees.txt\n")
file_config.write("genes_file deco_genes_reconciledTrees\n")
file_config.write("species_file species_tree.tree\n")
file_config.write("adjacencies_file deco_adjacences_reconciledTrees\n")
file_config.write("exp_name deco_"+name+"\n")
file_config.write("directory RESULTS\n")
file_config.write("ReconcilDone false\n")
file_config.write("INPUT_FORMAT 0\n")
file_config.write("OUTPUT_FORMAT 1\n")
file_config.write("sep |\n")
file_config.write("Adj_percentage 0\n")
file_config.write("Gain 2\n")
file_config.write("Break 1\n\n")
file_config.close()


name="ensemblTrees"
file_config = open("deco_config_ensemblTrees","w")
file_config.write("trees_file Compara.73.protein.nhx.emf_arbres\n")
file_config.write("genes_file deco_genes_reconciledTrees\n")
file_config.write("species_file species_tree.tree\n")
file_config.write("adjacencies_file deco_adjacences_reconciledTrees\n")
file_config.write("exp_name deco_"+name+"\n")
file_config.write("directory RESULTS\n")
file_config.write("ReconcilDone false\n")
file_config.write("INPUT_FORMAT 0\n")
file_config.write("OUTPUT_FORMAT 1\n")
file_config.write("sep |\n")
file_config.write("Adj_percentage 0\n")
file_config.write("Gain 2\n")
file_config.write("Break 1\n\n")
file_config.close()



print(">>>> DECO ON RECONCILED TREES\n")

launchString = "/panhome/bigot/deco_eric " + os.getcwd() +"/deco_config_reconciledTrees"
print(launchString)
os.system(launchString)

print("\n>>>>> END OF Deco on Reconciled Trees\n\n")

print(">>>> DECO ON STARTING TREES\n")

launchString = "/panhome/bigot/deco_eric " + os.getcwd() +"/deco_config_startingTrees"
print(launchString)
os.system(launchString)

print("\n>>>> END OF DeCo on Starting Trees\n\n")


print(">>>> DECO ON ENSEMBL TREES \n")

launchString = "/panhome/bigot/deco_eric " + os.getcwd() +"/deco_config_ensemblTrees"
print(launchString)
os.system(launchString)

print("\n>>>> END OF DeCo on Ensembl Trees\n\n")


os.chdir("RESULTS")

especes_actuelles = {}
degres = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
file_genes = open("../deco_genes_reconciledTrees").readlines()
file_adj = open("../deco_adjacences_reconciledTrees").readlines()
graphe = {}
especes_actuelles = {}
for line in file_genes:
  mots = line.split()
  if especes_actuelles.has_key(mots[0]):
    especes_actuelles[mots[0]] = especes_actuelles[mots[0]] + 1
  else:
    especes_actuelles[mots[0]] = 1
  graphe[mots[1]] = []
for line in file_adj:
  mots = line.split()
  graphe[mots[0]].append(mots[1])
  graphe[mots[1]].append(mots[0])
for v in graphe.keys():
  degres[len(graphe[v])] = degres[len(graphe[v])] + 1
sortie = open("summary_content_extant","w")
for e in especes_actuelles.keys():
	sortie.write(e+" "+str(especes_actuelles[e])+"\n")
sortie = open("summary_degres_extant","w")
for d in range(len(degres)):
	sortie.write(str(d)+" "+str(degres[d])+"\n")


#for name in ["deco_ProfileNJ_from_Ensembl"]:
for name in ["deco_ensemblTrees","deco_reconciledTrees","deco_startingTrees"]:
	print name
	fichier_especes = open(name+"_OUTPUT_species","r").readlines()

	especes_ancestrales = {}
	especes_actuelles = {}

	for line in fichier_especes:
		if len(line.split()) > 2:
			especes_ancestrales[line.split()[0]] = []
		else:
			especes_actuelles[line.split()[0]] = []
			#print line.split()[0]
			
	total = 0
	fichier_genes = open(name+"_OUTPUT_genes","r").readlines()
	for line in fichier_genes:
		words = line.split()
		if not especes_actuelles.has_key(words[0]):
			especes_ancestrales[words[0]].append(words[1])
			total = total + 1.0
		else:
			especes_actuelles[words[0]].append(words[1])

	sortie = open("summary_content_"+name,"w")
	for e in especes_ancestrales.keys():
		sortie.write(e+" "+str(len(especes_ancestrales[e]))+"\n")
	#~ sortie = open("summary_content_extant","w")
	#~ for e in especes_actuelles.keys():
		#~ sortie.write(str(len(especes_actuelles[e]))+"\n")
		
	fichier_adjacences = open(name+"_OUTPUT_adjacencies","r").readlines()
	graph = {}
	for e in especes_ancestrales.keys():
		graph[e] = {}
		for g in especes_ancestrales[e]:
			graph[e][g] = []
	for line in fichier_adjacences:
		words = line.split()
		graph[words[0]][words[1]].append(words[2])
		graph[words[0]][words[2]].append(words[1])
	
	degres = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	for e in especes_ancestrales.keys():
		for g in especes_ancestrales[e]:
			if len(degres) > len(graph[e][g]):
				degres[len(graph[e][g])] = degres[len(graph[e][g])]  + 1
	sortie = open("summary_degres_"+name,"w")
	for d in range(len(degres)):
		sortie.write(str(d)+" "+str(degres[d])+"\n")
		
sortie.close()
	
os.system("R --vanilla < ../script_graphes.R")

os.system("rm -f "+pathOfTrees+"/*jpg")
os.system("rm -f "+pathOfTrees+"/*pdf")
os.system("mv *pdf "+pathOfTrees)
os.system("mv ../*csv "+pathOfTrees)
