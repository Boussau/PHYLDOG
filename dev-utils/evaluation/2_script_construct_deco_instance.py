import script_tree,sys

print("This step will prepare deco execution: generating experiments with a Ensembl tree collection.")

if not len(sys.argv) == 4:
  print("Usage :\n "+sys.argv[0]+" basepath experimentName Compara.XX.protein.nhx.emf")
  sys.exit()

basepath = sys.argv[1]
experimentName = sys.argv[2]
comparaFile = sys.argv[3]

analysisDir = basepath + "/analysis/"+experimentName+"/"

for name in ("starting","reconciled"):
  treeFileName = analysisDir+name+"Trees.txt"
  tree_file = open(treeFileName,"r").readlines()

  print "READING COMPARA FILE"
  fichier_ensembl = open(comparaFile,"r").readlines()

  genes = {}
  l = 0
  while l < len(fichier_ensembl):
          if fichier_ensembl[l][:3] == "SEQ":
                  mots = fichier_ensembl[l].split()
                  species = mots[1].split("_")[0].title()+"_"+mots[1].split("_")[1]
                  genes[mots[2]] = [species,mots[3],int(mots[4]),int(mots[5]),mots[6]] 
          l = l + 1
          
  adjacence_file = open(analysisDir+"/deco_adjacences_"+name,"w")
  gene_file = open(analysisDir+"/deco_genes_"+name,"w")

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
      
      
      
      
      
  file_config = open(analysisDir + "/deco_config_"+name,"w")
  file_config.write("trees_file "+treeFileName+"\n")
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
