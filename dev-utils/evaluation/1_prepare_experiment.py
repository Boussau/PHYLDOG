import glob,sys,re

print("This step will collect likelihood and concatenate trees to be used with Deco in the next step.")

if not len(sys.argv) == 3:
  print("Usage :\n "+sys.argv[0]+" basePath experimentName")
  sys.exit()
basePath = sys.argv[1]
experimentName = sys.argv[2]
destDirectory = basePath + "/analysis/" + experimentName + "/"

def concatenateTrees(path,suffix,resultFile):
  filesToConcatenate = glob.glob(path+"/*"+suffix)
  numberOfFiles = len(filesToConcatenate)
  numberDone = 0
  with open(resultFile, 'w') as outfile:
    for currTree in filesToConcatenate:
      numberDone += 1
      print("\r"+str(int(float(numberDone)/numberOfFiles*100))+" %"),
      with open(currTree) as infile:
	outfile.write(infile.read())
      sys.stdout.flush()
  print("\nRead "+ str(len(filesToConcatenate)) + " files")
  if len(filesToConcatenate) == 0:
    print("No trees in " + path + " with suffix ``" + suffix + "''. I give up this experiment.\n____________________\n\n")
    sys.exit(0)
    
def collectLk_beforeAndAfter(dest_file,run_path):
  print("Collecting likelihood in "+run_path)
  regexp_lk = re.compile("/([^\/.]+).opt total logLk: -([0-9\.]+) ; scenario loglk: -([0-9\.]+)")
  regexp_initlk = re.compile("/([^\/.]+).opt ; Initial sequence and reconciliation likelihood: -([0-9\.]+)")
  lkfile = open(dest_file,"w")
  clientFiles = glob.glob(run_path+"/Client_*.out")
  dict_families_final_llk = dict()
  dict_families_init_loglk = dict()
  for currClientFile in clientFiles:
    with open(currClientFile, 'r') as currClient:
      for clientLine in currClient:
	found = regexp_lk.search(clientLine)
	if found:
	  family = found.group(1)
	  lk_total = found.group(2)
	  lk_scenario = found.group(3)
	  lk_sequences = str(float(lk_total) - float(lk_scenario))
	  dict_families_final_llk[family] = (lk_total,lk_sequences,lk_scenario)
	found = regexp_initlk.search(clientLine)
	if found:
	  family = found.group(1)
	  lk = found.group(2)
	  dict_families_init_loglk[family] = lk

	  
	  
  for family in dict_families_final_llk.keys():
    lkfile.write(family + "," + dict_families_init_loglk[family] + "," + dict_families_final_llk[family][0]+ "," + dict_families_final_llk[family][1]+ "," + dict_families_final_llk[family][2] +"\n")


concatenateTrees(basePath+"/result/"+experimentName,"StartingTree",destDirectory+"/startingTrees.txt")
concatenateTrees(basePath+"/result/"+experimentName,"ReconciledTree",destDirectory+"/reconciledTrees.txt")
collectLk_beforeAndAfter(destDirectory+"/collected_lk.csv",basePath+"/run/"+experimentName)

