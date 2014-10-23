import os
import textwrap

#recursively builds a list "liste" of all the files in a directory "nomdir"
def explore(nomdir, liste):
		listefiles=os.listdir(nomdir)
		for i in listefiles:
				if (os.path.isdir(os.path.join(nomdir,i))):
						explore(os.path.join(nomdir,i), liste)
				elif not i.startswith("."):
						liste.append(os.path.abspath(os.path.join(nomdir,i)))
		return()



if __name__ == '__main__':
	print ("\n\n\tWelcome to the PHYLDOG file preparation script.\n\n")
	for l in textwrap.wrap("This script will ask you for two input directories, one containing only alignment files, and one containing only link files. Link files relate sequences in the alignments to species. They are formatted as follows:\n\nSpecies_A:gene_A1\nSpecies_B:gene_B1;gene_B2\nSpecies_A:gene_A2\nSpecies_C:gene_C1\n...\n", replace_whitespace=False):
		print l
	for l in textwrap.wrap ("The script will ask you for a folder that will contain option files, and another that will contain results output by PHYLDOG."):
		print l
	for l in textwrap.wrap ("This script does not intend to cover all the options that PHYLDOG can accept. Rather, it is meant to be a convenient way to quickly generate a list of reasonable option files, with default values that should work in the general case. In particular, this script is pretty dumb about the models of sequence evolution it can set up, and it sets up the same model for all the gene families. PHYLDOG can work with a much larger array of models, including different models for different gene families, but to use them you will need to edit the option files using a text editor or custom scripts.\n\n", replace_whitespace=False):
		print l
	print("For comments and questions, please contact Bastien Boussau:\n\n\tboussau@gmail.com\n\n") 
	print ("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

	#Alignment folder
	listAlns=list()
	SEQDIR = raw_input("Absolute PATH to the folder containing all the alignment files?\t")
	while ( not os.path.isdir(SEQDIR) ) :
		print (SEQDIR + " does not exist.\n")
		SEQDIR = raw_input("Absolute PATH to the folder containing all the alignment files?\t") 
	explore(SEQDIR, listAlns)
	print (SEQDIR + " exists and contains "+str(len(listAlns))+" files, presumably all alignments.\n")
	
	#Data type
	DATATYPE = raw_input("What alphabet? DNA or RNA or CODON or PROTEIN?\t").upper()
	if (not ( "DNA".startswith(DATATYPE) or "RNA".startswith(DATATYPE) or "CODON".startswith(DATATYPE) or "PROTEIN".startswith(DATATYPE) ) ): 
		while (not ( "DNA".startswith(DATATYPE) or "RNA".startswith(DATATYPE) or "CODON".startswith(DATATYPE) or "PROTEIN".startswith(DATATYPE) ) ) :
			print (DATATYPE + " is not an option I understand for an alphabet.\n")
			DATATYPE = raw_input("What alphabet? DNA or CODON or PROTEIN?\t").upper()
	if ("DNA".startswith(DATATYPE) or "RNA".startswith(DATATYPE) ) :
		if ("DNA".startswith(DATATYPE)):
			print ("Alphabet: DNA\n")
		else:
			print ("Alphabet: RNA\n")	
		print ("We are going to be setting a GTR+Gamma4 model, by default, in all gene family-specific option files.\n")
		print ("This can be changed directly in the family-specific option files along with many other options.\n")
	elif ("CODON".startswith(DATATYPE)):
		print ("Alphabet: CODON\n")
		print ("We are going to be setting a YN98 model with a universal code, by default, in all gene family-specific option files.\n")
		print ("This can be changed directly in the family-specific option files along with many other options.\n")
	else :
		print ("Alphabet: PROTEIN\n")
		print ("We are going to be setting a LG08+Gamma4 model, by default, in all gene family-specific option files.\n")
		print ("This can be changed directly in the family-specific option files along with many other options.\n")
	
	#File format
	DATAFORMAT = raw_input("What alignment format? FASTA or MASE or PHYLIP or CLUSTAL or DCSE or NEXUS?\t").upper()
	while (not ( "FASTA".startswith(DATAFORMAT) or "MASE".startswith(DATAFORMAT) or "PHYLIP".startswith(DATAFORMAT) or "CLUSTAL".startswith(DATAFORMAT) or "DCSE".startswith(DATAFORMAT) or "NEXUS".startswith(DATAFORMAT) ) ) :
		print (DATAFORMAT + " is not an option I understand for an alignment format.\n")
		DATAFORMAT = raw_input("What alignment format? FASTA or MASE or PHYLIP or CLUSTAL or DCSE or NEXUS?\t").upper()
	if ("FASTA".startswith(DATAFORMAT) ):  
		print ("Alignment format: FASTA\n")
	elif ("MASE".startswith(DATAFORMAT)):
		print ("Alignment format: MASE\n")
	elif ("PHYLIP".startswith(DATAFORMAT)):
		print ("Alignment format: PHYLIP\n")
	elif ("CLUSTAL".startswith(DATAFORMAT)):
		print ("Alignment format: CLUSTAL\n")
	elif ("DCSE".startswith(DATAFORMAT)):
		print ("Alignment format: DCSE\n")
	elif ("NEXUS".startswith(DATAFORMAT)):
		print ("Alignment format: NEXUS\n")

	#Filtering of the data ?
        FILTERING = raw_input("Do you want to filter the alignments and trees using PHYLDOG or do you trust your data is good to go?\t").upper()
        while (not ( "YES".startswith(FILTERING) or "NO".startswith(FILTERING) ) ) :
                print (FILTERING + " is not an option I understand for whether you want to filter the data or not.\n")
                FILTERING = raw_input("Do you want to filter the data?\t").upper()
        if ("YES".startswith(FILTERING)):
		print ("Filtering the data\n.")
	else:
                print ("Not filtering the data\n.")

	#link file folder
	listLinks = list()
	LINKDIR = raw_input("Absolute PATH to the folder containing all the link files?\t")
	pathExists = os.path.isdir(LINKDIR)
	if ( not pathExists) : 
		while ( not pathExists ) :
			print (LINKDIR + " does not exist.")
			LINKDIR = raw_input("Absolute PATH to the folder containing all the link files?\t") 
			pathExists = os.path.isdir(LINKDIR)
	explore(LINKDIR, listLinks)
	print (LINKDIR + " exists and contains "+ str(len(listLinks)) +" files, presumably all link files.\n")

	#Option file folder
	OPTDIR = raw_input("Absolute PATH to the folder that will contain all the option files?\t")
	pathExists = os.path.isdir(OPTDIR)
	if ( not pathExists) : 
		print ("\nCreating the folder "+OPTDIR+"\n")
		try:
			os.makedirs(OPTDIR)
		except OSError:
			if not os.path.isdir(OPTDIR):
				raise
	print (OPTDIR + " exists and will contain the option files.\n")

	#Result file folder
	RESDIR = raw_input("Absolute PATH to the folder that will contain all the result files?\t")
	pathExists = os.path.isdir(RESDIR)
	if ( not pathExists) : 
		print ("\nCreating the folder "+RESDIR+"\n")
		try:
			os.makedirs(RESDIR)
		except OSError:
			if not os.path.isdir(RESDIR):
				raise
	print (RESDIR + " exists and will contain the result files.\n")

	#Species tree
	TREEFILEGIVEN = raw_input("Do you want to provide a starting species tree (Newick)?\t").upper()
	while (not ( "YES".startswith(TREEFILEGIVEN) or "NO".startswith(TREEFILEGIVEN) ) ) :
		print (TREEFILEGIVEN + " is not an option I understand for whether you give a starting species tree or not.\n")
		TREEFILEGIVEN = raw_input("Do you want to provide a starting species tree (Newick)?\t").upper()
	if ( "YES".startswith(TREEFILEGIVEN) ):
		TREEFILE = raw_input("Absolute PATH to the starting species tree file (Newick)?\t")
		pathExists = os.path.isfile(TREEFILE)
		if ( not pathExists) : 
			while ( not pathExists ) :
				print (TREEFILE + " does not exist.\n")
				TREEFILE = raw_input("Absolute PATH to the starting species tree file (Newick)?\t") 
				pathExists = os.path.isfile(TREEFILE)
		print ("Starting from "+TREEFILE + " \n")
	else : 
		STARTINGTREE = raw_input("What kind of species tree? RANDOM or MRP?\t").upper()
		if (not ( "RANDOM".startswith(STARTINGTREE) or "MRP".startswith(STARTINGTREE) ) ): 
			while (not ( "RANDOM".startswith(STARTINGTREE) or "MRP".startswith(STARTINGTREE) ) ) :
				print (STARTINGTREE + " is not an option I understand for a starting species tree.\n")
				STARTINGTREE = raw_input("What kind of species tree? RANDOM or MRP?\t").upper()
		if ("RANDOM".startswith(STARTINGTREE)):
			print ("Starting species tree: RANDOM\n")
		else:
			print ("Starting species tree: MRP\n")
		

	#Do we optimize the species tree?
	TOPOSPECIES = raw_input("Do you want to optimize the species tree topology?\t").upper()
	while (not ( "YES".startswith(TOPOSPECIES) or "NO".startswith(TOPOSPECIES) ) ) :
		print (TOPOSPECIES + " is not an option I understand for whether you want to optimize the species tree topology or not.\n")
		TOPOSPECIES = raw_input("Do you want to optimize the species tree topology?\t").upper()
	if ("YES".startswith(TOPOSPECIES) ):
		print ("Optimizing the species tree topology.\n")
	else:
		print ("NOT optimizing the species tree topology.\n")

	#Do we optimize the duplication and loss parameters, and how?
	DLPARAM = raw_input("Do you want to optimize the duplication and loss parameters?\t").upper()
	while (not ( "YES".startswith(DLPARAM) or "NO".startswith(DLPARAM) ) ) :
		print (DLPARAM + " is not an option I understand for whether you want to optimize the duplication and loss parameters or not.\n")
		DLPARAM = raw_input("Do you want to optimize the duplication and loss parameters?\t").upper()
	if ("YES".startswith(DLPARAM) ):
			DLOPT = raw_input("How do you want to optimize the duplication and loss parameters? Branchwise or average ?\t").upper()
			while (not("BRANCHWISE".startswith(DLOPT) or "AVERAGE".startswith(DLOPT) ) ) :
				print ( DLOPT +" is not an option I understand for the type of duplication and loss parameters optimization. Should be branchwise or average.\n")
				DLOPT = raw_input("How do you want to optimize the duplication and loss parameters? Branchwise or average ?\t").upper()			
			if ("BRANCHWISE".startswith(DLOPT)):
				print ("Optimizing the duplication and loss parameters, using BRANCHWISE.\n")
			else:
				print ("Optimizing the duplication and loss parameters, using AVERAGE.\n")
	else:
		print ("NOT optimizing the duplication and loss parameters.\n")

	#Do we want to assume that all species have roughly the same number of genes, 
	#and use that for computing an expected amount of missing data?
	EQUGENOMES = raw_input("Do you want to assume that all genomes, if they had no missing data, would have the same number of genes?\t").upper()
	while (not ( "YES".startswith(EQUGENOMES) or "NO".startswith(EQUGENOMES) ) ) :
		print (EQUGENOMES + " is not an option I understand for whether you want to assume that all genomes, save missing data, would have the same number of genes.\n")
		EQUGENOMES = raw_input("Do you want to assume that all genomes, if they had no missing data, would have the same number of genes?\t").upper()
	if ("YES".startswith(EQUGENOMES) ):
		print ("Assuming that all the genomes have the same number of genes, so that all missing genes are artifactual.\n")
	else:
		print ("NOT assuming that all the genomes have the same number of genes.\n")
		
	#Do we want to provide gene trees?
	INPUTGENETREES = raw_input("Do you want to provide input gene trees?\t").upper()
	listTrees = list()
        while (not ( "YES".startswith(INPUTGENETREES) or "NO".startswith(INPUTGENETREES) ) ) :
                print (INPUTGENETREES + " is not an option I understand for whether you want to provide starting gene trees or not.\n")
                INPUTGENETREES = raw_input("Do you want to provide input gene trees?\t").upper()
        if ("YES".startswith(INPUTGENETREES) ):
                print ("User-provided gene trees. \n")
                INPUTGENETREEDIR = raw_input("ABSOLUTE PATH to the folder containing starting gene trees?\t")
		while ( not os.path.isdir(INPUTGENETREEDIR) ) :
			print (INPUTGENETREEDIR + " does not exist.\n")
			INPUTGENETREEDIR = raw_input("Absolute PATH to the folder containing all the alignment files?\t") 
		explore(INPUTGENETREEDIR, listTrees)
		print (INPUTGENETREEDIR + " exists and contains "+str(len(listTrees))+" files, presumably all gene trees.\n")
        else:
                print ("Initial gene trees constructed using BioNJ.\n")


	#Do we optimize the gene trees?
	TOPOGENE = raw_input("Do you want to optimize the gene tree topologies?\t").upper()
	while (not ( "YES".startswith(TOPOGENE) or "NO".startswith(TOPOGENE) ) ) :
		print (TOPOGENE + " is not an option I understand for whether you want to optimize the gene tree topologies or not.\n")
		TOPOGENE = raw_input("Do you want to optimize the gene tree topologies?\t").upper()
	if ("YES".startswith(TOPOGENE) ):
		print ("Optimizing the gene tree topologies.\n")
	else:
		print ("NOT optimizing the gene tree topologies (but they will still be rerooted).\n")
	
	#Time limit
	TIMELIMIT = raw_input("What time limit do you want to set, in hours (minimum 2h)?\t")
	while (not TIMELIMIT.isdigit()):
		print (TIMELIMIT+" does not look like a number of hours (decimal numbers are not allowed).")
		TIMELIMIT = raw_input("What time limit do you want to set, in hours (minimum 2h)?\t")
	print ("We will stop at the latest about 1 hour before "+ TIMELIMIT +".\n")

	###########################################
	###########################################
	#Creating gene family-specific option files
	listSpecies = list()
	listOptionFiles = list()
	listSizes = list()
	dictLinks=dict()
	print ("Now we create one option file per gene family, but we do so only if we find an alignment file and a link file that correspond to each other.\n")
	print ("Correspondance is based on the radical of the file names (file name without .extension).\n")
	listAlns.sort()
	listLinks.sort()
	dictTrees = dict()
	if ("YES".startswith(INPUTGENETREES) ):
		listTrees.sort()
		for treef in listTrees:
			dictTrees[os.path.basename(treef).split('.')[0]] = treef
	for link in listLinks:
		dictLinks[os.path.basename(link).split('.')[0]] = link
	for aln in listAlns:
		radical = os.path.basename(aln).split('.')[0]
		if (dictLinks.__contains__(radical)):
			#First, open the link file in order to build a list of species
			try:
				fin=open(dictLinks[radical], 'r')
			except IOError, e:
				print "Unknown file: ",dictLinks[radical]
				sys.exit()
			for l in fin:
				num=1
				sp = l.split(":")[0]
				if "," in l:
					num=len(l.split(";"))
				for i in range(num):
					listSpecies.append(sp)
			fin.close()
			try:
				fopt=open(os.path.join(OPTDIR,radical)+'.opt', 'w')
			except IOError, e:
				print "Unknown file: ",os.path.join(OPTDIR,radical)+'.opt'
				sys.exit()
			fopt.write("\n######## First, data files ########\n")
			fopt.write("PATH="+os.path.join(SEQDIR,"")+"\n")
			fopt.write("RESULT="+os.path.join(RESDIR,"")+"\n")
			fopt.write("DATA="+radical+"\n")
			fopt.write("taxaseq.file="+dictLinks[radical]+"\n")
			fopt.write("input.sequence.file="+aln+"\n")
			if ("FASTA".startswith(DATAFORMAT) ):  
				fopt.write("input.sequence.format=Fasta\n")
			elif ("MASE".startswith(DATAFORMAT)):
				fopt.write("input.sequence.format=Fasta\n")
			elif ("PHYLIP".startswith(DATAFORMAT)):
				fopt.write("input.sequence.format=Phylip\n")
			elif ("CLUSTAL".startswith(DATAFORMAT)):
				fopt.write("input.sequence.format=Clustal\n")
			elif ("DCSE".startswith(DATAFORMAT)):
				fopt.write("input.sequence.format=Dcse\n")
			elif ("NEXUS".startswith(DATAFORMAT)):
				fopt.write("input.sequence.format=Nexus\n")
			listSizes.append(str(os.stat( aln )[6])	)		
			fopt.write("output.reconciled.tree.file=$(RESULT)$(DATA).ReconciledTree\n")
			fopt.write("output.duplications.tree.file=$(RESULT)$(DATA).DuplicationTree\n")
			fopt.write("output.losses.tree.file=$(RESULT)$(DATA).LossTree\n")
			fopt.write("output.numbered.tree.file=$(RESULT)$(DATA).NumberedTree\n")
			fopt.write("input.sequence.sites_to_use=all\n")
			fopt.write("input.sequence.max_gap_allowed=100%\n")
			if ("YES".startswith(INPUTGENETREES)):
				fopt.write("init.gene.tree=user\n") #the program will use a user-provided input gene tree
				fopt.write("gene.tree.file="+dictTrees[radical]+"\n")
                        else:
				fopt.write("init.gene.tree=bionj\n") #the program will build a starting user gene tree
				fopt.write("gene.tree.file=$(RESULT)$(DATA).GeneTree\n")
			fopt.write("output.starting.gene.tree.file=$(RESULT)$(DATA).StartingTree\n")
			fopt.write("\n######## Second, model options ########\n")
			if ("DNA".startswith(DATATYPE) ):
				fopt.write("alphabet=DNA\n")
				fopt.write("model=GTR(a=1.17322, b=0.27717, c=0.279888, d=0.41831, e=0.344783, theta=0.523374, theta1=0.542411, theta2=0.499195)\n")
				fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
			elif ("RNA".startswith(DATATYPE) ):
				fopt.write("alphabet=RNA\n")
				fopt.write("model=GTR( initFreqs=observed )\n")
				fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
			elif ("CODON".startswith(DATATYPE) ):
				fopt.write("alphabet=Codon(letter=DNA)\n")
				fopt.write("input.sequence.remove_stop_codons = yes\n")
				fopt.write("genetic_code=Standard\n")
				fopt.write("model=YN98( kappa=1, omega=1.0, initFreqs=observed )\n")
			elif ("PROTEIN".startswith(DATATYPE) ):
				fopt.write("alphabet=Protein\n")
				fopt.write("model=LG08+F(initFreqs=observed )\n")
				fopt.write("rate_distribution=Gamma(n=4,alpha=1)\n")
			#fopt.write("optimization.ignore_parameter=dist_Gamma.alpha, GTR.a, GTR.b, GTR.c, GTR.d, GTR.e, GTR.theta, GTR.theta1, GTR.theta2\n")
			#fopt.write("\n######## Then, algorithm options ########\n")
			#fopt.write("heuristics.level=0\n")
			#fopt.write("species.id.limit.for.root.position=3 #Useless unless heuristics.level=1\n") 
			fopt.write("\n######## Finally, optimization options ########\n")
			if ("YES".startswith(TOPOGENE) ):
				fopt.write("optimization.topology=yes\n")
			else:
				fopt.write("optimization.topology=no\n")			 
			fopt.write("optimization.topology.algorithm_nni.method=fast\n")
			fopt.write("optimization.tolerance=0.01\n")
			fopt.write("optimization.method_DB.nstep=0\n")
			fopt.write("optimization.topology.numfirst=false\n")
			fopt.write("optimization.topology.tolerance.before=100\n")
			fopt.write("optimization.topology.tolerance.during=100\n")
			fopt.write("optimization.max_number_f_eval=1000000\n")
			fopt.write("optimization.final=none\n")
			fopt.write("optimization.verbose=0\n")
			fopt.write("optimization.message_handler=none\n")
			fopt.write("optimization.profiler=none\n")
			fopt.write("optimization.reparametrization=no\n")
			fopt.close()
			listOptionFiles.append(os.path.join(OPTDIR,radical)+'.opt')
		else:
			print ("Alignment "+ aln+ " does not have a corresponding link file. Skipping it.\n")
	#Now the gene option files have been created

	#Creating the list of gene option files
	try:
		fout=open(os.path.join(OPTDIR,"listGenes.txt"), 'w')
	except IOError, e:
		print "Unknown file: ", os.path.join(OPTDIR,"listGenes.txt")
		sys.exit()
	for i in range(len(listOptionFiles)):
		fout.write(listOptionFiles[i]+":"+ listSizes[i] +"\n")
	fout.close()

	#Creating the species list file.
	dictSpecies = dict()
	for i in listSpecies:
		if (dictSpecies.__contains__(i)):
			dictSpecies[i] = dictSpecies[i]+1
		else:
			dictSpecies[i]=1
			
	#We have some data, why not output it?
	print ("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	print ("\n\tTaxonomic distribution in the gene families for which we created option files:\n")
	max = 0 
	min = 1000000
	sum = 0
	try:
	   	fsp=open(os.path.join(OPTDIR,"listSpecies.txt"), 'w')
   	except IOError, e:
	   	print "Unknown file: ",os.path.join(OPTDIR,"listSpecies.txt")
	   	sys.exit()
	for (k, v) in dictSpecies.items():
		print (k + " : " + str(v) + " genes.")
		fsp.write(k+"\n")
		sum = sum +v
		if (v>max):
			max=v
		if (v<min):
			min=v
	fsp.close()
	print ("\n\tAverage number of genes: "+str(sum/len(dictSpecies))+" ; Maximum: "+str(max)+" ; Minimum: "+str(min)+"\n")


	#Creating a taxonomic coverage file	
	if ("YES".startswith(EQUGENOMES) ):
		try:
			fsp=open(os.path.join(OPTDIR,"listSpeciesWithSequenceCoverage.txt"), 'w')
		except IOError, e:
			print "Unknown file: ",os.path.join(OPTDIR,"listSpeciesWithSequenceCoverage.txt")
			sys.exit()
		for k, v in dictSpecies.items():
			fsp.write (k + " : " + str(100*v/max) + " \n")
		fsp.close()

	
	#Now, the general options.
	try:
	   	fopt=open(os.path.join(OPTDIR,"GeneralOptions.txt"), 'w')
   	except IOError, e:
	   	print "Unknown file: ",os.path.join(OPTDIR,"GeneralOptions.txt")
	   	sys.exit()
	fopt.write("\n######## First, data files ########\n")
	fopt.write("OPT="+os.path.join(OPTDIR,"")+"\n")
	fopt.write("RESULT="+os.path.join(RESDIR,"")+"\n")
	fopt.write("DATA="+radical+"\n")
	if ("YES".startswith(TREEFILEGIVEN)):
		fopt.write("init.species.tree=user\n")
		fopt.write("species.tree.file="+TREEFILE+"\n") 
	elif ("RANDOM".startswith(STARTINGTREE)):
		fopt.write("init.species.tree=random #user\n")
	else:
		fopt.write("init.species.tree=mrp\n")
	fopt.write("species.names.file=$(OPT)listSpecies.txt\n")
	fopt.write("starting.tree.file=$(RESULT)StartingTree.tree\n")
	fopt.write("output.tree.file=$(RESULT)OutputSpeciesTree.tree\n")
	#There follows a list of genes for which a tree needs to be built.
	fopt.write("genelist.file=$(OPT)listGenes.txt\n")
	fopt.write("output.duplications.tree.file=$(RESULT)OutputSpeciesTree_ConsensusDuplications.tree\n")
	fopt.write("output.losses.tree.file=$(RESULT)OutputSpeciesTree_ConsensusLosses.tree\n")
	fopt.write("output.numbered.tree.file=$(RESULT)OutputSpeciesTree_ConsensusNumbered.tree\n")
	fopt.write("\n######## Second, options ########\n")
	if ("YES".startswith(TOPOSPECIES)):
		fopt.write("optimization.topology=yes\n")
	else:
		fopt.write("optimization.topology=no\n")
	#fopt.write("species.id.limit.for.root.position=3\n")
	if ("YES".startswith(DLPARAM)):
		if ("BRANCHWISE".startswith(DLOPT) ):
			fopt.write("branchProbabilities.optimization=average_then_branchwise\n")
		else:
			fopt.write("branchProbabilities.optimization=average\n")		
	else:
		fopt.write("branchProbabilities.optimization=no")
	if ("YES".startswith(EQUGENOMES) ):
		fopt.write("genome.coverage.file=$(PATH)HomolensSpeciesSequenceCoverage\n")
	fopt.write("spr.limit=5\n")
	fopt.write("time.limit="+TIMELIMIT+"\n")
	fopt.close()
	
	print ("\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	print ("\n\n\n\tTo launch PHYLDOG, you will need to run something like:\n\n")
	print ("\tmpirun -np NUMBER_OF_PROCESSES PATH_TO_PHYLDOG/phyldog param="+ os.path.join(OPTDIR,"GeneralOptions.txt") +"\n\n")

	print ("\n\tThank you for using this PHYLDOG file preparation script!\n\n")

	
		
