import sys

especes_actuelles = {}
degres = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
file_genes = open("./deco_genes_reconciled").readlines()
file_adj = open("./deco_adjacences_reconciled").readlines()
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

for name in ["deco_reconciled","deco_starting"]:
	print name
	fichier_especes = open("RESULTS/" + name+"_OUTPUT_species","r").readlines()

	especes_ancestrales = {}
	especes_actuelles = {}

	for line in fichier_especes:
		if len(line.split()) > 2:
			especes_ancestrales[line.split()[0]] = []
		else:
			especes_actuelles[line.split()[0]] = []
			#print line.split()[0]
			
	total = 0
	fichier_genes = open("RESULTS/" + name+"_OUTPUT_genes","r").readlines()
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
		
	fichier_adjacences = open("RESULTS/" + name+"_OUTPUT_adjacencies","r").readlines()
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
	degre_par_espece = dict()
	for e in especes_ancestrales.keys():
		degre_par_espece[e] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
		for g in especes_ancestrales[e]:
			if len(degres) > len(graph[e][g]):
				degres[len(graph[e][g])] = degres[len(graph[e][g])]  + 1
				degre_par_espece[e][len(graph[e][g])] += 1
	sortie = open("summary_degres_"+name,"w")
	for d in range(len(degres)):
		sortie.write(str(d)+" "+str(degres[d])+"\n")
	sortie.close()
	sortie = open("summary_degrees_per_ancspecies_"+name,"w")
        for e in degre_par_espece.keys():
                sortie.write(str(e) + " " +str(degre_par_espece[e][2]/float(sum(degre_par_espece[e])))+"\n")
	sortie.close()
	
