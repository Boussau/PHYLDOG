/*
 *  GeneTreeAlgorithms.cpp
 *  ReconcileDuplications.proj
 *
 *  Created by boussau on 29/06/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */

#include "GeneTreeAlgorithms.h"

/**************************************************************************
 * This function creates a sequence tree from a species tree and a std::map 
 * containing the link between the species and their sequences.
 **************************************************************************/

TreeTemplate<Node> * buildARandomSequenceTreeFromASpeciesTree (std::map <std::string, 
                                                               std::deque<std::string> > & spSeqs, 
                                                               TreeTemplate<Node> & tree, 
                                                               std::map <std::string, std::string> & spSelectedSeq) {
  
	TreeTemplate<Node> * seqTree ;
	seqTree = new TreeTemplate<Node>(tree);
	//seqTree = tree;
	spSelectedSeq.clear();
  
	for(std::map<std::string, std::deque<std::string> >::iterator it = spSeqs.begin(); it != spSeqs.end(); it++){
		int numSeq = (it->second).size();
		//    std::cout <<"Numseqs : "<<numSeq<<std::endl;
		if (numSeq > 0) {
      std::string chosenSequence = (it->second).at(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(numSeq));
			int leafId = seqTree->getLeafId(it->first);
			seqTree->setNodeName(leafId, chosenSequence);
			spSelectedSeq.insert( make_pair(it->first,chosenSequence));
		}
		else {
			int leafId = seqTree->getLeafId(it->first);
			seqTree->getNode(leafId)->getFather()->removeSon(seqTree->getNode(leafId));
		}
	}
	return seqTree;
  
}


/**************************************************************************
 * This function creates a sequence tree from a species tree and the std::map 
 * containing the link between the species and the putative orthologous sequence.
 **************************************************************************/

TreeTemplate<Node> * buildASequenceTreeFromASpeciesTreeAndCorrespondanceMap (TreeTemplate<Node> & tree, 
                                                                             std::map <std::string, 
                                                                             std::string> & spSelectedSeq) {
  
	TreeTemplate<Node> * seqTree ;
	seqTree = new TreeTemplate<Node>(tree);
  
	std::vector <std::string> names = seqTree->getLeavesNames();
	for(std::vector < std::string >::iterator it = names.begin(); it != names.end(); it++){
		if (spSelectedSeq.count(*it)>0) {
      std::map < std::string, std::string >::iterator iter = spSelectedSeq.find(*it);
      std::string chosenSequence = iter->second;
			int leafId = seqTree->getLeafId(*it);
			seqTree->setNodeName(leafId, chosenSequence);
		}
		else {
			int leafId = seqTree->getLeafId(*it);
			seqTree->getNode(leafId)->getFather()->removeSon(seqTree->getNode(leafId));
		}
	}
  
	return seqTree;
  
}

/******************************************************************************/
// This function refines branch lengths of a gene tree.
/******************************************************************************/

void refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                             TreeTemplate<Node>  *unrootedGeneTree, 
                                                             VectorSiteContainer * sites, 
                                                             SubstitutionModel* model, 
                                                             DiscreteDistribution* rDist, 
                                                             string file, Alphabet *alphabet)
{ 
  std::string backupParamnumEval = params[ std::string("optimization.max_number_f_eval")];
  std::string backupParamOptTol = params[ std::string("optimization.tolerance")];
  std::string backupParamOptTopol = params[ std::string("optimization.topology")];
  params[ std::string("optimization.tolerance")] = "0.00001";
  params[ std::string("optimization.max_number_f_eval")] = "100000";
  params[ std::string("optimization.topology")] = "no";
  DiscreteRatesAcrossSitesTreeLikelihood* tl;
  //DRHomogeneousTreeLikelihood * 
  tl = new DRHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, true, true);
  tl->initialize();//Only initializes the parameter list, and computes the likelihood
  double logL = tl->getValue();
  if (std::isinf(logL))
    {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for (unsigned int i = 0; i < pl.size(); i++)
      {
      if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
      }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
    }
  ApplicationTools::displayResult("Initial log likelihood for family", file+": "+TextTools::toString(-logL, 15)); 
  if (std::isinf(logL))
    {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    CodonAlphabet *pca=dynamic_cast<CodonAlphabet*>(alphabet);
    if (pca){
      bool f=false;
      unsigned int  s;
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++){
        if (std::isinf(tl->getLogLikelihoodForASite(i))){
          const Site& site=sites->getSite(i);
          s=site.size();
          for (unsigned int j=0;j<s;j++)
            if (pca->isStop(site.getValue(j))){
              (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
              f=true;
            }
        }
      }
      if (f)
        exit(-1);
    }
    ApplicationTools::displayError("!!! Looking at each site:");
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
    }
  tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
                                                             PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params));
  
  params[ std::string("optimization.tolerance")] = backupParamOptTol;
  params[ std::string("optimization.max_number_f_eval")] = backupParamnumEval;
  params[ std::string("optimization.topology")] = backupParamOptTopol;
  //delete unrootedGeneTree;
  unrootedGeneTree = new TreeTemplate<Node> ( tl->getTree() );
  delete tl;
}

/******************************************************************************/
// This function builds a bionj tree
/******************************************************************************/

TreeTemplate<Node>  * buildBioNJTree (std::map<std::string, std::string> & params, 
                                      SiteContainer * sites, 
                                      SubstitutionModel* model, 
                                      DiscreteDistribution* rDist, 
                                      string file, 
                                      Alphabet *alphabet) {
  TreeTemplate<Node>  *unrootedGeneTree = 0;
  
  DistanceEstimation distEstimation(model, rDist, sites, 1, false);
  
  BioNJ * bionj = new BioNJ();     
  
  bionj->outputPositiveLengths(true);  
  
  std::string type = ApplicationTools::getStringParameter("bionj.optimization.method", params, "init");
  if(type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
  else if(type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
  else if(type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
  else throw Exception("Unknown parameter estimation procedure '" + type + "'.");
  // Should I ignore some parameters?
  ParameterList allParameters = model->getParameters();
  allParameters.addParameters(rDist->getParameters());
  ParameterList parametersToIgnore;
  std::string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", "", true, false);
  bool ignoreBrLen = false;
  StringTokenizer st(paramListDesc, ",");
  
  while(st.hasMoreToken())
    {
    try
      {
      std::string param = st.nextToken();
      if(param == "BrLen")
        ignoreBrLen = true;
      else
        {
        if (allParameters.hasParameter(param))
          {
          Parameter* p = &allParameters.getParameter(param);
          parametersToIgnore.addParameter(*p);
          }
        else ApplicationTools::displayWarning("Parameter '" + param + "' not found."); 
        }
      } 
    catch(ParameterNotFoundException pnfe)
      {
      ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
      }
    }
  double tolerance = ApplicationTools::getDoubleParameter("bionj.optimization.tolerance", params, .000001);
  
  unrootedGeneTree = OptimizationTools::buildDistanceTree(distEstimation, *bionj, parametersToIgnore, !ignoreBrLen, false, type, tolerance);
  std::vector<Node*> nodes = unrootedGeneTree->getNodes();
  
  for(unsigned int k = 0; k < nodes.size(); k++)
    {
    if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
    if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() >= 5) throw Exception("Found a very large branch length in a gene tree; will avoid this family.");
    }
  
  delete bionj;  
  return unrootedGeneTree;
}


/******************************************************************************/
// This function refines a gene tree topology and branch lengths using the PhyML 
// algorithm.
/******************************************************************************/

void refineGeneTreeUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                TreeTemplate<Node>  *unrootedGeneTree, 
                                                VectorSiteContainer * sites, 
                                                SubstitutionModel* model, 
                                                DiscreteDistribution* rDist, 
                                                string file, 
                                                Alphabet *alphabet)
{ 
  std::string backupParamOptTopo = params[ std::string("optimization.topology.algorithm_nni.method")];
  std::string backupParamnumEval = params[ std::string("optimization.max_number_f_eval")];
  std::string backupParamOptTol = params[ std::string("optimization.tolerance")];
  params[ std::string("optimization.topology.algorithm_nni.method")] = "phyml";
  params[ std::string("optimization.tolerance")] = "0.00001";
  params[ std::string("optimization.max_number_f_eval")] = "100000";
  NNIHomogeneousTreeLikelihood * tl = new NNIHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, true, true);
  tl->initialize();//Only initializes the parameter list, and computes the likelihood
  double logL = tl->getValue();
  if (std::isinf(logL))
    {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for (unsigned int i = 0; i < pl.size(); i++)
      {
      if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
      }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
    }
  ApplicationTools::displayResult("Initial log likelihood for family", file+": "+TextTools::toString(-logL, 15)); 
  if (std::isinf(logL))
    {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    CodonAlphabet *pca=dynamic_cast<CodonAlphabet*>(alphabet);
    if (pca){
      bool f=false;
      unsigned int  s;
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++){
        if (std::isinf(tl->getLogLikelihoodForASite(i))){
          const Site& site=sites->getSite(i);
          s=site.size();
          for (unsigned int j=0;j<s;j++)
            if (pca->isStop(site.getValue(j))){
              (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
              f=true;
            }
        }
      }
      if (f)
        exit(-1);
    }
    ApplicationTools::displayError("!!! Looking at each site:");
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
    }
  
  tl = dynamic_cast<NNIHomogeneousTreeLikelihood*>(PhylogeneticsApplicationTools::optimizeParameters(tl, 
                                                                                                     tl->getParameters(), 
                                                                                                     params));
  
  params[ std::string("optimization.topology.algorithm_nni.method")] = backupParamOptTopo ;
  params[ std::string("optimization.tolerance")] = backupParamOptTol;
  params[ std::string("optimization.max_number_f_eval")] = backupParamnumEval;
  //delete unrootedGeneTree;
  unrootedGeneTree = new TreeTemplate<Node> ( tl->getTree() );
  delete tl;
}



