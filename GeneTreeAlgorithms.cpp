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

double refineGeneTreeBranchLengthsUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                             TreeTemplate<Node>  *& unrootedGeneTree, 
                                                             VectorSiteContainer * sites, 
                                                             SubstitutionModel* model, 
                                                             DiscreteDistribution* rDist, 
                                                             string file, Alphabet *alphabet, bool mapping)
{ 
    DiscreteRatesAcrossSitesTreeLikelihood* tl;
    //DRHomogeneousTreeLikelihood * DRHomogeneousTreeLikelihood
    // std::cout << TreeTemplateTools::treeToParenthesis(*unrootedGeneTree) <<std::endl;

    tl = new DRHomogeneousTreeLikelihood(*unrootedGeneTree, *sites, model, rDist, true);  
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
          if (f) {
              MPI::COMM_WORLD.Abort(1);
          exit(-1);
          }
      }
      ApplicationTools::displayError("!!! Looking at each site:");
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
        {
        (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
        }
      ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
          MPI::COMM_WORLD.Abort(1);
      exit(-1);
      }

  if (!mapping) 
    {
    
    tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
                                                               PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), tl->getParameters(), params));

    }
  else {
    optimizeBLMapping(dynamic_cast<DRTreeLikelihood*>(tl), 10);
    tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
                                                             PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), tl->getParameters(), params));

    
  }

//  std::cout << "Sequence likelihood: "<< -logL << "; Optimized sequence log likelihood "<< -tl->getValue()  << std::endl;
  logL = -tl->getValue();
  if (unrootedGeneTree)
    delete unrootedGeneTree;
  unrootedGeneTree = new TreeTemplate<Node> ( tl->getTree() );
  delete tl;   
  return logL;
}


/******************************************************************************/
// This function maps all types of substitutions in a gene tree.
/******************************************************************************/

vector< vector<unsigned int> > getCountsPerBranch(
                                                  DRTreeLikelihood& drtl,
                                                  const vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg,
                                                  SubstitutionCount *count,
                                                  bool stationarity,
                                                  double threshold)
{
  //auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  
  //SubstitutionCount* count = new SimpleSubstitutionCount(reg);
/*  const CategorySubstitutionRegister* creg = 0;
  if (!stationarity) {
    try {
      creg = &dynamic_cast<const CategorySubstitutionRegister&>(reg);
    } catch (Exception& ex) {
      throw Exception("The stationarity option can only be used with a category substitution register.");
    }
  }

  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));
*/
  vector< vector<unsigned int> > counts(ids.size());/*
  unsigned int nbSites = mapping->getNumberOfSites();
  unsigned int nbTypes = mapping->getNumberOfSubstitutionTypes();
  for (size_t k = 0; k < ids.size(); ++k) {
    //vector<double> countsf = SubstitutionMappingTools::computeSumForBranch(*mapping, mapping->getNodeIndex(ids[i]));
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    unsigned int nbIgnored = 0;
    bool error = false;
    for (unsigned int i = 0; !error && i < nbSites; ++i) {
      double s = 0;
      for (unsigned int t = 0; t < nbTypes; ++t) {
        tmp[t] = (*mapping)(k, i, t);
        error = isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
    }

  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (unsigned int t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }

      if (!stationarity) {
        vector<double> freqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(drtl, ids[k]);
        //Compute frequencies for types:
        vector<double> freqsTypes(creg->getNumberOfCategories());
        for (size_t i = 0; i < freqs.size(); ++i) {
          unsigned int c = creg->getCategory(static_cast<int>(i));
          freqsTypes[c - 1] += freqs[i];
        }
        //We divide the counts by the frequencies and rescale:
        double s = VectorTools::sum(countsf);
        for (unsigned int t = 0; t < nbTypes; ++t) {
          countsf[t] /= freqsTypes[creg->getCategoryFrom(t + 1) - 1];
        }
        double s2 = VectorTools::sum(countsf);
        //Scale:
        countsf = (countsf / s2) * s;
      }
    }

    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j) {
      counts[k][j] = static_cast<unsigned int>(floor(countsf[j] + 0.5)); //Round counts
    }
  }*/
  return counts;
}


/******************************************************************************/
// This function maps total numbers of substitutions per branch in a gene tree.
/******************************************************************************/

vector< vector<unsigned int> > getTotalCountsOfSubstitutionsPerBranch(
                                                  DRTreeLikelihood& drtl,
                                                  const vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg,
                                                  SubstitutionCount *count,
                                                  double threshold)
{
  //auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  
  //SubstitutionCount* count = new SimpleSubstitutionCount(reg);  
/*  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, *count, false));
 */ 
  vector< vector<unsigned int> > counts(ids.size());
 /* unsigned int nbSites = mapping->getNumberOfSites();
  unsigned int nbTypes = mapping->getNumberOfSubstitutionTypes();
  for (size_t k = 0; k < ids.size(); ++k) {
    //vector<double> countsf = SubstitutionMappingTools::computeSumForBranch(*mapping, mapping->getNodeIndex(ids[i]));
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    unsigned int nbIgnored = 0;
    bool error = false;
    for (unsigned int i = 0; !error && i < nbSites; ++i) {
      double s = 0;
      for (unsigned int t = 0; t < nbTypes; ++t) {
        tmp[t] = (*mapping)(k, i, t);
        error = isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
    }
    
  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (unsigned int t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
    }
    
    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j) {
      counts[k][j] = static_cast<unsigned int>(floor(countsf[j] + 0.5)); //Round counts
    }
  }*/
  return counts;
}


/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/

void optimizeBLMapping(
                       DRTreeLikelihood* tl,
                       double precision)
{
 //   std::cout<< "BEFORE: "<<TreeTemplateTools::treeToParenthesis(tl->getTree(), true)<< std::endl;
/*
  double currentValue = tl->getValue();
  double newValue = currentValue - 2* precision;
  ParameterList bls = tl->getBranchLengthsParameters () ;
  ParameterList newBls = bls;
  vector< vector<unsigned int> > counts;
  double numberOfSites = (double) tl->getNumberOfSites();
  vector<int> ids = tl->getTree().getNodesId();
  ids.pop_back(); //remove root id.
  SubstitutionRegister* reg = 0;  
  //Counting all substitutions
  reg = new TotalSubstitutionRegister(tl->getAlphabet());
  bool first = true;
    SubstitutionCount *count = new SimpleSubstitutionCount( reg);  // new UniformizationSubstitutionCount(tl->getSubstitutionModel(0,0), reg);   
  while (currentValue > newValue + precision) {
    if (first)
      first=false;
    else {
      currentValue = newValue;
    }
    //Perform the mapping:
    counts = getTotalCountsOfSubstitutionsPerBranch(*tl, ids, tl->getSubstitutionModel(0,0), *reg, count, -1);
    double value;
    string name;
    for (unsigned int i = 0 ; i < counts.size() ; i ++) {
      name = "BrLen" + TextTools::toString(i);
      value = double(VectorTools::sum(counts[i])) / (double)numberOfSites;
      newBls.setParameterValue(name, newBls.getParameter(name).getConstraint()->getAcceptedLimit (value));
    }
    
    tl->matchParametersValues(newBls);
    
    newValue = tl->getValue();
  // std::cout <<"before: "<< currentValue<<" AFTER: "<< newValue <<std::endl;
      if (currentValue > newValue + precision) { //Significant improvement
      bls = newBls;
        tl->setParameters(bls);
 //       std::cout<< "AFTER: "<<TreeTemplateTools::treeToParenthesis(tl->getTree(), true)<< std::endl;

    }
    else { 
      if (currentValue < newValue) //new state worse, getting back to the former state
      tl->matchParametersValues(bls);
    }
  }
  */
  
  /*
  //Counting all substitutions
  reg = new ComprehensiveSubstitutionRegister(tl->getAlphabet(), false);
  bool first = true;
  SubstitutionCount *count = new UniformizationSubstitutionCount(tl->getSubstitutionModel(0,0), reg->clone());
  while (currentValue > newValue + precision) {
    if (first)
      first=false;
      else {
        currentValue = newValue;
      }
    //Perform the mapping:
    counts = getCountsPerBranch(*tl, ids, tl->getSubstitutionModel(0,0), *reg, count, stationarity, -1);
    double value;
    string name;
    for (unsigned int i = 0 ; i < counts.size() ; i ++) {
      name = "BrLen" + TextTools::toString(i);
      value = double(VectorTools::sum(counts[i])) / (double)numberOfSites;
      newBls.setParameterValue(name, newBls.getParameter(name).getConstraint()->getAcceptedLimit (value));
    }

    tl->matchParametersValues(newBls);

    newValue = tl->getValue();
    std::cout <<"before: "<< currentValue<<" AFTER: "<< newValue <<std::endl;
    if (currentValue > newValue + precision) { //Significant improvement
      bls = newBls;
    }
    else { //getting back to the former state
      tl->matchParametersValues(bls);
    }
  }*/
}


/******************************************************************************/
// This function optimizes branch lengths in a gene tree using substitution mapping
/******************************************************************************/

void optimizeBLMappingForSPRs(
                       DRTreeLikelihood* tl,
                       double precision, map<string, string> params)
{    
 /*   if (!tl->isInitialized () ) {
        tl->initialize();
    }
    double currentValue = tl->getValue();
    double newValue = currentValue - 2* precision;
    ParameterList bls = tl->getBranchLengthsParameters () ;
    ParameterList newBls = bls;
    vector< vector<unsigned int> > counts;
    double numberOfSites = (double) tl->getNumberOfSites();
    vector<int> ids = tl->getTree().getNodesId();
    ids.pop_back(); //remove root id.
    SubstitutionRegister* reg = 0;  

    //Counting all substitutions
    reg = new TotalSubstitutionRegister(tl->getAlphabet());
    
    //Then, normal optimization.
    
    
    SubstitutionCount *count = new SimpleSubstitutionCount( reg);  // new UniformizationSubstitutionCount(tl->getSubstitutionModel(0,0), reg);   
   //We do the mapping based thing only once:
        //Perform the mapping:
    
        counts = getTotalCountsOfSubstitutionsPerBranch(*tl, ids, tl->getSubstitutionModel(0,0), *reg, count, -1);

        double value;
        string name;
        for (unsigned int i = 0 ; i < counts.size() ; i ++) {
          //  if (abs(bls[i].getValue() - 0.1) < 0.000001) {
            value = double(VectorTools::sum(counts[i])) / (double)numberOfSites;
             // LITTLE ATTEMPT TO TRY TO AVOID HUGE BRANCHES THAT MAY OCCUR FOR VERY BAD ALIGNMENTS
            if (value > 1.0) value = 1.0;
            name = "BrLen" + TextTools::toString(i);
            newBls.setParameterValue(name, newBls.getParameter(name).getConstraint()->getAcceptedLimit (value));
        }

       // tl->matchParametersValues(newBls);
    tl->setParametersValues(newBls);
    
        newValue = tl->getValue();

        if (currentValue > newValue + precision) { //Significant improvement
            bls = newBls;
            tl->setParametersValues(bls);     
        }
        else { 
            if (currentValue < newValue) //new state worse, getting back to the former state
                tl->matchParametersValues(bls);
            // tl->getValue();

        }

  //  }
    //Then, normal optimization.

//ATTEMPT WITHOUT FULL OPTIMIZATION 16 10 2011   
    //But with few rounds of optimization
    
    int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
    {
        params[ std::string("optimization.max_number_f_eval")] = 10;
    }
        
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params, "", true, false);
    params[ std::string("optimization.max_number_f_eval")] = backup;
     */
    
}

/******************************************************************************/
// This function optimizes branch lengths in a gene tree without substitution mapping
/******************************************************************************/

void optimizeBLForSPRs(
                              DRTreeLikelihood* tl,
                              double precision, map<string, string> params)
{    
    //But with few rounds of optimization
    int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
    {
        params[ std::string("optimization.max_number_f_eval")] = 10;
    }
    
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), params, "", true, false);
    params[ std::string("optimization.max_number_f_eval")] = backup;    
    
}


/******************************************************************************/
// This function optimizes branch lengths in a gene tree without substitution mapping
/******************************************************************************/

void optimizeBranchLengthsFast( DRHomogeneousTreeLikelihood* tl,
                               const ParameterList& parameters,
                               double tolerance)
throw (Exception)
{
    TreeTemplate<Node> tree (tl->getTree());
    std::vector<Node*> nodes = TreeTemplateTools::getNodes(*(tree.getRootNode()) );
    double oldLk = tl->getValue();
    double Lk = oldLk / 2;
    const DRASDRTreeLikelihoodNodeData* fatherData;
    const DRASDRTreeLikelihoodNodeData* sonData;
     VVVdouble                   fatherArray ;
     VVVdouble                   sonArray ;
    BranchLikelihood * brLikFunction= new BranchLikelihood(tl->getLikelihoodData()->getWeights());
    BrentOneDimension * brentOptimizer = new BrentOneDimension();

    while (oldLk - Lk > tolerance) {
        for (unsigned int i = 0 ; i < nodes.size(); i++) {
            if (nodes[i]->hasFather() ) {
                //Retrieving arrays of interest:
                 fatherData = &(tl->getLikelihoodData()->getNodeData(nodes[i]->getFather()->getId()));
                sonData = &(tl->getLikelihoodData()->getNodeData(nodes[i]->getId()));
                fatherArray   = fatherData->getLikelihoodArrayForNeighbor(nodes[i]->getId());
				vector<double> rootFreqs = tl->getRootFrequencies(0);
                if (nodes[i]->getFather() == tree.getRootNode() ) {
                    //This is the root node, we have to account for the ancestral frequencies:
                    for(unsigned int k = 0; k < fatherArray.size(); k++) //Nb of sites
                    {
                        for(unsigned int j = 0; j < tl->getNumberOfClasses(); j++)
                        {
                            for(unsigned int x = 0; x < tl->getNumberOfStates(); x++)
                                fatherArray[k][j][x] *= rootFreqs[x];
                        }
                    }
                }
                sonArray   = sonData->getLikelihoodArrayForNeighbor(nodes[i]->getFather()->getId());
                //Now we have the two arrays, we can build a branch likelihood object.
                
                //Initialize BranchLikelihood:
                brLikFunction->initModel(tl->getSubstitutionModel(), tl->getRateDistribution());
                brLikFunction->initLikelihoods(&fatherArray, &sonArray);
                ParameterList parameters;
                tl->getParameters().printParameters(std::cout);
                std::string name = "BrLen" + TextTools::toString(nodes[i]->getId());
                Parameter brLen = tl->getParameter(name);
                brLen.setName(name);
                parameters.addParameter(brLen);
                brLikFunction->setParameters(parameters);
                
                //Re-estimate branch length:
                brentOptimizer->setFunction(brLikFunction);
                brentOptimizer->getStopCondition()->setTolerance(0.1);
                brentOptimizer->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
                brentOptimizer->init(parameters);
                brentOptimizer->optimize();
                //brLenNNIValues_[nodeId] = brLikFunction_->getParameterValue("BrLen");
                double length = brentOptimizer->getParameters().getParameter(name).getValue();
                //brLenParameters_.setParameterValue(name, length);
                tl->setParameterValue(name, length);
                nodes[i]->setDistanceToFather(length);

                brLikFunction->resetLikelihoods(); //Array1 and Array2 will be destroyed after this function call.
                //We should not keep pointers towards them...

                
                
            }
        }
    }
    
    
    
}


/******************************************************************************/
// This function builds a bionj tree
/******************************************************************************/

TreeTemplate<Node>  * buildBioNJTree (std::map<std::string, std::string> & params, 
                                      SiteContainer * sites, 
                                      SubstitutionModel* model, 
                                      DiscreteDistribution* rDist, 
                                      Alphabet *alphabet) {
    TreeTemplate<Node>  *unrootedGeneTree = 0;
    
    //We don't want to use rate heterogeneity across sites with bionj, it seems to behave weirdly (e.g. very long branches)
    auto_ptr<DiscreteDistribution>   rDist2 (new ConstantDistribution(1.) );
    DistanceEstimation distEstimation(model->clone(), rDist2.release(), sites, 1, false);
    
    //DistanceEstimation distEstimation(model, rDist, sites, 1, false);
    	
    BioNJ * bionj = new BioNJ();     

    bionj->outputPositiveLengths(true);  

    std::string type = ApplicationTools::getStringParameter("bionj.optimization.method", params, "init");
    if(type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
    else if(type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
    else if(type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
    else throw Exception("Unknown parameter estimation procedure '" + type + "'.");
    // Should I ignore some parameters?
    ParameterList allParameters = model->getParameters();
    // allParameters.addParameters(rDist->getParameters());
    allParameters.addParameters(distEstimation.getRateDistribution().getParameters());

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

    unrootedGeneTree = OptimizationTools::buildDistanceTree(distEstimation, *bionj, parametersToIgnore, !ignoreBrLen, type, tolerance);
    std::vector<Node*> nodes = unrootedGeneTree->getNodes();
    
    for(unsigned int k = 0; k < nodes.size(); k++)
    {
        if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() < 0.000001) nodes[k]->setDistanceToFather(0.000001);
        if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() >= 5) throw Exception("Found a very large branch length in a gene tree; will avoid this family.");
        //Remove long branches, in case it helps
        //if(nodes[k]->hasDistanceToFather() && nodes[k]->getDistanceToFather() > 0.5) nodes[k]->setDistanceToFather(0.05);        
    }

    //delete rDist2;
    delete bionj;  
    return unrootedGeneTree;
}


/******************************************************************************/
// This function refines a gene tree topology and branch lengths using the PhyML 
// algorithm.
/******************************************************************************/

void refineGeneTreeUsingSequenceLikelihoodOnly (std::map<std::string, std::string> & params, 
                                                TreeTemplate<Node>  *& unrootedGeneTree, 
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
        if (f) {
            MPI::COMM_WORLD.Abort(1);
        exit(-1);
        }
    }
    ApplicationTools::displayError("!!! Looking at each site:");
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
        MPI::COMM_WORLD.Abort(1);
    exit(-1);
    }
  
  tl = dynamic_cast<NNIHomogeneousTreeLikelihood*>(PhylogeneticsApplicationTools::optimizeParameters(tl, 
                                                                                                     tl->getParameters(), 
                                                                                                     params));
  std::cout << "Sequence likelihood: "<< -logL << "; Optimized sequence log likelihood "<< -tl->getValue() <<" for family: " << file << std::endl; 
  params[ std::string("optimization.topology.algorithm_nni.method")] = backupParamOptTopo ;
  params[ std::string("optimization.tolerance")] = backupParamOptTol;
  params[ std::string("optimization.max_number_f_eval")] = backupParamnumEval;
  //delete unrootedGeneTree;
  unrootedGeneTree = new TreeTemplate<Node> ( tl->getTree() );
  delete tl;
}




/*
string parenthesisWithSpeciesNamesToGeneTree (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp ) {  
  //,                                             std::vector<string> &geneNames) {
  
}
*/

/**************************************************************************
 * This function produces a string version of a gene tree, 
 * with gene names replaced by species names. 
 **************************************************************************/
string geneTreeToParenthesisWithSpeciesNames (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp ) {
  TreeTemplate<Node> * geneTreeWithSpNames = geneTree->clone();
  vector <Node*> leaves = geneTreeWithSpNames->getLeaves();
  std::map<std::string, std::string >::const_iterator seqtosp;
  for (unsigned int i=0 ; i<leaves.size() ; i++ ) {
    seqtosp=seqSp.find(leaves[i]->getName());
    if (seqtosp!=seqSp.end()){
      leaves[i]->setName(seqtosp->second);
    }
    else {
      std::cout <<"Error in assignSpeciesIdToLeaf: "<< leaves[i]->getName() <<" not found in std::map seqSp"<<std::endl;
        MPI::COMM_WORLD.Abort(1);
      exit(-1);
    }
  }
  return (TreeTemplateTools::treeToParenthesis(*geneTreeWithSpNames, false) );  
}

/**************************************************************************
 * This function produces a string version of a gene tree, 
 * with gene names changed to include species names. 
 **************************************************************************/
string geneTreeToParenthesisPlusSpeciesNames (TreeTemplate<Node> * geneTree,
                                              std::map<std::string, std::string > seqSp ) {
  TreeTemplate<Node> * geneTreeWithSpNames = geneTree->clone();
  vector <Node*> leaves = geneTreeWithSpNames->getLeaves();
  std::map<std::string, std::string >::const_iterator seqtosp;
  for (unsigned int i=0 ; i<leaves.size() ; i++ ) {
    seqtosp=seqSp.find(leaves[i]->getName());
    if (seqtosp!=seqSp.end()){
      leaves[i]->setName(seqtosp->second + "%" + leaves[i]->getName());
    }
    else {
      std::cout <<"Error in assignSpeciesIdToLeaf: "<< leaves[i]->getName() <<" not found in std::map seqSp"<<std::endl;
      MPI::COMM_WORLD.Abort(1);
      exit(-1);
    }
  }
  return (TreeTemplateTools::treeToParenthesis(*geneTreeWithSpNames, false) );  
}

/**************************************************************************
 * This function produces a gene tree from a string version in which 
 * gene names have been changed to include species names. 
 **************************************************************************/
TreeTemplate<Node> * parenthesisPlusSpeciesNamesToGeneTree (string geneTreeStr) {
  TreeTemplate<Node> * geneTree = TreeTemplateTools::parenthesisToTree(geneTreeStr);
  vector <Node*> leaves = geneTree->getLeaves();
  for (vector<Node *>::iterator it = leaves.begin(); it!=leaves.end(); it++ )
    {
    StringTokenizer stk ((*it)->getName(), "%");
//    Tokenize((*it)->getName(),tokens,"%");
    (*it)->setName(stk.getToken(1));
    (*it)->setNodeProperty("S", BppString(stk.getToken(0)));
  }
  return geneTree;  
}


/**************************************************************************
 * This function produces a gene tree with leaves annotated with species names.
 **************************************************************************/
void annotateGeneTreeWithSpeciesNames (TreeTemplate<Node> * geneTree,
                                       std::map<std::string, std::string > seqSp ) {
  vector <Node*> leaves = geneTree->getLeaves();
  std::map<std::string, std::string >::const_iterator seqtosp;
  for (unsigned int i=0 ; i<leaves.size() ; i++ ) {
    seqtosp=seqSp.find(leaves[i]->getName());
    if (seqtosp!=seqSp.end()){
      leaves[i]->setNodeProperty("S", BppString(seqtosp->second));
    }
    else {
      std::cout <<"Error in assignSpeciesIdToLeaf: "<< leaves[i]->getName() <<" not found in std::map seqSp"<<std::endl;
      MPI::COMM_WORLD.Abort(1);
      exit(-1);
    }
  }
  return;  
}




/**************************************************************************
 * This function optimizes a gene tree based on the reconciliation score only.
 * It uses SPRs and NNIs, and calls findMLReconciliationDR to compute the likelihood.
 **************************************************************************/

double refineGeneTreeDLOnly (TreeTemplate<Node> * spTree, 
                             TreeTemplate<Node> *& geneTree, 
                             std::map<std::string, std::string > seqSp,
                             std::map<std::string, int > spID,
                             std::vector< double> &lossExpectedNumbers, 
                             std::vector < double> &duplicationExpectedNumbers, 
                             int & MLindex, 
                             std::vector <int> &num0lineages, 
                             std::vector <int> &num1lineages, 
                             std::vector <int> &num2lineages, 
                             std::set <int> &nodesToTryInNNISearch)
{
    TreeTemplate<Node> *tree = 0;
    TreeTemplate<Node> *bestTree = 0;
    TreeTemplate<Node> *currentTree = 0;
    currentTree = geneTree->clone();
    breadthFirstreNumber (*currentTree);//, duplicationExpectedNumbers, lossExpectedNumbers);
    int sprLimit = 10; //Arbitrary.
    std::vector <int> nodeIdsToRegraft; 
    double bestlogL;
    double logL;
    bool betterTree;
    unsigned int numIterationsWithoutImprovement = 0;
    double startingML = findMLReconciliationDR (spTree, 
                                                currentTree, 
                                                seqSp,
                                                spID,
                                                lossExpectedNumbers, 
                                                duplicationExpectedNumbers, 
                                                MLindex, 
                                                num0lineages, 
                                                num1lineages, 
                                                num2lineages, 
                                                nodesToTryInNNISearch, false);
    tree = currentTree->clone();
    bestTree = currentTree->clone();
    bestlogL = startingML;
    
    while (numIterationsWithoutImprovement + 1 < geneTree->getNumberOfNodes()) {
        for (unsigned int nodeForSPR=geneTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) {
            betterTree = false;
            tree = currentTree->clone();
            buildVectorOfRegraftingNodesLimitedDistance(*tree, nodeForSPR, sprLimit, nodeIdsToRegraft);
            
            for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) {
                if (tree) {
                    delete tree;
                }
                tree = currentTree->clone();
                makeSPR(*tree, nodeForSPR, nodeIdsToRegraft[i], false);
                logL = findMLReconciliationDR (spTree, 
                                               tree, 
                                               seqSp,
                                               spID,
                                               lossExpectedNumbers, 
                                               duplicationExpectedNumbers, 
                                               MLindex, 
                                               num0lineages, 
                                               num1lineages, 
                                               num2lineages, 
                                               nodesToTryInNNISearch, false);
                if (logL-0.01>bestlogL) {
                    betterTree = true;
                    bestlogL =logL;
                    if (bestTree)
                        delete bestTree;
                    bestTree = tree->clone();  
                    /*std::cout << "Gene tree SPR: Better candidate tree likelihood : "<<bestlogL<< std::endl;
                     std::cout << TreeTemplateTools::treeToParenthesis(*tree, true)<< std::endl;*/
                }
            }
            if (betterTree) {
                logL = bestlogL; 
                if (currentTree)
                    delete currentTree;
                currentTree = bestTree->clone();
                breadthFirstreNumber (*currentTree);//, duplicationExpectedNumbers, lossExpectedNumbers);
                //std::cout <<"NEW BETTER TREE: \n"<< TreeTemplateTools::treeToParenthesis(*currentTree, true)<< std::endl;
                numIterationsWithoutImprovement = 0;
            }
            else {
                logL = bestlogL; 
                if (currentTree)
                    delete currentTree;
                currentTree = bestTree->clone(); 
                numIterationsWithoutImprovement++;
            }
        }
    }
    if (geneTree)
        delete geneTree;
    geneTree = bestTree->clone();
    if (tree) delete tree;
    if (bestTree) delete bestTree;
    if (currentTree) delete currentTree;
   // std::cout << "DL initial likelihood: "<< startingML << "; Optimized DL log likelihood "<< bestlogL <<"." << std::endl; 
    return bestlogL;
}




/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/


void buildVectorOfRegraftingNodesGeneTree(TreeTemplate<Node> &spTree, 
                                          TreeTemplate<Node> &tree, 
                                          int nodeForSPR, 
                                          int distance, 
                                          std::vector <int> & nodeIdsToRegraft) {
    
    Node * N = tree.getNode(nodeForSPR);
    Node * brother = 0;
    if (N->hasFather()) {
        if (N->getFather()->getSon(0) != N) {
            brother = N->getFather()->getSon(0);
        }
        else {
            brother = N->getFather()->getSon(1);
        }
    }
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(N));
    
    // std::vector <int> allNodeIds = tree.getNodesId();
    std::vector <int> allNodeIds;
    std::vector <int> allProximalNodes;
  //TEST12 11 2011
    getNeighboringNodesIdLimitedDistance(tree, nodeForSPR, distance, allProximalNodes);
    //Here, all nodes below some threshold have been added.
    
    //Another idea: try all regrafting points such that the new father node has a species id between the pruned node species id and the species id of the father of the pruned node.
    
    //Now we add nodes that may be good regrafting points based on their species ID.
    
//    if (tree.getNode(nodeForSPR)->hasFather() && tree.getNode(nodeForSPR)->getFather()->hasNodeProperty("S")) {
    if (N->hasNodeProperty("S")) {

        std::vector <std::string> spIds;
        int spId = TextTools::toInt ( (dynamic_cast<const BppString *>(N->getNodeProperty("S")))->toSTL());
        int fatherSpId ;
        Node * node = spTree.getNode(spId);
        if (N->hasFather() && (N->getFather()->hasNodeProperty("S"))) {

            fatherSpId = TextTools::toInt ( (dynamic_cast<const BppString *>(N->getFather()->getNodeProperty("S")))->toSTL());
            //we can try all nodes in the tree downstream from the pruned subtree, with species Ids that are between spId and fatherSpId
            Node * n = node->getFather();
            vector <int> PotentialSpeciesIds ;
            PotentialSpeciesIds.push_back(spId);
            PotentialSpeciesIds.push_back(fatherSpId);


            while (n->hasFather() && n->getId() > fatherSpId) {
                PotentialSpeciesIds.push_back(n->getId());
                n = n->getFather();
            }

            for (unsigned int i = 0 ; i < PotentialSpeciesIds.size() ; i++) {
                getSonsOfNodesWithSimilarSpeciesIds(tree.getRootNode(), TextTools::toString(PotentialSpeciesIds[i]), allNodeIds);
                //TEST getSonsOfNodesWithSimilarSpeciesIds(brother, TextTools::toString(PotentialSpeciesIds[i]), allNodeIds);
            }
           if (N->hasFather() && N->getFather()->hasNodeProperty("S")) {
                getNodesWithSimilarSpeciesIds(N->getFather(), (dynamic_cast<const BppString *>(N->getFather()->getNodeProperty("S")))->toSTL(), allNodeIds);
            }
            getNodesWithSimilarSpeciesIdsUpstream(N, (dynamic_cast<const BppString *>(N->getFather()->getNodeProperty("S")))->toSTL(), allNodeIds);
        }
    /*    spIds.push_back(TextTools::toString(spId));
        //Now we look for the brother species, which we put in spIds
        if (node->hasFather()) {
            if (node->getFather()->getSon(0) != node) {
                spIds.push_back( TextTools::toString(node->getFather()->getSon(0)->getId()));
            }
            else {
                spIds.push_back( TextTools::toString(node->getFather()->getSon(1)->getId()));
            }
        }
        getAllCandidateBranchingPointsFromSpeciesID (tree, spIds, allNodeIds);*/
    }
            //We also add nodes with the same species id, downstream from the father node in the gene tree
    
    
            
    //std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)->getFather()));
    
    forbiddenIds.push_back(tree.getRootNode()->getId());
    
    //We don't want to reinsert the pruned subtree in the same place!
    forbiddenIds.push_back(brother->getId());
    
    //We remove the nodes that are not appropriate for regrafting
    std::vector <int> toRemove;
    for (unsigned int i = 0 ; i< allNodeIds.size() ; i++) {
        if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
            toRemove.push_back(i);
        }
    }
    
    /*  std:: cout <<"nodeForSPR: "<< nodeForSPR <<"; FORBIDDEN IDS: "<<std::endl;
     VectorTools::print(forbiddenIds);*/
    sort(toRemove.begin(), toRemove.end(), cmp);
    /*VectorTools::print(forbiddenIds);
     sort(allNodeIds.begin(), allNodeIds.end(), anticmp);*/
    for (unsigned int i = 0 ; i< toRemove.size() ; i++) {
        std::vector<int>::iterator vi = allNodeIds.begin();
        allNodeIds.erase(vi+toRemove[i]);
    }
    
    allNodeIds = VectorTools::unique (allNodeIds);
    //Now we only consider those nodes not far from the pruned node:
    allNodeIds = VectorTools::vectorIntersection(allNodeIds, allProximalNodes);
    
    //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
    nodeIdsToRegraft = allNodeIds;
}


/**************************************************************************
 * This function returns a vector of branching points that may diminish the number of duplications/losses.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/
void getAllCandidateBranchingPointsFromSpeciesID (TreeTemplate<Node> &tree, 
                                                  std::vector <std::string> spIds, 
                                                  std::vector <int> & allNodeIds) {
    for (unsigned int k = 0 ; k < spIds.size(); k++) {
        vector <Node *> nodes = tree.getNodes();
        for (unsigned int i = 0 ; i < nodes.size() ; i++ ) {
            if (nodes[i]->hasNodeProperty("S") && (dynamic_cast<const BppString *>(nodes[i]->getNodeProperty("S")))->toSTL() == spIds[k]) {
                allNodeIds.push_back(nodes[i]->getId()); //necessary?
              /*  for (unsigned int j = 0 ; j < nodes[i]->getNumberOfSons() ; j++) {
                    allNodeIds.push_back(nodes[i]->getSon(j)->getId());
                }*/
            }
        }
    }
    allNodeIds = VectorTools::unique (allNodeIds);
    return;
}

/**************************************************************************
 * This recursive function returns node ids with a given species id, in a subtree starting in node.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIds(Node * node, string spId, std::vector <int> & allNodeIds) {
    for (unsigned int i = 0 ; i < node->getNumberOfSons() ; i++ ) {
        if (node->getSon(i)->hasNodeProperty("S") && (dynamic_cast<const BppString *>(node->getSon(i)->getNodeProperty("S")))->toSTL() == spId)
        {
            allNodeIds.push_back(node->getSon(i)->getId());
            getNodesWithSimilarSpeciesIds(node->getSon(i), spId, allNodeIds);
        }
    }
    return;
}

/**************************************************************************
 * This recursive function returns node ids with a given species id, upstream from Node node.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getNodesWithSimilarSpeciesIdsUpstream(Node * node, string spId, std::vector <int> & allNodeIds) {
        if (node->hasFather() && node->getFather()->hasNodeProperty("S") && (dynamic_cast<const BppString *>(node->getFather()->getNodeProperty("S")))->toSTL() == spId)
        {
            if (node->getFather()->getSon(0) == node) {
                allNodeIds.push_back(node->getFather()->getSon(1)->getId());
            }
            else {
                allNodeIds.push_back(node->getFather()->getSon(0)->getId());
            }
            getNodesWithSimilarSpeciesIdsUpstream(node->getFather(), spId, allNodeIds);
        }
    return;
}


/**************************************************************************
 * This recursive function returns node ids of sons of nodes with a given species id, in a subtree starting in node.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/

void getSonsOfNodesWithSimilarSpeciesIds(Node * node, string spId, std::vector <int> & allNodeIds) {
    for (unsigned int i = 0 ; i < node->getNumberOfSons() ; i++ ) {
        if (node->hasNodeProperty("S") && (dynamic_cast<const BppString *>(node->getNodeProperty("S")))->toSTL() == spId)
        {
            allNodeIds.push_back(node->getSon(i)->getId());
        }
        getSonsOfNodesWithSimilarSpeciesIds(node->getSon(i), spId, allNodeIds);
    }
    return;
}

/**************************************************************************
 * This function returns a vector of branching points for gene tree SPR in the coalescent framework.
 * The gene tree has to be rooted and annotated with species numbers.
 **************************************************************************/


void buildVectorOfRegraftingNodesCoalGeneTree(TreeTemplate<Node> &spTree, 
											  TreeTemplate<Node> &tree, 
											  int nodeForSPR, 
											  int distance, 
											  std::vector <int> & nodeIdsToRegraft) {
    
    Node * N = tree.getNode(nodeForSPR);
    Node * brother = 0;
    if (N->hasFather()) {
        if (N->getFather()->getSon(0) != N) {
            brother = N->getFather()->getSon(0);
        }
        else {
            brother = N->getFather()->getSon(1);
        }
    }
    std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(N));
    
    // std::vector <int> allNodeIds = tree.getNodesId();
    std::vector <int> allNodeIds;
    std::vector <int> allProximalNodes;
	//TEST12 11 2011
    getNeighboringNodesIdLimitedDistance(tree, nodeForSPR, distance, allProximalNodes);
    //Here, all nodes below some threshold have been added.
    
    //Another idea: try all regrafting points such that the new father node has a species id between the pruned node species id and the species id of the father of the pruned node.
    
    //Now we add nodes that may be good regrafting points based on their species ID.
    
    
	
    //std::vector <int> forbiddenIds = TreeTemplateTools::getNodesId(*(tree.getNode(nodeForSPR)->getFather()));
    
    forbiddenIds.push_back(tree.getRootNode()->getId());
    
    //We don't want to reinsert the pruned subtree in the same place!
    forbiddenIds.push_back(brother->getId());
    
    //We remove the nodes that are not appropriate for regrafting
    std::vector <int> toRemove;
    for (unsigned int i = 0 ; i< allNodeIds.size() ; i++) {
        if (VectorTools::contains(forbiddenIds, allNodeIds[i])) {
            toRemove.push_back(i);
        }
    }
    
    /*  std:: cout <<"nodeForSPR: "<< nodeForSPR <<"; FORBIDDEN IDS: "<<std::endl;
     VectorTools::print(forbiddenIds);*/
    sort(toRemove.begin(), toRemove.end(), cmp);
    /*VectorTools::print(forbiddenIds);
     sort(allNodeIds.begin(), allNodeIds.end(), anticmp);*/
    for (unsigned int i = 0 ; i< toRemove.size() ; i++) {
        std::vector<int>::iterator vi = allNodeIds.begin();
        allNodeIds.erase(vi+toRemove[i]);
    }
    
    allNodeIds = VectorTools::unique (allNodeIds);
    //Now we only consider those nodes not far from the pruned node:
    allNodeIds = VectorTools::vectorIntersection(allNodeIds, allProximalNodes);
    
    //Now allNodeIds contains all the Ids of nodes where the subtree can be reattached.
    nodeIdsToRegraft = allNodeIds;
}



/******************************************************************************/

unsigned int optimizeBranchLengthsParameters(
											 DiscreteRatesAcrossSitesTreeLikelihood* tl,
											 const ParameterList& parameters,
											 double target,
											 OptimizationListener* listener,
											 double tolerance,
											 unsigned int tlEvalMax,
											 OutputStream* messageHandler,
											 OutputStream* profiler,
											 unsigned int verbose,
											 const std::string& optMethodDeriv)
throw (Exception)
{
	// Build optimizer:
	Optimizer* optimizer = 0;
	if (optMethodDeriv == OptimizationTools::OPTIMIZATION_GRADIENT)
	{
		tl->enableFirstOrderDerivatives(true);
		tl->enableSecondOrderDerivatives(false);
		optimizer = new ConjugateGradientMultiDimensions(tl);
	}
	else if (optMethodDeriv == OptimizationTools::OPTIMIZATION_NEWTON)
	{
		tl->enableFirstOrderDerivatives(true);
		tl->enableSecondOrderDerivatives(true);
		optimizer = new PseudoNewtonOptimizer(tl);
	}
	else if (optMethodDeriv == OptimizationTools::OPTIMIZATION_BFGS)
	{
		tl->enableFirstOrderDerivatives(true);
		tl->enableSecondOrderDerivatives(false);
		optimizer = new BfgsMultiDimensions(tl);
	}
	else
		throw Exception("OptimizationTools::optimizeBranchLengthsParameters. Unknown optimization method: " + optMethodDeriv);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->getStopCondition()->setTolerance(tolerance);
	// Optimize TreeLikelihood function:
	ParameterList pl = parameters.getCommonParametersWith(tl->getBranchLengthsParameters());
	optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
	if (listener)
		optimizer->addOptimizationListener(listener);
	optimizer->init(pl);

	//We do the optimization one iteration at a time
	optimizer->setMaximumNumberOfEvaluations( 2 ) ; //tlEvalMax);
	std::vector <double> values;
	values.push_back(tl->getValue());
	bool stop = false;
	while (!stop) {
		optimizer->optimize();
		values.push_back(optimizer->getFunction()->getValue());
		/*std::cout << "Target: "<< target << "; currentlk: "<< optimizer->getFunction()->getValue()<<std::endl;
		std::cout << "values.size(): "<<values.size()<<"; values[values.size()-1]: "<<values[values.size()-1]<< "; values[values.size()-2]: "<<values[values.size()-2] << std::endl;*/
		if ( ( values[values.size()-1] + (30-values.size()) * (values[values.size()-1] - values[values.size()-2]) > target /*+ 10*/) || (values[values.size()-1] < target) ){
			stop = true;
		}
		if (verbose > 0){
			std::cout << "Iteration "<< values.size()-1 <<"; New value: "<< optimizer->getFunction()->getValue() <<std::endl;
			ApplicationTools::displayMessage("\n");
		}
	}
	if ( (values[values.size()-1] < target) && ((values[values.size()-2] - values[values.size()-1]) > 0.1) ) {
		//optimizer->setMaximumNumberOfEvaluations( 2 ) ; //tlEvalMax);
		optimizer->optimize();
	}
	
	/*std::cout << "Values ";
	VectorTools::print(values);*/
	// We're done.
	unsigned int n = values.size()-1 ; //optimizer->getNumberOfEvaluations();
	delete optimizer;
	return n;
}

/******************************************************************************/

