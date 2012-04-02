//
//  CoalTools.cpp
//  phyldog
//
//  Created by Bastien Boussau on 07/03/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include <iostream>
#include "CoalTools.h"


//mpic++ -g -lbpp-core -lbpp-seq -lbpp-phyl -lboost_serialization -lboost_mpi CoalTools.cpp ReconciliationTools.cpp -o CoalTools

//Test main function
int main(int args, char ** argv)
{
   std::string spStr = "(((A,B),C),(D,E));";
   std::string gStr = "((A,(B,C)),(D,E));";
//    std::string gStr = "(((A,D), (B,C)),E);";
   TreeTemplate <Node>* spTree = TreeTemplateTools::parenthesisToTree(spStr);
      TreeTemplate <Node>* geneTree = TreeTemplateTools::parenthesisToTree(gStr);
   
    std::map<std::string, std::string > seqSp;
    
  /*  seqSp.insert( pair<std::string, std::string >("taxon1","taxon1") );
    seqSp.insert( pair<std::string, std::string >("taxon2","taxon2") );
    seqSp.insert( pair<std::string, std::string >("taxon3","taxon3") );
    seqSp.insert( pair<std::string, std::string >("taxon4","taxon4") );
    seqSp.insert( pair<std::string, std::string >("taxon5","taxon5") );
    seqSp.insert( pair<std::string, std::string >("taxon6","taxon6") );
    seqSp.insert( pair<std::string, std::string >("taxon7","taxon7") );
    seqSp.insert( pair<std::string, std::string >("taxon8","taxon8") );*/

    
    
   seqSp.insert( pair<std::string, std::string >("A","A") );
    seqSp.insert( pair<std::string, std::string >("B","B") );
    seqSp.insert( pair<std::string, std::string >("C","C") );
    seqSp.insert( pair<std::string, std::string >("D","D") );
    seqSp.insert( pair<std::string, std::string >("E","E") );/*
    seqSp.insert( pair<std::string, std::string >("F","F") );
    seqSp.insert( pair<std::string, std::string >("G","G") );
    seqSp.insert( pair<std::string, std::string >("H","H") );
*/
    
    
    
    
    //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
    std::vector < std::vector< std::vector< std::vector< unsigned int > > > > coalCounts;
    std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
    std::vector< std::vector<unsigned int> > coalCounts3;
    std::vector< unsigned int > coalCounts4;
    //speciesIDs: vector of genetreenbnodes vectors of 3 (3 directions) ints
    std::vector <std::vector<unsigned int> > speciesIDs;
    std::vector < unsigned int > speciesIDs2;
    for (unsigned int i = 0 ; i < 2 ; i++ ) {
        coalCounts4.push_back(0);
    }
    for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
        coalCounts3.push_back(coalCounts4);
    }
    for (unsigned int i = 0 ; i < 3 ; i++ ) {
        coalCounts2.push_back(coalCounts3);
        speciesIDs2.push_back(0);
    }
    for (unsigned int i = 0 ; i < geneTree->getNumberOfNodes() ; i++ ) {
        coalCounts.push_back(coalCounts2);
        speciesIDs.push_back(speciesIDs2);
    }
    
    breadthFirstreNumber (*spTree);
    std::map<std::string, int > spID = computeSpeciesNamesToIdsMap(*spTree);
    
    for(std::map<std::string, int >::iterator it = spID.begin(); it != spID.end(); it++){
        std::cout <<"it->first: "<< it->first << " : "<< it->second <<std::endl;
    }

    computeSubtreeCoalCountsPostorder(*spTree, 
                                      *geneTree, 
                                      geneTree->getRootNode(), 
                                      seqSp, 
                                      spID, 
                                      coalCounts,
                                      speciesIDs);
    //Add the starting lineage at the root
    for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
        if (coalCounts[geneTree->getRootNode()->getId()][0][i][0]==0 && coalCounts[geneTree->getRootNode()->getId()][0][i][1] !=0)
        {
            //coalCounts[geneTree->getRootNode()->getId()][0][0][0] = 1;
            coalCounts[geneTree->getRootNode()->getId()][0][i][0]=1; 
            break;
            
        }
    }

    std::cout <<"Species Tree: \n" <<    TreeTools::treeToParenthesis(*spTree, true) << std::endl;
    std::cout <<"Gene Tree: \n" <<    TreeTools::treeToParenthesis(*geneTree, true) << std::endl;

      for (unsigned int i = 0 ; i < coalCounts[geneTree->getRootNode()->getId()][0].size() ; i++) {
          std::cout << "Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] << " Num coal out: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][1]<<std::endl;
     }
    
    std::vector<double> bls (spTree->getNumberOfNodes(), 2.0);

    double initialLikelihood = computeCoalLikelihood ( coalCounts[geneTree->getRootNode()->getId()][0], bls ) ;
    
    std::cout << "Initial Likelihood: "<< initialLikelihood <<std::endl;
    
    std::map <double, Node*> LksToNodes;
    //Now doing the preorder tree traversal
    Node * geneRoot = geneTree->getRootNode();
    std::vector <Node *> sons = geneRoot->getSons();
    if (sons.size()!=2) {
        std::cerr <<"Error: "<<sons.size()<< "sons at the root!"<<std::endl; 
    }

    LksToNodes[initialLikelihood]=sons[0];
    //We fill the likelihood and species ID data for the root node.
    //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
    //with other nodes.
    coalCounts[geneRoot->getId()][1] = coalCounts[geneRoot->getSon(1)->getId()][0];
    coalCounts[geneRoot->getId()][2] = coalCounts[geneRoot->getSon(0)->getId()][0];
    speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon(1)->getId()][0];
    speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon(0)->getId()][0];

    
    
    
    for (unsigned int i = 0; i< sons.size(); i++){
        for (unsigned int j =0; j<sons[i]->getNumberOfSons(); j++) {
            computeSubtreeCoalCountsPreorder(*spTree, 
                                             *geneTree, 
                                             sons[i], 
                                             seqSp, 
                                             spID, 
                                             coalCounts,
                                             bls, 
                                             speciesIDs, j, LksToNodes);
        }
    }
    for (unsigned int j = 0 ; j < geneTree->getNumberOfNodes() ; j++) {
        std::cout << "Node j: " <<std::endl;
        for (unsigned int i = 0 ; i < coalCounts[j][0].size() ; i++) {
            std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][0][i][0] << " Num coal out: "<< coalCounts[j][0][i][1]<<std::endl;
            std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][1][i][0] << " Num coal out: "<< coalCounts[j][1][i][1]<<std::endl;
            std::cout << "\tSp Branch "<<i<<" Num coal in: "<< coalCounts[j][2][i][0] << " Num coal out: "<< coalCounts[j][2][i][1]<<std::endl;
                                                                                                            
        }
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    Newick *treeReader = new Newick(true);
    vector<Tree *> gTrees;
    //treeReader->read("testGeneTrees", gTrees);
//   treeReader->read("geneTrees7", gTrees);
//treeReader->read("5000SimulatedTrees8Taxa.trees", gTrees);
   // treeReader->read("mesquite500SimCoalTreesNe10000_Depth100000.proj.phy", gTrees);
  //  treeReader->read("mesquite500SimCoalTreesNe10000_Depth100000.proj.phy", gTrees);
    treeReader->read("40Sp_500GeneTrees.phy", gTrees);

    
    vector<Tree *> spTrees;
//    treeReader->read("spTree", spTrees);
//    treeReader->read("spTree2.0", spTrees);
//    treeReader->read("spTree7", spTrees);
  //  treeReader->read("SpTree8.tree", spTrees);
//    treeReader->read("1SimulatedSpTree8Taxa_Depth100000.tree", spTrees);
  treeReader->read("40Sp_First.phy", spTrees);
    
    delete treeReader;
//TEMP    TreeTemplate <Node>* spTree = dynamic_cast<TreeTemplate<Node>*>(spTrees[0]);
    vector<string> leaves = spTree->getLeavesNames();
    for (unsigned int i = 0 ; i < leaves.size() ; i++) {
        seqSp.insert( pair<std::string, std::string >(leaves[i],leaves[i]) );
    }
    
    breadthFirstreNumber (*spTree);
 //TEMP   std::map<std::string, int > spID = computeSpeciesNamesToIdsMap(*spTree);
    
    for(std::map<std::string, int >::iterator it = spID.begin(); it != spID.end(); it++){
        std::cout <<"it->first: "<< it->first << " : "<< it->second <<std::endl;
    }

    
    
    /* Checking the likelihood computation against results from Rosenberg 2002: it works.
    //11
    vector<unsigned int > vec(2, 1);
    double loglk = computeCoalLikelihood ( vec, 1.0 ) ;
    std::cout << "11: "<< loglk<< std::endl;
    
    vec[1] = 2;
     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
    std::cout << "12: "<<loglk<< std::endl;

    vec[0] = 2;
     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
    std::cout << "22: "<<loglk<< std::endl;

    vec[1] = 3;
     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
    std::cout << "23: "<<loglk<< std::endl;

    vec[0] = 3;
     loglk = computeCoalLikelihood ( vec, 1.0 ) ;
    std::cout << "33: "<<loglk<< std::endl;
*/
    
    
    
    
  //TEMP  std::vector<double> bls (spTree->getNumberOfNodes(), 2.0);
/*    for (unsigned int i = 0 ; i < bls.size() ; i++) {
        if ( spTree->getNode(i)->hasFather() )
            bls[i] = spTree->getNode(i)->getDistanceToFather();
    }
  */  
    std::cout <<"Species Tree: \n" <<    TreeTools::treeToParenthesis(*spTree, true) << std::endl;
//TEMP   TreeTemplate <Node>* geneTree = 0;
    
    //allGeneCounts: nbSpeciesBranches vectors of gTrees.size() vectors of 2 ints
    std::vector < std::vector < std::vector<unsigned int> > >  allGeneCounts;
    std::vector< std::vector< unsigned int > > allGeneCounts2;
    std::vector< unsigned int > allGeneCounts3;
    for (unsigned int i = 0 ; i < 2 ; i++ ) {
        allGeneCounts3.push_back(0);
    }
    for (unsigned int i = 0 ; i < gTrees.size() ; i++ ) {
        allGeneCounts2.push_back(allGeneCounts3);
    }
    for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
        allGeneCounts.push_back(allGeneCounts2);
    }
    
    
    
    for (unsigned int j = 0 ; j < gTrees.size() ; j++) {
        geneTree= dynamic_cast<TreeTemplate<Node>*>(gTrees[j]);
        //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
        std::vector < std::vector< std::vector< std::vector< unsigned int > > > > coalCounts;
        std::vector< std::vector< std::vector< unsigned int > > > coalCounts2;
        std::vector< std::vector<unsigned int> > coalCounts3;
        std::vector< unsigned int > coalCounts4;
        //speciesIDs: vector of genetreenbnodes vectors of 3 (3 directions) ints
        std::vector <std::vector<unsigned int> > speciesIDs;
        std::vector < unsigned int > speciesIDs2;
        for (unsigned int i = 0 ; i < 2 ; i++ ) {
            coalCounts4.push_back(0);
        }
        for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
            coalCounts3.push_back(coalCounts4);
        }
        for (unsigned int i = 0 ; i < 3 ; i++ ) {
            coalCounts2.push_back(coalCounts3);
            speciesIDs2.push_back(0);
        }
        for (unsigned int i = 0 ; i < geneTree->getNumberOfNodes() ; i++ ) {
            coalCounts.push_back(coalCounts2);
            speciesIDs.push_back(speciesIDs2);
        }
        
        computeSubtreeCoalCountsPostorder(*spTree, 
                                          *geneTree, 
                                          geneTree->getRootNode(), 
                                          seqSp, 
                                          spID, 
                                          coalCounts,
                                          speciesIDs);
        //Add the starting lineage at the root
        for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
            if (coalCounts[geneTree->getRootNode()->getId()][0][i][0]==0 && coalCounts[geneTree->getRootNode()->getId()][0][i][1] !=0)
            {
        //coalCounts[geneTree->getRootNode()->getId()][0][0][0] = 1;
                coalCounts[geneTree->getRootNode()->getId()][0][i][0]=1; 
                break;

            }
        }

         /*   for (unsigned int i = 0 ; i < coalCounts[geneTree->getRootNode()->getId()][0].size() ; i++) {
         std::cout << "Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][0] << " Num coal out: "<< coalCounts[geneTree->getRootNode()->getId()][0][i][1]<<std::endl;
         }*/
        //Now we compute the likelihood of the rooted gene tree:
        double loglk = computeCoalLikelihood ( coalCounts[geneTree->getRootNode()->getId()][0], bls ) ;
      //  std::cout << "lk: "<<exp(loglk)<<std::endl;
        for (unsigned int i = 0 ; i < spTree->getNumberOfNodes() ; i++ ) {
            allGeneCounts[i][j] = coalCounts[geneTree->getRootNode()->getId()][0][i];
        }
    }
   
    //If we want to do numerical optimization
    // BrentOneDimension *brentOptimizer = new BrentOneDimension();
    vector<double> blAnalytical = bls;
    for (unsigned int i = 0 ; i < allGeneCounts.size() ; i++) {
        CoalBranchLikelihood *brLikFunction = new CoalBranchLikelihood(allGeneCounts[i]);
        //Initialize BranchLikelihood:
        blAnalytical[i] = brLikFunction->initModel();

        //If we want to do numerical optimization
        /*
        ParameterList parameters;
        Parameter brLen = Parameter("BrLen", blAnalytical[i]);//, ExcludingPositiveReal (0) );
        parameters.addParameter(brLen);
        brLikFunction->setParameters(parameters);
        
        //Re-estimate branch length:
        brentOptimizer->setVerbose(0);
        brentOptimizer->setFunction(brLikFunction);
        brentOptimizer->getStopCondition()->setTolerance(0.000001);
        brentOptimizer->setInitialInterval(0, brLen.getValue()+2);
        brentOptimizer->init(parameters);
        brentOptimizer->optimize();
        bls[i] =brentOptimizer->getParameters().getParameter("BrLen").getValue();
        std::cout <<"Value for Sp branch "<<i<<": "<< bls[i] << std::endl;
        */
    }
    for (unsigned int i = 0 ; i < bls.size() ; i++) {
        if ( spTree->getNode(i)->hasFather() )
            spTree->getNode(i)->setDistanceToFather(bls[i]);
    }
    std::cout <<"Species Tree optimized: \n" << TreeTools::treeToParenthesis(*spTree, true) << std::endl;
    std::cout <<"Species Tree optimized: \n" << TreeTools::treeToParenthesis(*spTree, false) << std::endl;

    for (unsigned int i = 0 ; i < bls.size() ; i++) {
        if ( spTree->getNode(i)->hasFather() )
            spTree->getNode(i)->setDistanceToFather(blAnalytical[i]);
    }
    std::cout <<"Species Tree analytical: \n" << TreeTools::treeToParenthesis(*spTree, true) << std::endl;
    std::cout <<"Species Tree analytical: \n" << TreeTools::treeToParenthesis(*spTree, false) << std::endl;

    return 0;

}




/*****************************************************************************
 * This function performs a postorder tree traversal in order to  
 * fill up vectors of counts of coalescence events for rootings. 
 * For each branch of the species tree, we need to record how many lineages got in,
 * and how many lineages got out.
 * Thus, for each branch of the species tree, we have 2 ints: 
 * vec[0]: number of incoming lineages
 * vec[1]: number of outgoing lineages
 * When followed by the preorder tree traversal function, 
 * vectors of counts for all rootings are computed.
 * coalCounts contains all lower counts for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeCoalCountsPostorder(TreeTemplate<Node> & spTree, 
                                 TreeTemplate<Node> & geneTree, 
                                 Node * node, 
                                  std::map<std::string, std::string > & seqSp, 
                                  std::map<std::string, int > & spID, 
                                 std::vector < std::vector< std::vector< std::vector< unsigned int > > > > & coalCounts,
                                 std::vector <std::vector<unsigned int> > & speciesIDs) {
	int id=node->getId();
 	if (node->isLeaf()) { //In principle this should be done only once at the beginning of the algorithm, if node ids are kept
        //Fill all leaf count vectors with 0s
        initializeCountVectors(coalCounts[id]);
        speciesIDs[id][0] = assignSpeciesIdToLeaf(node, seqSp, spID);
//        int i = spID[ seqSp[ node->getName()] ];
        incrementOutCount(coalCounts[id][0], speciesIDs[id][0]);
        return;
    }
    else {
        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            computeSubtreeCoalCountsPostorder(spTree, geneTree, 
                                        sons[i], seqSp, 
                                        spID, coalCounts, 
                                        speciesIDs);
        }
        
        int idSon0 = sons[0]->getId();
        int idSon1 = sons[1]->getId();
        unsigned int directionSon0, directionSon1;
        std::vector <Node *> neighbors = sons[0]->getNeighbors();
        for (unsigned int i=0; i<neighbors.size(); i++) {
            if (neighbors[i]==node) {
                directionSon0 = i;
            }
        }
        neighbors = sons[1]->getNeighbors();
        for (unsigned int i=0; i<neighbors.size(); i++) {
            if (neighbors[i]==node) {
                directionSon1 = i;
            }
        }
        computeCoalCountsFromSons (spTree, sons, 
                                   speciesIDs[id][0], 
                                   speciesIDs[idSon0][directionSon0], 
                                   speciesIDs[idSon1][directionSon1],
                                   coalCounts[id][0],
                                   coalCounts[idSon0][directionSon0], 
                                   coalCounts[idSon1][directionSon1]
                                    );                
        return;
	}	
}


/*****************************************************************************
 * Utilitary functions to initialize vectors of counts at leaves.
 ****************************************************************************/
void initializeCountVectors(std::vector< std::vector< std::vector<unsigned int> > > &vec) {
    int size = vec[0].size();
    for (unsigned int i = 0 ; i < 3 ; i++ ) { //I should not be necessary to initialize all vectors, i==0 should be enough
        initializeCountVector(vec[i]);
    }
    return;
}


void initializeCountVector(std::vector<std::vector<unsigned int> >  &vec) {
    for (unsigned int j = 0 ; j < vec.size() ; j++ ) {
         for (unsigned int i = 0 ; i < 2 ; i++ ) {
             vec[j][i] = 0;
         }
    }
    return;
}

/*****************************************************************************
 * Utilitary functions to increment vectors of counts.
 ****************************************************************************/

void incrementOutCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos) {
    vec[pos][1] = vec[pos][1] +1;
    return;
}

void incrementInCount(std::vector< std::vector<unsigned int> > & vec, const unsigned int pos) {
    vec[pos][0] = vec[pos][0] +1;
    return;
}


/*****************************************************************************
 * Computes a vector of counts from two son vectors, and assigns species ID to the
 * father node.
 ****************************************************************************/

void computeCoalCountsFromSons (TreeTemplate<Node> & tree, std::vector <Node *> sons, 
                                unsigned int & rootSpId, 
                                const unsigned int & son0SpId,
                                const unsigned int & son1SpId,
                                std::vector< std::vector<unsigned int> > & coalCountsFather,
                                std::vector< std::vector<unsigned int> > & coalCountsSon0,
                                std::vector< std::vector<unsigned int> > & coalCountsSon1)
{
    
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
    for (unsigned int i= 0 ; i < coalCountsFather.size() ; i++) { //for each branch of the sp tree
        for (unsigned int j= 0 ; j < 2 ; j++) { //for incoming and outgoing counts
            coalCountsFather[i][j] = coalCountsSon0[i][j] + coalCountsSon1[i][j];
        }
        //std::cout << "coalCountsFather[i][j] For Sp branch "<<i<<" Num coal in: "<< coalCountsFather[i][0] << " Num coal out: "<< coalCountsFather[i][1]<<std::endl;
    }
    
    
    
    Node temp0 = *(tree.getNode(son0SpId));
    Node temp1 = *(tree.getNode(son1SpId));
    
    while (a!=b) { //There is ILS!
        if (a>b) {
            recoverILS(temp0, a, olda, coalCountsFather);
        }
        else {
            recoverILS(temp1, b, oldb, coalCountsFather);
        }
    }
    rootSpId = a;
    
    return;    
}


/*****************************************************************************
 * This function recovers ILS by comparing a subtree in a gene tree to
 * a species tree.
 * WARNING: MAY NEED TO CHANGE: NEED TO TAKE THE LOWEST POSSIBLE BRANCH FOR DOING THE COALESCENCE
 ****************************************************************************/
void recoverILS(Node & node, int & a, int & olda, 
                std::vector <std::vector < unsigned int > > &vec) {
    olda=a;
    Node* nodeA;
    if (node.hasFather()) {
        nodeA = node.getFather();
    }
    else {
        std::cout <<"Problem in recoverILS , nodeA has no father"<<std::endl; 
    }
    a = nodeA->getId();
    node = *nodeA;
    
    incrementInCount(vec, olda);
    incrementOutCount(vec, a);
    
    return;
}

/*****************************************************************************
 * Computes the likelihood using our coalescence model, 
 * given a vector of vector giving, for each branch of the species tree,
 * the number of incoming lineages, and the number of outgoing lineages.
 * Formula from Degnan and Salter (2005), Evolution 59(1), pp. 24-37.
 * 3 versions of the function:
 * - working for lots of gene families at once
 * - working for one gene family
 * - working for one branch of one gene family
 ****************************************************************************/


double computeCoalLikelihood (std::vector < std::vector<std::vector<unsigned int> > > vec, std::vector < double > CoalBl ) 
{
    double logLk = 0;
    for (unsigned int i = 0 ; i < vec.size() ; i++ ) 
    {
        logLk += computeCoalLikelihood (vec[i],  CoalBl );
    }
    return logLk;
}


double computeCoalLikelihood (std::vector < std::vector<unsigned int> > vec, std::vector < double > CoalBl ) 
{
    double logLk = 0;
    for (unsigned int i = 0 ; i < vec.size() ; i++ ) 
    {
        logLk += computeCoalLikelihood (vec[i],  CoalBl[i] );
    }
    return logLk;
}

double computeCoalLikelihood ( std::vector<unsigned int>  vec, double CoalBl ) 
{
    
    double logLk = 0;
    double prod;
    int kminus1;

    double v = (double) vec[0];
    double u = (double) vec[1];

    
    for (unsigned int k = vec[0] ; k <= vec[1] ; k++ ) 
    {
        prod = 1.0;
        kminus1 = k-1;
        double dy;
        for (unsigned int y = 0 ; y <= kminus1 ; y++ ) {
            dy = double(y);
            prod *= (v + dy) * (u - dy) / ( u + dy ) ;
        }
        if (prod == 0) {
            prod = 1;
        }
        double dk = (double)k;
        double dkminus1 = (double)kminus1;
        double dkminusv = dk - v ;
       /* std::cout << "CoalBl "<<CoalBl<<std::endl;
        std::cout << "prod "<< prod <<std::endl;
        std::cout <<"fact "<<NumTools::fact(vec[0])<<std::endl;
        std::cout <<"fact2 "<<NumTools::fact(k - vec[0])<<std::endl;
        std::cout <<"exp (-k * (k-1) * CoalBl/2 ) "<<exp (-(double)k * ((double)k-1.0) * CoalBl/2.0  )<<std::endl;
        std::cout <<"(2*k -1) "<<(2.0*(double)k -1)<<std::endl;
        std::cout <<"( (-1)^(k-vec[0]) ) "<<( pow(-1.0, ( (double)k-(double)vec[0]) ) )<<std::endl;
        std::cout <<"( NumTools::fact(vec[0]) * NumTools::fact(k - vec[0]) * (vec[0] + kminus1 ) ) "<<( (double)NumTools::fact(vec[0]) * (double)NumTools::fact(k - vec[0]) * ((double)vec[0] + kminus1 ) )<<std::endl;*/
        
        logLk += ( exp (-dk * dkminus1 * CoalBl/2.0 ) ) * (2.0*dk -1.0) * ( pow(-1.0, dkminusv ) ) / ( NumTools::fact(v) * NumTools::fact(dkminusv) * (v + dkminus1 ) ) * prod;
        //std::cout << "logLk: " << logLk<<std::endl;
    }
    return log(logLk);
}







/*****************************************************************************
 * This function performs a preorder tree traversal in order to fill vectors of counts. 
 * When used after the postorder tree traversal function, counts for all rootings are computed.
 * coalCounts contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * node: node of the gene tree used as son of the root.
 ****************************************************************************/

void computeSubtreeCoalCountsPreorder(TreeTemplate<Node> & spTree, 
                                      TreeTemplate<Node> & geneTree, 
                                      Node * node, 
                                      const std::map<std::string, std::string > & seqSp, 
                                      const std::map<std::string, int > & spID, 
                                      std::vector < std::vector< std::vector< std::vector<unsigned  int > > > > & coalCounts, 
                                      std::vector<double> & bls, 
                                      std::vector <std::vector<unsigned int> > & speciesIDs, 
                                      int sonNumber, 
                                      std::map <double, Node*> & LksToNodes) {
    computeRootingCoalCounts(spTree, node, 
                             coalCounts, bls, speciesIDs, 
                             sonNumber, LksToNodes);
    if (node->isLeaf()) {
        return; 
    }
    Node * son;
    if (sonNumber==1) {
        son= node->getSon(1);
    }
    else {
        son= node->getSon(0);
    }
    //  for (int i = 0; i< sons.size(); i++){
    for (unsigned int j =0; j<son->getNumberOfSons(); j++) {
        computeSubtreeCoalCountsPreorder(spTree, geneTree, 
                                         son, seqSp, spID, 
                                         coalCounts, 
                                         bls, 
                                         speciesIDs, j, LksToNodes);
    }
    //  }
	return;
	
}




/*****************************************************************************
 * This function computes the Coalescent counts of a rooting. 
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 * node is a node of the gene tree. 
 ****************************************************************************/
void computeRootingCoalCounts(TreeTemplate<Node> & spTree, 
                              Node * node, 
                              std::vector < std::vector< std::vector< std::vector<unsigned  int > > > > & coalCounts, 
                              const std::vector< double> & bls, 
                              std::vector <std::vector<unsigned int> > & speciesIDs, 
                              int sonNumber, 
                              std::map <double, Node*> & LksToNodes) {
    int geneNodeId = node->getId();
    int directionForFather;
    std::vector <Node*> nodes;

    nodes.push_back(node->getFather());

    if (sonNumber==0) { //If sonNumber==0, the subtree we're interested in is composed of son 1 and father of node.
        nodes.push_back(node->getSon(1));
    }
    else { //If sonNumber==1, the subtree we're interested in is composed of son 0 and father of node.
        nodes.push_back(node->getSon(0));
    }
    if (node->getFather()->getSon(0)==node) {
        directionForFather = 1; //node #1 is son 0: from node, the useful subtree is in 1
    }
    else {
        directionForFather = 2; //node #2 is son 1: from node, the useful subtree is in 2
    }

    int idNode0, idNode1;
    idNode0 = nodes[0]->getId();
    idNode1 = nodes[1]->getId();
    unsigned int directionNode0, directionNode1;
    directionNode0 = directionForFather;
    directionNode1 = 0;


    unsigned int directionSon0, directionSon1;
    directionSon0 = sonNumber+1;
    directionSon1 = 0;
    
    computeCoalCountsFromSons (spTree, nodes, 
                               speciesIDs[geneNodeId][directionSon0], 
                               speciesIDs[idNode0][directionNode0], 
                               speciesIDs[idNode1][directionNode1],
                               coalCounts[geneNodeId][directionSon0],
                               coalCounts[idNode0][directionNode0], 
                               coalCounts[idNode1][directionNode1]
                               );  
    
    //Now we have the counts for the upper subtree, 
    //as well as the counts of the lower subtree (which we already had)
    //We can thus compute the total likelihood of the rooting.
    
    std::vector <Node*> sons;
    sons.push_back(node);
    
    sons.push_back(node->getSon(sonNumber));
    int idSon0 = geneNodeId;
    int idSon1 = sons[1]->getId();
    
    double rootLikelihood = 0.0;
    unsigned int rootSpId;
    int rootDupData = 0;
    std::vector< std::vector<unsigned  int > > rootCounts = coalCounts[geneNodeId][directionSon0];
    

    /*
    for (unsigned int i = 0 ; i < coalCounts[geneNodeId][directionSon0].size() ; i++) {
        std::cout << "BEFORE For branch "<< geneNodeId<<"; Sp Branch "<<i<<" Num coal in: "<< coalCounts[geneNodeId][directionSon0][i][0] << " Num coal out: "<< coalCounts[geneNodeId][directionSon0][i][1]<<std::endl;
    }
*/
    computeCoalCountsFromSons (spTree, sons, 
                               rootSpId, 
                               speciesIDs[geneNodeId][directionSon0], 
                               speciesIDs[idSon1][0],
                               rootCounts,
                               coalCounts[geneNodeId][directionSon0], 
                               coalCounts[idSon1][0]
                               );  

    
    
    //Adding the two tables, from the upper subtree, and from the lower subtree
 /*   for (unsigned int i = 0 ; i < spTree.getNumberOfNodes() ; i++ ) {
        coalCounts[geneNodeId][directionSon0][i][0] = coalCounts[geneNodeId][directionSon0][i][0] + coalCounts[idSon1][0][i][0];
        coalCounts[geneNodeId][directionSon0][i][1] = coalCounts[geneNodeId][directionSon0][i][1] + coalCounts[idSon1][0][i][1];
    }*/
    
    //Add the starting lineage at the root
    for (unsigned int i = 0 ; i < spTree.getNumberOfNodes() ; i++ ) {
        if (rootCounts[i][0]==0 && rootCounts[i][1] !=0)
        {
            rootCounts[i][0]=1; 
            break;
        }
    }
    /*
    for (unsigned int i = 0 ; i < coalCounts[geneNodeId][directionSon0].size() ; i++) {
        std::cout << "For branch "<< geneNodeId<<"; Sp Branch "<<i<<" Num coal in: "<< rootCounts[i][0] << " Num coal out: "<< rootCounts[i][1]<<std::endl;
    }
     */

    //What to put?
    rootLikelihood = computeCoalLikelihood ( rootCounts, bls ) ;

 /*   
    computeConditionalLikelihoodAndAssignSpId(spTree, sons, 
                                              rootLikelihood, 
                                              likelihoodData[idSon0][directionSon0], 
                                              likelihoodData[idSon1][directionSon1], 
                                              lossRates, duplicationRates, 
                                              rootSpId, speciesIDs[idSon0][directionSon0], 
                                              speciesIDs[idSon1][directionSon1], 
                                              rootDupData, dupData[idSon0][directionSon0], 
                                              dupData[idSon1][directionSon1], true);*/
    // std::cout <<"LK FOUND "<<rootLikelihood<<std::endl;
    while (LksToNodes.find(rootLikelihood)!=LksToNodes.end()) {
        // std::cout <<"changing rootLikelihood !!!!!!!!!!!!!!!!!!!"<<std::endl;
        rootLikelihood+=SMALLPROBA;
    }
    LksToNodes[rootLikelihood]=node->getSon(sonNumber);
    return;
}









/*****************************************************************************
 * This function aims at finding the most likely coalescent reconciliation, 
 * using a double recursive tree traversal. 
 * The first traversal is post-order, and then the second traversal is pre-order.
 * This is a modification of an algorithm quickly explained in 
 * Chen, Durand, Farach-Colton, J. Comp. Biol. pp429-447, 2000.
 * Conditional likelihoods are recorded in a table. 
 * This table has (number of nodes) elements, and for each node, 
 * contains three conditional likelihoods. 
 * The table is thus (number of nodes)*3 cells. For each node i, 
 * likelihoodData[i][j] contains the conditional likelihood of the subtree 
 * having its root in subtree opposite neighbour j of node i.
 * Node species IDs are also recorded in a (number of nodes)*3 cells table.
 * The boolean "fillTables" is here to tell whether we want to update the vectors num*lineages.
 ****************************************************************************/

double findMLCoalReconciliationDR (TreeTemplate<Node> * spTree, 
                               TreeTemplate<Node> * geneTree, 
                               std::map<std::string, std::string > seqSp,
                               std::map<std::string, int > spID,
                               std::vector< double> coalBl, 
                               int & MLindex, 
                               std::vector < std::vector<std::vector<unsigned int> > > coalCounts,
                               std::set <int> &nodesToTryInNNISearch, 
                               bool fillTables)
{
	if (!geneTree->isRooted()) {
		std::cout << TreeTools::treeToParenthesis (*geneTree, true)<<std::endl;
		std::cout <<"!!!!!!gene tree is not rooted in findMLReconciliationDR !!!!!!"<<std::endl;
        MPI::COMM_WORLD.Abort(1);
		exit(-1);
    }		
	std::vector <double> nodeData(3, 0.0);
	std::vector <std::vector<double> > likelihoodData(geneTree->getNumberOfNodes(), nodeData);
    
    std::vector <int> nodeSpId(3, 0);
    std::vector <std::vector<int> > speciesIDs(geneTree->getNumberOfNodes(), nodeSpId);
    std::vector <std::vector<int> > dupData = speciesIDs;
    
    double initialLikelihood;
    //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
    std::map <double, Node*> LksToNodes;
    
    /* std::cout << TreeTools::treeToParenthesis (*geneTree, true)<<std::endl;*/
    //  std::cout << "IN findMLReconciliationDR "<<TreeTools::treeToParenthesis (*spTree, true)<<std::endl;
    
    //	std::cout <<"CLOCKBEFOREPOSTORDER "<<clock()<<std::endl;
	Node * geneRoot = geneTree->getRootNode();
    /*  std::cout <<"root number :"<<geneRoot->getId()<<std::endl;
     std::cout << TreeTools::treeToParenthesis (*geneTree, true)<<std::endl;*/
    initialLikelihood = computeSubtreeLikelihoodPostorder(*spTree, *geneTree, 
                                                          geneRoot, seqSp, spID, 
                                                          likelihoodData, lossRates, 
                                                          duplicationRates, speciesIDs, dupData);
    /*    VectorTools::print(duplicationRates);
     std::cout << TreeTools::treeToParenthesis (*spTree, true)<<std::endl;
     std::cout << TreeTools::treeToParenthesis (*geneTree, true)<<std::endl;
     std::cout << "initialLikelihood: "<< initialLikelihood << std::endl;*/
    //  std::cout <<"CLOCKAFTERPOSTORDER "<<clock()<<std::endl;
    //std::cout <<"Postorder tree traversal over: Initial likelihood: "<<initialLikelihood<<std::endl;
	//computeSubtreeLikelihoodPrefix then computes the other conditional likelihoods, and also returns the best rooting.
    std::vector <Node *> sons = geneRoot->getSons();
    
    if (sons.size()!=2) {
        std::cerr <<"Error: "<<sons.size()<< "sons at the root!"<<std::endl; 
    }
    LksToNodes[initialLikelihood]=sons[0];
    //We fill the likelihood and species ID data for the root node.
    //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
    //with other nodes.
    likelihoodData[geneRoot->getId()][1] = likelihoodData[geneRoot->getSon(1)->getId()][0];
    likelihoodData[geneRoot->getId()][2] = likelihoodData[geneRoot->getSon(0)->getId()][0];
    speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon(1)->getId()][0];
    speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon(0)->getId()][0];
    dupData[geneRoot->getId()][1] = dupData[geneRoot->getSon(1)->getId()][0];
    dupData[geneRoot->getId()][2] = dupData[geneRoot->getSon(0)->getId()][0];
    
    /*  std::cout <<"likelihoodData[geneRoot->getId()][1]"<<likelihoodData[geneRoot->getId()][1]<<std::endl;
     std::cout <<"likelihoodData[geneRoot->getId()][2]"<<likelihoodData[geneRoot->getId()][2]<<std::endl;
     
     for (int i=0 ; i<geneTree->getNumberOfNodes(); i++) {
     std::cout <<"ID "<<i<<"likelihoodData[ID][0]"<<likelihoodData[i][0]<<std::endl;
     }
     */
    
    for (unsigned int i = 0; i< sons.size(); i++){
        for (unsigned int j =0; j<sons[i]->getNumberOfSons(); j++) {
            computeSubtreeLikelihoodPreorder(*spTree, *geneTree, 
                                             sons[i], seqSp, spID, 
                                             likelihoodData, 
                                             lossRates, duplicationRates, 
                                             speciesIDs, dupData, j, LksToNodes);
        }
    }
    
    
    /* 
     std::cout <<"Printing all rooting likelihoods as found by the DR tree traversal"<<std::endl;
     std::map<double, Node*>::iterator it;
     
     for ( it=LksToNodes.begin() ; it != LksToNodes.end(); it++ )
     std::cout << (*it).second->getId() << " => " << (*it).first << std::endl;
     */
    //  std::cout << "IN findMLReconciliationDR: "<<TreeTools::treeToParenthesis (*geneTree, false)<<std::endl;
    
    vector<Node*> nodes = geneTree->getNodes();
    for (unsigned int i = 0 ; i < nodes.size() ; i++ ) {
        if (nodes[i]->hasNodeProperty("outgroupNode") ) {
            nodes[i]->deleteNodeProperty("outgroupNode");
            break;
        }
    }
    //    std::cout << TreeTools::treeToParenthesis (*geneTree, false)<<std::endl;
    
    LksToNodes.rbegin()->second->setNodeProperty("outgroupNode", BppString("here") );
    //    std::cout << TreeTools::treeToParenthesis (*geneTree, false)<<std::endl;
    
    geneTree->newOutGroup(LksToNodes.rbegin()->second); //uncomment that if you want to keep gene family trees fixed except for the root
    //    std::cout << TreeTools::treeToParenthesis (*geneTree, false)<<std::endl;
    
    
    if (fillTables) {
        
        //Now the best root has been found. I can thus run a function with this best root to fill all the needed tables. This additional tree traversal could be avoided.
        //To this end, the needed tables should be filled by the postfix and prefix traversals. This has not been done yet.
        //resetting
        speciesIDs = std::vector<std::vector<int> > (geneTree->getNumberOfNodes(), nodeSpId);
        
        dupData = speciesIDs;
        // Getting a well-rooted tree
        TreeTemplate<Node > * tree = geneTree->clone();
        
        //tree->newOutGroup(LksToNodes.rbegin()->second->getId());
        
        
        //std::cout << TreeTools::treeToParenthesis (*tree, true)<<std::endl;
        
        nodesToTryInNNISearch.clear();
        
        //Resetting numLineages std::vectors
        resetVector(num0lineages);
        resetVector(num1lineages);
        resetVector(num2lineages);
        
        
        // std::cout <<"HERE_rooted_tree "<<TreeTools::treeToParenthesis (*tree, true)<<std::endl;
        //  std::cout <<"HERE_rooted_tree2 "<<TreeTools::treeToParenthesis (*geneTree, true)<<std::endl;
        //  std::cout <<"HERE_SP_tree "<<TreeTools::treeToParenthesis (*spTree, true)<<std::endl;
        computeNumbersOfLineagesFromRoot(spTree, tree, 
                                         tree->getRootNode(), 
                                         seqSp, spID, 
                                         num0lineages, num1lineages, 
                                         num2lineages, speciesIDs, 
                                         dupData, nodesToTryInNNISearch);
        delete tree;
        
    }

    //We return the best likelihood
    MLindex = LksToNodes.rbegin()->second->getId();
    // std::cout <<"Bestlikelihood"<< LksToNodes.rbegin()->first<<std::endl;
    
	return LksToNodes.rbegin()->first;
    
    
    
}


























/*******************************************************************************/

//Creating the compressed vector, and the map giving the weights for the "site" patterns
double CoalBranchLikelihood::initModel()
{
    string pattern = "";
    compressedVec_.clear();
    std::vector <unsigned int> v (2, 1);
    for (unsigned int i = 0 ; i < vec_.size() ; i++) {
        pattern = TextTools::toString(vec_[i][0]) + TextTools::toString(vec_[i][1]);
        if ( patternToWeights_.find(pattern)!= patternToWeights_.end() ) {
            patternToWeights_[pattern] += 1; 
        }
        else {
            patternToWeights_.insert( pair<std::string, unsigned int >(pattern, 1) );
            v[0] = vec_[i][0];
            v[1] = vec_[i][1];
            compressedVec_.push_back(v);
        }
    }
    lks_.resize(compressedVec_.size(), 0.0);
    for(std::map<std::string, unsigned int >::iterator it = patternToWeights_.begin(); it != patternToWeights_.end(); it++){
        std::cout <<"\t\tit->first: "<< it->first << " : "<< it->second <<std::endl;
    }
    double n12 = (double) patternToWeights_["12"];
    double n22 = (double) patternToWeights_["22"];
    double estimate;
    std::cout <<"n12 "<< n12 <<" n22 "<< n22 <<std::endl;
    if (n22 ==0) {
        n22 = 1;
    }
    if (n12 ==0) {
        n12 = 1;
    }
    estimate = log((n12+n22)/n22);
    
    std::cout <<"\t\tAnalytical estimate: "<< estimate<<std::endl;
    
    return estimate;
}

/*******************************************************************************/
/*
void CoalBranchLikelihood::computeAllTransitionProbabilities()
{
    double l = getParameterValue("BrLen"); 
    
    //Computes all pxy once for all:
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
        VVdouble * pxy__c = & pxy_[c];
        RowMatrix<double> Q = _model->getPij_t(l * _rDist->getCategory(c));
        for(unsigned int x = 0; x < nbStates_; x++)
        {
            Vdouble * pxy__c_x = & (* pxy__c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
                (* pxy__c_x)[y] = Q(x, y);
            }
        }
    } 
}
*/
/*******************************************************************************/

void CoalBranchLikelihood::computeLogLikelihood()
{    
    double l = getParameterValue("BrLen"); 
    lnL_ = 0;   
    string pattern = "";
    for (unsigned int i = 0 ; i < compressedVec_.size() ; i++) {
        pattern = TextTools::toString(vec_[i][0]) + TextTools::toString(vec_[i][1]);
        lks_[i] = computeCoalLikelihood(compressedVec_[i], l);
        lnL_ -= patternToWeights_[pattern] * lks_[i];
    }

    return;
    /*
    _arrayTmp = *_array1;
    lnL_ = 0;
    for(unsigned int i = 0; i < _array1->size(); i++)
    {
        VVdouble * arrayTmp_i = & _arrayTmp[i];
        const VVdouble * array2_i = & (*_array2)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
            Vdouble * arrayTmp_i_c = & (*arrayTmp_i)[c];
            const Vdouble * array2_i_c = & (*array2_i)[c];
            VVdouble *pxy__c = & pxy_[c];
            for(unsigned int x = 0; x < nbStates_; x++)
            {
                Vdouble *pxy__c_x = & (*pxy__c)[x];
                double likelihood = 0;
                for(unsigned int y = 0; y < nbStates_; y++)
                {
                    likelihood += (*pxy__c_x)[y] * (*array2_i_c)[y];
                }
                (*arrayTmp_i_c)[x] *= likelihood;
            }
        }
    }
    
    vector<double> la(_array1->size());
    for(unsigned int i = 0; i < _array1->size(); i++)
    {
        VVdouble * arrayTmp_i = & _arrayTmp[i];
        double Li = 0;
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
            Vdouble * arrayTmp_i_c = & (*arrayTmp_i)[c];
            double rc = _rDist->getProbability(c);
            for(unsigned int x = 0; x < nbStates_; x++)
            {
                //Li += rc * _model->freq(x) * (* arrayTmp_i_c)[x];
                //freq is already accounted in the array
                Li += rc * (* arrayTmp_i_c)[x];
            }
        }
        la[i] = _weights[i] * log(Li);
    }
    sort(la.begin(), la.end());
    for(unsigned int i = _array1->size(); i > 0; i--)
        lnL_ -= la[i-1];*/
    
}


/*******************************************************************************/





