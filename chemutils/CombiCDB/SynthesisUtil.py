"""Assorted support classes and utilities for the synthesis design modules"""

import sys, os;
from optparse import OptionParser;
import random;
#from sets import Set;
from openeye.oechem import OECreateIsoSmiString, OECreateCanSmiString;
from openeye.oechem import OEDetermineRingSystems;

from chemutils.Common.Const import SMILES_MOL_DELIM;
from chemutils.Common.Util import ProgressDots, molBySmiles, generateFingerprint;
from chemutils.Common.Util import splitCompositeMol, splitCompositeMolToSmilesList;
from chemutils.CombiCDB.ReactionModel import REACTANTS, REAGENT, PRODUCTS;

from chemutils.CombiCDB.Const import MAX_TRIES;
from chemutils.CombiCDB.Util import log;

class BaseReactantPool:
    """Base class defining what a reactant pool must do to satisfy the
    needs of the RetroSynthesis class.
    """
    def search(self, queryMol, resultsLimit):
        """Given a query molecule object, return the top scored / similar results
        out of the reactant pool, up to resultsLimit items.  The results should be
        2-ples each consisting of a result molecule object and the similarity score.
        """
        raise NotImplementedError();

    def __contains__(self, queryMol):
        """Do an exact search / lookup.  Return true if and only if the queryMol is found within the reactant pool.
        Don't just take the first result from the "search" method, because it turns out multiple results
        may yield 1.0 similarity, even if they're not quite exactly what you were looking for.
        """
        raise NotImplementedError();

class SimpleReactantPool(BaseReactantPool):
    """Simple reactant pool just initialized in memory by a list of SMILES strings,
    using simple Tanimoto / Tversky similarity measure.  
    This wouldn't scale well for very large pools, but should be convenient for small tests sets.
    """
    # List of 3-ples, each consisting of a reactant 
    #   - isomeric SMILES string, 
    #   - molecule object and 
    #   = fingerprint representation
    reactantPool = None;
    
    def __init__(self, reactantSmilesList):
        from FINGER import tversky; # Only import as needed

        self.reactantPool = [];
        for reactantSmiles in reactantSmilesList:
            reactantMol = molBySmiles(reactantSmiles);
            fingerprint = generateFingerprint(reactantSmiles);
            self.reactantPool.append( (OECreateIsoSmiString(reactantMol), reactantMol, fingerprint) );

    def search(self, queryMol, resultsLimit):
        querySmi = OECreateIsoSmiString(queryMol);
        #print >> sys.stderr, querySmi
        queryFP = generateFingerprint(querySmi);
        scoredPool = [];
        for reactantSmi, reactantMol, fingerprint in self.reactantPool:
            score = tversky( queryFP, fingerprint, 1.0, 1.0); # Set to 0.1, 0.9 to favor query being a super-structure of reactants
            scoredPool.append( (score, reactantMol) );
        scoredPool.sort();  # Sort by score
        scoredPool.reverse();   # Want it in desc. order
        return scoredPool[0 : resultsLimit];

    def __contains__(self, queryMol):
        querySmi = OECreateIsoSmiString(queryMol);
        for reactantSmi, reactantMol, fingerprint in self.reactantPool:
            if reactantSmi == querySmi:
                return True;
        return False;

class ChemDBReactantPool(BaseReactantPool):
    """Use the ChemDB database as the reactant pool, 
    using similarity search servers and the database to identify similar compounds.
    To regulate performance, should initialize with a minimum similarity threshold 
    for returned results.
    """
    searchModel = None;
    chemicalSearch = None;
    
    def __init__(self, similarityThreshold):
        # Only import external modules if going to use them
        from chemutils.Common.ChemicalSearch import ChemicalSearch;
        from chemutils.Common.Model import ChemicalSearchModel;

        self.searchModel = ChemicalSearchModel();
        self.searchModel.similarityThreshold = similarityThreshold;
        self.includeSourceData = False; # Save on performance, don't want this information
        self.chemicalSearch = ChemicalSearch();
    
    def search(self, queryMol, resultsLimit):
        print >> sys.stderr, "Search: ", OECreateIsoSmiString(queryMol);

        self.searchModel.maxResults = resultsLimit;
        
        self.searchModel.addSimilarMol(queryMol);
        chemicalModels = self.chemicalSearch.findChemicals( self.searchModel );
        self.searchModel.similarMols.pop(); # Revert to original state

        scoredMols = [];
        for chemical in chemicalModels:
            mol = molBySmiles( chemical["can_smiles"] );
            scoredMols.append( (chemical["similarityScore"], mol) );

        return scoredMols;        

    def __contains__(self, queryMol):
        querySmi = OECreateIsoSmiString(queryMol);
        self.searchModel.addDiscreteCriteria( "can_smiles", [querySmi] );
        chemicalModels = self.chemicalSearch.findChemicals( self.searchModel );
        self.searchModel.popDiscreteCriteria(); # Revert to original state
        
        return len(chemicalModels) > 0;

class BaseScoreAggregator:
    """Base class for deriving some aggrevate score from a list of scores"""
    def __call__(self, scoreList):
        raise NotImplementedError();

class MinimumScoreAggregator(BaseScoreAggregator):
    """Score aggregator that seeks to minimize the maximum loss.
    Work by simply taking the minimum score in the list as the total score
    (list is only as good as the weakest component).
    """
    def __call__(self, scoreList):
        return min(scoreList);

class AverageScoreAggregator(BaseScoreAggregator):
    """Score aggregator that summarizes the list of scores by taking the average"""
    def __call__(self, scoreList):
        sum = 0;
        n = len(scoreList);
        for score in scoreList:
            sum += score;
        return sum / n;

class BaseRetroProductSelector:
    """Base class to define how a retro product selector should behave.
    Given a collection of retro-reagents and a target molecule, it will apply the retro-reagents
    to the target molecule.  Any of the retro-products (precursors) it deems
    worthy will be returned, along with a reference to the retro-reagent that created it.
    """
    def __init__(self):
        pass
    
    def __call__(self, retroReagentList, targetMol):
        """Should return a list of 3-ples, consisting of the accepted / selected retro product information
            - retroProduct molecule object
            - retroReagent object reference
            - iReagent index of the retro reagent in the list of reagents
        """
        raise NotImplementedError();

class GreedyRetroProductSelector(BaseRetroProductSelector):
    """Generates all posible retro-products (1-level deep), but only returns the top K results,
    with precedence defined by a respective similarity search measure on available
    reactants in the reactant pool.  Also needs a score aggregator in case a retro-reaction
    generates multiple products (combinatorial reactions) and it needs to generate
    a single summary score for that collection.
    """
    reactantPool = None;
    scoreAggregator = None;
    returnCount = None;
    
    def __init__(self, reactantPool, scoreAggregator, returnCount=1):
        BaseRetroProductSelector.__init__(self);
        self.reactantPool = reactantPool;
        self.scoreAggregator = scoreAggregator;
        self.returnCount = returnCount;

    def __call__(self, retroReagentList, targetMol):
        scoredRetroProductsByReagent = [];

        for iReagent, retroReagent in enumerate(retroReagentList):
            #print >> sys.stderr, "Reagent", iReagent, retroReagent["depict_smiles"],;
            retroProducts = retroReagent( [targetMol] );
            #print >> sys.stderr, len(retroProducts);

            for retroProduct in retroProducts:
                componentProducts = splitCompositeMol(retroProduct);
                componentScores = [];
                for componentProduct in componentProducts:
                    topResults = self.reactantPool.search( componentProduct, 1 );
                    topScore = 0.0; # 0 score when nothing found
                    if len(topResults) > 0:
                        (topScore, topMol) = topResults[0];
                        #log.debug( OECreateIsoSmiString(componentProduct), topScore, OECreateIsoSmiString(topMol) );
                    componentScores.append( topScore );
                aggregateScore = self.scoreAggregator( componentScores );
                scoredRetroProductsByReagent.append( (aggregateScore, (retroProduct, retroReagent, iReagent)) );
        scoredRetroProductsByReagent.sort();
        scoredRetroProductsByReagent.reverse(); # Want desc. order by score
        
        topRetroProductsByReagent = [];
        for iProduct in xrange( min(self.returnCount, len(scoredRetroProductsByReagent)) ):
            topRetroProductsByReagent.append( scoredRetroProductsByReagent[iProduct][1] );    # Copy only the top results, without scores

        return topRetroProductsByReagent;

class ExhaustiveRetroProductSelector(BaseRetroProductSelector):
    """Generates all posible retro-products (1-level deep), and just returns all of them
    to indicate an exhaustive exploration of the search space.
    """
    def __init__(self):
        BaseRetroProductSelector.__init__(self);

    def __call__(self, retroReagentList, targetMol):
        retroProductsByReagent = [];

        for iReagent, retroReagent in enumerate(retroReagentList):
            #print >> sys.stderr, "Reagent", iReagent, retroReagent["depict_smiles"],;
            retroProducts = retroReagent( [targetMol] );
            #print >> sys.stderr, len(retroProducts);

            for retroProduct in retroProducts:
                retroProductsByReagent.append( (retroProduct, retroReagent, iReagent) );

        return retroProductsByReagent;


def synthesisTreeTableFromReactionSteps(reactionSteps):
    """Given a collection of reaction steps 
    (3-ples with a list of reactants, reagent, list of products), construct
    a tree table based on the implied series of reactants -> products.
    
    Presumably the "reactants" and "products" are SMILES strings, but that
    is not strictly encorced.  Just assume they're something that can be
    used to key dictionary objects.
    """
    # Organize the reaction steps by their targetSmiles.
    # To make the most sense, each step should only have one product 
    #   (to yield a tree, and not a hypergraph), but that can still be compensated for
    #   with redundant tree entries
    reactionStepsByProduct = dict();
    for reactionStep in reactionSteps:
        for product in reactionStep[PRODUCTS]:
            reactionStepsByProduct[product] = reactionStep;
    
    # Look for "root" nodes of the synthesis tree, by finding target / products
    #   which are not precursor reactants for any other steps
    rootProductSet = reactionStepsByProduct.keys();
    for reactionStep in reactionSteps:
        for reactant in reactionStep[REACTANTS]:
            if reactant in rootProductSet:
                rootProductSet.remove(reactant);

    # First construct the tree as a variable width 2-D table object
    treeTable = [];
    rootProductList = list(rootProductSet);
    rootProductList.sort(); # Produce in a consistent order
    for rootProduct in rootProductList:
        addSynthesisTreeTableSubRows( rootProduct, treeTable, reactionStepsByProduct );

    return treeTable;

def addSynthesisTreeTableSubRows( product, treeTable, reactionStepsByProduct, depth=0 ):
    """Recursive function to do most of the work of building the synthesis tree table
    """
    row = [None]*depth; # Extend row to current recursion depth
    row.append( product );
    treeTable.append(row);  # Make sure this comes before recursion, for pre-traversal result

    if product in reactionStepsByProduct:
        reactionStep = reactionStepsByProduct[product];
        reagent = reactionStep[REAGENT];
        row.append(reagent);    # Record what reagent can be used to deconstruct this product

        # Recurse upon all of the precursor reactants for this product
        for reactant in reactionStep[REACTANTS]:
            addSynthesisTreeTableSubRows( reactant, treeTable, reactionStepsByProduct, depth+1 );
    else:
        # Base case, got to a leaf node, no more precursors for the given product to recurse upon
        pass;


def synthesisStartingMaterials( reactionSteps ):
    """Given a list of ordered reaction steps, representing a synthesis plan,
    figure out which reactant molecules are needed as starting materials
    (as opposed to molecules carried over from products of previous reactions).
    Actually returns a 2-ple, the other element being a set of SMILES strings
    representing every other molecule that appears in the synthesis which is
    not a starting material (must be a final, intermediate or side product).
    """
    startingMaterials = set();
    productSmiSet = set();
    for reactionStep in reactionSteps:
        for reactant in reactionStep[REACTANTS]:
            reactantSmi = OECreateIsoSmiString(reactant);
            if reactantSmi not in productSmiSet:
                startingMaterials.add(reactantSmi);
        for product in reactionStep[PRODUCTS]:
            componentProducts = splitCompositeMol(product, retainCounterIons=True);
            for component in componentProducts:
                productSmiSet.add( OECreateIsoSmiString(component) );
    return (startingMaterials, productSmiSet);

def analyzeProductMixture(productList, target=None):
    """Given a list of product molecules, figure out how many
    different products there are, both considering and not considering stereoisomers as different.
    If an intended target product is provided, also indicate how many unintended side products
    not matching the target were observed.  Considering complete different products and those
    that are just stereoisomers.
    
    Return these stats as a tuple:
    (distinctProducts, distinctIsmProducts, unintendedProducts, unintendedStereoisomers)
    
    """
    targetCanSmi = None;
    targetIsmSmi = None;
    if target is not None:
        targetCanSmi = OECreateCanSmiString(target);
        targetIsmSmi = OECreateIsoSmiString(target);
    
    unintendedProducts = 0;
    unintendedStereoisomers = 0;
    
    productCanSmiSet = set();   # Non-stereoisomer specific canonical, not isomer SMILES
    productIsmSmiSet = set();
    for product in productList:
        productCanSmi = OECreateCanSmiString(product);
        productIsmSmi = OECreateIsoSmiString(product);
        
        if targetCanSmi != productCanSmi and productCanSmi not in productCanSmiSet:
            # Product mismatch that was not already counted
            unintendedProducts += 1;

        if productCanSmi == targetCanSmi and productIsmSmi != targetIsmSmi:
            # Product match, but stereoisomer mismatch
            unintendedStereoisomers += 1;
        productCanSmiSet.add( productCanSmi );
        productIsmSmiSet.add( productIsmSmi );
        
    return (len(productCanSmiSet), len(productIsmSmiSet), unintendedProducts, unintendedStereoisomers);
        
def precursorStatistics(precursorList):
    """Collect some statistics regarding properties of the proposed precursors 
    (or any collection of molecules).  These would be the kind of statistics of interest
    for suggesting favorable retrosynthetic pathways.
    Returned as a dictionary of named value pairs.
    """
    statDict = dict();
    statDict["maxCarbonCount"] = 0;
    statDict["maxAtomCount"] = 0;
    statDict["maxRingSystemCount"] = 0;
    statDict["chiralAtoms"] = 0;
    statDict["chiralBonds"] = 0;
    
    for precursor in precursorList:
        (ringSystemCount, atomRingMap) = OEDetermineRingSystems(precursor);
        statDict["maxRingSystemCount"] = max(ringSystemCount, statDict["maxRingSystemCount"]);
    
        heavyAtoms = 0;
        carbonCount = 0;
        for atom in precursor.GetAtoms():
            if atom.HasStereoSpecified():
                statDict["chiralAtoms"] += 1;
            if atom.IsCarbon():
                carbonCount += 1;
            if not atom.IsHydrogen():
                heavyAtoms += 1;

        for bond in precursor.GetBonds():
            if bond.HasStereoSpecified():
                statDict["chiralBonds"] += 1;

        statDict["maxAtomCount"] = max(heavyAtoms, statDict["maxAtomCount"]);
        statDict["maxCarbonCount"] = max(carbonCount, statDict["maxCarbonCount"]);
    return statDict;

if __name__=="__main__":
    main(sys.argv)
