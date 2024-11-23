"""
    This class maybe used to generate a series of reactions.
    Using a list of reagents and a list of reactants it will return all products that were created.
 
    Currently this class only uses one reagent per reaction.
    
    1 Cycle is one rotation through all unique reactant with reagent pairs.
    Previous product[s] are the new reactant[s] in the next cycle.
    Reagents are used until maximum number of cycles are met, reusing the reagents if necessary.
    
    2 reactants and 2 reagents will result in 6 reactions after one cycle.
    Each individual reactant , and unique reactant pair, is used with each reagent.
    
    Not sure where this file should be, probably will be moved to another location when the full purpose of this class is realized (eg. Who will use it and why).
"""
import string

from openeye.oechem import OEMolecularFormula
from chemutils.Common.Util import molBySmiles
from chemutils.CombiCDB.ReactionModel import SMIRKSReagent;
from chemutils.CombiCDB.ReactionModel import SupplementalDataModel, reactionStepStr
from chemutils.CombiCDB.ReagentManager import ReagentManager, ReagentSearchModel;
from chemutils.CombiCDB.test.ReagentIdentifiers import *
from chemutils.Common.Util import createStandardSmiString, standardizeSmiles

class CycleReactions():
    def __init__( self, reactants=["CC=C"], reagentIDs=[1], maxCycles=1 ):
        """Initialize CycleReactions class variables.
        """
        #Grab reagents based on IDs provided 
        reagentManager = ReagentManager();
        reagentQuery = ReagentSearchModel()
        reagentQuery.reagentIDs = reagentIDs
        self.reagents = reagentManager.loadDBInstances( reagentQuery )

        #Set cycles
        self.currentCycle = 0
        self.maxCycles = maxCycles
        
        #These initial variables are only saved so that they maybe used in __str__()
        self.reagentIDs = reagentIDs
        self.reactants = reactants

        #Construct initial data containers
        self.oldReactants = set([])
        self.productDetails = {}
        for reactant in reactants:
            self.productDetails.setdefault( reactant, None )
        
    def __call__( self ):
        """This will call the main recursive function.
        """
        return self.react( set(self.productDetails.keys()) );
    
    def __str__( self ):
        """This function is used to display useful internal info about an instance of this class.
           Usually Python only gives the name and address when trying to print an instance.
           Hope this helps for debugging.
        """
        #Format reactants and products lists into a strings
        reactants = string.join( self.reactants, ', ' )
        products = string.join( self.productDetails.keys(), ', ' )
        
        #Format reagentIDs list into a string
        reagentIDs = ''
        for ID in self.reagentIDs:
            reagentIDs += ', ' + str(ID)
        reagentIDs = reagentIDs[2:]        
            
        objInfoStr =  """
ERROR:: Debugging? You tried to print, or str(), a CycleReactions object
Object Details::
    Class: CycleReactions
    File: chemutils/CombiCDB/CycleReactions.py
    
    Function( parameters )
        init( reactants=["CC=C"], reagentIDs=[1], maxCycles=1 )
        __call__()
        products(format="smiles")
        productDetails()
    
    Variable: value
        init Reactants: %s
        Reagent IDs: %s
        Current Cycle: %d
        Max Cycles: %d
        Current Products: %s"""% (reactants, reagentIDs, self.currentCycle, self.maxCycles, products) 
        
        return objInfoStr
    
    def productsAsSmiles(self):
        """Return all current unique products as SMILES strings.
           Order of occurrence is not preserved.
        """
        #Return keys, this list of products are already in smiles format
        return self.productDetails.keys()
    
    def productsAsFormulae(self):
        """Return all current unique products as molecular formulae.
            Order of occurrence is not preserved.
        """
        #Change from smiles format to molecular formula format
        formulas = []
        for smiles in self.productDetails.keys():
            mol = molBySmiles(smiles)
            formulas.append( OEMolecularFormula(mol) )
        return formulas
    
    def productDetails(self):
        """Return all current product data:
            productDetails[product] = (reactants, reagentID) }
        """
        return self.productDetails
    
    def react( self, reactants):
        self.currentCycle += 1
        #print "Starting cycle:", self.currentCycle, "out of", self.maxCycles
        #print "New Reactant[s]:", reactants
        #print "Total Product[s]:", self.productDetails.keys()
        
        #Find the the difference between these so that we can prepend them in the proper order
        #Prepending is done within uniquePairs()
        oldProducts = set(self.productDetails.keys()) - reactants 
        
        #Generate a smiles list of all unique combinations, including singular and pair combos
        reactantPairs = self.uniquePairs( list(reactants), list(oldProducts) )
        
        #Clear reactants, but keep 1 and 2 cycles back as oldReactants
        #This way we can compare with new reactants in order to break cyclic reactions
        if self.currentCycle % 3 == 0:
            self.oldReactants = set([])
        self.oldReactants.union( reactants )
        reactants = set([])
        
        # Go through each pair of reactants
        for reactantPair in reactantPairs:                                               
            #Print out which reactants where used
            #print "Choose:"
            #print reactantPair
                
            #Loop through each reagent and react with the pair of reactants
            productsList = {}
            for reagent in self.reagents:
                # Not interested in inherent products
                reagent.removeInherentProducts = True
                # Optional object to receive extra information
                supplementalData = SupplementalDataModel()
                # Plug reactants into reagent, like a function to generate products
                if len(reactantPair) > 1:
                    #Only plug in two reactants if the reagent can take two
                    if reagent["expected_reactants"] > 1:
                        productList = reagent( [molBySmiles(reactantPair[0]), molBySmiles(reactantPair[1])], supplementalData )
                else:
                    productList = reagent( [molBySmiles(reactantPair[0])], supplementalData )
                
                #print "\tReagent:\n\t\t", reagent["description"]
                #print "\tProduct[s]:"
                
                for products in productList:         
                    #Convert, so saved copy is in smiles format
                    products = createStandardSmiString(products)
                    #Split on the '.' char, so we can make sure that all products are unique
                    #In smiles 2 molecules may be produced, 'C=CC.C=CC', separated by a '.', but we only care about 'C=CC'  
                    products = products.split('.')
                        
                    for product in products:
                        #Save, reactant[s] and reagent (values) that created this product (key) 
                        self.productDetails.setdefault( product, ({'Reactants':reactantPair}, {'ReagentID':reagent["reagent_id"]}) )
                        
                        #Add new products, we need to keep track of which are new
                        reactants.add(product)

                        #Print product
                        #print "\t\t", product
                        #print
        
        #Removing cyclic reactions by calculating only the different reactants from the last cycle, and current cycle 
        reactanats = reactants - self.oldReactants
        
        #Check to see if we have reached the last cycle or not        
        if self.maxCycles > self.currentCycle:
            #Send off new reactants into the next cycle, ie:recursion/iteration
            self.react( reactants )
    
    def uniquePairs(self, reactants, products):
        """This function will create a unique list of combination pair sublists, it is the algorithm that drives the recursive function.
           All past products must be passed and only new reactants for the current cycle should be passed.
           With these conditions met, no duplications will be produced.
           
           EX:
                #All products
                products = ["A","B","C"]
                
                #Only new reactants
                reactants  = ["X","Y","Z"]  
                
                pairs = uniquePairs(products, reactants)
                print pairs
                
                [['X'], ['X', 'Y'], ['X', 'Z'], ['X', 'A'], ['X', 'B'], ['X', 'C'], 
                 ['Y'], ['Y', 'Z'], ['Y', 'A'], ['Y', 'B'], ['Y', 'C'], 
                 ['Z'], ['Z', 'A'], ['Z', 'B'], ['Z', 'C']]
                 
                Note that combination pair "A, B" is not listed above, that is because it is assumed to have been calculated on the previous cycle.
                The next example would be that previous cycle, which coincidentally is first cycle.
                No reactants have been used yet, thus the empty list for products. 
            EX:
                #All products
                products = []
                
                #Only new reactants
                reactants  = ["A","B","C"]  
                
                pairs = uniquePairs(products, reactants)
                print pairs
                
                [['A'], ['A', 'B'], ['A', 'C'], 
                 ['B'], ['B', 'C'], 
                 ['C']]
        """
        pairs=[]
        #Prepend new reactants on to products so that they are in the same order in both lists, reactants and products
        #Order matters in this case, the first element in the list is poped off, and it must be the same in both lists.
        #This will also make this algorithm run faster, it only needs to parse the list once.
        products = reactants + products
        
        for reactant in reactants:
            for product in products:
                if reactant != product:
                    pairs.append([reactant, product])
                else:
                    pairs.append([reactant])
            products.pop(0)
        
        return pairs

"""
##This is a test case, in which you may test the functionality of the CycleReactions class.##
#If you want a simpler example, try reducing reactants/reagents/maxcycles.

#Build instance and run reactions
reactions = CycleReactions(["C[C@@H]1C=C[C@H]1C","CC=C"], [HYDROBROMINATION,PERICYCLIC_THERMAL], 2)
reactions()

#This will output SMILES format
print "\nFinal cycle produced:", reactions.productsAsSmiles(), "\n"

#This will output molecular formula format
print "\nFinal cycle produced:", reactions.productsAsFormulae(), "\n"

#When giving a product, in smiles format, it will return the reactant[s] and reagent that created it
print reactions.productDetails["CC/C=C/[C@H](C)Br"], "\n"

#This last function is more for debugging, or if you just want to know what this object is, you can print it
#Although it may not be useful for automated purposes, the result is a huge string, at least its human readable.
#print reactions
"""
