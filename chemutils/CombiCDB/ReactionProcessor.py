#!/usr/bin/env python
from openeye.oechem import oemolistream, oemolostream, OEWriteMolecule, OEFormat_ISM;
from openeye.oechem import OELibraryGen, OEUniMolecularRxn;
from openeye.oechem import OERxnRole_Product, OEAddMols, OEGraphMol;
from openeye.oechem import OECreateIsoSmiString, OEAssignImplicitHydrogens;
from openeye.oechem import OEBondStereo_CisTrans, OEBondStereo_Trans, OEDetermineRingSystems;
from chemutils.Common.Util import ProgressDots;
from chemutils.Common.Util import molToSmiList, createStandardSmiString;
from chemutils.Common.MolExt import enumerateRacemicMixture;
from chemutils.Common.MolExt import determineRingSystemInfo, isBridgehead;
from chemutils.CombiCDB.Const import COMMENT_TAG, EST_INPUT, REACTION_LABEL, REACTION_LABEL_END, REACTANT_LIST_LABEL;
from chemutils.CombiCDB.Const import SENTINEL_NEUTRAL_CHARGE, SENTINEL_REJECT_CHARGE, SENTINEL_REJECT_IMMEDIATE_CHARGE;
from chemutils.CombiCDB.Util import readIDFile, log;
import sys, os
#from sets import Set;
from optparse import OptionParser
from chemutils.Common.Util import stdOpen, isStdFile, virtual_oemolistream, virtual_oemolostream;
from chemutils.Common import IteratorFactory;

def main(argv):
    """Command-line main method"""
    usageStr =  "usage: %prog [options] <reactantFile> <smirksFile> <productFile>\n"+\
                "   For more information, run the following under Python\n"+\
                "   >>> import ReactionProcessor\n"+\
                "   >>> help(ReactionProcessor)\n"
    parser = OptionParser(usage=usageStr)
    parser.add_option("-r", "--incReactants", action="store_true", help="When generating products, set this flag if you wish to have the reactants listed along side the products, for easy depiction / viewing of the whole reaction.")
    parser.add_option("-s", "--ignoreSelfReactions", action="store_true", help="Don't count products generated from a reactant reacting with another copy of itself");
    parser.add_option("-d", "--DBPrepare", dest="DBPrepare",metavar="<DBInputFile>", help="Generate a results file that is easily inserted into the database.  Interpret the 3 command line arguments as <reactantIDFile> <smirksIDFile> and <productIDFile> instead, where each has one line per respective object that ends with the database ID for the respective molecule / reaction.  Use with -p to specify the productFile already generated by this module, whose molecule labels include the information to piece together which IDs go with which.  Note it is imperative that the ID files' order match the source molecule / SMIRKS files exactly since the labels match them by index position.")
    parser.add_option("-p", "--productFile", dest="productFile",metavar="<productFile>", help="For use with -d option.  The product file created by a normal run of this module.  The label generated for each product molecule includes the information needed by the database preparation script.")
    (options, args) = parser.parse_args(argv[1:])

    instance = ReactionProcessor()
    instance.setIncludeReactants(options.incReactants)
    instance.ignoreSelfReactions = options.ignoreSelfReactions;

    if options.DBPrepare != None and options.productFile != None and len(args) >= 3:
        instance.formatDBFileByFilename( options.productFile, args[0], args[1], args[2], options.DBPrepare )
    elif len(args) >= 3:
        instance.generateProductsByFilename(args[0],args[1],args[2])
    else:
        parser.print_help()
        sys.exit(-1)

class ReactionProcessor:
    """Given a set of SMIRKS reactions and reactant molecules, generates as many 
    combinatorial products as possible by running every reactant permutation 
    through each reaction

    Also includes a script to generate the output in a format easily inserted
    into the application database.  Assuming starting with some reactant
    and SMIRKS files that have NOT been inserted to the database, a complete
    run, including inserting the product info into the database could be
    accomplished with the following from the command line:

    ===========================================================================
    python ReactionProcessor.py reactant.smi example.smirks product.smi
    python DBUtil.py -ireactant.smi     -tMOLECULE -oreactant.smi.id    CAN_SMILES LABEL
    python DBUtil.py -iexample.smirks   -tREACTION -oexample.smirks.id  SMIRKS LABEL
    python DBUtil.py -iproduct.smi      -tMOLECULE -oproduct.smi.id     CAN_SMILES LABEL
    python ReactionProcessor.py -dsynthesis.txt -pproduct.smi reactant.smi.id example.smirks.id product.smi.id
    python DBUtil.py -isynthesis.txt    -tSYNTHESIS -osynthesis.id      PRODUCT_ID  REACTION_ID REACTANT_ID REACTANT_POSITION
    ===========================================================================

    Alternatively, if you wish to use reactants and SMIRKS from the database, something like this:

    ===========================================================================
    python DBUtil.py "select CAN_SMILES, LABEL, MOLECULE_ID from MOLECULE"  reactant.smi
    python DBUtil.py "select SMIRKS, LABEL, REACTION_ID from REACTION"      example.smirks
    python ReactionProcessor.py reactant.smi example.smirks product.smi
    python DBUtil.py -iproduct.smi      -tMOLECULE -oproduct.smi.id     CAN_SMILES LABEL
    python ReactionProcessor.py -dsynthesis -pproduct.smi reactant.smi example.smirks product.smi.id
    python DBUtil.py -isynthesis        -tSYNTHESIS -osynthesis.id      PRODUCT_ID  REACTION_ID REACTANT_ID REACTANT_POSITION
    ===========================================================================

    Input: 
    - Reacant molecule file
        Can be any format understandable by oemolistream, assuming a properly 
        named extension.  For example, "molecules.smi" for SMILES format.

    - SMIRKS reaction file
        File containing one SMIRKS reaction string per line that will 
        be used to process the reactants

    Either of the above can take stdin as their source by specifying the 
    filename "-" or ".smi" or something similar.  See documentation of 
    oemolistream for more information

    Output:
    - Product molecule file
        Outputs all possible products generated from the SMIRKS reactions
        out of the reactant molecules.  Again, redirection to stdout possible 
        by specifying the filename "-".  Each product SMILES will be followed 
        by a molecule "title" of the format "SMIRKS[A]Reactants[X,Y,Z,etc.]" where
            A = Index / position in the SMIRKS reaction file of the reaction 
                used to generate this product.  Index is zero-based
            X,Y,Z,etc. = Index / position in the Reactant molecule file of 
                the respective reactant used
    """

    """Set to ignore products generated by self-reactions.  
    That is, any reactions which take more than one reactant, 
    don't allow the same reactant twice.
    """
    ignoreSelfReactions = True;

    """
    Primarily use if wish to include the reactants that yield each
    product as part a composite product molecule.
    """
    includeReactants = None;

    """If true, any unused reactants from the reactantOEISFactory
    will be reported as additional (unreacted) products.
    """
    includeUnusedReactants = False;

    """Normalization transformations to check for and correct certain error patterns 
    in generated products.  These are especially foudn to occur when the rules break aromatic bonds,
    and they incorrectly assume where double vs. single bonds would be in the Kekule format
    """
    normalizationTransforms = None;

    def __init__(self):
        self.normalizationTransforms = list();
        #self.normalizationTransforms.append(OEUniMolecularRxn("[CX3+0:1]1-[*:2]=[*:3]-[*:4]=[CX3+:5]-[*:6]1>>[CX3:1]1=[*:2]-[*:3]=[*:4]-[CX3+:5]-[*:6]1"));
        #self.normalizationTransforms.append(OEUniMolecularRxn("[CX3+0:1]([*:10])([*:11])-[CX3+0:2]([*:12])([*:13])>>[CX3+0:1]([*:10])([*:11])=[CX3+0:2]([*:12])([*:13])"));
        #self.normalizationTransforms.append(OEUniMolecularRxn("[CX3+1:1]1=[*:2][*:3]=[CX4+0:4][*:5]=[*:6]1>>[C+1:1]1[*:2]=[*:3][C:4][*:5]=[*:6]1"));
        #self.normalizationTransforms.append(OEUniMolecularRxn("[CX3+1:1]1=[*:2][*:3]=[*:4][*:5]=[CX4+0:6]1>>[CX3+1:1]1-[*:2]=[*:3]-[*:4]=[*:5]-[CX4+0:6]1"));

    def getIncludeReactants(self):
        return self.includeReactants;

    def setIncludeReactants(self,value):
        self.includeReactants = value;

    def generateProductsByFilename( self, reactantFilename, smirksFilename, productFilename ):
        """Opens files with respective names and delegates most work to "generateProducts"
        """
        reactantOEISFactory = IteratorFactory.oemolistreamFactory(reactantFilename);
        smirksFile = stdOpen(smirksFilename,"r",sys.stdin)
        productOEOS = oemolostream( productFilename )

        #reactantOEIS = oemolistream()
        #reactantOEIS.openstring(open(reactantFilename,"r").read())
        #reactantOEIS.SetFormat(1)

        try:
            self.generateProducts( reactantOEISFactory, smirksFile, productOEOS )
        finally:
            smirksFile.close()
            productOEOS.close()

    def generateProducts( self, reactantOEISFactory, smirksFile, productOEOS ):
        """Primary method, reads the source files to generate products to 
        the output file.  See module documentation for more information.

        Note:  This method takes actual File objects, oemolistreams and oemolostreams,
        not filenames, to allow the caller to pass "virtual Files" for the purpose 
        of testing and interfacing.  Use the "main" method to have the module take 
        care of opening files from filenames.

        Note that the reactantOEISFactory is not a simple oemolistream either,
        but a factory object that can generate oemolistreams over the list of
        reactants.  This is necessary as nested loops iterating over the
        reactants simultaneously is required.
        """

        smirksList = self.readSMIRKSFile(smirksFile)
        log.debug(smirksList)

        for iSmirks, smirks in enumerate(smirksList):
            libgen = OELibraryGen(smirks)
            # libgen.SetValenceCorrection(True)
            # libgen.SetExplicitHydrogens(True)

            nProducts = self.applyReaction( libgen, reactantOEISFactory, productOEOS, iSmirks )
            log.info("%d products generated from reaction %s",nProducts,smirks)


    def readSMIRKSFile(self, smirksFile):
        """Read the contents of the file as a list of SMIRKS strings.
        Comment lines prefixed with "#" will be ignored.  
        Expects one SMIRKS string per line of the file.  Each SMIRKS string can be followed
            by any title / comment, etc. separated by whitespace.  These will be ignored.
        """
        smirksList = list()
        for line in smirksFile:
            if not line.startswith(COMMENT_TAG):
                chunks = line.split()
                if len(chunks) > 0:
                    smirks = chunks[0]
                    smirksList.append(smirks)
        smirksFile.close()

        return smirksList

    def applyReaction( self, libgen, reactantOEISFactory, productOEOS, reactionIndex=0, currReactantIndexes=None, reactantList=None, rejectProductSmiSet=None ):
        """Recursive function to apply the reaction in libgen to all possible
        permutations of reactants from the reactantOEIS and outputting the results
        to the productOEOS.  Returns the number of products added for this function call.

        libgen = OELibraryGen initialized with a SMIRKS (or other reaction) string
                and any number of reactants upto libgen.NumReactants()
        reactantOEISFactory = IteratorFactory object that can generate
                oemolistreams over the reactant molecules to feed into libgen.
        productOEOS = Oemolostream to write output products from reaction processing to
        reactionIndex = Index indicating what reaction was used in this libgen.
                Just for labelling purposes of output.
        currReactantIndexes = List of indexes of reactants that have already been set
                on the libgen.  Length of list indicates how many have already been set
                (i.e. the current depth of recursion)
                and the actual indexes are again useful for labelling the output.

        If reactions always had 2 reactants, this design would not be necessary.
        A simple doubly nested loop could enumerate all permutations.  However,
        since an arbitrary reaction may have n reactants, an "n-leveled" nested
        loop would be required, which cannot be determined until runtime.  Thus,
        this recursive approach is used instead.
        """
        if reactantList is None:
            reactantList = []; # List of reactants added so far.

        if currReactantIndexes is None: 
            currReactantIndexes = list();
        reactantPosition = len(currReactantIndexes)

    
        if rejectProductSmiSet is None:
            # When generate products, if any produced should be specially rejected,
            #   store them here in case the caller needs to be aware of them.
            rejectProductSmiSet = set();

        productCount = 0
        if reactantPosition < libgen.NumReactants():
            # Core recursion.  Attempt to set each source molecule to the current reactant position.
            # If the molecule is accepted (fits the respective SMARTS pattern), then recursively
            # call this function, but with one less reactant to fill in

            progressDots = None # Track progress at top iteration level
            if reactantPosition == 0: progressDots = ProgressDots(EST_INPUT/10,EST_INPUT/200,"1st position reactants")

            reactantOEISIter = iter(reactantOEISFactory);
            for iReactant, reactant in enumerate(reactantOEISIter):
                if not self.ignoreSelfReactions or iReactant not in currReactantIndexes:
                    reactantFits = libgen.SetStartingMaterial(reactant, reactantPosition, False)
                    #log.debug("%d%d%d",reactantPosition,iReactant,reactantFits)
                    if progressDots != None: progressDots.Update()
                    if reactantFits:
                        # Adjust parameters for next recursion level
                        reactantList.append(reactant)
                        currReactantIndexes.append(iReactant) # Set last index in list to that of current reactant

                        productCount += self.applyReaction( libgen, reactantOEISFactory, productOEOS, reactionIndex, currReactantIndexes, reactantList, rejectProductSmiSet );

                        # Reset parameters to pre-recursion state
                        currReactantIndexes.pop() # Restore index list back to original state
                        reactantList.pop()
                    libgen.ClearStartingMaterial(reactantPosition)
        else:
            # Base case.  No more reactants to add, just generate product output
            # See module documentation for description of label format and contents
            productLabel = "%s%d%s%s%s" % (REACTION_LABEL, reactionIndex, REACTION_LABEL_END, REACTANT_LIST_LABEL, str(currReactantIndexes).replace(" ",""))

            for product in libgen.GetProducts():
                if self.getIncludeReactants():
                    self.addReactants(product, libgen, reactantList)
                product.SetTitle(productLabel)
                validProduct = self.productPostProcessing(product, reactantOEISFactory, reactantList);
                if validProduct:
                    OEWriteMolecule(productOEOS, product)
                    productCount += 1
                else:
                    # Invalid product for some reason.  
                    # Store in special variable so caller can be aware of rejecting these
                    rejectSmi = createStandardSmiString(product);   # Don't just use OECreateIsoSmiString,
                    rejectProductSmiSet.add( rejectSmi );           #   sometimes that has parsing errors on modified molecules
                    
        return productCount

    def applyReactionBySmirks( self, smirks, reactantList, uniqueOnly=True ):
        """Convenience method.  Parse out the SMIRKS string for the caller and
        collect product results in a list, rather than requiring an OE output stream.
        This instantiates separate copies of every product list, thus being less
        efficient in memory usage than a streaming process.  Should only be used
        for convenience.
        
        Reactant list parameter is expected to be a list of molecule objects.
        
        If uniqueOnly, will only return non-redundant results
        """
        libgen = OELibraryGen(smirks);
        productOEOS = virtual_oemolostream(OEFormat_ISM);
        nProducts = self.applyReaction( libgen, reactantList, productOEOS );
        
        productList = [];
        productSmiSet = set();
        productOEIS = virtual_oemolistream( productOEOS.GetString() );
        for product in productOEIS.GetOEGraphMols():
            productSmi = OECreateIsoSmiString(product);
            if productSmi not in productSmiSet:
                productCopy = OEGraphMol(product);
                productList.append( productCopy );
                productSmiSet.add(productSmi);
        return productList;
    
    def addReactants(self,product,libgen, reactantList):
        """Given a reaction product molecule (OEMolBase) and the library generator
        that created it (OELibraryGen), find all of the reactants (starting materials)
        from the libGen and add them as part of the product molecule such that
        the product molecule will instead represent the whole reaction, reactants included.
        
        Note that this method assumes that there is only one starting material per
        reactant position.  For other library generation applications, this assumption
        may be true.
        """
        # First specify the product as a reaction molecule with product atoms
        product.SetRxn(True)
        for atom in product.GetAtoms():
            atom.SetRxnRole(OERxnRole_Product)
        # Now find a reactant to add for each position
        for reactant in reactantList:
            OEAddMols(product,reactant)

    def productPostProcessing(self, product, reactantOEISFactory, reactantList ):
        """Post-processing of product just before it is finally written to the output.
        Return True if everything is okay.  Return False if this product has an error
        and should be rejected from final output
        """
        ringSystemInfo = determineRingSystemInfo(product);

        for bond in product.GetBonds():
            acceptableBond = True;
            
            if bond.GetOrder() == 2 and bond.IsInRing():

                if bond.GetBgn().GetDegree() == 2+1 and bond.GetEnd().GetDegree() == 2+1:
                    if isBridgehead(bond.GetBgn(), ringSystemInfo) or isBridgehead(bond.GetEnd(), ringSystemInfo):
                        # Anti-bredt olefin, unacceptable
                        acceptableBond = False;

                if bond.HasStereoSpecified():
                    # Should not have trans double bonds internal to a ring
                    (ringCount, atomRingMap) = OEDetermineRingSystems(product);
                    bondRingIdx = atomRingMap[ bond.GetBgn().GetIdx() ];
                    # Find bond neighbors part of the same ring
                    neighbors = [];
                    for atom in bond.GetBgn().GetAtoms():
                        if atom.GetIdx() != bond.GetEnd().GetIdx() and atomRingMap[atom.GetIdx()] == bondRingIdx:
                            neighbors.append(atom);
                            break;
                    for atom in bond.GetEnd().GetAtoms():
                        if atom.GetIdx() != bond.GetBgn().GetIdx() and atomRingMap[atom.GetIdx()] == bondRingIdx:
                            neighbors.append(atom);
                            break;
                    if bond.GetStereo(neighbors, OEBondStereo_CisTrans) == OEBondStereo_Trans:
                        # Intra-ring trans pi bond.  Unacceptable
                        acceptableBond = False;
            
            if not acceptableBond:
                # Reject product, but must be careful later.  
                # If just return nothing, subsequent methods may not even realize anything was tried at all
                # Isolate the component that has the bad trans ring bond, deleting any other components
                visitedAtomIndexes = self.__visitBondedAtoms( bond.GetBgn() );
                for atom in product.GetAtoms():
                    if atom.GetIdx() not in visitedAtomIndexes:
                        product.DeleteAtom(atom);
                return False;


        # Check for certain bond error patterns that occur when the rules break aromatic bonds,
        #   and they incorrectly assume where double vs. single bonds would be in the Kekule format
        for i, transform in enumerate(self.normalizationTransforms):
            preC = createStandardSmiString(product)
            changed = transform(product);
            postC = createStandardSmiString(product)
            
            if changed:
                print >> sys.stderr, preC, postC
                pass
                 

        # Correct any sentinel value charges meant to represent a neutral charge
        for atom in product.GetAtoms():
            if atom.GetFormalCharge() == SENTINEL_NEUTRAL_CHARGE:
                atom.SetFormalCharge(0);
        if self.includeUnusedReactants:
            usedReactantSmiList = molToSmiList(reactantList);
            for reactant in reactantOEISFactory:
                reactantSmi = OECreateIsoSmiString(reactant);
                if reactantSmi in usedReactantSmiList:
                    usedReactantSmiList.remove(reactantSmi);
                else:
                    # Unused reactant, add it to the product as an unused "agent / catalyst"
                    if self.getIncludeReactants():
                        # If including reaction side of the equation, first add it there
                        OEAddMols(product, reactant);
                        # Use temp label to distinguish which parts to label as product vs. reactant
                        for atom in product.GetAtoms():
                            atom.SetBoolData("RxnRoleLabeled",True);
                    OEAddMols(product, reactant);
                    if self.getIncludeReactants():
                        # Must label the product side if separating reactant vs. product side
                        for atom in product.GetAtoms():
                            if not atom.HasBoolData("RxnRoleLabeled"):
                                atom.SetRxnRole(OERxnRole_Product);
        return True;


    def formatDBFileByFilename( self, productFilename, reactantIDFilename, smirksIDFilename, productIDFilename, dbFilename ):
        """Opens files with respective names and delegates most work to "formatDBFile"
        """
        productFile     = oemolistream(productFilename)
        reactantIDFile  = stdOpen(reactantIDFilename,"r",sys.stdin)
        smirksIDFile    = stdOpen(smirksIDFilename,"r",sys.stdin)
        productIDFile   = stdOpen(productIDFilename,"r",sys.stdin)
        dbFile          = stdOpen(dbFilename,"w",sys.stdout)

        try:
            self.formatDBFile( productFile, reactantIDFile, smirksIDFile, productIDFile, dbFile )
        finally:
            productFile.close()
            reactantIDFile.close()
            smirksIDFile.close()
            productIDFile.close()
            dbFile.close()


    def formatDBFile( self, productOEIS, reactantIDFile, smirksIDFile, productIDFile, dbFile ):
        """Given the database IDs of reactants, reactions (smirks) and products and 
        information indicating how they are all related, generate a simple text file
        that should be very easy to import into the database to persist that
        association information.

        Each line should correspond to a row in the SYNTHESIS table, with values
        to insert respective to PRODUCT_ID, REACTION_ID, REACTANT_ID and REACTANT_POSITION
        """
        reactantIDList  = readIDFile( reactantIDFile )
        smirksIDList    = readIDFile( smirksIDFile )
        productIDList   = readIDFile( productIDFile )

        dbFile.write("# SYNTHESIS\n")
        dbFile.write("# PRODUCT_ID\tREACTION_ID\tREACTANT_ID\tREACTANT_POSITION\n")

        for iProduct, mol in enumerate(productOEIS.GetOEGraphMols()):
            label = mol.GetTitle().split()[0]   # Take first token as functional label.  Dump any extra items
            # Parse it out, assuming a label structure like Smirks[x]Reactants[a,b,c,etc.]
            iSmirks, iReactantList = label.split(REACTION_LABEL)[1].split(REACTION_LABEL_END+REACTANT_LIST_LABEL)  # Remove labels
            iSmirks = int(iSmirks)
            iReactantList = iReactantList[1:].split("]")[0].split(",")  # Clip off brackets and split by commas
            for i in range(len(iReactantList)):
                iReactantList[i] = int(iReactantList[i])

            # Values all parsed, now just right out the IDs upon lookup
            for reactantPosition, iReactant in enumerate(iReactantList):
                dbFile.write("%s\t%s\t%s\t%d\n" % (productIDList[iProduct],smirksIDList[iSmirks],reactantIDList[iReactant],reactantPosition) )

    def __visitBondedAtoms( self, atom, visitedAtomIndexes=None ):
        """Starting from the given atom, add every visited atom's index to the 
        visistedAtomIndexes Set.  Recursively visit all bonded atoms
        """
        if visitedAtomIndexes is None:
            visitedAtomIndexes = set();
        visitedAtomIndexes.add( atom.GetIdx() );
        for neighbor in atom.GetAtoms():    # Recurse on neighbors
            if neighbor.GetIdx() not in visitedAtomIndexes:
                self.__visitBondedAtoms( neighbor, visitedAtomIndexes );
        return visitedAtomIndexes;

if __name__=="__main__":
    main(sys.argv)
    