"""Models for electron flow arrows representing reaction mechanism diagrams.
Key functions that external callers will probably be interested in is

    generateMechanismPairs:
        to take apply mechanism arrow information from a reaction_profile object to a set of reactant molecules.
    moveAtomElectrons:
    moveBondElectrons:
        Use to rearrange electron configuration around molecules.
        Should map directly to arrows in a mechanism diagram.
"""
import sys
#from sets import Set;

from openeye.oechem import OEGraphMol, OECreateCanSmiString
from openeye.oechem import OERxnRole_Product

from chemutils.Common.MolStdValue import atomElectronegativity
from chemutils.Common.Util import createStandardSmiString, molBySmiles, standardizeSmiles
from chemutils.Common.Util import createAtomMapSmiString
from chemutils.Common.OrbitalModel import orbitalInfo
from chemutils.Common.MolExt import clearAtomMaps, clearMolStereo
from chemutils.CombiCDB.ReactionProcessor import ReactionProcessor
from chemutils.CombiCDB.Const import NORMAL_CHARGE_THRESHOLD, SENTINEL_LABEL_CHARGE, SENTINEL_RADICAL_CHARGE, SENTINEL_CHARGE, SENTINEL_OPTIONAL_STEP_CHARGE, SENTINEL_REJECT_IMMEDIATE_CHARGE, SENTINEL_REJECT_CHARGE
from chemutils.CombiCDB.Const import OK, CAUTION, WARNING, ERROR, DISFAVORED
from chemutils.CombiCDB.Util import log

"""Mechanism arrow code formatting characters"""
ATOM_DELIM = ","
ARROW_DELIM = ";"
SINGLE_ARROW = "-"
DOUBLE_ARROW = "="

class ElectronArrow:
    """Class / struct containing electron flow information for reaction mechanisms.
    Should map directly to a single curved arrow for an arrow-pushing reaction mechanism.
    """
    mol = None;
    sourceIndexes = None;
    sourceAtoms = None;
    nElectrons = 0;
    targetIndexes = None;
    targetAtoms = None;
    
    def __init__(self, mol, sourceIndexes, nElectrons, targetIndexes):
        """Initialization constructor given the (composite) molecule object
        the arrows would act upon.  The source and target indexes correspond
        to atom mapping indexes to indicate where the arrow comes from and goes to.
        """
        # Ensure the source and target indexes are consistently represented tuples
        sourceIndexes = list(sourceIndexes);
        sourceIndexes.sort();
        sourceIndexes = tuple(sourceIndexes);

        targetIndexes = list(targetIndexes);
        targetIndexes.sort();
        targetIndexes = tuple(targetIndexes);
        
                
        self.mol = mol;
        self.sourceIndexes = sourceIndexes;
        self.targetIndexes = targetIndexes;
        self.sourceAtoms = [];
        self.targetAtoms = [];
        
        self.__rebuildAtomLists()
        
        self.nElectrons = nElectrons;
    
    def __rebuildAtomLists(self):
        """Simple function to rebuild the atom list if the indices have changed"""
        # Track molecule atoms by map index for rapid retrieval
        atomByMapIdx = dict();
        for atom in self.mol.GetAtoms():
            if atom.GetMapIdx() >= 0:
                atomByMapIdx[atom.GetMapIdx()] = atom;
        
        self.sourceAtoms = [];
        self.targetAtoms = [];
        
        for sourceIndex in self.sourceIndexes:           
            sourceAtom = atomByMapIdx[sourceIndex];
            self.sourceAtoms.append( sourceAtom );    
        
        for targetIndex in self.targetIndexes:           
            targetAtom = atomByMapIdx[targetIndex];
            self.targetAtoms.append( targetAtom );
        
    def isBondDissociationArrow(self):
        """Simple method to see if the movement matches x,y->y pattern"""
        return len(self.sourceIndexes) == 2 and len(self.targetIndexes) == 1 and\
                self.targetIndexes[0] in self.sourceIndexes;
    
    def isBondPivotArrow(self):
        """Simple method to see if arrow matches x,y->y,z pattern"""
        return len(self.sourceIndexes) == 2 and len(self.targetIndexes) == 2 and \
                len(set(self.sourceIndexes).intersection(set(self.targetIndexes))) == 1;
    
    def isBondFormArrow(self):
        """Simple method to check if matches x->y"""
        return len(self.sourceIndexes) == 1 and len(self.targetIndexes) == 1;
    
    
    def __str__(self):
        """Format into arrow code"""
        sourceStr = [];
        for sourceIndex in self.sourceIndexes:
            sourceStr.append(str(sourceIndex));
        sourceStr = str.join(ATOM_DELIM, sourceStr );

        electronCode = "";
        if self.nElectrons == 1:
            electronCode = SINGLE_ARROW;
        else:
            electronCode = DOUBLE_ARROW;

        targetStr = [];
        for targetIndex in self.targetIndexes:
            targetStr.append(str(targetIndex));
        targetStr = str.join(ATOM_DELIM, targetStr );

        return "%s%s%s" % (sourceStr, electronCode, targetStr);

    def apply( self, mol, halfBondIndexes=None ):
        """Apply the electron arrow pushing to the respective molecule to actually
        "move" the electrons indicated, producing an (intermediate) product.
        """
        log.debug("%s %s %s" % (createStandardSmiString(mol), str(self), halfBondIndexes) );
        
        if len(self.sourceAtoms) > 1:
            # Source must be a bond
            bond = mol.GetBond( self.sourceAtoms[0], self.sourceAtoms[1] );

            pivotAtom = None;
            farTargetAtom = None;
            for targetAtom in self.targetAtoms:
                targetIndex = targetAtom.GetMapIdx();
                if targetIndex in self.sourceIndexes:
                    # When a source atom is used as a target, means it's the "pivot" to remain connected to
                    pivotAtom = targetAtom;
                else:
                    farTargetAtom = targetAtom;

            moveBondElectrons( bond, farTargetAtom, self.nElectrons, pivotAtom, halfBondIndexes=halfBondIndexes );
        else:
            # Source must be an atom
            sourceAtom = self.sourceAtoms[0];

            # Should only be 1 target. but just in case user includes the source atom as part of a "forming bond" target, ignore it.
            farTargetAtom = None;
            for targetAtom in self.targetAtoms:
                targetIndex = targetAtom.GetMapIdx();
                if targetIndex not in self.sourceIndexes:
                    farTargetAtom = targetAtom;

            # Going to be a problem for radical reactions.  
            #   Cannot create a new bond until have 2nd radical meeting the 1st
            moveAtomElectrons( sourceAtom, farTargetAtom, self.nElectrons, halfBondIndexes=halfBondIndexes );
    
    def applyAll( mol, arrowObjList ):
        """Apply all of the electron arrow objects in the list to the respective
        (composite) molecule to actually "move" the electrons and produce a product.

        For now just process each arrow object sequentially.  Works for most cases,
        but will cause problems with radical reactions where often 2 simultaneous
        movements / arrows are needed to complete or break a bond.  Will eventually
        need some more wholistic approach that keeps some tracking variables to deal
        with such cases.
        """
        halfBondIndexes = set();
        for arrowObj in arrowObjList:
            arrowObj.apply( mol, halfBondIndexes );
    applyAll = staticmethod(applyAll);

    def parseArrowCode( arrowCode ):
        """Convenience function to parse out arrow code strings.  Syntax described elsewhere, including
        http://baldigwiki.ics.uci.edu/index.php/Web_Development#Details_on_.27arrow_description.27_syntax    
            
        Returns a 3-ple consisting of 
        - The atom map indexes of the source atoms
        - The number of electrons moving with the arrow (1 or 2)
        - The atom map indexes of the target atoms
        """
        nElectrons = 0; # Number of electrons moving according to the arrow

        electronCode = None;
        if SINGLE_ARROW in arrowCode:
            nElectrons = 1;
            electronCode = SINGLE_ARROW;
        elif DOUBLE_ARROW in arrowCode:
            nElectrons = 2;
            electronCode = DOUBLE_ARROW;

        arrowTerminusCodes = arrowCode.split(electronCode);
        sourceIndexes = arrowTerminusCodes[0].split(ATOM_DELIM);
        targetIndexes = arrowTerminusCodes[1].split(ATOM_DELIM);

        # Convert index strings to actual integers
        for i, indexStr in enumerate(sourceIndexes):
            sourceIndexes[i] = int(indexStr);
        for i, indexStr in enumerate(targetIndexes):
            targetIndexes[i] = int(indexStr);

        # Sort the lists (may just be single item if representing atoms anyway, but maybe 2 if representing bond electrons)
        #   and convert to constant tuples for consistent ordering and lookup and subsequent indexing
        sourceIndexes.sort();
        sourceIndexes = tuple(sourceIndexes);
        targetIndexes.sort();
        targetIndexes = tuple(targetIndexes);

        arrowData = (sourceIndexes, nElectrons, targetIndexes);
        return arrowData;
    parseArrowCode = staticmethod(parseArrowCode);

    def makeCanonicalArrowList(arrowList):
        """In order to assist in the arrow->orbital mapping, rewrite the source/target idx lists so common atoms are 
        inside."""
        [arr.makeCanonicalArrow() for arr in arrowList];
        return arrowList;
    makeCanonicalArrowList = staticmethod(makeCanonicalArrowList)
    
    def makeCanonicalArrow(self):
        """Rewrite source/target idx lists so common atoms are inside."""
        if len(self.sourceIndexes) == 2:
            commonElem = list(set(self.sourceIndexes).intersection(set(self.targetIndexes)));
            if len(commonElem) == 1:
                commonElem = commonElem[0];
                if commonElem in self.sourceIndexes:
                    self.sourceIndexes = tuple([self.sourceIndexes[1-list(self.sourceIndexes).index(commonElem)], commonElem]);
                if len(self.targetIndexes) == 2 and commonElem in self.targetIndexes:
                    self.targetIndexes = tuple([commonElem, self.targetIndexes[1-list(self.targetIndexes).index(commonElem)]]);
        elif len(self.sourceIndexes) == 1 and len(self.targetIndexes) == 2:
            # Source is a single atom but target is multiple atoms.
            #   This doesn't make sense, unless the source is one of the targets,
            #   but that causes redundant confusion, so remove it.
            reducedTargetSet = set(self.targetIndexes) - set(self.sourceIndexes);
            self.targetIndexes = tuple(reducedTargetSet);
        
        self.__rebuildAtomLists();
        
        return;
    

    def parseArrowCodes( mol, arrowCodes ):
        """Further convenience function.  Given a composite string representing multiple
        electron flow arrows and the (composite) molecule they are intended to act upon,
        parse out all of them and return a list of respective ElectronArrow object instances.
        """
        if arrowCodes.strip() == "":
            # Beware of blank / empty code
            return [];
            
        arrowCodeList = arrowCodes.split(ARROW_DELIM);
        arrowObjList = [];
        for arrowCode in arrowCodeList:
            arrowData = ElectronArrow.parseArrowCode(arrowCode.strip());
            (sourceIndexes, nElectrons, targetIndexes) = arrowData;
            arrowObjList.append( ElectronArrow(mol, sourceIndexes, nElectrons, targetIndexes) );
        return arrowObjList;
    parseArrowCodes = staticmethod( parseArrowCodes );

    def prepareArrowCodes( arrowList ):
        """Given a list of ElectronArrow objects, return an arrow code string representation
        of the combined list.  However, first look for opportunities to reduce or simplify
        the list, such as removing redundancies and combining arrows.
        """
        # Eliminate redundancies by keying through a dictionary
        arrowsByCode = dict();
        # Track by sources and targets to look for opportunities to combine
        arrowListsBySourceIndexes = dict();
        for arrow in arrowList:
            arrowCode = str(arrow);
            arrowsByCode[arrowCode] = arrow;
            if arrow.sourceIndexes not in arrowListsBySourceIndexes:
                arrowListsBySourceIndexes[arrow.sourceIndexes] = list();
            arrowListsBySourceIndexes[arrow.sourceIndexes].append(arrow);

        # Look for arrows that could be combinined with a subsequent one.
        # In particular, bond electrons ending on one of the bond's own atoms,
        #   which could then be followed by the electrons on the atom moving on to another target
        comboArrowList = [];
        removeArrowCodes = set();
        for arrowCode, arrow in arrowsByCode.items():
            comboArrow = None;
            if len(arrow.sourceIndexes) == 2 and len(arrow.targetIndexes) == 1 and arrow.targetIndexes[0] in arrow.sourceIndexes:
                # Possibility for combination, look for a target match
                pivotAtomIndex = arrow.targetIndexes[0];
                if arrow.targetIndexes in arrowListsBySourceIndexes:
                    possibleLinkArrows = arrowListsBySourceIndexes[arrow.targetIndexes];
                    for linkArrow in possibleLinkArrows:
                        if arrow.nElectrons == linkArrow.nElectrons:
                            # Looks like a complete match, but make sure this arrow has not already been used
                            linkArrowCode = str(linkArrow);
                            if linkArrowCode not in removeArrowCodes:
                                # Make a record to add this combo
                                sourceIndexes = arrow.sourceIndexes;
                                targetIndexes = set(linkArrow.targetIndexes);
                                # Pivot atom must be in the target list to ensure the single combo arrow
                                #   unambiguously represents the proper electron transfer
                                targetIndexes.add(pivotAtomIndex);  
                                comboArrow = ElectronArrow( arrow.mol, sourceIndexes, arrow.nElectrons, targetIndexes);

                                # Record which arrows can be removed in favor of the combo arrow
                                removeArrowCodes.add( arrowCode );
                                removeArrowCodes.add( linkArrowCode );

                                # Don't bother looking for anymore
                                break;
            
            if comboArrow is not None:
                # Do not need to keep looking for matches for this arrow
                comboArrowList.append( comboArrow );
                continue;
        
        # Replace the component arrows with new combo arrows
        for arrowCode in removeArrowCodes:
            del arrowsByCode[arrowCode];
        for comboArrow in comboArrowList:
            arrowCode = str(comboArrow);
            arrowsByCode[arrowCode] = comboArrow;
        
        # Finally combine all of the code representations of the reduced arrow list into a single code string
        arrowCodeList = arrowsByCode.keys();
        arrowCodeList = sorted(arrowCodeList)
        #arrowCodeList.sort();   # Try to achieve some consistent order
        totalArrowCodeStr = str.join( ARROW_DELIM, arrowCodeList );
        return totalArrowCodeStr;
        
    prepareArrowCodes = staticmethod( prepareArrowCodes );

def prepLabeledReactant( labeledReactant ):
    """Assume reactant came out of a mechanism labeling run, in which case many atoms
    should have bogus formal charges in the 1000s range that should be converted to atom map indexes.
    Any hydrogens participating in the mechanism must be assigned an artificial isotope number so
    that generation of an isomeric SMILES string will not just ignore them as implicit hydrogens.
    """
    # Keep track of electron orbital statistics
    orbitalInfoByAtomIdx = dict();  
    
    for atom in labeledReactant.GetAtoms():
        # Translate extreme charge (in increments of SENTINEL_LABEL_CHARGE (1000)) to atom mapping index
        atomMapIndex = 0;
        while atom.GetFormalCharge() > NORMAL_CHARGE_THRESHOLD:
            # Clear out any label charges that may cause confusion
            if atom.GetFormalCharge() in (SENTINEL_RADICAL_CHARGE, SENTINEL_CHARGE, SENTINEL_OPTIONAL_STEP_CHARGE, SENTINEL_REJECT_IMMEDIATE_CHARGE, SENTINEL_REJECT_CHARGE):
                atom.SetFormalCharge(0);
            else:
                atom.SetFormalCharge( atom.GetFormalCharge()-SENTINEL_LABEL_CHARGE );
                atomMapIndex += 1;
            
        atom.SetMapIdx(atomMapIndex);

        if atom.GetAtomicNum() == 1 and atomMapIndex > 0:   
            # Hydrogen atom involved in mechanism.  Artificially set isotope, 
            #   otherwise will be lost as "implicit hydrogen" upon canonization
            atom.SetIsotope(1);

        try:
            orbitalInfoByAtomIdx[atom.GetIdx()] = orbitalInfo( atom );
        except:
            # Just in case don't know how to calculate orbital information for some atoms (transition metals)
            orbitalInfoByAtomIdx[atom.GetIdx()] = None;

    # Do an extra pass to see if any atom valences were corrupted, usually a result of
    #   aromatic atoms and bonds getting mock charges assigned, which confuses the bond order assignment
    for atom in labeledReactant.GetAtoms():
        orbInfo = orbitalInfoByAtomIdx[atom.GetIdx()];
        if orbInfo is not None and orbInfo["nRadicals"] > 0:
            # See if any adjacent atoms have have a radical as well.  
            #   Highly unlikely, suggests bond order corruption
            for neighbor in atom.GetAtoms():
                neighborOrbInfo = orbitalInfoByAtomIdx[neighbor.GetIdx()];
                if neighborOrbInfo is not None and neighborOrbInfo["nRadicals"] > 0:
                    # Neighbor is also a radical.  Unlikely.  Terminate the radicals by increasing bond order
                    bond = atom.GetBond(neighbor);
                    bond.SetOrder( bond.GetOrder()+1 );
                    # Recalculate orbital info based on bond revision
                    orbitalInfoByAtomIdx[neighbor.GetIdx()] = orbitalInfo( neighbor );
                    orbitalInfoByAtomIdx[atom.GetIdx()] = orbitalInfo( atom );
                    # Don't bother checking anymore neighbors
                    break;

def createNonIsoSmi( mol ):
    """Generate a canonical SMILES representation of the molecule
    which does not include isomer (stereochemistry) information,
    but does retain atom-mapping needed for mechanism diagrams.
    
    Many of the mechanism steps done so far do not support stereochemistry,
    so better not to display stereochemistry information erroneously.
    """
    mol = OEGraphMol(mol);  # Create a copy we're free to modify
    
    # Temporary hack.  Can only support canonical, not isomer SMILES strings for now,
    #   but then canonical SMILES misses hydrogens with atom mapping, even if they have
    #   an artificial isotope set.  Hack 'em as mock placeholder atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0 and atom.GetFormalCharge() == 0:
            atom.SetAtomicNum(0);
    smi = OECreateCanSmiString(mol);   # Stereochemistry not yet supported for arrow diagrams
    smi = standardizeSmiles(smi);   # Strange step.  Some cases with labeled hydrogens and "fake" stereocenters, canonical smi algorithm messes up
    smi = smi.replace("[R","[H:");
    return smi;

def prepLabeledSmi( labeledSmi ):
    """Back hack alterations made by the prepLabeledReactant function.  In particular,
    those added an artificial isotope label to hydrogen atoms to force their appearance
    in the generated SMILES string.  Once the SMILES string is ready, no need to keep
    the artificial label.
    """
    return labeledSmi.replace("[1H","[H"); # Get rid of the artificial 1H isotopes at this point

def reactionProfileSupportsMechanism( reactionProfile ):
    return \
        reactionProfile["mechanism_arrow_codes"] is not None and \
        reactionProfile["mechanism_arrow_codes"] != "" and \
        reactionProfile["warning_level_id"] not in (ERROR, DISFAVORED);

def generateMechanismPairs( reactants, reagent, reactionProfile, reactionProcessor=None, clearLabels=False ):
    """Return pairs (2-ples) of molecules representing reactants and products
    after the mechanism arrows of the reactionProfile are applied to the given reactants.
    
    reactionProcessor: Caller can supply their own reactionProcessor with modified attributes,
        otherise a default one will be used
    clearLabels: If set, returned molecules will not have atom map and isotope labels
        which would normally be used to communicate to image renderer on where to draw arrows
    """
    if reactionProcessor is None:
        reactionProcessor = ReactionProcessor();
    
    extReactants = list( reactants );   # Make list copy
    extReactants.extend( reagent.inherentReactants );
    label_smirks = reactionProfile["mechanism_label_smirks"];
    arrow_codes = reactionProfile["mechanism_arrow_codes"]

    resultPairs = [];

    labeledReactantList = reactionProcessor.applyReactionBySmirks( label_smirks, extReactants );
    for labeledReactant in labeledReactantList:
        # Process any mock charge labels and other issues
        prepLabeledReactant( labeledReactant );

        # Use the arrow information to "complete" the mechanism step, 
        #   telling us what (intermediate) product would be produced
        arrowObjList = ElectronArrow.parseArrowCodes( labeledReactant, arrow_codes );

        # Create a copy of the reactant molecule to get a cleared out version
        reactantCopy = OEGraphMol(labeledReactant);
        if clearLabels:
            clearMechanismLabels(reactantCopy);

        # Actually "move" the electrons to generate a product copy molecule
        labeledProduct = labeledReactant; #OEGraphMol(labeledReactant);

        # Should now have all of the arrow information in arrowObjList
        # Apply these to the labeled molecule to produce the target
        ElectronArrow.applyAll( labeledProduct, arrowObjList );

        # Clear out the product molecule in the same way
        productCopy = OEGraphMol(labeledProduct);
        if clearLabels:
            clearMechanismLabels(productCopy);
    
        resultPairs.append( (reactantCopy, productCopy) );
    
    return resultPairs;

def clearMechanismLabels(mol):
    """Clear any labels used for mechanism specifications.
    In particular, clear (set to 0) all isotope and atom map indexes.
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0);
        atom.SetMapIdx(0);

def clearProductMechanismLabels(mol):
    """Clear mechanism labels, but only for molecules on the product
    side of a reaction equation.
    """
    for atom in mol.GetAtoms():
        if atom.GetRxnRole() == OERxnRole_Product:
            atom.SetIsotope(0);
            atom.SetMapIdx(0);

def equivalentMechanismMol(mol1, mol2):
    """Determine if the two molecules are equivalent for the purposes of mechanism
    diagrams.  This means that stereochemistry and atom map labels should be ignored.
    """
    clearMechanismLabels(mol1);
    clearMechanismLabels(mol2);
    
    mechSmi1 = createNonIsoSmi(mol1);
    mechSmi2 = createNonIsoSmi(mol2);
    
    return (mechSmi1 == mechSmi2);

def moveAtomElectrons( sourceAtom, targetAtom, electronCount=2, arrowObjList=None, halfBondIndexes=None ):
    """Move electrons from sourceAtom to targetAtom.  If these are bonded already,
    just change the bond order and adjust the formal charges.
    If they are not already bonded, then create a new bond and adjust
    the formal charges.
    
    If moving single electrons (radical reactions) better provide halfBonds list.
    Used to keep track of which bonds have only been "half-formed" by a radical movement,
    and are waiting for a complementary one to complete or eliminate the bond.
    Note that the bond could be "1 1/2" or "2 1/2."  The point is that it's just
    not a whole number bond yet.
    
    Return reference to new bond object if it was created.
    
    If the arrowObjList is provided, then an item will be added for the
    ElectronArrow object representing the electron movement.
    Make this optional as excessive ElectronArrow object instantiation can be 
    relatively expensive, and is often unused.
    """

    if halfBondIndexes is None:
        halfBondIndexes = set();
    
    mol = sourceAtom.GetParent();
    bond = sourceAtom.GetBond(targetAtom);
    if bond is not None:
        # Atoms already bonded
        if electronCount == 1:
            # Check if it's an existing "half-bond"
            if bond.GetIdx() in halfBondIndexes:
                # Is a half-bond.  Adding this one more electron just completes the bond
                halfBondIndexes.remove(bond.GetIdx());
            else:
                # Was not a half-bond.  Only adding one electron, so can't complete.  Add it to half-bond list
                bond.SetOrder( bond.GetOrder()+1 );
                halfBondIndexes.add( bond.GetIdx() );
        else:
            # Typical case, moving 2 electrons
            bond.SetOrder( bond.GetOrder()+1 );
    else:
        # Atoms not yet bonded, create a new one
        if sourceAtom.GetParent() != targetAtom.GetParent():
            # Beware that the atoms must exist in the same composite molecule graph object
            #   to avoid object reference inconsistencies.
            # Should really be using the "is" and "not is" comparators for molecule object
            #   reference comparison, but this does not seem to work right for OEGraphMols
            raise Exception("Attempting to move electrons between atoms from separate molecule graph objects");
        
        bond = mol.NewBond( sourceAtom, targetAtom, 1 ); 
        if electronCount == 1:
            halfBondIndexes.add( bond.GetIdx() );
    if electronCount == 2:
        sourceAtom.SetFormalCharge( sourceAtom.GetFormalCharge()+1 );
        targetAtom.SetFormalCharge( targetAtom.GetFormalCharge()-1 );

    if arrowObjList is not None:
        sourceIndexes = [sourceAtom.GetMapIdx()];
        targetIndexes = [targetAtom.GetMapIdx()];
        arrow = ElectronArrow( mol, sourceIndexes, electronCount, targetIndexes );
        
        
        if electronCount == 1:
            # CASE 1:
            # Have to do a fix similar to that done in moveBondElectrons below.
            # If freeRadical, then look for something with a target of the curr source/target
            # and with a source that includes the curr target.
            # If we find this, correct the current one to have a target of curr source/target.
            haveFoundMatch = False;
            matchTarget = [sourceIndexes[0], targetIndexes[0]]
            matchTarget.sort();
            matchTarget = tuple(matchTarget)
            matchSource = targetIndexes[0]
            for arrowObj in arrowObjList:
                if len(arrowObj.sourceIndexes) == 2 and matchSource in arrowObj.sourceIndexes and matchTarget == arrowObj.targetIndexes:
                    haveFoundMatch = True;
                    break;
            
            if haveFoundMatch:
                arrow = ElectronArrow(mol, sourceIndexes, electronCount, [sourceIndexes[0], targetIndexes[0]])
            else:
                # CASE 2:
                # Have a x -> y, and could have a pre-existing y -> x, should make the target of both 
                # to be (x,y)
                idxToFix = -1;
                for i, arrowObj in enumerate(arrowObjList):
                    if len(arrowObj.sourceIndexes) == 1 and len(arrowObj.targetIndexes) == 1 \
                            and arrowObj.sourceIndexes[0] == targetIndexes[0] and \
                            arrowObj.targetIndexes[0] == sourceIndexes[0]:
                        idxToFix = i;
                        break;
                
                if idxToFix != -1:
                    oldArrow = arrowObjList.pop(idxToFix);
                    targetIndexes = [sourceIndexes[0], targetIndexes[0]]
                    arrow = ElectronArrow(mol, sourceIndexes, electronCount, targetIndexes);
                    newOldArr = ElectronArrow(mol, oldArrow.sourceIndexes, electronCount, targetIndexes);
                    
                    arrowObjList.append(newOldArr);
        arrowObjList.append(arrow); 
        
    return bond;

def moveBondElectrons( bond, atom, electronCount=2, pivotAtom=None, arrowObjList=None, halfBondIndexes=None ):
    """Move electrons from the bond to the target atom.
    If the atom is already part of the bond, then just down shift the bond
    order (maybe deleting the bond altogether) and adjust atom formal charges.
    
    Otherwise, down shift this bond order and create a new one to the atom,
    again adjusting atom formal charges.  Note that this only works if the
    "pivotAtom" is specified.  That is, the atom which appears in both
    the source and the target bond.  Otherwise, without the polarity of
    the source bond specified, can't tell which end to remain connected to.
    
    If the arrowObjList is provided, then an item will be added for the
    ElectronArrow object representing the electron movement.
    Make this optional as excessive ElectronArrow object instantiation can be 
    relatively expensive, and is often unused.
    """
    if bond is None:
        raise BondDoesNotExistException(atom, pivotAtom, electronCount);

    if atom is None and pivotAtom is not None:
        # No target specified, but there is a pivotAtom.  The pivot itself must be the target then, just shifting bond electrons to one side
        atom = pivotAtom;
        pivotAtom = None;

    if halfBondIndexes is None:
        halfBondIndexes = set();
    
    mol = bond.GetParent();
    neighborAtom = bond.GetNbr(atom);
    if neighborAtom is not None:
        # Atom already a part of the bond, just shift electrons
        if electronCount == 2:
            atom.SetFormalCharge( atom.GetFormalCharge()-1 );
            neighborAtom.SetFormalCharge( neighborAtom.GetFormalCharge()+1 );
            bond.SetOrder( bond.GetOrder()-1 );
        else:
            # Must be moving only a single electron.  Check if it's a half bond already, then break, otherwise make it a half
            if bond.GetIdx() in halfBondIndexes:
                # Was a half bond, break it down with one more electron out
                bond.SetOrder( bond.GetOrder()-1 );
                halfBondIndexes.remove( bond.GetIdx() );
            else:
                # Was not a half bond, but only moving 1 electron so cannot remove bond entirely
                halfBondIndexes.add( bond.GetIdx() );
        if bond.GetOrder() < 1:
            mol.DeleteBond(bond);
    elif pivotAtom is not None:
        # Just defer to 2-step process
        moveBondElectrons( bond, pivotAtom, electronCount, halfBondIndexes=halfBondIndexes );
        moveAtomElectrons( pivotAtom, atom, electronCount, halfBondIndexes=halfBondIndexes );
    else:
        raise Exception("Ambiguous scenario, don't know which source atom of the bond to leave the new bond attached to.  Call moveBondElectrons to one of the attached atoms first.\n%s %s,%s %s" % (createAtomMapSmiString(mol), bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx(), atom.GetMapIdx()) );
    
    if arrowObjList is not None:
        sourceIndexes = [bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx()];
        targetIndexes = [atom.GetMapIdx()];
        if pivotAtom is not None:
            targetIndexes.append( pivotAtom.GetMapIdx() );
        arrow = ElectronArrow( mol, sourceIndexes, electronCount, targetIndexes );
        
        #Possibly Fix the arrow codes if we have a free radical
        # This will occur if have a source lone radical and a target bond orbital
        # then the movement from both of these 'arrows' should be 'towards' the bond
        # btwn the sourceOrb.atom and the targetOrb.atom (from the moverOrbitalElectrons data structs) 
        if electronCount == 1:
            idxOfFix = -1;
            for i, arrowObj in enumerate(arrowObjList):
                if len(arrowObj.sourceIndexes) == 1 and len(arrowObj.targetIndexes) == 1 and \
                        arrowObj.sourceIndexes[0] == atom.GetMapIdx() and arrowObj.targetIndexes[0] == pivotAtom.GetMapIdx():
                    idxOfFix = i;
                    break;
            if idxOfFix != -1:
                wrongArrow = arrowObjList.pop(idxOfFix);
                fixedArrow = ElectronArrow(mol, wrongArrow.sourceIndexes, electronCount, targetIndexes)
                arrowObjList.append(fixedArrow)
        
        arrowObjList.append( arrow );




class BondDoesNotExistException(Exception):
    """Specialized exception representing an attempt to move bond electrons
    from a bond that no longer exists, probably because user attempted
    to move the bond electrons with more than one electron flow arrow.
    """
    def __str__(self, atom=None, pivotAtom=None, electronCount=None):
        atomStr = str(atom)
        pAtomStr = str(pivotAtom)
        if atom is not None:
            atomStr = atom.GetMapIdx();
        if pivotAtom is not None:
            pAtomStr = pivotAtom.GetMapIdx();
        eCountStr = str(electronCount)
        return \
            """Attempted to move electrons from a bond that no longer exists.  
            This occurs when multiple electron flow 
            arrows are drawn from a single bond.
            atom : %s   pAtom : %s  electrons : %s
            """ % (atomStr, pAtomStr, eCountStr);
            
            
def applyAndCanonicalizeReaction(oldMol, arrowStr=None, removeStereo=False):
    """Convenience to apply arrows and return a can Smi Str."""
    mol = OEGraphMol(oldMol);
    if arrowStr is not None:
        eArrowList = ElectronArrow.parseArrowCodes(mol, arrowStr);
        ElectronArrow.applyAll(mol, eArrowList);
    if removeStereo:
        clearMolStereo(mol);
    clearAtomMaps(mol);
    smi = OECreateCanSmiString(mol);
    del mol;    
    return smi;

def determineTransitionArrows( sourceMol, targetMol, arrowObjList=None, sourceAtomByMapIdx=None, targetAtomByMapIdx=None, sourceIndexSetsUsed=None ):
    """Given 2 molecule objects that should represent alternate forms of the same
    contents and have ***matching non-zero atom map indexes***, 
    determine what electron flow arrows are necessary to
    convert the sourceMol into the targetMol.  
    Return a list of ElectronFlow arrow objects to represent this.
    
    This method largely relies on the caller to use it in good faith.
    Not much pre-testing for valid matching structures, and
    will only check atoms and bonds with map indexes applied.
    
    Currently only supports ionic structures (i.e., movements of 2 electrons, not radicals)
    """
    #log.debug("%s\t%s" % (createStandardSmiString(sourceMol), createStandardSmiString(targetMol)) );
    
    origSourceMol = None;
    
    if arrowObjList is None:
        # Must be initial call.  Initialize tracking data
    
        arrowObjList = list();
        # Create a copy of the source molecule the first time, as we will be modifying it as we go
        origSourceMol = sourceMol;
        sourceMol = OEGraphMol(sourceMol);
    
        sourceAtomByMapIdx = dict();
        for sourceAtom in sourceMol.GetAtoms():
            if sourceAtom.GetMapIdx() > 0:
                sourceAtomByMapIdx[sourceAtom.GetMapIdx()] = sourceAtom;

        targetAtomByMapIdx = dict();
        for targetAtom in targetMol.GetAtoms():
            if targetAtom.GetMapIdx() > 0:
                targetAtomByMapIdx[targetAtom.GetMapIdx()] = targetAtom;
    
    
        # Keep track of arrow source index sets to avoid repeats (thus avoiding infinite recursion)
        sourceIndexSetsUsed = set();
    
    chargeDiffFound = False;
    arrowFound = False;
    for mapIdx, sourceAtom in sourceAtomByMapIdx.items():
        targetAtom = targetAtomByMapIdx[mapIdx];
        
        if arrowFound:
            # Don't continue to loop once found an arrow.  Leave it to recursive calls instead
            break;
        
        # See if a formal charge changed on this atom
        diffCharge = targetAtom.GetFormalCharge() - sourceAtom.GetFormalCharge();
        if (diffCharge != 0):
            chargeDiffFound = True;
        if diffCharge == +1:
            # Look to the atom neighbors and bonds to figure out what to do
            sourceBondByNeighborMapIdx = dict();
            for sourceNeighbor in sourceAtom.GetAtoms():
                if sourceNeighbor.GetMapIdx() > 0:
                    sourceBond = sourceAtom.GetBond( sourceNeighbor );
                    sourceBondByNeighborMapIdx[sourceNeighbor.GetMapIdx()] = sourceBond;

            targetBondByNeighborMapIdx = dict();
            for targetNeighbor in targetAtom.GetAtoms():
                if targetNeighbor.GetMapIdx() > 0:
                    targetBond = targetAtom.GetBond( targetNeighbor );
                    targetBondByNeighborMapIdx[targetNeighbor.GetMapIdx()] = targetBond;

            # Charge increased by +1, possibly a bond was pulled away
            for neighborMapIdx, sourceBond in sourceBondByNeighborMapIdx.items():
                targetBond = None;
                if neighborMapIdx in targetBondByNeighborMapIdx:
                    targetBond = targetBondByNeighborMapIdx[neighborMapIdx];

                if targetBond is None or targetBond.GetOrder() < sourceBond.GetOrder():
                    # A bond was lost in the source that must have been pulled away
                    arrowObj = ElectronArrow(sourceMol, [mapIdx, neighborMapIdx], 2, [neighborMapIdx])
                    isNewSource = applyAddArrow( targetMol, arrowObj, arrowObjList, sourceIndexSetsUsed );
                    if isNewSource:
                        # Found a valid new arrow to apply.  target structure has changed, need to reevaluate from the top.  Can stop the current loop
                        arrowObjList.extend( determineTransitionArrows(sourceMol, targetMol, arrowObjList, sourceAtomByMapIdx, targetAtomByMapIdx, sourceIndexSetsUsed) );
                        arrowFound = True;
                        break;

            # Charge increased by +1 could also be a lone pair pulled away to make a new bond
            for neighborMapIdx, targetBond in targetBondByNeighborMapIdx.items():
                sourceBond = None;
                if neighborMapIdx in sourceBondByNeighborMapIdx:
                    sourceBond = sourceBondByNeighborMapIdx[neighborMapIdx];

                if sourceBond is None or sourceBond.GetOrder() < targetBond.GetOrder():
                    # A new bond was created in the target that must have pulled away the lone pair / charge equivalent
                    arrowObj = ElectronArrow(sourceMol, [mapIdx], 2, [neighborMapIdx])
                    isNewSource = applyAddArrow( sourceMol, arrowObj, arrowObjList, sourceIndexSetsUsed );
                    if isNewSource:
                        # Found a valid new arrow to apply.  Source structure has changed, need to reevaluate from the top.  Can stop the current loop
                        arrowObjList.extend( determineTransitionArrows(sourceMol, targetMol, arrowObjList, sourceAtomByMapIdx, targetAtomByMapIdx, sourceIndexSetsUsed) );
                        arrowFound = True;
                        break;

        elif diffCharge == -1:
            # Charge decreased by one, a bond must have been pulled in, either directly attached or from a neighbor
            # Do not have to do anything here, because a decrease in charge here must correspond to an increase
            #   elsewhere, so the above cases will cover them
            pass;    

    if not chargeDiffFound:
        # Iteration through atoms found no formal charge differences to guide arrow selection.
        # Look for changes in bond orders then
        for mapIdx, sourceAtom in sourceAtomByMapIdx.items():
            targetAtom = targetAtomByMapIdx[mapIdx];

            if arrowFound:
                # Don't continue to loop once found an arrow.  Leave it to recursive calls instead
                break;

            # Look to the atom neighbors and bonds to figure out what to do
            sourceBondByNeighborMapIdx = dict();
            for sourceNeighbor in sourceAtom.GetAtoms():
                if sourceNeighbor.GetMapIdx() > 0:
                    sourceBond = sourceAtom.GetBond( sourceNeighbor );
                    sourceBondByNeighborMapIdx[sourceNeighbor.GetMapIdx()] = sourceBond;

            targetBondByNeighborMapIdx = dict();
            for targetNeighbor in targetAtom.GetAtoms():
                if targetNeighbor.GetMapIdx() > 0:
                    targetBond = targetAtom.GetBond( targetNeighbor );
                    targetBondByNeighborMapIdx[targetNeighbor.GetMapIdx()] = targetBond;

            for neighborMapIdx, sourceBond in sourceBondByNeighborMapIdx.items():
                targetBond = None;
                if neighborMapIdx in targetBondByNeighborMapIdx:
                    targetBond = targetBondByNeighborMapIdx[neighborMapIdx];

                if targetBond is None or targetBond.GetOrder() < sourceBond.GetOrder():
                    # A bond was lost in the source that must have been pulled away
                    
                    # Generate an arrow to heterolytically break the bond, direction based on any prevailing formal charges / atom electronegativity
                    sinkMapIdx = neighborMapIdx;
                    sourceNeighbor = sourceBond.GetNbr(sourceAtom);
                    if sourceAtom.GetFormalCharge() > sourceNeighbor.GetFormalCharge():
                        # Move arrow in direction of more positive formal charge
                        sinkMapIdx = mapIdx;
                    elif sourceNeighbor.GetFormalCharge() == sourceAtom.GetFormalCharge():
                        # Same formal charge, so move in direction of general electronegativity
                        if atomElectronegativity(sourceAtom.GetAtomicNum()) > atomElectronegativity(sourceNeighbor.GetAtomicNum()):
                            sinkMapIdx = mapIdx;
                    
                    arrowObj = ElectronArrow(sourceMol, [mapIdx, neighborMapIdx], 2, [sinkMapIdx])
                    isNewSource = applyAddArrow( targetMol, arrowObj, arrowObjList, sourceIndexSetsUsed );
                    if isNewSource:
                        # Found a valid new arrow to apply.  target structure has changed, need to reevaluate from the top.  Can stop the current loop
                        arrowObjList.extend( determineTransitionArrows(sourceMol, targetMol, arrowObjList, sourceAtomByMapIdx, targetAtomByMapIdx, sourceIndexSetsUsed) );
                        arrowFound = True;
                        break;

    # By the end of the loop / recursion, all available arrow opportunities should have been found
    #   and the sourceMol should now match the targetMol.  
    # If not, then these must not be valid / matching structures
    sourceSmi = createStandardSmiString(sourceMol);
    targetSmi = createStandardSmiString(targetMol);
    
    if sourceSmi != targetSmi:
        raise Exception("Unable to find arrow to convert source into target structure: %s vs. %s" % (sourceSmi, targetSmi) );

    if origSourceMol is not None:
        # The arrow objects were actually created against a modifiable copy of the source molecule.
        # Recreate a set of arrow objects based on the original molecule before returning them.
        arrowCodes = ElectronArrow.prepareArrowCodes( arrowObjList );   # Added benefit of merging extra arrows that pass through pivot atoms
        arrowObjList = ElectronArrow.parseArrowCodes( origSourceMol, arrowCodes );
    
    return arrowObjList;

def applyAddArrow( mol, arrowObj, arrowObjList, sourceIndexSetsUsed ):
    """Apply the ElectronArrow object to the mol, add a copy of the arrow to the arrowObjList,
    and the source indexes to the sourceIndexSetsUsed.  First check if the sourceIndexSetsUsed
    already included this source however, and skip if this was a redundant source 
    and return True/False as to whether this was a new source or not.
    """
    isNewSource = (arrowObj.sourceIndexes not in sourceIndexSetsUsed);
    if isNewSource:
        arrowObj.apply( mol );
        arrowObjList.append( arrowObj );
        sourceIndexSetsUsed.add( arrowObj.sourceIndexes );
    return isNewSource;

