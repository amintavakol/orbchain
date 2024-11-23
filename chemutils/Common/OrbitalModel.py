"""Atomic and molecular orbital class.

Key functions for external callers:

    orbitalIter:
        Conceptually should be a member function of the atom class (OEAtomBase), 
        but since we don't own the atom class ourselves, this function can be called on
        one instead to return a list of orbital objects that are related to that atom.
    orbitalInfo:
        Again, should really be an atom member function.  Provides some basic 
        information about the atom to help determine orbital information, such as
        number of lone pairs, empty orbitals, etc.

"""
from pprint import pformat
#from sets import Set;
from openeye.oechem import OEGetAtomicSymbol, OECreateIsoSmiString, OEGraphMol, OEParseSmiles, OEGetHybridization
from openeye.oechem import OEBondStereo_CisTrans, OEBondStereo_Cis, OEBondStereo_Trans
from openeye.oechem import OEAtomStereo_Tetrahedral, OEAtomStereo_Left, OEAtomStereo_Right, OEAtomStereo_Undefined
from openeye.oechem import OEAssignImplicitHydrogens, OEAssignFormalCharges, OEDetermineComponents
from chemutils.Common.MolStdValue import stableValenceShell, nStdValenceElectrons, nMaxHyperValenceElectrons, STD_OCTECT
from chemutils.Common.MolStdValue import PERIODIC_ROW
from chemutils.Common.Const import SMILES_MOL_DELIM
from chemutils.Common.Util import createAtomMapSmiString, log;
from chemutils.Common.Util import atomByMapIdx, getMaxMapIdx;
import string;

#EXPECTED_HYBRIDIZATIONS = set([OEHybridization_sp, OEHybridization_sp2, OEHybridization_sp3])
#DEFAULT_HYBRIDIZATION = OEHybridization_sp3;
"""Expected possible orbital types.
Note sp == sp1, s = sp0
Note sp4, sp5 does not really make sense or exist, but can show up for 
mock pivot orbitals while the application is in the middle of moving electrons 
over a multi-step progression in a hyper-valent system (e.g., sulfuric acid).
"""
ORBITAL_TYPES = set(["sigma","pi","s","sp0","sp","sp1","sp2","sp3","p", "sp3d", "d"]);
MOCK_ORBITAL_TYPES = set(["sp4","sp5"]);
BOND_ORBITAL_TYPES = set(["sigma","pi"]);

ORBITAL_DELIM = ";";
ORBITAL_FACE_DELIM = ","
INTERACTION_SYMBOL = "&gt;";

# Default label increment when laying out atom mapped reaction SMILES with multiple orbitals to label
ORBITAL_LABEL_INCREMENT = 10;

class Orbital:
    mol = None
    atom = None
    type = None
    neighbor = None    # Neighbor atom if this is a bond orbital.
    #Don't store actual bond reference, because that could be lost
    electrons = None
    score = None
    flag = None

    # Dictionary to hold any additional arbitrary data or annotations
    dataTags = None

    # Support serial chaining (extension) to further orbitals (e.g., enolates and E2 eliminations)
    #   with a linked list style structure
    extOrbital = None
    
    # Support having a preferred reactive sp2 atoms.  This is only used in atoms that are 
    # sp2 hybridized to determine on which side should an addition take place.
    # Here store a sequence of the atoms attached such that when looking from the reactive side, these 
    # atoms are in clockwise order.
    # Also used to support a specification of direction of attack for empty p orbitals in sp2 hybridized atoms  
    atomFaceList = None
    nbrFaceList = None
    
    # Here we need to track if an orbital movement involves syn addition.  If so, then we need to flip the source 
    # face when constructing an inverse orb.  Inverse orbs could possibly have both the source and the 
    __flagIsSynAddition = False;
    __flagIsNbrSynAddition = False;
    
    # Keep a flag here to remove face list info when it is used simply for breadcrumbs to be able to move 
    # backwards.  If this is True, the move/undoMoveOrbitalElectron functions should kill facelists when
    # done with the movement.  Basically a flag to remove the breadcrumb when done.
    isFaceListOnlyBreadCrumb = True;
    
    def simpleCopy(self):
        """Simply copy an orb object without extOrbitals  Useful for making extended orbitals.
        
        Handle extOrbitals now in a recursive fashion
        """
        bond = None;
        if self.neighbor is not None:
            bond = self.atom.GetBond(self.neighbor)
            
        nOrb = Orbital(self.atom, self.type, bond, self.electrons, self.mol, self.atomFaceList, self.nbrFaceList)
        if self.extOrbital is not None:
            copyExtOrbital = self.extOrbital.simpleCopy()
            nOrb.extOrbital = copyExtOrbital;
        
        return nOrb
    
    
    def __init__(self, atom, type, bond, electrons, mol=None, atomFaceList=None, nbrFaceList=None):
        if type is None:
            # Try to guess the type if the caller did not specify
            if bond is not None:
                # Bond orbital
                type = "sigma"
                if bond.GetOrder() > 1:
                    type = "pi"
            else:
                # Non-bonding orbital, either some kind of lone pair, radical or empty
                if electrons == 0:
                    ## Check if there really is 
                    type = "p"
                else:
                    #type = "sp%d" % atomHybridization(atom);
                    hybridNum = atomHybridization(atom);
                    if hybridNum < 4:
                        type = "sp%d" % atomHybridization(atom);
                    else:
                        spHyb = hybridNum % 10
                        if spHyb < 4:
                            type = 'sp%d' % (spHyb) ## Then lone pair and is sp3 
                        else:
                            type = 'sp3d'  # This is only case (if there truly is a lone pair)
                            #type = "sp%dd%d" % ((hydridNum % 10), int(hydridNum/10));
        #if type == 'sigma' and bond is None:
        #    raise Exception('Tried to init orbital as sigma with NO BOND!')

        if type == "sp1": type = "sp"
        if type == "sp0": type = "s"
        
        if mol is None:
            # Keep a reference to the mol object, otherwise the wrapped C reference (atom.GetParent())
            #   may become lost, and garbage collection will inadvertently recycle the mol object,
            #   but even this seems imperfect.  The GetParent reference seems non-identical to the 
            #   original source.  For assurance, caller should directly set the mol reference.
            mol = atom.GetParent()
        self.mol = mol
        self.atom = atom
        self.type = type
        if bond is not None:
            self.neighbor = bond.GetNbr(atom)
        self.electrons = electrons
        self.score = None
        self.flag = None
        if type not in ORBITAL_TYPES:
            if type in MOCK_ORBITAL_TYPES:
                log.warning("Using mock orbital type: %s - %s (%s)" % (type, createAtomMapSmiString(mol), atom.GetMapIdx()) )
            else:
                # Not even a mock orbital, who knows what the caller is trying to produce
                raise Exception("Unexpected orbital type: %s - %s (%s)" % (type, createAtomMapSmiString(mol), atom.GetMapIdx()) );

        self.dataTags = dict();
        
        self.extOrbital = None;
        
        self.atomFaceList = atomFaceList;
        self.nbrFaceList = nbrFaceList;
        self.__flagIsSynAddition = False;
        if self.atomFaceList is not None or self.nbrFaceList is not None:
            self.isFaceListOnlyBreadCrumb = False;
    
    
    def removeFaceBreadCrumbs(self):
        """Simple recursive method to remove face information for self and extOrbitals if isFaceListOnlyBreadCrumb is set"""
        if self.isFaceListOnlyBreadCrumb:
            self.atomFaceList = None;
            self.nbrFaceList = None;
        if self.extOrbital is not None:
            self.extOrbital.removeFaceBreadCrumbs();
    
    
    def isPerceivedFreeRadical(self):
        """Simple recursive method that checks if any of the orbs in chain have a single electron.  
        
        Used to help 'guess' if an interaction will involve a free radical."""
        if self.electrons == 1:
            return True;
        if self.extOrbital is None:
            return False;
            
        return self.extOrbital.isPerceivedFreeRadical();
    
    
    def setSynAdditionFlags(self):
        """Function to determine if an interaction with the bond orbital represents a syn addition.
        
        Wrap this, so the flag itself is a 'private' var of the orbital, caller should only access this flag
        through this function or the accessor below.
        
        There are two cases:
            -Tetrahedral - the face specified on the atom is the opposite of what would be 
                            if the bond was simply broken.
                            The atom or neighbor could fit this case (and the respective flag should be set.)
            -planar - both faces of the pi bond are the same.  *ie, both up or both down. 
                            This will be the case if the first atoms in the face lists of 
                            the source/neighbor (arranged such that the respective bond nbr
                            atom is last) are trans to each other.
        """
        if not self.isBondOrbital() or (self.atomFaceList is None and self.nbrFaceList is None):
            self.__flagIsSynAddition = False;
            self.__flagIsNbrSynAddition = False;
            return;
        
        if self.type == 'sigma':
            if not self.atom.HasStereoSpecified() or self.atomFaceList is None:
                self.__flagIsSynAddition = False;
            else:
                tempAtomList = self.atomFaceList[:];
                tempAtomList.insert(0, self.neighbor);
                self.__flagIsSynAddition = self.atom.GetStereo(tempAtomList, OEAtomStereo_Tetrahedral) == OEAtomStereo_Left;
                
            if not self.neighbor.HasStereoSpecified() or self.nbrFaceList is None:
                self.__flagIsNbrSynAddition = False;
            else:
                tempAtomList = self.nbrFaceList[:];
                tempAtomList.insert(0, self.atom);
                self.__flagIsNbrSynAddition = self.neighbor.GetStereo(tempAtomList, OEAtomStereo_Tetrahedral) == OEAtomStereo_Left;
                
        elif self.type =='pi':
            # No matter what the neighbor should have no corrections.
            self.__flagIsNbrSynAddition = False;
            
            if self.atomFaceList is None or self.nbrFaceList is None or not self.getBond().HasStereoSpecified():
                self.__flagIsSynAddition = False;
                return;
            bondAtomIdx = -1;
            for idx, atm in enumerate(self.atomFaceList):
                if atm == self.neighbor:
                    bondAtomIdx = idx;
                    break;
            atomFaceList = self.atomFaceList[bondAtomIdx+1:];
            atomFaceList.extend(self.atomFaceList[:bondAtomIdx+1]);
        
            bondAtomIdx = -1;
            for idx, atm in enumerate(self.nbrFaceList):
                if atm == self.atom:
                    bondAtomIdx = idx;
                    break;
            nbrFaceList = self.nbrFaceList[bondAtomIdx+1:];
            nbrFaceList.extend(self.nbrFaceList[:bondAtomIdx+1]);
        
            self.__flagIsSynAddition = self.getBond().GetStereo([atomFaceList[0], nbrFaceList[0]], OEBondStereo_CisTrans) == OEBondStereo_Trans;
        
    
    def isSynAddition(self):
        """Simple accessor for the synAddition flag"""
        return self.__flagIsSynAddition;
    
    def isNbrSynAddition(self):
        """Simple accessor for nbrSynAddition flag"""
        return self.__flagIsNbrSynAddition
    
    def flipFaceList(faceList):
        """Method to flip the face of the atom (or neighbor if specified).  This is used to correctly make inverse orbitals"""
        if faceList is not None and len(faceList) >= 3:
            tempList = [faceList[1], faceList[0]]
            tempList.extend(faceList[2:]);
            return tempList;
        
        return faceList;
    flipFaceList = staticmethod(flipFaceList)
    
    def coveredOrbitalAtomIdx(self ):
        """Returns a set of all atom idx covered by this orbital (including extOrbitals)"""
        atomIdx = set([self.atom.GetIdx()]);
        if self.neighbor is not None:
            atomIdx.add(self.neighbor.GetIdx());
        if self.extOrbital is not None:
            atomIdx.update(self.extOrbital.coveredOrbitalAtomIdx());
        return atomIdx; 
    
    
    def labelOrbitalAtoms(self, baseLabel):
        """Set atom map index labels on the molecule object
        this orbital is based upon.  Start with the given baseLabel
        and increment for the neighboring orbital atoms.
        
        Atoms involved in the specification of the reactive pi face are labeled here.
        Note:  Because of the way we label these, the actual numbers are really dependent on the order in the smi
        string.  Also, any 'face' labeled atom will have a different number if involved in an extended orbital.
        """
        self.atom.SetMapIdx(baseLabel);
        baseLabel += 1;

        if self.neighbor is not None:
            self.neighbor.SetMapIdx(baseLabel);
            baseLabel += 1;
            
        # Also label atoms connected to pi-bond atom if faces are specified.
        if self.atomFaceList is not None:
            for nbr in self.atom.GetAtoms():
                if self.neighbor is None or nbr.GetIdx() != self.neighbor.GetIdx():
                    nbr.SetMapIdx(baseLabel);
                    baseLabel += 1;
            if self.neighbor is not None and self.nbrFaceList is not None:
                for nbr in self.neighbor.GetAtoms():
                    if nbr.GetIdx() != self.atom.GetIdx():
                        nbr.SetMapIdx(baseLabel);
                        baseLabel += 1;
        
        if self.extOrbital is not None:
            self.extOrbital.labelOrbitalAtoms(baseLabel);
        
    
    def fromLabeledMol(mol, atomLabel, orbType, neighborLabel, electrons, atomFaceLabelList=None, nbrFaceLabelList=None):
        """Factory method to produce an orbital object based on the
        contents of the given mol.  Instead of specifying an actual atom (and bond)
        object, instead specify the atom map index label that will be found on the
        atom and, optionally, the (bond) neighbor atom.

        For example, given the molecule derived from this SMILES string
            [CH2:21]=[CH:20][O-:10].[Na+]
        
        Specific orbital models can be extracted by the following calls:
        
            # Lone pair on the oxygen
            Orbital.fromLabeledMol( mol, 10, "sp3", None, 2 );
        
            # Carbon-carbon pi bond
            Orbital.fromLabeledMol( mol, 20, "pi", 21, 2 );
            
        Also have an optional faceAtomLabelList which specify the 'side' of a pi-bonded atom which is
        reactive.  These are specified such that looking from reactive face, the atoms are in clockwise order.
        """

        orbAtom = None
        orbBond = None
        atomFaceAtomList = None
        nbrFaceAtomList = None

        for atom in mol.GetAtoms():
            if atom.GetMapIdx() == atomLabel:
                orbAtom = atom
                
                # Look for face specification index labels and translate into atom references
                atomFaceAtomList = None
                if atomFaceLabelList is not None:
                    mapIdxDict = {}
                    for nbr in orbAtom.GetAtoms():
                        mapIdxDict[nbr.GetMapIdx()] = nbr
                    atomFaceAtomList = list()
                    for label in atomFaceLabelList:
                        if label in mapIdxDict:
                            atomFaceAtomList.append( mapIdxDict[label] )
                        else:
                            # Unable to find an atom matching the label, must be a bogus specification, ignore it then
                            log.warning("Unable to find an face atom with label: %s" % label );
                            atomFaceAtomList = None;
                            break;
                
                if neighborLabel is not None:
                    # Check neighbors if this is a bond orbital
                    for neighbor in atom.GetAtoms():
                        if neighbor.GetMapIdx() == neighborLabel:
                            orbBond = orbAtom.GetBond(neighbor);
                            
                            if nbrFaceLabelList is not None:
                                # Find references to these atoms
                                mapIdxDict = {};
                                for nbr in neighbor.GetAtoms():
                                    mapIdxDict[nbr.GetMapIdx()] = nbr

                                nbrFaceAtomList = list()
                                for label in nbrFaceLabelList:
                                    if label in mapIdxDict:
                                        nbrFaceAtomList.append( mapIdxDict[label] )
                                    else:
                                        # Unable to find an atom matching the label, must be a bogus specification, ignore it then
                                        log.warning("Unable to find an face atom with label: %s" % label )
                                        nbrFaceAtomList = None
                                        break
                                
                            break  # No need to look further
                
                break # Not need to look further

        orbital = None

        if orbAtom is not None:
            orbital = Orbital(orbAtom, orbType, orbBond, electrons, mol, 
                                atomFaceList=atomFaceAtomList, nbrFaceList=nbrFaceAtomList)

        return orbital

    fromLabeledMol = staticmethod(fromLabeledMol)
    
    def involvedAtomList(self):
        """Function to return a list of all the atom/neighbors in self and all extOrbitals.
        
        NOTE: There is danger of an infinite loop here if the extOrbitals are linked.  Use with caution.
        """
        aList = set([self.atom])
        if self.neighbor is not None:
            if self.neighbor not in aList:
                aList.add(self.neighbor);
            else: 
                return aList;
            
        if self.extOrbital is not None:
            if self.extOrbital.atom not in aList:
                aList.update(self.extOrbital.involvedAtomList())
        
        return aList;
    
    
    def isMultipleComponent(self):
        """Function to test if the orbtial spans mutliple disconnected components.
        
        This will be the case in at least one orbital in pericyclic or pseudo-pericyclic 
        (bromination of an alkene) reactions
        """
        count, parts = OEDetermineComponents(self.mol)
        involvedAtomParts = set([parts[atm.GetIdx()] for atm in self.involvedAtomList()])
        if len(involvedAtomParts) > 1:
            return True
        
        return False
    
    
    def fromLabeledMolAndInfoStr(mol, orbInfoStr):
        """Even more general method to facilitate conversions 
        to / from simple strings.  Even supports the option
        of extended orbital chains with the orbital info strings
        representing a delimited list.
        
        Allow specification of the reactive face of a pi-bond orbital.  
        """
        headOrb = None
        prevOrb = None
        
        if orbInfoStr is not None:
            orbInfoStrList = orbInfoStr.split(ORBITAL_DELIM)
            for componentOrbInfoStr in orbInfoStrList:
                tokens = componentOrbInfoStr.split()    # Assume simple white-space separation

                atomLabel = int(tokens[0])
                orbType = tokens[1]
                neighborLabel = None
                if tokens[2] != str(None):
                    neighborLabel = int(tokens[2])
                electrons = int(tokens[3])
                
                #If the information is given, label the faces
                atomFaceLabelList = None
                nbrFaceLabelList = None
                if len(tokens) >= 5:
                    if tokens[4] != str(None):
                        rawLabels = tokens[4].split(',')
                        atomFaceLabelList = [int(label) for label in rawLabels]
                if len(tokens) >= 6:
                    if tokens[5] != str(None):
                        rawLabels = tokens[5].split(',')
                        nbrFaceLabelList = [int(label) for label in rawLabels]
                
                orb = Orbital.fromLabeledMol(mol, atomLabel, orbType, neighborLabel, electrons, atomFaceLabelList, nbrFaceLabelList)

                if headOrb is None:
                    # First orbital encountered, store as the head orbital
                    headOrb = orb
                    prevOrb = orb
                else:
                    # Have previous orbitals, this must be further down an extended chain
                    prevOrb.extOrbital = orb
                    prevOrb = orb
        else:
            # Valid option to return a None orbital, as a sentinel object
            #   to represent the "filled orb" of a bond dissociation reaction.
            pass
        return headOrb
    
    fromLabeledMolAndInfoStr = staticmethod(fromLabeledMolAndInfoStr)

    def fromLabeledMolAndOrbital(mol, baseAtomLabel, templateOrb):
        """Likely used if reproducing the template orbital, just
        in a different molecule object copy.
        
        The labels for reactive pi face are reproduced here.  However, we are using actual mapIdx from 
        templateOrb instead of incrementing baseAtomLabel to not run into problems when some face atoms 
        not truly involved in the orb mess up the baseAtom label numbering.
        """
        neighborLabel = None;
        if templateOrb.neighbor is not None:
            neighborLabel = templateOrb.neighbor.GetMapIdx();
        
        # If we have it, grab the face information.
        atomFaceLabelList = None;
        nbrFaceLabelList = None;
        if templateOrb.atomFaceList is not None:
            atomFaceLabelList = [atm.GetMapIdx() for atm in templateOrb.atomFaceList];
        if templateOrb.nbrFaceList is not None:
            nbrFaceLabelList = [atm.GetMapIdx() for atm in templateOrb.nbrFaceList];
        
        
        newOrb = Orbital.fromLabeledMol(mol, baseAtomLabel, templateOrb.type, baseAtomLabel + 1, templateOrb.electrons,
                                                atomFaceLabelList, nbrFaceLabelList);
        
        if templateOrb.extOrbital is not None:
            newOrb.extOrbital = Orbital.fromLabeledMolAndOrbital( mol, templateOrb.extOrbital.atom.GetMapIdx(), 
                                                        templateOrb.extOrbital );
        
        return newOrb;
    fromLabeledMolAndOrbital = staticmethod(fromLabeledMolAndOrbital);

    def toLabeledSmiAndInfoStr(self, baseLabel=None):
        """Generate a SMILES string representing the (composite) molecule
        that this orbital is based on, with atom map labels to identify
        where the orbital applies to.  This should then be usable to
        reconstruct the orbital information with the fromLabeledMol function.
        
        The primary orbital atom will be given the atom map index specified
        by baseLabel.  All other orbital related atoms (i.e., the neighbor atom
        in a bond orbital) will be mapped with progressively incremented indexes:
        baseLabel+1, baseLabel+2, etc..
        
        Alternatively, if the baseLabel is unspecified, will just use whatever
        atom map labels already exist on the molecule.  Beware that this is
        risky as inconsistencies will arrive if run on an orbital whose 
        participating atoms have no atom map labels on them.
        """
        # Cannot prepare a modifiable copy of the orbital's molecule object
        #   because the orbital object would still be pointing at the original.
        mol = self.mol
        
        # Instead, need to keep track of all changes made so they can be undone,
        #   reverting the parent molecule back to its original state.
        originalAtomMaps = None

        if baseLabel is not None:
            originalAtomMaps = dict()
            # Wipe out any existing atom map index labels, but keep track so they can be restored
            for atom in mol.GetAtoms():
                originalAtomMaps[atom.GetIdx()] = atom.GetMapIdx()
                atom.SetMapIdx(0)

            # Apply labels based on the current orbital
            self.labelOrbitalAtoms(baseLabel)
        
        labeledSmi = createAtomMapSmiString(mol)

        infoStr = self.formatInfoStr();

        if baseLabel is not None:
            # Revert the molecule back to its original state before leaving
            for atom in mol.GetAtoms():
                atom.SetMapIdx( originalAtomMaps[atom.GetIdx()] )

        return (labeledSmi, infoStr)

    def formatInfoStr(self):
        """Return basic information string on the orbital which
        can be parsed back into an object with an accompanying labeled molecule.
        """
        infoList = [];
        infoList.append( str(self.atom.GetMapIdx()) )
        infoList.append( self.type );
        if self.neighbor is not None:
            infoList.append( str(self.neighbor.GetMapIdx()) );
        else:
            infoList.append( str(None) );
        infoList.append( str(self.electrons) );
        
        # Add on the face data if we have it
        if self.atomFaceList is not None or self.nbrFaceList is not None:
            if self.atomFaceList is not None:
                atomFaceLabelList = [str(atm.GetMapIdx()) for atm in self.atomFaceList] 
                infoList.append( str.join(ORBITAL_FACE_DELIM, atomFaceLabelList));
            else:
                infoList.append(str(None))
            if self.nbrFaceList is not None:
                nbrFaceLabelList = [str(atm.GetMapIdx()) for atm in self.nbrFaceList]
                infoList.append( str.join(ORBITAL_FACE_DELIM, nbrFaceLabelList));
            else:
                infoList.append( str(None));
        
        infoStr = str.join("\t", infoList );
        
        if self.extOrbital is not None:
            extInfoStr = self.extOrbital.formatInfoStr();
            infoStr = "%s%s %s" % (infoStr, ORBITAL_DELIM, extInfoStr);
        
        return infoStr;
    
    def cloneOrbitalsWithCommonMol( orbitalList ):
        """Assuming all of the orbitals in the list are based on
        the same molecule object.  Create a clone of the molecule,
        but also a matching set of cloned orbital objects
        which correspond to the original orbitals.
        """
        originalMol = orbitalList[0].mol
        
        copyMol = OEGraphMol(originalMol)
        
        # Enable atom retrieval by index
        copyAtomsByIndex = dict()
        for atom in copyMol.GetAtoms():
            copyAtomsByIndex[atom.GetIdx()] = atom
        
        copyOrbitalList = []
        for originalOrb in orbitalList:
            # Assume that the mol copy was an exact copy, including atom index positions
            copyOrb = Orbital.cloneOrbital(originalOrb, copyMol, copyAtomsByIndex)

            copyOrbitalList.append( copyOrb )
            
        # Does not seem necessary to return copyMol since that is implicitly
        #   attached under the orbital.mol reference, but getting
        #   errors without it, seems like incorrect treating the molecule
        #   object as going out of scope and thus garbage collected.
        return (copyOrbitalList, copyMol)
    cloneOrbitalsWithCommonMol = staticmethod(cloneOrbitalsWithCommonMol)

    def cloneOrbital( originalOrb, copyMol, copyAtomsByIndex ):
        """Utility function for cloneOrbitalWithcommonMol,
        implemented recursively to handle extended orbital chains.
        """
        copyBaseAtom = copyAtomsByIndex[ originalOrb.atom.GetIdx() ];   
        copyBond = None;
        if originalOrb.neighbor is not None:
            copyNeighbor = copyAtomsByIndex[ originalOrb.neighbor.GetIdx() ];
            copyBond = copyBaseAtom.GetBond( copyNeighbor );
            
        #Carry over the face information
        copyAtomFaceList = None;
        copyNbrFaceList = None;
        if originalOrb.atomFaceList is not None:
            copyAtomFaceList = [copyAtomsByIndex[atm.GetIdx()] for atm in originalOrb.atomFaceList];
        if originalOrb.nbrFaceList is not None:
            copyNbrFaceList = [copyAtomsByIndex[atm.GetIdx()] for atm in originalOrb.nbrFaceList]
            
        copyOrb = Orbital(copyBaseAtom, originalOrb.type, copyBond, originalOrb.electrons, mol=copyMol,
                                atomFaceList=copyAtomFaceList, nbrFaceList=copyNbrFaceList);

        if originalOrb.extOrbital is not None:
            # Extended orbital chain, recursively continue cloning
            copyExt = Orbital.cloneOrbital( originalOrb.extOrbital, copyMol, copyAtomsByIndex );
            copyOrb.extOrbital = copyExt;
        return copyOrb;            
    cloneOrbital = staticmethod(cloneOrbital);

    def compositeOrbitalsFromDistinctMols( orbitalList ):
        """Assuming all of the orbitals in the list are based on 
        distinct molecule objects.  Create a single composite molecule object 
        which contains clones of all the original molecules,
        with a matching set of cloned orbital objects
        which correspond to the original orbitals, but pointint towards the composite.
        """
        orbSmiList = [];
        orbInfoStrList = [];
        for iOrb, orb in enumerate(orbitalList):
            orbBaseLabel = ORBITAL_LABEL_INCREMENT * (iOrb+1);
            (orbSmi, orbInfoStr) = orb.toLabeledSmiAndInfoStr(orbBaseLabel);
            log.debug('%d: orbSmi: %s info :%s' %(iOrb, orbSmi, str(orbInfoStr)))
            orbSmiList.append(orbSmi);
            orbInfoStrList.append(orbInfoStr);

        compositeSmi = str.join(SMILES_MOL_DELIM, orbSmiList);
        log.debug(compositeSmi)
        log.debug('Comp Info Str List : %s' % str(orbInfoStrList))
        compositeMol = OEGraphMol();
        OEParseSmiles(compositeMol, compositeSmi);
        
        
        compositeOrbList = [];
        for iOrb, componentOrb in enumerate(orbitalList):
            orbBaseLabel = ORBITAL_LABEL_INCREMENT * (iOrb+1);
            log.debug('orbBaseLabel : %d, componentOrb : %s' % (orbBaseLabel, str(componentOrb.toLabeledSmiAndInfoStr())))
            #compositeOrb = Orbital.fromLabeledMolAndOrbital( compositeMol, orbBaseLabel, componentOrb );
            compositeOrb = Orbital.fromLabeledMolAndInfoStr(compositeMol, orbInfoStrList[iOrb])
            compositeOrbList.append( compositeOrb );
            log.debug('After return compositeOrb : %s' % str(compositeOrb.toLabeledSmiAndInfoStr()))

        return (compositeOrbList, compositeMol);

    compositeOrbitalsFromDistinctMols = staticmethod(compositeOrbitalsFromDistinctMols);

    def getBond(self):
        """If a neighbor atom was specified, this should be a bond orbital.
        Get a reference to the bond object indirectly.  Don't depend on a direct
        reference because "moveOrbitalElectrons" functions can delete bonds,
        causing such references to become invalid.
        """
        if self.neighbor is not None:
            return self.atom.GetBond(self.neighbor);
        return None;

    def __eq__(self, other):
        return self.equals(other, True);

    def __ne__(self, other):
        return not self.__eq__(other);

    def equals(self, other, matchReferences=True):
        eq = True;
        eq = eq and self.type == other.type;
        eq = eq and self.electrons == other.electrons;

        if matchReferences:
            eq = eq and (self.atom is None) == (other.atom is None);
            if self.atom is not None:
                eq = eq and self.atom.GetIdx() == other.atom.GetIdx();
            
            eq = eq and (self.neighbor is None) == (other.neighbor is None);
            if self.neighbor is not None:
                eq = eq and self.neighbor.GetIdx() == other.neighbor.GetIdx();

        eq = eq and (self.extOrbital is None) == (other.extOrbital is None);
        if self.extOrbital is not None:
            eq = eq and self.extOrbital.equals( other.extOrbital, matchReferences );
        
        return eq;

    def __str__(self):
        atomStr = "None";
        if self.atom is not None:
            atomStr = OEGetAtomicSymbol(self.atom.GetAtomicNum())
            if self.atom.GetMapIdx() > 0:
                atomStr += ":%d" % self.atom.GetMapIdx();
        bondStr = "None";
        bond = self.getBond();
        if bond is not None:
            neighborAtom = bond.GetNbr(self.atom);
            bondStr = "-%s-%s" % ( bond.GetOrder(), OEGetAtomicSymbol(neighborAtom.GetAtomicNum()) );
            if neighborAtom.GetMapIdx() > 0:
                bondStr += ":%d" % neighborAtom.GetMapIdx();
        elif self.neighbor is not None:
            bondStr = "-?-%s" % OEGetAtomicSymbol(self.neighbor.GetAtomicNum());
            if self.neighbor.GetMapIdx() > 0:
                bondStr += ":%d" % self.neighbor.GetMapIdx();
        orbStr = "%s %s %s %s" % ( atomStr, self.type, bondStr, self.electrons );
        
        if self.atomFaceList is not None or self.nbrFaceList is not None:
            # Add in face info.., 
            atomFaceStr = "None";
            nbrFaceStr = "None";
            if self.atomFaceList is not None:
                atomFaceLabelList = [str(atm.GetMapIdx()) for atm in self.atomFaceList];
                atomFaceStr = str.join(ORBITAL_FACE_DELIM, atomFaceLabelList);
            if self.nbrFaceList is not None:
                nbrFaceLabelList = [str(atm.GetMapIdx()) for atm in self.nbrFaceList];
                nbrFaceStr = str.join(ORBITAL_FACE_DELIM, nbrFaceLabelList);
            orbStr = "%s %s %s" % (orbStr, atomFaceStr, nbrFaceStr);
        
        if self.extOrbital is not None:
            orbStr = "%s%s %s" % ( orbStr, ORBITAL_DELIM, str(self.extOrbital) ); 
        return orbStr;

    def htmlLabel( self, occupied=True ):
        """Prepare an HTML string label to represent the orbital.
        
        For now, face information is not represented.
        """
        htmlList = [];

        if self.extOrbital is not None and occupied:
            # Render extended source orbitals before (post-traversal)
            htmlList.append( self.extOrbital.htmlLabel(occupied) );

        typeStr = "n"
        if self.isBondOrbital():
            typeStr = "&%s;" % self.type;
            if not occupied:
                typeStr += "<sup>*</sup>"
        elif self.isEmpty():
            typeStr = "p";

        atomStr = "None";
        if self.atom is not None:
            atomStr = OEGetAtomicSymbol(self.atom.GetAtomicNum())

        bondStr = "";
        if self.neighbor is not None:
            bondSymbol = "-";    # Default to generic dash
            bond = self.getBond();
            if bond is not None:
                # Try to make bond symbol more specific if possible
                if bond.IsAromatic():
                    bondSymbol = ":";
                elif bond.GetOrder() == 2:
                    bondSymbol = "=";
                elif bond.GetOrder() == 3:
                    bondSymbol = "#";
            bondStr = "%s%s" % ( bondSymbol, OEGetAtomicSymbol(self.neighbor.GetAtomicNum()) );

        htmlList.append("%s<sub>%s%s</sub>" % (typeStr, atomStr, bondStr) );
        
        if self.extOrbital is not None and not occupied:
            # Render extended target orbitals after (pre-traversal)
            htmlList.append( self.extOrbital.htmlLabel(occupied) );

        return str.join(" %s " % INTERACTION_SYMBOL, htmlList );

    def isBondOrbital(self):
        return self.type in BOND_ORBITAL_TYPES;

    def isEmpty(self):
        return self.electrons < 1;

    def isLonePair(self):
        return not self.isBondOrbital() and self.electrons == 2;

    def isLoneRadical(self):
        return not self.isBondOrbital() and self.electrons == 1;

    def inPiSystem(self):
        result = self.type in ("p","pi");
        # Assume sp3 arranged lone pairs could be p arranged (sp2 total) partially
        result = result or self.type == "sp3"; 
        # sp2 not as common, especially if geometrically constrained, 
        #   but could be the case for separated triple bonds 
        #   (nitriles, carbon monoxide, acylium ions, etc.)
        result = result or (self.type == "sp2" and not self.atom.IsInRing());
        return result;

    def pCharacter(self):
        """Intended to be called on non-bond orbitals only, gives the p character
        of the orbitals hybridization.  Essentially equal to the hybridization
        integer value (3 for sp3, 2 for sp2, etc.)
        """
        if self.type == "sp3":  return 3;
        if self.type == "sp2":  return 2;
        if self.type in ("sp","sp1"):   return 1;
        return 0;

    def getNetFormalCharge(self):
        """Return the net formal charge.  Same as the atom formal charge for atomic orbitals,
        otherwise the some of atom formal charges for bond / molecular orbitals.
        """
        netCharge = self.atom.GetFormalCharge();
        if self.isBondOrbital():
            if self.neighbor is not None:
                # If was None, assume means an implicit hydrogen (uncharged)
                netCharge += self.neighbor.GetFormalCharge();
        return netCharge

    def flippedBondOrbital(self):
        """For bond orbitals, return a copy of the orbital but with the polarity swapped.
        That is, the primary atom and neighbor atom assignments will be exchanged.
        """
        if not self.isBondOrbital():
            raise Exception("Trying to flip a non-bond orbital: %s" % str(self) );
        flipOrb = Orbital( self.neighbor, self.type, None, self.electrons, self.mol, self.nbrFaceList, self.atomFaceList);
        flipOrb.neighbor = self.atom;    # Add after the fact, because cannot count on a bond existing in the current state
        
        # Syn Addition flags should propagate down.., 
        flipOrb.__flagIsSynAddition = self.__flagIsNbrSynAddition;
        flipOrb.__flagIsNbrSynAddition = self.__flagIsSynAddition;
        
        return flipOrb
        

def orbitalIter( atom, includeRedundant=False, useSymmetries=False):
    """Iterate over all of the apparent orbitals in the current atom.
    These will include sigma and pi molecular orbitals or possible sp, p, or d atomic orbitals.
    For each, will have a bond associated with it, or a number of unbonded electrons.
    
    includeRedundant:
        If set, will include all orbitals.  Otherwise, normally
        avoids yielding functionally redundant or uninteresting orbitals.  
        In particular:
            - Atoms with multiple lone pairs, only yield 1
            - Atoms with multiple implict hydrogens, only yield 1
            - Double and triple bonds, only yield the pi bond component, not the sigma bond
            - Triple bonds, only yield one pi bond, not both
            - for atoms with hypervalent capability, only yield a max of one empty d orb
            
    
    useSymmetries:
        - assumes that OEPercieveSymmetry was properly called on mol beforehand.
        - will only propose a single version of a bond orbital over symmetries (unless 
         includeRedundant)
    
    """
    spCount = 0;
    pCount = 0;
    dCount = 0;
    seenHydrogen = False;
    seenSymClassIdx = set([])
    nbrSymClass = 0
    for bond in atom.GetBonds():
        if useSymmetries:
            nbrSymClass = bond.GetNbr(atom).GetSymmetryClass()
        # Each bond will include a sigma component
        if includeRedundant or bond.GetOrder() < 2:
            ## But only yield if truly redundant or we haven't seen this yet
            if includeRedundant or (not useSymmetries) or (nbrSymClass not in seenSymClassIdx):
                orb = Orbital(atom, "sigma", bond, 2);
                yield orb;
        spCount += 1;    
        
        # Higher order bonds consist of additional pi components
        for iPi in range(1, bond.GetOrder()):
            if includeRedundant or iPi < 1+1:
                ## Same thing, only yield if not seen before
                if includeRedundant or (not useSymmetries) or (nbrSymClass not in seenSymClassIdx):
                    orb = Orbital(atom, "pi", bond, 2);
                    yield orb;
            pCount += 1;
            
        # Keep track to see if we have seen a sigma bond to a hydrogen.
        seenHydrogen = seenHydrogen or bond.GetNbr(atom).GetAtomicNum() == 1;
        if useSymmetries:
            seenSymClassIdx.add(nbrSymClass)
    
    # Implicit hydrogens
    for iHydrogen in range(atom.GetImplicitHCount()):
        if includeRedundant or (iHydrogen < 1 and not seenHydrogen):
            orb = Orbital(atom, "sigma", None, 2);
            log.debug('About to yield an implicit Hydrogen orbital: %s ' % str(orb.toLabeledSmiAndInfoStr()));
            yield orb;
        spCount += 1;
    
    # Collect statistics to figure out implicit orbitals
    orbInfo = orbitalInfo(atom);
    hybridization = atomHybridization(atom);   # Guess the atom's hybridization state
    stableShell = stableValenceShell( atom.GetAtomicNum() );    # How many electrons needed to make a stable outer shell?
    maxHyperValElectrons = nMaxHyperValenceElectrons(atom.GetAtomicNum())
    
    
    expectedOrbitals = stableShell / 2; 
    dHybrid = hybridization / 10
    spHybrid = hybridization % 10
    expectedSP  = spHybrid + 1; # Assuming OEChem indexing, expected sp orbitals is the index + 1;  e.g., sp3 should have 4
    expectedSP  = min(expectedSP, expectedOrbitals);    # Hydrogen won't be sp3 really, just base on expected orbitals
    expectedP   = expectedOrbitals - expectedSP;   # Assuming 4 orbitals for 8 stable valence electron shell
    expectedD   = maxHyperValElectrons / 2;    # Are there potential d orbitals?
    
    # Non-bonded orbitals
    for iLonePair in range(int(orbInfo["nLonePairs"])):
        orbType = None;
        if dHybrid == 0:
            if spCount < expectedSP:
                orbType = "sp%d" % spHybrid;
                spCount += 1;
            else:
                orbType = "p";
                pCount += 1;
        else:
            ## Means we have a lone pair, AND a d orbital is already involved in bonding
            ## The possible cases are: sulfoxide or tetra-coordinate neutral Sulfur, or
            ##  neg charge penta-coord sulfur
            if spCount < expectedSP:
                orbType = "sp%d" % (spCount);
                spCount += 1;
            else:
                orbType = "sp3d";
                dCount += 1;
            
        if includeRedundant or iLonePair < 1:
            orb = Orbital( atom, orbType, None, 2 );
            yield orb;
    
    ## Right now, we do NOT handle hypervalent radicals.
    for iRadical in range(orbInfo["nRadicals"]):
        orbType = None;
        if spCount < expectedSP:
            orbType = "sp%d" % hybridization;
            spCount += 1;
        else:
            orbType = "p";
            pCount += 1;
        
        orb = Orbital( atom, orbType, None, 1 );
        yield orb;
    
    # Account for any empty orbitals as well
    while spCount < expectedSP:
        orb = Orbital(atom, "sp%d" % hybridization, None, 0);
        yield orb;
        spCount += 1;
    
    while pCount < expectedP:
        orb = Orbital(atom, "p", None, 0);
        yield orb;
        pCount += 1;
    
    ## Need to fix up the dCounts:
    if spCount + pCount > expectedOrbitals:
        ## Captures any pi bonds incorrectly counted as p orbitals
        dCount += spCount + pCount - expectedOrbitals;
    
    nEmptyDOrbs = 0
    while dCount < expectedD:
        if includeRedundant or nEmptyDOrbs < 1:
            orb = Orbital(atom, "d", None, 0);
            yield orb;
            nEmptyDOrbs += 1
        dCount += 1;


def radicalOrbitalIter(atom, includeRedundant=False, useSymmetries=False):
    """Simple method for iterating over possible single electron (radical) orbitals of an atom

    The behavior of includeRedundant and useSymmetries are similar to the behavior in function above.

    This method will give bond orbitals (with a single electron), and single radical orbitals.
    The one tricky thing this method does is to handle carbenes and nitrenes, treating the lone pair as
    up to two potential single radical orbitals.
    """
    seenHydrogen = False
    seenSymClassIdx = set([])
    nbrSymClass = 0
    
    for bond in atom.GetBonds():
        if useSymmetries:
            nbrSymClass = bond.GetNbr(atom).GetSymmetryClass()
        ## Each bond has a sigma component
        if includeRedundant or bond.GetOrder() < 2:
            if includeRedundant or (not useSymmetries) or (nbrSymClass not in seenSymClassIdx):
                orb = Orbital(atom, 'sigma', bond, 1)
                yield orb

        ## Then pi bonds
        for iPi in range(1, bond.GetOrder()):
            if includeRedundant or iPi < 2:
                if includeRedundant or (not useSymmetries) or (nbrSymClass not in seenSymClassIdx):
                    orb = Orbital(atom, 'pi', bond, 1)
                    yield orb

        ## Keep track of hydrogens and the symmClass
        seenHydrogen = seenHydrogen or bond.GetNbr(atom).GetAtomicNum() == 1
        if useSymmetries:
            seenSymClassIdx.add(nbrSymClass)

    ## Handle any implicit hydrogens
    for iHydrogen in range(atom.GetImplicitHCount()):
        if includeRedundant or (iHydrogen < 1 and not seenHydrogen):
            orb = Orbital(atom,'sigma', None, 1)
            log.debug('About to yield an implicit Hydrogen orbital: %s ' % str(orb.toLabeledSmiAndInfoStr()));
            yield orb;
        
    ## Now its all about the implicit orbitals.
    orbInfo = orbitalInfo(atom)
    hybridization = atomHybridization(atom);   # Guess the atom's hybridization state
    spHybrid = hybridization % 10;
    
    ## Then loop through and iterate over radicals
    for iRadical in range(orbInfo['nRadicals']):
        if includeRedundant or iRadical < 1:
            orb = Orbital(atom, 'sp%d' % spHybrid, None, 1)
            yield orb

    ## And a final check - do we have a carbene/nitrene.
    isSingleSpinOrbs = orbInfo['nLonePairs'] > 0 and orbInfo['nEmptyOrbitals'] > 0
    if isSingleSpinOrbs:
        orb = Orbital(atom, 'sp%d' % spHybrid, None, 1)
        yield orb
        if includeRedundant:
            ## Simply a copy (two single electron orbitals)
            orb = Orbital(atom, 'sp%d' % spHybrid, None, 1)
            yield orb


def orbitalIterByMapIdx(mol, mapIdx):
    """Convenience function to grab an orb iter from a mol and mapping index"""
    atm = atomByMapIdx(mol, mapIdx);
    if atm is None:
        return None;
    return orbitalIter(atm);
    
def bondOrbitalByMapIdx(mol, sourceIdx, neighborIdx, piFirst=True, giveDefaultOnMissing=True):
    """Convenience function to get a bond orbital given mapping idx.
    
    If a bond orbital between these two atoms does not exist, create a sigma one as a dummy,
    this is going to be an intermediate orbital.
    """
    orbIter = orbitalIterByMapIdx(mol, sourceIdx);
    orb = [orb for orb in orbIter if orb.neighbor is not None and orb.neighbor.GetMapIdx() == neighborIdx];
    
    if orb == []:
        if giveDefaultOnMissing:
            atom = atomByMapIdx(mol, sourceIdx)
            orb = Orbital(atom, 'sigma', None, 2, mol)
            orb.neighbor = atomByMapIdx(mol, neighborIdx);
            return orb;
        else:
            return None;
    elif piFirst:
        return orb[-1];
    else:
        return orb[0];
        
def atomOrbitalByMapIdx(mol, atomIdx, numElectrons, giveDefaultOnMissing=True):
    """Convenience function to get an atomic orbital by map idx
    
    If an atom orbital does not currently exist, then make a dummy one that is empty.
    """
    orbIter = orbitalIterByMapIdx(mol, atomIdx);
    orb = [orb for orb in orbIter if orb.type not in BOND_ORBITAL_TYPES and \
                orb.electrons == numElectrons];
    if orb == []:
        if giveDefaultOnMissing:
            atom = atomByMapIdx(mol, atomIdx);
            orb = Orbital(atom, 'p', None, 0, mol);
            return orb;
        else:
            return None;
    else:
        return orb[0];
    

def orbitalInfo( atom ):
    """Returns various statistics related to the atom's non-bonded orbitals.
    This is a dictionary with keys:
        - nUnbondedElectrons
        - nLonePairs
        - nRadicals
        - nEmptyOrbitals
        - nPossDOrbitals - Is there hypervalent possibility?

    Given an atom object, figure out how many electrons are in its outer
    valence shell that are NOT part of a bond.  Based on standard valence electron number.
    """
    infoDict = dict();

    expectedOrbitals = stableValenceShell(atom.GetAtomicNum()) / 2
    nElectrons = nStdValenceElectrons( atom.GetAtomicNum() )
    nElectrons -= atom.GetFormalCharge()
    for bond in atom.GetBonds():
        expectedOrbitals -= bond.GetOrder()
        nElectrons -= bond.GetOrder()
    expectedOrbitals -= atom.GetImplicitHCount()
    nElectrons -= atom.GetImplicitHCount()

    nLonePairs = nElectrons / 2    # Integer division
    nRadicals = nElectrons % 2     # Modulo.  Assume pairs will always try to stick together for lowest energy orbitals
    
    expectedOrbitals -= nLonePairs
    expectedOrbitals -= nRadicals
    
    nPossDOrbitals = 0
    if PERIODIC_ROW[atom.GetAtomicNum()] > 2:
        nPossDOrbitals = 2
    #max(PERIODIC_ROW[atom.GetAtomicNum()] - 2, 0)
    
    infoDict["nUnbondedElectrons"] = nElectrons
    infoDict["nLonePairs"] = nLonePairs
    infoDict["nRadicals"] = nRadicals
    infoDict["nEmptyOrbitals"] = expectedOrbitals  # Anything left must be empty
    infoDict['nPossDOrbitals'] = nPossDOrbitals

    return infoDict
    
def isOrbitalOverlap(sourceOrb, sinkOrb):
    """Convenience function that there are not atom overlaps.
    """
    sourceIdx = sourceOrb.coveredOrbitalAtomIdx();
    sinkIdx = sinkOrb.coveredOrbitalAtomIdx();
    return sinkIdx.intersection(sourceIdx) != set([])

def clearNonInvolvedAtomMaps(orbitalList):
    """Clear the atom maps of any atoms not involved in any of orbitalList.
    
    Assumes all same mol.  Uses orbitalList[0].mol
    """
    mol = orbitalList[0].mol
    atomSet = orbitalList[0].involvedAtomList()
    for orb in orbitalList[1:]:
        atomSet.update(orb.involvedAtomList())
    
    for atom in mol.GetAtoms():
        if atom not in atomSet:
            atom.SetMapIdx(0)

def oe_atomHybridization(atom):
    return OEGetHybridization(atom)

def atomHybridization( atom, visitedAtomIndexes=None ):
    """OEChem has OEGetHybridization method, and this works similarly.
    
    Return value is the p character of the hybridization, for example 2 for sp2
    
    Simple rule based now, that should work for atoms up to the second row at least.
    Easy to tell if this atom has more than one pi bond.  Otherwise, assume is sp3
    unless have an empty orbital then assume is sp2.  Or, if have a lone pair
    in an aromatic system, must be sp2.  In cases like oxygen of C=CO, is probably
    somewhere in between.
    
    Hydrogen is special case.
    
    For hyper-valent atoms, add 10 for each d-orbital involved.  
    Penta-coordinated atoms will be 14 (1 d + 4 p)
    Hexa-coordinated atoms will be 25 (2 d + 5 p)
    This doesn't make complete sense because there are only max 3 p orbitals, but is a 
    sufficient marker for the fact that they have non-standard coordination.
    """
    if visitedAtomIndexes is None:
        visitedAtomIndexes = set()
    visitedAtomIndexes.add( atom.GetIdx() )
    nValenceEle = stableValenceShell(atom.GetAtomicNum())

    if nValenceEle == 2:    # H, Li, Na, K, etc. that are 1 electron more than an empty outer shell
        return 0   # 1s orbital, no p character
    
    pOrbitals = 0
    nElectrons = 0
    
    orbInfo = orbitalInfo(atom)
    nElectrons += orbInfo['nUnbondedElectrons']
    pOrbitals += orbInfo["nEmptyOrbitals"]    # Empty orbital, assume must be p, not spX
    
    for bond in atom.GetBonds():
        pOrbitals += bond.GetOrder()-1
        nElectrons += 2 * bond.GetOrder()
    
    if orbInfo["nUnbondedElectrons"] > 0:
        # Potential to be p orbital lone pair.  Assume so if aromatic lone pair, otherwise
        #   stick to sp3 and let conjugated pi system perceiver deal with it.
        if atom.IsAromatic() and pOrbitals < 1: # Don't count if already have pi bonds
            pOrbitals += 1
    
    dOrbitals = 0
    if nMaxHyperValenceElectrons(atom.GetAtomicNum()) > 0:
        ## There is some possible hypervalancy going on.
        ## Check how many electrons are in the system
        if nElectrons > STD_OCTECT:
            dOrbitals = (nElectrons - STD_OCTECT)/2
    
    return 3 - pOrbitals + 10 * dOrbitals   # Default sp3. but discount any p-orbital only bonds / electrons
            

def findClosureAtoms(srcOrb, sinkOrb):
    """Given a reaction determine the closure atoms (which atoms are part of new/broken sigma bond)."""
    atomSet = set([])
    if srcOrb.atom.GetBond(sinkOrb.atom) is None:
        atomSet.add(srcOrb.atom)
        atomSet.add(sinkOrb.atom)
    ## Idea is to step through the orbitals.  If any are sigma bonds then, these need to be added.
    ## If no bond exists between adjacent chained orbs, then break them.
    ## 
    orb = srcOrb
    prevOrb = None
    while orb is not None:
        if orb.type == 'sigma':
            atomSet.add(orb.atom)
            atomSet.add(orb.neighbor)
        if prevOrb and prevOrb.neighbor and orb.atom.GetBond(prevOrb.neighbor) is None:
            atomSet.add(orb.atom)
            atomSet.add(prevOrb.neighbor)
        prevOrb = orb
        orb = orb.extOrbital
    
    orb = sinkOrb
    prevOrb = None
    while orb is not None:
        if orb.type == 'sigma':
            atomSet.add(orb.atom)
            atomSet.add(orb.neighbor)
        if prevOrb and prevOrb.neighbor and orb.atom.GetBond(prevOrb.neighbor) is None:
            atomSet.add(orb.atom)
            atomSet.add(prevOrb.neighbor)
        prevOrb = orb
        orb = orb.extOrbital

    #log.critical('Found atom Maps: %s' % pformat([a.GetMapIdx() for a in atomSet]))
    return list(atomSet)        
