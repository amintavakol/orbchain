import sys,os
from atom.utils import *

from chemutils.Common.Util import splitCompositeSmilesToList, molBySmiles
from chemutils.Common.Util import joinSmilesListToCompositeSmiles
from openeye.oechem import OEAssignAromaticFlags, OEAroModelMMFF, OEClearAromaticFlags
from openeye.oechem import OEPerceiveChiral, OEPerceiveSymmetry
from openeye.oechem import OEFindRingAtomsAndBonds
from chemutils.Common.Util import standardizeSmiles
from chemutils.Common.Util import createStdAtomMapSmiString, clearAtomMapsSmiStr
from chemutils.Common.CanonicalAtomMapSmiles import canonicalizeAtomMapSmiString, createCanonicalAtomMapSmiString
from chemutils.Common.MolExt import clearAtomMaps, removeNonsenseStereo
from chemutils.Common.MolExt import setSingleExplicitHydrogens, getAtmMappedSmilesFromCompositeSmiles


class SimpleAtomObject(object):
    """Simple object to capture a representation of an atom.
    
    Basically, given an atom object, store a copy of the atom object, full smiles, and the smiles of 
    only connected component containing self.
    """
    def __init__(self, oeAtom, oeMol, symmetricAtoms):
        """Set the basic params"""
        self.atom = oeAtom
        self.mol = oeMol
        
        self.symmetricAtoms = symmetricAtoms
        
        self.fullSmiles = None
        self.connectedSmiles = None
        self.connectedNonMappedSmiles = None
        
        clearAtomMaps(self.mol)
        self.atom.SetMapIdx(1)
        self.fullSmiles = createCanonicalAtomMapSmiString(self.mol)
        self.atom.SetMapIdx(0)
        self.connectedSmiles = getAtmMappedSmilesFromCompositeSmiles(self.fullSmiles, clearMaps=False)
        self.connectedSmiles = canonicalizeAtomMapSmiString(self.connectedSmiles)
        self.connectedNonMappedSmiles = clearAtomMapsSmiStr(self.connectedSmiles);
        self.connectedKekuleSmiles = canonicalKekule(self.connectedSmiles)
        
        # Eventual ML stuff
        self.rawFeatDict = None
        self.normFeatDict = None
        self.fill_pred_value = None
        self.unfill_pred_value = None
        self.fillPredDecision = None
        self.unfillPredDecision = None
    
    def __str__(self):
        """Just set this to return the connectedSmiles"""
        return self.connectedSmiles
    
    def __repr__(self):
        return "'%s'" % self.__str__();
    
    @staticmethod
    def atomObjFromReactantSmi(reactantSmiles):
        """Convenience method to make a list of SimpleAtomObjects"""
        """Given a single reactant (connected component), return the list of atoms."""
        # First, ensure a single copy!
        #reactantSmiles = canonicalKekule(reactantSmiles)
        
        reactantSmiles = clearAtomMapsSmiStr(reactantSmiles)
        mol = molBySmiles(reactantSmiles)
        OEAssignAromaticFlags(mol, OEAroModelMMFF)
        reactantSmiles = createCanonicalAtomMapSmiString(mol)
        smilesSet = set(splitCompositeSmilesToList(reactantSmiles))
        newReactantSmiles = joinSmilesListToCompositeSmiles(list(smilesSet))
        newReactantSmiles = canonicalizeAtomMapSmiString(newReactantSmiles)
        
        mol = molBySmiles(newReactantSmiles);
        setSingleExplicitHydrogens(mol)
        OEAssignAromaticFlags(mol, OEAroModelMMFF)
        clearAtomMaps(mol)
        OEPerceiveChiral(mol)
        removeNonsenseStereo(mol)
        OEPerceiveSymmetry(mol)
        OEFindRingAtomsAndBonds(mol)
        
        resBySymmetryClass = {}
        seenSymClassSet = set([])
        for atm in mol.GetAtoms():
            theSymClass = atm.GetSymmetryClass()
            if theSymClass not in resBySymmetryClass:
                resBySymmetryClass[theSymClass] = SimpleAtomObject(atm, mol, [atm])
            else:
                resBySymmetryClass[theSymClass].symmetricAtoms.append(atm)
        return resBySymmetryClass.values()
    
    
    @staticmethod
    def combineAtomObjListToSmiles(atomObjList):
        """Given some 'filtered' atomObjs in a list, return a composite smiles with all of these atoms labeled."""
        if len(atomObjList) == 0:
            return None
        mol = atomObjList[0].mol
        clearAtomMaps(mol)
        for atomObj in atomObjList:
            [atm.SetMapIdx(1) for atm in atomObj.symmetricAtoms]
            #atomObj.atom.SetMapIdx(1)
        smi = createStdAtomMapSmiString(mol)
        clearAtomMaps(mol)
        return smi

    @staticmethod
    def combineAtomObjListToSmiles_diff_indices(atomObjList):
        """same as combineAtomObjListToSmiles but assign differenet numbers to each simpleAtomObj"""
        if len(atomObjList) == 0:
            return None
        mol = atomObjList[0].mol
        clearAtomMaps(mol)
        i = 1
        for atomObj in atomObjList:
            [atm.SetMapIdx(i) for atm in atomObj.symmetricAtoms]
            #atomObj.atom.SetMapIdx(1)
            i+=1
        smi = createStdAtomMapSmiString(mol)
        clearAtomMaps(mol)
        return smi
