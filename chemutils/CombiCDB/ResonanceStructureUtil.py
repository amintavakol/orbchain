"""Extensions for molecule object descriptors to support enumerating the different
resonance structure forms for a molecule.

Because of the exponential worst case of enumerating all the resonance structures 
in mols with a number of different pi systems, we work with a slightly odd 
representation.  The standard resonance structure is a set of sets of major 
resonance of each pi system fragment.  
Major resonance structures are either neutral, aromatic, or try to apply formal
negative charge to the most electronegative atoms.

To check if two molecules are resonance structs of each other.  You should first check
the molSkeleton (chemutils.Common.MolExt).  Then calc the standard resonance structure.  
Only if both are the same are the mols equivalent resonance structs.

Another timesaving technique used is to NOT try to move around electrons if we find a 
neutral and non-charge separated pi system.  If this is aromatic, flags will be set
to canonicalize (so we don't have to worry about equiv res structs with diff kekule forms)


Major functions include
    resonanceStructureIter:
        Rearrange pi system orbitals in molecule in sensible ways to show
        the different possible resonance structure representations of the molecule.
    
    standardizeResonanceStructure:
        Given an input molecule, rearrange the pi system orbitals until
        an "optimal" configuration is found, preferably one that is most representative
        of the true structure.  More important that this produces a consistent
        result given equivalent inputs (essentially looking for a canonical format).
        Instead of returning *a* standard molecule, returns a set of smiles string that depict
        equivalent major resonance structures of the compound.
    
    equivalentResonanceStructures:
        Given two input molecules, determine if they represent equivalent
        resonance structures of the same "true" structure.
        Could just run standardizeResonanceStructure on each and compare
        results, but some simple steps up front could rule out equivalency first.
    
    fragmentPiSystems:
        Given an input mol, return a list of mol objects that represent the disconnected 
        pi systems of the input mol.
    

"""

from openeye.oechem import OECreateIsoSmiString, OEGraphMol
from openeye.oechem import OEClearAromaticFlags, OEAddExplicitHydrogens
from openeye.oechem import OEMolecularFormula, OENetCharge, OESubsetMol
from openeye.oechem import OEAssignAromaticFlags

from chemutils.Common.MolStdValue import atomElectronegativity
from chemutils.Common.Util import (
    molBySmiles,
    createStandardSmiString,
    splitCompositeMolToSmilesList,
)
from chemutils.Common.Util import splitCompositeMol
from chemutils.Common.Util import splitCompositeSmilesToList
from chemutils.Common.MolExt import grossFormalCharges, createStandardMol, molSkeleton
from chemutils.Common.MolExt import atomIsInPiSystem, wrapAsOEUnaryAtomPred
from chemutils.Common.OrbitalModel import orbitalIter, orbitalInfo
from chemutils.CombiCDB.OrbitalInteraction import (
    moveOrbitalElectrons,
    undoMoveOrbitalElectrons,
)
from chemutils.CombiCDB.OrbitalInteraction import compatibleOrbitalPair

from chemutils.CombiCDB.Util import log

# from chemutils.score.molecule.FreeEnergyNoResonanceScore import FreeEnergyNoResonance;

# Adjust the molSkeleton function to be sure to retain all the singleton Hs
# def molSkeleton(mol):
#    """Wrapper around origMolSkeleton to ensure retainSingleHs is set True"""
#    return origMolSkeleton(mol, retainSingleHs=True);

# Set a maxDepth of 3.  This should be enough for an overlap of things up to napthalene level
MAX_DEPTH = 3


# Atomic number of carbon, the basis for organic chemistry
#   whose properties define borderlines between electronegative vs. electropositive elements, etc.
CARBON = 6
carbonEN = atomElectronegativity(CARBON)


def fragmentPiSystems(mol):
    """Return a tuple with the disconnected pi systems, and the disconnected NON-pi systems

    Using OESubsetMol and atomIsInPiSystem
    """
    # First make a OEUnaryAtomPredObject that will check for inclusion in a pi system.
    IsInPiClass = wrapAsOEUnaryAtomPred(atomIsInPiSystem)
    IsNotInPiClass = wrapAsOEUnaryAtomPred(atomIsInPiSystem, doNot=True)
    isInPiObj = IsInPiClass()
    isNotInPiObj = IsNotInPiClass()

    # Then can subset the mol
    # mol = createStandardMol(mol)
    newMol = OEGraphMol()
    OESubsetMol(newMol, mol, isInPiObj, True)

    # Make a set to have just non-redundant pi systems
    piSmiSet = set(splitCompositeMolToSmilesList(newMol))
    if piSmiSet == set([""]):
        piSmiList = []
    else:
        piSmiList = [smi for smi in piSmiSet]

    # Do the non-pi
    newMol = OEGraphMol()
    OESubsetMol(newMol, mol, isNotInPiObj, True)

    # Make a set to have just non-redundant non-pi systems (but leave this as a set)
    nonPiSmiSet = frozenset(splitCompositeMolToSmilesList(newMol))
    if nonPiSmiSet == frozenset([""]):
        nonPiSmiSet = frozenset([])

    return (piSmiList, nonPiSmiSet)


def standardizeResonanceStructure(mol, seenResStructDict={}):
    """Given an input molecule, find a canonical representation of the major resonance structures

    The standardization involves a tuple of the (molSkeletonSmi, [set(resStructures)])
    where the list is over each disparate pi system of the molecule, and the sets are over
    the major resonance contributors

    Note:
    """
    # First standardize to remove any spurious stereo values
    mol = createStandardMol(mol)

    # Then calc the mol skeleton:
    molSkel = molSkeleton(mol)
    molSkelSmi = createStandardSmiString(molSkel)

    # Then break the mol into separate pi systems

    piSystemSmiList, nonPiSmiSet = fragmentPiSystems(mol)

    repResStructSetList = []
    for componentSmi in piSystemSmiList:
        if componentSmi in seenResStructDict:
            repResStructSet = seenResStructDict[componentSmi]
        else:
            componentMol = molBySmiles(componentSmi)
            repResStructSet = set()
            # First check do we have any kind of charge separation
            grossPosCharge, grossNegCharge = grossFormalCharges(componentMol)
            if abs(grossPosCharge) + abs(grossNegCharge) > 0:
                # Have to enumerate res structs
                for structure in resonanceStructureIter(
                    componentMol,
                    includeMinor=False,
                    breakNeutralMol=False,
                    breakByDepth=True,
                ):
                    copyStruct = OEGraphMol(structure)
                    OEAssignAromaticFlags(copyStruct)
                    structureSmi = createStandardSmiString(copyStruct)

                    # Before adding, check if we see the charge separation go away.
                    # If so, this should be the only thing in the set and we can stop looking
                    grossPosCharge, grossNegCharge = grossFormalCharges(copyStruct)
                    if abs(grossPosCharge) + abs(grossNegCharge) == 0:
                        repResStructSet = set([structureSmi])
                        break
                    else:
                        repResStructSet.add(structureSmi)
            else:
                # Then simply add in this smi...
                OEAssignAromaticFlags(componentMol)
                repResStructSet.add(createStandardSmiString(componentMol))
        repResStructSetList.append(repResStructSet)
        seenResStructDict[componentSmi] = repResStructSet

    # print stdRep
    return (molSkelSmi, repResStructSetList, nonPiSmiSet)


def standardizeResonanceStructureSmi(smi, seenResStructDict={}):
    """Given an input smi, return a set of smiles with all the major res structs
    Wrapper around standardizeResonanceStructure.

    >>> test = ("[H-].C.C.Br.[Br+]","[H+].C.C.Br.[Br-]")
    >>> standardizeResonanceStructureSmi(test[0]);
    ('[C].[C].[Br].[Br]', [], frozenset(['C', '[H-]', '[Br+]', 'Br']))
    >>> standardizeResonanceStructureSmi(test[1]);
    ('[C].[C].[Br].[Br]', [], frozenset(['[H+]', 'C', '[Br-]', 'Br']))
    >>> test = ("[C-]#[O+]","[C]=O")
    >>> standardizeResonanceStructureSmi(test[0]);
    ('[C][O]', [set(['[C]=O'])], frozenset([]))
    >>> standardizeResonanceStructureSmi(test[1]);
    ('[C][O]', [set(['[C]=O'])], frozenset([]))
    >>> test = ('[N-]1[CH+]C=CC=C1','[N-]1C=C[CH+]C=C1', )
    >>> standardizeResonanceStructureSmi(test[0])
    ('[C]1[C][C][N][C][C]1', [set(['c1ccncc1'])], frozenset([]))
    >>> standardizeResonanceStructureSmi(test[1])
    ('[C]1[C][C][N][C][C]1', [set(['c1ccncc1'])], frozenset([]))
    >>> standardizeResonanceStructureSmi('CCO[C@@](C)([O-])Cl')
    ('[C][C]OC([C])([O])Cl', [], frozenset(['CCO[C@@](C)([O-])Cl']))
    >>> standardizeResonanceStructureSmi('CC[O@@H+][C@@H](C)C')
    ('[C][C]O[C]([C])[C]', [], frozenset(['CC[OH+]C(C)C']))
    >>> standardizeResonanceStructureSmi('CC[OH+][C@@](=O)C') # Should remove spurious stereo
    ('[C][C]O[C]([C])[O]', [set(['C(=O)[OH2+]'])], frozenset(['CC', 'C']))
    >>> standardizeResonanceStructureSmi('CC[O@@H+][C@@](=O)C') # this should remove the spurious onium stereo
    ('[C][C]O[C]([C])[O]', [set(['C(=O)[OH2+]'])], frozenset(['CC', 'C']))
    >>> standardizeResonanceStructureSmi('CCO.C=C(c1ccccc1Cl)[O-].[Na+]')
    ('[C][C][O].[C][C]([C]1[C][C][C][C][C]1Cl)[O].[Na]', [set(['[CH2-]C(=O)c1ccccc1Cl', 'C=C(c1ccccc1Cl)[O-]'])], frozenset(['CCO', '[Na+]']))
    >>> standardizeResonanceStructureSmi('CC(=O)C1=C(C=CO1)Br.Br')
    ('[C][C]([C]1[C]([C][C]O1)Br)[O].[Br]', [set(['c1coc(c1Br)C=O'])], frozenset(['C', 'Br']))
    >>> standardizeResonanceStructureSmi('C(CO)C(C#N)C#N.O')
    ('[C]([C][O])[C]([C][N])[C][N].[O]', [set(['C#N'])], frozenset(['CCCO', 'O']))
    >>> standardizeResonanceStructureSmi('c1ccc(cc1)[N+]#N')
    ('[C]1[C][C][C]([C][C]1)[N][N]', [set(['c1ccc(cc1)[N+]#N'])], frozenset([]))
    >>> standardizeResonanceStructureSmi('c1ccc(cc1)N=N')
    ('[C]1[C][C][C]([C][C]1)[N][N]', [set(['c1ccc(cc1)N=N'])], frozenset([]))
    >>> standardizeResonanceStructureSmi('C=C[C+](C=C)N=NC1=C(C=CC=C1)O')
    ('[C][C][C]([C][C])[N][N][C]1[C][C][C][C][C]1[O]', [set(['C=C[C+](C=C)N=Nc1ccccc1O', 'C=CC(=C[CH2+])N=Nc1ccccc1O', 'C=CC(=[N+]=Nc1ccccc1O)C=C'])], frozenset([]))
    >>> standardizeResonanceStructureSmi('[CH-]1C=CC=NC1N')
    ('[C]1[C][C][N][C]([C]1)[N]', [set(['[CH2-]C=CC=N', 'C=C[CH-]C=N', 'C=CC=C[NH-]'])], frozenset(['CN']))
    >>> standardizeResonanceStructureSmi('C1=CC([N-]C=C1)N')
    ('[C]1[C][C][N][C]([C]1)[N]', [set(['[CH2-]C=CC=N', 'C=C[CH-]C=N', 'C=CC=C[NH-]'])], frozenset(['CN']))
    """
    mol = molBySmiles(smi)
    mol = createStandardMol(mol)
    return standardizeResonanceStructure(mol, seenResStructDict)


def combineSmiles(componentsSetsList, combinedSoFar, lengthSoFar=0):
    """From a list of component *sets* of SMILES, return a set of as many corresponding
    SMILES strings composed of one SMILES from each set.
    Recursive.
    """
    if lengthSoFar == len(componentsSetsList):
        combinedStrings = set()
        for combined in combinedSoFar:
            combined.sort()  # alphabetic sort ensures canonization
            combinedStrings.add(".".join(combined))
        return combinedStrings

    if lengthSoFar == 0:
        # initialize combinedSoFar
        newCombinedSoFar = []
        for representer in componentsSetsList[lengthSoFar]:
            newCombinedSoFar.append([representer])
    else:
        # update all elements from combinedSoFar by joining them a new compound
        newCombinedSoFar = []
        for representer in combinedSoFar:
            for compoundToAdd in componentsSetsList[lengthSoFar]:
                representerCompounds = representer.append(compoundToAdd)
                newCombinedSoFar.append(representer)
    return combineSmiles(componentsSetsList, newCombinedSoFar, lengthSoFar + 1)


def assignAtomMaps(mol, undoDict=None):
    """
    Assign atom idx to every atom of the given molecule that does not already have an idx.
    Store in undoDict whatever values already exist, to reset via undoAssignAtomMaps.
    """
    # Add hydrogens for atom mapping, if they are already explicit
    # OEAssignImplicitHydrogens(mol)
    OEAddExplicitHydrogens(mol)

    usedIndexes = []
    if undoDict is not None:
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() != 0:
                usedIndexes.append(atom.GetMapIdx())
                undoDict[atom.GetIdx()] = atom.GetMapIdx()

    if usedIndexes == []:
        idx = 1
    else:
        idx = max(usedIndexes) + 1
    for atom in mol.GetAtoms():
        if atom.GetMapIdx() == 0:
            atom.SetMapIdx(idx)
            idx += 1
        # print "xxx atom %s with label %s" % (atom.GetAtomicNum(), atom.GetMapIdx())


def undoAssignAtomMaps(mol, undoDict=None):
    """
    Erase atom idx of every atom of the given molecule except for those in undoDict.
    """
    for atom in mol.GetAtoms():
        if undoDict is not None:
            if not undoDict.has_key(atom.GetMapIdx()):
                atom.SetMapIdx(0)
        else:
            atom.SetMapIdx(0)


def equivalentResonanceStructuresSmi(smi1, smi2, seenResStructDict={}):
    """Given two input smiles strings, determine if they represent equiv resonance structures

    Basically a wrapper around the equivalentResonanceStructures function.

    >>> test = ("[C-]#[O+]","[C]=O")
    >>> equivalentResonanceStructuresSmi(test[0], test[1]);
    True
    >>> test = ('[N-]1[CH+]C=CC=C1','[N-]1C=C[CH+]C=C1', "N1=CC=CC=C1")
    >>> equivalentResonanceStructuresSmi(test[0], test[1])
    True
    >>> equivalentResonanceStructuresSmi(test[0], test[2])
    True
    >>> equivalentResonanceStructuresSmi(test[1], test[2])
    True
    >>> testNotEquiv = ('C[N+](=O)[O-]', 'N1=CC=CC=C1', 'CC(=O)O')
    >>> equivalentResonanceStructuresSmi(testNotEquiv[0], testNotEquiv[1])
    False
    >>> equivalentResonanceStructuresSmi(testNotEquiv[0], testNotEquiv[2])
    False
    >>> equivalentResonanceStructuresSmi(testNotEquiv[1], testNotEquiv[2])
    False
    >>> testCloseButNotEquiv = ('CC(=O)[O-]', '[CH2-]C(=O)O')
    >>> equivalentResonanceStructuresSmi(testCloseButNotEquiv[0], testCloseButNotEquiv[1])
    False
    >>> testCloseButNotEquiv = ("[H-].C.C.Br.[Br+]","[H+].C.C.Br.[Br-]")
    >>> equivalentResonanceStructuresSmi(testCloseButNotEquiv[0], testCloseButNotEquiv[1])
    False
    >>> test = ("CCC[CH+][O-]","CCCC=O")
    >>> equivalentResonanceStructuresSmi(test[0], test[1])
    True
    >>> test = ("[CH2+]CC[CH+][O-]","CCCC=O")
    >>> equivalentResonanceStructuresSmi(test[0], test[1])
    False
    >>> test = ("[CH2-]CC[CH+][O-]","CCCC=O")
    >>> equivalentResonanceStructuresSmi(test[0], test[1])
    False
    >>> test = ("[CH2-]CCC=O","CCCC=O")
    >>> equivalentResonanceStructuresSmi(test[0], test[1])
    False
    >>> # These might cause a runtime exception
    >>> test =['CCOC(CC(=O)Cc1ccccc1)([O-])OC(=O)c2ccccc2', 'CCOC(C[C+]([O-])Cc1ccccc1)([O-])OC(=O)c2ccccc2']
    >>> equivalentResonanceStructuresSmi(test[0], test[1]);
    True
    >>> test =['CCOC(CC(=O)Cc1ccccc1)([O-])OC(=O)c2ccccc2', 'CCOC(C[C+]([O-])Cc1ccccc1)([O-])O[C+]([O-])c2ccccc2']
    >>> equivalentResonanceStructuresSmi(test[0], test[1]);
    True
    >>> test = ['[CH-]1C=CC=NC1N.C1=CC([N-]C=C1)N', 'C1=CC([N-]C=C1)N.C1=CC([N-]C=C1)N']
    >>> equivalentResonanceStructuresSmi(test[0], test[1]);
    True

    """
    if smi1 == smi2:
        return True
    mol1 = molBySmiles(smi1)
    mol2 = molBySmiles(smi2)
    return equivalentResonanceStructures(mol1, mol2, seenResStructDict)


def equivalentResonanceStructures(mol1, mol2, seenResStructDict={}):
    """Given two input molecules, determine if they represent equivalent
    resonance structures of the same "true" structure.
    Could just run standardizeResonanceStructure on each and compare
    results, but some simple steps up front could rule out equivalency faster.
    """
    # Check mol formula
    molecularFormula1 = OEMolecularFormula(mol1)
    molecularFormula2 = OEMolecularFormula(mol2)
    if molecularFormula1 != molecularFormula2:
        # print "Different molecular formulas"
        return False

    # Check net charge
    if OENetCharge(mol1) != OENetCharge(mol2):
        # print "Different net charges"
        return False

    # Check mol skeleton (chemutils.Common.MolExt.molSkeleton)
    # Copy the mol objects first
    copyMol1 = OEGraphMol(mol1)
    copyMol2 = OEGraphMol(mol2)
    skeletonSmi1 = OECreateIsoSmiString(molSkeleton(copyMol1))
    skeletonSmi2 = OECreateIsoSmiString(molSkeleton(copyMol2))
    if skeletonSmi1 != skeletonSmi2:
        # print "Different molecular skeletons"
        return False

    # Check standardized structures
    stdRep1 = standardizeResonanceStructure(mol1, seenResStructDict)
    stdRep2 = standardizeResonanceStructure(mol2, seenResStructDict)

    return equivalentStandardResonanceRepresentations(stdRep1, stdRep2)


def equivalentResSmiStandardRep(
    smi, stdRep, returnResStructForSmi=False, seenResStructDict={}
):
    """Given a smiles and a stdRep of a resonance structure, return if they're equivalent.

    Can do some simple things by doing a mol skel check first before calling standardizeResonanceStructureSmi.

    >>> test = ("[C-]#[O+]","[C]=O")
    >>> stdRep = standardizeResonanceStructureSmi(test[1]);
    >>> equivalentResSmiStandardRep(test[0], stdRep);
    True
    >>> equivalentResSmiStandardRep(test[1], stdRep);
    True
    >>> test = ('[N-]1[CH+]C=CC=C1','[N-]1C=C[CH+]C=C1', )
    >>> stdRep = standardizeResonanceStructureSmi(test[0])
    >>> equivalentResSmiStandardRep(test[0], stdRep)
    True
    >>> equivalentResSmiStandardRep(test[1], stdRep)
    True
    >>> testNotEquiv = ('C[N+](=O)[O-]', 'N1=CC=CC=C1', 'CC(=O)O')
    >>> equivalentResSmiStandardRep(testNotEquiv[0], standardizeResonanceStructureSmi(testNotEquiv[1]))
    False
    >>> equivalentResSmiStandardRep(testNotEquiv[0], standardizeResonanceStructureSmi(testNotEquiv[2]))
    False
    >>> equivalentResSmiStandardRep(testNotEquiv[1], standardizeResonanceStructureSmi(testNotEquiv[2]))
    False
    >>> testCloseButNotEquiv = ('CC(=O)[O-]', '[CH2-]C(=O)O')
    >>> equivalentResSmiStandardRep(testCloseButNotEquiv[0], standardizeResonanceStructureSmi(testCloseButNotEquiv[1]))
    False
    >>> equivalentResSmiStandardRep(testCloseButNotEquiv[1], standardizeResonanceStructureSmi(testCloseButNotEquiv[0]))
    False
    """
    mol = molBySmiles(smi)
    skeletonSmi1 = OECreateIsoSmiString(molSkeleton(mol))
    skeletonSmi2 = stdRep[0]
    if skeletonSmi1 != skeletonSmi2:
        if returnResStructForSmi:
            return False, None
        else:
            return False

    newStdRep = standardizeResonanceStructure(molBySmiles(smi), seenResStructDict)
    if returnResStructForSmi:
        return (
            equivalentStandardResonanceRepresentations(newStdRep, stdRep),
            newStdRep,
        )
    else:
        return equivalentStandardResonanceRepresentations(newStdRep, stdRep)


#
def isSubsetResonanceStructuresSmi(smi1, smi2, seenResStructDict={}):
    """Given two smiles, determine if smi1 is completely a subset of the other.

    >>> test = ("[C-]#[O+]","[C]=O")
    >>> isSubsetResonanceStructuresSmi(test[0], test[1]);
    True
    >>> test = ("[C-]#[O+]","[C]=O.CCCO")
    >>> isSubsetResonanceStructuresSmi(test[0], test[1]);
    True
    >>> test = ("[C-]#[O+].CCCO","[C]=O")
    >>> isSubsetResonanceStructuresSmi(test[0], test[1]);
    False
    """
    if smi1 == smi2:
        return True
    mol1 = molBySmiles(smi1)
    mol2 = molBySmiles(smi2)
    return isSubsetResonanceStructureMol(mol1, mol2, seenResStructDict)


def isSubsetResonanceStructureMol(mol1, mol2, seenResStructDict={}):
    """Given two mols, is mol1 a subset of resonance structure representations of mol2.

    Can basically call the isSubsetResonanceRepresentation but can do a simple test before hand as well.
    """
    # Check mol skeleton (chemutils.Common.MolExt.molSkeleton)
    # Copy the mol objects first
    copyMol1 = OEGraphMol(mol1)
    copyMol2 = OEGraphMol(mol2)
    skeletonSmi1 = set(
        splitCompositeSmilesToList(OECreateIsoSmiString(molSkeleton(copyMol1)))
    )
    skeletonSmi2 = set(
        splitCompositeSmilesToList(OECreateIsoSmiString(molSkeleton(copyMol2)))
    )
    if not skeletonSmi1.issubset(skeletonSmi2):
        # print "Different molecular skeletons"
        return False

    # Check standardized structures
    stdRep1 = standardizeResonanceStructure(mol1, seenResStructDict)
    stdRep2 = standardizeResonanceStructure(mol2, seenResStructDict)

    return isSubsetResonanceRepresentation(stdRep1, stdRep2)


def isSubsetResonanceRepresentation(stdRepA, stdRepB):
    """Test that the stdRep is a subset of stdRepB.

    The idea being I have some reactants (stdRepB) and then some new products from a reaction, where
    I have removed any unchanged reactants.  Then I want to see that the set of molSkeletonA disconnected
    components is a subset of the molSkeletonB disconnected components, and subset of the nonPi and the
    rearrangedPi Sets.
    """
    molSkelSmiA = set(splitCompositeSmilesToList(stdRepA[0]))
    molSkelSmiB = set(splitCompositeSmilesToList(stdRepB[0]))
    if not molSkelSmiA.issubset(molSkelSmiB):
        return False

    nonPiSetA = stdRepA[2]
    nonPiSetB = stdRepB[2]
    if not nonPiSetA.issubset(nonPiSetB):
        return False

    #
    molAResSetList = stdRepA[1]
    molBResSetList = stdRepB[1]

    intersectArr = []
    for aSet in molAResSetList:
        subIntList = [len(aSet.intersection(bSet)) > 0 for bSet in molBResSetList]
        intersectArr.append(subIntList)

    allACovered = allBoolean([anyBoolean(subList) for subList in intersectArr])
    # allBCovered = allBoolean([anyBoolean([intersectArr[i][j] for i in range(len(intersectArr))]) for j in range(len(molBResSetList))])

    return allACovered


def equivalentStandardResonanceRepresentations(stdRepA, stdRepB):
    """Helper function to encapsulating testing the equivalence of two standard representations.

    The representation is (molSkelSmi, [set(resStructs)]) as described in the standardizeResonanceStructure
    docstrings.  Test for equality here.
    """
    molSkelSmiA = stdRepA[0]
    molSkelSmiB = stdRepB[0]

    if molSkelSmiB != molSkelSmiA:
        return False

    nonPiSetA = stdRepA[2]
    nonPiSetB = stdRepB[2]
    if nonPiSetA != nonPiSetB:
        return False

    molAResSetList = stdRepA[1]
    molBResSetList = stdRepB[1]

    intersectArr = []
    for aSet in molAResSetList:
        subIntList = [len(aSet.intersection(bSet)) > 0 for bSet in molBResSetList]
        intersectArr.append(subIntList)

    allACovered = allBoolean([anyBoolean(subList) for subList in intersectArr])
    allBCovered = allBoolean(
        [
            anyBoolean([intersectArr[i][j] for i in range(len(intersectArr))])
            for j in range(len(molBResSetList))
        ]
    )

    return allACovered and allBCovered


def resonanceStructureIter(
    mol,
    includeMinor=False,
    observedSmilesSet=None,
    depth=0,
    arrowList=None,
    breakNeutralMol=True,
    breakByDepth=False,
):
    """For the given molecule or intermediate, enumerate all of the reasonable
    resonance structures by shifting electrons around available pi system atomic
    and molecular orbitals.  For efficiency, will always be the same molecule object
    yielded, just modified.  Caller should make sure that it comes back in the same state.

    includeMinor:
        Option whether or not to return minor resonance structures like
        charge separated versions of neutral molecules (e.g., [CH2+][O-])

    Strategy: Look for any 2 adjacent atoms hybridized with p or pi orbitals
    and rearrange electrons in a chemically natural manner.

    When includeMinor is False:
    Filter to only consider 2 adjacent atoms where at least one of them has a formal
    charge.  This will not create *new* charge separation, but will allow charges to
    *move* around.

    After any case, recursively look for more possibilities.  Use observedSmilesSet
    to avoid recursion loops.

    returnArrow:
        Option whether or not to return a tuple (mol, arrowList) where arrowList is
        the transformation that led to this particular resonance structure or not.

    """
    if observedSmilesSet is None:
        # OEClearAromaticFlags(mol);
        observedSmilesSet = set()

    # OEAssignAromaticFlags(mol)
    smiles = OECreateIsoSmiString(mol)
    if smiles in observedSmilesSet:
        # Base case, already tried this resonance structure
        # print "Already tried that."
        # log.info('Res Iter returning with back to old smi');
        return
    elif depth >= MAX_DEPTH and breakByDepth:
        # log.info('Res Iter returning at depth: %d' % depth);
        return
    else:
        observedSmilesSet.add(smiles)

    if screenResonanceStruct(mol, includeMinor):
        # print "Keeping: ",  smiles
        if arrowList is not None:
            yield (mol, arrowList)
        else:
            yield mol

    # log.critical('At depth %d, with mol: %s' % (depth, smiles));

    grossPos, grossNeg = grossFormalCharges(mol)
    # Only allow a completely neutral mol to break up all pi bonds.
    completelyNeutralMol = (
        grossNeg == 0 and grossPos == 0 and depth == 0 and breakNeutralMol
    )

    piPairs = []
    # Check every pair of bonded atoms (filter so at least one has a formal charge)
    for bond in mol.GetBonds():
        bgnAtm = bond.GetBgn()
        endAtm = bond.GetEnd()

        # Only explore all the bonds if we meet one of these conditions
        keepBond = (
            completelyNeutralMol
            or includeMinor
            or (bgnAtm.GetFormalCharge() != 0 or endAtm.GetFormalCharge() != 0)
        )

        if not keepBond:
            continue
        else:
            # log.critical('At depth %d, and keeping a bond' % depth)
            pass

        # print "Begin", bond.GetBgn().GetAtomicNum()
        beginPOrbs = []
        for beginOrb in orbitalIter(bond.GetBgn()):
            if beginOrb.inPiSystem():
                beginPOrbs.append(beginOrb)
                # print str(beginOrb)

        # print "End", bond.GetEnd().GetAtomicNum()
        endPOrbs = []
        for endOrb in orbitalIter(bond.GetEnd()):
            if endOrb.inPiSystem():
                endPOrbs.append(endOrb)
                # print str(endOrb)

        if len(beginPOrbs) > 0 and len(endPOrbs) > 0:
            piPairs.append([beginPOrbs, endPOrbs])

    for [beginPOrbs, endPOrbs] in piPairs:
        for resonanceStructure, newArrowList in rearrangeOrbitalElectrons(
            mol, bond, beginPOrbs, endPOrbs, includeMinor
        ):
            if arrowList is not None:
                newArrowList.extend(arrowList)
                # Recursive call
                for recursiveStruct, newerArrowList in resonanceStructureIter(
                    resonanceStructure,
                    includeMinor,
                    observedSmilesSet,
                    depth + 1,
                    newArrowList,
                ):
                    yield (recursiveStruct, newerArrowList)
            else:
                for recursiveStruct in resonanceStructureIter(
                    resonanceStructure,
                    includeMinor,
                    observedSmilesSet,
                    depth + 1,
                    arrowList,
                    breakByDepth=breakByDepth,
                ):
                    yield recursiveStruct

        # Try the inverse direction too, since order shouldn't matter
        for resonanceStructure, newArrowList in rearrangeOrbitalElectrons(
            mol, bond, endPOrbs, beginPOrbs, includeMinor
        ):
            if arrowList is not None:
                newArrowList.extend(arrowList)
                # Recursive call
                for recursiveStruct, newerArrowList in resonanceStructureIter(
                    resonanceStructure,
                    includeMinor,
                    observedSmilesSet,
                    depth + 1,
                    newArrowList,
                    breakByDepth=breakByDepth,
                ):
                    yield (recursiveStruct, newerArrowList)
            else:
                # Recursive call
                for recursiveStruct in resonanceStructureIter(
                    resonanceStructure,
                    includeMinor,
                    observedSmilesSet,
                    depth + 1,
                    arrowList,
                    breakByDepth=breakByDepth,
                ):
                    yield recursiveStruct


def uniqueAtomResonanceStructureIter(mol, includeMinor=False, undoLabels=True):
    """Simple wrapper around resonance structureIter method.  Before generating
    the resonance structures however, will give a unique atom map label to
    every atom in the molecule with the effect that symmetric resonance structures
    will each be counted individually, instead of being skipped over as if identical.

    undoLabels: Specifies whether to clear out these fake labels before
                    each result is yielded
    """
    undoDict = dict()
    labelUniqueAtomMaps(mol, undoDict)

    for result in resonanceStructureIter(mol, includeMinor):
        if undoLabels:
            undoLabelUniqueAtomMaps(mol, undoDict)

        yield result

        if undoLabels:
            # Relabel to state before finding more resonance structs
            labelUniqueAtomMaps(mol, undoDict)

    if undoLabels:
        undoLabelUniqueAtomMaps(mol, undoDict)


def labelUniqueAtomMaps(mol, undoDict=None):
    """Helper to label the atom map attribute of all of the atoms of the molecule
    equal to their index position in the molecule, which ensures they will all
    receive a unique value.

    If undoDict provided, then store the original values of the atom map indexes
    so they can later be restored by the undoLabelUniqueAtomMaps function.
    """
    # Add hydrogens
    OEAddExplicitHydrogens(mol)

    for atom in mol.GetAtoms():
        if undoDict is not None:
            undoDict[atom.GetIdx()] = atom.GetMapIdx()
        atom.SetMapIdx(atom.GetIdx())


def undoLabelUniqueAtomMaps(mol, undoDict=None):
    """Undo the effects of the labelUniqueAtomMaps function.  If no specific undoDict
    is provided, then just reset all of the atom map indexes to "0" indicating standard values.
    Otherwise, use the specific values stored in the undoDict.
    """
    origMapIdx = 0
    for atom in mol.GetAtoms():
        if undoDict is not None:
            origMapIdx = undoDict[atom.GetIdx()]
        atom.SetMapIdx(origMapIdx)


def rearrangeOrbitalElectrons(mol, bond, beginPOrbs, endPOrbs, includeMinor=False):
    """Given just the p and pi orbitals for the given bond,
    enumerate each natural way of rearranging the electrons within them.
    Don't have to consider mirror cases, leave that as caller's responsibility
    by calling this function twice, but swapping the begin and end atoms and pOrbs.

    Cases of interest (not counting radicals for now):
        - Both atoms share a pi molecular orbital.  That is, the bond includes a pi bond.
            . Separate the pi molecular orbitals into one empty orbital (+ charge) and
                one lone pair orbital (- charge).
                Then caller can recurse on this intermediate which should fit the above rules.
        - One atom has lone pair, other has empty p orbital (artificial case representing charge separation modeling)
            . Move lone pair to empty p orbital to form new pi bond here
        - One atom has pi molecular orbital other has empty p orbital
            . Move the pi electrons to form a new pi bond here
        - One atom has a lone pair in p orbital, other has a pi molecular orbital
            . Move the lone pair into a new pi bond
            . Move the pi electrons to a lone pair over the far (3rd) atom

    Additional case where atoms each have a pi molecular orbital, but not the same one,
    like the central bond in C=C-C=C.  Do not necessarily have to do
    anything here.  Let another case separate
    the pi orbital of one of the adjacent double bonds, and then it will look like
    a case we can work with (lone pair or empty orbital next to pi orbital).

    Note that this keeps yielding the same (but modified) molecule object.
    Thus, in implementing, after each yield usage, the atom and bond states must
    be returned to how they were.

    Can almost just use the moveOrbitalElectrons functions from every beginPOrb
    to every endPOrb and let that sort out the details.  This would work, except
    that we need to undo the changes at each step to revert the molecule to
    its original state, which is hard to do if we don't know the specific steps.

    Also returns the arrowList leading to the corresponding rearrangement.

    Note:  Before yielding, some check are done to cull out structures that don't lead to
    good resonance structures.  Calls exploreResonanceStructure:
    """
    # Keep track of whether either end of the bond includes empty orbitals or non-bonded electrons
    nonBondedBeginOrbs = False
    emptyBeginOrbs = False
    emptyEndOrbs = False

    for beginOrb in beginPOrbs:
        if beginOrb.isEmpty():
            emptyBeginOrbs = True
        elif not beginOrb.isBondOrbital():
            nonBondedBeginOrbs = True
    for endOrb in endPOrbs:
        if endOrb.isEmpty():
            emptyEndOrbs = True
            break
        elif endOrb.isBondOrbital():
            # Check if the other end is has an empty orbital
            neighborInfo = orbitalInfo(endOrb.neighbor)
            if neighborInfo["nEmptyOrbitals"] > 0:
                emptyEndOrbs = True
                break

    # print OECreateIsoSmiString(mol), nonBondedBeginOrbs, emptyBeginOrbs, emptyEndOrbs

    if not emptyBeginOrbs:
        # If empty orbitals exist on source atom, no combination will make sense
        for beginPOrb in beginPOrbs:
            for endPOrb in endPOrbs:
                if not emptyEndOrbs or endPOrb.isEmpty():
                    # If sink atom has any empty orbitals, the only combinations that make sense
                    #   is if the empty orbital itself is the target

                    # Check if the source and sink orbital are from the same bond,
                    #   to reflect a bond dissociation
                    sameBond = (
                        beginPOrb.isBondOrbital()
                        and endPOrb.isBondOrbital()
                        and (
                            beginPOrb.coveredOrbitalAtomIdx()
                            == endPOrb.coveredOrbitalAtomIdx()
                        )
                    )

                    if (
                        not nonBondedBeginOrbs
                        or not beginPOrb.isBondOrbital()
                        or sameBond
                    ):
                        # If source has any non-bonded electrons, then only combination
                        #   that make sense is if the non-bonded orbital is itself the source
                        #   or if they are the same bond, representing a bond dissociation

                        if compatibleOrbitalPair(beginPOrb, endPOrb):
                            arrowList = []
                            # print "source: %s" % str(beginPOrb.toLabeledSmiAndInfoStr())
                            # print "sink: %s" % str(endPOrb.toLabeledSmiAndInfoStr())
                            # print "BEFORE: %s " % OECreateIsoSmiString(mol);
                            # Before yielding anything first see if we get an exception:
                            try:
                                moveOrbitalElectrons(beginPOrb, endPOrb, 2, arrowList)
                                undoMoveOrbitalElectrons(beginPOrb, endPOrb)
                                # Revert state
                                moveOrbitalElectrons(beginPOrb, endPOrb, 2, arrowList)
                            except:
                                log.info(
                                    "Getting an exception on moveOrbitalElectrons in resIter.rearrange.  Skipping"
                                )
                                continue
                            # print "AFTER: %s " % OECreateIsoSmiString(mol);

                            OEClearAromaticFlags(mol)
                            # Needed to ensure...
                            oldSmi = OECreateIsoSmiString(mol)

                            # Test to cull double formal charge on any atom.
                            shouldYield = exploreResonanceStructure(
                                mol, beginPOrb, endPOrb, includeMinor
                            )

                            OEClearAromaticFlags(mol)
                            # ...SMILES string comes out consistently
                            newSmi = OECreateIsoSmiString(mol)
                            if newSmi != oldSmi:
                                # raise Exception("Resonance structure SMILES string failed to retain state.  Expected %s, got %s" % (oldSmi, newSmi) );
                                log.info(
                                    "Resonance structure SMILES string failed to retain state. Skipping  Expected %s, got %s"
                                    % (oldSmi, newSmi)
                                )
                                continue

                            if shouldYield:
                                yield (mol, arrowList)
                            else:
                                # log.critical('Successfully culled a res struct with a double charge!')
                                pass

                            undoMoveOrbitalElectrons(beginPOrb, endPOrb)
                            # Revert state


def exploreResonanceStructure(mol, beginPOrb, endPOrb, includeMinor=False):
    """Function to test whether the molecule is sufficiently bad enough to not only NOT return as
    a representative res structure, but also bad enough we should stop exploring anything obtainable from
    here as well.

    Cases to consider:
        - Double formal charge on the orbs involved
        - Absolute Gross formal charges >= 3 total

    """
    MAX_GROSS_CHARGE_MINOR = 5
    MAX_GROSS_CHARGE_MAJOR = 3
    maxGrossCharge = MAX_GROSS_CHARGE_MAJOR
    if includeMinor:
        maxGrossCharge = MAX_GROSS_CHARGE_MINOR

    # Testing double formal charge
    shouldYield = True
    shouldYield = shouldYield and abs(beginPOrb.atom.GetFormalCharge()) < 2
    shouldYield = shouldYield and abs(endPOrb.atom.GetFormalCharge()) < 2
    if beginPOrb.neighbor is not None:
        shouldYield = shouldYield and abs(beginPOrb.neighbor.GetFormalCharge()) < 2
    if endPOrb.neighbor is not None:
        shouldYield = shouldYield and abs(endPOrb.neighbor.GetFormalCharge()) < 2

    # return shouldYield;

    if shouldYield is False:
        return shouldYield

    # Then test grossFormalCharges.
    grossPos, grossNeg = grossFormalCharges(mol)
    if abs(grossPos) + abs(grossNeg) > maxGrossCharge:
        return False

    return True


def screenResonanceStruct(mol, includeMinor=False):
    """The molecule should represent a proposed resonance structure.
    Decide whether it should be accepted as a valid / meaningful one.

    For example:
        - Empty orbitals are generally bad, particularly if a pi system can be traced
            from the empty orbital to an available lone pair.
            May allow carbocations as minor contributors here if the none of the potential
            feeder lone pairs are negatively charged or if the negative charge
            rests on an electronegative atom

        - Separated free radicals would not make sense

    Inefficient implementation for now, have to rescan through all atoms
    for each resonance structure considered.  If analyzing multiple
    resonance structures representing the same molecule, should really only
    need to do the thorough analysis once and just reuse the information
    when considering the alternatives.

    includeMinor:
        Option to indicate how strict the screening should be
    """
    for componentMol in splitCompositeMol(mol):
        # print "screening: ", OECreateIsoSmiString(componentMol)
        for atom in componentMol.GetAtoms():
            if abs(atom.GetFormalCharge()) > 1:
                # Formal charges should not exceed 1 magnitude
                # print "Formal charges exceed 1"
                return False

            orbInfo = orbitalInfo(atom)
            if orbInfo["nEmptyOrbitals"] > 0:
                if atomElectronegativity(atom.GetAtomicNum()) > carbonEN:
                    # Electronegative atom with empty orbital is almost never acceptable
                    # Only exception maybe nitrenes.  Will have to adjust for those cases later
                    # print "EN atom with empty orb"
                    return False

                complementaryAtoms = findComplementaryResonanceAtoms(atom, 2)
                if len(complementaryAtoms) > 0:
                    # Resonance alternatives exist.  Case dependent whether this is acceptable
                    if atomElectronegativity(atom.GetAtomicNum()) > carbonEN:
                        # Electronegative atom with empty orbital is never acceptable,
                        #    given alternatives
                        # print "EN complementary atom with empty orb"
                        return False
                    else:
                        # Electropositive atom (probably a carbocation)
                        for compAtom in complementaryAtoms:
                            if compAtom.GetFormalCharge() < 0:
                                # Complementary atom exists with a negative charge
                                #   (indicates charge separated structure).
                                # This could only be acceptable if the complementary atom
                                #   is more electronegative, and even then, only as a
                                #   minor structure
                                if includeMinor:
                                    if (
                                        not atomElectronegativity(
                                            compAtom.GetAtomicNum()
                                        )
                                        > carbonEN
                                    ):
                                        return False
                                else:
                                    # print "Complementary atom with negative charge"
                                    return False
                            else:
                                # Complementary atom without a negative charge.
                                # Make sure no extended neighbors have negative charge
                                #    (lone pairs) either, which can happen for structures
                                #    like nitro groups
                                for extNeighbor in compAtom.GetAtoms():
                                    if extNeighbor.GetFormalCharge() < 0:
                                        # print "Extended neighbor with negative charge"
                                        return False
                else:
                    # No resonance structure alternatives available, so we must accept this
                    pass

            if atom.GetFormalCharge() > 0:
                # Positively charged atom, see if it is part of a pi bond
                #   that could accept electrons
                for bond in atom.GetBonds():
                    if bond.GetOrder() > 1:
                        # More than a sigma bond, must have pi bond potential
                        neighbor = bond.GetNbr(atom)

                        # Look for complementary atoms that could feed into this one by resonance
                        visitedAtoms = set([atom.GetIdx()])
                        complementaryAtoms = findComplementaryResonanceAtoms(
                            neighbor, 2, visitedAtoms
                        )

                        for compAtom in complementaryAtoms:
                            if compAtom.GetFormalCharge() < 0:
                                # Complementary atom exists with a negative charge
                                #   (indicates charge separated structure).
                                # This could only be acceptable only as a minor structure
                                if not includeMinor:
                                    # print "complementary atom has neg charge"
                                    return False

    # Completed iteration without failure, then accept as a pass
    return True


def findComplementaryResonanceAtoms(atom, electrons, visitedAtomIndexes=None):
    """Recursive search for atoms in the same pi system as the parameter atom which
    include a non-bonded orbital containing the specified number of electrons.

    Just look for immediately neighboring atoms.  If those neighboring atoms
    are pi bonded to other atoms, recursively look beyond them for farther
    reaching resonance structure equivalents.
    """
    if visitedAtomIndexes is None:
        # Breadcrumbs to avoid cyclic traversals
        visitedAtomIndexes = set()

    visitedAtomIndexes.add(atom.GetIdx())

    compAtoms = []
    for neighbor in atom.GetAtoms():
        if neighbor.GetIdx() not in visitedAtomIndexes:
            neighborInfo = orbitalInfo(neighbor)
            if electrons == 2 and neighborInfo["nLonePairs"] > 0:
                compAtoms.append(neighbor)
            elif electrons == 1 and neighborInf["nRadicals"] > 0:
                compAtoms.append(neighbor)
            else:
                # Neighbor is not a valid complementary atom.
                # Check for pi bonds though, for extended systems
                for bond in neighbor.GetBonds():
                    if bond.GetOrder() > 1:
                        # More than a sigma bond, means pi bonding potential
                        extNeighbor = bond.GetNbr(neighbor)
                        if extNeighbor.GetIdx() not in visitedAtomIndexes:
                            visitedAtomIndexes.add(neighbor.GetIdx())
                            extCompAtoms = findComplementaryResonanceAtoms(
                                extNeighbor, electrons, visitedAtomIndexes
                            )
                            compAtoms.extend(extCompAtoms)

    return compAtoms


def atomsAltered(mol1, mol2):
    """Given 2 molecules, that should be based on a common one, determine how many
    atoms in the molecule have been altered by electron rearrangement.  One of the
    molecules should probably have been a copy of the other, because this method
    depends on the mol.GetAtoms iterator returning the atoms in the same order.

    If that's the case, will go through each matching pair of atoms and check whether
    they share the same bond orders and formal charge.  If not, assume the electron
    arrangement has changed.
    """
    total = 0
    atomDict1 = dict()
    atomDict2 = dict()
    undoDict1 = {}
    assignAtomMaps(mol1, undoDict1)
    # log.debug("mol1: %s" % OECreateIsoSmiString(mol1));
    undoAssignAtomMaps(mol1, undoDict1)
    undoDict2 = {}
    assignAtomMaps(mol2, undoDict2)
    # log.debug("mol2: %s" % OECreateIsoSmiString(mol2));
    undoAssignAtomMaps(mol2, undoDict2)
    for atom1, atom2 in zip(mol1.GetAtoms(), mol2.GetAtoms()):
        # Generate dictionaries to describe the electron arrangement around each atom
        #   then check whether they are equal
        atomDict1.clear()
        atomDict1["index"] = atom1.GetIdx()
        atomDict1["charge"] = atom1.GetFormalCharge()
        for neighbor in atom1.GetAtoms():
            bond = atom1.GetBond(neighbor)
            atomDict1[neighbor.GetIdx()] = bond.GetOrder()

        atomDict2.clear()
        atomDict2["index"] = atom2.GetIdx()
        atomDict2["charge"] = atom2.GetFormalCharge()
        for neighbor in atom2.GetAtoms():
            bond = atom2.GetBond(neighbor)
            atomDict2[neighbor.GetIdx()] = bond.GetOrder()

        if atomDict1 != atomDict2:
            total += 1
    return total


def allBoolean(iterable):
    """Returns True if *all* in iterable are True.
    Copy of all function implemented in Python v2.5.  Added here to ensure backwards compatibility.
    """
    for element in iterable:
        if not element:
            return False
    return True


def anyBoolean(iterable):
    """Returns True if *any* in iterable are True.
    Copy of all function implemented in Python v2.5.  Added here to ensure backwards compatibility.
    """
    for element in iterable:
        if element:
            return True
    return False
