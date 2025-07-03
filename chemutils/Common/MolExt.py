"""Extensions for molecule object descriptors.

Major functions include
    enumerateRacemicMixture:
        Enumerating specific stereoisomers given a stereochemistry unspecified
        (or incompletely specified) molecule.
"""

# from sets import Set;
# from math import inf;
# from openeye.oechem import OEGetHybridization, OEHybridization_sp, OEHybridization_sp2, OEHybridization_sp3;
from openeye.oechem import (
    OEGetAtomicSymbol,
    OECreateIsoSmiString,
    OEClearAromaticFlags,
    OEGraphMol,
    OEParseSmiles,
)
from openeye.oechem import OEPerceiveChiral, OEGetAverageWeight
from openeye.oechem import (
    OEAtomStereo_Tetrahedral,
    OEAtomStereo_Left,
    OEAtomStereo_Right,
    OEAtomStereo_Undefined,
)
from openeye.oechem import (
    OEBondStereo_CisTrans,
    OEBondStereo_Cis,
    OEBondStereo_Trans,
    OEBondStereo_Undefined,
)
from openeye.oechem import OEDetermineRingSystems
from openeye.oechem import OEGasteigerPartialCharges
from openeye.oechem import OEAssignAromaticFlags, OEKekulize

from chemutils.Common.Util import (
    molBySmiles,
    createStandardSmiString,
    createAtomMapSmiString,
)
from chemutils.Common.Util import (
    getMaxMapIdx,
    splitCompositeSmilesToList,
    joinSmilesListToCompositeSmiles,
)
from chemutils.Common.Util import createStdAtomMapSmiString

from chemutils.Common.OrbitalModel import orbitalInfo

from chemutils.Common.MolStdValue import METAL_ATOMIC_NUMS
from chemutils.Common.MolStdValue import STD_ATOMIC_RADIUS

KEKULIZATION_FLAG = "isKek"


def labelBondToSmi(mol, bond):
    oldBgnMap = bond.GetBgn().GetMapIdx()
    oldEndMap = bond.GetEnd().GetMapIdx()
    bond.GetBgn().SetMapIdx(10)
    bond.GetEnd().SetMapIdx(20)
    smi = createStdAtomMapSmiString(mol)
    bond.GetBgn().SetMapIdx(oldBgnMap)
    bond.GetEnd().SetMapIdx(oldEndMap)
    return smi


def labelAtomToSmi(mol, atom):
    oldMap = atom.GetMapIdx()
    atom.SetMapIdx(10)
    smi = createStdAtomMapSmiString(mol)
    atom.SetMapIdx(oldMap)
    return smi


def atomLocalKekuleIterator(mol, atom):
    """Given a particular atom (part of mol), iterate over kekule structs resulting from flipping each aro bond to be double.

    If the BoolData of KEKULIZATION_FLAG is set, then should not alter for the first run, and should return the mol to the
    original kekulization when all done.

    We do multiple other things in here to ensure proper kekulizations and that future calls to this iterator with diff atoms
    does not alter the mol for the current state of the iterator.  This is handled with adding/removing flags on the bonds
    affected, and the use of a check that the set of bonds that is perceived as aromatic is invariant.
    """
    aroBondTypeMap = {}
    if not mol.GetBoolData(KEKULIZATION_FLAG):
        # log.info('First time popping into kek iter!')
        mol = kekulizeMol(mol)

    # log.info('Popping into kek iter, smi: %s' % createStdAtomMapSmiString(mol))
    aroBondTypeMap = createAroBondTypeMap(mol)
    mol.SetBoolData(KEKULIZATION_FLAG, True)

    OEAssignAromaticFlags(mol)
    aroBonds = []
    currBond = None
    for bond in atom.GetBonds():
        if bond.IsAromatic():
            if bond.GetOrder() < 2:
                aroBonds.append(bond)
            else:
                bond.SetBoolData(KEKULIZATION_FLAG, True)
                bond.SetIntType(2)
                currBond = bond

    # Yield up the first one:
    OEClearAromaticFlags(mol)
    yield mol
    if currBond is not None:
        currBond.SetBoolData(KEKULIZATION_FLAG, False)

    ## Then to check whether a particular kekulization makes sense.  Let's grab the set of aromatic bond indices.
    ## This should be invariant after a particular kekulization
    OEAssignAromaticFlags(mol)
    aroBondIndSet = set([b.GetIdx() for b in mol.GetBonds() if b.IsAromatic()])
    OEClearAromaticFlags(mol)

    # If we get here, then go through and set inttype on each of these to be two, then kekulize, then yield
    for bond in aroBonds:
        # log.info('About to look at new bond.  smi: %s' % labelBondToSmi(mol, bond))
        currAroBondTypeMap = createAroBondTypeMap(mol)
        OEAssignAromaticFlags(mol)
        for aBond in mol.GetBonds():
            if aBond.IsAromatic() and not aBond.GetBoolData(KEKULIZATION_FLAG):
                aBond.SetIntType(5)

        bond.SetIntType(2)
        OEClearAromaticFlags(mol)
        success = False
        try:
            success = OEKekulize(mol)
        except:
            success = False
        ## If we ever get a situation where this bond assignment does not work, then should simply skip this.
        ## Putatively looks good, but double check that the aromatization is the same
        if success:
            OEAssignAromaticFlags(mol)
            newAroBondIndSet = set(
                [b.GetIdx() for b in mol.GetBonds() if b.IsAromatic()]
            )
            success = success and newAroBondIndSet == aroBondIndSet

        if success:
            OEClearAromaticFlags(mol)
            bond.SetBoolData(KEKULIZATION_FLAG, True)
            yield mol
            bond.SetBoolData(KEKULIZATION_FLAG, False)
        else:
            # log.info('Have to reset kek.  OldSmi: %s' % labelAtomToSmi(mol, atom))
            kekulizeMol(mol, currAroBondTypeMap)
            # log.info('Had to reset kek. Smi %s' % labelAtomToSmi(mol, atom))
            # log.info('Had to reset kek. currAroBondTypeMap %s' % pformat(currAroBondTypeMap))

    mol = kekulizeMol(mol, aroBondTypeMap)
    # log.info('Popping out of kekIter, smiles: %s' % createStdAtomMapSmiString(mol))


def deKekulizeMol(mol):
    """Remove kekulization info about aro bonds"""
    OEAssignAromaticFlags(mol)
    for bond in mol.GetBonds():
        if bond.IsAromatic():
            bond.SetIntType(5)
        elif bond.GetOrder() != 0:
            bond.SetIntType(bond.GetOrder())
        else:
            bond.SetIntType(1)
    return mol


def findBondByAtomMapIdx(mol, mapTuple):
    """Convenience to find bonds given a tuple of map indices.  Assumes only two atoms have these"""
    aList = []
    for a in mol.GetAtoms():
        if a.GetMapIdx() in mapTuple:
            aList.append(a)
    if len(aList) != 2:
        return None
    return aList[0].GetBond(aList[1])


def kekulizeMol(mol, aroBondTypeMap=None):
    """Set up a reasonable kekule form.  If we have an aroBondTypeMap, then set it to this."""
    OEAssignAromaticFlags(mol)

    for bond in mol.GetBonds():
        if aroBondTypeMap is not None:
            bond.SetIntType(aroBondTypeMap[bond.GetIdx()])
        elif bond.IsAromatic():
            bond.SetIntType(5)
        elif bond.GetOrder() != 0:
            bond.SetIntType(bond.GetOrder())
        else:
            bond.SetIntType(1)

    OEClearAromaticFlags(mol)
    OEKekulize(mol)
    return mol


def setNonAroBondIntType(mol):
    """Go through a molecule and setup the bond IntType properties for all non-aromatic bonds"""
    OEAssignAromaticFlags(mol)
    for bond in mol.GetBonds():
        if not bond.IsAromatic():
            if bond.GetOrder() != 0:
                bond.SetIntType(bond.GetOrder())
            else:
                bond.SetIntType(1)
    # OEClearAromaticFlags()


def createAroBondTypeMap(mol):
    """Go through the mol and construct a dictionary of bond.idx -> GetIntType for aroBonds"""
    resDict = {}
    OEAssignAromaticFlags(mol)
    for bond in mol.GetBonds():
        # if bond.IsAromatic():
        # resDict[bond.GetIdx()] = bond.GetOrder()
        resDict[bond.GetIdx()] = bond.GetIntType()
    return resDict


def aromatizeMol(mol):
    """Simply call OEAssignAromaticFlags on mol"""
    OEAssignAromaticFlags(mol)


def grossFormalCharges(mol):
    """Figure out the total gross positive and negative charge on the molecule,
    which may be non-zero even if the net charge on the molecule is zero.

    For example:
    >>> from CHEM.Common.Util import molBySmiles;
    >>> print grossFormalCharges(molBySmiles( "CC(=O)Cl" ));
    (0, 0)
    >>> print grossFormalCharges(molBySmiles( "[OH3+]" ));
    (1, 0)
    >>> print grossFormalCharges(molBySmiles( "CC[N+](=O)[O-]" ));
    (1, -1)
    >>> print grossFormalCharges(molBySmiles( "[N-]=[N+]=[N-]" ));
    (1, -2)
    """
    grossPositive = 0
    grossNegative = 0
    for atom in mol.GetAtoms():
        atomCharge = atom.GetFormalCharge()
        if atomCharge > 0:
            grossPositive += atomCharge
        if atomCharge < 0:
            grossNegative += atomCharge
    return (grossPositive, grossNegative)


def atomIsChiral(atom):
    """Extension of OpenEye OEAtomBase.IsChiral method.
    Manual corrections for
        - Organometallics do not have configurationally stable stereochemistry
        - Carbanions should be readily invertible (similar to nitrogen) so achiral
        - onium ions should be considered achiral as well.

    """
    chiral = atom.IsChiral()
    for neighbor in atom.GetAtoms():
        if neighbor.GetAtomicNum() in METAL_ATOMIC_NUMS:
            chiral = False
            break
            # No point in checking more
    if atom.IsCarbon() and atom.GetFormalCharge() == -1:
        chiral = False
    if atom.IsOxygen() and atom.GetFormalCharge() == +1:
        chiral = False
    return chiral


def bondIsChiral(bond):
    """Extension of OpenEye OEBondBase.IsChiral method.
    Manual corrections for
        - Vinylic carbocations should not have stereochemistry, should be linear, sp hybridized
    """
    chiral = bond.IsChiral()
    chiral = chiral and not (
        bond.GetBgn().GetDegree() < 3 and bond.GetBgn().GetFormalCharge() != 0
    )
    chiral = chiral and not (
        bond.GetEnd().GetDegree() < 3 and bond.GetEnd().GetFormalCharge() != 0
    )
    return chiral


def removeNonsenseStereo(mol):
    """Method to remove stereo information about any non-chiral atoms or bonds"""
    OEPerceiveChiral(mol)
    for atom in mol.GetAtoms():
        if not atomIsChiral(atom):
            clearAtomStereo(atom)
    for bond in mol.GetBonds():
        if not bondIsChiral(bond):
            clearBondStereo(bond)


def removeNonsenseStereoSmi(smi):
    """Convenience to run the removeNonsenseStereo on a smiles string"""
    mol = molBySmiles(smi)
    removeNonsenseStereo(mol)
    return createAtomMapSmiString(mol)


def wrapAsOEUnaryAtomPred(func, doNot=False):
    """Function that makes an arbitrary single arg function into a class that extends OEUnaryAtomPred
    And thus is usable in OESubsetMol."""
    from openeye.oechem import OEUnaryAtomPred

    class klass(OEUnaryAtomPred):
        def __call__(self, atom):
            if doNot:
                return not func(atom)
            return func(atom)

    return klass


def atomIsInPiSystem(atom):
    """Predicate function.  Is the atom part of a pi system?
    Cases:
        has a bond with order > 1
        has a lone pair and
            next to pi bond
            next to atom with empty orb
        has an empty orb and
            next to pi bond
            next to atom with lone pair
    """
    atomInfo = orbitalInfo(atom)
    hasLonePair = atomInfo["nLonePairs"] > 0
    hasEmptyOrb = atomInfo["nEmptyOrbitals"] > 0
    for nbr in atom.GetAtoms():
        bond = atom.GetBond(nbr)
        if bond.GetOrder() > 1:
            # Is a pi system
            return True
        nbrInfo = orbitalInfo(nbr)
        if hasLonePair and nbrInfo["nEmptyOrbitals"] > 0:
            return True
        if hasEmptyOrb and nbrInfo["nLonePairs"] > 0:
            return True

        # Finally look for the adjacent pi bonds
        if hasLonePair or hasEmptyOrb:
            for extBond in nbr.GetBonds():
                if extBond.GetOrder() > 1:
                    return True
    return False


def createStandardMol(mol):
    """Reparse the molecule SMILES string to ensure any unusual properties are cleared out.
    Happens often when doing direct molecule graph manipulations or reaction processing.

    In addition, check molecule bonds and atoms and clear out any bogus stereoassignments.
    """
    stdSmi = createStandardSmiString(mol)

    mol = molBySmiles(stdSmi)
    # New instance / copy

    atomClearCount = 1000
    # Number of atoms with bogus stereoassignment cleared.  Start with sentinel value
    while (
        atomClearCount > 0
    ):  # Loop structure because clearing one stereocenter can affect the perception of another stereocenter
        atomClearCount = 0
        stdSmi = createStandardSmiString(mol)
        # Have to reparse standard SMILES again to ensure changes picked up
        mol.Clear()
        OEParseSmiles(mol, stdSmi)
        OEPerceiveChiral(mol)
        for atom in mol.GetAtoms():
            if not atomIsChiral(atom) and atom.HasStereoSpecified():
                # False chiral center, clear the stereo flag
                clearAtomStereo(atom)
                atomClearCount += 1

    bondClearCount = 1000
    # Number of bonds with bogus stereoassignment cleared.  Start with sentinel value
    while (
        bondClearCount > 0
    ):  # Loop structure because clearing one stereocenter can affect the perception of another
        bondClearCount = 0
        stdSmi = createStandardSmiString(mol)
        mol.Clear()
        OEParseSmiles(mol, stdSmi)
        OEPerceiveChiral(mol)
        for bond in mol.GetBonds():
            if not bondIsChiral(bond) and bond.HasStereoSpecified():
                # False chiral center, clear the stereo flag
                clearBondStereo(bond)
                bondClearCount += 1

    return mol


def enumerateRacemicMixture(
    mol, unspecifiedAtoms=None, unspecifiedBonds=None, allCisTrans=False
):
    """Assume the molecule object represents a racemic mixture of otherwise identical isomers.
    Go through every unassigned stereocenter and generate all possible specific isomers
    by assiging all possible stereo configurations.

    allCisTrans: If true, will enumerate all cis-trans isomers for unspecified chiral bonds.
        Otherwise, will only generate one with the bulkiest (largest by molecular weight)
        groups trans to each other
    """
    if unspecifiedAtoms is None or unspecifiedBonds is None:
        # Initial call.  Standardize the molecule first to ensure nothing funny happens
        mol = createStandardMol(mol)

    if unspecifiedAtoms is None:
        unspecifiedAtoms = []
        for atom in mol.GetAtoms():
            if atomIsChiral(atom) and not atom.HasStereoSpecified():
                unspecifiedAtoms.append(atom)

    if unspecifiedBonds is None:
        OEPerceiveChiral(mol)
        unspecifiedBonds = []
        for bond in mol.GetBonds():
            if bondIsChiral(bond) and not bond.HasStereoSpecified():
                unspecifiedBonds.append(bond)

    # Pop one feature at a time and enumerate both possibilities, then recurse
    if len(unspecifiedAtoms) > 0:
        chiralAtom = unspecifiedAtoms.pop()
        neighbors = []
        for atom in chiralAtom.GetAtoms():
            neighbors.append(atom)

        chiralAtom.SetStereo(neighbors, OEAtomStereo_Tetrahedral, OEAtomStereo_Left)
        for isomer in enumerateRacemicMixture(mol, unspecifiedAtoms, unspecifiedBonds):
            yield isomer

        chiralAtom.SetStereo(neighbors, OEAtomStereo_Tetrahedral, OEAtomStereo_Right)
        for isomer in enumerateRacemicMixture(mol, unspecifiedAtoms, unspecifiedBonds):
            yield isomer

        unspecifiedAtoms.append(chiralAtom)
        # Restore state

    elif len(unspecifiedBonds) > 0:
        chiralBond = unspecifiedBonds.pop()
        neighbors = []
        # Find bulkiest adjacent neighbor atom, other than the same bonded atom
        bulkiestNeighbor = None
        bulkiestNeighborWeight = None
        visitedAtomIndexes = set()

        for neighbor in chiralBond.GetBgn().GetAtoms():
            if neighbor.GetIdx() != chiralBond.GetEnd().GetIdx():
                visitedAtomIndexes.clear()
                visitedAtomIndexes.add(chiralBond.GetBgn().GetIdx())
                neighborWeight = neighborTreeWeight(neighbor, visitedAtomIndexes)
                if bulkiestNeighbor is None or neighborWeight > bulkiestNeighborWeight:
                    bulkiestNeighbor = neighbor
                    bulkiestNeighborWeight = neighborWeight
        neighbors.append(bulkiestNeighbor)

        bulkiestNeighbor = None
        for neighbor in chiralBond.GetEnd().GetAtoms():
            if neighbor.GetIdx() != chiralBond.GetBgn().GetIdx():
                visitedAtomIndexes.clear()
                visitedAtomIndexes.add(chiralBond.GetEnd().GetIdx())
                neighborWeight = neighborTreeWeight(neighbor, visitedAtomIndexes)
                if bulkiestNeighbor is None or neighborWeight > bulkiestNeighborWeight:
                    bulkiestNeighbor = neighbor
                    bulkiestNeighborWeight = neighborWeight
        neighbors.append(bulkiestNeighbor)

        chiralBond.SetStereo(neighbors, OEBondStereo_CisTrans, OEBondStereo_Trans)
        for isomer in enumerateRacemicMixture(mol, unspecifiedAtoms, unspecifiedBonds):
            yield isomer

        if allCisTrans:
            chiralBond.SetStereo(neighbors, OEBondStereo_CisTrans, OEBondStereo_Cis)
            for isomer in enumerateRacemicMixture(
                mol, unspecifiedAtoms, unspecifiedBonds
            ):
                yield isomer

        unspecifiedBonds.append(chiralBond)
        # Restore state

    else:
        # Base case, no more unspecified atoms.  Just yield a copy of the current molecule
        yield OEGraphMol(mol)


def clearAtomStereo(atom):
    if atom.HasStereoSpecified():
        neighbors = []
        for nbr in atom.GetAtoms():
            neighbors.append(nbr)
        atom.SetStereo(neighbors, OEAtomStereo_Tetrahedral, OEAtomStereo_Undefined)


def clearBondStereo(bond):
    if bond.HasStereoSpecified():
        neighbors = []
        # Find any adjacent neighbor atom, other than the same bonded atom
        for atom in bond.GetBgn().GetAtoms():
            if atom.GetIdx() != bond.GetEnd().GetIdx():
                neighbors.append(atom)
                break
        for atom in bond.GetEnd().GetAtoms():
            if atom.GetIdx() != bond.GetBgn().GetIdx():
                neighbors.append(atom)
                break

        bond.SetStereo(neighbors, OEBondStereo_CisTrans, OEBondStereo_Undefined)


def countNumAtmMappedSingleSmi(compositeSmi):
    """Method to count the number of individual single smi in a compositeSmi where there are atom maps.
    Example:

    >>> smi = 'C[CH+:1]C.[Br-:2].[Br-]'
    >>> countNumAtmMappedSingleSmi(smi);
    2
    >>> smi = 'C[CH+:1]C.[Br-].[Br-]'
    >>> countNumAtmMappedSingleSmi(smi);
    1
    >>> smi = 'C[CH+]C.[Br-].[Br-]'
    >>> countNumAtmMappedSingleSmi(smi);
    0
    """
    numAtmMapped = 0
    mol = OEGraphMol()
    for smi in splitCompositeSmilesToList(compositeSmi):
        mol.Clear()
        OEParseSmiles(mol, smi)
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() > 0:
                numAtmMapped += 1
                break

    return numAtmMapped


def hasAtomMaps(mol):
    """Predicate to determine if there are atom maps in this mol"""
    for atom in mol.GetAtoms():
        if atom.GetMapIdx() != 0:
            return True
    return False


def getAtmMappedSmilesFromCompositeSmiles(compositeSmi, clearMaps=True):
    """Method to return a smile string that has all of the mols in the compositeSmi that are atom mapped

    >>> smi = 'C[CH+:1]C.[Br-:2].[Br-]'
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, True);
    'C[CH+]C.[Br-]'
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, False);
    'C[CH+:1]C.[Br-:2]'
    >>> smi = 'C[CH+]C.[Br-].[Br-]';
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, False);
    ''
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, True);
    ''
    >>> smi = 'C[CH+]C.[Br-:2].[Br-]'
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, True);
    '[Br-]'
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, False);
    '[Br-:2]'
    >>> smi = 'C[CH+:20]C1=CC=CC=C1.C[CH+]C1=CC=CC=C1.C=C[C:13]1=[CH:12][CH:11]=[CH:10]C=C1.C=CC1=CC=CC=C1'
    >>> getAtmMappedSmilesFromCompositeSmiles(smi, False);
    'C[CH+:20]C1=CC=CC=C1.C=C[C:13]1=[CH:12][CH:11]=[CH:10]C=C1'
    """
    atmMappedSmiList = []
    mol = OEGraphMol()
    for smi in splitCompositeSmilesToList(compositeSmi):
        mol.Clear()
        OEParseSmiles(mol, smi)
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() > 0:
                atmMappedSmiList.append(smi)
                break

    totalSmi = joinSmilesListToCompositeSmiles(atmMappedSmiList)
    mol.Clear()
    OEParseSmiles(mol, totalSmi)

    if clearMaps:
        clearAtomMaps(mol)

    return createStdAtomMapSmiString(mol)


#
# def canonicalKekule(smi):
#    """Function to make sure that the kekule structure is canonical
#
#    >>> smi = 'c1ccc(cc1)[N+](=O)[O-]'
#    >>> canonicalKekule(smi)
#    'C1=CC=C(C=C1)[N+](=O)[O-]'
#    >>> smi = '
#    """
#    mol = value_to_oemol(smi)
#    removeNonsenseStereo(mol)
#    OEAssignAromaticFlags(mol)
#    smi = createStdAtomMapSmiString(mol)
#    mol = value_to_oemol(smi)
#    OEClearAromaticFlags(mol)
#    return createStdAtomMapSmiString(mol)


def setSingleExplicitHydrogens(mol, setAtomMaps=False):
    """Method to make sure a single implicit hydrogen on each atom is made explicit.

    if setAtomMaps - then the new hydrogens will be given mapIdx's based on the maxMapIdx of
    the mol.  Example:

    >>> from CHEM.Common.Util import molBySmiles, createAtomMapSmiString;
    >>> mol = molBySmiles('[BH3:1]')
    >>> setSingleExplicitHydrogens(mol, True)
    >>> createAtomMapSmiString(mol)
    '[H:2][BH2:1]'
    >>> setSingleExplicitHydrogens(mol, True)
    >>> createAtomMapSmiString(mol)
    '[H:2][BH2:1]'
    >>> mol = molBySmiles('O')
    >>> len([atm for atm in mol.GetAtoms()])
    1
    >>> setSingleExplicitHydrogens(mol)
    >>> len([atm for atm in mol.GetAtoms()])
    2
    >>> setSingleExplicitHydrogens(mol)
    >>> len([atm for atm in mol.GetAtoms()])
    2
    """
    maxMapIdx = getMaxMapIdx(mol)

    for atom in mol.GetAtoms():
        implicitHCount = atom.GetImplicitHCount()
        if implicitHCount > 0:
            seenHydrogen = False
            for nbr in atom.GetAtoms():
                if nbr.GetAtomicNum() == 1:
                    seenHydrogen = True
                    break

            if not seenHydrogen:
                newAtom = mol.NewAtom(1)
                newBond = mol.NewBond(atom, newAtom, 1)
                atom.SetImplicitHCount(implicitHCount - 1)

                if setAtomMaps:
                    newAtom.SetMapIdx(maxMapIdx + 1)
                    maxMapIdx += 1


def clearMolStereo(mol):
    """Convenience to clear stereo on all atoms and bonds in mol"""
    [clearBondStereo(bond) for bond in mol.GetBonds()]
    [clearAtomStereo(atom) for atom in mol.GetAtoms()]


def clearAtomMaps(mol):
    """Convenience to clear the atom mapping on a whole mol
    Removes isotope settings on the atom mapped hydrogens as well."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0:
            # Hydrogen atom with atom mapping.  Artificially set isotope, to zero
            atom.SetIsotope(0)
    [atm.SetMapIdx(0) for atm in mol.GetAtoms()]


def createUnAtomMapIsoSmiString(mol):
    """Convenience to take a mol and return an unatommapped iso smi string in non-destructive manner"""
    copyMol = OEGraphMol(mol)
    clearAtomMaps(copyMol)
    return OECreateIsoSmiString(copyMol)


def isSmiStereoSpecified(smi):
    """Wrapper to call isMolStereoSpecified with a smile string"""
    mol = molBySmiles(smi)
    return isMolStereoSpecified(mol)


def isMolStereoSpecified(mol):
    """Boolean function returns True if any atom or bond in a mol has stereo specified."""
    stereo = reduce(
        lambda x, y: x or y, [atm.HasStereoSpecified() for atm in mol.GetAtoms()], False
    )
    return stereo or reduce(
        lambda x, y: x or y, [bnd.HasStereoSpecified() for bnd in mol.GetBonds()], False
    )


def atomMapAll(mol):
    """Convenience to assign mapIDx to all atoms that dont already have one"""
    currIdx = [atm.GetMapIdx() for atm in mol.GetAtoms() if atm.GetMapIdx() != 0]
    newIdx = [i + 1 for i in range(mol.NumAtoms()) if i + 1 not in currIdx]
    i = 0
    for atm in mol.GetAtoms():
        if atm.GetMapIdx() == 0:
            atm.SetMapIdx(newIdx[i])
            i += 1


def enumerateRacemicMixtureList(molList, currStereoList=None):
    """Given a list of molecule objects, call enumerateRacemicMixture on each to
    yield a list of all possible isomer combinations of the molecules.

    currStereoList is the current list of fully-stereospecified molecules.
    One it reaches the same length of the original molList, one complete set
    of stereoisomers has been produced and can be returned.  Use to track depth in recursion
    """
    if currStereoList is None:
        currStereoList = []

    if len(currStereoList) >= len(molList):
        # Base case, no further to recurse
        yield currStereoList
    else:
        # Need to add more
        nextPosition = len(currStereoList)
        for stereoisomer in enumerateRacemicMixture(molList[nextPosition]):
            currStereoList.append(stereoisomer)
            for resultList in enumerateRacemicMixtureList(molList, currStereoList):
                yield resultList
                # Recursive generation for subsequent positions
            currStereoList.pop()
            # Restore state


def neighborTreeWeight(atom, visitedAtomIndexes):
    """Determine the molecular weight of all atoms connected up to the atom,
    that are not counted in the visitedAtomIndexes.  Ignore implicit hydrogens,
    which should be fine for relative comparisons as long as they are
    consistently ignored throughout the molecule.
    """
    weight = OEGetAverageWeight(atom.GetAtomicNum())
    visitedAtomIndexes.add(atom.GetIdx())
    for childAtom in atom.GetAtoms():
        if childAtom.GetIdx() not in visitedAtomIndexes:
            weight += neighborTreeWeight(childAtom, visitedAtomIndexes)
    return weight


def assignPartialCharges(mol):
    """If non-organic atoms are in the mol, OEGasteigerPartialCharges will set partial charges to be 'nan'
    This happens with cations or Grignard reagents.  A fix is if this happens, simply assign all
    the partial charges on the atoms to be the formal charges.
    """
    doesPass = OEGasteigerPartialCharges(mol)
    # atom = mol.GetAtoms().next();

    # This is a way to test for nan's (they will not be equivalent to self.)
    # This is fixed in python 2.6 (with math.isnan)
    if not doesPass:
        # atom.GetPartialCharge() != atom.GetPartialCharge() or atom.GetPartialCharge() == float('inf') or atom.GetPartialCharge() == -float('inf') or atom.GetPartialCharge() == float('nan'):
        for atom in mol.GetAtoms():
            atom.SetPartialCharge(atom.GetFormalCharge())


def atomSize(atomicNum):
    """Abstracted measure of atom size for assessing steric bulk"""
    return STD_ATOMIC_RADIUS[atomicNum]


def substituentStericBulk(atom, visitedAtomIndexes=None, depth=0, baseAtom=None):
    """Similar to neighborTreeWeight.
    Use to estimate the steric bulk of substituent groups.
    Ideally, would quantitatively predict A-values,
    but at least get a coarse approximation for now.

    TO DO:
    #** Need something to approximate A values:
    #       MolBranchAValue: (adapt from MolExt.neighborTreeWeight function)
    #       Input orbital (no, use general idea below),
    #           accumulate atom sizes with decreasing effect over distance?
    #           base on atomic radii?
    #       Input root atom (and excluded atoms, probably the opposite end of bond
    #           of interest), and accumulate from there.
    #       Try using A values to calibrate scale


    Basic strategy, recursively accumulate sizes of substituent atoms
    (measured by atomic radii), but contribution becomes progressively less by distance
    (recursive depth) to reflect reduced proximity to the relevant space.

    atom:
        Root atom of the substituent branch to be considered
    visitedAtomIndexes:
        Set of atom indexes that have already been traversed in recursion to prevent
        back-recursion.  The base atom of the structure the substituent root atom
        is attached to should start in this Set to ensure only the substituent atoms
        are counted.
    depth:
        Recursion depth.  Limit recursions to fixed number, as effects
        should quickly diminish with distance / recursion.
    baseAtom:
        Can provide this directly instead of through the visitedAtomIndexes.
        Just a convenience to have the function take care of adding it to the set.
    """
    MAX_RECURSIONS = 2
    SCALE_PER_RECURSION = 2

    if visitedAtomIndexes is None:
        visitedAtomIndexes = set()
    if baseAtom is not None:
        visitedAtomIndexes.add(baseAtom.GetIdx())

    bulk = atomSize(atom.GetAtomicNum())
    visitedAtomIndexes.add(atom.GetIdx())
    if depth < MAX_RECURSIONS:
        for childAtom in atom.GetAtoms():
            if childAtom.GetIdx() not in visitedAtomIndexes:
                bulk += (
                    substituentStericBulk(childAtom, visitedAtomIndexes, depth + 1)
                    / SCALE_PER_RECURSION
                )

        # Account for implicit hydrogens
        bulk += atomSize(1) * atom.GetImplicitHCount() / SCALE_PER_RECURSION

    return bulk


def molSkeleton(
    mol,
    clearChiralFlags=True,
    copyMol=True,
    carbonSkeleton=False,
    retainSingleHs=False,
    retainHs=False,
):
    """Reduces the molecule down to just its (carbon) skeleton.
    - All implicit and explicit hydrogen atoms are removed (unless they have incorrectly have >1 bond)
    - All bond orders are set to 1
    - All formal charges are neutralized
    - All isotopes are reset
    - All atom mapping indexes are reset

    If carbonSkeleton is True
        - All non-carbon atoms that are not within a ring or chain are removed
    If retainSingleHs is True, then only delete these if they are connected to some besides hydrogen.
    """

    if copyMol:
        mol = OEGraphMol(mol)

    # First wipe out all of the unnecessary hydrogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and not atom.GetDegree() > 1:
            if retainSingleHs and (
                atom.GetDegree() == 0
                or all([nbr.GetAtomicNum() == 1 for nbr in atom.GetAtoms()])
            ):
                continue
            else:
                if atom.GetExplicitDegree() == 1 and retainHs:
                    nbr = [a for a in atom.GetAtoms()][0]
                    nbr.SetImplicitHCount(nbr.GetImplicitHCount() + 1)
                mol.DeleteAtom(atom)

    # OEAssignImplicitHydrogens(mol)

    # Separate pass for heavy atoms that won't get confused by hydrogens
    for atom in mol.GetAtoms():
        if (
            carbonSkeleton
            and atom.GetAtomicNum() != 6
            and not atom.IsInRing()
            and not (atom.GetDegree() - atom.GetImplicitHCount()) > 1
        ):
            mol.DeleteAtom(atom)
        else:
            if not retainHs:
                atom.SetImplicitHCount(0)
            atom.SetFormalCharge(0)
            atom.SetIsotope(0)
            atom.SetMapIdx(0)
            atom.SetAromatic(False)
            if clearChiralFlags:
                clearAtomStereo(atom)

    for bond in mol.GetBonds():
        bond.SetOrder(1)
        if clearChiralFlags:
            clearBondStereo(bond)

    return mol


def clearMol(mol):
    """Delete all of the atoms in the molecule.  Don't just use mol.Clear(),
    if want to retain other information, such as SD tagged or other name-value pair data
    """
    for atom in mol.GetAtoms():
        mol.DeleteAtom(atom)


def isCarbonRadical(atom):
    """Determine if the atom represents a carbon free radical"""
    return atom.IsCarbon() and atom.GetFormalCharge() == 0 and atom.GetValence() == 3


def determineRingSystemInfo(mol):
    """Extension of OEDetermineRingSystems to not only report
    the number of ring systems and map which atoms belong to which ring systems,
    but provide additional information on ring branching points which can be
    used to infer spirocycle centers and bridgehead carbons.

    Specifically returns an information dictionary with the following keyed items
        ringSystemCount:
            Number of ring systems
        ringSystemAtomMap:
            Actually a list, but works like a map from atom index numbers to
            the index of the ring system the atom belongs to
        ringBranchingAtomMap:
            Another atom property list / map.  Values equal the ring branching
            number for each atom.  That is, for an atom in a ring system,
            return the number of neighbor atoms that belong to the same ring system.
            Values >= 3 indicate branch points, values == 4 indicate spirocenters
        ringNeighborBranchesAtomMap:
            Another atom property list / map.  Values equal the number of
            neighboring atoms that are branch points in the same ring system.
            Values >= 1 indicate a ring fusion point.
            Atoms where (neighborBranches - ringBranching) >= 3 indicate a bridgehead.

    Beware: Above defintion of "bridgehead" is probably too narrow.  It should work for
    the vast majority of cases, but consider a molecule like this SMILES:

        C1CC2CC[CH:1]3CCC(C1)[CH:2]2[CH2:3]3

    Carbon 2 is almost certain a bridgehead (bridges to branch point 1 via bridge carbon 2),
    but it also happens to be a part of fused ring system components (branched neighbors)
    which make it less obvious.  Stricter criteria would look for all pairs of branch points
    and trace every path between each branch point to find non-fused paths, but this would
    be a fairly expensive operation, so suffice with the current approximation.
    """
    (ringSystemCount, ringSystemAtomMap) = OEDetermineRingSystems(mol)

    infoDict = dict()
    infoDict["ringSystemCount"] = ringSystemCount
    infoDict["ringSystemAtomMap"] = ringSystemAtomMap

    infoDict["ringBranchingAtomMap"] = [0] * mol.NumAtoms()
    infoDict["ringNeighborBranchesAtomMap"] = [0] * mol.NumAtoms()

    # First pass to figure out branching numbers
    for atom in mol.GetAtoms():
        if atom.IsInRing():
            atomRingIndex = ringSystemAtomMap[atom.GetIdx()]
            branchingNumber = 0
            for neighbor in atom.GetAtoms():
                neighborRingIndex = ringSystemAtomMap[neighbor.GetIdx()]
                if neighborRingIndex == atomRingIndex:
                    branchingNumber += 1
            infoDict["ringBranchingAtomMap"][atom.GetIdx()] = branchingNumber

    # Second pass to figure out which atoms have neighboring branch points (branching number >= 3)
    for atom in mol.GetAtoms():
        neighborBranches = 0
        for neighbor in atom.GetAtoms():
            neighborBranchingNumber = infoDict["ringBranchingAtomMap"][
                neighbor.GetIdx()
            ]
            isNeighborBranch = neighborBranchingNumber >= 3
            if isNeighborBranch:
                neighborBranches += 1
        infoDict["ringNeighborBranchesAtomMap"][atom.GetIdx()] = neighborBranches

    return infoDict


def isBridgehead(atom, ringSystemInfo):
    """Determine if the given atom is a "bridgehead" atom
    of a ring system, using the ringSystemInfo
    to perceive which atoms and neighbors are part of the same ring system
    and possible branch points.

    Bridgeheads exist when there are all (3) multiple paths from one branch point
    to another branch point with intervening "bridge" atoms (non-branch points)
    """
    ringBranches = ringSystemInfo["ringBranchingAtomMap"][atom.GetIdx()]
    branchNeighbors = ringSystemInfo["ringNeighborBranchesAtomMap"][atom.GetIdx()]

    isBranchPoint = ringBranches >= 3
    nonBranchNeighbors = ringBranches - branchNeighbors
    # Should equal number of "bridge" neighbors

    return isBranchPoint and (nonBranchNeighbors >= 3)


def molContentDiff(mol1, mol2):
    """Return a dictionary to catalog the difference in number of
    specific atoms and electrons that exist in the two input molecules.
    Fair comparisons between molecular energy states should
    satisfy having zero "diff" in this measure between the (composite) molecules.
    """
    # First look at the atoms.
    diffDict = {}
    diffDict["e-"] = 0
    for atom in mol2.GetAtoms():
        atomSymb = OEGetAtomicSymbol(atom.GetAtomicNum())
        if atomSymb not in diffDict:
            diffDict[atomSymb] = 0
        diffDict[atomSymb] += 1

        # Check for implicit hydrogens
        atomSymb = OEGetAtomicSymbol(1)
        if atomSymb not in diffDict:
            diffDict[atomSymb] = 0
        diffDict[atomSymb] += atom.GetImplicitHCount()

        # Look at electrons
        diffDict["e-"] += orbitalInfo(atom)["nUnbondedElectrons"]
        diffDict["e-"] += atom.GetImplicitHCount() * 2
        # Assume a single bond to each implict hydrogen

    for atom in mol1.GetAtoms():
        atomSymb = OEGetAtomicSymbol(atom.GetAtomicNum())
        if atomSymb not in diffDict:
            diffDict[atomSymb] = 0
        diffDict[atomSymb] -= 1

        # Check for implicit hydrogens
        atomSymb = OEGetAtomicSymbol(1)
        if atomSymb not in diffDict:
            diffDict[atomSymb] = 0
        diffDict[atomSymb] -= atom.GetImplicitHCount()

        # Look at electrons
        diffDict["e-"] -= orbitalInfo(atom)["nUnbondedElectrons"]
        diffDict["e-"] -= atom.GetImplicitHCount() * 2
        # Assume a single bond to each implict hydrogen

    # Traverse the bonds
    numE = 0
    for bond in mol2.GetBonds():
        numE += bond.GetOrder() * 2
    for bond in mol1.GetBonds():
        numE -= bond.GetOrder() * 2
    diffDict["e-"] += numE

    keyList = diffDict.keys()
    for key in keyList:
        if diffDict[key] == 0:
            del diffDict[key]

    return diffDict
