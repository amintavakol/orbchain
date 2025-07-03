from openeye.oechem import OEDetermineComponents

from chemutils.Common.OrbitalModel import Orbital
from chemutils.Common.MolExt import atomLocalKekuleIterator, kekulizeMol

from chemutils.CombiCDB.proposers.Util import compatibleOrbitalPair, labelOrbPair


class AtomBasedReactionProposer(object):
    """Propose orbital interactions based on restricted atoms.

    Assumes we are given two lists of SimpleAtomObjects, one for potential source, one for sinks
    """

    def __init__(
        self,
        srcAtomList,
        sinkAtomList,
        orbitalPairFilters=None,
        checkResonance=True,
        checkIntra=True,
        checkRadical=False,
    ):
        """Constructor"""
        self.srcAtomList = srcAtomList
        self.sinkAtomList = sinkAtomList
        self.orbitalPairFilters = orbitalPairFilters
        self.mol = None
        self.seenOrbPairSet = set([])
        self.componentParts = None
        self.checkRadical = checkRadical
        self.checkIntra = checkIntra
        self.interOnly = False

        if len(self.srcAtomList) > 0:
            self.mol = list(self.srcAtomList)[0].mol
            self.numParts, self.partsMap = OEDetermineComponents(self.mol)
            self.mol = kekulizeMol(self.mol)

            # big errors in original code here. used to be called "intraOnly" which meant the opposite.
            # set interOnly to false, meaning we want intra considered no matter how many reactants.
            # self.interOnly = self.numParts > 1 and self.checkIntra
            self.interOnly = False

        # self.checkResonance = checkResonance
        # self.resonanceStructUtil = None
        # if self.checkResonance and self.mol is not None:
        #    self.resonanceStructUtil = ResonanceStructUtil(self.mol, checkRadical=checkRadical)

        self.srcOrbFactory = None
        self.sinkOrbFactory = None
        self.setOrbitalFactories()

        # self.inter = inter

    def setOrbitalFactories(self):
        """Abstract - Must implement"""
        raise NotImplementedError()

    def passOrbitalPairFilters(self, src, sink):
        strRep = labelOrbPair(src, sink)
        if strRep in self.seenOrbPairSet:
            return False
        self.seenOrbPairSet.add(strRep)

        # Check if any of the orbital pair filters do NOT pass
        for filterFunc in self.orbitalPairFilters:
            if not filterFunc(src, sink):
                return False

        # if self.checkResonance:
        #    if not self.resonanceStructUtil.checkOrbPair(src, sink):
        #        return False

        return True

    def proposeOrbitalPairs(self):
        """Main work of proposing.  Loop over the factories and make inter/intra orb pairs.

        NOTE: To get full coverage of the kekule structures when we have potential intra-molecular
        reactions within the same aromatic system, we need to NOT pre-compute the sinkOrbList.
        """

        ## Pre-make the sinkOrbs
        # sinkOrbList = [orb for orb in self.sinkOrbFactory]

        for srcOrb, srcAtom in self.srcOrbFactory:
            for sinkOrb, sinkAtom in self.sinkOrbFactory:  # sinkOrbList:
                # First do the inter check.  May want to skip if not from different components
                if self.interOnly and (
                    self.partsMap[srcOrb.atom.GetIdx()]
                    == self.partsMap[sinkOrb.atom.GetIdx()]
                ):
                    continue

                # First intra:
                if compatibleOrbitalPair(srcOrb, sinkOrb):
                    (
                        (clonedSrcOrb, clonedSinkOrb),
                        copyMol,
                    ) = Orbital.cloneOrbitalsWithCommonMol([srcOrb, sinkOrb])
                    yield (clonedSrcOrb, clonedSinkOrb, srcAtom, sinkAtom)

                    # NOTE: this appears to control whether duplicate copies of species are used for src/sink
                    # Then inter:
                    ## Only do this if from same mol
            # else:
            # if self.partsMap[srcOrb.atom.GetIdx()] == self.partsMap[sinkOrb.atom.GetIdx()]:
            # ( (compositeSrcOrb, compositeSinkOrb), compositeMol ) =  \
            #                Orbital.compositeOrbitalsFromDistinctMols( [srcOrb, sinkOrb] );
            # yield (compositeSrcOrb, compositeSinkOrb, srcAtom, sinkAtom)

    def __iter__(self):
        """Iterator to yield the orbitals.  Propose and check if we pass orbital filters"""
        for srcOrb, sinkOrb, srcAtom, sinkAtom in self.proposeOrbitalPairs():
            if self.passOrbitalPairFilters(srcOrb, sinkOrb):
                yield (srcOrb, sinkOrb, srcAtom, sinkAtom)


class AtomBasedOrbitalFactory(object):
    """Abstract Virtual class defining the expected interface.

    Handles how to run the atomIterator.
    User must implement iter
    """

    def __init__(self, atomObjList, allSymmetry=True, allKekule=True):
        """Basic Constructor"""
        self.atomObjList = atomObjList
        self.allSymmetry = allSymmetry
        self.allKekule = allKekule

    def atomIterator(self):
        """Convenience to handle how we iterate over the atoms in the"""
        for atomObj in self.atomObjList:
            atomList = [atomObj.atom]
            if self.allSymmetry:
                atomList = atomObj.symmetricAtoms

            for oeAtom in atomList:
                # If we need to iterate over the kekule structs, do so:
                mol = atomObj.mol
                # log.info('Starting with mol with smiles: %s' % (createCanonicalAtomMapSmiString(mol)))
                if self.allKekule:
                    for iMol, nMol in enumerate(atomLocalKekuleIterator(mol, oeAtom)):
                        oeAtom.SetMapIdx(0)
                        yield oeAtom, atomObj
                else:
                    yield oeAtom, atomObj

    def __iter__(self):
        """Needs to be implemented.

        Expected to return the orbital and the current atomObj
        """
        raise NotImplementedError()
