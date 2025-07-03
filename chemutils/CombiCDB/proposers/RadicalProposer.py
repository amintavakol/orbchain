from openeye.oechem import OEBondIsInRingSize

from chemutils.Common.MolStdValue import PERIODIC_ROW
from chemutils.CombiCDB.proposers.AtomBasedReactionProposer import (
    AtomBasedReactionProposer,
    AtomBasedOrbitalFactory,
)
from chemutils.Common.OrbitalModel import radicalOrbitalIter, orbitalInfo
from chemutils.score.orbital.ReactivityScore import ReactiveOrbitalFactory
from chemutils.CombiCDB.proposers import filters

defaultFilters = [
    filters.NonReasonableOrbitalsFilter(),
    filters.OverlapAtomsFilter(),
    filters.ExtremeChargeFilter(),
    filters.ExtremeBondFilter(),
    filters.PericyclicFilter(),
]


class RadicalProposer(AtomBasedReactionProposer):
    """Exhaustive radical reaction proposer (given list of possible atoms)

    Take into account a flag for photoexcitation
    Unlike the Polar Proposer, a simple set of rules are used for proposals:

    Making new radical species
    - Sigma bond homolysis will only occur for the following:
      - O-O, C-Heavy, X-X, C-N=N bonds, Strained Other (3-rings)
      - Can occur with photoexcitation or heat, so allow always
    - Pi bond homolysis will only occur
      - with photoexcitation
    - Single electron transfers
      - Are handled as single arrows (not multiple as normally done.)
        - Need to work out the details in the moveOrbitalElectrons methods
      - Li, Na (anytime) (could include SmI2 in here, but don't)
      - From a lone pair (only with photoexcitation)
      - Other types not handled:
        - Pb(OAc)4, Mn(OAc)3, CAN, DDQ, chloranil for oxidations
    - Cycloaromatizations - not handled here

    Reactions with existing radicals:
    - addition to Pi bond
    - Fragmentation -
    - atom abstraction - Abstraction of C DOES NOT occur
    - Radical-radical combination
    - Disproportionation
    - 1,2 radical shifts do not occur, but can be done for unsaturated species via
      addition-Fragmentation through a cyclopropane intermediate

    """

    def __init__(
        self,
        srcAtomList,
        sinkAtomList,
        orbitalPairFilters=defaultFilters,
        photo=True,
        checkIntra=True,
    ):
        """Basic constructor - call super and set the photo attribute"""
        self.photo = photo
        super(RadicalProposer, self).__init__(
            srcAtomList,
            sinkAtomList,
            orbitalPairFilters,
            checkIntra=checkIntra,
            checkRadical=True,
        )

    def setOrbitalFactories(self):
        """Specify the orbital factories to use"""
        self.srcOrbFactory = RadicalOrbitalFactory(
            self.srcAtomList, occupied=True, allSymmetry=True, photo=self.photo
        )
        self.sinkOrbFactory = RadicalOrbitalFactory(
            self.sinkAtomList, occupied=False, allSymmetry=False, photo=self.photo
        )


class RadicalOrbitalFactory(AtomBasedOrbitalFactory):
    """Class to handle the basic generation of radical based orbitals.

    Super class handles kekulization and symmetry in terms of the atomiteration through
    its implementation of the atomIterator method.
    """

    def __init__(
        self,
        atomObjList,
        occupied=False,
        allSymmetry=False,
        allKekule=True,
        photo=False,
    ):
        """Basic constructor call super and set the occupied and photo attributes"""
        super(RadicalOrbitalFactory, self).__init__(atomObjList, allSymmetry, allKekule)
        self.occupied = occupied
        self.photo = photo

    @staticmethod
    def prepareChainOrbital(orbital, parentOrb, occupied, photo):
        """Prepare a chain orbital.  Evaluate possibility and handle bookkeeping"""
        chainOrbital = None

        hasExtendedPotential = (
            parentOrb.type == "sigma" and (not occupied) and orbital.isLoneRadical()
        )

        if hasExtendedPotential:
            parentOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(parentOrb)
            candidateOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(orbital)
            if len(parentOrbIdxSet.intersection(candidateOrbIdxSet)) == 0:
                ## No overlap, looks good.
                chainOrbital = parentOrb.simpleCopy()
                chainOrbital.extOrbital = orbital
        return chainOrbital

    def proposeChainedOrbitals(self, orb):
        """Method to handle details of proposing chained orbitals

        We need to have chained orbitals for Disproportionation reactions
        (For now only for sink orbitals)
        """
        seenSymSet = set([])
        if not self.occupied and orb.type == "sigma":
            for atom in orb.neighbor.GetAtoms():
                if atom == orb.atom or atom.GetSymmetryClass() in seenSymSet:
                    continue
                seenSymSet.add(atom.GetSymmetryClass())
                for potChainOrb in radicalOrbitalIter(atom, useSymmetries=False):
                    chainOrb = self.prepareChainOrbital(
                        potChainOrb, orb, self.occupied, self.photo
                    )
                    if chainOrb is not None:
                        yield chainOrb
        return

    @staticmethod
    def isValidOrbitalOccupation(orbital, occupied, photo):
        """Verify that the orbital is actually a valid orbital for Radical Proposal purposes .

        Here is where the basic rules outlined in the module level docstrings come into play.

        Basic Ideas:
        If occupied:
         - Sigma only on non-C,H atoms, in 3-rings, or adjacent to a radical
         - pi only if photoexcited.
         - If is a free radical (or a carbene, then let go)
        else:
         - Sigma only on non-C,H atoms, in 3-rings, or adjacent to a radical
         - pi all the time
         - or free radical
        """
        accept = True
        if occupied:
            if orbital.type == "pi" and not photo:
                accept = False
            elif orbital.type == "sigma":
                accept = (
                    (
                        orbital.atom.GetAtomicNum() not in set([1, 6])
                        and orbital.neighbor.GetAtomicNum() not in set([1])
                    )
                    or PERIODIC_ROW[orbital.atom.GetAtomicNum()] > 2
                    or orbitalInfo(orbital.neighbor)["nRadicals"] > 0
                    or OEBondIsInRingSize(orbital.atom.GetBond(orbital.neighbor), 3)
                )
        else:
            if orbital.type == "sigma":
                ## Allow H-based abstractions, but not C in general
                accept = (
                    (
                        orbital.atom.GetAtomicNum() not in set([6])
                        and orbital.neighbor.GetAtomicNum() not in set([1])
                    )
                    or PERIODIC_ROW[orbital.atom.GetAtomicNum()] > 2
                    or orbitalInfo(orbital.neighbor)["nRadicals"] > 0
                    or OEBondIsInRingSize(orbital.atom.GetBond(orbital.neighbor), 3)
                )
                ## Handle one last case - Allow any sigma bond if any neighbor of the atom
                ## has a radical
                if not accept:
                    for atom in orbital.atom.GetAtoms():
                        if orbitalInfo(atom)["nRadicals"] > 0:
                            accept = True
                            break
        return accept

    def __iter__(self):
        """Workhorse. Iterate through atoms, then orbitals, check validity, prep any chains, yield"""
        for oeAtom, atomObj in self.atomIterator():
            # log.info('Looking at atom: %s' % atomObj)
            for orb in radicalOrbitalIter(oeAtom, useSymmetries=(not self.allSymmetry)):
                # log.info('Potential Orbital :%s' % pformat(orb.toLabeledSmiAndInfoStr(10)))
                if self.isValidOrbitalOccupation(orb, self.occupied, self.photo):
                    yield orb, atomObj
                    ## Look for chains (None for radicals?)
                    for chainOrb in self.proposeChainedOrbitals(orb):
                        yield chainOrb, atomObj
