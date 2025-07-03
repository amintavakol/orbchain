from chemutils.CombiCDB.proposers.AtomBasedReactionProposer import (
    AtomBasedReactionProposer,
)
from chemutils.Common.OrbitalModel import Orbital

from chemutils.CombiCDB.proposers import filters

defaultFilters = [
    filters.NonReasonableOrbitalsFilter(),
    filters.ExtremeChargeFilter(),
    filters.SubstitutionHybridizationFilter(),
    filters.ExtremeBondFilter(),
    filters.PiBondBondDissociationFilter(),
    filters.PossiblePericyclicFilter(),
]

from chemutils.CombiCDB.proposers.PeriPattern import defaultPeriPatterns


class PericyclicProposer(AtomBasedReactionProposer):
    """Reaction proposer"""

    def __init__(
        self, closureList, orbitalPairFilters=defaultFilters, periPatternObjs=None
    ):
        """Constructor.  Setup the necessary pieces to propose"""
        self.periPatterns = periPatternObjs
        self.closureList = closureList
        super(PericyclicProposer, self).__init__(
            closureList, closureList, orbitalPairFilters, checkResonance=False
        )

    def setOrbitalFactories(self):
        """Must implement for the base class, use this to setup the default periPatterns"""
        if self.periPatterns is None:
            self.periPatterns = defaultPeriPatterns

    def proposeOrbitalPairs(self):
        """Complete re-implementation of this function to handle special case of pericyclics
        Loop over the patterns, loop over orbs from patterns.

        Instead of atoms, we will yield the pattern that matched.
        """
        for pattObj in self.periPatterns:
            for srcOrb, sinkOrb in pattObj.proposeOrbPairs(self.closureList):
                ## We assume here that the patterns themselves handle the inter/intra
                ## differences.
                (
                    (clonedSrcOrb, clonedSinkOrb),
                    copyMol,
                ) = Orbital.cloneOrbitalsWithCommonMol([srcOrb, sinkOrb])
                yield (clonedSrcOrb, clonedSinkOrb, pattObj, None)
