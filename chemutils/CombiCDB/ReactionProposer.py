from openeye.oechem import OEGraphMol, OEParseSmiles
from openeye.oechem import OEAddMols
from openeye.oechem import OEAddExplicitHydrogens
from openeye.oechem import OESubSearch

from chemutils.Common.Util import (
    createStandardSmiString,
)
from chemutils.Common.OrbitalModel import Orbital, orbitalIter
from chemutils.Common.OrbitalModel import ORBITAL_LABEL_INCREMENT
from chemutils.Common.OrbitalModel import atomHybridization
from chemutils.CombiCDB.OrbitalInteraction import (
    moveOrbitalElectrons,
    undoMoveOrbitalElectrons,
)
from chemutils.CombiCDB.OrbitalInteraction import (
    compatibleOrbitalPair,
)
from chemutils.score.orbital.ReactivityScore import ReactiveOrbitalFactory
from chemutils.score.orbital.ReactivityScore import (
    OrbitalReactivityProbability,
    AcceptAllOrbitalReactivityProbability,
)
from chemutils.CombiCDB.MechanismModel import clearMechanismLabels



class BaseReactionProposer:
    """Proposes possible reaction steps in a simulation, acting as the
    "hypothesis generation" phase.  Given some a starting reaction mixture and
    conditions, yield all possible consistent orbital interaction combinations that
    could be treated as possible elementary reaction steps.
    It will be up to a separate module to further score and filter down these
    proposals to the reaction steps that are actually most favorable and likely.
    """

    # List of filters to eliminate undesired orbital pair proposals.
    # Set to an empty list to not filter anything.  If set to None, a default set will be used.
    orbitalPairFilters = None
    reactionProposerParams = None

    def __init__(self, reactionProposerParams=None):
        """Default constructor, but this is supposed to be an abstract base class"""
        self.reactionProposerParams = reactionProposerParams
        pass

    def getInstance(reactionProposerParams=None):
        instance = ReactionProposer(reactionProposerParams)
        return instance

    getInstance = staticmethod(getInstance)

    def getOrbitalPairFilters(self):
        """Get the filters to exclude certain orbital pair proposals,
        or a default set if none are specified.
        Sub-classes can override to alter the default set
        """
        if self.orbitalPairFilters is None:
            # Compose a default set of filters.
            self.orbitalPairFilters = list()
            # Default to empty list
        return self.orbitalPairFilters

    def __call__(self, reactionMixture, analyzedMoleculeSmiSet=None):
        """General use iterator function to yield
        the proposed orbital pairs.  Additional function
        of looking for any filters to eliminated undesired proposals.

        # TO DO:  Add extra step to standardize any proposals (for consistency against resonance structures)
        """
        for filledOrb, unfilledOrb in self.proposeOrbitalPairs(
            reactionMixture, analyzedMoleculeSmiSet
        ):
            if self.passOrbitalPairFilters(filledOrb, unfilledOrb):
                yield (filledOrb, unfilledOrb)

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        """Primary abstract method that subclasses should override
        to generate the collection of orbital pair reaction proposals.
        """
        raise NotImplementedError()

    def passOrbitalPairFilters(self, filledOrb, unfilledOrb):
        """Determine whether the proposed orbital
        pair passes all of the filters set on this object.
        """
        passFilters = True
        for orbitalPairFilter in self.getOrbitalPairFilters():
            # log.debug('*******FILTER*******: About to Call %s' % orbitalPairFilter.__class__.__name__ )
            if not orbitalPairFilter(filledOrb, unfilledOrb):
                passFilters = False
                break
                # Don't need to check anymore if one already failed
        return passFilters

    def checkAnalysisHistory(reactionMixture):
        """Check the elementary step history for reactants that have been observed already,
        to avoid repeating previous orbital pair analysis.
        Assume that all of the reactive combinations that include these have been pre-analyzed.
        """
        analyzedMoleculeSmiSet = set()
        for (
            reactionSmiles,
            elementaryStepData,
        ) in reactionMixture.elementaryStepHistory.iteritems():
            for reactantSmi in elementaryStepData["reactantSmilesList"]:
                analyzedMoleculeSmiSet.add(reactantSmi)

        # If there are no new species in the mixture, not worth checking anything, must have all been pre-analyzed
        newSpeciesFound = False
        for smiles in reactionMixture:
            if smiles not in analyzedMoleculeSmiSet:
                newSpeciesFound = True

        return (analyzedMoleculeSmiSet, newSpeciesFound)

    checkAnalysisHistory = staticmethod(checkAnalysisHistory)

    def isNewReactant(mol, analyzedMoleculeSmiSet):
        """Check that we have not already analyzed this molecule.
        The analyzedMoleculeSmiSet should contain the set of SMILES string representing
        the reactants that have already been analyzed before.
        """
        molCopy = OEGraphMol(mol)
        # Make a copy we're free to modify
        clearMechanismLabels(molCopy)
        molSmi = createStandardSmiString(molCopy)
        molCopy.Clear()

        previouslyAnalyzed = molSmi in analyzedMoleculeSmiSet

        return not previouslyAnalyzed

    isNewReactant = staticmethod(isNewReactant)

    def createStandardOrbitalPair(filledOrb, unfilledOrb):
        """Given a pair of orbitals that represents the interaction in
        an elementary reaction step, generate another orbital pair
        for an equivalent but standard representation.

        This is of particular relevance when the product can be drawn
        with multiple equivalent resonance structures, we want to
        consistently represent different equivalent orbital interactions
        with a reproducible standard representation.  This need not be
        the "best" representation, but it must be consistently reproduced
        for every possible variable but equivalent input.

        Assumes the starting reactant molecule(s) the input orbitals are
        based upon are already in a standard resonance structure form,
        perhaps by chemutils.CombiCDB.ResonanceStructureUtil.standardizeResonanceStructure.
        Thus, we're only worried about standardizing the resonance structure
        form of the product.
        """
        return (filledOrb, unfilledOrb)

    createStandardOrbitalPair = staticmethod(createStandardOrbitalPair)


class CompositeReactionProposer(BaseReactionProposer):
    """General proposer object that just composites the results
    of many component proposers.
    """

    reactionProposers = None

    def __init__(self, reactionProposerParams=None):
        """Default constructor."""
        BaseReactionProposer.__init__(self, reactionProposerParams)
        self.reactionProposers = None
        # Start with empty list of possibilities

    def addReactionProposer(self, proposer):
        if self.reactionProposers is None:
            self.reactionProposers = list()
        self.reactionProposers.append(proposer)

    def setupDefaultReactionProposers(self):
        """Initialize the composite proposer
        based on default set of standard ones.
        Should be overriden by sub-classes
        """
        raise NotImplementedError()

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        """General which just sequentially compiles the results
        of all the sub proposers.
        """
        # Make sure we at least have some default proposers setup
        if self.reactionProposers is None:
            self.setupDefaultReactionProposers()

        # Check the elementary step history to see if we actually have any new molecules to analyze.
        (analyzedMoleculeSmiSet, newSpeciesFound) = self.checkAnalysisHistory(
            reactionMixture
        )

        if newSpeciesFound:
            for reactionProposer in self.reactionProposers:
                for orbitalPair in reactionProposer(
                    reactionMixture, analyzedMoleculeSmiSet
                ):
                    yield orbitalPair


class ReactionProposer(CompositeReactionProposer):
    """Primary proposer object for actual use,
    composites results from several different proposer objects to generate
    many different possibilities.

    # TO DO: Add special cases / patterns of pairs of orbitals that should be considered,
    # Bond dissociations (Sn1 / E1, free radical initiation)
    # Brominium ion formation / cyclization, similar to carbene cyclopropanation
    # Pericyclic reactions
    #   Cycloadditions (+retro)
    #       [4+2] mostly
    #       [3+2] dipolar
    #       [2+2] some, but uncommon, may need photochemical activation?
    #   Electrocyclizations (+retro), 6 and 4 ring
    #   Sigmatropic rearrangements (mostly look for [3,3] pattern *=****=*)
    #       1,5 hydrogen sigmatropic rearrangement
    #       1,3 alkyl sigmatropic rearrangement (very minor, probably not worth it)
    #       [2,3] sigmatropic rearrangement
    # Peroxyacid alkene epoxidation?
    # Sigma chains.  Retro-Aldol, Grob fragmentation, etc.  Like an E2 elimination, but the
    #   start of the target orbital chain is a C-C (or maybe Si-C) bond, rather than H-C.
    #   Generally only effective if the source filled orbital is directly adjacent to the target,
    #   otherwise requires too many concerted bond angles.
    """

    def getOrbitalPairFilters(self):
        """Get the filters to exclude certain orbital pair proposals,
        or a default set if none are specified.
        """
        if self.orbitalPairFilters is None:
            # Compose a default set of filters.
            self.orbitalPairFilters = list()
            # Default to empty list
            # Actually add on one about non-overlap.
            self.orbitalPairFilters.append(OverlapAtomsFilter())
            self.orbitalPairFilters.append(ExtremeChargeFilter())
            self.orbitalPairFilters.append(SubstitutionHybridizationFilter())
            self.orbitalPairFilters.append(ExtremeBondFilter())
        return self.orbitalPairFilters

    def setupDefaultReactionProposers(self):
        """Initialize the composite proposer
        based on default set of standard ones.
        """
        self.reactionProposers = list()

        proposer = OrbitalPairCombinationProposer(self.reactionProposerParams)
        self.reactionProposers.append(proposer)

        proposer = LocalInteractionProposer(self.reactionProposerParams)
        # self.reactionProposers.append( proposer );

        proposer = BondDissociationProposer(self.reactionProposerParams)
        self.reactionProposers.append(proposer)

        # proposer = PericyclicReactionProposer();
        # self.reactionProposers.append( proposer );


class BaseFilledUnfilledProposer(BaseReactionProposer):
    """Abstract class for reaction proposers that are
    based on combinations involving a pre-screened
    collection of possible filled or unfilled orbitals.
    """

    def __init__(self, reactionProposerParams=None):
        """Default constructor."""
        BaseReactionProposer.__init__(self, reactionProposerParams)

        # Option whether to accept all orbitals in mixtures as possibilities, or pre-filter them by a reactivity ranking
        self.acceptAllOrbitals = False
        # Generators of possibly reactive filled and unfilled orbitals.
        # Default to None so the class will fill in some standard ones, but allow caller to override here
        self.filledOrbGenerator = None
        self.unfilledOrbGenerator = None

        if self.reactionProposerParams is not None:
            self.acceptAllOrbitals = self.reactionProposerParams.acceptAllOrbitals

        # self.restrictSymmetryFilled = True;
        # self.restrictSymmetryUnfilled = False;

    def getFilledOrbitalGenerator(self, reactionMixture, restrictSymmetry=True):
        """Get the (reactive) filled orbital generator set in the object,
        or prepare a standard, default one if no specific one is set
        If acceptAll option is set, the default one generated will simply accept
        all generated orbitals instead of filtering down to those that appear to be "reactive."
        """
        if self.filledOrbGenerator is None:
            filledOrbReactivityCalc = None
            if self.acceptAllOrbitals:
                filledOrbReactivityCalc = AcceptAllOrbitalReactivityProbability(
                    True, reactionMixture
                )
            else:
                filledOrbReactivityCalc = OrbitalReactivityProbability(
                    True, reactionMixture
                )
            self.filledOrbGenerator = ReactionMixtureOrbitalFactory(
                reactionMixture, filledOrbReactivityCalc, restrictSymmetry
            )
        return self.filledOrbGenerator

    def getUnfilledOrbitalGenerator(self, reactionMixture, restrictSymmetry=False):
        """Get the (reactive) unfilled orbital generator set in the object,
        or prepare a standard, default one if no specific one is set
        If acceptAll option is set, the default one generated will simply accept
        all generated orbitals instead of filtering down to those that appear to be "reactive."
        """
        if self.unfilledOrbGenerator is None:
            unfilledOrbReactivityCalc = None
            if self.acceptAllOrbitals:
                unfilledOrbReactivityCalc = AcceptAllOrbitalReactivityProbability(
                    False, reactionMixture
                )
            else:
                unfilledOrbReactivityCalc = OrbitalReactivityProbability(
                    False, reactionMixture
                )
            self.unfilledOrbGenerator = ReactionMixtureOrbitalFactory(
                reactionMixture, unfilledOrbReactivityCalc, restrictSymmetry
            )
        return self.unfilledOrbGenerator


class OrbitalPairCombinationProposer(BaseFilledUnfilledProposer):
    """Major reaction proposer that takes a collection
    of possibly reactive filled and unfilled orbitals
    and proposed all pair-wise combinations.

    Needs to take into consideration
    the possibility of the filled and unfilled orbitals coming from
    different molecule objects, in which case, need to merge the molecule
    objects into a common composite molecule that a fresh set of
    orbitals will properly refer to.

    If the filled and unfilled orbitals come from the same molecule,
    should be ready to interpret this in two ways.  Indeed as a single
    molecule for an intramolecular reaction, or reaction with a second
    instance of the molecule for an intermolecular reaction.
    """

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        """Given a list (or generator) of filled and unfilled orbitals,
        return all pair-wise combinations.
        """
        if analyzedMoleculeSmiSet is None:
            (analyzedMoleculeSmiSet, newSpeciesFound) = self.checkAnalysisHistory(
                reactionMixture
            )

        # Prepare generators for all of the separate filled and unfilled orbitals we want to consider in combination
        filledOrbGenerator = self.getFilledOrbitalGenerator(
            reactionMixture, restrictSymmetry=False
        )
        unfilledOrbGenerator = self.getUnfilledOrbitalGenerator(
            reactionMixture, restrictSymmetry=True
        )

        # Pre-iterate through the orbitals and record them in memory in a list, otherwise
        #   will keep regenerating them multiple times excessively in a nested loop
        unfilledOrbList = []
        for unfilledOrb in unfilledOrbGenerator:
            unfilledOrb.dataTags["fromNewReactant"] = self.isNewReactant(
                unfilledOrb.mol, analyzedMoleculeSmiSet
            )

            unfilledOrbList.append(unfilledOrb)

        for filledOrb in filledOrbGenerator:
            filledOrb.dataTags["fromNewReactant"] = self.isNewReactant(
                filledOrb.mol, analyzedMoleculeSmiSet
            )

            for unfilledOrb in unfilledOrbList:
                isUnfilledOrbNew = self.isNewReactant(
                    unfilledOrb.mol, analyzedMoleculeSmiSet
                )

                if (
                    not filledOrb.dataTags["fromNewReactant"]
                    and not unfilledOrb.dataTags["fromNewReactant"]
                ):
                    # Neither orbital comes from a "new" reactant.  Don't bother analyzing,
                    #   as this must be redundant to a previous simulation iteration's analysis.
                    continue

                if filledOrb.mol == unfilledOrb.mol:
                    # Orbitals from the same molecule object.
                    # Should be using the "is" comparator, but that is returning False
                    #   even when working with the same molecule object for some reason.
                    # Interpret this in 2 ways, as indeed an intra-molecular reaction,
                    #   or prepare separate molecule copies for an inter-molecular reaction.
                    # Should already work as intra-molecular example, but prepare a copy
                    #   so caller is free to manipulate it without disrupting the object contents

                    # Beware, if on the same molecule, these might be based on overlapping
                    #   atoms, which would not make sense to "react" them together.

                    if compatibleOrbitalPair(filledOrb, unfilledOrb):
                        (
                            (cloneFilledOrb, cloneUnfilledOrb),
                            copyMol,
                        ) = Orbital.cloneOrbitalsWithCommonMol([filledOrb, unfilledOrb])
                        cloneFilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 1)
                        cloneUnfilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 2)

                        # If we have gotten to this point, and the cloned orbs are still exactly the same,
                        # then it doesn't make sense to 'react' them together
                        if (
                            cloneFilledOrb.toLabeledSmiAndInfoStr()
                            != cloneUnfilledOrb.toLabeledSmiAndInfoStr()
                        ):
                            yield (cloneFilledOrb, cloneUnfilledOrb)

                    # Prepare copies to yield inter-molecular example
                    # Don't need to do it here, just reuse the usual 2 molecule case,
                    #   which will work out equivalently
                    pass

                # Orbitals from distinct molecule objects.  Prepare a composite
                #   molecule object that both orbitals will be based upon
                (
                    (compositeFilledOrb, compositeUnfilledOrb),
                    compositeMol,
                ) = Orbital.compositeOrbitalsFromDistinctMols([filledOrb, unfilledOrb])

                yield (compositeFilledOrb, compositeUnfilledOrb)


class LocalInteractionProposer(BaseFilledUnfilledProposer):
    """Should mostly produce similar results as the Orbital Combination proposer,
    but this specifically looks at orbitals combinations that are "local" / directly
    adjacent to each other.  These can be relevant when one of the interacting
    orbitals alone may seem to unreactive to be important, but the fact that it's
    so close to a complementary reactive orbital makes a reaction possible.

    For example:  Carbocation rearrangements, collapse of carbonyl tetrahedral intermediates,
    retro-aldol, Grob fragementations.
    """

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        """General use function to propose orbital interaction pairs
        as possible elementary steps from the given reaction mixture,
        using default orbital generators and selection criteria.
        """
        if analyzedMoleculeSmiSet is None:
            (analyzedMoleculeSmiSet, newSpeciesFound) = self.checkAnalysisHistory(
                reactionMixture
            )

        # Prepare generators for all of the separate filled and unfilled orbitals we want to consider in combination
        filledOrbGenerator = self.getFilledOrbitalGenerator(
            reactionMixture, restrictSymmetry=True
        )
        unfilledOrbGenerator = self.getUnfilledOrbitalGenerator(
            reactionMixture, restrictSymmetry=True
        )

        # Look for local interactions based on the "reactive" filled orbitals
        for filledOrb in filledOrbGenerator:
            filledOrb.dataTags["fromNewReactant"] = self.isNewReactant(
                filledOrb.mol, analyzedMoleculeSmiSet
            )

            if filledOrb.dataTags["fromNewReactant"]:
                # Check for local interactions with this orbital that may be missed in the main loop
                for (
                    localFilledOrb,
                    localUnfilledOrb,
                ) in self.prepareFilledLocalOrbitalPairs(filledOrb):
                    yield (localFilledOrb, localUnfilledOrb)

        # Look for local interactions based on the "reactive" unfilled orbitals
        for unfilledOrb in unfilledOrbGenerator:
            unfilledOrb.dataTags["fromNewReactant"] = self.isNewReactant(
                unfilledOrb.mol, analyzedMoleculeSmiSet
            )

            if unfilledOrb.dataTags["fromNewReactant"]:
                # Check for local interactions with this orbital that may be missed in the main loop
                for (
                    localFilledOrb,
                    localUnfilledOrb,
                ) in self.prepareUnfilledLocalOrbitalPairs(unfilledOrb):
                    yield (localFilledOrb, localUnfilledOrb)

    def prepareFilledLocalOrbitalPairs(self, filledOrb):
        """Given a filled orbital, prepare potential orbital interaction pairs with immediately
        adjacent (local) unfilled orbitals.  This is to help identify orbital interaction
        pairs that may be missed because the unfilled orbitals alone are rated as unlikely
        to react, when in fact they may react, but only if that happen to be immediately
        adjacent to a reactive filled orbital.

        e.g., Retro-aldol, Grob fragmentations?  These are actually even tricker because
        the "local" unfilled orbital should also be an extended chain.

        e.g., Elimination step of nucleophilic acylation to collapse a tetrahedral intermediate.
        Arguably the leaving group is at a tertiary center, making it unreactive, but actually
        this should collapse very easily.
        """

        # Ensure that none of the proposed orbital atoms overlap
        filledOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(filledOrb)

        for neighbor in filledOrb.atom.GetAtoms():
            for neighborOrb in orbitalIter(neighbor):
                hasInteractionPotential = not neighborOrb.isLonePair()
                # Anything but a lone pair could be interpreted as an unfilled orbital

                unfilledOrb = neighborOrb

                if hasInteractionPotential:
                    unfilledOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(
                        unfilledOrb
                    )

                    if len(filledOrbIdxSet.intersection(unfilledOrbIdxSet)) > 0:
                        # Orbitals have overlapping atoms, does not make sense for interaction
                        #   Beware of confounding cases with pericyclic reactions and multi-bond forming cycles
                        #   (e.g., bromonium ion formation).
                        pass
                    else:
                        # Yield clones so the caller is free to modify them without causing consistency issues
                        (
                            (cloneFilledOrb, cloneUnfilledOrb),
                            copyMol,
                        ) = Orbital.cloneOrbitalsWithCommonMol([filledOrb, unfilledOrb])

                        cloneFilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 1)
                        cloneUnfilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 2)

                        yield (cloneFilledOrb, cloneUnfilledOrb)

    def prepareUnfilledLocalOrbitalPairs(self, unfilledOrb):
        """Given an unfilled orbital, prepare potential orbital interaction pairs with immediately
        adjacent (local) filled orbitals.  This is to help identify orbital interaction
        pairs that may be missed because the filled orbitals alone are rated as unlikely
        to react, when in fact they may react, but only if that happen to be immediately
        adjacent to a reactive unfilled orbital.

        e.g., All 1,2 rearrangements (e.g., carbocations).  The C-H or C-C migrating bond
            can be a source when a good sink like an empty p orb is directly adjacent,
            but otherwise these  C-H and C-C bonds should be essentially unreactive.
        """
        # Ensure that none of the proposed orbital atoms overlap
        unfilledOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(unfilledOrb)

        for neighbor in unfilledOrb.atom.GetAtoms():
            for neighborOrb in orbitalIter(neighbor):
                hasInteractionPotential = not neighborOrb.isEmpty()
                # Anything but a an empty orbital could be interpreted as a filled orbital

                filledOrb = neighborOrb

                if hasInteractionPotential:
                    filledOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(
                        filledOrb
                    )

                    if len(filledOrbIdxSet.intersection(unfilledOrbIdxSet)) > 0:
                        # Orbitals have overlapping atoms, does not make sense for interaction
                        #   Beware of confounding cases with pericyclic reactions and multi-bond forming cycles
                        #   (e.g., bromonium ion formation).
                        pass
                    else:
                        # Flip bond orbitals to represent typical 1,2 rearrangements
                        if filledOrb.isBondOrbital():
                            filledOrb = filledOrb.flippedBondOrbital()

                        # Yield clones so the caller is free to modify them without causing consistency issues
                        (
                            (cloneFilledOrb, cloneUnfilledOrb),
                            copyMol,
                        ) = Orbital.cloneOrbitalsWithCommonMol([filledOrb, unfilledOrb])

                        cloneFilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 1)
                        cloneUnfilledOrb.labelOrbitalAtoms(ORBITAL_LABEL_INCREMENT * 2)

                        yield (cloneFilledOrb, cloneUnfilledOrb)


class BondDissociationProposer(BaseReactionProposer):
    """Proposal for bond-dissociation reactions that will be missed
    by the general orbital combination proposers.  These are unusual
    in that the the source filled orbital could be considered as identical
    to the target unfilled orbital, or perhaps that there is no filled / unfilled
    orbital at all.

    TO DO: Implement this for both hetero and homo-lytic bond dissociation
    """

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        raise StopIteration


class PericyclicReactionProposer(CompositeReactionProposer):
    """Proposal for pericyclic reactions that will be missed
    by the general orbital combination proposers.  These are unusual
    in that the electrons will move in a closed loop, so it makes no
    sense, or at least is completely arbitrary, what to define
    as the filled vs. unfilled orbital halves of the reaction.

    Further break this down into pericylcic reactions based on
    single or multiple (two) molecules, as the latter will
    affect how the combinatorial possibilities must be assessed.
    """

    def setupDefaultReactionProposers(self):
        """Initialize the composite proposer
        based on default set of standard ones.
        """
        self.reactionProposers = list()

        proposer = PericyclicUnimolecularProposer()
        self.reactionProposers.append(proposer)

        proposer = PericyclicBimolecularProposer()
        self.reactionProposers.append(proposer)


class PericyclicUnimolecularProposer(BaseReactionProposer):
    """Propose uni-molecular pericyclic reactions.
    Should make it simpler to not have to look for combinatorial
    possibilities, just particular rearrangement patterns.

    Knowledge based on patternMatcherLists.
       Contains sub-lists, containing tuples, each containing
        - SMARTS string
        - OESubSearch matcher object
        - Source orbital info string
        - Target orbital info string
        - Priority Value
       Multiple tuples per sub-list (patternMatcherList) to reflect varying
       degrees of specificity for the same kind of reactions.  Should be listed
       in order of decreasing specificity so the less specific ones will
       only be used if a more specific match was not already found.
       Ordering specified by the priority value, with option of equal priority
       values to reflect pattern variants of equal likelihood.
    """

    patternMatcherLists = None

    def addPatternList(self, patternInfoTuples):
        """Add a (related) list of pattern info to identify reaction proposers.
        List of pattern info tuples containing 3-ple of a smartsPattern, sourceOrbInfoStr and targetOrbInfoStr.
        """
        if self.patternMatcherLists is None:
            self.patternMatcherLists = list()

        patternMatcherList = list()
        for (
            priority,
            smartsPattern,
            sourceOrbInfoStr,
            targetOrbInfoStr,
        ) in patternInfoTuples:
            # Only want to initialize the pattern matchers once.  Relatively expensive process.
            matcher = OESubSearch()
            matcher.Init(smartsPattern)

            patternMatcherList.append(
                (priority, smartsPattern, matcher, sourceOrbInfoStr, targetOrbInfoStr)
            )

        # Ensure descending order by priority value
        patternMatcherList.sort()
        patternMatcherList.reverse()

        self.patternMatcherLists.append(patternMatcherList)

    def setupDefaultPatternMatchers(self):
        """Initialize a default set of pattern matchers."""
        # Reset to an empty list
        self.patternMatcherLists = list()

        # May look strange to terminate pericyclic reactions on a p orbital instead of a pi orbital,
        #   but part of the artificial way we have to manage pericyclics since the pi orbital you'd want
        #   to target would have already moved on by the opposite end of the cyclic chain of electron flow

        # Electrocyclization
        # 6 pi
        self.addPatternList(
            [
                (
                    100,
                    "[*:1](/[*:10])([*:11])=[*:2]/[*:3]=[*:4]\[*:5]=[*:6]([*:60])(/[*:61])",
                    "1 pi 2 2 2,10,11; 3 pi 4 2; 5 pi 6 2 None 5,60,61",
                    "6 p None 0 5,60,61",
                ),
                (
                    0,
                    "[*:1]=[*:2][*:3]=[*:4][*:5]=[*:6]",
                    "1 pi 2 2; 3 pi 4 2; 5 pi 6 2",
                    "6 p None 0",
                ),
            ]
        )

        # Sigmatropic rearrangements
        # [3,3]
        self.addPatternList(
            [
                (
                    0,
                    "[*:1]=[*:2][*:3]-[*:4][*:5]=[*:6]",
                    "1 pi 2 2; 3 sigma 4 2; 5 pi 6 2",
                    "6 p None 0",
                )
            ]
        )

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        # If caller has not custom-overriden the matcher list, setup defaults
        if self.patternMatcherLists is None:
            self.setupDefaultPatternMatchers()

        for smiles, mixtureComponent in reactionMixture.iteritems():
            # Need a fresh mol object for back / parent pointers on orbitals
            mol = OEGraphMol()
            OEParseSmiles(mol, smiles)

            # SMARTS patterns looking for hydrogens (evan as wildcard "*" atoms)
            #   will only match explicit hydrogen atoms, not implicit ones
            OEAddExplicitHydrogens(mol)

            if mixtureComponent.currentQuantity() > 0 and self.isNewReactant(
                mol, analyzedMoleculeSmiSet
            ):
                # Look for pattern matching on any of the mixture components
                #   but track redundancies to preferentially use higher priority (more specific) patterns
                priorityByRedundancyKey = dict()
                for patternMatcherList in self.patternMatcherLists:
                    for (
                        priority,
                        smartsPattern,
                        matcher,
                        sourceOrbInfoStr,
                        targetOrbInfoStr,
                    ) in patternMatcherList:
                        # Every match for a given mixture component
                        for match in matcher.Match(mol, False):
                            # For every atom that was matched, apply an atom map label so we can build orbitals off of it
                            for atomMatchPair in match.GetAtoms():
                                atomMatchPair.target.SetMapIdx(
                                    atomMatchPair.pattern.GetMapIdx()
                                )

                            sourceOrb = Orbital.fromLabeledMolAndInfoStr(
                                mol, sourceOrbInfoStr
                            )
                            targetOrb = Orbital.fromLabeledMolAndInfoStr(
                                mol, targetOrbInfoStr
                            )

                            # Wipe out atom map indexes to revert the molecule back to its original state
                            for atomMatchPair in match.GetAtoms():
                                atomMatchPair.target.SetMapIdx(0)

                            # Make sure we have not already proposed a more specific version of this proposal
                            redundancyKey = self.orbitalInteractionRedundancyKey(
                                sourceOrb, targetOrb
                            )

                            redundantPriority = None
                            if redundancyKey in priorityByRedundancyKey:
                                redundantPriority = priorityByRedundancyKey[
                                    redundancyKey
                                ]

                            if (
                                redundantPriority is None
                                or priority >= redundantPriority
                            ):
                                # Pattern is not redundant or it is priority ranked as high or higher than any previous match
                                priorityByRedundancyKey[redundancyKey] = priority

                                # Create clones of the orbitals and molecule so the caller is free to modify them
                                (
                                    (cloneSourceOrb, cloneTargetOrb),
                                    copyMol,
                                ) = Orbital.cloneOrbitalsWithCommonMol(
                                    [sourceOrb, targetOrb]
                                )

                                yield (cloneSourceOrb, cloneTargetOrb)

    def orbitalInteractionRedundancyKey(self, sourceOrb, targetOrb):
        """Come up with a unique key that represents the relevant information in
        proposed orbital pair combination that should be used to differentiate
        redundant proposals that likely just represent more or less specific
        versions of the same reaction.

        Usually the difference is stereochemistry, so we don't want to let orbital
        facial selectivity affect this redundancy key, and we can't count on consistent
        atom mapping indexes, so mainly just the use underlying atom index numbers
        on the molecule object that should remain constant across pattern matches.
        """
        sourceAtomIndexes = list()
        currentOrb = sourceOrb
        while currentOrb is not None:
            sourceAtomIndexes.append(currentOrb.atom.GetIdx())
            if currentOrb.neighbor is not None:
                sourceAtomIndexes.append(currentOrb.neighbor.GetIdx())
            currentOrb = currentOrb.extOrbital

        targetAtomIndexes = list()
        currentOrb = targetOrb
        while currentOrb is not None:
            targetAtomIndexes.append(currentOrb.atom.GetIdx())
            if currentOrb.neighbor is not None:
                targetAtomIndexes.append(currentOrb.neighbor.GetIdx())
            currentOrb = currentOrb.extOrbital

        redundancyKey = (tuple(sourceAtomIndexes), tuple(targetAtomIndexes))

        return redundancyKey


class PericyclicBimolecularProposer(BaseReactionProposer):
    """Propose bi-molecular pericyclic reactions, cycloadditions in particular.
    Must be handled similar to orbital combinations proposer in that we must
    look for reactive site combinations on multiple molecules, multiple instances
    of the same molecule and also intramolecular reactions.
    Knowledge based on patternMatcherListPairs with tuples.
    Pairs of lists to reflect separate items for source vs. target orbitals.
    Lists contain multiple tuples, each with 4 items to describe the orbital pattern:
    - SMARTS string
    - OESubSearch matcher object
    - Orbital info string
    - Priority Value
    Multiple tuples per sub-list to reflect varying
    degrees of specificity for the same kind of reactions.  Should be listed
    in order of decreasing specificity so the less specific ones will
    only be used if a more specific match was not already found.
    Ordering specified by the priority value, with option of equal priority
    values to reflect pattern variants of equal likelihood.

    Separate lists should indicate which map indexes in the source and target
    patterns represent the key reacting atoms that should be counted
    when looking for redundant matches.
    """

    patternMatcherListPairs = None

    def addPatternListPair(
        self, sourceTuples, targetTuples, sourceKeyMapIndexes, targetKeyMapIndexes
    ):
        """Add a pair of pattern lists to identify reaction proposers."""
        if self.patternMatcherListPairs is None:
            self.patternMatcherListPairs = list()

        # Setup source and target matcher lists, but do so only once per pattern.
        sourcePatternMatcherList = list()
        for priority, smartsPattern, orbInfoStr in sourceTuples:
            matcher = OESubSearch()
            matcher.Init(smartsPattern)
            sourcePatternMatcherList.append(
                (priority, smartsPattern, matcher, orbInfoStr)
            )

        targetPatternMatcherList = list()
        for priority, smartsPattern, orbInfoStr in targetTuples:
            matcher = OESubSearch()
            matcher.Init(smartsPattern)
            targetPatternMatcherList.append(
                (priority, smartsPattern, matcher, orbInfoStr)
            )

        # Ensure descending order by priority values
        sourcePatternMatcherList.sort()
        sourcePatternMatcherList.reverse()

        targetPatternMatcherList.sort()
        targetPatternMatcherList.reverse()

        # Add completed matcher lists to the master object attribute
        self.patternMatcherListPairs.append(
            (
                sourcePatternMatcherList,
                targetPatternMatcherList,
                sourceKeyMapIndexes,
                targetKeyMapIndexes,
            )
        )

    def setupDefaultPatternMatchers(self):
        """Initialize a default set of pattern matchers."""
        # Reset to an empty list
        self.patternMatcherListPairs = list()

        # May look strange to terminate pericyclic reactions on a p orbital instead of a pi orbital,
        #   but part of the artificial way we have to manage pericyclics since the pi orbital you'd want
        #   to target would have already moved on by the opposite end of the cyclic chain of electron flow

        # Cycloaddition [4+2]
        # 6 pi
        self.addPatternListPair(
            [
                (
                    100,
                    "[*:1]([*:10])(/[*:11])=[*:2]/[*:3]=[*:4]([*:12])(/[*:13])",
                    "1 pi 2 2 2,10,11 None; 3 pi 4 2 None 3,12,13",
                ),
                (0, "[*:1]=[*:2][*:3]=[*:4]", "1 pi 2 2; 3 pi 4 2"),
            ],
            # Weird, have to specify a target atom in the other (source) pattern,
            #   including potential facial selectivity data.  Orbital must be robust enough
            #   to ignore these if they do not apply.
            [
                (
                    100,
                    "[*:5](/[*:20])([*:21])=[*:6](\[*:22])([*:23])",
                    "6 pi 5 2 5,22,23 6,20,21; 4 p None 0",
                ),
                (0, "[*:5]=[*:6]", "6 pi 5 2; 4 p None 0"),
            ],
            (1, 2, 3, 4),
            (5, 6),
        )

    def proposeOrbitalPairs(self, reactionMixture, analyzedMoleculeSmiSet=None):
        # If caller has not custom-overriden the matcher list, setup defaults
        if self.patternMatcherListPairs is None:
            self.setupDefaultPatternMatchers()

        # Setup iterators to go through every pattern match through every reaction mixture component
        for (
            sourcePatternMatcherList,
            targetPatternMatcherList,
            sourceKeyMapIndexes,
            targetKeyMapIndexes,
        ) in self.patternMatcherListPairs:
            sourceMatchesIter = ReactionMixturePatternMatchFactory(
                reactionMixture, sourcePatternMatcherList, sourceKeyMapIndexes
            )
            for sourceMatch, sourceOrbInfoStr in sourceMatchesIter:
                sourceMol = None
                for sourceAtom in sourceMatch.GetTargetAtoms():
                    sourceMol = sourceAtom.GetParent()
                    break
                sourceIsNew = self.isNewReactant(sourceMol, analyzedMoleculeSmiSet)

                targetMatchesIter = ReactionMixturePatternMatchFactory(
                    reactionMixture, targetPatternMatcherList, targetKeyMapIndexes
                )
                for targetMatch, targetOrbInfoStr in targetMatchesIter:
                    targetMol = None
                    for targetAtom in targetMatch.GetTargetAtoms():
                        targetMol = targetAtom.GetParent()
                        break
                    targetIsNew = self.isNewReactant(targetMol, analyzedMoleculeSmiSet)

                    if not sourceIsNew and not targetIsNew:
                        # Neither the source nor target molecuel proposed comes from a "new" reactant.
                        #   Don't bother analyzing, as this must be redundant to a previous analysis.
                        continue

                    if sourceMol == targetMol:
                        # Patterns are from the same molecule object.
                        # Should be using the "is" comparator, but that is returning False
                        #   even when working with the same molecule object for some reason.
                        # Interpret this in 2 ways, as indeed an intra-molecular reaction,
                        #   or prepare separate molecule copies for an inter-molecular reaction.
                        # Should already work as intra-molecular example, but prepare a copy
                        #   so caller is free to manipulate it without disrupting the object contents

                        if self.compatiblePatternMatches(
                            sourceMatch,
                            targetMatch,
                            sourceKeyMapIndexes,
                            targetKeyMapIndexes,
                        ):
                            # Beware, if on the same molecule, these might be based on overlapping
                            #   atoms, which would not make sense to "react" them together.
                            #   Only allow the intra-molecular case if the patterns are compatible (do not overlap)

                            # For every atom that was matched, apply an atom map label so we can build orbitals off of it
                            for atomMatchPair in sourceMatch.GetAtoms():
                                atomMatchPair.target.SetMapIdx(
                                    atomMatchPair.pattern.GetMapIdx()
                                )

                            for atomMatchPair in targetMatch.GetAtoms():
                                atomMatchPair.target.SetMapIdx(
                                    atomMatchPair.pattern.GetMapIdx()
                                )

                            # Generate the orbital objects based on the labeled molecule and info strings.
                            sourceOrb = Orbital.fromLabeledMolAndInfoStr(
                                sourceMol, sourceOrbInfoStr
                            )
                            targetOrb = Orbital.fromLabeledMolAndInfoStr(
                                targetMol, targetOrbInfoStr
                            )

                            # Wipe out atom map indexes to revert the molecule back to its original state
                            for atomMatchPair in sourceMatch.GetAtoms():
                                atomMatchPair.target.SetMapIdx(0)

                            for atomMatchPair in targetMatch.GetAtoms():
                                atomMatchPair.target.SetMapIdx(0)

                            # Create clones of the orbitals and molecule so the caller is free to modify them
                            (
                                (cloneSourceOrb, cloneTargetOrb),
                                copyMol,
                            ) = Orbital.cloneOrbitalsWithCommonMol(
                                [sourceOrb, targetOrb]
                            )

                            yield (cloneSourceOrb, cloneTargetOrb)

                        # Prepare copies to yield inter-molecular example
                        # Don't need to do it here, just reuse the usual 2 molecule case,
                        #   which will work out equivalently
                        pass

                    # Patterns from distinct molecule objects.  Prepare a composite
                    #   molecule object that both orbitals will be based upon
                    compositeMol = OEGraphMol()

                    # For every atom that was matched, apply an atom map label so we can build orbitals off of it
                    for atomMatchPair in sourceMatch.GetAtoms():
                        atomMatchPair.target.SetMapIdx(
                            atomMatchPair.pattern.GetMapIdx()
                        )
                    OEAddMols(compositeMol, sourceMol)
                    # Wipe out atom map indexes to revert the original molecule back to its original state
                    for atomMatchPair in sourceMatch.GetAtoms():
                        atomMatchPair.target.SetMapIdx(0)

                    # For every atom that was matched, apply an atom map label so we can build orbitals off of it
                    for atomMatchPair in targetMatch.GetAtoms():
                        atomMatchPair.target.SetMapIdx(
                            atomMatchPair.pattern.GetMapIdx()
                        )
                    OEAddMols(compositeMol, targetMol)
                    # Wipe out atom map indexes to revert the original molecule back to its original state
                    for atomMatchPair in targetMatch.GetAtoms():
                        atomMatchPair.target.SetMapIdx(0)

                    # Generate the orbital objects based on the labeled molecule and info strings.
                    sourceOrb = Orbital.fromLabeledMolAndInfoStr(
                        compositeMol, sourceOrbInfoStr
                    )
                    targetOrb = Orbital.fromLabeledMolAndInfoStr(
                        compositeMol, targetOrbInfoStr
                    )

                    # Wipe out atom map indexes on the composite mol to a clean state
                    for atom in compositeMol.GetAtoms():
                        atom.SetMapIdx(0)

                    yield (sourceOrb, targetOrb)

    def compatiblePatternMatches(
        self, sourceMatch, targetMatch, sourceKeyMapIndexes, targetKeyMapIndexes
    ):
        """Determine if the source and target matches are compatible.
        Assuming they are on different molecule objects, it's definitely fine,
        but if they are from the same molecule, they must not overlap.

        Could check just the key map indexes, but actually, if they overlap at all,
        this can mess up subsequent orbital perception from information strings, so disallow any overlap right now
        """
        sourceMol = None
        targetMol = None

        sourceIndexes = set()
        targetIndexes = set()

        for atomMatchPair in sourceMatch.GetAtoms():
            if sourceMol is None:
                sourceMol = atomMatchPair.target.GetParent()
            sourceIndexes.add(atomMatchPair.target.GetIdx())

        for atomMatchPair in targetMatch.GetAtoms():
            if targetMol is None:
                targetMol = atomMatchPair.target.GetParent()
            targetIndexes.add(atomMatchPair.target.GetIdx())

        isCompatible = True
        if sourceMol == targetMol:
            # Should be using "is" operator, but that doesn't seem to work for OEGraphMol objects
            # If from the same molecule object, only compatible if there is no overlap in the atoms covered
            isCompatible = len(sourceIndexes.intersection(targetIndexes)) < 1
        return isCompatible


class ReactionMixturePatternMatchFactory:
    """Factory object that yields pattern matches (OEMatchBase)
    for all of the molecules in a reaction mixture
    based on a prioritized list of OESubSearch matcher objects.
    """

    reactionMixture = None
    patternMatcherList = None
    keyMapIndexes = None

    def __init__(self, reactionMixture, patternMatcherList, keyMapIndexes):
        """Initialize with the reaction mixture to iterate through,
        and pre-initialized lists of OESubSearch matchers object to find the patterns.
        """
        self.reactionMixture = reactionMixture
        self.patternMatcherList = patternMatcherList
        self.keyMapIndexes = keyMapIndexes

    def __iter__(self):
        for smiles, mixtureComponent in self.reactionMixture.iteritems():
            # Need a mol object for back / parent pointers on orbitals,
            #   but check if one has already been prepared.
            # Want to use the same ones so that intramolecular reactions can be recognized.
            if mixtureComponent.mol is None:
                mixtureComponent.mol = OEGraphMol()
                OEParseSmiles(mixtureComponent.mol, smiles)
            mol = mixtureComponent.mol

            if mixtureComponent.currentQuantity() > 0:
                # SMARTS patterns looking for hydrogens (evan as wildcard "*" atoms)
                #   will only match explicit hydrogen atoms, not implicit ones
                OEAddExplicitHydrogens(mol)

                # Look for pattern matching on any of the mixture components,
                #   but preferentially use the highest priority one if applying to redundant atoms
                maxPriorityByRedundancyKey = dict()

                for (
                    priority,
                    smartsPattern,
                    matcher,
                    orbInfoStr,
                ) in self.patternMatcherList:
                    # Every match for a given mixture component
                    for match in matcher.Match(mol, False):
                        # Only accept matches if this pattern has a priority ranking
                        #   as high or higher than any redundant matches found
                        redundancyKey = self.patternRedundancyKey(match)
                        if (
                            redundancyKey not in maxPriorityByRedundancyKey
                            or priority >= maxPriorityByRedundancyKey[redundancyKey]
                        ):
                            yield (match, orbInfoStr)
                            maxPriorityByRedundancyKey[redundancyKey] = priority

    def patternRedundancyKey(self, match):
        """Come up with a unique key that represents the relevant information in
        the pattern match that should be used to differentiate
        redundant proposals that likely just represent more or less specific
        versions of the same reaction.

        Usually the difference is stereochemistry, so we don't want to let orbital
        facial selectivity affect this redundancy key, and we can't count on consistent
        atom mapping indexes, so mainly just the use underlying atom index numbers
        on the molecule object that should remain constant across pattern matches.
        """
        atomMatchPairByMapIndex = dict()
        for atomMatchPair in match.GetAtoms():
            patternMapIndex = atomMatchPair.pattern.GetMapIdx()
            atomMatchPairByMapIndex[patternMapIndex] = atomMatchPair

        redundancyKey = list()

        for keyMapIndex in self.keyMapIndexes:
            atomMatchPair = atomMatchPairByMapIndex[keyMapIndex]
            redundancyKey.append(atomMatchPair.target.GetIdx())

        return tuple(redundancyKey)


class ReactionMixtureOrbitalFactory:
    """Factory object that "yields" orbitals for all of the molecules
    in the reaction mixture.  Has the same effect as returning an actual list
    of all of the relevant orbitals from the reaction mixture,
    but using this generator factory approach can save temporary memory usage
    since only one orbital ever needs to be referenced at a time.
    """

    reactionMixture = None
    reactiveOrbitalFactory = None

    def __init__(self, reactionMixture, orbitalReactivityCalc, restrictSymmetry=True):
        self.reactionMixture = reactionMixture
        self.reactiveOrbitalFactory = ReactiveOrbitalFactory(
            orbitalReactivityCalc, restrictSymmetry
        )

    def __iter__(self):
        """Return all of the filled (occupied) or unfilled (unoccupied)
        orbitals for the molecules in the reactionMixture, but only those
        which pass the given orbital reactivity probability estimator.
        Skip orbitals from components that have 0 quantity in the reaction mixture.
        """
        for orb in self.scoreIter():
            if orb.dataTags["passThreshold"]:
                yield orb

    def scoreIter(self):
        """More detailed version of the iterator.  Yields all
        orbitals, not just the ones that "passThreshold" for reactivity.
        Scoring information can be found as elements in each orbitals'
        dataTags dictionary: (reactivityScore, relevanceProb, passThreshold).
        """
        for smiles, mixtureComponent in self.reactionMixture.iteritems():
            # Need a mol object for back / parent pointers on orbitals,
            #   but check if one has already been prepared.
            # Want to use the same ones so that intramolecular reactions can be recognized.
            if mixtureComponent.mol is None:
                mixtureComponent.mol = OEGraphMol()
                OEParseSmiles(mixtureComponent.mol, smiles)
            mol = mixtureComponent.mol

            if mixtureComponent.currentQuantity() > 0:
                for orb in self.reactiveOrbitalFactory.scoreIter(mol):
                    yield orb


class BaseOrbitalPairFilter:
    """Simple base class for filters that analyze a proposed orbital
    pair and return True/False whether the proposal should pass the filter or not.
    """

    def __call__(self, filledOrb, unfilledOrb):
        raise NotImplementedError()


class OverlapAtomsFilter(BaseOrbitalPairFilter):
    """Disallow any orbital pairs where the atoms in one orbital are in another orbital"""

    def __call__(self, filledOrb, unfilledOrb):
        return compatibleOrbitalPair(filledOrb, unfilledOrb)


class ExtremeChargeFilter(BaseOrbitalPairFilter):
    """Disallow any orbital pairs that would result in a
    product with extreme formal charges on any single atom.
    """

    def __call__(self, filledOrb, unfilledOrb):
        moveOrbitalElectrons(filledOrb, unfilledOrb)

        violationFound = False
        for atom in filledOrb.mol.GetAtoms():
            if not self.validAtomCharge(atom):
                # Already found a violation, don't need to look further
                violationFound = True
                break

        undoMoveOrbitalElectrons(filledOrb, unfilledOrb)

        return not violationFound

    def validAtomCharge(atom):
        valid = True
        valid = (
            valid
            and atom.GetFormalCharge()
            >= ExtremeChargeFilter.minAllowedCharge(atom.GetAtomicNum())
        )
        valid = (
            valid
            and atom.GetFormalCharge()
            <= ExtremeChargeFilter.maxAllowedCharge(atom.GetAtomicNum())
        )
        return valid

    validAtomCharge = staticmethod(validAtomCharge)

    def minAllowedCharge(atomicNum):
        return -1

    minAllowedCharge = staticmethod(minAllowedCharge)

    def maxAllowedCharge(atomicNum):
        # Weird, but third row elements like S and P could be written in +2 form like sulfuric acid
        if atomicNum in (
            15,
            16,
            12,
        ):
            return +2
        return +1

    maxAllowedCharge = staticmethod(maxAllowedCharge)


class ExtremeBondFilter(BaseOrbitalPairFilter):
    """Disallow orbital pairs that create bonds with order >=4"""

    MAX_BOND_ORDER = 3

    def __call__(self, filledOrb, unfilledOrb):
        moveOrbitalElectrons(filledOrb, unfilledOrb)

        violationFound = False
        for bond in filledOrb.mol.GetBonds():
            if not self.validBondOrder(bond):
                # Already found a violation, don't need to look further
                violationFound = True
                break

        undoMoveOrbitalElectrons(filledOrb, unfilledOrb)

        return not violationFound

    def validBondOrder(self, bond):
        """Test that a bond is within an acceptable range"""
        return bond.GetOrder() <= self.MAX_BOND_ORDER


class SubstitutionHybridizationFilter(BaseOrbitalPairFilter):
    """Disallow any substitution reactions (sigma* targets)
    at sp2 or sp hybridized centers.  Should only be
    reasonable at sp3 hybridized centers.

    This IS reasonable in a halogenation of an alkyne.  One step here is the
    substitution to break the halonium ion cycle at an sp2 center.

    Allow a special case where there is a bridged halonium.

    """

    def __call__(self, filledOrb, unfilledOrb):
        if unfilledOrb.type == "sigma" and atomHybridization(unfilledOrb.atom) in (
            1,
            2,
        ):
            # Looks bad with substitution (sigma*) target hybridization,
            #   but could be okay if this is a local orbital interaction between directly neighboring atoms
            #   or could be ok if the substitution is to break up a bridged halonium ion.
            isLocal = filledOrb.atom.GetBond(unfilledOrb.atom) is not None
            is3Ring = any(
                [
                    (unfilledOrb.atom.GetBond(atom) is not None)
                    for atom in unfilledOrb.neighbor.GetAtoms()
                    if atom is not unfilledOrb.atom
                ]
            )
            # log.info('fOrb: %s, uOrb: %s, isLocal:%s, is3Ring:%s' % \
            #    (pformat(filledOrb.toLabeledSmiAndInfoStr(10)), pformat(unfilledOrb.toLabeledSmiAndInfoStr(20)),
            #    str(isLocal), str(is3Ring)))
            if not (isLocal or is3Ring):
                # Atoms are not direct neighbors, or in a 3-ring
                # looks like a bad substitution then
                return False
        return True


class ReactionProposerParameters:
    """Simple struct to encapsulate parameters that could be set at the UI level."""

    useLocalOrbitalPairs = False
    acceptAllOrbitals = False

    def __init__(self, useLocalOrbitalPairs=False, acceptAllOrbitals=False):
        """Basic Constructor"""
        self.useLocalOrbitalPairs = useLocalOrbitalPairs
        self.acceptAllOrbitals = acceptAllOrbitals
