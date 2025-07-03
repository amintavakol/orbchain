"""Miscellaneous utility functions used across the application

"""
from chemutils.CombiCDB.proposers import Const
import sys
import logging

from chemutils.Common.CanonicalAtomMapSmiles import (
    canonicalizeAtomMapSmiString,
)
from chemutils.CombiCDB.OrbitalInteraction import isPericyclic
from chemutils.Common.MolExt import clearAtomMaps

## Set up logging
log = logging.getLogger(Const.APPLICATION_NAME)
log.setLevel(Const.LOGGER_LEVEL)

handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter(Const.LOGGER_FORMAT)

handler.setFormatter(formatter)
log.addHandler(handler)


## Utils - Input for all of these are assumed to be CHEM.Common.OrbitalModel.Orbital objs
def isBondDissociation(sourceOrb, targetOrb):
    """Determine whether the given filled and unfilled orbital combination
    represents a bond dissociation reaction.

    Note:  Bond dissociation should have exactly the same orbitals
    """
    bondDissociation = (
        sourceOrb is None and targetOrb and targetOrb.neighbor is not None
    )
    bondDissociation = bondDissociation or (
        sourceOrb.neighbor is not None
        and targetOrb.neighbor is not None
        and (sourceOrb.atom.GetIdx(), sourceOrb.neighbor.GetIdx())
        == (targetOrb.atom.GetIdx(), targetOrb.neighbor.GetIdx())
    )
    return bondDissociation


def isReasonableOrbital(orb):
    """Simple predicate test to check that an orbital truly is reasonable.  This is basically a check
    to be sure that pi bond orbitals reference bonds with order > 1, and that sigma bond orbitals
    reference bonds with order == 1."""
    if orb.neighbor is not None:
        bondOrder = orb.atom.GetBond(orb.neighbor).GetOrder()
        if orb.type == "pi" and bondOrder < 2:
            return False
        if orb.type == "sigma" and bondOrder > 1:
            return False
        if orb.extOrbital is not None:
            return isReasonableOrbital(orb.extOrbital)
    return True


def compatibleOrbitalPair(filledOrb, unfilledOrb):
    """Assuming the orbitals are based on the same molecule object,
    ensure that it makes sense to combine / react these together.
    In particular, ensure that they do not overlap in terms of the
    atoms they are based on, which would make any attempt to
    "react" them lead to an inconsistent molecular state.
    """
    compatible = True

    # None of the orbitals in the filled chain can be empty
    # Iterate through chain has slight risk of infinite loop (cyclic linked list)
    currentOrb = filledOrb
    while currentOrb is not None:
        compatible = compatible and not currentOrb.isEmpty()
        currentOrb = currentOrb.extOrbital

    # None of the orbitals in the unfilled chain can be full (lone pair)
    currentOrb = unfilledOrb
    while currentOrb is not None:
        compatible = compatible and not currentOrb.isLonePair()
        currentOrb = currentOrb.extOrbital

    if compatible:
        # Incompatible if any atoms in common, unless a
        #   special case representing a bond dissociation or pericyclic reaction
        # Specialized reaction types where overlapping orbital proposals can make sense
        if isBondDissociation(filledOrb, unfilledOrb):
            # log.debug('IS BOND DISSOCIATION (compatibleOrbitalPair)')
            # Bond Dissocations should NOT be chained.
            compatible = filledOrb is None or (
                filledOrb is not None and filledOrb.extOrbital is None
            )
            compatible = compatible and unfilledOrb.extOrbital is None
        elif isPericyclic(filledOrb, unfilledOrb):
            compatible = True
        else:
            # General case, disallow any proposals where the orbitals have overlapping atoms
            filledAtomIndexes = filledOrb.coveredOrbitalAtomIdx()
            unfilledAtomIndexes = unfilledOrb.coveredOrbitalAtomIdx()

            indexIntersection = filledAtomIndexes.intersection(unfilledAtomIndexes)

            compatible = len(indexIntersection) < 1

    return compatible


def labelOrbPair(src, sink, arom=True):
    """Convenience to turn an orb pair into a string tuple"""
    clearAtomMaps(src.mol)
    if not isBondDissociation(src, sink):
        src.labelOrbitalAtoms(10)
    sink.labelOrbitalAtoms(20)

    nSrcInfo = src.toLabeledSmiAndInfoStr()
    nSinkInfo = sink.toLabeledSmiAndInfoStr()

    return (
        canonicalizeAtomMapSmiString(nSrcInfo[0], arom=arom),
        nSrcInfo[1],
        nSinkInfo[1],
    )
