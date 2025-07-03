#!/usr/bin/env python
# encoding: utf-8
"""
Classes for proposing polar reactions

File: PolarProposer.py
Created 2012-01-06
Copyright 2012. Institute of Genomics and Bioinformatics

Subclassing the AtomBasedReactionProposer
"""


from chemutils.score.orbital.ReactivityScore import ReactiveOrbitalFactory
from chemutils.Common.OrbitalModel import orbitalIter

from chemutils.CombiCDB.proposers.AtomBasedReactionProposer import (
    AtomBasedReactionProposer,
    AtomBasedOrbitalFactory,
)

from chemutils.CombiCDB.proposers import filters

defaultFilters = [
    filters.NonReasonableOrbitalsFilter(),
    filters.OverlapAtomsFilter(),
    filters.ExtremeChargeFilter(),
    filters.SubstitutionHybridizationFilter(),
    filters.ExtremeBondFilter(),
    filters.PericyclicFilter(),
    filters.PiBondBondDissociationFilter(),
    filters.AromaticSingleBondSourceFilter(),
    filters.SourceSinkDistanceFilter(),
]


class PolarProposer(AtomBasedReactionProposer):
    """Propose polar 2-electron source to sink reactions, based on possible atom sources

    Sub-class of AtomBasedReactionProposer - set the attributes.
    """

    def __init__(
        self,
        srcAtomList,
        sinkAtomList,
        orbitalPairFilters=defaultFilters,
        checkIntra=True,
    ):
        """Set specifics for Polar reactions"""
        super(PolarProposer, self).__init__(
            srcAtomList, sinkAtomList, orbitalPairFilters, checkIntra=checkIntra
        )

    def setOrbitalFactories(self):
        """Specify what the factory objects will be.  Called by super's __init__"""
        self.srcOrbFactory = PolarOrbitalFactory(
            self.srcAtomList, occupied=True, allSymmetry=True
        )
        self.sinkOrbFactory = PolarOrbitalFactory(
            self.sinkAtomList, occupied=False, allSymmetry=False
        )


class PolarOrbitalFactory(AtomBasedOrbitalFactory):
    """Class to handle atom based 2-electron, polar orbital proposals"""

    def __init__(self, atomObjList, occupied=False, allSymmetry=True, allKekule=True):
        """Constructor setup attributes"""
        super(PolarOrbitalFactory, self).__init__(atomObjList, allSymmetry, allKekule)
        self.occupied = occupied

    @staticmethod
    def isValidOrbitalOccupation(orbital, occupied):
        """Verify that the orbital is actually filled / unfilled
        No probability that an empty orbital is a good filled orbital or
        that a lone pair is a good unfilled orbital.
        """
        return not (
            (occupied and orbital.isEmpty()) or (not occupied and orbital.isLonePair())
        )

    @staticmethod
    def prepareChainOrbital(orbital, parentOrb, occupied):
        """Prepare a chain orbital.

        This evaluates an orbital based on its suitability for chaining.1

        For a filled orbital, this should be an adjacent pi bond or a lone pair.
          or C-Metal sigma bond
        For an unfilled orbital, this should be an adjacent pi bond,
          or if the atom is a hydrogen, an adjacent empty orbital, or sigma bond
          or parentOrb is a pi bond, orb is sigma, and neighbor is not a carbon or hydrogen

        If the orbital is a valid candidate, return it, otherwise, return None
        """
        chainOrbital = None

        hasExtendedPotential = orbital.type == "pi"
        hasExtendedPotential = hasExtendedPotential or (
            occupied and orbital.isLonePair()
        )
        hasExtendedPotential = hasExtendedPotential or (
            occupied
            and orbital.type == "sigma"
            and orbital.atom.GetAtomicNum() == 6
            and orbital.neighbor.GetAtomicNum() in (3, 12)
        )
        hasExtendedPotential = hasExtendedPotential or (
            not occupied and orbital.isEmpty()
        )
        hasExtendedPotential = hasExtendedPotential or (
            not occupied
            and parentOrb.atom.GetAtomicNum() == 1
            and orbital.neighbor is not None
        )
        hasExtendedPotential = hasExtendedPotential or (
            not occupied
            and parentOrb.type == "pi"
            and orbital.type == "sigma"
            and orbital.neighbor.GetAtomicNum() not in [1, 6]
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
        """Method to handle the details of proposing the chained orbitals.

        Doing a couple of things here:
        - First, using symmetries to make reasonable culls of equivalent proposals
        - Second, proposing doubly-chained proposals if particularly favorable.
            e.g., if the ultimate sink is an empty orbital.
        """
        if not orb.isBondOrbital():
            return
        seenSymSet = set([])

        for atom in orb.neighbor.GetAtoms():
            if atom == orb.atom or atom.GetSymmetryClass() in seenSymSet:
                continue
            seenSymSet.add(atom.GetSymmetryClass())
            seenNeighborSymClass = set([])
            for potChainOrb in orbitalIter(atom, useSymmetries=False):
                ## Check if we are looking back to ourselves or at symmetrically equiv atom
                if potChainOrb.neighbor is not None:
                    if (
                        potChainOrb.neighbor == orb.neighbor
                        or potChainOrb.neighbor.GetSymmetryClass()
                        in seenNeighborSymClass
                    ):
                        continue
                    seenNeighborSymClass.add(potChainOrb.neighbor.GetSymmetryClass())

                chainOrb = self.prepareChainOrbital(potChainOrb, orb, self.occupied)
                if chainOrb is not None:
                    yield chainOrb
                    for doubleChainedOrb in self.findPotentialDoubleChainedOrbital(
                        chainOrb
                    ):
                        yield doubleChainedOrb

    def findPotentialDoubleChainedOrbital(self, chainOrb):
        """Given a chained orb.  Is there a nice orbital close by to extend the chain."""
        if chainOrb.extOrbital.type != "pi":
            return

        ## Look for a potential extra cation
        for potEmptyOrbAtom in chainOrb.extOrbital.neighbor.GetAtoms():
            potential = potEmptyOrbAtom is not chainOrb.extOrbital.atom
            orbList = []
            if self.occupied:
                ## Looking for a lone pair
                orbList = [
                    aOrb for aOrb in orbitalIter(potEmptyOrbAtom) if aOrb.isLonePair()
                ]
            else:
                ## Looking for an empty orb (or a non-carbon neighbor double bond.)
                orbList = [
                    aOrb
                    for aOrb in orbitalIter(potEmptyOrbAtom)
                    if aOrb.isEmpty()
                    or (aOrb.type == "pi" and aOrb.neighbor.GetAtomicNum() != 6)
                ]
            for potOrb in orbList:
                chainOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(chainOrb)
                candidateOrbIdxSet = ReactiveOrbitalFactory.orbitalAtomIdxSet(potOrb)
                if len(chainOrbIdxSet.intersection(candidateOrbIdxSet)) == 0:
                    ## No overlap, looks good.
                    nChainOrbital = chainOrb.simpleCopy()
                    nChainOrbital.extOrbital.extOrbital = potOrb
                    yield nChainOrbital

    def __iter__(self):
        """Generator for the orbitals"""
        for oeAtom, atomObj in self.atomIterator():
            for orb in orbitalIter(oeAtom, useSymmetries=(not self.allSymmetry)):
                if self.isValidOrbitalOccupation(orb, self.occupied):
                    yield orb, atomObj
                    # Look for neighbor chained orbitals
                    for chainOrb in self.proposeChainedOrbitals(orb):
                        yield chainOrb, atomObj
