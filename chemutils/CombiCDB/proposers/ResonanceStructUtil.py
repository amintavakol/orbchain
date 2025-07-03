from openeye.oechem import *
from chemutils.Common.Util import molBySmiles

from chemutils.Common.MolExt import molSkeleton
from chemutils.Common.Const import REACTION_COMPONENT_DELIM, SMILES_MOL_DELIM
from chemutils.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair
from chemutils.Common.OrbitalModel import orbitalInfo
from chemutils.CombiCDB.Util import log


class ResonanceStructUtil(object):
    """Class that handles telling whether two mols are resonance equivalent.

    The key idea is that two equivalent resonance structures will have the same:
    - Molecular formula
    - Molecular skeleton
    - Net Formal Charge
    - Number of radicals
    """

    def __init__(self, reactantMol=None, checkRadical=False):
        """Constructor, if reactantMol is set, then set up some data structures to determine
        future reaction equivalence.
        """
        self.reactantMol = reactantMol
        self.reactantIdentTuple = None
        self.checkRadical = checkRadical
        if self.reactantMol is not None:
            self.reactantIdentTuple = self.identifyingTupleFromMol(self.reactantMol)

        self.seenIdentSet = set()

    def identifyingTupleFromMol(self, mol):
        """Turn a mol into the identifyingTuple  (if self.checkRadical) include number of radicals"""
        smiStr = OECreateIsoSmiString(mol)
        chunks = smiStr.split(SMILES_MOL_DELIM)

        forms = [OEMolecularFormula(molBySmiles(ch)) for ch in chunks]
        # mol = OEGraphMol()
        # forms = []
        # for c in chunks:
        #    mol.Clear()
        #    OEParseSmiles(mol, c)
        #    forms.append(OEMolecularFormula(mol))
        forms.sort()

        if self.checkRadical:
            nradicals = sum([orbitalInfo(a)["nRadicals"] for a in mol.GetAtoms()])
            return (
                tuple(forms),
                OECreateIsoSmiString(molSkeleton(mol, retainHs=True)),
                OENetCharge(mol),
                nradicals,
            )
        return (
            tuple(forms),
            OECreateIsoSmiString(molSkeleton(mol, retainHs=True)),
            OENetCharge(mol),
        )

    def identifyingTupleFromSmi(self, smi):
        """Turn smiles into identifyingTuple"""
        #        mol = molBySmiles(smi)
        mol = OEGraphMol()
        try:
            OEParseSmiles(mol, smi)
        except:
            log.debug("Exception! mol: %s, smi: %s" % (mol, smi))
            log.debug("mol: %s" % mol)
            log.debug("smi: %s" % smi)

        return self.identifyingTupleFromMol(mol)

    def identifyingTupleFromOrbPair(self, src, sink):
        """Turn a src/sink orb pair's product into an identifyingTuple"""
        rsmiles = reactionSmilesFromOrbitalPair(src, sink)
        chunks = rsmiles.split(REACTION_COMPONENT_DELIM)
        prodSmi = chunks[-1]
        return self.identifyingTupleFromSmi(prodSmi)

    def checkOrbPair(self, src, sink):
        """Returns False if the src, sink combo should be filtered. Given a src and sink, check whether we've seen before or equiv to reactants.
        If not, add to seen set."""
        # identTuple = self.identifyingTupleFromOrbPair(src, sink)
        # return identTuple != self.reactantIdentTuple

        # Want to allow resonance structs of *reactants* to prevent dead-ends
        identTuple = self.identifyingTupleFromOrbPair(src, sink)
        if identTuple not in self.seenIdentSet:
            # for reactant resonance structs
            # if identTuple == self.reactantIdentTuple:
            #    return True

            self.seenIdentSet.add(identTuple)
            return True
        else:
            return False
