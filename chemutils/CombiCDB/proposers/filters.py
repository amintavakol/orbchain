from chemutils.CombiCDB.proposers.Util import log
from chemutils.CombiCDB.proposers.Util import isBondDissociation, isReasonableOrbital, compatibleOrbitalPair
from chemutils.CombiCDB.proposers.Util import labelOrbPair
from chemutils.CombiCDB.OrbitalInteraction import isPericyclic
from chemutils.CombiCDB.ReactionProposer import ExtremeChargeFilter, BaseOrbitalPairFilter
from chemutils.CombiCDB.ReactionProposer import SubstitutionHybridizationFilter, ExtremeBondFilter

from openeye.oechem import OEDetermineComponents, OEAssignAromaticFlags, OEClearAromaticFlags, OEAroModelMMFF, OEMolToSmiles
from openeye.oechem import *

class PossiblePericyclicFilter(BaseOrbitalPairFilter):
    """Simple predicate to test if this orb pair is possibly a """
    def __call__(self, srcOrb, sinkOrb):
        """Head and tail atoms have to be the same or next to each other.

        Bond Dissociations aren't pericyclic, and there must be at least one extOrbital
        """
        if isBondDissociation(srcOrb, sinkOrb):
            return False
        if not srcOrb.extOrbital and not sinkOrb.extOrbital:
            return False
        
        srcHeadAtom = self.findEndAtomOrb(srcOrb)
        sinkHeadAtom = self.findEndAtomOrb(sinkOrb)

        if srcHeadAtom == sinkHeadAtom:
            return True
        if srcHeadAtom.GetBond(sinkHeadAtom) is not None:
            return True
        return False

    def findEndAtomOrb(self, orb):
        """Recursive function to find the end atom of an orbital"""
        if orb.extOrbital is not None:
            return self.findEndAtomOrb(orb.extOrbital)
        elif orb.neighbor is not None:
            return orb.neighbor
        return orb.atom
        
class AromaticSingleBondSourceFilter(BaseOrbitalPairFilter):
    """Make sure we don't use the single bond of a Kekulized aromatic ring to attack and open the ring."""
    def __call__(self, srcOrb, sinkOrb):
        
        for o in self.chained_orb_iter(srcOrb):
            #log.debug('orb: %s' % o)
            if o.type == 'sigma':
                if o.atom.IsInRing() and o.neighbor is not None and o.neighbor.IsInRing():
                    # check if atom is part of aromatic ring
                    mol = o.atom.GetParent()
                    o.atom.SetMapIdx(1)
                    #log.info('smi: %s' % OEMolToSmiles(mol))
                    mol_smi = OEMolToSmiles(mol)
                    new_mol = OEGraphMol()
                    OESmilesToMol(new_mol, mol_smi)
                    OEAssignAromaticFlags(new_mol, OEAroModelMMFF)
                    for a in new_mol.GetAtoms():
                        if a.GetMapIdx() == 1 and a.IsAromatic():
                            return False

        return True

    def chained_orb_iter(self, orb):
        """Recursive function to yield all orbs in extended orb"""
        if orb.extOrbital is not None:
            self.chained_orb_iter(orb.extOrbital)
        yield orb

class PericyclicFilter(BaseOrbitalPairFilter):
    """Disallow pericyclic reactions"""
    def __call__(self, srcOrb, sinkOrb):
        """Simple filter that returns true if a reaction is pericyclic"""
        return not isPericyclic(srcOrb, sinkOrb);

class OverlapAtomsFilter(BaseOrbitalPairFilter):
    """Disallow any orbital pairs where the atoms in one orbital are in another orbital"""
    def __call__(self, filledOrb, unfilledOrb):
        return compatibleOrbitalPair(filledOrb, unfilledOrb);

class PiBondBondDissociationFilter(BaseOrbitalPairFilter):
    """Disallow a common source of resonance errors, a pi bond dissociating"""
    def __call__(self, filledOrb, unfilledOrb):
        return not (isBondDissociation(filledOrb, unfilledOrb) and unfilledOrb.type == 'pi')


class NonReasonableOrbitalsFilter(BaseOrbitalPairFilter):
    """Disallow any interactions where the orbs no longer make sense"""
    def __call__(self, filledOrb, unfilledOrb):
        reasonableFilled = isReasonableOrbital(filledOrb)
        reasonableUnfilled = isReasonableOrbital(unfilledOrb)
        
        if not reasonableUnfilled: 
            #log.info('Going to filter as unreasonable unfilledOrb: %s' %
            # pformat(unfilledOrb.toLabeledSmiAndInfoStr(10)))
            pass
        if not reasonableFilled:
            #log.info('Going to filter as unreasonable filledOrb: %s' %
            # pformat(filledOrb.toLabeledSmiAndInfoStr(20)))
            pass
        
        return (reasonableFilled and reasonableUnfilled)

class SourceSinkDistanceFilter(BaseOrbitalPairFilter):
    """Per Aaron request (see plot), do not allow intramolecular attack beyond N atoms
	NOTE this was originally implemented in ~/rxnpred/tools/Filters.py. 
	Trying here to see if it improves performance."""

    def __init__(self):
        self.MAX_SOURCE_SINK_DISTANCE = 8

    def __call__(self, filled_orb, unfilled_orb):
        distance = OEGetPathLength(filled_orb.atom, unfilled_orb.atom)
        if distance > self.MAX_SOURCE_SINK_DISTANCE:
            #log.error('filled_orb.atom: %s' % filled_orb.atom)
            #log.error('unfilled_orb.atom: %s' % unfilled_orb.atom)
            #log.error('SrcSinkDistance was: %d' % distance)
            return False

        return True
