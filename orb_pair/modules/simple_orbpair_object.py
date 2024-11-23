from chemutils.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair, ElementaryStepData
from chemutils.CombiCDB.proposers.Util import isBondDissociation
from chemutils.CombiCDB.MechanismModel import ElectronArrow
from chemutils.CombiCDB.Const import REACTION_DELIM
from chemutils.Common.MolExt import clearAtomMaps
from chemutils.Common.Util import exact_mass

from chemutils.CombiCDB.proposers.PolarProposer import PolarProposer
from chemutils.CombiCDB.proposers.RadicalProposer import RadicalProposer
from chemutils.CombiCDB.proposers.PericyclicProposer import PericyclicProposer

from openeye.oechem import *  # for mass calculation added 02.16.16
from orb_pair.modules.Filters import Filters
import random


class SimpleOrbPairObject:

    # Filtering
    pass_filters = Filters()

    """Class to represent a basic Orb Pair"""
    def __init__(self, srcOrb, sinkOrb, srcAtom, sinkAtom, radical=False, pericyclic=False):
        """Constructor"""
        self.srcOrb = srcOrb
        self.sinkOrb = sinkOrb
        self.srcAtom = srcAtom
        self.sinkAtom = sinkAtom

        self.radical = radical
        self.pericyclic = pericyclic

        self.labelOrbPair(self.srcOrb, self.sinkOrb)

        self.connectedNonMappedSmiles = self.srcAtom.connectedNonMappedSmiles

        self.arrowObjectList = []
        self.reactionSmiles = reactionSmilesFromOrbitalPair(self.srcOrb, self.sinkOrb, False, self.arrowObjectList)
        self.reactantSmiles, self.productSmiles = self.reactionSmiles.split(REACTION_DELIM)
        # self.reactant_masses = []
        self.product_masses = self.get_masses(self.productSmiles)
        # self.product_mass_string = self.masses_to_string(self.product_masses)
        self.arrowCodes = ElectronArrow.prepareArrowCodes(self.arrowObjectList)

        # For eventual ML
        self.rawFeatDict = None
        self.normFeatDict = None
        self.predValue = None
        self.rawPredValues = None
        self.atom_pred_vals_sum = None

        # self.filtered_orbpairs = {}

    def toElementaryStepData(self):
        """Make an elementary step object version of self"""
        eStep = ElementaryStepData()
        eStep['compositeMol'] = self.srcOrb.mol
        eStep['filledOrb'] = self.srcOrb
        eStep['unfilledOrb'] = self.sinkOrb
        return eStep
    
    def __str__(self):
        """Show something"""
        return str((self.reactionSmiles, self.arrowCodes))
    
    def __repr__(self):
        return "'%s'" % self.__str__()
    
    @classmethod
    def orbPairObjectsFromAtoms(cls, srcAtomList, sinkAtomList, photo=False, radical=False, pericyclic=False,
                                checkIntra=True, rand_one: bool = False):
        """Given some atoms, run the proposer and make a list of SimpleOrbPairObjects"""
        from utils.timing_utils import Times
        import time
        if radical:
            proposer = RadicalProposer(srcAtomList, sinkAtomList, photo=photo, checkIntra=checkIntra)
        elif pericyclic:
            proposer = PericyclicProposer(srcAtomList)
        else:
            proposer = PolarProposer(srcAtomList, sinkAtomList, checkIntra=checkIntra)
        orbPairObjectList = []

        orbParisProp = [x for x in proposer.proposeOrbitalPairs()]
        if rand_one:
            random.shuffle(orbParisProp)

        times = Times()
        for srcOrb, sinkOrb, srcAtom, sinkAtom in orbParisProp:
            start_time = time.time()
            pass_filters1: bool = proposer.passOrbitalPairFilters(srcOrb, sinkOrb)
            times.record_time("passFilters1", time.time() - start_time)
            if not pass_filters1:
                continue

            start_time = time.time()
            if pericyclic:
                new_op = cls(srcOrb, sinkOrb, srcAtomList[0], None, radical, pericyclic)
            else:
                new_op = cls(srcOrb, sinkOrb, srcAtom, sinkAtom, radical, pericyclic)
            times.record_time("init", time.time() - start_time)

            start_time = time.time()
            pass_filters2: bool = cls.pass_filters(new_op)
            times.record_time("passFilters2", time.time() - start_time)
            if pass_filters2:
                if rand_one:
                    return [new_op]
                orbPairObjectList.append(new_op)

        # print(times.get_time_str())
        return orbPairObjectList
    
    @staticmethod
    def labelOrbPair(src, sink):
        """Convenience to turn an orb pair into a string tuple"""
        clearAtomMaps(src.mol)
        if not isBondDissociation(src, sink):
            src.labelOrbitalAtoms(10)   
        sink.labelOrbitalAtoms(20)

    def get_masses(self, smi_str):
        mass_list = []
        smi_list = smi_str.split('.')
        mol = OEGraphMol()
        for smi in smi_list:
            OEParseSmiles(mol, smi)
            mass = exact_mass(mol)
            mass_list.append( round(mass, 1) )
            mol.Clear()

        return mass_list

    def masses_to_string(self, mass_list):
        # map elements to string type and join them
        return ", ".join(map(str, mass_list))

    def sum_atom_pred_values(self):
        if not self.atom_pred_vals_sum:
            self.atom_pred_vals_sum = self.srcAtom.fill_pred_value + self.sinkAtom.unfill_pred_value

        return self.atom_pred_vals_sum

