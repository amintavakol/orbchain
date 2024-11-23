"""
Module for tools and methods for a BaseReaction
"""
import logging as log
log.basicConfig(level=log.DEBUG)

from chemutils.Common.MolExt import getAtmMappedSmilesFromCompositeSmiles, clearAtomMaps
from chemutils.Common.MolExt import clearMolStereo
from chemutils.Common.Util import molBySmiles, splitCompositeSmilesToList
from chemutils.Common.Util import createAtomMapSmiString, clearAtomMapsSmiStr
from chemutils.Common.Util import standardizeSmiles, createStdAtomMapSmiString, smi_to_unique_smi_fast 
from chemutils.Common.CanonicalAtomMapSmiles import canonicalizeAtomMapSmiString
from chemutils.CombiCDB.OrbitalInteraction import moveOrbitalElectrons, undoMoveOrbitalElectrons
from chemutils.CombiCDB.OrbitalInteraction import perceiveInvertedOrbitals
from chemutils.CombiCDB.OrbitalInteraction import ElementaryStepData, reactionSmilesFromOrbitalPair
from chemutils.CombiCDB.OrbitalInteraction import isBondDissociation
from chemutils.CombiCDB.MechanismModel import ElectronArrow
from chemutils.Common.OrbitalModel import Orbital as chemutilsOrbital, clearNonInvolvedAtomMaps
from chemutils.CombiCDB.Const import REACTION_DELIM
from chemutils.CombiCDB.MarvinArrowLookup import arrowToOrbital, ArrowConversionError

from openeye.oechem import *
from random import shuffle

START_ORB_LABEL = 10
START_SINK_ORB_LABEL = 20

class Reaction():
    """Abstract class for a reaction.

    This is largely based on a pair of orbitals (Source, Sink).
    Future versions of this could be more complex (ie, pericyclic couldTest
    have a 3rd 'glue' orbital)

    When looking at reactions to compare, we need to be sure to grab all
    reactions that are from all the types of similar MolPairs.  For
    example, if we have a reaction with a MolPair(a, b), then we need to
    compare all reactions from MolPair's (a, b), (a,None), (a,a), (b,
    b), (b, None). However, the logic of this should be handled
    ELSEWHERE.

    In this class, we will consider it to have a MolPair of (a, None) if
    intramolecular, a MolPair of (a, a) if homo-intermolecular, and (a,
    b) if hetero-intermolecular.  We will only store atom mapped smiles
    of what is actually involved in the reaction.    
    """
    def __init__(self, smirks, arrows, *args, **kwargs):
        if ">>" in smirks:
            log.info("good format! no reagent labeled.")
        else:
            raise ValueError("do not label reagents.")
            log.error("If the reaction is in correct format, you probably are including the reagents.")
    
        reactants, products = smirks.split(">>")
 
        self.atm_map_reactant_smiles = None
        self.arrwos = None
        
        self.srcOrb = None
        self.sinkOrb = None
        self.num_electrons = None
        self.is_bond_dissociation = None

        self.srcatom = None
        self.sinkatom = None

        self.is_intermolecular = None

        self.spectators, _, _ = self.find_spectators(reactants, products)
        
        products_smiles = smi_to_unique_smi_fast(products).split(".")
        reactive_products = [i for i in products_smiles if i not in self.spectators.split(".")]
        self.products_smiles = smi_to_unique_smi_fast('.'.join(reactive_products))

        
        # at this time, we do not process different reaction conditions.
        self.reaction_condition = {
                              'temperature': 298,
                              'anion_solvation_potential': 1.0,
                              'cation_solvation_potential': 1.0,
                              'photochemical': 0.0,
                              'lewis_acid_catalyst': 0.0
                              }
        
        self.category='unknown' # for now

        self.setup_reaction(reactants, arrows)
        self.arrows = self.arrows()

    @property
    def src_orbital_str(self):
        """Simple property accessor"""
        if self.src_orbital:
            return self.src_orbital.orbital_str
        return ''

    @property
    def sink_orbital_str(self):
        """Simple property accessor"""
        if self.sink_orbital:
            if not self.is_bond_dissociation:
                return self.sink_orbital.orbital_str_as_sink
            else:
                return self.sink_orbital.orbital_str
        return ''

    @property
    def aromatic_atm_map_smiles(self):
        """The atm_map_reactant_smiles aromatized.  This is necessary for checking equivalence of reactions"""
        return clean_atom_map_smiles(self.atm_map_reactant_smiles, doKekule=False, doArom=True)
    

    def create_full_smirks(self):
        reactants = self.atm_map_reactant_smiles
        products = self.products_smiles
        spectators = self.spectators 
        arrows = self.arrows
        left = reactants + "." + spectators
        left = left.split(".")
        right = products + "." + spectators
        right = right.split(".")
        shuffle(right)
        shuffle(left)
        smirks = ".".join(left) + ".".join(right) + " " + arrows

        return smirks
    
    
    @staticmethod
    def find_spectators(reactants, products):
        canon_reactants = [smi_to_unique_smi_fast(i) for i in reactants.split(".")]
        canon_products = [smi_to_unique_smi_fast(i) for i in products.split(".")]
        
        specs = list(set(canon_reactants) & set(canon_products))
        
        reactants = [i for i in canon_reactants if i not in specs]
        products = [i for i in canon_products if i not in specs]

        reactant_smiles = '.'.join(reactants)
        product_smiles = '.'.join(products)
        specs_smiles = '.'.join(specs)
        
        return smi_to_unique_smi_fast(specs_smiles), smi_to_unique_smi_fast(reactant_smiles), smi_to_unique_smi_fast(product_smiles)

    def _toElementaryStepData(self):
        """Function to convert to chemutilsOrbital format
        """
        oemol = molBySmiles(str(self.atm_map_reactant_smiles))

        elemStepData = ElementaryStepData()
        elemStepData['compositeMol'] = oemol
        elemStepData['filledOrb'] = self.srcOrb
        elemStepData['unfilledOrb'] = self.sinkOrb
        return elemStepData

    def toArrowCodes(self):
        """Convenience to return arrow codes for this reaction"""
        elemStepData = self._toElementaryStepData()
        arrowObjList = []
        bond = moveOrbitalElectrons(elemStepData['filledOrb'], elemStepData['unfilledOrb'], arrowObjList=arrowObjList)
        return ElectronArrow.prepareArrowCodes(arrowObjList)

    def arrows(self):
        """Better name. For use in django templates etc."""
        return self.toArrowCodes()
    

    def equivalentReaction(self, rxn):
        """Predicate to see if the rxn is equivalent to self.

        We are not overriding __eq__ because we want to maintain address and dbid level equality checks for other
        purposes.  This is mainly to see, is a new BaseReaction functionally equivalent to self.

        This is determined by equality of:  mol_pair, reaction_conditions, iso_product_smiles,
        src_orbital_str, sink_orbital_str, aromatic_atm_map_smiles
        """
        equiv = self.reaction_condition == rxn.reaction_condition
        equiv = equiv and self.src_orbital_str == rxn.src_orbital_str
        equiv = equiv and self.sink_orbital_str == rxn.sink_orbital_str
        equiv = equiv and self.aromatic_atm_map_smiles == rxn.aromatic_atm_map_smiles
        return equiv
    
    def __unicode__(self):
        """Make a unicode representation of a reaction

        Note: handle the case when we have some blank things
        """
        complete = hasattr(self, 'src_orbital') and hasattr(self, 'sink_orbital') and self.iso_product_smiles;
        
        if not complete:
            return 'INCOMPLETE: %s' % self.atm_map_reactant_smiles;

        return '%s>>%s (%s, %s) (%s)' % (self.atm_map_reactant_smiles,
                                      self.product_smiles,
                                      self.src_orbital_str,
                                      self.sink_orbital_str,
                                      self.reaction_condition)
    
    def init_cond(self):
        """
        Currently we are reading the initial condition from the file.
        ltimately, this must be implemented here so we extract the
        intial condition from the Reaction object
        """
        init_cond = "unknown"
        self.init_cond = init_cond
        return init_cond


    def cat1(self):
        """
        Currently we are reading the category 1 from the file.
        ltimately, this must be implemented here so we extract the
        category 1 from the Reaction object
        """
        cat1 = "unknown"
        self.cat1 = cat1
        return cat1

    def cat2(self):
        """
        Currently we are reading the categor 2 from the file.
        ltimately, this must be implemented here so we extract the
        category 2 from the Reaction object
        """
        cat2 = "unknown"
        self.cat2 = cat2
        return cat2

    @staticmethod
    def productSmiFromArrowCodes(atm_map_reactant_smiles, arrowCodes):
        """Convenience function to make the actual output smiles from a given input smiles and arrowCodes

        Returns a tuple with the canonical input and canonical output smiles
        """
        input_smi = canonicalizeAtomMapSmiString(atm_map_reactant_smiles)
        mol = molBySmiles(input_smi)
        arrowList = ElectronArrow.parseArrowCodes(mol, arrowCodes)
        ElectronArrow.applyAll(mol, arrowList)
        output_smi = createStdAtomMapSmiString(mol)
        output_smi = canonicalizeAtomMapSmiString(output_smi)
        return (input_smi, output_smi)
    
    def setup_reaction(self, atm_map_reactant_smiles, arrowCodes, *args, **kwargs):
        """Process a given reaction.

        Note: Right now, this does not account for pericyclics.

        Note as well that this does NOT handle lookup if the particular
        reaction already exists in the DB.  Caller is expected to know
        what to do with the returned instance.

        There is the option for NOT saving the underlying components.
        This is meant as a way to save on hitting the db.  Be careful though,
        if you choose not to save, but then use like it is saved, the results
        can't be determined.

        NOTE that e.g. srcOrb actually contains the entire composite mol
        """
        
        srcOrbStr, sinkOrbStr, srcOrb, sinkOrb, nElectrons = arrowToOrbital(atm_map_reactant_smiles, arrowCodes)
       
        #self.nElectrons = nElectrons

        is_bond_dissociation = isBondDissociation(srcOrb, sinkOrb)
        self.is_bond_dissociation = is_bond_dissociation
        clearMolStereo(srcOrb.mol)
        clearAtomMaps(srcOrb.mol)

        #store here the set of disconnected smi components, to later check if we have spectators. A.B.C to set(A,B,C)
        component_mols_str_set = set(OEMolToSmiles(srcOrb.mol).split('.'))

        if is_bond_dissociation:
            sinkOrb.labelOrbitalAtoms(START_ORB_LABEL)
        else:
            srcOrb.labelOrbitalAtoms(START_ORB_LABEL)
            sinkOrb.labelOrbitalAtoms(START_SINK_ORB_LABEL)
        
        initSmi = createStdAtomMapSmiString(srcOrb.mol)
        atm_map_reactant_smiles = getAtmMappedSmilesFromCompositeSmiles(initSmi, clearMaps=False)
        smiList = splitCompositeSmilesToList(atm_map_reactant_smiles)
        isInterMolecular = len(smiList) > 1

        self.is_intermolecular = isInterMolecular

        self.srcOrb = srcOrb
        self.sinkOrb = sinkOrb
        
        # Make a new copy of the orbs to ensure only working with atm
        # mapped components
        newOEMol = molBySmiles(atm_map_reactant_smiles)
        nSrc = chemutilsOrbital.fromLabeledMolAndInfoStr(newOEMol,
                                                      srcOrb.formatInfoStr())
        nSink = chemutilsOrbital.fromLabeledMolAndInfoStr(newOEMol,
                                                     sinkOrb.formatInfoStr())
        
        ## Then make sure we have a canonincal version of the atom mapped smiles
        clearNonInvolvedAtomMaps([nSrc, nSink])
        atm_map_reactant_smiles = createStdAtomMapSmiString(newOEMol)
        
        self.atm_map_reactant_smiles = atm_map_reactant_smiles

        # Next figure out the orbitals
        clearAtomMaps(newOEMol)
        nSrc.labelOrbitalAtoms(START_ORB_LABEL)
        srcAtmSmi, srcInfoStr = nSrc.toLabeledSmiAndInfoStr()
        srcAtmSmi = getAtmMappedSmilesFromCompositeSmiles(srcAtmSmi, clearMaps=False)  # here is where the non-connected components are dropped. part w map is saved
        
        srcOEMol = molBySmiles(srcAtmSmi)
        singleMolSrcOrb = chemutilsOrbital.fromLabeledMolAndInfoStr(srcOEMol, srcInfoStr)
        if singleMolSrcOrb.isPerceivedFreeRadical():
            num_elec_src = 1
        else:
            num_elec_src = 2

        self.srcAtom = self.clean_other_non_reactive_atoms(srcAtmSmi)

        
        clearAtomMaps(newOEMol)
        nSink.labelOrbitalAtoms(START_ORB_LABEL)
        sinkAtmSmi, sinkInfoStr = nSink.toLabeledSmiAndInfoStr()
        sinkAtmSmi = getAtmMappedSmilesFromCompositeSmiles(sinkAtmSmi, clearMaps=False)
       
        sinkOEMol = molBySmiles(sinkAtmSmi)
        singleMolSinkOrb = chemutilsOrbital.fromLabeledMolAndInfoStr(sinkOEMol, sinkInfoStr)
        if singleMolSinkOrb.isPerceivedFreeRadical():
            num_elec_sink = 1
        else:
            num_elec_sink = 2

        self.sinkAtom = self.clean_other_non_reactive_atoms(sinkAtmSmi)

        # Then figure out the products
        reactionSmi = reactionSmilesFromOrbitalPair(nSrc, nSink)
        reactantSmi, productSmi = reactionSmi.split(REACTION_DELIM)
        iso_product_smiles = standardizeSmiles(clearAtomMapsSmiStr(productSmi), aromatize=True)

        self.iso_product_smiles = iso_product_smiles

        self.num_electrons = min(num_elec_src, num_elec_sink)

    def clean_other_non_reactive_atoms(self, smi):
        mol = molBySmiles(smi)
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() == START_ORB_LABEL:
                pass
            else:
                atom.SetMapIdx(0)
        
        return OEMolToSmiles(mol)

    @classmethod
    def get_or_create_from_object(cls, obj):
        """Convenience to lookup or insert based on some common attr
        
        """
        obj.clean()
        return cls.objects.get_or_create(atm_map_reactant_smiles=obj.atm_map_reactant_smiles,
                                         iso_product_smiles=obj.iso_product_smiles,
                                         reaction_conditions=obj.reaction_conditions,
                                         src_orbital=obj.src_orbital,
                                         sink_orbital=obj.sink_orbital,
                                         mol_pair=obj.mol_pair,
                                         is_bond_dissociation=obj.is_bond_dissociation,
                                         num_electrons=obj.num_electrons)

    @classmethod
    def get_from_object(cls, obj):
        """Convenience to lookup or insert based on some common attr
        
        """
        obj.clean()
        return cls.objects.get(atm_map_reactant_smiles=obj.atm_map_reactant_smiles,
                                         iso_product_smiles=obj.iso_product_smiles,
                                         reaction_conditions=obj.reaction_conditions,
                                         src_orbital=obj.src_orbital,
                                         sink_orbital=obj.sink_orbital,
                                         mol_pair=obj.mol_pair,
                                         is_bond_dissociation=obj.is_bond_dissociation,
                                         num_electrons=obj.num_electrons)
    
