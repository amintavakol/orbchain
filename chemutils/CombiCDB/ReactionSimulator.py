import sys, os;
from optparse import OptionParser;
from sets import Set;
import math;
from openeye.oechem import OEGraphMol, OEParseSmiles;
from openeye.oechem import OEAddExplicitHydrogens;

from chemutils.Common.Const import SMILES_MOL_DELIM, REACTION_COMPONENT_DELIM;
from chemutils.Common.Util import molBySmiles;
from chemutils.Common.Util import ProgressDots, stdOpen;
from chemutils.Common.Util import EyringFormula;
from chemutils.Common.Util import createAtomMapSmiString, createStandardSmiString, standardizeSmiles;
from chemutils.Common.IteratorFactory import oemolistreamFactory;
from chemutils.Common.MolExt import grossFormalCharges;
from chemutils.Common.OrbitalModel import Orbital, orbitalIter;
from chemutils.Common.OrbitalModel import ORBITAL_LABEL_INCREMENT;
from chemutils.CombiCDB.MechanismModel import ElectronArrow;
from chemutils.CombiCDB.OrbitalInteraction import moveOrbitalElectrons, undoMoveOrbitalElectrons;
from chemutils.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair, compatibleOrbitalPair;
from chemutils.CombiCDB.OrbitalInteraction import ElementaryStepData;
from chemutils.score.orbital.TransitionStateEnergyScore import TransitionStateEnergy, FreeEnergyChange;
from chemutils.score.orbital.ReactivityScore import ReactiveOrbitalFactory;
from chemutils.score.orbital.ReactivityScore import OrbitalReactivityProbability, AcceptAllOrbitalReactivityProbability;
from MechanismModel import clearMechanismLabels;
from ReactionProposer import ReactionProposer;
from chemutils.Common.MolStdValue import ROOM_TEMPERATURE;

from chemutils.CombiCDB.Const import VALUE_DELIM;
from chemutils.CombiCDB.Util import log;

DEFAULT_MAX_RECURSIONS = 1;

"""Maximum fraction of a reaction mixture component that 
can be removed in a single iteration.
"""
SIMULATION_ITERATION_MAX_FRACTION = 0.5;
"""Minimum unit of a reaction mixture component that
will be removed in a single iteration, if it is going to 
react at all.
"""
SIMULATION_ITERATION_MIN_UNIT = 5;

"""A default volume scalar to divide reaction component quantities by
when calculating rate values.  In terms of relative rate values,
this will cancel out when comparing any reactions of the same molecularity.
The primary effect it has is to make unimolecular reactions much faster
than bimolecular reactions (since bimolecular reactions have to divide their
rate by this scalar a second time).
"""
DEFAULT_VOLUME_SCALAR = 100;

class ReactionSimulator:
    """Driver class to run a reaction simulation.  Given some a starting reaction mixture and 
    conditions, try to predict the resulting products of the reaction, including intermediates,
    mechanisms and relative energy levels.
    """
    maxRecursions = None;
    freeEnergyChangeCalc = None;
    transitionEnergyCalc = None;
    
    def __init__(self):
        """Default constructor."""
        self.maxRecursions = DEFAULT_MAX_RECURSIONS;
        # These calculation modules can be set to customized objects,
        #   otherwise will just setup some default ones
        self.freeEnergyChangeCalc = None;
        self.transitionEnergyCalc = None;
    
    def runSimulation(self, reactionMixture, depth=0, reactionProposerParams=None):
        """Primary driver function to carry out a general simulation,
        but mostly defers to other functions and calculator objects
        to do the real work.
        """
        
        ## Step (0) - set up the energy calculations
        if self.freeEnergyChangeCalc is None:
            self.freeEnergyChangeCalc = FreeEnergyChange.getInstance();
        
        if self.transitionEnergyCalc is None:
            self.transitionEnergyCalc = TransitionStateEnergy.getInstance();    # Should propagate reactionMixture.temperature parameter?
        
        # Step (1,2) - Prepare reaction proposer and get generator
        # (2) Prepare all possible pair-wise reactive orbital possibilities,
        #       adjusting for orbitals on different molecules, considering
        #       intramolecular and intermolecular possibilities
        #       Check previous history as well to avoid repeating previous orbital pair analysis
        # NOTE:  This step is wrapped up in the reactionProposer.__call__ below.., 
        reactionProposer = ReactionProposer.getInstance(reactionProposerParams);
        
        
        # (3) Score all (new) pair-wise reactive orbital possibilities
        for (filledOrb, unfilledOrb) in reactionProposer(reactionMixture):
            # Check if the filled orb / unfilled orb pair is already
            #   represented in the elementaryStepHistory.  
            #   If so, no need to reanalyze the pair nor overwrite the change history
            reactionSmiles = reactionSmilesFromOrbitalPair(filledOrb, unfilledOrb);

            if reactionSmiles not in reactionMixture.elementaryStepHistory:
                elementaryStep = self.analyzeOrbitalPair( filledOrb, unfilledOrb, reactionMixture );
                if elementaryStep.isProductive():
                    reactionMixture.elementaryStepHistory[reactionSmiles] = elementaryStep;
        
        # (4) Apply a reaction "iteration" by altering the contents and quantities in the
        #       reaction mixture based on the scored orbital interactions
        reactionMixture = self.applyElementaryStepsToReactionMixture( reactionMixture, reactionMixture.elementaryStepHistory );
        
        # (5) Evaluate simulation stopping criteria (reaction mixture stabilized, or max recursions reached)
        stoppingCriteriaMet = self.simulationStoppingCriteriaMet( reactionMixture, depth );
        
        # (6) Recursively continue simulation if no stopping criteria met.
        if not stoppingCriteriaMet:
            self.runSimulation( reactionMixture, depth+1, reactionProposerParams=reactionProposerParams );

    def analyzeOrbitalPair( self, filledOrb, unfilledOrb, reactionMixture ):
        """Given a filled and unfilled orbital and a scoring function,
        calculate the interaction score for the pair (transition state energy
        and fundamental reaction rate constant), as well as the SMILES
        string representations of the (composite) molecules before and
        after the orbital interaction to capture all of the core information
        of interest regarding this elementary step.
        The reaction SMILES string should contain labels for the orbital locations,
        so that the collection of all this data can be used to key a 
        unique string representing the orbital interaction.
        """
        elementaryStep = ElementaryStepData();

        # Assume the filledOrb and unfilledOrb are based on the same compositeMol object
        elementaryStep["compositeMol"] = filledOrb.mol;
        elementaryStep["filledOrb"] = filledOrb;
        elementaryStep["unfilledOrb"] = unfilledOrb;
        
        arrowObjList = [];
        reactionSmiles = reactionSmilesFromOrbitalPair(filledOrb, unfilledOrb, arrowObjList=arrowObjList);
        elementaryStep["reactionSmiles"] = reactionSmiles;
        elementaryStep["arrowStr"] = ElectronArrow.prepareArrowCodes(arrowObjList);

        reactionMol = OEGraphMol();    
        OEParseSmiles(reactionMol, elementaryStep["reactionSmiles"] );
        clearMechanismLabels( reactionMol );

        # Separate reaction molecule into components
        unlabeledRxnSmi = createStandardSmiString(reactionMol);
        tokens = unlabeledRxnSmi.split(REACTION_COMPONENT_DELIM);
        reactantsSmi = tokens[0];
        reagentSmi = tokens[1];
        productsSmi = tokens[2];
        elementaryStep["reactantSmilesList"] = reactantsSmi.split(SMILES_MOL_DELIM);
        elementaryStep["productSmilesList"] = productsSmi.split(SMILES_MOL_DELIM);

        # Copy over conditions parameters for scoring function
        elementaryStep["temperature"] = reactionMixture.temperature;
        elementaryStep["externalLight"] = reactionMixture.externalLight;
        elementaryStep["cationSolvationPotential"] = reactionMixture.solventModel.cationSolvationPotential;
        elementaryStep["anionSolvationPotential"] = reactionMixture.solventModel.anionSolvationPotential;
        
        elementaryStep["freeEnergyChange"] = self.freeEnergyChangeCalc(elementaryStep);
        elementaryStep["transitionEnergy"] = self.transitionEnergyCalc(elementaryStep);

        historyLength = 0;
        if reactionMixture.currentIteration is not None:
            historyLength = reactionMixture.currentIteration;

        elementaryStep["changeHistory"] = [0] * historyLength;

        return elementaryStep;

    def applyElementaryStepsToReactionMixture( self, reactionMixture, elementaryStepHistory ):
        """Apply a reaction "iteration" by altering the contents and quantities in the
        reaction mixture based on the scored orbital interactions.
        """
        
        # Increment a current slot in reaction mixture history
        # Startup current slot as a copy of the previous slot
        for molSmi in reactionMixture:
            currentQuantity = reactionMixture[molSmi].currentQuantity();
            reactionMixture[molSmi].appendQuantity( currentQuantity );
        reactionMixture.currentIteration += 1;

        # Figure out how much to modify the reaction mixture components by.
        simulationScale = self.calculateSimulationScale( reactionMixture, elementaryStepHistory );
        reactionMixture.scaleHistory.append( simulationScale );
        
        # Calculate all simulations deltas before actually applying any of them
        #   to ensure that all changes effectively happen simultaneously
        simDeltaByReactionSmiles = dict();
        for reactionSmiles, elementaryStep in elementaryStepHistory.iteritems():
            # Pre-calculated simulation scale to achieve an interesting (not too big / small)
            #   simulation iteration, and prevent molecule quantities dropping below 0.
            simulationDelta = self.calculateSimulationDelta( elementaryStep, reactionMixture, simulationScale );

            simDeltaByReactionSmiles[reactionSmiles] = simulationDelta;
        
        for reactionSmiles, simulationDelta in simDeltaByReactionSmiles.iteritems():
            elementaryStep = elementaryStepHistory[reactionSmiles];
            
            if simulationDelta != 0:
                for reactantSmi in elementaryStep["reactantSmilesList"]:
                    currentQuantity = reactionMixture[reactantSmi].currentQuantity();
                    reactionMixture[reactantSmi].updateQuantity( currentQuantity - simulationDelta );

                for productSmi in elementaryStep["productSmilesList"]:
                    # Standardize the smiles here. (ensureAtomMaps)
                    productSmi = standardizeSmiles(productSmi, ensureAtomMaps=True);
                    if productSmi not in reactionMixture:
                        # Newly generated product molecule,  
                        # Initialize it's history to all 0
                        reactionMixture.addMoleculeQuantity( productSmi, 0 );
                    currentQuantity = reactionMixture[productSmi].currentQuantity();
                    reactionMixture[productSmi].updateQuantity( currentQuantity + simulationDelta );

            # Record deltas so we can reconstruct the reaction pathway history later
            elementaryStep["changeHistory"].append( simulationDelta );
        
        return reactionMixture;

    def calculateSimulationScale( self, reactionMixture, elementaryStepHistory ):
        """Calculate a scalar to translate reaction rates into actual
        simulation delta values to apply to the reaction mixture contents.
        
        Reviews the current state of the reaction mixture and the 
        calculated orbital pair analysis for all of the relevant possible
        elementary reaction steps to figure out a globally useful scale.
        
        Scale should not be too large or else important simulation
        dynamics may be missed and should not be too small, or else
        too many iterations will be wasted on unproductive 
        (minimally productive) steps.  Steps also cannot be too large
        or else they may remove too many of a reacting molecule,
        resulting in an unreal negative quantity.
        """
        # First just assume no scaling factor and sum up the net reaction rates
        #   affecting every mixture component.  This has the important effect
        #   of accounting for dynamic equilibrium.  If one reaction would increase
        #   a component quantity, but another one would decrease it, only the
        #   net change really matters.
        netRatePerComponentSmi = dict();
        for reactionSmiles, elementaryStep in elementaryStepHistory.iteritems():
            rate = self.calculateSimulationDelta( elementaryStep, reactionMixture );
            for reactantSmi in elementaryStep["reactantSmilesList"]:
                if reactantSmi not in netRatePerComponentSmi:
                    netRatePerComponentSmi[reactantSmi] = 0;
                netRatePerComponentSmi[reactantSmi] -= rate;

            for productSmi in elementaryStep["productSmilesList"]:
                if productSmi not in netRatePerComponentSmi:
                    netRatePerComponentSmi[productSmi] = 0;
                netRatePerComponentSmi[productSmi] += rate;

        # Now figure out which is the limiting reactant.  This would be
        #   one where the net rate is negative (the reactant is being consumed),
        #   with the largest (magnitude) ratio of net rate / current quantity.
        limitingRatio = None;
        limitingComponentSmi = None;
        
        for componentSmi, netRate in netRatePerComponentSmi.iteritems():
            if netRate < 0:
                currentQuantity = reactionMixture[componentSmi].currentQuantity();
                if currentQuantity > 0:
                    # If none of the component currently exists, there will be no reaction that uses it
                    componentRatio = abs(netRate / currentQuantity);
                    if limitingRatio is None or componentRatio > limitingRatio:
                        limitingRatio = componentRatio;
                        limitingComponentSmi = componentSmi;

        simulationScale = 1.0;
        if limitingComponentSmi is not None:
            # Adjust simulation scale such that the limiting component is removed
            #   by up to 50% of its current quantity (or other specified fraction),
            #   but at least 1 unit (or other specified amount).
            # Actual simulation iteration delta is based on product of the net rate
            #   and current quantity.
            limitRate = netRatePerComponentSmi[limitingComponentSmi];
            limitQuantity = reactionMixture[limitingComponentSmi].currentQuantity();
            standardDelta = limitRate * limitQuantity;

            if (limitQuantity * SIMULATION_ITERATION_MAX_FRACTION) < SIMULATION_ITERATION_MIN_UNIT:
                # The current quantity is small enough that changing the max fraction is
                #   still less than the minimum unit, so just change by the minimum unit then
                # Set scalar to satisfy equation:
                #   scledDelta = scalar * limitRate = SIMULATION_ITERATION_MIN_UNIT;
                # But, ensure we don't remove more than is available
                minimumQuantity = float(min(limitQuantity,SIMULATION_ITERATION_MIN_UNIT));
                simulationScale = minimumQuantity / limitRate;
            else:
                # The current quantity is large enough that we can safely change
                #   by the max fraction.  Set simulation scalar to satisfy equation:
                #   scaledDelta = limitQuantity * MAX_FRACTION = scalar * limitRate
                simulationScale = SIMULATION_ITERATION_MAX_FRACTION * limitQuantity / limitRate;
        
            #raise Exception(simulationScale, limitQuantity*limitRate)

        simulationScale = abs(simulationScale);
        return simulationScale;
    
    def calculateSimulationDelta( self, elementaryStep, reactionMixture, simulationScale=1.0 ):
        """Given the reactivity analysis of an orbital pair and the overall reaction mixture,
        determine how many units the reaction mixture should be modified by this 
        orbital interaction in the next simulation iteration.
        
        Scale simulation so will not alter too much of the content at one time, otherwise
        may overshoot reasonable simulation intervals.  Do not want to do all small
        simulations either because of wasted computation / simulation cycles.
        
        Mass action kinetic rate laws are normally based on molecule concentrations,
        not absolute amounts.  This qualitatively corresponds to the probability
        of encountering any particular type of molecule in the mixture.  
        For this simulation then, simulation deltas may need to be scaled by the molecule
        amounts over the total amount of all molecules in the mixture.
        """
        # Eyring equation extension of Arrhenius equation as a 
        #   formal theoretical basis for translating transition state 
        #   energies into rate constants.
        energy2rateConst = EyringFormula(reactionMixture.temperature);
        rateConstant = energy2rateConst( elementaryStep["transitionEnergy"] );
        
        # Law of mass action takes product of rate constant and concentrations 
        #   of all reacting molecules to produce an overall reaction rate.
        # Note that if the current quantity of a reactant is 0 (no more left),
        #   this will naturally supress the reaction rate to 0.
        rate = rateConstant;
        for reactantSmi in elementaryStep["reactantSmilesList"]:
            currentQuantity = 0;
            if reactantSmi in reactionMixture:
                currentQuantity = reactionMixture[reactantSmi].currentQuantity();

            rate *= currentQuantity;
            # Really should be concentration, not quantity, so divide by a fixed
            #   volume scalar to help ensure that intramolecular reactions are favored
            rate /= reactionMixture.volumeScalar;
            
        
        # Scale the reaction rate by the simulation scale being used.
        # The idea is to translate reaction rates into actual delta values
        #   to be applied to reaction mixture contents.
        simDelta = (rate * simulationScale);
        simDelta = int(round(simDelta));    # Round to the nearest integer
        
        return simDelta;
        
    def simulationStoppingCriteriaMet( self, reactionMixture, depth ):
        """Evaluate simulation stopping criteria (reaction mixture stabilized, or max recursions reached)

        # TO DO: Implement this to check if simulation has stabilized
        """
        stoppingCriteriaMet = False;
        if (depth+1) >= self.maxRecursions:
            stoppingCriteriaMet = True;

        return stoppingCriteriaMet;

class ReactionMixture(dict):
    """Representation of reaction mixture to simulate with.  Extension of dict object,
    keyed by SMILES strings of mixture components, values are a list of counts or other
    population / concentration of the respective component.  Multiple values in the list
    to represent a time series of population data as the reaction mixture changes.
    
    Addition information can be stored on object such as reaction conditions:
    temperature, solvents, etc.
    """
    
    temperature = None;
    externalLight = None;
    volumeScalar = None;
    solventModel = None;
    
    currentIteration = None;
    scaleHistory = None;
    
    # Dictionary keyed by reaction smiles representing elementary step
    #   with data and information on the elementary step,
    #   including a running history of the net simulation deltas
    #   the steps have applied to the reaction mixture contents
    elementaryStepHistory = None;
    
    def __init__(self, initReactionMixture=None):
        """Constructor, with copy option if supply another
        ReactionMixture instance.
        """
        self.temperature = ROOM_TEMPERATURE;
        self.externalLight = False;
        self.volumeScalar = 1.0;
        self.solventModel = SolventModel();
        self.currentIteration = 0;  # Keep track of how much simulation history has been recorded
        self.scaleHistory = []; # Keep track of relative step size of each change history
        self.elementaryStepHistory = dict();

        if initReactionMixture is not None:
            dict.__init__(self, initReactionMixture);
        
            self.temperature = initReactionMixture.temperature;
            self.externalLight = initReactionMixture.externalLight;
            self.volumeScalar = initReactionMixture.volumeScalar;
            self.solventModel = initReactionMixture.solventModel;

    def addMoleculeQuantity( self, smiles, quantity, pastHistory=None ):
        """Add a record in the mixture for the molecule represented
        by the SMILES string and the quantity specified.
        If an existing (history) record for the SMILES already 
        exists, then just update the the last quanity measure to the molecule's history.
        Allow possibility of pastHistory to be in a string list format.
        
        Added in a standardization where the smiles are standardized keeping any existing 
        atom maps and kekulizing if aromatic.
        
        # TO DO: Cope better with redundancy checks.
        #   Separate functions for add molecule vs. update molecule?
        """
        smiles = standardizeSmiles(smiles, ensureAtomMaps=True);
        if smiles not in self:
            # New component, never seen before
            componentData = ReactionMixtureComponent();
            componentData["smiles"] = smiles;
            # Possibly several simulation iterations later, so fill in all needed spots
            componentData["quantityHistory"] = [0] * self.currentIteration; 
            mol = OEGraphMol();
            OEParseSmiles(mol, smiles);
            componentData.mol = mol;
            self[smiles] = componentData;

        componentData = self[smiles];

        if pastHistory is not None:
            # Use this to fill in or replace any current history,
            # Input may be list of ints already, but also accept string list
            if isinstance( pastHistory, str ):
                pastHistoryStr = pastHistory;
                pastHistory = pastHistoryStr.split(VALUE_DELIM);
                for iHistory, historyStr in enumerate(pastHistory):
                    pastHistory[iHistory] = int(historyStr);

            # Consistency check for history length and iteration count
            if len(pastHistory) != self.currentIteration:
                raise Exception('Expected %d elements in reaction quantity history but instead found %s' % (self.currentIteration, pastHistory) );
            componentData["quantityHistory"] = pastHistory;

        componentData.appendQuantity( quantity );

    def listSortedMixtureSmiles( self ):
        """Go through the reaction mixture contents and return
        the list of SMILES string keys sorted in some meaningful order, 
        such as by last quantity, carbon content, etc.
        """
        dataTuples = [];
        for smiles, mixtureElement in self.iteritems():
            # Sort by the total carbon content (number carbons * latest quantity value)
            carbonCount = 0;
            for atom in mixtureElement.mol.GetAtoms():
                if atom.IsCarbon():
                    carbonCount += 1;
            sortCol = ( carbonCount * mixtureElement.currentQuantity() );
            dataTuple = (sortCol, smiles);
            dataTuples.append( dataTuple );
        dataTuples.sort();
        dataTuples.reverse();   # Descending order
        
        orderedSmilesList = [];
        for (sortCol, smiles) in dataTuples:
            orderedSmilesList.append( smiles );

        return orderedSmilesList;

    def scaleHistoryStr(self):
        """List all of the quantity scaling history elements),
        in a single delimited string format.
        """
        scaleHistoryStrList = [];
        for value in self.scaleHistory:
            scaleHistoryStrList.append( str(value) );
        historyStr = str.join(VALUE_DELIM, scaleHistoryStrList);
        return historyStr;

    def parseScaleHistory(self, historyStr):
        """Parse a delimited string representing the quantity scaling
        history and store the results as a list attribute.
        """
        if historyStr is None or historyStr == "":
            self.scaleHistory = [];
        else:
            self.scaleHistory = historyStr.split(VALUE_DELIM);
            for iHistory, historyStr in enumerate(self.scaleHistory):
                self.scaleHistory[iHistory] = float(historyStr);

        # Consistency check for history length and iteration count
        if len(self.scaleHistory) != self.currentIteration:
            raise Exception('Expected %d elements in reaction simulation scaling history but instead found %s' % (self.currentIteration, self.scaleHistory) );

    def clearHistory(self):
        """Keep the current reaction mixture contents, but clear the past histories.
        """
        for componentData in self.itervalues():
            componentData["quantityHistory"] = componentData["quantityHistory"][-1:];
        for elementaryStepData in self.elementaryStepHistory.itervalues():
            elementaryStepData["changeHistory"] = [];
        self.currentIteration = 0;        
        self.scaleHistory = [];

class ReactionMixtureComponent(dict):
    """Capture all of the core information of interest 
    regarding an element of a reaction mixture.
    For the most part, this would just be the SMILES string and
    perhaps mol object representing the molecule element,
    and a quantity value indicating the relative amount in the mixture.
    
    For more advanced tracking featuers, may want information like
    a historical series of quantity values to track the progress
    of the reaction mixture over the course of a multi-step simulation.
    """
    mol = None;
    
    def __init__(self):
        """Default constructor"""
        self["smiles"] = None;
        self["quantityHistory"] = [];
        self.mol = None;

    def currentQuantity(self):
        """Last value in quantity history should be the current value"""
        return self["quantityHistory"][-1];

    def cumulativeQuantity(self):
        """Sum of all values in quantity history.
        Useful to check if all 0, meaning this component really never participated.
        """
        return sum(self["quantityHistory"]);

    def pastHistoryStr(self):
        """List all of the quantity history elements from the past,
        (excludes the current value), in a single delimited string format.
        """
        pastHistoryStrList = [];
        for value in self["quantityHistory"][:-1]:  # Skip the current value as not being a *past* value
            pastHistoryStrList.append( str(value) );
        historyStr = str.join(VALUE_DELIM, pastHistoryStrList ); # Past values
        return historyStr;

    def appendQuantity(self, newQuantity ):
        """Add a new quantity value to the history for this component,
        which will then represent the current quantity.
        """
        self["quantityHistory"].append( newQuantity );

    def updateQuantity(self, newQuantity):
        """Update the value of the current (latest) quantity.
        """
        self["quantityHistory"][-1] = newQuantity;


class SolventModel:
    """Object to contain any solvent model parameters.
    Pretty simple struct for now.
    """
    cationSolvationPotential = None;
    anionSolvationPotential = None;
    
    def __init__(self):
        self.cationSolvationPotential = 1.0;
        self.anionSolvationPotential = 1.0;
    
    

def main(argv):
    """Main method, callable from command line"""
    usageStr =  "usage: %prog [options] <reactantFile> <productFile>\n"+\
                "   For more information, run the following under Python\n"+\
                "   >>> import ReactionSimulator\n"+\
                "   >>> help(ReactionSimulator)\n"
    parser = OptionParser(usage=usageStr)
    (options, args) = parser.parse_args(argv[1:])

    instance = ReactionSimulator()

    if len(args) >= 2:
        molstreamFactory    = oemolistreamFactory(args[0]);
        productFile         = stdOpen(args[1],"w",sys.stdout);
        instance.predictReactionsStream( molstreamFactory, productFile );
    else:
        parser.print_help()
        sys.exit(-1)
    
if __name__=="__main__":
    main(sys.argv)
