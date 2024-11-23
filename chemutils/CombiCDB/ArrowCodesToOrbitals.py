import sys
import os
from optparse import OptionParser;
from pprint import pformat;

from chemutils.Common.OrbitalModel import Orbital, orbitalIter, BOND_ORBITAL_TYPES, bondOrbitalByMapIdx, \
                    atomOrbitalByMapIdx, atomHybridization;
from chemutils.Common.MolExt import clearAtomStereo;
from chemutils.Common.Util import  createAtomMapSmiString
from chemutils.CombiCDB.MechanismModel import ATOM_DELIM, ARROW_DELIM, SINGLE_ARROW, DOUBLE_ARROW;
from chemutils.CombiCDB.MechanismModel import moveAtomElectrons, moveBondElectrons, ElectronArrow;
from chemutils.CombiCDB.MechanismModel import BondDoesNotExistException
from chemutils.CombiCDB.OrbitalInteraction import isBondDissociation

from chemutils.CombiCDB.Util import log

#TODO -- set logging level appropriately for production/testing
#import logging
#log.setLevel(logging.DEBUG)

class ArrowConversionError(Exception):
    """Simple class to capture an error on arrow conversion"""
    def __init__(self, msg, mol=None, arrowList=None, orbTupleList=None):
        """Several ways in which this could be called."""
        #Always have the msg, everything else is optional
        self.msg = msg;
        
        self.mol = mol
        self.arrowList = arrowList
        self.orbTupleList = orbTupleList

        self.extraInfoStr = ''

        ## if we have the mol and arrowList, 
        if self.arrowList is not None:
            self.extraInfoStr = "Arrows: %s" % ElectronArrow.prepareArrowCodes(arrowList)
        elif self.mol is None and self.orbTupleList is not None:
            self.mol = self.orbTupleList[0][0].mol
            self.extraInfoStr =  'orbTupleList : %s' % pformat(
                [[orb.toLabeledSmiAndInfoStr() for orb in theTuple] for theTuple in self.orbTupleList]);
        self.smi = ''
        try:
            self.smi = createAtomMapSmiString(self.mol)
        except Exception:
            ## If mal-formed mol, let pass.
            pass;

    def __str__(self):
        return """ArrowConversionError - %s - %s - %s""" % (self.smi, self.extraInfoStr, self.msg)


def orbitalPairFromArrowStr(mol, arrowStr):
    """General method to return tuple of source/sink orbitals given a mol and arrowStr"""
    arrowList = ElectronArrow.parseArrowCodes( mol, arrowStr )
    return orbitalPairFromArrowList(mol, arrowList);

def orbitalPairFromArrowList(mol, arrowList):
    """General method to return tuple of source/sink orbs given a list of ElectronArrow objs."""
    #First check for single electrons.., and other non-implemented things
    isFreeRadical = False
    for arr in arrowList:
        if arr.nElectrons != 2:
            isFreeRadical = True
            break
    
    #First, canonicalize the arrowList (make the common source/target atoms be in expected places)
    arrowList = ElectronArrow.makeCanonicalArrowList(arrowList);
    log.debug('After Canonical the arrowList is %s ' % str(ElectronArrow.prepareArrowCodes(arrowList)));
    
    if not isFreeRadical:
        sourceOrb, sinkOrb = orbitalPairFromTwoElectronArrowList(mol, arrowList)
        return (sourceOrb, sinkOrb, 2)
    else:
        sourceOrb, sinkOrb = orbitalPairFromRadicalElectronArrowList(mol, arrowList)
        return (sourceOrb, sinkOrb, 1)



def orbitalPairFromRadicalElectronArrowList(mol, arrowList):
    """General method to return tuple of src/sink orbs given list of ElectronArrow objs involving 1e movement.
    
    Idea:  First try to chain together bond pivots.  Then attach bond dissociation or bond formation.
            In the case of NO bond pivots, then either homolytic bond dissociation or dimerization.
    """
    sourceOrb = None;
    targetOrb = None
    if len(arrowList) == 2:
        #homolytic bond dis
        if arrowList[0].isBondDissociationArrow() and arrowList[1].isBondDissociationArrow():
            atom = arrowList[0].sourceAtoms[0];
            bond = atom.GetBond(arrowList[0].sourceAtoms[1])
            sourceOrb = Orbital(atom, None, bond, 1, mol);
            targetOrb = sourceOrb;
        # dimerization
        elif arrowList[0].isBondFormArrow() and arrowList[1].isBondFormArrow():
            sourceOrb = Orbital(arrowList[0].sourceAtoms[0], None, None, electrons=1, mol=mol);
            targetOrb = Orbital(arrowList[1].sourceAtoms[0], None, None, electrons=1, mol=mol);
        else:
            raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Ambigous radical arrows.')
        
    elif len(arrowList) >= 3:
        orderedArrowList = orderRadicalElectronArrowList(mol, arrowList);
        
        #First will tell source
        currArr = orderedArrowList.pop(0);
        atom = currArr.sourceAtoms[0]
        currTargetAtom= None;
        if currArr.isBondFormArrow():
            # Source radical
            sourceOrb = Orbital(atom, None, None, electrons=1, mol=mol);
            currTargetAtom  = currArr.targetAtoms[0] 
        else:
            # Is a source bond
            bond = atom.GetBond(currArr.sourceAtoms[1]);
            sourceOrb = Orbital(atom, None, bond, electrons=1, mol=mol);
            # and the next one is just a completion of the partial form
            currArr = orderedArrowList.pop(0)
            #And we know that the next target.atom will be the Z atom here
            currTargetAtom = currArr.targetAtoms[1]
        
        currArr = orderedArrowList.pop(0);
        currArr = orderedArrowList.pop(0);
        # Here, we may know what the target already is.
        if currTargetAtom is not None:
            atom = currTargetAtom;
            nbrAtom = currArr.sourceAtoms[0]
            if nbrAtom == atom:
                nbrAtom = currArr.sourceAtoms[1]
            bond = atom.GetBond(nbrAtom);
        else:
            atom = currArr.sourceAtoms[1]
            bond = atom.GetBond(currArr.sourceAtoms[0]);
        
        currOrb = Orbital(atom, None, bond, electrons=1, mol=mol)
        targetOrb = currOrb;
        
        #Now the rest
        prevOrb = currOrb;
        while orderedArrowList != []:
            currArr = orderedArrowList.pop(0);
            if currArr.isBondDissociationArrow() or currArr.isBondFormArrow():
                #Target is an atom orb
                atom = currArr.sourceAtoms[0]
                currOrb = Orbital(atom, None, None, electrons=1, mol=mol)
            else:
                currArr = orderedArrowList.pop(0);
                
                atom = currArr.sourceAtoms[0]
                bond = atom.GetBond(currArr.sourceAtoms[1]);
                
                currOrb = Orbital(atom, None, bond, electrons=1, mol=mol)
            prevOrb.extOrbital = currOrb;
            prevOrb = currOrb;
    else:
        raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Not enough radical arrows')
    sourceOrb, targetOrb = canonicalOrbitalPair(sourceOrb, targetOrb, 1);
    return sourceOrb, targetOrb;

def orderRadicalElectronArrowList(mol, arrowList):
    """
    Idea:  First try to chain together bond pivots.  Then attach bond dissociation or bond formation.
            In the case of NO bond pivots, then either homolytic bond dissociation or dimerization.
    """
    pivotArrowList = [arr for arr in arrowList if arr.isBondPivotArrow()];
    bondFormArrowList = [arr for arr in arrowList if arr.isBondFormArrow()];
    bondDisArrowList = [arr for arr in arrowList if arr.isBondDissociationArrow()];
        
    # Should have some chain where we either have a bond form/pivots/bond dis or
    # Bond dis/pivots/bond dis, so sum of these lists should be len 2.
    if len(bondFormArrowList) + len(bondDisArrowList) != 2:
        raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Ambiguous free radical arrows involving bond pivots!');
    
    orderedArrowList = [pivotArrowList.pop()];
    foundOneBack = True
    while pivotArrowList != [] and foundOneBack:
        currYAtom = orderedArrowList[0].sourceIndexes[1];
        foundOneBack = False;
        for i, pivotArr in enumerate(pivotArrowList):
            #Match if the  X or Z matches
            if pivotArr.sourceIndexes[0] == currYAtom or pivotArr.targetIndexes[1] == currYAtom:
                orderedArrowList.insert(0, pivotArrowList.pop(i));
                foundOneBack = True;
                break;
    
    foundOneForward = True;
    while pivotArrowList != [] and foundOneForward:
        currYAtom = orderedArrowList[-1].sourceIndexes[1];
        foundOneForward = False;
        for i, pivotArr in enumerate(pivotArrowList):
            #Match if the  X or Z matches
            if pivotArr.sourceIndexes[0] == currYAtom or pivotArr.targetIndexes[1] == currYAtom:
                orderedArrowList.append(pivotArrowList.pop(i));
                foundOneForward = True;
                break;
    
    if pivotArrowList != []:
        raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Unable to chain radical bond pivot arrows. %s, %s' % \
                    (str([str(p) for p in pivotArrowList]), str([str(p) for p in orderedArrowList])));
    
    # Check that the form and dis actually chain up correctly.
    if len(bondDisArrowList) == 1:
        bondFormArrow = bondFormArrowList[0];
        bondDisArrow = bondDisArrowList[0];
                
        if orderedArrowList[0].sourceIndexes[1] != bondFormArrow.targetIndexes[0] or  \
                orderedArrowList[-1].sourceIndexes[1] != bondDisArrow.sourceIndexes[0]:
            # This could be b/c the ordered list is flipped.
            orderedArrowList.reverse();
            if orderedArrowList[0].sourceIndexes[1] != bondFormArrow.targetIndexes[0] or  \
                    orderedArrowList[-1].sourceIndexes[1] != bondDisArrow.sourceIndexes[0]:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Non-Matching dissociation and formation radical arrows!')
        
        orderedArrowList.insert(0, bondFormArrow);
        orderedArrowList.append(bondDisArrow);
    elif len(bondDisArrowList) == 2:
        # two bond disocciations figure out which goes at the fornt and the back
        bondDisArrow = bondDisArrowList.pop();
        
        if bondDisArrow.sourceIndexes[1] == orderedArrowList[0].sourceIndexes[0]:
            orderedArrowList.insert(0, bondDisArrow);
            bondDisArrow = bondDisArrowList.pop();
            if bondDisArrow.sourceIndexes[0] == orderedArrowList[-1].sourceIndexes[1]:
                orderedArrowList.append(bondDisArrow);
            else:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond dissociation arrow')
        elif bondDisArrow.sourceIndexes[0] == orderedArrowList[-1].sourceIndexes[1]:
            orderedArrowList.append(bondDisArrow);
            bondDisArrow = bondDisArrowList.pop();
            if bondDisArrow.sourceIndexes[1] == orderedArrowList[0].sourceIndexes[0]:
                orderedArrowList.insert(0, bondDisArrow);
            else:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond dissociation arrow')    
        else:
            raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond dissociation arrow')
    elif len(bondFormArrowList) == 2:
        # two bond formations figure out which goes at the front and the back
        bondFormArrow = bondFormArrowList.pop();
        
        if bondFormArrow.targetIndexes[0] == orderedArrowList[0].sourceIndexes[1]:
            orderedArrowList.insert(0, bondFormArrow);
            bondFormArrow = bondFormArrowList.pop();
            if bondFormArrow.sourceIndexes[0] == orderedArrowList[-1].targetIndexes[1]:
                orderedArrowList.append(bondFormArrow);
            else:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond form arrow')
        elif bondFormArrow.sourceIndexes[0] == orderedArrowList[-1].targetIndexes[1]:
            orderedArrowList.append(bondFormArrow);
            bondFormArrow = bondFormArrowList.pop();
            if bondFormArrow.targetIndexes[0] == orderedArrowList[0].sourceIndexes[1]:
                orderedArrowList.insert(0, bondFormArrow);
            else:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond form arrow')    
        else:
            raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous radical bond form arrow')
    
    
    return orderedArrowList;
    
def orbitalPairFromTwoElectronArrowList(mol, arrowList):
    """General method to return a tuple of source/sink orbs given a list of ElectronArrow objs involving 2 e movement.
    """
    log.debug('Before start the arrowList is %s ' % str(ElectronArrow.prepareArrowCodes(arrowList)))
    #log.debug('Len arrowList is %d' % len(arrowList));
    
    #Begin by making pairs of source/sink orbs
    orbTupleList = []
    while arrowList != []:
        arr = arrowList.pop(0)
        log.debug('Just popped this: %s' % str(arr))
        log.debug('After pop, len(arrowList) = %d' % len(arrowList))
        sourceIdx = arr.sourceIndexes
        targetIdx = arr.targetIndexes
        log.debug('arr.sourceIndices %s' % str(arr.sourceIndexes))
        log.debug('arr.targetIndices %s' % str(arr.targetIndexes))
        sourceOrb = None
        targetOrb = None
        if len(sourceIdx) == 1 and len(targetIdx) == 1:
            log.debug('Have a pattern of a1=a2, looking for potential a2,a3=a3.')
            #have a pattern of a1=a2
            #Look for other arrows of pattern a2,a3=a3
            log.debug('ArrowList is %s' % ElectronArrow.prepareArrowCodes(arrowList))
            log.debug(' len(arrowList) = %d' % len(arrowList))

            potIdx = [i for i in xrange(len(arrowList))\
                            if len(arrowList[i].sourceIndexes) == 2 and \
                            len(arrowList[i].targetIndexes) == 1 and \
                            arrowList[i].sourceIndexes[0] == targetIdx[0] \
                            and arrowList[i].sourceIndexes[1] == arrowList[i].targetIndexes[0]];
            log.debug('Got these for potIdx %s' % str(potIdx))
            if len(potIdx) > 1:
                raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous arrow codes!!');
        
            #No matter what the sourceOrb should be a non-bonded orbital.., (any of them should work)
            log.debug('Making source orb of non-bonded orbital on source atom')
            sourceOrb = atomOrbitalByMapIdx(mol, sourceIdx[0], 2);   
            if len(potIdx) == 1:
                log.debug('Have a a2,a3=a3 match, popping and handling')
                #Pop this match
                extArrow = arrowList.pop(potIdx[0])
            
                #orb(a2, a3)
                #Here grab the orbital that is the bond from targetIdx[0] and potIdx[0].sourceIndexes[1]
                log.debug('Finding the bonded orbital from a2 to a3')
                targetOrb = bondOrbitalByMapIdx(mol, targetIdx[0], extArrow.sourceIndexes[1]);
                
            else:
                #orb(a2)
                log.debug('No matching a2,a3=a3 patterns, just the singleton arrow movement.  Make non-bonded target orb')
                targetOrb = atomOrbitalByMapIdx(mol, targetIdx[0], 0)
                
        else:
            #2 source atms..
            # some logic is missing here for Pd...
            if sourceIdx[1] == targetIdx[0] and len(targetIdx) == 1:
                log.debug('have a1,a2=a2 pattern!  Look for a3=a1 or a3,a4=a4,a1 or nothing');
                targetOrb = bondOrbitalByMapIdx(mol, sourceIdx[0], sourceIdx[1]);
                
                
                #First a3=a1 case
                potIdx = [i for i in xrange(len(arrowList)) \
                                if arrowList[i].targetIndexes[0] == sourceIdx[0] and \
                                len(arrowList[i].targetIndexes) == 1];
                
                if len(potIdx) > 1:
                    raise ArrowConversionError(mol=mol, arrowList=arrowList, msg= 'Ambiguous arrow codes!!');
                elif len(potIdx) == 1:
                    log.debug('Have the a3=a1 pattern followed by a1,a2=a2');
                    extArrow = arrowList.pop(potIdx[0])
                    
                    sourceOrb = atomOrbitalByMapIdx(mol, extArrow.sourceIndexes[0], 2);
                   
                else:
                    #Need to look for a3,a4=a4,a1
                    potIdx = [i for i in xrange(len(arrowList)) \
                                if len(arrowList[i].targetIndexes) == 2 and \
                                arrowList[i].targetIndexes[1] == sourceIdx[0] and \
                                arrowList[i].sourceIndexes[1] == arrowList[i].targetIndexes[0] ]
                    
                    if len(potIdx) > 1:
                        raise NotImplementedError('Ambiguous Arrow codes!');
                    elif len(potIdx) == 1:
                        log.debug('Have the a3,a4=a4,a1 followed by a1,a2=a2');
                        extArrow = arrowList.pop(potIdx[0])
                        sourceOrb = bondOrbitalByMapIdx(mol, extArrow.sourceIndexes[1], extArrow.sourceIndexes[0])
                        
                    else:
                        log.debug('Only a1,a2=a2 pattern.  Bond dissociation.')
                        sourceOrb = targetOrb;
            else:
                #a1,a2=a2,a3
                log.debug('Have a1,a2=a2,a3 pattern!')
                log.debug('Looking for the a3,a4=a4 pattern')
                log.debug('len targetIdx %d' % len(targetIdx))
                log.debug('len arrowList %d' % len(arrowList))
                potIdx = [i for i in xrange(len(arrowList)) \
                            if arrowList[i].sourceIndexes[0] == targetIdx[1] and \
                            len(arrowList[i].sourceIndexes) == 2 and \
                            len(arrowList[i].targetIndexes) == 1 and \
                            arrowList[i].sourceIndexes[1] == arrowList[i].targetIndexes[0]];
                log.debug('potIdx %s' % str(potIdx))
                if len(potIdx) > 1:
                    raise ArrowConversionError(mol=mol, arrowList=arrowList, msg='Ambiguous arrow codes!!');
                
                #Find the orbs around a2#
                sourceOrb = bondOrbitalByMapIdx(mol, sourceIdx[1], sourceIdx[0])
                
                if len(potIdx) == 1:
                    log.debug('Have the a3,a4=a4 pattern following a1,a2=a2,a3')
                    extArrow = arrowList.pop(potIdx[0])
                    #Grab the atm and orbIter for a3
                    targetOrb = bondOrbitalByMapIdx(mol, extArrow.sourceIndexes[0], extArrow.targetIndexes[0]);
                else:
                    log.debug('Have a1,a2=a2,a3 on its own.')
                    #Grab the second target atom now.., 
                    targetOrb = atomOrbitalByMapIdx(mol, targetIdx[1], 0)
        
        orbTupleList.append((sourceOrb, targetOrb))
   
    log.debug('sourceOrb is %s' % str(sourceOrb))
    log.debug('targetOrb is %s' % str(targetOrb))

    log.debug('The curr orbTupleList is : %s' % str([(str(pair[0]), str(pair[1])) for pair in orbTupleList]))
    orbTuple = collapseOrbitalPairs(orbTupleList)
    
    actSourceOrbInfoStr = orbTuple[0].formatInfoStr();
    actSourceOrbInfoStr = actSourceOrbInfoStr.replace('\t', ' ')
    actTargetOrbInfoStr = orbTuple[1].formatInfoStr();
    actTargetOrbInfoStr = actTargetOrbInfoStr.replace('\t', ' ')
    
    log.debug('Source After first :%s' % actSourceOrbInfoStr);
    log.debug('Target After first :%s' % actTargetOrbInfoStr);
    
    
    return canonicalOrbitalPair( orbTuple[0], orbTuple[1], 2);


def collapseOrbitalPairs(orbTupleList, numElectrons=2):
    """Function given a list of source/target tuples, collapse into single source/target.
    
    Idea:
        Find the break point, when source.atom/target.atom not bonded.
        Create list of sources/targets
        Add to targets as much as possible,
        Add to sources.
        If anything left over, then Exception.
    """
    if len(orbTupleList) == 1:
        return orbTupleList[0];
    
    #Find the break point, by simply looking for a bond form. :
    potIdx = [];
    for i, orbTuple in enumerate(orbTupleList):
        if orbTuple[0].atom.GetBond(orbTuple[1].atom) is None:
            # Bond form if a bond between orbs has not been made yet.
            potIdx.append(i);
    
    # if found none, then look for a bond break
    if len(potIdx) == 0:
        for i, orbTuple in enumerate(orbTupleList):
            if orbTuple[1].type == 'sigma':
                potIdx.append(i);
    
    # If none found at this point.., then its simply a resonance rearrangement..., Doesn't really matter what the starting point is.
    # log a warning.., then just make it as the first.
    if len(potIdx) == 0:
        log.warning('Cant find a good break point, just using first tuple.  orbTupleList: %s'% 
                            str([[orb.toLabeledSmiAndInfoStr() for orb in theTuple] for theTuple in orbTupleList])); 
        potIdx = [0];
    
    
    # Then decide which of these to use
    if len(potIdx) == 1:
        potIdx = potIdx[0];
    else:
        #More than one, grab the one with the lowest source atomidx
        atomIdxDict = {};
        for i in potIdx:
            atomIdxDict[orbTupleList[i][0].atom.GetMapIdx()]=i;
        sortedKeys = atomIdxDict.keys()
        sortedKeys.sort();
        potIdx = atomIdxDict[sortedKeys[0]];
    
    #Make the source, target lists
    breakOrbTuple = orbTupleList.pop(potIdx);
    log.debug('SourceBreak:%s, and TargetBreak:%s' % (str(breakOrbTuple[0]), str(breakOrbTuple[1])))
    sourceList = [breakOrbTuple[0]];
    targetList = [breakOrbTuple[1]];
    
    contLookingForTargets = targetList[-1].type not in BOND_ORBITAL_TYPES;
    while orbTupleList != [] and contLookingForTargets:
        foundSomething = False;
        currTargetAtmIdx = targetList[-1].atom.GetMapIdx();
        potIdx = [i for i in xrange(len(orbTupleList))\
                    if orbTupleList[i][0].neighbor is not None and \
                    orbTupleList[i][0].neighbor.GetMapIdx() == currTargetAtmIdx];
        if len(potIdx) > 0:
            log.debug('Found a match to chain a target')
            #Found a match, replace last target with flipped source, chain with new target.
            foundSomething = True;
            newOrbTuple = orbTupleList.pop(potIdx[0])
            targetList[-1] = newOrbTuple[0].flippedBondOrbital();
            targetList.append(newOrbTuple[1])
        contLookingForTargets = foundSomething and targetList[-1].type not in BOND_ORBITAL_TYPES;
    
    contLookingForSources = sourceList[-1].type in BOND_ORBITAL_TYPES;
    while orbTupleList != [] and contLookingForSources:
        foundSomething = False;
        currSourceNeigIdx = sourceList[-1].neighbor.GetMapIdx();
        potIdx = [i for i in xrange(len(orbTupleList)) \
                    if orbTupleList[i][1].atom.GetMapIdx() == currSourceNeigIdx and \
                    orbTupleList[i][1].type not in BOND_ORBITAL_TYPES];
        if len(potIdx) > 0:
            log.debug('Found a match to chain a source.')
            foundSomething = True;
            newOrbTuple = orbTupleList.pop(potIdx[0]);
            sourceList.append(newOrbTuple[0]);
        
        contLookingForSources = foundSomething and sourceList[-1].type in BOND_ORBITAL_TYPES;
        
    #Check that nothing is left in the orbTupleList.., 
    if orbTupleList != []:
        ##log.critical('orbTupleList is: %s' % pformat(orbTupleList))
        raise ArrowConversionError(msg='Unable to chain the orbitals. Left over orbs. \n sourceOrbList : %s \n targetList: %s'
                                   % ( str([sourceOrb.toLabeledSmiAndInfoStr() for sourceOrb in sourceList]),
                                       str([targOrb.toLabeledSmiAndInfoStr() for targOrb in targetList])),
                                       orbTupleList = orbTupleList)

       
    sourceOrb = sourceList.pop(0)
    currOrb = sourceOrb; 
    for orb in sourceList:
        currOrb.extOrbital = orb;
        currOrb = orb;
    currOrb.extOrbital = None;
    targetOrb = targetList.pop(0)
    currOrb = targetOrb;
    for orb in targetList:
        currOrb.extOrbital = orb;
        currOrb = orb;
    currOrb.extOrbital = None;
    
    return (sourceOrb, targetOrb);
    
    
def canonicalOrbitalPair(sourceOrb, targetOrb, nElectrons=None):
    """Function to make a canonical pair of orbitals.  
    
    The canonicalization is to take the bond forming point with the lowest source.atom.mapidx
    as the first source/target.  If there are no bond forming points, use a bond break point.
    Then build targets as much as possible, then add sources.
    Cyclic reactions will have many targets and a single source.

    nElectrons is guessed
    """
    if nElectrons is None:
        if (sourceOrb is not None and sourceOrb.isPerceivedFreeRadical()) \
            or targetOrb.isPerceivedFreeRadical():
            nElectrons = 1
        else:
            nElectrons = 2
    #log.critical('nElectrons = %d' % nElectrons)
    #log.critical('sourceOrb = %s' % pformat(sourceOrb.toLabeledSmiAndInfoStr(10)))
    #log.critical('targetOrb = %s' % pformat(targetOrb.toLabeledSmiAndInfoStr(10)))
    #For ease, flatten linked lists
    if sourceOrb is None and nElectrons == 2:
        return sourceOrb, targetOrb;
    # For bond dissociation, the atom should be lower mapIdx
    if nElectrons == 1 and isBondDissociation(sourceOrb, targetOrb):
        idxTuple = (targetOrb.atom.GetMapIdx(), targetOrb.neighbor.GetMapIdx())
        if idxTuple[0] > idxTuple[1]:
            targetOrb = targetOrb.flippedBondOrbital()
            return (targetOrb, targetOrb);
    
    sourceOrbList = [sourceOrb]
    currOrb = sourceOrb
    while currOrb.extOrbital is not None:
        sourceOrbList.append(currOrb.extOrbital);
        currOrb = currOrb.extOrbital;
        
    targetOrbList = [targetOrb];
    currOrb = targetOrb;
    while currOrb.extOrbital is not None:
        targetOrbList.append(currOrb.extOrbital)
        currOrb = currOrb.extOrbital;
    
    #if there is no chaining..., then can just return
    if len(sourceOrbList) == 1 and len(targetOrbList) == 1:
        # With free radical reactions, set the atom orbital to be the source (if one exists).
        if targetOrb.isPerceivedFreeRadical() and sourceOrb.isBondOrbital() and not targetOrb.isBondOrbital():
            return targetOrb, sourceOrb;
        # Dimerization, one with lowest mapIdx should be first
        if targetOrb.isPerceivedFreeRadical() and sourceOrb.isPerceivedFreeRadical() \
                and not sourceOrb.isBondOrbital() and not targetOrb.isBondOrbital():
            atomIdxList = [sourceOrb.atom.GetMapIdx(), targetOrb.atom.GetMapIdx()];
            if atomIdxList[0] > atomIdxList[1]:
                return targetOrb, sourceOrb;
        return sourceOrb, targetOrb;
    
    #First identify bond making points/bond breaking points.
    potBreakSources = [];
    potBreakTargets = [];
    potFormSources = [];
    potFormTargets = [];
    
    #Is the first a break/form:
    ## Breaks are simply sigma bonds, forms are no bonds between two atoms adjacent
    ## in a planned movement.
    if sourceOrb.atom.GetBond(targetOrb.atom) is None:
        potFormSources.append(0);
    if sourceOrb.type == 'sigma':
        potBreakSources.append(0)
    
    for i in xrange(1,len(sourceOrbList)):
        if sourceOrbList[i].atom.GetBond(sourceOrbList[i-1].neighbor) is None:
            potFormSources.append(i);
        if sourceOrbList[i].type == 'sigma':
            potBreakSources.append(i)
    for i in xrange(1, len(targetOrbList)):
        if targetOrbList[i-1].neighbor.GetBond(targetOrbList[i].atom) is None:
            potFormTargets.append(i-1);
        if targetOrbList[i-1].type == 'sigma':
            potBreakTargets.append(i-1);
    
    # Look for the lowest atomMapIdx of a break orbital atom.
    atmIdxDict = {};
    for i in potFormSources:
        atmIdxDict[sourceOrbList[i].atom.GetMapIdx()] = (0, i);
    for i in potFormTargets:
        atmIdxDict[targetOrbList[i].atom.GetMapIdx()] = (1, i);
    
    sortedKeys = atmIdxDict.keys()
    sortedKeys.sort();
    
    # Check that we actually have a new bond being formed.  If NOT, need to 
    # try the breaks.  
    if len(sortedKeys) < 1:
        atmIdxDict = {};
        for i in potBreakSources:
            atmIdxDict[sourceOrbList[i].atom.GetMapIdx()] = (0, i);
        for i in potBreakTargets:
            atmIdxDict[targetOrbList[i].atom.GetMapIdx()] = (1, i);
        sortedKeys = atmIdxDict.keys()
        sortedKeys.sort();
        
    ## if sortedKeys is still empty, then no breaks/forms to be found.
    ## log a warning and return as is
    if len(sortedKeys) < 1:
        log.warning("IN Cononicalization, can't find breaks/forms in %s, %s" % \
                (pformat(sourceOrb.toLabeledSmiAndInfoStr()), pformat(targetOrb.toLabeledSmiAndInfoStr())))
        return sourceOrb, targetOrb
    
    typeId, breakIdx = atmIdxDict[sortedKeys[0]];
    
    if typeId == 0:
        #A source.., 
        newSourceList = sourceOrbList[breakIdx:]
        newTargetList = sourceOrbList[:breakIdx];
        newTargetList.reverse();
        newTargetList = [orb.flippedBondOrbital() for orb in newTargetList];
        newTargetList.extend(targetOrbList);
    else:
        #A target
        newTargetList = targetOrbList[breakIdx + 1:]
        newSourceList = targetOrbList[:breakIdx + 1];
        newSourceList.reverse()
        newSourceList = [orb.flippedBondOrbital() for orb in newSourceList];
        newSourceList.extend(sourceOrbList);
    
    sourceOrbList = newSourceList;
    targetOrbList = newTargetList;
    
    #Fix sources to target if necessary.., 
    if sourceOrbList[-1].neighbor is not None and \
            targetOrbList[-1].type not in BOND_ORBITAL_TYPES and \
            sourceOrbList[-1].neighbor.GetMapIdx() == targetOrbList[-1].atom.GetMapIdx():
        #Here, if the final source can be chained to the final target, then 
        #   can chain all sources to start on end of targets.,  
        newTargets = sourceOrbList[1:];
        newTargets.reverse()
        newTargets = [orb.flippedBondOrbital() for orb in newTargets];
        #Drop original last target.
        targetOrbList = targetOrbList[:-1];
        targetOrbList.extend(newTargets);
        sourceOrbList = sourceOrbList[:1];
        #Make a new 'empty' orbital (really the source of the start)
        targetOrbList.append(Orbital(sourceOrbList[0].neighbor, 'p', None, 0))
    
    # For radicals, source could be target.  Arbitrary, but have the src be the one with the lower atom.MapIdx
    if nElectrons == 1:
        sourceAtomMapIdx = sourceOrbList[0].atom.GetMapIdx();
        targetAtomMapIdx = targetOrbList[0].atom.GetMapIdx();
        if targetAtomMapIdx < sourceAtomMapIdx:
            temp = sourceOrbList;
            sourceOrbList = targetOrbList;
            targetOrbList = temp;
    
    sourceOrb = sourceOrbList.pop(0)
    currOrb = sourceOrb;
    for orb in sourceOrbList:
        currOrb.extOrbital = orb;
        currOrb = orb;
    currOrb.extOrbital = None;
    targetOrb = targetOrbList.pop(0)
    currOrb = targetOrb;
    for orb in targetOrbList:
        currOrb.extOrbital = orb;
        currOrb = orb;
    currOrb.extOrbital = None;
    
    return (sourceOrb, targetOrb);


