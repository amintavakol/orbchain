from pprint import pformat

from openeye.oechem import OESubSearch, OEMatchPairAtom
from openeye.oechem import OEAssignAromaticFlags, OEGraphMol
from chemutils.Common.MolExt import kekulizeMol, deKekulizeMol, findBondByAtomMapIdx
from chemutils.Common.MolExt import createAroBondTypeMap
from chemutils.Common.OrbitalModel import Orbital
from chemutils.Common.CanonicalAtomMapSmiles import createCanonicalAtomMapSmiString
from chemutils.Common.Util import clearAtomMapsSmiStr

from chemutils.CombiCDB.OrbitalInteraction import reactionSmilesFromOrbitalPair
from chemutils.Common.Const import REACTION_COMPONENT_DELIM

from chemutils.CombiCDB.Util import log


class PeriPattern(object):
    """Class to capture proposal of reactions for a single pattern."""

    def __init__(
        self,
        smarts,
        srcOrbStr,
        sinkOrbStr,
        breakAtomIds,
        aroFirst=False,
        orderNoMatter=True,
        kekBondMaps=[],
    ):
        """Constructor - Setup everything we need."""
        self.smarts = smarts
        self.srcOrbStr = srcOrbStr
        self.sinkOrbStr = sinkOrbStr
        self.breakAtomIds = breakAtomIds

        self.aroFirst = aroFirst

        self.orderNoMatter = orderNoMatter
        self.kekBondMaps = kekBondMaps

        ## Pre-process to get the subsearch object and references
        ## to the important pattern atoms
        self.ss = OESubSearch(self.smarts)
        self.breakAtomsDict = {}
        for atom in self.ss.GetPattern().GetAtoms():
            if atom.GetMapIdx() in self.breakAtomIds:
                self.breakAtomsDict[atom.GetMapIdx()] = atom
        self.breakAtomIds.sort()
        self.breakAtoms = [self.breakAtomsDict[i] for i in self.breakAtomIds]

        ## Statistics
        self.numIniConstrFailed = 0
        self.numOtherConstrFailed = 0
        self.seenLabeledSmi = set([])
        self.seenProductSet = set([])
        self.numDuplicateLabels = 0
        self.numDuplicateProds = 0

    def proposeOrbPairs(self, closureList):
        """Actually do the matching, and check against the closures.!

        ## First a simple check.  If the number of possible closure atoms is
        less than the number needed, can return right away.

        Then for each atom in break list:
        - setup a constraint
        - find matches
        - check rest of constraints are satisfied
        - turn into a orb pair and yield

        ## Finally, do a check have we seen this exact reaction before?
        ## By exact, it means we have used all the same atoms, and the products are
        ## the same.
        """
        self.numIniConstrFailed = 0
        self.numOtherConstrFailed = 0

        ## Also include a way to return only unique labeling
        self.seenLabeledSmi = set([])
        self.numDuplicateLabels = 0
        self.seenProductSet = set([])
        self.numDuplicateProds = 0

        ## Initial check:
        completePossBreakSet = set([])
        for atomObj in closureList:
            completePossBreakSet.update(atomObj.symmetricAtoms)
        if len(self.breakAtomIds) > len(completePossBreakSet):
            return
        ## If we pass, then need to go through and do the actual proposal
        mol = closureList[0].mol
        if self.aroFirst:
            OEAssignAromaticFlags(mol)
            deKekulizeMol(mol)
            # log.info('After dekek, smi:  %s' % createCanonicalAtomMapSmiString(mol))

        ## Then for each atom in the closur eList, try setting as a constraint for atom id 1.
        for atomObj in closureList:
            self.ss.ClearConstraints()
            constraint = OEMatchPairAtom(self.breakAtoms[0], atomObj.atom)
            constrOK = self.ss.AddConstraint(constraint)
            if not constrOK:
                log.info("Problem loading a constraint: %s" % pformat(atomObj))
                self.numIniConstrFailed += 1
                continue

            for i, match in enumerate(self.ss.Match(mol, self.orderNoMatter)):
                # log.info('Have a match: %d!' % i)
                [a.SetMapIdx(0) for a in mol.GetAtoms()]
                if self.labelMolFromMatch(match, completePossBreakSet):
                    if self.checkUniqueLabel(mol):
                        orbPair = self.orbPairFromLabeledMol(mol)
                        if self.checkUniqueProducts(orbPair):
                            yield orbPair
                        else:
                            self.numDuplicateProds += 1
                    else:
                        self.numDuplicateLabels += 1
                else:
                    self.numOtherConstrFailed += 1

    def orbPairFromLabeledMol(self, mol):
        """Given a properly labeled mol, makes copies, ensure kekulization, return the orbPair"""
        nMol = OEGraphMol(mol)
        self.ensureKekulization(nMol)
        srcOrb = Orbital.fromLabeledMolAndInfoStr(nMol, self.srcOrbStr)
        sinkOrb = Orbital.fromLabeledMolAndInfoStr(nMol, self.sinkOrbStr)
        return (srcOrb, sinkOrb)

    def ensureKekulization(self, mol):
        """Ensure proper kekulization if specified"""
        if self.kekBondMaps == []:
            return

        ## We have some bonds that need to be set to be kekulized in a particular way
        aromMap = {}
        for mapTuple in self.kekBondMaps:
            bond = findBondByAtomMapIdx(mol, mapTuple)
            if bond.IsAromatic():
                aromMap[bond.GetIdx()] = 2

        if len(aromMap) > 0:
            defAroMap = createAroBondTypeMap(mol)
            defAroMap.update(aromMap)
            kekulizeMol(mol, defAroMap)

    def checkUniqueProducts(self, orbPair):
        """Check that the products have not yet been seen!"""
        srcOrb, sinkOrb = orbPair
        rxnSmiles = reactionSmilesFromOrbitalPair(srcOrb, sinkOrb)
        prodsmi = clearAtomMapsSmiStr(rxnSmiles.split(REACTION_COMPONENT_DELIM)[-1])
        if prodsmi not in self.seenProductSet:
            # log.info('Proposing prodsmi : %s' % prodsmi)
            self.seenProductSet.add(prodsmi)
            return True
        return False

    def checkUniqueLabel(self, mol):
        """Given a labeled mol, check that we haven't seen it before"""
        smi = createCanonicalAtomMapSmiString(mol)
        if smi in self.seenLabeledSmi:
            self.numDuplicateLabels += 1
            return False

        self.seenLabeledSmi.add(smi)
        return True

    def labelMolFromMatch(self, match, possBreakSet):
        """Given a match and possBreakSet, check constraints, label properly then return

        Assumes that the mapidxs for the mol have been reset.
        Returns whether constraints were met

        NOTE: does in-place change of input mol
        """
        retValue = True
        for ma in match.GetAtoms():
            ## Check that constraints are satified.
            if ma.pattern in self.breakAtoms and ma.target not in possBreakSet:
                retValue = False
                # log.info('Failed on ma.pattern.GetMapIdx() : %d' % ma.pattern.GetMapIdx())
            ma.target.SetMapIdx(ma.pattern.GetMapIdx())

        return retValue


defaultPeriPatterns = [  ## Electrocylic:
    PeriPattern(
        "[C:1]=[C:2][C:3]=[C:4][C:5]=[C:6]",
        "1 pi 2 2",
        "6 pi 5 2; 4 pi 3 2; 2 p None 0",
        [1, 6],
        aroFirst=True,
    ),  ## EC 6pi close
    PeriPattern(
        "[C:1]1[C:2]=[C:3][C:4]=[C:5][C:6]1",
        "3 pi 2 2",
        "4 pi 5 2; 6 sigma 1 2; 2 p None 0",
        [1, 6],
        aroFirst=True,
    ),  ## EC 6pi open
    PeriPattern(
        "[C;X4:1]1[#6:2]=,:[#6:3][C;X4:4]1",
        "3 pi 2 2",
        "4 sigma 1 2; 2 p None 0",
        [1, 4],
        aroFirst=True,
        kekBondMaps=[(2, 3)],
    ),  ## EC 4pi open
    PeriPattern(
        "[C:1]=[C:2][C:3]=[C:4]",
        "2 pi 1 2",
        "3 pi 4 2; 1 p None 0",
        [1, 4],
        aroFirst=True,
    ),  ## EC 4pi close
    PeriPattern(
        "[C;+1:1]1[C:2]=[C:3][C:4][C:5]1",
        "3 pi 2 2",
        "4 sigma 5 2; 1 p None 0",
        [4, 5],
        aroFirst=True,
    ),  ## EC 4pi open cation
    PeriPattern(
        "[C:1]1=[C:2][C:3][C:4]1", "2 pi 1 2", "3 sigma 4 2; 1 p None 0", [3, 4]
    ),
    PeriPattern(
        "[C;!$(C1=C[C;+1]C=C1):1]=[C:2][C;+1:3][C:4]=[C:5]",
        "1 pi 2 2",
        "5 pi 4 2; 3 p None 0",
        [1, 5],
        aroFirst=True,
    ),  ## EC 4pi close cation
    ## Cycloaddition:
    PeriPattern(
        "[C,N,O:1]=[C,N:2][C,N:3]=[C:4].[C,O,N:10]#,=[C,N,S:11]",
        "10 pi 11 2",
        "1 pi 2 2; 3 pi 4 2; 11 p None 0",
        [1, 4, 10, 11],
        aroFirst=True,
    ),  ## CA 4+2
    PeriPattern(
        "[C:1]=[C,O,N:2].[C,N;+0:3]=[P,C&$(C(=O)):4]",
        "1 pi 2 2",
        "3 pi 4 2; 2 p None 0",
        [1, 2, 3, 4],
        aroFirst=True,
    ),  # CA 2+2
    PeriPattern(
        "[C;+0:10]#,=[C,O;+0:11].[O,N;-1:1]=,-[O,N;+1:2]=,#[C,O,N;+0:3]",
        "1 sp3 None 2",
        "10 pi 11 2; 3 pi 2 2",
        [1, 3, 10, 11],
        aroFirst=True,
        orderNoMatter=False,
    ),  ## CA 3+2 dipolar
    PeriPattern(
        "[O;+0:1]=[Mn,Os:2]=[O;+0:3].[C:4]=,#[C:5]",
        "1 pi 2 2",
        "4 pi 5 2; 3 pi 2 2",
        [1, 3, 4, 5],
        aroFirst=True,
    ),  ## CA 3+2 Metal (non-dipolar)
    PeriPattern(
        "[C;v2;+0:1].[C:2]=[C:3]",
        "1 sp2 None 2",
        "2 pi 3 2; 1 p None 0",
        [1, 2, 3],
        aroFirst=True,
    ),  ## CA 2+1 Cheletropic
    PeriPattern(
        "[S;$(S(=O)=O):10].[C:1]=[C:2][C:3]=[C:4]",
        "10 sp3 None 2",
        "1 pi 2 2; 3 pi 4 2; 10 p None 0",
        [10, 1, 4],
        aroFirst=True,
    ),  ## CA 4+1 Cheletropic
    PeriPattern(
        "[C:1]1=[C:2][C:3][C&$(C=O),S&$(S(=O)=O):4][C:5]1",
        "5 sigma 4 2",
        "1 pi 2 2; 3 sigma 4 2",
        [3, 4, 5],
        aroFirst=True,
    ),  ## Retro-CA 4+1
    PeriPattern(
        "[C,N:1]1[C,N:2]=[C,N:3][C:4][C,N,S:11]=,-[C,O,N:10]1",
        "2 pi 3 2",
        "1 sigma 10 2; 11 sigma 4 2; 3 p None 0",
        [1, 10, 11, 4],
        aroFirst=True,
    ),  ## Retro-CA 4+2
    PeriPattern(
        "[P,N:1]1[O:2][C:3][C:4]1",
        "4 sigma 1 2",
        "3 sigma 2 2; 1 p None 0",
        [1, 2, 3, 4],
        aroFirst=True,
    ),  ## Retro-CA 2+2 (very restricted.)
    PeriPattern(
        "[O:1]1[O:2][O:3][C:4][C:5]1",
        "3 sp3 None 2",
        "4 sigma 5 2; 1 sigma 2 2",
        [1, 2, 4, 5],
        aroFirst=True,
    ),  ## Retro-CA 3+2 (very restricted)
    ## Sigmatropic Rearrangements
    PeriPattern(
        "[#6,O&!$([#6,O]1=,:[#6,#7][C,O]CC=C1):1]=,:[#6,#7:2][C,O:3][C:4][C:5]=[C:6]",
        "1 pi 2 2",
        "6 pi 5 2; 4 sigma 3 2; 2 p None 0",
        [1, 6, 4, 3],
        aroFirst=True,
        kekBondMaps=[(1, 2)],
    ),  # SigmaR 3,3
    PeriPattern(
        "[-1;!$(*1[+1,O]C[#6]=,:[#6]*1):1][+1,O:2][C,N:3][#6:4]=,:[#6:5]",
        "1 sp3 None 2",
        "5 pi 4 2; 3 sigma 2 2",
        [5, 1, 3, 2],
        aroFirst=True,
        kekBondMaps=[(4, 5)],
    ),  ## SigmaR 2,3
    PeriPattern(
        "[#1,C:1][C:2]1[C:3]=[C:4][C:5]=[C:6]1",
        "1 sigma 2 2",
        "3 pi 4 2; 5 pi 6 2; 2 p None 0",
        [1, 2, 3],
        orderNoMatter=False,
        aroFirst=True,
    ),  ## SigmaR 1,5
    ## Ene Reactions
    PeriPattern(
        "[#1,Mg:1][*;+0:2][*;+0:3]=[*;+0:4].[*:5]#,=[*;+0:6]",
        "1 sigma 2 2",
        "6 pi 5 2; 4 pi 3 2; 2 p None 0",
        [1, 2, 6, 5, 4],
        aroFirst=True,
        orderNoMatter=False,
    ),  ## basic Ene
    PeriPattern(
        "[#1:1][*;+0:2][*;+0:3][*;+0:4][*;+0:5]=[*;+0:6]",
        "1 sigma 2 2",
        "6 pi 5 2; 4 sigma 3 2; 2 p None 0",
        [1, 6, 2, 4, 3],
        aroFirst=True,
    ),  ## Retro-ene
    PeriPattern(
        "[O,C;-1:1][!C;+1:2][*;+0:3][C;+0:4][#1:5]",
        "1 sp3 None 2",
        "5 sigma 4 2; 3 sigma 2 2",
        [1, 5, 4, 3, 2],
        aroFirst=True,
    ),  ## Retro-hetero-ene
    ## Concerted Peroxyacid reactions:
    PeriPattern(
        "[O:1]([H:2])[O:3][C:4]=[O:5].[C:6]=[C:7]",
        "5 pi 4 2",
        "2 sigma 1 2; 6 pi 7 2; 1 sigma 3 2; 4 sp3 None 0",
        [1, 3, 5, 6, 7],
        aroFirst=True,
    ),
    ## Halogenation
    PeriPattern(
        "[Cl,Br,I:1][Cl,Br,I:2].[#6:3]=,#[#6:4]",
        "1 sp3 None 2",
        "3 pi 4 2; 1 sigma 2 2",
        [1, 2, 3, 4],
        aroFirst=True,
    ),  ## Simple halonium ion formation
    ## Boro-hydride reactions
    PeriPattern(
        "[B;+0:2][H:1].[C:3]=,#[C:4]",
        "1 sigma 2 2",
        "3 pi 4 2; 2 sp3 None 2",
        [1, 2, 3, 4],
        aroFirst=True,
        orderNoMatter=False,
    ),  ## Boro-hydride
]
