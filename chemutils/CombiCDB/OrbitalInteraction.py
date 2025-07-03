"""Simple interaction models for fundamental reactions of
orbital overlap interactions.

Key functions for external callers:

    moveOrbitalElectrons:
        Fundamental reaction unit, moves the electrons from a given filled orbital
        to an unfilled orbital, rearranging the electron configuration (bonds, formal charges)
        accordingly.
    undoMoveOrbitalElectrons:
        After a call to moveOrbitalElectrons, this should undo the changes made
        in electron configuration, reverting the composite molecule to its prior state.
    perceiveInvertedOrbitals:
        After a call to moveOrbitalElectrons, perceive the orbitals in the
        new composite molecule state that when applied would effectively invert
        the transformation and revert the molecule to its prior state.
    reactionSmilesFromOrbitalPair:
        Apply a filled-unfilled orbital pair and produce the reactions SMILES
        string that represents the change.

    ElementaryStepData:
        Class / struct to manage core information describing an orbital interaction

"""
# from sets import Set;
from openeye.oechem import (
    OEGraphMol,
    OEAtomStereo_Tetrahedral,
    OEAtomStereo_RightHanded,
    OEAtomStereo_LeftHanded,
    OEBondStereo_Cis,
    OEBondStereo_CisTrans,
)

from openeye.oechem import OEParseSmiles
from openeye.oechem import OEClearAromaticFlags
from chemutils.Common.Const import REACTION_COMPONENT_DELIM, SMILES_MOL_DELIM
from chemutils.Common.Util import (
    createStandardSmiString,
)
from chemutils.Common.Util import createStdAtomMapSmiString
from chemutils.Common.MolExt import clearAtomStereo
from chemutils.Common.OrbitalModel import (
    Orbital,
    atomHybridization,
)
from chemutils.CombiCDB.MechanismModel import (
    moveAtomElectrons,
    moveBondElectrons,
)
from chemutils.CombiCDB.MechanismModel import clearMechanismLabels
from chemutils.CombiCDB.ReactionFormatter import aReactionCanonizer
from chemutils.CombiCDB.Util import log
from chemutils.CombiCDB.Const import (
    VALUE_DELIM,
    EXT_TYPE_UNDEFINED,
)

# Some constant to denote sp3 hybridization
SP3_HYB = 3

FOUR_BONDS = 4

# Constants for Bond order
SIGMA_BOND = 1
DOUBLE_BOND = 2
TRIPLE_BOND = 3

# Constant for Types
SIGMA_BOND_TYPE = "sigma"
PI_BOND_TYPE = "pi"
P_ORB_TYPE = "p"


def moveOrbitalElectrons(
    sourceOrb, targetOrb, electronCount=None, arrowObjList=None, halfBondIndexes=None
):
    """Wrapper around the recursive moveOrbitalElectrons_recursive which handles aromaticity."""
    mol = targetOrb.mol
    OEClearAromaticFlags(mol)

    log.debug("Enter moveOrbitalElectrons (non-recursive)")
    if sourceOrb is not None:
        log.debug("Src Orb: %s" % str(sourceOrb.toLabeledSmiAndInfoStr()))
    if targetOrb is not None:
        log.debug("Sink Orb: %s" % str(targetOrb.toLabeledSmiAndInfoStr()))

    bond = moveOrbitalElectrons_recursive(
        sourceOrb, targetOrb, electronCount, arrowObjList, halfBondIndexes
    )

    # OEAssignAromaticFlags(mol)

    return bond


def moveOrbitalElectrons_recursive(
    sourceOrb, targetOrb, electronCount=None, arrowObjList=None, halfBondIndexes=None
):
    """Most general half-reaction, movement of electrons from
    one molecular orbital to another.  Function should figure out
    all of the substeps necessary to properly represent the change,
    including consideration for whether the orbitals are atomic (non-bond) or molecular (bond).

    Assume this is called in a sensible manner however.
    Cannot use any arbitary orbitals as the arguments.
    The source orbital (sourceOrb) must contain electrons to move, no empty orbitals allowed.
    The target orbital (targetOrb) must have space for electrons to move in.
    In general this means it is either an empty orbital or a bond orbital which can
    displace the existing bond electrons to the neighbor atom.  No lone pairs allowed.

    Additionally support a special case for bond dissociation reactions.
    The source and target orbital should be based on the same bond, with polarity
    of electron movement based on the target orbital polarity.
    Alternatively, specifying "None" for the sourceOrb as a sentinel object
    will be interpreted as a bond dissociation reaction for the targetOrb bond.

    Returns the bond just created or added to by the move.

    If the arrowObjList is provided, then items will be added for the
    collection of ElectronArrow objects representing the respective "arrow-pushing"
    mechanism diagram.  These will only make sense if the atoms involved have
    non-zero atom mapping index labels.
    Make this optional as excessive ElectronArrow object instantiation can be
    relatively expensive, and is often unused.

    To handle aromatic compounds correctly, we clear the flags upon entry and then try to reinstantiate upon exit.
    """
    log.debug("ENTER moveOrbitalElectrons")
    if sourceOrb is not None:
        log.debug("sourceOrb : %s" % str(sourceOrb.toLabeledSmiAndInfoStr()))
    log.debug("targetOrb : %s" % str(targetOrb.toLabeledSmiAndInfoStr()))

    # Percieve free radical reactions
    if electronCount is None:
        # Should percieve here.
        isFreeRadical = targetOrb.isPerceivedFreeRadical() or (
            sourceOrb is not None and sourceOrb.isPerceivedFreeRadical()
        )
        if isFreeRadical:
            electronCount = 1
        else:
            electronCount = 2
    else:
        isFreeRadical = electronCount == 1

    if isFreeRadical and halfBondIndexes is None:
        halfBondIndexes = set()

    # The existing or new bond creating by moving the orbital electrons
    bond = None

    # Before any movements, set these flags.
    preMovementSetStereoFlags(sourceOrb, targetOrb)

    if isBondDissociation(sourceOrb, targetOrb):
        # Special case.  Source and target orbitals are based on the same bond.
        #   Already bonded and just want to move the bond electrons to one side.
        #   Base the direction of electron movement on the targetOrb.
        moveBondElectrons(
            targetOrb.getBond(),
            targetOrb.neighbor,
            electronCount=electronCount,
            arrowObjList=arrowObjList,
            halfBondIndexes=halfBondIndexes,
        )
        if isFreeRadical:
            moveBondElectrons(
                targetOrb.getBond(),
                targetOrb.atom,
                electronCount=electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )
    else:
        if not sourceOrb.isBondOrbital():
            # Atom orbital (lone pair or radical), simple electron movement
            bond = moveAtomElectrons(
                sourceOrb.atom,
                targetOrb.atom,
                electronCount=electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )
        else:
            # Direct the bond electrons the the local atom to represent the polarity of the bond electrons
            moveBondElectrons(
                sourceOrb.getBond(),
                targetOrb.atom,
                pivotAtom=sourceOrb.atom,
                electronCount=electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )

            # log.debug('arrowObjList is %s' % str([str(arr) for arr in arrowObjList]))
            # Hard to retrieve generated bond directly from function, have to get it indirectly
            bond = sourceOrb.atom.GetBond(targetOrb.atom)

            if isFreeRadical:
                moveBondElectrons(
                    sourceOrb.getBond(),
                    sourceOrb.neighbor,
                    electronCount=electronCount,
                    arrowObjList=arrowObjList,
                    halfBondIndexes=halfBondIndexes,
                )
        if targetOrb.isBondOrbital():
            # Displace the electrons in the target atom
            moveBondElectrons(
                targetOrb.getBond(),
                targetOrb.neighbor,
                electronCount=electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )
            if isFreeRadical:
                # An electron from this bond will contribute to half-bond made with the source atom
                moveBondElectrons(
                    targetOrb.getBond(),
                    sourceOrb.atom,
                    pivotAtom=targetOrb.atom,
                    electronCount=electronCount,
                    arrowObjList=arrowObjList,
                    halfBondIndexes=halfBondIndexes,
                )
        else:
            # For a non-bond target orbital, and there was a free radical, then one electron from the target will add to existing
            # half-bond
            if isFreeRadical:
                moveAtomElectrons(
                    targetOrb.atom,
                    sourceOrb.atom,
                    electronCount=electronCount,
                    arrowObjList=arrowObjList,
                    halfBondIndexes=halfBondIndexes,
                )
    #
    # This should cover all of the possible stereo corrections.  (And all necessary information
    # is encoded in the faceLists on the orbitals )
    postMovementDoStereoCorrection(sourceOrb, targetOrb, bond)

    if sourceOrb is not None and sourceOrb.extOrbital is not None:
        if sourceOrb.isBondOrbital():
            # After the original orbital electron move, the source neighbor atom should be left
            #   with an empty orbital that the extended chain can move into
            pivotOrb = Orbital(
                sourceOrb.neighbor, "p", None, 0, sourceOrb.mol, sourceOrb.nbrFaceList
            )
            pivotBond = moveOrbitalElectrons_recursive(
                sourceOrb.extOrbital,
                pivotOrb,
                electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )
        else:
            raise NotImplementedError("Only bond orbitals support extended chaining.")

    if targetOrb.extOrbital is not None:
        if targetOrb.isBondOrbital():
            # After the original orbital electron move, the target neighbor atom should be left
            #   with a lone pair that can move into the extended chain
            pivotAtom = targetOrb.neighbor

            pivotOrb = Orbital(
                pivotAtom,
                None,
                None,
                electronCount,
                sourceOrb.mol,
                targetOrb.nbrFaceList,
            )
            pivotBond = moveOrbitalElectrons_recursive(
                pivotOrb,
                targetOrb.extOrbital,
                electronCount,
                arrowObjList=arrowObjList,
                halfBondIndexes=halfBondIndexes,
            )
        else:
            raise NotImplementedError("Only bond orbitals support extended chaining.")

    # Set any flags necessary to create the inverse..,
    # postRecursionSetInvStereoFlags(sourceOrb, targetOrb);

    return bond


def postRecursionSetInvStereoFlags(sourceOrb, targetOrb):
    """Function to set any stereo flags necessary to create the inverses."""
    # Here if we have a pi-bond (source or sink) that just took part in a syn addition, then simply flip the
    # atom face list so the inverse is able to correctly figure out what is going on..,
    log.debug("ENTER postRecursionSetInvStereoFlags")
    if sourceOrb is not None and sourceOrb.isSynAddition():
        log.debug("SRC SYN ADDITION CORRECTION")
        sourceOrb.flipFacesForSynAddition()
    if targetOrb.isSynAddition():
        log.debug("TARGET SYN ADDITION CORRECTION")
        targetOrb.flipFacesForSynAddition()


def preMovementSetStereoFlags(sourceOrb, targetOrb):
    """Function to do pre-processing of orbitals to help with stereo correction"""
    # Propagate face information if isPericyclic.
    if isPericyclic(sourceOrb, targetOrb):
        doPericyclicFaceCorrection(sourceOrb, targetOrb)

    # Flags about what stereo info to track.
    targetAtomHasStereo = targetOrb.atom.HasStereoSpecified()
    targetNbrHasStereo = (
        targetOrb.neighbor is not None and targetOrb.neighbor.HasStereoSpecified()
    )
    sourceAtomHasStereo = sourceOrb is not None and sourceOrb.atom.HasStereoSpecified()
    sourceNbrHasStereo = (
        sourceOrb is not None
        and sourceOrb.neighbor is not None
        and sourceOrb.neighbor.HasStereoSpecified()
    )

    # have to track some info about the target atom before any bonds are made to it..,
    # Track the hybridization of the targetAtom before a bond has been added;
    if targetAtomHasStereo:
        targetAtomIsSp3Hyb = atomHybridization(targetOrb.atom) == SP3_HYB
        targetAtomAtomFaceList = targetOrb.atomFaceList
        if (
            targetAtomAtomFaceList is None
            and targetOrb.isBondOrbital()
            and targetAtomIsSp3Hyb
        ):
            targetAtomAtomFaceList = determineAtomFaceList(
                targetOrb.atom, targetOrb.neighbor
            )

    # Also, need to track whether a pibond break is syn or anti addition.
    # simply add a flag to the orbital if this represents a syn additions
    if sourceOrb is not None:
        sourceOrb.setSynAdditionFlags()
        log.debug("sourceOrb.isSynAddition() : %s" % str(sourceOrb.isSynAddition()))
    targetOrb.setSynAdditionFlags()
    log.debug("targetOrb.isSynAddition() : %s" % str(targetOrb.isSynAddition()))

    # Bond Dissociation
    if isBondDissociation(sourceOrb, targetOrb):
        # Need to save info about the targetAtom to ensure inverse returns stereo
        # pi-bonds will figure out inverse faces
        if targetAtomHasStereo and targetAtomIsSp3Hyb:
            targetOrb.atomFaceList = targetAtomAtomFaceList
        clearAtomStereo(targetOrb.atom)
        if (
            targetOrb.nbrFaceList is None
            and targetNbrHasStereo
            and atomHybridization(targetOrb.neighbor) == SP3_HYB
        ):
            targetOrb.nbrFaceList = determineAtomFaceList(
                targetOrb.neighbor, targetOrb.atom
            )
        clearAtomStereo(targetOrb.neighbor)
        return

    # Here we are breaking a bond, (source.neighbor to source.atom) so need to leave a bread crumb for backwards movement
    if sourceOrb.type == SIGMA_BOND_TYPE:
        if (
            sourceOrb.atomFaceList is None
            and sourceAtomHasStereo
            and atomHybridization(sourceOrb.atom) == SP3_HYB
        ):
            sourceOrb.atomFaceList = determineAtomFaceList(
                sourceOrb.atom, sourceOrb.neighbor
            )
        clearAtomStereo(sourceOrb.atom)
        if (
            sourceOrb.nbrFaceList is None
            and sourceNbrHasStereo
            and atomHybridization(sourceOrb.neighbor) == SP3_HYB
        ):
            sourceOrb.nbrFaceList = determineAtomFaceList(
                sourceOrb.neighbor, sourceOrb.atom
            )
        clearAtomStereo(sourceOrb.neighbor)

    if targetOrb.isBondOrbital():
        # Breaking bond, so save breadcrumb
        if (
            targetOrb.atomFaceList is None
            and targetAtomHasStereo
            and targetAtomIsSp3Hyb
        ):
            targetOrb.atomFaceList = targetAtomAtomFaceList
        clearAtomStereo(targetOrb.atom)
        if (
            targetOrb.nbrFaceList is None
            and targetNbrHasStereo
            and atomHybridization(targetOrb.neighbor) == SP3_HYB
        ):
            targetOrb.nbrFaceList = determineAtomFaceList(
                targetOrb.neighbor, targetOrb.atom
            )
        clearAtomStereo(targetOrb.neighbor)


def postMovementDoStereoCorrection(sourceOrb, targetOrb, bond):
    """Actually correct the stereochemistry of the reaction using face information in the orbitals

    What the correction is will depend on the bond formed.
    """
    anyHasFace = (
        targetOrb.atomFaceList is not None
        or targetOrb.nbrFaceList is not None
        or (
            sourceOrb is not None
            and (
                sourceOrb.atomFaceList is not None or sourceOrb.nbrFaceList is not None
            )
        )
    )
    if not anyHasFace or bond is None:
        return

    log.debug("ENTER doStereoCorrection;")
    if sourceOrb.atomFaceList is not None:
        log.debug(
            "Source.atomFaceList is %s"
            % str([atm.GetMapIdx() for atm in sourceOrb.atomFaceList])
        )
    if sourceOrb.nbrFaceList is not None:
        log.debug(
            "Source.nbrFaceList is %s"
            % str([atm.GetMapIdx() for atm in sourceOrb.nbrFaceList])
        )
    if targetOrb.atomFaceList is not None:
        log.debug(
            "Target.atomFaceList is %s"
            % str([atm.GetMapIdx() for atm in targetOrb.atomFaceList])
        )
    if targetOrb.nbrFaceList is not None:
        log.debug(
            "targetOrb.nbrFaceList is %s"
            % str([atm.GetMapIdx() for atm in targetOrb.nbrFaceList])
        )

    if bond.GetOrder() == SIGMA_BOND:
        if targetOrb.type == SIGMA_BOND_TYPE:
            if targetOrb.atomFaceList is not None:  # and targetOrb.atom.IsChiral():
                log.debug("GOING TO CORRECT SN2 Chemistry")
                # Sn2:Invert previous (facelist is defined to be 'RightHanded' from previous bond..., so flip)
                targetOrb.atomFaceList.insert(0, sourceOrb.atom)
                targetOrb.atom.SetStereo(
                    targetOrb.atomFaceList,
                    OEAtomStereo_Tetrahedral,
                    OEAtomStereo_LeftHanded,
                )
                # Here we clear the faceList because this could propagate the wrong face info.
                targetOrb.atomFaceList = None
        elif targetOrb.type == PI_BOND_TYPE or targetOrb.type == P_ORB_TYPE:
            if (
                targetOrb.atomFaceList is not None
                and targetOrb.atom.GetDegree() == FOUR_BONDS
            ):
                # Then backwards E2 with initial pi bond break, (or extension after pi bond break)
                # set the new bond to be RightHanded to the face
                log.debug("GOING TO CORRECT FOR BACK E2 LIKE Movement")
                targetOrb.atomFaceList.insert(0, sourceOrb.atom)
                targetOrb.atom.SetStereo(
                    targetOrb.atomFaceList,
                    OEAtomStereo_Tetrahedral,
                    OEAtomStereo_RightHanded,
                )
                targetOrb.atomFaceList.pop(0)

        # No matter what the target is, we have formed a new bond from source -> target, correct source stereo here if
        # face info is given and the source atom now has 4 bonds.
        if (
            sourceOrb.atomFaceList is not None
            and sourceOrb.atom.GetDegree() == FOUR_BONDS
        ):
            # Backwards E2 second step (from initially displaced electrons on beta carbon)
            # From the breadcrumb, know this should be RightHanded to latest bond
            log.debug("CORRECTING SOURCE.ATOM FACE INFO ON SIGMA BOND FORM")
            sourceOrb.atomFaceList.insert(0, targetOrb.atom)
            sourceOrb.atom.SetStereo(
                sourceOrb.atomFaceList,
                OEAtomStereo_Tetrahedral,
                OEAtomStereo_RightHanded,
            )
            sourceOrb.atomFaceList.pop(0)

    elif bond.GetOrder() == DOUBLE_BOND:
        # If we have faces specified for both the targetAtom and sourceAtom, should be able to figure out C/S
        if targetOrb.atomFaceList is not None and sourceOrb.atomFaceList is not None:
            # First arrange both of the face lists such that the bond atom is last
            bondAtomIdx = -1
            for idx, atm in enumerate(targetOrb.atomFaceList):
                if atm == sourceOrb.atom:
                    bondAtomIdx = idx
                    break
            targetAtomFaceList = targetOrb.atomFaceList[bondAtomIdx + 1 :]
            targetAtomFaceList.extend(targetOrb.atomFaceList[: bondAtomIdx + 1])

            bondAtomIdx = -1
            for idx, atm in enumerate(sourceOrb.atomFaceList):
                if atm == targetOrb.atom:
                    bondAtomIdx = idx
                    break
            sourceAtomFaceList = sourceOrb.atomFaceList[bondAtomIdx + 1 :]
            sourceAtomFaceList.extend(sourceOrb.atomFaceList[: bondAtomIdx + 1])

            # Now the first in each of these lists will be cis to each other (opposite faces)
            # NOTE: could define this however we want, but this will follow from the E2 (allowing the
            # same data structure to recover chirality in inverse)
            log.debug(
                "SETTING C,S of double bond %s"
                % str(
                    [
                        sourceAtomFaceList[0].GetMapIdx(),
                        targetAtomFaceList[0].GetMapIdx(),
                    ]
                )
            )
            bond.SetStereo(
                [sourceAtomFaceList[0], targetAtomFaceList[0]],
                OEBondStereo_CisTrans,
                OEBondStereo_Cis,
            )

    # Then no mattter what the new bond is, could have some face info from the displacement of electrons from a sourceBond
    # Could have broken a double bond here, and need to specify the tetrahedral stereo
    shouldCorrectSrcNbr = (
        sourceOrb.nbrFaceList is not None
        and sourceOrb.neighbor is not None
        and sourceOrb.atom.GetBond(sourceOrb.neighbor) is not None
        and sourceOrb.neighbor.GetDegree() == FOUR_BONDS
    )

    if shouldCorrectSrcNbr:
        log.debug("CORRECTING SOURCE.NBR FACE INFO ON SIGMA BOND FORM")
        sourceOrb.nbrFaceList.insert(0, sourceOrb.atom)
        sourceOrb.neighbor.SetStereo(
            sourceOrb.nbrFaceList, OEAtomStereo_Tetrahedral, OEAtomStereo_RightHanded
        )
        sourceOrb.nbrFaceList.pop(0)


def determineAtomFaceList(atom, removedAtom):
    """Simple function to return a list of atoms to denote reactive face on atom when
    a bond btwn atom and removedAtom is broken

    Note: only makes sense when atom is sp3 hybridized and OEAtomStereo_Tetrahedral is specified.

    The face is defined as a list of atoms such that when looking from the removed bond,
    the list is in right-handed (clockwise) order.
    """
    if not atom.HasStereoSpecified() or atom.GetBond(removedAtom) is None:
        return None
    nbrList = [atm for atm in atom.GetAtoms()]
    removedAtomIdx = -1
    for i, atm in enumerate(nbrList):
        if atm.GetIdx() == removedAtom.GetIdx():
            removedAtomIdx = i
            break
    nbrList[removedAtomIdx : removedAtomIdx + 1] = []
    nbrList.insert(0, removedAtom)
    stereo = atom.GetStereo(nbrList, OEAtomStereo_Tetrahedral)
    nbrList.pop(0)

    if stereo != OEAtomStereo_RightHanded:
        nbrList.reverse()
    return nbrList


def undoMoveOrbitalElectrons(
    sourceOrb, targetOrb, electronCount=None, arrowObjList=None
):
    """Assume the moveOrbitalElectrons function was just called for sourceOrb and targetOrb.
    Revert the molecule state to where it was before then based on the info in the orbitals.

    Make a good guess about the electronCount if None
    """
    if electronCount is None:
        if (
            sourceOrb is not None and sourceOrb.isPerceivedFreeRadical()
        ) or targetOrb.isPerceivedFreeRadical():
            electronCount = 1
        else:
            electronCount = 2

    (invSourceOrb, invTargetOrb) = perceiveInvertedOrbitals(
        sourceOrb, targetOrb, electronCount
    )
    log.debug("perceived invsrc : %s" % str(invSourceOrb.toLabeledSmiAndInfoStr()))
    log.debug("perceived invsink : %s" % str(invTargetOrb.toLabeledSmiAndInfoStr()))
    invBond = moveOrbitalElectrons(
        invSourceOrb, invTargetOrb, electronCount, arrowObjList=arrowObjList
    )
    if sourceOrb is not None:
        sourceOrb.removeFaceBreadCrumbs()
    if targetOrb is not None:
        targetOrb.removeFaceBreadCrumbs()

    # Finally, apply stereo corrections on the

    return invBond


def perceiveInvertedOrbitals(sourceOrb, targetOrb, electronCount=None):
    """Assume the moveOrbitalElectrons function was just called for sourceOrb and targetOrb.
    Produce a new pair of "inverted" source and sink orbitals such that, when these
    are applied to each other, the underlying composite molecule will revert back to
    the state it was in before the original source and target orbitals interacted.
    """
    log.debug("Enter perceiveInvertedOrbitals")
    if sourceOrb is not None:
        log.debug("sourceOrb: %s" % str(sourceOrb.toLabeledSmiAndInfoStr()))
    log.debug("targetOrb: %s" % str(targetOrb.toLabeledSmiAndInfoStr()))

    # Need to 'perceive' the electron count if it is None
    if electronCount is None:
        # Should percieve here.
        isFreeRadical = targetOrb.isPerceivedFreeRadical() or (
            sourceOrb is not None and sourceOrb.isPerceivedFreeRadical()
        )
        if isFreeRadical:
            electronCount = 1
        else:
            electronCount = 2

    mol = targetOrb.mol

    invSourceOrb = None
    invTargetOrb = None

    if isBondDissociation(sourceOrb, targetOrb):
        # Special case.  Was already bonded together and just wanted to move the
        #   electrons in the bond from one side to the other (bond dissociation)
        invSourceOrb = Orbital(
            targetOrb.neighbor,
            None,
            None,
            targetOrb.electrons,
            mol,
            targetOrb.nbrFaceList,
        )
        invTargetOrb = Orbital(
            targetOrb.atom, None, None, 2 - electronCount, mol, targetOrb.atomFaceList
        )

    else:
        # Unless there was a single bond dissociation, the moveOrbitalElectrons should have
        #   resulted in a bond formation between the atoms of orbitals 1 and 2
        bond = sourceOrb.atom.GetBond(targetOrb.atom)

        # Grab and flip if necessary the face lists
        atomFaceList = targetOrb.atomFaceList
        nbrFaceList = sourceOrb.atomFaceList
        if targetOrb.isSynAddition():
            atomFaceList = Orbital.flipFaceList(atomFaceList)
        if sourceOrb.isSynAddition():
            nbrFaceList = Orbital.flipFaceList(nbrFaceList)

        # General cases, separate source and target orbitals
        # Invert source orbital into inversion target orbital.
        if not sourceOrb.isBondOrbital():
            # Non-bonded electron source, probably a lone pair or free radical
            # Polarity on bond atom choice may seem strange, but
            #   works for consistency with electron movement patterns for anti-bonding orbital targets
            invTargetOrb = Orbital(
                targetOrb.atom,
                None,
                bond,
                sourceOrb.electrons,
                mol,
                atomFaceList,
                nbrFaceList,
            )
        else:
            # Bond electron source, will need to track extended orbitals
            #   similar to an E2 mechanism for the inverted reaction
            invTargetOrb = Orbital(
                targetOrb.atom,
                None,
                bond,
                sourceOrb.electrons,
                mol,
                atomFaceList,
                nbrFaceList,
            )
            # Must extend target orbital
            if sourceOrb.extOrbital is None:
                # No extension known, then should have just been an empty orbital leftover (or a single free radical)
                # Might have to flip the face here tho.
                atomFaceList = sourceOrb.nbrFaceList
                if sourceOrb.isNbrSynAddition():
                    atomFaceList = Orbital.flipFaceList(atomFaceList)
                invTargetOrb.extOrbital = Orbital(
                    sourceOrb.neighbor, "p", None, 2 - electronCount, mol, atomFaceList
                )
            else:
                # Came from an extended orbital, recursively chain through together then
                # Extension will need to consider the chaining orbital with a different
                #   polarity as a target orbital vs. source orbital, so provide a flipped version
                flipSourceOrb = sourceOrb.flippedBondOrbital()
                (extInvSourceOrb, extInvTargetOrb) = perceiveInvertedOrbitals(
                    sourceOrb.extOrbital, flipSourceOrb, electronCount
                )
                invTargetOrb.extOrbital = extInvTargetOrb
                # extInvSourceOrb can be discarded as it represents an unused imaginary orbital
                #   intermediary between the chaining

        # Invert target orbital into inversion source orbital.
        if not targetOrb.isBondOrbital():
            # Non-bond target, must be an empty p orbital

            if not sourceOrb.isBondOrbital():
                # Lone pair / radical > empty p orbital
                # Direct bond formation, inverse is just bond dissociation
                # Sentinel None object to represent inverted reaction is bond dissociation would work too
                invSourceOrb = invTargetOrb.flippedBondOrbital()

                # For a free radical, should set that one of these has only a sinlge electron
                invSourceOrb.electrons = electronCount
            else:
                # Bond > empty p orbital
                # Inversion reactions needs one less orbital, taken from the target side
                invSourceOrb = invTargetOrb.flippedBondOrbital()
                invTargetOrb = invTargetOrb.extOrbital
        else:
            # Attacking an anti-bonding orbital
            if targetOrb.extOrbital is None:
                # Simple case with no extensions, electrons would have ended up in non-bonding orbital
                # Might need to flip this face
                atomFaceList = targetOrb.nbrFaceList
                if targetOrb.isNbrSynAddition():
                    atomFaceList = Orbital.flipFaceList(atomFaceList)
                invSourceOrb = Orbital(
                    targetOrb.neighbor,
                    None,
                    None,
                    targetOrb.electrons,
                    mol,
                    atomFaceList,
                )
            else:
                # Extended orbital to invert by recursively chaining through together
                # Extension will need to consider the chaining orbital with a different
                #   polarity as a target orbital vs. source orbital, so provide a flipped version
                flipTargetOrb = targetOrb.flippedBondOrbital()
                (extInvSourceOrb, extInvTargetOrb) = perceiveInvertedOrbitals(
                    flipTargetOrb, targetOrb.extOrbital, electronCount
                )

                if not extInvTargetOrb.isBondOrbital():
                    # Empty orbital intervening, no point in recording, just skip past it
                    invSourceOrb = extInvSourceOrb
                else:
                    # More general case with orbital chaining with an intervening bond
                    # Flip the extended target bond to make the correct source orb
                    invSourceOrb = extInvTargetOrb.flippedBondOrbital()
                    invSourceOrb.extOrbital = extInvSourceOrb

    ## Make sure to pass along the flag about whether the face info is for backwards compatibitily only
    if invSourceOrb is not None and sourceOrb is not None:
        invSourceOrb.isFaceListOnlyBreadCrumb = sourceOrb.isFaceListOnlyBreadCrumb
    if invTargetOrb is not None and targetOrb is not None:
        invTargetOrb.isFaceListOnlyBreadCrumb = targetOrb.isFaceListOnlyBreadCrumb

    return (invSourceOrb, invTargetOrb)


def isBondDissociation(sourceOrb, targetOrb):
    """Determine whether the given filled and unfilled orbital combination
    represents a bond dissociation reaction.
    """
    bondDissociation = sourceOrb is None and targetOrb.neighbor is not None
    bondDissociation = bondDissociation or (
        sourceOrb.neighbor is not None
        and targetOrb.neighbor is not None
        and set([sourceOrb.atom.GetIdx(), sourceOrb.neighbor.GetIdx()])
        == set([targetOrb.atom.GetIdx(), targetOrb.neighbor.GetIdx()])
    )
    return bondDissociation


def isPericyclic(sourceOrb, targetOrb):
    """Rough identification of pericyclic reactions is when there exists a
    chain of filled or unfilled orbitals and the extreme head and tail of
    these chains meet at a common atom.

    Not perfectly accurate, since 1,3 dipolar reactions don't appear to have
    any overlap, but should not affect the point of identifying valid reactions
    with orbitals that "overlap" common atoms.
    """
    result = False

    if sourceOrb is not None and (
        sourceOrb.extOrbital is not None or targetOrb.extOrbital is not None
    ):
        # Must have at least one extended orbital chain to make a pericyclic reaction

        sourceSet = set(sourceOrb.coveredOrbitalAtomIdx())
        targetSet = set(targetOrb.coveredOrbitalAtomIdx())
        overlapSet = sourceSet.intersection(targetSet)

        # If it is a pericyclic proposal by our definition,
        #   it must overlap at the head / tail by exactly one atom
        if len(overlapSet) == 1:
            # Iterate through any orbital chains until we find the extreme head (source) and tail (target)
            # Slight danger here of an infinite loop (cyclic linked list)
            extremeSourceOrb = sourceOrb
            while extremeSourceOrb.extOrbital is not None:
                extremeSourceOrb = extremeSourceOrb.extOrbital

            extremeTargetOrb = targetOrb
            while extremeTargetOrb.extOrbital is not None:
                extremeTargetOrb = extremeTargetOrb.extOrbital

            # Identify the index of the apparent "linking" atom
            extremeSourceIndex = extremeSourceOrb.atom.GetIdx()
            if extremeSourceOrb.neighbor is not None:
                extremeSourceIndex = extremeSourceOrb.neighbor.GetIdx()

            extremeTargetIndex = extremeTargetOrb.atom.GetIdx()
            if extremeTargetOrb.neighbor is not None:
                extremeTargetIndex = extremeTargetOrb.neighbor.GetIdx()

            result = extremeSourceIndex == extremeTargetIndex

    return result


def doPericyclicFaceCorrection(sourceOrb, targetOrb):
    """The pericyclic reactions have extreme source equal to extreme target.  If face info is specified
    for one, then should copy to the other.
    ASSUMES:  isPericyclic(sourceOrb, targetOrb) == True"""
    extremeSourceOrb = sourceOrb
    while extremeSourceOrb.extOrbital is not None:
        extremeSourceOrb = extremeSourceOrb.extOrbital

    extremeTargetOrb = targetOrb
    while extremeTargetOrb.extOrbital is not None:
        extremeTargetOrb = extremeTargetOrb.extOrbital

    # Ugly if setup, but no easy elegant way
    if extremeSourceOrb.neighbor is not None:
        if extremeTargetOrb.neighbor is not None:
            if (
                extremeSourceOrb.nbrFaceList is not None
                and extremeTargetOrb.nbrFaceList is None
            ):
                extremeTargetOrb.nbrFaceList = extremeSourceOrb.nbrFaceList
            elif (
                extremeSourceOrb.nbrFaceList is None
                and extremeTargetOrb.nbrFaceList is not None
            ):
                extremeSourceOrb.nbrFaceList = extremeTargetOrb.nbrFaceList
        else:
            if (
                extremeSourceOrb.nbrFaceList is not None
                and extremeTargetOrb.atomFaceList is None
            ):
                extremeTargetOrb.atomFaceList = extremeSourceOrb.nbrFaceList
            elif (
                extremeSourceOrb.nbrFaceList is None
                and extremeTargetOrb.atomFaceList is not None
            ):
                extremeSourceOrb.nbrFaceList = extremeTargetOrb.atomFaceList
    else:
        if extremeTargetOrb.neighbor is not None:
            if (
                extremeSourceOrb.atomFaceList is not None
                and extremeTargetOrb.nbrFaceList is None
            ):
                extremeTargetOrb.nbrFaceList = extremeSourceOrb.atomFaceList
            elif (
                extremeSourceOrb.atomFaceList is None
                and extremeTargetOrb.nbrFaceList is not None
            ):
                extremeSourceOrb.atomFaceList = extremeTargetOrb.nbrFaceList
        else:
            if (
                extremeSourceOrb.atomFaceList is not None
                and extremeTargetOrb.atomFaceList is None
            ):
                extremeTargetOrb.atomFaceList = extremeSourceOrb.atomFaceList
            elif (
                extremeSourceOrb.atomFaceList is None
                and extremeTargetOrb.atomFaceList is not None
            ):
                extremeSourceOrb.atomFaceList = extremeTargetOrb.atomFaceList


def reactionSmilesFromOrbitalPair(
    filledOrb, unfilledOrb, canonize=False, arrowObjList=None
):
    """Given a filled and unfilled orbital pair,
    return the reaction SMILES string representing the orbital overlap.
    If the arrowObjList is provided, will also add respective curved arrow objects.

    WARNING: Function assumes that the composite mol object the orbitals are
    based upon already has atom mapping labels on the relevant orbital atoms.

    canonize:
        Option to relabel any atom mapping labels in a canonical way
        so the same reaction SMILES will be generated, even if different
        initial labels were used on the orbitals.
    """
    compositeMol = filledOrb.mol

    reactantsSmi = createStdAtomMapSmiString(compositeMol)
    bond = moveOrbitalElectrons(filledOrb, unfilledOrb, arrowObjList=arrowObjList)
    # Apply the electron movement
    productsSmi = createStdAtomMapSmiString(compositeMol)
    undoMoveOrbitalElectrons(filledOrb, unfilledOrb)
    # Revert the compositeMol back to its original (reactant) state

    # Extra standardization check because always a risk when move atoms and bonds around
    #   that the canonization algorithm gets confused
    # reactantsSmi = standardizeSmiles(reactantsSmi, True);
    # productsSmi = standardizeSmiles(productsSmi, True);

    reagentSmi = ""  # No recording of reagents for now

    reactionSmiles = str.join(
        REACTION_COMPONENT_DELIM, [reactantsSmi, reagentSmi, productsSmi]
    )
    # print 'reactionSmiles in OrbitalInteraction', reactionSmiles  #TODO remove this
    if canonize:
        reactionSmiles = aReactionCanonizer.canonizeReactionSmiles(reactionSmiles)
    return reactionSmiles


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
            log.debug("IS BOND DISSOCIATION (compatibleOrbitalPair)")
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


class ElementaryStepData(dict):
    """Capture all of the core information of interest
    regarding an elementary step (orbital interaction).
    """

    def __init__(self):
        """Default constructor"""
        self["reactionSmiles"] = None
        self["reactantSmilesList"] = []
        self["productSmilesList"] = []
        self["arrowStr"] = None
        self["freeEnergyChange"] = None
        self["transitionEnergy"] = None
        self["compositeMol"] = None
        # Should represent the reactants and act as basis for...
        self["filledOrb"] = None
        #   orbital objects that when applied to each other
        self["unfilledOrb"] = None
        #   should convert the reactants in the compositeMol into products
        self["changeHistory"] = []

        # Supplementary string versions of other data to facilitate temporary persistence as text data
        self["filledOrbInfoStr"] = None
        self["unfilledOrbInfoStr"] = None
        self["changeHistoryStr"] = None

    def invert(originalStep):
        """Return an inverted version of this elementary step, representing
        the inverse / retro reaction step.  Copy over data elements,
        except for ones which should not be applicable in the inverse case,
        but which may be relevant to represent the reaction conditions.
        """
        filledOrb = originalStep["filledOrb"]
        unfilledOrb = originalStep["unfilledOrb"]

        # Setup the inverse reaction to get the inverse transition energy
        if filledOrb.mol != unfilledOrb.mol:
            # Orbitals from distinct molecule objects.  Prepare a composite
            #   molecule object copy that both orbitals will be based upon,
            #   otherwise will cause problems if ever try moving electrons between the orbitals
            (
                (compositeFilledOrb, compositeUnfilledOrb),
                compositeMol,
            ) = Orbital.compositeOrbitalsFromDistinctMols([filledOrb, unfilledOrb])

            filledOrb = compositeFilledOrb
            unfilledOrb = compositeUnfilledOrb
        else:
            # Comes from the same molecule.  Make copies we're free to manipulate
            (
                (cloneFilledOrb, cloneUnfilledOrb),
                copyMol,
            ) = Orbital.cloneOrbitalsWithCommonMol([filledOrb, unfilledOrb])

            filledOrb = cloneFilledOrb
            unfilledOrb = cloneUnfilledOrb

        # Invert the reaction under consideration
        moveOrbitalElectrons(filledOrb, unfilledOrb)
        # Convert to  product state
        (invSourceOrb, invTargetOrb) = perceiveInvertedOrbitals(filledOrb, unfilledOrb)

        invStep = ElementaryStepData()

        # Copy over data elements from the original step that may represent reaction condition parameters,
        #   but exclude most that are specific to the forward direction version
        specificKeys = set(
            [
                "reactionSmiles",
                "reactantSmilesList",
                "productSmilesList",
                "arrowStr",
                "freeEnergyChange",
                "transitionEnergy",
                "compositeMol",
                "filledOrb",
                "unfilledOrb",
                "changeHistory",
                "filledOrbInfoStr",
                "unfilledOrbInfoStr",
                "changeHistoryStr",
            ]
        )
        for key, value in originalStep.iteritems():
            if key not in specificKeys:
                invStep[key] = value

        # Initialize with key data of relevance
        invStep["filledOrb"] = invSourceOrb
        invStep["unfilledOrb"] = invTargetOrb
        invStep["compositeMol"] = invSourceOrb.mol
        invStep["reactantSmilesList"] = originalStep["productSmilesList"]
        invStep["productSmilesList"] = originalStep["reactantSmilesList"]
        if originalStep["freeEnergyChange"] is not None:
            invStep["freeEnergyChange"] = -originalStep["freeEnergyChange"]

        return invStep

    invert = staticmethod(invert)

    def loadTextProperties(self):
        """Restore object data (composite molecule and orbitals) from text attributes.
        Useful for retrieving based on database data or
        persisting temporarily in HTML form data, etc.
        Other attributes can be restored more or less directly, just have to
        do format conversions.
        """
        reactionMol = OEGraphMol()
        OEParseSmiles(reactionMol, self["reactionSmiles"])
        clearMechanismLabels(reactionMol)
        # Separate reaction molecule into components
        unlabeledRxnSmi = createStandardSmiString(reactionMol)
        tokens = unlabeledRxnSmi.split(REACTION_COMPONENT_DELIM)
        reactantsSmi = tokens[0]
        reagentSmi = tokens[1]
        productsSmi = tokens[2]
        self["reactantSmilesList"] = reactantsSmi.split(SMILES_MOL_DELIM)
        self["productSmilesList"] = productsSmi.split(SMILES_MOL_DELIM)

        # Fresh copy, retaining atom map labels and only reactant side for orbital object basis
        reactionMol.Clear()
        reactantsSmi = self["reactionSmiles"].split(REACTION_COMPONENT_DELIM)[0]
        OEParseSmiles(reactionMol, reactantsSmi)
        self["compositeMol"] = reactionMol
        self["filledOrb"] = Orbital.fromLabeledMolAndInfoStr(
            reactionMol, self["filledOrbInfoStr"]
        )
        self["unfilledOrb"] = Orbital.fromLabeledMolAndInfoStr(
            reactionMol, self["unfilledOrbInfoStr"]
        )

        if self["changeHistoryStr"] is None or self["changeHistoryStr"] == "":
            self["changeHistory"] = []
        else:
            historyList = self["changeHistoryStr"].split(VALUE_DELIM)
            for iHistory, strValue in enumerate(historyList):
                historyList[iHistory] = int(strValue)
            self["changeHistory"] = historyList

        self["freeEnergyChange"] = float(self["freeEnergyChange"])
        self["transitionEnergy"] = float(self["transitionEnergy"])

    def changeHistoryStr(self):
        """Return the change history in the format of a single string."""
        strList = []
        for value in self["changeHistory"]:
            strList.append(str(value))
        return str.join(VALUE_DELIM, strList)

    def cumulativeChangeHistory(self):
        """Return the sum of all change history elements,
        useful for comparing against 0, which would indicate a step
        that was never actually used, or for getting a sense
        as to which steps have been the most important
        in affecting the reaction outcome.
        """
        return sum(self["changeHistory"])

    def isProductive(self):
        """Check whether this is a productive reaction or not.
        An unproductive reaction is one where the products and reactants
        are the same or equivalent.  There is no value in simulating
        such reactions.
        """
        changeDict = self.reactionChangeCounts(
            self["reactantSmilesList"], self["productSmilesList"]
        )

        for smiles, change in changeDict.iteritems():
            if change != 0:
                return True

        # Did not find any with a net change, means unproductive reaction
        return False

    def reactionChangeCounts(self, reactantSmiList, productSmiList):
        """Returns a dictionary keyed by the reactant and product
        SMILES strings lists provided.  Values equals the net change
        to the component that the reaction represents (-1 for each
        reactant occurence, +1 for each product occurence).
        """
        countDict = dict()

        for smiles in reactantSmiList:
            if smiles not in countDict:
                countDict[smiles] = 0
            countDict[smiles] -= 1

        for smiles in productSmiList:
            if smiles not in countDict:
                countDict[smiles] = 0
            countDict[smiles] += 1

        return countDict


class SupplementalStereoDataModel:
    """Simple struct to store stereochemical information about a previous electron movement"""

    def __init__(self, extType=EXT_TYPE_UNDEFINED):
        """Constructor, set some defaults"""
        self.extType = extType

        # If a sigma bond is broken, (and electrons are displaced to a chiral atom),
        #   the following should be set.
        # A list of the previous stereo
        self.oldStereoList = None
        self.removedAtom = None

        # If a
