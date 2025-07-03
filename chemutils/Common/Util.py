"""Miscellaneous utility functions used across the application
"""
from chemutils.Common import Const, Env
import sys
import os
import logging
import urllib
import time
import math
from queue import Queue, Empty, Full

try:
    from openeye.oechem import (
        oemolistream,
        oemolostream,
        OEFormat_SMI,
        OEFormat_ISM,
        OEFormat_CAN,
    )
    from openeye.oechem import (
        OEOFlavor_ISM_Default,
        OEOFlavor_ISM_AtomMaps,
        OEOFlavor_ISM_BondStereo,
        OEOFlavor_ISM_AtomStereo,
    )
    from openeye.oechem import (
        OEGraphMol,
        OEParseSmiles,
        OECreateIsoSmiString,
        OECreateSmiString,
        OEGetDefaultMass,
    )
    from openeye.oechem import OENetCharge, OEMolecularFormula
    from openeye.oechem import OEGetAtomicNum, OEGetAtomicSymbol
    from openeye.oechem import OEClearAromaticFlags, OEAssignAromaticFlags
    from openeye.oechem import *
except:
    # Looks like OpenEye dependency not installed.  Some functions won't work,
    #   but at least try to continue without crashing
    print(
        sys.stderr,
        "OEChem dependency does not appear to be installed.  Attempting to continue anyway.",
    )
    # Insert a couple fake values so that later calls won't crash completely (though they won't work properly either)
    oemolistream = list
    oemolostream = list
    OEFormat_SMI = -1

from chemutils.Common.Const import (
    FP_MAX,
    FP_MIN,
    FP_SIZE,
    SMILES_MOL_DELIM,
)
from chemutils.Common.MolStdValue import (
    BOLTZMANN_CONST,
    PLANCKS_CONST,
    ATOM_PER_MOL,
    ROOM_TEMPERATURE,
)

log = logging.getLogger(Const.APPLICATION_NAME)
log.setLevel(Const.LOGGER_LEVEL)

handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter(Const.LOGGER_FORMAT)

handler.setFormatter(formatter)
log.addHandler(handler)


def smi_to_unique_smi_fast(smi):
    """Creates a unique smiles with aromaticity perceived and (atom maps, atom stereo, bond stereo) omitted.
    Based on latest BP re: OEMolToSmiles()"""
    mol = OEGraphMol()
    OEParseSmiles(mol, smi)
    OEAssignAromaticFlags(mol)
    return OECreateSmiString(
        mol,
        OEOFlavor_ISM_Default
        - OEOFlavor_ISM_AtomMaps
        - OEOFlavor_ISM_AtomStereo
        - OEOFlavor_ISM_BondStereo,
    )


def smi_to_unique_smi_map(smi):
    """Creates a unique smiles with aromaticity perceived and (atom stereo, bond stereo) omitted.
    Based on latest BP re: OEMolToSmiles()"""
    mol = OEGraphMol()
    OEParseSmiles(mol, smi)
    OEAssignAromaticFlags(mol)
    return OECreateSmiString(
        mol, OEOFlavor_ISM_Default - OEOFlavor_ISM_AtomStereo - OEOFlavor_ISM_BondStereo
    )


def mol_to_unique_smi_fast(mol):
    """Creates a unique smiles with aromaticity perceived and atom maps omitted.
    Based on latest BP re: OEMolToSmiles()"""
    OEAssignAromaticFlags(mol)
    return OECreateSmiString(
        mol,
        OEOFlavor_ISM_Default
        - OEOFlavor_ISM_AtomMaps
        - OEOFlavor_ISM_AtomStereo
        - OEOFlavor_ISM_BondStereo,
    )


def to_unique_smi_fast(input):
    """Expects input to be either a smiles string, or an OEGraphMol().
    Returns OECreateSmiString(mol, OEOFlavor_ISM_Default - OEOFlavor_ISM_AtomMaps - OEOFlavor_ISM_AtomStereo - OEOFlavor_ISM_BondStereo)
    """
    if isinstance(input, (str, unicode)):
        return smi_to_unique_smi_fast(str(input))
    else:
        return mol_to_unique_smi_fast(input)


def exact_mass(mol):
    """Calculate the exact mass assuming all isotopes are present in their most abundant variations."""
    implicitH = 0
    mass = 0.0
    for atom in mol.GetAtoms():
        elemno = atom.GetAtomicNum()
        # most abundant isotopic mass
        default_isotopic_mass = OEGetIsotopicWeight(elemno, OEGetDefaultMass(elemno))
        implicitH += atom.GetImplicitHCount()
        mass += default_isotopic_mass
    mass += implicitH * OEGetIsotopicWeight(OEElemNo_H, OEGetDefaultMass(OEElemNo_H))
    return mass


def isStdFile(filename):
    """Given a filename, determine if it is meant to represent
    sys.stdin / sys.stdout or just a regular file based on the
    convention of the filename "-"
    """
    sourceDir, sourceBase = os.path.split(filename)
    return filename == Const.STD_FILE or sourceBase in Const.STD_MOL_EXT.values()


def stdOpen(filename, mode="r", stdFile=None):
    """Wrapper around basic file "open" method.  If specified filename
    is the recognized as representing stdin or stdout (see "isStdFile")
    then return the designated stdFile object (presumably sys.stdin or sys.stdout).

    Also provide automatic unzipping of gzipped files.  If find default
    GZip file extension, will automatically open the file using gzip opener
    for transparent access.
    """
    if stdFile == None:
        if mode.startswith("w"):
            stdFile = sys.stdout
        else:
            stdFile = sys.stdin

    if isStdFile(filename):
        return stdFile
    else:
        if filename.endswith(Const.GZIP_EXT):
            import gzip

            return gzip.open(filename, mode)
        else:
            return open(filename, mode)


def getMaxMapIdx(mol):
    """Simple function to determine the maximum current map idx on a molecule"""
    maxMapIdx = 0
    for atm in mol.GetAtoms():
        currMapIdx = atm.GetMapIdx()
        if currMapIdx > maxMapIdx:
            maxMapIdx = currMapIdx

    return maxMapIdx


def molBySmiles(smiles):
    """Return an OEGraphMol by smiles"""
    mol = OEGraphMol()
    OEParseSmiles(mol, smiles)
    return mol


def molToSmiList(molList):
    """Convert a list of molecule objects into a list of respective SMILES strings."""
    smiList = []
    for mol in molList:
        smiList.append(createStandardSmiString(mol))
    return smiList


def splitCompositeMolToSmilesList(
    compositeMol, retainCounterIons=False, retainAtomMaps=False
):
    """Given a (composite) molecule, return a list of SMILES strings,
    with one element per component molecule of the composite.
    """
    if retainAtomMaps:
        compositeSmi = createAtomMapSmiString(compositeMol)
    else:
        compositeSmi = createStandardSmiString(compositeMol)
    return splitCompositeSmilesToList(compositeSmi, retainCounterIons)


def splitCompositeSmilesToList(compositeSmi, retainCounterIons=False):
    """retainCounterIons option will try to not separate out components
    if they have opposing (and neutralizing) formal charges.
    """
    smilesList = compositeSmi.split(SMILES_MOL_DELIM)
    if retainCounterIons:
        # Sort out by net formal charges on each component
        chargeSmiList = []
        mol = OEGraphMol()
        for smi in smilesList:
            OEParseSmiles(mol, smi)
            chargeSmiList.append((OENetCharge(mol), smi))
            mol.Clear()
        chargeSmiList.sort()

        # Constants to extract tuple components of the chargeSmiList items
        CHARGE = 0
        SMILES = 1

        # Now put these back to the smilesList, but combine any charged components
        #   from either end of the sorted list
        smilesList = []
        iFirst = 0
        iLast = len(chargeSmiList) - 1

        currComponent = []
        currNetCharge = 0
        while iFirst <= iLast:
            if len(currComponent) < 1:
                # Haven't added any components yet, just add the first one
                currComponent.append(chargeSmiList[iFirst][SMILES])
                currNetCharge += chargeSmiList[iFirst][CHARGE]
                iFirst += 1
            elif currNetCharge < 0:
                # Negative charge, look for a positive addition at end of the list to neutralize
                if chargeSmiList[iLast][CHARGE] > 0:
                    currComponent.append(chargeSmiList[iLast][SMILES])
                    currNetCharge += chargeSmiList[iLast][CHARGE]
                    iLast -= 1
                else:
                    # No positive charged components to add for neutralization, have to just accept the lone ion then
                    smilesList.append(str.join(SMILES_MOL_DELIM, currComponent))
                    currComponent = []
                    currNetCharge = 0
            elif currNetCharge > 0:
                # Positive charge, look for a negative addition at front of the list to neutralize
                if chargeSmiList[iFirst][CHARGE] < 0:
                    currComponent.append(chargeSmiList[iFirst][SMILES])
                    currNetCharge += chargeSmiList[iFirst][CHARGE]
                    iFirst += 1
                else:
                    # No negative charged components to add for neutralization, have to just accept the lone ion then
                    smilesList.append(str.join(SMILES_MOL_DELIM, currComponent))
                    currComponent = []
                    currNetCharge = 0
            elif currNetCharge == 0:
                # Found enough components to get a neutral component
                smilesList.append(str.join(SMILES_MOL_DELIM, currComponent))
                currComponent = []
                currNetCharge = 0
            else:
                raise Exception("Sanity check, this shouldn't happen")
        # Whatever's left as a component, have to accept it
        smilesList.append(str.join(SMILES_MOL_DELIM, currComponent))
    return smilesList


def splitCompositeMol(compositeMol, retainCounterIons=False, retainAtomMaps=False):
    """Given a (composite) molecule, return a list of molecule objects,
    with one element per component molecule of the composite.
    Based on isomeric SMILES string delimiters, so no 3D or other annotation data is retained.
    """
    componentMolList = []
    smilesList = splitCompositeMolToSmilesList(
        compositeMol, retainCounterIons, retainAtomMaps
    )
    for smiles in smilesList:
        componentMolList.append(molBySmiles(smiles))
    return componentMolList


def joinSmilesListToCompositeSmiles(smiList):
    """Simple convenience to join together a set of smiles"""
    return SMILES_MOL_DELIM.join(smiList)


def joinMolListToCompositeSmiles(molList, retainAtomMaps=False):
    """Given a list of mols, return a string with the smiles of a composite mol."""
    if retainAtomMaps:
        return SMILES_MOL_DELIM.join([createAtomMapSmiString(mol) for mol in molList])
    else:
        return SMILES_MOL_DELIM.join([createStandardSmiString(mol) for mol in molList])


def atomSetFromSmiles(smi):
    """Convenience to get a set of atomic numbers for all of the atoms in a particular set"""
    mol = molBySmiles(smi)
    return set([atm.GetAtomicNum() for atm in mol.GetAtoms()])


def joinMolListToCompositeMol(molList, retainAtomMaps=False):
    """Given a list of mols, return a new mol obj with all as a composite mol"""
    compSmiStr = joinMolListToCompositeSmiles(molList, retainAtomMaps)
    return molBySmiles(compSmiStr)


def isSmilesFormat(formatCode):
    """Check if the supplied formatCode is one of the OpenEye SMILES formats."""
    return formatCode in (OEFormat_SMI, OEFormat_ISM, OEFormat_CAN)


class virtual_oemolistream(oemolistream):
    """Similar to using StringIO to simulate a virtual file object.
    Use this to simulate a virtual oemolistream object that
    reads from a string (molstring) instead of from a file.

    Since no file names is specified, can't automatically guess
    the molecule format from the file extension, thus it
    must be explicitly specified (or assumed to be SMILES).
    """

    emptyStream = False

    def __init__(self, molstring, format=OEFormat_SMI):
        oemolistream.__init__(self)
        if molstring.strip() == "":
            self.emptyStream = True
        else:
            self.openstring(molstring)
            self.SetFormat(format)

    def GetOEGraphMols(self):
        """Overridden method to account for special behavior if input is an empty string"""
        if self.emptyStream:
            return []
        else:
            return oemolistream.GetOEGraphMols(self)


class virtual_oemolostream(oemolostream):
    """Similar to using StringIO to simulate a virtual file object.
    Use this to simulate a virtual oemolostream object that
    writes to a string instead of to a file.

    Since no file names is specified, can't automatically guess
    the molecule format from the file extension, thus it
    must be explicitly specified (or assumed to be SMILES).

    To retrieve the contents of the virtual file output after
    done with the oemolostream, use the method oemolostream.GetString().
    """

    def __init__(self, format=OEFormat_SMI):
        oemolostream.__init__(self)
        self.SetString("")
        # Not sure if this is best method.  OpenEye doesn't document this API
        self.SetFormat(format)


def generateFingerprint(smiles, fpMax=FP_MAX, fpMin=FP_MIN, fpSize=FP_SIZE):
    """Returns a standard fingerprint object for a given (SMILES) string.
    Just defer to basic FINGER.Finger(...) constructor, but check to see if
    want to use the web service instead.  Useful for Windows clients
    that cannot import the FINGER module.
    """
    fingerprint = None
    if not Env.FP_USE_WEB_SERVICE:
        try:
            # Try to generate fingerprint directly
            from FINGER import Encoder

            enc = Encoder(length=fpSize, maxdepth=fpMax, mindepth=fpMin)
            fingerprint = enc(smiles)
        except:
            fingerprint = None

    if fingerprint is None:
        # Either option set to specifically use the web service,
        #   or the alternative failed, probably because FINGER module not compiled
        #   or installed on local system, try pulling it from the
        #   http based Fingerprint web service.
        #   Note that this isn't quite as good, because it only gets the
        #   fingerprint string, not object.  Good enough for most purposes though.
        query = urllib.quote_plus(
            "ENCODE %s %d %d %d\n" % (smiles, fpMax, fpMin, fpSize)
        )
        url = "%s%s" % (Env.FP_WEB_SERVICE, query)
        fingerprint = urllib.urlopen(url).read().strip()

        # Only got a string, not object version (probably on a Windows machine)
        # Use Python Finger object to parse it
        from CHEM.datatype.pyFINGER import Finger

        fingerprint = Finger(fingerprint)

    return fingerprint


def clearAtomMapsSmiStr(smi):
    """Convenience to clean up a smi string.  (Basically remove atom maps)

    Because the setting of hydrogen isotopes to 1 is used to ensure atom maps elsewhere,
    remove all isotopes of 1 on the hydrogens here.

    Note: this does mess up Hydrogens that are truly an isotope."""
    from chemutils.Common.MolExt import clearAtomMaps

    mol = molBySmiles(smi)

    # This clears the isotopes as well.
    clearAtomMaps(mol)
    return createStandardSmiString(mol)


def clearAromaticFlagsSmiStr(smi):
    """Convenience to clear out aromatic flags of a smi string, and return a smi string"""
    mol = molBySmiles(smi)
    OEClearAromaticFlags(mol)
    return createAtomMapSmiString(mol)


def assignAromaticFlagsSmiStr(smi):
    """Convenience to assign aromatic flags of a smi string, and return a smi string"""
    mol = molBySmiles(smi)
    OEAssignAromaticFlags(mol)
    return createAtomMapSmiString(mol)


def standardizeSmiles(smi, ensureAtomMaps=False, kekulize=False, aromatize=False):
    """Given a SMILES string, generate a standardized isomeric SMILES version.
    Option to use createAtomMapSmiString to ensure any hydrogens with atom mappings
    are retained in the generated SMILES string.
    Also add option to kekulize.
    """
    mol = molBySmiles(smi)
    if kekulize:
        OEClearAromaticFlags(mol)

    if aromatize:
        OEAssignAromaticFlags(mol)

    stdSmi = None
    if ensureAtomMaps:
        stdSmi = createAtomMapSmiString(mol)
    else:
        stdSmi = OECreateIsoSmiString(mol)
    return stdSmi


def createStandardSmiString(mol):
    """Strange.  Should be perfectly redundant with OECreateIsoSmiString, but in some cases,
    the direct call messes up, such as alkene stereochemistry giving C\\C=C\\C instead of C/C=C/C.
    Functionally the same, but distinct when doing unique SMILES checks.
    """
    if mol is None:
        return ""

    smi = OECreateIsoSmiString(mol)
    return standardizeSmiles(smi)


def makeNonMappedHydrogensImplicit(mol):
    """To ensure common format, we need to make nonMapped Hydrogens to be implicit"""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() == 0:
            bonds = [b for b in atom.GetBonds()]
            if len(bonds) == 1:
                nbr = bonds[0].GetNbr(atom)
                implHCount = nbr.GetImplicitHCount()
                mol.DeleteAtom(atom)
                nbr.SetImplicitHCount(implHCount + 1)


def createStdAtomMapSmiString(mol, semiImplicitH=True):
    """Just a wrapper to ensure proper flags are set to handle hydrogens.

    This is a replacement for the createAtomMapSmiString function (that should work better)
    ie, it should NOT alter aromaticity or kekulization.
    """
    from openeye import oechem as oe

    nMol = oe.OEGraphMol(mol)
    if semiImplicitH:
        makeNonMappedHydrogensImplicit(nMol)
    flavor = (
        oe.OESMILESFlag_DEFAULT
        | oe.OESMILESFlag_Hydrogens
        | oe.OESMILESFlag_AtomStereo
        | oe.OESMILESFlag_BondStereo
    )
    return oe.OECreateSmiString(nMol, flavor)


def createAtomMapSmiString(mol):
    """Should be the same as createStandardSmiString, but beware of molecules
    that should generate a SMILES string like [H:1][Br:2].
    The standard canonization function will ignore
    "implicit hydrogens" and thus yield [BrH:2] as the result, losing the atom
    mapping information on the hydrogen.  Can hack around this, though beware
    that doing so will mess up any hydrogens specified with an isotope of 1
    like [1H]Br.
    """
    originalIsotopes = dict()
    for atom in mol.GetAtoms():
        originalIsotopes[atom.GetIdx()] = atom.GetIsotope()
        if atom.GetAtomicNum() == 1 and atom.GetMapIdx() > 0:
            # Hydrogen atom with atom mapping.  Artificially set isotope,
            #   otherwise will be lost as "implicit hydrogen" upon canonization
            atom.SetIsotope(1)

    smi = createStandardSmiString(mol)
    smi = smi.replace("[1H", "[H")
    # Get rid of the artificial 1H isotopes at this point

    # Revert the molecule object back to its original state
    for atom in mol.GetAtoms():
        atom.SetIsotope(originalIsotopes[atom.GetIdx()])

    return smi


def atomByMapIdx(mol, mapIdx):
    """Convenience function to grab an atom by mapping index."""
    atm = [atm for atm in mol.GetAtoms() if atm.GetMapIdx() == mapIdx]
    if atm == []:
        atm = None
    else:
        atm = atm[0]
    return atm


def atomCount(mol, atomicNumSet=None):
    """Return number of atoms in the molecule, including hydrogens.
    If atomicNumSet parameter is provided, then only count atoms
    that have an atomic number which is in the set.
    """
    total = 0
    for atom in mol.GetAtoms():
        if atomicNumSet is None or atom.GetAtomicNum() in atomicNumSet:
            total += 1
        if atomicNumSet is None or 1 in atomicNumSet:
            # 1 is atomic number of hydrogen, check if we should add these
            total += atom.GetImplicitHCount()
    return total


# Characters representing formal charges
CHARGE_CHARS = ("-", "+")


def recordAtomCount(currentAtom, currentCount, atomCountDict):
    """Utility function used by standardMolecularFormula"""
    if currentAtom != "":
        if currentAtom not in CHARGE_CHARS:
            atomicNum = OEGetAtomicNum(currentAtom)
            # Re-convert through atomic number
            currentAtom = OEGetAtomicSymbol(atomicNum)
            #   to ensure standard representation (esp. capitalization)
        if currentCount != "":
            currentCount = int(currentCount)
        else:
            # No count number specified, default to 1
            currentCount = 1
        atomCountDict[currentAtom] = currentCount


def standardizeMolecularFormula(molecularFormula):
    """Read in a string representing a molecular formula and extract
    out the atom contents.  Feed these into a temporary molecule
    (no meaningful connectivity information required) and then generate
    the standard molecular formula representation.
    """

    # Dictionary to keep track of how many times each atom type appears
    atomCountDict = dict()

    currentAtom = ""
    currentCount = ""
    for char in molecularFormula:
        if char.isupper() or char in CHARGE_CHARS:
            # Upper-case character or formal charge designation, must be starting a new atom type
            recordAtomCount(currentAtom, currentCount, atomCountDict)
            # Reset tracking variables to prepare for next token
            currentAtom = ""
            currentCount = ""
        elif char.isalpha() and currentCount != "":
            # Starting a new alphabetic token when in the middle of a numeric, must be starting new atom type
            recordAtomCount(currentAtom, currentCount, atomCountDict)
            # Reset tracking variables to prepare for next token
            currentAtom = ""
            currentCount = ""

        if char.isalpha() or char in CHARGE_CHARS:
            currentAtom += char
        elif char.isdigit():
            currentCount += char

    # Record final token
    recordAtomCount(currentAtom, currentCount, atomCountDict)

    # Construct a mock SMILES string representing a molecule of the given molecular formula
    mockSmi = []
    for atomSymbol, atomCount in atomCountDict.iteritems():
        if atomSymbol not in CHARGE_CHARS:
            smiToken = "[%s]" % atomSymbol
            mockSmi.append(smiToken * atomCount)
    mockSmi = str.join("", mockSmi)

    mol = OEGraphMol()
    OEParseSmiles(mol, mockSmi)
    # Check if we need to add formal charges
    for chargeChar in CHARGE_CHARS:
        if chargeChar in atomCountDict:
            # Select a random atom from the molecule and assign the respective formal charge
            # Assumes only one net charge type will be found.  Shouldn't see + and - charges in net formula
            for atom in mol.GetAtoms():
                chargeCount = atomCountDict[chargeChar]
                formalCharge = int("%s%s" % (chargeChar, chargeCount))
                # Convert charge string into an integer
                atom.SetFormalCharge(formalCharge)
                break
                # Only apply to one atom

    stdMolFormula = OEMolecularFormula(mol)
    return stdMolFormula


class ProgressDots:
    """Clone of OEChem OEDots class, to add progress indicator to long processes,
    without actually requiring OEChem as a dependency.
    """

    def __init__(
        self,
        big=Const.PROG_BIG,
        small=Const.PROG_SMALL,
        name="items",
        stream=sys.stderr,
    ):
        """Constructor.
        big - Number of updates before completing a progress output line.
        small - Number of updates before outputting a progress dot.
        name - Name of the items being processed.
        stream - Stream to send progress output to.  Defaults to sys.stderr.
        """
        self.big = big
        self.small = small
        self.name = name
        self.stream = sys.stderr
        self.count = 0
        self.start = time.time()

    def Update(self, step=1):
        """Update the progress counter by an increment of size step (default=1).
        Well output progress dots or line information to the stream if
        reached an appropriate big or small increment.
        """
        self.count += step
        if self.small > 0 and self.count % self.small == 0:
            self.stream.write(".")
        if self.big > 0 and self.count % self.big == 0:
            self.PrintStatus()

    def GetCounts(self):
        """Get the current count of updates"""
        return self.count

    def GetTime(self):
        """Get the time (in seconds) since the progressindicator was created."""
        return time.time() - self.start

    def PrintStatus(self):
        # Assume 0 values means no feedback desired
        if self.small > 0 and self.big > 0:
            print(
                self.stream,
                "%d %s processed after %d seconds."
                % (self.count, self.name, self.GetTime()),
            )


def loadObjectsFromQueuesById(objectsById, objectQueuesById, neededObjectIds):
    """
    objectsById:
        Dictionary to be populated with (id: object) instances pulled out of the
        respective objectQueues from objectQueuesById
    objectQueuesById:
        The "cache" of objects, keyed by ID.  Actually store a Queue of possibly multiple objects
        for each ID, and just retrieve the first available to store in objectsById.
    neededObjectIds:
        Set of IDs of the objects to look for and load from the objectQueues
    return:
        Returns the set of "newObjectIds."  That is, the IDs of objects which were needed but were NOT
        found to be available in the objectQueues.  Thus, caller will have to make own arrangements
        to fill in these blanks.
    """
    newObjectIds = set()
    for objectId in neededObjectIds:
        if objectId not in objectQueuesById:
            objectQueuesById[objectId] = Queue()
        objectQueue = objectQueuesById[objectId]
        try:
            cacheObject = objectQueue.get_nowait()
            objectsById[objectId] = cacheObject
        except Empty:
            # No copy of this object is available pre-cached.
            # Will ultimately have to find for a fresh one
            newObjectIds.add(objectId)

    return newObjectIds


def saveObjectsToQueuesById(objectsById, objectQueuesById):
    """
    objectsById:
        Dictionary of (id: object) pairs to deposit in the
        respective objectQueues of objectQueuesById
    objectQueuesById:
        The "cache" of objects, keyed by ID.  Actually store a Queue of possibly multiple objects
        for each ID, and just add any available in objectsById.
    """
    for objectId in objectsById.keys():
        # Deposit objects back to cache queues
        cacheObject = objectsById.pop(objectId)
        # Retrieve AND remove object from dictionary, to prevent usage while it's sitting in cache
        objectQueue = objectQueuesById[objectId]
        try:
            objectQueue.put_nowait(cacheObject)
        except Full:
            # If object queue is full, forget it then, just discard this instance
            del cacheObject


class Codec:
    """Simple class to "encrypt" and "decrypt" strings since the rotor
    module got deprecated after Python 2.3.
    Pretty simple now, but just to lightly obfuscate the text
    to deter spurious disclosures.
    """

    DEFAULT_DELIM_CHAR = "%"
    CUSTOM_DELIM_CHAR = "@"
    ENCRYPTION_TAG = "@@@"
    # Tag expected to begin and end encrypted strings

    # String of characters that will be used in the encryption maps
    BASE_CHARS = "-=~!@#$%^&*()[]{}|\;':,./<>?abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

    def __init__(self, key):
        """Constructor.
        Caller specifies the encryption key.
        """
        self.key = key
        self.encryptMap = None
        self.decryptMap = None

        (self.encryptMap, self.decryptMap) = self.generateTranslationMaps(self.key)

    def generateTranslationMaps(self, keyStr):
        """Generate an encoding and decoding translation map,
        adjusted by the given key string.
        """
        encryptMap = dict()
        decryptMap = dict()
        currentLength = 0

        # Setup initial set of mappings based on key string.  Look out for redundancies.
        for keyChar in keyStr:
            if keyChar not in decryptMap:
                sourceChar = self.BASE_CHARS[currentLength]
                decryptMap[keyChar] = sourceChar
                encryptMap[sourceChar] = keyChar
                currentLength += 1

        # Fill in the rest of the alphabet not covered in the key string
        for targetChar in self.BASE_CHARS:
            # Check if we have not already covered this in the key string
            if targetChar not in decryptMap:
                sourceChar = self.BASE_CHARS[currentLength]
                decryptMap[targetChar] = sourceChar
                encryptMap[sourceChar] = targetChar
                currentLength += 1

        return (encryptMap, decryptMap)

    def isEncrypted(self, source):
        return (
            source is None
            or source.startswith(self.ENCRYPTION_TAG)
            and source.endswith(self.ENCRYPTION_TAG)
        )

    def applyTranslation(self, sourceStr, translateMap):
        """Simple character replacement encryption based on the provided translation map"""
        targetStrList = []
        # Build up as an array of characters
        for sourceChar in sourceStr:
            targetChar = sourceChar
            # Default to source char unless a translation is found
            if sourceChar in translateMap:
                targetChar = translateMap[sourceChar]
            targetStrList.append(targetChar)

        return str.join("", targetStrList)

    def reverseStr(self, source):
        """Extra flourish to obscure text by reversal of string"""
        target = list(source)
        target.reverse()
        target = str.join("", target)
        return target

    def encrypt(self, source):
        target = None
        if source is not None:
            target = self.applyTranslation(source, self.encryptMap)
            target = urllib.quote(target)
            # Extra flourishes to obscure symbols
            target = self.reverseStr(target)
            targetList = [self.ENCRYPTION_TAG, target, self.ENCRYPTION_TAG]
            target = str.join("", targetList)
        return target

    def decrypt(self, source):
        target = None
        if source is not None:
            if not self.isEncrypted(source):
                raise Exception(
                    "This string does not appear to be encrypted: '%s'" % source
                )

            tagLength = len(self.ENCRYPTION_TAG)
            detagSource = source[tagLength:-tagLength]

            target = self.reverseStr(detagSource)
            target = urllib.unquote(target)
            target = self.applyTranslation(target, self.decryptMap)
        return target


class BoltzmannProbabilityWeight:
    """Functor which, given an energy value, should provide
    a proportional probability weight based on Boltzmann probability distributions.
    Basically just e^(-E/N kB T)

    e: Base of the natural logarithm
    E: Input energy value, expected to be in units of kcal / mol
    N: Avogadro's number to translate kcal / mol to simple kcal (per particle)
    kB: Boltzmann Constant (kcal / Kelvin)
    T: Absolute temperature (Kelvin)

    http://en.wikipedia.org/wiki/Boltzmann_distribution.

    Use to assign relative probabilities of existence of different species
    with relative energy values.
    """

    # Temperature at which to evaluate the probability
    temperature = None
    # Total exponential scalar
    scalar = None

    def __init__(self, temperature=ROOM_TEMPERATURE):
        """Initialization constructor expects the temperature at which
        to calculate the function.  If none specified, assume room temperature.
        """
        self.setTemperature(temperature)

    def setTemperature(self, temperature):
        self.temperature = temperature
        self.scalar = -1.0 / (ATOM_PER_MOL * BOLTZMANN_CONST * self.temperature)

    def __call__(self, energy):
        # print energy
        # print energy * self.scalar
        return math.exp(energy * self.scalar)


class EyringFormula:
    """Functor which, given a transition state energy value,
    provides an exact theoretical estimate for the reaction rate
    of an elementary reaction step.

    e: Base of the natural logarithm
    E: Input energy value, expected to be in units of kcal / mol
    N: Avogadro's number to translate kcal / mol to simple kcal (per particle)
    kB: Boltzmann Constant (kcal / Kelvin)
    h: Planck's Constant (kcal * sec)
    T: Absolute temperature (Kelvin)

    rate = (kB T / h) * e^(-E/N kB T)

    http://en.wikipedia.org/wiki/Eyring_equation
    http://www.chemie.uni-regensburg.de/Organische_Chemie/Didaktik/Keusch/eyr-e.htm

    Units for rate?  Works out to be Hz (1/sec), but how does this map to rate constants (k)
    having different units depending on the order of the reaction?  Have to eventually multiply
    rate constant by reactant concentrations for an overall reaction rate, since the concentration
    levels will influence the likelihood of reaction.  First order rate constants normally are
    in units of 1/sec, but second order is usually 1/M*s so that the overall reaction rate in
    both cases ends up being M/s (molarity per second).
    Reference above outlines the derivation and implies how the the rate constant should still work out.
    There is a step in the derivation translating a Keq to dG which is fair, but ends up losing
    a M unit in the Keq (assuming a bimolecular reaction), indicating where the rate unit was lost.
    This means the formula is actually derived based on assumptions for a bimolecular reaction.
    Should be able to work for any molecularity reaction then, but roll in entropy factors into
    transition state energy to account for the reactivity differences.

    You may note that the exponential term is just the Boltzmann probability distribution
    function and that the overall formula takes the form of the Arrhenius equation,
    except that the A coefficient is explicitly defined in terms of theoretical
    components instead of empirically determined.
    """

    # Temperature at which to evaluate the probability
    temperature = None
    # Based on Boltzmann probability distribution function
    boltzmannProbCalc = None
    # Total front coefficient, only need to calculate once
    coefficient = None

    def __init__(self, temperature=ROOM_TEMPERATURE):
        """Initialization constructor expects the temperature at which
        to calculate the function.  If none specified, assume room temperature.
        """
        self.temperature = temperature
        self.boltzmannProbCalc = BoltzmannProbabilityWeight(temperature)
        self.coefficient = BOLTZMANN_CONST * temperature / PLANCKS_CONST

    def __call__(self, energy):
        return self.coefficient * self.boltzmannProbCalc(energy)


def main(argv):
    """Main method, callable from command line"""
    pass


if __name__ == "__main__":
    main(sys.argv)
