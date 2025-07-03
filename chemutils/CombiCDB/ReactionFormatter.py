import sys
from optparse import OptionParser
from chemutils.Common.Util import stdOpen, virtual_oemolistream, createStandardSmiString
from openeye.oechem import (
    oemolistream,
    oemolostream,
    OEReadMolecule,
    OEWriteMolecule,
    OEFormat_SDF,
)
from openeye.oechem import OEGraphMol, OEParseSmiles
from openeye.oechem import OECanonicalOrderAtoms
from openeye.oechem import OEAddMols, OESetSDData
from openeye.oechem import OERxnRole_None, OERxnRole_Reactant, OERxnRole_Product


def main(argv):
    """Main method, callable from command line"""

    usageStr = (
        "usage: %prog [options] <inputFile> <outputFile> \n"
        + "   <inputFile>     File containing reaction information, such as SMIRKS / SMILES or RDF / SDF\n"
        + "   <outputFile>    File for same reaction information but in different file format (.smi or .sdf, .rdf not currently supported for output).\n"
    )

    parser = OptionParser(usage=usageStr)
    parser.add_option(
        "-c",
        "--canonizeAtomMapping",
        action="store_true",
        dest="canonizeAtomMapping",
        help="If set, before outputting reaction molecuel data, will canonize atom mapping labels and atom ordering",
    )
    (options, args) = parser.parse_args(argv[1:])

    if len(args) >= 2:
        inputFilename = args[0]
        outputFilename = args[1]

        inputMolStream = reactionmolistream(inputFilename)
        outputMolStream = oemolostream(outputFilename)

        if not options.canonizeAtomMapping:
            # Just do simple read / write / translate
            for mol in inputMolStream.GetOEGraphMols():
                OEWriteMolecule(outputMolStream, mol)
        else:
            # Don't just read / write / translate the molecule.  Canonize reaction information as needed
            reactionCanonizer = ReactionCanonizer()
            reactionCanonizer.runByMolStream(inputMolStream, outputMolStream)

        inputMolStream.close()
        outputMolStream.close()
    else:
        parser.print_help()
        sys.exit(-1)


def reactionmolistream(inputFilename):
    """Simple function to mimic oemolistream, but checks for possibility
    of RDF input, then used specialized RDF reader
    """
    inputMolStream = None
    lowerInputFilename = inputFilename.lower()
    if lowerInputFilename.endswith(".rdf") or lowerInputFilename.endswith(".rdf.gz"):
        # Looks like an RDF input file.  Need specialized reader.
        # Standard OEChem oemolistream doesn't handle as would like
        if lowerInputFilename.endswith(".gz"):
            # Looks like it's gzipped, open it accordingly
            import gzip

            inputMolStream = RDFReader(gzip.open(inputFilename))
        elif lowerInputFilename == ".rdf":
            # Looks like special indication of using std input, but as RDF format
            inputMolStream = RDFReader(sys.stdin)
        else:
            # "Normal" RDF file, use the RDFReader to encase it in an oemolistream equivalent
            inputMolStream = RDFReader(stdOpen(inputFilename))
    else:
        # Not an RDF file, just use standard oemolistream
        return oemolistream(inputFilename)
    return inputMolStream


class RDFReader:
    """Reader for RDF reaction files.
    OEChem supposedly supports reading these already, in terms of listing
    out all of the molecules from the reactions, but it loses important information.
    - It misses the very first molecule of the first reaction completely,
    - It doesn't keep the reaction components together, so it's unclear
        where one reaction begins and another ends, as well as
        being unclear whether each molecule encountered is a reactant or product.
    - Data annotations applied to the reaction as a whole are missed.

    Given an RDF file as text input, will generate an iterator over
    the reactions in the file with a single composite reaction molecule
    object representing each reaction, not each component.  Furthermore,
    each of these "reaction molecules" will include all of the relevant
    annotations as SD data pairs.

    Probably don't need to instantiate this directly, use reactionmolistream instead,
    which takes a filename and figures out whether to use this wrapper or just
    a standard oemolistream based on the filename extension.
    """

    # Constants used to parse out the RDF contents
    TAG_PREFIX = "$"
    REACTION_BEGIN_TAG = "$RXN"
    MOLECULE_BEGIN_TAG = "$MOL"
    MOLECULE_END_TAG = "M  END"
    REACTANT_PRODUCT_COUNT_LINE = 4
    ANNOTATION_NAME_TAG = "$DTYPE"
    ANNOTATION_VALUE_TAG = "$DATUM"
    ANNOTATION_MOL_TAG = "$MFMT"
    ANNOTATION_SMI_TAG = "$SMI"
    ANNOTATION_WRAP_COLUMN = 80
    ANNOTATION_WRAP_CHAR = "+"

    inputFile = None

    def __init__(self, inputFile):
        """Initialization constructor, taking the text input file stream
        to read from.
        """
        self.inputFile = inputFile

    def __iter__(self):
        """Primary method, produce an iterator over the reaction contents of the file."""

        iReaction = 0
        reactionMol = None
        # The reaction molecule object that is "yielded" each time
        iReactionLine = None
        # Line number within the current reaction being parsed
        nReactants = None
        nProducts = None

        annotationName = None
        annotationValue = None
        incompleteName = False
        incompleteValue = False
        molAnnotation = False
        # Tracks whether the annotation actually represents a formatted molecule

        # Molecule text and object for each component of
        #   a reaction parsed out (reactant or product, etc.)
        componentText = None
        componentMol = OEGraphMol()

        # Molecule text and object for "$MFMT" annotations
        annotationMolText = None
        annotationMol = OEGraphMol()

        # Parse through the text file to idenfity how to delimit and repackage the data
        for line in self.inputFile:
            lineWithoutNewline = line[:-1]
            # Strip the ending new line "\n" character

            if (
                annotationValue is not None
                and not incompleteValue
                and annotationValue != self.ANNOTATION_MOL_TAG
            ):
                # Already have an annotation value, and doesn't look like a special, long, "wrapped" line or molecule format case
                # Check the next (current) line.  If it looks like starting a new tag, then
                #   assume this annotation is completed and save it.
                #   Otherwise, assume this is a continued part of a multi-line annotation, separated by a newline "\n"
                if line.startswith(self.TAG_PREFIX) or len(line.strip()) < 1:
                    # Starting a new tag, or reached the end of the reaction definition
                    OESetSDData(reactionMol, annotationName, annotationValue)

                    # Reset annotation fields
                    annotationName = None
                    annotationValue = None
                    incompleteName = False
                    incompleteValue = False
                else:
                    # Looks like continuing a multi-line case
                    annotationValue += "\n" + lineWithoutNewline
                    incompleteValue = len(
                        lineWithoutNewline
                    ) > self.ANNOTATION_WRAP_COLUMN and annotationValue.endswith(
                        self.ANNOTATION_WRAP_CHAR
                    )
                    if incompleteValue:
                        # Strip off the wrap / join character
                        annotationValue = annotationValue[
                            : -len(self.ANNOTATION_WRAP_CHAR)
                        ]

            if (
                annotationValue is None
                and annotationName is not None
                and not incompleteName
            ):
                # Same idea as annotationValue.  Beware of multi-line values
                if line.startswith(self.TAG_PREFIX):
                    # Starting a new tag, probably the annotation value.  Just let the annotation name be completed then
                    pass
                else:
                    # Looks like continuing a multi-line case
                    annotationName += "\n" + lineWithoutNewline
                    incompleteName = len(
                        lineWithoutNewline
                    ) > self.ANNOTATION_WRAP_COLUMN and annotationName.endswith(
                        self.ANNOTATION_WRAP_CHAR
                    )
                    if incompleteName:
                        # Strip off the wrap / join character
                        annotationName = annotationName[
                            : -len(self.ANNOTATION_WRAP_CHAR)
                        ]

            if line.startswith(self.REACTION_BEGIN_TAG):
                # Starting up a new reaction.  Prepare a new molecule object
                #   to store the data, and if there was a previous one, "yield" it
                #   to the iterator.
                if reactionMol is not None:
                    yield reactionMol
                    reactionMol.Clear()
                    # Blank for a new, fresh object without requiring another instantiation
                    iReaction += 1

                    # Verify that all tracking variables are cleared and didn't leave anything hanging
                    if annotationName is not None or incompleteName:
                        raise Exception(
                            'Starting new reaction %d before completing annotation: "%s"'
                            % (iReaction, annotationName)
                        )
                    if annotationValue is not None or incompleteValue:
                        raise Exception(
                            'Starting new reaction %d before completing annotation value: "%s"'
                            % (iReaction, annotationValue)
                        )
                    if componentText is not None:
                        raise Exception(
                            'Starting new reaction %d before completing component molecule: "%s"'
                            % (iReaction, str.join("", componentText))
                        )
                else:
                    # No previous reaction, must be the first one.  Instantiate the object once
                    reactionMol = OEGraphMol()
                iReactionLine = 0

            elif iReactionLine == self.REACTANT_PRODUCT_COUNT_LINE:
                # Use this line to figure out how many reactants and products to look for
                chunks = line.split()
                nReactants = int(chunks[0])
                nProducts = int(chunks[1])

            elif componentText is not None:
                # Expecting more text to add to represent the molecule
                componentText.append(line)
                if line.startswith(self.MOLECULE_END_TAG):
                    # This was the last line of text representing the molecule
                    # Load it up as a molecule object
                    componentText = str.join("", componentText)
                    # Convert list of lines into a single multi-line string
                    oeis = virtual_oemolistream(componentText, OEFormat_SDF)
                    OEReadMolecule(oeis, componentMol)
                    oeis.close()

                    # Figure out what part of the reaction this component represents
                    reactionMol.SetRxn(True)
                    # If don't set, will get strange formatting for reaction SMILES strings
                    reactionRole = OERxnRole_None
                    if nReactants > 0:
                        reactionRole = OERxnRole_Reactant
                        nReactants -= 1
                    elif nProducts > 0:
                        reactionRole = OERxnRole_Product
                        nProducts -= 1
                    self.applyReactionRole(componentMol, reactionRole)

                    OEAddMols(reactionMol, componentMol)
                    componentMol.Clear()
                    componentText = None
                    # Reset to blank to indicate molecule finished

            elif incompleteName:
                # An annotation name field must have been started already, but looks like a long line spanning multi-line text.
                annotationName += lineWithoutNewline
                # Check if there are still more lines to look for
                incompleteName = len(
                    lineWithoutNewline
                ) > self.ANNOTATION_WRAP_COLUMN and annotationName.endswith(
                    self.ANNOTATION_WRAP_CHAR
                )
                if incompleteName:
                    # Strip off the wrap / join character
                    annotationName = annotationName[: -len(self.ANNOTATION_WRAP_CHAR)]

            elif incompleteValue:
                # An annotation Value field must have been started already, but looks like a long line multi-line text.
                # Typical multi-line text data
                annotationValue += lineWithoutNewline
                # Check if there are still more lines to look for
                incompleteValue = len(
                    lineWithoutNewline
                ) > self.ANNOTATION_WRAP_COLUMN and annotationValue.endswith(
                    self.ANNOTATION_WRAP_CHAR
                )
                if incompleteValue:
                    # Strip off the wrap / join character
                    annotationValue = annotationValue[: -len(self.ANNOTATION_WRAP_CHAR)]

            elif annotationMolText is not None:
                # Multi-line data representing formatted molecule text
                annotationMolText.append(line)
                if line.startswith(self.MOLECULE_END_TAG):
                    # This was the last line of text representing the molecule
                    # Load it up as a molecule object
                    annotationMolText = str.join("", annotationMolText)
                    # Convert list of lines into a single multi-line string
                    oeis = virtual_oemolistream(annotationMolText, OEFormat_SDF)
                    OEReadMolecule(oeis, annotationMol)
                    oeis.close()

                    # Alter the annotation value to instead just be a single-line notation (SMILES) string representation of the molecule
                    annotationValue = "%s %s" % (
                        self.ANNOTATION_SMI_TAG,
                        createStandardSmiString(annotationMol),
                    )

                    annotationMol.Clear()
                    annotationMolText = None
                    # Reset to blank to indicate molecule finished

            elif line.startswith(self.MOLECULE_BEGIN_TAG):
                # The subsequent lines are expected to represent a reaction component molecule
                if componentText is not None:
                    raise Exception(
                        'Starting new component molecule before completing previous one: "%s"'
                        % (str.join("", componentText))
                    )
                componentText = list()

            elif line.startswith(self.ANNOTATION_NAME_TAG):
                if annotationName is not None or incompleteName:
                    raise Exception(
                        'Starting new annotation before completing previous: "%s"'
                        % (annotationName)
                    )

                # Extract annotation name from the line by excising the tag header and stripping any flanking whitespace
                annotationName = line[len(self.ANNOTATION_NAME_TAG) :].strip()
                # Check if looks like a multi-line text value.  Next text line iteration should look to extend it if so
                incompleteName = len(
                    lineWithoutNewline
                ) > self.ANNOTATION_WRAP_COLUMN and annotationName.endswith(
                    self.ANNOTATION_WRAP_CHAR
                )
                if incompleteName:
                    # Strip off the wrap / join character
                    annotationName = annotationName[: -len(self.ANNOTATION_WRAP_CHAR)]

            elif line.startswith(self.ANNOTATION_VALUE_TAG):
                if annotationValue is not None or incompleteValue:
                    raise Exception(
                        'Starting new annotation value before completing previous: "%s"'
                        % (annotationValue)
                    )
                if annotationMolText is not None:
                    raise Exception(
                        'Starting new annotation before completing prior annotation molecule: "%s"'
                        % (str.join("", annotationMolText))
                    )

                # Extract annotation value from the line by excising the tag header and stripping any flanking whitespace
                annotationValue = line[len(self.ANNOTATION_VALUE_TAG) :].strip()
                # Check if looks like a multi-line text value.  Next text line iteration should look to extend it if so
                incompleteValue = len(
                    lineWithoutNewline
                ) > self.ANNOTATION_WRAP_COLUMN and annotationValue.endswith(
                    self.ANNOTATION_WRAP_CHAR
                )
                if incompleteValue:
                    # Strip off the wrap / join character
                    annotationValue = annotationValue[: -len(self.ANNOTATION_WRAP_CHAR)]

                # Special case of molecule formatted annotations, like sub-SDF files
                if annotationValue == self.ANNOTATION_MOL_TAG:
                    # Prepare to read subsequent lines in as molecule text
                    annotationMolText = list()

            if iReactionLine is not None:
                iReactionLine += 1

        # Check if there is one last annotation to add
        if annotationName is not None and annotationValue is not None:
            OESetSDData(reactionMol, annotationName, annotationValue)

            # Reset annotation fields
            annotationName = None
            annotationValue = None
            incompleteName = False
            incompleteValue = False

        # After final line, should have one last molecule to yield
        yield reactionMol
        reactionMol.Clear()
        # Blank for a new, fresh object without requiring another instantiation
        iReaction += 1

        # Verify that all tracking variables are cleared and didn't leave anything hanging
        if annotationName is not None or incompleteName:
            raise Exception(
                'Final reaction up to %d before completing annotation: "%s"'
                % (iReaction, annotationName)
            )
        if annotationValue is not None or incompleteValue:
            raise Exception(
                'Final reaction up to %d before completing annotation value: "%s"'
                % (iReaction, annotationValue)
            )
        if componentText is not None:
            raise Exception(
                'Final reaction up to %d before completing component molecule: "%s"'
                % (iReaction, str.join("", componentText))
            )

    def GetOEGraphMols(self):
        """Mirrors the OEChem oemolistream interface."""
        return iter(self)

    def close(self):
        """Mirrors the OEChem oemolistream interface."""
        self.inputFile.close()

    def applyReactionRole(mol, reactionRole):
        """Convenience method.  Set all of the atoms in the molecule
        as belonging to the respective reaction role
        (i.e., reactant, reagent, product).
        """
        for atom in mol.GetAtoms():
            atom.SetRxnRole(reactionRole)

    applyReactionRole = staticmethod(applyReactionRole)


class ReactionCanonizer:
    """Given a reaction representation, encoded through one of OpenEye's
    standard parsers (SMILES for SMIRKS, SDF for RDF), produce and canonize
    a SMIRKS string output for each.

    Primarily based on canonical SMILES algorithm, but extra steps to ensure
    consistent atom mapping numbers and hydrogen atom mapping.

    Known Issues:
    -   This actually only works well for "Reaction SMILES" not for full-fledged
        SMIRKS strings.  The difference is, SMIRKS strings can contain SMARTS
        pattern strings that are not directly manageable by OEChem's molecule objects.
    -   Symmetric molecules may not yield consistent SMILES when atom-mappings are used.
        If they're symmetric, usually it doesn't matter which branch is chosen,
        since the pieces are equivalent.  With atom-mapping however, such branches
        may be non-equivalent, but OEChem doesn't use atom-mapping to break these ties.
    """

    tempMol = None

    def __init__(self):
        """Constructor."""
        self.tempMol = OEGraphMol()

    def canonizeReactionMol(self, reactionMol):
        """Given an OEMolBase object representing a reaction,
        modify it's atom mapping contents and number to ensure
        canonical labeling.
        """
        OECanonicalOrderAtoms(reactionMol)

        oldToNewAtomMapping = dict()
        # If there were old / existing atom-mapping numbers, keep track of what they were replaced by
        for iAtom, atom in enumerate(reactionMol.GetAtoms()):
            standardMapIdx = iAtom + 1
            # +1 because 0 map index has special "null" meaning

            if atom.GetMapIdx() > 0:  # Only apply map numbers if started with one
                if atom.GetMapIdx() in oldToNewAtomMapping:
                    atom.SetMapIdx(oldToNewAtomMapping[atom.GetMapIdx()])
                else:
                    # Set a standard atom mapping number and record it
                    oldToNewAtomMapping[atom.GetMapIdx()] = standardMapIdx
                    atom.SetMapIdx(standardMapIdx)

                if atom.GetAtomicNum() == 1 and atom.GetIsotope() < 1:
                    # Hydrogen special handling, or else atom-mapping won't output
                    atom.SetIsotope(1)

        return reactionMol

    def canonizeReactionSmiles(self, reactionSmi):
        """Just create a molecule object representation of the SMILES string
        and call canonizeReactionMol on it, then return a SMILES representation back.
        """
        OEParseSmiles(self.tempMol, reactionSmi)
        self.tempMol = self.canonizeReactionMol(self.tempMol)
        canonizedSmi = createStandardSmiString(self.tempMol)
        self.tempMol.Clear()

        # Correct wildcard atoms (gets translated into "R-groups," but SMIRKS doesn't accept this format)
        canonizedSmi = canonizedSmi.replace("R", "*:")

        return canonizedSmi

    def runByMolStream(self, inputMolStream, outputMolStream):
        """Read through contents of the inputFile, canonize each reactionMol
        and then write to the outputFile.
        """
        for reactionMol in inputMolStream.GetOEGraphMols():
            self.canonizeReactionMol(reactionMol)
            OEWriteMolecule(outputMolStream, reactionMol)

    def runByFilename(self, inputFilename, outputFilename):
        """Read through contents of the inputFile, canonize each reactionMol
        and then write to the outputfile.
        """
        ifs = oemolistream(inputFilename)
        ofs = oemolostream(outputFilename)

        self.runByMolStream(ifs, ofs)


"""Static instance that should be good enough for normal users.
Beware that this is NOT thread-safe by default because of the
instance variable "tempMol."
"""
aReactionCanonizer = ReactionCanonizer()

if __name__ == "__main__":
    main(sys.argv)
