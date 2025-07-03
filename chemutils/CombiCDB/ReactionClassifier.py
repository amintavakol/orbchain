import sys
from optparse import OptionParser
from chemutils.Common.Util import virtual_oemolistream, createStandardSmiString
from openeye.oechem import (
    oemolistream,
)
from openeye.oechem import OEGraphMol, OEParseSmiles
from openeye.oechem import OERxnRole_None, OERxnRole_Reactant, OERxnRole_Product
from chemutils.CombiCDB.ReactionProcessor import ReactionProcessor
from chemutils.CombiCDB.ReactionFormatter import RDFReader
from chemutils.Common.Util import *
from openeye.oechem import *





def main(argv):
    """Main method, callable from command line"""

    usageStr = (
        "usage: %prog [options] <inputFile> <outputFile> \n"
        + "   <reactionDatabaseFile>     File containing reaction information, such as SMIRKS / SMILES or RDF / SDF\n"
        + "   <reactionPatternFile>    File containing SMIRKS reaction patterns for different name reaction.\n"
    )

    parser = OptionParser(usage=usageStr)
    (options, args) = parser.parse_args(argv[1:])

    if len(args) >= 2:
        reactionDatabaseFile = args[0]
        reactionPatternFile = args[1]

        reactionDatabaseStream = reactionDatabaseistream(reactionDatabaseFile)
        # Read records from a reaction Database
        for molReaction in reactionDatabaseStream.GetOEGraphMols():
            molReactionString = OECreateIsoSmiString(molReaction)
            # reactionPatternStream =reactionPatternistream(reactionPatternFile);
            # for molPattern in reactionPatternStream.GetOEGraphMols():
            ## for reaction patterns

            print(molReactionString)

            rpf = open(reactionPatternFile)
            for linePattern in rpf:
                molPatternString, molPatternName = linePattern.split(" ", 1)
                # judge if two reactions are the same or not
                reactionMatcher = ReactionMatcher()
                # molPatternt = OECreateIsoSmiString(molPattern);
                print("run by one pattern")
                if (
                    reactionMatcher.runByReaction(molReactionString, molPatternString)
                    == 1
                ):
                    # log.debug("-------------------------------reactions matched;-)");
                    print("reaction matched by one pattern")
                    print(molPatternString)
                else:
                    rpfSecond = open(reactionPatternFile)
                    # rpfSecond.close();
                    print("run two patterns")
                    for linePatternSecond in rpfSecond:
                        (
                            molPatternStringSecond,
                            molPatternNameSecond,
                        ) = linePatternSecond.split(" ", 1)
                        # judge if two reactions are the same or not
                        if (
                            reactionMatcher.runByTwoReactionPatterns(
                                molReactionString,
                                molPatternString,
                                molPatternStringSecond,
                            )
                            == True
                        ):
                            # log.debug("--------------reactions matched by Two patterns;-)");
                            print("reaction matched by Two patterns")
                        else:
                            rpfThird = open(reactionPatternFile)
                            # rpfThird.close();
                            print("run three patterns")

                            for linePatternThird in rpfThird:
                                (
                                    molPatternStringThird,
                                    molPatternNameThird,
                                ) = linePatternThird.split(" ", 1)
                                if (
                                    reactionMatcher.runByThreeReactionPatterns(
                                        molReactionString,
                                        molPatternString,
                                        molPatternStringSecond,
                                        molPatternStringThird,
                                    )
                                    == True
                                ):
                                    print("reaction matched by Three patterns")
                            rpfThird.close()
                    rpfSecond.close()

            rpf.close()

        reactionDatabaseStream.close()
    else:
        parser.print_help()
        sys.exit(-1)


class ReactionMatcher:
    """
    Given a reaction (rdf, sdf, smi) and a reaction pattern (smirks), compare if they are the same or not.
    """

    def __init__(self):
        """Constructor."""
        pass

    def runByReaction(self, molReaction, molPattern):
        """
        Compare two records from reaction inputFile and pattern inputFile. To output yes(1) or no(0)
        """
        # extract reactants
        # print( "      reading Reactants   ")
        reactantMolList = extractReactants(molReaction)

        # generate products
        generatedProductOEOS = virtual_oemolostream(OEFormat_ISM)
        # print("      generating Products  ");
        generatedProductOEOS = generateProducts(molPattern, reactantMolList)
        # Read output back through a virtual input stream
        generatedProductIS = virtual_oemolistream(
            generatedProductOEOS.GetString(), OEFormat_ISM
        )
        # extract products

        # for generatedProductMol in generatedProductIS.GetOEGraphMols():
        #    print "run by one pattern "+createStandardSmiString(generatedProductMol)
        # print "run by one pattern "
        expectedProductMol = extractProducts(molReaction)

        # compare extracted and generated products
        # print("      comparing Products   ")
        isMatched = compareProducts(
            expectedProductMol, generatedProductIS.GetOEGraphMols(), reactantMolList
        )
        return isMatched

    def runByTwoReactionPatterns(self, molReaction, molPattern1, molPattern2):
        """
        Match a reaction by two patterns"""
        isMatched = False
        reactantMolList1 = extractReactants(molReaction)
        generatedProductOEOS1 = generateProducts(molPattern1, reactantMolList1)
        generatedProductIS1 = virtual_oemolistream(
            generatedProductOEOS1.GetString(), OEFormat_ISM
        )
        expectedProductMol = extractProducts(molReaction)

        for generatedProductMol1 in generatedProductIS1.GetOEGraphMols():
            canonizedGeneratedSmi = OECreateCanSmiString(generatedProductMol1)
            molReaction2 = canonizedGeneratedSmi + "." + molReaction
            reactantMolList2 = extractReactants(molReaction2)

            libgen2 = OELibraryGen(molPattern2, True)
            productOEOS2 = virtual_oemolostream(OEFormat_ISM)

            reactionProcessor = ReactionProcessor()
            reactionProcessor.applyReaction(libgen2, reactantMolList2, productOEOS2, 0)

            # generatedProductOEOS2 = generateProducts(molPattern2,reactantMolList2);
            generatedProductIS2 = virtual_oemolistream(
                productOEOS2.GetString(), OEFormat_ISM
            )
            # print "run by two patterns "
            # for generatedProductMol2 in generatedProductIS2.GetOEGraphMols():
            #    print "run by two patterns "+createStandardSmiString(generatedProductMol2)
            if (
                compareProducts(
                    expectedProductMol,
                    generatedProductIS2.GetOEGraphMols(),
                    reactantMolList2,
                )
                == True
            ):
                isMatched = True
        return isMatched

    def runByThreeReactionPatterns(
        self, molReaction, molPattern1, molPattern2, molPattern3
    ):
        """
        Match a reaction by three patterns. To output yes(1) or no(0)
        """
        # extact reactants
        isMatched = False
        reactantMolList1 = extractReactants(molReaction)
        generatedProductOEOS1 = generateProducts(molPattern1, reactantMolList1)
        generatedProductIS1 = virtual_oemolistream(
            generatedProductOEOS1.GetString(), OEFormat_ISM
        )
        expectedProductMol = extractProducts(molReaction)

        for generatedProductMol1 in generatedProductIS1.GetOEGraphMols():
            canonizedGeneratedSmi = OECreateCanSmiString(generatedProductMol1)
            molReaction2 = canonizedGeneratedSmi + "." + molReaction
            reactantMolList2 = extractReactants(molReaction2)
            generatedProductOEOS2 = generateProducts(molPattern2, reactantMolList2)
            generatedProductIS2 = virtual_oemolistream(
                generatedProductOEOS2.GetString(), OEFormat_ISM
            )

            for generatedProductMol2 in generatedProductIS2.GetOEGraphMols():
                canonizedGeneratedSmi2 = OECreateCanSmiString(generatedProductMol2)
                # print "canonizedGeneratedSmi2::: "+canonizedGeneratedSmi2
                molReaction3 = canonizedGeneratedSmi2 + "." + molReaction
                reactantMolList3 = extractReactants(molReaction3)
                generatedProductOEOS3 = generateProducts(molPattern3, reactantMolList3)
                generatedProductIS3 = virtual_oemolistream(
                    generatedProductOEOS3.GetString(), OEFormat_ISM
                )
                # for generatedProductMol3 in generatedProductIS3.GetOEGraphMols():
                # print "run by three patterns "+createStandardSmiString(generatedProductMol3)
                # print "run by three patterns "
                if (
                    compareProducts(
                        expectedProductMol,
                        generatedProductIS3.GetOEGraphMols(),
                        reactantMolList3,
                    )
                    == True
                ):
                    isMatched = True
        return isMatched

    def runByReagent(self, reactionSmi, reagent):
        """
        Return whether the reaction SMILES fits under the reagent model.
        That is, will applying the reagent model to the reactant side of
        the reaction equation yield the same product side?
        """
        reactantMolList = extractReactants(reactionSmi)
        try:
            generatedProductMolList = reagent(reactantMolList)
            # print "generatedProductMolList ", generatedProductMolList

            expectedProductMol = extractProducts(reactionSmi)
            # print "expectedProductMol ",expectedProductMol

            matchFound = compareProducts(
                expectedProductMol, generatedProductMolList, reactantMolList
            )

            return matchFound
        except:
            pass
        return False


# extract reactant molecules from a reaction record
def extractReactants(molReaction):
    # print "molReaction"+molReaction
    reactionMol = OEGraphMol()
    OEParseSmiles(reactionMol, molReaction)
    # Create copy of the whole reaction, then delete all non-reactants
    reactantsMol = OEGraphMol(reactionMol)
    for atom in reactantsMol.GetAtoms():
        if atom.GetRxnRole() != OERxnRole_Reactant:
            reactantsMol.DeleteAtom(atom)
        else:
            # Change "reactants" to just ordinary molecules
            atom.SetRxnRole(OERxnRole_None)
    reactantsMol.SetRxn(False)
    # Treat as "regular" composite molecule now, not a reaction
    return reactantsMol


# generate product molecules by using reactants and one reation pattern
def generateProducts(smirksPattern, reactantsMol):
    reactantMolList = splitCompositeMol(reactantsMol)
    libgen = OELibraryGen(smirksPattern, True)
    productOEOS = virtual_oemolostream(OEFormat_ISM)

    reactionProcessor = ReactionProcessor()
    reactionProcessor.applyReaction(libgen, reactantMolList, productOEOS, 2)
    productIS = virtual_oemolistream(productOEOS.GetString(), OEFormat_ISM)
    return productOEOS


# extract product molecules from a reaction record
def extractProducts(molReaction):
    reactionMol = OEGraphMol()
    OEParseSmiles(reactionMol, molReaction)
    productsMol = OEGraphMol(reactionMol)
    for atom in productsMol.GetAtoms():
        if atom.GetRxnRole() != OERxnRole_Product:
            productsMol.DeleteAtom(atom)
        else:
            # Change "products" to just ordinary molecules
            atom.SetRxnRole(OERxnRole_None)
    productsMol.SetRxn(False)
    # Treat as "regular" composite molecule now, not a reaction

    return productsMol


#
def compareProducts(expectedProductMol, generatedProductIter, reactantMolList):
    """
    Go through every product in the generated product iterator
    and return True if any of them match the expected product.
    Match determined by standard (isomeric) SMILES representation.
    """

    OEAssignImplicitHydrogens(expectedProductMol)
    OEAddExplicitHydrogens(expectedProductMol)
    OEAssignAromaticFlags(expectedProductMol)
    OEAssignMDLHydrogens(expectedProductMol)

    canonizedExpectedSmi = createStandardSmiString(expectedProductMol)
    # print "expected "+canonizedExpectedSmi

    OEAssignImplicitHydrogens(reactantMolList)
    OEAddExplicitHydrogens(reactantMolList)
    OEAssignAromaticFlags(reactantMolList)
    OEAssignMDLHydrogens(reactantMolList)
    reactantsMol = createStandardSmiString(reactantMolList)
    # print tt

    for generatedProductMol in generatedProductIter:
        try:
            OEAssignImplicitHydrogens(generatedProductMol)
            OEAddExplicitHydrogens(generatedProductMol)
            OEAssignAromaticFlags(generatedProductMol)
            # OESuppressHydrogens(generatedProductMol)
            OEAssignMDLHydrogens(generatedProductMol)
            canonizedGeneratedSmi = OECreateCanSmiString(generatedProductMol)

            # print "generateded "+canonizedGeneratedSmi
            # to get individual generated product molecule
            countGenerated = canonizedGeneratedSmi.count(".")
            countExpected = canonizedExpectedSmi.count(".")

            for canonizedGeneratedSmiIndividual in canonizedGeneratedSmi.split("."):
                # print "canonizedGeneratedSmiIndividual"+canonizedGeneratedSmiIndividual
                for canonizedExpectedSmiIndividual in canonizedExpectedSmi.split("."):
                    # print "canonizedExpectedSmiIndividual"+canonizedExpectedSmiIndividual
                    # and compare with individual extracted product molecule
                    # reactionMoltt = OEGraphMol();
                    # reactionMol(reactantMolList)
                    # print OECreateCanSmiString(reactantMolList)
                    if (
                        inReactantsMolList(
                            reactantsMol, canonizedGeneratedSmiIndividual
                        )
                        == False
                    ):
                        if (
                            canonizedGeneratedSmiIndividual
                            == canonizedExpectedSmiIndividual
                        ):
                            # print canonizedGeneratedSmiIndividual
                            # print "canonizedGeneratedSmiIndividual"
                            return True
        except:
            pass
    return False


def inReactantsMolList(reactantsMol, canonizedGeneratedSmiIndividual):
    countReactantMol = reactantsMol.count(".")
    for reactantMol in reactantsMol.split(".", countReactantMol + 1):
        # print "reactntMol                     "+reactantMol
        # print "canonizedGeneratedSmiIndividual"+canonizedGeneratedSmiIndividual
        if canonizedGeneratedSmiIndividual == reactantMol:
            # print "true"
            return True
    # print "false"
    return False


def reactionDatabaseistream(reactionDatabaseFile):
    """read reaction records"""
    reactionDatabaseStream = None
    lowerreactionDatabaseFile = reactionDatabaseFile.lower()
    if lowerreactionDatabaseFile.endswith(".rdf") or lowerreactionDatabaseFile.endswith(
        ".rdf.gz"
    ):
        # Looks like an RDF input file.  Need specialized reader.
        # Standard OEChem oemolistream doesn't handle as would like
        if lowerreactionDatabaseFile.endswith(".gz"):
            # Looks like it's gzipped, open it accordingly
            import gzip

            reactionDatabaseStream = RDFReader(gzip.open(reactionDatabaseFile))
        else:
            # "Normal" RDF file, use the RDFReader to encase it in an oemolistream equivalent
            reactionDatabaseStream = RDFReader(open(reactionDatabaseFile))
    else:
        # Not an RDF file, just use standard oemolistream
        return oemolistream(reactionDatabaseFile)
    # reactionDatabaseStream = oemolistream(reactionDatabaseFile);
    # return reactionDatabaseStream;
    return reactionDatabaseStream


def reactionPatternistream(reactionPatternFile):
    """read reaction patterns"""
    reactionPatternStream = oemolistream(reactionPatternFile)
    return reactionPatternStream


if __name__ == "__main__":
    main(sys.argv)
