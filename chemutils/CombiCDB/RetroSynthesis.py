import sys
from optparse import OptionParser
from openeye.oechem import OECreateIsoSmiString, OESmilesAtomCount

from chemutils.Common import DBUtil
from chemutils.Common.Util import (
    ProgressDots,
    molBySmiles,
    stdOpen,
)
from chemutils.Common.Util import splitCompositeMolToSmilesList
from chemutils.Common.Util import createStandardSmiString
from chemutils.CombiCDB.ReactionModel import SupplementalDataModel
from chemutils.CombiCDB.ReactionModel import (
    PRODUCTS,
    reactionStepStr,
)
from chemutils.CombiCDB.ReagentManager import ReagentManager, ReagentSearchModel
from chemutils.CombiCDB.SynthesisGenerator import SynthesisGenerator
from chemutils.CombiCDB.SynthesisUtil import SimpleReactantPool, ChemDBReactantPool
from chemutils.CombiCDB.SynthesisUtil import (
    GreedyRetroProductSelector,
    ExhaustiveRetroProductSelector,
)

from chemutils.CombiCDB.Const import SIMILARITY_THRESHOLD


class RetroSynthesis:
    """Initialized with a set of reagent objects, generates a retro version
    of each reagent.  Also needs a link to some pool of reactants including
    a similarity measure / score between intermediate products found to items
    in the pool.

    Given a desired target structure, will attempt to find a path from starting
    components in the reactant pool to the target via a series of reagent reactions,
    by searching backwards.

    General algorithm:
        While the target mol(s) are not found in the reactant pool (closest item has < 1.0 similarity)
        and we haven't exceeded the maximum number of steps for consideration:
            Apply each retro reagent to the target mol and see if it yields any retro products (source reactants).
            Pick one or more of these retro steps based on how some scoring function of the retro products.
            Recurse with the retro products being target molecules themselves.

    Limitations:
        - Assumes all retro-reactions expect a single (non-inherent) reactant
        - If no exact match is found, just returns the last one tried, rather than sum integration over all or even a ranked list
    """

    # Maximum number of backwards reaction steps to consider
    # Final number of steps may be greater than this if the pathway has branch points,
    #   this limits the depth of the pathway tree.
    maxDepth = 0

    reagentList = None
    retroReagentList = None
    reactantPool = None
    retroProductSelector = None

    def __init__(self, reagentList, reactantPool, retroProductSelector):
        """Constructor"""
        self.reactantPool = reactantPool
        self.retroProductSelector = retroProductSelector
        self.reagentList = reagentList
        self.retroReagentList = []
        for reagent in reagentList:
            retroReagent = reagent.retro()
            retroReagent.countIntermediates = True
            # Reactant may have already been half-processed (e.g., aldehyde, not alcohol, oxidized to carboxylate)
            retroReagent.removeInherentProducts = True
            # Not interested in these
            self.retroReagentList.append(retroReagent)

    def __call__(self, targetMol, depth=0):
        """Generate a series of reaction steps to produce the targetMol.
        The starting reactants for the steps should be as similar as possible
        (if not exactly matching) the reactants in the asssociated reactantPool.
        Return a 2-ple, the first item being the reaction steps and the second
        being a boolean as to whether an exact match was successfully found.
        """

        targetSmi = OECreateIsoSmiString(targetMol)

        print >> sys.stderr, targetSmi, depth

        # Check if the targetMol is already in reactant pool
        if targetMol in self.reactantPool:
            return ([], True)
            # Base case, return blank synthesis pathway, target is already in the reactant pool

        # Recursion base case
        if depth >= self.maxDepth:
            return ([], False)

        # Try applying each retro reagent to the target mol and select a subset of them to pursue
        # For the time being, always assume the retroReagents expect 1 molecule
        # Expect list of 3-ples: (retroProducts, retroReagent, iReagent)
        retroProductsByReagent = self.retroProductSelector(
            self.retroReagentList, targetMol
        )

        reactionSteps = []
        exactMatch = False

        for retroProduct, retroReagent, iReagent in retroProductsByReagent:
            # Recursively find synthesis pathway for retroProduct(s) (if any)
            exactMatch = True
            reactionSteps = []
            componentProductSmiList = splitCompositeMolToSmilesList(retroProduct)
            for componentProductSmi in componentProductSmiList:
                componentMol = molBySmiles(componentProductSmi)
                (precursorSteps, precursorExactMatch) = self(componentMol, depth + 1)
                # Recursive call
                reactionSteps.extend(precursorSteps)
                exactMatch = exactMatch and precursorExactMatch
            # Add to the end the step we just used to create the retro prorduct in the first place
            reactionSteps.append(
                [componentProductSmiList, self.reagentList[iReagent], [targetSmi]]
            )

            if exactMatch:
                # Found an exact match, successful search, return the result
                return (reactionSteps, exactMatch)
            else:
                # Not an exact match, dump this pathway and keep looking among the other possibilities
                pass
        # If reached this point, means no exact solution was found, just return the last one tried
        return (reactionSteps, exactMatch)

    def main(argv):
        """Main method, callable from command line"""
        usageStr = """usage: %prog [options] <productSmiles>
            
            Looks for synthesis pathways to generate the productSmiles based on the
            known database of starting material reactants and reagents that can apply reactions to them.
            
            productSmiles:  SMILES string of the target product molecule to find a synthesis pathway for.
                            If omitted, will auto-generate a random one with the SynthesisGenerator
            """
        parser = OptionParser(usage=usageStr)
        parser.add_option(
            "-d",
            "--maxDepth",
            dest="maxDepth",
            metavar="<maxDepth>",
            default="4",
            help="Limit the retro-synthesis search depth to only this many reaction steps.",
        )
        parser.add_option(
            "-s",
            "--synthSteps",
            dest="synthSteps",
            metavar="<synthSteps>",
            default="4",
            help="If no product is specified, will auto-generate a random one with the SynthesisGenerator, using up to this many reaction steps",
        )
        parser.add_option(
            "-i",
            "--iterations",
            dest="iterations",
            metavar="<iterations>",
            default="1",
            help="If no product is specified, will test by auto-generating this many random products.",
        )
        parser.add_option(
            "-C",
            "--ChemDB",
            dest="ChemDB",
            action="store_true",
            help="If set, will use the entire linked ChemDB as the available reactant pool instead of just the sample test reactant pool.",
        )
        parser.add_option(
            "-g",
            "--greedyChoices",
            dest="greedyChoices",
            metavar="<greedyChoices>",
            default="1",
            help="When doing a greedy search, number of paths to consider at any branch point.  Set to 0 to do an exhaustive search.",
        )
        parser.add_option(
            "-o",
            "--outputFile",
            dest="outputFile",
            metavar="<outputFile>",
            default="-",
            help="Name of file to output results to.  Defaults to stdout",
        )
        (options, args) = parser.parse_args(argv[1:])

        # Load list of test reactants to work with
        reactantSmilesList = []
        reactantTable = DBUtil.execute("select smiles from reactant")
        for row in reactantTable:
            reactantSmilesList.append(row[0])

        # Load known reagents
        reagentManager = ReagentManager()
        reagentQuery = ReagentSearchModel()
        reagentQuery.enabled = True
        reagentList = reagentManager.loadDBInstances(reagentQuery)

        # Check if product is specified
        origProductSmiles = None
        if len(args) > 0:
            origProductSmiles = args[0]

        # Prepare product generator if doing a random test
        generator = SynthesisGenerator(options.synthSteps)

        # Choose which reactant pool to use for "available starting materials"
        reactantPool = None
        if options.ChemDB:
            reactantPool = ChemDBReactantPool(SIMILARITY_THRESHOLD)
        else:
            reactantPool = SimpleReactantPool(reactantSmilesList)

        # Select which product selector will be used to guide the retro-synthetic search
        productSelector = None
        greedyChoices = int(options.greedyChoices)
        if greedyChoices > 0:
            productSelector = GreedyRetroProductSelector(
                reactantPool, AverageScoreAggregator(), greedyChoices
            )
        else:
            productSelector = ExhaustiveRetroProductSelector()

        # Output file to write results to
        outfile = stdOpen(options.outputFile, "w", sys.stdout)

        # Generate an actual instance to use
        retroSynth = RetroSynthesis(reagentList, reactantPool, productSelector)
        retroSynth.maxDepth = int(options.maxDepth)

        iterations = 1
        if origProductSmiles is None:
            iterations = int(options.iterations)

        nSuccess = 0
        progress = ProgressDots(1, 1)
        for i in xrange(iterations):
            productSmiles = origProductSmiles
            if productSmiles is None:
                # No product was specified, randomly generate one here
                reactionSteps = generator(reactantSmilesList, reagentList)
                finalProductSmilesList = reactionSteps[-1][PRODUCTS]
                maxProductAtoms = 0
                finalProductSmiles = None
                for productSmiles in finalProductSmilesList:
                    productAtoms = OESmilesAtomCount(productSmiles)
                    if productAtoms > maxProductAtoms:
                        maxProductAtoms = productAtoms
                        finalProductSmiles = productSmiles
                productSmiles = finalProductSmiles
            productMol = molBySmiles(productSmiles)

            print >> outfile, i, productSmiles
            (retroReactionSteps, exactMatch) = retroSynth(productMol)
            print >> outfile, exactMatch

            for sampleStep in retroReactionSteps:
                print >> outfile, reactionStepStr(sampleStep)

            if exactMatch:
                nSuccess += 1

            progress.Update()

        progress.PrintStatus()
        print >> sys.stderr
        print >> sys.stderr, nSuccess, "Successes"

    main = staticmethod(main)


class ReactionPlanner:
    """Looks for a pathway to produce a target product using a
    limited available reactant and reagent pool.
    Further limited by the maximum number of steps that can be used in the synthesis.

    Similar purpose to retro-synthesis, but searching is done in the forward direction,
    by considering all possible combinatorial reaction pathways from the starting materials
    until the target product is found, or all avenues are exhausted.

    Instead of an exhaustive traversal of all possibilities, could also use some kind of
    similarity heuristics to guide the search, but right now this is just for
    reconstructing apparent reaction tutorial problem solutions which should have very
    limited pool of starting materials (only the exact ones used) anyway.
    """

    def __init__(self, reactionCategoryIDs=None):
        # Instantiate a SynthesisGenerator object to reuse some screening functions
        self.synthesisGenerator = SynthesisGenerator()
        self.reactionCategoryIdSet = None
        if reactionCategoryIDs is not None:
            # Allow reaction category filters to help restrict the set of reagent usage that is appropriate in the plan
            self.reactionCategoryIdSet = set(reactionCategoryIDs)

    def __call__(self, startingMaterialList, reagentList, targetMol, maxSteps):
        """Return an iterator (generator) over series of reaction steps to produce the
        targetMol from the available starting material and reagent list, using up to
        the specified number of max steps.

        Beware, this yields an iterator over the SAME list representing the reaction steps,
        so after each is yielded, the previous contents are discarded.
        Make copies of the yielded reactionSteps list objects if you want to retain them for longer.
        """
        return self.addReactionSteps(
            startingMaterialList, reagentList, targetMol, maxSteps
        )

    def addReactionSteps(
        self,
        reactantList,
        reagentList,
        targetMol,
        maxSteps,
        reactionSteps=None,
        usedReactantSmiList=None,
    ):
        """Primary recursive function to drive the general call method.
        Progressively add to the reactionSteps list with all combinations from the available
        reactant and reagent lists.  If eventually encounter the targetMol, yield the accumulated
        reaction steps that produce the target.
        """
        targetSmi = createStandardSmiString(targetMol)
        reactantSmiSet = set()
        for reactant in reactantSmiSet:
            reactantSmiSet.add(createStandardSmiString(reactant))

        if reactionSteps is None:
            reactionSteps = list()
            # Check if the targetMol is already in reactant pool, nothing really to do then
            if targetSmi in reactantSmiSet:
                yield reactionSteps
                return

        if usedReactantSmiList is None:
            # Keep track of which reactants have been used over the course of the current steps
            #   to ensure original reactants are not re-used multiple times
            usedReactantSmiList = list()

        # Recursion base case
        if len(reactionSteps) >= maxSteps:
            # Reached or exceeded the maximum number of steps, apparently without having
            #   found the intended target product.  Return to indicate failure.
            return

        # List of molecules we should consider as potential (required) reactants for the current step
        keyReactantList = list()
        if len(reactionSteps) < 1:
            # No prior steps to take products from, just pull from the original reactant list
            keyReactantList.extend(reactantList)
        else:
            # Use carryover products from the last reaction step
            carryoverProductList = reactionSteps[-1][PRODUCTS]
            keyReactantList.extend(
                SynthesisGenerator.selectBestProducts(carryoverProductList)
            )

        # Try starting from each key reactant
        for keyReactant in keyReactantList:
            keyReactantSmi = createStandardSmiString(keyReactant)

            currentReactants = list()
            currentReactants.append(keyReactant)
            usedReactantSmiList.append(keyReactantSmi)

            # Try each reagent on the current reactant set
            for reagent in reagentList:
                for matchSteps in self.applyReagentStep(
                    currentReactants,
                    reagent,
                    reactantList,
                    reagentList,
                    targetMol,
                    maxSteps,
                    reactionSteps,
                    usedReactantSmiList,
                ):
                    yield matchSteps

                if reagent["expected_reactants"] > 1:
                    # No matches or maybe no products generated when used lone reactant,
                    #   but this reagent may use multiple reactants.
                    # Try adding combos from the original reactant list.
                    for originalReactant in reactantList:
                        originalReactantSmi = createStandardSmiString(originalReactant)

                        if originalReactantSmi in usedReactantSmiList:
                            # Avoid reusing original reactants
                            continue

                        currentReactants.append(originalReactant)
                        usedReactantSmiList.append(originalReactantSmi)

                        for matchSteps in self.applyReagentStep(
                            currentReactants,
                            reagent,
                            reactantList,
                            reagentList,
                            targetMol,
                            maxSteps,
                            reactionSteps,
                            usedReactantSmiList,
                        ):
                            yield matchSteps

                        # Revert to previous state
                        currentReactants.pop()
                        usedReactantSmiList.pop()

            # Revert to previous state
            # currentReactants.pop();   # Not necessary because start with new list each time.
            #                           # Don't want to modify since this could be part of a previously yield result object.
            usedReactantSmiList.pop()

    def applyReagentStep(
        self,
        currentReactants,
        reagent,
        reactantList,
        reagentList,
        targetMol,
        maxSteps,
        reactionSteps,
        usedReactantSmiList,
    ):
        """Apply the reagent model to the current reactants.
        Check if this generates the targetMol, then yield
        the completed reactionSteps synthesis list.
        Otherwise (assuming the reagent generates anything at all)
        recursively call the addReactionSteps function.
        """
        reagent.removeInherentProducts = True
        reagent.ignoreSelfReactions = True
        supplementalData = SupplementalDataModel()
        currentProducts = reagent(currentReactants, supplementalData)

        acceptableStep = False

        if len(currentProducts) > 0:
            acceptableStep = self.synthesisGenerator.acceptableReactionStep(
                currentReactants,
                reagent,
                currentProducts,
                supplementalData,
                self.reactionCategoryIdSet,
            )

            if acceptableStep:
                # Make a copy of the reactant list or else it can be modified by future iterations
                currentStep = (list(currentReactants), reagent, currentProducts)
                reactionSteps.append(currentStep)

                # print >> sys.stderr, reactionStepStr(currentStep), len(reactionSteps);

                # See if any of the generated products match the target
                targetSmi = createStandardSmiString(targetMol)
                for product in currentProducts:
                    productSmi = createStandardSmiString(product)
                    if productSmi == targetSmi:
                        # Found a match return the steps used to indicate success
                        yield reactionSteps

                # No matches, but maybe just needs more steps
                for subMatch in self.addReactionSteps(
                    reactantList,
                    reagentList,
                    targetMol,
                    maxSteps,
                    reactionSteps,
                    usedReactantSmiList,
                ):
                    yield subMatch

                # Revert to previous state
                reactionSteps.pop()

        # print reactionStepStr( (currentReactants, reagent, currentProducts) ), acceptableStep


if __name__ == "__main__":
    RetroSynthesis.main(sys.argv)
