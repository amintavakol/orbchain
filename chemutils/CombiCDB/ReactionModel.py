"""Common objects / base classes used to support reaction processing. 

Primarily revolves around Reagent object models to predict
products from input reactants.  The example code below illustrates
how to load existing reagents from the database using the
ReagentSearchModel and ReagentManager.loadDBInstances method.
Once these are ready, and any settings are made (such as
the removal of inherent products), the reagent acts as a functor
which can take in a list of reactant molecule objects and
generate predicted product molecule objects.


>>> from CHEM.Common.Util import molBySmiles;
>>> from CHEM.CombiCDB.ReagentManager import ReagentManager, ReagentSearchModel;
>>> from CHEM.CombiCDB.ReactionModel import SMIRKSReagent;
>>> from CHEM.CombiCDB.ReactionModel import SupplementalDataModel, reactionStepStr;
>>> from CHEM.CombiCDB.test.ReagentIdentifiers import HYDROBROMINATION;
>>> reagentQuery = ReagentSearchModel();
>>> reagentQuery.reagentIDs = [HYDROBROMINATION];
>>> manager = ReagentManager();
>>> reagents = manager.loadDBInstances( reagentQuery );
>>> reagent = reagents[0];  # Expect only the one hydrobromination reagent
>>> reagent.removeInherentProducts = True;  # Not interested in inherent products
>>> print reagent["reagent_id"], reagent["description"];    # Reagent works as a dictionary
1 Hydrobromination
>>> supplementalData = SupplementalDataModel(); # Optional object to receive extra information
>>> reactants = [molBySmiles("CC=C")];  # List of molecule objects as reactant input
>>> products = reagent( reactants, supplementalData );    # Use reagent as function to generate products
>>> reactionStep = ( reactants, reagent, products );    # Combined elements of a complete reaction
>>> print reactionStepStr( reactionStep ).replace("\\t"," ");   # Tab characters mess up doctest
CC=C>>CC(C)Br Hydrobromination
>>> for elementaryStep in supplementalData.filteredElementarySteps( reactants, products ):
...     (detailReactants, reactionProfile, detailProducts) = elementaryStep;
...     print "%(reaction_profile_id)s: %(description)s" % reactionProfile;
47: Alkene, Protic Acid Addition, Secondary
1414: Carbocation, Halide Addition

Alternative example showing how to more easily retrieve and recycle reagent objects
without always having to query from the database (a relatively expensive operation).
This only works if you know the exact reagentIDs you're interested in.
>>> manager.objectCache = dict();   # This will probably be a persistent dictionary you keep somewhere to cache the objects
>>> neededReagentIds = [HYDROBROMINATION];
>>> reagentsById = manager.retrieveReagentsById( neededReagentIds );
>>> reagent = reagentsById[HYDROBROMINATION];
>>> print reagent["reagent_id"], reagent["description"];
1 Hydrobromination
>>> manager.returnReagentsById( reagentsById );

"""
# from sets import Set;
from threading import Lock
from openeye.oechem import OEGraphMol, OEParseSmiles, OEAddMols
from openeye.oechem import OELibraryGen
from openeye.oechem import OENetCharge
from openeye.oechem import OEFormat_ISM
from openeye.oechem import OEAssignAromaticFlags

from chemutils.Common.Const import SMILES_MOL_DELIM, NULL_STRING, EMPTY_SET
from chemutils.Common.Util import (
    virtual_oemolistream,
    virtual_oemolostream,
    molBySmiles,
)
from chemutils.Common.Util import createStandardSmiString
from chemutils.Common.Util import (
    splitCompositeMol,
    splitCompositeMolToSmilesList,
    splitCompositeSmilesToList,
    molToSmiList,
    atomCount,
)
from chemutils.Common.Util import Codec
from chemutils.Common.Model import (
    RowItemModel,
)
from chemutils.Common.MolExt import (
    enumerateRacemicMixture,
    enumerateRacemicMixtureList,
    createStandardMol,
    isCarbonRadical,
)
from ReactionProcessor import ReactionProcessor
from Util import log
from chemutils.CombiCDB.Const import (
    REACTION_DELIM,
    MOL_LIST_DELIM,
    REACTION_PROFILE_DELIM,
    VALUE_DELIM,
)
from chemutils.CombiCDB.Const import (
    SENTINEL_CHARGE,
    SENTINEL_OPTIONAL_STEP_CHARGE,
    NORMAL_CHARGE_THRESHOLD,
    SENTINEL_REJECT_CHARGE,
    SENTINEL_REJECT_IMMEDIATE_CHARGE,
)
from chemutils.CombiCDB.Const import SENTINEL_CONJUGATE_BASE, SENTINEL_CONJUGATE_ACID
from chemutils.CombiCDB.Const import REQUIRED_ATOMIC_NUMS
from chemutils.CombiCDB.Const import DISFAVORED, RETRO_ONLY
from chemutils.CombiCDB.Const import ENCRYPTION_KEY

"""Indexes for identifying components of a reaction step"""
REACTANTS = 0
REAGENT = 1
PRODUCTS = 2
SUPPLEMENTAL = 3


class BaseReagent(RowItemModel):
    """Generic object class to model reaction reagents.
    Should be able to provide a description of what the reagent is and,
    given some reactant molecules, produce legitimate product molecules
    that such a reagent would result in.
    """

    # Reactant molecules inherent to using this reagent
    inherentReactants = None

    # Product molecules inherently expected when using this reagent (e.g., H2O in condensation reactions)
    inherentProducts = None
    inherentProductSmiSet = None

    # When a reaction can be recursively applied multiple times, decide
    #   whether or not to count the intermediate products along the way
    #   as final products to be returned in the end.
    countIntermediates = False

    # When reporting the final product list from the reagent, indicate whether or not
    #   to filter out expected inherent products (e.g., H2O in condensation reactions)
    removeInherentProducts = True

    # Pass option to reaction processor whether to ignore self-reaction products
    ignoreSelfReactions = True

    # Link back to original reagent if this is a retro-reagent
    originalReagent = None

    # If a retroReagent, when generating "retro-products" (i.e., precursors)
    #   optionally reapply these through the forward forward reagent
    #   to verify that the original target product ("retro-reactant") is reproduced
    # If "enforcement" is on, do not accept retro-products that cannot reproduce the retro-reactants
    enforceRetroReactions = False

    def __init__(self, initData=None, dataKeys=None):
        """Constructor sets up expected attributes.
        Python dictionary is a superclass, so mostly expects
        initialization data in the form of key-value string pairs
        that will then be parsed out by the loadTextProperties method.
        """
        self["inherent_reactant_smiles"] = " "
        self["inherent_product_smiles"] = " "

        # Number of times the reaction should be applied to the reactants
        #   recursively, assuming the reactants can use this same
        #   reaction multiple times.  Default to just 1 step
        self["recursions"] = 1

        # Template string for the reagent description HTML.  Will use
        #   self as a template dictionary to fill in any parameters
        self["html_template"] = " "

        # Reagent depictions
        self["depict_smiles"] = " "
        self["depict"] = " "

        # Javascript function / action to take if image is clicked.  Default to no action
        self["action"] = " "

        self.inherentReactants = []

        # Load in any initialization data
        self.update(initData, dataKeys)

        self.loadTextProperties()

    def loadTextProperties(self):
        """Look for any text-valued properties that should be
        parsed into some object or list for persistence in the model
        """
        if self["inherent_reactant_smiles"] is not None:
            self.inherentReactants = []
            smilesList = splitCompositeSmilesToList(
                self["inherent_reactant_smiles"], retainCounterIons=False
            )
            for smiles in smilesList:
                if len(smiles.strip()) > 0:
                    # log.critical('Loading these SMILES : "%s"' % smiles)
                    # log.critical('SMILES Type : %s' % str(type(smiles)))
                    self.inherentReactants.append(molBySmiles(smiles.strip()))

        if self["inherent_product_smiles"] is not None:
            self.inherentProducts = []
            self.inherentProductSmiSet = set()
            smilesList = splitCompositeSmilesToList(
                self["inherent_product_smiles"], retainCounterIons=False
            )
            for smiles in smilesList:
                if len(smiles.strip()) > 0:
                    mol = molBySmiles(smiles)
                    self.inherentProducts.append(mol)
                    self.inherentProductSmiSet.add(createStandardSmiString(mol))

        if self["recursions"] is not None:
            self["recursions"] = int(self["recursions"])

    def html(self):
        """Return an HTML representation of the reagent suitable for
        presentation on the Web.  Should be based on the reagent's
        html_template property.  If '%(depict)s' is present in the
        template, the reagent's "depict" property should be set
        to a suitable IMG tag first.
        """
        return self["html_template"] % self

    def __str__(self):
        """Return a simple string representation of the reagent,
        suitable for list and debugging displays.
        Usually just the reagent's description, unless that's
        unavailable, then fall back on reagent_id.
        """
        if "description" in self:
            return self["description"]
        else:
            return str(self["reagent_id"])

    def isRetro(self):
        """Determine whether or not this represents a retroReagent"""
        return self.originalReagent is not None

    def __call__(self, reactants, supplementalData=None, enumerateRacemicMixtures=True):
        """Function call (parantheses) operator overload for standard
        call to generate products based off this reagent and any
        starting reactant materials.

        supplementalData: If provided, should be an instance of SupplementalDataModel.
            Besides the standard product return values, additional data will be stored
            in this model for the caller to inspect.
        enumerateRacemicMixtures:  Will check all of the provided reactants
            and ensure all of their stereocenters are specified one way or another.
            If they are not specified, all possible isomer configurations will be tried.
        """
        if supplementalData is None:
            # Even if user is not interested, just prepare a container to catch and discard
            #   the data rather than having to do a bunch of null / None pointer checks later
            supplementalData = SupplementalDataModel()

        # Pre-build info on reactants so can check against them for retro-reagents and post-processing
        reactantSmiList = []
        for reactant in reactants:
            reactantSmiList.extend(
                splitCompositeMolToSmilesList(reactant, retainCounterIons=False)
            )
        reactantSmiSet = set(reactantSmiList)

        inherentReactantSmiSet = set()
        for inherentReactant in self.inherentReactants:
            inherentReactantSmiSet.add(createStandardSmiString(inherentReactant))

        products = []
        if enumerateRacemicMixtures:
            # First check if need to enumerate any reactant mixtures
            for stereoSpecifiedReactants in enumerateRacemicMixtureList(reactants):
                # Recursive call afer mixtures enumerated
                products.extend(
                    self(
                        stereoSpecifiedReactants,
                        supplementalData,
                        enumerateRacemicMixtures=False,
                    )
                )
        else:
            # Standard generation run
            products = self.generateProducts(reactants, supplementalData)

            if self.isRetro() and (
                self.enforceRetroReactions or supplementalData is not None
            ):
                # This must be a retro reagent.  Make sure reapplying the forward reactions
                #   can actually regenerate the (retro) reactant
                if supplementalData is None:
                    # Don't actually want to return supplemental data, but do want to enforce retro-reactions
                    # Produce a placeholder supplemental data container in the meantime
                    supplementalData = SupplementalDataModel()
                if supplementalData.retroValidationProductsByPrecursorSmi is None:
                    # Initialize, but only if currently blank.  Maybe reused over multiple function calls (as in enumerating racemic mixtures)
                    supplementalData.retroReactantsReproduced = False
                    supplementalData.retroValidationProductsByPrecursorSmi = dict()
                for iProduct in xrange(len(products) - 1, -1, -1):
                    forwardReactantComposite = products[iProduct]
                    forwardReactantCompositeSmi = createStandardSmiString(
                        forwardReactantComposite
                    )

                    # Reapply these after dumping inherent reactants, and separating out composites
                    forwardReactantSmiList = splitCompositeMolToSmilesList(
                        forwardReactantComposite, retainCounterIons=False
                    )
                    forwardReactantList = []
                    for forwardReactantSmi in forwardReactantSmiList:
                        if forwardReactantSmi not in inherentReactantSmiSet:
                            forwardReactant = molBySmiles(forwardReactantSmi)
                            forwardReactantList.append(forwardReactant)
                            log.debug(
                                "Forward Reactant: "
                                + createStandardSmiString(forwardReactant)
                            )

                    # Generate the products.
                    forwardProducts = self.originalReagent(forwardReactantList)
                    supplementalData.retroValidationProductsByPrecursorSmi[
                        forwardReactantCompositeSmi
                    ] = forwardProducts
                    # Link a copy accessible to the caller as a member variable

                    # Our original (retro) reactants had better be amongst these
                    retroReactantsFound = False
                    for forwardProduct in forwardProducts:
                        forwardProductSmiList = splitCompositeMolToSmilesList(
                            forwardProduct, retainCounterIons=False
                        )

                        log.debug("Forward Product: " + str(forwardProductSmiList))

                        for forwardProductSmi in forwardProductSmiList:
                            if forwardProductSmi in reactantSmiSet:
                                reactantSmiSet.remove(forwardProductSmi)

                        if len(reactantSmiSet) < 1:
                            # All retro reactants accounted for, this was a valid retro reaction
                            retroReactantsFound = True
                            break

                    # Make accessible as a member attribute for caller to discover
                    supplementalData.retroReactantsReproduced = (
                        supplementalData.retroReactantsReproduced or retroReactantsFound
                    )

                    log.debug("RetroReactantsFound: " + str(retroReactantsFound))
                    if not retroReactantsFound and self.enforceRetroReactions:
                        # This does not appear to be a valid retro-reaction then,
                        #   since reapplying the forward reaction doesn't reproduce the same results
                        del products[iProduct]

        # Do post-processing run, even if did "enumerateRacemicMixtures,
        #   as this step will check for redundancies across all recursions
        self.postProcessing(products, reactantSmiSet, supplementalData)

        return products

    def postProcessing(self, products, reactantSmiSet, supplementalData):
        """Post-processing to clean up product information.
        The product generation / prediction may have left behind
        sentinel value formal charge codes on some atoms and
        other similar issues that the end user would be uninterested in.
        """

        productSmiSet = set()
        # Check for redundancies
        extraProducts = []
        # Any extra products generated, perhaps by chirality enumeration

        for iProduct in xrange(
            len(products) - 1, -1, -1
        ):  # Reverse iteration in case need to delete any items
            product = products[iProduct]

            product = self.postProcess(product, extraProducts, reactantSmiSet)
            if product is None:
                # Must have been completely rejected
                products.pop(iProduct)
                continue
            else:
                products[iProduct] = product

            productSmi = createStandardSmiString(product)

            # Check to make sure this wasn't identified as a "disfavored" product
            if productSmi in supplementalData.disfavoredProductDict:
                products.pop(iProduct)
                continue

            # No other modifications done.  Add to unique set, or remove if it is redundant
            if productSmi in productSmiSet:
                products.pop(iProduct)
                # Already in product set, remove redundant copy
                continue
            else:
                productSmiSet.add(productSmi)

        # In case any extra were produced
        for product in extraProducts:
            productSmi = createStandardSmiString(product)
            if productSmi not in productSmiSet:
                productSmiSet.add(productSmi)
                products.append(product)

    def postProcess(
        self, product, extraProducts=None, reactantSmiSet=None, allowIntermediate=None
    ):
        """Do post-processing cleanup for an individual product molecule.
        extraProducts: Optional list to store any additional products generated during the post-processing
        reactantSmiSet: Set of reactant SMILES to avoid counting these as products
        """
        if extraProducts is None:
            extraProducts = list()
        if reactantSmiSet is None:
            reactantSmiSet = set()
        if allowIntermediate is None:
            allowIntermediate = self.countIntermediates

        # Clear out dev information labels
        changesMade1 = clearSentinelCharges(product, SENTINEL_CHARGE)
        changesMade2 = clearSentinelCharges(product, SENTINEL_OPTIONAL_STEP_CHARGE)
        changesMade = changesMade1 or changesMade2
        # Must separate OR logic, otherwise may short-circuit logic and skip the method call

        if changesMade:
            # See if any apparent carbon radicals immediately adjacent to each other
            # Doesn't make sense, probably a label corrupted bond (especially aromatics)
            #   combine the radicals into an extra bond order
            for atom in product.GetAtoms():
                if isCarbonRadical(atom):
                    # Look for carbon radical neighbors
                    for neighbor in atom.GetAtoms():
                        if isCarbonRadical(neighbor):
                            # Found 2 carbon radicals adjacent to each other
                            # Highly unlikely, should just be a bond between them
                            bond = atom.GetBond(neighbor)
                            bond.SetOrder(bond.GetOrder() + 1)
                            # This may have reestablished aromaticity, allow perception
                            OEAssignAromaticFlags(product)

        # Disallow intermediate / transition state products,
        #   crudely identified as those with net formal charges
        # Later should add checks to eliminate free radical products as well
        # Counter ions should be retained to allow for charged components
        netCharge = OENetCharge(product)
        componentProductSmiList = splitCompositeMolToSmilesList(
            product, retainCounterIons=True
        )
        componentsRemoved = False

        mol = OEGraphMol()
        for iComponent in xrange(len(componentProductSmiList) - 1, -1, -1):
            componentSmi = componentProductSmiList[iComponent]
            OEParseSmiles(mol, componentSmi)
            rejectComponent = False

            componentCharge = OENetCharge(mol)
            if not rejectComponent and componentSmi in reactantSmiSet:
                # Don't list any original reactants in the products since it's kind of pointless
                rejectComponent = True

            if not rejectComponent:
                # If have a special "reject charge" value explicitly set, do so
                for atom in mol.GetAtoms():
                    if atom.GetFormalCharge() == SENTINEL_REJECT_CHARGE:
                        rejectComponent = True
                        break

            if rejectComponent:
                netCharge -= componentCharge
                componentProductSmiList.pop(iComponent)
                componentsRemoved = True

            mol.Clear()

        # Do a separate round for checking unbalanced charges, otherwise sentinel charges
        #   like the SENTINEL_REJECT_CHARGE will mess up the net charge calculations
        for iComponent in xrange(len(componentProductSmiList) - 1, -1, -1):
            componentSmi = componentProductSmiList[iComponent]
            OEParseSmiles(mol, componentSmi)
            rejectComponent = False

            componentCharge = OENetCharge(mol)
            if (
                not rejectComponent
                and not allowIntermediate
                and componentCharge * netCharge > 0
            ):
                # Net charge is non-zero and this component is contributing to it.  Dump it
                rejectComponent = True

            if rejectComponent:
                netCharge -= componentCharge
                componentProductSmiList.pop(iComponent)
                componentsRemoved = True

            mol.Clear()

        if componentsRemoved:
            if len(componentProductSmiList) < 1:
                # Nothing left, dump this product
                product = None
            else:
                newProduct = molBySmiles(
                    str.join(SMILES_MOL_DELIM, componentProductSmiList)
                )
                orig_warning_level_id = product.GetIntData("warning_level_id")
                orig_warning_message = product.GetStringData("warning_message")

                # Likewise do final chirality check after sentinel charges were cleared
                for iCopy, copyProduct in enumerate(
                    enumerateRacemicMixture(newProduct)
                ):
                    copyProduct.SetIntData("warning_level_id", orig_warning_level_id)
                    copyProduct.SetStringData("warning_message", orig_warning_message)
                    if iCopy == 0:
                        product = copyProduct
                    else:
                        extraProducts.append(copyProduct)

        return product

    def generateProducts(self, reactants, supplementalData=None):
        """Main method to be overriden by subclasses"""
        raise NotImplementedError("generateProducts should be overriden in subclass")

    def retro(self):
        """Return a retro version of the current reagent.  Not quite a perfect inverse,
        but roughly speaking, applying the retro reagent to products of the normal
        reagent, should yield the normal reactants as retro products.
        """
        raise NotImplementedError("retro should be overriden in subclass")


class SMIRKSReagent(BaseReagent):
    """Specific implementation of BaseReagent that uses
    SMIRKS string processing to execute the transformations.
    """

    # Store reaction profiles as dictionary objects (RowIteModels) with entries for
    #   smirks, priority, warning_level, warning_message,
    #   reaction_category_ids (+ reactionCategoryIDSet) and an actual OELibraryGen object
    reactionProfileList = None
    reactionProcessor = None

    # Extra check, maybe unnecessary.  Accepted product molecules must contain
    #   at at least one atom from this list (usually carbon).
    requiredAtomicNums = REQUIRED_ATOMIC_NUMS

    # Whether the reagent has specific intended uses specified
    intendedUsageSpecified = None

    """Default encryption / decryption for SMIRKS strings in the database"""
    smirksCodec = Codec(ENCRYPTION_KEY)

    def reactionModelCompositeKey(reactionModel):
        """Return a tuple representing the composite key to identify
        the core features of a reaction profile model.
        """
        return (
            reactionModel["reagent_id"],
            reactionModel["reaction_profile_id"],
            reactionModel["priority"],
        )

    reactionModelCompositeKey = staticmethod(reactionModelCompositeKey)

    def __init__(self, initData=None, dataKeys=None):
        """Constructor, based on same dictionary arrangement
        of text-based key-value properties which are parsed
        by the loadTextProperties method.
        """

        self["reaction_profile_id_list"] = None
        self["reagent_reaction_profile_id_list"] = None
        self["smirks_list"] = " "
        self["description_list"] = " "
        self["mechanism_label_smirks_list"] = " "
        self["mechanism_arrow_codes_list"] = " "
        self["priority_list"] = None
        self["warning_level_id_list"] = None
        self["warning_message_list"] = None
        self["reaction_category_ids_list"] = None
        self["pre_status_threshold_list"] = None
        self["pre_status_minimum_list"] = None
        self["post_status_list"] = None
        self["retro_ready_list"] = None
        self["intended_usage_list"] = None
        self["expected_reactants"] = None
        BaseReagent.__init__(self, initData, dataKeys)

    def loadTextProperties(self):
        """Parse text-based key-value properties
        into actual objects, numbers, lists, etc.
        """

        BaseReagent.loadTextProperties(self)

        # Record reaction_profile_id for different SMIRKS
        reaction_profile_idList = None
        if self["reaction_profile_id_list"] is not None:
            reaction_profile_idList = []
            for reaction_profile_idStr in self["reaction_profile_id_list"].split(
                REACTION_PROFILE_DELIM
            ):
                reaction_profile_id = None
                try:
                    reaction_profile_id = int(reaction_profile_idStr)
                except:
                    reaction_profile_id = 0
                    # Default value if none specified
                reaction_profile_idList.append(reaction_profile_id)

        # Record reagent_reaction_profile_id for different SMIRKS
        reagent_reaction_profile_idList = None
        if self["reagent_reaction_profile_id_list"] is not None:
            reagent_reaction_profile_idList = []
            for reagent_reaction_profile_idStr in self[
                "reagent_reaction_profile_id_list"
            ].split(REACTION_PROFILE_DELIM):
                reagent_reaction_profile_id = None
                try:
                    reagent_reaction_profile_id = int(reagent_reaction_profile_idStr)
                except:
                    reagent_reaction_profile_id = 0
                    # Default value if none specified
                reagent_reaction_profile_idList.append(reagent_reaction_profile_id)

        # Record priority ratings for different SMIRKS
        priorityList = None
        if self["priority_list"] is not None:
            priorityList = []
            for priorityStr in self["priority_list"].split(REACTION_PROFILE_DELIM):
                priority = None
                try:
                    priority = int(priorityStr)
                except:
                    priority = 0
                    # Default value if none specified
                priorityList.append(priority)

        # Record warning levels and messages
        warningLevelList = None
        if self["warning_level_id_list"] is not None:
            warningLevelList = []
            for warningLevelStr in self["warning_level_id_list"].split(
                REACTION_PROFILE_DELIM
            ):
                warning_level_id = None
                try:
                    warning_level_id = int(warningLevelStr)
                except:
                    warning_level_id = 0
                    # Default value if none specified
                warningLevelList.append(warning_level_id)

        warningMessageList = None
        if self["warning_message_list"] is not None:
            warningMessageList = []
            for warning_message in self["warning_message_list"].split(
                REACTION_PROFILE_DELIM
            ):
                warningMessageList.append(str(warning_message))

        # Record reaction category criteria
        reactionCategoryIDsStrList = None
        reactionCategoryIDSetList = None
        if self["reaction_category_ids_list"] is not None:
            reactionCategoryIDsStrList = []
            reactionCategoryIDSetList = []
            for reactionCategoryIDsStr in self["reaction_category_ids_list"].split(
                REACTION_PROFILE_DELIM
            ):
                reactionCategoryIDSet = None
                try:
                    idList = reactionCategoryIDsStr.split(VALUE_DELIM)
                    for i, idStr in enumerate(idList):
                        idList[i] = int(idStr)
                    reactionCategoryIDSet = set(idList)
                except:
                    # Default value if none specified
                    reactionCategoryIDsStr = ""
                    reactionCategoryIDSet = EMPTY_SET
                reactionCategoryIDsStrList.append(reactionCategoryIDsStr)
                reactionCategoryIDSetList.append(reactionCategoryIDSet)

        # Record post status modifications and pre status checks
        preStatusThresholdList = None
        if self["pre_status_threshold_list"] is not None:
            preStatusThresholdList = []
            for preStatusThresholdStr in self["pre_status_threshold_list"].split(
                REACTION_PROFILE_DELIM
            ):
                pre_status_threshold = None
                try:
                    pre_status_threshold = int(preStatusThresholdStr)
                except:
                    pre_status_threshold = 0
                    # Default value if none specified
                preStatusThresholdList.append(pre_status_threshold)

        preStatusMinimumList = None
        if self["pre_status_minimum_list"] is not None:
            preStatusMinimumList = []
            for preStatusMinimumStr in self["pre_status_minimum_list"].split(
                REACTION_PROFILE_DELIM
            ):
                pre_status_minimum = None
                try:
                    pre_status_minimum = int(preStatusMinimumStr)
                except:
                    pre_status_minimum = 0
                    # Default value if none specified
                preStatusMinimumList.append(pre_status_minimum)

        postStatusList = None
        if self["post_status_list"] is not None:
            postStatusList = []
            for postStatusStr in self["post_status_list"].split(REACTION_PROFILE_DELIM):
                post_status = None
                try:
                    post_status = int(postStatusStr)
                except:
                    post_status = 0
                    # Default value if none specified
                postStatusList.append(post_status)

        retroReadyList = None
        if self["retro_ready_list"] is not None:
            retroReadyList = []
            for retroReadyStr in self["retro_ready_list"].split(REACTION_PROFILE_DELIM):
                retro_ready = retroReadyStr in ("True", "t", "1")
                retroReadyList.append(retro_ready)

        intendedUsageList = None
        if self["intended_usage_list"] is not None:
            intendedUsageList = []
            for intendedUsageStr in self["intended_usage_list"].split(
                REACTION_PROFILE_DELIM
            ):
                intended_usage = intendedUsageStr in ("True", "t", "1")
                intendedUsageList.append(intended_usage)

        descriptionList = None
        if self["description_list"] is not None:
            descriptionList = self["description_list"].split(REACTION_PROFILE_DELIM)

        mechanism_label_smirksList = None
        if (
            self["mechanism_label_smirks_list"] is not None
            and len(self["mechanism_label_smirks_list"].strip()) > 0
        ):
            mechanism_label_smirksList = self["mechanism_label_smirks_list"].split(
                REACTION_PROFILE_DELIM
            )

        mechanism_arrow_codesList = None
        if (
            self["mechanism_arrow_codes_list"] is not None
            and len(self["mechanism_arrow_codes_list"].strip()) > 0
        ):
            mechanism_arrow_codesList = self["mechanism_arrow_codes_list"].split(
                REACTION_PROFILE_DELIM
            )

        #####################################
        #   Mkayala
        #####################################

        # Now parse the actual SMIRKS reaction profiles
        if self["smirks_list"] is not None:
            smirksList = self["smirks_list"].split(REACTION_PROFILE_DELIM)
            sortedList = []
            self.intendedUsageSpecified = False
            for i, smirks in enumerate(smirksList):
                reactionProfile = RowItemModel()
                reactionProfile["smirks"] = smirks

                description = ""
                if descriptionList is not None:
                    description = descriptionList[i]
                reactionProfile["description"] = description

                mechanism_label_smirks = ""
                if mechanism_label_smirksList is not None:
                    mechanism_label_smirks = mechanism_label_smirksList[i]
                reactionProfile["mechanism_label_smirks"] = mechanism_label_smirks

                mechanism_arrow_codes = ""
                if mechanism_arrow_codesList is not None:
                    mechanism_arrow_codes = mechanism_arrow_codesList[i]
                reactionProfile["mechanism_arrow_codes"] = mechanism_arrow_codes

                reaction_profile_id = 0
                if reaction_profile_idList is not None:
                    reaction_profile_id = reaction_profile_idList[i]
                reactionProfile["reaction_profile_id"] = reaction_profile_id

                reagent_reaction_profile_id = 0
                if reagent_reaction_profile_idList is not None:
                    reagent_reaction_profile_id = reagent_reaction_profile_idList[i]
                reactionProfile[
                    "reagent_reaction_profile_id"
                ] = reagent_reaction_profile_id

                priority = 0
                if priorityList is not None:
                    priority = priorityList[i]
                reactionProfile["priority"] = priority

                warning_level_id = 0
                if warningLevelList is not None:
                    warning_level_id = warningLevelList[i]
                reactionProfile["warning_level_id"] = warning_level_id

                warning_message = ""
                if warningMessageList is not None:
                    warning_message = warningMessageList[i]
                reactionProfile["warning_message"] = warning_message

                reactionCategoryIDsStr = ""
                if reactionCategoryIDsStrList is not None:
                    reactionCategoryIDsStr = reactionCategoryIDsStrList[i]
                reactionProfile["reaction_category_ids"] = reactionCategoryIDsStr

                reactionCategoryIDSet = None
                if reactionCategoryIDSetList is not None:
                    reactionCategoryIDSet = reactionCategoryIDSetList[i]
                else:
                    reactionCategoryIDSet = EMPTY_SET
                reactionProfile["reactionCategoryIDSet"] = reactionCategoryIDSet

                pre_status_threshold = 0
                if preStatusThresholdList is not None:
                    pre_status_threshold = preStatusThresholdList[i]
                reactionProfile["pre_status_threshold"] = pre_status_threshold

                pre_status_minimum = 0
                if preStatusMinimumList is not None:
                    pre_status_minimum = preStatusMinimumList[i]
                reactionProfile["pre_status_minimum"] = pre_status_minimum

                post_status = 0
                if postStatusList is not None:
                    post_status = postStatusList[i]
                reactionProfile["post_status"] = post_status

                retro_ready = False
                if retroReadyList is not None:
                    retro_ready = retroReadyList[i]
                reactionProfile["retro_ready"] = retro_ready

                intended_usage = False
                if intendedUsageList is not None:
                    intended_usage = intendedUsageList[i]
                reactionProfile["intended_usage"] = intended_usage

                self.intendedUsageSpecified = (
                    self.intendedUsageSpecified or intended_usage
                )

                # The actual library generator object.
                libgen = OELibraryGen(smirks)
                reactionProfile["libgen"] = libgen

                # Only want to instantiate libgens once (initialization is relatively expensive), but
                #   these are NOT thread-safe.  They communicate through the setting and retrieval
                #   of member variables.  Have a threading.Lock object for each to track each ones isolated usage
                reactionProfile["threadLock"] = Lock()

                sortedList.append((priority, reactionProfile))

            # Sort the reaction profiles by priority score (descending)
            sortedList.sort()
            sortedList.reverse()
            # Desc order

            # Now extract back out only the reaction profile objects
            self.reactionProfileList = []
            for priority, reactionProfile in sortedList:
                self.reactionProfileList.append(reactionProfile)

        self.reactionProcessor = ReactionProcessor()
        self.reactionProcessor.includeUnusedReactants = True

    def generateProducts(self, reactants, supplementalData=None):
        """Main method for generating / predicting products from some input reactants.

        Although this can be called directly, it is advised to use the __call__ construct
        instead.  That is a wrapper method with additional functions such as
        post-processing of generated products and automatic enumeration of reactant
        stereoisomers.
        """

        # Fresh model for tracking data per reaction
        trackingData = TrackingDataModel()

        # Keep record of reactants as one composite mol, in case need to recurse
        reactantSmiList = []
        for reactant in reactants:
            reactantSmiList.append(createStandardSmiString(reactant))
        trackingData.originalReactantMolList.extend(reactants)

        self.reactionProcessor.ignoreSelfReactions = self.ignoreSelfReactions

        finalProductList = self.__recursiveGeneration(
            reactants, trackingData=trackingData, supplementalData=supplementalData
        )

        return finalProductList

    def retro(self):
        """Return a new Reagent model which is the retro / inverse
        version of the current one.  Will only produce meaningful results
        if the starting reagent is properly setup with "retro_ready"
        reaction profiles.
        """
        retroProps = dict(self)
        # Copy over initialization properties

        # Reverse all of the SMIRKS reaction profiles that are "retro_ready," ignore the rest
        retroSmirksList = []
        retroDescriptionList = []
        retroPriorityList = []
        retroPreStatusThresholdList = []
        retroPreStatusMinimumList = []
        retroPostStatusList = []
        for reactionProfile in self.reactionProfileList:
            if reactionProfile["retro_ready"]:
                # Flip arrangement of SMIRKS string
                smirks = reactionProfile["smirks"]
                smirksComponents = smirks.split(REACTION_DELIM)
                smirksComponents.reverse()
                retroSmirks = str.join(REACTION_DELIM, smirksComponents)
                retroSmirksList.append(retroSmirks)

                description = reactionProfile["description"]
                retroDescriptionList.append("(RETRO) " + description)

                # Actually retain priority list, cause most important reactions should still be most important?
                priority = reactionProfile["priority"]
                retroPriorityList.append(str(priority))

                # Invert pre / post-status conditions to work in reverse
                # Awkard when have a range for pre-status minimum and maximum
                #   but a single post-status
                preStatusThreshold = reactionProfile["pre_status_threshold"]
                preStatusMinimum = reactionProfile["pre_status_minimum"]
                postStatus = reactionProfile["post_status"]

                retroPreStatusThresholdList.append(str(postStatus))
                retroPreStatusMinimumList.append(str(0))
                # Necessary to allow reaction chain to initiate, when assume start with status = 0
                retroPostStatusList.append(
                    str((preStatusThreshold + preStatusMinimum) / 2)
                )

        retroProps["smirks_list"] = str.join(REACTION_PROFILE_DELIM, retroSmirksList)
        retroProps["description_list"] = str.join(
            REACTION_PROFILE_DELIM, retroDescriptionList
        )
        retroProps["priority_list"] = str.join(
            REACTION_PROFILE_DELIM, retroPriorityList
        )
        retroProps["pre_status_threshold_list"] = str.join(
            REACTION_PROFILE_DELIM, retroPreStatusThresholdList
        )
        retroProps["pre_status_minimum_list"] = str.join(
            REACTION_PROFILE_DELIM, retroPreStatusMinimumList
        )
        retroProps["post_status_list"] = str.join(
            REACTION_PROFILE_DELIM, retroPostStatusList
        )

        retroProps["inherent_reactant_smiles"] = self["inherent_product_smiles"]
        retroProps["inherent_product_smiles"] = self["inherent_reactant_smiles"]

        retroProps["description"] = "(RETRO) " + self["description"]

        retroReagent = SMIRKSReagent(retroProps)

        # Count intermediates since forward reaction may have already been acting on a halfway product
        retroReagent.countIntermediates = True

        # Link back to original
        retroReagent.originalReagent = self

        return retroReagent

    def __recursiveGeneration(
        self,
        reactants,
        depth=0,
        requiredReactant=None,
        trackingData=None,
        supplementalData=None,
    ):
        """Returns any products generated in this step AND all subsequent recursions

        requiredReactant:
            If specified, should match one of the provided reactants.  Only accept
            reaction products where this reactant was used (i.e., where it does
            NOT appear UNmodified in the product).
        """
        # Recursion base cases
        if depth >= self["recursions"]:
            return []

        # Check if we've already seen these reactants
        reactantSmiList = molToSmiList(reactants)
        smiListStr = str.join(MOL_LIST_DELIM, reactantSmiList)
        if smiListStr not in trackingData.recursiveReactantSmiDict:
            log.debug("Recursion depth %d, Begin" % depth)
            for iReactant, reactant in enumerate(reactants):
                log.debug(
                    "Reactant %d: %s" % (iReactant, createStandardSmiString(reactant))
                )

            # Main action
            (alreadyGenerated, productList) = self.__singleStepGeneration(
                reactants,
                requiredReactant,
                trackingData=trackingData,
                supplementalData=supplementalData,
            )

            # See if we should do any more recursions, and filter product list to acceptable ones
            acceptedProductList = []
            # Weird, should be able to use original productList and pop unwanted items, but bug seems to result in empty lists unexpectedly
            for product in productList:
                productSmiles = createStandardSmiString(product)

                # Resolve into component products if mixture
                productSmilesList = splitCompositeSmilesToList(
                    productSmiles, retainCounterIons=False
                )

                # Recursive generation,
                #   Only count recursive products that actually use the new components
                #   and only continue if no error reaction messages
                allSubProducts = []

                reactants = []
                # Next recursion reactant container
                componentSubProductsByComponentSmi = dict()
                for componentProductSmiles in productSmilesList:
                    if componentProductSmiles in reactantSmiList:
                        # Skip this component, looks just like an unused reactant carried over
                        continue

                    componentProduct = molBySmiles(componentProductSmiles)
                    componentProduct.SetIntData(
                        "warning_level_id", product.GetIntData("warning_level_id")
                    )
                    # Propagate warning data
                    componentProduct.SetStringData(
                        "warning_message", product.GetStringData("warning_message")
                    )
                    # Propagate warning data
                    componentProduct.SetIntData(
                        "post_status", product.GetIntData("post_status")
                    )
                    componentProduct.SetIntData(
                        "priority", product.GetIntData("priority")
                    )
                    reactants.append(componentProduct)
                    # Use each of last recursion's products as next reactant

                    componentSubProducts = []

                    if self["possible_reactants"] < 2:
                        componentSubProducts.extend(
                            self.__recursiveGeneration(
                                reactants,
                                depth + 1,
                                trackingData=trackingData,
                                supplementalData=supplementalData,
                            )
                        )
                    elif self["possible_reactants"] >= 2:
                        # Additional spots expected
                        # First check if any unimolecular (just prior product) reactions can occur
                        preReactionOptionalCharges = countAtomCharges(
                            reactants, SENTINEL_OPTIONAL_STEP_CHARGE
                        )
                        unimolecularProducts = self.__recursiveGeneration(
                            reactants,
                            depth + 1,
                            trackingData=trackingData,
                            supplementalData=supplementalData,
                        )
                        postReactionOptionalCharges = countAtomCharges(
                            unimolecularProducts, SENTINEL_OPTIONAL_STEP_CHARGE
                        )
                        componentSubProducts.extend(unimolecularProducts)

                        if (
                            len(unimolecularProducts) < 1
                            or postReactionOptionalCharges - preReactionOptionalCharges
                            > 0
                        ):
                            # No unimolecular products or unimolecular reaction was just some kind of labeling or
                            #   other optional step that left sentinel charge labels (i.e., allylic resonance)
                            # Try to fill in additional positions with original reactants (one at a time)
                            # At this point, only consider up to binary combinations
                            # If ternary or higher, can't just add 1 item from the original reactant list,
                            #   would need to add all possible (n choose k) combinations, where n = number of original reactants, k = number of spots to fill
                            #   but that's too complicated to be worthwhile

                            # Don't accept all sub products from different original reactant combinations,
                            #   only those from the highest reaction profile priority
                            topSubPriority = -1
                            originalReactantByPriority = dict()
                            for (
                                originalReactant
                            ) in trackingData.originalReactantMolList:
                                reactants.append(originalReactant)

                                # Look-ahead one recursion to see which will yield the next highest reaction profile priorty.  Only follow those branches
                                # Should not have to worry much about redundant calculation since the results are cached in the trackingData object
                                subPriority = -1
                                (
                                    alreadyGenerated,
                                    subProductList,
                                ) = self.__singleStepGeneration(
                                    reactants,
                                    componentProduct,
                                    trackingData=trackingData,
                                    supplementalData=supplementalData,
                                )
                                for subProduct in subProductList:
                                    subPriority = max(
                                        subPriority, subProduct.GetIntData("priority")
                                    )
                                if subPriority not in originalReactantByPriority:
                                    originalReactantByPriority[subPriority] = list()
                                originalReactantByPriority[subPriority].append(
                                    originalReactant
                                )
                                topSubPriority = max(topSubPriority, subPriority)

                                reactants.pop()
                                # Revert state
                            # Now collect / accept only those original reactant branches based on the top sub priority
                            if topSubPriority in originalReactantByPriority:
                                for originalReactant in originalReactantByPriority[
                                    topSubPriority
                                ]:
                                    reactants.append(originalReactant)
                                    componentSubProducts.extend(
                                        self.__recursiveGeneration(
                                            reactants,
                                            depth + 1,
                                            componentProduct,
                                            trackingData=trackingData,
                                            supplementalData=supplementalData,
                                        )
                                    )
                                    reactants.pop()
                                    # Revert state
                    reactants.pop()
                    # Restore state

                    componentSubProductsByComponentSmi[
                        componentProductSmiles
                    ] = componentSubProducts
                    allSubProducts.extend(componentSubProducts)

                if self.countIntermediates:
                    # Record the mixture product to the final list if we're counting intermediates
                    self.__addToFinalProductList(
                        product,
                        productSmilesList,
                        acceptedProductList,
                        depth,
                        trackingData=trackingData,
                        supplementalData=supplementalData,
                    )

                if len(allSubProducts) > 0:
                    # Record the mixture of all replaced sub-products as well, like modified versions of the products,
                    #   becomes important for reactions that yield more products than starting materials (i.e., retro-combinatorial reactions)
                    componentSubProductsList = []
                    for componentProductSmiles in productSmilesList:
                        componentSubProducts = [molBySmiles(componentProductSmiles)]
                        if componentProductSmiles in componentSubProductsByComponentSmi:
                            lookupProducts = componentSubProductsByComponentSmi[
                                componentProductSmiles
                            ]
                            if len(lookupProducts) > 0:
                                componentSubProducts = lookupProducts
                        componentSubProductsList.append(componentSubProducts)

                    for compositeSubProduct in self.__allComponentCombinations(
                        componentSubProductsList
                    ):
                        subProductSmiList = splitCompositeMolToSmilesList(
                            compositeSubProduct, retainCounterIons=False
                        )
                        self.__addToFinalProductList(
                            compositeSubProduct,
                            subProductSmiList,
                            acceptedProductList,
                            depth,
                            trackingData=trackingData,
                            supplementalData=supplementalData,
                        )

                if len(allSubProducts) < 1:
                    # No subProducts generated, this must be a terminal path.  Add to final product list,
                    if not self.countIntermediates:
                        # If counting intermediates, this product would have already been counted
                        self.__addToFinalProductList(
                            product,
                            productSmilesList,
                            acceptedProductList,
                            depth,
                            trackingData=trackingData,
                            supplementalData=supplementalData,
                        )

            log.debug("Recursion depth %d, End" % depth)
            trackingData.recursiveReactantSmiDict[smiListStr] = acceptedProductList
        return trackingData.recursiveReactantSmiDict[smiListStr]

    def __singleStepGeneration(
        self, reactants, requiredReactant=None, trackingData=None, supplementalData=None
    ):
        """Single reaction step vs. general case which may have multiple recursions

        requiredReactant:
            If specified, should match one of the provided reactants.  Only accept
            reaction products where this reactant was used (i.e., where it does
            NOT appear UNmodified in the product).
        """

        requiredReactantSmi = None
        if requiredReactant is not None:
            requiredReactantSmi = createStandardSmiString(requiredReactant)

        # Check for the worst warning messages to propagate
        # Check post status of reactants before proceeding past pre status threshold of subsequent reactions
        extremeReactantPostStatus = 0
        reactant_warning_level_id = 0
        reactant_warning_message = ""
        for reactant in reactants:
            (
                reactant_warning_level_id,
                reactant_warning_message,
            ) = self.__worstWarningLevel(
                reactant_warning_level_id,
                reactant_warning_message,
                reactant.GetIntData("warning_level_id"),
                reactant.GetStringData("warning_message"),
            )
            # Check for most extreme (farthest from 0) post status
            # Don't just use max because retro-reagents will use negative values
            reactantPostStatus = reactant.GetIntData("post_status")
            if abs(extremeReactantPostStatus) < abs(reactantPostStatus):
                extremeReactantPostStatus = reactantPostStatus

        # Check if we've already cached the result for these reactant parameters
        #   If not, then record it, so we don't repeat upon recursions
        reactantSmiList = molToSmiList(reactants)
        smiListStr = str.join(MOL_LIST_DELIM, reactantSmiList)
        smiListStr += "%s%s%s%s" % (
            MOL_LIST_DELIM,
            extremeReactantPostStatus,
            MOL_LIST_DELIM,
            requiredReactantSmi,
        )

        alreadyGenerated = smiListStr in trackingData.reactantSmilesDict

        if not alreadyGenerated:
            # Create a copy list, and add inherent reactants
            extendedReactants = list(reactants)
            extendedReactants.extend(self.inherentReactants)
            # Unmodified copy
            reactantsCopy = list(reactants)

            # Product SMILES that should be rejected if ever encountered
            rejectProductSmiSet = set()

            # Try using every available libgen reaction profile,
            #   but desc. order of priority.  Only move down in priority if the higher priority ones produce no product
            productList = []
            currentProductSmiSet = set()
            topPriorityWithProducts = None
            disfavoredProductsPriority = -1
            # Priority of the last reaction profile found which yielded "disfavored" products
            for reactionIndex, reactionProfile in enumerate(self.reactionProfileList):
                priority = reactionProfile["priority"]
                libgen = reactionProfile["libgen"]
                threadLock = reactionProfile["threadLock"]

                acceptableReactionProfile = (
                    topPriorityWithProducts is None
                    or priority >= topPriorityWithProducts
                )
                acceptableReactionProfile = (
                    acceptableReactionProfile
                    and extremeReactantPostStatus
                    <= reactionProfile["pre_status_threshold"]
                )
                acceptableReactionProfile = (
                    acceptableReactionProfile
                    and extremeReactantPostStatus
                    >= reactionProfile["pre_status_minimum"]
                )
                acceptableReactionProfile = acceptableReactionProfile and not (
                    reactionProfile["warning_level_id"] == DISFAVORED
                    and reactionProfile["priority"] < disfavoredProductsPriority
                )
                # If already found some disfavored products, don't keep looking for more down the priority list
                acceptableReactionProfile = acceptableReactionProfile and (
                    reactionProfile["warning_level_id"] != RETRO_ONLY or self.isRetro()
                )

                if acceptableReactionProfile:
                    productString = None
                    productCount = None

                    lockAcquired = threadLock.acquire(False)
                    # Verify the libgen is available for use and not locked by another thread
                    try:
                        if not lockAcquired:
                            # Another thread is using this libgen, fall back on new instance (slower, but at least will work)
                            libgen = OELibraryGen(reactionProfile["smirks"])

                        productOEOS = virtual_oemolostream(OEFormat_ISM)
                        productCount = self.reactionProcessor.applyReaction(
                            libgen,
                            extendedReactants,
                            productOEOS,
                            reactionIndex,
                            rejectProductSmiSet=rejectProductSmiSet,
                        )
                        productString = productOEOS.GetString()

                    finally:
                        if lockAcquired:
                            threadLock.release()

                    # Now look for valid products
                    if productCount > 0:
                        (
                            worst_warning_level_id,
                            worst_warning_message,
                        ) = self.__worstWarningLevel(
                            reactant_warning_level_id,
                            reactant_warning_message,
                            reactionProfile["warning_level_id"],
                            reactionProfile["warning_message"],
                        )

                        productOEIS = virtual_oemolistream(productString, OEFormat_ISM)
                        componentMol = OEGraphMol()
                        componentCopy = OEGraphMol()
                        productsAccepted = 0
                        for product in productOEIS.GetOEGraphMols():
                            # Dump any inherent products or otherwise rejected items from the composite mol
                            componentProductSmiList = splitCompositeMolToSmilesList(
                                product, retainCounterIons=False
                            )

                            if requiredReactantSmi in componentProductSmiList:
                                # Required reactant left in product pool
                                #   Don't accept this product, skip past it
                                continue

                            anyComponentsRejected = False
                            rejectAll = False

                            for iComponent in xrange(
                                len(componentProductSmiList) - 1, -1, -1
                            ):  # Reverse iteration for safety since may be deleting items
                                componentProductSmi = componentProductSmiList[
                                    iComponent
                                ]
                                OEParseSmiles(componentMol, componentProductSmi)
                                rejectComponent = False

                                if (
                                    not rejectComponent
                                    and self.removeInherentProducts
                                    and self.inherentProducts is not None
                                ):
                                    # Discard inherent products
                                    if (
                                        componentProductSmi
                                        in self.inherentProductSmiSet
                                    ):
                                        rejectComponent = True
                                    else:
                                        # Check if the post-processed version is there
                                        OEParseSmiles(
                                            componentCopy, componentProductSmi
                                        )
                                        componentProcessed = self.postProcess(
                                            componentCopy, allowIntermediate=True
                                        )
                                        componentCopySmi = createStandardSmiString(
                                            componentProcessed
                                        )
                                        if (
                                            componentCopySmi
                                            in self.inherentProductSmiSet
                                        ):
                                            rejectComponent = True
                                        componentCopy.Clear()

                                if (
                                    not rejectComponent
                                    and componentProductSmi in rejectProductSmiSet
                                ):
                                    # Eliminate components specially recognized for rejection
                                    # In fact, reject entire product since this reaction should never have been able to happen
                                    # Necessary for special handling of non-trans diaxial eliminations in ring systems
                                    rejectComponent = True
                                    rejectAll = True

                                # Do not immediately reject disfavored products.  Keep a record of them so
                                #   subsequent steps can know to avoid them
                                # if not rejectComponent and componentProductSmi in supplementalData.disfavoredProductDict:
                                #    # Likewise for products and intermediates identified as disfavored
                                #    rejectComponent = True;
                                #    rejectAll;

                                if not rejectComponent:
                                    # Eliminate components specially labeled for rejection
                                    for atom in componentMol.GetAtoms():
                                        if (
                                            atom.GetFormalCharge()
                                            == SENTINEL_REJECT_IMMEDIATE_CHARGE
                                        ):
                                            rejectComponent = True

                                if rejectComponent:
                                    componentProductSmiList.pop(iComponent)
                                    anyComponentsRejected = True

                                componentMol.Clear()

                            if rejectAll:
                                product = None
                            elif anyComponentsRejected:
                                # If one was rejected, then modify the final product returned
                                productSmi = str.join(
                                    SMILES_MOL_DELIM, componentProductSmiList
                                )
                                if productSmi != "":
                                    product = molBySmiles(productSmi)
                                    # Check if eliminated everything interesting (organic compounds with carbons),
                                    #   and maybe only left with pointless counter ions ([Na+].[Cl-] or something).
                                    requiredAtomFound = False
                                    for atom in product.GetAtoms():
                                        if (
                                            atom.GetAtomicNum()
                                            in self.requiredAtomicNums
                                        ):
                                            requiredAtomFound = True
                                            break
                                    if not requiredAtomFound:
                                        product = None
                                        # No useful components left, reject product
                                else:
                                    product = None
                                    # All components rejected, no product left

                            if product is not None:
                                # Products just from this one reaction profile.
                                # The primary productList could contain products from multiple reaction profiles
                                #   (though all among them should have the same priority ranking)
                                reactionProfileProducts = []
                                for productCopy in enumerateRacemicMixture(product):
                                    productCopySmi = createStandardSmiString(
                                        productCopy
                                    )
                                    if (
                                        productCopySmi not in currentProductSmiSet
                                    ):  # Don't bother adding redundant products
                                        productCopy.SetIntData(
                                            "warning_level_id", worst_warning_level_id
                                        )
                                        productCopy.SetStringData(
                                            "warning_message", worst_warning_message
                                        )
                                        if reactionProfile["post_status"] >= 0:
                                            productCopy.SetIntData(
                                                "post_status",
                                                reactionProfile["post_status"],
                                            )
                                        else:
                                            # Negative post_status doesn't make sense.  Sentinel value to just propagate last status
                                            productCopy.SetIntData(
                                                "post_status", extremeReactantPostStatus
                                            )
                                        productCopy.SetIntData(
                                            "priority", reactionProfile["priority"]
                                        )

                                        if worst_warning_level_id == DISFAVORED:
                                            # Disfavored product.  Don't accept into the normal product pool,
                                            #   but keep a record of it in case caller wants information on why these were rejected
                                            # Post-process it as well to ensure robust comparison later
                                            productCopy = self.postProcess(
                                                productCopy, allowIntermediate=True
                                            )
                                            productCopySmi = createStandardSmiString(
                                                productCopy
                                            )
                                            if (
                                                productCopy is not None
                                                and productCopySmi
                                                not in supplementalData.disfavoredProductDict
                                            ):
                                                supplementalData.disfavoredProductDict[
                                                    productCopySmi
                                                ] = productCopy
                                                elementaryStep = (
                                                    reactantsCopy,
                                                    reactionProfile,
                                                    [productCopy],
                                                )
                                                supplementalData.elementarySteps.append(
                                                    elementaryStep
                                                )

                                                # Second copy which does not allow intermediate (i.e., net charged or radical structures)
                                                nonIntermedCopy = OEGraphMol(
                                                    productCopy
                                                )
                                                nonIntermedCopy = self.postProcess(
                                                    nonIntermedCopy,
                                                    allowIntermediate=False,
                                                )
                                                nonIntermedCopySmi = (
                                                    createStandardSmiString(
                                                        nonIntermedCopy
                                                    )
                                                )
                                                if (
                                                    nonIntermedCopy is not None
                                                    and nonIntermedCopySmi
                                                    not in supplementalData.disfavoredProductDict
                                                ):
                                                    supplementalData.disfavoredProductDict[
                                                        nonIntermedCopySmi
                                                    ] = nonIntermedCopy
                                                    elementaryStep = (
                                                        reactantsCopy,
                                                        reactionProfile,
                                                        [nonIntermedCopy],
                                                    )
                                                    supplementalData.elementarySteps.append(
                                                        elementaryStep
                                                    )

                                            disfavoredProductsPriority = priority
                                        else:
                                            # Looks like a normal product to report, one way or another
                                            currentProductSmiSet.add(productCopySmi)
                                            reactionProfileProducts.append(productCopy)
                                            productsAccepted += 1
                                        if productCopy is not None:
                                            log.debug(
                                                "Product Generated: %s, by SMIRKS %s %s %s"
                                                % (
                                                    createStandardSmiString(
                                                        productCopy
                                                    ),
                                                    reactionProfile["description"],
                                                    reactionProfile["priority"],
                                                    productCopy.GetIntData(
                                                        "post_status"
                                                    ),
                                                )
                                            )
                                productList.extend(reactionProfileProducts)

                                # Keep a record of the elementary steps used in case caller wants more details
                                if len(reactionProfileProducts) > 0:
                                    elementaryStep = (
                                        reactantsCopy,
                                        reactionProfile,
                                        reactionProfileProducts,
                                    )
                                    supplementalData.elementarySteps.append(
                                        elementaryStep
                                    )

                        if productsAccepted > 0 and topPriorityWithProducts is None:
                            topPriorityWithProducts = priority

            trackingData.reactantSmilesDict[smiListStr] = productList
        else:
            productList = trackingData.reactantSmilesDict[smiListStr]
            productSmiList = molToSmiList(productList)
            productSmiListStr = str.join(MOL_LIST_DELIM, productSmiList)
            log.debug(
                "Products already generated for reactants: %s >> %s"
                % (smiListStr, productSmiListStr)
            )

        return (alreadyGenerated, trackingData.reactantSmilesDict[smiListStr])

    def __addToFinalProductList(
        self,
        product,
        productSmilesList,
        acceptedProductList,
        depth,
        trackingData=None,
        supplementalData=None,
    ):
        productSmi = createStandardSmiString(product)
        if productSmi not in trackingData.compositeProductSmilesSet:
            log.debug(
                "Adding to accepted product list at depth %d: %s %d"
                % (
                    depth,
                    str(productSmilesList),
                    product.GetIntData("warning_level_id"),
                )
            )
            trackingData.compositeProductSmilesSet.add(productSmi)
        trackingData.productSmilesSet.update(productSmilesList)
        acceptedProductList.append(product)

    def __allComponentCombinations(self, componentSubProductsList, currentList=None):
        """componentSubProductsList is a list of sub product lists.
        Pick one from each list, in every combination, to yield a composite product.
        """
        if currentList is None:
            currentList = []
        currentPosition = len(currentList)
        if currentPosition >= len(componentSubProductsList):
            # Base case, no more component positions to fill in
            compositeProduct = OEGraphMol()
            warning_level_id = 0
            warning_message = ""
            for componentProduct in currentList:
                (warning_level_id, warning_message) = self.__worstWarningLevel(
                    warning_level_id,
                    warning_message,
                    componentProduct.GetIntData("warning_level_id"),
                    componentProduct.GetStringData("warning_message"),
                )
                componentSmi = createStandardSmiString(componentProduct)
                OEAddMols(compositeProduct, componentProduct)
            compositeProduct.SetIntData("warning_level_id", warning_level_id)
            compositeProduct.SetStringData("warning_message", warning_message)
            yield compositeProduct
        else:
            # Still more positions to fill, try adding all possible combinations
            nextComponentSubProducts = componentSubProductsList[currentPosition]
            for componentProduct in nextComponentSubProducts:
                currentList.append(componentProduct)
                for compositeProduct in self.__allComponentCombinations(
                    componentSubProductsList, currentList
                ):
                    yield compositeProduct
                currentList.pop()
                # Restore state

    def __worstWarningLevel(self, level1, message1, level2, message2):
        """Given a pair of warning_level_ids with their respective warning_messages,
        return the warning of greater severity, and combine the warning messages.
        """
        if message1 is None or message1 == NULL_STRING:
            message1 = ""
        if message2 is None or message2 == NULL_STRING:
            message2 = ""
        comboMessage = ""
        if level2 > level1:
            comboMessage = message2
            if message1 not in message2:
                # If the next message has not already been observed, then append it
                if len(comboMessage) > 0 and len(message1) > 0:
                    comboMessage += "<br>" + message1
                else:
                    comboMessage += message1
                    # One or the other is blank, no point in adding <br> line break
            return (level2, comboMessage)
        else:
            comboMessage = message1
            if message2 not in message1:
                # If the next message has not already been observed, then append it
                if len(comboMessage) > 0 and len(message2) > 0:
                    comboMessage += "<br>" + message2
                else:
                    comboMessage += message2
                    # One or the other is blank, no point in adding <br> line break
            return (level1, comboMessage)


class TrackingDataModel:
    """Simple struct to store tracking variables while processing reagents.
    Better than using member variables as it keeps the reagent thread-safe this way.
    """

    # Keep track of reactants and products observed to avoid redundancy.
    reactantSmilesDict = None
    # Keyed by reactant list string, value is list of product molecules that came out of it
    recursiveReactantSmiDict = None
    # Same as above but for recursive generation calls, not just single step
    productSmilesSet = None
    compositeProductSmilesSet = None
    originalReactantMolList = None

    def __init__(self):
        """Instantiate member variables in constructor.
        Don't do it at declaration time, otherwise will be treated as a class variable
        shared by all instances of the class.
        """
        self.reactantSmilesDict = dict()
        # Keyed by reactant list string, value is list of product molecules that came out of it
        self.recursiveReactantSmiDict = dict()
        # Same as above but for recursive generation calls, not just single step
        self.productSmilesSet = set()
        self.compositeProductSmilesSet = set()
        self.originalReactantMolList = []


class SupplementalDataModel:
    """Simple struct to store supplemental data that can be returned after
    reagent processing if the caller is interested.
    """

    # If a retro-reagent, report on whether the retro-reactants were reproduced by the retro-products
    retroReactantsReproduced = None
    # Store any that are found in the dict attribute below,
    #   to report on side-products if any
    retroValidationProductsByPrecursorSmi = None

    # After each reaction call, collect all of the sub (elementary) steps that were used
    # Caller may then access this member variable for additional step-wise information,
    #   instead of just the final product result
    # Data should be a tuple including:
    #   - Reactant molecules used in the step as a list of molecule objects
    #   - Reaction Profile object / dictionary used to execute the step (including information like SMIRKS, description, warning levels, etc.)
    #   - Product molecules that were produced in the step as a list of molecule objects
    elementarySteps = None

    # Keep track of products that were rejected as disfavored, but which may still be of
    #   interest to the caller as to why they were rejected.
    # Dictionary keyed by disfavored product SMILES string, value is the molecule object, including warning message
    disfavoredProductDict = None

    def __init__(self):
        """Instantiate member variables in constructor.
        Don't do it at declaration time, otherwise will be treated as a class variable
        shared by all instances of the class.
        """
        self.elementarySteps = []
        self.disfavoredProductDict = dict()

    def filteredProductSmiSet(self, productList):
        """Take the (final) product list for a reaction
        and break it down into component molecules and only counting
        those with relevant / interesting atoms (not boring counter-ions, etc.).
        Collect these as a set of SMILES strings to quickly identify the key products
        of interest in a reaction.
        """
        productSmiSet = set()
        for product in productList:
            componentProducts = splitCompositeMol(product)
            for componentProduct in componentProducts:
                clearSentinelCharges(componentProduct)
                # Only count products which have relevant atoms (and not just boring counter-ions)
                if atomCount(componentProduct, REQUIRED_ATOMIC_NUMS) > 0:
                    componentSmi = createStandardSmiString(componentProduct)
                    productSmiSet.add(componentSmi)
        return productSmiSet

    def filteredElementarySteps(self, reactantList, productList):
        """Given the list of elementarySteps
        (3-ples consisting of (detailReactants, reactionProfile, detailProducts)) for steps
        used to generate the productList from the reactantList,
        return a subset of the list of elementarySteps to only include those that are
        relevant for showing the production of the productList and not any side reactions
        or disfavored steps.
        """
        reactantSmiSet = set()
        productSmiSet = set()

        for reactant in reactantList:
            reactantSmiSet.add(createStandardSmiString(reactant))
        productSmiSet = self.filteredProductSmiSet(productList)

        acceptedStepIndexes = set()

        moreStepsFound = True
        while moreStepsFound:
            moreStepsFound = False

            for iStep, elementaryStep in enumerate(self.elementarySteps):
                (detailReactants, reactionProfile, detailProducts) = elementaryStep
                acceptStep = True

                for detailProduct in detailProducts:
                    detailProductSmi = createStandardSmiString(detailProduct)

                    if detailProductSmi in self.disfavoredProductDict:
                        # Don't include steps that yield disfavored products
                        acceptStep = False
                    else:
                        # Check if any component product is in the desired set
                        detailProductCopy = OEGraphMol(detailProduct)
                        clearSentinelCharges(detailProductCopy)
                        detailProductCopy = createStandardMol(detailProductCopy)
                        componentSmiSet = set(
                            splitCompositeMolToSmilesList(detailProductCopy)
                        )
                        if len(componentSmiSet.intersection(productSmiSet)) < 1:
                            acceptStep = False

                if acceptStep and iStep not in acceptedStepIndexes:
                    # Looks like a valid step to include and have not already counted it
                    moreStepsFound = True
                    acceptedStepIndexes.add(iStep)
                    # Add this steps detail reactants as "products" to traceback further steps,
                    #   but don't traceback to original starting materials, or else may accept some circular sequences
                    for detailReactant in detailReactants:
                        detailReactantCopy = OEGraphMol(detailReactant)
                        clearSentinelCharges(detailReactantCopy)
                        detailReactantCopy = createStandardMol(detailReactantCopy)
                        detailReactantSmi = createStandardSmiString(detailReactantCopy)
                        if detailReactantSmi not in reactantSmiSet:
                            productSmiSet.add(detailReactantSmi)

        """        
        print >> sys.stderr, reactantSmiSet
        print >> sys.stderr, productSmiSet
        for iStep, elementaryStep in enumerate(self.elementarySteps):
            (detailReactants, reactionProfile, detailProducts) = elementaryStep;
            for detailReactant in detailReactants:
                detailReactantSmi = createStandardSmiString(detailReactant);
                print >> sys.stderr, detailReactantSmi,
            print >> sys.stderr, '>>',
            for detailProduct in detailProducts:
                detailProductSmi = createStandardSmiString(detailProduct);
                print >> sys.stderr, detailProductSmi,
            print >> sys.stderr
        """

        # Now copy over the accepted steps
        filteredSteps = []
        for iStep, elementaryStep in enumerate(self.elementarySteps):
            if iStep in acceptedStepIndexes:
                filteredSteps.append(elementaryStep)
        return filteredSteps


def countAtomCharges(molList, chargeValue):
    """Go through every atom of every molecule in the list.
    Count up the number of them whose formal charge is equal
    to the specified chargeValue
    """
    charges = 0
    for mol in molList:
        for atom in mol.GetAtoms():
            correctedCharge = atom.GetFormalCharge() - chargeValue
            if abs(correctedCharge) <= NORMAL_CHARGE_THRESHOLD:
                charges += 1
    return charges


def clearSentinelCharges(mol, chargeValue=None):
    """Look for atoms in the molecule with the specified sentinel charge value
    and clear them out to normal values.  Returns whether any such changes were made.
    If no specific chargeValue is specified, will just set to formal charge to 0
    for any atom whose formal charge is outside of the NORMAL_CHARGE_THRESHOLD.

    This excludes specific sentinel charges for acid-base stand-ins that are meant
        to have specific charge states:
    SENTINEL_CONJUGATE_BASE     #  A-
    SENTINEL_ACID               # HA
    SENTINEL_BASE               #  B
    SENTINEL_CONJUGATE_ACID     # HB+
    """
    changesMade = False
    for atom in mol.GetAtoms():
        if chargeValue is not None:
            correctedCharge = atom.GetFormalCharge() - chargeValue
            if abs(correctedCharge) <= NORMAL_CHARGE_THRESHOLD:
                atom.SetFormalCharge(correctedCharge)
                changesMade = True
        elif atom.GetFormalCharge() == SENTINEL_CONJUGATE_BASE:
            atom.SetFormalCharge(-1)
            changesMade = True
        elif atom.GetFormalCharge() == SENTINEL_CONJUGATE_ACID:
            atom.SetFormalCharge(+1)
            changesMade = True
        elif abs(atom.GetFormalCharge()) > NORMAL_CHARGE_THRESHOLD:
            atom.SetFormalCharge(0)
            changesMade = True
    return changesMade


class ReactionStep(list):
    """Class to represent a single step in a reaction synthesis pathway.
    Actually this is mostly unused since it's just a simple list structure,
    but looking at the implementation here indicates how that list is organized.

    List with 3 cells:
        [0] = List of reactant SMILES strings
        [1] = Reagent object
        [2] = List of product SMILES strings
    Can use the constants REACTANTS, REAGENT and PRODUCTS to reference the cells
    """

    def __init__(self):
        list.__init__(self)
        self.append(["*"])
        # Reactant list with placeholder * molecule
        self.append(None)
        # Null reagent object
        self.append(["*"])
        # Product list with placeholder * molecule

    def __str__(self):
        return reactionStepStr(self)


def reactionStepStr(reactionStep):
    """Convenience function for converting a reaction step into a comprehensible string.
    Reaction step model is currently just a 3-ple consisting of
        - List of reactant SMILES strings (or molecule objects)
        - Reagent object
        - List of product SMILES strings (or molecule objects)

    Optionally may be a 4-ple with the last position holding a SupplementalDataModel
    """
    reactantList = reactionStep[REACTANTS]
    if len(reactantList) > 0 and not isinstance(reactantList[0], str):
        # Assume is a list of molecule objects.  Convert to SMILES strings
        reactantList = []
        for mol in reactionStep[REACTANTS]:
            reactantList.append(createStandardSmiString(mol))
    reactantStr = str.join(MOL_LIST_DELIM, reactantList)

    productList = reactionStep[PRODUCTS]
    if len(productList) > 0 and not isinstance(productList[0], str):
        # Assume is a list of molecule objects.  Convert to SMILES strings
        productList = []
        for mol in reactionStep[PRODUCTS]:
            productList.append(createStandardSmiString(mol))
    productStr = str.join(MOL_LIST_DELIM, productList)

    try:
        if "description" in reactionStep[REAGENT]:
            reagentStr = reactionStep[REAGENT]["description"]
        else:
            reagentStr = reactionStep[REAGENT]["depict_smiles"]
    except:
        # Failover to default string representation
        reagentStr = str(reactionStep[REAGENT])
    return "%s>>%s\t%s" % (reactantStr, productStr, reagentStr)


def reactionStep_GetIntData(reactionStep, dataName):
    """Convenience function to extract a molecule integer data
    associated with a reaction step (assuming it has one),
    that will be stored as a data attribute on the reaction mols.
    """
    for reactantMol in reactionStep[REACTANTS]:
        dataValue = reactantMol.GetIntData(dataName)
        # Just need to find it in the first reactant, the rest should be redundant
        return dataValue
