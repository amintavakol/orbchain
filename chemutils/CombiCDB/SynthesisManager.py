import sys

# from sets import Set;

# from mx.DateTime import now;
from datetime import datetime
from chemutils.Common.Util import molToSmiList
from chemutils.Common.Util import molBySmiles
from chemutils.Common.Util import createStandardSmiString
from chemutils.CombiCDB.ReactionModel import REACTANTS, REAGENT, PRODUCTS
from chemutils.CombiCDB.ReactionModel import (
    SupplementalDataModel,
)
from chemutils.CombiCDB.ReagentManager import ReagentManager, ReagentSearchModel
from chemutils.Common import DBUtil
from chemutils.Common.Model import SQLQuery, modelListFromTable, modelDictFromList

from chemutils.CombiCDB.Const import MOL_LIST_DELIM, REACTION_PROFILE_DELIM
from chemutils.CombiCDB.Const import (
    MAX_TRIES,
    CACHE_CODE,
    DISABLED_CODE,
    SAMPLE_REACTION_CLASS,
)


class SynthesisManager:
    """Utility object to manage general saving, loading, and analysis
    functions for reaction synthesis steps.
    """

    # If a dictionary is provided, when doing reagent lookups, etc.
    #   will first look for them in the cache, instead of querying the database.
    #   will subsequently store any queried reagents in the cache as well, to save future query times.
    objectCache = None

    # Can specify a source to produce DB connection objects other than the default
    connFactory = None

    def __init__(self):
        """Initialization constructor.  Parameters to specify"""
        self.objectCache = None
        self.connFactory = None

    def loadReactionSynthesisIDs(self, synthesisRequest, cacheOnly=False):
        """Give the reaction synthesis request, just return the reaction_synthesis_ids
        for any matching syntheses.
        """
        # Query for reaction_synthesis_ids based upon the request
        query = SQLQuery()
        query.addSelect("synth.reaction_synthesis_id")
        query.addFrom("reaction_synthesis as synth")

        if cacheOnly:
            query.addWhereEqual("code", CACHE_CODE)
        if synthesisRequest.maxSteps is not None:
            query.addWhereEqual("max_steps", synthesisRequest.maxSteps)
        if len(synthesisRequest.reactionSynthesisIDs) > 0:
            query.addWhereIn(
                "synth.reaction_synthesis_id", synthesisRequest.reactionSynthesisIDs
            )

        if synthesisRequest.exampleReagentId is not None:
            query.addWhereEqual("synth.class", SAMPLE_REACTION_CLASS)
            # query.addWhereEqual("step.reagent_id", synthesisRequest.exampleReagentId); # Leave this for loadReactionSyntheses query to avoid bogging down an unoptimized query here
            query.addWhere("generation_date is null")
            # Do not include auto-generated reactions in demonstration examples

        if (
            synthesisRequest.reactionStepId is not None
            or synthesisRequest.reagentId is not None
            or synthesisRequest.reactantSmi is not None
        ):
            query.addFrom("reaction_step as step")
            query.addWhere("synth.reaction_synthesis_id = step.reaction_synthesis_id")
            if synthesisRequest.reactionStepId is not None:
                query.addWhereEqual(
                    "step.reaction_step_id", synthesisRequest.reactionStepId
                )
            if synthesisRequest.reactantSmi is not None:
                query.addWhereLike("step.reactant_smiles", synthesisRequest.reactantSmi)
            elif synthesisRequest.reagentId is not None:
                query.addWhereEqual("step.reagent_id", synthesisRequest.reagentId)

        if len(synthesisRequest.reactionCategoryIDs) > 0:
            if (
                not synthesisRequest.allowCategorySubsets
                and not synthesisRequest.allowCategorySupersets
            ):
                # Subsets not allowed, only an exact match will be accepted
                query.addWhereEqual(
                    "reaction_category_ids", synthesisRequest.reactionCategoryIDsStr()
                )
            else:
                # Subsets or supersets allowed.  With n-categories, could be 2^n possibilites to check for,
                #   but in practice, very few of these will actually exist in the database.
                #   Query from the direction of what's already in the database to determine
                #   which subsets are worth querying for.
                requestIdSet = synthesisRequest.getReactionCategoryIDSet()

                databaseIdSetData = self.loadSynthesisCategoryIdSetData()

                # To be DB safe, need to convert ID integers into strings to match reaction_category_ids text column
                idStrList = [synthesisRequest.reactionCategoryIDsStr()]
                # Start with exact matches to ensure a non-empty list
                for idSetStr, idSet in databaseIdSetData.iteritems():
                    if (
                        synthesisRequest.allowCategorySubsets
                        and idSet.issubset(requestIdSet)
                    ) or (
                        synthesisRequest.allowCategorySupersets
                        and idSet.issuperset(requestIdSet)
                    ):
                        idStrList.append(idSetStr)
                query.addWhereIn("reaction_category_ids", idStrList)

        reactionSynthesisIDTable = DBUtil.execute(query, connFactory=self.connFactory)

        reactionSynthesisIDs = []
        for row in reactionSynthesisIDTable:
            reactionSynthesisIDs.append(row[0])

        return reactionSynthesisIDs

    def loadReactionSyntheses(self, synthesisRequest, cacheOnly=False):
        """Given the reaction synthesis request,
        return the syntheses and the reaction steps for the respective syntheses.
        """
        reactionSynthesisIDs = self.loadReactionSynthesisIDs(
            synthesisRequest, cacheOnly
        )
        if len(reactionSynthesisIDs) < 1:
            return []
            # No IDs to search, nothing to return

        query = SQLQuery()
        query.addSelect("reaction_step_id")
        query.addSelect("reaction_synthesis_id")
        query.addSelect("reactant_smiles")
        query.addSelect("reagent_id")
        query.addSelect("product_smiles")
        query.addSelect("warning_level_ids")
        query.addSelect("warning_messages")
        query.addFrom("reaction_step as step")
        query.addWhereIn("reaction_synthesis_id", reactionSynthesisIDs)
        query.addOrderBy("reaction_synthesis_id")
        # Must order by this first for subsequent loop to work correctly
        query.addOrderBy("position")
        query.addOrderBy("reaction_step_id")
        # Ensure a consistent ordering
        if synthesisRequest.exampleReagentId is not None:
            query.addWhereEqual("reagent_id", synthesisRequest.exampleReagentId)

        reactionStepTable = DBUtil.execute(
            query, includeColumnNames=True, connFactory=self.connFactory
        )
        reactionStepModels = modelListFromTable(reactionStepTable)

        # Need to replace reagent_ids with actual reagent objects
        reagentIDs = []
        for reactionStepModel in reactionStepModels:
            reagentIDs.append(reactionStepModel["reagent_id"])

        reagentList = []
        if len(reagentIDs) > 0:
            reagentManager = ReagentManager()
            reagentManager.connFactory = self.connFactory

            reagentQuery = ReagentSearchModel()
            reagentQuery.reagentIDs = reagentIDs
            reagentQuery.shallowReagents = synthesisRequest.shallowReagents

            reagentList = reagentManager.loadDBInstances(reagentQuery)

        reagentById = dict()
        for reagent in reagentList:
            reagent_id = reagent["reagent_id"]
            reagentById[reagent_id] = reagent

        # Build up reactions steps for each distinct synthesis
        synthesisList = []
        lastSynthesisId = None
        reactionSteps = None
        # Swap in reagent objects, and convert composite SMILES into lists of molecule objects
        #   with reaction_synthesis_id data tags on the reactant molecules
        for iStep, reactionStepModel in enumerate(reactionStepModels):
            reactants = reactionStepModel["reactant_smiles"].split(MOL_LIST_DELIM)
            for iReactant, reactant in enumerate(reactants):
                reactantMol = molBySmiles(reactant)
                # Tag objects with IDs to facilitate subsequent traceback
                reactantMol.SetIntData(
                    "reaction_step_id", reactionStepModel["reaction_step_id"]
                )
                reactantMol.SetIntData(
                    "reaction_synthesis_id", reactionStepModel["reaction_synthesis_id"]
                )
                reactants[iReactant] = reactantMol

            reagent = reagentById[reactionStepModel["reagent_id"]]

            # Load products, but also tack on any warning messages
            products = reactionStepModel["product_smiles"].split(MOL_LIST_DELIM)
            for iProduct, product in enumerate(products):
                products[iProduct] = molBySmiles(product)
            if (
                reactionStepModel["warning_level_ids"] is not None
                and len(reactionStepModel["warning_level_ids"]) > 0
            ):
                warning_level_ids = reactionStepModel["warning_level_ids"].split(
                    REACTION_PROFILE_DELIM
                )
                warning_messages = reactionStepModel["warning_messages"].split(
                    REACTION_PROFILE_DELIM
                )

                for product, warning_level_id, warning_message in zip(
                    products, warning_level_ids, warning_messages
                ):
                    product.SetIntData("warning_level_id", int(warning_level_id))
                    product.SetStringData("warning_message", warning_message)

            reactionStep = [reactants, reagent, products]

            if lastSynthesisId is None:
                # First synthesis, prepare a blank set of reactionSteps to fill in
                reactionSteps = []
            elif reactionStepModel["reaction_synthesis_id"] != lastSynthesisId:
                # New set of synthesis steps, add the previous set to the final list
                # and prepare a new blank set of reaction steps
                synthesisList.append(reactionSteps)
                reactionSteps = []

            lastSynthesisId = reactionStepModel["reaction_synthesis_id"]
            reactionSteps.append(reactionStep)

        if reactionSteps is not None:
            synthesisList.append(reactionSteps)

        # If userID is specified in request parameters, check which products
        #   the user has already seen and avoid reproducing them
        if synthesisRequest.userID is not None:
            productSmilesSeenBeforeSet = set()
            alreadySeenTable = DBUtil.execute(
                """SELECT product_smiles 
                    FROM problem_record 
                    WHERE user_id = %(placeholder)s
                    AND completion_date is not null
                    """
                % {"placeholder": DBUtil.SQL_PLACEHOLDER},
                (synthesisRequest.userID,),
                connFactory=self.connFactory,
            )
            for row in alreadySeenTable:
                productSmilesSeenBeforeSet.add(row[0])

            filteredList = []
            for preCachedSynthesis in synthesisList:
                reactionSteps = preCachedSynthesis
                productSeenBefore = False
                for productMol in reactionSteps[-1][
                    PRODUCTS
                ]:  # Check final product list for products already seen / solved for
                    productSmiles = createStandardSmiString(productMol)
                    productSeenBefore = (
                        productSeenBefore or productSmiles in productSmilesSeenBeforeSet
                    )
                if not productSeenBefore:
                    filteredList.append(reactionSteps)
            if not synthesisRequest.limitedUserFilter or len(filteredList) > 0:
                # Replace the primary list with the user filtered list.  But if there is nothing left in the
                #   filtered list and we are using the limited option, skip this step to use the complete available list
                synthesisList = filteredList

        return synthesisList

    def retrieveReactionStepsBySynthesisId(self, reaction_synthesis_id):
        """Look for reaction steps (synthesis) data by the given ID in the database.
        If have an objectCache, then look there first and store any
        new copies there to save time later.
        """
        reactionSteps = None
        reactionStepsBySynthesisId = None

        # Check for availability in cache
        if self.objectCache is not None:
            if "reactionStepsBySynthesisId" not in self.objectCache:
                self.objectCache["reactionStepsBySynthesisId"] = dict()
            reactionStepsBySynthesisId = self.objectCache["reactionStepsBySynthesisId"]

            if reaction_synthesis_id in reactionStepsBySynthesisId:
                # Already exists in cache, just return it directly
                reactionSteps = reactionStepsBySynthesisId[reaction_synthesis_id]

        # Load from the database if not available in the cache
        if reactionSteps is None:
            synthesisRequest = SynthesisRequestModel()
            synthesisRequest.reactionSynthesisIDs = [reaction_synthesis_id]
            reactionStepsList = self.loadReactionSyntheses(synthesisRequest)
            reactionSteps = []
            if len(reactionStepsList) > 0:
                reactionSteps = reactionStepsList[0]

            if self.objectCache is not None:
                # Store the newly loaded copy in the cache while we're here to save time later
                reactionStepsBySynthesisId[reaction_synthesis_id] = reactionSteps

        return reactionSteps

    def saveReactionSynthesis(
        self,
        reactionSteps,
        synthesisRequest=None,
        code=None,
        problemClass=None,
        conn=None,
    ):
        """Save the given reactions steps as a reaction synthesis in the database.
        If the synthesisRequest is supplied, then also save the generation parameters
        for subsequent retrieval filter.  Return the reaction_synthesis_id generated.
        """
        # Use one connection for transactional integrity
        extConn = conn is not None
        if not extConn:
            if self.connFactory is not None:
                conn = self.connFactory.connection()
            else:
                conn = DBUtil.connection()

        try:
            reaction_category_ids = None
            max_steps = None

            # Setup default values if none specified
            if code is None:
                code = CACHE_CODE

            if synthesisRequest is not None:
                max_steps = synthesisRequest.maxSteps
                reaction_category_ids = synthesisRequest.reactionCategoryIDsStr()
                if reaction_category_ids == "":
                    reaction_category_ids = None

            generation_date = datetime.now()

            DBUtil.execute(
                """INSERT INTO reaction_synthesis
                (code, reaction_category_ids, max_steps, generation_date, class) 
                VALUES
                (%(placeholder)s, %(placeholder)s, %(placeholder)s, %(placeholder)s, %(placeholder)s)
                """
                % {"placeholder": DBUtil.SQL_PLACEHOLDER},
                (code, reaction_category_ids, max_steps, generation_date, problemClass),
                conn=conn,
            )

            reaction_synthesis_id = DBUtil.execute(
                DBUtil.identityQuery("reaction_synthesis"), conn=conn
            )[0][0]

            for iStep, reactionStep in enumerate(reactionSteps):
                reactant_smiles = str.join(
                    MOL_LIST_DELIM, molToSmiList(reactionStep[REACTANTS])
                )
                reagent_id = reactionStep[REAGENT]["reagent_id"]
                product_smiles = str.join(
                    MOL_LIST_DELIM, molToSmiList(reactionStep[PRODUCTS])
                )

                warning_level_ids = []
                warning_messages = []
                for product in reactionStep[PRODUCTS]:
                    warning_level_ids.append(
                        str(product.GetIntData("warning_level_id"))
                    )
                    warning_messages.append(product.GetStringData("warning_message"))
                warning_level_ids = str.join(REACTION_PROFILE_DELIM, warning_level_ids)
                warning_messages = str.join(REACTION_PROFILE_DELIM, warning_messages)

                DBUtil.execute(
                    """INSERT INTO reaction_step
                    (   reaction_synthesis_id, 
                        reactant_smiles, 
                        reagent_id, 
                        product_smiles, 
                        warning_level_ids,
                        warning_messages,
                        position
                    )
                    VALUES
                    (   %(placeholder)s, %(placeholder)s, %(placeholder)s, 
                        %(placeholder)s, %(placeholder)s, %(placeholder)s, %(placeholder)s
                    )
                    """
                    % {"placeholder": DBUtil.SQL_PLACEHOLDER},
                    (
                        reaction_synthesis_id,
                        reactant_smiles,
                        reagent_id,
                        product_smiles,
                        warning_level_ids,
                        warning_messages,
                        iStep,
                    ),
                    conn=conn,
                )

                reaction_step_id = DBUtil.execute(
                    DBUtil.identityQuery("reaction_step"), conn=conn
                )[0][0]

                # Reconstruct and store reaction_step_reaction_profile records
                reactionProfileIdSet = self.assignStepReactionProfileIds(
                    reactant_smiles, reagent_id, conn=conn
                )
                self.updateReactionStepReactionProfile(
                    reaction_step_id, reactionProfileIdSet, conn=conn
                )

            conn.commit()

            try:
                # Let any object caches know that the synthesis table data has changed
                self.objectCache.recordDataUpdate("synthesisReactionCategoryIdSetData")
            except:
                # Apparently this is not a DatabaseTimestampCache.  Forget it then, nothing to record.
                pass

        finally:
            if not extConn:
                # Done with this connection object
                conn.close()
        return reaction_synthesis_id

    def deleteReactionSynthesis(self, clearIDs, conn=None, softDelete=False):
        """Delete the identified reaction synthesis records
        including any dependent child records.

        softDelete:  If true, will not actual remove the DB cache records, will just
                    relabel the cached problem codes to DISABLED_CODE.
        """
        if len(clearIDs) < 1:
            # Nothing to do
            return

        # Use a single connection to ensure transactional integrity
        if conn is None:
            connFactory = self.connFactory
            if connFactory is None:
                # Default connection factory if none specified
                connFactory = DBUtil.ConnectionFactory()
            conn = connFactory.connection()

        cursor = conn.cursor()
        try:
            clearCount = len(clearIDs)
            placeholders = str.join(",", [DBUtil.SQL_PLACEHOLDER] * clearCount)
            clearIDs = tuple(clearIDs)
            # So it will work with DB parameter replacement

            if not softDelete:
                # Do a real DB delete to remove the records
                query = (
                    """DELETE FROM reaction_step_reaction_profile 
                    WHERE reaction_step_id IN
                    (   SELECT reaction_step_id
                        FROM reaction_step
                        WHERE reaction_synthesis_id IN (%s)
                    )"""
                    % placeholders
                )
                cursor.execute(query, clearIDs)

                query = (
                    """DELETE FROM reaction_step WHERE reaction_synthesis_id IN (%s)"""
                    % placeholders
                )
                cursor.execute(query, clearIDs)

                query = (
                    """DELETE FROM reaction_synthesis WHERE reaction_synthesis_id IN (%s)"""
                    % placeholders
                )
                cursor.execute(query, clearIDs)
            else:
                # Just "soft delete" the problems by relabeling them
                query = (
                    """UPDATE reaction_synthesis SET code = '%s' WHERE reaction_synthesis_id IN (%s)"""
                    % (DISABLED_CODE, placeholders)
                )
                cursor.execute(query, clearIDs)

            conn.commit()
            # Only commit after both deletes were completed together as a single transaction
        finally:
            cursor.close()
            conn.close()

    def cloneReactionSyntheses(self, synthesisIds):
        """Retrieve the synthesis / reaction lists identified by the given IDs.
        Persist copies of these records with revised code / title to indicate the difference.
        """
        # Load the detailed reaction step data
        synthesisRequest = SynthesisRequestModel()
        synthesisRequest.reactionSynthesisIDs = synthesisIds
        synthesisRequest.shallowReagents = True
        synthesisList = self.loadReactionSyntheses(synthesisRequest)

        # Common connection for transactional integrity and efficiency
        conn = None
        if self.connFactory is not None:
            conn = self.connFactory.connection()
        else:
            conn = DBUtil.connection()

        try:
            # Load the descriptive / metadata for the synthesis / reaction lists
            query = SQLQuery()
            query.addSelect("*")
            query.addFrom("reaction_synthesis")
            query.addWhereIn("reaction_synthesis_id", synthesisIds)
            resultTable = DBUtil.execute(query, includeColumnNames=True, conn=conn)
            resultModels = modelListFromTable(resultTable)
            modelsById = modelDictFromList(resultModels, "reaction_synthesis_id")

            for reactionSteps in synthesisList:
                if len(reactionSteps) > 0 and len(reactionSteps[0][REACTANTS]) > 0:
                    oldSynthesisId = reactionSteps[0][REACTANTS][0].GetIntData(
                        "reaction_synthesis_id"
                    )
                    oldSynthesisModel = modelsById[oldSynthesisId]

                    # Clone a copy of the detailed reaction steps
                    newSynthesisId = self.saveReactionSynthesis(
                        reactionSteps, conn=conn
                    )

                    # Prepare a revised version of the problem data code / title for the clone
                    newSynthesisModel = dict(oldSynthesisModel)
                    del newSynthesisModel["reaction_synthesis_id"]
                    # Keep the new ID
                    del newSynthesisModel["generation_date"]
                    # Keep the new date
                    newSynthesisModel["code"] = (
                        "%(code)s (CLONED from #%(reaction_synthesis_id)s)"
                        % oldSynthesisModel
                    )

                    # Update the the clone with respective descriptive metadata
                    DBUtil.updateRow(
                        "reaction_synthesis",
                        newSynthesisModel,
                        newSynthesisId,
                        conn=conn,
                    )

            conn.commit()
        finally:
            conn.close()

    def loadSynthesisCategoryIdSetData(self):
        """Return a dictionary listing all of the reaction_category_ids strings and sets
        that currently existing in the database of reaction_synthesis records.

        Save this result to the objectCache if available so subsequent requests can
        be fulfilled without hitting the database.  Possibility of cached data list
        becoming stale, so methods which could alter the table contents should use the
        DatabaseTimestampCache.recordDataUpdate method to signal changes.

        Dictionary items keyed by the reaction_category_ids string found in the database,
        value is a set of int representations of each string.
        """
        objectCache = self.objectCache
        if objectCache is None:
            # Just build a temporary cache dictionary to populate for now then
            objectCache = dict()

        if "synthesisReactionCategoryIdSetData" not in objectCache:
            # Query for reaction_category_ids actually used in existing synthesis records
            query = SQLQuery()
            query.addSelect("reaction_category_ids")
            query.addFrom("reaction_synthesis")
            query.addWhere("reaction_category_ids is not null")
            query.addGroupBy("reaction_category_ids")
            categoryIdSetTable = DBUtil.execute(query, connFactory=self.connFactory)

            # Use a synthesis request object to consistently parse out the id set string
            synthesisRequest = SynthesisRequestModel()

            idSetData = dict()
            for (idSetStr,) in categoryIdSetTable:
                synthesisRequest.reactionCategoryIDs = idSetStr.split(MOL_LIST_DELIM)
                idSet = synthesisRequest.getReactionCategoryIDSet()

                idSetData[idSetStr] = idSet

            objectCache["synthesisReactionCategoryIdSetData"] = idSetData

        return objectCache["synthesisReactionCategoryIdSetData"]

    def moveStepPosition(self, reactionStepId, direction):
        """Find all of the synthesis / reaction list the given reaction step is under.
        Swap the values of the reaction step's position field with the next step adjacent to it,
        to effectively change the step order in the list.
        The direction parameter is just a positive or negative value, indicating whether this
        step should be moved / swapped further or sooner in the list, respectively.
        """
        # Create a single connection for transactional integrity
        connFactory = self.connFactory
        if connFactory is None:
            connFactory = DBUtil.ConnectionFactory()
            # Use default connection source

        conn = connFactory.connection()
        try:
            # Query the database for the step order / position information
            stepDataTable = DBUtil.execute(
                """select reaction_step_id, position
                    from reaction_step
                    where reaction_synthesis_id = 
                    (   select reaction_synthesis_id
                        from reaction_step
                        where reaction_step_id = %s
                    )
                    order by position
                    """
                % DBUtil.SQL_PLACEHOLDER,
                (reactionStepId,),
                conn=conn,
            )

            # Find the specified step in the list
            nSteps = len(stepDataTable)
            for iStep, (stepId, position) in enumerate(stepDataTable):
                if stepId == reactionStepId:
                    # Specify which steps to swap positions
                    (stepIdA, positionA) = (None, None)
                    (stepIdB, positionB) = (None, None)

                    if direction < 0 and iStep > 0:
                        # Move the step sooner in the list
                        (stepIdA, positionA) = tuple(stepDataTable[iStep - 1])
                        (stepIdB, positionB) = tuple(stepDataTable[iStep + 0])

                    elif direction > 0 and iStep < nSteps:
                        # Move the step farther along the list
                        (stepIdA, positionA) = tuple(stepDataTable[iStep + 0])
                        (stepIdB, positionB) = tuple(stepDataTable[iStep + 1])

                    if stepIdA is not None:
                        # Execute the swap by updating the database
                        updateQuery = (
                            """update reaction_step set position = %(p)s where reaction_step_id = %(p)s"""
                            % {"p": DBUtil.SQL_PLACEHOLDER}
                        )
                        DBUtil.execute(updateQuery, (positionB, stepIdA), conn=conn)
                        DBUtil.execute(updateQuery, (positionA, stepIdB), conn=conn)
                        conn.commit()

                    break
                    # Do not need to look further
        finally:
            conn.close()

    def reassignSteps(self, reactionStepIdList, reassignListId):
        """Reassign the selected reaction step records to another reaction_synthesis (reaction list).
        Update the step position numbers to be unique and greater than any steps under the existing list
        to avoid duplicates.
        """
        # Create a single connection for transactional integrity
        connFactory = self.connFactory
        if connFactory is None:
            connFactory = DBUtil.ConnectionFactory()
            # Use default connection source

        conn = connFactory.connection()
        try:
            # Query the database for position under the existing synthesis / list
            currentMaxPosition = DBUtil.execute(
                """select max(position)
                    from reaction_step
                    where reaction_synthesis_id = %s
                    """
                % DBUtil.SQL_PLACEHOLDER,
                (reassignListId,),
                conn=conn,
            )[0][0]

            if currentMaxPosition is None:
                currentMaxPosition = 0

            updateQuery = """update reaction_step
                set reaction_synthesis_id = %(p)s,
                    position = %(p)s
                where reaction_step_id = %(p)s
                """ % {
                "p": DBUtil.SQL_PLACEHOLDER
            }

            completedStepIds = set()
            # Keep track of completion to avoid repeats
            for reactionStepId in reactionStepIdList:
                if reactionStepId not in completedStepIds:
                    currentMaxPosition += 10
                    # Arbitrary spacer number to expand out position values
                    DBUtil.execute(
                        updateQuery,
                        (reassignListId, currentMaxPosition, reactionStepId),
                        conn=conn,
                    )
                    completedStepIds.add(reactionStepId)

            conn.commit()
        finally:
            conn.close()

    def minimizeSynthesisId(self, reaction_synthesis_id):
        """Given the ID for a reaction synthesis record stored in the database,
        alter the ID of the synthesis record (and any underlying reaction_step records)
        to the smallest (positive) value currently unused in the database.

        Beware that this will only work for a new, unused synthesis problem.
        Once users start recording problem records against it or classes start adding
        it to their assignments, dependent records will be created that need to maintain
        their foreign key constraints.  Could auto-correct all of those child records
        as well, but may be more trouble than it is worth.  Bad idea to change ID
        values ones an object has entered regular use anyway.
        """
        # Create a single connection for transactional integrity
        connFactory = self.connFactory
        if connFactory is None:
            connFactory = DBUtil.ConnectionFactory()
            # Use default connection source

        minAvailableId = None

        conn = connFactory.connection()
        try:
            minAvailableId = DBUtil.smallestAvailableValue(
                "reaction_synthesis", "reaction_synthesis_id", conn=conn
            )

            # Temporarily move child records (reaction_steps) under another synthesis record
            # Necessary step to adhere to foreign key constraint.
            # Need to keep track of records moved
            reactionStepIdTable = DBUtil.execute(
                """select reaction_step_id
                        from reaction_step
                        where reaction_synthesis_id = %s
                        """
                % DBUtil.SQL_PLACEHOLDER,
                (reaction_synthesis_id,),
                conn=conn,
            )
            reactionStepIds = list()
            for row in reactionStepIdTable:
                stepId = row[0]
                reactionStepIds.append(str(stepId))

            if len(reactionStepIds) > 0:
                # We're assuming a placeholder synthesis with ID (-1) exists.
                placeholderSynthesisId = -1
                DBUtil.execute(
                    """update reaction_step
                        set reaction_synthesis_id = %s
                        where reaction_step_id in (%s)
                        """
                    % (DBUtil.SQL_PLACEHOLDER, str.join(",", reactionStepIds)),
                    (placeholderSynthesisId,),
                    conn=conn,
                )

            # Update the primary synthesis record to the available ID value
            DBUtil.execute(
                """update reaction_synthesis
                    set reaction_synthesis_id = %s
                    where reaction_synthesis_id = %s
                    """
                % (DBUtil.SQL_PLACEHOLDER, DBUtil.SQL_PLACEHOLDER),
                (minAvailableId, reaction_synthesis_id),
                conn=conn,
            )

            # Reconnect the child records to the primary synthesis record
            if len(reactionStepIds) > 0:
                DBUtil.execute(
                    """update reaction_step
                        set reaction_synthesis_id = %s
                        where reaction_step_id in (%s)
                        """
                    % (DBUtil.SQL_PLACEHOLDER, str.join(",", reactionStepIds)),
                    (minAvailableId,),
                    conn=conn,
                )

            # Now minimize all of the child record IDs
            for reactionStepIdStr in reactionStepIds:
                reactionStepId = int(reactionStepIdStr)
                self.minimizeStepId(reactionStepId, conn=conn)

        finally:
            conn.close()

        return minAvailableId

    def minimizeStepId(self, reaction_step_id, conn):
        """Called within minimizeSynthesisId for a comparable effect,
        just on the reaction_step table.
        Require external connection object supplied so caller will be responsible
        for maintaining transactional integrity.
        """
        minAvailableId = DBUtil.smallestAvailableValue(
            "reaction_step", "reaction_step_id", conn=conn
        )

        # Temporarily move child records under another record
        # Necessary step to adhere to foreign key constraint.
        # Need to keep track of records moved
        linkProfileIdTable = DBUtil.execute(
            """select reaction_step_reaction_profile_id
                    from reaction_step_reaction_profile
                    where reaction_step_id = %s
                    """
            % DBUtil.SQL_PLACEHOLDER,
            (reaction_step_id,),
            conn=conn,
        )
        linkProfileIds = list()
        for row in linkProfileIdTable:
            linkProfileId = row[0]
            linkProfileIds.append(str(linkProfileId))

        if len(linkProfileIds) > 0:
            # Must ensure a placeholder reaction step (with ID -1) exists
            placeholderId = -1
            placeholderModel = {
                "reaction_step_id": placeholderId,
                "reaction_synthesis_id": placeholderId,
            }
            DBUtil.findOrInsertItem("reaction_step", placeholderModel, conn=conn)

            DBUtil.execute(
                """update reaction_step_reaction_profile
                    set reaction_step_id = %s
                    where reaction_step_reaction_profile_id in (%s)
                    """
                % (DBUtil.SQL_PLACEHOLDER, str.join(",", linkProfileIds)),
                (placeholderId,),
                conn=conn,
            )

        # Update the primary record to the available ID value
        DBUtil.execute(
            """update reaction_step
                set reaction_step_id = %s
                where reaction_step_id = %s
                """
            % (DBUtil.SQL_PLACEHOLDER, DBUtil.SQL_PLACEHOLDER),
            (minAvailableId, reaction_step_id),
            conn=conn,
        )

        # Reconnect the child records to the primary record
        if len(linkProfileIds) > 0:
            DBUtil.execute(
                """update reaction_step_reaction_profile
                    set reaction_step_id = %s
                    where reaction_step_reaction_profile_id in (%s)
                    """
                % (DBUtil.SQL_PLACEHOLDER, str.join(",", linkProfileIds)),
                (minAvailableId,),
                conn=conn,
            )

        return minAvailableId

    def assignStepReactionProfileIds(self, reactantSmiles, reagentId, conn=None):
        """Assign / perceive the reaction profiles used for an individual
        reaction step by applying a designated reagent model to
        a list of input starting materials.
        Reactant SMILES input expected with MOL_LIST_DELIM to separate multiple inputs.
        """
        reactionProfileIdSet = set()

        reagentManager = ReagentManager()
        reagentManager.objectCache = self.objectCache
        reagentManager.connFactory = self.connFactory

        reagentsById = reagentManager.retrieveReagentsById([reagentId], conn=conn)

        reagent = reagentsById[reagentId]
        reagent.removeInherentProducts = True
        # Not interested in inherent products
        reagent.ignoreSelfReactions = True
        supplementalData = SupplementalDataModel()
        # Optional object to receive extra information
        reactants = [molBySmiles(smi) for smi in reactantSmiles.split(MOL_LIST_DELIM)]
        # List of molecule objects as reactant input
        products = reagent(reactants, supplementalData)
        # Use reagent as function to generate products
        reactionStep = (reactants, reagent, products)
        # Combined elements of a complete reaction
        for elementaryStep in supplementalData.filteredElementarySteps(
            reactants, products
        ):
            (detailReactants, reactionProfile, detailProducts) = elementaryStep
            reactionProfileIdSet.add(reactionProfile["reaction_profile_id"])

        reagentManager.returnReagentsById(reagentsById)
        # Return used reagent objects to cache for later reuse

        return reactionProfileIdSet

    def assignSynthesisReactionProfileIds(
        self, reactionSynthesisId, checkExisting=True, conn=None
    ):
        """Assign the reaction profiles used in the course of the
        identified synthesis.  This may be a random or non-random synthesis,
        as long as a respective set of database records exists for it.
        Returns a list of sets, one set per step used in the synthesis.

        checkExisting:
            Option to indicate whether to first check if records already
            exist under reaction_step_reaction_profile.  If so, just pick those up
            instead of reconstructing them.
            If no such records exist, or this option is turned off, force the
            system to reconstruct the results.
        """
        reactionProfileIdSetList = list()

        if checkExisting:
            # Optimistic search, the synthesis may already have the linked reaction profiles pre-computed
            grandChildQuery = (
                """select 
                    rs.reaction_step_id, 
                    rs.position,
                    rsrp.reaction_profile_id
                from 
                    reaction_step_reaction_profile rsrp,
                    reaction_step as rs
                where 
                    rsrp.reaction_step_id = rs.reaction_step_id and
                    rs.reaction_synthesis_id = %s
                order by
                    rs.position,
                    rsrp.reaction_profile_id
                """
                % DBUtil.SQL_PLACEHOLDER
            )
            linkedData = DBUtil.execute(
                grandChildQuery,
                (reactionSynthesisId,),
                conn=conn,
                connFactory=self.connFactory,
            )

            # Build up lists for each reaction step
            stepIdsFound = set()
            for stepId, position, ruleId in linkedData:
                if stepId not in stepIdsFound:
                    # Found a new reaction step.  Prepare a new, empty Id set
                    stepIdsFound.add(stepId)
                    reactionProfileIdSetList.append(set())

                ruleIdSet = reactionProfileIdSetList[-1]
                ruleIdSet.add(ruleId)

        if len(reactionProfileIdSetList) < 1:
            # Have not found any data yet, try to reconstruct it then
            reactionStepData = DBUtil.execute(
                """SELECT reaction_step_id, reagent_id, reactant_smiles 
                    FROM reaction_step 
                    WHERE reaction_synthesis_id = %d
                    ORDER BY position
                    """
                % reactionSynthesisId,
                conn=conn,
                connFactory=self.connFactory,
            )

            for stepId, reagentId, reactantSmiles in reactionStepData:
                stepReactionProfileIds = self.assignStepReactionProfileIds(
                    reactantSmiles, reagentId, conn
                )
                reactionProfileIdSetList.append(stepReactionProfileIds)

        return reactionProfileIdSetList

    def updateReactionSynthesisReactionProfile(
        self, reactionSynthesisId, reactionProfileIdSetList, conn=None
    ):
        """Given reactionSynthesisId and a list of sets of reaction profiles,
        add in the relevant rows to the reaction_step_reaction_profile table.
        Assumes that the list of ID sets is in the same order as the
        reaction step records under the synthesis (ordered by the "position" column).
        """
        reactionStepData = DBUtil.execute(
            """SELECT reaction_step_id 
                FROM reaction_step 
                WHERE reaction_synthesis_id = %d
                ORDER BY position
                """
            % reactionSynthesisId,
            conn=conn,
            connFactory=self.connFactory,
        )
        reactionStepIdList = list()
        for (stepId,) in reactionStepData:
            reactionStepIdList.append(stepId)

        for stepId, reactionProfileIdSet in zip(
            reactionStepIdList, reactionProfileIdSetList
        ):
            self.updateReactionStepReactionProfile(
                stepId, reactionProfileIdSet, conn=conn
            )

        # Commit changes to database so we don't leave too much hanging in the buffer
        if conn is not None:
            conn.commit()

    def updateReactionStepReactionProfile(
        self, reactionStepId, reactionProfileIdSet, conn=None
    ):
        """Given reactionStepId and a set of reaction profiles,
        add in the relevant rows to the reaction_step_reaction_profile table.
        """
        cols = [
            "reaction_step_id",
            "reaction_profile_id",
        ]

        insertQuery = DBUtil.buildInsertQuery("reaction_step_reaction_profile", cols)

        for reactionProfileId in reactionProfileIdSet:
            params = (reactionStepId, reactionProfileId)
            DBUtil.execute(insertQuery, params, conn=conn, connFactory=self.connFactory)

    def clearProblemCacheByReactants(self, reactionCategoryId, reactantSmiSet):
        """Find any problems in the saved DB cache that include this reaction category
        and use any of the included reactant SMILES and remove / delete them from the cache.
        This is important after removing reactants from a reaction category, as otherwise
        the cache will still have "stale" problems that include the use of those reactants.
        """
        synthesisRequest = SynthesisRequestModel()
        synthesisRequest.reactionCategoryIDs = [reactionCategoryId]
        synthesisRequest.allowCategorySupersets = True
        # Any problem which includes this category is suspect

        deleteIds = set()
        for reactantSmi in reactantSmiSet:
            # Add flanking wildcards to accomodate substring search
            # This may be too broad, with some false positive deletes, but that's okay,
            #   these are just temporary cache problems.  Safer to be more aggressive in their removal,
            #   rather than leaving false negatives (stale problems) behind.
            synthesisRequest.reactantSmi = "%" + reactantSmi + "%"

            deleteIds.update(
                self.loadReactionSynthesisIDs(synthesisRequest, cacheOnly=True)
            )

        self.deleteReactionSynthesis(deleteIds, softDelete=True)

        return deleteIds

    def clearProblemCacheByReagents(self, reactionCategoryId, reagentIds):
        """Find any problems in the saved DB cache that include this reaction category
        and use any of the included reagents and remove / delete them from the cache.
        This is important after removing reagents from a reaction category, as otherwise
        the cache will still have "stale" problems that include the use of those reagents.
        """
        synthesisRequest = SynthesisRequestModel()
        synthesisRequest.reactionCategoryIDs = [reactionCategoryId]
        synthesisRequest.allowCategorySupersets = True
        # Any problem which includes this category is suspect

        deleteIds = set()
        for reagentId in reagentIds:
            synthesisRequest.reagentId = reagentId

            deleteIds.update(
                self.loadReactionSynthesisIDs(synthesisRequest, cacheOnly=True)
            )

        self.deleteReactionSynthesis(deleteIds, softDelete=True)

        return deleteIds

    def clearProblemCacheByCategoryId(self, reactionCategoryId):
        """Find any problems in the saved DB cache that include this reaction category
        and remove / delete them from the cache.
        This can be important after modifying reaction category contents, as otherwise
        the cache will still have "stale" problems that include the use of the previous content.
        """
        synthesisRequest = SynthesisRequestModel()
        synthesisRequest.reactionCategoryIDs = [reactionCategoryId]
        synthesisRequest.allowCategorySupersets = True
        # Any problem which includes this category is suspect

        deleteIds = set()
        deleteIds.update(
            self.loadReactionSynthesisIDs(synthesisRequest, cacheOnly=True)
        )

        self.deleteReactionSynthesis(deleteIds, softDelete=True)

        return deleteIds


class SynthesisRequestModel:
    """Simple struct to contain query parameters indicating
    which reaction synthesis to load from the database.
    """

    # List of specific reaction synthesis IDs to load
    reactionSynthesisIDs = None

    # List of reaction category IDs under which the synthesis was generated
    reactionCategoryIDs = None

    # ID of a specific reaction step whose overall synthesis is being searched for
    reactionStepId = None

    # ID of a specific reagent to look for example reactions for
    exampleReagentId = None

    # ID of specific reagent to look for any syntheses that include use of the reagent
    reagentId = None

    # SMILES strings to look for syntheses that include one of these in
    #   one of their steps' reactant / starting materials.
    # Searched with the "like" operator, so % wildcards can be added for substring searches.
    reactantSmi = None

    # Maximum number of steps parameter used when generating the synthesis.
    # The actual number of steps in the synthesis could be less
    maxSteps = None

    # Maximum number of tries to make an optimal synthesis if end up
    #   needing to generate some dynamically.
    maxTries = MAX_TRIES

    # If specified, avoid returning a synthesis whose product this user has
    #   already solved in a past problem.
    userID = None

    # Limited filtering by user ID.  If user has somehow solved every possible problem available,
    #   don't return an empty list (will result in costly problem generation time), just return what is available.
    limitedUserFilter = True

    # If specified with a user ID, when retrieving problems, will filter the possible results down to this
    #   size, based on those deemed closest to the "fringe" of the user's knowledge state.
    knowledgeMapLimit = None

    # If set, will not load all of the reagents' reaction profiles, only the top-level shallow attributes
    # Efficient if you're only interested in loading the reagent for display purposes
    shallowReagents = True

    # If set, when request is used to retrieve existing (cached) syntheses,
    #   will not require exact category set matches, but allow subsets as well.
    # For example, if request a synthesis that covers reaction categories X,Y,Z
    #   accept any synthesis based on a non-empty subset of categories (X,Y,Z).
    allowCategorySubsets = False

    # Same as above, but for supersets
    allowCategorySupersets = False

    def __init__(self):
        self.reactionSynthesisIDs = list()
        self.reactionCategoryIDs = list()
        self.reactantSmiSet = set()

    def copyModiferAttributes(self, otherRequest):
        """Copy over the search parameter attributes from the other SynthesisRequestModel,
        that do not necessarily modify the SQL query that will be used to find the
        synthesis models, but which are used as modifiers in the loadReactionSyntheses
        method to alter (filter) a queried list of syntheses.
        """
        self.shallowReagents = otherRequest.shallowReagents
        self.exampleReagentId = otherRequest.exampleReagentId
        self.userID = otherRequest.userID
        self.limitedUserFilter = otherRequest.limitedUserFilter

    def reactionCategoryIDsStr(self):
        """Translate ID list into a consistent string representation for
        storage / lookup in the database.  Should eliminate redundancies
        down to a unique set, presented in sorted numerical order.
        """

        # Extract unique set of IDs as numerical (integer) values
        idSet = self.getReactionCategoryIDSet()

        # Copy to a list for sorting
        idList = list(idSet)
        idList.sort()

        # Convert to a list of strings for easy concatenation
        idStrList = list()
        for id in idList:
            idStrList.append(str(id))
        compositeIdStr = str.join(MOL_LIST_DELIM, idStrList)
        return compositeIdStr

    def getReactionCategoryIDSet(self):
        """Return reaction category IDs as a unique Set,
        and ensure that all are converted to ints
        """
        idSet = set()
        for reactionCategoryID in self.reactionCategoryIDs:
            idSet.add(int(reactionCategoryID))
        return idSet


if __name__ == "__main__":
    instance = SynthesisManager()
    instance.main(sys.argv)
