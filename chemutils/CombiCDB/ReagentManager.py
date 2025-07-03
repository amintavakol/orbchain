
# from mx.DateTime import now;
from chemutils.Common.Util import molToSmiList
from chemutils.Common.Util import loadObjectsFromQueuesById, saveObjectsToQueuesById
from chemutils.CombiCDB.ReactionModel import (
    SMIRKSReagent,
    SupplementalDataModel,
)
from chemutils.Common.Model import SQLQuery, modelListFromTable, modelDictFromList
from chemutils.Common import DBUtil

from chemutils.CombiCDB.Const import (
    MOL_LIST_DELIM,
    REACTION_PROFILE_DELIM,
)
from chemutils.CombiCDB.Const import LIST_SUFFIX

from chemutils.CombiCDB.Util import log


class ReagentSearchModel:
    """Simple struct containing query parameters indicating
    which reagents to pull out of the database.
    """

    # List of reagent IDs to select by
    reagentIDs = None

    # List of reaction category IDs under which the reagents to get must be listed
    reactionCategoryIDs = None

    # Whether to retrieve reagents marked as "retro ready" or not.  That is, those where the "retro" function is known to work properly
    retroReady = None

    # Whether to retrieve reagents marked as enabled for predictive reactions
    enabled = None

    # If set, will not load all of the reaction profiles, only the top-level shallow attributes
    # Efficient if you're only interested in loading the reagent for display purposes
    shallowReagents = False

    def __init__(self):
        """Constructor, setup default parameters"""
        self.reagentIDs = list()
        self.reactionCategoryIDs = list()


class ReagentManager:
    """Utility object to facilitate loading of reagent model objects,
    either from the database or from an in-memory cache.
    """

    # If a dictionary is provided, when doing reagent lookups, etc.
    #   will first look for them in the cache, instead of querying the database.
    #   will subsequently store any queried reagents in the cache as well, to save future query times.
    objectCache = None

    # If provided, will save any applied reagent results here, to reserve if get
    #   repeat requests.  Should probably use a TimedCache so an expiration time on entries
    #   so don't accumulate indefinitely.
    resultCache = None

    # Can specify a source to produce DB connection objects other than the default
    connFactory = None

    def __init__(self):
        """Initialization constructor.  Default components."""
        self.objectCache = dict()
        # Replace with external one as appropriate
        self.resultCache = None
        self.connFactory = DBUtil.ConnectionFactory()

    def loadDBInstances(self, reagentSearchModel=None, conn=None):
        """Return a list of reagent objects based on data loaded from the database.

        reagentSearchModel
            Simple struct containing various query parameters
        """
        extConn = conn is not None
        if not extConn:
            if self.connFactory is not None:
                conn = self.connFactory.connection()
            else:
                conn = DBUtil.connection()
                # Default connection source

        if reagentSearchModel is None:
            reagentSearchModel = ReagentSearchModel()
            # Just generate a new one with default options

        reagentIDs = reagentSearchModel.reagentIDs
        reactionCategoryIDs = reagentSearchModel.reactionCategoryIDs

        try:
            # Load all of the basic reagent data, filtered by the search criteria
            query = SQLQuery()
            query.addSelect("distinct r.*")
            query.addFrom("reagent as r")
            query.addOrderBy("r.position, r.description")
            if len(reagentIDs) > 0:
                query.addWhereIn("r.reagent_id", reagentIDs)
            if len(reactionCategoryIDs) > 0:
                query.addFrom("reaction_category_reagent as rcr")
                query.addWhere("r.reagent_id = rcr.reagent_id")
                query.addWhereIn("rcr.reaction_category_id", reactionCategoryIDs)
            if reagentSearchModel.enabled is not None:
                query.addWhereEqual("r.enabled", reagentSearchModel.enabled)
            if reagentSearchModel.retroReady is not None:
                query.addWhereEqual("r.retro_ready", reagentSearchModel.retroReady)

            reagentTable = DBUtil.execute(query, includeColumnNames=True, conn=conn)
            reagentModels = modelListFromTable(reagentTable)
            reagentsById = modelDictFromList(reagentModels, "reagent_id")

            if len(reagentModels) > 0 and not reagentSearchModel.shallowReagents:
                # Load the SMIRKS reaction profiles
                query = SQLQuery()
                query.addSelect("rp.*, rrp.*")
                query.addFrom("reagent_reaction_profile as rrp")
                query.addFrom("reaction_profile as rp")
                query.addWhere("rrp.reaction_profile_id = rp.reaction_profile_id")
                query.addOrderBy("rrp.reagent_id")
                query.addWhereIn("rrp.reagent_id", reagentsById.keys())

                reactionTable = DBUtil.execute(
                    query, includeColumnNames=True, conn=conn
                )
                reactionModels = modelListFromTable(reactionTable)
                nReactions = len(reactionModels)

                # Decrypt SMIRKS strings from database as necessary
                for reactionModel in reactionModels:
                    if SMIRKSReagent.smirksCodec.isEncrypted(reactionModel["smirks"]):
                        reactionModel["smirks"] = SMIRKSReagent.smirksCodec.decrypt(
                            reactionModel["smirks"]
                        )
                    if SMIRKSReagent.smirksCodec.isEncrypted(
                        reactionModel["mechanism_label_smirks"]
                    ):
                        reactionModel[
                            "mechanism_label_smirks"
                        ] = SMIRKSReagent.smirksCodec.decrypt(
                            reactionModel["mechanism_label_smirks"]
                        )

                # Track reaction models by identifying features for lookup later
                reactionModelsByKey = dict()
                for reactionModel in reactionModels:
                    reactionKey = SMIRKSReagent.reactionModelCompositeKey(reactionModel)
                    reactionModelsByKey[reactionKey] = reactionModel

                # Before loading the reaction profiles onto the primary reagent object
                #   recursively check for inherited reaction profiles
                # First find IDs of any inherited reagents
                inheritanceQuery = SQLQuery()
                inheritanceQuery.addSelect("reagent_id, inherited_reagent_id")
                inheritanceQuery.addFrom("reagent_inheritance")
                inheritanceQuery.addWhereIn("reagent_id", reagentsById.keys())
                inheritedReagentTable = DBUtil.execute(inheritanceQuery, conn=conn)
                if len(inheritedReagentTable) > 0:
                    # Now recursively load these reagent objects
                    inheritedReagentSearchModel = ReagentSearchModel()
                    inheritedReagentSearchModel.reagentIDs = []
                    for row in inheritedReagentTable:
                        inheritedReagentSearchModel.reagentIDs.append(row[1])
                    inheritedReagents = self.loadDBInstances(
                        inheritedReagentSearchModel, conn
                    )
                    inheritedReagentsById = modelDictFromList(
                        inheritedReagents, "reagent_id"
                    )
                    # Now add the reaction profiles from the inherited reagents to the primary reagents
                    for row in inheritedReagentTable:
                        primary_reagent_id = row[0]
                        inherited_reagent_id = row[1]
                        reagent = reagentsById[primary_reagent_id]
                        inherited_reagent = inheritedReagentsById[inherited_reagent_id]
                        for reactionProfile in inherited_reagent.reactionProfileList:
                            copyReactionProfile = dict(reactionProfile)
                            # Create a copy we can modify
                            copyReactionProfile["reagent_id"] = primary_reagent_id
                            # Alter the ownership of the copy from the inherited reagent to the primary (inheriting) reagent

                            # Check if an equivalent reactionProfile already exists in the primary reagent.
                            # If so, that one should "override" any inherited ones
                            inheritedReactionKey = (
                                SMIRKSReagent.reactionModelCompositeKey(
                                    copyReactionProfile
                                )
                            )
                            if inheritedReactionKey not in reactionModelsByKey:
                                reactionModels.append(copyReactionProfile)

                # Initialize reagents with blank data
                for reagent in reagentModels:
                    reagent["reaction_profile_id_list"] = []
                    reagent["reagent_reaction_profile_id_list"] = []
                    reagent["smirks_list"] = []
                    reagent["description_list"] = []
                    reagent["mechanism_label_smirks_list"] = []
                    reagent["mechanism_arrow_codes_list"] = []
                    reagent["priority_list"] = []
                    reagent["warning_level_id_list"] = []
                    reagent["warning_message_list"] = []
                    reagent["reaction_category_ids_list"] = []
                    reagent["pre_status_threshold_list"] = []
                    reagent["pre_status_minimum_list"] = []
                    reagent["post_status_list"] = []
                    reagent["retro_ready_list"] = []
                    reagent["intended_usage_list"] = []

                # Go through reaction profiles and match up to reagents to build smirks_list
                for reactionModel in reactionModels:
                    reagent_id = reactionModel["reagent_id"]

                    # Default values for optional parameters
                    if reactionModel["description"] is None:
                        reactionModel["description"] = ""
                    if reactionModel["mechanism_label_smirks"] is None:
                        reactionModel["mechanism_label_smirks"] = ""
                    if reactionModel["mechanism_arrow_codes"] is None:
                        reactionModel["mechanism_arrow_codes"] = ""

                    # Accumulate list attributes as strings
                    if reagent_id in reagentsById:
                        reagent = reagentsById[reagent_id]
                        for key in reactionModel.iterkeys():
                            listKey = key + LIST_SUFFIX
                            if listKey in reagent:
                                reagent[listKey].append(str(reactionModel[key]))

            # Transform the reagent RowItemModels into functional SMIRKSReagent objects
            smirksReagents = []
            for reagentModel in reagentModels:
                # Load sub-reaction profile data if it exists (not a shallow query)
                if "smirks_list" in reagentModel:
                    for key in reagentModel.iterkeys():
                        if key.endswith(LIST_SUFFIX):
                            reagentModel[key] = str.join(
                                REACTION_PROFILE_DELIM, reagentModel[key]
                            )
                smirksReagents.append(SMIRKSReagent(reagentModel))

            return smirksReagents

        finally:
            if not extConn:
                conn.close()

    def retrieveReagentsById(self, neededReagentIds, conn=None):
        """Given a list of reagent IDs, return a dictionary
        of reagent models, keyed by those IDs.  In general this could
        just be just like the loadDBInstances function from the database,
        but also have the possibility of retrieving from
        the provided object cache to reduce load times.
        """
        reagentsById = dict()

        # See what's available in the object cache, otherwise load a fresh copies from the database
        objectCache = self.objectCache
        if objectCache is None:
            objectCache = dict()
        if "reagentsById" not in objectCache:
            objectCache["reagentsById"] = dict()
        reagentsById.update(objectCache["reagentsById"])

        # Determine which reagents are not already in the cache
        neededIds = set(neededReagentIds)
        newReagentIds = neededIds.difference(reagentsById.keys())

        if len(newReagentIds) > 0:
            # Reagent was unavailable in the cache, will have to load it fresh from the database
            query = ReagentSearchModel()
            query.reagentIDs = newReagentIds
            query.shallowReagents = False
            reagents = self.loadDBInstances(query, conn=conn)
            for reagent in reagents:
                reagentsById[reagent["reagent_id"]] = reagent
        return reagentsById

    def retrieveRetroReagentsById(self, neededRetroReagentIds, conn=None):
        """Given a list of reagent IDs, return a dictionary
        of "retrified" reagent models, keyed by those IDs.  In general this could
        just be just like the loadDBInstances function from the database,
        but also have the possibility of retrieving from
        the provided object cache to reduce load times.

        NOTE:  These are keyed by the original reagent ids
        """
        log.debug("Attempting to load retro reagents from cache.")
        retroReagentsById = dict()
        # Container to track loaded object

        # See if it's available in the object cache, otherwise load a fresh copy from the database
        objectCache = self.objectCache
        if objectCache is None:
            objectCache = dict()
        if "retroReagentQueuesById" not in objectCache:
            objectCache["retroReagentQueuesById"] = dict()
        retroReagentQueuesById = objectCache["retroReagentQueuesById"]

        # Determine which reagents are already available and which are needed right now
        newRetroReagentIds = loadObjectsFromQueuesById(
            retroReagentsById, retroReagentQueuesById, neededRetroReagentIds
        )

        log.debug(
            "Loaded %d Retro Reagents From cache"
            % (len(neededRetroReagentIds) - len(newRetroReagentIds))
        )

        if len(newRetroReagentIds) > 0:
            log.debug("Loading %d retro reagents from the DB" % len(newRetroReagentIds))
            # Retro Reagent was unavailable in the cache, will have to load it fresh
            # from the database
            query = ReagentSearchModel()
            query.reagentIDs = newRetroReagentIds
            query.shallowReagents = False
            reagents = self.loadDBInstances(query, conn=conn)
            for reagent in reagents:
                retroReagentsById[reagent["reagent_id"]] = reagent.retro()

        return retroReagentsById

    def retrieveShallowReagentsByReactionCategoryIds(self, reactionCategoryIDs):
        """Load shallow reagent models that fall under the give reaction categories.
        Load from the database if necessary, but check in the objectCache first.
        """
        if "reagentIdsByReactionCategoryId" not in self.objectCache:
            # Not already pre-cached.  Get it from the database then, doing a full
            #   query for all data, don't even bother filtering by category.
            query = SQLQuery()
            query.addSelect("rcr.reagent_id, rcr.reaction_category_id")
            query.addFrom("reaction_category_reagent as rcr")
            reagentData = DBUtil.execute(query, connFactory=self.connFactory)

            reagentIdsByReactionCategoryId = dict()
            for row in reagentData:
                reagent_id = row[0]
                reactionCategoryID = row[1]

                if reactionCategoryID not in reagentIdsByReactionCategoryId:
                    reagentIdsByReactionCategoryId[reactionCategoryID] = set()
                    # Initialize empty set for this category

                reagentIdsByReactionCategoryId[reactionCategoryID].add(reagent_id)

            self.objectCache[
                "reagentIdsByReactionCategoryId"
            ] = reagentIdsByReactionCategoryId

        reagentIdsByReactionCategoryId = self.objectCache[
            "reagentIdsByReactionCategoryId"
        ]
        # Track which reagents belong to which categories (caching reaction_category_reagent table)

        expectedReagentIds = set()
        if reactionCategoryIDs is not None:
            for reactionCategoryID in reactionCategoryIDs:
                reactionCategoryID = int(reactionCategoryID)
                # Ensure proper type for dictionary lookup
                if reactionCategoryID in reagentIdsByReactionCategoryId:
                    expectedReagentIds.update(
                        reagentIdsByReactionCategoryId[reactionCategoryID]
                    )
        else:
            # Assume all are expected
            for reagentIds in reagentIdsByReactionCategoryId.itervalues():
                expectedReagentIds.update(reagentIds)

        shallowReagents = self.retrieveShallowReagents()

        # Extract just the expected shallow reagents from the list now
        reagentList = []
        for reagent in shallowReagents:
            if reagent["reagent_id"] in expectedReagentIds:
                reagentList.append(reagent)

        return reagentList

    def retrieveShallowReagents(self):
        """Load all available shallow reagent models.
        Load from the database if necessary, but check in the objectCache first.
        """
        if "shallowReagents" not in self.objectCache:
            # None pre-cached, load them all now (but just shallow representations for view-only)
            reagentQuery = ReagentSearchModel()
            reagentQuery.enabled = True
            reagentQuery.shallowReagents = True
            # Just for display and selection purposes here

            shallowReagents = self.loadDBInstances(reagentQuery)

            self.objectCache["shallowReagents"] = shallowReagents
        shallowReagents = self.objectCache["shallowReagents"]

        return shallowReagents

    def retrieveShallowReagentsById(self):
        """Similar to retrieveShallowReagents, except return as
        a dictionary of objects, keyed by reagentId.
        """
        if "shallowReagentsById" not in self.objectCache:
            # None pre-cached, load them all now (but just shallow representations for view-only)
            shallowReagents = self.retrieveShallowReagents()
            shallowReagentsById = modelDictFromList(shallowReagents, "reagent_id")

            self.objectCache["shallowReagentsById"] = shallowReagentsById
        shallowReagentsById = self.objectCache["shallowReagentsById"]

        return shallowReagentsById

    def returnReagentsById(self, reagentsById):
        """Should be called to followup and cleanup retrieveReagentsById.
        Any reagents that were loaded into the reagentsById dictionary
        should be stored back in the objectCache for future reuse.
        """
        objectCache = self.objectCache
        if objectCache is not None:
            if "reagentsById" not in objectCache:
                objectCache["reagentsById"] = dict()
            objectCache["reagentsById"].update(reagentsById)

    def returnRetroReagentsById(self, retroReagentsById):
        """Should be called to followup and cleanup retrieveRetroReagentsById.
        Any retro reagents that were loaded into the retroReagentsById dictionary
        should be stored back in the objectCache for future reuse.
        """
        objectCache = self.objectCache
        if objectCache is not None and "retroReagentQueuesById" in objectCache:
            retroReagentQueuesById = objectCache["retroReagentQueuesById"]
            saveObjectsToQueuesById(retroReagentsById, retroReagentQueuesById)

    def applyReagent(self, reactantList, reagent):
        """May be as simple as just running the reagent object against the reactantList,
        but include option to check for timedCache presence, in which case,
        check for previously cached results to just return instead of having to calculate
        fresh every time.
        """
        # Construct key string out of input items
        key = "%s>%s>" % (
            str.join(MOL_LIST_DELIM, molToSmiList(reactantList)),
            reagent["reagent_id"],
        )

        if self.resultCache is not None and key in self.resultCache:
            # Previously cached result, just return directly
            return self.resultCache[key]
        else:
            # No prior cached result, calculate from scratch then
            supplementalData = SupplementalDataModel()
            productList = reagent(reactantList, supplementalData)
            result = (productList, supplementalData)
            if self.resultCache is not None:
                self.resultCache[key] = result
                # Store in result cache if available for future use
            return result
