"""Common objects / base classes used to support DB interactions.
"""
from chemutils.Common.Env import SQL_PLACEHOLDER
from chemutils.Common.Util import generatePlaceholders
from chemutils.Common.Const import ALPHA, BETA, WEIGHT

from openeye.oechem import OEGraphMol, OEParseSmiles


class RowItemModel(dict):
    """Generic object class to model rows from database tables.
    Basically is just a dictionary mapping column names to values, but
    extra abstraction for additional data and functions.
    """

    """Name of the table this is supposed to represent a row of"""
    tableName = None

    """Simple boolean field indicating whether this item is newly inserted
    this run, or if it was retrieved as an existing row.
    """
    isNew = None

    """Molecule object associated with this record"""
    mol = None

    """If generated as a series of items, can record the index position in the original list here"""
    index = None
    """Can also count up the number of new items found, not just the total count"""
    newCount = None

    """Arbitrary child data if needed, perhaps to precompute some child data before actual traversal to children"""
    childData = None

    def __init__(self, initData=None, dataKeys=None):
        """Initialization constructor.  If no parameters provided, will
        just create an empty model object.  If a single dictionary
        parameter is provided, will copy over the contents into the row
        item model.  If two parameters are provided, expect these to be
        a list of data values and a list of data names / keys to be
        added in the row item model.
        """
        dict.__init__(self)
        self.update(initData, dataKeys)

    def update(self, initData=None, dataKeys=None):
        """Same like the constructor but can do at any time to update (overwrite)
        or extend the data already in the model.
        """
        if initData != None:
            if dataKeys != None:
                # Have both initData and dataKeys.  Expect these to be lists of data and names / keys
                for key, value in zip(dataKeys, initData):
                    self[key] = value
            else:
                # Only have initData, expect this to be a dictionary.  Copy over contents
                for key, value in initData.iteritems():
                    self[key] = value

    def valuesByName(self, columnNames):
        """Return the values in the dictionary as a list.  Unlike the basic
        dict.values() method, can provide a list of columnNames to only return
        the values keyed by the names provided in that list, and in that order.
        """
        values = []
        for col in columnNames:
            values.append(self[col])
        return values


def modelListFromTable(results, columnNames=None):
    """Given a table (2D array) of result items, where each row represents
    an row / item / object from the database, create an equivalent list of
    RowItemModels representing those items.

    If the columnNames are not supplied, then requires that the first row
    of the results table actually be the names
    of the columns to key the items of the RowItemModels by.
    If used the DBUtil.execute method, this is easy as you just need to set
    the "includeColumnNames" parameter.
    """
    if columnNames is None:
        columnNames = results.pop(0)

    modelList = []
    for row in results:
        modelList.append(RowItemModel(row, columnNames))
    return modelList


def modelDictFromList(modelList, columnName, listValues=False):
    """Given a list of model objects, create a dictionary keyed by
    the models' values of the columnName attribute with values
    equal to the corresponding model objects.

    If listValues is True, rather the dictionary values being the
    model objects themselves, they will instead be list objects
    with the models appended into the lists.  This would be uninteresting
    for "unique key" columnNames as it would just create dictionary
    where every value is a list of size 1.  However, this is important
    for non-unique key columnNames to make sure that model objects
    from the original list are not overwritten / lost.
    """
    modelDict = {}
    for model in modelList:
        key = model[columnName]
        if listValues:
            if not modelDict.has_key(key):
                modelDict[key] = []
                # Create a new list
            modelDict[key].append(model)
        else:
            modelDict[key] = model
    return modelDict


def columnFromModelList(modelList, columnName):
    """Given a list of model objects and the name of a
    "column" / attribute, return a list representing
    a "column" of those values, one from each model object.
    """
    column = []
    for model in modelList:
        column.append(model[columnName])
    return column


class RowItemFieldComparator:
    """Comparator object to compare pairs of RowItemModel objects,
    though it will work on any dictionary object.
    Takes a field name (column name / key) as an initialization parameter.
    Compares the total RowItemModel objects based on the each objects
    field value namd by the given key.
    """

    fieldName = None

    def __init__(self, fieldName):
        self.fieldName = fieldName

    def __call__(self, item1, item2):
        value1 = item1[self.fieldName]
        value2 = item2[self.fieldName]
        return cmp(value1, value2)


class ChemicalSearchModel:
    """Model object to build to define search criteria for chemicals.

    For similarity searches, should have OEBaseMol objects in the similarMols list.
    Furthermore, similarity search options should be set as "FloatData" on the
    molecules.  For example, to set the alpha / beta parameters for super vs. sub-structure
    searches, do the following first, or else default values will be assumed:

    from Const import ALPHA, BETA, WEIGHT;
    mol.SetFloatData(ALPHA, 0.9);
    mol.SetFloatData(BETA,  0.1);

    Rather than having multiple options for search by can_smiles, chemical_id, etc.
    just build a list of "discreteCriteria" and have the user supply the
    name of the column to compare against and a list of values to expect.
    Likewise, rather than separate options for searching by logp, num_heavy_atoms, etc.,
    just build a list of "rangeCriteria" and have the user supply the
    name of the column and the boundary values to expect.
    """

    def __init__(self):
        # Public criteria that may be accessed and set directly
        self.similarityThreshold = None
        # Minimum similarity a target must have to satisfy similarity search
        self.sourceIds = None
        # ID numbers of the sources to filter by
        self.strictSubstructure = False
        # If true and have similarMols to search by, apply them as SMARTS patterns to match strict substructures
        self.includeSourceData = False
        # Whether to include info on which sources provide this chemical.  Disable to improve performance
        self.maxResults = None
        # Maximum number of results to return
        self.resultsStart = None
        # Offset of first result to display, usually used in combination with maxResults

        # Special criteria that should added through respective class methods
        self.similarMols = []
        # Molecules to find similarity to
        self.__textCriteria = None
        # Text criteria to search by (e.g., chemical names)
        self.__textField = None
        #   Which field to do the text search by
        self.__discreteCriteria = []
        # Discrete criteria 2-ples to search by.(col,value) pairs
        self.__rangeCriteria = []
        # Range criteria 3-ples to search by.   (col,minValue,maxValue) triples

        # Internal use attributes to facilitate searches.  Public users trying to set them will achieve no effect
        self.substructureFilters = None
        # SMARTS based substructure filters when need to do strict searching
        self.similarQuery = None
        # Actual text query string to send similarity search server

        # Output attributes
        self.maxAdvancedResultsExceeded = False
        # Results indicator if similarity search results went farther than acceptable

    def addSimilarMol(self, mol, alpha=1.0, beta=1.0, weight=1.0):
        """Convenience function to specify similar molecule parameters"""
        mol.SetFloatData(ALPHA, alpha)
        mol.SetFloatData(BETA, beta)
        mol.SetFloatData(WEIGHT, weight)
        self.similarMols.append(mol)

    def addSimilarMolBySmiles(self, smiles, alpha=1.0, beta=1.0, weight=1.0):
        """Convenience function to specify molecule by SMILES"""
        mol = OEGraphMol()
        valid = OEParseSmiles(mol, smiles)
        if valid:
            self.addSimilarMol(mol, alpha, beta, weight)
            return True
        else:
            return False

    #            raise Exception \
    #            (   """Unable to parse SMILES string: %s
    #                <ul>
    #                    <li>Try consulting the <a href="http://www.daylight.com/dayhtml_tutorials/languages/smiles/index.html">Daylight SMILES Tutorial</a>
    #                    <li>If you are trying to input an alternative molecular file format (e.g., SDF) try using the <a href="BabelWeb.py">Babel Conversion Tool</a>.
    #                    <li>If you are trying to input a molecular formula, instead use the "Molecular Formula" option under "Basic Filters"
    #                    <li>If you are trying to input a chemical name or other identifier, instead use the "Text Search" field.
    #                </ul>
    #                """ % smiles
    #            );

    def setTextCriteria(self, textValue, field):
        """Search by text criteria."""
        if textValue is not None:
            textValue = textValue.strip()
            # Dump flanking whitespace
            # textValue = textValue.lower();  # Make case-insensitive
        self.__textCriteria = textValue
        self.__textField = field

    def getTextCriteria(self):
        return (self.__textCriteria, self.__textField)

    def addDiscreteCriteria(self, column, valueList):
        """Tell search to find rows where the named column
        has values in the provided value list.

        Note, for extra safety, the column name should include
        the table name as well, like "chemical.can_smiles" not just
        "can_smiles" in case the query joins to other tables with the
        same column name, an ambiguous column reference will break the query.
        """
        # Add a 2-ple of the column and the value list
        self.__discreteCriteria.append((column, valueList))

    def popDiscreteCriteria(self):
        """Remove the last discrete criteria specified, and return the
        respective tuple.  If none exist, return None.
        """
        if len(self.__discreteCriteria) > 0:
            return self.__discreteCriteria.pop()
        else:
            return None

    def iterDiscreteCriteria(self):
        """Return iterator over discrete criteria tuples"""
        return iter(self.__discreteCriteria)

    def addRangeCriteria(self, column, minValue, maxValue):
        """Tell search to find rows where the named column
        has values >= the minValue and <= the maxValue.

        Note, for extra safety, the column name should include
        the table name as well, like "chemical.molecular_weight" not just
        "molecular_weight".
        """
        # Add a 3-ple of the search criteria
        self.__rangeCriteria.append((column, minValue, maxValue))

    def iterRangeCriteria(self):
        """Return iterator over range criteria tuples"""
        return iter(self.__rangeCriteria)


class SQLQuery:
    """Helper class to dynamically build a SQL statement.
    Unlike just appending to a string, allows for addition
    of select, from, where, etc. clauses in any order,
    then generate them in normal order at the end.
    """

    def __init__(self):
        self.delete = False
        # If set, will ignore the select list and make a delete query instead
        self.select = []
        self.from_ = []
        self.where = []
        self.groupBy = []
        self.having = []
        self.orderBy = []
        self.limit = -1
        self.offset = -1
        self.params = []

    def addSelect(self, aSelect):
        self.select.append(aSelect)

    def addFrom(self, aFrom):
        self.from_.append(aFrom)

    def addWhere(self, aWhere):
        self.where.append(aWhere)

    def addWhereEqual(self, field, param):
        self.where.append("%s = %s" % (field, SQL_PLACEHOLDER))
        self.params.append(param)

    def addWhereLike(self, field, param):
        self.where.append("%s like %s" % (field, SQL_PLACEHOLDER))
        self.params.append(param)

    def addWhereOp(self, field, op, param):
        self.where.append("%s %s %s" % (field, op, SQL_PLACEHOLDER))
        self.params.append(param)

    def addWhereIn(self, field, paramList):
        placeholders = generatePlaceholders(len(paramList))
        self.where.append("%s in (%s)" % (field, placeholders))
        self.params.extend(paramList)

    def addWhereNotIn(self, field, paramList):
        placeholders = generatePlaceholders(len(paramList))
        self.where.append("%s not in (%s)" % (field, placeholders))
        self.params.extend(paramList)

    def addGroupBy(self, aGroupBy):
        self.groupBy.append(aGroupBy)

    def addHaving(self, aHaving):
        self.having.append(aHaving)

    def addOrderBy(self, aOrderBy, dir=None):
        if dir != None:
            aOrderBy = "%s %s" % (aOrderBy, dir)
        self.orderBy.append(aOrderBy)

    def setLimit(self, limit):
        self.limit = limit

    def setOffset(self, offset):
        self.offset = offset

    def addParam(self, aParam):
        self.params.append(aParam)

    def getParams(self):
        return self.params

    def __str__(self):
        query = []  # Build list of query components, then join into a string at the end

        if self.delete:
            query.append("delete")
        else:
            query.append("select")
            for item in self.select:
                query.append(item)
                query.append(",")
            query.pop()  # Remove the last extraneous comma (",")

        query.append("from")
        for item in self.from_:
            query.append(item)
            query.append(",")
        query.pop()

        if len(self.where) > 0:
            query.append("where")
            for item in self.where:
                query.append(item)
                query.append("and")
            query.pop()

        if len(self.groupBy) > 0:
            query.append("group by")
            for item in self.groupBy:
                query.append(item)
                query.append(",")
            query.pop()

        if len(self.having) > 0:
            query.append("having")
            for item in self.having:
                query.append(item)
                query.append("and")
            query.pop()

        if len(self.orderBy) > 0:
            query.append("order by")
            for item in self.orderBy:
                query.append(item)
                query.append(",")
            query.pop()

        if self.limit > -1:
            query.append("limit")
            query.append(str(self.limit))
        if self.offset > -1:
            query.append("offset")
            query.append(str(self.offset))

        query = str.join(
            " ", query
        )  # Join list items into a single string, space-separated
        return query
