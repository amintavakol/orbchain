import sys, os
from optparse import OptionParser
#from sets import Set;

from chemutils.Common import Const, Util, Env
from chemutils.Common.Util import stdOpen, virtual_oemolistream, virtual_oemolostream, ProgressDots;
from chemutils.Common.ResultsFormatter import TextResultsFormatter;
from chemutils.Common import DBUtil;
from chemutils.Common.Env import chemutils_DB_PARAM;

from chemutils.Common.Const import PROG_BIG, PROG_SMALL;
from chemutils.Common.Env import SQL_PLACEHOLDER
from chemutils.Common.Model import RowItemModel, SQLQuery, ChemicalSearchModel;
from chemutils.Common.Model import modelListFromTable, modelDictFromList, columnFromModelList
from chemutils.Common.ChemicalSearch import ChemicalSearch;

from openeye.oechem import oemolistream, oemolostream, OEReadMolecule, OEWriteMolecule, OEFormat_SDF; 
from openeye.oechem import OEGraphMol, OEParseSmiles, OESetSDData, OEDeleteSDData, OEGetSDDataPairs;

def main(argv):
    """Main method, callable from command line"""
    usageStr =  "usage: %prog [options] <outputFile> <searchString>\n"+\
                "   <outputFile>    File to output to.  Specify \".sdf\" to send to stdout.\n"+\
                "   <searchString>  SMILES string to search for chemical details by.\n"+\
                "                   Omitting this param when doing \n"  +\
                "                   an SDF download indicates to download the ENTIRE database.\n"  +\
                "                   Multiple search values can be specified, separated by spaces.\n" 
                
    parser = OptionParser(usage=usageStr)
    parser.add_option("-a", "--annotation", dest="annotation",  action="store_true",help="If set, returns the annotation table for this chemical.")
    parser.add_option("-s", "--sdf",        dest="sdf",         action="store_true",help="If set, returns the SDF 3D file contents for this chemical.")
    (options, args) = parser.parse_args(argv[1:])

    if len(args) >= 1:
        instance = ChemicalDetail()

        if options.annotation and len(args) >= 2:
            outputFile = stdOpen(args[0],"w",sys.stdout)
            results = instance.annotationDictsBySmiles( args[1] )
            formatter = TextResultsFormatter(outputFile)
            formatter.formatResultDicts(results,ChemicalDetail.ANNOTATION_COLS)
        else:
            if len(args) < 2:   smilesSearch = None;
            else:               smilesSearch = args[1:];
            instance.retrieveIsomersBySmiles( smilesSearch, oemolostream(args[0]) );
    else:
        parser.print_help()
        sys.exit(-1)

class ChemicalDetail:
    """Various utility methods for retrieving details for some
    chemical from the database, usually keyed by the canonical 
    SMILES string.  See individual methods for details.
    """

    """Columns to expect in annotation results"""
    ANNOTATION_COLS = ["iso_smiles","abbrev","website","external_id","name","svalue","fvalue"]

    def __init__(self):
        """Constructor"""
        pass

    def chemicalDetails(self, searchColumn, searchValue):
        """Find an item in the chemical table based on a specified searchColumn
        and value.  Will only return the first chemical result so typically only
        want to search by a unique identifier like "iso_smiles" or "chemical_id."
        
        Returns a RowItemModel dictionary representing all of the data found 
        in that row.  In addition, count up the number of isomer3d records 
        under the chemical and store that in an entry named count_isomer3d.
        """
        # Usually a bad idea to do select *, better to explicitly specify columns,
        #   but since storing to a RowItemModel object that can accept any number
        #   of attributes in any order (dictionary based), it should be good
        query = "select * from chemical where %s = %s" % (searchColumn,SQL_PLACEHOLDER);
        conn = DBUtil.connection( chemutils_DB_PARAM );
        results = DBUtil.execute(query,(searchValue,),includeColumnNames=True,conn=conn)        

        chemicalModel = None;
        # Since set to "includeColumnNames," the first row should be the column names
        #   to key the data by, and the second row should be the actual data.
        if len(results) > 1:
            colNames = results[0];
            data = results[1];

            chemicalModel = RowItemModel( data, colNames );
        
            # Extra query to get isomer count.  Do it separately for "outer join" effect
            query = "select sdf, isomer from isomer3d where chemical_id = %s" % SQL_PLACEHOLDER
            results = DBUtil.execute(query,(chemicalModel["chemical_id"],))

            chemicalModel["count_isomer3d"] = len(results);
            for row in results: # Sample the 1st predicted 3D coordinates in SDF format
                if row[1] == 1:
                    chemicalModel["isomer3d.sdf"] = row[0];
        
        return chemicalModel

    def annotationDictsBySmiles(self, smiles, limit=None):
        """Find an item in the chemical table based on the given SMILES string.
        Then trace through the whole DB schema
        (chemical - mixturecomponent - chemicalmix - source2chemicalmix - source - annotation)
        to find all of the chemicalmix and sources this chemical belongs to.
        Furthermore, find all of the names and values of annotations applied to such
        chemical mixes, ordered by the source that specified the annotation.

        Do an outer-join.  That is, if a chemical exists under some
        source, but no annotations exist, still at least report a row for that
        source with blank annotation values.
        
        The results should be a list of RowItemModel dictionaries, at least containing
        the following "columns" in each model object.
        
            chemicalmix.iso_smiles
            source.abbrev
            source2chemicalmix.external_id
            annotation.name
            annotation.svalue
            annotation.fvalue
        """
        # The "left join ... on" syntax does an outer-join.
        #   This is DB specific (PostgreSQL) syntax and would
        #   have to be changed if the app DB ever does.
        #
        # Actually, break this into multiple queries, cause
        #   this isn't working optimally.  Causes sequence scans
        #   when all should be doable with index scans.
        #
        """
        annotationQuery = \
            ""
            select 
               cm.iso_smiles, 
               s.abbrev, 
               s.website,
               s2c.external_id,
               a.name, 
               a.svalue, 
               a.fvalue
            from 
               chemical as c, 
               mixturecomponent as mc,
               chemicalmix as cm,
               source as s,

               source2chemicalmix as s2c left join annotation as a
               on s2c.source2chemicalmix_id = a.source2chemicalmix_id 
            where 
               c.chemical_id = mc.chemical_id  and 
               mc.chemicalmix_id = cm.chemicalmix_id and 
               s2c.chemicalmix_id = cm.chemicalmix_id and 
               s2c.source_id = s.source_id and 
               c.iso_smiles = %s
            order by
               iso_smiles,
               abbrev,
               name
            "" % SQL_PLACEHOLDER
        """             
        
        chemicalmixQuery = SQLQuery();
        chemicalmixQuery.addSelect("cm.chemicalmix_id");
        chemicalmixQuery.addFrom("chemical as c");
        chemicalmixQuery.addFrom("mixturecomponent as mc");
        chemicalmixQuery.addFrom("chemicalmix as cm");
        chemicalmixQuery.addWhere("c.chemical_id = mc.chemical_id");
        chemicalmixQuery.addWhere("mc.chemicalmix_id = cm.chemicalmix_id");
        chemicalmixQuery.addWhereEqual("c.iso_smiles",smiles);
        if isinstance( limit, int ):
            chemicalmixQuery.setLimit(limit);

        chemicalmix_ids = DBUtil.execute(chemicalmixQuery);
        for i in xrange(len(chemicalmix_ids)):
            chemicalmix_ids[i] = chemicalmix_ids[i][0]; # Convert lists of size 1, to single items
        if len(chemicalmix_ids) < 1:
            chemicalmix_ids.append(-1); # Append dummy value so future queries don't break

        s2cQuery = SQLQuery();
        s2cQuery.addSelect("s2c.source2chemicalmix_id");
        s2cQuery.addSelect("s2c.external_id");
        s2cQuery.addSelect("cm.iso_smiles");
        s2cQuery.addSelect("s.abbrev");
        s2cQuery.addSelect("s.website");
        s2cQuery.addFrom("chemicalmix as cm")
        s2cQuery.addFrom("source2chemicalmix as s2c");
        s2cQuery.addFrom("source as s")
        s2cQuery.addWhere("cm.chemicalmix_id = s2c.chemicalmix_id");
        s2cQuery.addWhere("s.source_id = s2c.source_id");
        s2cQuery.addWhereIn("cm.chemicalmix_id",chemicalmix_ids);

        s2cTable = DBUtil.execute(s2cQuery,includeColumnNames=True);
        s2cModels= modelListFromTable(s2cTable);
        s2cById  = modelDictFromList(s2cModels,"source2chemicalmix_id");
        
        annotationQuery = SQLQuery();
        annotationQuery.addSelect("s2c.source2chemicalmix_id");
        annotationQuery.addSelect("a.name");
        annotationQuery.addSelect("a.svalue");
        annotationQuery.addSelect("a.fvalue");
        annotationQuery.addFrom("source2chemicalmix as s2c left join annotation as a on s2c.source2chemicalmix_id = a.source2chemicalmix_id");
        annotationQuery.addWhereIn("s2c.chemicalmix_id",chemicalmix_ids);
        annotationQuery.addOrderBy("s2c.chemicalmix_id");
        annotationQuery.addOrderBy("s2c.source2chemicalmix_id");
        annotationQuery.addOrderBy("a.name");
        
        annotationTable = DBUtil.execute(annotationQuery,includeColumnNames=True);
        annotationModels = modelListFromTable(annotationTable);
        # "Join" with info from s2cQuery
        for annotation in annotationModels:
            s2cId = annotation["source2chemicalmix_id"];
            s2cModel = s2cById[s2cId];
            annotation.update(s2cModel);
        
        return annotationModels;

    def separateAnnotationDictsByMetadata( self, annotationDicts, metadataColumn ):
        """Given a list of dictionaries (RowItemModels) representing Annotation rows
        from the annotationDictsBySmiles method, see if any of the rows match
        items in the annotation_metadata table whose specified metadataColumn is True.
        
        If so, remove those items from the original list, and return them in a new,
        separate list of annotation dictionary models.
        """
        # Find the metadata rows that satisfy the property
        query = \
            """select s.abbrev, am.name 
            from annotation_metadata as am, source as s 
            where s.source_id = am.source_id 
            and %s is true
            """ % metadataColumn;
        metadataResults = DBUtil.execute(query,includeColumnNames=True);
        metadataSet = set();    # Convert into a Set of tuples for rapid lookup
        for row in metadataResults:
            metadataSet.add( tuple(row) );
        
        filteredDicts = [];
        for i in xrange(len(annotationDicts)-1,-1,-1):  # Iterate in reverse order since we may be deleting items
            annotationModel = annotationDicts[i];
            searchTuple = (annotationModel["abbrev"], annotationModel["name"]);
            if searchTuple in metadataSet:
                annotationDicts.pop(i);
                filteredDicts.append(annotationModel);
        return filteredDicts;
        
    def annotatedChemicalMols(self, smilesList, idInfoOnly=False, streamOutput=True):
        """Given one or a list of chemical iso_smiles values, return an iterator over OEMolBase
        objects representing the chemicals.  Extra value of running the
        annotationDictsBySmiles method for each one and recording these values as
        SD Tags on the molecule object.  If the molecule is subsequently output in SDF format,
        the caller can view these in file output.
        
        Note that this returns an iterator over the same molecule object, it just gets cleared
        and rewritten at each iteration.  If for some reason you want a separate copy for each,
        you'll need to do that yourself by doing something like "copy = OEGraphMol(mol)"
        
        Since a chemical can belong to multiple chemical mixes, this is a little
        sloppy situation.  May only want to count when the chemical is the
        "first largest component" of the mix or something.  In the meantime, only record 
        sources and external IDs this is under.
        
        idInfoOnly: If set to True, just load the found chemical_ids and molecule titles.
            Don't spend time searching for 3D structures or source information.
        """
        # If only provided one value, convert into a list of size 1
        if isinstance(smilesList,str):
            smilesList = [smilesList];
        
        # Forget all the annotations for now.  Just record the sources and external_id value.
        # Will enter a lot of redundant data with this loop structure, but it's okay, will just
        #   keep overwriting the same source, external_id value
        """
        mol = OEGraphMol();
        for smiles in smilesList:
            OEParseSmiles(mol,smiles);
            chemicalModel   = self.chemicalDetails("iso_smiles",smiles);
            mol.SetTitle( str(chemicalModel["chemical_id"]) );
            annotationDicts = self.annotationDictsBySmiles(smiles);
            for annotationModel in annotationDicts:
                sourceAbbrevTag = annotationModel["abbrev"]+Const.EXTERNAL_ID_SUFFIX;
                OESetSDData( mol, sourceAbbrevTag, annotationModel["external_id"] );
            yield mol;
            mol.Clear();        
        """    

        searchModel = ChemicalSearchModel();
        searchModel.addDiscreteCriteria("iso_smiles",smilesList);
        searchModel.includeSourceData = not idInfoOnly;
        aChemicalSearch = ChemicalSearch();
        chemicalModels = aChemicalSearch.findChemicals( searchModel );

        if streamOutput:
            mol = OEGraphMol();
            for chemicalModel in chemicalModels:
                if not idInfoOnly:
                    # Try to get more chemical detail information.  A sample 3D SDF file in particular
                    chemicalDetailModel = self.chemicalDetails("iso_smiles",chemicalModel["iso_smiles"]);
                    if "isomer3d.sdf" in chemicalDetailModel:
                        # An SDF formatted isomer with 3D coordinates is available.
                        self.__parseSDFString( mol, chemicalDetailModel["isomer3d.sdf"] );
                    else:
                        # No 3D structures available.  At least get the connectivity information from the SMILES string
                        OEParseSmiles( mol, chemicalModel["iso_smiles"] );
                else:
                    # Not interested in 3D structure.  Just get the connectivity information from the SMILES string
                    OEParseSmiles( mol, chemicalModel["iso_smiles"] );
                    
                mol.SetTitle( str(chemicalModel["chemical_id"]) );
                
                if "source.abbrev" in chemicalModel and "source2chemicalmix.external_id" in chemicalModel:
                    for sourceAbbrev, external_id in zip(chemicalModel["source.abbrev"],chemicalModel["source2chemicalmix.external_id"]):
                        OESetSDData( mol, sourceAbbrev+Const.EXTERNAL_ID_SUFFIX, external_id );
    
                yield mol;
                
                # Reset state
                mol.Clear();    # Doesn't clear SD Tags
                for dataPair in OEGetSDDataPairs(mol):
                    OEDeleteSDData(mol, dataPair.GetTag());
        else:
            yield chemicalModels        
   
    def retrieveIsomersBySmiles(self,smiles,targetOEOS=None,streamOutput=True):
        """Find an item in the chemical table based on the given SMILES string.
        Then find respective records in the isomer3d table and any SDF 3D molecular
        data files in particular.  If so, then return an oemolistream of OEMolBase
        objects to represent the contents of these SD Files.
        
        If the smiles provided is actually a list, will return all SDFs under 
        those chemical SMILES.

        Furthermore, include SD annotation tags for:
           - All primary chemical annotations / columns (except mixturecomponent_id)
           - Isomer #
           - Isomeric SMILES
        
        If a targetOEOS (oemolostream) is provided, will write the output 
        directly there instead of managing it in memory.  No OEIS will be returned then.
        This is in case the query size is huge (i.e., the whole database) and
        should be processed as a stream of data.
        """
        directOutput = True;
        if targetOEOS is None: 
            directOutput = False;
            targetOEOS = virtual_oemolostream(OEFormat_SDF);

        # Convert to a list as necessary
        if isinstance(smiles,str): 
            smiles = [smiles];
        
        query = SQLQuery();
        query.addSelect("i.sdf, i.isomer, i.iso_smiles, c.*");
        query.addFrom("chemical as c, isomer3d as i");
        query.addWhere("c.chemical_id = i.chemical_id");
        if smiles is not None:
            query.addOrderBy("c.iso_smiles, i.isomer");
            query.addWhereIn("c.iso_smiles",smiles); 

        conn = DBUtil.connection();
        cursor = conn.cursor();
        cursor.execute(str(query),tuple(query.getParams()));
#        cursor.execute("""SELECT i.sdf, i.isomer, i.iso_smiles, c.*
#                          FROM chemical AS c, isomer3d AS i
#                          WHERE c.chemical_id = i.chemical_id
#                          AND c.chemical_id IN (4115654)
#                          ORDER BY c.iso_smiles, i.isomer""");
        colNames = DBUtil.columnNamesFromCursor(cursor);

        progress = ProgressDots(PROG_BIG,PROG_SMALL,"Isomers");
        mol = OEGraphMol();
        row = cursor.fetchone();
        isomerModels = []
        while row is not None:
            isomerModel = RowItemModel(row,colNames);
            mol.Clear();
            self.__parseSDFString( mol, isomerModel["sdf"] );
            self.__applySDTags(mol,isomerModel);
            OEWriteMolecule(targetOEOS,mol);

            row = cursor.fetchone();
            progress.Update();
            
            if isomerModel["sdf"] != None:
                isomerModel["sdf"] = "yes";
            else:
                isomerModel["sdf"] = "no";
            isomerModels.append(isomerModel)
        targetOEOS.close();

        if not directOutput:
            if streamOutput:
                return virtual_oemolistream(targetOEOS.GetString(),OEFormat_SDF);
            else:
                return isomerModels;
        return None;
    
    def __parseSDFString(mol,sdfString):
        """Given an SDF string, parse it into the given OEMolBase molecule object
        as if the string represented an SDF file.  Will only get the first molecule
        in the SDF if it contains multiple.
        """
        sourceOEIS = virtual_oemolistream( sdfString, OEFormat_SDF );

        success = OEReadMolecule(sourceOEIS,mol);
        if success < 0: 
            # HACK!  Seems like some SDF files in DB don't have a label / title
            #   line to start with.  If so, oemolistream can't read it right.
            #   If this happens, try to compensate by just droping in a blank title line
            sourceOEIS = virtual_oemolistream("\n"+sdfString, OEFormat_SDF ); # Restart, but try adding blank title line
            OEReadMolecule(sourceOEIS,mol);
    __parseSDFString = staticmethod(__parseSDFString);

    def __applySDTags(self,isomerMol,isomerModel):
        """Given an isomerMol OEMolBase object, add SD annotation tags on it
        based on the contents of the isomerModel RowItemModel (dictionary)
        representing column values from the database.
        """
        # Just add every column known in the isomerModel except for a specific ones
        #   we want to exclude.  This includes the SDF (don't want the SDF text to
        #   be an annotation of the SDF molecule itself) and the mixturecomponent_id
        #   which is higher level than should be cared about.
        for col, value in isomerModel.iteritems():
            if not col in Const.EXCLUDED_SD_TAGS:
                OESetSDData(isomerMol,col,str(value));
        # Modify the isomer title to include parent chemical ID information
        isomerMol.SetTitle( str(isomerModel["chemical_id"])+isomerMol.GetTitle() );

if __name__ == "__main__":
    main(sys.argv)
