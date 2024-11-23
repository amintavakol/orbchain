"""Various constants for use by the application modules"""

import logging
from chemutils.Common import Env
#from sets import Set, ImmutableSet;

"""See oemolistream documentation for more info on this pattern"""
STD_FILE = "-"
STD_MOL_EXT = {}
try:
    from openeye.oechem import OEFormat_SMI,OEFormat_MDL,OEFormat_PDB,OEFormat_MOL2,OEFormat_ISM,OEFormat_MOL2H,OEFormat_SDF,OEFormat_CAN,OEFormat_MF,OEFormat_XYZ,OEFormat_FASTA,OEFormat_FASTA,OEFormat_MOPAC,OEFormat_OEB
    STD_MOL_EXT[OEFormat_SMI]   = '.smi'    # 1
    STD_MOL_EXT[OEFormat_MDL]   = '.mdl'    # 2
    STD_MOL_EXT[OEFormat_PDB]   = '.pdb'    # 3
    STD_MOL_EXT[OEFormat_MOL2]  = '.mol2'   # 4
    #STD_MOL_EXT[OEFormat_BIN]   = '.bin'    # 5 // no longer in OEChem 
    STD_MOL_EXT[OEFormat_ISM]   = '.ism'    # 7
    STD_MOL_EXT[OEFormat_MOL2H] = '.mol2h'  # 8
    STD_MOL_EXT[OEFormat_SDF]   = '.sdf'    # 9
    STD_MOL_EXT[OEFormat_CAN]   = '.can'    # 10
    STD_MOL_EXT[OEFormat_MF]    = '.mf'     # 11
    STD_MOL_EXT[OEFormat_XYZ]   = '.xyz'    # 12
    STD_MOL_EXT[OEFormat_FASTA] = '.fasta'  # 13
    STD_MOL_EXT[OEFormat_MOPAC] = '.mopac'  # 14
    STD_MOL_EXT[OEFormat_OEB]   = '.oeb'    # 15
except Exception:
    # OpenEye dependency is probably not installed, will just have to skip / ignore this function then
    import warnings
    warnings.warn("OEChem dependency does not appear to be installed.  Attempting to continue anyway.") #, ImportWarning)
    
"""File extension to identify GZipped files"""
GZIP_EXT = ".gz";

"""SMILES separator to indicate distinct molecules"""
SMILES_MOL_DELIM = "."

"""Reaction SMILES component separator (reactant>reagent>product)"""
REACTION_COMPONENT_DELIM = ">";

"""Image tag construction formats"""
FORMAT_OGHAM = 0;
FORMAT_DAYLIGHT = 1;
FORMAT_JME = 2;
FORMAT_OGHAM_ENCRYPT = 3;
FORMAT_CHEMAXON_MARVIN = 4;
FORMAT_CHEMAXON_MARVIN_APPLET = 5;

"""Construction size parameters for generating fingerprints"""
FP_SIZE = Env.FP_SIZE;
FP_MAX = Env.FP_MAX;
FP_MIN = Env.FP_MIN;
FP_COLUMN = Env.FP_COLUMN;

"""Characters marking the search server results"""
FP_RESULTS_START = "{";
FP_RESULTS_END = "}";
FP_RESULTS_ERROR = "***";

TEXT_RESULTS_START = "searchName";
TEXT_RESULTS_COL_DELIM = "|";
TEXT_RESULTS_LINE_DELIM = "#!#";

"""Character to add to end of string to inform Lucene to do fuzzy name search"""
FUZZY_NAME_SUFFIX = "~";

"""Default max result set size.  Must have some value for similarity searches.  
Also have a limit on how big the results set size can be.  Don't let users abuse
or otherwise swamp the system with queries for thousands of records at a time."""
DEFAULT_MAX_RESULTS = 10;
MAX_RESULTS_LIMIT   = 100;

"""Max advanced (similarity or name) search results.  When doing advanced searches, 
don't look anymore past the top X matches, where X is this number.  
Important for limiting the depth of the search
when doing integrated advanced + basic searches."""
ADVANCED_RESULTS_LIMIT = 2500;

"""Tag to mark comments in processed files"""
COMMENT_TAG = "#"

"""Delimiter of SQL commands in a DB script file"""
SQL_DELIM = ";"

"""Tag to indicate a multi-part or multi-line parameter counts as a single token"""
TOKEN_END = '"'

"""Null string used to represent DB null value"""
NULL_STRING = str(None);

"""Null tag used to represent DB null value"""
NULL_TAG = "<NULL>";

"""Constants to key similarity search criteria"""
ALPHA = "alpha";
BETA  = "beta";
WEIGHT= "weight";

"""Suffix to add after a source abbreviation to identify respective SD tags"""
EXTERNAL_ID_SUFFIX = ".external_id";

"""Columns to be EXCLUDED from SD annotations for isomer downloads"""
EXCLUDED_SD_TAGS = set();
EXCLUDED_SD_TAGS.add("sdf");
EXCLUDED_SD_TAGS.add("mixturecomponent_id");
EXCLUDED_SD_TAGS.add("dev_flag");
EXCLUDED_SD_TAGS.add("count_isomer3d");
#EXCLUDED_SD_TAGS.add("fingerprint_1");
#EXCLUDED_SD_TAGS.add("fingerprint512");
#EXCLUDED_SD_TAGS.add("oepassed_filter");

"""Default suffix for database table id fields, particularly the primary key, but also
foreign key columns that reference other tables."""
DEFAULT_ID_COL_SUFFIX = "_id";

"""Default name of the code / name column of a lookup code table"""
CODE_COL = "code";

"""Default name of the code table column to order by"""
CODE_ORDER = "position, code";

"""Wildcard string used in SQL queries"""
SQL_WILDCARD = "%";

"""Single common instance of an empty set to save on instantiation costs when have many of them"""
EMPTY_SET = frozenset();

"""Not functionally important","just indicates an estimate of the number
of lines in input files to hint a proper scale for the progress indicators.
If you do not wish to see those dot indicators","set this value to 0."""
EST_INPUT = Env.EST_INPUT

"""Updates to process before reporting progress"""
PROG_BIG = 10000
PROG_SMALL = 500

"""Application name","for example to identify a common logger object"""
APPLICATION_NAME = "CHEM.Common.app"

"""Default level for application logging.  Modify these for different scenarios.  See Python logging package documentation for more information"""
LOGGER_LEVEL = Env.LOGGER_LEVEL

"""Default format of logger output"""
LOGGER_FORMAT = "[%(asctime)s %(levelname)s] %(message)s"
