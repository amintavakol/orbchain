"""Various constants for use by the application modules,
but these can / should be changed depending on the 
platform / environment where they are installed.
"""

import logging
import socket

"""Default level for application logging.  Modify these for different scenarios.  
See Python logging package documentation for more information"""
# LOGGER_LEVEL = logging.DEBUG
# LOGGER_LEVEL = logging.INFO
# LOGGER_LEVEL = logging.WARNING
LOGGER_LEVEL = logging.ERROR
# LOGGER_LEVEL = logging.CRITICAL

"""Not functionally important, just indicates an estimate of the number
of lines in input files to hint a proper scale for the progress indicators.
If you do not wish to see those dot indicators, set this value to 0.

Several usages divide this value by 200 to indicate the desired number
of total dots ~200.  Thus, do not set the value < 200 or else there
will effectively be no results, unless that is your desired behavior (set to 0)
"""
EST_INPUT = 10000

"""Parameter placeholder for SQL queries.  Depends on the DB-API module
used.  Can tell which by the module's paramstyle attribute.
Read the DB-API 2.0 specs for more info. (http://www.python.org/peps/pep-0249.html)"""
# SQL_PLACEHOLDER = "?"   # "qmark"
SQL_PLACEHOLDER = "%s"  # "format" and "pyFormat"

"""Strings to use for boolean parameters."""
BOOLEAN_STR = dict()
# pyPgSQL
BOOLEAN_STR[True] = str(True)
BOOLEAN_STR[False] = str(False)
# MS Access
# BOOLEAN_STR[True] = str(-1);
# BOOLEAN_STR[False]= str(0);

"""Construction size parameters for generating fingerprints"""
FP_SIZE = 1024
FP_MAX = 8
FP_MIN = 0
FP_COLUMN = "fingerprint_1"

"""Local host name and IP address which can be translated by DNS.
On ICS computers, has annoying habit of just giving back the prefix (e.g., "contact6",
instead of the entire name "contact6.ics.uci.edu", so may have to do some extra work here.
"""
try:
    HOSTNAME = socket.gethostname()
    IP_ADDR = socket.gethostbyname(HOSTNAME)
except:
    HOSTNAME = "localhost"
    IP_ADDR = "127.0.0.1"

"""Hostname for Smi2Depict services"""
# DEPICT_HOST = "c1ccccc1-34.ics.uci.edu";
DEPICT_HOST = "new-c1ccccc1-34.ics.uci.edu:8081"
# JAVA_DEPICT_HOST = "cdb.ics.uci.edu:8080";
# JAVA_DEPICT_HOST = "c1ccccc1-34.ics.uci.edu:8081";
JAVA_DEPICT_HOST = "new-c1ccccc1-34.ics.uci.edu:8081"
# JAVA_DEPICT_HOST = "localhost:8080";

"""Directory of Marvin applets relative to cgibin dir.
Use to specifiy and distinguish specific version installations,
otherwise client browsers may get confused on which version to use
when doing applet upgrades.
"""
MARVIN_APPLETS_REL_DIR = "../resource/marvin-5.4.1.2"

"""Connection parameters for fingerprint server"""
FP_SERVER = "c1ccccc1-34.ics.uci.edu"
FP_PORT = 5555
# FP_SERVER = "finger.ics.uci.edu"
# FP_PORT = 6002
FP_WEB_SERVICE_BASE = "http://c1ccccc1-34.ics.uci.edu/CHEM/Web/cgibin/admin/FingerServerWeb.py?FingerServerWeb=1&fpHost=%s&fpPort=%s&outputOnly=1&query="
FP_WEB_SERVICE = FP_WEB_SERVICE_BASE % (FP_SERVER, FP_PORT)
FP_USE_WEB_SERVICE = False

"""Connection parameters for name search server"""
TEXT_WEB_SERVICE = "http://c1ccccc1-34.ics.uci.edu:8081/chemdb/results.jsp?query=%s&field=%s&maxresults=%s"
# TEXT_WEB_SERVICE = "http://%s.ics.uci.edu:8081/chemdb/results.jsp" % HOSTNAME + "?query=%s&field=%s&maxresults=%s";
# TEXT_WEB_SERVICE = "http://cdb.ics.uci.edu:8081/chemdb/results.jsp?query=%s&field=%s&maxresults=%s";
# TEXT_WEB_SERVICE = "http://mine1.ics.uci.edu:8081/chemdb/results.jsp?query=%s&field=%s&maxresults=%s";

"""Parameters needed to open a connection to the database.  
Dependent upon particular connection interface and database implementation
"""
DB_PARAM = {}
DB_PARAM["HOST"] = "easyv1.ics.uci.edu"
DB_PARAM["DSN"] = "cdbv2"
DB_PARAM["UID"] = "cdbweb"
DB_PARAM["PWD"] = "FxnGrp#2004"

"""Connection parameters for reaction tutorial data.
"""
TUTORIAL_DB_PARAM = dict()
TUTORIAL_DB_PARAM["HOST"] = "easyv1.ics.uci.edu"
TUTORIAL_DB_PARAM["DSN"] = "reaction_tutorial"
TUTORIAL_DB_PARAM["UID"] = "cdbweb"
TUTORIAL_DB_PARAM["PWD"] = "FxnGrp#2004"
TUTORIAL_DB_PARAM["PORT"] = 6002

"""Connection parameters for reaction applications, relating to
'ReactDB,' classification, retro-synthesis, etc. 
"""
REACTION_DB_PARAM = dict()
REACTION_DB_PARAM["HOST"] = "easyv1.ics.uci.edu"
REACTION_DB_PARAM["DSN"] = "reaction_dev"
REACTION_DB_PARAM["UID"] = "cdbweb"
REACTION_DB_PARAM["PWD"] = "FxnGrp#2004"
REACTION_DB_PARAM["PORT"] = 5999

"""Connection parameters for ChemDB searches
"""
CHEM_DB_PARAM = dict()
CHEM_DB_PARAM["HOST"] = "breezyv5.ics.uci.edu"
CHEM_DB_PARAM["DSN"] = "cdbv2"
CHEM_DB_PARAM["UID"] = "cdbweb"
CHEM_DB_PARAM["PWD"] = "FxnGrp#2004"

"""Connection parameters for OrbPair Scoring database work
"""
ORB_DB_PARAM = dict()
ORB_DB_PARAM["HOST"] = "db-igb.ics.uci.edu"
ORB_DB_PARAM["DSN"] = "orb_pair_score_dev"
ORB_DB_PARAM["UID"] = "devtest"
ORB_DB_PARAM["PWD"] = None
ORB_DB_PARAM["PORT"] = 5999

"""Failover backup for production DB.  Should essentially be a clone, but data inserts
may be running here, whereas production should remain untouched."""
FAILOVER_DB_PARAM = {}
FAILOVER_DB_PARAM["HOST"] = "chemdb.ics.uci.edu"
FAILOVER_DB_PARAM["DSN"] = "cdbtest"
FAILOVER_DB_PARAM["UID"] = "devtest"
FAILOVER_DB_PARAM[
    "PWD"
] = None  # Comment out this line to force user to enter password at runtime


"""Connection parameters for development database.  Any inserts, updates, etc. should
generally only apply to the development database.  The production database should essentially
be read-only making it an easy process to resync by just dropping the production database
and reimporting from the dev server.
"""
DEV_DB_PARAM = {}
DEV_DB_PARAM["HOST"] = "chemdb.ics.uci.edu"
DEV_DB_PARAM["DSN"] = "cdbtest"
DEV_DB_PARAM["UID"] = "devtest"
DEV_DB_PARAM["PWD"] = None
# Comment out this line to force user to enter password at runtime

"""Connection parameters for unit tests.  Should use these DB and fingerprint server connection 
parameters to avoid messing the "real" stuff up"""
TEST_DB_PARAM = {}
TEST_DB_PARAM["HOST"] = "breezyv5.ics.uci.edu"
TEST_DB_PARAM["DSN"] = "cdbv2"
TEST_DB_PARAM["UID"] = "cdbweb"
TEST_DB_PARAM["PWD"] = "FxnGrp#2004"

""" TEST Connection parameters for TEST OrbPair Scoring database work
"""
TEST_ORB_DB_PARAM = dict()
TEST_ORB_DB_PARAM["HOST"] = "db-igb.ics.uci.edu"
TEST_ORB_DB_PARAM["DSN"] = "test_orb_pair_score_dev"
TEST_ORB_DB_PARAM["UID"] = "devtest"
TEST_ORB_DB_PARAM["PWD"] = None
TEST_ORB_DB_PARAM["PORT"] = 5999


TEST_FP_SERVER = "finger.ics.uci.edu"
TEST_FP_PORT = 6000

TEST_TEXT_WEB_SERVICE = "No Text Search Service exists for Test Data!"


def formatDBConnectString(dbParamDict):
    connStr = ""
    for key, value in dbParamDict.iteritems():
        connStr += "%s=%s;" % (key, value)
    return connStr
