from urllib import urlopen, urlencode
from cgi import parse_qsl

## Some chemutils specific entries:
from chemutils.Common.Env import JAVA_DEPICT_HOST

from chemutils.Common.Util import molBySmiles
from chemutils.CombiCDB.ArrowCodesToOrbitals import (
    orbitalPairFromArrowStr,
)

CHEMAXON_MARVIN_URL = "http://%s/arrow-webapp/ArrowWebService" % JAVA_DEPICT_HOST

def arrowToOrbital(smiles, arrowCodes):
    """Convert arrow representations into orbitals format"""
    mol = molBySmiles(smiles)
    orbInfoTuple = orbitalPairFromArrowStr(mol, arrowCodes)

    srcOrb, sinkOrb, nElectrons = orbInfoTuple

    srcInfoStr = srcOrb.formatInfoStr()
    sinkInfoStr = sinkOrb.formatInfoStr()

    return (srcInfoStr, sinkInfoStr)


def marvinXMLToDict(marvinXML):
    """Call the marvin webservice and get the respective smiles and arrow codes

    Arguments:
    - `marvinXML`:

    Returns a dict with:
    - smiles
    - arrowdesc
    - warnings
    """
    queryDict = {"action": "mrv2smi", "mrvxml": marvinXML}
    queryStr = urlencode(queryDict)

    theURL = CHEMAXON_MARVIN_URL + "?" + queryStr
    dataString = urlopen(theURL).read()

    res = dict(parse_qsl(dataString))

    if "smiles" not in res:
        res["smiles"] = ""
    if "arrowdesc" not in res:
        res["arrowdesc"] = ""
    if "warnings" not in res:
        res["warnings"] = ""

    return res
