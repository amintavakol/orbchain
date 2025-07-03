"""Various constants for use by the reaction processing modules"""

from chemutils.CombiCDB import Env

"""Not functionally important, just indicates an estimate of the number
of lines in input files to hint a proper scale for the progress indicators.
If you do not wish to see those dot indicators, set this value to 0."""
EST_INPUT = Env.EST_INPUT

"""Tag to mark comments in processed files"""
COMMENT_TAG = "#"

"""SMIRKS string representing a separation of reactants and products"""
REACTION_DELIM = ">>"

"""SMIRKS string representing a separation of reactants, reagents, products"""
SUB_REACTION_DELIM = ">"

"""Delimiter to indicate separate molecule list components.
Maybe this should be the same as SMILES_MOL_DELIM, or it should be a non-SMILES character."""
MOL_LIST_DELIM = ","

"""Delimiter for items in a list of reaction profiles"""
REACTION_PROFILE_DELIM = "\n"

"""Default delimiter for parameter value strings"""
VALUE_DELIM = ","

"""Label for reaction index"""
REACTION_LABEL = "SMIRKS["
REACTION_LABEL_END = "]"

"""Label for reactant index list"""
REACTANT_LIST_LABEL = "Reactants"

"""Application name, for example to identify a common logger object"""
APPLICATION_NAME = "CHEM.CombiCDB.app"

"""Default level for application logging.  Modify these for different scenarios.  See Python logging package documentation for more information"""
LOGGER_LEVEL = Env.LOGGER_LEVEL

"""Default format of logger output"""
LOGGER_FORMAT = "[%(asctime)s %(levelname)s] %(message)s"

"""Similarity threshold for finding related products in a retro-synthesis search.
Not actively used right now.
"""
SIMILARITY_THRESHOLD = 0.8

"""Sentinel value to change a charged atom into an uncharged neutral atom.  
Change the atom into one with a charge of this sentinel value
and we'll specially override it.  Actually, further testing indicates 
this is unncessary, can specify neutralization of charge in a SMIRKS
string by explcitly stating +0"""
SENTINEL_NEUTRAL_CHARGE = +100

"""Similarly, sentinel value to represent radicals"""
SENTINEL_RADICAL_CHARGE = +200

"""Other sentinel values to be used for whatever developer purpose"""
SENTINEL_CHARGE = +50
SENTINEL_LEAVING_GROUP_CHARGE = +300

"""Sentinel charges for proton transfer (acid-base) reagents.
Should really only need the acid (HA) and conjugate acid (HB+)
indicators, since the others can be derived by just "deprotonation."
"""
SENTINEL_CONJUGATE_BASE = +299
#  A-
SENTINEL_ACID = +298
# HA
SENTINEL_BASE = +301
#  B
SENTINEL_CONJUGATE_ACID = +302
# HB+


"""Special sentinel value, if appears in product, 
means a unimolecular step that was optional.  So far, specifically designed for allylic resonance.
"""
SENTINEL_OPTIONAL_STEP_CHARGE = +75

"""Sentinel charge value label for molecules not to be presented as part of the product.
For example, used to label leaving groups"""
SENTINEL_REJECT_CHARGE = +400

"""Sentinel charge value label for molecules to be removed from any products
as soon as they are encountered, not just at the end when products are collected.
For example, deprotonating bases"""
SENTINEL_REJECT_IMMEDIATE_CHARGE = +350

"""Mock atom-mapping labels.  Use extreme charges instead"""
SENTINEL_LABEL_CHARGE = +1000

"""Threshold for "normally" expected absolute charges"""
NORMAL_CHARGE_THRESHOLD = 3

"""Default maximum number of tries to regenerate a synthesis pathway if 
to try and find one with as many good / interesting steps as possible."""
MAX_TRIES = 2

"""For standard reagents, assume looking at organic chemistry and thus products
must contain at least one carbon atom to be of interest."""
REQUIRED_ATOMIC_NUMS = [6]

"""Warning level constants regarding misuse / unexpected consequences of reactions"""
OK = 0
# Can use anywhere, no restrictions
CAUTION = 1
# Can use anywhere, but synthesis generator will avoid if possible
WARNING = 2
# Can use in experiments, but synthesis generator won't use it
ERROR = 3
# Can't use in synthesis generator or experiments.  Invalidates any other products.
DISFAVORED = 4
# Can't use in synthesis generator or experiments.  Other products may still be valid, this could be presented informationally to explain why "No Products Predicted"
RETRO_ONLY = 5
# Completely ignore unless flagged as a retro-reagent

"""Maximum number of syntheses for each request pattern that should be cached"""
SYNTHESIS_CACHE_MAX = 100

"""Minimum number of syntheses for each request pattern that would be happy with
before start throttling down new synthesis generation.
"""
SYNTHESIS_CACHE_MIN = 3

"""Number of old records to delete when synthesis cache exceeds maximum size"""
SYNTHESIS_CACHE_SCALEBACK = 20


"""ID of Garbage / Trash Can synthesis / reaction list to discard unwanted / unused reactions"""
GARBAGE_ID = -2

"""Sentinel value for a reaction step index to indicate no specific step"""
SENTINEL_STEP_INDEX = -1

"""Reaction synthesis code if is just a temporarily cached problem"""
CACHE_CODE = "Cache"

"""Reaction synthesis code if the synthesis was reconstructed from a completed problem,
which leaves some question as to it's accuracy.
"""
RECONSTRUCTED_CODE = "Reconstructed"

"""Reaction synthesis code if the synthesis is to be disabled, like a soft delete.
"""
DISABLED_CODE = "Disabled"

"""Class label for example reaction synthesis problems."""
SAMPLE_SYNTHESIS_CLASS = "Synthesis"

"""Class label for reaction "syntheses" that are just containers for many sample reactions"""
SAMPLE_REACTION_CLASS = "Example"

"""Class label for reaction mechanism "synthesis" that are just containers for many reactions with sample mechanism problems"""
SAMPLE_MECHANISM_CLASS = "Mechanism"

"""Class label for reaction lists custom authored by users"""
USER_REACTION_CLASS = "UserRecipe"

"""Suffix appended to reagent reaction profile attributes to indicate a list of these occur under one reagent"""
LIST_SUFFIX = "_list"

"""Reaction category position scaling factors.
For example, if the "position" for a category is 2107, and
the scaling factor is 100.0, then this will actually be
presented to the user as 21.07 (like chapter 21, section 7).
"""
REACTION_CATEGORY_POSITION_SCALE = 100.0
REACTION_CATEGORY_POSITION_FORMAT = "%5.2f"


"""Some Constants to Denote whether stereo info is for extSource or extTarget"""
EXT_TYPE_UNDEFINED = None
EXT_TYPE_SOURCE = "source"
EXT_TYPE_TARGET = "target"

"""Decryption key for possibly encrypted reaction mdoel SMIRKS strings"""
ENCRYPTION_KEY = "ReactionExplorer2010"
