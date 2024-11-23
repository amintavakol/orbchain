"""Various constants for use by the reaction processing modules"""

import logging
import chemutils.CombiCDB.Env as Env

"""Not functionally important, just indicates an estimate of the number
of lines in input files to hint a proper scale for the progress indicators.
If you do not wish to see those dot indicators, set this value to 0."""
EST_INPUT = Env.EST_INPUT

"""Application name, for example to identify a common logger object"""
APPLICATION_NAME = "CHEM.CombiCDB.proposers.app"

"""Default level for application logging.  Modify these for different scenarios.  See Python logging package documentation for more information"""
LOGGER_LEVEL = Env.LOGGER_LEVEL

"""Default format of logger output"""
LOGGER_FORMAT = "[%(asctime)s %(levelname)s] %(message)s"

