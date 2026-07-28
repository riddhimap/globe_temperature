"""
globe_temperature
-----------------

A Python package for calculating Globe Temperature (Tg), 
Natural Wet Bulb Temperature (Tnwb), and Wet Bulb Globe Temperature (WBGT) 
from climate model data.

"""

from .wbgt import Tg_10mwind, Tnwb_10mwind, WBGT

__version__ = "1.0.1"
__author__ = "Riddhima Puri"

__all__ = [
    "Tg_10mwind",
    "Tnwb_10mwind",
    "WBGT",
]
