"""Standard values for different molecular and atomic properties
and some accessor functions.
"""


def nStdValenceElectrons(atomicNum):
    """Return the standard number of valence electrons
    for an atom type specified by atomic number.
    Based on lookup table for now, can come up with more generalized scheme later.
    """
    return STD_VALENCE[atomicNum]


def nMaxHyperValenceElectrons(atomicNum):
    """Return the maximum number of hyper electrons
    for an atom type specified by atomic number.

    The basics are.  If the PERIODIC_ROW is great enough, return 4.
    otherwise, return 0.

    This allows for hexa-coordinated species.  ()

    """
    if PERIODIC_ROW[atomicNum] >= 3:
        # if PERIODIC_ROW[atomicNum] >= 3 or atomicNum == 3: Li may require special handling per DVV requests
        return 4
    else:
        return 0


def atomIonizationEnergy(atom, formalCharge=None):
    """Return a reference value for the atom's ionization energy (kcal/mol).
    If the atom has a formal charge, should return an adjusted value
    to estimate 2nd ionization energies, etc.
    If a formal charge parameter is specified use that value instead of
    any set on the actual atom object.
    """
    atomicNum = None
    if type(atom) == int:
        # Just specified an atomic number
        atomicNum = atom
        if formalCharge is None:
            formalCharge = 0
    else:
        # Assume is an atom object
        atomicNum = atom.GetAtomicNum()
        if formalCharge is None:
            formalCharge = atom.GetFormalCharge()

    # Find the standard value for the neutral atom
    stdValue = STD_IONIZATION_ENERGY[atomicNum]

    # Heuristic, double or halve the value for each formal charge.
    #   For each positive charge, should be even harder to ionize the next electron
    #   For each negative charge, should be easier to ionize the "extra" electron
    value = stdValue * 2 ** (formalCharge)
    return value


def atomElectronAffinity(atom, formalCharge=None):
    """Return a reference value for the atom's electron affinity (kcal/mol).
    If the atom has a formal charge, should return an adjusted value
    to estimate 2nd electron affinities, etc.
    If a formal charge parameter is specified use that value instead of
    any set on the actual atom object.
    """
    atomicNum = None
    if type(atom) == int:
        # Just specified an atomic number
        atomicNum = atom
        if formalCharge is None:
            formalCharge = 0
    else:
        # Assume is an atom object
        atomicNum = atom.GetAtomicNum()
        if formalCharge is None:
            formalCharge = atom.GetFormalCharge()

    # Find the standard value for the neutral atom
    stdValue = STD_ELECTRON_AFFINITY[atomicNum]

    # Heuristic, double or halve the value for each formal charge.
    #   For each positive charge, should have greater affinity to replace an electron
    #   For each negative charge, should have less affinity for more electrons
    # This is weird or unreliable because electron affinity
    #   can attain both positive and negative values.
    value = stdValue * 2 ** (formalCharge)
    # value = stdValue;

    return value


def atomElectronegativity(atom, formalCharge=None):
    """Return a reference value for the electronegativity
    of a given atom.  This uses a general Mulliken scale where
    EN = (IE - EA) / 2,
    (IE = Ionization Energy, EA = Electron Affinity) which can in turn
    be found by standard values and inferred for charged atoms.
    """
    atomicNum = None
    if type(atom) == int:
        # Just specified an atomic number
        atomicNum = atom
        if formalCharge is None:
            formalCharge = 0
    else:
        # Assume is an atom object
        atomicNum = atom.GetAtomicNum()
        if formalCharge is None:
            formalCharge = atom.GetFormalCharge()

    return (
        atomIonizationEnergy(atomicNum, formalCharge)
        - atomElectronAffinity(atomicNum, formalCharge)
    ) / 2

    """    
    stdValue = STD_ELECTRONEGATIVITY[atom.GetAtomicNum()];
    # Arbitrary modifier to make charge count for a lot
    stdValue += atom.GetFormalCharge()*100;
    return stdValue
    """


def stableValenceShell(atomicNum):
    """How many electrons should be in the outer shell for this atom type to be stable?
    Most commonly should just be 8 to make the octet rule, except for hydrogen.
    Note also that essentially any can become stable with "zero" outer shell electrons,
    since that just means the next complete inner shell becomes the outer shell.

    Stable valence shell is the same for hypervalent possible compounds.
    """
    # if atomicNum in (1, 3, 11, 19):
    if atomicNum == 1:  # removed Li to allow add'l bonds per DVV
        # H
        return 2
    elif atomicNum in (3, 11, 19):
        # Li, Na, K
        return 6
    else:
        return 8


#####################################################################
#################### Begin Constants Section ########################
#####################################################################


# Permitivity of free space, electric constant.  (Units: C^2 / N m^2)
EPSILON_0 = 8.854e-12

# Elementary charge.  Magnitude of charge for one proton or electron (Units: C = Coulomb)
E_CHARGE = 1.60217653e-19

# Conversion factor: 1 kJ/mol = 0.010364 eV/atom
EV_PER_ATOM_VS_KJ_PER_MOL = 0.010364
# Conversion factor: 4.1868 kJ / 1 kcal
KJ_PER_KCAL = 4.1868

# Boltzmann constant to relate temperature to energy
#   note that this is proportional to the universal gas constant (R), just different units
BOLTZMANN_CONST_KJ = 1.3806505e-026
# kJ / Kelvin
BOLTZMANN_CONST_KCAL = BOLTZMANN_CONST_KJ / KJ_PER_KCAL
BOLTZMANN_CONST = BOLTZMANN_CONST_KCAL

# Planck's constant (h)
PLANCKS_CONST_J = 6.6260689633e-034
# Joules * sec
PLANCKS_CONST_KJ = PLANCKS_CONST_J / 1000
PLANCKS_CONST_KCAL = PLANCKS_CONST_KJ / KJ_PER_KCAL
PLANCKS_CONST = PLANCKS_CONST_KCAL


# Avogadro's number
ATOM_PER_MOL = 6.0221415e23

# Standard temperatures (absolute temperature in kelvin)
ROOM_TEMPERATURE = 298
BODY_TEMPERATURE = 310
STD_TEMPERATURE = 273  # Freezing point of water
DRY_ICE_BATH_TEMPERATURE = 195  # Sublimation of dry ice

# standard octet number which is 8
STD_OCTECT = 8

# Valence and Orbital Electron Functions
STD_VALENCE = {
    0: 1,  # *
    1: 1,  # H
    2: 2,  # He
    3: 1,  # Li
    4: 2,  # Be
    5: 3,  # B
    6: 4,  # C
    7: 5,  # N
    8: 6,  # O
    9: 7,  # F
    10: 8,  # Ne
    11: 1,  # Na
    12: 2,  # Mg
    13: 3,  # Al
    14: 4,  # Si
    15: 5,  # P
    16: 6,  # S
    17: 7,  # Cl
    18: 8,  # Ar
    19: 1,  # K
    20: 2,  # Ca
    21: 2,  # Sc
    22: 4,  # Ti
    25: 7,  # Mn?
    26: 3,  # Fe?
    29: 1,  # Cu?
    30: 2,  # Zn?
    33: 5,  # As?
    34: 6,  # Se
    35: 7,  # Br
    46: 2,  # Pd?
    47: 1,  # Ag
    50: 4,  # Sn
    55: 1,  # Cs
    53: 7,  # I
    76: 8,  # Os?
}


# Just as STD_VALENCE maps to columns of the periodic table, use this to identify trends down rows of the table
PERIODIC_ROW = {
    1: 1,  # H
    2: 1,  # He
    3: 2,  # Li
    4: 2,  # Be
    5: 2,  # B
    6: 2,  # C
    7: 2,  # N
    8: 2,  # O
    9: 2,  # F
    10: 2,  # Ne
    11: 3,  # Na
    12: 3,  # Mg
    13: 3,  # Al
    14: 3,  # Si
    15: 3,  # P
    16: 3,  # S
    17: 3,  # Cl
    18: 3,  # Ar
    19: 4,  # K
    20: 4,  # Ca
    21: 4,  # Sc
    22: 4,  # Ti
    25: 4,  # Mn?
    26: 4,  # Fe?
    29: 4,  # Cu?
    30: 4,  # Zn?
    33: 4,  # As?
    34: 4,  # Se
    35: 4,  # Br
    46: 5,  # Pd
    50: 5,  # Sn
    53: 5,  # I
    55: 6,  # Cs
    76: 6,  # Os?
}

# Standard Pauling scale electronegativies, taken from http://en.wikipedia.org/wiki/Electronegativity
STD_ELECTRONEGATIVITY = {
    1: 2.20,  # H
    2: 5.50,  # He
    3: 0.98,  # Li
    4: 1.57,  # Be
    5: 2.04,  # B
    6: 2.55,  # C
    7: 3.04,  # N
    8: 3.44,  # O
    9: 3.98,  # F
    11: 0.93,  # Na
    12: 1.31,  # Mg
    13: 1.61,  # Al
    14: 1.90,  # Si
    15: 2.19,  # P
    16: 2.58,  # S
    17: 3.16,  # Cl
    18: 3.20,  # Ar
    19: 0.82,  # K
    20: 1.00,  # Ca
    21: 1.36,  # Sc
    22: 1.54,  # Ti
    23: 1.63,  # V
    24: 1.66,  # Cr
    25: 1.55,  # Mn
    26: 1.83,  # Fe
    27: 1.88,  # Co
    28: 1.91,  # Ni
    29: 1.90,  # Cu
    30: 1.65,  # Zn
    31: 1.81,  # Ga
    32: 2.01,  # Ge
    33: 2.18,  # As
    34: 2.55,  # Se
    35: 2.96,  # Br
    36: 3.00,  # Kr
    37: 0.82,  # Rb
    38: 0.95,  # Sr
    39: 1.22,  # Y
    40: 1.33,  # Zr
    41: 1.60,  # Nb
    42: 2.16,  # Mo
    43: 1.90,  # Tc
    44: 2.20,  # Ru
    45: 2.28,  # Rh
    46: 2.20,  # Pd
    47: 1.93,  # Ag
    48: 1.69,  # Cd
    49: 1.78,  # In
    50: 1.96,  # Sn
    51: 2.05,  # Sb
    52: 2.10,  # Te
    53: 2.66,  # I
    54: 2.60,  # Xe
    55: 0.79,  # Cs
    76: 2.2,  # Os
}

IONIC_EN_DIFF = 1.7
POLAR_EN_DIFF = 0.4

# Standard (1st) Ionization Energy values (kJ/mol) adapted from http://hyperphysics.phy-astr.gsu.edu/hbase/chemical/ionize.html
## 1 kJ/mol = 0.010364 eV/atom
STD_IONIZATION_ENERGY_KJ = {
    1: 1192.0,  # H     # Hacked value!  Original was 1312.  Trying to make EN(H) just less than EN(C)
    2: 2372.3,  # He
    3: 520.26,  # Li
    4: 899.46,  # Be
    5: 800.66,  # B
    6: 1166.5,  # C     # Hacked value!  See above.  Original was 1086.5
    7: 1402.4,  # N
    8: 1314.0,  # O
    9: 1681.0,  # F
    10: 2080.7,  # Ne
    11: 495.85,  # Na
    12: 737.75,  # Mg
    13: 577.58,  # Al
    14: 786.47,  # Si
    15: 1011.8,  # P
    16: 999.61,  # S
    17: 1251.2,  # Cl
    18: 1520.6,  # Ar
    19: 418.85,  # K
    20: 589.83,  # Ca
    21: 633.1,  # Sc
    22: 658.8,  # Ti
    26: 759.36,  # Fe
    29: 745.45,  # Cu
    30: 906.4,  # Zn
    33: 947.0,  # As
    34: 941.0,  # Se
    35: 1140.0,  # Br
    47: 731.0,  # Ag
    50: 708.61,  # Sn
    53: 1008.4,  # I
    55: 375.7,  # Cs
    76: 839.4,  # Os
}
# Recalculate into kcal / mol units
STD_IONIZATION_ENERGY_KCAL = dict()
for key, value in STD_IONIZATION_ENERGY_KJ.items():
    STD_IONIZATION_ENERGY_KCAL[key] = value / KJ_PER_KCAL

STD_IONIZATION_ENERGY = STD_IONIZATION_ENERGY_KCAL

# Standard Electron Affinity values (kJ/mol) taken from http://en.wikipedia.org/wiki/Electron_affinity
STD_ELECTRON_AFFINITY_KJ = {
    1: -73,  # H
    2: 21,  # He
    3: -60,  # Li
    4: 19,  # Be
    5: -27,  # B
    6: -122,  # C
    7: 7,  # N
    8: -141,  # O
    9: -328,  # F
    10: 29,  # Ne
    11: -53,  # Na
    12: 19,  # Mg
    13: -43,  # Al
    14: -134,  # Si
    15: -72,  # P
    16: -200,  # S
    17: -349,  # Cl
    18: 35,  # Ar
    19: -48,  # K
    20: 10,  # Ca
    21: -18,  # Sc
    22: -8,  # Ti
    23: -51,  # V
    24: -64,  # Cr
    26: -16,  # Fe
    27: -64,  # Co
    28: -112,  # Ni
    29: -118,  # Cu
    30: 47,  # Zn
    31: -29,  # Ga
    32: -116,  # Ge
    33: -78,  # As
    34: -195,  # Se
    35: -325,  # Br
    36: 39,  # Kr
    37: -47,  # Rb
    39: -30,  # Y
    40: -41,  # Zr
    41: -86,  # Nb
    42: -72,  # Mo
    43: -53,  # Tc
    44: -101,  # Ru
    45: -110,  # Rh
    46: -54,  # Pd
    47: -126,  # Ag
    48: 32,  # Cd
    49: -29,  # In
    50: -116,  # Sn
    51: -103,  # Sb
    52: -190,  # Te
    53: -295,  # I
    54: 41,  # Xe
    55: -45,  # Cs
    73: -31,  # Ta
    74: -79,  # W
    75: -14,  # Re
    76: -106,  # Os
    77: -151,  # Ir
    78: -205,  # Pt
    79: -223,  # Au
    80: 61,  # Hg
    81: -20,  # Tl
    82: -35,  # Pb
    83: -91,  # Bi
    84: -183,  # Po
    85: -270,  # At
    86: 41,  # Rn
}
# Recalculate into kcal / mol units
STD_ELECTRON_AFFINITY_KCAL = dict()
for key, value in STD_ELECTRON_AFFINITY_KJ.items():
    STD_ELECTRON_AFFINITY_KCAL[key] = value / KJ_PER_KCAL
STD_ELECTRON_AFFINITY = STD_ELECTRON_AFFINITY_KCAL

# Calculated atomic radii, expressed in Angstroms (10^-10 m)
#   http://en.wikipedia.org/wiki/Atomic_radius
#   E. Clementi, D.L.Raimondi, and W.P. Reinhardt, J. Chem. Phys. 1967, 47, 1300
STD_ATOMIC_RADIUS = {
    1: 0.23,  # H   ?
    2: 0.31,  # He
    3: 1.67,  # Li
    4: 1.12,  # Be
    5: 0.87,  # B
    6: 0.67,  # C
    7: 0.56,  # N
    8: 0.48,  # O
    9: 0.42,  # F
    10: 0.38,  # Ne
    11: 1.90,  # Na
    12: 1.45,  # Mg
    13: 1.18,  # Al
    14: 1.11,  # Si
    15: 0.98,  # P
    16: 0.88,  # S
    17: 0.79,  # Cl
    18: 0.71,  # Ar
    19: 2.43,  # K
    20: 1.94,  # Ca
    21: 1.84,  # Sc
    22: 1.76,  # Ti
    23: 1.71,  # V
    24: 1.66,  # Cr
    25: 1.61,  # Mn
    26: 1.56,  # Fe
    27: 1.52,  # Co
    28: 1.49,  # Ni
    29: 1.45,  # Cu
    30: 1.42,  # Zn
    31: 1.36,  # Ga
    32: 1.25,  # Ge
    33: 1.14,  # As
    34: 1.03,  # Se
    35: 0.94,  # Br
    36: 0.88,  # Kr
    37: 2.65,  # Rb
    38: 2.19,  # Sr
    39: 2.12,  # Y
    40: 2.06,  # Zr
    41: 1.98,  # Nb
    42: 1.90,  # Mo
    43: 1.83,  # Tc
    44: 1.78,  # Ru
    45: 1.73,  # Rh
    46: 1.69,  # Pd
    47: 1.65,  # Ag
    48: 1.61,  # Cd
    49: 1.56,  # In
    50: 1.45,  # Sn
    51: 1.33,  # Sb
    52: 1.23,  # Te
    53: 1.15,  # I
    54: 1.08,  # Xe
    55: 2.98,  # Cs
    56: 2.53,  # Ba
    72: 2.08,  # Hf
    73: 2.00,  # Ta
    74: 1.93,  # W
    75: 1.88,  # Re
    76: 1.85,  # Os
    77: 1.80,  # Ir
    78: 1.77,  # Pt
    79: 1.74,  # Au
    80: 1.71,  # Hg
    81: 1.56,  # Tl
    82: 1.54,  # Pb
    83: 1.43,  # Bi
    84: 1.35,  # Po
    86: 1.20,  # Rn
}

# Organometallic forming atoms which are not configurationally stable
METAL_ATOMIC_NUMS = {
    3: True,  # Li
    11: True,  # Na
    12: True,  # Mg
    19: True,  # K
    20: True,  # Ca
    29: True,  # Cu
}
