"""
Some simple fitlers to remove some nonsense reactions which get put up on top.

"""

from chemutils.Common.Util import molBySmiles
from openeye.oechem import *

from chemutils.Common.Util import smi_to_unique_smi_fast



class Filters(object):
    """A filter to throw out nonsense reactions"""

    def __init__(self, filterList=None):
        """Setup default filterList"""
        self.filterList = filterList
        if self.filterList is None:
            self.setupDefaultFilters()

    def setupDefaultFilters(self):
        """Set some reasonable defaults"""
        self.filterList = [
            BadProductFilter(),
            LoneHydrogenFilter(),
        ]

    def __call__(self, op):
        """Run through each of the filters, if any return false, return false"""
        for f in self.filterList:
            val = f(op)
            if not val:
                return val
        return True


class BaseOPFilter(object):
    """Simple structure for how one of the filters should look"""

    def __init__(self):
        """Constructor"""

    def __call__(self, op):
        """Subclass must implement"""
        raise NotImplementedError()


class LoneHydrogenFilter(BaseOPFilter):
    """Return false if there is a lone hydrogen in the products"""

    def __call__(self, op):
        """Like it says above"""
        mol = molBySmiles(op.productSmiles)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1 and atom.GetDegree() == 0:
                return False

        return True


class BadProductFilter(BaseOPFilter):
    """Re: DVV request. Return false if there is a single unwanted product by formula"""

    def __init__(self):
        self.bad_product_formulae = set(["CH3+", "CH3-", "H+", "H-"])
        self.failed_orbpair = []

    def __call__(self, op):
        mol_parts = op.productSmiles.split(".")
        mol = OEGraphMol()
        for part in mol_parts:
            OEParseSmiles(mol, part)
            if OEMolecularFormula(mol) in self.bad_product_formulae:
                return False
                self.failed_orbpair.append(op)
            mol.Clear()

        return True


class OrbPairRules(object):
    """A filter to throw out nonsense reactions"""

    def __init__(self, rule_list=None):
        """Setup default filterList"""
        self.rule_list = rule_list
        if self.rule_list is None:
            self.setup_default_rules()

    def setup_default_rules(self):
        """Set some reasonable defaults"""
        self.rule_list = [
            TransitionState(),
            HydrogenSink(),
            SourceOrbitals(),
            TetrahedralIntermediate(),
            Pka_Pkb(),
            GeometryRule(),
            NegativeClBrI(),
            NeutralOxygen(),
        ]

    def __call__(self, op):
        """Run through each of the filters, if any return false, return false"""
        if (
            op.sinkAtom and op.srcAtom
        ):  # for pathway search we need to check this before applying the rules. Otherwise you will get into trouble!
            for f in self.rule_list:
                val = f(op)
                if not val:
                    return val
            return True


class TransitionState(object):
    """Transition State filter."""

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []
        self.MAX_SOURCE_SINK_DISTANCE = 6

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.src_sink_distance(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def src_sink_distance(self, op):
        """Only allow sources and sinks to combine
        if they form transition states of less than 7
        or 8 atoms. Maybe make it a checkbox such as
        allow 7 or 8 member transition states"""
        distance = OEGetPathLength(op.srcAtom.atom, op.sinkAtom.atom)
        if distance > self.MAX_SOURCE_SINK_DISTANCE:
            return True

        return False


class HydrogenSink(object):
    """Hydrogen Sink filter."""

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []
        self.ALLOWED_TRANSITION_MEMBER = [1, 2, 3, 6]

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.hydrogen_sink_transition_member(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def hydrogen_sink_transition_member(self, op):
        """only allow sources to combine with hydrogen sinks
        on the same molecule that form 5 and 6 member
        transition states.
        """
        if op.sinkAtom.atom.GetAtomicNum() == 1:
            if OEMolToSmiles(op.srcAtom.atom.GetParent()) == OEMolToSmiles(
                op.sinkAtom.atom.GetParent()
            ):
                distance = OEGetPathLength(op.srcAtom.atom, op.sinkAtom.atom)
                if distance in self.ALLOWED_TRANSITION_MEMBER:
                    return True

        return False


class SourceOrbitals(object):
    """Source Orbitals filter."""

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []
        self.atom_list = [7, 8, 9, 15, 16, 17, 35, 53]  # [N, O, F, P, S, Cl, Br, I]
        self.ALLOWED_ORBITALS = ["sp", "sp2", "sp3"]

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.srcOrbital_check(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def srcOrbital_check(self, op):
        """allow only sp3, sp2, and sp orbitals of
        O, S, N, P, F, Cl, Br and, I for SOURCES."""
        if op.srcAtom.atom.GetAtomicNum() in self.atom_list:
            if op.srcOrb.type == "sigma":
                distance = OEGetPathLength(op.srcAtom.atom, op.sinkAtom.atom)
                # if distance != 1 and op.srcAtom.atom != op.sinkAtom.atom:
                #    return True
                if distance == 1:
                    return False
                elif op.srcAtom.atom == op.sinkAtom.atom:
                    return False
                else:
                    return True

            elif op.srcOrb.type not in self.ALLOWED_ORBITALS:
                return True

        return False


class TetrahedralIntermediate(object):
    """Source Orbitals filter."""

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.is_tetrahedral_intermediate(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def is_tetrahedral_intermediate(self, op):
        """Only allow carbon tetrahedral intermediates sinks to combine with sources
        1 atom away (or themselves)"""
        if op.sinkAtom.atom.GetAtomicNum() == 6:
            carbon = op.sinkAtom.atom
            if self.is_carbon_tetrahedral_intermediate(carbon):
                distance = OEGetPathLength(op.srcAtom.atom, op.sinkAtom.atom)
                if distance == 1:
                    return False
                elif op.srcAtom.atom == op.sinkAtom.atom:
                    return False
                else:
                    return True

        return False

    def is_carbon_tetrahedral_intermediate(self, atom):
        """given an OEAtomBase object, returns true if it's
        a carbon tetrahedral intermediate"""
        # TODO: Now we check the tetrahedral intermediate structure
        # just fo CHNO. Later it has to cover all other structures.
        neighbors = {}
        # carbon degree must be 4:
        if atom.GetDegree() == 4:
            for bond in atom.GetBonds():
                nbr = bond.GetNbr(atom)
                neighbors[
                    nbr.GetAtomicNum(), nbr.GetFormalCharge(), nbr.GetDegree()
                ] = nbr
        # there are 6 substructures which should be checked for carbon
        # tetrahedral intermediate.

        # 1: O+ and N are neighbors. They are attached to some other groups.
        if (7, 0, 3) in neighbors and (8, 1, 3) in neighbors:
            # log.debug("First condition for tetrahedral intermediate is met!")
            return True

        # 2: O- is one of the neighbors and it has no other attachments (it's a leaf in mol tree).
        elif (8, -1, 1) in neighbors:
            # log.debug("Second condition for tetrahedral intermediate is met!")
            return True

        # 3: O+ and neutral O are the neighbors with degrees: 3, 2.
        elif (8, 0, 2) in neighbors and (8, 1, 3) in neighbors:
            # log.debug("Third condition for tetrahedral intermediate is met!")
            return True

        # 4: N+ and N are two of the neighbors and they are attached to some other groups.
        elif (7, 0, 3) in neighbors and (7, 1, 4) in neighbors:
            # log.debug("Fourth condition of tetrahedral intermediate is met!")
            return True

        # 5: Neutral N and neutral O connecting to other groups are the neighbors of this carbon.
        elif (7, 0, 3) in neighbors and (8, 0, 2) in neighbors:
            # log.debug("Fifth condition of tetrahedral intermediate is met!")
            return True

        # 6: Neutral O and N+ are neighbors. They are attached to other groups.
        elif (7, 1, 4) in neighbors and (8, 0, 2) in neighbors:
            # log.debug("Sixth condition of tetrahedral intermediate is met!")
            return True

        else:
            return False


class Pka_Pkb(object):
    """Filter based on the reactivity of different substructures."""

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []
        self.pka_dict_original = {
            "I[H:1]": -10,
            "Br[H:1]": -9,
            "Cl[H:1]": -8,
            "[OH2+][H:1]": -1.7,
            "O[H:1]": 15.7,
            "F[H:1]": 3.17,
            "[H:1][N+]1=CC=CC=C1": 5.25,
            "[H:1]C#N": 9.4,
            "[H:1]N=[N+]=[N-]": 4.27,
            "[H:1]SC#N": 4,
            "O=S(O[H:1])O": 1.9,
            "O=S([O-])O[H:1]": 7.21,
            "O=S(O[H:1])(O)=O": -3,
            "O=S([O-])(O[H:1])=O": 1.99,
            "O=P(O[H:1])(O)O": 2.12,
            "O=P([O-])(O)O[H:1]": 7.21,
            "O=P([O-])(O[H:1])[O-]": 12.32,
            "O=[N+]([O-])O[H:1]": -1.3,
            "O=NO[H:1]": 3.29,
            "O=S(O[H:1])(C)=O": -2.6,
            "O=S(O[H:1])(C(F)(F)F)=O": -14,
            "OO[H:1]": 11.6,
            "O=S(O[H:1])(C1=CC=C(C)C=C1)=O": -2.8,
            "O=C(O[H:1])C": 4.76,
            "O=C(O[H:1])C(Cl)(Cl)Cl": 0.65,
            "O=C(O[H:1])C(F)(F)F": -0.25,
            "O=CO[H:1]": 3.77,
            "O=C(O)O[H:1]": 3.6,
            "O=C(O[H:1])[O-]": 10.3,
            "O=C(C1=CC=CC=C1)O[H:1]": 4.2,
            "O=C(O)C(O[H:1])=O": 1.25,
            "O=C(O[H:1])C([O-])=O": 4.14,
            "CO[H:1]": 15.5,
            "CCO[H:1]": 15.9,
            "CC(C)O[H:1]": 16.5,
            "CC(C)(C)O[H:1]": 17,
            "[H:1]OC1=CC=CC=C1": 9.95,
            "COO[H:1]": 11.5,
            "O=C(OO[H:1])C": 8.2,
            "O=C(C1=CC(Cl)=CC=C1)OO[H:1]": 7.57,
            "CC[O+]([H:1])CC": -3.6,
            "[H:1][O+]1CCCC1": -2.05,
            "C[OH+][H:1]": -2.2,
            "CC[OH+][H:1]": -2.4,
            "[H:1][OH+]C(C)C": -2.8,
            "[H:1][OH+]C(C)(C)C": -3,
            "[H:1]O[N+]1=CC=CC=C1": 0.79,
            "[NH3+][H:1]": 9.2,
            "CC[NH2+][H:1]": 10.6,
            "C[NH2+][H:1]": 10.64,
            "CC[NH+](CC)[H:1]": 11.09,
            "C[NH+]([H:1])C": 10.73,
            "CC([NH+]([H:1])C(C)C)C": 11.05,
            "CC[N+:1](CC)([H])CC": 10.75,
            "CC[N+](C(C)C)([H:1])C(C)C": 10.75,
            "[H:1][NH2+]C1=CC=CC=C1": 4.6,
            "C[N+](C)([H:1])C1=CC=CC=C1": 5.2,
            "[H:1][NH+](C1=CC=CC=C1)C2=CC=CC=C2": 0.78,
            "[H:1][NH2+]C1=CC2=CC=CC=C2C=C1": 4.16,
            "N[NH2+][H:1]": 8.12,
            "O[NH2+][H:1]": 5.96,
            "[H:1][N+]12CCC(CC2)CC1": 11,
            "[H:1][NH+]1CCOCC1": 8.36,
            "[H:1][N+]1(C)CCOCC1": 7.38,
            "[H:1][N+]12CC[NH+](CC2)CC1": 2.97,
            "[H:1][N+]12CCN(CC2)CC1": 8.93,
            "[NH3+]CC[NH2+][H:1]": 6.9,
            "[H]N([H])CC[NH2+][H:1]": 9.95,
            "O=C1N([H:1])C(C2=CC=CC=C21)=O": 8.3,
            "O=C1N([H:1])C(CC1)=O": 9.5,
            "O=C1N([H:1])C(C=C1)=O": 8.8,
            "[H:1][N+]1=C2N(CCCCC2)CCC1": 13.5,
            "CN(C)C1=CC=[N+]([H:1])C=C1": 9.2,
            "[H:1][N+]1=CNC=C1": 6.95,
            "[H:1][N+]1=C(C(C)(C)C)C=CC=C1C(C)(C)C": 4.95,
            "CC(C)(C)C1=[N+]([H:1])C(C(C)(C)C)=CC(C)=C1": 4.6,
            "CC1=[N+]([H:1])C(C)=CC=C1": 6.75,
            "O=C(CC)O[H:1]": 4.88,
            "[SH2+][H:1]": 1.1,
            "[H]S[H:1]": 7,
            "O=C(O[H:1])C(C)(C)C": 5.01,
            "O=C(C(CC1=CC=CC=C1)[NH3+])O[H:1]": 1.83,
            "O=C(C(CC1=CC=CC=C1)[NH2+][H:1])[O-]": 9.13,
            "[H:1][NH+]1CCCC1": 11.27,
            "[H:1][NH+]1CCCCC1": 11.22,
            "CC([N+]([Li])([H:1])C(C)C)C": 36,
            "C[Si](C)(C)[N+]([Li])([H:1])[Si](C)(C)C": 26,
            "CC(C)(C)[O+]([H:1])[Li]": 16.99,
            "CC(C)(C)[O+]([H:1])[Na]": 16.99,
            "CC(C)(C)[O+]([H:1])[K]": 16.99,
            "CC(C)[O+]([H:1])[Li]": 16.49,
            "CC(C)[O+]([H:1])[Na]": 16.49,
            "CC(C)[O+]([H:1])[K]": 16.49,
            "CC[O+]([H:1])[Li]": 15.89,
            "CC[O+]([H:1])[Na]": 15.89,
            "CC[O+]([H:1])[K]": 15.89,
            "C[O+]([H:1])[Li]": 15.49,
            "C[O+]([H:1])[Na]": 15.49,
            "C[O+]([H:1])[K]": 15.49,
            "[Li][OH+][H:1]": 15.69,
            "[Na][OH+][H:1]": 15.69,
            "[K][OH+][H:1]": 15.69,
            "[H:1][N]([H])[H]": 32.2,
            "[H:1][N]([H])CC": 33.6,
            "[H:1][N]([H])C": 33.64,
            "CC[N](CC)[H:1]": 34.09,
            "C[N]([H:1])C": 33.73,
            "CC(C)[N](C(C)C)[H:1]": 34.05,
            "[H:1][N]([H])C1=CC=CC=C1": 27.6,
            "[H:1][N](C1=CC=CC=C1)C2=CC=CC=C2": 23.78,
            "[H:1][N]([H])C1=CC2=CC=CC=C2C=C1": 27.16,
            "N[N]([H:1])[H]": 31.12,
            "O[N]([H:1])[H]": 28.96,
            "[H:1][N]1CCOCC1": 31.36,
            "NCC[N]([H:1])[H]": 32.95,
            "[H:1][N]1CCCC1": 34.27,
            "[H:1][N]1CCCCC1": 34.22,
        }
        self.pkb_dict_original = {
            "[I-:1]": 24,
            "[Br-:1]": 23,
            "[Cl-:1]": 22,
            "[H][O:1][H]": 15.7,
            "[O-:1][H]": -1.7,
            "[F-:1]": 10.83,
            "C1=CC=CC=[N:1]1": 8.75,
            "[C-:1]#N": 4.6,
            "[N-:1]=[N+]=[N-]": 9.73,
            "[S-:1]C#N": 10,
            "O=S([O-:1])O": 12.1,
            "O=S([O-])[O-:1]": 6.79,
            "O=S([O-:1])(O)=O": 17,
            "O=S([O-])([O-:1])=O": 12.01,
            "O=P([O-:1])(O)O": 11.88,
            "O=P([O-:1])(O)[O-]": 6.79,
            "O=P([O-])([O-])[O-:1]": 1.68,
            "O=[N+]([O-])[O-:1]": 15.3,
            "O=N[O-:1]": 10.71,
            "O=S([O-:1])(C)=O": 16.6,
            "O=S([O-:1])(C(F)(F)F)=O": 28,
            "O[O-:1]": 2.4,
            "O=S([O-:1])(C1=CC=C(C)C=C1)=O": 16.8,
            "O=C([O-:1])C": 9.24,
            "O=C([O-:1])C(Cl)(Cl)Cl": 13.35,
            "O=C([O-:1])C(F)(F)F": 14.25,
            "O=C[O-:1]": 10.23,
            "O=C(O)[O-:1]": 10.4,
            "O=C([O-])[O-:1]": 3.7,
            "O=C(C1=CC=CC=C1)[O-:1]": 9.8,
            "O=C(O)C([O-:1])=O": 12.75,
            "O=C([O-])C([O-:1])=O": 9.86,
            "C[O-:1]": -1.5,
            "CC[O-:1]": -1.9,
            "CC(C)[O-:1]": -2.5,
            "CC(C)(C)[O-:1]": -3,
            "[O-:1]C1=CC=CC=C1": 4.05,
            "CO[O-:1]": 2.5,
            "O=C(O[O-:1])C": 5.8,
            "O=C(C1=CC(Cl)=CC=C1)O[O-:1]": 6.43,
            "CC[O:1]CC": 17.6,
            "C1CC[O:1]C1": 16.05,
            "C[O:1][H]": 16.2,
            "CC[O:1][H]": 16.4,
            "[H][O:1]C(C)C": 16.8,
            "[H][O:1]C(C)(C)C": 17,
            "[O-:1][N+]1=CC=CC=C1": 13.21,
            "[H][N:1]([H])[H]": 4.8,
            "[H][N:1]([H])CC": 3.4,
            "[H][N:1]([H])C": 3.36,
            "CC[N:1](CC)[H]": 2.91,
            "C[N:1]([H])C": 3.27,
            "CC(C)[N:1](C(C)C)[H]": 2.95,
            "CC[N:1](CC)CC": 3.25,
            "CC[N:1](C(C)C)C(C)C": 3.25,
            "[H][N:1]([H])C1=CC=CC=C1": 9.4,
            "C[N:1](C)C1=CC=CC=C1": 8.8,
            "[H][N:1](C1=CC=CC=C1)C2=CC=CC=C2": 13.22,
            "[H][N:1]([H])C1=CC2=CC=CC=C2C=C1": 9.84,
            "N[N:1]([H])[H]": 5.88,
            "O[N:1]([H])[H]": 8.04,
            "C1(CC2)CC[N:1]2CC1": 3,
            "[H][N:1]1CCOCC1": 5.64,
            "C[N:1]1CCOCC1": 6.62,
            "[H][N+]12CC[N:1](CC2)CC1": 11.03,
            "N12CC[N:1](CC2)CC1": 5.07,
            "[NH3+]CC[N:1]([H])[H]": 7.1,
            "NCC[N:1]([H])[H]": 4.05,
            "O=C1[N-:1]C(C2=CC=CC=C21)=O": 5.7,
            "O=C1[N-:1]C(CC1)=O": 4.5,
            "O=C1[N-:1]C(C=C1)=O": 5.2,
            "N1(CCCCC2)C2=[N:1]CCC1": 0.5,
            "CN(C)C1=CC=[N:1]C=C1": 4.8,
            "C1=C[N:1]=CN1": 7.05,
            "CC(C)(C)C1=[N:1]C(C(C)(C)C)=CC=C1": 9.05,
            "CC(C)(C)C1=[N:1]C(C(C)(C)C)=CC(C)=C1": 9.4,
            "CC1=[N:1]C(C)=CC=C1": 7.25,
            "O=C(CC)[O-:1]": 9.12,
            "CCC[C:1]([H])([H])[Li]": -36,
            "C[C:1](C)(C)[Li]": -37,
            "[H][S:1][H]": 12.9,
            "[S-:1][H]": 7,
            "O=C([O-:1])C(C)(C)C": 8.99,
            "O=C(C(CC1=CC=CC=C1)[NH3+])[O-:1]": 12.17,
            "O=C(C(CC1=CC=CC=C1)[N:1]([H])[H])[O-]": 4.87,
            "O=C(C[N:1]([H])[H])O": 5.4,
            "[H][N:1]1CCCC1": 2.73,
            "[H][N:1]1CCCCC1": 2.78,
            "CC([N:1]([Li])C(C)C)C": -22,
            "C[Si](C)(C)[N:1]([Li])[Si](C)(C)C": -12,
            "[Li][H:1]": -21,
            "[Na][H:1]": -21,
            "[K][H:1]": -21,
            "CC(C)(C)[O:1][Li]": -3.01,
            "CC(C)(C)[O:1][Na]": -3.01,
            "CC(C)(C)[O:1][K]": -3.01,
            "CC(C)[O:1][Li]": -2.51,
            "CC(C)[O:1][Na]": -2.51,
            "CC(C)[O:1][K]": -2.51,
            "CC[O:1][Li]": -1.91,
            "CC[O:1][Na]": -1.91,
            "CC[O:1][K]": -1.91,
            "C[O:1][Li]": -2.51,
            "C[O:1][Na]": -2.51,
            "C[O:1][K]": -2.51,
            "[Li][O:1][H]": -1.71,
            "[Na][O:1][H]": -1.71,
            "[K][O:1][H]": -1.71,
            "O=C(C(CC1=CC=CC=C1)[N:1]([H])[H])O": 5.87,
            "O=C(C(C)[NH3+])[O-:1]": 11.68,
            "O=C(C(C)[N:1]([H])[H])[O-]": 4.38,
            "O=C(C(C)[N:1]([H])[H])O": 5.38,
            "O=C(C1[NH2+]CCC1)[O-:1]": 12.01,
            "O=C(C(CCC1)[N:1]1[H])[O-]": 3.04,
            "O=C(C(CCC1)[N:1]1[H])O": 4.04,
            "O=C(C(C(C)C)[NH3+])[O-:1]": 11.65,
            "O=C(C(C(C)C)[N:1]([H])[H])[O-]": 4.31,
            "O=C(C(C(C)C)[N:1]([H])[H])O": 5.31,
            "O=C(C[NH3+])[O-:1]": 11.66,
            "O=C(C[N:1]([H])[H])[O-]": 4.4,
        }
        self.pka_dict = dict(
            zip(
                map(
                    lambda x: self.smi_to_unique_smi_map(x),
                    self.pka_dict_original.keys(),
                ),
                self.pka_dict_original.values(),
            )
        )
        self.pkb_dict = dict(
            zip(
                map(
                    lambda x: self.smi_to_unique_smi_map(x),
                    self.pkb_dict_original.keys(),
                ),
                self.pkb_dict_original.values(),
            )
        )
        self.canonicalized_pka = map(
            lambda x: smi_to_unique_smi_fast(x), self.pka_dict.keys()
        )
        self.canonicalized_pkb = map(
            lambda x: smi_to_unique_smi_fast(x), self.pkb_dict.keys()
        )
        self.canon_pka_dict = dict(zip(self.canonicalized_pka, self.pka_dict.values()))
        self.canon_pkb_dict = dict(zip(self.canonicalized_pkb, self.pkb_dict.values()))

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.pka_pkb_check(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def pka_pkb_check(self, op):
        """Among all hydrogen sinks, the one within the molecule with the lowest pka is favourable.
        Among all the sources attacking to a hydrogen sink, the one woth the lowest pkb is favourable.
        """
        if op.sinkAtom.atom.GetAtomicNum() == 1:
            clean_sink_mol = smi_to_unique_smi_fast(
                self.take_out_active_mol(op.sinkAtom.atom)
            )
            clean_source_mol = smi_to_unique_smi_fast(
                self.take_out_active_mol(op.srcAtom.atom)
            )
            sink_mol = self.take_out_active_mol(op.sinkAtom.atom)
            source_mol = self.take_out_active_mol(op.srcAtom.atom)
            sink_reactants = set((smi_to_unique_smi_fast(op.reactantSmiles).split(".")))
            sink_reactants.remove(clean_source_mol)
            # log.debug("clean source is removed")
            source_reactants = set(
                (smi_to_unique_smi_fast(op.reactantSmiles).split("."))
            )
            source_reactants.remove(clean_sink_mol)
            # log.debug("clean sink is removed")
            if (
                clean_source_mol in self.canonicalized_pkb
                and clean_sink_mol in self.canonicalized_pka
            ):
                return False
            elif (
                len(sink_reactants.intersection(set(self.canonicalized_pka))) >= 2
                and sink_mol in self.pka_dict.keys()
            ):
                # log.debug("more than one for sinks")
                # log.debug("sink:%s"%(sink_mol))
                sink_pka_scores = {}
                sink_pka = float(self.pka_dict[sink_mol])
                for mol in sink_reactants.intersection(set(self.canonicalized_pka)):
                    sink_pka_scores[mol] = self.canon_pka_dict[mol]
                if (
                    min(sink_pka_scores, key=lambda k: sink_pka_scores[k])
                    == clean_sink_mol
                ):
                    return False
                else:
                    return True
            elif (
                len(source_reactants.intersection(set(self.canonicalized_pkb))) >= 2
                and source_mol in self.pkb_dict.keys()
            ):
                # log.debug("more than one for sinks")
                # log.debug("sink:%s"%(sink_mol))
                source_pkb_scores = {}
                source_pkb = float(self.pkb_dict[source_mol])
                for mol in source_reactants.intersection(set(self.canonicalized_pkb)):
                    source_pkb_scores[mol] = self.canon_pkb_dict[mol]
                if (
                    min(source_pkb_scores, key=lambda k: source_pkb_scores[k])
                    == clean_source_mol
                ):
                    return False
                else:
                    return True

        return False

    def take_out_active_mol(self, active_atom):
        """A stupid method to take out the mol
        containing src or sink atom from an op object"""
        mol = active_atom.GetParent()
        for atom in mol.GetAtoms():
            if atom == active_atom:
                atom.SetMapIdx(1)
            else:
                atom.SetMapIdx(0)
        smiles = OEMolToSmiles(mol).split(".")
        for smi in smiles:
            small_mol = OEGraphMol()
            OESmilesToMol(small_mol, smi)
            for atom in small_mol.GetAtoms():
                if atom.GetMapIdx() == 1:
                    active_mol = smi
                    return active_mol

    def smi_to_unique_smi_map(self, smi):
        """Creates a unique smiles with aromaticity perceived and (atom stereo, bond stereo) omitted.
        Based on latest BP re: OEMolToSmiles()"""
        mol = OEGraphMol()
        OEParseSmiles(mol, smi)
        OEAssignAromaticFlags(mol)
        return OECreateSmiString(
            mol,
            OEOFlavor_ISM_Default - OEOFlavor_ISM_AtomStereo - OEOFlavor_ISM_BondStereo,
        )


class GeometryRule(object):
    def __init__(self):
        self.failed_orbpair = []
        self.funcs = [
            self.bredt_check,
            self.geometry_check_1,
            self.geometry_check_2,
            self.geometry_check_3,
        ]

    def __call__(self, op):
        if OEGetPathLength(op.srcAtom.atom, op.sinkAtom.atom):
            if self.geometry_check(op):
                self.failed_orbpair.append(op)
                return False
            else:
                return True
        else:
            return True

    def geometry_check(self, op):
        source = op.srcAtom.atom
        sink = op.sinkAtom.atom
        atoms, bonds = self.make_path_and_bonds(source, sink)
        funcs = [
            self.bredt_check,
            self.geometry_check_1,
            self.geometry_check_2,
            self.geometry_check_3,
        ]
        # funcs = [self.geometry_check_1]
        for fun in funcs:
            if fun(atoms, bonds):
                return True
        atoms, bonds = self.make_path_and_bonds(sink, source)
        for fun in self.funcs:
            if fun(atoms, bonds):
                return True
        return False

    def make_path_and_bonds(self, atom_1, atom_2):
        # atoms = [atom_1]
        atoms = []
        # index = 3
        # for atom in atom_1.GetParent().GetAtoms():
        #    if atom == atom_1:
        #        atom.SetMapIdx(1)
        #    elif atom == atom_2:
        #        atom.SetMapIdx(2)
        for atom in OEShortestPath(atom_1, atom_2):
            #    if atom.GetMapIdx() == 0:
            #        atom.SetMapIdx(index)
            #        index+=1
            atoms.append(atom)
        # atoms.append(atom_2) #we need source and sink to be at the first and end of the list
        bonds = [atoms[i].GetBond(atoms[i + 1]) for i in range(len(atoms) - 1)]

        return atoms, bonds

    def bredt_check(self, atoms, bonds):
        "return true for the substructures not obeying the Bredt's rule"
        if len(atoms) > 3:  # distance 2
            if OEAtomIsInRingSize(atoms[0], 6):
                if OEAtomIsInRingSize(atoms[1], 6):
                    # if OEAtomIsInRingSize(atoms[2], 6) and bonds[0].IsAromatic() and bonds[1].IsAromatic(): #aromatic ring
                    if (
                        OEAtomIsInRingSize(atoms[2], 6)
                        and OEBondIsInRingSize(bonds[0], 6)
                        and OEBondIsInRingSize(bonds[1], 6)
                    ):  # aromatic ring
                        return True
                    else:
                        return False
                else:
                    return False
            elif (
                OEAtomIsInRingSize(atoms[1], 6)
                and OEAtomIsInRingSize(atoms[2], 6)
                and bonds[0].GetOrder() == 1
                and OEBondIsInRingSize(bonds[1], 6)
                and OEBondIsInRingSize(bonds[2], 6)
            ):
                return True
            else:
                return False
        elif len(atoms) == 3:
            if (
                OEAtomIsInRingSize(atoms[0], 6)
                and OEAtomIsInRingSize(atoms[1], 6)
                and OEBondIsInRingSize(bonds[0], 6)
                and OEBondIsInRingSize(bonds[1], 6)
            ):
                return True
            else:
                return False
        else:
            return False

    def geometry_check_1(self, atoms, bonds):
        if len(atoms) == 3:
            if (
                OEAtomIsInRingSize(atoms[0], 6)
                and OEAtomIsInRingSize(atoms[1], 6)
                and OEBondIsInRingSize(bonds[0], 6)
                and bonds[1].GetOrder() == 1
            ):
                return True
            else:
                return False
        else:
            return False

    def geometry_check_2(self, atoms, bonds):
        if len(atoms) == 3:
            # for i in range(2):
            #    if atoms[i].IsInRing():
            #        return False
            if bonds[0].GetOrder() == 3 and bonds[1].GetOrder() == 1:
                return True
            else:
                return False
        elif len(atoms) == 4:
            # for i in range(2):
            #    if atoms[i].IsInRing():
            #        return False
            if bonds[0].GetOrder() == 3 and bonds[1].GetOrder() == 1:
                return True
            else:
                return False
        else:
            return False

    def geometry_check_3(self, atoms, bonds):
        if len(atoms) == 4:
            if (
                OEAtomIsInRingSize(atoms[1], 6)
                and OEAtomIsInRingSize(atoms[2], 6)
                and OEBondIsInRingSize(bonds[1], 6)
            ):
                return True
            else:
                return False
        else:
            return False


class NegativeClBrI(object):
    """
    Negative ClBrI attacks carbon.
    """

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []
        self.special_atoms = [17, 35, 53]  # Cl, Br, I

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.main(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def main(self, op):
        """Do not let Cl-, Br-, I-, attack carbon on neutral
        molecules if it makes less molecules in the product"""
        if (
            op.srcAtom.atom.GetAtomicNum() in self.special_atoms
            and op.srcAtom.atom.GetFormalCharge() == -1
            and op.sinkAtom.atom.GetAtomicNum() == 6
        ):
            net_charge = 0
            mol = op.sinkAtom.atom.GetParent()
            OEAddExplicitHydrogens(mol)
            for atom in mol.GetAtoms():
                net_charge += atom.GetFormalCharge()
            if (
                len(op.reactantSmiles.split(".")) > len(op.productSmiles.split("."))
                and net_charge == 0
            ):
                return True
            else:
                return False
        else:
            return False


class NeutralOxygen(object):
    """
    Neutral Oxygen attacks carbon in a non-positive molecule.
    """

    def __init__(self):
        """constructor"""
        self.failed_orbpair = []

    def __call__(self, op):
        """Returns True if the orbpair passes the filter"""
        if self.main(op):
            self.failed_orbpair.append(op)
            return False
        else:
            return True

    def main(self, op):
        """Do not let neutral oxygen attack carbon on neutral
        or negative molecules if it makes less molecules in
        the product.
        (NOTE:we might need to make exceptions for anhydrides
        later on but we can have this for now.)"""
        if (
            op.srcAtom.atom.GetAtomicNum() == 8
            and op.srcAtom.atom.GetFormalCharge() == 0
            and op.sinkAtom.atom.GetAtomicNum() == 6
        ):
            net_charge = 0
            mol = op.sinkAtom.atom.GetParent()
            OEAddExplicitHydrogens(mol)
            for atom in mol.GetAtoms():
                net_charge += atom.GetFormalCharge()
            if (
                len(op.reactantSmiles.split(".")) > len(op.productSmiles.split("."))
                and net_charge <= 0
            ):
                return True
            else:
                return False
        else:
            return False
