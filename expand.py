import sys
import argparse

from orb_pair.utils import *
from atom.utils import *
from orb_pair.modules.simple_orbpair_object import SimpleOrbPairObject as SOO
from atom.modules.simple_atom_object import SimpleAtomObject as SAO
from chemutils.Common.Util import smi_to_unique_smi_fast, smi_to_unique_smi_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_reactants", type=str, help="the reactants for which you want to generate all the "
                                                          "products.")
    parser.add_argument("reaction_type", help="polar or radical")

    args = parser.parse_args()
    if args.input_reactants:
        reactants = args.input_reactants
    else:
        reactants = "CCC[CH2].CCCCCC"

    prods = expand(reactants, args.reaction_type)

    for item in prods:
        reaction = reactants + "" + items[0] + " " + items[1]
        print(reaction)


def expand(reactants, reaction_type="polar"):
    """
    Given a set of reactant molecules, extract all
    source and sink atoms and generate all the plausible
    reactions.
    NOTE: works for polar mechanisms now.

    Args:
        reactants (str): SMILES string delim with ".".
        reaction_type (str): polar or radical.

    returns:
        (list): list of all the generated SMILES strings.
    """
    atoms = SAO.atomObjFromReactantSmi(reactants)
    if reaction_type.lower() == "polar":
        ops = SOO.orbPairObjectsFromAtoms(atoms, atoms)
    elif reaction_type.lower() == "radical":
        ops = SOO.orbPairObjectsFromAtoms(atoms, atoms, radical=True)
    else:
        sys.exit("reaction type is not supported.")

    return [(op.productSmiles, op.arrowCodes) for op in ops]


if __name__ == "__main__":
   sys.exit(main()) 
