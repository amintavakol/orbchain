import argparse

from ranker.utils import *
from atom.utils import *
from ranker.modules.simple_orbpair_object import SimpleOrbPairObject as SOO
from atom.modules.simple_atom_object import SimpleAtomObject as SAO
from chemutils.Common.Util import smi_to_unique_smi_fast, smi_to_unique_smi_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_reactants", type=str, help="the reactants for which you want to generate all the "
                                                          "products.")
    parser.add_argument('--rand_one', action='store_true', default=False, help="")
    # parser.add_argument("reaction_type", help="either polar or radical")

    args = parser.parse_args()
    if args.input_reactants:
        reactants = args.input_reactants
    else:
        reactants = "[NH2:1]CCCC(NC)=O.C[CH+:2]C1=CC=C(C(OC)=O)C=C1"

    start_time = time.time()
    prods = expand(reactants, args.rand_one)
    # prods = expand(args.input_reactants, reaction_type)

    # print(prods[random.randint(1, len(prods))])
    # print(prods)
    print("Number of products: %i" % len(prods))
    print("Time: %s" % (time.time() - start_time))


def expand(reactants, rand_one: bool):
    """
    Given a set of reactant molecules, extract all
    source and sink atoms and generate all the plausible
    reactions.
    NOTE: works for polar mechanisms now.

    Args:
        reactants (str): SMILES string delim with ".".
        rand_one (bool): if true will just randomly get one

    returns:
        (list): list of all the generated SMILES strings.
    """
    atoms = SAO.atomObjFromReactantSmi(reactants)
    ops = SOO.orbPairObjectsFromAtoms(atoms, atoms, rand_one=rand_one)

    return [op.productSmiles for op in ops]


if __name__ == "__main__":
   sys.exit(main()) 
