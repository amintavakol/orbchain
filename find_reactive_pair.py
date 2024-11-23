import sys
import argparse

from orb_pair.Reaction import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_reactants", type=str, help="the reaction must be in the form of reactants>>product arrows")

    args = parser.parse_args()
    reaction = args.input_reactants
    try:
        smirks, arrows = reaction.split(" ")
    except:
        sys.exit("the reaction must be in the form of reactants>>product arrows")
    

    rxn = Reaction(smirks, arrows)
    print("reactive atom 1: {}".format(rxn.srcAtom))
    print("reactive MO 1: {}".format(rxn.srcOrb))
    print("reactive atom 2: {}".format(rxn.sinkAtom))
    print("reactive MO 2: {}".format(rxn.sinkOrb))

if __name__ == "__main__":
    sys.exit(main)
