import rdkit
import numpy as np
from collections import Counter
import random

from rdkit import Chem

import logging
import sys

from chemutils.Common.Util import  molBySmiles
from openeye.oechem import OEAssignAromaticFlags, OEAroModelMMFF, OEClearAromaticFlags, OEKekulize
from openeye.oechem import OECanonicalOrderAtoms, OECanonicalOrderBonds
from chemutils.Common.CanonicalAtomMapSmiles import createCanonicalAtomMapSmiString
from chemutils.Common.MolExt import removeNonsenseStereo


BOS={
     Chem.BondType.SINGLE:1.0,
     Chem.BondType.DOUBLE:2.0,
     Chem.BondType.TRIPLE:3.0,
     Chem.BondType.AROMATIC:1.5,
     Chem.BondType.UNSPECIFIED:0.0
     }
    
ELEMENTS = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 
             'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 
			 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb', 'W', 
			 'Ru', 'Nb', 'Re', 'Te', 'Rh', 'Tc', 'Ba', 'Bi', 'Hf', 'Mo', 'U', 'Sm', 'Os', 'Ir', 
			 'Ce','Gd','Ga','Cs', 'unknown']

def bond_to_string(bond):
    """
    bond is an rdkit bond.
    """
    strings = {1.0: "-", 2.0: "=", 3.0: "#", 1.5: "a", 0.0: "u"}
    order = BOS[bond.GetBondType()]
    return strings[order]

def mol_with_hydrogens(smi):
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    mol = Chem.MolFromSmiles(smi, ps)
    mol = Chem.AddHs(mol)
    return mol

def prepare_atom_indicator(m, idx1):
    a = m.GetAtomWithIdx(idx1)
    if a.GetFormalCharge()==0:
        if a.IsInRing():
            a_name = "(00r%s)"%(str(a.GetSymbol()))
        else:
            a_name = "(00n%s)"%(str(a.GetSymbol()))
    else:
        if a.IsInRing():
            a_name = "(%sr%s)"%(str(a.GetFormalCharge()), str(a.GetSymbol()))
        else:
            a_name = "(%sn%s)"%(str(a.GetFormalCharge()), str(a.GetSymbol()))

    return a_name

def prepare_path_names(m, path):
    """
    path is a list of extracted paths
    """
    out = []
    for p in path:
        new_p = []
        for i in range(len(p)):
            if i==0:
                continue    # we don't wand to add the atom itself in the path chain
            elif i%2==0:
                new_p.append(prepare_atom_indicator(m, p[i]))
            else:
                new_p.append(p[i])

        s = ""
        for item in new_p:
            s+=str(item)
        out.append(s)
    
    return out

def rdkit_atom_from_labeled_smi(smi):
    mol = mol_with_hydrogens(smi)
    for a in mol.GetAtoms():
        if a.GetAtomMapNum()==1:
            return a

def shuffle_two_lists(a, b):
    c = list(zip(a, b))
    random.shuffle(c)
    a, b = zip(*c)
    return a, b

def label_each_atom(smi):
    mol = mol_with_hydrogens(smi)
    out = []
    idxs = []
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        idxs.append(atom.GetIdx())
    for i, atom in enumerate(mol.GetAtoms()):
        if i>0:
            a = mol.GetAtomWithIdx(idxs[i])
            a.SetAtomMapNum(1)
            a_prev = mol.GetAtomWithIdx(idxs[i-1])
            a_prev.SetAtomMapNum(0)
        else:
            a = mol.GetAtomWithIdx(idxs[i])
            a.SetAtomMapNum(1)
        out.append(Chem.MolToSmiles(mol))

    return out

def onek_encoding_unk(x, items):
    if x not in items:
        x = items[-1]
    return list(map(lambda s: x == s, items))

def atom_features(atom):
    attributes = onek_encoding_unk(atom.GetSymbol(), ELEMENTS) \
            + onek_encoding_unk(atom.GetDegree(), [0,1,2,3,4,5]) \
            + onek_encoding_unk(atom.GetExplicitValence(), [1,2,3,4,5,6]) \
            + [atom.GetIsAromatic()] \
            + [atom.GetNumRadicalElectrons()] \
            + [atom.GetIsAromatic() == False and any([neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors()])] \
            + [atom.IsInRing()] \
            + [atom.GetAtomicNum() in [9, 17, 35, 53, 85, 117]] \
            + [atom.GetAtomicNum() in [8, 16, 34, 52, 84, 116]] \
            + [atom.GetAtomicNum() in [7, 15, 33, 51, 83]] \
            + [atom.GetAtomicNum() in [3, 11, 19, 37, 55, 87]] \
            + [atom.GetAtomicNum() in [4, 12, 20, 38, 56, 88]] \
            + [atom.GetAtomicNum() in [13, 22, 24, 25, 26, 27, 28, 29, 30, 33, 42, 44, 45, 46, 47, 48, 49, 50, 78, 80, 82]]
    attributes = np.array(attributes, dtype=np.float32)
    attributes[np.isnan(attributes)] = 0.0 # filter nan
    attributes[np.isinf(attributes)] = 0.0 # filter inf
    return attributes

def path_list_to_count_dict(paths):
    return Counter(paths)   


def canonicalKekule(smi):
    """Simple method to canonicalize the kekule form"""
    mol = molBySmiles(smi);
    removeNonsenseStereo(mol)
    OEAssignAromaticFlags(mol, OEAroModelMMFF)
    for bond in mol.GetBonds():
        if bond.IsAromatic():
            bond.SetIntType(5)
    OECanonicalOrderAtoms(mol)
    OECanonicalOrderBonds(mol)
    OEClearAromaticFlags(mol)
    OEKekulize(mol)
    return createCanonicalAtomMapSmiString(mol)


