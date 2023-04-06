from biotite.structure import get_residue_count
from biotite.structure.residues import get_residues
from biotite.structure.io import pdb
from biotite.sequence import ProteinSequence

def add_metric(metrics_dict, protein_name, metric_name, value):
  if protein_name not in metrics_dict:
    metrics_dict[protein_name] = dict()
  metrics_dict[protein_name][metric_name] = value

def get_pdb_sequence(pdb_path):
    with open(pdb_path) as f:
        pdb_file = pdb.PDBFile.read(pdb_path)
        atoms  = pdb_file.get_structure()
        residues = get_residues(atoms)[1]
    return ''.join([ProteinSequence.convert_letter_3to1(r) for r in residues])

def residues_in_pdb(pdb_path):
    with open(pdb_path) as f:
        pdb_file = pdb.PDBFile.read(pdb_path)
        atoms  = pdb_file.get_structure()
    return get_residue_count(atoms)