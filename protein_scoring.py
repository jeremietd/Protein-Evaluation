device = 'cuda:0'

from glob import glob
from pgen.utils import parse_fasta
import pandas as pd
import os
import argparse

from scoring_metrics import structure_metrics as st_metrics
from scoring_metrics import single_sequence_metrics as ss_metrics
from scoring_metrics import alignment_based_metrics as ab_metrics

#Reset calculated metrics (creates a new datastructure to store results, clearing any existing results)
results = dict()

# Default directories
file_dir = os.path.dirname(os.path.realpath(__file__))
default_pdb_dir = os.path.join(file_dir, "pdbs")
default_reference_dir = os.path.join(file_dir, "reference_seqs")
default_target_dir = os.path.join(file_dir, "target_seqs")

os.makedirs(default_pdb_dir) if not os.path.exists(default_pdb_dir) else None
os.makedirs(default_reference_dir) if not os.path.exists(default_reference_dir) else None
os.makedirs(default_target_dir) if not os.path.exists(default_target_dir) else None

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--pdb_dir", type=str, default=default_pdb_dir, help="Directory containing pdb files")
parser.add_argument("--reference_dir", type=str, default=default_reference_dir, help="Directory containing reference fasta files")
parser.add_argument("--target_dir", type=str, default=default_target_dir, help="Directory containing target fasta files")
parser.add_argument("--sub_matrix", type=str, choices=["blosum62", "pfasum15"], default="blosum62", help="Substitution matrix to use for alignment-based metrics")
parser.add_argument("--score_mean", type=bool, action="store_false", help="Whether to score the mean of the scores for mutated sequences")
parser.add_argument("--identity", type=bool, action="store_false", help="Whether to score the identity of the mutated sequence to the closest reference sequence")
parser.add_argument("--sub_gap_open", type=int, default=10, help="Gap open penalty for alignment-based metrics")
parser.add_argument("--sub_gap_extend", type=int, default=2, help="Gap extend penalty for alignment-based metrics")
args = parser.parse_args()

# Check that the required directories exist
pdb_dir = os.path.abspath(args.pdb_dir)
reference_dir = os.path.abspath(args.reference_dir)
target_dir = os.path.abspath(args.target_dir)

assert os.path.exists(pdb_dir), f"PDB directory {pdb_dir} does not exist"
assert os.path.exists(reference_dir), f"Reference directory {reference_dir} does not exist"
assert os.path.exists(target_dir), f"Target directory {target_dir} does not exist"

# Check that the required files exist
pdb_files = glob(pdb_dir + "/*.pdb")
reference_files = glob(reference_dir + "/*.fasta")
target_files = glob(target_dir + "/*.fasta")

assert len(pdb_files) > 0, f"No pdb files found in {pdb_dir}"
assert len(reference_files) > 0, f"No reference fasta files found in {reference_dir}"
assert len(target_files) > 0, f"No target fasta files found in {target_dir}"


# Structure metrics
# ESM-IF, ProteinMPNN, MIF-ST, AlphaFold2 pLDDT
st_metrics.ESM_IF(pdb_files, results)
st_metrics.ProteinMPNN(pdb_files, results)
st_metrics.MIF_ST(pdb_files, results, device)
st_metrics.AlphaFold2_pLDDT(pdb_files, results)

# Single sequence metrics
# ESM-1v, ESM-1v-mask6, CARP-640m-logp, Repeat-1, Repeat-2, Repeat-3, Repeat-4
target_seqs_file = "/tmp/target_seqs.fasta"
with open(target_seqs_file,"w") as fh:
  for target_fasta in target_files:
    for name, seq in zip(*parse_fasta(target_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

ss_metrics.CARP_640m_logp(target_seqs_file, results, device)
ss_metrics.ESM_1v(target_files, results, device)
ss_metrics.ESM_1v_mask6(target_files, results, device)
ss_metrics.Repeat(target_files, results)

# Alignment-based metrics
# ESM-MSA, Identity to closest reference, Subtitution matix (BLOSUM62 or PFASUM15) score mean of mutated positions
sub_matrix = args.sub_matrix.upper()
score_mean = args.sub_score_mean
identity = args.identity
sub_gap_open = args.sub_gap_open
sub_gap_extend = args.sub_gap_extend

#concatenate reference sequences
# Reference sequences
reference_seqs_file = "/tmp/reference_seqs.fasta"
with open(reference_seqs_file,"w") as fh:
  for reference_fasta in reference_files:
    for name, seq in zip(*parse_fasta(reference_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

# Target sequences
target_seqs_file = "/tmp/target_seqs.fasta"
with open(target_seqs_file,"w") as fh:
  for target_fasta in reference_files:
    for name, seq in zip(*parse_fasta(target_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

ab_metrics.ESM_MSA(target_seqs_file, reference_seqs_file, results)
ab_metrics.substitution_score(target_seqs_file, reference_seqs_file,
                              substitution_matrix=sub_matrix, 
                              Substitution_matrix_score_mean_of_mutated_positions=score_mean, 
                              Identity_to_closest_reference=identity,
                              results=results,
                              gap_open=sub_gap_open,
                              gap_extend=sub_gap_extend,)

# Download results
df = pd.DataFrame.from_dict(results, orient="index")
save_path = "calculated_metrics.csv"
df.to_csv(save_path)
print(f"Results saved to {save_path}")

