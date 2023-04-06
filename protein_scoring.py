device = 'cuda:0'

from glob import glob
from pgen.utils import parse_fasta
import pandas as pd

#Reset calculated metrics (creates a new datastructure to store results, clearing any existing results)
results = dict()

# Structure metrics
# ESM-IF, ProteinMPNN, MIF-ST, AlphaFold2 pLDDT
from scoring_metrics import structure_metrics as st_metrics

st_metrics.ESM_IF(results)
st_metrics.ProteinMPNN(results)
st_metrics.MIF_ST(results, device)
st_metrics.AlphaFold2_pLDDT(results)

# Single sequence metrics
# ESM-1v, ESM-1v-mask6, CARP-640m-logp, Repeat-1, Repeat-2, Repeat-3, Repeat-4

target_seqs_file = "/tmp/target_seqs.fasta"
with open(target_seqs_file,"w") as fh:
  for target_fasta in glob("/target_seqs/*"):
    for name, seq in zip(*parse_fasta(target_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

from scoring_metrics import single_sequence_metrics as ss_metrics

ss_metrics.CARP_640m_logp(target_seqs_file, results, device)
ss_metrics.ESM_1v(results, device)
ss_metrics.ESM_1v_mask6(results, device)
ss_metrics.Repeat(results)

# Alignment-based metrics
# ESM-MSA, Identity to closest reference, Subtitution matix (BLOSUM62 or PFASUM15) score mean of mutated positions

# substitution_matrix = "BLOSUM62" # ["BLOSUM62", "PFASUM15"]

#concatenate reference sequences
reference_seqs_file = "/tmp/reference_seqs.fasta"
with open(reference_seqs_file,"w") as fh:
  for reference_fasta in glob("/reference_seqs/*"):
    for name, seq in zip(*parse_fasta(reference_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

target_seqs_file = "/tmp/target_seqs.fasta"
with open(target_seqs_file,"w") as fh:
  for target_fasta in glob("/target_seqs/*"):
    for name, seq in zip(*parse_fasta(target_fasta, return_names=True, clean="unalign")):
      print(f">{name}\n{seq}", file=fh)

from scoring_metrics import alignment_based_metrics as ab_metrics
ab_metrics.ESM_MSA(results)
ab_metrics.substitution_score(target_seqs_file, reference_seqs_file,
                              substitution_matrix="BLOSUM62", 
                              Substitution_matrix_score_mean_of_mutated_positions=True, 
                              Identity_to_closest_reference=True,
                              results=results,
                              gap_open=10,
                              gap_extend=2,)

# Download results
df = pd.DataFrame.from_dict(results, orient="index")
df.to_csv(f"/tmp/calculated_metrics.csv")

