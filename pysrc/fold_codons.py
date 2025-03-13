"""Code for interfacing with executables from Python"""
from typing import List
from dataclasses import dataclass, field
import subprocess
import tempfile
import RNA  # ViennaRNA Python bindings


@dataclass
class FoldCodonResult:
    """Datastructure to store results from calling the executable"""
    # RNA sequence
    primary: str
    # Dot bracket notation
    secondary:  str
    # The score is MFE - ln(cai)*cai_lambda
    # This is the objective function during mRNA folding
    # Suboptimal samples are ranked by this score
    score: float
    # Codon adaptation index
    cai: float
    # Ensemble free energy
    efe: float


@dataclass
class BannedCodon:
    """Defines a banned codon at a specific position"""
    # Codon, e.g., "AUC" or "GGU"
    # No check is performed to ensure that the codon is valid
    # May crash if the codon is not valid
    codon: str
    # Zero-indexed position of the codon in the aa sequence
    pos: int


@dataclass
class ForcedCodon:
    """Defines a forced codon at a specific position"""
    # Codon, e.g., "AUC" or "GGU"
    # No check is performed to ensure that the codon is valid
    # May crash if the codon is not valid
    codon: str
    # Zero-indexed position of the codon in the aa sequence
    pos: int


@dataclass
class FoldCodonsConfig:
    """Stores configuration for calling the executable"""
    # Amino acid sequence to fold. E.g., "MVSKGEELFTG"
    aa_seq: str
    # Path to the fold_codons executable
    exe_path: str = "../build/exe/fold_codon_graph"
    # True if multi-core parallelization is enabled
    # Note that this will only work if the executable was compiled
    # with the Threaded Building Blocks (TBB) library
    parallel: bool = False
    # List of codons that must be used at the specified positions
    forced_codons: List[ForcedCodon] = field(default_factory=list)
    # List of codons that are banned at the specified positions
    banned_codons: List[BannedCodon] = field(default_factory=list)
    # Weighting factor of codon adaptation index
    # The score is MFE - ln(cai)*cai_lambda
    cai_lambda: float = 1.0
    # Amount of noise in the suboptimal sampling
    # Useful for increasing the diversity of the samples
    # This is a percentage, and should be between 0 and 1
    # 0.1 would mean that 10% of the score of a traceback is contributed by noise
    # The score is MFE - ln(cai)*cai_lambda + noise
    # It is recommended to use either this or keep_chance, but not both
    # It is recommended to use subopt_randomness instead of keep_chance
    # This is because subopt_randomness is easier to interpret
    subopt_randomness: float = 0.0
    # Number of suboptimal traces to sample
    # Setting this to a lower number will make the algorithm faster
    # It will also use less memory
    # This is because the algorithm only keeps the top-k partial traces
    # in memory during traceback generation
    num_subopt_traces: int = 100
    # Probability that a traceback is kept during suboptimal sampling
    # Useful for increasing the diversity of the samples
    # Suboptimal sampling will be done in a way equivalent to the following process:
    # All tracebacks are generated and put into a list
    # Tracebacks in the list are randomly deleted such that the chance of keeping
    # a traceback is keep_chance
    # Note that a much more efficient algorithm is used, but the result is the same
    # Increasing this is a good way of increasing the diversity of the samples
    # It is recommended to use either this or subopt_randomness, but not both
    # Should be between 0 and 1
    keep_chance: float = 1.0
    # The 5' UTR. This sequence is prepended to the codon sequence
    # The UTRs are used in the optimization but are never modified
    # Should be all caps e.g., AUUACUGU
    utr_5p: str = ""
    # The 3' UTR. This sequence is appended to the codon sequence
    # The UTRs are used in the optimization but are never modified
    # Should be all caps e.g., AUUACUGU
    utr_3p: str = ""
    # A span of nucleotides (zero-index) that are encouraged to be unapaired
    # This is done by not counting the contribution of any pairs using a nt in the span
    # This is not a hard constraint, but a gentle encouragement
    # It is recommended to combine this with suboptimal sampling
    # to ensure a good sequence is found
    # Also, increasing sample diversity is recommended
    # (e.g., via subopt_randomness or keep_chance, and num_subopt_rounds)
    encourage_unpaired: List[tuple[int, int]] = field(default_factory=list)
    # A list of pairs of nucleotides that are discouraged from being paired
    # This is done similarly to encourage_unpaired. See above for details
    discouraged_pairs: List[tuple[int, int]] = field(default_factory=list)
    # A list of pairs of nucleotides that are encouraged to be paired
    # Nucleotides are zero-indexed
    # This is done by adding a bonus to the score for each pair
    # This is not a hard constraint, but a gentle encouragement
    encouraged_pairs: List[tuple[int, int]] = field(default_factory=list)
    # The number of rounds of suboptimal sampling to run
    # The total number of samples will be num_subopt_traces * num_subopt_rounds
    # Note that rounds get executed in parallel if parallel is set to True
    # Each round will get a different random seed
    # This is useful for getting a more diverse set of samples
    # It is also useful for making use of multiple cores during sampling
    # It is recommended to use this in combination with either subopt_randomness or keep_chance
    # to increase the diversity of the samples
    num_subopt_rounds: int = 1
    # Seed to use for the random number generator
    # This is important for keep_chance and subopt_randomness
    rng_seed: int = 0
    # The organism to use for the codon frequency table
    # Currently supports "human" and "mouse"
    # Uses a built-in table from the Kazusa Codon Usage Database
    # See "cpp_src/mwmrna_lib/src/codon_table_strings.hpp"
    organism: str | None = None
    # Path to codon frequency table
    # Format is assumed to be from Kazusa Codon Usage Database
    # Specifically, the "style like CodonFrequency output in GCG Wisconsin PackageTM"
    # Either organism or codon_freq_path should be set, but not both
    # See codon_tables/ for examples
    codon_freq_path: str | None = None


class FoldCodonsError(Exception):
    """Base class for exceptions in this module"""

# Given the configuration, call the fold_codons executable
# Returns a list containing the suboptimal samples


def call_fold_codons(cfg: FoldCodonsConfig) -> List[FoldCodonResult]:
    """
    Calls the fold_codons executable with the given configuration
    and returns the results as a list of suboptimal samples
    """
    tmp_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    tmp_file.write(f"aa_seq {cfg.aa_seq}\n")
    tmp_file.write(f"parallel {str(cfg.parallel).lower()}\n")
    tmp_file.write(f"lambda {cfg.cai_lambda}\n")
    tmp_file.write(f"subopt_randomness {cfg.subopt_randomness}\n")
    tmp_file.write(f"num_subopt_traces {cfg.num_subopt_traces}\n")
    tmp_file.write(f"num_subopt_rounds {cfg.num_subopt_rounds}\n")
    tmp_file.write(f"rng_seed {cfg.rng_seed}\n")
    if cfg.organism is not None:
        tmp_file.write(f"organism {cfg.organism}\n")
    if cfg.codon_freq_path is not None:
        tmp_file.write(f"codon_freq_path {cfg.codon_freq_path}\n")
    for fc in cfg.forced_codons:
        tmp_file.write(f"forced_codon {fc.pos} {fc.codon}\n")
    for bc in cfg.banned_codons:
        tmp_file.write(f"banned_codon {bc.pos} {bc.codon}\n")
    tmp_file.write(f"keep_chance {cfg.keep_chance}\n")
    if cfg.utr_5p != "":
        tmp_file.write(f"utr_5p {cfg.utr_5p}\n")
    if cfg.utr_3p != "":
        tmp_file.write(f"utr_3p {cfg.utr_3p}\n")
    for span in cfg.encourage_unpaired:
        tmp_file.write(f"encourage_unpaired {span[0]},{span[1]}\n")
    for pair in cfg.discouraged_pairs:
        tmp_file.write(f"discourage_pair {pair[0]},{pair[1]}\n")
    for pair in cfg.encouraged_pairs:
        tmp_file.write(f"encourage_pair {pair[0]},{pair[1]}\n")
    tmp_file.close()
    cmd = [cfg.exe_path, tmp_file.name]
    results = []
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.terminate()
    if p.returncode != 0:
        stderr_str = stderr.decode() if stderr is not None else "None"
        raise FoldCodonsError(f"fold_codons failed with return code: "
                              f"{p.returncode} and stderr: {stderr_str}")
    lines = stdout.decode().split('\n')
    if lines[-1] == "":
        lines = lines[:-1]
    for i in range(0, len(lines), 5):
        score = float(lines[i+2].split(": ")[1])
        efe = float(lines[i+3].split(": ")[1])
        cai = float(lines[i+4].split(": ")[1])
        results.append(FoldCodonResult(lines[i], lines[i+1], score, cai, efe))
    return results


def main():
    """Demonstrates how to use call_fold_codons"""

    def get_fc(rna: str, dangles: int = 0):
        # Default to dangles=0 since that is what mRNAfold uses
        md = RNA.md()
        md.dangles = dangles
        fc = RNA.fold_compound(rna, md)
        return fc

    config = FoldCodonsConfig(
        # miniGFP example
        aa_seq="MEKSFVITDPWLPDYPIISASDGFLELTEYSREEIMGRNARFLQGPETDQATVQKIRDAIRDRRPTTVQLINYTKSGKKFWNLLHLQPVFDGKGGLQYFIGVQLVGSDHV",
        num_subopt_traces=20,
        num_subopt_rounds=4,
        parallel=True,
        cai_lambda=1.0,
        utr_3p="GGGAGAACCCA",
        utr_5p="GGAGCGCAUUACGAUGCACUGAUGCUGAGCGCAUUACGAUGCUGAUGCUG",
        subopt_randomness=0.05,
        rng_seed=1234,
    )

    res = call_fold_codons(config)
    for r in res:
        print(f"Sample Primary:   {r.primary}")
        print(f"Sample Secondary: {r.secondary}")
        print(f"Sample Score: {r.score}")
        print(f"Sample CAI: {r.cai}")
        print(f"Sample EFE: {r.efe}")
        fc = get_fc(r.primary)
        print(f"Sample Vienna EFE: {fc.pf()[1]}")
        db, mfe = fc.mfe()
        print(f"Sample Vienna MFE: {mfe}")
        print(f"Sample Vienna Secondary: {db}")
        print("---")


if __name__ == '__main__':
    main()
