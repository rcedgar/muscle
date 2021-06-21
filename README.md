# muscle
Multiple sequence alignment with MSA ensemble bootstrap

## Downloads

Binary files are self-contained, no dependencies.

[muscle5.0.1278_linux64](https://github.com/rcedgar/muscle/raw/main/binaries/muscle5.0.1278_linux64)   
[muscle5.0.1278_win64.exe](https://github.com/rcedgar/muscle/raw/main/binaries/muscle5.0.1278_win64.exe)   

## Usage

    Amino acid sequences, runs MPC or Super5 algorithm:
      muscle -mpc seqs.fa -out aln.fa
      muscle -super5 seqs.fa -out aln.fa
    
    Nucleotide or amino acid sequences, runs muscle v3 algorithm:
      muscle -in seqs.fa -out aln.fa [-maxiters 2]
    
    Input is FASTA, output is aligned FASTA.
    
    For amino acid sequences, use MPC for up to around 500 sequences, for larger
    sets, use Super5. The muscle v3 algorithm is less accurate, not recommended.
    
    For nucleotide sequences, use -maxiters 2 if the default is too slow, which
    typically happens with a few hundred sequences.
    
    mpc and super5:
      -perturb n           Integer seed for randomized HMM perturbations (default 0).
                             To make an ensemble, run several times with different seeds.
                             Zero (default) disables perturbing HMM parameters.
      -perm type           Guide tree permutation none, abc, acb, or bca (default none).
      -hmmin filename      Read HMM parameters from file.
      -hmmout filename     Save HMM parameters to file.
      -consiters n         Nr consistency iterations (default 2).
      -refineiters n       Nr refinement iterations (default 100).
      -threads n           Use n threads. Default is min(number of CPU cores, 20).
    
    super5:
      -maxcoarse n         Max size for coarse clusters (default 500).
      -paircount n         Max pairs for profile alignment (default 2000).
    
    muscle v3:
      -maxiters n          Max iterations. 2 is progressive, >2 does refinement.
                           Default is continue until convergence.
    
    All algorithms:
      -quiet               Turn off progress messages.
      -log filename        Write time, memory use and progress messages to file.
