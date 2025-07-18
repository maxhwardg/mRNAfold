# mRNAfold

## Table of Contents
- [mRNAfold](#mrnafold)
  - [Table of Contents](#table-of-contents)
  - [Basic Information](#basic-information)
  - [Citation](#citation)
  - [Build Guide](#build-guide)
  - [Example Config Files](#example-config-files)
  - [Python Script](#python-script)

## Basic Information
This repository contains the "mRNAfold" software package. This package comprises code for mRNA folding algorithms written by Dr Max Ward. These algorithms can be used to design optimised mRNA sequences quickly and easily. The algorithms are written in C++ and a Python runner is included.

The phrase "mRNA folding" refers to a specific algorithmic strategy that extends RNA folding algorithms to designing optimised mRNA sequences. Similar algorithms include [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign), [CDSfold](https://github.com/gterai/CDSfold), and [DERNA](https://github.com/elkebir-group/derna).

This software package has several advantages over existing packages. It supports the following notable features:

- Suboptimal folding allowing users to get multiple mRNA sequences that are near optimal. At time of writing no other software has this feature. Novel algorithms were developed to enable this, which are soon to be published
- Techniques to that increase the diversity of suboptimal samples. At time of writing no other software has this feature. Novel algorithms were developed to enable this, which are soon to be published
- Supports MFE and CAI optimisation (like LinearDesign and DERNA)
- Multithreaded algorithms (utilizing the C++ standard library)
- Very efficient code. mRNAfold is competitive with LinearDesign despite mRNAfold employing an exact algorithm and LinearDesign using an approximation. In addition, with multithreading enabled mRNAfold is extremely fast. On a 16-core AMD 7950X with a protein comprising 1500 amino acids LinearDesign took 302s and mRNAfold took 46s (DERNA and CDSfold are much slower than both).
- Can incoporate the UTRs into the optimisation
- Can force or ban codons at specific positions in the sequence

## Citation

For now please cite as in the following BibTeX entry.

```
@software{mrnafold2024,
  author = {Ward, Max},
  year = {2024},
  month = {12},
  title = {{mRNAfold Software Package}},
  url = {https://github.com/maxhwardg/mRNAfold},
  version = {0.0.0},
}
```

When a paper is published it will be added here.

  

## Build Guide
Several steps must be followed exactly to build the program. These steps should work on up-to-date Linux images. These steps have been tested to work on a clean Ubuntu 24.04 install.

Start by using update/upgrade the ensure our system is up-to-date.
`sudo apt update`
`sudo apt upgrade`

Make sure the GNU C++ compiler is installed: `sudo apt install g++`

Alternatively Clang can be used: `sudo apt install clang`

Next, let's check the compiler version.

`g++ --version` or `clang++ --version`

The code was developed and tested using clang 18.1.3. The code has also been tested with g++ 13.3.0.

Since the code is written using modern C++20, we need up-to-date compiler and toolchain support.  Older version of Ubuntu and other Linux distros may require manually installing a newer toolchain. You may need to configure CMake to use a specific compiler if your system default is insufficient.

We also use several libraries and tools including CMake and Intel TBB. We will install these next.
`sudo apt install cmake libtbb-dev`

Now we are ready to build the program. We use out of source builds with CMake. Let's start by making the build folder. First, make sure we are the root directory of this project.

Run: `ls`

This should show: `configs  cpp_src  LICENSE  pysrc  README.md`

Next, create the build directory: `mkdir build`

Next, configure build/ directory in release mode with cmake: `cmake -S cpp_src/ -B build/ -DCMAKE_BUILD_TYPE=Release`

This should produce something like:
```
-- The CXX compiler identification is GNU 13.3.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Configuring done (0.3s)
-- Generating done (0.0s)
-- Build files have been written to: /pathto/mRNAfold/build
```

Next, we must go to the build directory and build it. First, let's go to there: `cd build`

Next, we use `make` to build the code: `make all`

This may produce some output. This can be ignored so long as it does not error or crash. No errors should occur and the process should complete with the following final line: `[100%] Built target fold_codon_graph`

Now we can run the fold_codons_graph program. It is best to do this from the project directory. Let's get back there: `cd ../`. It is called the `fold_codon_graph` program since it does mRNA folding using a "codon graph", which is a technique related to those used by LinearDesign and CDSfold.

Now, we can run fold_codons: `./build/exe/fold_codon_graph`

This will print out usage instructions:

```
Usage: ./build/exe/fold_codon_graph <config_file_path>
```

The way fold_codons works is to take a single argument, which is the path to a configuration file.
This file sets all the parameters for the program to run. Some example config files are given in the
`configs` directory. Running `ls configs` should show something like the following:

```
egfp_example.config  epo_example.config  ffluc_example.config  ...
```

Looking at these files should provide details about how to make your own config files.
In short, a config must provide a number of arguments each on its own line. Each line contains the name of the argument, then the value separated by a space.
These lines can be in any order, and all of them are optional except `aa_seq`. If they are not given, a default value is used.

We can run an example with `./build/exe/fold_codon_graph ./configs/minigfp_example2.config`

This should produce many suboptimal folds that look like the following.

```
AUGGAAAAGUCCUUUGUGAUUACUGACCCAUGGCUGCCCGACUAUCCUAUCAUCAGCGCCAGCGACGGCUUUCUGGAGCUGACUGAGUAUAGUCGGGAGGAGAUCAUGGGCAGGAACGCCCGGUUUCUGCAGGGGCCUGAGACCGAUCAGGCCACCGUGCAGAAGAUCCGGGAUGCUAUACGGGAUCGGAGGCCCACUACUGUGCAGCUGAUCAACUACACCAAGAGCGGGAAGAAGUUCUGGAAUUUGUUGCAUCUGCAGCCCGUGUUUGACGGAAAAGGCGGGCUGCAGUAUUUCAUCGGCGUGCAGCUGGUGGGCUCUGAUCACGUG
(.((.(...).))...)......((.(((((((((.(((((((((....((((((((.(((((....)))....)).))))).)))..))))))))).))...)))))))))((....(((((.((((((((((((((((......))))))).)).)))))))...)))))...)).(((((.((((((.((((((((((((((.((((..((((....((.(((((........))))))).....))))...(((((((((((.(((......))).)))))))))))........)))))))))).))))))))))))))).))))
Score: -113.5977
EFE: -161.1318
CAI: 0.9048
AUGGAAAAAUCCUUCGUGAUUACUGACCCAUGGCUGCCCGACUAUCCUAUCAUCAGCGCCAGCGAUGGCUUCCUGGAGCUGACUGAGUAUAGUCGGGAGGAGAUCAUGGGCAGGAACGCCCGGUUUCUGCAGGGGCCUGAGACCGAUCAGGCCACCGUGCAGAAGAUCCGGGAUGCUAUACGGGAUCGGAGGCCCACUACUGUGCAGCUGAUCAACUACACCAAGAGCGGGAAGAAGUUCUGGAAUUUGUUGCAUCUGCAGCCCGUGUUUGACGGAAAAGGCGGGCUGCAGUAUUUCAUCGGCGUGCAGCUGGUGGGCUCUGAUCACGUG
(((((.......)))))......((.(((((((((.(((((((((....((((((((.((((.((....)).)))).))))).)))..))))))))).))...)))))))))((....(((((.((((((((((((((((......))))))).)).)))))))...)))))...)).(((((.((((((.((((((((((((((.((((..((((....((.(((((........))))))).....))))...(((((((((((.(((......))).)))))))))))........)))))))))).))))))))))))))).))))
Score: -122.9577
EFE: -163.9839
CAI: 0.9038
```

## Example Config Files

There are several useful features demonstrated in the example config files. Let us begin with a description of some of the most useful parameters.

- `aa_seq MVSKGEELFTG`: Amino acid sequence to fold. E.g., "MVSKGEELFTG". This is a required parameter
- `parallel false`: Set to `true` if a parallel algorithm should be used and `false` otherwise. Defaults to `false` if nothing is provided.
- `lambda 1.0`: The weight given to Codon Adaptation Index (CAI). The program will produce solutions in order of their score, which is `MFE - ln(cai)*lambda`. A lambda of 1.0 weights CAI and MFE roughly equally, while 0.0 ignores CAI and lambda=infinity ignores MFE. Note that MFE means Minimum Free Energy, and is a measure of structural stability. This defaults to 1.0 if nothing is provided.
- `subopt_randomness 0.0`: Amount of noise in the suboptimal sampling. Useful for increasing the diversity of the samples. This is a percentage, and should be between 0 and 1 where 0.1 would mean that 10% of the score of a traceback is contributed by noise. The score with this parameter included is MFE - ln(cai)*cai_lambda + noise. This defaults to 0.0 if not provided.
- `num_subopt_traces 100`: The number of results to generate. If not provided defaults to 100.
- `num_subopt_rounds 1`: The number of rounds of suboptimal sampling to run. The total number of samples will be num_subopt_traces * num_subopt_rounds. Note that rounds get executed in parallel if parallel is set to True. Each round will get a different random seed derived from the root seed value. This is useful for getting a more diverse set of samples. It is also useful for making use of multiple cores during sampling. It is recommended to use this in combination with either subopt_randomness otherwise each round will generate the same results. If not provided defaults to 1.
- `rng_seed 0`: The root seed value to use for randomness (e.g., in subopt_randomness). Defaults to 0 if not provided.

More details of the parameters can be found by looking at the `FoldCodonsConfig` documentation in `pysrc/fold_codons.py`. In addition, `configs/` contains examples.

The following files are good examples of the basic parmeters: `minigfp_example1.config`, `minigfp_example2.config`, and `egfp_example.config`. In addition, `ffluc_example.config` and `lacz_example.config` are similarly good examples using longer sequences, but are much slower to run.


See `minigfp_example3.config`. It demonstrates the usage of utr_5p, utr_3p, and
banned_motifs, which can be used to add UTRs and banned sequences respectively. The UTRs will not be modified by the algorithm but will be considered when calculating MFE.

The `minigfp_example4.config` shows an example of the encourage_unpaired parameter. This encourages certain regions to be unpaired by using a penalty. Note that this setting is subtle. It will find samples where the sequence-structure pair does not have pairs in the given region, but that does not enforce that the MFE structure for that sequence does not have pairs there. Nonetheless, it is a good way to encourage certain regions to have fewer pairs. It works well in combination with suboptimal folding to find sequences with unpaired regions.

See `ffluc_example2.config` for keep_chance, which is a good way of increasing sampling diversity. This is an alternative method to subopt_randomness. It is conceptually the same as taking all possible samples and deleting them at random until the remaining proportion is keep_chance, then sampling the top num_subopt_traces from there.

One can use num_subopt_rounds to generate suboptimal samples in parallel. This is useful, since the folding algorithm is already parallel, but sampling is not. This can be achieved by running multiple sampling rounds in different threads. This is designed to work with subopt_randomness (or keep_chance), and each thread will use a different random seed and therefore sample differently. See `egfp_example.config` for an example usage.

See `egfp_example2.config`. It shows how to use a custom codon frequency table. Note that it assumes the user is running the program from the project root. E.g., `./build/exe/fold_codon_graph configs/egfp_example2.config`. By default, a human codon table is used, which is baked into the executable (no need to pass a file path).

See `epo_example.config`. It shows how to force/ban codons at specific positions. Note that these positions are zero-indexed.

## Python Script

Once you have build the binaries, there is a python script that works as a convient interface. It calls the compiled binaries and assumes there are are stored in the folder given in the above instructions. This script contains documentation for each of the possible arguments and is useful to read for crafting the .config files mentioned above.

A example python program that gives a python function interface to the binary is provided in `pysrc`.
We can run it by going to the pysrc directory `cd pysrc` then executing the script `python3 fold_codons.py`. The main function in the script shows an example use. The function and its parameters are documented. It can be used in other scripts the usual way via `import fold_codons`.

Be aware that the `fold_codons.py` script has a dependency. It uses the ViennaRNA package. You can install this using pip via `pip install ViennaRNA`. This is used only to show to evaluate the results against ViennaRNA as a reference.
