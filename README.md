# *DeepRho(v2.0)*
DeepRho: software accompanyment for "DeepRho: Accurate Estimation of Recombination Rate from Inferred Genealogies using Deep Learning", Haotian Zhang and Yufeng Wu, manuscript, 2021.

DeepRho constructs images from population genetic data and takes advantage of the power of convolutional neural network (CNN) in image classification to etstimate recombination rate. The key idea of DeepRho is generating genetics-informative images based on inferred gene geneaologies and linkage disequilibrium from population genetic data.

# Code
`deeprho` is an open-source software developed for per-base recombination rate estimation from inferred genealogies using deep learning. `deeprho` makes estimates based on LD patterns and local genealogical trees inferred by [*RENT+*](https://github.com/SajadMirzaei/RentPlus).

---
### Prerequisite
- OS: Linux, Windows, MacOS
- Software: [Conda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/download.html)
- Device: [CUDA-Enabled GPU](https://developer.nvidia.com/cuda-gpus) (optional, CPU-only settings see below)

### Installation
1. Clone from GitHub: `git clone https://github.com/haotianzh/deeprho_v2.git` or Download & unzip the file to your local directory.
2. Enter the directory: `cd deeprho_v2`
3. Create a virtual environment through conda: `conda create -n deeprho python=3.7`
4. activate conda environment: `conda activate deeprho`
5. install [msprime](https://tskit.dev/msprime/docs/stable/installation.html): `conda install msprime -c conda-forge`
6. build deeprho: `pip install .`
7. inspect correct installation: `deeprho -v`


### Input formats
- ms-formatted input (the first line is position (seperated by space) followed by haplotype sequences, check `examples/data.ms` for details)
- VCF file (check `examples/data.vcf`)

### Usage (Examples)
- #### Estimate
    ```jupyterpython
    # given a .vcf file, use default models and plot the estimated recombination map. 
    deeprho estiamte --file examples/data.vcf --length 1e5 --ne 1e5 --ploidy 2 --m1 models/model_fine.hdf5 --m2 models/model_large.hdf5 --plot
  ```
- #### simuate
    ```jupyterpython
    # sample 10 genomes with recombination rates uniformly spaced between 1e-9 and 1e-8, and use 8 cpus for parallel speeding.
    deeprho simulate --nsam 10 --npop 50 --mutation-rate 1e-8 --ne 1e5 --rmin 1e-9 --rmax 1e-8 --num-thread 8 
    ```

### Docs
- #### Estimate
  ```python
    deeprho estimate [-h] [--file FILE] [--length LENGTH] [--ne NE] [--ploidy PLOIDY] [--res RES] \
                      [--threshold THRESHOLD] [--gws GWS] [--ws WS] [--ss SS] [--m1 MODEL_FINE] \
                      [--m2 MODEL_LARGE] [--plot] [--savenp] [--num-thread NUM_THREAD]  
  ```
    | Arguments                   | Descriptions                                                      |
    |------------------------------|-------------------------------------------------------------------|
    | `--file <FILE>`             | Input file                                                        |
    | `--length <LENGTH>`         | Length of chromosome                                              |
    | `--m1 <MODEL_FINE>`         | Path of fine model                                                |
    | `--m2 <MODEL_LARGE>`        | Path of large model                                               |
    | `--ploidy <PLOIDY>`         | Ploidy (default 1)                                                |
    | `--ne <NE>`                 | Effective population size (default 10<sup>5</sup>)                |
    | `--res <RES>`               | Resolution of estimation (default 10<sup>4</sup>)                 |
    | `--threshold <THRESHOLD>`   | Threshold of recombination Hotspot (default 5x10<sup>-8</sup>)    |
    | `--gws <GWS>`               | Window size for inferring genealogy (default 10<sup>3</sup> SNPs) |
    | `--ws <WS>`                 | Window size for performing `deeprho` (default 50 SNPs)            |
    | `--ss <SS>`                 | Step size for performing `deeprho` (default same as `--ws`)       |
    | `--savenp`                  | Save estimated rates as numpy ndarray (saved as `<FILE>.out.npy`) |
    | `--plot`                    | Plot recombination map (saved as `<FILE>.out.png`)                |
    | `--num_thread <NUM_THREAD>` | Specify number of workers for parallel (default 4)                |
    | `--help, -h`                | Show usage                                                        |

    - `<LENGTH>` can be either explicitly specified or inferred from input, if the latter, `<LENGTH>`= S<sub>n</sub>-S<sub>1</sub>, 
       where S<sub>n</sub> is physical position of the last SNP site, S<sub>1</sub> is the position of the first SNP site. 
    - `<MODEL_FINE>, <MODEL_LARGE>` are two pretrained-models, `deeprho` takes two-stages strategies to estimate recombination rate, 
       `<MODEL_FINE>` is applied for estimating recombination background regions while `<MODEL_LARGE>` is used to fine-tune hotspot regions.
        two default models are included in this repo, users are also allowed to train their own models through following sections.
    - `<THRESHOLD>` define a threshold above which a region can be regarded as a hotspot. 5x10<sup>-8</sup> is set as default.
    - `<GWS>` guides how large region the genealogies are inferred from. As our test, 1000 is a great choice to include as much information as possible
       for improving local genealogical inference.
    - `<WS>` specified window size during estimation, `deeprho` is a sliding-window based algorithm, too large window will weaken its ability to detect
       hotspot, however, too small region may undermine its performance. 50 is highly recommended in practice.
  

- #### Simulate
    ```python
        deeprho simulate [-h] [--nsam NSAM] [--npop NPOP] [--ne NE] [--ploidy PLOIDY] \
                          [--mutation-rate MUTATION_RATE] [--demography DEMOGRAPHY] \
                          [--rmin RMIN] [--rmax RMAX] [--num-thread NUM_THREAD] [--out OUT]
    ```
    | Arguments                         | Descriptions                                                                          |
    |---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
    | `--nsam <NSAM>`                   | Number of different recombination rates sampled from \<RMIN> to \<RMAX> (default 200) |
    | `--ndraw <NDRAW>`                 | Number of random draws at each different recombination rate (default 5)               |
    | `--npop <NPOP>`                   | Number of individuals (default 100)                                                   |
    | `--mutation-rate <MUTATION_RATE>` | Mutation rate (default 2.5x10<sup>-8</sup>)                                           |
    | `--ploidy <PLOIDY>`               | Ploidy (default 1)                                                                    |
    | `--ne <NE>`                       | Effective population size (default 10<sup>5</sup>)                                    |       
    | `--demography <DEMOGRAPHY_FILE>`      | Demographic history (see details)                                                     |
    | `--num_thread <NUM_THREAD>`       | Specify number of workers for parallel (default 4)                                    |
    | `--rmin <RMIN>`                   | Minimum recombination rate (default 10<sup>-9</sup>)                                  |
    | `--rmax <RMAX>`                   | Maximum recombination rate (default 5x10<sup>-7</sup>)                                |
    | `--help, -h`                      | Show usage                                                                            |

    - `<NSAM>` means how many **uniform interpolation points** we are expecting to have between `<RMIN>` to `<RMAX>`, suppose `<RMIN>=10^-9, <RMAX>=10^-8, <NSAM>=3`, 
       then three genomes with recombination rates as [10<sup>-9</sup>, 5.5x10<sup>-9</sup>, 10<sup>-8</sup>] will be simulated respectively.
       Usually, larger `<NSAM>` results in more accurate estimation, but more time-consuming meanwhile.
    - `<NDRAW>` means how many windows we randomly draw from a single genome simulated as describe above, larger `<NDRAW>` results in more data points sampled at the same recombination rate.
    - `<NPOP>` should be considered together with `<PLOIDY>`, there will be `<NPOP>*<PLOIDY>` haplotypes in a single genome. For example, `<NPOP>=50, <PLOIDY>=2` means 100 haplotypes in total.
    - `<NE>, <DEMOGRAPHY_FILE>` these are important parameters to define effective population size changes, set either `<NE>` or `<DEMOGRAPHY_FILE>`. Simply setting `<NE>=10^5` refers to a constant population size, 
       however, if you have confidence with a more detailed demographic history, a `.csv` file can also be provided instead, see example. If both `<NE>, <DEMOGRAPHY_FILE>` are specified, `deeprho` keeps `<DEMOGRAPHY_FILE>` while `<NE>` is ignored.
  

