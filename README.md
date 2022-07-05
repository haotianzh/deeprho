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

### Data formats
- ms-formatted input (the first line is position (seperated by space) followed by haplotype sequences, check `examples/data.ms` for details)
- VCF file (check `examples/data.vcf`)

### Usage
- Estimate
  - 