# DSA-learn
This repository provides tools to generate datasets for dynamic security assessment using our proposed method and baseline benchmarks. Below is a guide to help you navigate the repository and generate your own dataset.

For more information, read our paper 'A Dataset Generation Toolbox for Dynamic Security Assessment: On the Role of the Security Boundary'. If you use this work, please cite it as follows:

```bibtex
@article{giraud2025datasets,
  author = {Bastien Giraud, Lola Charles, Agnes Marjorie Nakiganda, Johanna Vorwerk and Spyros Chatzivasileiadis},
  title = {A Dataset Generation Toolbox for Dynamic Security Assessment: On the Role of the Security Boundary},
  journal = {Submitted to IREP 2025 (under review)},
  year = {2025},
  url = {(http://arxiv.org/abs/2501.09513)},
}
```

## How to Use

1. **Run the Main Script**:  
   The main script is `main.jl`. Execute this file to generate your dataset.

2. **Configure Parameters**:  
   All parameters required for dataset generation are specified in `init.jl`. This includes configurations for both the proposed method and the benchmark methods.

3. **Generate Benchmark Datasets**:  
   To generate datasets using one of the baseline benchmarks, execute the corresponding code in the `baseline` folder.

## Repository Structure

Here’s an overview of where to find everything in this repository:

- **Julia Packages**:  
  The required Julia packages and their specific versions are listed in `functions/env.jl`.

- **Toolbox Functions**:  
  All functions used throughout the toolbox are located in the `functions` folder.

- **Case and Generator Data**:  
  Case data and generator parameters are stored in the `cases` folder.

- **Baseline Benchmarks**:  
  The benchmarks used for comparison with the proposed method are found in the `baseline` folder.

- **Dependencies and Environments**:  
  Use the `Manifest.toml` and `Project.toml` files to import all necessary packages and their dependencies.


