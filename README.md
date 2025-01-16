# DSA-learn
A dataset generation toolbox for data-driven dynamic security assessment: DSA-learn.jl.

For more information, read our paper 'A Dataset Generation Toolbox for Dynamic Security Assessment: On the Role of the Security Boundary'. If you like the paper or the toolbox, please cite our work:

@article{giraud2025datasets,
  author = {Bastien Giraud, Lola Charles, Agnes Marjorie Nakiganda, Johanna Vorwerk and Spyros Chatzivasileiadis},
  title = {A Dataset Generation Toolbox for Dynamic Security Assessment: On the Role of the Security Boundary},
  journal = {Submitted to IREP 2025 (under review)},
  year = {2025},
  url = {ArXiV_url},
}


# Repository

The main.jl file is the main file that needs to be run to generate your own dataset. You can specify your parameters in the init.jl file. This file contains all the parameters to generate datasets for the proposed method and for the benchmarks. To generate a dataset using one of the benchmarks, you have to run the code for the benchmarks located in the 'baseline' folder.

A quick run through of where everything is located in the repository:

- You can find the necessary Julia packages and their versions under functions/env.jl.
- The functions used throughout the whole toolbox can be found in the 'functions' folder.
- The case- and generator data can be found in the 'cases' folder.
- The benchmarks used to compare our proposed method can be found in the 'baseline' folder.
- Use the Manifest.toml and Project.toml files to import all the used packages and their dependencies.


