bioscrape - Biological Stochastic Simulation of Single Cell Reactions and Parameter Estimation
==============================================================================================

|Build Status|

Bioscrape is a Chemical Reaction Network Simulator written in cython for speed and python compatability. It can be used for deterministic, stochastic, or single cell simulation and can also be used for sensitivity analysis. Bioscrape has built-in plug-and-play style parameter inference features. Refer to the Wiki page of the repository for a detailed description of the features.

To install Bioscrape, in your terminal run the command:
   
   pip install bioscrape

Please note that Bioscrape is a Cython extension module and requires a C++ compiler to be set up on your computer for installation.

Try online without installing:
1. [Bioscrape core features](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Basic%20Examples%20-%20START%20HERE.ipynb#scrollTo=Jmm8mTPfhMMS) (modeling and simulation)
2. [Bioscrape sensitivity analysis](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb#scrollTo=Gu1r4H4ti_z7)
3. [Bioscrape inference](https://colab.research.google.com/github/biocircuits/bioscrape/blob/colab-ipynb/inference%20examples/Bioscrape%20Inference%20-%20Getting%20Started.ipynb#scrollTo=yvYVliBgjyzF)

For more information, please visit https://github.com/biocircuits/bioscrape/wiki. If you face any issues, please feel free to raise an issue on this Github.

.. |Build Status| image:: https://github.com/biocircuits/bioscrape/actions/workflows/bioscrape.yml/badge.svg
   :target: https://github.com/biocircuits/bioscrape
