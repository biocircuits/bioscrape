---
title: 'Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape'
tags:
  - Python
  - synthetic biology
  - systems biology
  - deterministic and stochastic simulations
  - parameter inference
authors:
  - name: Anandh Swaminathan^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0001-9935-6530
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: William Poole^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0002-2958-6776
    affiliation: "2" # (Multiple affiliations must be quoted)
  - name: Ayush Pandey^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0003-3590-4459
    affiliation: "3" # (Multiple affiliations must be quoted)
  - name: Victoria Hsiao # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0001-9297-1522
    affiliation: "4" # (Multiple affiliations must be quoted)
  - name: Richard M Murray # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0002-5785-7481
    affiliation: "5" # (Multiple affiliations must be quoted)
  # - name: Author Without ORCID^[Co-first author] # note this makes a footnote saying 'Co-first author'
  #   affiliation: 2
  # - name: Author with no affiliation^[Corresponding author]
  #   affiliation: 3
affiliations:
 - name: Ghost Locomotion, Mountain View, CA, USA
   index: 1
 - name: Computation and Neural Systems, California Institute of Technology, Pasadena, CA, USA
   index: 2
 - name: Control and Dynamical Systems, California Institute of Technology, Pasadena, CA, USA
   index: 3
 - name: Amyris, Emeryville, CA, USA
   index: 4
 - name: Control and Dynamical Systems and Biology and Biological Engineering, California Institute of Technology, Pasadena, CA, USA
   index: 5
date: 17 May 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
In systems and synthetic biology, it is common to build chemical reaction network (CRN) models of biochemical circuits and networks. Although automation and other high-throughput techniques have led to an abundance of data enabling data-driven quantitative modeling and parameter estimation, the intense amount of simulation needed for these methods still frequently results in a computational bottleneck. Here we present bioscrape (Bio-circuit Stochastic Single-cell Reaction Analysis and Parameter Estimation) - a Python package for fast and flexible modeling and simulation of highly customizable chemical reaction networks. Specifically, bioscrape supports deterministic and stochastic simulations, which can incorporate delay, cell growth, and cell division. All functionalities - reaction models, simulation algorithms, cell growth models, partioning models, and Bayesian inference - are implemented as interfaces in an easily extensible and modular object-oriented framework. Models can be constructed via Systems Biology Markup Language (SBML) or specified programmatically via a Python API. Simulation run times obtained with the package are comparable to those obtained using C code - this is particularly advantageous for computationally expensive applications such as Bayesian inference or simulation of cell lineages. We first show the package's simulation capabilities on a variety of example simulations of stochastic gene expression. We then further demonstrate the package by using it to do parameter inference on a model of integrase enzyme-mediated DNA recombination dynamics with experimental data. The bioscrape package is publicly available online (https://github.com/biocircuits/bioscrape) along with more detailed documentation and examples.



# Statement of need

In the fields of systems and synthetic biology, it has become increasingly common to build mathematical models of biochemical networks. In principle, such models allow for quantitative predictions of the behavior of complex biological systems and efficient testing of hypotheses regarding how real biological networks function. Such predictions would transform the way in which we design and debug synthetic engineered biological circuits. 


Biological circuits can often be noisy~\cite{elowitz_stochastic_2002,eldar_functional_2010}, especially in single cells with low molecular copy numbers~\cite{paulsson2005models}. In these cases, a stochastic model is often necessary to capture the noise characteristics of a circuit.

Stochastic simulation also allows for the inclusion of delay into chemical reactions. Processes like protein production are not instantaneous, and there is often a significant delay between when transcription of a gene is initiated and when a mature protein is produced. This type of delay can lead to non-trivial behavior such as oscillations~\cite{stricker}, and thus it is often important to incorporate delay into the modeling framework. 

Cell growth and division are also critical aspects of biological circuits that operate in single cells. Typically, a dilution term in the model accounts for cell growth. However, in stochastic models, modeling the continuous dilution process with a stochastic and discrete degradation reaction might not be accurate. Another source of noise is the partitioning of molecules between daughter cells at cell division, which can be difficult to distinguish from other forms of noise~\cite{paulsson_partition}. Therefore, modeling cell growth as well as division and partitioning is important for investigating noise in gene expression across a lineage of cells.

Regardless of simulation framework, it is necessary to first specify the values of the parameters of each propensity function in the model along with the initial levels of the model species. In some cases, these parameters and initial conditions are experimentally known. Often, however, they have to be inferred from from biological data via a process known as parameter inference, parameter estimation, or parameter identification~\cite{sun2012parameter}. Bayesian inference~\cite{golightly2011bayesian,komorowski_bayesian_2009} is one of the most rigorous methods of parameter identification. It provides a posterior distribution over the parameter space so that the stochastic effects from the experimental data are modeled by the parameter distributions instead of a fixed optimal point. This gives insight into the accuracy and identifiability of the model. Also, such an approach allows for an easy comparison between different model classes using the model evidence. The drawback of these approaches is that their implementation is computationally expensive and is based on repeated forward simulations of the model within the framework of Markov chain Monte Carlo (MCMC)~\cite{golightly2011bayesian}. Therefore, it is important to have the underlying simulations running as fast as possible in order to speed up computation time.

Once a given model is fully specified, it is then important to validate the model against additional biological data. In this workflow, it is often necessary to add or remove reactions from the model or to perform a different type of simulation. For example, one might decide that a circuit behaves too noisily for deterministic simulations and want to switch to a stochastic simulation framework. If delays are playing a significant role in the dynamics, one might want to incorporate previously unmodeled delays into the model.

The result is that a very large amount of data is needed to first parameterize and then validate models. The use of technologies for lab automation makes this data collection increasingly accessible and economical. For deterministic models, this may include data collected at many different operating conditions which can be achieved with high throughput measurement techniques involving liquid handling automation~\cite{freemont_echo}. For stochastic models this may include large sample sizes of single cell cell measurements such as flow cytometry~\cite{sachs_causal_2005,zechner2012moment} and tracking single cell lineages with fluorescent microscopy~\cite{kretzschmar2012lineage}. 

# Summary of features

This paper presents bioscrape (Bio-circuit Stochastic Single-cell Reaction Analysis and Parameter Estimation), which is a Python package for fast and flexible modeling and simulation of biological circuits. The bioscrape package uses Cython~\cite{cython}, an extension for Python that compiles code using a C compiler to vastly increase speed. This helps assuage the computational time issues that arise in parameter estimation and stochastic simulation. Bioscrape provides an object oriented framework which allows for easily customizable models that can be simulated in many different ways including deterministically, stochastically, or as growing and dividing lineages of single cells. Flexible easy-to-use wrapper and a Python API make it straightforward for a researcher to change their model and try simulations under diverse conditions. Some popular software packages that do somewhat similar tasks as the bioscrape package are MATLAB's SimBiology toolbox~\cite{MATLAB_2016} and Stochpy~\cite{stochpy}. However, the bioscrape package is faster, supports fully general propensity functions, and allows more kinds of simulation than these alternatives making it more flexible and more efficient than alternative packages.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

The figure \autoref{model_simulation} shows an example.

![(a) A simple model of gene expression with transcription, translation, mRNA degradation, and protein degradation. The quantity of the gene encoding for mRNA is considered constant and absorbed into the transcription rate $\beta$. (b) Example Python code to construct a CRN model of gene expression using Bioscrape. (c) Models constructed via SBML or the Python API can be easily simulated with results returned as a Pandas Dataframe~\cite{mckinney-proc-scipy-2010}. (d) Deterministic and stochastic simulations (with and without delays) using Bioscrape.The empirical probability distribution and the autocorrelation function for mRNA in the stochastic simulation matches the theoretical Poisson and exponential curve respectively.\label{model_simulation}](examples/joss_figure.pdf)

# Acknowledgements

AS, AP, and VH were supported by the Defense Advanced Research Projects Agency (Agreement HR0011-17-2-0008). The content of the information does not necessarily reflect the position or the policy of the Government, and no official endorsement should be inferred. AS was also supported by AFOSR grant FA9550-14-1-0060. AP was also supported by the NSF grant CBET-1903477.
WP was supported by an NSF Graduate Research Fellowship (No.2017246618).

The authors acknowledge members of the Murray lab at Caltech for assistance with experiments and helpful feedback and also acknowledge all the members of the scientific community at large who have used and provided feedback on bioscrape.


# References
