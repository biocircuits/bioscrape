---
title: 'Fast and flexible simulation and parameter estimation for synthetic biology using bioscrape'
tags:
  - Python
  - synthetic biology
  - systems biology
  - deterministic and stochastic simulations
  - parameter inference
authors:
  - name: Ayush Pandey^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0003-3590-4459
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: William Poole^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0002-2958-6776
    affiliation: "2" # (Multiple affiliations must be quoted)
  - name: Anandh Swaminathan^[Co-first author] # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0001-9935-6530
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
 - name: Control and Dynamical Systems, California Institute of Technology, Pasadena, CA, USA
   index: 1
 - name: Altos Labs, San Francisco, CA, USA
   index: 2
 - name: Ghost Locomotion, Mountain View, CA, USA
   index: 3
 - name: Amyris, Emeryville, CA, USA
   index: 4
 - name: Control and Dynamical Systems and Biology and Biological Engineering, California Institute of Technology, Pasadena, CA, USA
   index: 5
date: 26 Jan 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

In systems and synthetic biology, it is common to build chemical reaction network (CRN) models of biochemical circuits and networks. Although automation and other high-throughput techniques have led to an abundance of data enabling data-driven quantitative modeling and parameter estimation, the intense amount of simulation needed for these methods still frequently results in a computational bottleneck. Here we present bioscrape (Bio-circuit Stochastic Single-cell Reaction Analysis and Parameter Estimation) - a Python package for fast and flexible modeling and simulation of highly customizable chemical reaction networks. Specifically, bioscrape supports deterministic and stochastic simulations, which can incorporate delay, cell growth, and cell division. All functionalities - reaction models, simulation algorithms, cell growth models, partitioning models, and Bayesian inference - are implemented as interfaces in an easily extensible and modular object-oriented framework. Models can be constructed via Systems Biology Markup Language (SBML) or specified programmatically via a Python API. Simulation run times obtained with the package are comparable to those obtained using C code - this is particularly advantageous for computationally expensive applications such as Bayesian inference or simulation of cell lineages. We show the package's simulation capabilities on a variety of example simulations of stochastic gene expression. We also demonstrate the package by using it to do parameter inference on a model of integrase enzyme-mediated DNA recombination dynamics with experimental data. The bioscrape package is publicly available online [@github_bioscrape] along with more detailed documentation and examples.



# Statement of need

A central theme of research in systems and synthetic biology is the quantitative predictions of the behavior of complex biological systems and efficient testing of hypotheses. Mathematical modeling and analysis tools play an integral role in such predictions and can transform the way in which we design and debug synthetic engineered biological circuits. 

Cell growth and division are critical aspects of biological circuits which are typically represented as a dilution term in the model. However, in stochastic models, modeling the continuous dilution process with a stochastic and discrete degradation reaction might not be accurate. Moreover, the partitioning of molecules between daughter cells at cell division may introduce noise that is difficult to distinguish from other forms of noise [@paulsson_partition]. Therefore, modeling cell growth as well as division and partitioning is important for investigating noise in gene expression across a lineage of cells.

Regardless of simulation framework, it is necessary to first specify the values of the parameters of each propensity function in the model along with the initial levels of the model species. In some cases, these parameters and initial conditions are experimentally known. Often, however, they have to be inferred from biological data via a process known as parameter inference, parameter estimation, or parameter identification [@sun2012parameter]. Bayesian inference [@golightly2011bayesian;@komorowski_bayesian_2009] is one of the most rigorous methods of parameter identification. It provides a posterior distribution over the parameter space so that the stochastic effects from the experimental data are modeled by the parameter distributions instead of a fixed optimal point. This gives insight into the accuracy and identifiability of the model. Also, such an approach allows for an easy comparison between different model classes using the model evidence. The drawback of these approaches is that their implementation is computationally expensive and is based on repeated forward simulations of the model within the framework of Markov chain Monte Carlo (MCMC) [@golightly2011bayesian]. Therefore, it is important to have the underlying simulations running as fast as possible in order to speed up computation time.

Once a given model is fully specified, it is then important to validate the model against additional biological data. In this workflow, it is often necessary to add or remove reactions from the model or to perform a different type of simulation. For example, one might decide that a circuit behaves too noisily for deterministic simulations and want to switch to a stochastic simulation framework. If delays are playing a significant role in the dynamics, one might want to incorporate previously unmodeled delays into the model. 

![(a) A simple model of gene expression with transcription, translation, mRNA degradation, and protein degradation. The quantity of the gene encoding for mRNA is considered constant and absorbed into the transcription rate $\beta$. (b) Example Python code to construct a CRN model of gene expression using Bioscrape. (c) Models constructed via SBML or the Python API can be easily simulated with results returned as a Pandas Dataframe [@mckinney-proc-scipy-2010]. (d) Deterministic and stochastic simulations (with and without delays) using Bioscrape.The empirical probability distribution and the autocorrelation function for mRNA in the stochastic simulation matches the theoretical Poisson and exponential curve respectively
\label{fig:model_simulation}](examples/joss_figure.pdf)

The result is that a very large amount of data is needed to first parameterize and then validate models. The use of technologies for lab automation makes this data collection increasingly accessible and economical. For deterministic models, this may include data collected at many different operating conditions which can be achieved with high throughput measurement techniques involving liquid handling automation [@freemont_echo]. For stochastic models this may include large sample sizes of single cell measurements such as flow cytometry [@sachs_causal_2005;@zechner2012moment] and tracking single cell lineages with fluorescent microscopy [@kretzschmar2012lineage]. The Python API, simulation tools, and lineage module in bioscrape provide an ideal platform for such applications.

Some popular software packages that do somewhat similar tasks as the bioscrape package are MATLAB's SimBiology toolbox [@MATLAB_2016], Stochpy [@stochpy], COPASI [@copasi], and Tellurium [@tellurium]. The SBML simulator libRoadRunner [@libroadrunner1;@libroadrunner2] is the state-of-the-art in deterministic and stochastic simulations of SBML models. Bioscrape simulation performance is of the same order as libRoadRunner and only around 20-30% slower while being an order of magnitude faster than other GUI-based simulators such as MATLAB and COPASI. However, the bioscrape package provides features beyond SBML simulations as it supports fully general propensity functions, provides easy-to-use parameter identification interfaces, and allows simulation of delays, and cell populations. The target audience for bioscrape includes researchers from diverse fields such as systems biology, synthetic biology, and chemical engineering. It is also aimed as an educational tool for classes on mathematical and computational biology.

# Summary of features

\autoref{fig:model_simulation} shows an example of gene expression model created and simulated stochastically and deterministically using Bioscrape. 

We conclude with a list of Bioscrape features:

1. Bioscrape provides a Cython [@cython] based simulator that compiles code using a C compiler to vastly increase speed. This helps assuage the computational time issues that arise in parameter estimation and stochastic simulation. 
2. Kinds of possible simulations include: deterministic, stochastic, growing and dividing lineages of single cells, and stochastic simulation of delayed chemical reaction networks. A flexible easy-to-use wrapper and a Python API make it straightforward for a researcher to change their model and try simulations under diverse conditions. 
3. Markov Chain Monte Carlo (MCMC) sampler based inference tools to identify parameter distributions of biological circuit models using experimental data. Bioscrape provides interfaces to easily use common biological data types such as time-series fluorescence data and flow cytometry data. The MCMC sampler is a wrapper around Python emcee [@emcee].
4. Bioscrape can be used to perform local sensitivity analysis of a model to study the sensitivities of each parameter with time.

[@pandey2022characterization] demonstrates Bioscrape's features for quantification and predictive modeling of an engineered biological system.

# Acknowledgements

AP, AS, and VH were supported by the Defense Advanced Research Projects Agency (Agreement HR0011-17-2-0008). The content of the information does not necessarily reflect the position or the policy of the Government, and no official endorsement should be inferred. AP was also supported by the NSF Grant CBET-1903477 and AFOSR MURI Grant FA9550-22-1-0316. WP was supported by an NSF Graduate Research Fellowship (No.2017246618). AS was also supported by AFOSR Grant FA9550-14-1-0060. 

The authors acknowledge members of the Murray Biocircuits lab at Caltech for assistance with experiments and helpful feedback and also acknowledge all the members of the scientific community at large who have used and provided feedback on bioscrape.


# References
