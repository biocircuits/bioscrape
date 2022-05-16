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

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

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

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
