# Contributing to Bioscrape 

Thank you for your interest in contributing to Bioscrape!

In this file you will find detailed instructions on how you can start making contributions to the package. Bioscrape is hosted on the [Biocircuits](https://github.com/biocircuits) organization page on GitHub. Note that the legacy version of the package was hosted on the lead developer, ananswam's account and is not actively maintained. Use the [biocircuits/bioscrape](https://github.com/biocircuits/bioscrape/) for all purposes. For more information on getting started with the package, refer to the README file on the [home page](https://github.com/biocircuits/bioscrape/) and the tutorial style example jupyter notebooks under the [examples](https://github.com/biocircuits/bioscrape/tree/master/examples) directory. The software documentation is available on the Github Wiki page [here](https://github.com/biocircuits/bioscrape/wiki). 

## How to contribute?
To get started, set up your Bioscrape fork under your account on Github - [raise an issue](https://github.com/biocircuits/bioscrape/issues/new) if you are not sure how to set it up. All contributions to Bioscrape should be made as a Github pull request to the **dev** branch [here](https://github.com/biocircuits/bioscrape/tree/dev). 

### Reporting Bugs/Asking for help

Use the Github issues page on Bioscrape to report a bug or to ask for help with running Bioscrape. The Github issues have labels that you can use so that the issues can be filtered easily:

* If you are unsure where to begin contributing to Bioscrape, you can start by looking through the issues with the label `beginner` or `help-wanted`. A beginner issue usually requires only changing a few lines of code to fix something or add a new enhancement. The `help-wanted` issues are slightly more involved and may require an understanding of the Bioscrape modules.

* If you have a particular feature idea in mind, feel free to suggest that as an `enhancement` or a `feature-request` tagged issue. If you would like to get in touch with the developers working on the package to discuss your contribution ideas, you can also join our Slack channel (details at the end of this page).

### Pull Requests

All pull requests should be made to the `dev` branch of Bioscrape. To maintain code readability and validity, we encourage you to document your pull requests using the following ways:
* Add a detailed comment when creating the pull request that summarizes the changes, features and/or bugs fixed.
* We recommend that all new functions and classes have [docstrings](https://www.python.org/dev/peps/pep-0257/) and other documentation, if needed.
* If possible, we encourage you to add test functions in the Tests directory that validate the code contributions. 
* If your pull request is adding a new feature to the package, we also highly recommend a jupyter notebook example that goes along with that feature that discusses the use case. 

## Styleguides

* Following the PEP8 guideline, limit the first line to 72 characters or less
* Reference issues and pull requests in your pull request comment 

### Docstrings

Docstrings in `bioscrape/*.py`, `*.pyx`, `*.pxd`, and `lineage/*.pyx`
follow the numpydoc convention (see the [numpydoc style guide](
https://numpydoc.readthedocs.io/en/latest/format.html)), prioritizing
what reads well as plain text via `help()`/`obj.__doc__`, since this
package isn't always built alongside its Sphinx documentation.

* Wrap docstring *text* (not code) to ~78 characters, greedily. Don't
  break a wrap point inside a single- or double-backtick span.
* The one-line summary (first line) must fit on a single physical
  line, ≤75 characters, start with a capital letter, and end with a
  period. If it doesn't fit, shorten it and move detail into the
  extended summary or Parameters/Returns instead of wrapping it.
* Avoid LaTeX/math notation (`$...$`) in docstrings -- it doesn't
  render in a terminal. Write formulas in plain text, e.g.
  `rate = k*(x1/K)^n / (1 + (x1/K)^n)` or `d(x_i)/d(p_j)`.
* Use single backticks around Python objects and parameter names
  (e.g. `` `Model` ``, `` `param_name` ``); double backticks for
  inline code/short fragments. Quote plain strings with single quotes
  (e.g. `'uniform'`), not double-backtick-wrapped. Write True/False/
  None with no backticks, capitalized.
* Class docstrings host the constructor's Parameters, not
  `__init__`'s. Give `__init__` a one-line stub,
  `"""See class docstring."""` -- leaving it with no docstring at all
  is unsafe if it subclasses another documented class, since some
  tools then fall back to the parent's `__init__` docs and silently
  describe the wrong constructor. Classes should not have a "Returns"
  section.
* Every parameter and keyword argument must be documented.
  Functions/attributes not meant for user access start with an
  underscore (or, for Cython, are `cdef` rather than `def`/`cpdef`).

## Have any questions?

If you have questions or would like to connect to the Bioscrape team on a regular basis, you can join our Slack channel. Bioscrape is a channel under the Synthetic Biology Modeling and Analysis Tools (SBTools) slack team.

* [Join the SBTools Slack](https://join.slack.com/t/sbtools/shared_invite/zt-g82qjmvm-GAsNFLjyXGPlRBapqGDgFg)
    * Use the `#bioscrape` channel for general questions or discussion about Bioscrape
    * Use the `#general` channel for general questions about SBTools
    * There are many other channels available for other synthetic biology modeling and analysis tools, check the channel list
