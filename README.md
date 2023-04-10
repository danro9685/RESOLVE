RESOLVE
=======

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/danro9685/RESOLVE/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/danro9685/RESOLVE/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/danro9685/RESOLVE/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/danro9685/RESOLVE/actions/workflows/check-bioc.yml) |

Cancer is a genetic disease caused by somatic mutations in genes controlling key biological functions such as cellular growth and division. Such mutations may arise both through cell-intrinsic and exogenous processes, generating characteristic mutational patterns over the genome named mutational signatures. The study of mutational signatures have become a standard component of modern genomics studies, since it can reveal which (environmental and endogenous) mutagenic processes are active in a tumor, and may highlight markers for therapeutic response.

Mutational signatures computational analysis presents many pitfalls. First, the task of determining the number of signatures is very complex and depends on heuristics. Second, several signatures have no clear etiology, casting doubt on them being computational artifacts rather than due to mutagenic processes. Last, approaches for signatures assignment are greatly influenced by the set of signatures used for the analysis. To overcome these limitations, we developed RESOLVE (Robust EStimation Of mutationaL signatures Via rEgularization), a framework that allows the efficient extraction and assignment of mutational signatures.

RESOLVE implements a novel algorithm that enables (i) the efficient extraction, (ii) exposure estimation, and (iii) confidence assessment during the computational inference of mutational signatures.
