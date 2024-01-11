# gsmmutils


## Description

This is a collection of utilities for working with Genome-Scale Metabolic Models and associated experimental data.
It offers a set of functions for (some of them are executed with docker containers):
- Reading and writing GSMMs in SBML format
- Reading and writing experimental data in tabular format
- Performing genome analysis like InterProScan searches and BUSCO analysis
- Performing dynamic flux balance analysis (dFBA) with dFBA (https://gitlab.com/davidtourigny/dynamic-fba)
- Integrating omics data into GSMMs with _Troppo_ and MEWpy


## Installation

To install gsmmutils, clone the repository and run the setup.py script:

```
git clone  https://github.com/ecunha1996/gsmmutils.git
cd gsmmutils
python setup.py install
```

## Credits and License
Developed at:
-  Centre of Biological Engineering, University of Minho (2020-)

