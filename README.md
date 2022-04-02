DLSuR: Dynamic light scattering microrheology in Python
===========================================

DLSuR is a data analysis package for analyzing the scattering intensity from a dynamic light scattering instrument and deriving the microrheology spectrum in the Python programming language.

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

To use DLSuR, you need to:
* have data from a dynamic light scattering instrument,
* save the data in the specific format that is listed in [this paper](https://pubs.rsc.org/en/content/articlelanding/2021/sm/d0sm01597k), and
* be sure to collect data following the methods listed in [this paper](https://pubs.rsc.org/en/content/articlelanding/2021/sm/d0sm01597k)

# The DLSuR environment

## Easy Implementation

The DLSuR method is simple to implement, utilizing just the scattering autocorrelation of embedded particles in a given soft material sample. The methods are split into different ways to analyze and visualize one's data.

By using only the scattering autocorrelation, the methodology of analyzing the mean-squared displacement of embedded particles to derive the frequency-dependent complex modulus becomes much simpler than other microrheology techniques such as video particle tracking.

## Large Range of Rheological Behavior

DLSuR has the capability of measuring up to six decades in rheological behavior without using time-temperature superposition. This is a major advantage over state-of-the-art rheological techniques such as oscillatory rheometers. 

## How to cite

If you use this package, please cite the following paper:

Cai P. C., Krajina B. A., Kratochvil M. J., Zou L., Zhu A., Burgener E. B., Bollyky P. L., Milla C. E., Webber M. J., Spakowitz A. J., Heilshorn S. C. (2021). Dynamic light scattering microrheology for soft and living materials. *Soft Matter, 17*(7), 1929-1939.

# Installation


## Dependencies

DLSuR requires:

* Python (>= 3.7)
* SciPy 
* NumPy
* Matplotlib
* Pandas
* Seaborn
* Sphinx (>=1.4)


## Standard installation (on CPU hardware)
We strongly recommend running DLSuR in an Anaconda environment, because this simplifies the installation of other
dependencies. The first step is to create a new Anaconda environment:

```
conda create -n myenv
```

This creates an environment called `myenv` (replace the bolded word with whatever you want to name your environment) that you can enter by doing:

```
conda activate myenv
```

Next, you can install the latest version of DLSuR using the package manager `pip`, which will automatically download
DLSuR from the Python Package Index (PyPI):

```
pip install DLSuR
```

Windows, Linux, and macOS are the officially supported operating systems. NOTE: sometimes there will be an error requiring the package `keyring` (version >=15.1).


## Installation from source

Assuming the DLSuR source has been downloaded, you may install it by running

```
pip install -r requirements.txt
python setup.py install
```

#  Analysis Using DLSuR

Once installed, DLSuR can be imported and utilized by using the following at the top of your Python scripts:

```
import dlsmicro
```

Alternatively, you can import the functions, such as `analyze_conditions` using the syntax:

```
from dlsmicro import analyze_conditions
```

Or, sometimes the path to the function is incomplete if using just the above, so you can also try importing the functions using:

```
from dlsmicro.analyze_conditions import analyze_conditions
```

A great place to start is by looking at the file `test_new.py` in the package. This file contains examples of how to use each function included in the package to analyze the data in the `example_data` folder. The data structure within the `example_data` folder also gives you an example of how to set up the file structure for your data in order for the functions within this package to be able to complete the analysis efficiently.


# Documentation

The documentation of DLSuR is officially hosted on the [DLSuR](https://dlsur.readthedocs.io/) website.


## Online resources

* [GitHub repository](https://github.com/PamCai/DLSuR)
* [GitHub issue tracker](https://github.com/PamCai/DLSuR/issues)
* [BSD-3-Clause license](https://github.com/PamCai/DLSuR/blob/master/LICENSE.md)


## Building the documentation from source
The documentation can also be found in the `doc/` subfolder of the GitHub repository.
To build the documentation locally, please clone this repository and run

```
pip install -r requirements_optional.txt
cd doc; make clean; make html
```

## Support

We wish to thank Stanford University, National Science Foundation, Stanford Bio-X Initiative for their financial support.