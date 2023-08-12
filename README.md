UW Department of Atmospheric Sciences research sample code
==========================
Austin Ahlstrom

This is a fork of [http://github.com/darothen/pyrcel](), a package that we used as a point of reference for experimental computational research during the time I spent doing mentored graduate research with [Dr. Peter Blossey](https://atmos.washington.edu/~bloss/) and [Dr. Robert Wood](https://faculty.washington.edu/robwood2/wordpress/) of University of Washington's [Department of Atmospheric Sciences](https://atmos.uw.edu/).

This fork doesn't particularly add anything or attempt to contribute to the Pyrcel package itself; it simply makes use of it as a baseline for comparison with other methods of water vapor and cloud formation modeling that we were looking into. Forking the repo just seemed like the easiest way to reference it.

Some of the experimentation we were doing was based on a 2007 paper by Morrison and Grabowski: [https://journals.ametsoc.org/view/journals/atsc/65/3/2007jas2374.1.xml](https://journals.ametsoc.org/view/journals/atsc/65/3/2007jas2374.1.xml), specifically in the appendix when it describes the computational method used for modeling. A description of the simplified version we were using for comparison is included in a PDF file within this root directory.

I was using [Anaconda](https://www.anaconda.com/) for running this code locally. The environment settings I was using have been exported into the file environment.yaml.

The file [supersaturation_compare_methods.py](https://github.com/aaustin9/atmo/blob/master/supersaturation_compare_methods.py) contains an example of the kind of code I was working on over the course of my research, comparing a simple Euler's method approach to Pyrcel's more sophisticated parcel model. To run this file, after cloning this directory in a local location enabled with Anaconda:

* Import the conda environment: `conda env create -n ENVNAME --file environment.yml`, using whatever environment name ENVNAME you want
* Run the setup file for Pyrcel: `python setup.py build`
* Run the file with my code: `python supersaturation_compare_methods.py`

The script `supersaturation_compare_methods.py` displays several graphs and prints out some console output. There are some slight discrepancies between the results of using the two different methods, which was essentially what we were looking into when I graduated from the master's program at UW and was no longer an eligible student research assistant. I was going to format this work as a Jupyter notebook, but my local Jupyter notebook is apparently currently messed up and I never got it working with the dependency libraries for this work 😅