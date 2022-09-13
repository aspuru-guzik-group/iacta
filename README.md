# Automatic Discovery of Chemical Reactions Using Imposed Activation
This repository contains code for the paper: [Automatic Discovery of Chemical Reactions Using Imposed Activation](https://doi.org/10.26434/chemrxiv.13008500.v2). 

Original authors: Cyrille Lavigne, Gabe Gomes, Robert Pollice, Alán Aspuru-Guzik 


## Prerequsites: 

`iacta` needs [`python >= 3.7`](https://www.python.org/downloads/). Before being able to use the repository, you will need to install a few additional programs and packages. First, common packages that are required are `numpy` and `pandas`.

Second, you will need to install [`xtb`](https://xtb-docs.readthedocs.io/en/latest/contents.html) version >= 6.3.0. Importantly, by default the `xtb` and `crest` binaries are assumed to be in `$PATH`.

Third, you will need to install [openbabel](https://open-babel.readthedocs.io/en/latest/Installation/install.html) version >= 3.1.1 and the associated python bindings. The easiest is to do this is via conda:

```
conda install -c conda-forge openbabel
```

Finally, you also need to install `yaml` and `pyaml`, again easiest to do via conda:

```
conda install -c conda-forge pyyaml
```

Afterwards, simply clone the GitHub repository and installation is complete.

## How to run: 

There are two ways to run `iacta`. The easiest is via `yaml` input files. There are a few example input files for the case studies shown in the paper in the subdirectories of the `test-set` and `other-examples` directories. These input files are used with `rsearch-restart.py` for `iacta` simulations as follows:

```
python3 /path/to/ts-search/rsearch-restart.py /path/to/user.yaml -o path/to/output/
```

The second possibility is to use `rsearch.py`.

## Questions, problems?
Make a github issue. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (gabegomes[AT]cmu[DOT]edu, r[DOT]pollice[AT]rug[DOT]nl)


## License

iacta is provided under the MIT license.
