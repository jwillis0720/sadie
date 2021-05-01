<h1 align="center">
  <br>
  <img src="images/Social3.png" alt="SADIE" style="width:80%">
</h1>

<div class="flex-container" align="center">
    <a href="https://github.com/jwillis0720/sadie/commits/master">
    <img src="https://img.shields.io/github/commit-activity/y/jwillis0720/sadie?style=flat-square"
         alt="GitHub commits">
    <a href="https://github.com/jwillis0720/sadie/workflows/Linux%20Build%20and%20Test/badge.svg">
    <img src="https://github.com/jwillis0720/sadie/workflows/Linux%20Build%20and%20Test/badge.svg"
         alt="Linux Build">
    <a href="https://github.com/jwillis0720/sadie/workflows/MacOS%20Build%20and%20Test/badge.svg">
    <img src="https://github.com/jwillis0720/sadie/workflows/MacOS%20Build%20and%20Test/badge.svg"
         alt="Mac Build">
    <a href="https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9-blue">
    <img src="https://img.shields.io/badge/Python-3.7%7C3.8%7C3.9-blue"
        alt="Python Version">
    <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg"
        alt="Format Version">
    <a href="https://codecov.io/gh/jwillis0720/sadie">
    <img src="https://codecov.io/gh/jwillis0720/sadie/branch/master/graph/badge.svg?token=EH9QEX4ZMP"
        alt="Code Coverage">
    <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white"
        alt="pre commit">
    <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/pypi/v/sadie-antibody?color=blue"
        alr='pypi'>
</div>

<p align="center" style="color:green">
  <a href="#about">About</a> •
  <a href="#installation">Installation</a> •
  <a href="#usage">Usage Quickstart</a> •
  <!-- <a href="#contributing">Contributing</a> • -->
  <!-- <a href="#credits">Credits</a> • -->
  <!-- <a href="#support">Support</a> • -->
  <a href="#license">License</a>
</p>

## About

 SADIE is the **S**equencing **A**nalysis and **D**ata library for **I**mmunoinformatics **E**xploration. The goals of the SADIE project is to:

- Provide pre-built command line apps for popular immunoformatics applications.
- Provide a low level API framework for immunoformatics developers to build higher level tools.
- Provide testable and reusable library that WORKS!
- Maintain data formats consistent with standards governed by the [AIRR community](https://docs.airr-community.org/en/stable/#table-of-contents)

SADIE is billed as a "complete antibody library" because it contains both low, mid and high level functionality for immunoformatics tools and workflows. You can use SADIE as a framework to develop your own tools, use many of the prebuilt contributed tools, or run it in a notebook to enable data exploration. In addition, SADIE aims to port all code to python because relies heavily on the [Pandas](www.pandas.org) library, the workhorse of the data science/machine learning age.

## Installation

Installation can be handled with the python package installer `pip`

```
# Install with pip
pip install sadie-antibody


#Or Install and develop
pip install -r requirements.txt
pip install -e .
```

## Usage

SADIE is divided into modules depending on the task.

### SADIE.AIRR

At the heart of every workflow is the need to annotate a nucleotide sequence. This is accomplished using `sadie.airr`.

#### Command line usage
```
# annotate antibody sequences only from functional human imgt antibodies
# sequences will output to a gzipped csv
airr -q my_sequecnes.fasta --functional -s human -d imgt -z gzip
```

#### API

```
# define a single sequence
pg9_seq = """
    CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
    CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
    GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
    AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
    TTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
    TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC""".replace(
    "\n", ""
)

# initialize the api
air_api = Airr("human")

# run single sequence
airr_table = air_api.run_single("PG9", pg9_seq)

# or run file
airr_table = air_api.run_file("myfile.fasta")
```

## License

[![License](https://img.shields.io/github/license/jwillis0720/sadie)](https://opensource.org/licenses/MIT)

- Copyright © Jordan R. Willis
