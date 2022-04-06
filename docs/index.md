<h2 align="center" style="font-family:verdana;font-size:150%"> <b>S</b>equencing <b>A</b>nalysis and <b>D</b>ata Library for <b>I</b>mmunoinformatics <b>E</b>xploration</h2>
<div align="center">
  <img src="https://sadiestaticcrm.s3.us-west-2.amazonaws.com/Sadie.svg" alt="SADIE" style="margin:0.5em;width:50%">
</div>

<div class="flex-container" align="center">
    <!-- <a href="https://github.com/jwillis0720/sadie/commits/master">
    <img src="https://img.shields.io/github/commit-activity/y/jwillis0720/sadie?style=flat-square"
         alt="GitHub commits"> -->
    <a href="https://github.com/jwillis0720/sadie/workflows/Linux%20Build%20and%20Test/badge.svg">
    <img src="https://github.com/jwillis0720/sadie/workflows/Linux%20Build%20and%20Test/badge.svg"
         alt="Linux Build">
    <a href="https://github.com/jwillis0720/sadie/workflows/MacOS%20Build%20and%20Test/badge.svg">
    <img src="https://github.com/jwillis0720/sadie/workflows/MacOS%20Build%20and%20Test/badge.svg"
         alt="Mac Build">
    <a href="https://github.com/jwillis0720/sadie/actions/workflows/mypy.yml/badge.svg">
    <img src="https://github.com/jwillis0720/sadie/actions/workflows/mypy.yml/badge.svg"
         alt="Static Type">
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
    <a href=https://pypi.org/project/sadie-antibody/">
    <img src="https://img.shields.io/pypi/v/sadie-antibody?color=blue"
        alt='pypi'>
    <a href="https://sadie.jordanrwillis.com" >
    <img src="https://api.netlify.com/api/v1/badges/59ff956c-82d9-4900-83c7-758ed21ccb34/deploy-status"
        alt="Documentation">
    </a>
</div>

## About

---

<!-- use a href so you can use _blank to open new tab -->

**Documentation**: <a href="https://sadie.jordanrwillis.com" target="_blank">https://sadie.jordanrwillis.com</a>

**Source Code**: <a href="https://github.com/jwillis0720/sadie" target="_blank">https://github.com/jwillis0720/sadie</a>

---

SADIE is the **S**equencing **A**nalysis and **D**ata library for **I**mmunoinformatics **E**xploration. The key feautures include:

- Provide pre-built **command line apps** for popular immunoinformatics applications.

- Provide a **low-level API framework** for immunoinformatics developers to build higher level tools.

- Provide a **testable** and **reusable** library that WORKS!

- Provide a **customizable** and **verified** germline reference library.

- Maintain data formats consistent with standards governed by the [**AIRR community**](https://docs.airr-community.org/en/stable/#table-of-contents)

- **Portability** ready to use out the box.

SADIE is billed as a "**complete antibody library**", not because it aims to do everything, but because it aims to meet the needs of all immunoinformatics users. SADIE contains both low, mid and high level functionality for immunoinformatics tools and workflows. You can use SADIE as a framework to develop your own tools, use many of the prebuilt contributed tools, or run it in a notebook to enable data exploration. In addition, SADIE aims to port all code to python because relies heavily on the [Pandas](https://www.pandas.org) library, the workhorse of the data science/machine learning age.

## Installation

---

Installation is handled using the python package installer `pip`

<div class="termy">

```console
$ pip install sadie-antibody
---> 100%
```

</div>

### Development installation.

!!! info
Pull requests are highly encouraged [here](https://github.com/jwillis0720/sadie/pulls). The development installation uses [pre-commit](https://pre-commit.com/), [flake8](https://flake8.pycqa.org/en/latest/) linting and [black](https://github.com/psf/black) style formatting to maintain code readability and reausability.

<div class="termy">

A development installation is available by running the following command:

```console
$ git clone git@github.com/jwillis0720/sadie.git
$ pip install -e .[dev]
```

</div>

## The Littlest Usage

Consult the [documentation](https://sadie.jordanrwillis.com) for complete usage

### Command Line Usage

<!-- get these icons through icon search https://squidfunk.github.io/mkdocs-material/reference/icons-emojis/#search -->

=== ":material-console-line: Command Line Usage"
Annotate antibody sequences only from functional human imgt antibodies to a gzip output

    <div class="termy">

    ```console
    $ sadie airr -s human -db-type imgt my_sequences.fasta output.csv
    ```

    </div>

=== " :material-api: Python"
Use the SADIE library to annotate seqeuences

    ```Python
    {!> docs_src/annotation/tutorial002.py!}
    ```

## License

[![License](https://img.shields.io/github/license/jwillis0720/sadie)](https://opensource.org/licenses/MIT)

- Copyright © Jordan R. Willis and Troy Sincomb

<!-- Page last revised on: {{ git_revision_date }} -->
