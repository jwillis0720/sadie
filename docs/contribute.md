Please contribute to *SADIE*!

!!! Key Notes
- SADIE does not support he M1/M2 Mac chip set architecture.
  - SADIE uses [pyhmmer](https://github.com/althonos/pyhmmer), an amazing library for Cython bindings to the sequence aligner [HMMER3](http://hmmer.org/). Due to HMMER3 not supporting the M1/M2 chip sets, SADIE does not support them either. HMMER3 will support it with it's next merge from their develop branch. When that happens we will update SADIE to support it as well.

## Issues

Questions, feature requests and bug reports are all welcome as [discussions or issues](https://github.com/jwillis0720/sadie/issues/new/choose). **However, to report a security
vulnerability, please see our [security policy](https://github.com/jwillis0720/sadie/security/policy).**

To make it as simple as possible for us to help you, please include the version output in your issue:

```bash
sadie --version
```

Please try to always include the above unless you're unable to install *SADIE* or **know** it's not relevant
to your question or feature request.

## Pull Requests

*SADIE* has an automated release. This means that if you submit a pull request, it will be released as soon as it is accepted. This is to ensure that the latest version of *SADIE* is always available to the community.

You'll need to have a version between **Python 3.8 and 3.10**, **poetry**, and **git** installed.

1. Clone your fork from Github and cd into your repo directory

    <div class="termy">
        ```console
        $ git clone git@github.com:YOUR_USERNAME/sadie.git
        $ cd sadie
        ```
    </div>

1. Set up a poetry for running tests
    <div class="termy">
    ```console
    $ pip install poetry
    ---> 100%
    ```
    </div>

    !!! Info
        Currently poetry does not support python 3.11

2. Install sadie, dependencies, test dependencies and doc dependencies
    <div class="termy">
    ```console
    $ poetry install --with dev
    ---> 100%
    ```
    </div>

3. Checkout a new branch and make your changes
    <div class="termy">
    ```console
    $ git checkout -b my-new-feature-branch
    ```
    </div>

4. Fix formatting and imports
    <div class="termy">
    ```console
    $ pre-commit run --all-files
    ```
    </div>

    !!! Info
        SADIE uses [black](https://github.com/psf/black) to enforce formatting, [isort](https://github.com/PyCQA/isort) to fix imports, and [pyright](https://github.com/microsoft/pyright) for type checking.

5. Run tests and linting
    <div class="termy">
    ```console
    $ poetry run pytest tests
    ```
    </div>

6. Build documentation
    <div class="termy">
    ```console
    $ mkdocs build
    $ mkdocs serve
    INFO     -  Building documentation...
    INFO     -  [14:27:11] Serving on http://127.0.0.1:8000/
    ```
    <div>

    !!! Info
        Our netlify.toml is used to create our documentation site.
        This is not needed for a pull request, but is useful for checking your changes locally.

7. Commit your changes and submit a pull request to the `development` branch

    ... add, commit, push, and create your pull request point to our development branch. Thank you in advance!
