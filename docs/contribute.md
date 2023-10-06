Please contribute to _SADIE_!

## Issues

Questions, feature requests, and bug reports are welcome as [discussions or issues](https://github.com/jwillis0720/sadie/issues/new/choose). **However, to report a security
vulnerability, please see our [security policy](https://github.com/jwillis0720/sadie/security/policy).**

To make it as simple as possible for us to help you, please include the version output in your issue:

```bash
sadie --version
```

Please try to always include the above unless you're unable to install _SADIE_ or **know** it's not relevant
to your question or feature request.

## Pull Requests

_SADIE_ has an automated release. If you submit a pull request, it will be released as soon as it is accepted. This ensures that the latest version of _SADIE_ is always available to the community.

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
   Currently, poetry does not support Python 3.11

1. Install sadie, dependencies, test dependencies, and doc dependencies
<div class="termy">
```console
$ poetry install --with dev
---> 100%
```
</div>

1. Checkout a new branch and make your changes
<div class="termy">
```console
$ git checkout -b my-new-feature-branch
```
</div>

1. Fix formatting and imports
   <div class="termy">
   ```console
   $ pre-commit run --all-files
   ```
   </div>

   !!! Info
   SADIE uses [black](https://github.com/psf/black) to enforce formatting, [isort](https://github.com/PyCQA/isort) to fix imports, and [pyright](https://github.com/microsoft/pyright) for type checking.

1. Run tests and linting
<div class="termy">
```console
$ poetry run pytest tests
```
</div>

1. Build documentation
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
   This is not needed for a pull request but is useful for checking your changes locally.

1. Commit your changes and submit a pull request to the `development` branch

   ... add, commit, push, and create your pull request point to our development branch. Thank you in advance!
