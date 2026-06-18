# Reference Module

The SADIE reference module abstracts the underlying reference data used by the [AIRR](annotation.md) and [Numbering](renumbering.md) modules. Both of these modules use external database files. Their organization (particularly by AIRR, which ports [IGBlast](https://www.ncbi.nlm.nih.gov/igblast/)) can be extremely complicated. Making a new reference database is a tedious and time-consuming task. This module provides a simple interface for making your own reference databases.

!!! Abstract "Builtin reference"
SADIE ships with a reference database that contains the most common species along with functional genes. The average user will not need to use this module as the database is comprehensive. You can see each entry by looking either directly at the paths used `src/sadie/airr/data/` for AIRR and `src/sadie/anarci/data` for the renumbering module. Another convenient way to look at the reference database is to view the [reference.yml](https://github.com/jwillis0720/sadie/blob/master/src/sadie/reference/data/reference.yml). More on how that file is structured will be [provided](#the-reference-yaml).

## Bundled Germline Gene Data

New germline gene segments are being discovered at a rapid pace. To meet the needs of this changing landscape, SADIE ships bundled germline gene records from [IMGT](https://www.imgt.org) and custom genes annotated by programs such as [IGDiscover](http://docs.igdiscover.se/en/stable/). Reference generation reads these package resources offline.

### Example: read bundled IMGT genes offline

The example below reads `src/sadie/reference/data/imgt-g3.json.gz` through `importlib.resources` and writes the first 5 human V-gene segments in IMGT notation.

```Python
{!> docs_src/reference/tutorial001.py!}
```

The output will be a JSON file containing the V-gene segment records and all relevant information needed by SADIE to write out databases used by the AIRR and Numbering modules.

??? example "human_v.json"

    <div id='json_block_div'>
    ```json
        {!> docs_src/reference/human_v.json!}
    ```
    </div>

!!! Tip "Offline reference data"

    No network call is needed. SADIE resolves reference genes from bundled package data; this direct gzip read is only for inspecting the records used by the docs.

## Generating AIRR Reference Database

=== ":material-console-line: Command Line Usage"

    <div class="termy">

    ```console
    $ {!> docs_src/reference/tutorial002.bash!}
    ```

    </div>

=== " :material-api: Python"

    ```Python
    {!> docs_src/reference/tutorial002.py!}
    ```

## The reference YAML

The reference YAML file is a simple YAML file that takes the following structure.

```yaml
name:
  database:
    species:
    -gene1
    -gene2
    species2:
    -gene3
    -gene4
```

| Field      | Description                                                                        | Example                 |
| ---------- | ---------------------------------------------------------------------------------- | ----------------------- |
| `name`     | :material-check: The name that this reference will be called in SADIE              | `human`, `mouse`, `clk` |
| `database` | :material-check: The database that the gene comes from                             | `IMGT` or `custom`      |
| `species`  | :material-check: The name of the species that will be used in the annotation table | `human`, `mouse`        |
| `gene`     | :material-check: The full gene name                                                | `IGHV3-23*01`           |

!!! Info "Why do we allow multiple species?"

    Most of the time the name and species will be the same thing. i.e.

    ```yaml
    human
        imgt:
            human:
                -IGHV3-23*01
                -IGHD3-3*01
                -IGHJ6*01
    ```

    However, sometimes, you may work with chimeric models where a transgene is knocked into a model species. Consider the HuGL mouse models from [Deli et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32873644/)

    ```yaml
    hugl18:
        imgt:
            human:
            - IGHV4-59*01
            - IGHD3-3*01
            - IGHJ3*02
            mouse:
            - IGHV1-11*01
            - IGHV1-12*01
            - IGHV1-13*01
            - IGHV1-14*01
        ...
    ```

    The HuGL18 model will have the full mouse background and three gene segments knocked-in from a human.

Again, a full list of built-in databases, species and genes can be found in the bundled [reference.yml](https://github.com/jwillis0720/sadie/blob/master/src/sadie/reference/data/reference.yml) and the package data under `src/sadie/reference/data/`.

## Generating AIRR database with Reference Class

Rather than generate a pre-configured database, SADIE can also generate a reference file on the fly. This is useful for procedural analysis, where you generate custom genes for multiple species.

```Python
{!> docs_src/reference/tutorial003.py!}
```

or we can use the YAML file as a template to add more genes

```Python
{!> docs_src/reference/tutorial004.py!}
```

- Copyright © Jordan R. Willis, Troy Sincomb, and Caleb K. Kibet
