# Reference Module

The SADIE reference module abstracts the underlying reference data used by the [AIRR](annotation.md) and [Numbering](renumbering.md) modules. Both of these modules use external database files. Their organization (particularly by AIRR, which ports [IGBlast](https://www.ncbi.nlm.nih.gov/igblast/)) can be extremely complicated. Making a new reference database is a tedious and time-consuming task. This module provides a simple interface for making your own reference databases.

!!! Abstract "Builtin reference"
SADIE ships with a reference database that contains the most common species along with functional genes. The average user will not need to use this module as the database is comprehensive. You can see each entry by looking either directly at the paths used `src/sadie/airr/data/` for AIRR and `src/sadie/anarci/data` for the renumbering module. Another convenient way to look at the reference database is to view the [reference.yml](https://github.com/jwillis0720/sadie/blob/master/src/sadie/reference/data/reference.yml). More on how that file is structured will be [provided](#the-reference-yaml).

## Germline Data Sources

SADIE supports multiple germline data sources for maximum flexibility:

| Source | Description | Use Case |
|--------|-------------|----------|
| **imgt** | [IMGT](https://www.imgt.org) reference database | Default, most comprehensive |
| **ogrdb** | [OGRDB](https://ogrdb.airr-community.org/) curated alleles | Novel/validated alleles |
| **vdjbase** | [VDJbase](https://www.vdjbase.org/) germlines | Population-specific alleles |
| **custom** | Your own sequences | Novel discoveries, proprietary |

!!! tip "New in v1.2"
    SADIE now supports **OGRDB** and **VDJbase** as first-class data sources. You can mix sources in a single reference database to combine curated alleles from different providers.

### Germline Gene Gateway (Legacy)

For backwards compatibility, SADIE also supports the [Germline Gene Gateway](https://g3.jordanrwillis.com/docs/) API which provides IMGT and custom genes. This RESTful API conforms to the [OpenAPI 3.0](https://swagger.io/specification/) specification.

### Examples of how to use the G3 API

The following examples show how to pull genes programmatically using the command line utilities `curl`, `wget` and the `requests` library in Python. It will fetch the first 5 V-Gene segments in IMGT notation.

=== ":material-console-line: curl"

    <div class="termy">

    ```console
    $ curl -X 'GET' 'https://g3.jordanrwillis.com/api/v1/genes?source=imgt&segment=V&common=human&limit=5' -H 'accept: application/json' -o 'human_v.json'
    ```
    </div>

=== ":material-console-line: wget"

    <div class="termy">
    ```console
    $ wget 'https://g3.jordanrwillis.com/api/v1/genes?source=imgt&segment=V&common=human&limit=5' -O human_v.json
    ```
    </div>

=== " :material-api: Python"

    ```Python
    {!> docs_src/reference/tutorial001.py!}
    ```

The output will be a JSON file containing the V-Gene segment and all relevant information needed by SADIE to write out databases needed by the AIRR and Numbering modules.

??? example "human_v.json"

    <div id='json_block_div'>
    ```json
        {!> docs_src/reference/human_v.json!}
    ```
    </div>

!!! Tip "G3 API"

    The G3 API can be explored live through the documentation. Go to the [G3 API Documentation](https://g3.jordanrwillis.com/docs/) to do so. It is a clean non-redundant dataset that can be used for any project programatically. To learn more, [explore the source code](https://github.com/jwillis0720/g3). SADIE abstracts most connections with G3, so you should not have to interact with the API directly.

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
| `source`   | :material-check: The data source for the genes                                     | `imgt`, `ogrdb`, `vdjbase`, `custom` |
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

### Multi-Source Example (v1.2+)

You can now combine genes from multiple sources in a single reference:

```yaml
human-comprehensive:
  # Baseline from IMGT
  imgt:
    human:
      - IGHV1-2*02
      - IGHV3-23*01
      - IGHD3-10*01
      - IGHJ4*02

  # Add curated alleles from OGRDB
  ogrdb:
    human:
      - IGHV1-69*13    # Novel allele validated by OGRDB
      - IGHV4-34*09    # Population-specific variant

  # Include VDJbase alleles
  vdjbase:
    human:
      - IGHV3-30*18    # Curated from VDJbase
```

This allows you to:

- Start with a comprehensive IMGT baseline
- Add novel alleles from OGRDB that have been experimentally validated
- Include population-specific alleles from VDJbase

For more details on multi-source configurations, see [Custom Reference Databases](reference-workflow.md).

Again, a full list of databases, species and genes can be found by exploring the [G3 API](https://g3.jordanrwillis.com/docs#/G3/find_genes_api_v1_genes_get), click the `Try it out` button.

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
