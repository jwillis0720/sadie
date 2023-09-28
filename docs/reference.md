# Reference Module

The SADIE reference module abstracts the underlying reference data used by the [AIRR](annotation.md) and [Numbering](renumbering.md) modules. Both of these modules use external database files. Their organization (particularly by AIRR, which ports [IGBlast](https://www.ncbi.nlm.nih.gov/igblast/)) can be extremely complicated. Making a new reference database is a tedious and time-consuming task. This module provides a simple interface for making your own reference databases.

!!! Abstract "Builtin reference"
SADIE ships with a reference database that contains the most common species along with functional genes. The average user will not need to use this module as the database is comprehensive. You can see each entry by looking either directly at the paths used `src/sadie/airr/data/` for AIRR and `src/sadie/anarci/data` for the renumbering module. Another convenient way to look at the reference database is to view the [reference.yml](https://github.com/jwillis0720/sadie/blob/master/src/sadie/reference/data/reference.yml). More on how that file is structured will be [provided](#the-reference-yaml).

## Germline Gene Gateway

New germline gene segments are being discovered at a rapid pace. To meet the needs of this changing landscape, SADIE gets all of the germline gene info from a programmatic API called the [Germline Gene Gateway](https://g3.jordanrwillis.com/docs/). This API is hosted as a free service. It consists of germline genes from [IMGT](https://www.imgt.org) as well as custom genes that have been annotated and cataloged by programs such as [IGDiscover](http://docs.igdiscover.se/en/stable/). To explore the API, visit the [Germline Gene Gateway](https://g3.jordanrwillis.com/docs/). This RESTful API conforms to the [OpenAPI 3.0](https://swagger.io/specification/) specification.

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

- Copyright © Jordan R. Willis and Troy Sincomb
