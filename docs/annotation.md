
#AIRR Annotation

Annotation is the bedrock of all immunoformatics workflows. It is the process of identifying CDRs/frameworks, levels of somatic mutation, locus use, productive rearragements, and other features that describe the B cell receptor or T cell recptor (BCR/TCR). In the descritiption of a BCR/TCR, how can we use the data file output from one data pipeline can be compared to another? In other words, what if the description of a reperotire has different fields and datatypes that describe a repertoire or even a single BCR/TCR? Fear not! [The AIRR community to the rescue](https://docs.airr-community.org/en/stable/)!

---

"_AIRR Data Representations are versioned specifications that consist of a file format and a well-defined schema[...] The schema defines the data model, field names, data types, and encodings for AIRR standard objects. Strict typing enables interoperability and data sharing between different AIRR-seq analysis tools and repositories[...]_"

<a href='https://docs.airr-community.org/en/stable/datarep/overview.html'><div style="text-align: right; margin-right: 10%;"> The AIRR Standards 1.3 documentation </div></a>

---

SADIE leverages the AIRR to provide a standardized data representation for BCRs. You can read all the fields and values in the AIRR Rearrangment shema standard [here](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields)

## Single Sequence Annotation

```Python
{!docs_src/annotation/tutorial001.py!}

```

The output will contain `<class 'sadie.airr.airrtable.airrtable.AirrTable'>` and shows that the output is an instance of the `AirrTable` class.

!!! info
    Running an AIRR method generates an AIRR table object. The AIRR table is a subclass of a [pandas dataframe](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) and thus can be used by any pandas method. Pandas is the workhorse of the SADIE library so we highly encourage some rudimentary knowledge of pandas to get maximize SAIDIE functionality.

### Writing Files

#### AIRR Rearrangment File

To output an AIRR file, we can use the `AirrTable.to_airr()` method.

```Python hl_lines="13 16-17"
{!docs_src/annotation/tutorial002.py!}
```

The tsv file `PG9 AIRR.tsv` generated will be a tabular datafile that will resemble the following:

{!docs_output/annotation/tutorial001.md!}

This `.tsv` file is a [Rearrangement Schema compliant AIRR table](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#file-format-specification). These files have a certain specification, including a `.tsv` file suffix. Since they are AIRR compliant, they can be used by other [AIRR compliant software.](https://docs.airr-community.org/en/stable/resources/rearrangement_support.html). For instance, we could use the output `.tsv ` in any module in the [immcantation portal](https://immcantation.readthedocs.io/en/stable/).

#### Other Output Formats

While the `.tsv` AIRR table is the recognized standard for AIRR, you can also output for any other formats that [pandas supports](https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html).

```Python hl_lines="13 16-17"
{!docs_src/annotation/tutorial003.py!}
```
