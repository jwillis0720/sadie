#AIRR Annotation

Annotation is the bedrock of all immunoformatics workflows. It is the process of identifying CDRs/frameworks, levels of somatic mutation, locus use, productive rearrangements, and other features that describe the B cell receptor or T cell receptor (BCR/TCR). In the description of a BCR/TCR, how can we compare the data file output from one data pipeline to another? In other words, what if the description of a repertoire has different fields and datatypes that describe a repertoire or even a single BCR/TCR? Fear not! [The AIRR community to the rescue](https://docs.airr-community.org/en/stable/)!

---

"_AIRR Data Representations are versioned specifications that consist of a file format and a well-defined schema[...] The schema defines the data model, field names, data types, and encodings for AIRR standard objects. Strict typing enables interoperability and data sharing between different AIRR-seq analysis tools and repositories[...]_"

<a href='https://docs.airr-community.org/en/stable/datarep/overview.html'><div style="text-align: right; margin-right: 10%;"> The AIRR Standards 1.3 documentation </div></a>

---

SADIE leverages the AIRR to provide a standardized data representation for BCRs. You can read all the fields and values in the AIRR Rearrangement schema standard [here](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields)

## Single Sequence Annotation

```Python
{!docs_src/annotation/tutorial001.py!}

```

The output will contain `<class 'sadie.airr.airrtable.airrtable.AirrTable'>` and shows that the output is an instance of the `AirrTable` class.

!!! info

    Running an AIRR method generates an AIRR table object. The AIRR table is a subclass of a [pandas dataframe](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html) and thus can be used by any pandas method. Pandas is the workhorse of the SADIE library, so we highly encourage some rudimentary knowledge of pandas to get maximize SAIDIE functionality.

### Writing Files

#### AIRR Rearrangement File

To output an AIRR file, we can use the `AirrTable.to_airr()` method.

```Python
{!docs_src/annotation/tutorial002.py!}
```

The tsv file `PG9 AIRR.tsv` generated will be a tabular datafile that will resemble the following:

{!docs_output/annotation/tutorial001.md!}

This `.tsv` file is a [Rearrangement Schema compliant AIRR table](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#file-format-specification). These files have certain specifications, including a `.tsv` file suffix. Since they are AIRR compliant, they can be used by other [AIRR compliant software.](https://docs.airr-community.org/en/stable/resources/rearrangement_support.html). For instance, we could use the output `.tsv ` in any module in the [immcantation portal](https://immcantation.readthedocs.io/en/stable/).

#### Other Output Formats

While the `.tsv` AIRR table is the recognized standard for AIRR, you can also output to any other formats that [pandas supports](https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html).

```Python
{!docs_src/annotation/tutorial003.py!}
```

!!! attention

    Because `AirrTable` is a subclass of `pandas.DataFrame`, you can use any pandas IO methods to write to a file of your choosing. However, it must be noted that these are not official [Rearrangement Schema compliant AIRR tables](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#file-format-specification). They may only be read in by software that reads those file types or be read back in by **SADIE** and probably will not work in other software that supports the AIRR standard. But, these file formats are extremely useful for much larger files.

### Reading Files

To read in an AIRR file, we have to create an `AirrTable` object.

#### Reading an AIRR.tsv

You can read official AIRR.tsv using the `AirrTable.from_airr()` method or with pandas and casting to an `AirrTable` object.

```Python
{!docs_src/annotation/tutorial004.py!}
```

Outputs:

```output
<class 'sadie.airr.airrtable.airrtable.AirrTable'> True
<class 'sadie.airr.airrtable.airrtable.AirrTable'> True
True # The airr tables are equal
```

#### Reading other file formats

Any other file formats that are readable by [pandas IO](https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html) can be read in by passing them to AirrTable.

```Python
{!docs_src/annotation/tutorial005.py!}
```
