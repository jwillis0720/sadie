#Renumbering Antibody Sequences

One of the hardest parts of working with antibody sequences is that have different definitions of numbering. There is also different definitions of where the framework and CDR regions are depending on the scheme. SADIE provides a simple interface to renumber antibody sequences to a common numbering scheme. We borrow heavily from the Antigen receptor Numbering and Receptor Classification ([ANARCI](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/anarci/))

---

## Single Sequence Annotation

```Python
{!docs_src/renumbering/simple_usage.py!}

```

The output will contain `<sadie.renumbering.result.NumberingResults'>` and shows that the output is an instance of the `DataFrame` class.
