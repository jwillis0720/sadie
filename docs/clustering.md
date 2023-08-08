# Clustering

Clustering is handled by using an input of [AirrTables](annotation.md#single-sequence-annotation). Consider the following example.

```Python
{!docs_src/cluster/tutorial001.py!}
```

will print out:

```
cluster name: 35 contains ['PCT64-18B', 'PCT64-18C', 'PCT64-18D', 'PCT64-18F', 'PCT64-24F', 'PCT64-24G', 'PCT64-24H', 'PCT64-35B', 'PCT64-35C', 'PCT64-35D', 'PCT64-35E', 'PCT64-35F', 'PCT64-35G', 'PCT64-35H', 'PCT64-35I', 'PCT64-35K', 'PCT64-35N', 'PCT64-35O', 'PCT64-35S']

cluster name: 15 contains ['VRC26.08', 'VRC26.09', 'VRC26.20', 'VRC26.21', 'VRC26.22', 'VRC26.26', 'VRC26.27', 'VRC26.28', 'VRC26.29', 'VRC26.30', 'VRC26.31']

cluster name: 0 contains ['DH270.1', 'DH270.2', 'DH270.3', 'DH270.4', 'DH270.5', 'DH270.6', 'DH270.IA1', 'DH270.IA2', 'DH270.IA3', 'DH270.IA4']

cluster name: 1 contains ['VRC38.02', 'VRC38.03', 'VRC38.04', 'VRC38.05', 'VRC38.08', 'VRC38.09', 'VRC38.10', 'VRC38.11']

cluster name: 34 contains ['10J4', '10M6', '13I10', '2N5', '35O22', '4O20', '7B9', '7K3']

cluster name: 2 contains ['PGT151', 'PGT152', 'PGT153', 'PGT154', 'PGT155', 'PGT156', 'PGT157', 'PGT158']

...
```

So what happened here.

- We read in a `LinkedAirrTable` which has the heavy and light table paired together. This is also a pandas dataframe so we can use any pandas method on it.

- Instantiate a Cluster object that contains the following parameters.

| Parameter     | Description                                                                                                                                             |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `linkage`     | The hierarchical clustering linkage method. [single, average or complete linkage](https://en.wikipedia.org/wiki/Hierarchical_clustering)                |
| `groupby`     | Do we pre-group by any fields. e.g. v_call_heavy will only cluster things that are the same v_call_heavy.                                               |
| `lookup`      | What fields to take the Levenshtein distance to make a cluster distance.                                                                                |
| `pad_somatic` | Any somatic mutations that are present in both sequences will be subtracted from the total distance. This is useful for somatic hypermutation analysis. |

- Run `cluseter_api.cluster(distace)` where `distance` is a distance cutoff in hierarchical clustering. This will return an airrtable with the a field called cluster.
