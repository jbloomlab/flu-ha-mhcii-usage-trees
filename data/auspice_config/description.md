## Overview
These are trees of influenza {tree} hemagglutinin that can be colored (using the *Color By* dropdown) by their predicted usage of MHC class II based on summing the effects of how their constituent mutations affect either binding to MHCII or entry into cells expressing MHCII instead of sialic acid (compared to normal sialic acid expressing cells).

The coloring can be done according to the total effect of all mutations, the single mutation with the maximum magnitude effect, or the identity of the mutation with the maximum magnitude effect.

You can also color sites by amino-acid identity (using the *Genotype* option in *Color By* dropdown).
The site numbering is sequential 1, 2, ... numbering of each of the three protein regions (*SigPep*, *HA1*, and *HA2*).
For a mapping of that site numbering to H3 numbering, compare the *protein* and *protein_site* columns to the *H3_site* column in the TSV at [https://github.com/jbloomlab/flu-ha-trees-mhcII-binding/blob/main/data/trees/{tree}/protein_sites.tsv](https://github.com/jbloomlab/flu-ha-trees-mhcII-binding/blob/main/data/trees/{tree}/protein_sites.tsv).

For the numerical data on the mutation effects mapped on this tree, see [https://github.com/jbloomlab/flu-ha-trees-mhcII-binding/blob/main/data/trees/{tree}/mutation_effects.tsv](https://github.com/jbloomlab/flu-ha-trees-mhcII-binding/blob/main/data/trees/{tree}/mutation_effects.tsv).

See [https://github.com/jbloomlab/flu-ha-trees-mhcII-binding](https://github.com/jbloomlab/flu-ha-trees-mhcII-binding) for the code used to generate these trees.

