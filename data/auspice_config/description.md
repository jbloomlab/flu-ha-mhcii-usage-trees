## Overview
These are trees of influenza {tree} hemagglutinin that can be colored (using the *Color By* dropdown) by their predicted or measured usage of MHC class II for cell entry.
Specifically, you can color:

  - The experimentally measured effects of how a strain's constituent mutations affect either binding to MHCII or entry into cells expressing MHCII instead of sialic acid (compared to normal sialic acid expressing cells); note that these values are still only predicted as they assume mutation effects can be extrapolated additively across genetic backgrounds. The coloring can be done according to the total effect of all mutations, the single mutation with the maximum magnitude effect, or the identity of the mutation with the maximum magnitude effect.

  - In some cases, direct measurements of the log10 titer change of entry into cells expressing different MHCII proteins, sialic acid, or no sialic acid.

You can also color sites by amino-acid identity (using the *Genotype* option in *Color By* dropdown).
The site numbering is sequential 1, 2, ... numbering of each of the three protein regions (*SigPep*, *HA1*, and *HA2*).
For a mapping of that site numbering to H3 numbering, compare the *protein* and *protein_site* columns to the *H3_site* column in the TSV at [https://github.com/jbloomlab/flu-ha-mhcii-usage-trees/blob/main/data/trees/{tree}/protein_sites.tsv](https://github.com/jbloomlab/flu-ha-mhcii-usage-trees/blob/main/data/trees/{tree}/protein_sites.tsv).

See [https://github.com/jbloomlab/flu-ha-mhcii-usage-trees](https://github.com/jbloomlab/flu-ha-mhcii-usage-trees) for the code used to generate these trees.
