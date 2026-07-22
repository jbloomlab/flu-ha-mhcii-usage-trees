## Overview
These are trees of influenza {tree} hemagglutinin that can be colored (using the *Color By* dropdown) by their predicted or measured usage of MHC class II for cell entry as measured by [Dadonaite et al (2026)](https://doi.org/10.64898/2026.07.17.738765) using pseudovirus deep mutational scanning or measurements of pseudovirus cell entry in different cell variants.

Specifically, you can color the tree in various ways using the *Color By* dropdown:

  - The experimentally measured effects of how a strain's constituent mutations affect either binding to MHCII or entry into cells expressing MHCII instead of sialic acid (compared to normal sialic acid expressing cells); note that these values are _only predicted as they assume mutation effects can be extrapolated additively across genetic backgrounds_. The coloring can be done according to the total effect of all mutations, the single mutation with the maximum magnitude effect, or the identity of the mutation with the maximum magnitude effect. These colorings have names of the form *MHCII entry increase (H5) total effect* to indicate coloring by the total (summed) effect of all mutations in the strain based on deep mutational scanning of a H5 HA. There are similar colorings for the max effect of any constituent mutation and the mutation with max effect, and mutation effects measured in H5 or H7 deep mutational scanning.

  - Some strains have been directly measured for the titer of pseuodoviruses expressing that HA on 293 cells expressing sialic acid, no sialic acid, or different MHC II proteins. Such strains have the *has_titers* key set to *yes* (so they can be filtered by that), and then the colorings have names like *log10_293-noSA-tufted-duck-MHCII_over_293-noSA_titer* for a coloring givin the log10 ration of the titer on 293-noSA-tufted-duck-MHCII cells to 293-noSA cells.

  - HA amino-acid identity (using the *Genotype* option in *Color By* dropdown).
    The site numbering is sequential 1, 2, ... numbering of each of the three protein regions (*SigPep*, *HA1*, and *HA2*).
    For a mapping of that site numbering to H3 numbering, compare the *protein* and *protein_site* columns to the *H3_site* column in the TSV at [https://github.com/jbloomlab/flu-ha-mhcii-usage-trees/blob/master/data/trees/{tree}/protein_sites.tsv](https://github.com/jbloomlab/flu-ha-mhcii-usage-trees/blob/master/data/trees/{tree}/protein_sites.tsv).

See [https://github.com/jbloomlab/flu-ha-mhcii-usage-trees](https://github.com/jbloomlab/flu-ha-mhcii-usage-trees) for the code used to generate these trees, and [Dadonaite et al (2026)](https://doi.org/10.64898/2026.07.17.738765) for more study details.
