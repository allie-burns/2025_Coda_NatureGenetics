# 2024_Coda

This is the code used to perform the ChIP analyses in the manuscript, "Locus-specific epigenetic editing in engram cells regulates memory retention" (Coda et al., 2024).

### Abstract
With their Janus-faced property of being at once dynamic and stable, epigenetic mechanisms have for long been proposed to act as potential molecular mnemonics, but a cell type-specific, locus-restricted and temporally controllable interrogation thereof has thus far been lacking. Over the past years, accumulating evidence has shown that memories are likely encoded in sparse populations of brain cells, so-called engrams, which have, with few exceptions, received little molecular attention. Here, we combine c-Fos driven engram tagging with CRISPR-based epigenetic editing in vivo to assess whether and the extent to which the epigenetic regulation of a single locus within engram cells can contribute to memory formation and storage. Focusing on the promoter region of Arc, a master regulator of learning and synaptic plasticity, we find that its temporally restricted epigenetic regulation within engram cells in the mouse dentate gyrus is both necessary and sufficient for memory retention after learning as well as for memory maintenance after recall. Furthermore, such epigenetic editing and its behavioral consequence is reversible, and capable of bidirectionally altering memory capacities even outside the labile phase of memory consolidation. Together, these findings suggest that epigenetic mechanisms are causally involved in regulating memory expression.

### Analysis Information
This analysis was performed in order to determine off-target effects of Cas9.

Fastq files were aligned and processed using Bowtie2, samtools and picard (`/code/1_runAlignment.sh`).  The R tool, csaw, was used to call peaks (`/code/2_PeakCalling_csaw.R`).  And Cas9 targets were determined by performing a differential enrichment analysis between IP and input samples (`/code/3_DEanalysis_edgeR.R`)
