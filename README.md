[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5610348.svg)](https://doi.org/10.5281/zenodo.5610348)

**Haplotype-aware inference of human chromosome abnormalities**

<img src="https://rajivmccoy.files.wordpress.com/2021/05/method.png" alt="drawing" width="400"/>

Extra or missing chromosomes—a phenomenon termed aneuploidy—frequently arises during human meiosis and embryonic mitosis and is the leading cause of pregnancy loss, including in the context of *in vitro* fertilization (IVF). While meiotic aneuploidies affect all cells and are deleterious, mitotic errors generate mosaicism, which may be compatible with healthy live birth. Large-scale abnormalities such as triploidy and haploidy also contribute to adverse pregnancy outcomes, but remain hidden from standard sequencing-based approaches to preimplantation genetic testing (PGT-A). The ability to reliably distinguish meiotic and mitotic aneuploidies, as well as abnormalities in genome-wide ploidy may thus prove valuable for enhancing IVF outcomes. Here, we describe a statistical method for distinguishing these forms of aneuploidy based on analysis of low-coverage whole-genome sequencing data, which is the current standard in the field. Our approach overcomes the sparse nature of the data by leveraging allele frequencies and linkage disequilibrium (LD) measured in a population reference panel. The method, which we term LD-informed PGT-A (LD-PGTA), retains high accuracy down to coverage as low as 0.05x and at higher coverage can also distinguish between meiosis I and meiosis II errors based on signatures spanning the centromeres. LD-PGTA provides fundamental insight into the origins of human chromosome abnormalities, as well as a practical tool with the potential to improve genetic testing during IVF.

Embryo aneuploidy data and results described in the manuscript can be found here: https://github.com/mccoy-lab/LD-PGTA/tree/master/data

A tutorial and description of the software can be found here: https://github.com/mccoy-lab/LD-PGTA/wiki/LD-PGTA-tutorial

© 2021 The Johns Hopkins University 
