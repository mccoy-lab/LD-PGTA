### Aneuploidy data from PGT-A of 8153 embryos

`bluefuse_by_chrom.txt`: Inference of chromosome gains and losses for all autosomes, produced by Bluefuse Multi. 

Description of fields:
- `sample_code`: unique identifier of the embryo biopsy
- `cycle_code`: unique identifier of the IVF cycle to which the embryo belongs
- `patient_age`: age of the female patient; "NA" indicates cycles for which this metadata was not available
- `chrom`: chromosome identifier (only autosomes were analyzed)
- `bluefuse`: inferred ploidy status of the corresponding chromosome using coverage-based analysis with Bluefuse Multi. "G" indicates chromosome gains, "L" indicates chromosome losses, and "NA" indicates normal copy number
- `dlr`: derivative log ratio, measuring variability in copy number signal across adjacent bins; see Bluefuse Multi manual for description
- `sample_coverage`: mean depth of coverage of aligned reads across the genome


`ldpgta_by_chrom.txt`: Inference of haplotype patterns for all autosomes, produced by LD-PGTA.

Description of fields:
- `sample_code`: unique identifier of the embryo biopsy
- `cycle_code`: unique identifier of the IVF cycle to which the embryo belongs
- `patient_age`: age of the female patient; "NA" indicates cycles for which this metadata was not available
- `chrom`: chromosome identifier (only autosomes were analyzed)
- `mean_llr_bph_sph`: mean of the log-likelihood ratio of comparison of the BPH and SPH trisomy hypotheses for the chromosome, computed by LD-PGTA
- `sd_llr_bph_sph`: standard deviation of the log-likelihood ratio of comparison of the BPH and SPH trisomy hypotheses for the chromosome, computed by LD-PGTA
- `mean_llr_monosomy_disomy`: mean of the log-likelihood ratio of comparison of the monosomy and disomy hypotheses for the chromosome, computed by LD-PGTA
- `sd_llr_monosomy_disomy`: standard deviation of the log-likelihood ratio of comparison of the monosomy and disomy hypotheses for the chromosome, computed by LD-PGTA
- `mean_llr_bph_disomy`: mean of the log-likelihood ratio of comparison of the BPH trisomy and disomy hypotheses for the chromosome, computed by LD-PGTA
- `sd_llr_bph_disomy`: standard deviation of the log-likelihood ratio of comparison of the BPH trisomy and disomy hypotheses for the chromosome, computed by LD-PGTA
