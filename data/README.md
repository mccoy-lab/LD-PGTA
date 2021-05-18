### Aneuploidy data from PGT-A of 8153 embryos

`bluefuse_by_chrom.txt`: Automated inference of chromosome gains and losses for all autosomes, produced by Bluefuse Multi. 

Description of fields:
- `sample_code`: unique identifier of the embryo biopsy
- `cycle_code`: unique identifier of the IVF cycle to which the embryo belongs
- `patient_age`: age of the female patient. "NA" indicates cycles for which this metadata was not available.
- `chrom`: chromosome identifier (only autosomes were analyzed)
- `bluefuse`: inferred ploidy status of the corresponding chromosome using coverage-based analysis with Bluefuse Multi. "G" indicates chromosome gains, "L" indicates chromosome losses, and "NA" indicates normal copy number.
- `dlr`: derivative log ratio, measuring variability in copy number signal across adjacent bins. See Bluefuse Multi manual for description.
- `sample_coverage`: mean depth of coverage of aligned reads across the genome

