# gVCFCreator

A graphical interface (Tkinter) for generating genomic VCFs (`.g.vcf.gz`) from BAM/CRAM files with GATK HaplotypeCaller (`-ERC GVCF`).

## Features

- Batch processing of multiple BAM/CRAM files, one gVCF per sample
- Optional restriction to target regions (BED / interval_list) with configurable padding
- Configurable quality filters: minimum MAPQ, minimum base quality, calling confidence threshold
- Configurable ploidy, PairHMM threads and Java memory
- Optional inclusion of reads marked as duplicates
- Automatic creation of the reference `.fai` and `.dict` if missing
- Handling of BAMs containing several samples (different `SM` tags):
  - one gVCF per sample (`--sample-name`), or
  - merge into a single sample (`AddOrReplaceReadGroups`)
- BAMs without read groups are handled automatically (a read group is added)
- Input validation, live log, progress bar, and the ability to stop a running analysis
- Existing outputs are skipped unless overwriting is enabled

## Requirements

- Python 3.8+ with Tkinter
- GATK 4.x (`gatk` launcher or `gatk-package-*.jar`) and Java 17
- samtools (optional: needed only to create the reference `.fai` and to read CRAM headers)

No third-party Python packages are required.

## Usage

```bash
python3 gvcfcreator.py
```

1. Set the GATK path (auto-detected if `gatk` is in `PATH`) and the Java memory.
2. Select the reference genome (FASTA) and the output folder.
3. Optionally select a regions file and padding.
4. Add the BAM/CRAM files. The index (`.bai`/`.csi` or `.crai`) must be located next to each file.
5. Adjust the quality filters and parameters if needed.
6. Click **Validate inputs**, then **Run**.

## Output

For each sample, in the output folder:

| File | Description |
|------|-------------|
| `<sample>.g.vcf.gz` | genomic VCF |
| `<sample>.g.vcf.gz.tbi` | tabix index |
| `<sample>.haplotypecaller.log` | full GATK log |
| `<sample>.addrg.log` | AddOrReplaceReadGroups log (only when read groups are rewritten) |

For BAMs with multiple samples processed separately, outputs are named `<bam name>.<SM>.g.vcf.gz`.

## Author

Serena Spartano

## License

Copyright (c) 2026 Serena Spartano. All rights reserved. See [LICENSE](LICENSE).
