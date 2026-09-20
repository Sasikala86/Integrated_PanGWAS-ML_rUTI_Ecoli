# Workflow carried out for MOB_suite
                 
       MOB-suite Workflow

MOB-suite v3.1.9 was applied to all 491 E. coli genome assemblies to identify chromosome- and plasmid-associated sequences and reconstruct putative plasmids.

```text
491 E. coli genomes
        │
        ▼
   MOB-suite mob_recon
        │
        ▼
mob_output/<isolate>/
        │
   ┌────┴──────────────┐
   │                   │
   ▼                   ▼
Chromosome          Plasmids
                       │
                       ▼
                Collect & rename
                       │
                       ▼
              Plasmid QC analysis
                       │
              ┌────────┴─────────┐
              │                  │
              ▼                  ▼
        PlasmidFinder       Plasmid Panaroo
        validation
```
Scripts
01_mobsuite_batch.sh — runs mob_recon on all genome assemblies.
02_collect_plasmids.sh — collects plasmid FASTA files from individual MOB-suite output directories.
03_rename_plasmids.sh — renames plasmid FASTA files to retain the isolate identifier.
04_MOB_QC_analysis.py — generates MOB-suite contig- and plasmid-level QC results.
05_run_plasmidfinder_batch.py — performs independent replicon-level assessment of MOB-suite plasmid reconstructions using PlasmidFinder.

Note: MOB-suite was run once on the complete set of 491 genomes. UTI_ and rUTI_ prefixes were used only for isolate/phenotype identification and did not represent separate MOB-suite analyses.
