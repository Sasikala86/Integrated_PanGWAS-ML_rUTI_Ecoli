# Workflow carried out for MOB_suite
                 
                  491 E. coli genomes
                         │
                         ▼
                MOB-suite mob_recon
                         │
                         ▼
              mob_output/<isolate>/
                         │
              ┌──────────┴──────────┐
              │                     │
              ▼                     ▼
         Chromosome             Plasmid
         sequences             sequences
              │                     │
              │                     ▼
              │              Collect plasmids
              │                     │
              │                     ▼
              │              Rename with isolate ID
              │                     │
              │                     ▼
              │              Downstream plasmid
              │               Panaroo analysis
              │
              │
              ▼
       contig_report.txt
              │
              ▼
       MOB-suite QC Python
              │
              ▼
       QC tables and reports
