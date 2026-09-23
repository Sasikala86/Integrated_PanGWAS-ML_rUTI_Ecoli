panaroo \
-i /home/sastra/Sasikala/Pangenome/panaroo_analysis/*.gff \
-o /home/sastra/Sasikala/Pangenome/panaroo_analysis/panaroo_alignment_results \
--clean-mode strict \
-a core \
--aligner mafft \
--core_threshold 0.99 \
-t 18
