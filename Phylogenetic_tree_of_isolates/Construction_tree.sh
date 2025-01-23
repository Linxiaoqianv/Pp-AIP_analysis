mafft --auto latest2.total.16s.fa > latest2.16s_all_aligned.fasta
trimal -in latest2.16s_all_aligned.fasta -out latest2.16s_all_aligned_filtered.fasta -automated1 -keepheader
fasttree -quote -nt latest2.16s_all_aligned_filtered.fasta > latest2.16s_all_aligned_filtered.tree
