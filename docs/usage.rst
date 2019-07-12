========
Usage
========

.. highlight:: none

You can use refchooser to select a good reference from a list of assemblies. You will need either
a directory of assemblies or a file containing the paths to the assemblies. Refchooser prints a
table of metrics for each assembly.

The captured metrics are:

* N50
* N90
* Number of contigs
* Assembly length
* Mean mash distance to all other assemblies

The results can be sorted by any metric you choose. By default, the assemblies are sorted by a simple score
which is the N50/Distance ratio.

To print the table of assemblies sorted by N50::

    # Choose the top 10 from a collection of 900 assemblies
    refchooser metrics --sort N50 --top 10 assembly_paths.txt sketches

    Assembly   N50    N90    Contigs Length  Mean_Distance Path                    Score
    SRR5868281 791291 119061 30      4845891 7.045173e-04  fasta/SRR5868281.fasta  1.123168e+09
    SRR7439260 775386 146629 24      4815254 4.927176e-04  fasta/SRR7439260.fasta  1.573693e+09
    SRR7906469 775033 432700 18      4779049 6.352714e-04  fasta/SRR7906469.fasta  1.220003e+09
    SRR6949545 774499 146519 21      4814882 5.308503e-04  fasta/SRR6949545.fasta  1.458978e+09
    SRR6949610 774132 105140 33      4888983 8.929775e-04  fasta/SRR6949610.fasta  8.669110e+08
    SRR7426190 774120 146629 30      4820457 5.317999e-04  fasta/SRR7426190.fasta  1.455660e+09
    SRR7426155 774120 146449 29      4775352 6.618484e-04  fasta/SRR7426155.fasta  1.169633e+09
    SRR7441818 774120 146519 25      4797608 5.506614e-04  fasta/SRR7441818.fasta  1.405800e+09
    SRR7439259 774120 146629 25      4815750 4.911346e-04  fasta/SRR7439259.fasta  1.576187e+09
    SRR7439242 774120 146519 32      4803747 5.681594e-04  fasta/SRR7439242.fasta  1.362505e+09


To print the table of assemblies sorted by mean mash distance::

    # Choose the top 10 from a collection of 900 assemblies
    refchooser metrics --sort Mean_Distance --top 10 assembly_paths.txt sketches

    Assembly   N50    N90    Contigs Length  Mean_Distance Path                    Score
    SRR1645597 226490 55522  55      4803421 4.611227e-04  fasta/SRR1645597.fasta  4.911708e+08
    SRR1965968 237440 55508  59      4804728 4.614244e-04  fasta/SRR1965968.fasta  5.145805e+08
    SRR1963305 166064 47353  61      4804588 4.618624e-04  fasta/SRR1963305.fasta  3.595530e+08
    SRR1646405 226711 56774  58      4800826 4.629222e-04  fasta/SRR1646405.fasta  4.897389e+08
    SRR1967694 287598 63327  54      4800251 4.637451e-04  fasta/SRR1967694.fasta  6.201639e+08
    SRR7458586 216691 68846  48      4802064 4.642351e-04  fasta/SRR7458586.fasta  4.667700e+08
    SRR7439539 333102 76943  45      4796679 4.646953e-04  fasta/SRR7439539.fasta  7.168180e+08
    SRR5584738 216691 54594  62      4797573 4.649960e-04  fasta/SRR5584738.fasta  4.660061e+08
    SRR7439240 774109 146629 34      4814374 4.658887e-04  fasta/SRR7439240.fasta  1.661575e+09
    SRR8691682 216324 75764  47      4795530 4.659022e-04  fasta/SRR8691682.fasta  4.643121e+08


To print the table of assemblies sorted by the N50/Mean_Distance ratio score::

    # Choose the top 10 from a collection of 900 assemblies
    refchooser metrics --top 10 assembly_paths.txt sketches

    Assembly   N50    N90    Contigs Length  Mean_Distance Path                    Score
    SRR7439240 774109 146629 34      4814374 4.658887e-04  fasta/SRR7439240.fasta  1.661575e+09
    SRR7439252 774092 146519 26      4810069 4.722259e-04  fasta/SRR7439252.fasta  1.639241e+09
    SRR5237981 749843 146968 30      4811975 4.681156e-04  fasta/SRR5237981.fasta  1.601833e+09
    SRR7439259 774120 146629 25      4815750 4.911346e-04  fasta/SRR7439259.fasta  1.576187e+09
    SRR7439260 775386 146629 24      4815254 4.927176e-04  fasta/SRR7439260.fasta  1.573693e+09
    SRR7140222 773996 105140 27      4810707 4.939575e-04  fasta/SRR7140222.fasta  1.566928e+09
    SRR7426191 774120 146629 28      4813941 5.066112e-04  fasta/SRR7426191.fasta  1.528036e+09
    SRR7347002 774109 146519 34      4822419 5.137192e-04  fasta/SRR7347002.fasta  1.506872e+09
    SRR1793292 774006 119046 56      4833765 5.247949e-04  fasta/SRR1793292.fasta  1.474873e+09
    SRR6945020 774120 146519 21      4813673 5.303144e-04  fasta/SRR6945020.fasta  1.459738e+09
