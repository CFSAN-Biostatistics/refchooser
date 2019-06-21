========
Usage
========

.. highlight:: bash


You can use refchooser to select a good reference from a list of assemblies. The recommended approach is
to first prepare a list of assemblies with low Mash distance to all the others, then examine those assemblies
to verify the number of contigs and total size are acceptable.


To find the assemblies with the lowest average Mash distance to all the other assemblies::

    refchooser sketch assembly_paths.txt sketch_directory
    refchooser distance --top 25 sketch_paths.txt


To find the assemblies with the fewest contigs::

    refchooser contigs --top 25 assembly_paths.txt
