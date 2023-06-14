SurfaceMatch
============

(c) 1993 Andrew C.R. Martin
---------------------------

SurfaceMatch is a simple pair of programs for trying to match two
surfaces for guiding docking.  It was written while self-employed and
visiting DKfz.

`surface` takes a PDB file and generates a list of surface features.
Optionally it can take a file containing residue ranges to restrict
the analysis to part of the protein surface.

`match` takes two output files from `surface`.  The first is a pattern
from one structure which is searched across the pattern of the second
structure.  The first pattern will be smaller than the second.  For
example, the first might be a peptide or ligand being searched against
a protein.  Alternatively, the first might be likely contact residues
in the CDRs of an antibody being searched against an antigen.

### Compilation

You must install BiopLib (http://www.bioinf.org.uk/software/bioplib)
and can then simply run the `compile.sh` script.

