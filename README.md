# frag2lead

Set of scripts to (1) find specific superstructures in chemical libraries, (2) retrieve molecules that engage with a protein through specific interactions, and (3) retrieve molecules that retain a specific binding mode relative to a given input molecule and chemical SMARTS pattern.

## Chemical pattern filter
```
python pattern-filter.py -i input_smiles.ism -p $(cat smarts_pattern.txt) -o output_smiles.ism
```

## Interaction filter
```
python interaction-filter.py -p protein.pdb -i docked_poses.mol.gz -s "760H" -o interacting_compounds.txt
```

## Binding mode filter
```
python rmsd-filter.py -r reference.mol2 -c docked_poses.mol2.gz -t 2.0 -o filtered_poses.mol2 -s $(cat smarts_pattern.txt)
```

# Troubleshooting
If you have any problems, questions or suggestions please open a GitHub Issue/Discussion.

# Citing our work
If you use the UniverseGenerator in your research, please cite our work as follows: [Luttens, A., *et al* Virtual Fragment Screening for DNA repair inhibitors in Vast Chemical Space. Nat Commun 16, 1741 (2025). [https://doi.org/10.1038/s41467-025-56893-9](https://doi.org/10.1038/s41467-025-56893-9).


Thanks for your interest in our work!

Happy modeling! ~ Andreas
