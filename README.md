# UniProt-PDB-mapper

## Quick mapping of UniProt sequences to PDB structures

This script maps residue numbering and chains between a UniProt ID and a PDB structure. It can also find all structures or the best 
structure for a protein. I wrote this because I was fed up with the UniProt website and API, and parsing SIFTS files for 6h every time
I wanted to play with a (real) structure. 

I wrote it to be simple - it doesn't require overcomplicated nextflow and pytorch quantum computing, 
just python and a basic ability to parse json. But I'm sure you can handle that (or you can ask chatgpt, I don't care much). 

## It offers these functionalities:
1. If the user provides a UniProt ID and PDB ID: It Maps residues between the UniProt protein and the specified PDB structure.
2. Automatic mode: Given only a UniProt ID, it identifies the best PDB model based on sequence coverage, 
   resolution and R-factor, retrieves its biological assembly (or asymmetric unit if specified), and performs 
   the residue and chain mapping.
3. Option to list all PDB entries mapped to the UniProt ID with metadata (experimental technique, resolution, 
   sequence coverage, and R-factor).

## Fancy features:
- Supports mapping to either the biological assembly or the asymmetric unit (user-specified).
- Outputs results in a JSON format.

## Usage:
    python uniprot_pdb_mapper.py --uniprot_id <UniProt_ID> [--pdb_id <PDB_ID>] [--assembly <biological|asymmetric>] [--list_all]

## Arguments:
- `--uniprot_id`: Required. The UniProt ID of the protein to map.
- `--pdb_id`: Optional. The PDB ID of the structure to use for mapping.
- `--list_all`: Optional. List all PDB entries mapped to the UniProt ID with metadata.
- `--outfile`: Optional. Output file to write results in .json.

## Dependencies:
- Python 3
- Requests library (install with `pip install requests`)

## Example Usage:
1. Provide both UniProt ID and PDB ID:
    python uniprot_pdb_mapper.py --uniprot_id P69905 --pdb_id 1BZ0 --assembly biological

2. Automatic mode:
    python uniprot_pdb_mapper.py --uniprot_id P69905 

3. List all PDBs mapped to the UniProt ID:
    python uniprot_pdb_mapper.py --uniprot_id P69905 --list_all
