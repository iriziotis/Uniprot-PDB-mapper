# UniProt-PDB mapper

## Quick mapping of UniProt sequences to PDB structures

This script maps residue numbering and chains between a UniProt ID and a PDB structure. It can also find all structures or the best 
structure for a protein. I wrote this because I was fed up with the UniProt website and API, and parsing SIFTS files for 6h every time
I wanted to play with a (real) structure. 

I wrote it to be simple - it doesn't require overcomplicated nextflow and pytorch quantum computing, 
just python and a basic ability to parse json. But I'm sure you can handle that (or you can ask chatgpt, I don't care much). 

## It offers these functionalities:
- If the user provides a UniProt ID and PDB ID: It Maps residues between the UniProt protein and the specified PDB structure.
- Automatic mode: Given only a UniProt ID, it identifies the best PDB model based on sequence coverage and resolution, retrieves
  its biological assembly (or asymmetric unit if specified), and performs the residue and chain mapping.
- Option to list all PDB entries mapped to the UniProt ID with metadata (experimental technique, resolution and coverage).
- Outputs mapping in a JSON format.
- Can automatically download structure files for assemblies, assymetric units or both (it will not unzip them though, this still burdens the user)

## Usage:
    python uniprot_pdb_mapper.py --uniprot_id <UniProt_ID> [--pdb_id <PDB_ID>] [--list_all] [--outfile] [--get-assymetric] [--get-assembly]

## Arguments:
- `--uniprot-id <UniProt_ID>`: Required. The UniProt ID of the protein to map.
- `--pdb-id <PDB_ID>`: Optional. The PDB ID of the structure to use for mapping.
- `--list-all`: Optional. List all PDB entries mapped to the UniProt ID with metadata.
- `--outfile`: Optional. Output file to write results in .json. Defaults to stdout.
- `--get-assymetric`: Optional. Get an unzipped mmCIF with the raw PDB coordinates
- `--get-assembly`: Optional. Get an unzipped mmCIF with the preferred assembly coordinates

## Dependencies:
- Python 3
- Requests library (install with `pip install requests`)

## Example Usage:
1. Provide both UniProt ID and PDB ID:
   
    `python uniprot_pdb_mapper.py --uniprot-id P69905 --pdb-id 1BZ0`

3. Automatic mode:
   
    `python uniprot_pdb_mapper.py --uniprot-id P69905`

5. List all PDBs mapped to the UniProt ID:
   
    `python uniprot_pdb_mapper.py --uniprot-id P69905 --list-all`

7. List all PDBs mapped to the UniProt ID and get their assembly coordinate files
   
    `python uniprot_pdb_mapper.py --uniprot-id Q00526 --list-all --get-assembly`

## Contribution
Please report bugs and request features either by pull request or by emailing me: <ioannis.riziotis@crick.ac.uk>
