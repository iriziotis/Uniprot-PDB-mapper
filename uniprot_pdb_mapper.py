#!/usr/bin/env python
"""
Residue and Chain Mapper for UniProt and PDB

This script maps residue numbering and chains between a UniProt ID and a PDB structure. 
It offers the following functionalities:
1. User-provided UniProt ID and PDB ID: Maps residues between the UniProt protein and the specified PDB structure.
2. Automatic mode: Given only a UniProt ID, the script identifies the best PDB model based on sequence coverage, 
   resolution, and R-factor, retrieves its biological assembly (or asymmetric unit if specified), and performs 
   the residue and chain mapping.
3. Option to list all PDB entries mapped to the UniProt ID with metadata (experimental technique, resolution, 
   sequence coverage, and R-factor).

Features:
- Supports mapping to either the biological assembly or the asymmetric unit (user-specified).
- Outputs results in a JSON format.

Usage:
    python mapper.py --uniprot_id <UniProt_ID> [--pdb_id <PDB_ID>] [--assembly <biological|asymmetric>] [--list_all]

Arguments:
- `--uniprot_id`: Required. The UniProt ID of the protein to map.
- `--pdb_id`: Optional. The PDB ID of the structure to use for mapping.
- `--list_all`: Optional. List all PDB entries mapped to the UniProt ID with metadata.
- `--outfile`: Optional. Output file to write results in .json.

Dependencies:
- Python 3
- Requests library (install with `pip install requests`)

Example Usage:
1. Provide both UniProt ID and PDB ID:
    python mapper.py --uniprot_id P69905 --pdb_id 1BZ0 --assembly biological

2. Automatic mode:
    python mapper.py --uniprot_id P69905 

3. List all PDBs mapped to the UniProt ID:
    python mapper.py --uniprot_id P69905 --list_all
"""

import requests
import argparse
import json
import sys


def fetch_pdb_entries(uniprot_id):
    """Fetch all PDB entries for a given UniProt ID."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json().get(uniprot_id, [])
    else:
        raise ValueError(f"Failed to fetch PDB entries for UniProt ID: {uniprot_id}")


def select_best_pdb_model(pdb_entries):
    """Select the best PDB model based on sequence coverage, resolution, and R-factor."""
    if not pdb_entries:
        raise ValueError("No PDB entries found for this UniProt ID.")

    sorted_entries = sorted(
        pdb_entries,
        key=lambda entry: (
            -entry["coverage"],  # Higher coverage is better
            entry.get("resolution", float("inf")),  # Lower resolution is better
            entry.get("r_factor", float("inf")),  # Lower R-factor is better
        ),
    )
    return sorted_entries[0]  # Best entry


def fetch_sifts_mapping(pdb_id):
    """Fetch SIFTS residue and chain mapping for the given PDB ID."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(f"Failed to fetch SIFTS mapping for PDB ID: {pdb_id}")


def fetch_biological_assembly(pdb_id):
    """Fetch the preferred biological assembly for a PDB ID."""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        assembly_data = response.json().get(pdb_id, [{}])[0].get('assemblies', [])
        if not assembly_data:
            raise ValueError(f"No biological assembly found for PDB ID: {pdb_id}")
        return next((assembly for assembly in assembly_data if assembly["preferred"] == True), assembly_data[0])
    else:
        raise ValueError(f"Failed to fetch biological assembly for PDB ID: {pdb_id}")


def map_residues(uniprot_id, pdb_id):
    """Map residue numbering and chains between UniProt and PDB."""
    sifts_data = fetch_sifts_mapping(pdb_id)

    if pdb_id not in sifts_data or "UniProt" not in sifts_data[pdb_id]:
        raise ValueError(f"No SIFTS mapping found for PDB ID: {pdb_id} and UniProt ID: {uniprot_id}")

    uniprot_mappings = sifts_data[pdb_id]["UniProt"]
    if uniprot_id not in uniprot_mappings:
        raise ValueError(f"No mapping found for UniProt ID: {uniprot_id} in PDB ID: {pdb_id}")

    chain_mappings = uniprot_mappings[uniprot_id]["mappings"]

    residue_map = {}
    for mapping in chain_mappings:
        pdb_chain = mapping["chain_id"]
        for i, unp_residue in enumerate(range(mapping["unp_start"], mapping["unp_end"] + 1)):
            try:
                pdb_residue = mapping["start"]['residue_number'] + i
            except TypeError:
                pdb_residue = None
            try:
                pdb_author_residue = mapping["start"]['author_residue_number'] + i
            except TypeError:
                pdb_author_residue = None
            residue_map[unp_residue] = {"pdb_residue": pdb_residue, 
                                        "pdb_author_residue": pdb_author_residue,  
                                        "pdb_chain": pdb_chain}
    return residue_map


def main():
    # Command line arguments
    parser = argparse.ArgumentParser(description="Residue and Chain Mapper for UniProt and PDB.")
    parser.add_argument("-u", "--uniprot_id", required=True, 
                        help="The UniProt ID of the protein.")
    parser.add_argument("-p", "--pdb_id", 
                        help="The PDB ID of the structure to map (optional).")
    parser.add_argument("--list_all", action="store_true",
                        help="List all PDB entries mapped to the UniProt ID with metadata.")
    parser.add_argument("-o", "--outfile", default=None,
                        help="File to write results in .json")

    # Parse arguments
    args = parser.parse_args()

    try:
        uniprot_id = args.uniprot_id
        pdb_id = args.pdb_id
        list_all = args.list_all
        outfile = args.outfile

        # Outfile handling
        if outfile:
            o = open(outfile, 'w')
        else:
            o = sys.stdout

        if list_all:
            # Mode: List all PDBs
            print(f"Listing all PDB models mapped to UniProt ID: {uniprot_id}.", file=sys.stderr)
            pdb_entries = fetch_pdb_entries(uniprot_id)
            if not pdb_entries:
                print(f"No PDB entries found for UniProt ID: {uniprot_id}", file=sys.stderr)
            else:
                pdb_list = []
                for entry in pdb_entries:
                    assembly = fetch_biological_assembly(entry['pdb_id'])
                    assembly_id = assembly["assembly_id"]
                    residue_map = map_residues(uniprot_id, entry["pdb_id"])
                    pdb_list.append({    
                        "pdb_id": entry["pdb_id"],
                        "experimental_method": entry["experimental_method"],
                        "resolution": entry.get("resolution"),
                        "sequence_coverage": entry["coverage"],
                        "r_factor": entry.get("r_factor"),
                        "preferred_assembly": assembly_id,
                        "residue_map": residue_map})
                print(json.dumps({"uniprot_id": uniprot_id, "pdb_entries": pdb_list}, indent=4), file=o)

        elif pdb_id:
            # Mode: UniProt ID and PDB ID provided
            print(f"Mapping residues for UniProt ID: {uniprot_id} and PDB ID: {pdb_id}.", file=sys.stderr)
            assembly = fetch_biological_assembly(pdb_id)
            assembly_id = assembly["assembly_id"]
            residue_map = map_residues(uniprot_id, pdb_id)
            output = {
                "uniprot_id": uniprot_id, 
                "pdb_id": pdb_id, 
                "preferred_assembly": assembly_id,
                "residue_map": residue_map}
            print(json.dumps(output, indent=4), file=o)

        else:
            # Mode: Automatic mapping
            pdb_entries = fetch_pdb_entries(uniprot_id)
            best_pdb = select_best_pdb_model(pdb_entries)
            pdb_id = best_pdb["pdb_id"]
            assembly = fetch_biological_assembly(pdb_id)
            assembly_id = assembly["assembly_id"]
            print(f"Selected PDB ID for {uniprot_id} [pdbid assembly]: {pdb_id} {assembly_id}", file=sys.stdout)
            residue_map = map_residues(uniprot_id, pdb_id)
            output = {
                "uniprot_id": uniprot_id,
                "pdb_id": pdb_id,
                "preferred_assembly": assembly_id,
                "residue_map": residue_map}
            print(json.dumps(output, indent=4), file=o)

    except ValueError as e:
        print(f"Error: {e}")

    # Close outfile
    o.close()

    # All done
    return

if __name__ == "__main__":
    main()

