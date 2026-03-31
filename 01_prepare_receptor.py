#!/usr/bin/env python3
"""
01_prepare_receptor.py — Download, fix, minimize, and prepare receptor structures.

Targets:
  - ROS1 G2032R (PDB 9QEK) — primary docking target
  - ROS1 WT     (PDB 7Z5X) — selectivity comparison
  - MET         (PDB 2WGJ) — dual-activity profiling

Outputs:
  - data/receptor_*.pdb, data/receptor_*.pdbqt
  - data/grid_*.txt (Vina grid box configs)
  - data/controls/zidesamtinib_crystal.pdb
"""

import os
import subprocess
import tempfile
from pathlib import Path

import gemmi
import numpy as np
import requests
from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import (
    NoCutoff,
    PDBFile,
    ForceField,
    Modeller,
    Simulation,
)
from openmm.unit import (
    nanometer,
    kelvin,
    picosecond,
    angstrom,
    kilocalorie_per_mole,
)

DATA_DIR = Path("data")
CONTROLS_DIR = DATA_DIR / "controls"

TARGETS = {
    "G2032R": {
        "pdb_id": "9QEK",
        "chain": "A",
        "catalytic_residues": [1980, 2027, 2102],  # K1980, E2027, D2102
    },
    "WT": {
        "pdb_id": "7Z5X",
        "chain": "A",
        "catalytic_residues": [1980, 2027, 2102],
    },
    "MET": {
        "pdb_id": "2WGJ",
        "chain": "A",
        "catalytic_residues": [1110, 1160, 1222],  # K1110, M1160, D1222
    },
}



def download_pdb(pdb_id: str) -> str:
    """Download PDB/mmCIF file from RCSB, return path. Converts CIF to PDB if needed."""
    out_path = DATA_DIR / f"{pdb_id}.pdb"
    if out_path.exists():
        print(f"  {pdb_id}.pdb already downloaded")
        return str(out_path)

    # Try PDB format first
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")
    resp = requests.get(url, timeout=30)

    if resp.status_code == 200:
        out_path.write_text(resp.text)
        print(f"  Saved {out_path} ({len(resp.text)} bytes)")
        return str(out_path)

    # Fall back to mmCIF and convert via PDBFixer
    print(f"  PDB format not available, downloading mmCIF...")
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    resp = requests.get(cif_url, timeout=30)
    resp.raise_for_status()

    cif_path = DATA_DIR / f"{pdb_id}.cif"
    cif_path.write_text(resp.text)
    print(f"  Saved {cif_path} ({len(resp.text)} bytes)")

    # Convert CIF to PDB using PDBFixer (reads CIF, writes PDB)
    fixer = PDBFixer(str(cif_path))
    with open(str(out_path), "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"  Converted to PDB: {out_path}")
    return str(out_path)


def extract_ligand(pdb_id: str, chain_id: str, output_path: str):
    """Extract co-crystallized ligand from structure for RMSD validation.
    Uses gemmi for robust CIF/PDB handling."""
    # Try CIF first (more reliable for newer structures), then PDB
    cif_path = DATA_DIR / f"{pdb_id}.cif"
    pdb_path = DATA_DIR / f"{pdb_id}.pdb"

    if cif_path.exists():
        structure = gemmi.read_structure(str(cif_path))
    else:
        structure = gemmi.read_structure(str(pdb_path))

    structure.remove_waters()

    # Standard amino acid names to exclude
    STANDARD_AA = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE",  # selenomethionine
    }
    SOLVENTS = {"HOH", "SO4", "GOL", "EDO", "PEG", "CL", "NA", "MG", "DMS", "ACT"}

    # Find ligands in the specified chain
    ligand_atoms = []
    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue
            for residue in chain:
                resname = residue.name.strip()
                if resname in STANDARD_AA or resname in SOLVENTS:
                    continue
                # Must be a non-standard residue (likely a ligand)
                if len(resname) > 3 or residue.het_flag != "\0":
                    print(f"  Found ligand: {resname} in chain {chain_id}")
                    ligand_atoms.append(residue)

    if not ligand_atoms:
        print(f"  WARNING: No ligands found in chain {chain_id}")
        return

    # Write ligand to PDB
    out_structure = gemmi.Structure()
    out_model = gemmi.Model("1")
    out_chain = gemmi.Chain(chain_id)
    for res in ligand_atoms:
        out_chain.add_residue(res)
    out_model.add_chain(out_chain)
    out_structure.add_model(out_model)
    out_structure.write_pdb(output_path)
    print(f"  Saved crystal ligand to {output_path}")


def extract_chain_and_fix(pdb_id: str, chain_id: str, output_path: str):
    """Extract chain, remove heterogens, fix missing atoms, add hydrogens, minimize."""
    # Step 1: Extract chain with gemmi (handles both PDB and CIF)
    cif_path = DATA_DIR / f"{pdb_id}.cif"
    pdb_path = DATA_DIR / f"{pdb_id}.pdb"

    if cif_path.exists():
        structure = gemmi.read_structure(str(cif_path))
    else:
        structure = gemmi.read_structure(str(pdb_path))

    # Keep only specified chain, remove heterogens and water
    structure.remove_waters()
    structure.remove_ligands_and_waters()

    # Write chain-only PDB for PDBFixer
    for model in structure:
        chains_to_remove = [c.name for c in model if c.name != chain_id]
        for cname in chains_to_remove:
            model.remove_chain(cname)

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as tmp:
        tmp_path = tmp.name
    structure.write_pdb(tmp_path)

    # Step 2: PDBFixer for missing residues/atoms and hydrogens
    fixer = PDBFixer(filename=tmp_path)
    os.unlink(tmp_path)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    # Step 3: Energy minimization with OpenMM
    print("  Running energy minimization (AMBER14, 1500 steps L-BFGS)...")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)

    # Add solvent isn't needed for vacuum minimization of a protein
    # Just use implicit solvent or vacuum
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )

    # Add position restraints on backbone heavy atoms
    from openmm import CustomExternalForce
    restraint = CustomExternalForce(
        "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    restraint.addGlobalParameter("k", 10.0 * kilocalorie_per_mole / angstrom**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    positions = modeller.positions
    for atom in modeller.topology.atoms():
        if atom.name in ("CA", "C", "N", "O") and atom.residue.name != "HOH":
            pos = positions[atom.index]
            restraint.addParticle(
                atom.index,
                [pos.x, pos.y, pos.z],
            )
    system.addForce(restraint)

    # Use a dummy integrator (we only minimize, don't run dynamics)
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1.0 / picosecond, 0.002 * picosecond)

    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # L-BFGS minimization
    simulation.minimizeEnergy(maxIterations=1500)

    # Get minimized positions
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    print(f"  Final energy: {state.getPotentialEnergy()}")

    # Save
    with open(output_path, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print(f"  Saved {output_path}")


def pdb_to_pdbqt(pdb_path: str, pdbqt_path: str):
    """Convert PDB to PDBQT using Meeko's mk_prepare_receptor."""
    # Strip .pdbqt extension for -o flag (meeko adds it)
    out_base = pdbqt_path.replace(".pdbqt", "")
    try:
        result = subprocess.run(
            [
                "mk_prepare_receptor.py",
                "--read_pdb", pdb_path,
                "-o", out_base,
                "-p",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode == 0:
            print(f"  Converted to PDBQT: {pdbqt_path}")
            return
        # If it fails due to a problematic residue, try deleting it
        stderr = result.stderr
        if "Expected" in stderr and "paddings" in stderr:
            # Extract residue from error message
            import re as _re
            match = _re.search(r"\(([A-Z]:\d+),", stderr)
            if match:
                bad_res = match.group(1)
                print(f"  Retrying with -d {bad_res}...")
                result = subprocess.run(
                    [
                        "mk_prepare_receptor.py",
                        "--read_pdb", pdb_path,
                        "-o", out_base,
                        "-p",
                        "-d", bad_res,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                if result.returncode == 0:
                    print(f"  Converted to PDBQT: {pdbqt_path}")
                    return
        print(f"  mk_prepare_receptor.py failed: {stderr[:200]}")
    except FileNotFoundError:
        print("  mk_prepare_receptor.py not found")



def calculate_grid_center(file_path: str, chain_id: str, residue_numbers: list) -> tuple:
    """Calculate grid box center as average CA position of catalytic residues.
    Works with PDB or CIF files via gemmi."""
    structure = gemmi.read_structure(file_path)

    ca_coords = []
    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue
            for residue in chain:
                seq_id = residue.seqid.num
                if seq_id in residue_numbers:
                    ca = residue.find_atom("CA", "*")
                    if ca:
                        ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])

    # If chain filter was too strict, try all chains
    if len(ca_coords) == 0:
        print(f"  WARNING: No catalytic residues in chain {chain_id}, scanning all chains...")
        for model in structure:
            for chain in model:
                for residue in chain:
                    seq_id = residue.seqid.num
                    if seq_id in residue_numbers:
                        ca = residue.find_atom("CA", "*")
                        if ca:
                            ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])

    if len(ca_coords) == 0:
        raise RuntimeError(
            f"Could not find any catalytic residues {residue_numbers} in {file_path}"
        )

    if len(ca_coords) != len(residue_numbers):
        print(f"  WARNING: Found {len(ca_coords)}/{len(residue_numbers)} catalytic residues")

    center = np.mean(ca_coords, axis=0)
    return tuple(center)


def write_grid_config(target_name: str, center: tuple, size: float = 25.0):
    """Write Vina grid box configuration file."""
    config_path = DATA_DIR / f"grid_{target_name}.txt"
    with open(config_path, "w") as f:
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size:.1f}\n")
        f.write(f"size_y = {size:.1f}\n")
        f.write(f"size_z = {size:.1f}\n")
    print(f"  Grid config: {config_path} (center: {center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    CONTROLS_DIR.mkdir(parents=True, exist_ok=True)

    # Extract zidesamtinib crystal pose from 9QEK before cleanup
    print("\n=== Extracting zidesamtinib crystal pose from 9QEK ===")
    download_pdb("9QEK")
    crystal_ligand_path = str(CONTROLS_DIR / "zidesamtinib_crystal.pdb")
    extract_ligand("9QEK", "A", crystal_ligand_path)

    for target_name, config in TARGETS.items():
        pdb_id = config["pdb_id"]
        print(f"\n=== Processing {target_name} ({pdb_id}) ===")

        # Download
        download_pdb(pdb_id)

        # Calculate grid center from the ORIGINAL structure
        # Use CIF if available (more reliable for newer entries)
        cif_path = DATA_DIR / f"{pdb_id}.cif"
        pdb_path = DATA_DIR / f"{pdb_id}.pdb"
        struct_path = str(cif_path) if cif_path.exists() else str(pdb_path)

        print("  Calculating grid center from original structure...")
        try:
            center = calculate_grid_center(
                struct_path, config["chain"], config["catalytic_residues"]
            )
        except RuntimeError as e:
            print(f"  ERROR: {e}")
            print("  Will retry after cleanup...")
            center = None

        # Extract chain, fix, minimize
        receptor_pdb = str(DATA_DIR / f"receptor_{target_name}.pdb")
        extract_chain_and_fix(pdb_id, config["chain"], receptor_pdb)

        # If grid center failed on original, try on cleaned structure
        if center is None:
            center = calculate_grid_center(
                receptor_pdb, "A", config["catalytic_residues"]
            )

        # Write grid config
        write_grid_config(target_name, center)

        # Convert to PDBQT
        receptor_pdbqt = str(DATA_DIR / f"receptor_{target_name}.pdbqt")
        pdb_to_pdbqt(receptor_pdb, receptor_pdbqt)

    print("\n=== Receptor preparation complete ===")
    print("Files created:")
    for f in sorted(DATA_DIR.glob("receptor_*")):
        print(f"  {f}")
    for f in sorted(DATA_DIR.glob("grid_*")):
        print(f"  {f}")


if __name__ == "__main__":
    main()
