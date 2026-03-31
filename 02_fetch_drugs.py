#!/usr/bin/env python3
"""
02_fetch_drugs.py — Load FDA-approved drugs from Selleckchem SDFs, prepare PDBQT files.

Sources:
  - Selleckchem L1300 FDA-approved Drug Library (~3,100 drugs)
  - Selleckchem L8000 FDA-approved Anticancer Drug Library (~1,700 drugs)
  - Positive control SMILES for known ROS1 TKIs

Outputs:
  - data/ligands/*.pdbqt — one PDBQT per drug
  - data/controls/*.pdbqt — known ROS1 TKI controls
  - data/drug_library.csv — drug metadata, properties, CNS MPO
"""

import re
import sys
import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Crippen,
    Lipinski,
    SaltRemover,
    inchi,
)
from meeko import MoleculePreparation, PDBQTWriterLegacy
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("data")
LIGANDS_DIR = DATA_DIR / "ligands"
CONTROLS_DIR = DATA_DIR / "controls"

# Selleckchem SDF files
SDF_FILES = [
    DATA_DIR / "20260107-L1300-FDA-approved-Drug-Library.SDF",
    DATA_DIR / "20251113-L8000-FDA-approved-Anticancer-Drug-Library.SDF",
]

# Known ROS1 TKIs — used as positive controls
KNOWN_ROS1_TKIS = {
    "crizotinib": "Clc1cc(OC2(N)CC2)c(Cl)cc1-c1ccc2[nH]nc(-c3ccc(F)cc3F)c2n1",
    "entrectinib": "O=C(NC1CC1)c1cccc(-n2ncc3cc(-c4cn(C5CCNCC5)nc4-c4cccc(F)c4)ccc32)c1",
    "repotrectinib": "CC1(c2cc(-c3cncc4[nH]nc(N)c34)ccn2)C(F)(F)C1",
    "taletrectinib": "CC1(O)CN(c2cc3c(Nc4ccc5[nH]ncc5c4)ncnc3s2)C(=O)C1(F)F",
    "lorlatinib": "O=C(Nc1cc(-c2cccc(F)c2)cn1C1CCC(NC(=O)C2(F)CC2)CC1)C1CC1",
    "ceritinib": "Cc1cc(Nc2ncc(Cl)c(Nc3ccccc3S(=O)(=O)C(C)C)n2)c(OC(C)C)cc1C1CCNCC1",
}

# Zidesamtinib SMILES from PubChem CID 137321816
ZIDESAMTINIB_SMILES = "CC(C)(O)c1cc(-c2ccc3[nH]nc(-c4cc(F)c(OC5CC(N)C5)c(F)c4)c3n2)ccn1"


def parse_selleckchem_sdf(sdf_path: Path, source_label: str) -> list[dict]:
    """Parse a Selleckchem SDF file. Returns list of drug dicts."""
    drugs = []
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)

    for mol in tqdm(supplier, desc=f"Parsing {sdf_path.name}"):
        if mol is None:
            continue

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        # Read properties individually to handle non-UTF8 metadata
        props = {}
        for pname in ["Name", "_Name", "Cat", "Target", "CAS", "Synonyms"]:
            try:
                if mol.HasProp(pname):
                    props[pname] = mol.GetProp(pname)
            except (UnicodeDecodeError, RuntimeError):
                pass

        name = props.get("Name", props.get("_Name", f"drug_{len(drugs)}"))
        cat = props.get("Cat", "")
        target = props.get("Target", "")
        cas = props.get("CAS", "")
        synonyms = props.get("Synonyms", "")

        smiles = Chem.MolToSmiles(mol)
        drugs.append({
            "drug_name": name,
            "catalog_id": cat,
            "smiles": smiles,
            "mol": mol,
            "source": source_label,
            "target": target,
            "cas": cas,
            "synonyms": synonyms,
        })

    log.info("Parsed %d drugs from %s", len(drugs), sdf_path.name)
    return drugs


def deduplicate_drugs(drugs: list[dict]) -> list[dict]:
    """Deduplicate drugs by InChIKey, preferring L1300 (broader library) entries."""
    seen = {}
    for drug in drugs:
        mol = drug["mol"]
        try:
            ik = inchi.MolToInchi(mol)
            if ik:
                ik_key = inchi.InchiToInchiKey(ik)
            else:
                ik_key = Chem.MolToSmiles(mol)
        except Exception:
            ik_key = drug["smiles"]

        drug["inchikey"] = ik_key

        if ik_key not in seen:
            seen[ik_key] = drug
        elif drug["source"] == "selleckchem_L1300" and seen[ik_key]["source"] != "selleckchem_L1300":
            # Merge target info from L8000 (anticancer) before overwriting
            seen[ik_key]["target"] = seen[ik_key].get("target", "") or drug.get("target", "")
            seen[ik_key] = {**drug, "target": seen[ik_key]["target"] or drug.get("target", "")}

    result = list(seen.values())
    log.info("Deduplicated: %d -> %d unique drugs", len(drugs), len(result))
    return result


def filter_drugs(drugs: list[dict]) -> list[dict]:
    """Filter out biologics, bad molecules, strip salts."""
    remover = SaltRemover.SaltRemover()
    filtered = []

    for drug in drugs:
        mol = drug["mol"]

        # Strip salts
        try:
            mol = remover.StripMol(mol)
            if mol is None or mol.GetNumAtoms() == 0:
                continue
        except Exception:
            pass

        # Skip molecules that fail sanitization
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            log.debug("Skipping %s: sanitization failed", drug["drug_name"])
            continue

        # Remove MW > 1000 Da (biologics/peptides)
        mw = Descriptors.MolWt(mol)
        if mw > 1000:
            log.debug("Skipping %s: MW %.0f > 1000", drug["drug_name"], mw)
            continue

        drug["mol"] = mol
        drug["smiles"] = Chem.MolToSmiles(mol)
        filtered.append(drug)

    log.info("Filtered: %d -> %d drugs (removed biologics, salts, bad molecules)",
             len(drugs), len(filtered))
    return filtered


def tag_known_tkis(drugs: list[dict]) -> list[dict]:
    """Tag known ROS1 TKIs by name/synonym matching."""
    tki_names_lower = {name.lower() for name in KNOWN_ROS1_TKIS}
    tki_names_lower.add("zidesamtinib")

    for drug in drugs:
        name_lower = drug["drug_name"].lower()
        synonyms_lower = drug.get("synonyms", "").lower()
        target_lower = drug.get("target", "").lower()

        # Match by name or synonyms
        name_match = any(tki in name_lower or tki in synonyms_lower
                         for tki in tki_names_lower)
        # Also check if the target field mentions ROS1
        target_match = "ros1" in target_lower

        drug["is_known_ros1_tki"] = name_match

    tagged = [d["drug_name"] for d in drugs if d["is_known_ros1_tki"]]
    log.info("Tagged %d known ROS1 TKIs: %s", len(tagged), tagged)
    return drugs


def add_control_smiles(drugs: list[dict]) -> list[dict]:
    """Add positive control TKIs if not already in the library."""
    existing_smiles = {d["smiles"] for d in drugs}

    controls = {**KNOWN_ROS1_TKIS, "zidesamtinib": ZIDESAMTINIB_SMILES}

    for name, smiles in controls.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log.warning("Failed to parse control SMILES for %s", name)
            continue

        canonical = Chem.MolToSmiles(mol)
        if canonical in existing_smiles:
            log.info("Control %s already in library", name)
            continue

        ik = inchi.MolToInchi(mol)
        ik_key = inchi.InchiToInchiKey(ik) if ik else canonical

        drugs.append({
            "drug_name": name,
            "catalog_id": "",
            "smiles": canonical,
            "mol": mol,
            "inchikey": ik_key,
            "source": "control",
            "target": "ROS1",
            "cas": "",
            "synonyms": "",
            "is_known_ros1_tki": True,
        })
        log.info("Added control: %s", name)

    return drugs


def generate_conformer(mol):
    """Generate 3D conformer with ETKDG v3 + MMFF94 optimization. Returns mol or None."""
    if mol.GetNumAtoms() == 0:
        return None

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        # Try with more iterations and random coords
        params.maxIterations = 500
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            return None

    # MMFF94 optimization
    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=500)
    except Exception:
        pass  # keep unoptimized conformer

    return mol


def mol_to_pdbqt(mol, output_path: str) -> bool:
    """Convert RDKit mol to PDBQT using Meeko. Returns True on success."""
    try:
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_path, "w") as f:
                    f.write(pdbqt_string)
                return True
            else:
                log.debug("Meeko PDBQT error: %s", error_msg)

    except Exception as e:
        log.debug("Meeko conversion failed: %s", str(e)[:100])

    return False


def compute_properties(mol) -> dict:
    """Compute molecular properties for a drug."""
    mol_noH = Chem.RemoveHs(mol) if mol.GetNumAtoms() > 0 else mol

    mw = Descriptors.MolWt(mol_noH)
    tpsa = Descriptors.TPSA(mol_noH)
    clogp = Crippen.MolLogP(mol_noH)
    hbd = Lipinski.NumHDonors(mol_noH)
    hba = Lipinski.NumHAcceptors(mol_noH)
    rot_bonds = Lipinski.NumRotatableBonds(mol_noH)

    # Pfizer CNS MPO desirability functions (0-1 each)
    def d_mw(x):
        if x <= 360: return 1.0
        if x >= 500: return 0.0
        return 1.0 - (x - 360) / 140

    def d_clogp(x):
        if x <= 3.0: return 1.0
        if x >= 5.0: return 0.0
        return 1.0 - (x - 3.0) / 2.0

    def d_tpsa(x):
        if 40 <= x <= 90: return 1.0
        if x < 20 or x > 120: return 0.0
        if x < 40: return (x - 20) / 20
        return 1.0 - (x - 90) / 30

    def d_hbd(x):
        if x <= 0.5: return 1.0
        if x >= 3.5: return 0.0
        return 1.0 - (x - 0.5) / 3.0

    # CNS MPO score (4 of 6 parameters available, scaled to 0-6)
    cns_mpo_raw = d_mw(mw) + d_clogp(clogp) + d_tpsa(tpsa) + d_hbd(hbd)
    cns_mpo = cns_mpo_raw * 1.5

    # Simple CNS-penetrant flag
    cns_penetrant = tpsa < 79 and hbd <= 2 and mw < 450

    return {
        "mw": round(mw, 1),
        "tpsa": round(tpsa, 1),
        "clogp": round(clogp, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot_bonds,
        "cns_mpo": round(cns_mpo, 2),
        "cns_penetrant": cns_penetrant,
    }


def main():
    LIGANDS_DIR.mkdir(parents=True, exist_ok=True)
    CONTROLS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load drugs from Selleckchem SDFs
    drugs = []
    source_labels = {
        "L1300": "selleckchem_L1300",
        "L8000": "selleckchem_L8000",
    }

    for sdf_path in SDF_FILES:
        if not sdf_path.exists():
            log.warning("SDF not found: %s", sdf_path)
            continue
        label = "selleckchem_L8000" if "L8000" in sdf_path.name else "selleckchem_L1300"
        drugs.extend(parse_selleckchem_sdf(sdf_path, label))

    if not drugs:
        log.error("No drugs found! Place Selleckchem SDF files in data/")
        sys.exit(1)

    # Step 2: Deduplicate
    drugs = deduplicate_drugs(drugs)

    # Step 3: Filter
    drugs = filter_drugs(drugs)

    # Step 4: Tag known TKIs
    drugs = tag_known_tkis(drugs)

    # Step 5: Add control SMILES if not in library
    drugs = add_control_smiles(drugs)

    # Step 6: Generate 3D conformers, convert to PDBQT, compute properties
    log.info("Generating 3D conformers and PDBQT files...")
    results = []
    failed_embed = 0
    failed_pdbqt = 0

    for drug in tqdm(drugs, desc="Processing drugs"):
        mol = drug["mol"]
        name = drug["drug_name"]
        safe_name = re.sub(r'[^\w\-.]', '_', name)

        # Generate 3D conformer
        try:
            mol_3d = generate_conformer(mol)
        except Exception as e:
            log.debug("Conformer error for %s: %s", name, str(e)[:80])
            failed_embed += 1
            continue
        if mol_3d is None:
            log.debug("Failed to embed %s", name)
            failed_embed += 1
            continue

        # Determine output directory
        is_control = drug.get("is_known_ros1_tki", False)
        out_dir = CONTROLS_DIR if is_control else LIGANDS_DIR
        pdbqt_path = str(out_dir / f"{safe_name}.pdbqt")

        # Convert to PDBQT
        if not mol_to_pdbqt(mol_3d, pdbqt_path):
            log.debug("Failed PDBQT conversion for %s", name)
            failed_pdbqt += 1
            continue

        # Compute properties
        props = compute_properties(mol)

        results.append({
            "drug_name": name,
            "catalog_id": drug.get("catalog_id", ""),
            "smiles": drug["smiles"],
            "inchikey": drug.get("inchikey", ""),
            "source": drug["source"],
            "target": drug.get("target", ""),
            "cas": drug.get("cas", ""),
            "is_known_ros1_tki": is_control,
            "pdbqt_path": pdbqt_path,
            **props,
        })

    log.info(
        "Processed %d drugs: %d succeeded, %d failed embedding, %d failed PDBQT",
        len(drugs), len(results), failed_embed, failed_pdbqt,
    )

    # Step 7: Save CSV
    csv_path = DATA_DIR / "drug_library.csv"
    df = pd.DataFrame(results)
    df.to_csv(csv_path, index=False)
    log.info("Saved %s (%d drugs)", csv_path, len(df))

    # Summary
    n_controls = df["is_known_ros1_tki"].sum()
    n_cns = df["cns_penetrant"].sum()
    print(f"\n=== Drug Library Summary ===")
    print(f"Total drugs: {len(df)}")
    print(f"Known ROS1 TKIs: {n_controls}")
    print(f"CNS-penetrant: {n_cns}")
    print(f"Sources: {df['source'].value_counts().to_dict()}")
    print(f"Saved to: {csv_path}")


if __name__ == "__main__":
    main()
