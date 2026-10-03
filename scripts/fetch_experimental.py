"""Fetch experimentally determined 3-D structures of flexible organic molecules.

Sources
-------
* RCSB PDB ligand instances (X-ray): best-resolution entry per chem-comp ID,
  coordinates from the ModelServer, template SMILES from the chem-comp
  dictionary, bond orders via ``AssignBondOrdersFromTemplate``.
* Crystallography Open Database (small-molecule crystal structures): the CIF
  is expanded to a 3x3x3 block of unit cells with gemmi, connectivity is
  inferred, and one complete molecule is mapped onto the template graph.
  This part needs gemmi, which is not a project dependency:

      uv run --with gemmi python scripts/fetch_experimental.py

Every structure is written as heavy-atom coordinates in the atom order of the
stored SMILES (``Chem.MolFromSmiles(smiles)``), which is the order used by
``TorsionSystem.heavy_idx``. Structures whose 3-D stereochemistry disagrees with
the template are dropped.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdForceFieldHelpers as ffh
from rdkit.Chem.MolStandardize import rdMolStandardize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from perlin_mole.chem import rotatable_torsions  # noqa: E402
from perlin_mole.experimental import canonical_with_coords, write_structure  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "experimental"

# Drug-like / organic PDB ligands (chem-comp ID -> short name)
PDB_LIGANDS = {
    "STI": "imatinib",
    "NIL": "nilotinib",
    "AQ4": "erlotinib",
    "IRE": "gefitinib",
    "1N1": "dasatinib",
    "B49": "sunitinib",
    "DB8": "bosutinib",
    "VIA": "sildenafil",
    "E20": "donepezil",
    "BRL": "rosiglitazone",
    "G39": "oseltamivir-carboxylate",
    "IBP": "ibuprofen",
    "FLP": "flurbiprofen",
    "IMN": "indomethacin",
    "CEL": "celecoxib",
    "BTN": "biotin",
    "RAL": "raloxifene",
    "DIF": "diclofenac",
    "MTX": "methotrexate",
    "TMQ": "trimetrexate",
    "0LI": "ponatinib",
}

# COD small-molecule crystal structures: name -> (search text, template SMILES). Entries are
# matched on exact molecular formula (excludes salts, co-crystals, solvates); the lowest R factor wins.
COD_TARGETS: dict[str, tuple[str, str]] = {
    "tolbutamide": ("tolbutamide", "CCCCNC(=O)NS(=O)(=O)c1ccc(C)cc1"),
    "mefenamic-acid": ("mefenamic", "Cc1cccc(Nc2ccccc2C(=O)O)c1C"),
    "aspirin": ("acetylsalicylic", "CC(=O)Oc1ccccc1C(=O)O"),
    "fenofibrate": ("fenofibrate", "CC(C)OC(=O)C(C)(C)Oc1ccc(cc1)C(=O)c1ccc(Cl)cc1"),
    "tamoxifen": ("tamoxifen", "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
    "ketoprofen": ("ketoprofen", "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1"),
    "ibuprofen-crystal": ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
}


def get(url: str, data: bytes | None = None, tries: int = 3) -> bytes:
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"} if data else {})
            with urllib.request.urlopen(req, timeout=60) as fh:
                return fh.read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))
    raise RuntimeError("unreachable")


def chemcomp_smiles(comp: str) -> str:
    d = json.loads(get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp}"))
    desc = {(x["type"], x["program"]): x["descriptor"] for x in d.get("pdbx_chem_comp_descriptor", [])}
    for key in (("SMILES_CANONICAL", "CACTVS"), ("SMILES_CANONICAL", "OpenEye OEToolkits"), ("SMILES", "CACTVS")):
        if key in desc and Chem.MolFromSmiles(desc[key]) is not None:
            return desc[key]
    raise ValueError(f"no usable SMILES for {comp}")


def best_entries(comp: str, n: int = 8) -> list[str]:
    q = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": comp,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "paginate": {"start": 0, "rows": n},
        },
    }
    d = json.loads(get("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q).encode()))
    return [r["identifier"] for r in d.get("result_set", [])]


def resolution(entry: str) -> float | None:
    d = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{entry}"))
    res = d.get("rcsb_entry_info", {}).get("resolution_combined")
    return float(res[0]) if res else None


def neutral_template(smiles: str) -> Chem.Mol | None:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    if any(a.GetFormalCharge() for a in m.GetAtoms()):
        m = rdMolStandardize.Uncharger().uncharge(m)
        Chem.SanitizeMol(m)
    if any(a.GetFormalCharge() for a in m.GetAtoms()):
        return None
    return m


def acceptable(template: Chem.Mol) -> tuple[bool, str]:
    if any(a.GetSymbol() not in ("C", "N", "O", "S", "F", "Cl", "Br", "P") for a in template.GetAtoms()):
        return False, "unsupported element"
    if len(Chem.GetMolFrags(template)) != 1:
        return False, "multi-fragment"
    k = len(rotatable_torsions(template))
    if not 3 <= k <= 10:
        return False, f"{k} torsions"
    if ffh.MMFFGetMoleculeProperties(Chem.AddHs(template)) is None:
        return False, "no MMFF94 parameters"
    return True, f"{k} torsions"


def pdb_instance(entry: str, comp: str, template: Chem.Mol) -> Chem.Mol | None:
    sdf = get(f"https://models.rcsb.org/v1/{entry.lower()}/ligand?label_comp_id={comp}&encoding=sdf").decode()
    m = Chem.MolFromMolBlock(sdf.split("$$$$")[0], sanitize=False, removeHs=False)
    if m is None:
        return None
    m = Chem.RemoveHs(m, sanitize=False)
    if m.GetNumAtoms() != template.GetNumAtoms():
        return None  # incomplete ligand or extra atoms
    xyz = m.GetConformer().GetPositions()
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1) + 10 * np.eye(len(xyz))
    if d.min() < 0.9:
        return None  # overlapping alternate locations
    try:
        for b in m.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
        for a in m.GetAtoms():
            a.SetFormalCharge(0)
            a.SetNoImplicit(False)
        m.UpdatePropertyCache(strict=False)
        return AllChem.AssignBondOrdersFromTemplate(template, m)
    except Exception:
        return None


def fetch_pdb(manifest: dict) -> None:
    for comp, name in PDB_LIGANDS.items():
        try:
            template = neutral_template(chemcomp_smiles(comp))
        except Exception as e:  # noqa: BLE001
            print(f"{comp:4s} {name:14s} skip: {e}")
            continue
        if template is None:
            print(f"{comp:4s} {name:14s} skip: charged/unparseable")
            continue
        ok, why = acceptable(template)
        if not ok:
            print(f"{comp:4s} {name:14s} skip: {why}")
            continue
        done = False
        for entry in best_entries(comp):
            res = resolution(entry)
            if res is None or res > 2.5:
                continue
            inst = pdb_instance(entry, comp, template)
            if inst is None:
                continue
            got = canonical_with_coords(template, inst)
            if got is None:
                print(f"{comp:4s} {name:14s} {entry}: stereo mismatch / mapping failed, trying next")
                continue
            smiles, mol = got
            path = write_structure(OUT / f"{name}.sdf", mol, name)
            manifest[name] = {
                "source": "PDB",
                "id": entry,
                "comp_id": comp,
                "resolution_A": res,
                "smiles": smiles,
                "n_torsions": len(rotatable_torsions(Chem.MolFromSmiles(smiles))),
                "n_heavy": mol.GetNumAtoms(),
                "file": path.name,
            }
            print(f"{comp:4s} {name:14s} {entry} {res:.2f} A  {why}")
            done = True
            break
        if not done:
            print(f"{comp:4s} {name:14s} no usable instance")


def _formula(text: str) -> str:
    return "".join(text.replace("-", " ").split())


def cod_candidates(text: str, formula: str) -> list[tuple[float, str]]:
    import urllib.parse

    d = json.loads(get(f"https://www.crystallography.net/cod/result?text={urllib.parse.quote(text)}&format=json"))
    out = []
    for x in d:
        r = x.get("Robs") or x.get("Rall")
        if _formula(x.get("formula") or "") == formula and r:
            out.append((float(r), str(x["file"])))
    return sorted(out)


def cod_instance(cif_text: str, template: Chem.Mol) -> tuple[Chem.Mol | None, str]:
    """Expand the cell, find one complete molecule, return it with template bond orders."""
    import gemmi
    from rdkit.Chem import rdDetermineBonds

    block = gemmi.cif.read_string(cif_text).sole_block()
    dis = [v for v in block.find_values("_atom_site_disorder_group") if v not in (".", "?")]
    occ = [v for v in block.find_values("_atom_site_occupancy") if v not in (".", "?")]
    if dis or any(abs(float(gemmi.cif.as_number(v)) - 1.0) > 1e-3 for v in occ):
        return None, "disordered"
    small = gemmi.make_small_structure_from_block(block)
    sites = [s for s in small.get_all_unit_cell_sites() if s.element.name != "H"]
    if not sites:
        return None, "no sites"
    rw = Chem.RWMol()
    pos = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                for s in sites:
                    f = s.fract
                    p = small.cell.orthogonalize(gemmi.Fractional(f.x + i, f.y + j, f.z + k))
                    rw.AddAtom(Chem.Atom(s.element.name.capitalize()))
                    pos.append((p.x, p.y, p.z))
    mol = rw.GetMol()
    from perlin_mole.experimental import attach_coords

    mol = attach_coords(mol, np.array(pos))
    rdDetermineBonds.DetermineConnectivity(mol, useHueckel=False)
    centre = small.cell.orthogonalize(gemmi.Fractional(0.5, 0.5, 0.5))
    centre = np.array([centre.x, centre.y, centre.z])
    xyz = np.array(pos)
    want = sorted(a.GetSymbol() for a in template.GetAtoms())
    frags = [f for f in Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False) if len(f) == template.GetNumAtoms()]
    frags = [f for f in frags if sorted(mol.GetAtomWithIdx(i).GetSymbol() for i in f) == want]
    if not frags:
        return None, "no complete molecule"
    frag = min(frags, key=lambda f: np.linalg.norm(xyz[list(f)].mean(0) - centre))
    idx = {a: n for n, a in enumerate(frag)}
    sub = Chem.RWMol()
    for a in frag:
        sub.AddAtom(Chem.Atom(mol.GetAtomWithIdx(a).GetSymbol()))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in idx and j in idx:
            sub.AddBond(idx[i], idx[j], Chem.BondType.SINGLE)
    sub = attach_coords(sub.GetMol(), xyz[list(frag)])
    for a in sub.GetAtoms():
        a.SetNoImplicit(False)
    sub.UpdatePropertyCache(strict=False)
    try:
        return AllChem.AssignBondOrdersFromTemplate(template, sub), "ok"
    except Exception as e:  # noqa: BLE001
        return None, f"template mismatch ({e.__class__.__name__})"


def fetch_cod(manifest: dict) -> None:
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    for name, (text, smi) in COD_TARGETS.items():
        template = Chem.MolFromSmiles(smi)
        ok, why = acceptable(template)
        if not ok:
            print(f"COD  {name:16s} skip: {why}")
            continue
        try:
            cands = cod_candidates(text, CalcMolFormula(template))
        except Exception as e:  # noqa: BLE001
            print(f"COD  {name:16s} search failed: {e}")
            continue
        done = False
        for r, cod_id in cands[:6]:
            try:
                cif = get(f"https://www.crystallography.net/cod/{cod_id}.cif").decode(errors="replace")
                inst, status = cod_instance(cif, template)
            except Exception as e:  # noqa: BLE001
                inst, status = None, f"error {e.__class__.__name__}: {e}"
            if inst is None:
                print(f"COD  {name:16s} {cod_id}: {status}")
                continue
            got = canonical_with_coords(template, inst)
            if got is None:
                print(f"COD  {name:16s} {cod_id}: stereo mismatch / mapping failed")
                continue
            smiles, mol = got
            path = write_structure(OUT / f"{name}.sdf", mol, name)
            manifest[name] = {
                "source": "COD",
                "id": cod_id,
                "r_factor": r,
                "resolution_A": None,
                "smiles": smiles,
                "n_torsions": len(rotatable_torsions(Chem.MolFromSmiles(smiles))),
                "n_heavy": mol.GetNumAtoms(),
                "file": path.name,
            }
            print(f"COD  {name:16s} {cod_id} R={r:.3f}  {why}")
            done = True
            break
        if not done:
            print(f"COD  {name:16s} no usable entry ({len(cands)} formula matches)")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_path = OUT / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if "--cod-only" not in sys.argv:
        fetch_pdb(manifest)
    if "--pdb-only" not in sys.argv:
        try:
            import gemmi  # noqa: F401
        except ImportError:
            print("COD fetch skipped: run with `uv run --with gemmi python scripts/fetch_experimental.py`")
        else:
            fetch_cod(manifest)
    manifest_path.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    print(f"{len(manifest)} structures in {manifest_path}")


if __name__ == "__main__":
    main()
