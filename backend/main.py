"""
main.py  —  Το BACKEND της web εφαρμογής (FastAPI)
==================================================
Κάθε endpoint (/api/...) είναι μια "πόρτα": ο frontend στέλνει δεδομένα, το
backend τα δίνει στον πυρήνα (bioinf_core), και επιστρέφει JSON.

ΕΠΙΠΕΔΟ 1 (robustness): όρια μεγέθους, validation, per-session μοντέλο, error handling.
ΕΠΙΠΕΔΟ 2 (χρησιμότητα):
  - /api/fetch_ncbi : κατεβάζει πραγματική ακολουθία από το NCBI με ένα accession ID.
  - /api/mutations  : υποστηρίζει ΚΑΙ alignment για άνισες ακολουθίες (indels).
  - /api/history    : ιστορικό αναλύσεων από τη βάση (SQLite, db.py).
  - Η ανάλυση επιστρέφει και reverse complement + πρωτεΐνη.

Τρέξιμο (από τον φάκελο backend):
    pip install -r ../requirements.txt
    uvicorn main:app --reload   ->  http://localhost:8000
"""

import uuid
from pathlib import Path
from io import StringIO

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from Bio import SeqIO, Entrez

import bioinf_core as core
import db

app = FastAPI(title="Sequence Toolkit API")
db.init_db()   # σιγουρεύει ότι υπάρχει ο πίνακας ιστορικού

FRONTEND_DIR = Path(__file__).resolve().parent.parent / "frontend"
app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")
FRONTEND = FRONTEND_DIR / "index.html"

MAX_CHARS = 5_000_000
MAX_MODELS = 50
MODELS = {}


# ---------- Σχήματα εισόδου ----------
class FastaIn(BaseModel):
    fasta: str

class MutationIn(BaseModel):
    reference: str
    variant: str
    align: bool = False        # αν True (ή άνισα μήκη) -> alignment αντί για Hamming

class PredictIn(BaseModel):
    sequence: str
    session_id: str = ""

class NcbiIn(BaseModel):
    accession: str
    email: str = "peri.oly123@gmail.com"   # το NCBI απαιτεί email


def too_big(text):
    if len(text) > MAX_CHARS:
        return True, f"Το αρχείο είναι πολύ μεγάλο (όριο {MAX_CHARS/1_000_000:.0f} MB)."
    return False, ""


# ---------- Frontend ----------
@app.get("/", response_class=HTMLResponse)
def home():
    return FRONTEND.read_text(encoding="utf-8")


# ---------- PROJECT 1: ανάλυση ----------
@app.post("/api/analyze")
def analyze(data: FastaIn):
    big, msg = too_big(data.fasta)
    if big:
        return {"error": msg}
    try:
        records = core.validate_records(core.read_fasta_text(data.fasta))
    except ValueError as e:
        return {"error": str(e)}
    except Exception:
        return {"error": "Δεν μπόρεσα να διαβάσω το FASTA. Έλεγξε τη μορφή του."}

    report = core.analyze_sequences(records)
    first_id, first_seq = records[0]
    first_orfs = core.find_orfs(first_seq)

    db.save_analysis("analyze", f"{len(records)} ακολουθίες",
                     {"ids": [r[0] for r in records]})
    return {
        "report": report.to_dict(orient="records"),
        "first_id": first_id,
        "orfs": [{k: o[k] for k in ("frame", "start", "end", "length_aa")}
                 for o in first_orfs],
        # ΕΠΙΠΕΔΟ 2: extra ανάλυση για την 1η ακολουθία
        "reverse_complement": core.reverse_complement(first_seq)[:120],
        "protein": core.translate_seq(first_seq)[:120],
    }


# ---------- PROJECT 2: μεταλλάξεις (με προαιρετικό alignment) ----------
@app.post("/api/mutations")
def mutations(data: MutationIn):
    big, msg = too_big(data.reference + data.variant)
    if big:
        return {"error": msg}
    try:
        ref = core.validate_sequence(data.reference)
        var = core.validate_sequence(data.variant)
    except ValueError as e:
        return {"error": str(e)}

    # Αν ζητηθεί alignment ή τα μήκη διαφέρουν -> alignment mode (χειρίζεται indels)
    if data.align or len(ref) != len(var):
        res = core.align_and_diff(ref, var)
        db.save_analysis("align", f"{res['n_mismatch']} mismatch, {res['n_gap']} gap")
        return {"mode": "align", **res}

    # Αλλιώς κλασικό SNP report
    report = core.variant_report(ref, var)
    db.save_analysis("mutations", f"{core.hamming_distance(ref, var)} SNPs")
    return {
        "mode": "snp",
        "hamming": core.hamming_distance(ref, var),
        "snps": report.to_dict(orient="records"),
    }


# ---------- PROJECT 3: εκπαίδευση + πρόβλεψη ----------
@app.post("/api/train")
def train(data: FastaIn):
    big, msg = too_big(data.fasta)
    if big:
        return {"error": msg}
    try:
        recs, labels = [], []
        for r in SeqIO.parse(StringIO(data.fasta), "fasta"):
            lab = next((t.split("=")[1] for t in r.description.split()
                        if t.startswith("class=")), None)
            if lab is not None:
                recs.append((r.id, str(r.seq).upper()))
                labels.append(lab)
        recs = core.validate_records(recs)
    except ValueError as e:
        return {"error": str(e)}
    except Exception:
        return {"error": "Δεν μπόρεσα να διαβάσω το FASTA εκπαίδευσης."}

    if len(set(labels)) < 2:
        return {"error": "Χρειάζονται ≥2 κατηγορίες (κεφαλίδες με class=A, class=B, ...)."}

    model, cols = core.train_classifier(recs, labels)
    if len(MODELS) >= MAX_MODELS:
        MODELS.pop(next(iter(MODELS)))
    session_id = uuid.uuid4().hex
    MODELS[session_id] = (model, cols)
    return {"trained": True, "session_id": session_id,
            "n_sequences": len(recs), "classes": sorted(set(labels))}


@app.post("/api/predict")
def predict(data: PredictIn):
    entry = MODELS.get(data.session_id)
    if entry is None:
        return {"error": "Δεν βρέθηκε εκπαιδευμένο μοντέλο. Κάνε πρώτα εκπαίδευση (Βήμα 1)."}
    try:
        seq = core.validate_sequence(data.sequence)
    except ValueError as e:
        return {"error": str(e)}
    model, cols = entry
    label, proba = core.predict_sequence(model, cols, seq)
    db.save_analysis("predict", f"κατηγορία {label}")
    return {"label": str(label),
            "probabilities": {str(k): round(float(v), 3) for k, v in proba.items()}}


# ---------- ΕΠΙΠΕΔΟ 2: NCBI fetch ----------
@app.post("/api/fetch_ncbi")
def fetch_ncbi(data: NcbiIn):
    """Κατεβάζει μια ακολουθία FASTA από το NCBI με βάση το accession ID
    (π.χ. NM_000410). Χρειάζεται σύνδεση στο internet."""
    acc = data.accession.strip()
    if not acc:
        return {"error": "Δώσε ένα NCBI accession ID (π.χ. NM_000410)."}
    Entrez.email = data.email or "anonymous@example.com"
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()
    except Exception:
        return {"error": f"Δεν βρέθηκε ή δεν κατέβηκε το «{acc}». Έλεγξε το ID και τη σύνδεση."}
    if not fasta.strip().startswith(">"):
        return {"error": f"Το «{acc}» δεν επέστρεψε έγκυρο FASTA."}
    return {"fasta": fasta}


# ---------- ΕΠΙΠΕΔΟ 2: Ιστορικό ----------
@app.get("/api/history")
def history():
    return {"history": db.recent_history(20)}
