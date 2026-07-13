# 🧬 Sequence Toolkit

**An end-to-end web application for DNA sequence analysis, mutation detection, and machine-learning classification — all in one place.**

🔗 **Live demo:** [https://sequence-toolkit.onrender.com](https://sequence-toolkit.onrender.com)

![Python](https://img.shields.io/badge/Python-3.11-blue)
![FastAPI](https://img.shields.io/badge/FastAPI-backend-009688)
![scikit-learn](https://img.shields.io/badge/scikit--learn-ML-F7931E)
![Biopython](https://img.shields.io/badge/Biopython-bioinformatics-2E7D32)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

> ⏱️ The live demo runs on a free tier that sleeps after inactivity — the first visit may take ~30 seconds to wake the server, then it responds instantly.

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Tech Stack](#tech-stack)
4. [Architecture](#architecture)
5. [Project Structure](#project-structure)
6. [Getting Started (Local)](#getting-started-local)
7. [API Reference](#api-reference)
8. [How It Works — The Science](#how-it-works--the-science)
9. [Deployment](#deployment)
10. [Limitations](#limitations)
11. [Roadmap](#roadmap)
12. [Acknowledgements](#acknowledgements)

---

## Overview

**Sequence Toolkit** turns three separate bioinformatics exercises into a single,
usable web tool. You paste or upload DNA sequences (in **FASTA** format), or fetch
them directly from the **NCBI** database, and the app returns statistics, detected
mutations, machine-learning predictions, and a searchable history — entirely through
the browser, with no code required.

It was built to demonstrate a complete, production-minded workflow: reusable analysis
code, a clean REST API, an interactive frontend, input validation, per-user isolation,
persistent storage, and a public deployment.

The analysis logic originates from a three-part bioinformatics curriculum:

| Project | Question it answers | Where it lives in the app |
|--------:|---------------------|---------------------------|
| **1 — DNA Analyzer** | *What is inside this sequence?* (length, GC content, base composition, ORFs) | Tab **Analysis** |
| **2 — Mutation Detection** | *What changed between two sequences, and does it matter?* (SNPs, transition/transversion, protein effect, alignment) | Tab **Mutations** |
| **3 — ML Classification** | *Which category does this sequence belong to?* (feature engineering + Random Forest) | Tab **Classification** |

---

## Features

### 🔬 Tab 1 — Sequence Analysis
- Upload a FASTA file (drag & drop) **or** paste text **or** fetch by **NCBI accession ID** (e.g. `NM_000410`).
- Per-sequence statistics table: length, GC fraction, base counts, number of `ATG`, longest ORF.
- **Interactive charts** (Chart.js): GC content per sequence and stacked base composition.
- Open Reading Frame (ORF) detection across the three forward frames.
- Extra outputs: **reverse complement** and **protein translation** of the first sequence.
- One-click CSV export of the results.

### 🧬 Tab 2 — Mutation Detection
- Compare a **reference** and a **variant** sequence.
- For equal-length sequences: **Hamming distance**, SNP table with **transition/transversion**
  classification and **synonymous / missense / nonsense** protein effect (with codon and amino-acid changes).
- For **unequal-length** sequences (indels): optional **Needleman–Wunsch alignment** that
  highlights mismatches (red) and gaps (yellow).
- Color-coded badges for every mutation type.

### 🤖 Tab 3 — Machine-Learning Classification
- Train a **Random Forest** on a labelled FASTA (headers contain `class=A`, `class=B`, …).
- Features combine hand-crafted statistics (GC, base fractions, homopolymer runs) with
  **k-mer frequencies**.
- Predict the class of any new sequence, with a probability bar for each category.
- Each user gets an isolated, per-session model — concurrent users never interfere.

### 🕘 Tab 4 — History
- Every analysis is persisted to a **SQLite** database and listed here, most-recent first.

### 🛡️ Robustness (built in)
- Input **validation and cleaning** (whitespace/case/digits stripped; only valid IUPAC DNA characters accepted).
- **Upload size limit** (5 MB) enforced on both client and server.
- Graceful, human-readable error messages instead of crashes.

---

## Tech Stack

| Layer | Technology |
|-------|-----------|
| **Backend** | Python 3.11, [FastAPI](https://fastapi.tiangolo.com/), Uvicorn |
| **Bioinformatics** | [Biopython](https://biopython.org/) (parsing, translation, alignment, NCBI Entrez) |
| **Machine Learning** | [scikit-learn](https://scikit-learn.org/) (Random Forest), pandas |
| **Database** | SQLite (standard library) |
| **Frontend** | Vanilla HTML / CSS / JavaScript, [Chart.js](https://www.chartjs.org/) |
| **Hosting** | [Render](https://render.com/) (free web service) |

No frontend build step, no npm — the interface is three static files served by the backend.

---

## Architecture

The app follows a strict **separation of concerns**, which keeps it easy to extend:

```
Browser (frontend)  ──HTTP/JSON──▶  FastAPI (backend)  ──▶  bioinf_core (analysis logic)
   index.html / style.css / app.js       main.py                gc, ORFs, SNPs, ML, alignment
                                            │
                                            ├──▶  db.py  (SQLite history)
                                            └──▶  Biopython Entrez (NCBI downloads)
```

- **`frontend/`** — the user interface. Sends requests with `fetch()`, renders tables/charts. Contains **no** analysis logic.
- **`backend/main.py`** — the API layer. Validates input, routes each request to the core, returns JSON. Contains **no** analysis logic itself.
- **`backend/bioinf_core.py`** — the pure analysis functions (identical to the original notebooks). The single source of truth for the science.
- **`backend/db.py`** — the persistence layer.

Change the science → edit the core. Change the look → edit the CSS. Change behaviour → edit the JS. Each concern is isolated.

---

## Project Structure

```
sequence_webapp/
├── backend/
│   ├── main.py          # FastAPI app: API endpoints + serves the frontend
│   ├── bioinf_core.py   # All analysis functions (Projects 1–3 + Level-2 extras)
│   ├── db.py            # SQLite helper (history)
│   └── history.db      # auto-generated at runtime (git-ignored)
├── frontend/
│   ├── index.html       # page structure
│   ├── style.css        # styling
│   └── app.js          # behaviour (API calls, charts, tabs)
├── requirements.txt     # Python dependencies
├── runtime.txt          # Python version for deployment
├── render.yaml          # Render deployment blueprint
├── .gitignore
└── README.md
```

---

## Getting Started (Local)

### Prerequisites
- Python 3.11+
- `pip`

### Installation & Run

```bash
# 1. clone the repository
git clone https://github.com/PeriStr/sequence-toolkit.git
cd sequence-toolkit

# 2. (recommended) create a virtual environment
python -m venv .venv
# Windows:  .venv\Scripts\activate
# macOS/Linux:  source .venv/bin/activate

# 3. install dependencies
pip install -r requirements.txt

# 4. start the server (from the backend folder)
cd backend
uvicorn main:app --reload
```

Then open **http://localhost:8000** in your browser.

Interactive API docs (auto-generated by FastAPI) are available at
**http://localhost:8000/docs**.

---

## API Reference

All endpoints accept/return JSON. The base URL is the site root.

| Method | Endpoint | Description | Body |
|-------:|----------|-------------|------|
| `GET`  | `/` | Serves the web interface | — |
| `POST` | `/api/analyze` | Sequence statistics, ORFs, reverse complement, protein | `{ "fasta": "..." }` |
| `POST` | `/api/mutations` | SNP report, or alignment for unequal lengths | `{ "reference": "...", "variant": "...", "align": false }` |
| `POST` | `/api/train` | Train a Random Forest, returns a `session_id` | `{ "fasta": "...(class= headers)" }` |
| `POST` | `/api/predict` | Predict a sequence's class using a trained model | `{ "sequence": "...", "session_id": "..." }` |
| `POST` | `/api/fetch_ncbi` | Download a FASTA record from NCBI | `{ "accession": "NM_000410" }` |
| `GET`  | `/api/history` | Recent analyses from the database | — |

**Example**

```bash
curl -X POST https://sequence-toolkit.onrender.com/api/analyze \
     -H "Content-Type: application/json" \
     -d '{"fasta": ">seq1\nATGCGTACGTTAGCATGCAA"}'
```

---

## How It Works — The Science

- **GC content** — the fraction of `G`+`C` bases; a basic fingerprint that varies between organisms and genome regions.
- **ORF (Open Reading Frame)** — a stretch from a start codon `ATG` to the next in-frame stop codon (`TAA`/`TAG`/`TGA`); a candidate gene. The app scans the three forward reading frames.
- **SNP classification** — a single-base change is a **transition** (purine↔purine or pyrimidine↔pyrimidine) or a **transversion** (across the two groups). Inside a coding region it is **synonymous** (same amino acid), **missense** (different amino acid), or **nonsense** (premature stop).
- **Alignment** — when two sequences differ in length (insertions/deletions), Hamming distance is meaningless, so a global **Needleman–Wunsch** alignment lines them up, distinguishing mismatches from gaps.
- **ML features** — each sequence is converted to a fixed-length numeric vector (statistics + k-mer frequencies) so a **Random Forest** can learn to separate classes.

---

## Deployment

The app deploys to **Render** via the included `render.yaml` blueprint:

- **Build:** `pip install -r requirements.txt`
- **Start:** `cd backend && uvicorn main:app --host 0.0.0.0 --port $PORT`

Any push to the `main` branch triggers an automatic rebuild and redeploy.

---

## Limitations

- **Free-tier sleep:** the hosted instance idles after ~15 minutes of inactivity; the next request wakes it (~30 s cold start).
- **Ephemeral history:** the SQLite database lives on Render's temporary disk and resets on each deploy. Persistent history would require a mounted disk or a managed Postgres database.
- **In-memory models:** trained ML models are held in memory per session, capped in number; they are not persisted across restarts.

These are acceptable trade-offs for a demonstration and are documented as future work below.

---

## Roadmap

- [ ] Persist trained models and history to a durable database (Postgres).
- [ ] Connect detected mutations to clinical significance (e.g. ClinVar lookups).
- [ ] Support 6-frame ORF search and richer visualizations in the UI.
- [ ] Add a deep-learning classifier (one-hot encoding + CNN) as an alternative model.
- [ ] User accounts so history and models are private and durable.

---

## Acknowledgements

Built as a self-initiated extension of a three-part bioinformatics curriculum
(DNA analysis, mutation detection, and ML sequence classification). The analysis
functions are powered by [Biopython](https://biopython.org/) and
[scikit-learn](https://scikit-learn.org/), served through
[FastAPI](https://fastapi.tiangolo.com/), and deployed on [Render](https://render.com/).

---

*Made with 🧬 and Python.*
