"""
db.py  —  Απλή βάση δεδομένων (SQLite) για ΙΣΤΟΡΙΚΟ αναλύσεων.
================================================================
Με απλά λόγια:
  - Μέχρι τώρα ό,τι έκανε ο χρήστης χανόταν μόλις έκλεινε τη σελίδα.
  - Εδώ αποθηκεύουμε κάθε ανάλυση σε ένα αρχείο βάσης (history.db), ώστε να
    υπάρχει ιστορικό που μπορεί να δει ξανά.

Γιατί SQLite: έρχεται μαζί με την Python (μηδέν εγκατάσταση), αποθηκεύει τα
πάντα σε ΕΝΑ αρχείο, και είναι ό,τι πρέπει για μικρές/μεσαίες εφαρμογές.

Πίνακας `analyses`:
  id          — αύξων αριθμός
  created_at  — πότε έγινε
  kind        — τι είδους (analyze / mutations / align / predict)
  summary     — μια σύντομη περιγραφή (π.χ. "3 ακολουθίες")
  detail      — πλήρη δεδομένα σε JSON (για μελλοντική χρήση)
"""

import json
import sqlite3
from pathlib import Path

# Το αρχείο της βάσης μπαίνει δίπλα σε αυτό το script.
DB_PATH = Path(__file__).resolve().parent / "history.db"


def _connect():
    con = sqlite3.connect(DB_PATH)
    con.row_factory = sqlite3.Row     # ώστε να παίρνουμε γραμμές σαν dict
    return con


def init_db():
    """Δημιουργεί τον πίνακα αν δεν υπάρχει. Καλείται μία φορά στην εκκίνηση."""
    with _connect() as con:
        con.execute("""
            CREATE TABLE IF NOT EXISTS analyses (
                id         INTEGER PRIMARY KEY AUTOINCREMENT,
                created_at TEXT DEFAULT (datetime('now','localtime')),
                kind       TEXT NOT NULL,
                summary    TEXT NOT NULL,
                detail     TEXT
            )
        """)


def save_analysis(kind, summary, detail=None):
    """Αποθηκεύει μία εγγραφή. detail = οποιοδήποτε dict (μετατρέπεται σε JSON)."""
    with _connect() as con:
        con.execute(
            "INSERT INTO analyses (kind, summary, detail) VALUES (?, ?, ?)",
            (kind, summary, json.dumps(detail or {}, ensure_ascii=False)),
        )


def recent_history(limit=20):
    """Επιστρέφει τις πιο πρόσφατες εγγραφές (λίστα από dicts)."""
    with _connect() as con:
        rows = con.execute(
            "SELECT id, created_at, kind, summary FROM analyses "
            "ORDER BY id DESC LIMIT ?", (limit,)
        ).fetchall()
    return [dict(r) for r in rows]


def clear_history():
    """Διαγράφει όλο το ιστορικό (χρήσιμο για δοκιμές)."""
    with _connect() as con:
        con.execute("DELETE FROM analyses")
