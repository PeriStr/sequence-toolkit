"""
bioinf_core.py
==============
Ο "πυρήνας" του εργαλείου: όλες οι συναρτήσεις ανάλυσης από τα 3 projects,
μαζεμένες σε ένα αρχείο ώστε να τις καλεί το app.py (το Streamlit dashboard).

Δεν υπάρχει τίποτα καινούριο εδώ σε σχέση με τα notebooks — απλώς οι ίδιες
συναρτήσεις, καθαρές και έτοιμες για επαναχρησιμοποίηση.

Ομάδες συναρτήσεων:
  1) PROJECT 1  — ανάλυση ακολουθίας (GC, βάσεις, ORFs)
  2) PROJECT 2  — ανίχνευση μεταλλάξεων (SNPs, transition/transversion, effect)
  3) PROJECT 3  — machine learning (features, εκπαίδευση, πρόβλεψη)
  4) I/O helpers — διάβασμα FASTA από κείμενο ή αρχείο
"""

from collections import Counter
from io import StringIO

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

STOP_CODONS = {"TAA", "TAG", "TGA"}
PURINES = set("AG")
PYRIMIDINES = set("CT")

# Έγκυροι χαρακτήρες DNA: A,T,G,C + N (άγνωστο) + U (RNA) + IUPAC ambiguity codes.
VALID_DNA = set("ACGTUNRYSWKMBDHV")


# ======================================================================
# 0) VALIDATION — έλεγχος & καθαρισμός εισόδου
# ======================================================================
def clean_sequence(seq):
    """Αφαιρεί κενά/newlines/αριθμούς και κάνει κεφαλαία. Επιστρέφει καθαρό string."""
    return "".join(ch for ch in seq.upper() if not ch.isspace() and not ch.isdigit())


def invalid_chars(seq):
    """Επιστρέφει ταξινομημένη λίστα με όσους χαρακτήρες ΔΕΝ είναι έγκυρο DNA."""
    return sorted(set(seq.upper()) - VALID_DNA)


def validate_sequence(seq, min_len=1):
    """Καθαρίζει & ελέγχει μία ακολουθία. Επιστρέφει καθαρή. Πετάει ValueError αν άκυρη."""
    s = clean_sequence(seq)
    if len(s) < min_len:
        raise ValueError("Η ακολουθία είναι κενή ή πολύ μικρή.")
    bad = invalid_chars(s)
    if bad:
        raise ValueError(f"Μη έγκυροι χαρακτήρες DNA: {', '.join(bad)}. "
                         "Επιτρέπονται μόνο A, T, G, C, N.")
    return s


def validate_records(records):
    """Ελέγχει λίστα (id, seq). Επιστρέφει καθαρή λίστα ή πετάει ValueError με σαφές μήνυμα."""
    if not records:
        raise ValueError("Δεν βρέθηκαν ακολουθίες. Είναι έγκυρο FASTA (γραμμές με '>' και ακολουθία);")
    cleaned = []
    for sid, seq in records:
        s = clean_sequence(seq)
        if not s:
            raise ValueError(f"Η ακολουθία «{sid}» είναι κενή.")
        bad = invalid_chars(s)
        if bad:
            raise ValueError(f"Η ακολουθία «{sid}» έχει μη έγκυρους χαρακτήρες: {', '.join(bad)}.")
        cleaned.append((sid, s))
    return cleaned


# ======================================================================
# 4) I/O — διάβασμα FASTA
# ======================================================================
def read_fasta_text(text):
    """Διαβάζει FASTA από ένα string (π.χ. ό,τι ανέβασε ο χρήστης).
    Επιστρέφει λίστα από (id, sequence)."""
    handle = StringIO(text)
    return [(rec.id, str(rec.seq).upper()) for rec in SeqIO.parse(handle, "fasta")]


def read_fasta_file(path):
    """Διαβάζει FASTA από αρχείο στο δίσκο. Επιστρέφει λίστα από (id, sequence)."""
    return [(rec.id, str(rec.seq).upper()) for rec in SeqIO.parse(path, "fasta")]


# ======================================================================
# 1) PROJECT 1 — ανάλυση ακολουθίας
# ======================================================================
def gc_fraction(seq):
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def base_counts(seq):
    seq = seq.upper()
    c = Counter(seq)
    out = {b.lower(): c.get(b, 0) for b in "ATGCN"}
    out["other"] = sum(v for k, v in c.items() if k not in "ATGCN")
    return out


def find_orfs(seq, min_aa_length=10):
    """Βρίσκει ORFs στα 3 εμπρός πλαίσια. Επιστρέφει λίστα από dicts."""
    seq = seq.upper().replace("U", "T")
    orfs = []
    for frame in (0, 1, 2):
        i = frame
        while i < len(seq) - 2:
            if seq[i:i+3] == "ATG":
                aas, j = [], i
                while j < len(seq) - 2:
                    c = seq[j:j+3]
                    if c in STOP_CODONS:
                        break
                    aas.append(str(Seq(c).translate()))
                    j += 3
                else:
                    i += 3
                    continue
                if len(aas) >= min_aa_length:
                    orfs.append({"frame": frame, "start": i, "end": j + 3,
                                 "length_aa": len(aas), "protein": "".join(aas)})
                i = j + 3
            else:
                i += 3
    return orfs


def find_orfs_6frames(seq, min_aa_length=10):
    """Επέκταση (stretch goal): ψάχνει και τους 2 κλώνους -> 6 πλαίσια."""
    seq = seq.upper().replace("U", "T")
    out = []
    for o in find_orfs(seq, min_aa_length):
        o["strand"] = "+"
        out.append(o)
    rc = str(Seq(seq).reverse_complement())
    for o in find_orfs(rc, min_aa_length):
        o["strand"] = "-"
        o["frame"] += 3
        out.append(o)
    return out


def longest_orf(seq, min_aa_length=10):
    orfs = find_orfs(seq, min_aa_length)
    return max(orfs, key=lambda o: o["length_aa"]) if orfs else None


def analyze_sequences(records):
    """Παίρνει λίστα (id, seq) και βγάζει τον πίνακα στατιστικών (DataFrame)."""
    rows = []
    for sid, seq in records:
        bc = base_counts(seq)
        lo = longest_orf(seq)
        rows.append({
            "id": sid,
            "length": len(seq),
            "gc_fraction": round(gc_fraction(seq), 3),
            **bc,
            "n_atg": seq.count("ATG"),
            "longest_orf_aa": lo["length_aa"] if lo else 0,
        })
    return pd.DataFrame(rows)


# ======================================================================
# 2) PROJECT 2 — ανίχνευση μεταλλάξεων
# ======================================================================
def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError(f"Different lengths: {len(s1)} vs {len(s2)}")
    return sum(a != b for a, b in zip(s1, s2))


def find_snps(reference, variant):
    if len(reference) != len(variant):
        raise ValueError("Sequences must be aligned (same length).")
    return [{"position": i, "ref": r, "alt": v}
            for i, (r, v) in enumerate(zip(reference.upper(), variant.upper()))
            if r != v]


def classify_substitution(ref_base, alt_base):
    if ref_base == alt_base:
        return "none"
    pair = {ref_base, alt_base}
    if pair <= PURINES or pair <= PYRIMIDINES:
        return "transition"
    return "transversion"


def classify_snp_effect(ref, alt, position, frame_offset=0):
    relative = position - frame_offset
    if relative < 0:
        return {"effect": "non-coding", "ref_codon": None, "alt_codon": None,
                "ref_aa": None, "alt_aa": None, "aa_position": None}
    codon_index = relative // 3
    pos_in_codon = relative % 3
    codon_start = frame_offset + codon_index * 3
    ref_codon = ref[codon_start:codon_start + 3]
    if len(ref_codon) < 3:
        return {"effect": "incomplete", "ref_codon": None, "alt_codon": None,
                "ref_aa": None, "alt_aa": None, "aa_position": None}
    alt_codon = list(ref_codon)
    alt_codon[pos_in_codon] = alt
    alt_codon = "".join(alt_codon)
    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())
    if ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == "*":
        effect = "nonsense"
    else:
        effect = "missense"
    return {"effect": effect, "ref_codon": ref_codon, "alt_codon": alt_codon,
            "ref_aa": ref_aa, "alt_aa": alt_aa, "aa_position": codon_index + 1}


def variant_report(reference, variant, frame_offset=0, var_id="variant"):
    rows = []
    for s in find_snps(reference, variant):
        eff = classify_snp_effect(reference, s["alt"], s["position"], frame_offset)
        rows.append({
            "variant": var_id,
            "position": s["position"],
            "ref_base": s["ref"],
            "alt_base": s["alt"],
            "ts_tv": classify_substitution(s["ref"], s["alt"]),
            "effect": eff["effect"],
            "ref_codon": eff["ref_codon"],
            "alt_codon": eff["alt_codon"],
            "ref_aa": eff["ref_aa"],
            "alt_aa": eff["alt_aa"],
            "aa_pos": eff["aa_position"],
        })
    return pd.DataFrame(rows)


def align_and_diff(seq1, seq2):
    """ΕΠΙΠΕΔΟ 2: ευθυγραμμίζει δύο ακολουθίες (ίσως ΑΝΙΣΟΥ μήκους) με Needleman-Wunsch
    και επιστρέφει τις διαφορές, ξεχωρίζοντας mismatches από gaps (indels).
    Χρησιμοποιείται όταν η Hamming distance δεν αρκεί (διαφορετικά μήκη)."""
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    aln = aligner.align(seq1.upper(), seq2.upper())[0]
    a1, a2 = str(aln[0]), str(aln[1])
    diffs, n_mismatch, n_gap = [], 0, 0
    for i, (x, y) in enumerate(zip(a1, a2)):
        if x == y:
            continue
        kind = "gap" if (x == "-" or y == "-") else "mismatch"
        if kind == "gap":
            n_gap += 1
        else:
            n_mismatch += 1
        diffs.append({"position": i, "ref": x, "alt": y, "kind": kind})
    return {"aligned_ref": a1, "aligned_var": a2,
            "n_mismatch": n_mismatch, "n_gap": n_gap, "diffs": diffs}


def reverse_complement(seq):
    """ΕΠΙΠΕΔΟ 2: ο ανάστροφος συμπληρωματικός κλώνος (ο 'άλλος' κλώνος του DNA)."""
    return str(Seq(clean_sequence(seq)).reverse_complement())


def translate_seq(seq, to_stop=True):
    """ΕΠΙΠΕΔΟ 2: μεταφράζει DNA -> πρωτεΐνη. to_stop=True σταματά στο 1ο stop codon."""
    return str(Seq(clean_sequence(seq)).translate(to_stop=to_stop))


# ======================================================================
# 3) PROJECT 3 — machine learning
# ======================================================================
def basic_features(seq):
    seq = seq.upper()
    n = len(seq)
    a, t, g, c = (seq.count(b) for b in "ATGC")
    return {"length": n, "gc": (g + c) / n if n else 0,
            "a_frac": a / n if n else 0, "t_frac": t / n if n else 0,
            "g_frac": g / n if n else 0, "c_frac": c / n if n else 0,
            "at_runs": seq.count("AAAA") + seq.count("TTTT"),
            "gc_runs": seq.count("GGGG") + seq.count("CCCC")}


def kmer_features(seq, k=3):
    seq = seq.upper()
    counts = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
    total = sum(counts.values())
    if total == 0:
        return {}
    return {f"k{k}_{km}": counts[km] / total for km in counts}


def build_feature_matrix(records, k=3):
    """Παίρνει λίστα (id, seq) και φτιάχνει τον πίνακα features X."""
    basic = pd.DataFrame([basic_features(s) for _, s in records])
    kmer = pd.DataFrame([kmer_features(s, k) for _, s in records]).fillna(0.0)
    return pd.concat([basic, kmer], axis=1)


def train_classifier(records, labels, k=3, seed=42):
    """Εκπαιδεύει Random Forest. Επιστρέφει (model, feature_columns)."""
    from sklearn.ensemble import RandomForestClassifier
    X = build_feature_matrix(records, k)
    model = RandomForestClassifier(n_estimators=200, random_state=seed)
    model.fit(X, labels)
    return model, list(X.columns)


def predict_sequence(model, feature_columns, seq, k=3):
    """Προβλέπει την κατηγορία μιας ΝΕΑΣ ακολουθίας. Επιστρέφει (label, prob_dict)."""
    feats = {}
    feats.update(basic_features(seq))
    feats.update(kmer_features(seq, k))
    row = pd.DataFrame([feats]).reindex(columns=feature_columns, fill_value=0.0)
    label = model.predict(row)[0]
    proba = dict(zip(model.classes_, model.predict_proba(row)[0]))
    return label, proba
