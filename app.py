import streamlit as st
import pandas as pd
import io
import base64
import matplotlib.pyplot as plt
from typing import List, Tuple

# Try to import Biopython SeqIO; fall back to simple parser
try:
    from Bio import SeqIO  # type: ignore
    BIOPYTHON_AVAILABLE = True
except Exception:
    BIOPYTHON_AVAILABLE = False

st.set_page_config(page_title="DNA GC Content Analyzer", layout="wide")

# ------------------ Utilities ------------------

def simple_fasta_parser(fasta_text: str) -> List[Tuple[str, str]]:
    """Parse FASTA from string. Returns list of (header, sequence)."""
    sequences = []
    header = None
    seq_lines = []
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                sequences.append((header, "".join(seq_lines).upper()))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        sequences.append((header, "".join(seq_lines).upper()))
    return sequences


def parse_fasta_file(uploaded_file) -> List[Tuple[str, str]]:
    content = uploaded_file.getvalue().decode("utf-8")
    if BIOPYTHON_AVAILABLE:
        # Use Biopython to ensure robust parsing
        records = list(SeqIO.parse(io.StringIO(content), "fasta"))
        return [(str(rec.id), str(rec.seq).upper()) for rec in records]
    else:
        return simple_fasta_parser(content)


def compute_gc(seq: str) -> float:
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    valid = sum(seq.count(n) for n in ["A", "T", "G", "C"])
    # Use length of sequence for denominator to include ambiguous bases; but also provide valid base count
    return 100.0 * (g + c) / len(seq) if len(seq) > 0 else 0.0


def sliding_window_gc(seq: str, window:int=100, step:int=10) -> Tuple[List[int], List[float]]:
    seq = seq.upper()
    gc_values = []
    positions = []
    if window <= 0:
        return positions, gc_values
    for start in range(0, max(len(seq)-window+1, 0), step):
        window_seq = seq[start:start+window]
        gc = compute_gc(window_seq)
        positions.append(start+1)  # 1-based position
        gc_values.append(gc)
    # Handle case where sequence shorter than window -> compute single window over full seq
    if not positions and len(seq) > 0:
        positions = [1]
        gc_values = [compute_gc(seq)]
    return positions, gc_values


def df_from_sequences(seqs: List[Tuple[str, str]]) -> pd.DataFrame:
    rows = []
    for header, seq in seqs:
        rows.append({
            "header": header,
            "length": len(seq),
            "gc_percent": round(compute_gc(seq), 3),
            "a_count": seq.count("A"),
            "t_count": seq.count("T"),
            "g_count": seq.count("G"),
            "c_count": seq.count("C"),
        })
    return pd.DataFrame(rows)


def get_table_download_link(df: pd.DataFrame, filename: str = "gc_results.csv") -> str:
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f"data:text/csv;base64,{b64}"
    return href

# ------------------ Streamlit UI ------------------

st.title("🧬 DNA GC Content Analyzer")
st.markdown(
    """
    Upload a FASTA file or paste DNA sequences to compute GC content per sequence and sliding-window GC profile.
    **Tips:**
    - For larger genomes, increase the window size.
    - Use the download button to save results as CSV for downstream analyses.
    """
)

col1, col2 = st.columns([2,1])

with col1:
    uploaded_file = st.file_uploader("Upload FASTA file", type=["fa","fasta","txt"]) 
    seq_text = st.text_area("Or paste FASTA / raw sequence(s) here (optional)", height=150)
    st.markdown("---")
    st.header("Sliding Window Settings")
    window = st.number_input("Window size (bp)", min_value=1, value=100, step=1)
    step = st.number_input("Step size (bp)", min_value=1, value=10, step=1)
    st.markdown("---")
    run_button = st.button("Analyze")

with col2:
    st.info("If you don't have Biopython installed, the app uses a built-in FASTA parser.")
    if BIOPYTHON_AVAILABLE:
        st.success("Biopython detected: SeqIO parsing enabled.")
    else:
        st.warning("Biopython not found: using a simple FASTA parser. Install biopython for robust parsing.")
    st.markdown("**Quick actions**")
    st.write("- Try the example sequences below to test the app.")
    if st.button("Load Example FASTA"):
        example = ">seq1\nATGCGCGATATTTGGCGCGCGAATTT\n>seq2\nATATATATATATATATATATATATATATATATATATATATAT\n>seq3\nGGGGCCCCAAAATTTTNNNNNNNNNN"
        seq_text = example
        st.experimental_set_query_params(example_loaded="1")

# Run analysis
if run_button:
    sequences = []
    if uploaded_file is not None:
        try:
            sequences = parse_fasta_file(uploaded_file)
        except Exception as e:
            st.error(f"Failed to parse uploaded file: {e}")
    if seq_text and not sequences:
        # If user pasted text, parse it
        sequences = simple_fasta_parser(seq_text)
        # If still empty but raw sequence (no header), accept it
        if not sequences and seq_text.strip():
            seq = seq_text.strip().replace("\n", "").upper()
            sequences = [("input_sequence", seq)]

    if not sequences:
        st.error("No sequences found. Please upload a FASTA file or paste sequence(s).")
    else:
        df = df_from_sequences(sequences)
        tab1, tab2 = st.tabs(["📊 Summary Table", "📈 Sliding Window Plots"])

        with tab1:
            st.subheader("Per-sequence Summary")
            st.dataframe(df.style.format({"gc_percent": "{:.3f}"}), use_container_width=True)
            csv_href = get_table_download_link(df, "gc_results.csv")
            st.markdown(f"[Download CSV]({csv_href})")

        with tab2:
            for header, seq in sequences:
                st.markdown(f"### {header} (len={len(seq)})")
                pos, gc_vals = sliding_window_gc(seq, window=window, step=step)
                if len(pos) == 0:
                    st.write("Sequence too short for the chosen window; showing whole-sequence GC.")
                    fig, ax = plt.subplots()
                    ax.bar([1], [compute_gc(seq)])
                    ax.set_xlabel("Window #")
                    ax.set_ylabel("GC%")
                    ax.set_ylim(0,100)
                    st.pyplot(fig)
                else:
                    fig, ax = plt.subplots(figsize=(8,3))
                    ax.plot(pos, gc_vals, marker="o")
                    ax.set_xlabel("Start position (1-based)")
                    ax.set_ylabel("GC%")
                    ax.set_title(f"Sliding-window GC (window={window}, step={step})")
                    ax.set_ylim(0,100)
                    st.pyplot(fig)

    st.success("Analysis complete.")


# Footer / About
st.markdown("---")
st.caption("App generated for educational/demo purposes. For heavy-duty genome analyses, consider using optimized bioinformatics pipelines and exact parsers (Biopython).")
