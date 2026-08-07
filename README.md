# 🧬 DNA GC Content Analyzer

A simple **Streamlit web app** to analyze DNA sequences for **GC content**.  
Upload FASTA files or paste raw DNA sequences and get **summary statistics** plus **sliding-window GC plots**.n

---

## 🚀 Features

- 📂 Upload FASTA files or paste raw sequences
    
- 🔬 Automatic parsing (uses **Biopython SeqIO** if available, else built-in parser)
    
- 📊 Per-sequence **summary table** (length, GC%, A/T/G/C counts)
    
- 📈 Interactive **sliding-window GC plots**
    
- 💾 Download results as CSV
    
- 🛠 Adjustable **window size** and **step size** for GC profiling  

---

## 📸 Demo

![GC Content Analyzer Screenshot](https://via.placeholder.com/800x400?text=Add+App+Screenshot+Here)

*(Replace with actual screenshot of your app!)*

---


📦 Requirements

Python 3.8+

Streamlit

Pandas

Matplotlib

(Optional) Biopython for robust FASTA parsing

Install dependencies with:

bash
Copy code

pip install streamlit pandas matplotlib biopython

📖 Usage
Upload a FASTA file (.fa, .fasta, .txt)

Or paste raw DNA/FASTA sequences into the text box

Adjust sliding window size & step

Click Analyze

View results in:

📊 Summary Table (GC%, length, base counts)

📈 Sliding Window Plots

Download results as CSV for further analysis

🧪 Example FASTA Input
shell
Copy code

>seq1
>
ATGCGCGATATTTGGCGCGCGAATTT

>seq2

ATATATATATATATATATATATATATATATATATATATATAT

>seq3

GGGGCCCCAAAATTTTNNNNNNNNNN

📂 Project Structure

bash

Copy code
├── app.py              # Streamlit app
├── requirements.txt    # Dependencies
├── README.md           # Documentation

✨ Future Improvements

Add AT/GC skew and CpG island detection

Option to handle very large genomes (chunked processing)

Integration with genomic file formats (e.g., GenBank, GFF)

📜 License

This project is for educational and demo purposes.
Feel free to fork and modify for your use.

🙌 Acknowledgments

Biopython for sequence parsing

Streamlit for quick web apps

## ⚙️ Installation

 Clone this repository:

```bash
git clone https://github.com/Bano733-code/gc-content-analyzer.git
cd gc-content-analyzer

Install required packages:
bash
Copy code
pip install -r requirements.txt
Run the app:
bash
Copy code
streamlit run app.py

yaml
Copy code
