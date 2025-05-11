import streamlit as st
import py3Dmol
import requests
from Bio.PDB import PDBParser
from io import StringIO

# Streamlit page settings
st.set_page_config(page_title="3D Protein Viewer", layout="wide")
st.title("🧬 3D Protein Structure Visualization")

# 🔎 Input PDB ID
pdb_id = st.text_input("Enter PDB ID (e.g., 1CRN, 6LU7):", "1CRN").upper()

# 📥 Function to fetch PDB file from RCSB
def fetch_pdb(pdb_id):
    """
    Fetches a PDB file from the RCSB PDB database using the PDB ID.
    
    Args:
        pdb_id (str): Protein PDB identifier.
        
    Returns:
        str: The PDB file content as text.
        None: If fetching failed.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        st.error(f"Failed to fetch data for ID: {pdb_id}")
        return None

# 🎨 Visualization settings in the sidebar
with st.sidebar:
    st.header("⚙️ Visualization Settings")
    style_option = st.selectbox("Visualization Style", ["cartoon", "stick", "surface"])
    color_option = st.selectbox("Color Scheme", ["spectrum", "chain", "residue"])
    show_ligands = st.checkbox("Show Ligands", value=True)

# 🧠 Function to generate protein visualization
def show_pdb(pdb_string, style="cartoon", color="spectrum", ligands=True):
    """
    Generates a 3D visualization of a protein structure from a PDB file.
    
    Args:
        pdb_string (str): Text representation of the PDB structure.
        style (str): Visualization style ("cartoon", "stick", "surface").
        color (str): Color scheme for the structure ("spectrum", "chain", "residue").
        ligands (bool): Flag to display ligands.
        
    Returns:
        view: A py3Dmol object containing the visualization.
    """
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_string, "pdb")
    view.setStyle({style: {"color": color}})
    if ligands:
        view.addStyle({"hetflag": True}, {"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    return view

# 🧬 Function to extract amino acid sequences from a PDB file
def extract_sequence(pdb_string):
    """
    Extracts the amino acid sequence from a PDB file.
    
    Args:
        pdb_string (str): Text representation of the PDB structure.
        
    Returns:
        dict: A dictionary containing sequences for each protein chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", StringIO(pdb_string))
    seqs = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            seq = ""
            for res in chain:
                if res.id[0] == " ":
                    resname = res.resname
                    one_letter = three_to_one.get(resname, "X")  # Convert 3-letter code to 1-letter code
                    seq += one_letter
            seqs[chain_id] = seq
    return seqs

# 🔄 Conversion of 3-letter code to 1-letter code (for protein sequences)
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O'
}

# 📊 Display the dashboard
pdb_data = fetch_pdb(pdb_id)
if pdb_data:
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader(f"🔬 Structure: {pdb_id}")
        st.markdown(f"[🔗 View on RCSB PDB](https://www.rcsb.org/structure/{pdb_id})")
        st.code(pdb_data[:300] + "\n...")

        # ➕ Amino acid sequence
        st.subheader("🧾 Amino Acid Sequence")
        sequences = extract_sequence(pdb_data)
        if sequences:
            for chain_id, seq in sequences.items():
                st.text_area(f"Chain {chain_id}", value=seq, height=150)
        else:
            st.write("No sequence available for display.")

    with col2:
        view = show_pdb(pdb_data, style_option, color_option, show_ligands)
        st.components.v1.html(view._make_html(), height=600, width=800)
