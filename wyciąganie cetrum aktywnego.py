from Bio.PDB import PDBList, PDBParser, NeighborSearch
import plotly.graph_objects as go

def download_pdb_file(id):
    """Downloads pdb file and saves it locally.

    Parameters
    ----------
    id : string
        PDB id of the target protein.

    Returns
    -------
    string
        Path to the downloaded file.
    """
    pdbl = PDBList()
    pdb_file = pdbl.retrieve_pdb_file(id, file_format="pdb")
    return pdb_file

def extract_active_site(pdb_file, cutoff=5.0):
    """Extracting protein's active site and its ligand (copper ion) from the PDB file.

    Parameters
    ----------
    pdb_file : string
        Path to the target PDB file.
    cutoff : float, optional
        Here we define active site as residues within a certain distance
        from the ligand. The distance is measured in angstroms, by default 5.0.

    Returns
    -------
    tuple
        copper_atoms : list of ndarray
            List of coordinates for copper atoms.
        active_site_atoms : list of ndarray
        
    """

    parser = PDBParser()
    structure = parser.get_structure(pdb_file, pdb_file)
    atom_list = list(structure.get_atoms())
    neighbors = NeighborSearch(atom_list)

    copper_ions = [atom.coord for atom in atom_list if atom.element == "CU"]
    active_site_atoms = []

    for copper in copper_ions:
        nearby_atoms = neighbors.search(copper, cutoff)
        for atom in nearby_atoms:
            active_site_atoms.append(atom.coord)
    

    return copper_ions, active_site_atoms

def visualize_structure(copper_ions, protein_atoms):
    """ Visualizes the protein structure, copper ions, and active site residues using Plotly.


    Parameters
    ----------
    copper_ions : list of ndarray
        Coordinates of the copper ion.
    protein_atoms : list of ndarray
        Coordinates of protein atoms.
    """
    copper_x, copper_y, copper_z = zip(*copper_ions) if copper_ions else ([], [], [])
    
    protein_x, protein_y, protein_z = zip(*protein_atoms) if protein_atoms else ([], [], [])

    # Create copper ion scatter plot
    copper_trace = go.Scatter3d(
        x=copper_x, y=copper_y, z=copper_z,
        mode='markers',
        marker=dict(size=10, color='orange'),
        name='Copper Ions'
    )
    
    # Create active site scatter plot
    protein_trace = go.Scatter3d(
        x=protein_x, y=protein_y, z=protein_z,
        mode='markers',
        marker=dict(size=5, color='blue'),
        name='Protein Atoms'
    )

    # Combine traces
    fig = go.Figure(data=[copper_trace, protein_trace])

    # Set the layout
    fig.update_layout(
        title="Protein-ligand active site",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z"
        ),
        margin=dict(r=10, l=10, b=10, t=10)
    )

    fig.show()


pdb_file = download_pdb_file("1AAC")
copper, protein_atoms = extract_active_site(pdb_file, 5)
visualize_structure(copper, protein_atoms)
print(pdb_file)
