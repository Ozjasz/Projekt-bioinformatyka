import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gaussian_blur import optimize_for_large_data, visualize_grids
from vdw_radii import VDW_RADII

def create_protein_example(n_atoms=200):
    """
    Creates a more realistic protein structure example with
    atoms distributed according to typical protein geometry.
    
    Returns a MultiIndex DataFrame with atom coordinates.
    """
    # Define atom types with their typical proportions in proteins
    atom_types = {
        'C': 0.35,   # Carbon (backbone + side chains)
        'N': 0.10,   # Nitrogen (backbone + side chains)
        'O': 0.10,   # Oxygen (backbone + carbonyl)
        'H': 0.40,   # Hydrogen (abundant)
        'S': 0.02,   # Sulfur (cysteine, methionine)
        'ligand': 0.03  # Small molecule ligand
    }
    
    # Calculate number of atoms for each type
    atom_counts = {atom: int(np.ceil(prop * n_atoms)) for atom, prop in atom_types.items()}
    
    # Adjust to match total
    total = sum(atom_counts.values())
    if total > n_atoms:
        diff = total - n_atoms
        # Remove from the most abundant first
        for atom in sorted(atom_counts.keys(), key=lambda k: atom_counts[k], reverse=True):
            if atom_counts[atom] > diff:
                atom_counts[atom] -= diff
                break
            elif atom_counts[atom] > 1:
                diff -= (atom_counts[atom] - 1)
                atom_counts[atom] = 1
    
    # Create data structure
    data = []
    indices = []
    
    # Create a spherical protein structure with radius ~7Å
    # Center the protein at (9, 9, 9) to fit in 18x18x18 cube
    center = np.array([9.0, 9.0, 9.0])
    
    # Function to generate atom positions in a spherical shell
    def generate_positions(n, min_r, max_r, center, spread=0.8):
        positions = []
        for _ in range(n):
            # Random direction (uniform on sphere)
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.arccos(2*np.random.random() - 1)
            
            # Random radius (slightly concentrated toward outer regions)
            r = np.random.uniform(min_r, max_r)
            
            # Convert to Cartesian coordinates
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)
            
            # Add some noise/spread
            noise = np.random.normal(0, spread, 3)
            pos = center + np.array([x, y, z]) + noise
            
            # Ensure position is within bounds (0-18Å)
            pos = np.clip(pos, 0.1, 17.9)
            
            positions.append(pos)
        return positions
    
    # Generate main protein atoms (carbon backbone forms the core)
    atom_id = 0
    for atom_type, count in atom_counts.items():
        if atom_type == 'ligand':
            # Place ligand in center
            positions = generate_positions(count, 0, 2.5, center, spread=0.5)
        elif atom_type == 'C':
            # Carbon atoms form the backbone (middle layer)
            positions = generate_positions(count, 3.5, 6.5, center, spread=1.0)
        elif atom_type == 'N' or atom_type == 'O':
            # N and O atoms (slightly outer layer, in functional groups)
            positions = generate_positions(count, 4.0, 7.0, center, spread=1.2)
        elif atom_type == 'S':
            # Sulfur atoms (sparse, can be anywhere)
            positions = generate_positions(count, 3.0, 7.0, center, spread=1.5)
        else:  # H atoms
            # Hydrogen atoms (mostly on the outside surface)
            positions = generate_positions(count, 5.0, 7.5, center, spread=0.9)
        
        # Add to data structure
        for pos in positions:
            data.append({'x': pos[0], 'y': pos[1], 'z': pos[2]})
            indices.append((atom_type, atom_id))
            atom_id += 1
    
    # Create MultiIndex
    index = pd.MultiIndex.from_tuples(indices, names=['atom_type', 'atom_id'])
    
    # Create dataframe
    return pd.DataFrame(data, index=index)

if __name__ == "__main__":
    # Create realistic protein example
    atom_coords_df = create_protein_example(n_atoms=200)
    
    # Print summary stats
    print(f"Total atoms: {len(atom_coords_df)}")
    type_counts = atom_coords_df.index.get_level_values('atom_type').value_counts()
    print("Atom type distribution:")
    for atom_type, count in type_counts.items():
        print(f"  {atom_type}: {count}")
    
    # Use Van der Waals radii from imported dictionary
    atom_types = type_counts.index.tolist()
    vdw_subset = {atom_type: VDW_RADII.get(atom_type, 1.5) for atom_type in atom_types}
    
    # Apply Gaussian blur
    print("Applying Gaussian blur...")
    original_grid, blurred_grid = optimize_for_large_data(atom_coords_df, vdw_subset)
    
    # Visualize
    print("Generating visualization...")
    fig = visualize_grids(original_grid, blurred_grid, type_counts.index.tolist())
    fig.write_html("protein_visualization.html")
    print("Visualization saved to protein_visualization.html") 