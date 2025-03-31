"""
Utility functions for protein visualization that can be shared across modules
to avoid circular imports.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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

def visualize_grids(original_grid, blurred_grid, atom_types, opacity_settings=None):
    """
    Visualize the original and blurred grids for each atom type.
    
    Parameters:
    -----------
    original_grid : numpy.ndarray
        4D array with original atom positions.
    blurred_grid : numpy.ndarray
        4D array with blurred atom positions.
    atom_types : list
        List of atom type names.
    opacity_settings : dict, optional
        Settings to control opacity of visualization. Dictionary with keys:
        - base_opacity: Base opacity value (default 0.5)
        - min_opacity: Minimum opacity for low values (default 0.05)
        - opacity_scale: Scale factor for opacity (default 0.5)
    """
    # Set default opacity settings if none provided
    if opacity_settings is None:
        opacity_settings = {
            'base_opacity': 0.5,
            'min_opacity': 0.05,
            'opacity_scale': 0.5
        }
    
    n_atom_types = original_grid.shape[0]
    
    # Create grid coordinates - actual physical coordinates in Angstroms
    grid_size = 18  # 18 x 18 x 18 Angstrom cube
    grid_resolution = grid_size / (original_grid.shape[1] - 1)  # Calculate actual resolution
    x = np.linspace(0, grid_size, original_grid.shape[1])
    y = np.linspace(0, grid_size, original_grid.shape[2])
    z = np.linspace(0, grid_size, original_grid.shape[3])
    
    # Create a figure with two subplots (original and blurred)
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=["Original Atom Positions", "Gaussian Blur Applied"],
        specs=[[{'type': 'scene'}], [{'type': 'scene'}]],
        vertical_spacing=0.05
    )
    
    # Create grid lines - use fewer lines to reduce file size
    # Just create boundary lines and every 3 Angstroms
    grid_lines = []
    for i in [0, 3, 6, 9, 12, 15, 18]:  # Reduced grid lines (every 3 Å)
        # X-axis grid lines (reduced)
        for j in [0, 6, 12, 18]:
            grid_lines.append(
                go.Scatter3d(
                    x=[i, i], y=[j, j], z=[0, 18],
                    mode='lines',
                    line=dict(color='rgba(150,150,150,0.2)', width=1),
                    showlegend=False,
                    visible=True
                )
            )
        # Y-axis grid lines (reduced)
        for j in [0, 6, 12, 18]:
            grid_lines.append(
                go.Scatter3d(
                    x=[i, i], y=[0, 18], z=[j, j],
                    mode='lines',
                    line=dict(color='rgba(150,150,150,0.2)', width=1),
                    showlegend=False,
                    visible=True
                )
            )
    
    # Define color scheme for each atom type (standard chemical colors)
    atom_colors = {
        'C': 'gray',          # Carbon - gray
        'N': 'blue',          # Nitrogen - blue
        'O': 'red',           # Oxygen - red
        'H': 'white',         # Hydrogen - white
        'S': 'yellow',        # Sulfur - yellow
        'ligand': 'green',    # Ligand - green (arbitrary)
        'P': 'orange',        # Phosphorus - orange
        'Cl': 'green',        # Chlorine - green
        'F': 'green',         # Fluorine - green
        'Br': 'brown',        # Bromine - brown
        'I': 'purple',        # Iodine - purple
        'Fe': 'orange',       # Iron - orange-brown
        'Mg': 'green',        # Magnesium - green
        'Ca': 'gray',         # Calcium - gray
        'Zn': 'gray',         # Zinc - gray
    }
    
    # For any atom type not explicitly defined, use a default color
    default_color = 'purple'
    
    # Create traces for each atom type
    original_traces = []
    blurred_traces = []
    
    # Store all traces with their atom type for filtering
    all_atom_type_data = {}
    
    # For each atom type, create traces but set visibility to False initially
    for i in range(n_atom_types):
        atom_type = atom_types[i]
        traces_for_type = []
        
        # Get the color for this atom type
        atom_color = atom_colors.get(atom_type, default_color)
        
        # ORIGINAL GRID VISUALIZATION - Use points instead of cuboids
        non_zero_indices = np.where(original_grid[i] > 0)
        
        if len(non_zero_indices[0]) > 0:
            # Extract coordinates of non-zero points
            x_points = x[non_zero_indices[0]]
            y_points = y[non_zero_indices[1]]
            z_points = z[non_zero_indices[2]]
            
            # Create a scatter plot for atom points
            original_scatter = go.Scatter3d(
                x=x_points,
                y=y_points,
                z=z_points,
                mode='markers',
                marker=dict(
                    size=5,
                    color=atom_color,
                    opacity=0.8,
                    symbol='circle'
                ),
                name=f"{atom_type} atoms",
                visible=False  # Hide initially
            )
            
            original_traces.append(original_scatter)
            traces_for_type.append(original_scatter)
        
        # BLURRED GRID VISUALIZATION - Use 1x1x1 fixed size cubes
        blur_values = blurred_grid[i].flatten()
        max_val = np.max(blur_values)
        
        if max_val > 0:
            # Print some statistics to help debug
            print(f"Atom type: {atom_type}, Max value: {max_val:.6f}")
            percentiles = np.percentile(blur_values[blur_values > 0], [1, 5, 10, 25, 50, 75, 90, 95, 99])
            print(f"Non-zero percentiles (1,5,10,25,50,75,90,95,99): {percentiles}")
            
            # Use a higher threshold to reduce the number of voxels and file size
            # Use 10th percentile instead of 1st to significantly reduce number of cubes
            threshold = max(0.0001, percentiles[2])  # Use 10th percentile
            
            # Get voxel coordinates and values above threshold
            voxel_indices = np.where(blurred_grid[i] > threshold)
            
            # Subsample the voxels if there are too many to reduce file size
            max_voxels = 5000  # Maximum voxels to display per atom type
            
            if len(voxel_indices[0]) > max_voxels:
                # Randomly select a subset of voxels
                subset_idx = np.random.choice(len(voxel_indices[0]), max_voxels, replace=False)
                vx = x[voxel_indices[0][subset_idx]]
                vy = y[voxel_indices[1][subset_idx]]
                vz = z[voxel_indices[2][subset_idx]]
                values = blurred_grid[i][voxel_indices[0][subset_idx], voxel_indices[1][subset_idx], voxel_indices[2][subset_idx]] / max_val
                print(f"Subsampled {atom_type} from {len(voxel_indices[0])} to {max_voxels} voxels")
            else:
                vx = x[voxel_indices[0]]
                vy = y[voxel_indices[1]]
                vz = z[voxel_indices[2]]
                values = blurred_grid[i][voxel_indices] / max_val
            
            # Calculate transparency values based on voxel values
            # Use exponential function: opacity = min_opacity + (base_opacity - min_opacity) * (1 - exp(-opacity_scale * normalized_value))
            base_opacity = opacity_settings['base_opacity']
            min_opacity = opacity_settings['min_opacity']
            opacity_scale = opacity_settings['opacity_scale']
            
            # Normalize values to [0,1] range for consistent opacity scale
            normalized_values = values / values.max()  # Already normalized if values are from blurred_grid/max_val
            
            # Calculate opacity with exponential function that gives more weight to smaller values
            opacity_values = min_opacity + (base_opacity - min_opacity) * (1 - np.exp(-opacity_scale * normalized_values))
            
            # Create cuboids with fixed 1x1x1 size for blurred visualization
            vertices, faces, colors = create_cuboids(
                vx, vy, vz, 
                size=1.0,  # Fixed 1x1x1 size
                color=atom_color,  # Use atom-specific color
                opacity=opacity_values
            )
            
            # Add as mesh3d trace with initial visibility set to False
            blurred_mesh = go.Mesh3d(
                x=vertices[:, 0],
                y=vertices[:, 1],
                z=vertices[:, 2],
                i=faces[:, 0],
                j=faces[:, 1],
                k=faces[:, 2],
                facecolor=colors,
                opacity=1.0,  # Opacity is already in the colors
                name=f"{atom_type} blur",
                visible=False  # Hide initially
            )
            
            blurred_traces.append(blurred_mesh)
            traces_for_type.append(blurred_mesh)
        
        # Store traces for this atom type
        all_atom_type_data[atom_type] = traces_for_type
    
    # Add grid lines to both plots
    for line in grid_lines:
        fig.add_trace(line, row=1, col=1)
        
        # Create a new line for the second plot
        line2 = go.Scatter3d(
            x=line.x, y=line.y, z=line.z,
            mode='lines',
            line=dict(color='rgba(150,150,150,0.2)', width=1),
            showlegend=False,
            visible=True
        )
        fig.add_trace(line2, row=2, col=1)
    
    # Add all traces to the figure with initial visibility set to False
    for trace in original_traces:
        fig.add_trace(trace, row=1, col=1)
    
    for trace in blurred_traces:
        fig.add_trace(trace, row=2, col=1)
    
    # Create buttons for filtering atom types
    buttons = []
    
    # Add "All" button
    all_button = dict(
        label="All Atoms",
        method="update",
        args=[{"visible": [True] * len(fig.data)}]
    )
    buttons.append(all_button)
    
    # Add button for each atom type
    grid_line_count = len(grid_lines) * 2  # Grid lines for both plots
    
    for atom_type in atom_types:
        # Create visibility list for this atom type
        visibility = [False] * len(fig.data)
        
        # Always show grid lines
        for i in range(grid_line_count):
            visibility[i] = True
        
        # Find indices of traces for this atom type
        for trace_idx, trace in enumerate(fig.data):
            if trace.name and atom_type in trace.name:
                visibility[trace_idx] = True
        
        # Create button for this atom type
        button = dict(
            label=atom_type,
            method="update",
            args=[{"visible": visibility}]
        )
        buttons.append(button)
    
    # Add buttons to layout
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                active=0,
                x=0.05,
                y=1.15,
                buttons=buttons
            )
        ]
    )
    
    # Initially show the "All" view
    for trace in fig.data:
        trace.visible = True
    
    # Update layout to take full screen and center in HTML
    fig.update_layout(
        title_text="Atom Distributions: Original vs Gaussian Blur<br><sub>Use buttons to filter atom types</sub>",
        height=950,
        width=1200,  # Reduced width for better loading
        autosize=True,
        margin=dict(l=20, r=20, t=100, b=20),
        scene=dict(
            xaxis=dict(
                range=[0, 18], 
                title="X (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,  # 3 Å ticks - reduced for better performance
                gridcolor='rgba(200,200,200,0.4)',
            ),
            yaxis=dict(
                range=[0, 18], 
                title="Y (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,
                gridcolor='rgba(200,200,200,0.4)',
            ),
            zaxis=dict(
                range=[0, 18], 
                title="Z (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,
                gridcolor='rgba(200,200,200,0.4)',
            ),
            aspectmode="cube"
        ),
        scene2=dict(
            xaxis=dict(
                range=[0, 18], 
                title="X (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,
                gridcolor='rgba(200,200,200,0.4)',
            ),
            yaxis=dict(
                range=[0, 18], 
                title="Y (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,
                gridcolor='rgba(200,200,200,0.4)',
            ),
            zaxis=dict(
                range=[0, 18], 
                title="Z (Å)",
                tickmode='linear',
                tick0=0,
                dtick=3,
                gridcolor='rgba(200,200,200,0.4)',
            ),
            aspectmode="cube"
        ),
        # Center plot in HTML container
        template="plotly",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    
    return fig

def create_cuboids(x, y, z, size=1.0, color='red', opacity=None):
    """
    Create vertices and faces for rendering cuboids at specified positions.
    
    Parameters:
    -----------
    x, y, z : array-like
        Coordinates of cuboid centers
    size : float
        Size of cuboids (can be a scalar or vector for different dimensions)
    color : str
        Color of cuboids
    opacity : float or array-like
        Opacity of each cuboid (can be a single value or an array matching x,y,z)
    
    Returns:
    --------
    vertices : ndarray
        Array of vertex coordinates
    faces : ndarray
        Array of face indices
    colors : list
        List of face colors including opacity
    """
    if isinstance(size, (int, float)):
        size = np.array([size, size, size])
    
    # Define the 8 vertices of a unit cube
    unit_vertices = np.array([
        [-0.5, -0.5, -0.5],  # 0
        [0.5, -0.5, -0.5],   # 1
        [0.5, 0.5, -0.5],    # 2
        [-0.5, 0.5, -0.5],   # 3
        [-0.5, -0.5, 0.5],   # 4
        [0.5, -0.5, 0.5],    # 5
        [0.5, 0.5, 0.5],     # 6
        [-0.5, 0.5, 0.5]     # 7
    ])
    
    # Define the 12 triangular faces of a cube (two triangles per square face)
    unit_faces = np.array([
        [0, 1, 2], [0, 2, 3],  # Bottom face
        [4, 5, 6], [4, 6, 7],  # Top face
        [0, 1, 5], [0, 5, 4],  # Front face
        [2, 3, 7], [2, 7, 6],  # Back face
        [0, 3, 7], [0, 7, 4],  # Left face
        [1, 2, 6], [1, 6, 5]   # Right face
    ])
    
    num_cubes = len(x)
    num_vertices_per_cube = len(unit_vertices)
    num_faces_per_cube = len(unit_faces)
    
    # Prepare arrays for all vertices and faces
    all_vertices = np.zeros((num_cubes * num_vertices_per_cube, 3))
    all_faces = np.zeros((num_cubes * num_faces_per_cube, 3), dtype=int)
    
    # Generate colors with varying opacity if provided
    if opacity is None:
        opacity = np.ones(num_cubes)
    elif isinstance(opacity, (int, float)):
        opacity = np.ones(num_cubes) * opacity
    
    all_colors = []
    
    # Define RGB values for common colors
    color_map = {
        'red': (255, 0, 0),
        'green': (0, 255, 0),
        'blue': (0, 0, 255),
        'white': (255, 255, 255),
        'gray': (128, 128, 128),
        'yellow': (255, 255, 0),
        'orange': (255, 165, 0),
        'purple': (128, 0, 128),
        'brown': (165, 42, 42),
        'black': (0, 0, 0),
    }
    
    # Get the RGB values for the specified color
    if color in color_map:
        r, g, b = color_map[color]
    else:
        # Default to gray if color not found
        r, g, b = 128, 128, 128
    
    # For each cube
    for i in range(num_cubes):
        # Calculate the vertices for this cube
        start_idx = i * num_vertices_per_cube
        cube_vertices = unit_vertices * size + np.array([x[i], y[i], z[i]])
        all_vertices[start_idx:start_idx + num_vertices_per_cube] = cube_vertices
        
        # Calculate the face indices for this cube
        start_face_idx = i * num_faces_per_cube
        face_offset = i * num_vertices_per_cube
        cube_faces = unit_faces + face_offset
        all_faces[start_face_idx:start_face_idx + num_faces_per_cube] = cube_faces
        
        # Set color with opacity for this cube's faces
        cube_color = f'rgba({r},{g},{b},{opacity[i]:.3f})'
        
        # Assign the color to all faces of this cube
        all_colors.extend([cube_color] * num_faces_per_cube)
    
    return all_vertices, all_faces, all_colors 