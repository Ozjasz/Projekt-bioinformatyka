# Protein Visualization with Gaussian Blur

A toolkit for visualizing protein structures using Gaussian blurring based on Van der Waals radii.

## Overview

This project provides tools to create interactive 3D visualizations of protein structures, applying Gaussian blur to represent the electron density around each atom. The visualization helps to better understand the spatial arrangement and types of atoms in protein structures.

## Project Structure

The codebase is organized into modular components:

- **gaussian_blur.py**: Main module with the core Gaussian blur algorithms
- **grid_utils.py**: Utility functions for grid creation and manipulation
- **atom_utils.py**: Functions for working with atom data
- **visualization_utils.py**: Tools for creating and saving visualizations
- **vdw_radii.py**: Standard Van der Waals radii values for different atom types
- **protein_utils.py**: Utilities for protein structure generation and visualization

## Key Features

- Generate 3D grid representations of protein structures
- Apply Gaussian blur based on Van der Waals radii
- Fast and optimized algorithm for large protein structures
- Interactive visualization with HTML output
- Beginner-friendly API with sensible defaults
- Customizable opacity settings for better visualization

## How to Use

For beginners, a simple way to create a visualization:

```python
from gaussian_blur import create_protein_visualization

# Assuming you have a DataFrame of atom coordinates with atom_type in the index
create_protein_visualization(atom_coords_df, output_filename="my_protein.html")
```

### Adjusting Opacity

You can make the blurred voxels less transparent for better visibility:

```python
from gaussian_blur import create_protein_visualization
from visualization_utils import adjust_opacity_settings

# Create settings with higher opacity
opacity_settings = adjust_opacity_settings(
    base_opacity=0.8,    # Higher values make voxels more opaque (0.0-1.0)
    min_opacity=0.3,     # Minimum opacity for low-value voxels
    opacity_scale=0.9    # How quickly opacity increases with value
)

# Use custom opacity settings
create_protein_visualization(
    atom_coords_df,
    output_filename="high_opacity_protein.html",
    opacity=opacity_settings
)
```

### Advanced Usage

For more advanced usage:

```python
from gaussian_blur import optimize_for_large_data
from atom_utils import get_atom_types, get_atom_radii
from visualization_utils import save_visualization, adjust_opacity_settings
from protein_utils import visualize_grids

# Get atom types and prepare radii dictionary
atom_types = get_atom_types(atom_coords_df)
vdw_radii_dict = get_atom_radii(atom_types)

# Apply Gaussian blur with optimized algorithm
original_grid, blurred_grid = optimize_for_large_data(atom_coords_df, vdw_radii_dict)

# Set custom opacity for better visibility
opacity = adjust_opacity_settings(base_opacity=0.8, opacity_scale=0.9)

# Visualize the results
fig = visualize_grids(original_grid, blurred_grid, atom_types, opacity_settings=opacity)

# Save to an HTML file
save_visualization(fig, "detailed_protein.html")
```

## Data Format

The expected input is a pandas DataFrame with:
- A MultiIndex containing 'atom_type' as the last level
- Columns for 'x', 'y', and 'z' coordinates

## Requirements

- NumPy
- pandas
- Plotly (for visualization) 