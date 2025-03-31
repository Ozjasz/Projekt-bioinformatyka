from protein_utils import create_protein_example, visualize_grids
from grid_utils import (
    create_grid,
    initialize_grids,
    coords_to_grid_indices,
    is_within_grid,
    get_neighbor_offsets,
    gaussian_3d_norm,
    gaussian_value,
)
from visualization_utils import save_visualization, adjust_opacity_settings
from atom_utils import get_atom_types, get_atom_radii, print_atom_stats

"""
Gaussian Blur Module for Protein Visualization

This module provides functions to apply Gaussian blur to protein atom coordinates,
creating smooth 3D representations of protein structures. The blur is based on 
the Van der Waals radii of different atom types.

Main functions:
- apply_gaussian_blur: Full implementation (slow but accurate)
- optimize_for_large_data: Fast implementation for larger datasets (recommended)
- create_protein_visualization: Simple wrapper for beginners

Example usage:
    from gaussian_blur import create_protein_visualization
    create_protein_visualization(atom_dataframe, output_filename="my_protein.html")
"""


def apply_gaussian_blur(atom_coords_df, vdw_radii_dict):
    """
    Apply Gaussian blur to atom coordinates to create a smooth 3D representation.

    This function takes atom coordinates and creates two 3D grid representations:
    1. An original grid where each atom is represented as a single point
    2. A blurred grid where each atom contributes a Gaussian distribution based on its Van der Waals radius

    Parameters:
    -----------
    atom_coords_df : pandas.DataFrame
        Dataframe with atom coordinates.
        Expected structure:
            - MultiIndex with the last level being 'atom_type'
            - Columns should include 'x', 'y', and 'z' coordinates

    vdw_radii_dict : dict
        Dictionary mapping atom types (str) to their Van der Waals radii (float).
        Example: {'C': 1.7, 'N': 1.55, 'O': 1.52}

    Returns:
    --------
    tuple of numpy.ndarray
        Two 4D arrays:
        - original_grid: Binary representation where atoms are single points (shape: [n_atom_types, grid_x, grid_y, grid_z])
        - blurred_grid: Smoothed representation with Gaussian blur (same shape)
    """
    # Get unique atom types from the dataframe
    atom_types = get_atom_types(atom_coords_df)

    # Set up the 3D grid
    x, y, z, n_points = create_grid(grid_size=18, grid_resolution=0.5)

    # Initialize empty grids
    original_grid, blurred_grid = initialize_grids(len(atom_types), n_points)

    # Process each atom type separately
    for i, atom_type in enumerate(atom_types):
        # Get coordinates for all atoms of current type
        type_coords = atom_coords_df.xs(atom_type, level="atom_type")

        # Get the Van der Waals radius for this atom type
        # Use default radius of 1.5 Angstroms if not found
        if atom_type in vdw_radii_dict:
            vdw_radius = vdw_radii_dict[atom_type]
        else:
            print(
                f"Warning: Van der Waals radius not found for {atom_type}. Using default value of 1.5 Angstroms."
            )
            vdw_radius = 1.5

        # Step 1: Create the original (point) representation
        for _, row in type_coords.iterrows():
            # Get atom coordinates
            x0, y0, z0 = row["x"], row["y"], row["z"]

            # Convert to grid indices
            ix, iy, iz = coords_to_grid_indices(x0, y0, z0, grid_resolution=0.5)

            # Ensure the atom is within the grid boundaries
            if is_within_grid(ix, iy, iz, n_points):
                # Mark this point as containing an atom
                original_grid[i, ix, iy, iz] = 1

        # Step 2: Create the blurred representation with Gaussian distribution
        # Calculate normalization once per atom type
        norm_factor = gaussian_3d_norm(vdw_radius)

        for _, row in type_coords.iterrows():
            # Get atom coordinates
            x0, y0, z0 = row["x"], row["y"], row["z"]

            # Loop through all grid points (this is computationally expensive!)
            for ix in range(n_points):
                for iy in range(n_points):
                    for iz in range(n_points):
                        # Get actual coordinates of this grid point
                        x_prime = x[ix]
                        y_prime = y[iy]
                        z_prime = z[iz]

                        # Calculate differences between atom and grid point
                        x_diff = x_prime - x0
                        y_diff = y_prime - y0
                        z_diff = z_prime - z0

                        # Calculate Gaussian value at this point
                        gauss_val = gaussian_value(
                            x_diff, y_diff, z_diff, vdw_radius, norm_factor
                        )

                        # Add this atom's contribution to the grid point
                        blurred_grid[i, ix, iy, iz] += gauss_val

    return original_grid, blurred_grid


def optimize_for_large_data(atom_coords_df, vdw_radii_dict, max_radius=5.0):
    """
    Fast version of Gaussian blur for large datasets (recommended for most uses).

    Instead of calculating Gaussian values for all grid points (which is very slow),
    this function only applies the blur to:
    1. The grid point containing each atom
    2. The 26 immediate neighboring grid points

    This optimization dramatically speeds up processing while maintaining acceptable
    visual quality for most applications.

    Parameters:
    -----------
    atom_coords_df : pandas.DataFrame
        Dataframe with atom coordinates.
        Expected structure:
            - MultiIndex with the last level being 'atom_type'
            - Columns should include 'x', 'y', and 'z' coordinates

    vdw_radii_dict : dict
        Dictionary mapping atom types (str) to their Van der Waals radii (float).
        Example: {'C': 1.7, 'N': 1.55, 'O': 1.52}

    max_radius : float
        Maximum radius to consider (not used in current implementation).
        Kept for future extensions.

    Returns:
    --------
    tuple of numpy.ndarray
        Two 4D arrays:
        - original_grid: Binary representation where atoms are single points
        - blurred_grid: Optimized smoothed representation with limited Gaussian blur
    """
    # Get unique atom types from the dataframe
    atom_types = get_atom_types(atom_coords_df)

    # Set up the 3D grid
    x, y, z, n_points = create_grid(grid_size=18, grid_resolution=0.5)

    # Initialize empty grids
    original_grid, blurred_grid = initialize_grids(len(atom_types), n_points)

    # Get the 26 neighboring offsets
    neighbor_offsets = get_neighbor_offsets()

    # Process each atom type separately
    for i, atom_type in enumerate(atom_types):
        # Get coordinates for all atoms of current type
        type_coords = atom_coords_df.xs(atom_type, level="atom_type")

        # Get the Van der Waals radius for this atom type
        # Use default radius of 1.5 Angstroms if not found
        if atom_type in vdw_radii_dict:
            vdw_radius = vdw_radii_dict[atom_type]
        else:
            print(
                f"Warning: Van der Waals radius not found for {atom_type}. Using default value of 1.5 Angstroms."
            )
            vdw_radius = 1.5

        # Step 1: Create the original (point) representation
        for _, row in type_coords.iterrows():
            # Get atom coordinates
            x0, y0, z0 = row["x"], row["y"], row["z"]

            # Convert to grid indices
            ix, iy, iz = coords_to_grid_indices(x0, y0, z0, grid_resolution=0.5)

            # Ensure the atom is within the grid boundaries
            if is_within_grid(ix, iy, iz, n_points):
                # Mark this point as containing an atom
                original_grid[i, ix, iy, iz] = 1

        # Calculate the normalizing factor for 3D Gaussian once per atom type
        norm_factor = gaussian_3d_norm(vdw_radius)

        # Step 2: Create the optimized blurred representation
        for _, row in type_coords.iterrows():
            # Get atom coordinates
            x0, y0, z0 = row["x"], row["y"], row["z"]

            # Convert to grid indices for the center point
            center_ix, center_iy, center_iz = coords_to_grid_indices(
                x0, y0, z0, grid_resolution=0.5
            )

            # Skip if atom center is outside grid boundaries
            if not is_within_grid(center_ix, center_iy, center_iz, n_points):
                continue

            # Step 2a: Apply Gaussian to the center point
            x_prime = x[center_ix]
            y_prime = y[center_iy]
            z_prime = z[center_iz]

            # Calculate differences between atom and grid point
            x_diff = x_prime - x0
            y_diff = y_prime - y0
            z_diff = z_prime - z0

            # Calculate Gaussian value at the center point
            gaussian_val = gaussian_value(
                x_diff, y_diff, z_diff, vdw_radius, norm_factor
            )

            # Set the value at this grid point (take maximum of current and new value)
            blurred_grid[i, center_ix, center_iy, center_iz] = max(
                blurred_grid[i, center_ix, center_iy, center_iz], gaussian_val
            )

            # Step 2b: Apply Gaussian to the 26 neighboring grid points
            for dx, dy, dz in neighbor_offsets:
                # Calculate neighbor coordinates
                nx, ny, nz = center_ix + dx, center_iy + dy, center_iz + dz

                # Check if neighbor is within grid boundaries
                if is_within_grid(nx, ny, nz, n_points):
                    # Get actual coordinates of this neighbor
                    x_prime = x[nx]
                    y_prime = y[ny]
                    z_prime = z[nz]

                    # Calculate differences between atom and neighbor point
                    x_diff = x_prime - x0
                    y_diff = y_prime - y0
                    z_diff = z_prime - z0

                    # Calculate Gaussian value at this neighbor
                    gaussian_val = gaussian_value(
                        x_diff, y_diff, z_diff, vdw_radius, norm_factor
                    )

                    # Set the value (take maximum of current and new value)
                    blurred_grid[i, nx, ny, nz] = max(
                        blurred_grid[i, nx, ny, nz], gaussian_val
                    )

    return original_grid, blurred_grid


def create_protein_visualization(
    atom_coords_df,
    vdw_radii_dict=None,
    use_optimized=True,
    output_filename="protein_visualization.html",
    opacity=None,
):
    """
    Simple function to create a protein visualization from atom coordinates.

    This is a beginner-friendly wrapper that handles the common steps needed to:
    1. Apply Gaussian blur to atom coordinates
    2. Create a 3D visualization
    3. Save the result as an interactive HTML file

    Parameters:
    -----------
    atom_coords_df : pandas.DataFrame
        Dataframe with atom coordinates.
        Must have a MultiIndex with 'atom_type' as the last level.

    vdw_radii_dict : dict, optional
        Dictionary of Van der Waals radii for each atom type.
        If None, uses the built-in VDW_RADII dictionary.

    use_optimized : bool, default=True
        If True, uses the faster optimized algorithm.
        If False, uses the more accurate but slower algorithm.

    output_filename : str, default="protein_visualization.html"
        Filename to save the visualization HTML.

    opacity : dict, optional
        Opacity settings for visualization.
        If None, uses settings that make voxels more transparent and less visible

    Returns:
    --------
    tuple
        The generated figure object and the path to the saved HTML file.
    """
    # Get atom types from the dataframe
    atom_types = get_atom_types(atom_coords_df)

    # If no radii dictionary provided, use the built-in one
    if vdw_radii_dict is None:
        vdw_radii_dict = get_atom_radii(atom_types)

    # Print basic information
    print(f"Processing {len(atom_coords_df)} atoms of {len(atom_types)} types")

    # Apply the appropriate Gaussian blur function
    print("Applying Gaussian blur...")
    if use_optimized:
        original_grid, blurred_grid = optimize_for_large_data(
            atom_coords_df, vdw_radii_dict
        )
    else:
        original_grid, blurred_grid = apply_gaussian_blur(
            atom_coords_df, vdw_radii_dict
        )

    # Set default opacity settings if none provided
    # Use settings that make voxels more transparent and less visible
    if opacity is None:
        opacity = adjust_opacity_settings(
            base_opacity=0.2,    # More transparent
            min_opacity=0.02,    # Lower minimum opacity
            opacity_scale=0.3    # More subtle scaling
        )

    # Generate the visualization with the opacity settings
    print("Creating visualization...")
    fig = visualize_grids(
        original_grid, blurred_grid, atom_types, opacity_settings=opacity
    )

    # Save the visualization
    output_path = save_visualization(fig, output_filename)

    return fig, output_path


# Example usage demonstration
if __name__ == "__main__":
    """
    This section demonstrates how to use the Gaussian blur functions.
    It creates a sample protein, applies the blur, and visualizes the results.
    """
    # Create example protein data with 300 atoms
    atom_coords_df = create_protein_example(n_atoms=300)

    # Print summary statistics
    print_atom_stats(atom_coords_df)

    # Create visualization with lower opacity for subtle visibility
    opacity_settings = adjust_opacity_settings(
        base_opacity=0.2,    # More transparent
        min_opacity=0.02,    # Lower minimum opacity
        opacity_scale=0.3    # More subtle scaling
    )

    # Use the simplified function with custom opacity
    create_protein_visualization(
        atom_coords_df,
        output_filename="protein_visualization_subtle.html",
        opacity=opacity_settings,
    )
