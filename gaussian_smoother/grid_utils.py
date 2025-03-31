import numpy as np

"""
Grid Utilities for Protein Visualization

This module provides functions for creating and manipulating 3D grids
used in protein structure visualization.
"""

def create_grid(grid_size=18, grid_resolution=0.5):
    """
    Create a 3D grid with specified size and resolution.
    
    Parameters:
    -----------
    grid_size : float, default=18
        Size of cubic space in Angstroms (e.g., 18 x 18 x 18)
    
    grid_resolution : float, default=0.5
        Distance between grid points in Angstroms
        
    Returns:
    --------
    tuple
        (x, y, z, n_points) where:
        - x, y, z are numpy arrays with grid coordinates
        - n_points is the number of grid points in each dimension
    """
    # Calculate number of grid points in each dimension
    n_points = int(grid_size / grid_resolution)
    
    # Create evenly spaced grid coordinates in each dimension
    x = np.linspace(0, grid_size, n_points)
    y = np.linspace(0, grid_size, n_points)
    z = np.linspace(0, grid_size, n_points)
    
    return x, y, z, n_points

def initialize_grids(n_atom_types, n_points):
    """
    Initialize empty grids for original and blurred representations.
    
    Parameters:
    -----------
    n_atom_types : int
        Number of different atom types to represent
    
    n_points : int
        Number of grid points in each dimension
        
    Returns:
    --------
    tuple of numpy.ndarray
        Two 4D arrays (original_grid, blurred_grid) initialized with zeros
    """
    original_grid = np.zeros((n_atom_types, n_points, n_points, n_points))
    blurred_grid = np.zeros((n_atom_types, n_points, n_points, n_points))
    
    return original_grid, blurred_grid

def coords_to_grid_indices(x0, y0, z0, grid_resolution):
    """
    Convert atom coordinates to grid indices.
    
    Parameters:
    -----------
    x0, y0, z0 : float
        Atom coordinates in Angstroms
    
    grid_resolution : float
        Distance between grid points in Angstroms
        
    Returns:
    --------
    tuple of int
        (ix, iy, iz) grid indices
    """
    ix = int(x0 / grid_resolution)
    iy = int(y0 / grid_resolution)
    iz = int(z0 / grid_resolution)
    
    return ix, iy, iz

def is_within_grid(ix, iy, iz, n_points):
    """
    Check if indices are within grid boundaries.
    
    Parameters:
    -----------
    ix, iy, iz : int
        Grid indices to check
    
    n_points : int
        Number of grid points in each dimension
        
    Returns:
    --------
    bool
        True if indices are within bounds, False otherwise
    """
    return 0 <= ix < n_points and 0 <= iy < n_points and 0 <= iz < n_points

def get_neighbor_offsets():
    """
    Get relative coordinates for the 26 neighbors in a 3x3x3 cube.
    
    Returns:
    --------
    list of tuples
        List of (dx, dy, dz) offsets excluding (0,0,0)
    """
    neighbor_offsets = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx != 0 or dy != 0 or dz != 0:  # Skip the center point
                    neighbor_offsets.append((dx, dy, dz))
    return neighbor_offsets

def gaussian_3d_norm(vdw_radius):
    """
    Calculate normalizing factor for 3D Gaussian distribution.
    
    Parameters:
    -----------
    vdw_radius : float
        Van der Waals radius (sigma in the Gaussian formula)
        
    Returns:
    --------
    float
        Normalization factor for 3D Gaussian
    """
    return 1 / ((2 * np.pi * vdw_radius**2)**(3/2))

def gaussian_value(x_diff, y_diff, z_diff, vdw_radius, norm_factor=None):
    """
    Calculate the value of a 3D Gaussian function at a given distance.
    
    Parameters:
    -----------
    x_diff, y_diff, z_diff : float
        Differences between atom coordinates and grid point
    
    vdw_radius : float
        Van der Waals radius (sigma in the Gaussian formula)
    
    norm_factor : float, optional
        Pre-calculated normalization factor. If None, it will be calculated.
        
    Returns:
    --------
    float
        Gaussian value at the specified point
    """
    # Calculate squared distance
    dist_sq = x_diff**2 + y_diff**2 + z_diff**2
    
    # Calculate or use provided normalization factor
    if norm_factor is None:
        norm_factor = gaussian_3d_norm(vdw_radius)
    
    # Apply Gaussian formula: G(r) = norm * exp(-r²/(2σ²))
    return norm_factor * np.exp(-dist_sq / (2 * vdw_radius**2)) 