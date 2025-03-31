import pandas as pd
from vdw_radii import VDW_RADII

"""
Atom Utilities for Protein Visualization

This module provides functions for working with atom data for protein visualization.
"""

def get_atom_types(atom_coords_df):
    """
    Extract unique atom types from the dataframe.
    
    Parameters:
    -----------
    atom_coords_df : pandas.DataFrame
        Dataframe with atom coordinates and atom_type in MultiIndex
        
    Returns:
    --------
    list
        List of unique atom types
    """
    return atom_coords_df.index.get_level_values('atom_type').unique().tolist()

def get_atom_radii(atom_types, custom_radii=None):
    """
    Get Van der Waals radii for specified atom types.
    
    Parameters:
    -----------
    atom_types : list
        List of atom types to get radii for
    
    custom_radii : dict, optional
        Custom radius values to use instead of defaults
        
    Returns:
    --------
    dict
        Dictionary mapping atom types to their Van der Waals radii
    """
    if custom_radii is None:
        custom_radii = {}
        
    # Use custom radii when provided, otherwise fall back to built-in values
    # Use 1.5 Angstroms as default if not found
    return {atom_type: custom_radii.get(atom_type, VDW_RADII.get(atom_type, 1.5)) 
            for atom_type in atom_types}

def print_atom_stats(atom_coords_df):
    """
    Print basic statistics about the atoms in the dataframe.
    
    Parameters:
    -----------
    atom_coords_df : pandas.DataFrame
        Dataframe with atom coordinates
    """
    print(f"Total atoms: {len(atom_coords_df)}")
    
    type_counts = atom_coords_df.index.get_level_values('atom_type').value_counts()
    print("Atom type distribution:")
    for atom_type, count in type_counts.items():
        print(f"  {atom_type}: {count}")
        
    return type_counts 