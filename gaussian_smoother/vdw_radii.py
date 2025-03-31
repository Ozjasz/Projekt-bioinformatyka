"""
Van der Waals radii for elements based on Bondi's values, which are widely used in molecular modeling.
Values are in Angstroms (Å).

For more information, see:
A. Bondi, "van der Waals Volumes and Radii," J. Phys. Chem., 1964, 68, 441-451.
"""

# Van der Waals radii in Angstroms (Å)
VDW_RADII = {
    # Period 1
    'H': 1.20,  # Hydrogen
    'He': 1.40,  # Helium
    
    # Period 2
    'Li': 1.82,  # Lithium
    'Be': 1.53,  # Beryllium
    'B': 1.92,  # Boron
    'C': 1.70,  # Carbon
    'N': 1.55,  # Nitrogen
    'O': 1.52,  # Oxygen
    'F': 1.47,  # Fluorine
    'Ne': 1.54,  # Neon
    
    # Period 3
    'Na': 2.27,  # Sodium
    'Mg': 1.73,  # Magnesium
    'Al': 1.84,  # Aluminum
    'Si': 2.10,  # Silicon
    'P': 1.80,  # Phosphorus
    'S': 1.80,  # Sulfur
    'Cl': 1.75,  # Chlorine
    'Ar': 1.88,  # Argon
    
    # Period 4
    'K': 2.75,  # Potassium
    'Ca': 2.31,  # Calcium
    'Sc': 2.15,  # Scandium
    'Ti': 2.11,  # Titanium
    'V': 2.07,  # Vanadium
    'Cr': 2.06,  # Chromium
    'Mn': 2.05,  # Manganese
    'Fe': 2.04,  # Iron
    'Co': 2.00,  # Cobalt
    'Ni': 1.97,  # Nickel
    'Cu': 1.96,  # Copper
    'Zn': 2.01,  # Zinc
    'Ga': 1.87,  # Gallium
    'Ge': 2.11,  # Germanium
    'As': 1.85,  # Arsenic
    'Se': 1.90,  # Selenium
    'Br': 1.85,  # Bromine
    'Kr': 2.02,  # Krypton
    
    # Period 5
    'Rb': 3.03,  # Rubidium
    'Sr': 2.49,  # Strontium
    'Y': 2.32,  # Yttrium
    'Zr': 2.23,  # Zirconium
    'Nb': 2.18,  # Niobium
    'Mo': 2.17,  # Molybdenum
    'Tc': 2.16,  # Technetium
    'Ru': 2.13,  # Ruthenium
    'Rh': 2.10,  # Rhodium
    'Pd': 2.10,  # Palladium
    'Ag': 2.11,  # Silver
    'Cd': 2.18,  # Cadmium
    'In': 1.93,  # Indium
    'Sn': 2.17,  # Tin
    'Sb': 2.06,  # Antimony
    'Te': 2.06,  # Tellurium
    'I': 1.98,  # Iodine
    'Xe': 2.16,  # Xenon
    
    # Period 6
    'Cs': 3.43,  # Cesium
    'Ba': 2.68,  # Barium
    'La': 2.43,  # Lanthanum
    'Ce': 2.42,  # Cerium
    'Pr': 2.40,  # Praseodymium
    'Nd': 2.39,  # Neodymium
    'Pm': 2.38,  # Promethium
    'Sm': 2.36,  # Samarium
    'Eu': 2.35,  # Europium
    'Gd': 2.34,  # Gadolinium
    'Tb': 2.33,  # Terbium
    'Dy': 2.31,  # Dysprosium
    'Ho': 2.30,  # Holmium
    'Er': 2.29,  # Erbium
    'Tm': 2.27,  # Thulium
    'Yb': 2.26,  # Ytterbium
    'Lu': 2.24,  # Lutetium
    'Hf': 2.23,  # Hafnium
    'Ta': 2.22,  # Tantalum
    'W': 2.18,  # Tungsten
    'Re': 2.16,  # Rhenium
    'Os': 2.16,  # Osmium
    'Ir': 2.13,  # Iridium
    'Pt': 2.13,  # Platinum
    'Au': 2.14,  # Gold
    'Hg': 2.23,  # Mercury
    'Tl': 1.96,  # Thallium
    'Pb': 2.02,  # Lead
    'Bi': 2.07,  # Bismuth
    'Po': 1.97,  # Polonium
    'At': 2.02,  # Astatine
    'Rn': 2.20,  # Radon
    
    # Period 7
    'Fr': 3.48,  # Francium
    'Ra': 2.83,  # Radium
    'Ac': 2.47,  # Actinium
    'Th': 2.45,  # Thorium
    'Pa': 2.43,  # Protactinium
    'U': 2.41,  # Uranium
    'Np': 2.39,  # Neptunium
    'Pu': 2.43,  # Plutonium
    'Am': 2.44,  # Americium
    
    # Common atom types in proteins
    'ligand': 1.80,  # Generic ligand (approximate value)
}

# Common amino acid-specific atom types
VDW_RADII_AA = {
    'CA': 1.70,  # Alpha carbon
    'CB': 1.70,  # Beta carbon
    'CD': 1.70,  # Delta carbon
    'CE': 1.70,  # Epsilon carbon
    'CZ': 1.70,  # Zeta carbon
    'NZ': 1.55,  # Zeta nitrogen (e.g., in Lysine)
    'OD': 1.52,  # Delta oxygen (e.g., in Aspartate)
    'OE': 1.52,  # Epsilon oxygen (e.g., in Glutamate)
    'SD': 1.80,  # Delta sulfur (e.g., in Methionine)
    'SG': 1.80,  # Gamma sulfur (e.g., in Cysteine)
}

if __name__ == "__main__":
    print(f"Van der Waals radii dictionary contains {len(VDW_RADII)} elements")
    print(f"Amino acid specific atom types dictionary contains {len(VDW_RADII_AA)} entries")
    
    # Print some common values
    common_elements = ['H', 'C', 'N', 'O', 'S', 'P', 'Fe', 'Zn', 'Mg', 'Ca']
    print("\nCommonly used VDW radii in biochemistry:")
    for elem in common_elements:
        if elem in VDW_RADII:
            print(f"{elem}: {VDW_RADII[elem]} Å") 