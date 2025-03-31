"""
Visualization Utilities for Protein Visualization

This module provides functions for creating and configuring visualizations
of protein structures.
"""

def get_plot_config(width=1900, height=950):
    """
    Create configuration for Plotly plots with sensible defaults.
    
    Parameters:
    -----------
    width : int, default=1900
        Width of the plot image in pixels
    
    height : int, default=950
        Height of the plot image in pixels
        
    Returns:
    --------
    dict
        Plotly configuration dictionary
    """
    return {
        'displayModeBar': True,
        'responsive': True,
        'toImageButtonOptions': {
            'format': 'png',
            'width': width,
            'height': height
        }
    }

def create_html_with_plot(fig, title="Protein Visualization"):
    """
    Create HTML string with centered plot and nice styling.
    
    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The plotly figure to embed
    
    title : str, default="Protein Visualization"
        Title for the HTML page
        
    Returns:
    --------
    str
        Complete HTML string with embedded plot
    """
    config = get_plot_config()
    
    html_string = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8" />
        <title>{title}</title>
        <style>
            body {{
                margin: 0;
                padding: 0;
                display: flex;
                justify-content: center;
                align-items: center;
                min-height: 100vh;
                background-color: #f5f5f5;
            }}
            .plot-container {{
                margin: 20px auto;
                max-width: 1200px;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                background: white;
                border-radius: 8px;
                overflow: hidden;
            }}
        </style>
    </head>
    <body>
        <div class="plot-container">
            {fig.to_html(config=config, include_plotlyjs='cdn', full_html=False)}
        </div>
    </body>
    </html>
    """
    
    return html_string

def adjust_opacity_settings(base_opacity=0.7, min_opacity=0.2, opacity_scale=0.8):
    """
    Create opacity settings for visualization.
    
    This function returns settings that control how transparent/opaque the
    blurred voxels appear in the visualization. Higher values make voxels
    more opaque and more visible.
    
    Parameters:
    -----------
    base_opacity : float, default=0.7
        Base opacity for voxels (0.0 = fully transparent, 1.0 = fully opaque)
        Recommend values between 0.5-0.8 for best visualization
    
    min_opacity : float, default=0.2
        Minimum opacity for low-value voxels
        
    opacity_scale : float, default=0.8
        How quickly opacity increases with voxel value
        Higher values make even dim voxels more visible
        
    Returns:
    --------
    dict
        Opacity settings that can be passed to visualize_grids
    """
    return {
        'base_opacity': base_opacity,
        'min_opacity': min_opacity,
        'opacity_scale': opacity_scale
    }

def save_visualization(fig, output_filename="protein_visualization.html", title="Protein Visualization"):
    """
    Save visualization as an HTML file.
    
    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The plotly figure to save
    
    output_filename : str, default="protein_visualization.html"
        Filename to save the HTML to
    
    title : str, default="Protein Visualization"
        Title for the HTML page
        
    Returns:
    --------
    str
        Path to the saved file
    """
    html_string = create_html_with_plot(fig, title)
    
    with open(output_filename, "w") as f:
        f.write(html_string)
    
    print(f"Visualization saved to {output_filename}")
    return output_filename 