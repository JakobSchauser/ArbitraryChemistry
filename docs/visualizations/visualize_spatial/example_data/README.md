# Example Spatial Visualization Data

This folder contains example lattice data files that demonstrate the CSV format used by the spatial visualization tool.

## Files

- `lattice_0.txt` - Initial state with Chemical 1 entering from the top
- `lattice_100.txt` - Diffusion spreading the chemical further
- `lattice_200.txt` - Chemical 2 begins to appear
- `lattice_300.txt` - More balanced distribution with Chemical 2
- `lattice_400.txt` - Chemical 3 appears, more complex mixing

## Format

Each file follows the standard format from the CUDA simulation:

```
# Lattice size: 10x10
# Molecules per site: 100
# Chemical range: 0-100 (0=empty)
# Format: x,y,chemical,count
x,y,chemical,count
[data rows...]
```

## Visualization Patterns

### Time Step 0
- **Chemical 1 (Red)**: Concentrated at top row (supply source)
- **Chemical 0 (Empty)**: Fills remaining space
- **Pattern**: Linear gradient from top to bottom

### Time Step 100
- **Chemical 1**: Spreading downward through diffusion
- **Pattern**: Wider gradient distribution

### Time Step 200
- **Chemical 1**: More evenly distributed
- **Chemical 2 (Green)**: Small amounts appear
- **Pattern**: Dual chemical mixing begins

### Time Step 300
- **Chemical 1**: Background levels
- **Chemical 2**: Increased presence
- **Pattern**: Green becomes more prominent

### Time Step 400
- **Chemical 1**: Low background
- **Chemical 2**: Medium levels
- **Chemical 3 (Blue)**: New chemical appears
- **Pattern**: Three-chemical RGB mixing

## How to Use

1. Open `spatial_visualization.html`
2. The example data path is pre-loaded: `docs/visualizations/visualize_spatial/example_data`
3. Click "Load Data"
4. Select chemicals for RGB channels:
   - Red: Chemical 1
   - Green: Chemical 2  
   - Blue: Chemical 3
5. Use time slider to see evolution
6. Hover over cells for detailed counts

## Expected Visualization

- **Early steps**: Red gradient from top to bottom
- **Middle steps**: Red-green mixing
- **Later steps**: Red-green-blue combinations
- **Patterns**: Diffusion and reaction create color gradients and mixing zones

This demonstrates how chemical concentrations can be visualized as RGB color mixing to reveal spatial patterns in reaction-diffusion systems.