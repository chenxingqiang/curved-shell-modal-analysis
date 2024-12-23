# Modal Analysis Package for Curved Surfaces

A comprehensive MATLAB package for modal analysis of curved surfaces, including spherical, ellipsoidal, cylindrical, and conical surfaces. The package supports advanced features such as multi-scale analysis, fluid-structure interaction, topology optimization, and progressive failure analysis.

## Features

### Modal Analysis

- Natural frequencies and mode shapes computation
- Mass normalization and orthogonality checks
- Frequency response analysis
- Random vibration analysis

### Surface Types

- Spherical shells
- Ellipsoidal shells
- Cylindrical shells
- Conical shells
- Custom parametric surfaces

### Analysis Types

- Modal analysis
- Dynamic response
- Thermal-structural coupling
- Buckling analysis
- Fluid-structure interaction
- Topology optimization

## Results Gallery

### Modal Analysis Results

#### Spherical Shell

| Mode 1 | Mode 2 | Mode 3 | Mode 4 |
|--------|--------|--------|--------|
| ![Mode 1](docs/images/spherical_mode1.png) | ![Mode 2](docs/images/spherical_mode2.png) | ![Mode 3](docs/images/spherical_mode3.png) | ![Mode 4](docs/images/spherical_mode4.png) |

#### Cylindrical Shell

| Mode 1 | Mode 2 | Mode 3 | Mode 4 |
|--------|--------|--------|--------|
| ![Mode 1](docs/images/cylindrical_mode1.png) | ![Mode 2](docs/images/cylindrical_mode2.png) | ![Mode 3](docs/images/cylindrical_mode3.png) | ![Mode 4](docs/images/cylindrical_mode4.png) |

### Fluid-Structure Interaction

| Pressure Field | Deformation |
|----------------|-------------|
| ![Pressure](docs/images/fsi_pressure.png) | ![Deformation](docs/images/fsi_deformation.png) |

### Topology Optimization

| Initial Design | Optimized Design | Convergence |
|----------------|------------------|-------------|
| ![Initial](docs/images/topo_initial.png) | ![Final](docs/images/topo_final.png) | ![Convergence](docs/images/topo_convergence.png) |

## Installation

1. Clone the repository:

```bash
git clone https://github.com/gxingqiang/modal-analysis.git
```

2. Add the package to your MATLAB path:

```matlab
addpath('/path/to/modal-analysis');
```

## Quick Start

```matlab
% Create a cylindrical shell
params = struct('R', 0.3, 'L', 0.6, 't', 0.002);
shell = CurvedShellAnalysis.CylindricalSurface(params);

% Set material properties (Steel)
material = struct('E', 210e9, 'nu', 0.3, 'rho', 7800);
shell.setMaterial(material);

% Perform modal analysis
modal = CurvedShellAnalysis.ModalAnalysis(shell, 10);
modal.analyze();

% Plot mode shapes
modal.plotModes(1:4);
```

## Examples

Check out the demos folder for comprehensive examples:

1. `demo_modal_analysis.m`: Basic modal analysis
   - Natural frequencies and mode shapes
   - Frequency response functions
   - Dynamic response analysis

2. `demo_composite_shell.m`: Composite shell analysis
   - Laminate properties
   - Thermal effects
   - Progressive failure

3. `demo_fsi_analysis.m`: Fluid-structure interaction
   - Pressure fields
   - Flow-induced vibration
   - Flutter analysis

4. `demo_topology_optimization.m`: Shell optimization
   - Compliance minimization
   - Frequency constraints
   - Manufacturing constraints

## Documentation

- [Theory Manual](docs/theory_manual.md)
- [Advanced Examples](docs/advanced_examples.md)
- [API Reference](docs/api_reference.md)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Xingqiang Chen

## Citation

If you use this package in your research, please cite:

```bibtex
@software{chen2024modal,
  author = {Chen, Xingqiang},
  title = {Modal Analysis Package for Curved Surfaces},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/gxingqiang/modal-analysis}
}
```

## Contact

For questions and feedback:

- Email: <gxingqiang@gmail.com>
- GitHub Issues: [https://github.com/gxingqiang/modal-analysis/issues](https://github.com/gxingqiang/modal-analysis/issues)
