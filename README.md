# Modal Analysis Package for Curved Surfaces

A comprehensive MATLAB package for modal analysis of curved surfaces, including spherical, ellipsoidal, cylindrical, and conical surfaces. The package supports advanced features such as multi-scale analysis, fluid-structure interaction, topology optimization, and progressive failure analysis.

## Results Gallery

### Modal Analysis
![Modal Analysis](docs/images/modal_analysis.png)

### Frequency Response Analysis
![Frequency Response](docs/images/frequency_response.png)

### Individual Mode Shapes
| Mode 1 (89.67 Hz) | Mode 2 (89.67 Hz) | Mode 3 (89.72 Hz) |
|-------------------|-------------------|-------------------|
| ![Mode 1](docs/images/mode_1_89667_hz.png) | ![Mode 2](docs/images/mode_2_89667_hz.png) | ![Mode 3](docs/images/mode_3_89720_hz.png) |

| Mode 4 (89.72 Hz) | Mode 5 (89.81 Hz) | Mode 6 (89.81 Hz) |
|-------------------|-------------------|-------------------|
| ![Mode 4](docs/images/mode_4_89723_hz.png) | ![Mode 5](docs/images/mode_5_89807_hz.png) | ![Mode 6](docs/images/mode_6_89813_hz.png) |

### Advanced Analysis Results

#### Comparative Modal Analysis
![Comparative Modes](docs/images/comparative_modes.png)
*Comparison of mode shapes between different surface types*

#### Composite Shell Analysis
![Composite Modes](docs/images/composite_modes.png)
*Modal analysis results for composite shell structures*

#### Dynamic Response Analysis
![Dynamic Response](docs/images/dynamic_response.png)
*Time-domain response under dynamic loading*

#### Frequency Analysis
| Frequency Response Function | MAC Matrix Analysis |
|----------------------------|---------------------|
| ![FRF Analysis](docs/images/FRF.png) | ![MAC Matrices](docs/images/mac_matrices.png) |

#### Strain Energy Distribution
![Strain Energy](docs/images/strain_energy.png)
*Strain energy distribution during modal deformation*

#### Thermal Effects
![Thermal Deformation](docs/images/thermal_deformation.png)
*Thermal-structural coupling analysis results*

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
- Email: chen.xingqiang@iechor.com
- GitHub Issues: [https://github.com/gxingqiang/modal-analysis/issues](https://github.com/gxingqiang/modal-analysis/issues)
