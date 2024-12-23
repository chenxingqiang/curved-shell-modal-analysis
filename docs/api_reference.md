# API Reference

## Surface Classes

### Surface
Base class for all surface types.

```matlab
surface = CurvedShellAnalysis.Surface(params)
```

#### Properties
- `Nodes`: Node coordinates (N×3 matrix)
- `Elements`: Element connectivity (M×4 matrix)
- `Material`: Material properties
- `Thickness`: Shell thickness
- `StiffnessMatrix`: Global stiffness matrix
- `MassMatrix`: Global mass matrix

#### Methods
- `calculateStiffnessMatrix()`: Compute element and global stiffness matrices
- `calculateMassMatrix()`: Compute element and global mass matrices
- `calculateStressStrain(displacement)`: Calculate stress and strain from displacement
- `plot()`: Visualize the surface

### SphericalSurface
```matlab
sphere = CurvedShellAnalysis.SphericalSurface(params)
```

Additional properties:
- `Radius`: Sphere radius
- `Angles`: Angular coverage [θ_start, θ_end, φ_start, φ_end]

### CylindricalSurface
```matlab
cylinder = CurvedShellAnalysis.CylindricalSurface(params)
```

Additional properties:
- `Radius`: Cylinder radius
- `Length`: Cylinder length
- `Angles`: Angular coverage [θ_start, θ_end]

## Analysis Classes

### ModalAnalysis
```matlab
analysis = CurvedShellAnalysis.ModalAnalysis(surface, num_modes)
```

#### Properties
- `Surface`: Surface object
- `NumModes`: Number of modes to compute
- `NaturalFrequencies`: Computed natural frequencies
- `ModeShapes`: Computed mode shapes

#### Methods
- `analyze()`: Perform modal analysis
- `plotMode(mode_num)`: Plot specific mode shape
- `animateMode(mode_num)`: Animate mode shape

### NonlinearAnalysis
```matlab
analysis = CurvedShellAnalysis.NonlinearAnalysis(surface)
```

#### Properties
- `Surface`: Surface object
- `LoadSteps`: Load step data
- `Tolerance`: Convergence tolerance
- `MaxIterations`: Maximum iterations per step

#### Methods
- `solveNonlinear(load_steps)`: Solve nonlinear problem
- `calculateTangentStiffness()`: Update tangent stiffness matrix
- `checkConvergence()`: Check solution convergence

### ThermalStructuralAnalysis
```matlab
analysis = CurvedShellAnalysis.ThermalStructuralAnalysis(surface)
```

#### Properties
- `Surface`: Surface object
- `Temperature`: Temperature field
- `ThermalLoads`: Thermal load vector
- `ThermalExpansion`: Thermal expansion coefficients

#### Methods
- `setTemperature(temp_field)`: Set temperature distribution
- `calculateThermalLoads()`: Calculate thermal load vector
- `solve()`: Solve coupled problem

## Material Classes

### AdvancedMaterial
```matlab
material = CurvedShellAnalysis.AdvancedMaterial(type, props)
```

#### Properties
- `Type`: Material type ('elastic', 'viscoelastic', 'damage', 'composite')
- `Properties`: Material properties
- `State`: Current material state

#### Methods
- `calculateStress(strain)`: Calculate stress response
- `updateState(strain, dt)`: Update internal state variables
- `getDamageTensor()`: Get current damage tensor

### FailureCriteria
```matlab
failure = CurvedShellAnalysis.FailureCriteria(type, props)
```

#### Properties
- `Type`: Failure criterion type
- `Properties`: Material properties
- `Options`: Analysis options

#### Methods
- `evaluateFailure(stress, strain)`: Evaluate failure criterion
- `tsaiWuCriterion(stress)`: Tsai-Wu criterion
- `puckCriterion(stress)`: Puck criterion
- `larcCriterion(stress, strain)`: LaRC criterion
- `multicontinuumFailure(stress, strain)`: MCT failure analysis

## Optimization Classes

### AdvancedOptimization
```matlab
opt = CurvedShellAnalysis.AdvancedOptimization(problem, algorithm, options)
```

#### Properties
- `Problem`: Optimization problem definition
- `Algorithm`: Optimization algorithm
- `Options`: Algorithm options
- `Results`: Optimization results

#### Methods
- `optimize()`: Run optimization
- `evaluateObjective(x)`: Evaluate objective function
- `checkConstraints(x)`: Check constraint satisfaction
- `updateSurrogate()`: Update surrogate model

## Parallel Computing

### ParallelManager
```matlab
pm = CurvedShellAnalysis.ParallelManager(num_workers, batch_size, memory_limit)
```

#### Properties
- `NumWorkers`: Number of parallel workers
- `BatchSize`: Batch size for parallel processing
- `MemoryLimit`: Memory limit per worker
- `GPUDevices`: Available GPU devices

#### Methods
- `initializeParallel()`: Initialize parallel environment
- `addTask(task_type, params)`: Add task to queue
- `executeTasks()`: Execute tasks in parallel
- `executeGPUBatch(batch_tasks)`: Execute batch on GPU
