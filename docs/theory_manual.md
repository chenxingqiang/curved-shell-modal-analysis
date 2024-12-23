# Theory Manual

## Shell Theory

### Kinematics

The displacement field of a shell is described by:

```
u(ξ¹,ξ²,ξ³) = v(ξ¹,ξ²) + ξ³θ(ξ¹,ξ²)
```

where:

- (ξ¹,ξ²,ξ³) are curvilinear coordinates
- v is the displacement of the middle surface
- θ is the rotation of the normal vector

### Strain-Displacement Relations

The Green-Lagrange strain tensor is:

```
E = ε + ξ³κ
```

where:

- ε is the membrane strain
- κ is the curvature change

### Constitutive Relations

For linear elastic materials:

```
σ = C : ε
```

For viscoelastic materials (Prony series):

```
σ(t) = E∞ε(t) + Σᵢ Eᵢexp(-t/τᵢ)∫₀ᵗ exp(s/τᵢ)(dε/ds)ds
```

### Element Formulation

MITC4 shell element with:

- 4 nodes
- 5 DOFs per node (3 displacements, 2 rotations)
- Mixed interpolation of tensorial components
- Assumed natural strain method

## Modal Analysis

### Eigenvalue Problem

```
(K - ω²M)φ = 0
```

where:

- K is the stiffness matrix
- M is the mass matrix
- ω are the natural frequencies
- φ are the mode shapes

### Solution Methods

1. Subspace iteration
2. Lanczos algorithm
3. Block Krylov methods

## Nonlinear Analysis

### Geometric Nonlinearity

Updated Lagrangian formulation with:

- Green-Lagrange strain
- Second Piola-Kirchhoff stress
- Total Lagrangian approach for large rotations

### Material Nonlinearity

1. Elastoplastic model:

   ```
   σ = Cᵉᶩ : (ε - εᵖ)
   ```

2. Damage model:

   ```
   σ = (1-D)C : ε
   ```

### Solution Methods

1. Newton-Raphson iteration:

   ```
   Kᵗδu = R
   ```

2. Arc-length method for limit points
3. Line search for improved convergence

## Thermal-Structural Coupling

### Governing Equations

```
Kᵐu + Kᵗδθ = F
Cθ̇ + Hθ = Q + Hᵉu̇
```

where:

- Kᵐ is mechanical stiffness
- Kᵗ is thermal stiffness
- H is conductivity matrix
- C is capacity matrix

### Solution Strategy

1. Staggered scheme
2. Monolithic coupling
3. Partitioned coupling

## Failure Criteria

### Tsai-Wu Criterion

```
F₁σ₁ + F₂σ₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂ = 1
```

### Puck Criterion

For fiber failure:

```
f_E = σ₁/X
```

For matrix failure:

```
f_E = √[(τₙₜ/(R^A_⊥⊥-μₙₜp^A_⊥⊥σₙ))² + (τₙₗ/R^A_⊥∥)²]
```

### LaRC Criterion

Matrix compression:

```
f_m = (τᵢ/S^L)² + (σₙ/Y_c)²
```

### Multi-Continuum Theory

1. Decompose composite stress into fiber and matrix stresses
2. Apply failure criteria to each constituent
3. Track damage evolution separately

## Optimization Methods

### Genetic Algorithm

1. Selection using tournament method
2. Crossover with adaptive rate
3. Mutation with dynamic probability
4. Elite preservation

### Particle Swarm Optimization

```
v_i(t+1) = wv_i(t) + c₁r₁(p_i-x_i(t)) + c₂r₂(g-x_i(t))
x_i(t+1) = x_i(t) + v_i(t+1)
```

### Surrogate-Based Optimization

1. Design of Experiments (Latin Hypercube)
2. Kriging metamodel construction
3. Expected Improvement criterion
4. Model update strategy

## Parallel Computing

### Domain Decomposition

1. Mesh partitioning using METIS
2. Interface handling with MPI
3. Load balancing strategies

### GPU Acceleration

1. Matrix assembly
2. Linear solver
3. Element calculations
4. Post-processing

### Task Parallelism

1. Task dependency graph
2. Dynamic scheduling
3. Load balancing
4. Memory management
