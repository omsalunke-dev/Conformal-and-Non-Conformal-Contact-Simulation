# Sphere–sphere contact study

## Geometry and material

- Sphere radius $R = 5.000\,\mathrm{mm}$, so $R^{\ast} = R/2 = 2.500\,\mathrm{mm}$.
- Cube side $L = 5.000\,\mathrm{mm}$, face area $A = 25.000\,\mathrm{mm}^2$.
- Equal-area punch radius $a = \sqrt{A/\pi} = 2.8209\,\mathrm{mm}$.
- $E = 70.00\,\mathrm{GPa}$, $\nu = 0.22$.

```math
E^{\ast} = \frac{E}{2\left(1-\nu^{2}\right)} = 36.780160\,\mathrm{GPa}.
```

## Contact laws

```math
F_{H} = \frac{4}{3} E^{\ast} \sqrt{R^{\ast}} \, \delta^{3/2}
```

```math
k_{H} = \frac{dF_{H}}{d\delta} = 2 E^{\ast} \sqrt{R^{\ast} \delta}
```

```math
F_{\mathrm{flat}} = 2 E^{\ast} a \, \delta
```

```math
k_{\mathrm{flat}} = 2 E^{\ast} a = 2.075098 \times 10^{8}\,\mathrm{N/m}
```

At $\delta = 2\,\mu\mathrm{m}$ the Hertz tangent stiffness is $k_{H} = 5.201500 \times 10^{6}\,\mathrm{N/m}$.

The stiffness ratio is $k_{\mathrm{flat}} / k_{H}(2\,\mu\mathrm{m}) = 39.8942$.

The rounded report values $5.20 \times 10^{6}\,\mathrm{N/m}$, $2.075 \times 10^{8}\,\mathrm{N/m}$, and a ratio of about $39.9$ are the same results.

## Linear spring-dashpot used in LAMMPS

The sphere spring is the Hertz tangent stiffness at $2\,\mu\mathrm{m}$:

```math
k_{n} = k_{H}(2\,\mu\mathrm{m}) = 5.201500 \times 10^{6}\,\mathrm{N/m}
```

This matches the slope of the Hertz curve at that overlap, not the Hertz force. At $2\,\mu\mathrm{m}$,

```math
F_{\mathrm{linear}} / F_{H} = 1.500000,
```

which is exactly $3/2$ for this calibration.

| Quantity | Value |
|---|---|
| Density | $2500.0\,\mathrm{kg/m^3}$ |
| Sphere mass | $1.308997 \times 10^{-3}\,\mathrm{kg}$ |
| Reduced mass | $6.544985 \times 10^{-4}\,\mathrm{kg}$ |
| Restitution $e$ | $0.90$ |
| Damping $\gamma_{n}$ | $3.911403\,\mathrm{N\,s/m}$ |

```math
\gamma_{n} = -\frac{2 \ln(e) \sqrt{m_{\mathrm{eff}} k_{n}}}{\sqrt{\pi^{2} + \left[\ln(e)\right]^{2}}}
```

LAMMPS uses `pair_style granular`, normal model hooke, `damping velocity`, and frictionless `linear_nohistory` tangential contact. With `damping velocity` the coefficient is $\gamma_{n}$ itself:

```math
F = k_{n} \delta + \gamma_{n} \dot{\delta}
```

during contact.

## LAMMPS tests

### Static overlap sweep

Overlaps run from $0.1\,\mu\mathrm{m}$ to $2.0\,\mu\mathrm{m}$ in steps of $0.1\,\mu\mathrm{m}$, with both particles at rest.

| Comparison | Maximum relative error | RMS relative error |
|---|---|---|
| Hooke versus $k_{n}\delta$ | $7.684 \times 10^{-15}$ | $3.810 \times 10^{-15}$ |
| `hertz/material` versus Hertz, $E = 70\,\mathrm{GPa}$, $\nu = 0.22$ | $2.783 \times 10^{-15}$ | $1.183 \times 10^{-15}$ |

At $2\,\mu\mathrm{m}$, $F_{\mathrm{linear}} / F_{H} = 1.500000$ from the two LAMMPS models.

### Head-on collision

Two free spheres start with a gap of $2\,\mu\mathrm{m}$.

| Quantity | Value |
|---|---|
| Prescribed approach speed | $0.187834\,\mathrm{m/s}$ |
| Measured incoming speed | $0.187834\,\mathrm{m/s}$ |
| Measured separating speed | $0.169001\,\mathrm{m/s}$ |
| Measured restitution | $0.899735$ (target $0.90$) |
| Peak overlap | $2.0008\,\mu\mathrm{m}$ |
| Peak repulsive force | $10.433828\,\mathrm{N}$ |
| Contact samples | 353 |

Force compared with $k_{n}\delta + \gamma_{n}\dot{\delta}$, using the half-step velocity at which LAMMPS evaluates damping:

- Maximum relative error $1.235 \times 10^{-12}$
- RMS relative error $1.032 \times 10^{-13}$

The same formula with the dumped end-of-step velocity is low by a uniform factor. The median relative offset is $2.989 \times 10^{-4}$, which matches $\gamma_{n} \Delta t / m = 2.988 \times 10^{-4}$.

## Modeling differences

- The LAMMPS hooke run reproduces the linear law. The largest relative force error on the static sweep is about $8 \times 10^{-15}$, which is floating-point roundoff.
- Hertz force grows as $\delta^{3/2}$. The linear spring calibrated to the tangent stiffness at $2\,\mu\mathrm{m}$ is stiffer, in the secant sense $F/\delta$, at every overlap in this range, and its force is $1.5$ times the Hertz force at $2\,\mu\mathrm{m}$.
- Loading and unloading in the collision do not follow a single force–overlap curve. Damping adds force while the overlap is increasing and subtracts force while the overlap is decreasing.
- LAMMPS writes the velocity after the second velocity-Verlet kick. The pair force was computed with the velocity from half a step earlier. Comparing the dumped velocity directly to $k_{n}\delta + \gamma_{n}\dot{\delta}$ therefore shifts every contact force by about $\gamma_{n}\Delta t/m$, here $0.030\%$. Removing that half-step leaves a residual at roundoff, so the contact law itself matches.
- The measured restitution, $0.8997$ against the target $0.90$, is the same timestep effect. The continuous-time damping formula is exact only as $\Delta t$ goes to zero. This run uses $10^{-7}\,\mathrm{s}$, about 350 steps of contact.
- Near the end of contact the damping term can exceed the spring term, so the pair force becomes tensile for a short interval. The input does not set `limit_damping`, which is the LAMMPS option that removes that tensile spike.
- The cube–cube flat-punch law is carried in the Python curves only. This LAMMPS model is a sphere–sphere pair; a square face is not represented.
- The flat-punch stiffness uses an equal-area circle of radius $2.8209\,\mathrm{mm}$ in place of the square face. That geometric substitution is unchanged from the weekly report.

## Outputs

- `figures/force_overlap.png`
- `figures/stiffness_overlap.png`
- `figures/stiffness_ratio.png`
- `figures/lammps_static_hooke.png`
- `figures/lammps_static_hertz.png`
- `figures/lammps_collision_force_time.png`
- `figures/lammps_collision_force_overlap.png`
- `data/theory_curves.csv`
- `data/parameters.csv`
- `data/static_comparison.csv`
- `data/collision_comparison.csv`
- `lammps/in.static_hooke.lmp`
- `lammps/in.static_hertz.lmp`
- `lammps/in.collision.lmp`
