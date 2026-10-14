# Conformal and non-conformal contact

Course project comparing the normal contact stiffness of round and flat-faced particles with one spring-dashpot framework.

Sphere–sphere contact is the non-conformal case and is described by Hertz theory. Cube–cube face-to-face contact is the conformal case and is described by a flat-punch (Boussinesq) model. Both laws are evaluated in Python. The sphere–sphere law is then checked in LAMMPS with a linear spring-dashpot.

## Geometry and material

Both particles are 5 mm in size and share the same material, so the stiffness difference comes from the contact geometry.

| Quantity | Sphere–sphere | Cube–cube |
|---|---|---|
| Size | Radius $R = 5\ \mathrm{mm}$ | Side $L = 5\ \mathrm{mm}$ |
| Contact measure | $R^{*} = R/2 = 2.5\ \mathrm{mm}$ | Equal-area punch radius $a = 2.821\ \mathrm{mm}$ |
| Modulus | $E = 70\ \mathrm{GPa}$, $\nu = 0.22$, $E^{*} = 36.78\ \mathrm{GPa}$ | same |
| Overlap range | $0$–$2\ \mu\mathrm{m}$ | same |
| Force | $F \propto \delta^{3/2}$ | $F \propto \delta$ |
| Stiffness | $k \propto \sqrt{\delta}$ | constant |
| Representative $k_n$ | $5.202 \times 10^{6}\ \mathrm{N/m}$ at $2\ \mu\mathrm{m}$ | $2.075 \times 10^{8}\ \mathrm{N/m}$ |

The effective modulus for two identical bodies is

$$
E^{*} = \frac{E}{2\left(1-\nu^{2}\right)}.
$$

## Contact laws

Hertz force and tangent stiffness for two spheres:

$$
F_H = \frac{4}{3} E^{*} \sqrt{R^{*}}\, \delta^{3/2},
\qquad
k_H = \frac{dF_H}{d\delta} = 2 E^{*} \sqrt{R^{*} \delta}.
$$

Flat punch on the cube face. The square face is replaced by a circle of the same area, $a = \sqrt{A/\pi}$ with $A = L^{2}$:

$$
F_{\mathrm{flat}} = 2 E^{*} a\, \delta,
\qquad
k_{\mathrm{flat}} = 2 E^{*} a.
$$

At $2\ \mu\mathrm{m}$, $k_{\mathrm{flat}}/k_H \approx 39.9$. The ratio is larger at smaller overlap because the Hertz stiffness falls toward zero as contact begins, while the flat-punch stiffness does not.

The DEM spring is linear,

$$
F_n = k_n \delta + \gamma_n \dot{\delta}.
$$

For the sphere, $k_n$ is the Hertz tangent stiffness at $2\ \mu\mathrm{m}$. That choice matches the slope of the Hertz curve at one overlap, not the Hertz force. At $2\ \mu\mathrm{m}$ the linear force is $3/2$ times the Hertz force. With density $2500\ \mathrm{kg/m^{3}}$ and restitution $e = 0.9$,

$$
\gamma_n = -\frac{2 \ln(e) \sqrt{m_{\mathrm{eff}} k_n}}{\sqrt{\pi^{2} + \left[\ln(e)\right]^{2}}} = 3.911\ \mathrm{N\,s/m}
$$

for the sphere pair.

## What is in this repository

| Path | Contents |
|---|---|
| `contact_study.py` | Theory curves, LAMMPS inputs, runs, and comparison |
| `lammps/` | Sphere–sphere inputs: static hooke sweep, static Hertz sweep, head-on collision |
| `figures/` | Force–overlap, stiffness–overlap, stiffness-ratio, and LAMMPS comparison plots |
| `data/` | Parameters, theory tables, LAMMPS logs, and comparison tables |
| `RESULTS.txt` | Numerical summary and modeling notes |

LAMMPS uses `pair_style granular` with a hooke normal spring, `damping velocity`, and frictionless `linear_nohistory` tangential contact. A second static sweep uses `hertz/material` so the nonlinear Hertz force can be checked directly.

## Reproduce

Python 3 with NumPy and Matplotlib, and a LAMMPS build that includes the GRANULAR package. Point `LAMMPS_EXE` in `contact_study.py` at `lmp` if it is not the Windows binary used for these runs, then:

```bash
python contact_study.py
```

The script rewrites the inputs, reruns the three LAMMPS cases, and refreshes `figures/`, `data/`, and `RESULTS.txt`.

## Agreement with LAMMPS

On the static overlap sweep ($0.1$–$2.0\ \mu\mathrm{m}$, particles at rest):

- Hooke force versus $k_n \delta$: maximum relative error about $8 \times 10^{-15}$.
- `hertz/material` versus the Hertz formula: maximum relative error about $3 \times 10^{-15}$.
- At $2\ \mu\mathrm{m}$ the two LAMMPS models still differ by the expected factor $F_{\mathrm{linear}}/F_H = 1.5$.

The head-on collision of two free spheres reaches a peak overlap of $2.001\ \mu\mathrm{m}$ and a peak force of $10.43\ \mathrm{N}$. After the velocity-Verlet half-step is accounted for, the spring-dashpot force matches LAMMPS to about $10^{-12}$ relative error. The measured restitution is $0.8997$ against the target $0.90$, with a timestep of $10^{-7}\ \mathrm{s}$ (about 350 steps in contact).

Loading and unloading do not share one force–overlap line: damping raises the force while the overlap grows and lowers it while the overlap shrinks. Near separation the damping term can exceed the spring term, so the force is briefly tensile. The inputs do not set `limit_damping`.

The cube–cube flat-punch law is evaluated in Python only. This LAMMPS model is a sphere pair and does not represent a square face.

## References

- Johnson, K. L. (1985). *Contact Mechanics*. Cambridge University Press.
- Lu, G., Third, J. R., & Müller, C. R. (2015). Discrete element models for non-spherical particle systems. *Chemical Engineering Science*, 127, 425–465.
