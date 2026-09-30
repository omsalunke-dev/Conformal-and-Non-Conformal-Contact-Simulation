"""Sphere-sphere and cube-cube contact study.

Theoretical Hertz and flat-punch laws, linear spring-dashpot parameters,
LAMMPS sphere-sphere checks, and comparison against those laws.

Geometry and material data follow the EN 649 weekly progress report:
equal spheres of radius 5 mm, cubes of side 5 mm, E = 70 GPa, nu = 0.22.
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
LMPDIR = ROOT / "lammps"
LAMMPS_EXE = ROOT.parent / "bin" / "lmp.exe"

# ---------------------------------------------------------------------------
# Geometry and material system (report sections 2.1-2.5)
# ---------------------------------------------------------------------------
E = 70.0e9  # Pa
NU = 0.22
R = 5.0e-3  # m, sphere radius
L = 5.0e-3  # m, cube side
RHO = 2500.0  # kg/m^3
RESTITUTION = 0.9
DELTA_REF = 2.0e-6  # m, overlap at which the linear sphere spring is calibrated
DELTA_MAX = 2.0e-6
N_STATIC = 20  # static samples from 0.1 um to 2.0 um


def contact_parameters() -> dict:
    e_star = E / (2.0 * (1.0 - NU**2))
    r_star = R / 2.0
    area = L**2
    punch_radius = math.sqrt(area / math.pi)
    k_flat = 2.0 * e_star * punch_radius
    k_hertz_ref = 2.0 * e_star * math.sqrt(r_star * DELTA_REF)
    volume_sphere = (4.0 / 3.0) * math.pi * R**3
    mass_sphere = RHO * volume_sphere
    mass_eff = 0.5 * mass_sphere
    log_e = math.log(RESTITUTION)
    gamma_n = (
        -2.0
        * log_e
        * math.sqrt(mass_eff * k_hertz_ref)
        / math.sqrt(math.pi**2 + log_e**2)
    )
    # Approach speed that places the damped peak overlap near DELTA_REF.
    omega0 = math.sqrt(k_hertz_ref / mass_eff)
    zeta = gamma_n / (2.0 * mass_eff)
    omega_d = math.sqrt(omega0**2 - zeta**2)
    t_peak = math.pi / (2.0 * omega_d)
    v_rel = DELTA_REF * omega_d / math.exp(-zeta * t_peak)
    return {
        "E_Pa": E,
        "nu": NU,
        "E_star_Pa": e_star,
        "R_m": R,
        "R_star_m": r_star,
        "L_m": L,
        "face_area_m2": area,
        "punch_radius_m": punch_radius,
        "k_flat_Npm": k_flat,
        "k_sphere_Npm": k_hertz_ref,
        "delta_ref_m": DELTA_REF,
        "rho": RHO,
        "restitution": RESTITUTION,
        "mass_sphere_kg": mass_sphere,
        "mass_eff_kg": mass_eff,
        "gamma_n": gamma_n,
        "v_rel_mps": v_rel,
        "omega0": omega0,
        "contact_half_period_s": math.pi / omega_d,
    }


def hertz_force(delta: np.ndarray, p: dict) -> np.ndarray:
    return (4.0 / 3.0) * p["E_star_Pa"] * np.sqrt(p["R_star_m"]) * delta**1.5


def hertz_stiffness(delta: np.ndarray, p: dict) -> np.ndarray:
    return 2.0 * p["E_star_Pa"] * np.sqrt(p["R_star_m"] * delta)


def flat_force(delta: np.ndarray, p: dict) -> np.ndarray:
    return p["k_flat_Npm"] * delta


def flat_stiffness(delta: np.ndarray, p: dict) -> np.ndarray:
    return np.full_like(delta, p["k_flat_Npm"], dtype=float)


def write_theory_tables(p: dict) -> None:
    delta = np.linspace(0.0, DELTA_MAX, 401)
    # Stiffness and the ratio are singular or zero at exact contact; keep a
    # copy that starts just above zero for the ratio curve.
    delta_pos = delta.copy()
    delta_pos[0] = delta[1] * 0.5

    f_h = hertz_force(delta, p)
    k_h = hertz_stiffness(delta_pos, p)
    k_h[0] = 0.0
    f_f = flat_force(delta, p)
    k_f = flat_stiffness(delta, p)
    f_lin_s = p["k_sphere_Npm"] * delta
    f_lin_c = p["k_flat_Npm"] * delta
    ratio = np.divide(k_f, k_h, out=np.full_like(k_h, np.nan), where=k_h > 0)

    curve_path = DATA / "theory_curves.csv"
    with curve_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "delta_m",
                "delta_um",
                "F_hertz_N",
                "F_flat_N",
                "F_linear_sphere_N",
                "F_linear_cube_N",
                "k_hertz_Npm",
                "k_flat_Npm",
                "k_linear_sphere_Npm",
                "k_linear_cube_Npm",
                "stiffness_ratio_flat_over_hertz",
            ]
        )
        for i in range(delta.size):
            writer.writerow(
                [
                    f"{delta[i]:.8e}",
                    f"{delta[i] * 1e6:.8e}",
                    f"{f_h[i]:.8e}",
                    f"{f_f[i]:.8e}",
                    f"{f_lin_s[i]:.8e}",
                    f"{f_lin_c[i]:.8e}",
                    f"{k_h[i]:.8e}",
                    f"{k_f[i]:.8e}",
                    f"{p['k_sphere_Npm']:.8e}",
                    f"{p['k_flat_Npm']:.8e}",
                    "" if not np.isfinite(ratio[i]) else f"{ratio[i]:.8e}",
                ]
            )

    param_path = DATA / "parameters.csv"
    with param_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["quantity", "value", "unit"])
        rows = [
            ("Young modulus E", E, "Pa"),
            ("Poisson ratio nu", NU, "-"),
            ("Effective modulus E*", p["E_star_Pa"], "Pa"),
            ("Sphere radius R", R, "m"),
            ("Effective radius R*", p["R_star_m"], "m"),
            ("Cube side L", L, "m"),
            ("Cube face area", p["face_area_m2"], "m^2"),
            ("Equal-area punch radius a", p["punch_radius_m"], "m"),
            ("Hertz stiffness at 2 um", p["k_sphere_Npm"], "N/m"),
            ("Flat-punch stiffness", p["k_flat_Npm"], "N/m"),
            ("Stiffness ratio at 2 um", p["k_flat_Npm"] / p["k_sphere_Npm"], "-"),
            ("Density", RHO, "kg/m^3"),
            ("Coefficient of restitution", RESTITUTION, "-"),
            ("Sphere mass", p["mass_sphere_kg"], "kg"),
            ("Sphere-sphere reduced mass", p["mass_eff_kg"], "kg"),
            ("Normal damping gamma_n", p["gamma_n"], "N s/m"),
            ("Collision approach speed", p["v_rel_mps"], "m/s"),
        ]
        for name, value, unit in rows:
            writer.writerow([name, f"{value:.10e}", unit])


def _style() -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "axes.grid": True,
            "grid.alpha": 0.35,
        }
    )


def make_theory_plots(p: dict) -> None:
    _style()
    delta_um = np.linspace(0.0, 2.0, 400)
    delta = delta_um * 1e-6
    delta_k = np.clip(delta, 1e-12, None)

    f_h = hertz_force(delta, p)
    f_f = flat_force(delta, p)
    f_ls = p["k_sphere_Npm"] * delta
    f_lc = p["k_flat_Npm"] * delta
    k_h = hertz_stiffness(delta_k, p)
    k_h[0] = 0.0
    k_f = flat_stiffness(delta, p)
    ratio = p["k_flat_Npm"] / hertz_stiffness(delta_k, p)

    # Force-overlap
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(delta_um, f_f, color="#c45911", lw=2.0, label="Flat punch, cube-cube")
    ax.plot(
        delta_um,
        f_lc,
        color="#ed7d31",
        lw=1.6,
        ls="--",
        label="Linear spring, cube $k_n$",
    )
    ax.plot(delta_um, f_h, color="#1f4e79", lw=2.0, label="Hertz, sphere-sphere")
    ax.plot(
        delta_um,
        f_ls,
        color="#5b9bd5",
        lw=1.6,
        ls="--",
        label="Linear spring, sphere $k_n$",
    )
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel("Normal force (N)")
    ax.set_title("Force-overlap: Hertz, flat punch, and linear springs")
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, None)
    ax.legend(loc="upper left")
    inset = ax.inset_axes([0.38, 0.18, 0.55, 0.38])
    inset.plot(delta_um, f_h, color="#1f4e79", lw=2.0)
    inset.plot(delta_um, f_ls, color="#5b9bd5", lw=1.6, ls="--")
    inset.set_xlim(0.0, 2.0)
    inset.set_ylim(0.0, max(f_ls) * 1.08)
    inset.set_title("Sphere curves", fontsize=9)
    inset.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(FIG / "force_overlap.png")
    plt.close(fig)

    # Stiffness-overlap
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(delta_um, k_f / 1e6, color="#c45911", lw=2.0, label="Flat punch, cube-cube")
    ax.plot(
        delta_um,
        np.full_like(delta_um, p["k_flat_Npm"]) / 1e6,
        color="#ed7d31",
        lw=1.6,
        ls="--",
        label="Linear spring, cube $k_n$",
    )
    ax.plot(delta_um, k_h / 1e6, color="#1f4e79", lw=2.0, label="Hertz tangent, sphere-sphere")
    ax.plot(
        delta_um,
        np.full_like(delta_um, p["k_sphere_Npm"]) / 1e6,
        color="#5b9bd5",
        lw=1.6,
        ls="--",
        label=r"Linear spring, sphere $k_n$ at 2 $\mu$m",
    )
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel(r"Normal stiffness ($10^6$ N/m)")
    ax.set_title("Stiffness-overlap: Hertz tangent stiffness and flat punch")
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, p["k_flat_Npm"] / 1e6 * 1.08)
    ax.legend(loc="center right")
    inset = ax.inset_axes([0.16, 0.45, 0.48, 0.42])
    inset.plot(delta_um, k_h / 1e6, color="#1f4e79", lw=2.0, label="Hertz")
    inset.axhline(p["k_sphere_Npm"] / 1e6, color="#5b9bd5", lw=1.6, ls="--", label="Linear $k_n$")
    inset.set_xlim(0.0, 2.0)
    inset.set_ylim(0.0, p["k_sphere_Npm"] / 1e6 * 1.15)
    inset.set_title("Sphere stiffness", fontsize=9)
    inset.grid(True, alpha=0.35)
    fig.tight_layout()
    fig.savefig(FIG / "stiffness_overlap.png")
    plt.close(fig)

    # Stiffness ratio
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    mask = delta_um >= 0.05
    ax.plot(
        delta_um[mask],
        ratio[mask],
        color="#1f4e79",
        lw=2.2,
        label=r"$k_\mathrm{flat} / k_\mathrm{Hertz}(\delta)$",
    )
    ratio_lin = p["k_flat_Npm"] / p["k_sphere_Npm"]
    ax.axhline(
        ratio_lin,
        color="#c45911",
        lw=1.6,
        ls="--",
        label=rf"Linear ratio $k_{{n,cube}}/k_{{n,sphere}} = {ratio_lin:.2f}$",
    )
    ax.scatter(
        [2.0],
        [ratio_lin],
        color="#c45911",
        zorder=3,
        label=rf"Report check at 2 $\mu$m: {ratio[-1]:.2f}",
    )
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel(r"Stiffness ratio $k_\mathrm{cube} / k_\mathrm{sphere}$")
    ax.set_title("Stiffness ratio of flat-faced contact to spherical contact")
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, None)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "stiffness_ratio.png")
    plt.close(fig)


def write_lammps_inputs(p: dict) -> None:
    kn = f"{p['k_sphere_Npm']:.12e}"
    gn = f"{p['gamma_n']:.12e}"
    v_half = f"{0.5 * p['v_rel_mps']:.12e}"
    diameter = f"{2.0 * R:.12e}"
    e_mod = f"{E:.12e}"

    static = f"""# Quasi-static sphere-sphere overlap sweep
# One particle pair, zero velocity. MODEL selects hooke or hertz.
units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
atom_modify     map array
newton          on

region          box block -0.02 0.04 -0.02 0.02 -0.02 0.02
create_box      1 box

variable        model string hooke
pair_style      granular
pair_coeff      * * hooke {kn} {gn} tangential linear_nohistory 0.0 0.0 damping velocity

neighbor        0.002 bin
neigh_modify    delay 0 every 1 check yes

create_atoms    1 single 0.0 0.0 0.0
create_atoms    1 single 0.009999 0.0 0.0
set             atom 1*2 diameter {diameter}
set             atom 1*2 density {RHO}
set             atom * vx 0.0 vy 0.0 vz 0.0

variable        i loop {N_STATIC}
variable        delta equal v_i*1.0e-7
variable        xpos equal 0.01-v_delta
variable        fx1 equal fx[1]
variable        x1 equal x[1]
variable        x2a equal x[2]

thermo          1
thermo_style    custom step pe fmax

label           loop
reset_timestep  $i
set             atom 1 x 0.0 y 0.0 z 0.0 vx 0.0 vy 0.0 vz 0.0
set             atom 2 x ${{xpos}} y 0.0 z 0.0 vx 0.0 vy 0.0 vz 0.0
run             0
print           "$i ${{delta}} ${{x1}} ${{x2a}} ${{fx1}}" append ../data/static_MODEL.dat
next            i
jump            SELF loop
"""
    # The print line must expand on each jump. MODEL is replaced per file.
    hooke = static.replace("static_MODEL.dat", "static_hooke.dat")
    hertz = static.replace("static_MODEL.dat", "static_hertz.dat")
    hertz = hertz.replace(
        f"pair_coeff      * * hooke {kn} {gn} tangential linear_nohistory 0.0 0.0 damping velocity",
        f"pair_coeff      * * hertz/material {e_mod} 0.0 {NU} tangential linear_nohistory 0.0 0.0 damping velocity",
    )
    (LMPDIR / "in.static_hooke.lmp").write_text(hooke, newline="\n")
    (LMPDIR / "in.static_hertz.lmp").write_text(hertz, newline="\n")

    collision = f"""# Head-on sphere-sphere collision with the linear spring-dashpot law
units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
atom_modify     map array
newton          on

region          box block -0.05 0.06 -0.02 0.02 -0.02 0.02
create_box      1 box

pair_style      granular
pair_coeff      * * hooke {kn} {gn} tangential linear_nohistory 0.0 0.0 damping velocity

neighbor        0.005 bin
neigh_modify    delay 0 every 1 check yes

# Initial gap of 2 um so contact begins after the neighbour list is active.
create_atoms    1 single 0.0 0.0 0.0
create_atoms    1 single 0.010002 0.0 0.0
set             atom 1*2 diameter {diameter}
set             atom 1*2 density {RHO}
set             atom 1 vx  {v_half} vy 0.0 vz 0.0
set             atom 2 vx -{v_half} vy 0.0 vz 0.0

fix             1 all nve/sphere
timestep        1.0e-7

dump            1 all custom 1 ../data/collision.dump id x vx fx radius
dump_modify     1 sort id format float %.16e

thermo_style    custom step ke pe etotal
thermo          50
run             1000
"""
    (LMPDIR / "in.collision.lmp").write_text(collision, newline="\n")


def run_lammps(script: str, log_name: str) -> None:
    if not LAMMPS_EXE.exists():
        raise FileNotFoundError(f"LAMMPS executable not found: {LAMMPS_EXE}")
    log_path = DATA / log_name
    cmd = [str(LAMMPS_EXE), "-in", script, "-log", str(log_path)]
    print("Running", " ".join(cmd))
    completed = subprocess.run(
        cmd,
        cwd=LMPDIR,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(f"LAMMPS failed on {script} (exit {completed.returncode})")


def load_static(path: Path) -> dict[str, np.ndarray]:
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        rows.append([float(x) for x in parts[:5]])
    arr = np.array(rows, dtype=float)
    return {
        "index": arr[:, 0],
        "delta_cmd": arr[:, 1],
        "x1": arr[:, 2],
        "x2": arr[:, 3],
        "fx1": arr[:, 4],
    }


def load_dump(path: Path) -> dict[str, np.ndarray]:
    steps = []
    records = []
    step = None
    natoms = 0
    reading = False
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if line == "ITEM: TIMESTEP":
                step = int(next(fh).strip())
                reading = False
                continue
            if line == "ITEM: NUMBER OF ATOMS":
                natoms = int(next(fh).strip())
                continue
            if line.startswith("ITEM: ATOMS"):
                block = []
                for _ in range(natoms):
                    block.append(next(fh).split())
                steps.append(step)
                records.append(block)
                reading = False
    # columns: id x vx fx radius
    n = len(records)
    out = {
        "step": np.array(steps, dtype=float),
        "time": np.array(steps, dtype=float) * 1.0e-7,
    }
    for key in ("x1", "x2", "vx1", "vx2", "fx1", "fx2", "r1", "r2"):
        out[key] = np.empty(n)
    for i, block in enumerate(records):
        atoms = {int(row[0]): row for row in block}
        a1, a2 = atoms[1], atoms[2]
        out["x1"][i] = float(a1[1])
        out["vx1"][i] = float(a1[2])
        out["fx1"][i] = float(a1[3])
        out["r1"][i] = float(a1[4])
        out["x2"][i] = float(a2[1])
        out["vx2"][i] = float(a2[2])
        out["fx2"][i] = float(a2[3])
        out["r2"][i] = float(a2[4])
    return out


def compare(p: dict) -> dict:
    _style()
    radius = R
    hooke = load_static(DATA / "static_hooke.dat")
    hertz = load_static(DATA / "static_hertz.dat")

    def overlap_from_x(data: dict) -> np.ndarray:
        return (2.0 * radius) - (data["x2"] - data["x1"])

    hooke_delta = overlap_from_x(hooke)
    hertz_delta = overlap_from_x(hertz)
    # Repulsive magnitude: left particle is pushed in -x.
    f_hooke = -hooke["fx1"]
    f_hertz = -hertz["fx1"]
    f_lin = p["k_sphere_Npm"] * hooke_delta
    f_hz_theory = hertz_force(hertz_delta, p)

    rel_hooke = np.abs(f_hooke - f_lin) / f_lin
    rel_hertz = np.abs(f_hertz - f_hz_theory) / f_hz_theory

    with (DATA / "static_comparison.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "delta_m",
                "F_lammps_hooke_N",
                "F_linear_N",
                "rel_err_hooke",
                "F_lammps_hertz_N",
                "F_hertz_theory_N",
                "rel_err_hertz",
                "F_linear_over_F_hertz",
            ]
        )
        for i in range(hooke_delta.size):
            writer.writerow(
                [
                    f"{hooke_delta[i]:.8e}",
                    f"{f_hooke[i]:.8e}",
                    f"{f_lin[i]:.8e}",
                    f"{rel_hooke[i]:.8e}",
                    f"{f_hertz[i]:.8e}",
                    f"{f_hz_theory[i]:.8e}",
                    f"{rel_hertz[i]:.8e}",
                    f"{(f_lin[i] / f_hz_theory[i]):.8e}",
                ]
            )

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    d_um = hooke_delta * 1e6
    ax.plot(d_um, f_lin, color="#1f4e79", lw=2.0, label="Linear theory $k_n\\delta$")
    ax.plot(
        d_um,
        f_hooke,
        color="#c45911",
        ls="none",
        marker="o",
        ms=5,
        label="LAMMPS hooke",
    )
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel("Repulsive force (N)")
    ax.set_title("Static sphere-sphere check: LAMMPS linear spring vs $k_n\\delta$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "lammps_static_hooke.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    d_um = hertz_delta * 1e6
    ax.plot(d_um, f_hz_theory, color="#1f4e79", lw=2.0, label="Hertz theory")
    ax.plot(d_um, p["k_sphere_Npm"] * hertz_delta, color="#5b9bd5", ls="--", lw=1.6, label="Linear $k_n\\delta$")
    ax.plot(d_um, f_hertz, color="#c45911", ls="none", marker="o", ms=5, label="LAMMPS hertz/material")
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel("Repulsive force (N)")
    ax.set_title("Static sphere-sphere check: Hertz law vs linear spring")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "lammps_static_hertz.png")
    plt.close(fig)

    col = load_dump(DATA / "collision.dump")
    dt = 1.0e-7
    delta = (col["r1"] + col["r2"]) - (col["x2"] - col["x1"])
    # Positive when the gap is closing. The dump stores the velocity after the
    # second velocity-Verlet kick. Pair damping was evaluated one half-kick
    # earlier, so the end-of-step rate is low by (F/m)*dt for equal masses.
    delta_dot_end = -(col["vx2"] - col["vx1"])
    contact = delta > 0.0
    f_meas = -col["fx1"]
    f_from_end = p["k_sphere_Npm"] * delta + p["gamma_n"] * delta_dot_end
    half_kick = p["gamma_n"] * dt / p["mass_sphere_kg"]
    f_pred = np.where(contact, f_from_end / (1.0 - half_kick), 0.0)
    delta_dot = delta_dot_end + np.where(contact, f_meas * dt / p["mass_sphere_kg"], 0.0)

    # Approach and separation speeds from the free-flight plateaus.
    incoming = col["time"] < col["time"][contact][0]
    outgoing_start = np.where(contact)[0][-1] + 1
    v_in = np.median(-(col["vx2"][incoming] - col["vx1"][incoming]))
    if outgoing_start < len(delta):
        v_out = np.median(col["vx2"][outgoing_start:] - col["vx1"][outgoing_start:])
    else:
        v_out = float("nan")
    e_meas = v_out / v_in if v_in else float("nan")

    rel_contact = np.abs(f_meas[contact] - f_pred[contact]) / np.maximum(np.abs(f_pred[contact]), 1e-12)
    rel_end = np.abs(f_meas[contact] - f_from_end[contact]) / np.maximum(
        np.abs(f_from_end[contact]), 1e-12
    )

    with (DATA / "collision_comparison.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "time_s",
                "delta_m",
                "delta_dot_mps",
                "F_lammps_N",
                "F_spring_dashpot_N",
                "in_contact",
            ]
        )
        for i in range(delta.size):
            writer.writerow(
                [
                    f"{col['time'][i]:.8e}",
                    f"{delta[i]:.8e}",
                    f"{delta_dot[i]:.8e}",
                    f"{f_meas[i]:.8e}",
                    f"{f_pred[i]:.8e}",
                    int(contact[i]),
                ]
            )

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(col["time"][contact] * 1e6, f_pred[contact], color="#1f4e79", lw=2.0, label="Spring-dashpot $k_n\\delta + \\gamma_n\\dot\\delta$")
    # Plot every 4th LAMMPS point so the line remains visible.
    idx = np.where(contact)[0][::4]
    ax.plot(
        col["time"][idx] * 1e6,
        f_meas[idx],
        color="#c45911",
        ls="none",
        marker="o",
        ms=3.5,
        label="LAMMPS",
    )
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel("Repulsive force (N)")
    ax.set_title("Head-on collision: LAMMPS force vs spring-dashpot")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "lammps_collision_force_time.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    loading = contact & (delta_dot_end > 0.0)
    unloading = contact & (delta_dot_end <= 0.0)
    ax.plot(
        delta[loading] * 1e6,
        f_meas[loading],
        color="#c45911",
        lw=2.0,
        label="LAMMPS, loading",
    )
    ax.plot(
        delta[unloading] * 1e6,
        f_meas[unloading],
        color="#1f4e79",
        lw=2.0,
        label="LAMMPS, unloading",
    )
    d_line = np.linspace(0.0, float(np.max(delta[contact])), 200)
    ax.plot(
        d_line * 1e6,
        p["k_sphere_Npm"] * d_line,
        color="#5b9bd5",
        ls="--",
        lw=1.5,
        label="Elastic part only, $k_n\\delta$",
    )
    ax.set_xlabel(r"Overlap $\delta$ ($\mu$m)")
    ax.set_ylabel("Repulsive force (N)")
    ax.set_title("Collision force-overlap: damping splits loading and unloading")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "lammps_collision_force_overlap.png")
    plt.close(fig)

    stats = {
        "hooke_max_rel": float(np.max(rel_hooke)),
        "hooke_rms_rel": float(np.sqrt(np.mean(rel_hooke**2))),
        "hertz_max_rel": float(np.max(rel_hertz)),
        "hertz_rms_rel": float(np.sqrt(np.mean(rel_hertz**2))),
        "collision_max_rel": float(np.max(rel_contact)),
        "collision_rms_rel": float(np.sqrt(np.mean(rel_contact**2))),
        "collision_end_step_rel": float(np.median(rel_end)),
        "half_kick": float(half_kick),
        "v_in": float(v_in),
        "v_out": float(v_out),
        "e_meas": float(e_meas),
        "delta_max": float(np.max(delta)),
        "f_max": float(np.max(f_meas)),
        "n_contact": int(np.count_nonzero(contact)),
        "force_ratio_at_2um": float(f_lin[-1] / f_hz_theory[-1]),
    }
    return stats


def write_summary(p: dict, stats: dict) -> None:
    ratio = p["k_flat_Npm"] / p["k_sphere_Npm"]
    f_h_ref = (4.0 / 3.0) * p["E_star_Pa"] * math.sqrt(p["R_star_m"]) * DELTA_REF**1.5
    f_lin_ref = p["k_sphere_Npm"] * DELTA_REF
    text = f"""Sphere-sphere contact study
==========================

Geometry and material
---------------------
Sphere radius R = {R*1e3:.3f} mm, so the Hertz effective radius is R* = R/2 = {p['R_star_m']*1e3:.3f} mm.
Cube side L = {L*1e3:.3f} mm, face area A = {p['face_area_m2']*1e6:.3f} mm^2.
Equal-area circular punch radius a = sqrt(A/pi) = {p['punch_radius_m']*1e3:.4f} mm.
E = {E/1e9:.2f} GPa, nu = {NU:.2f}.
E* = E / (2(1-nu^2)) = {p['E_star_Pa']/1e9:.6f} GPa.

Contact laws
------------
Hertz force:        F = (4/3) E* sqrt(R*) delta^(3/2)
Hertz stiffness:    k = dF/ddelta = 2 E* sqrt(R* delta)
Flat-punch force:   F = 2 E* a delta
Flat-punch stiffness is constant: k_flat = 2 E* a = {p['k_flat_Npm']:.6e} N/m.

At delta = 2 um the Hertz tangent stiffness is k_H = {p['k_sphere_Npm']:.6e} N/m.
The stiffness ratio k_flat / k_H(2 um) = {ratio:.4f}.
The report values 5.20e6 N/m, 2.075e8 N/m and a ratio of about 39.9 are the same results rounded.

Linear spring-dashpot used in LAMMPS
-------------------------------------
Sphere spring: k_n = k_H(2 um) = {p['k_sphere_Npm']:.6e} N/m.
This matches the tangent stiffness at 2 um. It does not match the Hertz force.
At 2 um, F_linear / F_Hertz = {f_lin_ref/f_h_ref:.6f} (exactly 3/2 for this calibration).
Density = {RHO:.1f} kg/m^3, sphere mass = {p['mass_sphere_kg']:.6e} kg,
reduced mass = {p['mass_eff_kg']:.6e} kg, e = {RESTITUTION:.2f}.
Damping coefficient gamma_n = {p['gamma_n']:.6e} N s/m,
from gamma_n = -2 ln(e) sqrt(m_eff k_n) / sqrt(pi^2 + ln(e)^2).
LAMMPS pair style: granular, normal model hooke, damping velocity,
tangential linear_nohistory with zero friction.
With damping velocity, the LAMMPS coefficient is gamma_n itself:
F = k_n delta + gamma_n * d(delta)/dt during contact.

LAMMPS tests
------------
1. Static overlap sweep, 0.1 um to 2.0 um in steps of 0.1 um, both particles at rest.
   Hooke model versus k_n * delta:
     max relative error = {stats['hooke_max_rel']:.3e}
     rms relative error = {stats['hooke_rms_rel']:.3e}
   Hertz/material model (E = 70 GPa, nu = 0.22) versus the Hertz formula:
     max relative error = {stats['hertz_max_rel']:.3e}
     rms relative error = {stats['hertz_rms_rel']:.3e}
   At 2 um, F_linear / F_Hertz from the two LAMMPS models = {stats['force_ratio_at_2um']:.6f}.

2. Head-on collision of two free spheres. Initial gap = 2 um.
   Prescribed approach speed = {p['v_rel_mps']:.6e} m/s.
   Measured incoming speed = {stats['v_in']:.6e} m/s.
   Measured separating speed = {stats['v_out']:.6e} m/s.
   Measured restitution = {stats['e_meas']:.6f} (target {RESTITUTION:.2f}).
   Peak overlap = {stats['delta_max']*1e6:.4f} um.
   Peak repulsive force = {stats['f_max']:.6f} N.
   Contact samples = {stats['n_contact']}.
   Force versus k_n*delta + gamma_n*delta_dot, using the half-step
   velocity at which LAMMPS evaluates damping:
     max relative error = {stats['collision_max_rel']:.3e}
     rms relative error = {stats['collision_rms_rel']:.3e}
   The same formula with the dumped end-of-step velocity is low by a
   uniform factor. Median relative offset = {stats['collision_end_step_rel']:.3e},
   which matches gamma_n*dt/m = {stats['half_kick']:.3e}.

Modeling differences
--------------------
- The LAMMPS hooke run reproduces the linear law. The largest relative force
  error on the static sweep is about 8e-15, which is floating-point roundoff.
- Hertz force grows as delta^(3/2). The linear spring calibrated to the tangent
  stiffness at 2 um is stiffer, in the secant sense F/delta, at every overlap
  in this range, and its force is 1.5 times the Hertz force at 2 um.
- Loading and unloading in the collision do not follow a single force-overlap
  curve. Damping adds force while the overlap is increasing and subtracts force
  while the overlap is decreasing.
- LAMMPS writes the velocity after the second velocity-Verlet kick. The pair
  force was computed with the velocity from half a step earlier. Comparing the
  dumped velocity directly to k_n*delta + gamma_n*delta_dot therefore shifts
  every contact force by about gamma_n*dt/m, here 0.030%. Removing that
  half-step leaves a residual at roundoff, so the contact law itself matches.
- The measured restitution, 0.8997 against the target 0.90, is the same
  timestep effect. The continuous-time damping formula is exact only as dt
  goes to zero. This run uses 1e-7 s, about 350 steps of contact.
- Near the end of contact the damping term can exceed the spring term, so the
  pair force becomes tensile for a short interval. The input does not set
  limit_damping, which is the LAMMPS option that removes that tensile spike.
- The cube-cube flat-punch law is carried in the Python curves only. This
  LAMMPS model is a sphere-sphere pair; a square face is not represented.
- The flat-punch stiffness uses an equal-area circle of radius {p['punch_radius_m']*1e3:.4f} mm
  in place of the square face. That geometric substitution is unchanged from
  the weekly report.

Outputs
-------
figures/force_overlap.png
figures/stiffness_overlap.png
figures/stiffness_ratio.png
figures/lammps_static_hooke.png
figures/lammps_static_hertz.png
figures/lammps_collision_force_time.png
figures/lammps_collision_force_overlap.png
data/theory_curves.csv
data/parameters.csv
data/static_comparison.csv
data/collision_comparison.csv
lammps/in.static_hooke.lmp
lammps/in.static_hertz.lmp
lammps/in.collision.lmp
"""
    (ROOT / "RESULTS.txt").write_text(text, newline="\n")
    print(text)


def main() -> None:
    DATA.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    LMPDIR.mkdir(exist_ok=True)
    for name in ("static_hooke.dat", "static_hertz.dat"):
        target = DATA / name
        if target.exists():
            target.unlink()

    p = contact_parameters()
    write_theory_tables(p)
    make_theory_plots(p)
    write_lammps_inputs(p)
    run_lammps("in.static_hooke.lmp", "log.static_hooke")
    run_lammps("in.static_hertz.lmp", "log.static_hertz")
    run_lammps("in.collision.lmp", "log.collision")
    stats = compare(p)
    write_summary(p, stats)


if __name__ == "__main__":
    main()
