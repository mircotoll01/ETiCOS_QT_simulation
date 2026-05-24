## CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Academic project for the *Quantum Transport* course at Politecnico di Torino (a.y. 2025-2026). It contains:

- A MATLAB simulation engine in [data/](data/) that solves 1D quantum transport in layered heterostructures (wells + barriers).
- A Beamer LaTeX report in [qt_report.tex](qt_report.tex) that includes per-topic frames from [sections/](sections/), each embedding a PNG produced by the MATLAB scripts in [images/](images/).

The two halves are decoupled: MATLAB scripts produce figures, the report includes them by filename. There is no build pipeline wiring them together — images are committed.

## Building the report

```powershell
pdflatex qt_report.tex
pdflatex qt_report.tex    # second pass for \tableofcontents
```

The preamble sets `\graphicspath{{images/}}`, so `\includegraphics{Foo.png}` resolves to `images/Foo.png`. New simulation figures must be saved into `images/` and referenced by bare filename in the corresponding `sections/*.tex`.

## Running simulations (MATLAB)

Open MATLAB in [data/](data/) and run a `Main_*.m` script directly. Each `Main_*` script is **self-contained**: it redefines the physical constants, unit conversions (`nm`, `eV`, …), geometry parameters (`Ubarr`, `twell`, `tbarr`, `nbarr`), and solver parameters (`NzPlot`, `NzSC`, `Vvet`, `Evet`) at the top, then calls into the shared `f_*` functions. Use [data/Main_0_Template.m](data/Main_0_Template.m) as the starting point when adding a new study — it sets up the standard preamble and the `U`/`L` geometry-vector construction loop without running any solver.

There are no tests, no linting, no package manager. The MATLAB engine has no external dependencies beyond core MATLAB.

## Engine architecture

The simulation engine works in **SI units** throughout. A heterostructure is described by two parallel vectors of equal length:

- `U` — potential energy of each layer (J)
- `L` — thickness of each layer (m)

By convention the geometry is built as alternating `well, barrier` pairs and terminated with a final well (or barrier — see [data/Main_8_BoundState_Search.m](data/Main_8_BoundState_Search.m) which inverts the pattern for bound-state search). The first and last entries represent the semi-infinite contacts.

The solver is a **transmission-line analogue** built on local reflection coefficients (not transfer matrices):

1. [f_EvalGamma.m](data/f_EvalGamma.m) — given energy `E` and geometry `(U, L)`, computes `kz` (wavevector per layer), `Zinf` (characteristic impedance), and `GammaRight` (right-looking reflection coefficient) via a backward recursion from the right contact. Returns `GammaLeft` analogously when called from the left.
2. [f_EvalTau.m](data/f_EvalTau.m) — propagates the transmission amplitude `tau` across layers using the same `GammaRight` recursion; `T = |tau|^2`.
3. [f_EvalLDOS.m](data/f_EvalLDOS.m) — evaluates the local density of states at `NzPlot` points by combining `GammaLeft` and `GammaRight` into a local impedance, then taking `1i*(V - conj(V))`. Requires that both `GammaRight` and `GammaLeft` have been computed for the same energy.
4. [f_EvalLoopGain.m](data/f_EvalLoopGain.m) — used by bound-state search (Main_8): zeros of `1 - loopGain` at energies below `Ubarr` give bound-state energies.
5. [f_ApplyField.m](data/f_ApplyField.m) — applies a constant electric field by linearly tilting the potential between contacts and re-discretizing into a staircase of `NzSC` layers. The output `(USC, LSC)` is in the same `(U, L)` format and feeds directly back into `f_EvalGamma`. This is the mechanism for the bias/field studies (Main_9–11).
6. [f_Geom2Plot.m](data/f_Geom2Plot.m) / [f_evalVelocity.m](data/f_evalVelocity.m) — geometry-to-plot conversion and drift-velocity post-processing.

When sweeping energy or voltage, the pattern in every `Main_*` script is the same: an outer `for indE = 1:length(Evet)` (and optionally an inner `for V = Vvet`) calls `f_EvalGamma` → `f_EvalTau` / `f_EvalLDOS` per energy.

## Report structure

[qt_report.tex](qt_report.tex) is a Beamer presentation that `\input`s exactly five section files in [sections/](sections/), each corresponding to one physical topic and pulling figures produced by the matching `Main_*` script:

| Section | Figures from |
|---|---|
| `tran_coeff.tex` | Main_3 |
| `LDOS.tex` | Main_4, Main_5 |
| `transmittance_spectrum.tex` | Main_5 (1, 2, 20 barriers) |
| `bound_states.tex` | Main_8 |
| `electric_field.tex` | Main_9, Main_10 |

The section frames currently contain `\includegraphics` with placeholder `ADD CAPTION` captions and no surrounding prose — expect to write that content rather than just edit it.
