# SimpleHarmonicOs

A small Python project that simulates a pendulum using a 4th-order Runge-Kutta (RK4) integrator.

This project was inspired by Walter Lewin's lectures:
https://www.youtube.com/watch?v=77ZF50ve6rs&t=24s

## What It Does

- Simulates angle and angular velocity over time
- Supports damping (`b`) so the pendulum can realistically settle
- Offers optional "stop when at rest" detection
- Plots:
  - Angle vs time
  - Angular velocity vs time
  - Mechanical energy vs time

## Physics Model

The equation solved is:

`theta'' + b * theta' + (g / l) * sin(theta) = 0`

where:
- `g` is gravity
- `l` is pendulum length
- `b` is linear damping coefficient

State vector:
- `theta` = angular displacement (rad)
- `omega` = angular velocity (rad/s)

## Requirements

- Python 3.10+
- `numpy`
- `matplotlib`
- Optional: Tkinter (only needed for `--gui`)

Install dependencies:

```bash
pip install numpy matplotlib
```

## Usage

Basic run:

```bash
python main.py
```

Run with custom parameters:

```bash
python main.py --theta-deg 30 --damping 0.08 --duration 25 --dt 0.005
```

Stop automatically when motion is effectively at rest:

```bash
python main.py --stop-when-rest --theta-deg 25 --damping 0.12
```

Use GUI prompt for initial angle:

```bash
python main.py --gui
```

Run without plotting (useful for checks/automation):

```bash
python main.py --no-plot
```

## CLI Parameters

- `--theta-deg`: initial angle in degrees (default: `20`)
- `--omega0`: initial angular velocity in rad/s (default: `0`)
- `--gravity`: gravitational acceleration in m/s^2 (default: `9.81`)
- `--length`: pendulum length in meters (default: `1.0`)
- `--mass`: pendulum bob mass in kg (default: `1.0`)
- `--damping`: linear damping coefficient in 1/s (default: `0.05`)
- `--dt`: integration step in seconds (default: `0.01`)
- `--duration`: max simulation time in seconds (default: `20`)
- `--stop-when-rest`: stop early if rest criteria are satisfied
- `--theta-tol-deg`: angular rest threshold in degrees (default: `0.05`)
- `--omega-tol`: angular velocity rest threshold in rad/s (default: `1e-3`)
- `--rest-window`: threshold hold time before stopping in seconds (default: `0.5`)
- `--gui`: asks initial angle via Tkinter dialog
- `--no-plot`: disables plotting

## Notes

- If `--damping` is `0`, total mechanical energy should remain nearly constant (numerical error aside).
- Smaller `--dt` generally improves accuracy but increases runtime.
