from __future__ import annotations

import argparse
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

try:
    import tkinter as tk
    from tkinter import simpledialog
except Exception:
    tk = None
    simpledialog = None


@dataclass
class PendulumParams:
    gravity: float = 9.81
    length: float = 1.0
    mass: float = 1.0
    damping: float = 0.05
    dt: float = 0.01
    duration: float = 20.0
    stop_when_rest: bool = False
    theta_tol: float = np.deg2rad(0.05)
    omega_tol: float = 1e-3
    rest_window: float = 0.5


def derivatives(state: np.ndarray, _t: float, params: PendulumParams) -> np.ndarray:
    theta, omega = state
    dtheta_dt = omega
    domega_dt = -(params.gravity / params.length) * np.sin(theta) - params.damping * omega
    return np.array([dtheta_dt, domega_dt], dtype=float)


def rk4_step(state: np.ndarray, t: float, dt: float, params: PendulumParams) -> np.ndarray:
    k1 = dt * derivatives(state, t, params)
    k2 = dt * derivatives(state + 0.5 * k1, t + 0.5 * dt, params)
    k3 = dt * derivatives(state + 0.5 * k2, t + 0.5 * dt, params)
    k4 = dt * derivatives(state + k3, t + dt, params)
    return state + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def mechanical_energy(theta: float, omega: float, params: PendulumParams) -> float:
    kinetic = 0.5 * params.mass * (params.length**2) * (omega**2)
    potential = params.mass * params.gravity * params.length * (1.0 - np.cos(theta))
    return kinetic + potential


def simulate(
    theta0: float, omega0: float, params: PendulumParams
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    max_steps = int(np.ceil(params.duration / params.dt)) + 1
    time = np.zeros(max_steps, dtype=float)
    theta = np.zeros(max_steps, dtype=float)
    omega = np.zeros(max_steps, dtype=float)
    energy = np.zeros(max_steps, dtype=float)
    state = np.array([theta0, omega0], dtype=float)

    rest_steps_required = max(1, int(np.ceil(params.rest_window / params.dt)))
    rest_counter = 0
    used_steps = max_steps

    for i in range(max_steps):
        current_time = i * params.dt
        time[i] = current_time
        theta[i], omega[i] = state
        energy[i] = mechanical_energy(theta[i], omega[i], params)

        is_rest_state = abs(theta[i]) < params.theta_tol and abs(omega[i]) < params.omega_tol
        if params.stop_when_rest and current_time > 0.0 and is_rest_state:
            rest_counter += 1
            if rest_counter >= rest_steps_required:
                used_steps = i + 1
                break
        else:
            rest_counter = 0

        if i < max_steps - 1:
            state = rk4_step(state, current_time, params.dt, params)

    return time[:used_steps], theta[:used_steps], omega[:used_steps], energy[:used_steps]


def plot_results(
    time: np.ndarray, theta: np.ndarray, omega: np.ndarray, energy: np.ndarray, params: PendulumParams
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(10, 9))

    axes[0].plot(time, theta, label="Angle (rad)")
    axes[0].set_title("Simple Harmonic Oscillator (Pendulum)")
    axes[0].set_ylabel("Theta (rad)")
    axes[0].legend(loc="upper right")
    axes[0].grid(alpha=0.3)

    axes[1].plot(time, omega, color="tab:orange", label="Angular velocity (rad/s)")
    axes[1].set_ylabel("Omega (rad/s)")
    axes[1].legend(loc="upper right")
    axes[1].grid(alpha=0.3)

    axes[2].plot(time, energy, color="tab:green", label="Mechanical energy (J)")
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Energy (J)")
    axes[2].legend(loc="upper right")
    axes[2].grid(alpha=0.3)

    fig.suptitle(
        f"g={params.gravity} m/s^2, l={params.length} m, damping={params.damping:.3f} 1/s, dt={params.dt}s",
        fontsize=10,
    )
    fig.tight_layout()
    plt.show()


def get_theta_from_gui(default_theta_deg: float) -> float:
    if tk is None or simpledialog is None:
        raise RuntimeError("Tkinter is not available in this Python environment.")

    root = tk.Tk()
    root.withdraw()
    try:
        value = simpledialog.askfloat(
            "Input",
            "Enter initial angle in degrees:",
            initialvalue=default_theta_deg,
        )
    finally:
        root.destroy()

    if value is None:
        raise ValueError("No input provided in GUI.")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Simulate a simple pendulum with RK4 integration.")
    parser.add_argument("--theta-deg", type=float, default=20.0, help="Initial angle in degrees.")
    parser.add_argument("--omega0", type=float, default=0.0, help="Initial angular velocity in rad/s.")
    parser.add_argument("--gravity", type=float, default=9.81, help="Gravity in m/s^2.")
    parser.add_argument("--length", type=float, default=1.0, help="Pendulum length in meters.")
    parser.add_argument("--mass", type=float, default=1.0, help="Pendulum bob mass in kilograms.")
    parser.add_argument("--damping", type=float, default=0.05, help="Linear damping coefficient in 1/s.")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step in seconds.")
    parser.add_argument("--duration", type=float, default=20.0, help="Maximum simulation duration in seconds.")
    parser.add_argument("--stop-when-rest", action="store_true", help="Stop once rest criteria are met.")
    parser.add_argument(
        "--theta-tol-deg",
        type=float,
        default=0.05,
        help="Rest threshold for |theta| in degrees (used with --stop-when-rest).",
    )
    parser.add_argument(
        "--omega-tol",
        type=float,
        default=1e-3,
        help="Rest threshold for |omega| in rad/s (used with --stop-when-rest).",
    )
    parser.add_argument(
        "--rest-window",
        type=float,
        default=0.5,
        help="How long thresholds must hold before stopping, in seconds.",
    )
    parser.add_argument("--gui", action="store_true", help="Prompt for initial angle with a Tkinter dialog.")
    parser.add_argument("--no-plot", action="store_true", help="Skip plotting (useful for automated runs).")
    return parser.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    if args.gravity <= 0:
        raise ValueError("gravity must be > 0")
    if args.length <= 0:
        raise ValueError("length must be > 0")
    if args.mass <= 0:
        raise ValueError("mass must be > 0")
    if args.dt <= 0:
        raise ValueError("dt must be > 0")
    if args.duration <= 0:
        raise ValueError("duration must be > 0")
    if args.damping < 0:
        raise ValueError("damping must be >= 0")
    if args.omega_tol <= 0:
        raise ValueError("omega-tol must be > 0")
    if args.theta_tol_deg <= 0:
        raise ValueError("theta-tol-deg must be > 0")
    if args.rest_window <= 0:
        raise ValueError("rest-window must be > 0")


def main() -> None:
    args = parse_args()
    validate_inputs(args)

    theta_deg = args.theta_deg
    if args.gui:
        theta_deg = get_theta_from_gui(theta_deg)

    params = PendulumParams(
        gravity=args.gravity,
        length=args.length,
        mass=args.mass,
        damping=args.damping,
        dt=args.dt,
        duration=args.duration,
        stop_when_rest=args.stop_when_rest,
        theta_tol=np.deg2rad(args.theta_tol_deg),
        omega_tol=args.omega_tol,
        rest_window=args.rest_window,
    )

    time, theta, omega, energy = simulate(np.deg2rad(theta_deg), args.omega0, params)
    final_time = float(time[-1])
    print(f"Simulated {final_time:.2f}s over {len(time)} steps.")
    print(f"Final state: theta={theta[-1]:.5f} rad, omega={omega[-1]:.5f} rad/s")

    if params.stop_when_rest:
        if final_time < params.duration - params.dt:
            print(f"Rest condition reached at ~{final_time:.2f}s.")
        else:
            print("Rest condition not reached within the selected duration.")

    if not args.no_plot:
        plot_results(time, theta, omega, energy, params)


if __name__ == "__main__":
    main()
