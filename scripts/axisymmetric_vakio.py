"""Approximate finite-volume heat conduction around an insulated wall sleeve.

Axisymmetry replaces the square insert by its inscribed circle. The cold air
boundary is fixed at outdoor temperature along the entire channel. The casing,
gasket, plaster finish, adhesive and anchors are not represented by this model.
"""

from dataclasses import dataclass
from math import pi

import numpy as np


@dataclass
class ThermalResult:
    radii_m: np.ndarray
    depths_m: np.ndarray
    temperatures_c: np.ndarray
    room_surface_c: np.ndarray
    iterations: int
    relative_residual: float
    relative_heat_imbalance: float

    def temperature_at_radius(self, radius_m: float) -> float:
        return float(np.interp(radius_m, self.radii_m, self.room_surface_c))


def solve(
    collar_lambda: float | None,
    step_m: float = 0.005,
    indoor_c: float = 25,
    outdoor_c: float = -20,
) -> ThermalResult:
    """Solve a 450 mm wall, Ø142 hole and optional 100 mm deep insert.

    Conductivities in W/(m·K): block .15, brick .70, inner insulation .04,
    pipe .17, external insulation .036, annular filling .04. These are input
    assumptions, not measured properties of the installed assembly.
    """
    radial_segments = (
        (0.0535, 0.0625, min(0.003, step_m)),
        (0.0625, 0.0635, min(0.001, step_m)),
        (0.0635, 0.0685, min(0.0025, step_m)),
        (0.0685, 0.071, min(0.00125, step_m)),
        (0.071, 0.10, step_m / 2),
        (0.10, 0.20, step_m),
        (0.20, 0.50, step_m * 2),
    )

    def edges(segments):
        parts = []
        for start, end, spacing in segments:
            count = int(np.ceil((end - start) / spacing))
            parts.extend(np.linspace(start, end, count + 1)[:-1])
        return np.array(parts + [segments[-1][1]])

    re = edges(radial_segments)
    ze = edges(((0, 0.10, step_m), (0.10, 0.33, step_m), (0.33, 0.45, step_m)))
    r, z = (re[:-1] + re[1:]) / 2, (ze[:-1] + ze[1:]) / 2
    dz = np.diff(ze)
    area = pi * (re[1:]**2 - re[:-1]**2)
    conductivity = np.where(z[None, :] < 0.33, 0.15, 0.70) * np.ones((len(r), 1))
    if collar_lambda is not None:
        conductivity = np.where((r[:, None] < 0.20) & (z[None, :] < 0.10), collar_lambda, conductivity)
    for start, end, value in (
        (0.0535, 0.0625, 0.04),
        (0.0625, 0.0635, 0.17),
        (0.0635, 0.0685, 0.036),
        (0.0685, 0.071, 0.04),
    ):
        conductivity = np.where((r[:, None] >= start) & (r[:, None] < end), value, conductivity)

    # Exact series resistance of two radial half-cells; annular axial areas.
    radial_g = 2 * pi * dz[None, :] / (
        np.log(re[1:-1] / r[:-1])[:, None] / conductivity[:-1]
        + np.log(r[1:] / re[1:-1])[:, None] / conductivity[1:]
    )
    axial_g = area[:, None] / (
        dz[:-1][None, :] / (2 * conductivity[:, :-1])
        + dz[1:][None, :] / (2 * conductivity[:, 1:])
    )
    cold_g = 2 * pi * dz / (np.log(r[0] / re[0]) / conductivity[0])
    room_g = area / (dz[0] / (2 * conductivity[:, 0]) + 0.25)
    outside_g = area / (dz[-1] / (2 * conductivity[:, -1]) + 0.04)
    diagonal = np.zeros_like(conductivity)
    diagonal[:-1] += radial_g
    diagonal[1:] += radial_g
    diagonal[:, :-1] += axial_g
    diagonal[:, 1:] += axial_g
    diagonal[0] += cold_g
    diagonal[:, 0] += room_g
    diagonal[:, -1] += outside_g
    rhs = np.zeros_like(conductivity)
    rhs[0] += cold_g * outdoor_c
    rhs[:, 0] += room_g * indoor_c
    rhs[:, -1] += outside_g * outdoor_c

    def multiply(temperature):
        value = diagonal * temperature
        value[:-1] -= radial_g * temperature[1:]
        value[1:] -= radial_g * temperature[:-1]
        value[:, :-1] -= axial_g * temperature[:, 1:]
        value[:, 1:] -= axial_g * temperature[:, :-1]
        return value

    # Jacobi-preconditioned conjugate gradients on a symmetric positive matrix.
    temperature = np.full_like(conductivity, (indoor_c + outdoor_c) / 2)
    residual = rhs - multiply(temperature)
    tolerance = 1e-9 * max(float(np.linalg.norm(rhs)), 1)
    direction = residual / diagonal
    product = float(np.sum(residual * direction))
    iterations = 0
    while float(np.linalg.norm(residual)) > tolerance:
        if iterations >= 10000:
            raise RuntimeError("Thermal solver did not converge")
        matrix_direction = multiply(direction)
        alpha = product / float(np.sum(direction * matrix_direction))
        temperature += alpha * direction
        residual -= alpha * matrix_direction
        preconditioned = residual / diagonal
        new_product = float(np.sum(residual * preconditioned))
        direction = preconditioned + (new_product / product) * direction
        product = new_product
        iterations += 1

    surface = indoor_c - (indoor_c - temperature[:, 0]) * 0.25 / (0.25 + dz[0] / (2 * conductivity[:, 0]))
    room_heat = float(np.sum(room_g * (indoor_c - temperature[:, 0])))
    outside_heat = float(np.sum(outside_g * (temperature[:, -1] - outdoor_c)))
    channel_heat = float(np.sum(cold_g * (temperature[0] - outdoor_c)))
    heat_scale = max(abs(room_heat), abs(outside_heat) + abs(channel_heat), 1)
    return ThermalResult(
        r, z, temperature, surface, iterations,
        float(np.linalg.norm(rhs - multiply(temperature))) / max(float(np.linalg.norm(rhs)), 1),
        abs(room_heat - outside_heat - channel_heat) / heat_scale,
    )


def summarize(result: ThermalResult) -> dict:
    return {
        "cells": list(result.temperatures_c.shape),
        "iterations": result.iterations,
        "relative_residual": result.relative_residual,
        "relative_heat_imbalance": result.relative_heat_imbalance,
        "room_surface_c_at_radius_mm": {
            str(radius): result.temperature_at_radius(radius / 1000)
            for radius in (71, 80, 89.5, 100, 150, 200, 300)
        },
    }
