"""Simple, steady-state temperature model for a two-layer wall."""

from dataclasses import dataclass
from math import log


TOTAL_THICKNESS_MM = 450
BRICK_THICKNESS_MM = 120
BLOCK_THICKNESS_MM = TOTAL_THICKNESS_MM - BRICK_THICKNESS_MM
BLOCK_CONDUCTIVITY = 0.15  # W/(m·K), illustrative value
BRICK_CONDUCTIVITY = 0.70  # W/(m·K), illustrative value
INDOOR_SURFACE_RESISTANCE = 0.13  # m²·K/W
OUTDOOR_SURFACE_RESISTANCE = 0.04  # m²·K/W


@dataclass(frozen=True)
class WallResult:
    dew_point_c: float
    positions_mm: tuple[float, float, float]
    temperatures_c: tuple[float, float, float]
    isotherm_mm: float | None


def dew_point(temperature_c: float, relative_humidity_percent: float) -> float:
    """Return the Magnus approximation of the room-air dew point in °C."""
    if not 0 < relative_humidity_percent <= 100:
        raise ValueError("Относительная влажность должна быть больше 0 и не выше 100 %")
    a, b = 17.502, 240.97
    gamma = log(relative_humidity_percent / 100) + a * temperature_c / (b + temperature_c)
    return b * gamma / (a - gamma)


def calculate_wall(
    indoor_c: float,
    outdoor_c: float,
    relative_humidity_percent: float,
    block_conductivity: float = BLOCK_CONDUCTIVITY,
    brick_conductivity: float = BRICK_CONDUCTIVITY,
) -> WallResult:
    """Calculate temperatures from room side to outdoors and the dew-point isotherm."""
    if block_conductivity <= 0 or brick_conductivity <= 0:
        raise ValueError("Теплопроводность должна быть больше нуля")

    block_resistance = BLOCK_THICKNESS_MM / 1000 / block_conductivity
    brick_resistance = BRICK_THICKNESS_MM / 1000 / brick_conductivity
    total_resistance = (
        INDOOR_SURFACE_RESISTANCE
        + block_resistance
        + brick_resistance
        + OUTDOOR_SURFACE_RESISTANCE
    )
    heat_flux = (indoor_c - outdoor_c) / total_resistance
    room_surface = indoor_c - heat_flux * INDOOR_SURFACE_RESISTANCE
    interface = room_surface - heat_flux * block_resistance
    outside_surface = interface - heat_flux * brick_resistance
    temperatures = (room_surface, interface, outside_surface)
    positions = (0.0, float(BLOCK_THICKNESS_MM), float(TOTAL_THICKNESS_MM))
    dew = dew_point(indoor_c, relative_humidity_percent)

    crossing = None
    for start, end in ((0, 1), (1, 2)):
        t0, t1 = temperatures[start], temperatures[end]
        if min(t0, t1) <= dew <= max(t0, t1):
            if t0 == t1:
                crossing = positions[start]
            else:
                crossing = positions[start] + (positions[end] - positions[start]) * (dew - t0) / (t1 - t0)
            break

    return WallResult(dew, positions, temperatures, crossing)
