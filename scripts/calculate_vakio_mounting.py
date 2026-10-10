"""Reproduce the preliminary VAKIO installation calculations and section drawing.

Run from the project root: python -m scripts.calculate_vakio_mounting
The thermal bridge is assessed approximately; the casing and moisture transport
are not simulated.
"""

import json
from math import cos, exp, log, pi, radians, tan
from pathlib import Path

from dewpoint.model import calculate_wall
from scripts.axisymmetric_vakio import solve, summarize


def saturation_pressure(temperature_c: float) -> float:
    """Magnus approximation over liquid water, in Pa."""
    return 613.65 * exp(17.502 * temperature_c / (240.97 + temperature_c))


def temperature_for_pressure(pressure_pa: float) -> float:
    gamma = log(pressure_pa / 613.65)
    return 240.97 * gamma / (17.502 - gamma)


def calculate() -> dict:
    indoor, outdoor, wall_mm = 25.0, -20.0, 450.0
    humidity_table = []
    for humidity in (40, 50, 60, 70):
        vapor_pressure = humidity / 100 * saturation_pressure(indoor)
        dew = temperature_for_pressure(vapor_pressure)
        surface_80 = temperature_for_pressure(vapor_pressure / 0.8)
        humidity_table.append({
            "indoor_rh_percent": humidity,
            "dew_point_c": dew,
            "surface_temperature_at_80_percent_rh_c": surface_80,
            "minimum_temperature_factor_for_80_percent": (surface_80 - outdoor) / (indoor - outdoor),
        })

    # Pipe dimensions are preliminary; measure the purchased pipe before drilling.
    bore_mm, pipe_wall_mm, outer_insulation_mm, seal_gap_mm = 125, 1, 5, 2.5
    pipe_od_mm = bore_mm + 2 * pipe_wall_mm
    insulated_od_mm = pipe_od_mm + 2 * outer_insulation_mm
    hole_mm = insulated_od_mm + 2 * seal_gap_mm
    slope_table = [{
        "angle_degrees": angle,
        "outward_drop_mm": wall_mm * tan(radians(angle)),
        "axial_length_mm": wall_mm / cos(radians(angle)),
    } for angle in (3, 4, 5)]

    room_pressure = 0.5 * saturation_pressure(indoor)
    insulation_lambda = 0.036
    planar_screening = []
    for surface_resistance in (0.13, 0.25):
        for thickness_mm in (5, 9, 19, 25, 32, 60):
            resistance = thickness_mm / 1000 / insulation_lambda
            temperature = indoor - (indoor - outdoor) * surface_resistance / (surface_resistance + resistance)
            planar_screening.append({
                "room_surface_resistance_m2k_w": surface_resistance,
                "insulation_mm": thickness_mm,
                "room_facing_surface_c": temperature,
                "calculated_surface_rh_percent": 100 * room_pressure / saturation_pressure(temperature),
            })

    screening_target_c = 19
    required_thickness = [{
        "room_surface_resistance_m2k_w": rsi,
        "insulation_mm_for_19_c_surface": 1000 * insulation_lambda * rsi * ((indoor - outdoor) / (indoor - screening_target_c) - 1),
    } for rsi in (0.13, 0.25)]

    wall = calculate_wall(indoor, outdoor, 50)
    collar_lambda, collar_thickness_m, sheltered_rsi = 0.036, 0.100, 0.25
    collar_surface = indoor - (indoor - outdoor) * sheltered_rsi / (sheltered_rsi + collar_thickness_m / collar_lambda)
    volume_flow_m3h = 60
    heating_load_w = 1.2 * 1005 * volume_flow_m3h / 3600 * (indoor - outdoor)
    return {
        "status": "preliminary; not a validated thermal-bridge or condensation simulation",
        "conditions": {"indoor_c": indoor, "outdoor_c": outdoor, "indoor_rh_percent": 50},
        "wall": {
            "total_mm": wall_mm,
            "block_mm_assumed": 330,
            "brick_mm_assumed": 120,
            "air_gap_mm_confirmed": 0,
            "internal_finish_confirmed": "plaster; thickness unknown",
            "plaster_thickness_mm_estimated": [5, 10],
            "block_wall_load_bearing_confirmed": False,
            "block_lambda_assumed_w_mk": 0.15,
            "brick_lambda_assumed_w_mk": 0.70,
            "undisturbed_wall_room_surface_c": wall.temperatures_c[0],
        },
        "humidity_thresholds": humidity_table,
        "candidate_dimensions_mm": {
            "pipe_internal": bore_mm,
            "pipe_wall": pipe_wall_mm,
            "pipe_external": pipe_od_mm,
            "outer_insulation": outer_insulation_mm,
            "insulated_external": insulated_od_mm,
            "radial_seal_gap": seal_gap_mm,
            "calculated_hole": hole_mm,
            "factory_internal_insulation_nominal_thickness": 9,
            "approximate_free_air_bore": bore_mm - 2 * 9,
        },
        "slopes": slope_table,
        "candidate_room_collar": {
            "material": "unfaced closed-cell cellular glass; grade to be selected",
            "status": "planar screening complemented by an approximate axisymmetric model; device and fixings excluded",
            "insulation_thickness_mm": 100,
            "plan_size_mm": [400, 400],
            "lambda_assumed_w_mk": collar_lambda,
            "room_surface_resistance_m2k_w": sheltered_rsi,
            "cold_back_face_assumed_c": outdoor,
            "screening_room_surface_c": collar_surface,
            "screening_surface_rh_percent": 100 * room_pressure / saturation_pressure(collar_surface),
        },
        "screening_planar_cold_back_face": {
            "warning": "Thicknesses below are not specifications for insulation around the sleeve or behind the device.",
            "cold_back_face_c": outdoor,
            "insulation_lambda_assumed_w_mk": insulation_lambda,
            "target_room_facing_surface_c": screening_target_c,
            "required_thickness": required_thickness,
            "cases": planar_screening,
        },
        "airflow": {
            "nominal_max_m3h": volume_flow_m3h,
            "assumed_reference_density_kg_m3": 1.2,
            "assumed_cp_j_kg_k": 1005,
            "room_heating_load_w_approx": heating_load_w,
            "velocity_in_125_mm_bore_m_s": volume_flow_m3h / 3600 / (pi * 0.125**2 / 4),
            "velocity_in_107_mm_bore_m_s_approx": volume_flow_m3h / 3600 / (pi * 0.107**2 / 4),
        },
    }


def save_figure(figure, path: Path) -> None:
    figure.savefig(path)
    if path.suffix.lower() == ".svg":
        content = path.read_text()
        path.write_text("\n".join(line.rstrip() for line in content.splitlines()) + "\n")


def draw_section(path: Path) -> None:
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle, Rectangle

    fig = Figure(figsize=(14, 6), layout="constrained")
    side, section = fig.subplots(1, 2, gridspec_kw={"width_ratios": (2.2, 1)})
    side.add_patch(Rectangle((0, -240), 330, 480, facecolor="#bed9ed"))
    side.add_patch(Rectangle((330, -240), 120, 480, facecolor="#e8c8a8"))
    side.add_patch(Rectangle((0, -200), 100, 400, facecolor="#c8b8d5"))
    xs = (0, 450)
    centers = (0, -450 * tan(radians(4)))
    # Band thickness in this elevation is adjusted for the 4-degree slope.
    for diameter, color in ((142, "#cccccc"), (137, "#b0c9aa"), (127, "#555555"), (125, "#eddca4"), (107, "#fafafa")):
        half_height = diameter / 2 / cos(radians(4))
        side.fill_between(xs, [y - half_height for y in centers], [y + half_height for y in centers], color=color)
    side.plot(xs, centers, linestyle="--", color="#175a82", linewidth=1)
    for x, color in ((0, "#bb4237"), (450, "#306c47")):
        y = -x * tan(radians(4))
        side.plot((x, x), (y - 73, y - 55), color=color, linewidth=4)
        side.plot((x, x), (y + 55, y + 73), color=color, linewidth=4)
    side.text(210, 210, "Полистиролбетон\n330 мм — принято", ha="center", fontsize=9)
    side.text(390, 210, "Силикатный\nкирпич 120 мм\n— принято", ha="center", fontsize=8)
    side.text(50, -150, "Вставка\n100 мм\nв зоне\n400×400 мм", ha="center", fontsize=9)
    side.text(-8, -140, "Помещение\n+25 °C", ha="right")
    side.text(458, -140, "Улица\n−20 °C", ha="left")
    side.annotate("Уклон 4° наружу\nперепад 31,5 мм", xy=(310, -22), xytext=(180, -170), arrowprops={"arrowstyle": "->"})
    side.annotate("Герметичное внутреннее\nпримыкание к вставке", xy=(0, 66), xytext=(5, 280), arrowprops={"arrowstyle": "->"}, fontsize=9)
    side.annotate("Наружное погодостойкое\nпримыкание; слив в канале открыт", xy=(450, 35), xytext=(270, 280), arrowprops={"arrowstyle": "->"}, fontsize=9)
    side.set_xlim(-60, 515)
    side.set_ylim(-265, 310)
    side.set_aspect("equal")
    side.set_xlabel("Толщина стены, мм")
    side.set_ylabel("Высота относительно внутреннего центра, мм")
    side.set_title("Продольный разрез прохода (корпус и решётка не показаны)")

    for radius, color, label in (
        (71, "#cccccc", "Монтажный шов: 2,5 мм"),
        (68.5, "#b0c9aa", "Наружная изоляция: 5 мм"),
        (63.5, "#555555", "Пластиковая гильза: 1 мм"),
        (62.5, "#eddca4", "Штатная внутренняя изоляция: ≈9 мм"),
        (53.5, "#fafafa", "Воздушный просвет: ≈107 мм"),
    ):
        section.add_patch(Circle((0, 0), radius, facecolor=color, edgecolor="#555555", linewidth=0.5, label=label))
    section.set_aspect("equal")
    section.set_xlim(-85, 85)
    section.set_ylim(-140, 85)
    section.axis("off")
    section.set_title("Поперечный разрез\nрасчётное отверстие Ø142 мм")
    section.legend(loc="lower center", fontsize=8)
    fig.suptitle("Предварительный узел: Ø142 и вставка из пеностекла; требуется проверка примыканий и крепежа", fontsize=12)
    save_figure(fig, path)


def draw_temperature_field(path: Path, result, comparisons: dict) -> None:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(13, 6), layout="constrained")
    field, surface = fig.subplots(1, 2)
    color = field.pcolormesh(result.depths_m * 1000, result.radii_m * 1000, result.temperatures_c, shading="auto", vmin=-20, vmax=25, cmap="coolwarm", rasterized=True)
    field.set_ylim(53.5, 300)
    field.set_xlabel("Глубина от комнатной поверхности, мм")
    field.set_ylabel("Радиус от оси канала, мм")
    field.set_title("Вставка 100 мм, λ = 0,036 Вт/(м·К)")
    field.axvline(100, color="black", linestyle=":", linewidth=1)
    field.axvline(330, color="black", linestyle=":", linewidth=1)
    fig.colorbar(color, ax=field, label="Температура, °C")
    for label, sample in comparisons.items():
        surface.plot(sample.radii_m * 1000, sample.room_surface_c, label=label)
    surface.set_xlim(71, 300)
    surface.set_ylim(10, 25)
    surface.axhline(19, color="#175a82", linestyle=":", label="Проектная цель +19 °C")
    threshold = temperature_for_pressure(0.5 * saturation_pressure(25) / 0.8)
    surface.axhline(threshold, color="#bb4237", linestyle=":", label="80% RH поверхности при 50% RH комнаты")
    surface.axvspan(71, 89.5, color="#cccccc", alpha=0.4)
    surface.set_xlabel("Радиус от оси канала, мм")
    surface.set_ylabel("Температура комнатной поверхности, °C")
    surface.set_title("Серый участок требует герметичного сопряжения")
    surface.legend(fontsize=8)
    surface.grid(alpha=0.3)
    fig.suptitle("Осесимметричная оценка +25/−20 °C: корпус, прокладка, отделка и крепёж не смоделированы", fontsize=11)
    save_figure(fig, path)


if __name__ == "__main__":
    docs = Path(__file__).resolve().parents[1] / "docs"
    docs.mkdir(exist_ok=True)
    report = docs / "vakio_openair_calculation.json"
    result = calculate()
    samples = {
        "Без вставки": solve(None),
        "Вставка: λ = 0,050": solve(0.050),
        "Вставка: λ = 0,036": solve(0.036, step_m=0.0025),
    }
    result["axisymmetric_screening"] = {
        "assumptions": "450 mm wall; 330/120 mm layers; circular insert radius 200 mm, depth 100 mm; channel face fixed at -20 C; Rsi .25, Rse .04; far radius .50 m adiabatic; materials as in solver; no casing, gasket, finish, glue or anchors; 4-degree slope omitted",
        "cases": {name: summarize(sample) for name, sample in samples.items()},
    }
    report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    drawing = docs / "vakio_openair_section.svg"
    draw_section(drawing)
    field = docs / "vakio_openair_thermal_field.svg"
    draw_temperature_field(field, samples["Вставка: λ = 0,036"], samples)
    print(report)
    print(drawing)
    print(field)
