"""Main application window and the wall diagram."""

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import (
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QMainWindow,
    QVBoxLayout,
    QWidget,
)

from dewpoint.model import (
    BLOCK_CONDUCTIVITY,
    BLOCK_THICKNESS_MM,
    BRICK_CONDUCTIVITY,
    BRICK_THICKNESS_MM,
    TOTAL_THICKNESS_MM,
    WallResult,
    calculate_wall,
)


class DiagramCanvas(FigureCanvasQTAgg):
    """Draw the wall section and steady-state temperature profile."""

    def __init__(self, parent: QWidget | None = None) -> None:
        self.figure = Figure(figsize=(9, 6), layout="constrained")
        super().__init__(self.figure)
        self.setParent(parent)

    def draw_result(self, result: WallResult, indoor_c: float, outdoor_c: float) -> None:
        self.figure.clear()
        wall, profile = self.figure.subplots(2, 1, sharex=True, height_ratios=(1, 2))

        for axes in (wall, profile):
            axes.axvspan(0, BLOCK_THICKNESS_MM, color="#bed9ed", alpha=0.8)
            axes.axvspan(BLOCK_THICKNESS_MM, TOTAL_THICKNESS_MM, color="#e8c8a8", alpha=0.8)
            axes.axvline(BLOCK_THICKNESS_MM, color="#555555", linewidth=1)
            if result.isotherm_mm is not None:
                axes.axvline(result.isotherm_mm, color="#bb4237", linestyle="--", linewidth=2)

        wall.set_ylim(0, 1)
        wall.set_yticks([])
        wall.text(BLOCK_THICKNESS_MM / 2, 0.5, f"Полистиролбетон\n{BLOCK_THICKNESS_MM} мм", ha="center", va="center")
        wall.text(BLOCK_THICKNESS_MM + BRICK_THICKNESS_MM / 2, 0.5, f"Силикатный\nкирпич\n{BRICK_THICKNESS_MM} мм", ha="center", va="center", fontsize=9)
        wall.set_title("Помещение  →  разрез стены  →  улица")

        profile.plot(result.positions_mm, result.temperatures_c, color="#175a82", marker="o", linewidth=2, label="Температура стены")
        profile.axhline(result.dew_point_c, color="#bb4237", linestyle=":", label=f"Точка росы воздуха: {result.dew_point_c:.1f} °C")
        profile.set_xlim(0, TOTAL_THICKNESS_MM)
        profile.set_xticks((0, BLOCK_THICKNESS_MM, TOTAL_THICKNESS_MM))
        profile.set_xlabel("Расстояние от внутренней поверхности, мм")
        profile.set_ylabel("Температура, °C")
        profile.grid(True, alpha=0.35)
        profile.legend(loc="best")
        profile.set_title(f"Воздух: внутри {indoor_c:g} °C, снаружи {outdoor_c:g} °C")
        self.draw_idle()


class MainWindow(QMainWindow):
    """Inputs and diagram for the fixed 450 mm wall."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Точка росы — стена дома")
        self.resize(1100, 800)
        self.setMinimumSize(720, 600)

        central = QWidget(self)
        layout = QVBoxLayout(central)
        self.setCentralWidget(central)

        inputs = QGroupBox("Исходные данные")
        form = QFormLayout(inputs)
        self.indoor = self._spin(-60, 80, 25, " °C", 1)
        self.outdoor = self._spin(-60, 80, -20, " °C", 1)
        self.humidity = self._spin(1, 100, 50, " %", 0)
        self.block_lambda = self._spin(0.01, 5, BLOCK_CONDUCTIVITY, " Вт/(м·К)", 3)
        self.brick_lambda = self._spin(0.01, 5, BRICK_CONDUCTIVITY, " Вт/(м·К)", 3)
        form.addRow("Температура внутри:", self.indoor)
        form.addRow("Температура снаружи:", self.outdoor)
        form.addRow("Влажность внутри:", self.humidity)
        form.addRow("λ полистиролбетона:", self.block_lambda)
        form.addRow("λ силикатного кирпича:", self.brick_lambda)
        form.addRow("Стена:", QLabel(f"{BLOCK_THICKNESS_MM} мм блок + {BRICK_THICKNESS_MM} мм кирпич = {TOTAL_THICKNESS_MM} мм"))
        layout.addWidget(inputs)

        self.diagram = DiagramCanvas(central)
        layout.addWidget(self.diagram, stretch=1)
        self.result_label = QLabel()
        layout.addWidget(self.result_label)
        layout.addWidget(QLabel("Условная изотерма: для оценки конденсации внутри стены нужен расчёт переноса водяного пара."))

        for spin in (self.indoor, self.outdoor, self.humidity, self.block_lambda, self.brick_lambda):
            spin.valueChanged.connect(self.update_diagram)
        self.update_diagram()

    @staticmethod
    def _spin(minimum: float, maximum: float, value: float, suffix: str, decimals: int) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setRange(minimum, maximum)
        spin.setDecimals(decimals)
        spin.setSingleStep(1 if decimals <= 1 else 0.01)
        spin.setSuffix(suffix)
        spin.setValue(value)
        return spin

    def update_diagram(self) -> None:
        result = calculate_wall(
            self.indoor.value(), self.outdoor.value(), self.humidity.value(),
            self.block_lambda.value(), self.brick_lambda.value(),
        )
        self.diagram.draw_result(result, self.indoor.value(), self.outdoor.value())
        if result.isotherm_mm is None:
            location = "В пределах стены изотерма точки росы отсутствует."
        else:
            material = "полистиролбетон" if result.isotherm_mm <= BLOCK_THICKNESS_MM else "силикатный кирпич"
            location = f"Условная изотерма: {result.isotherm_mm:.0f} мм от внутренней поверхности ({material})."
        if result.temperatures_c[0] <= result.dew_point_c:
            location += " Внутренняя поверхность не теплее точки росы: возможен поверхностный конденсат."
        self.result_label.setText(f"Точка росы комнатного воздуха: {result.dew_point_c:.1f} °C. {location}")
