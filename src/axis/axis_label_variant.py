from enum import Enum


class AxisLabelVariant(Enum):
    Radius = 0
    Mass = 1
    Bin = 2


class KernelErrorVariant(Enum):
    KgPerSecondPerDensity = 0
    KgPerCollision = 1
    PercentPerCollision = 2
