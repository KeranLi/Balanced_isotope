"""Compatibility entry point for the canonical siliciclastic Mg model.

The maintained implementation lives in :mod:`systems.mg.silicate`. This
module remains executable so older commands do not invoke the obsolete
Rayleigh equation that previously clipped retained Mg above 100%.
"""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from systems.mg.silicate import (
    SilicateMgSystem,
    SilicateModelResult,
    SilicateWeatheringModel,
    SilicateWeatheringParams,
)

__all__ = [
    "SilicateMgSystem",
    "SilicateModelResult",
    "SilicateWeatheringModel",
    "SilicateWeatheringParams",
]


def example_usage():
    model = SilicateMgSystem(basin="changjiang")
    print("Siliciclastic Mg isotope weathering model (Hu et al., 2023)")
    print("delta26Mg_clay  f_Mg   F_silicate (mol/yr)  status")
    for delta26 in (-0.15, -0.10, -0.05, 0.00):
        result = model.calculate_weathering_flux(delta26)
        print(
            f"{delta26:+.2f}           {result.f_Mg:.3f}  "
            f"{result.F_silicate:.3e}         {result.model_status}"
        )
    return model


if __name__ == "__main__":
    example_usage()
