#!/usr/bin/env python3
"""Original-NEC4 M9 ground-field routing and diagnostic policy.

This module selects the documented linear-table, model-fit, or asymptotic
region.  It deliberately does not evaluate a field: routing must remain
auditable independently of numerical implementations and fitted constants.
"""

from dataclasses import dataclass
import math

from nec4_m9_tables import near_table_contains, separation, table_scales


LINEAR_TABLE = "linear-table"
MODEL_FIT = "model-fit"
ASYMPTOTIC = "asymptotic"


@dataclass(frozen=True)
class Dispatch:
    method: str
    model_region: int | None
    grazing_risk: bool
    grazing_angle_degrees: float
    separation_wavelengths: float


def model_subregion(rho, observer_z, wavelength0, lg):
    """Return the lowest-numbered Figure 5-4 subregion on boundaries.

    The lower-left boundary follows z=min(rho,Lg) to rho=0.6 lambda0.
    Above/right of that corner, z=rho separates regions 1 and 3.  This is a
    direct transcription of the polygonal boundary in Figure 5-4; it is kept
    here rather than inferred from the sample lattices.
    """
    radial_corner = 0.6 * wavelength0
    if rho <= radial_corner:
        return 1 if observer_z >= min(rho, lg) else 2
    return 1 if observer_z >= rho else 3


def model_contains(rho, observer_z, source_z, wavelength0, lg):
    """Documented outer model extent following Figure 5-4 and its text."""
    return (rho >= 0 and observer_z >= 0 and source_z <= 0
            and rho <= 3.0 * wavelength0
            and observer_z <= 2.0 * wavelength0
            and abs(source_z) <= 0.8 * lg)


def grazing_risk(rho, observer_z, source_z, wavelength0, method):
    """Conservative diagnostic envelope for the original-GN2 limitation.

    This is a warning classifier, not a correction or claim of an exact NEC4
    implementation threshold.  It is restricted to the asymptotic path and
    requires the direct source-observer ray to lie within fifteen degrees of
    the interface.  The absolute endpoint heights are not the controlling
    variable: the published limitation is a grazing-angle effect.  M9
    qualification data freezes this conservative boundary against the M8
    direct oracle (severe through 11.31 degrees, materially improved at 21.80
    degrees for the ordinary second-order saddle term).
    """
    if method != ASYMPTOTIC or wavelength0 <= 0 or source_z > 0:
        return False, 90.0
    depth = abs(source_z)
    angle = math.degrees(math.atan2(observer_z + depth, rho)) if rho else 90.0
    risk = angle <= 15.0
    return risk, angle


def dispatch(rho, observer_z, source_z, k_lower, k_upper, wavelength0):
    if wavelength0 <= 0 or rho < 0 or observer_z < 0 or source_z > 0:
        raise ValueError("cross-interface dispatch requires rho,z>=0, z'<=0")
    scales = table_scales(k_lower, k_upper, wavelength0)
    if near_table_contains(rho, observer_z, source_z, scales):
        method, region = LINEAR_TABLE, None
    elif model_contains(rho, observer_z, source_z, wavelength0, scales.lg):
        method = MODEL_FIT
        region = model_subregion(rho, observer_z, wavelength0, scales.lg)
    else:
        method, region = ASYMPTOTIC, None
    risk, angle = grazing_risk(rho, observer_z, source_z,
                               wavelength0, method)
    return Dispatch(method, region, risk, angle,
                    separation(rho, observer_z, source_z) / wavelength0)
