#!/usr/bin/env python3
"""Shared centerline construction for the coated-wire bend envelope study."""

import math


WIRE_RADIUS_M = 0.002
COATING_RADIUS_M = 0.004
HALF_SPAN_M = 0.050
BASE_HEIGHT_M = 0.187


def case_name(chords, radius_ratio):
    if chords == 0:
        return "sharp"
    return f"arc-n{chords}-r{radius_ratio:g}b".replace(".", "p")


def bend_case(chords, radius_ratio=0.0):
    """Return a constant-length Z article with two rounded 90-degree bends."""
    if chords == 0:
        points = [(-HALF_SPAN_M, 0.0, -BASE_HEIGHT_M),
                  (-HALF_SPAN_M, 0.0, 0.0),
                  (HALF_SPAN_M, 0.0, 0.0),
                  (HALF_SPAN_M, 0.0, BASE_HEIGHT_M)]
        return {
            "name": "sharp", "chords_per_corner": 0,
            "radius_ratio": 0.0, "radius_m": 0.0,
            "maximum_turn_deg": 90.0, "minimum_chord_ratio": 0.0,
            "points_m": points,
        }
    radius = radius_ratio*COATING_RADIUS_M
    if not 0 < radius < HALF_SPAN_M:
        raise ValueError("bend radius must fit inside the U half-span")
    chord = 2*radius*math.sin(math.pi/(4*chords))
    arc_length = chords*chord
    height = BASE_HEIGHT_M + 2*radius-arc_length
    points = [(-HALF_SPAN_M, 0.0, -height)]
    left_center = (-HALF_SPAN_M+radius, -radius)
    for index in range(chords+1):
        angle = math.pi - index*(math.pi/2)/chords
        points.append((left_center[0]+radius*math.cos(angle), 0.0,
                       left_center[1]+radius*math.sin(angle)))
    right_center = (HALF_SPAN_M-radius, radius)
    # The straight bottom run is implicit between the two arc endpoint lists.
    for index in range(chords+1):
        angle = -math.pi/2 + index*(math.pi/2)/chords
        points.append((right_center[0]+radius*math.cos(angle), 0.0,
                       right_center[1]+radius*math.sin(angle)))
    points.append((HALF_SPAN_M, 0.0, height))
    # Remove duplicated points introduced at the vertical/arc boundaries.
    clean = []
    for point in points:
        if not clean or math.dist(point, clean[-1]) > 1e-12:
            clean.append(point)
    maximum_turn = 45.0 if chords == 1 else 90.0/chords
    return {
        "name": case_name(chords, radius_ratio),
        "chords_per_corner": chords,
        "radius_ratio": radius_ratio,
        "radius_m": radius,
        "maximum_turn_deg": maximum_turn,
        "minimum_chord_ratio": chord/COATING_RADIUS_M,
        "points_m": clean,
    }


def insert_feed_point(points):
    """Insert the origin on the horizontal run for surface-MoM feeding."""
    result = []
    for left, right in zip(points, points[1:]):
        result.append(left)
        if (abs(left[2]) < 1e-12 and abs(right[2]) < 1e-12 and
                left[0] < 0 < right[0]):
            result.append((0.0, 0.0, 0.0))
    result.append(points[-1])
    return result


def study_cases():
    cases = [bend_case(0)]
    for chords in (1, 2, 3, 4, 6):
        for ratio in (1, 2, 4, 6, 8, 12):
            cases.append(bend_case(chords, ratio))
    return cases
