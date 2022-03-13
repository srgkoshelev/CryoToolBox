"""Helper geometric functions."""
from math import pi


def circle_area(D):
    """Calculate area of a circle."""
    return pi * D**2 / 4


def cylinder_volume(D, L):
    """Calculate cylinder volume"""
    return circle_area(D) * L
