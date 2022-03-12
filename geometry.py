"""Helper geometric functions."""
from math import pi


def area_circle(D):
    """Calculate area of a circle."""
    return pi * D**2 / 4


def cylinder_volume(D, L):
    """Calculate cylinder volume"""
    return area_circle(D) * L
