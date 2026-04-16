"""Reusable geometry helpers for CryoToolBox.

These functions are intentionally small and unit-aware. They are most useful
when the same geometric relationship appears in several engineering
calculations and keeping the formula in one place improves readability and
consistency.
"""
from math import pi


def circle_area(D):
    """Calculate the area of a circle from its diameter.

    Parameters
    ----------
    D : ureg.Quantity {length: 1}
        The diameter of the circle.

    Returns
    -------
    ureg.Quantity {length: 2}
        The area of the circle.

    Notes
    -----
    The formula for calculating the area of a circle is:

    .. math:: A = \\frac{\\pi D^2}{4}

    where:
      - `A` is the area of the circle,
      - `D` is the diameter of the circle,
      - `pi` is a mathematical constant, approximately 3.14159.

    Examples
    --------
    >>> circle_area(4 * ureg.meter)
    <Quantity(12.566370614359172, 'meter ** 2')>

    >>> circle_area(2.5 * ureg.inch)
    <Quantity(4.908738521234052, 'inch ** 2')>

    """
    return pi * D**2 / 4


def cylinder_volume(D, H):
    """Calculate the volume of a cylinder from its diameter and height.

    Parameters
    ----------
    D : ureg.Quantity {length: 1}
        The diameter of the cylinder.
    H : ureg.Quantity {length: 1}
        The height of the cylinder.

    Returns
    -------
    ureg.Quantity {length: 3}
        The volume enclosed by the cylinder.

    Notes
    -----
    The formula for calculating the volume of a cylinder is:

    .. math:: V = \\frac{\\pi D^2 H}{4}

    where:
      - `V` is the volume of the cylinder,
      - `D` is the diameter of the cylinder,
      - `H` is the height of the cylinder,
      - `pi` is a mathematical constant, approximately 3.14159.

    Examples
    --------
    >>> cylinder_volume(4 * ureg.meter, 6 * ureg.meter)
    <Quantity(75.39822368615503, 'meter ** 3')>

    >>> cylinder_volume(2.5 * ureg.inch, 8 * ureg.inch)
    <Quantity(19.634954084936208, 'inch ** 3')>

    """
    return pi * D**2 * H / 4
