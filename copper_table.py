"""Copper pipe nominal dimensions from the Copper Tube Handbook by CDA.
"""
from . import ureg
COPPER_TABLE = {
    1/4: {
        'OD': .375 * ureg.inch,
        'Type K': .035 * ureg.inch,
    },
    3/8: {
        'OD': .5 * ureg.inch,
        'Type K': .049 * ureg.inch,
    },
    1/2: {
        'OD': .625 * ureg.inch,
        'Type K': .049 * ureg.inch,
    },
    5/8: {
        'OD': .75 * ureg.inch,
        'Type K': .049 * ureg.inch,
    },
    3/4: {
        'OD': .875 * ureg.inch,
        'Type K': .065 * ureg.inch,
    },
    1: {
        'OD': 1.125 * ureg.inch,
        'Type K': .065 * ureg.inch,
    },
    1.25: {
        'OD': 1.375 * ureg.inch,
        'Type K': .065 * ureg.inch,
    },
    1.5: {
        'OD': 1.625 * ureg.inch,
        'Type K': .072 * ureg.inch,
    },
    2: {
        'OD': 2.125 * ureg.inch,
        'Type K': .083 * ureg.inch,
    },
    2.5: {
        'OD': 2.625 * ureg.inch,
        'Type K': .095 * ureg.inch,
    },
    3: {
        'OD': 3.125 * ureg.inch,
        'Type K': .109 * ureg.inch,
    },
    3.5: {
        'OD': 3.625 * ureg.inch,
        'Type K': .120 * ureg.inch,
    },
    4: {
        'OD': 4.125 * ureg.inch,
        'Type K': .134 * ureg.inch,
    },
    5: {
        'OD': 5.125 * ureg.inch,
        'Type K': .160 * ureg.inch,
    },
    6: {
        'OD': 6.125 * ureg.inch,
        'Type K': .192 * ureg.inch,
    },
    8: {
        'OD': 8.125 * ureg.inch,
        'Type K': .271 * ureg.inch,
    },
    10: {
        'OD': 10.125 * ureg.inch,
        'Type K': .338 * ureg.inch,
    },
    12: {
        'OD': 12.125 * ureg.inch,
        'Type K': .405 * ureg.inch,
    },
}
