TABLE_2 = {
    'Globe, single port': {
        '3 V-port plug': (0.9, 0.7, 0.48),  # F_L, x_T, F_d
        '4 V-port plug': (0.9, 0.7, 0.41),
        '6 V-port plug': (0.9, 0.7, 0.3),
        'Contoured plug': {'open': (0.9, 0.72, 0.46),
                           'close': (0.8, 0.55, 1.0)},
        '60 equal diameter': (0.9, 0.68, 0.13),
        '120 equal diameter': (0.9, 0.68, 0.09),
        'Characterized cage, 4-port':{'outward': (0.9, 0.75, 0.41),
                                      'inward': (0.85, 0.7, 0.41)}
    },
    'Globe, double port': {
        'Ported plug': (0.9, 0.75, 0.28),
        'Contoured plug': (0.85, 0.7, 0.32)
    },
    'Globe, angle': {
        'Contoured plug': {
            'open': (0.9, 0.72, 0.46),
            'close': (0.8, 0.65, 1.00)
        },
        'Characterized cage, 4-port': {
            'outward': (0.9, 0.65, 0.41),
            'inward': (0.85, 0.6, 0.41)
        },
        'Venturi': (0.5, 0.2, 1.0),
    },
    'Globe, small flow trim': {
        'V-notch': (0.98, 0.84, 0.7),
        'Flat seat': (0.85, 0.7, 0.3),
        'Tapered needle': (0.95, 0.84, 'N/A')  # See ISA-75 for F_d
    },
    'Rotary': {
        'Eccentric spherical plug': {
            'open': (0.85, 0.6, 0.42),
            'close': (0.68, 0.4, 0.42)
        },
        'Eccentric conical plug': {
            'open': (0.77, 0.54, 0.44),
            'close': (0.79, 0.55, 0.44)
        }
    },
    'Butterfly, centered shaft': {
        'Swing-through 70 deg': (0.62, 0.35, 0.57),
        'Swing-through 60 deg': (0.70, 0.42, 0.5),
        'Fluted vane 70 deg': (0.67, 0.38, 0.3)
    },
    'Butterfly, eccentric shaft': {
        'Offset seat': (0.67, 0.35, 0.57)
    },
    'Ball': {
        'Full bore 70 deg': (0.74, 0.42, 0.99),
        'Segmented ball': (0.6, 0.3, 0.98)
    }
}
