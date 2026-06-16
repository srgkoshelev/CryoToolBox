"""Minimal average-property outlet-state propagation for CryoToolBox."""

from scipy.optimize import root

import CryoToolBox as ctb


u = ctb.ureg


def propagate_properties(m_dot, fluid_i, pipe, Q_dot=0 * u.W, *, tol=1e-8):
    """Return outlet state for one horizontal constant-area element.

    Solves outlet ``P,Hmass`` using midpoint properties for ``K(Re)`` and
    endpoint specific volumes for acceleration and kinetic-energy change:

    ``P1 - P2 = K_m G**2 (v1 + v2) / 4 + G**2 (v2 - v1)``

    ``h2 - h1 = Q_dot / m_dot - G**2 (v2**2 - v1**2) / 2``
    """

    # 0.40 is the common average-property pressure-drop limit from Rennels.
    # The volume and quality limits are step-size controls for property
    # curvature/two-phase resolution, not hard physical limits.
    max_dP_frac = 0.40
    max_dv_frac = 0.25
    max_dx = 0.05

    if m_dot < 0 * u.kg / u.s:
        raise ctb.piping.HydraulicError('m_dot must be non-negative.')
    if m_dot == 0 * u.kg / u.s:
        if Q_dot != 0 * u.W:
            raise ctb.piping.HydraulicError(
                'Cannot propagate nonzero heat input with m_dot = 0.')
        fluid_o = fluid_i.copy()
        fluid_o.update_kw(P=fluid_i.P, Hmass=fluid_i.Hmass)
        return fluid_o

    P1 = fluid_i.P
    h1 = fluid_i.Hmass
    v1 = 1 / fluid_i.Dmass
    G = m_dot / pipe.area
    q = Q_dot / m_dot

    dP_guess = pipe.K(ctb.Re(fluid_i, m_dot, pipe.ID, pipe.area)) * G**2 * v1 / 2
    dP_guess = min(dP_guess, 0.2 * P1)
    h_scale = max(abs(q), abs(G**2 * v1**2 / 2), 1 * u.J / u.kg)

    def state(P, h):
        fluid = fluid_i.copy()
        fluid.update_kw(P=P, Hmass=h)
        return fluid

    def residual(x):
        P2 = x[0] * u.Pa
        h2 = x[1] * u.J / u.kg
        if P2 <= 0 * u.Pa:
            return [1e6, 1e6]

        try:
            fluid_o = state(P2, h2)
            fluid_m = state((P1 + P2) / 2, (h1 + h2) / 2)
            v2 = 1 / fluid_o.Dmass
            K_m = pipe.K(ctb.Re(fluid_m, m_dot, pipe.ID, pipe.area))
            dP_loss = K_m * G**2 * (v1 + v2) / 4
            dP_acc = G**2 * (v2 - v1)
            dke = G**2 * (v2**2 - v1**2) / 2
            return [
                float(((P1 - P2) - dP_loss - dP_acc) / P1),
                float(((h2 - h1) - q + dke) / h_scale),
            ]
        except Exception:
            return [1e6, 1e6]

    solution = root(
        residual,
        [(P1 - dP_guess).m_as(u.Pa), (h1 + q).m_as(u.J / u.kg)],
        tol=tol)
    err = max(abs(r) for r in residual(solution.x))
    if not solution.success or err > tol:
        raise ctb.piping.HydraulicError(
            f'Outlet-state propagation did not converge: {solution.message!r}.')

    fluid_o = state(solution.x[0] * u.Pa, solution.x[1] * u.J / u.kg)
    dP = P1 - fluid_o.P
    dv_frac = abs(float((1 / fluid_o.Dmass - v1) / v1))

    if isinstance(pipe, ctb.piping.Tube):
        if abs(float(dP / P1)) > max_dP_frac:
            raise ctb.piping.HydraulicError(
                'Average-property pressure change is too large; split the element.')
        if dv_frac > max_dv_frac:
            raise ctb.piping.HydraulicError(
                'Average-property specific-volume change is too large; '
                'split the element.')
        if fluid_o.phase != fluid_i.phase:
            raise ctb.piping.HydraulicError(
                'Outlet state changed phase; split the element at the boundary.')
        if fluid_i.phase == 6 and abs(fluid_o.Q - fluid_i.Q) > max_dx:
            raise ctb.piping.HydraulicError(
                'Two-phase quality change is too large; split the element.')

    return fluid_o
