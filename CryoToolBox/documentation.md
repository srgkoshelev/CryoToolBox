CryoToolBox provides utility functions for variety of engineering
calculations, including fluid dynamic, heat transfer, and oxygen
deficiency analyses. The setup is simple:

``` python
>>> import CryoToolBox as ctb
>>> u = ctb.ureg  # Unit registry that handles unit conversions
```

Easy to create and get dimensions of pipes and tubes. Inputs in
different units is supported:

``` python
>>> D = ctb.Q_('2 inch')  # Nominal Diameter (here in inch)
>>> SCH = 10
>>> L = 10 * u.m  # Length (here in meter)
>>> pipe_2in = ctb.piping.Pipe(D, SCH, L)   # Creating the pipepipe_2in = ctb.piping.Pipe(D, SCH, L)   # Creating the pipe
>>> print(pipe_2in.ID)  # Pipe ID in inchprint(pipe_2in.ID)  # Pipe ID in inch
2.157 inch
>>> print(pipe_2in.ID.to(u.m))  # Pipe ID in m
0.0547878 meter
>>> print(pipe_2in.area)  # Pipe Area
3.6541819795329746 inch ** 2
>>> pipe_1in = ctb.piping.Pipe(1, SCH, L)  # Nominal diameter can be a number
>>> print(pipe_1in.OD)
1.315 inch
```

Tubes can have custom dimensions:

``` python
>>> tube = ctb.piping.Tube(2*u.inch, wall=0.035*u.inch)
>>> print(tube.ID)  # Tube ID
1.93 inch
```

Copper tubes are also supported:

``` python
>>> c_tube = ctb.piping.CopperTube(2, type_='K')
>>> print(c_tube.ID)
1.959 inch
```

Under the hood CoolProp calculates fluid properties, no matter what
units are used to define the thermodynamic state:

``` python
>>> fluid = ctb.ThermState('helium', T=300*u.K, P=3*u.bar)  # Fluid state (here helium at T=300K and P=3 bars)
>>> print(fluid.Dmass)  # Mass density
0.4807183970599733 kilogram / meter ** 3
>>> print(fluid.conductivity)  # Conductivity
0.1561263836300348 watt / kelvin / meter
>>> print(fluid.viscosity)  # Viscosity
1.993663632711359e-05 pascal * second
>>> print(fluid.Prandtl)  # Prandtl
0.6631570139483488
```

Additionally, pressure drop calculations for other piping elements,
e.g., elbows, tees, or valves, are supported:

``` python
>>> Cv = 10   # Cv of the valve
>>> valve = ctb.piping.Valve(D, Cv)
>>> print(valve.K())
142.46046037233918 dimensionless
>>> dP = ctb.piping.dP_adiab(m_dot, fluid, valve)  # Pressure drop for adiabatic compressible flow
>>> print(dP)
1282.9654742232524 pascal
```

The package also supports calculation of pressure drops for
incompressible flow, isothermal or adiabatic flows. Here\'s an example
of calculation for helium running through 1\" NPS pipe (both defined
above):

``` python
>>> m_dot = ctb.Q_('10 g/s')  # Mass flow (here in g/s)
>>> dP_isot = ctb.piping.dP_isot(m_dot, fluid, pipe_1in)  # Pressure drop for isothermal compressible flow
>>> print(dP_isot)
2896.1987828552374 pascal
>>> dP_adiab = ctb.piping.dP_adiab(m_dot, fluid, pipe_1in) # Pressure drop for adiabatic flow
>>> print(dP_adiab)
2632.691405886493 pascal
```

CryoToolBox can calculate some other useful quantities for manual like
checks, e.g., Reynolds number or Darcy friction factor:

``` python
>>> Re = ctb.Re(fluid, m_dot, pipe_1in.ID, pipe_1in.area)  # Reynolds number
>>> print(Re)
22920.172805808736
>>> eps = 0.000015 * u.m  # Pipe rugosity (here stainless steel)
>>> eps_r = eps/pipe_1in.ID
>>> f_Darcy = ctb.piping.f_Darcy(Re, eps_r)
>>> print(f_Darcy)
0.026363992431440168
```

The inverse problem of maximum flow rate calculation is also handled:

``` python
>>> P_out = 2.999 * u.bar
>>> m_dot = ctb.piping.m_dot_isot(fluid, pipe_1in, P_out)
>>> print(m_dot)
1.5228117452427699 gram / second
```

The package also supports pressure relief sizing, for example API-520
recommended direct integration method.

``` python
>>> G = ctb.piping.G_nozzle(fluid, P_out)  # Maximum mass flux with direct integration method (kg/sꞏm2)
>>> # Maximum mass flux with direct integration method (here the orifce is the pipe area) with Kb=Kc=Kd=K=1
>>> m_dot = G * pipe_1in.area
>>> print(m_dot.to(u.g/u.s))
5.978136976441524 gram / second
>>> A_orifice = m_dot/G   # Orifice size with Kb=Kc=Kd=Kv=1
>>> print(A_orifice)
0.9451552184159598 inch ** 2
```
