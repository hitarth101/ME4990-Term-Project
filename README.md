# Venus Orbiter Thermal Control Design

ME4550 Spacecraft Thermal Controls term project (Fall 2026).

Thermal-control design for a Venus orbiter that maps the surface with radar at 10 m resolution and measures the gravity field. The mission follows Magellan:
- elliptical polar orbit (294 x 8,540 km) for mapping;
- aerobraking into a 250 km circular orbit for the gravity survey.

## Team

| Member | Work package |
|---|---|
| Member 1 | Mission and spacecraft configuration |
| Member 2 | Component requirements and power |
| Member 3 | Environment, orbit geometry, design cases |
| Member 4 | Heat paths, thermal model, verification |
| Member 5 | Hot/cold cases, heaters, transients |
| Member 6 | Design trades, margins, remaining work |

## Repository

```
src/venus_env.py      orbit geometry, eclipse, solar/albedo/IR flux on each face
src/thermal_model.py  5-node thermal network (steady state + RK4 transient with thermostats)
src/run_cases.py      model checks, design cases, radiator/heater sizing, trades -> results/
src/make_figures.py   plots -> figures/
results/              summary.md (all tables), results.json, transient time histories
figures/              plots used in the slides
slides/               presentation decks
```

## How to run

```
pip install -r requirements.txt
cd src
python run_cases.py
python make_figures.py
```

`run_cases.py` takes about 2 minutes and rewrites everything in `results/`.

## Model

| Node | Represents | C (kJ/K) |
|---|---|---|
| EL | Electronics boxes bolted to the radiator panel | 144 |
| RAD | -Y radiator panel, quartz-mirror OSR | 13.5 |
| BAT | Li-ion battery and its own radiator patch | 16.5 |
| PROP | Hydrazine tank + propellant | 257 |
| BUS | Bus structure behind MLI | 162 |

- **Environment:**
  - Gilmore Eq. 2.9 solar flux (2,579-2,650 W/m2);
  - albedo 0.80 +/- 0.02;
  - planetary IR 138.5-178.4 W/m2 (Gilmore Tables 2.9, 2.10).
  - View factors and albedo factors are integrated numerically over the visible planet.
- **Attitude:** nadir-pointing, with a yaw flip at beta = 0 so the -Y radiator never sees the Sun.
- **Couplings:** box-to-radiator 20 W/K; battery and tank isolated (0.3, 0.4 W/K); MLI effective emittance 0.03.
- **Coatings:** OSR alpha 0.08 BOL / 0.20 EOL, emittance 0.80; Beta-cloth MLI cover alpha 0.32 / 0.55.
- **Margins:** 11 C on passively controlled limits; 25% heater capacity at the 24 V minimum bus.

## Checks

- **View factor:** nadir view factor matches (R/r)^2.
- **Eclipse:** duration matches the analytic value.
- **Reference sphere:** black reference sphere at 607 km matches Gilmore Table 2.11 within 1 C.
- **Energy balance:** residual over one orbit is about 2 J out of about 2 MJ.

## Current results

| Item | Value |
|---|---|
| Hot case | 250 km circular, beta = 30 deg, perihelion, EOL |
| Radiator | 1.20 m2 OSR on -Y + 0.05 m2 battery radiator |
| Worst cold case | Safe mode in cruise near 1 AU (211 W heater demand) |
| Installed heaters (24 V) | Electronics 230 W, battery 22 W, tank 16 W |
| Smallest margin | Battery hot side, +13.9 C |

Full tables: `results/summary.md`.

## References

- D. G. Gilmore (ed.), *Spacecraft Thermal Control Handbook, Vol. I*, 2nd ed., Aerospace Press/AIAA, 2002.
- NASA NSSDC Venus Fact Sheet; NSSDC Magellan Mission Plan.
- ME4550 course lectures and worked examples.
