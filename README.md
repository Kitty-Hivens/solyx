# Solyx

Computational thermodynamics library for the JVM.
CALPHAD-based phase equilibrium calculations for materials science and engineering.

> **Pre-release.** Core numerical methods and thermodynamic models are implemented.
> Results are not yet verified against experimental data. Do not use in production.

---

## What it does

Solyx computes thermodynamic equilibrium in multicomponent systems — which phases are
stable, their compositions, and their fractions — using the CALPHAD methodology.

Targets researchers and engineers who need thermodynamic calculations on the JVM
without depending on proprietary software (Thermo-Calc, Pandat) or Python.

---

## Installation

Not yet published to Maven Central. Currently available via source only.

```kotlin
// settings.gradle.kts
include(":solyx-core", ":solyx-thermo", ":solyx-calphad")
```

Maven Central publication is planned for `v1.0.0` after experimental verification.

---

## Modules

| Module          | Description                                 | Status      |
|-----------------|---------------------------------------------|-------------|
| `solyx-core`    | Numerical methods, units, error tracking    | Complete    |
| `solyx-thermo`  | Thermodynamic models                        | Complete    |
| `solyx-calphad` | TDB parser, Gibbs minimizer, phase diagrams | In progress |

---

## Quick start

```kotlin
import dev.hivens.solyx.thermo.*
import dev.hivens.solyx.core.units.*

// Binary Fe-C system at 1200°C
val state = SystemState(
    composition = mapOf(
        Element.Fe to MoleFraction(0.98),
        Element.C  to MoleFraction(0.02)
    ),
    temperature = Kelvin.fromCelsius(1200.0)
)

// Ideal mixing entropy
val mixing = IdealMixing()
val g = mixing.compute(state)
println(g) // -975.3 J/mol

// Regular solution with Fe-C interaction
val solution = RegularSolution(
    interactions = listOf(
        RedlichKisterParameter(
            a = Element.Fe,
            b = Element.C,
            coefficients = listOf(JoulePerMole(-10_000.0))
        )
    )
)
val gTotal = solution.compute(state)

// Magnetic contribution for BCC Fe
val gMag = MagneticContribution.BCC_FE.compute(Kelvin.fromCelsius(1200.0))
```

---

## Design

- `double` precision throughout — sufficient for thermodynamic data accuracy
- `Interval` arithmetic for propagated uncertainty — every result carries error bounds
- Strongly typed units — `Kelvin`, `JoulePerMole`, `MoleFraction` — physically incorrect
  input is rejected at compile time or construction time
- Zero external dependencies in `solyx-core`
- Java-compatible API via `@JvmStatic` and explicit overloads

---

## Verification

Before `v1.0.0`, results will be verified against:

- Binary Fe-C phase diagram (NIST data)
- Binary Al-Cu phase diagram
- Liquidus/solidus temperatures vs published experimental data

Uncertainty bounds are reported on all results via `Interval`. If `result.isReliable`
is false — the calculation did not converge to sufficient precision.

---

## License

MIT — no restrictions on commercial or academic use.

## Roadmap

See [ROADMAP.md](ROADMAP.md) for the full implementation plan and release schedule.
