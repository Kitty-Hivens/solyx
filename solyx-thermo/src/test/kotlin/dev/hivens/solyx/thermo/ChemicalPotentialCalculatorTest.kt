package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.PhysicalConstants.R
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.ln

class ChemicalPotentialCalculatorTest : DescribeSpec({

    val t1200 = Kelvin(1200.0)
    val ideal = IdealMixing()
    val calc  = ChemicalPotentialCalculator(ideal)

    // -------------------------------------------------------------------------
    // Analytical results — ideal mixing has known closed-form μ_i = R*T*ln(x_i)
    // -------------------------------------------------------------------------

    describe("analytical results — ideal mixing") {

        it("binary — μ_Fe = R*T*ln(x_Fe)") {
            val state = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            val mu = calc.compute(state).first { it.element == Element.Fe }
            val expected = R * 1200.0 * ln(0.3)
            mu.value.value.shouldBeWithinPercentageOf(expected, 0.01)
        }

        it("binary — μ_C = R*T*ln(x_C)") {
            val state = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            val mu = calc.compute(state).first { it.element == Element.C }
            val expected = R * 1200.0 * ln(0.7)
            mu.value.value.shouldBeWithinPercentageOf(expected, 0.01)
        }

        it("ternary — μ_i = R*T*ln(x_i) for each component") {
            val state = ternaryState(
                Element.Fe, 0.5,
                Element.Cr, 0.3,
                Element.Ni, 0.2,
                t1200
            )
            val potentials = calc.compute(state).associateBy { it.element }

            potentials[Element.Fe]!!.value.value
                .shouldBeWithinPercentageOf(R * 1200.0 * ln(0.5), 0.01)
            potentials[Element.Cr]!!.value.value
                .shouldBeWithinPercentageOf(R * 1200.0 * ln(0.3), 0.01)
            potentials[Element.Ni]!!.value.value
                .shouldBeWithinPercentageOf(R * 1200.0 * ln(0.2), 0.01)
        }

        it("quaternary — μ_i = R*T*ln(x_i) for each component") {
            val state = quaternaryState(t1200)
            val potentials = calc.compute(state).associateBy { it.element }

            mapOf(
                Element.Fe to 0.4,
                Element.Cr to 0.3,
                Element.Ni to 0.2,
                Element.Mo to 0.1
            ).forEach { (el, x) ->
                potentials[el]!!.value.value
                    .shouldBeWithinPercentageOf(R * 1200.0 * ln(x), 0.01)
            }
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants — must hold for ANY model, not just ideal mixing
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("Euler equation — G = Σ x_i * μ_i — binary") {
            val state = binaryState(Element.Fe, 0.4, Element.C, 0.6, t1200)
            val g = ideal.compute(state).value
            val potentials = calc.compute(state)
            val eulerSum = potentials.sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.01)
        }

        it("Euler equation — G = Σ x_i * μ_i — ternary") {
            val state = ternaryState(
                Element.Fe, 0.5, Element.Cr, 0.3, Element.Ni, 0.2, t1200
            )
            val g = ideal.compute(state).value
            val potentials = calc.compute(state)
            val eulerSum = potentials.sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.01)
        }

        it("Euler equation — G = Σ x_i * μ_i — quaternary") {
            val state = quaternaryState(t1200)
            val g = ideal.compute(state).value
            val potentials = calc.compute(state)
            val eulerSum = potentials.sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.01)
        }

        it("Euler equation holds for regular solution — binary") {
            val model = RegularSolution(
                interactions = listOf(
                    RedlichKisterParameter(
                        Element.Fe, Element.C,
                        listOf(JoulePerMole(-15_000.0), JoulePerMole(3_000.0))
                    )
                )
            )
            val calc2 = ChemicalPotentialCalculator(model)
            val state = binaryState(Element.Fe, 0.35, Element.C, 0.65, t1200)
            val g = model.compute(state).value
            val eulerSum = calc2.compute(state).sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.01)
        }

        it("Euler equation holds for regular solution — ternary") {
            val model = RegularSolution(
                interactions = listOf(
                    RedlichKisterParameter(
                        Element.Fe, Element.Cr,
                        listOf(JoulePerMole(-8_000.0))
                    ),
                    RedlichKisterParameter(
                        Element.Fe, Element.Ni,
                        listOf(JoulePerMole(-5_000.0))
                    )
                )
            )
            val calc3 = ChemicalPotentialCalculator(model)
            val state = ternaryState(
                Element.Fe, 0.6, Element.Cr, 0.25, Element.Ni, 0.15, t1200
            )
            val g = model.compute(state).value
            val eulerSum = calc3.compute(state).sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.01)
        }

        it("μ_i approaches -infinity as x_i -> 0 (Raoult's law)") {
            // For ideal mixing μ_i = R*T*ln(x_i) → -∞ as x_i → 0
            val dilute = SystemState(
                composition = mapOf(
                    Element.Fe to MoleFraction(1.0 - 1e-6),
                    Element.C  to MoleFraction(1e-6)
                ),
                temperature = t1200
            )
            val mu_C = calc.compute(dilute).first { it.element == Element.C }
            // Must be very negative
            (mu_C.value.value < -100_000.0) shouldBe true
        }

        it("μ_i is independent of labelling — swapping components gives same values") {
            val state1 = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            val state2 = binaryState(Element.C, 0.7, Element.Fe, 0.3, t1200)
            val p1 = calc.compute(state1).associateBy { it.element }
            val p2 = calc.compute(state2).associateBy { it.element }
            p1[Element.Fe]!!.value.value.shouldBeWithinPercentageOf(
                p2[Element.Fe]!!.value.value, 0.01
            )
            p1[Element.C]!!.value.value.shouldBeWithinPercentageOf(
                p2[Element.C]!!.value.value, 0.01
            )
        }

        it("result count matches component count") {
            val binary   = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            val ternary  = ternaryState(Element.Fe, 0.5, Element.Cr, 0.3, Element.Ni, 0.2, t1200)
            val quaternary = quaternaryState(t1200)

            calc.compute(binary).size    shouldBe 2
            calc.compute(ternary).size   shouldBe 3
            calc.compute(quaternary).size shouldBe 4
        }
    }

    // -------------------------------------------------------------------------
    // Numerical stability — extreme compositions (safe step-size fix)
    // -------------------------------------------------------------------------

    describe("numerical stability") {

        it("extreme composition x=0.999 — no crash, Euler holds") {
            val state = SystemState(
                composition = mapOf(
                    Element.Fe to MoleFraction(0.999),
                    Element.C  to MoleFraction(0.001)
                ),
                temperature = t1200
            )
            val g = ideal.compute(state).value
            val eulerSum = calc.compute(state).sumOf { mu ->
                (state.composition[mu.element]?.value ?: 0.0) * mu.value.value
            }
            eulerSum.shouldBeWithinPercentageOf(g, 0.1)
        }

        it("extreme composition x=0.001 — no crash") {
            val state = SystemState(
                composition = mapOf(
                    Element.Fe to MoleFraction(0.001),
                    Element.C  to MoleFraction(0.999)
                ),
                temperature = t1200
            )
            val potentials = calc.compute(state)
            potentials.forEach { mu ->
                mu.value.value.isNaN()      shouldBe false
                mu.value.value.isInfinite() shouldBe false
            }
        }

        it("all compositions across binary range — no NaN or Infinity") {
            (1..99).forEach { i ->
                val x = i / 100.0
                val state = binaryState(Element.Fe, x, Element.C, 1.0 - x, t1200)
                calc.compute(state).forEach { mu ->
                    mu.value.value.isNaN()      shouldBe false
                    mu.value.value.isInfinite() shouldBe false
                }
            }
        }
    }
})

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------

private fun binaryState(
    a: Element, xa: Double,
    b: Element, xb: Double,
    temperature: Kelvin
) = SystemState(
    composition = mapOf(a to MoleFraction(xa), b to MoleFraction(xb)),
    temperature = temperature
)

private fun ternaryState(
    a: Element, xa: Double,
    b: Element, xb: Double,
    c: Element, xc: Double,
    temperature: Kelvin
) = SystemState(
    composition = mapOf(
        a to MoleFraction(xa),
        b to MoleFraction(xb),
        c to MoleFraction(xc)
    ),
    temperature = temperature
)

private fun quaternaryState(temperature: Kelvin) = SystemState(
    composition = mapOf(
        Element.Fe to MoleFraction(0.4),
        Element.Cr to MoleFraction(0.3),
        Element.Ni to MoleFraction(0.2),
        Element.Mo to MoleFraction(0.1)
    ),
    temperature = temperature
)
