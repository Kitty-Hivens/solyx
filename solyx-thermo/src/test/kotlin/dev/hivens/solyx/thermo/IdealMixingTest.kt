package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.PhysicalConstants.R
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeLessThan
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.ln

class IdealMixingTest : DescribeSpec({

    val mixing = IdealMixing()
    val t1200  = Kelvin(1200.0)

    // -------------------------------------------------------------------------
    // Analytical results
    // -------------------------------------------------------------------------

    describe("analytical results") {

        it("binary 50/50 — maximum entropy, known value") {
            val state = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            val g = mixing.compute(state).value

            // G_mix = R * T * 2 * (0.5 * ln(0.5)) = R * T * ln(0.5)
            val expected = R * 1200.0 * ln(0.5)
            g.shouldBeWithinPercentageOf(expected, 1e-6)
        }

        it("binary 50/50 — symmetric, order of components does not matter") {
            val state1 = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            val state2 = binaryState(Element.C, 0.5, Element.Fe, 0.5, t1200)
            mixing.compute(state1).value.shouldBeWithinPercentageOf(
                mixing.compute(state2).value, 1e-10
            )
        }

        it("binary 80/20 — known value") {
            val state = binaryState(Element.Fe, 0.8, Element.C, 0.2, t1200)
            val g = mixing.compute(state).value
            val expected = R * 1200.0 * (0.8 * ln(0.8) + 0.2 * ln(0.2))
            g.shouldBeWithinPercentageOf(expected, 1e-6)
        }

        it("pure component x=1.0 — G_mix is zero") {
            val state = binaryState(Element.Fe, 1.0, Element.C, 0.0, t1200)
            mixing.compute(state).value.shouldBeWithinPercentageOf(0.0, 1e-6)
        }

        it("ternary equal fractions — known value") {
            val state = ternaryState(
                Element.Fe, 1.0 / 3.0,
                Element.Cr, 1.0 / 3.0,
                Element.Ni, 1.0 / 3.0,
                t1200
            )
            val g = mixing.compute(state).value
            // G_mix = R*T * 3 * (1/3 * ln(1/3)) = R*T*ln(1/3)
            val expected = R * 1200.0 * ln(1.0 / 3.0)
            g.shouldBeWithinPercentageOf(expected, 1e-6)
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("G_mix is always negative for mixed system") {
            listOf(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9).forEach { x ->
                val state = binaryState(Element.Fe, x, Element.C, 1.0 - x, t1200)
                mixing.compute(state).value shouldBeLessThan 0.0
            }
        }

        it("G_mix scales linearly with temperature") {
            val state = binaryState(Element.Fe, 0.5, Element.C, 0.5, Kelvin(600.0))
            val g600  = mixing.compute(state).value
            val g1200 = mixing.compute(
                binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            ).value
            // G_mix(1200) / G_mix(600) = 1200 / 600 = 2.0
            (g1200 / g600).shouldBeWithinPercentageOf(2.0, 1e-6)
        }

        it("magnitude increases with temperature") {
            val state600  = binaryState(Element.Fe, 0.3, Element.C, 0.7, Kelvin(600.0))
            val state1200 = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            val g600  = mixing.compute(state600).value
            val g1200 = mixing.compute(state1200).value
            // Both negative — higher T means more negative
            (g1200 < g600) shouldBe true
        }
    }

    // -------------------------------------------------------------------------
    // Edge cases
    // -------------------------------------------------------------------------

    describe("edge cases") {

        it("near-zero mole fraction — no NaN or Infinity") {
            val state = binaryState(Element.Fe, 1e-10, Element.C, 1.0 - 1e-10, t1200)
            val g = mixing.compute(state).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }

        it("very low temperature — G_mix approaches zero") {
            val state = binaryState(Element.Fe, 0.5, Element.C, 0.5, Kelvin(1.0))
            val g = mixing.compute(state).value
            // G_mix = R * 1.0 * ln(0.5) ≈ -5.76 J/mol — small but not zero
            g.shouldBeWithinPercentageOf(R * 1.0 * ln(0.5), 1e-6)
        }

        it("single component — returns zero without error") {
            val state = SystemState(
                composition = mapOf(Element.Fe to MoleFraction(1.0)),
                temperature = t1200
            )
            mixing.compute(state).value.shouldBeWithinPercentageOf(0.0, 1e-6)
        }
    }

    // -------------------------------------------------------------------------
    // Invalid input
    // -------------------------------------------------------------------------

    describe("invalid input") {

        it("composition not summing to 1.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                SystemState(
                    composition = mapOf(
                        Element.Fe to MoleFraction(0.5),
                        Element.C  to MoleFraction(0.3)
                    ),
                    temperature = t1200
                )
            }
        }

        it("negative temperature — throws") {
            shouldThrow<IllegalArgumentException> {
                Kelvin(-1.0)
            }
        }

        it("mole fraction above 1.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                MoleFraction(1.1)
            }
        }

        it("mole fraction below 0.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                MoleFraction(-0.1)
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
