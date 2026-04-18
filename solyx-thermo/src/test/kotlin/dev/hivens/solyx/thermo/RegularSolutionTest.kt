package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.PhysicalConstants.R
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.ln

class RegularSolutionTest : DescribeSpec({

    val t1200 = Kelvin(1200.0)

    // -------------------------------------------------------------------------
    // Analytical results
    // -------------------------------------------------------------------------

    describe("analytical results") {

        it("omega=0 — identical to IdealMixing") {
            val regular = RegularSolution(interactions = emptyList())
            val ideal   = IdealMixing()

            val state = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            regular.compute(state).value.shouldBeWithinPercentageOf(
                ideal.compute(state).value, 1e-6
            )
        }

        it("negative omega — G is more negative than ideal (attraction)") {
            val omega = -10_000.0
            val regular = regularWithOmega(omega)
            val ideal   = IdealMixing()

            val state = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            regular.compute(state).value shouldBeLessThan ideal.compute(state).value
        }

        it("positive omega — G is less negative than ideal (repulsion)") {
            val omega = 10_000.0
            val regular = regularWithOmega(omega)
            val ideal   = IdealMixing()

            val state = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            regular.compute(state).value shouldBeGreaterThan ideal.compute(state).value
        }

        it("binary 50/50 with L0 — known excess value") {
            val omega   = -10_000.0
            val regular = regularWithOmega(omega)
            val state   = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            val g       = regular.compute(state).value

            // G = R*T*ln(0.5) + 0.5*0.5*omega
            val expected = R * 1200.0 * ln(0.5) + 0.5 * 0.5 * omega
            g.shouldBeWithinPercentageOf(expected, 1e-6)
        }

        it("pure component — excess is zero") {
            val regular = regularWithOmega(-10_000.0)
            val state   = binaryState(Element.Fe, 1.0, Element.C, 0.0, t1200)
            // x_A * x_B = 1.0 * 0.0 = 0.0 — no excess
            regular.compute(state).value.shouldBeWithinPercentageOf(0.0, 1e-6)
        }

        it("excess is symmetric — G(x) = G(1-x) for L0 only") {
            val regular = regularWithOmega(-10_000.0)
            val state1  = binaryState(Element.Fe, 0.3, Element.C, 0.7, t1200)
            val state2  = binaryState(Element.Fe, 0.7, Element.C, 0.3, t1200)
            regular.compute(state1).value.shouldBeWithinPercentageOf(
                regular.compute(state2).value, 1e-6
            )
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("G is continuous across all compositions") {
            val regular = regularWithOmega(-10_000.0)
            val values  = (1..99).map { i ->
                val x = i / 100.0
                regular.compute(binaryState(Element.Fe, x, Element.C, 1.0 - x, t1200)).value
            }
            // No NaN or Infinity in the sequence
            values.forEach { g ->
                g.isNaN()      shouldBe false
                g.isInfinite() shouldBe false
            }
        }
    }

    // -------------------------------------------------------------------------
    // Edge cases
    // -------------------------------------------------------------------------

    describe("edge cases") {

        it("very large negative omega — no overflow") {
            val regular = regularWithOmega(-1_000_000.0)
            val state   = binaryState(Element.Fe, 0.5, Element.C, 0.5, t1200)
            val g = regular.compute(state).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }

        it("near-zero mole fraction — no NaN") {
            val regular = regularWithOmega(-10_000.0)
            val state   = binaryState(Element.Fe, 1e-9, Element.C, 1.0 - 1e-9, t1200)
            val g = regular.compute(state).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }
    }
})

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------

private fun regularWithOmega(omega: Double) = RegularSolution(
    interactions = listOf(
        RedlichKisterParameter(
            a = Element.Fe,
            b = Element.C,
            coefficients = listOf(JoulePerMole(omega))
        )
    )
)

private fun binaryState(
    a: Element, xa: Double,
    b: Element, xb: Double,
    temperature: Kelvin
) = SystemState(
    composition = mapOf(a to MoleFraction(xa), b to MoleFraction(xb)),
    temperature = temperature
)

private infix fun Double.shouldBeLessThan(other: Double) {
    (this < other) shouldBe true
}

private infix fun Double.shouldBeGreaterThan(other: Double) {
    (this > other) shouldBe true
}
