package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.PhysicalConstants.R
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.ln

class IonicLiquidTest : DescribeSpec({

    val t1500 = Kelvin(1500.0)

    val fe2 = IonicSpecies.Cation(Element.Fe, 2)
    val mn2 = IonicSpecies.Cation(Element.Mn, 2)
    val o2  = IonicSpecies.Anion(Element.O,  2)
    val va  = IonicSpecies.Vacancy

    // -------------------------------------------------------------------------
    // Type safety — compile-time and init-time checks
    // -------------------------------------------------------------------------

    describe("type safety") {

        it("AnionInteraction with Cation as a — throws IllegalArgumentException") {
            shouldThrow<IllegalArgumentException> {
                AnionInteraction(
                    a = fe2,  // Cation, not allowed
                    b = o2,
                    parameter = RedlichKisterParameter(Element.Fe, Element.O,
                        listOf(JoulePerMole(-5000.0)))
                )
            }
        }

        it("AnionInteraction with Cation as b — throws IllegalArgumentException") {
            shouldThrow<IllegalArgumentException> {
                AnionInteraction(
                    a = o2,
                    b = mn2,  // Cation, not allowed
                    parameter = RedlichKisterParameter(Element.O, Element.Mn,
                        listOf(JoulePerMole(-5000.0)))
                )
            }
        }

        it("AnionInteraction with Vacancy is valid") {
            // Should not throw
            AnionInteraction(
                a = o2,
                b = va,
                parameter = RedlichKisterParameter(Element.O, Element.Si,
                    listOf(JoulePerMole(-3000.0)))
            )
        }

        it("CationInteraction requires Cation typed at compile time — instantiates cleanly") {
            // Compiler enforces Cation types — just verify construction
            CationInteraction(
                a = fe2,
                b = mn2,
                parameter = RedlichKisterParameter(Element.Fe, Element.Mn,
                    listOf(JoulePerMole(-10000.0)))
            )
        }
    }

    // -------------------------------------------------------------------------
    // Analytical results
    // -------------------------------------------------------------------------

    describe("analytical results") {

        it("single cation single anion — G_ref only, no mixing") {
            val model = IonicLiquid(
                endmembers = listOf(
                    IonicEndmember(fe2, o2, JoulePerMole(-250_000.0))
                )
            )
            val state = IonicState(
                cationOccupancy = mapOf(fe2 to 1.0),
                anionOccupancy  = mapOf(o2  to 1.0),
                temperature = t1500
            )
            // G = 1*1*(-250000) + R*T*(P*0 + Q*0) = -250000
            model.compute(state).value
                .shouldBeWithinPercentageOf(-250_000.0, 1e-6)
        }

        it("two cations one anion — G_ref + G_mix on cation sublattice") {
            val model = IonicLiquid(
                endmembers = listOf(
                    IonicEndmember(fe2, o2, JoulePerMole(-250_000.0)),
                    IonicEndmember(mn2, o2, JoulePerMole(-290_000.0))
                )
            )
            val yFe = 0.4; val yCa = 0.6
            val state = IonicState(
                cationOccupancy = mapOf(fe2 to yFe, mn2 to yCa),
                anionOccupancy  = mapOf(o2  to 1.0),
                temperature = t1500
            )
            val g = model.compute(state).value

            // P = 1*2 = 2 (anion sublattice: one O2- with charge 2)
            // Q = yFe*2 + yCa*2 = 2
            val p = 2.0; val q = 2.0
            val gRef = yFe * (-250_000.0) + yCa * (-290_000.0)
            val gMix = R * 1500.0 * (
                    p * (yFe * ln(yFe) + yCa * ln(yCa)) +
                            q * 0.0  // single anion — no mixing
                    )
            g.shouldBeWithinPercentageOf(gRef + gMix, 1e-3)
        }

        it("single cation with anion and vacancy — G_mix on anion sublattice") {
            val model = IonicLiquid(
                endmembers = listOf(
                    IonicEndmember(fe2, o2, JoulePerMole(-250_000.0)),
                    IonicEndmember(fe2, va, JoulePerMole(-30_000.0))
                )
            )
            val yO = 0.7; val yVa = 0.3
            val state = IonicState(
                cationOccupancy = mapOf(fe2 to 1.0),
                anionOccupancy  = mapOf(o2 to yO, va to yVa),
                temperature = t1500
            )
            val g = model.compute(state).value

            // P = yO*2 + yVa*0 = 1.4
            // Q = 1.0*2 = 2.0
            val p = yO * 2.0
            val q = 2.0
            val gRef = 1.0 * yO * (-250_000.0) + 1.0 * yVa * (-30_000.0)
            val gMix = R * 1500.0 * (
                    p * 0.0 +  // single cation
                            q * (yO * ln(yO) + yVa * ln(yVa))
                    )
            g.shouldBeWithinPercentageOf(gRef + gMix, 1e-3)
        }

        it("P and Q electroneutrality constraint — P*Q_cation = Q*P_anion") {
            // For a neutral system: P = Σ y_j*|z_j|, Q = Σ y_i*z_i
            val yFe = 0.5; val yCa = 0.5
            val state = IonicState(
                cationOccupancy = mapOf(fe2 to yFe, mn2 to yCa),
                anionOccupancy  = mapOf(o2 to 1.0),
                temperature = t1500
            )
            state.p.shouldBeWithinPercentageOf(2.0, 1e-6)  // 1*2 = 2
            state.q.shouldBeWithinPercentageOf(2.0, 1e-6)  // 0.5*2 + 0.5*2 = 2
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("G is finite across cation composition range") {
            val model = IonicLiquid(
                endmembers = listOf(
                    IonicEndmember(fe2, o2, JoulePerMole(-250_000.0)),
                    IonicEndmember(mn2, o2, JoulePerMole(-290_000.0))
                )
            )
            (1..99).forEach { i ->
                val y = i / 100.0
                val state = IonicState(
                    cationOccupancy = mapOf(fe2 to y, mn2 to 1.0 - y),
                    anionOccupancy  = mapOf(o2 to 1.0),
                    temperature = t1500
                )
                val g = model.compute(state).value
                g.isNaN()      shouldBe false
                g.isInfinite() shouldBe false
            }
        }
    }

    // -------------------------------------------------------------------------
    // IonicState validation
    // -------------------------------------------------------------------------

    describe("IonicState validation") {

        it("cation fractions not summing to 1.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                IonicState(
                    cationOccupancy = mapOf(fe2 to 0.3, mn2 to 0.3), // sum=0.6
                    anionOccupancy  = mapOf(o2 to 1.0),
                    temperature = t1500
                )
            }
        }

        it("anion fractions not summing to 1.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                IonicState(
                    cationOccupancy = mapOf(fe2 to 1.0),
                    anionOccupancy  = mapOf(o2 to 0.5, va to 0.3), // sum=0.8
                    temperature = t1500
                )
            }
        }
    }
})
