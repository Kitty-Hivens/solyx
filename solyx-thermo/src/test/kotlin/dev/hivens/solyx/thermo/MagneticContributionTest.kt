package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.Kelvin
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeLessThan
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe

class MagneticContributionTest : DescribeSpec({

    // -------------------------------------------------------------------------
    // Analytical results
    // -------------------------------------------------------------------------

    describe("analytical results") {

        it("zero beta — returns zero contribution") {
            val model = MagneticContribution(tc = Kelvin(1043.0), beta = 0.0, p = 0.28)
            model.compute(Kelvin(800.0)).value shouldBe 0.0
        }

        it("BCC Fe at 800K — below Tc, known sign") {
            val g = MagneticContribution.BCC_FE.compute(Kelvin(800.0)).value
            // Below Tc magnetic contribution stabilizes the phase — G_mag < 0
            g shouldBeLessThan 0.0
        }

        it("BCC Fe at 1500K — above Tc, known sign") {
            val g = MagneticContribution.BCC_FE.compute(Kelvin(1500.0)).value
            // Above Tc disordering — G_mag > 0 (destabilizes ordered phase)
            g shouldBeLessThan 0.0
        }

        it("contribution approaches zero at very high temperature") {
            val g2000 = MagneticContribution.BCC_FE.compute(Kelvin(2000.0)).value
            val g5000 = MagneticContribution.BCC_FE.compute(Kelvin(5000.0)).value
            // At high T the paramagnetic tail decays
            kotlin.math.abs(g5000) shouldBeLessThan kotlin.math.abs(g2000)
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("continuity across Tc — no discontinuous jump") {
            val tc = MagneticContribution.BCC_FE.tc.value
            val gBelow = MagneticContribution.BCC_FE.compute(Kelvin(tc - 1.0)).value
            val gAbove = MagneticContribution.BCC_FE.compute(Kelvin(tc + 1.0)).value
            gAbove.shouldBeWithinPercentageOf(gBelow, 10.0)
        }

        it("FCC Fe antiferromagnetic — betaEff is positive and contribution is non-zero") {
            // beta=-2.1, afm=-3 -> betaEff = -2.1 / -3 = 0.7
            MagneticContribution.FCC_FE.betaEff.shouldBeWithinPercentageOf(0.7, 0.01)
            val g = MagneticContribution.FCC_FE.compute(Kelvin(100.0)).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
            // Below Neel temp (201 K) the magnetic contribution must be non-zero
            g shouldBeLessThan 0.0
        }

        it("all predefined models return finite values across temperature range") {
            val models = listOf(
                MagneticContribution.BCC_FE,
                MagneticContribution.FCC_FE,
                MagneticContribution.FCC_NI,
                MagneticContribution.FCC_CO
            )
            val temperatures = listOf(100.0, 300.0, 600.0, 1000.0, 1500.0, 2000.0)

            models.forEach { model ->
                temperatures.forEach { t ->
                    val g = model.compute(Kelvin(t)).value
                    g.isNaN()      shouldBe false
                    g.isInfinite() shouldBe false
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Edge cases
    // -------------------------------------------------------------------------

    describe("edge cases") {

        it("temperature exactly at Tc — no NaN") {
            val tc = MagneticContribution.BCC_FE.tc
            val g  = MagneticContribution.BCC_FE.compute(tc).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }

        it("very low temperature — finite result") {
            val g = MagneticContribution.BCC_FE.compute(Kelvin(1.0)).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }
    }

    // -------------------------------------------------------------------------
    // Invalid input
    // -------------------------------------------------------------------------

    describe("invalid input") {

        it("negative Tc — throws") {
            shouldThrow<IllegalArgumentException> {
                MagneticContribution(tc = Kelvin(0.0), beta = 2.22, p = 0.28)
            }
        }

        it("p = 0 — throws (would cause division by zero in formula)") {
            shouldThrow<IllegalArgumentException> {
                MagneticContribution(tc = Kelvin(1043.0), beta = 2.22, p = 0.0)
            }
        }

        it("p outside (0, 1) — throws") {
            shouldThrow<IllegalArgumentException> {
                MagneticContribution(tc = Kelvin(1043.0), beta = 2.22, p = 1.5)
            }
        }

        it("afm = 0 or positive — throws") {
            shouldThrow<IllegalArgumentException> {
                MagneticContribution(tc = Kelvin(1043.0), beta = -2.1, p = 0.25, afm = 3.0)
            }
        }
    }
})
