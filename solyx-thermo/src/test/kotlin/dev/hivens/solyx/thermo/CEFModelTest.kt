package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.PhysicalConstants.R
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.ln

class CEFModelTest : DescribeSpec({

    val t1200 = Kelvin(1200.0)

    val fe  = Species.of(Element.Fe)
    val cr  = Species.of(Element.Cr)
    val c   = Species.of(Element.C)
    val va  = Species.Vacancy

    // -------------------------------------------------------------------------
    // Analytical results
    // -------------------------------------------------------------------------

    describe("analytical results") {

        it("single endmember pure phase — G_ref = G°, G_mix = 0") {
            // (Fe)_1 (Va)_1 — only one configuration, pure endmember
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0))
                )
            )
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to 1.0)),
                Sublattice(1.0, mapOf(va to 1.0))
            )
            val g = model.compute(sublattices, t1200).value
            // G_ref = 1.0 * 1.0 * (-8000) = -8000
            // G_mix = R*T*(1*ln(1) + 1*ln(1)) = 0
            g.shouldBeWithinPercentageOf(-8000.0, 1e-6)
        }

        it("binary first sublattice — G_ref interpolates linearly") {
            // (Fe, Cr)_1 (Va)_1 — two endmembers, no mixing on second sublattice
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0)),
                    EndmemberParameter(listOf(cr, va), JoulePerMole(-6000.0))
                )
            )
            val yFe = 0.6; val yCr = 0.4
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to yFe, cr to yCr)),
                Sublattice(1.0, mapOf(va to 1.0))
            )
            val g = model.compute(sublattices, t1200).value

            // G_ref = 0.6*(-8000) + 0.4*(-6000) = -7200
            // G_mix = R*T*(0.6*ln(0.6) + 0.4*ln(0.4)) + R*T*(1*ln(1))
            val gRef = yFe * (-8000.0) + yCr * (-6000.0)
            val gMix = R * 1200.0 * (yFe * ln(yFe) + yCr * ln(yCr))
            g.shouldBeWithinPercentageOf(gRef + gMix, 1e-4)
        }

        it("interstitial solution — (Fe,Cr)_1 (C,Va)_1 known value") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, c),  JoulePerMole(-10000.0)),
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0)),
                    EndmemberParameter(listOf(cr, c),  JoulePerMole(-12000.0)),
                    EndmemberParameter(listOf(cr, va), JoulePerMole(-9000.0))
                )
            )
            val yFe = 0.7; val yCr = 0.3
            val yC  = 0.05; val yVa = 0.95
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to yFe, cr to yCr)),
                Sublattice(1.0, mapOf(c to yC, va to yVa))
            )
            val g = model.compute(sublattices, t1200).value

            val gRef = yFe * yC  * (-10000.0) +
                       yFe * yVa * (-8000.0)  +
                       yCr * yC  * (-12000.0) +
                       yCr * yVa * (-9000.0)
            val gMix1 = R * 1200.0 * (yFe * ln(yFe) + yCr * ln(yCr))
            val gMix2 = R * 1200.0 * (yC  * ln(yC)  + yVa * ln(yVa))
            g.shouldBeWithinPercentageOf(gRef + gMix1 + gMix2, 1e-3)
        }

        it("stoichiometry scaling — s=2 sublattice contributes 2x entropy") {
            // (Fe)_1 (C,Va)_2
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, c),  JoulePerMole(-10000.0)),
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0))
                )
            )
            val yC = 0.1; val yVa = 0.9
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to 1.0)),
                Sublattice(2.0, mapOf(c to yC, va to yVa))
            )
            val g = model.compute(sublattices, t1200).value

            val gRef = 1.0 * yC  * (-10000.0) + 1.0 * yVa * (-8000.0)
            val gMix = R * 1200.0 * 2.0 * (yC * ln(yC) + yVa * ln(yVa))
            g.shouldBeWithinPercentageOf(gRef + gMix, 1e-3)
        }

        it("pure endmember on both sublattices — no ideal mixing entropy") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0))
                )
            )
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to 1.0)),
                Sublattice(1.0, mapOf(va to 1.0))
            )
            // No mixing → G = G_ref exactly
            model.compute(sublattices, t1200).value
                .shouldBeWithinPercentageOf(-8000.0, 1e-6)
        }
    }

    // -------------------------------------------------------------------------
    // Physical invariants
    // -------------------------------------------------------------------------

    describe("physical invariants") {

        it("G is finite across full composition range") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, va), JoulePerMole(-8000.0)),
                    EndmemberParameter(listOf(cr, va), JoulePerMole(-6000.0))
                )
            )
            (1..99).forEach { i ->
                val y = i / 100.0
                val sublattices = listOf(
                    Sublattice(1.0, mapOf(fe to y, cr to 1.0 - y)),
                    Sublattice(1.0, mapOf(va to 1.0))
                )
                val g = model.compute(sublattices, t1200).value
                g.isNaN()      shouldBe false
                g.isInfinite() shouldBe false
            }
        }
    }

    // -------------------------------------------------------------------------
    // Guard conditions — empty config and OOB
    // -------------------------------------------------------------------------

    describe("guard conditions") {

        it("empty endmember configuration — throws IllegalArgumentException") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(emptyList(), JoulePerMole(-8000.0))
                )
            )
            val sublattices = listOf(Sublattice(1.0, mapOf(fe to 1.0)))
            shouldThrow<IllegalArgumentException> {
                model.compute(sublattices, t1200)
            }
        }

        it("endmember configuration exceeds sublattice count — throws IllegalArgumentException") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, va, c), JoulePerMole(-8000.0))
                )
            )
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to 1.0)),
                Sublattice(1.0, mapOf(va to 1.0))
                // missing 3rd sublattice
            )
            shouldThrow<IllegalArgumentException> {
                model.compute(sublattices, t1200)
            }
        }

        it("Sublattice fractions not summing to 1.0 — throws") {
            shouldThrow<IllegalArgumentException> {
                Sublattice(1.0, mapOf(fe to 0.4, cr to 0.4)) // sum = 0.8
            }
        }

        it("Sublattice negative site fraction — throws") {
            shouldThrow<IllegalArgumentException> {
                Sublattice(1.0, mapOf(fe to 1.1, cr to -0.1))
            }
        }

        it("Sublattice stoichiometry zero — throws") {
            shouldThrow<IllegalArgumentException> {
                Sublattice(0.0, mapOf(fe to 1.0))
            }
        }
    }

    // -------------------------------------------------------------------------
    // Edge cases
    // -------------------------------------------------------------------------

    describe("edge cases") {

        it("single sublattice single species — G_mix = 0") {
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe), JoulePerMole(-5000.0))
                )
            )
            val sublattices = listOf(Sublattice(1.0, mapOf(fe to 1.0)))
            model.compute(sublattices, t1200).value
                .shouldBeWithinPercentageOf(-5000.0, 1e-6)
        }

        it("three-sublattice model — computes without error") {
            val ni = Species.of(Element.Ni)
            val model = CEFModel(
                endmembers = listOf(
                    EndmemberParameter(listOf(fe, cr, va), JoulePerMole(-9000.0)),
                    EndmemberParameter(listOf(fe, ni, va), JoulePerMole(-7000.0)),
                    EndmemberParameter(listOf(cr, cr, va), JoulePerMole(-11000.0)),
                    EndmemberParameter(listOf(cr, ni, va), JoulePerMole(-8000.0))
                )
            )
            val sublattices = listOf(
                Sublattice(1.0, mapOf(fe to 0.6, cr to 0.4)),
                Sublattice(1.0, mapOf(cr to 0.5, ni to 0.5)),
                Sublattice(1.0, mapOf(va to 1.0))
            )
            val g = model.compute(sublattices, t1200).value
            g.isNaN()      shouldBe false
            g.isInfinite() shouldBe false
        }
    }
})
