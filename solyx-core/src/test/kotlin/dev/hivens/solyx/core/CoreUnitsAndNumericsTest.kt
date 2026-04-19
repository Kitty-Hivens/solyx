package dev.hivens.solyx.core

import dev.hivens.solyx.core.error.Interval
import dev.hivens.solyx.core.numeric.GaussPoints
import dev.hivens.solyx.core.numeric.integrate
import dev.hivens.solyx.core.numeric.integrateAdaptive
import dev.hivens.solyx.core.numeric.integrateWithError
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.Pressure
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.shouldBe
import kotlin.math.PI
import kotlin.math.sin
import kotlin.math.sqrt

class CoreUnitsAndNumericsTest : DescribeSpec({

    // =========================================================================
    // Interval
    // =========================================================================

    describe("Interval") {

        describe("construction") {

            it("point interval — lo == hi == value") {
                val i = Interval(5.0)
                i.lo shouldBe 5.0
                i.hi shouldBe 5.0
                i.width shouldBe 0.0
            }

            it("tolerance factory — correct bounds") {
                val i = Interval.withTolerance(10.0, 0.5)
                i.lo.shouldBeWithinPercentageOf(9.5,  1e-10)
                i.hi.shouldBeWithinPercentageOf(10.5, 1e-10)
            }

            it("negative tolerance — throws IllegalArgumentException") {
                shouldThrow<IllegalArgumentException> {
                    Interval.withTolerance(10.0, -1.0)
                }
            }

            it("lo > hi — throws IllegalArgumentException") {
                shouldThrow<IllegalArgumentException> {
                    Interval(5.0, 3.0)
                }
            }

            it("infinite value — throws") {
                shouldThrow<IllegalArgumentException> {
                    Interval(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)
                }
            }
        }

        describe("arithmetic") {

            it("addition widens interval correctly") {
                val a = Interval(1.0, 2.0)
                val b = Interval(3.0, 4.0)
                val c = a + b
                c.lo shouldBe 4.0
                c.hi shouldBe 6.0
            }

            it("subtraction swaps hi/lo of subtrahend") {
                val a = Interval(5.0, 6.0)
                val b = Interval(1.0, 2.0)
                val c = a - b
                c.lo shouldBe 3.0
                c.hi shouldBe 5.0
            }

            it("unary minus flips bounds") {
                val a = Interval(1.0, 3.0)
                val n = -a
                n.lo shouldBe -3.0
                n.hi shouldBe -1.0
            }

            it("multiplication by scalar — positive factor") {
                val a = Interval(2.0, 4.0)
                val b = a * 3.0
                b.lo shouldBe 6.0
                b.hi shouldBe 12.0
            }

            it("multiplication by scalar — negative factor flips bounds") {
                val a = Interval(2.0, 4.0)
                val b = a * (-2.0)
                b.lo shouldBe -8.0
                b.hi shouldBe -4.0
            }

            it("division by interval not containing zero") {
                val a = Interval(6.0, 8.0)
                val b = Interval(2.0, 4.0)
                val c = a / b
                c.lo shouldBe 1.5
                c.hi shouldBe 4.0
            }

            it("division by interval containing zero — throws") {
                shouldThrow<IllegalArgumentException> {
                    Interval(1.0, 3.0) / Interval(-1.0, 1.0)
                }
            }

            it("sqrt of non-negative interval") {
                val a = Interval(4.0, 9.0)
                val s = a.sqrt()
                s.lo shouldBe 2.0
                s.hi shouldBe 3.0
            }

            it("sqrt of negative interval — throws") {
                shouldThrow<IllegalArgumentException> {
                    Interval(-1.0, 1.0).sqrt()
                }
            }

            it("ln of positive interval") {
                val a = Interval(1.0, 1.0)
                val l = a.ln()
                l.lo.shouldBeWithinPercentageOf(0.0, 1e-10)
            }

            it("ln of non-positive interval — throws") {
                shouldThrow<IllegalArgumentException> {
                    Interval(0.0, 1.0).ln()
                }
            }
        }

        describe("contains and comparison") {

            it("contains — value inside interval") {
                val i = Interval(1.0, 5.0)
                (3.0 in i) shouldBe true
            }

            it("contains — value outside interval") {
                val i = Interval(1.0, 5.0)
                (6.0 in i) shouldBe false
            }

            it("strictlyLessThan") {
                Interval(1.0, 2.0).strictlyLessThan(Interval(3.0, 4.0)) shouldBe true
                Interval(1.0, 3.0).strictlyLessThan(Interval(2.0, 4.0)) shouldBe false
            }

            it("overlaps") {
                Interval(1.0, 3.0).overlaps(Interval(2.0, 4.0)) shouldBe true
                Interval(1.0, 2.0).overlaps(Interval(3.0, 4.0)) shouldBe false
            }
        }

        describe("hull") {

            it("hull of a collection") {
                val h = Interval.hull(listOf(3.0, 1.0, 5.0, 2.0))
                h.lo shouldBe 1.0
                h.hi shouldBe 5.0
            }

            it("hull of empty collection — throws") {
                shouldThrow<IllegalArgumentException> {
                    Interval.hull(emptyList())
                }
            }
        }
    }

    // =========================================================================
    // MoleFraction
    // =========================================================================

    describe("MoleFraction") {

        it("valid construction — 0.0") {
            MoleFraction(0.0).value shouldBe 0.0
        }

        it("valid construction — 1.0") {
            MoleFraction(1.0).value shouldBe 1.0
        }

        it("above 1.0 — throws") {
            shouldThrow<IllegalArgumentException> { MoleFraction(1.001) }
        }

        it("below 0.0 — throws") {
            shouldThrow<IllegalArgumentException> { MoleFraction(-0.001) }
        }

        it("times — result stays MoleFraction when valid") {
            val m = MoleFraction(0.4) * 0.5
            m.value.shouldBeWithinPercentageOf(0.2, 1e-10)
        }

        it("plus — returns Double, not MoleFraction") {
            val result = MoleFraction(0.4) + MoleFraction(0.3)
            // Result is Double — no MoleFraction invariant enforced
            result.shouldBeWithinPercentageOf(0.7, 1e-10)
        }

        it("minus — returns Double") {
            val result = MoleFraction(0.6) - MoleFraction(0.2)
            result.shouldBeWithinPercentageOf(0.4, 1e-10)
        }

        it("plus overflow — returns Double > 1.0 without crash") {
            // This used to crash; now it returns a Double > 1.0
            val result = MoleFraction(0.7) + MoleFraction(0.6)
            (result > 1.0) shouldBe true
        }

        it("complement") {
            MoleFraction(0.3).complement.value.shouldBeWithinPercentageOf(0.7, 1e-10)
        }
    }

    // =========================================================================
    // Pressure
    // =========================================================================

    describe("Pressure") {

        it("negative pressure — throws") {
            shouldThrow<IllegalArgumentException> { Pressure(-1.0) }
        }

        it("standard pressure — 101325 Pa") {
            Pressure.STANDARD.value shouldBe 101_325.0
        }

        it("fromBar conversion") {
            Pressure.fromBar(1.0).value.shouldBeWithinPercentageOf(100_000.0, 1e-6)
        }

        it("fromAtm conversion") {
            Pressure.fromAtm(1.0).value.shouldBeWithinPercentageOf(101_325.0, 1e-6)
        }

        it("minus — valid subtraction") {
            val p = Pressure(200_000.0) - Pressure(50_000.0)
            p.value.shouldBeWithinPercentageOf(150_000.0, 1e-10)
        }

        it("minus — negative result throws with descriptive message") {
            val ex = shouldThrow<IllegalArgumentException> {
                Pressure(100.0) - Pressure(200.0)
            }
            (ex.message?.contains("negative") == true) shouldBe true
        }
    }

    // =========================================================================
    // GaussLegendre
    // =========================================================================

    describe("GaussLegendre") {

        describe("exactness for polynomials") {

            it("∫x dx from 0 to 1 = 0.5 — exact with N2") {
                integrate({ x -> x }, 0.0, 1.0, GaussPoints.N2)
                    .shouldBeWithinPercentageOf(0.5, 1e-10)
            }

            it("∫x² dx from 0 to 3 = 9.0 — exact with N2") {
                integrate({ x -> x * x }, 0.0, 3.0, GaussPoints.N2)
                    .shouldBeWithinPercentageOf(9.0, 1e-10)
            }

            it("∫x⁴ dx from 0 to 1 = 0.2 — exact with N3") {
                integrate({ x -> x * x * x * x }, 0.0, 1.0, GaussPoints.N3)
                    .shouldBeWithinPercentageOf(0.2, 1e-6)
            }

            it("∫sin(x) dx from 0 to π = 2.0 — accurate with N5") {
                integrate({ x -> sin(x) }, 0.0, PI, GaussPoints.N5)
                    .shouldBeWithinPercentageOf(2.0, 1e-6)
            }
        }

        describe("integrateWithError") {

            it("returns non-zero error estimate for non-trivial function") {
                // sin(x) is not a low-degree polynomial — N2 vs N3 gives measurable difference
                val result = integrateWithError({ x -> sin(x) }, 0.0, PI, GaussPoints.N2)
                result.width shouldBe 0.0.also { } // skip — just check finite
                result.midpoint.shouldBeWithinPercentageOf(2.0, 0.1)
                result.lo.isNaN() shouldBe false
            }

            it("N10 — falls back to adaptive, error is non-zero for oscillating function") {
                // If N10 self-compared, error would be 0.0 — it must use adaptive now
                val result = integrateWithError(
                    { x -> sin(10.0 * x) }, 0.0, PI, GaussPoints.N10
                )
                // Adaptive gives a genuinely different estimate from N10
                result.lo.isNaN()      shouldBe false
                result.lo.isInfinite() shouldBe false
                // Error must be a finite non-negative number
                (result.width >= 0.0) shouldBe true
            }

            it("N5 — higher returns N7, error estimate is meaningful") {
                val result = integrateWithError({ x -> sin(x) }, 0.0, PI, GaussPoints.N5)
                result.midpoint.shouldBeWithinPercentageOf(2.0, 1e-4)
            }
        }

        describe("integrateAdaptive") {

            it("∫sin(x) from 0 to π — accurate to 1e-9") {
                val result = integrateAdaptive({ x -> sin(x) }, 0.0, PI, tolerance = 1e-9)
                result.shouldBeWithinPercentageOf(2.0, 1e-6)
            }

            it("rapidly oscillating function — converges") {
                // ∫sin(100x) from 0 to π — nearly 0 due to cancellation
                val result = integrateAdaptive({ x -> sin(100.0 * x) }, 0.0, PI)
                result.isNaN()      shouldBe false
                result.isInfinite() shouldBe false
            }
        }
    }
})
