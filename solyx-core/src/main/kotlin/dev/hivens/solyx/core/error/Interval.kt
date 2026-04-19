package dev.hivens.solyx.core.error

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/**
 * Interval arithmetic — tracks accumulated numerical error through calculations.
 *
 * Every operation widens the interval to guarantee the true result lies within [lo, hi].
 * When [width] grows too large, the result is considered unreliable.
 *
 * ```kotlin
 * val a = Interval(1.0)
 * val b = Interval(2.0, tolerance = 1e-10)
 * val c = a + b  // [2.9999999999, 3.0000000001]
 *
 * if (c.isReliable) println(c.midpoint)
 * else error("Result is unreliable: uncertainty = ${c.width}")
 * ```
 */
data class Interval(val lo: Double, val hi: Double) {

    init {
        require(lo <= hi) { "Invalid interval: lo=$lo > hi=$hi" }
        require(lo.isFinite() && hi.isFinite()) { "Interval must be finite: [$lo, $hi]" }
    }

    /** Midpoint — best estimate of the true value. */
    val midpoint: Double get() = (lo + hi) / 2.0

    /** Width — measure of accumulated uncertainty. */
    val width: Double get() = hi - lo

    /** Relative uncertainty with respect to midpoint. */
    val relativeWidth: Double get() = if (midpoint != 0.0) width / abs(midpoint) else width

    /**
     * Result is considered reliable when relative uncertainty is below this threshold.
     * For engineering calculations: 1e-9 is a reasonable default.
     */
    val isReliable: Boolean get() = relativeWidth < DEFAULT_TOLERANCE

    // -------------------------------------------------------------------------
    // Arithmetic — each operation widens the interval conservatively
    // -------------------------------------------------------------------------

    operator fun plus(other: Interval) = Interval(
        lo + other.lo,
        hi + other.hi
    )

    operator fun minus(other: Interval) = Interval(
        lo - other.hi,
        hi - other.lo
    )

    operator fun times(other: Interval): Interval {
        val products = doubleArrayOf(lo * other.lo, lo * other.hi, hi * other.lo, hi * other.hi)
        return Interval(products.min(), products.max())
    }

    operator fun div(other: Interval): Interval {
        require(0.0 !in other) { "Division by interval containing zero: $other" }
        return times(Interval(1.0 / other.hi, 1.0 / other.lo))
    }

    operator fun unaryMinus() = Interval(-hi, -lo)

    operator fun contains(value: Double) = value in lo..hi

    // -------------------------------------------------------------------------
    // Scalar operations — convenience
    // -------------------------------------------------------------------------

    operator fun plus(scalar: Double)  = plus(Interval(scalar))
    operator fun minus(scalar: Double) = minus(Interval(scalar))
    operator fun times(scalar: Double) = if (scalar >= 0)
        Interval(lo * scalar, hi * scalar)
    else
        Interval(hi * scalar, lo * scalar)
    operator fun div(scalar: Double)   = times(1.0 / scalar)

    // -------------------------------------------------------------------------
    // Mathematical functions — conservative bounds
    // -------------------------------------------------------------------------

    fun sqrt(): Interval {
        require(lo >= 0.0) { "sqrt of negative interval: $this" }
        return Interval(sqrt(lo), sqrt(hi))
    }

    fun ln(): Interval {
        require(lo > 0.0) { "ln of non-positive interval: $this" }
        return Interval(kotlin.math.ln(lo), kotlin.math.ln(hi))
    }

    fun pow(n: Int): Interval {
        if (n == 0) return Interval(1.0)
        if (n == 1) return this
        val p = pow(n / 2)
        return if (n % 2 == 0) p * p else p * p * this
    }

    // -------------------------------------------------------------------------
    // Comparison
    // -------------------------------------------------------------------------

    /** True if this interval is strictly less than [other]. */
    fun strictlyLessThan(other: Interval) = hi < other.lo

    /** True if this interval is strictly greater than [other]. */
    fun strictlyGreaterThan(other: Interval) = lo > other.hi

    /** True if intervals overlap. */
    fun overlaps(other: Interval) = lo <= other.hi && other.lo <= hi

    // -------------------------------------------------------------------------
    // Utilities
    // -------------------------------------------------------------------------

    /** Intersect two intervals — useful for constraint propagation. */
    fun intersect(other: Interval): Interval? {
        val newLo = max(lo, other.lo)
        val newHi = min(hi, other.hi)
        return if (newLo <= newHi) Interval(newLo, newHi) else null
    }

    /** Extend interval to include [value]. */
    fun extend(value: Double) = Interval(min(lo, value), max(hi, value))

    override fun toString() = "[${lo}, ${hi}] (±${width / 2})"

    companion object {
        const val DEFAULT_TOLERANCE = 1e-9

        /** Point interval — no uncertainty. */
        operator fun invoke(value: Double) = Interval(value, value)

        /**
         * Interval with explicit tolerance around a central value.
         *
         * ```kotlin
         * val x = Interval(1200.0, tolerance = 1e-10)
         * ```
         */
        operator fun invoke(value: Double, tolerance: Double) =
            Interval(value - tolerance, value + tolerance)

        /** Hull of a collection — smallest interval containing all values. */
        fun hull(values: Collection<Double>): Interval {
            require(values.isNotEmpty()) { "Cannot compute hull of empty collection" }
            return Interval(values.min(), values.max())
        }
    }
}

/** Wrap a plain Double as a point interval. */
val Double.interval get() = Interval(this)

/** Wrap with explicit tolerance. */
fun Double.interval(tolerance: Double) = Interval(this, tolerance)
