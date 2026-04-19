package dev.hivens.solyx.core.numeric

import dev.hivens.solyx.core.error.Interval
import kotlin.math.abs

/**
 * Gauss-Legendre quadrature — high-accuracy numerical integration.
 *
 * Exact for polynomials of degree up to 2n-1 using only n function evaluations.
 * Far more efficient than Simpson or trapezoidal rules for smooth functions.
 *
 * Used in Solyx for integrating thermodynamic quantities such as:
 * - heat capacity over a temperature range: ∫Cp(T)dT
 * - enthalpy increments: H(T2) - H(T1) = ∫Cp dT
 * - entropy increments: S(T2) - S(T1) = ∫(Cp/T) dT
 *
 * ```kotlin
 * // ∫x² dx from 0 to 3 = 9.0
 * val result = integrate({ x -> x * x }, 0.0, 3.0)
 *
 * // with error estimate
 * val result = integrateWithError({ x -> x * x }, 0.0, 3.0)
 * if (result.isReliable) println(result.midpoint)
 * ```
 */

/**
 * Integrate [f] over [[from], [to]] using Gauss-Legendre quadrature.
 *
 * @param f integrand
 * @param from lower bound
 * @param to upper bound
 * @param points number of quadrature points — higher means more accurate
 */
fun integrate(
    f: (Double) -> Double,
    from: Double,
    to: Double,
    points: GaussPoints = GaussPoints.N5
): Double {
    val (nodes, weights) = points.nodesAndWeights
    val mid = (to + from) / 2.0
    val half = (to - from) / 2.0

    return half * kahanSum(DoubleArray(nodes.size) { i ->
        weights[i] * f(mid + half * nodes[i])
    })
}

/**
 * Integrate [f] with an error estimate via Richardson extrapolation.
 *
 * Runs integration twice — with [points] and with the next higher order —
 * and returns the difference as the uncertainty bound.
 *
 * When [points] is already [GaussPoints.N10] (the highest fixed order),
 * falls back to [integrateAdaptive] for the fine estimate to avoid
 * reporting a spurious zero error from comparing a method with itself.
 */
fun integrateWithError(
    f: (Double) -> Double,
    from: Double,
    to: Double,
    points: GaussPoints = GaussPoints.N5
): Interval {
    val coarse = integrate(f, from, to, points)
    val fine = if (points == GaussPoints.N10) {
        integrateAdaptive(f, from, to, tolerance = 1e-12)
    } else {
        integrate(f, from, to, points.higher)
    }
    val error  = abs(fine - coarse)
    return Interval(fine, error / 2.0)
}

/**
 * Adaptive integration — subdivides interval where the integrand varies rapidly.
 *
 * Uses recursive bisection when local error estimate exceeds [tolerance].
 * More robust than fixed-order quadrature for functions with sharp features,
 * such as phase transition regions in thermodynamic calculations.
 *
 * @param f integrand
 * @param from lower bound
 * @param to upper bound
 * @param tolerance absolute error tolerance
 * @param maxDepth maximum recursion depth — prevents infinite subdivision
 */
fun integrateAdaptive(
    f: (Double) -> Double,
    from: Double,
    to: Double,
    tolerance: Double = 1e-9,
    maxDepth: Int = 50
): Double = adaptiveStep(f, from, to, tolerance, maxDepth, 0)

private fun adaptiveStep(
    f: (Double) -> Double,
    a: Double,
    b: Double,
    tolerance: Double,
    maxDepth: Int,
    depth: Int
): Double {
    val whole = integrate(f, a, b, GaussPoints.N5)
    val mid   = (a + b) / 2.0
    val left  = integrate(f, a, mid, GaussPoints.N5)
    val right = integrate(f, mid, b, GaussPoints.N5)
    val error = abs(left + right - whole)

    return if (error < tolerance * 15.0 || depth >= maxDepth) {
        left + right + (left + right - whole) / 15.0
    } else {
        adaptiveStep(f, a, mid, tolerance / 2.0, maxDepth, depth + 1) +
                adaptiveStep(f, mid, b, tolerance / 2.0, maxDepth, depth + 1)
    }
}

// -----------------------------------------------------------------------------
// Quadrature nodes and weights
// -----------------------------------------------------------------------------

/**
 * Pre-computed Gauss-Legendre nodes and weights on [-1, 1].
 *
 * Higher order = more accurate but more function evaluations.
 * N5 is sufficient for smooth thermodynamic functions.
 * N10 for functions with mild curvature near phase transitions.
 */
enum class GaussPoints {
    N2, N3, N5, N7, N10;

    val nodesAndWeights: Pair<DoubleArray, DoubleArray> get() = when (this) {
        N2 -> Pair(
            doubleArrayOf(-0.5773502691896258, 0.5773502691896258),
            doubleArrayOf(1.0, 1.0)
        )
        N3 -> Pair(
            doubleArrayOf(-0.7745966692414834, 0.0, 0.7745966692414834),
            doubleArrayOf(0.5555555555555556, 0.8888888888888888, 0.5555555555555556)
        )
        N5 -> Pair(
            doubleArrayOf(
                -0.9061798459386640, -0.5384693101056831, 0.0,
                0.5384693101056831,  0.9061798459386640
            ),
            doubleArrayOf(
                0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                0.4786286704993665, 0.2369268850561891
            )
        )
        N7 -> Pair(
            doubleArrayOf(
                -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0,
                0.4058451513773972,  0.7415311855993945,  0.9491079123427585
            ),
            doubleArrayOf(
                0.1294849661688697, 0.2797053914892767, 0.3818300505051189, 0.4179591836734694,
                0.3818300505051189, 0.2797053914892767, 0.1294849661688697
            )
        )
        N10 -> Pair(
            doubleArrayOf(
                -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                -0.4333953941292472, -0.1488743389816312,  0.1488743389816312,
                0.4333953941292472,  0.6794095682990244,  0.8650633666889845,
                0.9739065285171717
            ),
            doubleArrayOf(
                0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
                0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
                0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
                0.0666713443086881
            )
        )
    }

    /** Next higher order for error estimation. */
    val higher: GaussPoints get() = when (this) {
        N2  -> N3
        N3  -> N5
        N5  -> N7
        N7  -> N10
        N10 -> N10
    }
}
