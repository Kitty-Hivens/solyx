package dev.hivens.solyx.core.numeric

import dev.hivens.solyx.core.error.Interval
import kotlin.math.abs
import kotlin.math.cbrt

/**
 * Numerical differentiation via central difference.
 *
 * Central difference is used over forward difference because it has
 * O(h²) error instead of O(h), at the cost of one extra function evaluation.
 *
 * Step size h is chosen automatically based on machine epsilon to balance
 * truncation error (smaller h is better) against roundoff error (larger h is better).
 * Optimal h ≈ ε^(1/3) for central difference.
 *
 * ```kotlin
 * val df = derivative({ x -> x * x }, 3.0)  // ≈ 6.0
 *
 * val (value, error) = derivativeWithError({ x -> x * x }, 3.0)
 * ```
 */

private const val MACHINE_EPSILON = 2.220446049250313e-16

/**
 * Optimal step size for central difference: h = ε^(1/3) * max(1, |x|)
 */
private fun optimalStep(x: Double): Double =
    cbrt(MACHINE_EPSILON) * maxOf(1.0, abs(x))

/**
 * First derivative via central difference.
 *
 * f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
 *
 * @param f function to differentiate
 * @param x point at which to evaluate the derivative
 * @param h step size — auto-selected if null
 */
fun derivative(f: (Double) -> Double, x: Double, h: Double? = null): Double {
    val step = h ?: optimalStep(x)
    return (f(x + step) - f(x - step)) / (2.0 * step)
}

/**
 * First derivative with uncertainty estimate.
 *
 * Returns the derivative as an [Interval] where the width reflects
 * the estimated numerical error.
 *
 * @param f function to differentiate
 * @param x point at which to evaluate the derivative
 */
fun derivativeWithError(f: (Double) -> Double, x: Double): Interval {
    val h = optimalStep(x)
    val central = (f(x + h) - f(x - h)) / (2.0 * h)

    // Richardson extrapolation for error estimate
    val h2 = h * 2.0
    val coarse = (f(x + h2) - f(x - h2)) / (2.0 * h2)
    val error = abs(central - coarse) / 3.0

    return Interval(central, error)
}

/**
 * Second derivative via central difference.
 *
 * f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
 */
fun secondDerivative(f: (Double) -> Double, x: Double, h: Double? = null): Double {
    val step = h ?: optimalStep(x)
    return (f(x + step) - 2.0 * f(x) + f(x - step)) / (step * step)
}

/**
 * Partial derivative with respect to argument [index] of a multivariate function.
 *
 * ```kotlin
 * // f(x, y) = x² + y²
 * val f = { args: DoubleArray -> args[0].pow(2) + args[1].pow(2) }
 * val dfdx = partialDerivative(f, doubleArrayOf(3.0, 4.0), index = 0)  // ≈ 6.0
 * val dfdy = partialDerivative(f, doubleArrayOf(3.0, 4.0), index = 1)  // ≈ 8.0
 * ```
 */
fun partialDerivative(
    f: (DoubleArray) -> Double,
    point: DoubleArray,
    index: Int,
    h: Double? = null
): Double {
    val step = h ?: optimalStep(point[index])
    val forward = point.copyOf().also { it[index] += step }
    val backward = point.copyOf().also { it[index] -= step }
    return (f(forward) - f(backward)) / (2.0 * step)
}

/**
 * Gradient of a multivariate function at [point].
 *
 * Returns partial derivatives for all arguments.
 */
fun gradient(f: (DoubleArray) -> Double, point: DoubleArray): DoubleArray =
    DoubleArray(point.size) { i -> partialDerivative(f, point, i) }
