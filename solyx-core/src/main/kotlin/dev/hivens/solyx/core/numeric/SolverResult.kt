package dev.hivens.solyx.core.numeric

import dev.hivens.solyx.core.error.Interval
import kotlin.math.abs
import kotlin.math.sign

/**
 * Result of a root-finding operation.
 */
sealed class SolverResult {
    /** Root found within tolerance. */
    data class Success(
        val root: Double,
        val iterations: Int,
        val residual: Interval
    ) : SolverResult()

    /** Solver did not converge within [maxIterations]. */
    data class DidNotConverge(
        val bestEstimate: Double,
        val iterations: Int,
        val residual: Double
    ) : SolverResult()

    /** Precondition violated — invalid input. */
    data class InvalidInput(val reason: String) : SolverResult()
}

// -----------------------------------------------------------------------------
// Newton-Raphson solver
// -----------------------------------------------------------------------------

/**
 * Newton-Raphson root finder.
 *
 * Fast quadratic convergence near the root, but requires a good initial guess
 * and fails when the derivative is near zero. Use [brent] when robustness
 * matters more than speed.
 *
 * f'(x) is computed numerically via central difference if not provided.
 *
 * ```kotlin
 * val result = newton(
 *     f  = { x -> x * x - 4.0 },
 *     x0 = 1.0
 * )
 * when (result) {
 *     is SolverResult.Success -> println(result.root)  // ≈ 2.0
 *     else -> error("Did not converge")
 * }
 * ```
 *
 * @param f function whose root to find
 * @param df derivative of f — auto-computed if null
 * @param x0 initial guess
 * @param tolerance convergence criterion on |f(x)|
 * @param maxIterations maximum allowed iterations
 */
fun newton(
    f: (Double) -> Double,
    df: ((Double) -> Double)? = null,
    x0: Double,
    tolerance: Double = 1e-10,
    maxIterations: Int = 100
): SolverResult {
    val derivative = df ?: { x -> derivative(f, x) }
    var x = x0

    repeat(maxIterations) { iteration ->
        val fx = f(x)

        if (abs(fx) < tolerance) {
            return SolverResult.Success(
                root = x,
                iterations = iteration,
                residual = Interval(fx, tolerance)
            )
        }

        val dfx = derivative(x)
        if (abs(dfx) < 1e-14) {
            return SolverResult.DidNotConverge(x, iteration, abs(fx))
        }

        x -= fx / dfx
    }

    return SolverResult.DidNotConverge(x, maxIterations, abs(f(x)))
}

// -----------------------------------------------------------------------------
// Brent solver
// -----------------------------------------------------------------------------

/**
 * Brent's method — robust hybrid root finder.
 *
 * Combines bisection, secant, and inverse quadratic interpolation.
 * Guaranteed to converge if f(a) and f(b) have opposite signs.
 * Preferred over Newton when:
 * - derivative is unavailable or expensive
 * - initial guess quality is uncertain
 * - robustness is required over speed
 *
 * Used as the primary solver in the Gibbs minimizer for phase boundary detection.
 *
 * ```kotlin
 * val result = brent(
 *     f = { x -> x * x - 4.0 },
 *     a = 0.0,
 *     b = 3.0
 * )
 * ```
 *
 * @param f function whose root to find
 * @param a left bracket — f(a) and f(b) must have opposite signs
 * @param b right bracket
 * @param tolerance convergence criterion
 * @param maxIterations maximum allowed iterations
 */
fun brent(
    f: (Double) -> Double,
    a: Double,
    b: Double,
    tolerance: Double = 1e-10,
    maxIterations: Int = 100
): SolverResult {
    var fa = f(a)
    var fb = f(b)

    if (fa * fb > 0.0) {
        return SolverResult.InvalidInput(
            "f(a) and f(b) must have opposite signs: f($a)=$fa, f($b)=$fb"
        )
    }

    if (abs(fa) < tolerance) return SolverResult.Success(a, 0, Interval(fa, tolerance))
    if (abs(fb) < tolerance) return SolverResult.Success(b, 0, Interval(fb, tolerance))

    var lo = a; var hi = b
    var c = lo; var fc = fa
    var d = hi - lo; var e = d

    repeat(maxIterations) { iteration ->
        if (abs(fc) < abs(fb)) {
            lo = hi; fa = fb
            hi = c; fb = fc
            c = lo; fc = fa
        }

        val tol = 2.0 * 2.220446049250313e-16 * abs(hi) + tolerance / 2.0
        val mid = (c - hi) / 2.0

        if (abs(mid) <= tol || abs(fb) < tolerance) {
            return SolverResult.Success(
                root = hi,
                iterations = iteration,
                residual = Interval(fb, tolerance)
            )
        }

        if (abs(e) >= tol && abs(fa) > abs(fb)) {
            val s = fb / fa
            var p: Double
            var q: Double

            if (lo == c) {
                // Secant
                p = 2.0 * mid * s
                q = 1.0 - s
            } else {
                // Inverse quadratic interpolation
                val r = fb / fc
                q = fa / fc
                p = s * (2.0 * mid * q * (q - r) - (hi - lo) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            }

            if (p > 0.0) q = -q else p = -p

            if (2.0 * p < minOf(3.0 * mid * q - abs(tol * q), abs(e * q))) {
                e = d; d = p / q
            } else {
                d = mid; e = mid
            }
        } else {
            d = mid; e = mid
        }

        lo = hi; fa = fb
        hi += if (abs(d) > tol) d else tol * sign(mid)
        fb = f(hi)

        if (fb * fc > 0.0) { c = lo; fc = fa; e = hi - lo; d = e }
    }

    return SolverResult.DidNotConverge(hi, maxIterations, abs(fb))
}
