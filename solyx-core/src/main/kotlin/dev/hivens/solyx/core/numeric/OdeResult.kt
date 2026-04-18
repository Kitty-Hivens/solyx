package dev.hivens.solyx.core.numeric

import dev.hivens.solyx.core.error.Interval
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

/**
 * Result of an ODE integration.
 */
sealed class OdeResult {
    /**
     * Integration completed successfully.
     *
     * @param t time points
     * @param y state values at each time point
     * @param steps total number of steps taken
     */
    data class Success(
        val t: DoubleArray,
        val y: DoubleArray,
        val steps: Int
    ) : OdeResult() {
        /** Final state value. */
        val finalValue: Double get() = y.last()

        /** Final time. */
        val finalTime: Double get() = t.last()
        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (other !is Success) return false
            return t.contentEquals(other.t) &&
                    y.contentEquals(other.y) &&
                    steps == other.steps
        }

        override fun hashCode(): Int {
            var result = t.contentHashCode()
            result = 31 * result + y.contentHashCode()
            result = 31 * result + steps
            return result
        }
    }

    /** Integration failed — step size collapsed below minimum. */
    data class StepSizeTooSmall(
        val t: Double,
        val y: Double,
        val minStep: Double
    ) : OdeResult()

    /** Maximum number of steps reached before reaching [tEnd]. */
    data class MaxStepsReached(
        val t: Double,
        val y: Double,
        val steps: Int
    ) : OdeResult()
}

// -----------------------------------------------------------------------------
// RK4 — fixed step
// -----------------------------------------------------------------------------

/**
 * Classic 4th-order Runge-Kutta solver — fixed step size.
 *
 * Suitable when step size is known in advance and the solution is smooth.
 * For thermodynamic cooling/heating simulations with uniform time steps.
 *
 * ```kotlin
 * // dT/dt = -k * (T - T_ambient)
 * val result = rk4(
 *     f  = { t, y -> -0.1 * (y - 300.0) },
 *     y0 = 1500.0,
 *     t0 = 0.0,
 *     tEnd = 50.0,
 *     steps = 500
 * )
 * when (result) {
 *     is OdeResult.Success -> println(result.finalValue)
 *     else -> error("Integration failed")
 * }
 * ```
 *
 * @param f right-hand side: dy/dt = f(t, y)
 * @param y0 initial state
 * @param t0 initial time
 * @param tEnd final time
 * @param steps number of steps
 */
fun rk4(
    f: (t: Double, y: Double) -> Double,
    y0: Double,
    t0: Double,
    tEnd: Double,
    steps: Int
): OdeResult {
    require(steps > 0) { "steps must be positive: $steps" }
    require(tEnd > t0)  { "tEnd must be greater than t0: $tEnd <= $t0" }

    val dt = (tEnd - t0) / steps
    val tArr = DoubleArray(steps + 1)
    val yArr = DoubleArray(steps + 1)

    tArr[0] = t0
    yArr[0] = y0

    for (i in 0 until steps) {
        val t = tArr[i]
        val y = yArr[i]

        val k1 = f(t,            y)
        val k2 = f(t + dt / 2.0, y + dt * k1 / 2.0)
        val k3 = f(t + dt / 2.0, y + dt * k2 / 2.0)
        val k4 = f(t + dt,       y + dt * k3)

        tArr[i + 1] = t + dt
        yArr[i + 1] = y + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    }

    return OdeResult.Success(tArr, yArr, steps)
}

// -----------------------------------------------------------------------------
// RK45 — adaptive step (Dormand-Prince)
// -----------------------------------------------------------------------------

/**
 * Dormand-Prince RK45 — adaptive step ODE solver.
 *
 * Computes both a 4th and 5th order solution simultaneously using 6 function
 * evaluations. The difference between them estimates local error, which drives
 * automatic step size control:
 * - error too large → shrink step, retry
 * - error small     → accept step, grow step
 *
 * Preferred over [rk4] when:
 * - solution has regions of rapid change (phase transitions)
 * - optimal step size is not known in advance
 * - uniform accuracy across the integration interval is required
 *
 * This is the same algorithm as MATLAB ode45.
 *
 * ```kotlin
 * val result = rk45(
 *     f    = { t, y -> -0.1 * (y - 300.0) },
 *     y0   = 1500.0,
 *     t0   = 0.0,
 *     tEnd = 50.0
 * )
 * ```
 *
 * @param f right-hand side: dy/dt = f(t, y)
 * @param y0 initial state
 * @param t0 initial time
 * @param tEnd final time
 * @param tolerance absolute + relative error tolerance
 * @param initialStep initial step size — auto-selected if null
 * @param minStep minimum allowed step — fails with [OdeResult.StepSizeTooSmall] if reached
 * @param maxSteps maximum number of steps
 */
fun rk45(
    f: (t: Double, y: Double) -> Double,
    y0: Double,
    t0: Double,
    tEnd: Double,
    tolerance: Double = 1e-9,
    initialStep: Double? = null,
    minStep: Double = 1e-12,
    maxSteps: Int = 100_000
): OdeResult {
    require(tEnd > t0) { "tEnd must be greater than t0: $tEnd <= $t0" }

    // Dormand-Prince coefficients
    val c2 = 1.0/5.0;  val c3 = 3.0/10.0; val c4 = 4.0/5.0
    val c5 = 8.0/9.0

    val a21 = 1.0/5.0
    val a31 = 3.0/40.0;       val a32 = 9.0/40.0
    val a41 = 44.0/45.0;      val a42 = -56.0/15.0;     val a43 = 32.0/9.0
    val a51 = 19372.0/6561.0; val a52 = -25360.0/2187.0; val a53 = 64448.0/6561.0; val a54 = -212.0/729.0
    val a61 = 9017.0/3168.0;  val a62 = -355.0/33.0;     val a63 = 46732.0/5247.0; val a64 = 49.0/176.0;   val a65 = -5103.0/18656.0

    // 5th order weights
    val b1 = 35.0/384.0;  val b3 = 500.0/1113.0; val b4 = 125.0/192.0
    val b5 = -2187.0/6784.0; val b6 = 11.0/84.0

    // Error coefficients (difference between 4th and 5th order)
    val e1 =  71.0/57600.0;  val e3 = -71.0/16695.0; val e4 =  71.0/1920.0
    val e5 = -17253.0/339200.0; val e6 = 22.0/525.0; val e7 = -1.0/40.0

    val tList = mutableListOf(t0)
    val yList = mutableListOf(y0)

    var t = t0
    var y = y0
    var h = initialStep ?: ((tEnd - t0) / 100.0)
    var steps = 0

    while (t < tEnd) {
        if (steps >= maxSteps) return OdeResult.MaxStepsReached(t, y, steps)
        if (h < minStep)       return OdeResult.StepSizeTooSmall(t, y, minStep)

        h = min(h, tEnd - t)

        val k1 = f(t,           y)
        val k2 = f(t + c2 * h, y + h * a21 * k1)
        val k3 = f(t + c3 * h, y + h * (a31 * k1 + a32 * k2))
        val k4 = f(t + c4 * h, y + h * (a41 * k1 + a42 * k2 + a43 * k3))
        val k5 = f(t + c5 * h, y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
        val k6 = f(t + h,      y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))

        val yNext = y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)

        val k7 = f(t + h, yNext)

        // Local error estimate
        val err = abs(h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7))
        val scale = tolerance * (1.0 + max(abs(y), abs(yNext)))

        if (err <= scale) {
            // Accept step
            t += h
            y = yNext
            tList.add(t)
            yList.add(y)
            steps++

            // Grow step
            val factor = min(5.0, 0.9 * (scale / err).pow(0.2))
            h *= factor
        } else {
            // Reject step — shrink
            val factor = max(0.1, 0.9 * (scale / err).pow(0.25))
            h *= factor
        }
    }

    return OdeResult.Success(tList.toDoubleArray(), yList.toDoubleArray(), steps)
}
