package dev.hivens.solyx.core.numeric

/**
 * Kahan compensated summation — O(ε) error regardless of n.
 *
 * Standard summation accumulates floating-point error proportional to n.
 * Kahan's algorithm tracks the lost bits and compensates on each step,
 * keeping total error bounded by machine epsilon regardless of input size.
 *
 * Critical for thermodynamic calculations where millions of iterations
 * would otherwise produce unacceptable error accumulation.
 *
 * ```kotlin
 * val sum = kahanSum(doubleArrayOf(1.0, 1e-10, 1e-10, 1e-10))
 *
 * val accumulator = KahanAccumulator()
 * accumulator += 1.0
 * accumulator += 1e-10
 * val result = accumulator.sum
 * ```
 */

/**
 * Compute the Kahan sum of an array.
 *
 * @param values input values
 * @return compensated sum with O(ε) error
 */
fun kahanSum(values: DoubleArray): Double {
    var sum = 0.0
    var compensation = 0.0

    for (v in values) {
        val y = v - compensation
        val t = sum + y
        compensation = (t - sum) - y
        sum = t
    }

    return sum
}

/**
 * Compute the Kahan sum of a collection.
 */
fun kahanSum(values: Iterable<Double>): Double {
    var sum = 0.0
    var compensation = 0.0

    for (v in values) {
        val y = v - compensation
        val t = sum + y
        compensation = (t - sum) - y
        sum = t
    }

    return sum
}

/**
 * Compute the Kahan sum of a sequence with a transform.
 *
 * ```kotlin
 * val result = kahanSumOf(components) { (_, x) -> x * ln(x) }
 * ```
 */
inline fun <T> kahanSumOf(values: Iterable<T>, selector: (T) -> Double): Double {
    var sum = 0.0
    var compensation = 0.0

    for (v in values) {
        val y = selector(v) - compensation
        val t = sum + y
        compensation = (t - sum) - y
        sum = t
    }

    return sum
}

/**
 * Mutable Kahan accumulator — use when values arrive incrementally.
 *
 * ```kotlin
 * val acc = KahanAccumulator()
 * for (component in components) {
 *     acc += component.gibbsContribution()
 * }
 * val total = acc.sum
 * ```
 */
class KahanAccumulator {
    private var sum = 0.0
    private var compensation = 0.0

    val result: Double get() = sum

    operator fun plusAssign(value: Double) {
        val y = value - compensation
        val t = sum + y
        compensation = (t - sum) - y
        sum = t
    }

    operator fun minusAssign(value: Double) {
        plusAssign(-value)
    }

    fun reset() {
        sum = 0.0
        compensation = 0.0
    }
}
