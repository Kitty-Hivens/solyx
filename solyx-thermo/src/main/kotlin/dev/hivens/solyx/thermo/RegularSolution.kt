package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.numeric.kahanSumOf
import dev.hivens.solyx.core.units.JoulePerMole
import kotlin.math.pow

/**
 * Redlich-Kister interaction parameter between two components.
 *
 * Models the excess Gibbs energy between a pair of species on a sublattice:
 *
 * ```
 * G_excess = x_A * x_B * Σ L^n * (x_A - x_B)^n
 * ```
 *
 * @param a first component
 * @param b second component
 * @param coefficients L^0, L^1, L^2, ... polynomial coefficients in J/mol
 */
data class RedlichKisterParameter(
    val a: Element,
    val b: Element,
    val coefficients: List<JoulePerMole>
) {
    init {
        require(coefficients.isNotEmpty()) { "At least L^0 must be provided" }
    }

    /**
     * Evaluate excess contribution for given mole fractions [xA] and [xB].
     *
     * G_excess_AB = x_A * x_B * Σ L^n * (x_A - x_B)^n
     */
    fun evaluate(xA: Double, xB: Double): Double {
        val polynomial = kahanSumOf(coefficients.withIndex()) { (n, L) ->
            L.value * (xA - xB).pow(n)
        }
        return xA * xB * polynomial
    }

    /** True if this parameter covers the [first]-[second] pair in any order. */
    fun matches(first: Element, second: Element) =
        (a == first && b == second) || (a == second && b == first)
}

/**
 * Regular solution model — ideal mixing + Redlich-Kister excess energy.
 *
 * Full Gibbs energy of a phase:
 *
 * ```
 * G = G_ref + G_mix + G_excess
 *
 * G_mix    = R * T * Σ(x_i * ln(x_i))       [via IdealMixing]
 * G_excess = Σ x_A * x_B * Σ L^n * (x_A - x_B)^n
 * ```
 *
 * For a symmetric regular solution (L^0 only, L^1 = 0):
 * - L^0 < 0 — components attract, mixing is favorable
 * - L^0 > 0 — components repel, phase separation tendency
 * - L^0 = 0 — ideal solution, no excess energy
 *
 * ```kotlin
 * val feCModel = RegularSolution(
 *     idealMixing = IdealMixing(),
 *     interactions = listOf(
 *         RedlichKisterParameter(
 *             a = Element.Fe,
 *             b = Element.C,
 *             coefficients = listOf((-10000.0).jPerMol, 500.0.jPerMol)
 *         )
 *     )
 * )
 * val g = feCModel.compute(state)
 * ```
 *
 * @param idealMixing ideal mixing contribution — injected for composability
 * @param interactions Redlich-Kister parameters for all interacting pairs
 */
class RegularSolution(
    private val idealMixing: GibbsModel = IdealMixing(),
    private val interactions: List<RedlichKisterParameter> = emptyList()
) : GibbsModel {

    override fun compute(state: SystemState): JoulePerMole {
        val gMix = idealMixing.compute(state)
        val gExcess = computeExcess(state)
        return JoulePerMole(gMix.value + gExcess)
    }

    private fun computeExcess(state: SystemState): Double =
        kahanSumOf(interactions) { param ->
            val xA = state.composition[param.a]?.value ?: return@kahanSumOf 0.0
            val xB = state.composition[param.b]?.value ?: return@kahanSumOf 0.0
            param.evaluate(xA, xB)
        }
}
