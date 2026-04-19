package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.PhysicalConstants.R
import dev.hivens.solyx.core.numeric.kahanSumOf
import kotlin.math.ln

/**
 * State of a thermodynamic system at a given condition.
 *
 * Passed to [GibbsModel.compute] to evaluate the molar Gibbs energy.
 *
 * @param composition mole fractions of each component — must sum to 1.0
 * @param temperature temperature in Kelvin
 */
data class SystemState(
    val composition: Map<Element, MoleFraction>,
    val temperature: Kelvin
) {
    init {
        val sum = kahanSumOf(composition.values) { it.value }
        require(sum in 0.9999..1.0001) {
            "Mole fractions must sum to 1.0, got $sum"
        }
    }
}

/**
 * Base interface for all molar Gibbs energy models.
 *
 * Every phase model in Solyx implements this interface.
 * The total Gibbs energy of a phase is:
 *
 * ```
 * G = G_ref + G_mix + G_excess
 * ```
 *
 * Each model may implement one or more of these contributions.
 *
 * ```kotlin
 * val model: GibbsModel = RegularSolution(omega = -10000.0.jPerMol)
 * val g = model.compute(state)
 * ```
 */
fun interface GibbsModel {
    /**
     * Compute molar Gibbs energy at the given [state].
     *
     * @param state composition and temperature of the system
     * @return molar Gibbs energy in J/mol
     */
    fun compute(state: SystemState): JoulePerMole
}

/**
 * Ideal mixing contribution to Gibbs energy.
 *
 * Models the entropy gain from randomly mixing components:
 *
 * ```
 * G_mix = R * T * Σ(x_i * ln(x_i))
 * ```
 *
 * Always negative — mixing is always thermodynamically favorable in an ideal solution.
 * Components with x_i = 0 are skipped (0 * ln(0) → 0 by convention).
 *
 * This is the foundation of all solution models in CALPHAD.
 * Real solutions add excess terms on top of ideal mixing.
 *
 * ```kotlin
 * val mixing = IdealMixing()
 * val g = mixing.compute(state)
 * ```
 */
class IdealMixing : GibbsModel {

    override fun compute(state: SystemState): JoulePerMole {
        val entropyTerm = kahanSumOf(state.composition.values) { fraction ->
            val x = fraction.value
            if (x > 0.0) x * ln(x) else 0.0
        }
        return JoulePerMole(R * state.temperature.value * entropyTerm)
    }
}

/**
 * Combines multiple [GibbsModel] contributions additively.
 *
 * Used to build composite models from independent terms:
 *
 * ```kotlin
 * val model = CompositeGibbsModel(
 *     ReferenceEnergy(endmembers),
 *     IdealMixing(),
 *     RedlichKister(params)
 * )
 * ```
 */
class CompositeGibbsModel(private vararg val models: GibbsModel) : GibbsModel {

    override fun compute(state: SystemState): JoulePerMole {
        val total = kahanSumOf(models.toList()) { model ->
            model.compute(state).value
        }
        return JoulePerMole(total)
    }
}
