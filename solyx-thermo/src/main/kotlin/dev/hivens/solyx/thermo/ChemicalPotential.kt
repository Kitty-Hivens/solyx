package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.error.Interval
import dev.hivens.solyx.core.numeric.derivative
import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.MoleFraction
import dev.hivens.solyx.core.units.PhysicalConstants.R
import kotlin.math.exp
import kotlin.math.ln

/**
 * Chemical potential of a component in a phase.
 *
 * Defined as the partial molar Gibbs energy:
 * ```
 * μ_i = ∂G/∂n_i  (at constant T, P, n_j≠i)
 * ```
 *
 * Physical meaning — the change in total Gibbs energy when one mole
 * of component i is added to the system. At equilibrium, chemical
 * potential of each component is equal across all phases:
 * ```
 * μ_i(phase α) = μ_i(phase β) = μ_i(phase γ)
 * ```
 * This is the fundamental condition that [GibbsMinimizer] enforces.
 *
 * @param element component this potential belongs to
 * @param value   chemical potential in J/mol
 * @param uncertainty numerical uncertainty from differentiation
 */
data class ChemicalPotential(
    val element: Element,
    val value: JoulePerMole,
    val uncertainty: Interval = Interval(value.value)
) {
    override fun toString() = "μ(${element.symbol}) = ${value.value} J/mol"
}

/**
 * Computes chemical potentials for all components in a system.
 *
 * Uses numerical differentiation of the Gibbs energy model:
 * ```
 * μ_i = G + (1 - x_i) * ∂G/∂x_i   [binary tangent line construction]
 * ```
 *
 * For multicomponent systems, uses partial derivatives with respect
 * to each mole fraction, corrected for the constraint Σx_i = 1.
 *
 * ```kotlin
 * val calculator = ChemicalPotentialCalculator(regularSolution)
 * val potentials = calculator.compute(state)
 * potentials.forEach { println(it) }
 * ```
 *
 * @param model Gibbs energy model to differentiate
 */
class ChemicalPotentialCalculator(private val model: GibbsModel) {

    /**
     * Compute chemical potential for all components at [state].
     *
     * @return chemical potential for each element in the system
     */
    fun compute(state: SystemState): List<ChemicalPotential> {
        val elements = state.composition.keys.toList()
        val g = model.compute(state).value

        return elements.map { element ->
            val (dGdx, error) = partialGibbs(state, element)
            val xi = state.composition[element]?.value ?: 0.0

            // μ_i = G + (1 - x_i) * ∂G/∂x_i
            val mu = g + (1.0 - xi) * dGdx
            val uncertainty = Interval(mu, error * (1.0 - xi))

            ChemicalPotential(
                element = element,
                value = JoulePerMole(mu),
                uncertainty = uncertainty
            )
        }
    }

    /**
     * Compute partial derivative ∂G/∂x_i numerically.
     *
     * Perturbs mole fraction of [target] while adjusting a reference
     * component to maintain Σx_i = 1.
     *
     * @return Pair(derivative, error estimate)
     */
    private fun partialGibbs(state: SystemState, target: Element): Pair<Double, Double> {
        val elements = state.composition.keys.toList()
        val reference = elements.firstOrNull { it != target }
            ?: return Pair(0.0, 0.0)

        val xi = state.composition[target]?.value ?: return Pair(0.0, 0.0)

        val g = { dx: Double ->
            val newComposition = state.composition.toMutableMap()
            val newXi = (xi + dx).coerceIn(1e-10, 1.0 - 1e-10)
            val xRef  = state.composition[reference]?.value ?: 0.0
            val newXRef = (xRef - dx).coerceIn(1e-10, 1.0 - 1e-10)

            newComposition[target]    = MoleFraction(newXi)
            newComposition[reference] = MoleFraction(newXRef)

            val newState = state.copy(composition = newComposition)
            model.compute(newState).value
        }

        val h = 1e-6
        val central  = (g(h) - g(-h)) / (2.0 * h)
        val coarse   = (g(2 * h) - g(-2 * h)) / (4.0 * h)
        val error    = kotlin.math.abs(central - coarse) / 3.0

        return Pair(central, error)
    }
}

/**
 * Activity coefficient of a component — measure of deviation from ideal behavior.
 *
 * Defined from chemical potential:
 * ```
 * μ_i = G°_i + R*T*ln(a_i)
 * a_i = γ_i * x_i
 * ```
 *
 * Where:
 * - a_i — thermodynamic activity
 * - γ_i — activity coefficient
 * - x_i — mole fraction
 *
 * Interpretation:
 * - γ_i = 1.0 — ideal solution (Raoult's law)
 * - γ_i > 1.0 — positive deviation, components repel
 * - γ_i < 1.0 — negative deviation, components attract
 *
 * @param element   component
 * @param value     activity coefficient γ_i (dimensionless)
 * @param activity  thermodynamic activity a_i = γ_i * x_i
 */
data class ActivityCoefficient(
    val element: Element,
    val value: Double,
    val activity: Double
) {
    /** True if behavior is approximately ideal. */
    val isIdeal: Boolean get() = kotlin.math.abs(value - 1.0) < 1e-4

    override fun toString() = "γ(${element.symbol}) = $value, a = $activity"
}

/**
 * Computes activity coefficients from chemical potentials.
 *
 * Requires reference Gibbs energies G°_i for each pure component.
 *
 * ```kotlin
 * val calculator = ActivityCalculator(
 *     potentialCalculator = ChemicalPotentialCalculator(model),
 *     referenceEnergies   = mapOf(Element.Fe to (-8000.0).jPerMol)
 * )
 * val activities = calculator.compute(state)
 * ```
 */
class ActivityCalculator(
    private val potentialCalculator: ChemicalPotentialCalculator,
    private val referenceEnergies: Map<Element, JoulePerMole>
) {
    /**
     * Compute activity coefficients for all components at [state].
     */
    fun compute(state: SystemState): List<ActivityCoefficient> {
        val potentials = potentialCalculator.compute(state)

        return potentials.mapNotNull { potential ->
            val g0 = referenceEnergies[potential.element] ?: return@mapNotNull null
            val xi = state.composition[potential.element]?.value ?: return@mapNotNull null

            if (xi <= 0.0) return@mapNotNull null

            // ln(γ_i) = (μ_i - G°_i) / (R*T) - ln(x_i)
            val lnGamma = (potential.value.value - g0.value) /
                          (R * state.temperature.value) - ln(xi)
            val gamma    = exp(lnGamma)
            val activity = gamma * xi

            ActivityCoefficient(
                element  = potential.element,
                value    = gamma,
                activity = activity
            )
        }
    }
}
