package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.error.Interval
import dev.hivens.solyx.core.units.JoulePerMole
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
 * Uses numerical differentiation of the Gibbs energy model.
 *
 * For an N-component system with x_ref as the dependent variable
 * (x_ref = 1 - Σ x_i), the exact formula is:
 * ```
 * μ_i   = G + (1 - x_i) * ∂G/∂x_i - Σ_{j≠i,j≠ref} x_j * ∂G/∂x_j   [i ≠ ref]
 * μ_ref = G - Σ_{j≠ref} x_j * ∂G/∂x_j
 * ```
 *
 * For binary systems this reduces to the standard tangent-line construction:
 * `μ_1 = G + (1 - x_1) * dG/dx_1`
 *
 * All partial derivatives are computed with the same reference component
 * absorbing the perturbation, ensuring the mole fraction constraint
 * Σ x_i = 1 is maintained exactly throughout differentiation.
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

        // Use last element as the reference (x_ref = 1 - Σ x_i).
        // All partials are computed by perturbing x_i and compensating x_ref.
        val reference = elements.last()

        // Compute ∂G/∂x_i for every non-reference component.
        val partials: Map<Element, Pair<Double, Double>> = elements
            .filter { it != reference }
            .associateWith { partialGibbs(state, it, reference) }

        return elements.map { element ->
            if (element == reference) {
                // μ_ref = G - Σ_{j≠ref} x_j * ∂G/∂x_j
                val mu = g - partials.entries.sumOf { (el, pd) ->
                    (state.composition[el]?.value ?: 0.0) * pd.first
                }
                val errorBound = partials.values.sumOf { it.second }
                ChemicalPotential(element, JoulePerMole(mu), Interval.withTolerance(mu, errorBound))
            } else {
                val (dGdxi, err) = partials[element]!!
                val xi = state.composition[element]?.value ?: 0.0

                // μ_i = G + (1 - x_i)*∂G/∂x_i - Σ_{j≠i, j≠ref} x_j * ∂G/∂x_j
                val crossTerms = partials.entries
                    .filter { (el, _) -> el != element }
                    .sumOf { (el, pd) -> (state.composition[el]?.value ?: 0.0) * pd.first }

                val mu = g + (1.0 - xi) * dGdxi - crossTerms
                ChemicalPotential(element, JoulePerMole(mu), Interval.withTolerance(mu, err * (1.0 - xi)))
            }
        }
    }

    /**
     * Compute ∂G/∂x_target numerically, keeping x_reference as the compensating
     * component so that Σ x_i = 1 is preserved exactly at every evaluation point.
     *
     * The step size is chosen to be the largest value that keeps both
     * x_target and x_reference strictly inside (ε, 1-ε) for all four
     * evaluation points required by Richardson extrapolation (±h, ±2h).
     *
     * @return Pair(derivative, Richardson error estimate)
     */
    private fun partialGibbs(
        state: SystemState,
        target: Element,
        reference: Element
    ): Pair<Double, Double> {
        val xi   = state.composition[target]?.value    ?: return Pair(0.0, 0.0)
        val xRef = state.composition[reference]?.value ?: return Pair(0.0, 0.0)

        val eps = 1e-10
        // Maximum safe step: Richardson needs ±2h, so divide safe range by 2.
        val safeH = minOf(
            1e-6,
            (xi   - eps) / 2.0,
            (xRef - eps) / 2.0,
            (1.0 - eps - xi)   / 2.0,
            (1.0 - eps - xRef) / 2.0
        )
        if (safeH <= 0.0) return Pair(0.0, 0.0)

        val perturb = { dx: Double ->
            val newComposition = state.composition.toMutableMap()
            newComposition[target]    = MoleFraction(xi + dx)
            newComposition[reference] = MoleFraction(xRef - dx)
            model.compute(state.copy(composition = newComposition)).value
        }

        val central = (perturb(safeH)       - perturb(-safeH))       / (2.0 * safeH)
        val coarse  = (perturb(2.0 * safeH) - perturb(-2.0 * safeH)) / (4.0 * safeH)
        val error   = kotlin.math.abs(central - coarse) / 3.0

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
