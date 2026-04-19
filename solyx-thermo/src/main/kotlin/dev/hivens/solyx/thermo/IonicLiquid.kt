package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.numeric.kahanSumOf
import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.PhysicalConstants.R
import kotlin.math.ln

/**
 * Charged species for the ionic liquid model.
 *
 * Unlike [Species] used in CEF, ionic species carry an explicit charge
 * which is essential for maintaining electroneutrality in the liquid phase.
 *
 * ```kotlin
 * val fe2 = IonicSpecies.Cation(Element.Fe, charge = 2)
 * val o2  = IonicSpecies.Anion(Element.O,  charge = 2)
 * val va  = IonicSpecies.Vacancy
 * ```
 */
sealed class IonicSpecies {
    abstract val symbol: String

    /** Positively charged ion — occupies the cation sublattice. */
    data class Cation(val element: Element, val charge: Int) : IonicSpecies() {
        init { require(charge > 0) { "Cation charge must be positive: $charge" } }
        override val symbol get() = "${element.symbol}${charge}+"
        override fun toString() = symbol
    }

    /** Negatively charged ion — occupies the anion sublattice. */
    data class Anion(val element: Element, val charge: Int) : IonicSpecies() {
        init { require(charge > 0) { "Anion charge must be positive (sign is implicit): $charge" } }
        override val symbol get() = "${element.symbol}${charge}-"
        override fun toString() = symbol
    }

    /** Vacancy on the anion sublattice — represents metallic bonding character. */
    data object Vacancy : IonicSpecies() {
        override val symbol = "Va"
        override fun toString() = "Va"
    }
}

/**
 * State of the ionic liquid — site fractions on each sublattice.
 *
 * @param cationOccupancy site fractions y_i for each cation — must sum to 1.0
 * @param anionOccupancy  site fractions y_j for each anion/vacancy — must sum to 1.0
 */
data class IonicState(
    val cationOccupancy: Map<IonicSpecies.Cation, Double>,
    val anionOccupancy:  Map<IonicSpecies, Double>,       // Anion or Vacancy
    val temperature: Kelvin
) {
    init {
        val cSum = kahanSumOf(cationOccupancy.values) { it }
        require(cSum in 0.9999..1.0001) { "Cation fractions must sum to 1.0, got $cSum" }
        val aSum = kahanSumOf(anionOccupancy.values) { it }
        require(aSum in 0.9999..1.0001) { "Anion fractions must sum to 1.0, got $aSum" }
    }

    /**
     * P — effective size of the cation sublattice.
     * Computed as the mean charge of anions to maintain electroneutrality.
     *
     * P = Σ y_j * |z_j|   (sum over anions only, vacancies contribute 0)
     */
    val p: Double get() = kahanSumOf(anionOccupancy.entries) { (species, y) ->
        when (species) {
            is IonicSpecies.Anion   -> y * species.charge
            is IonicSpecies.Vacancy -> 0.0
            else                    -> 0.0
        }
    }

    /**
     * Q — effective size of the anion sublattice.
     * Computed as the mean charge of cations to maintain electroneutrality.
     *
     * Q = Σ y_i * z_i   (sum over cations)
     */
    val q: Double get() = kahanSumOf(cationOccupancy.entries) { (species, y) ->
        y * species.charge
    }
}

/**
 * Endmember parameter for the ionic liquid model.
 *
 * @param cation cation species
 * @param anion  anion or vacancy species
 * @param energy reference Gibbs energy G° in J/mol
 */
data class IonicEndmember(
    val cation: IonicSpecies.Cation,
    val anion:  IonicSpecies,
    val energy: JoulePerMole
)

/**
 * Two-sublattice ionic liquid model.
 *
 * Extends CEF for liquid phases where species carry electric charges.
 * The sublattice sizes P and Q are composition-dependent to enforce
 * electroneutrality at all compositions.
 *
 * Structure:
 * ```
 * (cations)_P  (anions, Va)_Q
 * ```
 *
 * Gibbs energy:
 * ```
 * G = G_ref + G_mix + G_excess
 *
 * G_ref = Σ_i Σ_j  y_i * y_j * G°(i:j)
 *
 * G_mix = R*T * [ P * Σ_i y_i * ln(y_i)      <- cation sublattice
 *               + Q * Σ_j y_j * ln(y_j) ]     <- anion sublattice
 *
 * G_excess = Σ interactions (Redlich-Kister on each sublattice)
 * ```
 *
 * Typical use — oxide and sulfide slags, ionic salt melts.
 *
 * ```kotlin
 * val slag = IonicLiquid(
 *     endmembers = listOf(
 *         IonicEndmember(IonicSpecies.Cation(Fe, 2), IonicSpecies.Anion(O, 2),  (-250000.0).jPerMol),
 *         IonicEndmember(IonicSpecies.Cation(Fe, 2), IonicSpecies.Vacancy,       (-30000.0).jPerMol),
 *         IonicEndmember(IonicSpecies.Cation(Ca, 2), IonicSpecies.Anion(O, 2),  (-310000.0).jPerMol),
 *         IonicEndmember(IonicSpecies.Cation(Ca, 2), IonicSpecies.Vacancy,       (-20000.0).jPerMol)
 *     )
 * )
 * val g = slag.compute(state)
 * ```
 */
class IonicLiquid(
    private val endmembers: List<IonicEndmember>,
    private val cationInteractions: List<CationInteraction> = emptyList(),
    private val anionInteractions:  List<AnionInteraction>  = emptyList()
) {
    /**
     * Compute molar Gibbs energy for the given ionic state.
     */
    fun compute(state: IonicState): JoulePerMole {
        val gRef    = computeReference(state)
        val gMix    = computeIdealMixing(state)
        val gExcess = computeExcess(state)
        return JoulePerMole(gRef + gMix + gExcess)
    }

    /** G_ref = Σ_i Σ_j y_i * y_j * G°(i:j) */
    private fun computeReference(state: IonicState): Double =
        kahanSumOf(endmembers) { em ->
            val yi = state.cationOccupancy[em.cation] ?: 0.0
            val yj = state.anionOccupancy[em.anion]   ?: 0.0
            yi * yj * em.energy.value
        }

    /**
     * G_mix = R*T * [ P * Σ y_i*ln(y_i)  +  Q * Σ y_j*ln(y_j) ]
     *
     * P and Q are composition-dependent — computed from [IonicState.p] and [IonicState.q].
     */
    private fun computeIdealMixing(state: IonicState): Double {
        val cationEntropy = kahanSumOf(state.cationOccupancy.values) { y ->
            if (y > 0.0) y * ln(y) else 0.0
        }
        val anionEntropy = kahanSumOf(state.anionOccupancy.values) { y ->
            if (y > 0.0) y * ln(y) else 0.0
        }
        return R * state.temperature.value * (state.p * cationEntropy + state.q * anionEntropy)
    }

    /** G_excess — Redlich-Kister on cation and anion sublattices independently. */
    private fun computeExcess(state: IonicState): Double {
        val cationExcess = kahanSumOf(cationInteractions) { interaction ->
            val yA = state.cationOccupancy[interaction.a] ?: 0.0
            val yB = state.cationOccupancy[interaction.b] ?: 0.0
            interaction.parameter.evaluate(yA, yB)
        }
        val anionExcess = kahanSumOf(anionInteractions) { interaction ->
            val yA = state.anionOccupancy[interaction.a] ?: 0.0
            val yB = state.anionOccupancy[interaction.b] ?: 0.0
            interaction.parameter.evaluate(yA, yB)
        }
        return cationExcess + anionExcess
    }
}

/**
 * Redlich-Kister interaction between two cations on the cation sublattice.
 *
 * Using a dedicated type instead of a generic `IonicInteraction(a: IonicSpecies, ...)`
 * prevents ClassCastException at runtime — the compiler enforces that only
 * cations can be passed here.
 */
data class CationInteraction(
    val a: IonicSpecies.Cation,
    val b: IonicSpecies.Cation,
    val parameter: RedlichKisterParameter
)

/**
 * Redlich-Kister interaction between two anion/vacancy species on the anion sublattice.
 */
data class AnionInteraction(
    val a: IonicSpecies,   // Anion or Vacancy
    val b: IonicSpecies,
    val parameter: RedlichKisterParameter
) {
    init {
        require(a is IonicSpecies.Anion || a is IonicSpecies.Vacancy) {
            "AnionInteraction.a must be Anion or Vacancy, got: $a"
        }
        require(b is IonicSpecies.Anion || b is IonicSpecies.Vacancy) {
            "AnionInteraction.b must be Anion or Vacancy, got: $b"
        }
    }
}
