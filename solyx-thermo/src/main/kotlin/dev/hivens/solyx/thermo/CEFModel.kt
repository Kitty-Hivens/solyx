package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.numeric.kahanSumOf
import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.PhysicalConstants.R
import kotlin.math.ln

/**
 * A single sublattice in the CEF model.
 *
 * Represents a set of crystallographically equivalent sites in the crystal
 * structure. Each sublattice has a fixed number of sites [stoichiometry] and
 * a set of species that can occupy those sites.
 *
 * ```kotlin
 * // Austenite: (Fe, Cr)_1 (C, Va)_1
 * val sublattice1 = Sublattice(
 *     stoichiometry = 1.0,
 *     occupancy = mapOf(Species.of(Element.Fe) to 0.9, Species.of(Element.Cr) to 0.1)
 * )
 * val sublattice2 = Sublattice(
 *     stoichiometry = 1.0,
 *     occupancy = mapOf(Species.of(Element.C) to 0.02, Species.Vacancy to 0.98)
 * )
 * ```
 *
 * @param stoichiometry number of sites per formula unit
 * @param occupancy site fractions y_i for each species — must sum to 1.0
 */
data class Sublattice(
    val stoichiometry: Double,
    val occupancy: Map<Species, Double>
) {
    init {
        require(stoichiometry > 0.0) { "Stoichiometry must be positive: $stoichiometry" }
        val sum = occupancy.values.sum()
        require(sum in 0.9999..1.0001) {
            "Site fractions must sum to 1.0, got $sum"
        }
        require(occupancy.values.all { it >= 0.0 }) {
            "Site fractions must be non-negative"
        }
    }

    /** Site fraction of [species] on this sublattice. */
    fun y(species: Species): Double = occupancy[species] ?: 0.0

    /** Ideal mixing entropy contribution for this sublattice. */
    fun idealEntropy(): Double = kahanSumOf(occupancy.values) { y ->
        if (y > 0.0) y * ln(y) else 0.0
    }
}

/**
 * Endmember parameter G° for a specific configuration across all sublattices.
 *
 * An endmember is a hypothetical state where each sublattice is occupied by
 * exactly one species. For (Fe,Cr)_1(C,Va)_1 there are four endmembers:
 * Fe:C, Fe:Va, Cr:C, Cr:Va.
 *
 * The colon notation `A:B` means species A on sublattice 1, B on sublattice 2.
 *
 * @param configuration one species per sublattice, in sublattice order
 * @param energy reference Gibbs energy G° for this configuration in J/mol
 */
data class EndmemberParameter(
    val configuration: List<Species>,
    val energy: JoulePerMole
)

/**
 * Compound Energy Formalism (CEF) — multi-sublattice Gibbs energy model.
 *
 * The total molar Gibbs energy consists of three contributions:
 *
 * ```
 * G = G_ref + G_mix + G_excess
 *
 * G_ref  = Σ (Π y_i^s) * G°(config)        — weighted endmember energies
 * G_mix  = R * T * Σ s_k * Σ y_i * ln(y_i) — ideal mixing per sublattice
 * G_excess = Σ interactions                  — Redlich-Kister on each sublattice
 * ```
 *
 * This is the standard model used in all modern CALPHAD databases (TDB files).
 *
 * ```kotlin
 * // Austenite (Fe,Cr)_1 (C,Va)_1
 * val austenite = CEFModel(
 *     endmembers = listOf(
 *         EndmemberParameter(listOf(Species.of(Fe), Species.of(C)),  (-10000.0).jPerMol),
 *         EndmemberParameter(listOf(Species.of(Fe), Species.Vacancy), (-8000.0).jPerMol),
 *         EndmemberParameter(listOf(Species.of(Cr), Species.of(C)),  (-12000.0).jPerMol),
 *         EndmemberParameter(listOf(Species.of(Cr), Species.Vacancy), (-9000.0).jPerMol)
 *     )
 * )
 * val g = austenite.compute(sublattices, temperature)
 * ```
 *
 * @param endmembers G° parameters for each configuration
 * @param interactions Redlich-Kister excess parameters per sublattice
 */
class CEFModel(
    private val endmembers: List<EndmemberParameter>,
    private val interactions: List<SublatticeInteraction> = emptyList()
) {
    /**
     * Compute molar Gibbs energy for the given sublattice occupancies.
     *
     * @param sublattices current site fractions for each sublattice
     * @param temperature temperature in Kelvin
     */
    fun compute(sublattices: List<Sublattice>, temperature: dev.hivens.solyx.core.units.Kelvin): JoulePerMole {
        val gRef    = computeReference(sublattices)
        val gMix    = computeIdealMixing(sublattices, temperature)
        val gExcess = computeExcess(sublattices)
        return JoulePerMole(gRef + gMix + gExcess)
    }

    /** G_ref = Σ (Π y_i^s) * G°(config) */
    private fun computeReference(sublattices: List<Sublattice>): Double =
        kahanSumOf(endmembers) { endmember ->
            require(endmember.configuration.isNotEmpty()) {
                "Endmember configuration must not be empty"
            }
            require(endmember.configuration.size <= sublattices.size) {
                "Endmember configuration has ${endmember.configuration.size} sublattices " +
                        "but model has only ${sublattices.size}"
            }
            val product = endmember.configuration.mapIndexed { s, species ->
                sublattices[s].y(species)
            }.fold(1.0, Double::times)
            product * endmember.energy.value
        }

    /** G_mix = R * T * Σ s_k * Σ y_i * ln(y_i) */
    private fun computeIdealMixing(
        sublattices: List<Sublattice>,
        temperature: dev.hivens.solyx.core.units.Kelvin
    ): Double {
        val entropySum = kahanSumOf(sublattices) { sub ->
            sub.stoichiometry * sub.idealEntropy()
        }
        return R * temperature.value * entropySum
    }

    /** G_excess = Σ interactions on each sublattice */
    private fun computeExcess(sublattices: List<Sublattice>): Double =
        kahanSumOf(interactions) { interaction ->
            val sub = sublattices[interaction.sublatticeIndex]
            val yA = sub.y(interaction.a)
            val yB = sub.y(interaction.b)
            interaction.parameter.evaluate(yA, yB)
        }
}

/**
 * Redlich-Kister excess interaction on a specific sublattice.
 *
 * @param sublatticeIndex which sublattice this interaction applies to
 * @param a first species
 * @param b second species
 * @param parameter Redlich-Kister coefficients
 */
data class SublatticeInteraction(
    val sublatticeIndex: Int,
    val a: Species,
    val b: Species,
    val parameter: RedlichKisterParameter
)
