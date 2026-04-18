package dev.hivens.solyx.thermo

/**
 * Chemical element with its fundamental properties.
 *
 * Elements are the building blocks of thermodynamic systems in CALPHAD.
 * Each element carries its atomic mass and reference state — the stable
 * physical form at 298.15 K and 101 325 Pa (SER — Standard Element Reference).
 *
 * ```kotlin
 * val fe = Element.Fe
 * val c  = Element("C", atomicMass = 12.011, referenceState = "GRAPHITE")
 * ```
 */
data class Element(
    val symbol: String,
    val atomicMass: Double,
    val referenceState: String = "SER"
) {
    init {
        require(symbol.isNotBlank())  { "Element symbol cannot be blank" }
        require(atomicMass > 0.0)     { "Atomic mass must be positive: $atomicMass" }
    }

    override fun toString() = symbol

    companion object {
        // Common elements in alloy thermodynamics
        val H  = Element("H",  1.008,   "GAS")
        val C  = Element("C",  12.011,  "GRAPHITE")
        val N  = Element("N",  14.007,  "GAS")
        val O  = Element("O",  15.999,  "GAS")
        val Al = Element("Al", 26.982,  "FCC_A1")
        val Si = Element("Si", 28.085,  "DIAMOND_A4")
        val Ti = Element("Ti", 47.867,  "HCP_A3")
        val Cr = Element("Cr", 51.996,  "BCC_A2")
        val Mn = Element("Mn", 54.938,  "CBCC_A12")
        val Fe = Element("Fe", 55.845,  "BCC_A2")
        val Co = Element("Co", 58.933,  "HCP_A3")
        val Ni = Element("Ni", 58.693,  "FCC_A1")
        val Cu = Element("Cu", 63.546,  "FCC_A1")
        val Zn = Element("Zn", 65.38,   "HCP_A3")
        val Nb = Element("Nb", 92.906,  "BCC_A2")
        val Mo = Element("Mo", 95.95,   "BCC_A2")
        val W  = Element("W",  183.84,  "BCC_A2")
        val Pb = Element("Pb", 207.2,   "FCC_A1")
    }
}

/**
 * Species — an occupant of a sublattice site in the CEF model.
 *
 * A species is either a chemical element or a vacancy (Va).
 * Vacancies represent empty sites and are essential for modeling
 * interstitial solutions such as carbon in austenite.
 *
 * ```kotlin
 * val fe  = Species.of(Element.Fe)
 * val c   = Species.of(Element.C)
 * val va  = Species.Vacancy
 * ```
 */
sealed class Species {
    abstract val symbol: String

    /** Species backed by a chemical element. */
    data class ElementSpecies(val element: Element) : Species() {
        override val symbol get() = element.symbol
        override fun toString() = symbol
    }

    /** Vacancy — empty sublattice site. */
    data object Vacancy : Species() {
        override val symbol = "Va"
        override fun toString() = "Va"
    }

    companion object {
        fun of(element: Element): Species = ElementSpecies(element)

        /** Parse species from string — "Va" becomes [Vacancy], anything else looks up in [known]. */
        fun parse(symbol: String, known: Map<String, Element>): Species {
            if (symbol == "Va" || symbol == "VA" || symbol == "va") return Vacancy
            val element = known[symbol]
                ?: error("Unknown element symbol: $symbol. Add it to the known elements map.")
            return ElementSpecies(element)
        }
    }
}
