package dev.hivens.solyx.calphad.tdb

import dev.hivens.solyx.core.units.Kelvin

/**
 * Parsed representation of a TDB thermodynamic database.
 *
 * This is the in-memory analogue of a .tdb file. All data is immutable
 * after parsing — modifications require a new database instance.
 */
data class TdbDatabase(
    val elements: Map<String, TdbElement>,
    val species:  Map<String, TdbSpecies>,
    val functions: Map<String, TdbFunction>,
    val phases:   Map<String, TdbPhase>,
    val parameters: List<TdbParameter>
) {
    /** Find all parameters for a given phase and type. */
    fun parametersFor(phaseName: String, type: ParameterType): List<TdbParameter> =
        parameters.filter { it.phaseName == phaseName && it.type == type }

    /** Find Gibbs energy parameter for a specific configuration. */
    fun gibbsParameter(phaseName: String, constituents: List<List<String>>): TdbParameter? =
        parameters.firstOrNull {
            it.type == ParameterType.G &&
            it.phaseName == phaseName &&
            it.constituentArray == constituents
        }
}

/**
 * Element declaration.
 *
 * ```
 * ELEMENT FE BCC_A2 5.5845E+01 4.4890E+03 2.7280E+01 !
 * ```
 */
data class TdbElement(
    val symbol:         String,
    val referenceState: String,
    val atomicMass:     Double,
    val h298MinusH0:    Double,    // H298 - H0 in J/mol
    val s298:           Double     // S298 in J/(mol·K)
)

/**
 * Species declaration — elemental or compound, possibly charged.
 *
 * ```
 * SPECIES FE+2  FE/+2 !
 * SPECIES O-2   O/-2  !
 * SPECIES SIO2  SI1O2 !
 * ```
 */
data class TdbSpecies(
    val name:    String,
    val formula: String,
    val charge:  Int = 0     // 0 for neutral, positive for cations, negative for anions
)

/**
 * Named Gibbs energy function with temperature ranges.
 *
 * ```
 * FUNCTION GHSERFE 298.15 1225.7+124.134*T-23.5143*T*LN(T); 1811 Y
 *  -25383.581+299.31255*T-46*T*LN(T); 6000 N !
 * ```
 */
data class TdbFunction(
    val name:   String,
    val ranges: List<TdbTemperatureRange>
) {
    /** Evaluate the function at [temperature]. Returns null if out of range. */
    fun evaluate(temperature: Kelvin): Double? {
        val range = ranges.firstOrNull {
            temperature.value >= it.tLow && temperature.value <= it.tHigh
        } ?: return null
        return range.expression.evaluate(temperature.value)
    }
}

/**
 * A single temperature range within a function or parameter.
 *
 * @param tLow lower temperature bound (inclusive)
 * @param tHigh upper temperature bound (inclusive)
 * @param expression mathematical expression for G(T)
 */
data class TdbTemperatureRange(
    val tLow:       Double,
    val tHigh:      Double,
    val expression: TdbExpression
)

/**
 * Phase declaration.
 *
 * ```
 * PHASE BCC_A2 %B 2 1 3 !
 * CONSTITUENT BCC_A2 :FE,CR : VA : !
 * ```
 *
 * @param name phase name
 * @param typeFlags type definition flags (e.g. 'B' for BCC magnetic)
 * @param sublattices stoichiometries of each sublattice
 * @param constituents species allowed on each sublattice
 * @param magneticAfm antiferromagnetic factor from TYPE_DEF (-1 BCC, -3 FCC/HCP)
 * @param magneticP structure constant p from TYPE_DEF (0.4 BCC, 0.28 FCC/HCP)
 */
data class TdbPhase(
    val name:         String,
    val typeFlags:    String = "",
    val sublattices:  List<Double>,
    val constituents: List<List<String>> = emptyList(),
    val magneticAfm:  Double? = null,
    val magneticP:    Double? = null
) {
    val sublatticeCount get() = sublattices.size
}

/**
 * Type of thermodynamic parameter.
 */
enum class ParameterType {
    G,      // Gibbs energy endmember
    L,      // Redlich-Kister interaction
    TC,     // Curie/Néel temperature
    BMAGN,  // Bohr magneton number
    V0,     // Molar volume
    UNKNOWN
}

/**
 * Thermodynamic parameter — G, L, TC, or BMAGN.
 *
 * ```
 * PARAMETER G(BCC_A2,FE:VA;0)    298.15 GHSERFE; 6000 N !
 * PARAMETER L(LIQUID,CU,MG;0)    298.15 -34000+4*T; 3000 N !
 * PARAMETER TC(BCC_A2,FE:VA;0)   298.15 1043; 6000 N !
 * PARAMETER BMAGN(BCC_A2,FE:VA;0) 298.15 2.22; 6000 N !
 * ```
 *
 * @param type parameter type
 * @param phaseName phase this parameter belongs to
 * @param constituentArray sublattice configuration — outer list = sublattices, inner = species
 * @param order Redlich-Kister polynomial order (0, 1, 2, ...)
 * @param ranges temperature-dependent expression ranges
 */
data class TdbParameter(
    val type:             ParameterType,
    val phaseName:        String,
    val constituentArray: List<List<String>>,
    val order:            Int = 0,
    val ranges:           List<TdbTemperatureRange>
) {
    /** Evaluate the parameter at [temperature]. */
    fun evaluate(temperature: Kelvin): Double? {
        val range = ranges.firstOrNull {
            temperature.value >= it.tLow && temperature.value <= it.tHigh
        } ?: return null
        return range.expression.evaluate(temperature.value)
    }
}
