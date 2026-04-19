package dev.hivens.solyx.calphad.tdb

import java.io.File
import java.io.Reader
import java.io.StringReader

/**
 * Parser for Thermo-Calc TDB thermodynamic database files.
 *
 * TDB is a text format with keyword-driven entries terminated by '!'.
 * This parser handles the core subset required for CALPHAD calculations:
 * ELEMENT, SPECIES, FUNCTION, TYPE_DEF, PHASE, CONSTITUENT, PARAMETER.
 *
 * Usage:
 * ```kotlin
 * val db = TdbParser.parse(File("steel.tdb"))
 * val db = TdbParser.parse(tdbString)
 * ```
 *
 * Known limitations:
 * - OPTIMIZATION blocks are ignored (used for parameter fitting, not calculation)
 * - ASSESSED_SYSTEMS is ignored
 * - Database references after N in parameters are ignored
 * - Maximum line length of 80 chars is not enforced (we handle continuation lines)
 */
object TdbParser {

    fun parse(file: File): TdbDatabase = parse(file.reader())
    fun parse(text: String): TdbDatabase = parse(StringReader(text))

    fun parse(reader: Reader): TdbDatabase {
        val raw = reader.readText()
        val entries = tokenize(raw)

        val elements   = mutableMapOf<String, TdbElement>()
        val species    = mutableMapOf<String, TdbSpecies>()
        val functions  = mutableMapOf<String, TdbFunction>()
        val typeDefs   = mutableMapOf<String, TdbTypeDef>()
        val phases     = mutableMapOf<String, TdbPhase>()
        val parameters = mutableListOf<TdbParameter>()

        for (entry in entries) {
            val tokens = entry.trim().split("\\s+".toRegex()).filter { it.isNotEmpty() }
            if (tokens.isEmpty()) continue

            when (tokens[0].uppercase()) {
                "ELEMENT"   -> parseElement(tokens)?.let { elements[it.symbol] = it }
                "SPECIES"   -> parseSpecies(tokens)?.let { species[it.name] = it }
                "FUNCTION"  -> parseFunction(entry, functions)
                "TYPE_DEF", "TYPE_DEFINITION" -> parseTypeDef(tokens)?.let { typeDefs[it.symbol] = it }
                "PHASE"     -> parsePhase(tokens, typeDefs)?.let { phases[it.name] = it }
                "CONSTITUENT" -> parseConstituent(tokens, phases)
                "PARAMETER" -> parseParameter(entry, functions)?.let { parameters.add(it) }
                // Ignored keywords
                "OPTIMIZATION", "ASSESSED_SYSTEMS", "ADD_REFERENCES",
                "DATABASE_INFO", "REFERENCE_FILE", "LIST_OF_REFERENCES" -> {}
                else -> {} // Unknown keyword — silently skip
            }
        }

        return TdbDatabase(elements, species, functions, phases, parameters)
    }

    // -------------------------------------------------------------------------
    // Tokenizer — splits raw TDB text into entries by '!'
    // -------------------------------------------------------------------------

    private fun tokenize(raw: String): List<String> {
        val withoutComments = raw.lines().joinToString("\n") { line ->
            // $ at start of line = comment
            if (line.trimStart().startsWith("$")) "" else line
        }

        return withoutComments
            .split("!")
            .map { it.trim() }
            .filter { it.isNotBlank() }
    }

    // -------------------------------------------------------------------------
    // ELEMENT parser
    // -------------------------------------------------------------------------

    /**
     * ELEMENT FE BCC_A2 5.5845E+01 4.4890E+03 2.7280E+01
     * ELEMENT VA VACUUM 0.0 0.0 0.0
     */
    private fun parseElement(tokens: List<String>): TdbElement? {
        if (tokens.size < 3) return null
        val symbol = tokens[1].uppercase()
        val refState = tokens[2].uppercase()
        val atomicMass  = tokens.getOrNull(3)?.toDoubleOrNull() ?: 0.0
        val h298MinusH0 = tokens.getOrNull(4)?.toDoubleOrNull() ?: 0.0
        val s298        = tokens.getOrNull(5)?.toDoubleOrNull() ?: 0.0
        return TdbElement(symbol, refState, atomicMass, h298MinusH0, s298)
    }

    // -------------------------------------------------------------------------
    // SPECIES parser
    // -------------------------------------------------------------------------

    /**
     * SPECIES FE+2  FE/+2
     * SPECIES O-2   O/-2
     * SPECIES SIO2  SI1O2
     */
    private fun parseSpecies(tokens: List<String>): TdbSpecies? {
        if (tokens.size < 3) return null
        val name    = tokens[1].uppercase()
        val formula = tokens[2].uppercase()
        val charge  = parseCharge(name)
        return TdbSpecies(name, formula, charge)
    }

    private fun parseCharge(name: String): Int {
        val plusIdx  = name.lastIndexOf('+')
        val minusIdx = name.lastIndexOf('-')
        return when {
            plusIdx > 0  -> name.substring(plusIdx + 1).toIntOrNull() ?: 1
            minusIdx > 0 -> -(name.substring(minusIdx + 1).toIntOrNull() ?: 1)
            else         -> 0
        }
    }

    // -------------------------------------------------------------------------
    // FUNCTION parser
    // -------------------------------------------------------------------------

    /**
     * FUNCTION GHSERFE 298.15 1225.7+124.134*T-23.5143*T*LN(T); 1811 Y
     *  -25383.581+299.31255*T-46*T*LN(T); 6000 N
     */
    private fun parseFunction(entry: String, functions: MutableMap<String, TdbFunction>) {
        val parts = entry.trim().split("\\s+".toRegex(), limit = 3)
        if (parts.size < 3) return
        val name = parts[1].uppercase()
        val rest = parts[2]
        val ranges = parseTemperatureRanges(rest, functions)
        if (ranges.isNotEmpty()) {
            functions[name] = TdbFunction(name, ranges)
        }
    }

    // -------------------------------------------------------------------------
    // TYPE_DEF / TYPE_DEFINITION parser
    // -------------------------------------------------------------------------

    /**
     * TYPE_DEF B GES A_P_D BCC_A2 MAGNETIC -1 0.400
     * TYPE_DEF A GES A_P_D FCC_A1 MAGNETIC -3 0.280
     */
    private fun parseTypeDef(tokens: List<String>): TdbTypeDef? {
        if (tokens.size < 2) return null
        val symbol = tokens[1]
        val magneticIdx = tokens.indexOfFirst { it.uppercase() == "MAGNETIC" }
        if (magneticIdx >= 0 && tokens.size > magneticIdx + 2) {
            val afm = tokens[magneticIdx + 1].toDoubleOrNull() ?: return null
            val p   = tokens[magneticIdx + 2].trimEnd(',').toDoubleOrNull() ?: return null
            return TdbTypeDef(symbol, afm, p)
        }
        return TdbTypeDef(symbol)
    }

    private data class TdbTypeDef(
        val symbol: String,
        val magneticAfm: Double? = null,
        val magneticP:   Double? = null
    )

    // -------------------------------------------------------------------------
    // PHASE parser
    // -------------------------------------------------------------------------

    /**
     * PHASE BCC_A2 %B 2 1 3
     * PHASE LIQUID % 1 1.0
     * PHASE SIGMA % 3 8 4 18
     */
    private fun parsePhase(tokens: List<String>, typeDefs: Map<String, TdbTypeDef>): TdbPhase? {
        if (tokens.size < 4) return null
        val name      = tokens[1].uppercase()
        val typeField = tokens[2]  // e.g. "%B", "%A", "%"
        val typeFlags = typeField.removePrefix("%")

        // Number of sublattices
        val nSublattices = tokens[3].toIntOrNull() ?: return null
        val stoich = (0 until nSublattices).mapNotNull {
            tokens.getOrNull(4 + it)?.toDoubleOrNull()
        }

        // Resolve magnetic parameters from TYPE_DEF
        var magneticAfm: Double? = null
        var magneticP:   Double? = null
        for (flag in typeFlags) {
            val td = typeDefs[flag.toString()]
            if (td?.magneticAfm != null) {
                magneticAfm = td.magneticAfm
                magneticP   = td.magneticP
            }
        }

        return TdbPhase(name, typeFlags, stoich, emptyList(), magneticAfm, magneticP)
    }

    // -------------------------------------------------------------------------
    // CONSTITUENT parser
    // -------------------------------------------------------------------------

    /**
     * CONSTITUENT BCC_A2 :FE,CR : VA :
     * CONSTITUENT FCC_A1 :AL,NI : VA :
     * CONSTITUENT LIQUID :CU,MG :
     */
    private fun parseConstituent(tokens: List<String>, phases: MutableMap<String, TdbPhase>) {
        if (tokens.size < 2) return
        val phaseName = tokens[1].uppercase()
        val phase = phases[phaseName] ?: return

        // Join remaining tokens and split by ':'
        val rest = tokens.drop(2).joinToString(" ")
        val sublattices = rest.split(":")
            .map { it.trim() }
            .filter { it.isNotBlank() }
            .map { sublattice ->
                sublattice.split(",")
                    .map { it.trim().trimEnd('%').uppercase() }
                    .filter { it.isNotBlank() }
            }

        phases[phaseName] = phase.copy(constituents = sublattices)
    }

    // -------------------------------------------------------------------------
    // PARAMETER parser
    // -------------------------------------------------------------------------

    /**
     * PARAMETER G(BCC_A2,FE:VA;0) 298.15 GHSERFE; 6000 N
     * PARAMETER L(LIQUID,CU,MG;0) 298.15 -34000+4*T; 3000 N
     * PARAMETER TC(BCC_A2,FE:VA;0) 298.15 1043; 6000 N
     * PARAMETER BMAGN(BCC_A2,FE:VA;0) 298.15 2.22; 6000 N
     */
    private fun parseParameter(entry: String, functions: Map<String, TdbFunction>): TdbParameter? {
        val parts = entry.trim().split("\\s+".toRegex(), limit = 3)
        if (parts.size < 3) return null

        val declaration = parts[1]  // e.g. G(BCC_A2,FE:VA;0)
        val rest        = parts[2]  // e.g. "298.15 GHSERFE; 6000 N"

        // Parse declaration
        val parenOpen  = declaration.indexOf('(')
        val parenClose = declaration.lastIndexOf(')')
        if (parenOpen < 0 || parenClose < 0) return null

        val typeStr = declaration.substring(0, parenOpen).uppercase()
        val type = when (typeStr) {
            "G"     -> ParameterType.G
            "L"     -> ParameterType.L
            "TC"    -> ParameterType.TC
            "BMAGN" -> ParameterType.BMAGN
            "V0", "VA", "VB", "VC", "VK" -> ParameterType.V0
            else    -> ParameterType.UNKNOWN
        }
        if (type == ParameterType.UNKNOWN) return null

        val inner = declaration.substring(parenOpen + 1, parenClose)

        // Split by ';' to get constituents and order
        val semiIdx = inner.lastIndexOf(';')
        val order = if (semiIdx >= 0) inner.substring(semiIdx + 1).toIntOrNull() ?: 0 else 0
        val constituentStr = if (semiIdx >= 0) inner.substring(0, semiIdx) else inner

        // First token is phase name, rest are constituents separated by ':'
        val constituentParts = constituentStr.split(",", limit = 2)
        val phaseName = constituentParts[0].uppercase()

        // For L parameters: L(LIQUID,CU,MG;0) — constituents on single sublattice
        // For G parameters: G(BCC_A2,FE:VA;0) — sublattices separated by ':'
        val constituentArray: List<List<String>> = if (constituentParts.size > 1) {
            val rest2 = constituentParts[1]
            if (rest2.contains(':')) {
                // CEF sublattice notation: FE:VA
                rest2.split(":").map { sub ->
                    sub.split(",").map { it.trim().uppercase() }.filter { it.isNotBlank() }
                }
            } else {
                // Single sublattice: CU,MG
                listOf(rest2.split(",").map { it.trim().uppercase() }.filter { it.isNotBlank() })
            }
        } else {
            emptyList()
        }

        val ranges = parseTemperatureRanges(rest, functions)
        if (ranges.isEmpty()) return null

        return TdbParameter(type, phaseName, constituentArray, order, ranges)
    }

    // -------------------------------------------------------------------------
    // Temperature range parser — shared between FUNCTION and PARAMETER
    // -------------------------------------------------------------------------

    /**
     * Parses: "298.15 expr1; 1357.77 Y expr2; 3200 N"
     * Into a list of TdbTemperatureRange.
     *
     * TDB layout after splitting by ';':
     *   segment 0 : "tLow expr1"
     *   segment 1 : "tHigh Y expr2"   <- upper bound, continuation flag, AND next expression
     *   segment 2 : "tHigh N"          <- upper bound of last range (no expression)
     *
     * The 'Y' means another range follows; 'N' means last range.
     * Reference labels after 'N' (e.g. "N REF123") are ignored.
     */
    private fun parseTemperatureRanges(
        input: String,
        functions: Map<String, TdbFunction>
    ): List<TdbTemperatureRange> {
        val ranges = mutableListOf<TdbTemperatureRange>()
        // Split on ';' — each segment is: tLow expression or tHigh Y/N
        val segments = input.split(";").map { it.trim() }.filter { it.isNotBlank() }
        if (segments.isEmpty()) return ranges

        // First segment: "tLow expr"
        var idx = 0
        val firstParts = segments[idx].split("\\s+".toRegex(), limit = 2)
        var tLow = firstParts[0].toDoubleOrNull() ?: return ranges
        var exprStr = firstParts.getOrNull(1) ?: return ranges
        idx++

        while (idx < segments.size) {
            // Segment: "tHigh Y/N [nextExpr]"
            // Split into at most 3 parts: tHigh, flag, rest-of-line (= expression for next range)
            val boundParts = segments[idx].trim().split("\\s+".toRegex(), limit = 3)
            val tHigh = boundParts[0].toDoubleOrNull() ?: break
            val continuation = boundParts.getOrNull(1)?.uppercase()
            idx++

            val expression = TdbExpressionParser.parse(exprStr).withFunctions(functions)
            ranges.add(TdbTemperatureRange(tLow, tHigh, expression))

            if (continuation != "Y") break

            // Y: expression for the next range is the remainder of this same segment (boundParts[2])
            tLow = tHigh
            exprStr = boundParts.getOrNull(2) ?: break
        }

        return ranges
    }
}
