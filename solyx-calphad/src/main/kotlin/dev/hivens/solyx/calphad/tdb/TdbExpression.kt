package dev.hivens.solyx.calphad.tdb

import kotlin.math.ln
import kotlin.math.log10
import kotlin.math.pow

/**
 * Parsed mathematical expression from a TDB file.
 *
 * TDB expressions are temperature-dependent polynomials of the form:
 * ```
 * A + B*T + C*T*LN(T) + D*T**2 + E*T**3 + F*T**(-1) + G*T**(-9)
 * ```
 *
 * Expressions can also reference named functions defined earlier in the TDB:
 * ```
 * +GHSERFE+500-0.3*T
 * ```
 *
 * The expression is stored as a list of terms which are evaluated at runtime.
 */
class TdbExpression(
    private val terms: List<Term>,
    private val functionRefs: Map<String, TdbFunction> = emptyMap()
) {
    /**
     * Evaluate expression at temperature [t] in Kelvin.
     */
    fun evaluate(t: Double): Double =
        terms.sumOf { it.evaluate(t, functionRefs) }

    /**
     * Return a copy with resolved function references.
     */
    fun withFunctions(functions: Map<String, TdbFunction>): TdbExpression =
        TdbExpression(terms, functions)

    override fun toString() = terms.joinToString(" + ")

    // -------------------------------------------------------------------------
    // Term types
    // -------------------------------------------------------------------------

    sealed class Term {
        abstract fun evaluate(t: Double, functions: Map<String, TdbFunction>): Double
    }

    /** Constant: A */
    data class Constant(val value: Double) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>) = value
    }

    /** T^n term: coefficient * T^exponent */
    data class PowerTerm(val coefficient: Double, val exponent: Double) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>) =
            coefficient * t.pow(exponent)
    }

    /** T*LN(T) term: coefficient * T * ln(T) */
    data class TLnT(val coefficient: Double) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>) =
            coefficient * t * ln(t)
    }

    /** LN(T) term: coefficient * ln(T) */
    data class LnT(val coefficient: Double) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>) =
            coefficient * ln(t)
    }

    /** LOG(T) term — same as LN in TDB context */
    data class LogT(val coefficient: Double) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>) =
            coefficient * log10(t)
    }

    /** Reference to a named FUNCTION: +GHSERFE, -GHSERMG */
    data class FunctionRef(val name: String, val sign: Double = 1.0) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>): Double {
            val fn = functions[name] ?: return 0.0
            return sign * (fn.evaluate(dev.hivens.solyx.core.units.Kelvin(t)) ?: 0.0)
        }
    }

    /** Scaled function reference: 2*GHSERFE */
    data class ScaledFunctionRef(val scale: Double, val name: String, val sign: Double = 1.0) : Term() {
        override fun evaluate(t: Double, functions: Map<String, TdbFunction>): Double {
            val fn = functions[name] ?: return 0.0
            return sign * scale * (fn.evaluate(dev.hivens.solyx.core.units.Kelvin(t)) ?: 0.0)
        }
    }
}

/**
 * Parser for TDB mathematical expressions.
 *
 * Parses strings like:
 * ```
 * -7770.458+130.485*T-24.112*T*LN(T)-.00265684*T**2+52478*T**(-1)
 * +GHSERCU+600+0.2*T
 * 2*GHSERFE+GHSERVA
 * ```
 */
object TdbExpressionParser {

    private val FUNCTION_NAME = Regex("[A-Z][A-Z0-9_#]*")
    private val NUMBER = Regex("(?:\\d+(?:\\.\\d+)?|\\.\\d+)(?:[Ee][+-]?\\d+)?")

    fun parse(input: String): TdbExpression {
        val terms = mutableListOf<TdbExpression.Term>()
        val cleaned = input.trim().replace("\\s+".toRegex(), "")

        var i = 0
        while (i < cleaned.length) {
            if (cleaned[i] !in "+-0123456789.TGHSERLNABCDEFIJKMOPQUVWXYZabcdefghijklmnopqrstuvwxyz") {
                i++
                continue
            }

            // Determine sign
            var sign = 1.0
            when {
                cleaned[i] == '+' -> { sign = 1.0; i++ }
                cleaned[i] == '-' -> { sign = -1.0; i++ }
            }

            if (i >= cleaned.length) break

            // Check if starts with a letter — could be function ref
            if (cleaned[i].isLetter()) {
                val match = FUNCTION_NAME.matchAt(cleaned, i) ?: break
                val name = match.value
                i += name.length

                // Check for multiplication: FUNCTION*something (unusual but possible)
                if (i < cleaned.length && cleaned[i] == '*') {
                    // skip, treat as function ref
                }
                terms.add(TdbExpression.FunctionRef(name, sign))
                continue
            }

            // Parse number
            val numMatch = NUMBER.matchAt(cleaned, i)
            val coefficient = if (numMatch != null && numMatch.range.first == i) {
                i += numMatch.value.length
                numMatch.value.toDouble() * sign
            } else {
                sign
            }

            // What follows?
            if (i >= cleaned.length || cleaned[i] == '+' || cleaned[i] == '-') {
                terms.add(TdbExpression.Constant(coefficient))
                continue
            }

            if (cleaned[i] == '*') {
                i++ // consume '*'

                // Check for function name after '*'
                if (i < cleaned.length && cleaned[i].isLetter()) {
                    val match = FUNCTION_NAME.matchAt(cleaned, i) ?: break
                    val name = match.value
                    i += name.length
                    terms.add(TdbExpression.ScaledFunctionRef(coefficient, name))
                    continue
                }

                // Must be T
                if (i < cleaned.length && cleaned[i] == 'T') {
                    i++ // consume 'T'

                    if (i >= cleaned.length || cleaned[i] == '+' || cleaned[i] == '-') {
                        // coefficient * T
                        terms.add(TdbExpression.PowerTerm(coefficient, 1.0))
                        continue
                    }

                    if (cleaned[i] == '*') {
                        i++ // consume '*'

                        // T*LN(T) or T**n
                        if (i < cleaned.length && cleaned[i] == 'L') {
                            // T*LN(T) or T*LOG(T)
                            val isLog = cleaned.startsWith("LOG(T)", i)
                            i += if (isLog) 6 else 5
                            terms.add(TdbExpression.TLnT(coefficient))
                            continue
                        }

                        if (cleaned[i] == '*') {
                            i++ // consume second '*' of '**'
                            // Parse exponent — may be negative in parens: T**(-9)
                            val expStr = if (cleaned[i] == '(') {
                                val end = cleaned.indexOf(')', i)
                                val s = cleaned.substring(i + 1, end)
                                i = end + 1
                                s
                            } else {
                                val m = NUMBER.matchAt(cleaned, i)!!
                                i += m.value.length
                                m.value
                            }
                            terms.add(TdbExpression.PowerTerm(coefficient, expStr.toDouble()))
                            continue
                        }
                    }

                    // Just T
                    terms.add(TdbExpression.PowerTerm(coefficient, 1.0))
                    continue
                }
            }

            // Lone T
            if (i < cleaned.length && cleaned[i] == 'T') {
                i++
                terms.add(TdbExpression.PowerTerm(coefficient, 1.0))
                continue
            }

            // Fallback — just a constant
            terms.add(TdbExpression.Constant(coefficient))
            i++ // prevent infinite loop on unrecognized character
        }

        return TdbExpression(terms)
    }
}
