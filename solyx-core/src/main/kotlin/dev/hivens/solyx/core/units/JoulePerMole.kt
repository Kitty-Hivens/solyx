package dev.hivens.solyx.core.units

/**
 * Energy in Joules per mole (J/mol) — standard unit for molar Gibbs energy.
 *
 * ```kotlin
 * val g = JoulePerMole(-8040.0)
 * val g = (-8040.0).jPerMol
 * ```
 */
@JvmInline
value class JoulePerMole(val value: Double) : Comparable<JoulePerMole> {

    operator fun plus(other: JoulePerMole)  = JoulePerMole(value + other.value)
    operator fun minus(other: JoulePerMole) = JoulePerMole(value - other.value)
    operator fun times(factor: Double)      = JoulePerMole(value * factor)
    operator fun div(factor: Double)        = JoulePerMole(value / factor)
    operator fun unaryMinus()               = JoulePerMole(-value)

    override fun compareTo(other: JoulePerMole) = value.compareTo(other.value)

    fun toKJPerMole() = value / 1000.0

    override fun toString() = "$value J/mol"
}

/**
 * Entropy in Joules per mole per Kelvin (J/mol·K).
 *
 * ```kotlin
 * val s = JoulePerMoleKelvin(2.703)
 * ```
 */
@JvmInline
value class JoulePerMoleKelvin(val value: Double) : Comparable<JoulePerMoleKelvin> {

    operator fun plus(other: JoulePerMoleKelvin)  = JoulePerMoleKelvin(value + other.value)
    operator fun minus(other: JoulePerMoleKelvin) = JoulePerMoleKelvin(value - other.value)
    operator fun times(factor: Double)            = JoulePerMoleKelvin(value * factor)
    operator fun unaryMinus()                     = JoulePerMoleKelvin(-value)

    override fun compareTo(other: JoulePerMoleKelvin) = value.compareTo(other.value)

    /** S * T -> G contribution */
    operator fun times(temperature: Kelvin) = JoulePerMole(value * temperature.value)

    override fun toString() = "$value J/(mol·K)"
}

/** Convenience extensions. */
val Double.jPerMol    get() = JoulePerMole(this)
val Double.jPerMolK   get() = JoulePerMoleKelvin(this)
val Double.kjPerMol   get() = JoulePerMole(this * 1000.0)
