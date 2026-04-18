package dev.hivens.solyx.core.units

/**
 * Temperature in Kelvin — the only unit used internally.
 *
 * All thermodynamic calculations require Kelvin. Input in other units
 * must be converted at the boundary via [fromCelsius] or [fromFahrenheit].
 *
 * ```kotlin
 * val t = Kelvin(1200.0)
 * val t = Kelvin.fromCelsius(926.85)
 * ```
 */
@JvmInline
value class Kelvin(val value: Double) : Comparable<Kelvin> {

    init {
        require(value >= 0.0) { "Temperature cannot be below absolute zero: $value K" }
    }

    operator fun plus(other: Kelvin)  = Kelvin(value + other.value)
    operator fun minus(other: Kelvin) = Kelvin(value - other.value)
    operator fun times(factor: Double) = Kelvin(value * factor)
    operator fun div(factor: Double)   = Kelvin(value / factor)

    override fun compareTo(other: Kelvin) = value.compareTo(other.value)

    fun toCelsius()    = value - 273.15
    fun toFahrenheit() = (value - 273.15) * 9.0 / 5.0 + 32.0

    override fun toString() = "$value K"

    companion object {
        val ABSOLUTE_ZERO = Kelvin(0.0)

        fun fromCelsius(celsius: Double)       = Kelvin(celsius + 273.15)
        fun fromFahrenheit(fahrenheit: Double) = Kelvin((fahrenheit - 32.0) * 5.0 / 9.0 + 273.15)
    }
}

/** Convenience extensions for readable DSL usage. */
val Double.kelvin  get() = Kelvin(this)
val Double.celsius get() = Kelvin.fromCelsius(this)
val Int.kelvin     get() = Kelvin(this.toDouble())
val Int.celsius    get() = Kelvin.fromCelsius(this.toDouble())
