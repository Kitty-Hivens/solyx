package dev.hivens.solyx.core.units

/**
 * Pressure in Pascals.
 *
 * Most CALPHAD calculations assume standard pressure [STANDARD].
 * Explicit pressure is required for high-pressure phase diagrams.
 *
 * ```kotlin
 * val p = Pressure.STANDARD
 * val p = Pressure.fromBar(1.0)
 * ```
 */
@JvmInline
value class Pressure(val value: Double) : Comparable<Pressure> {

    init {
        require(value >= 0.0) { "Pressure cannot be negative: $value Pa" }
    }

    operator fun plus(other: Pressure)  = Pressure(value + other.value)
    operator fun minus(other: Pressure) = Pressure(value - other.value)
    operator fun times(factor: Double)  = Pressure(value * factor)

    override fun compareTo(other: Pressure) = value.compareTo(other.value)

    fun toBar() = value / 100_000.0
    fun toAtm() = value / 101_325.0

    override fun toString() = "${value} Pa"

    companion object {
        /** 101 325 Pa — standard atmospheric pressure used as default in CALPHAD. */
        val STANDARD = Pressure(101_325.0)

        fun fromBar(bar: Double)   = Pressure(bar * 100_000.0)
        fun fromAtm(atm: Double)   = Pressure(atm * 101_325.0)
        fun fromKPa(kPa: Double)   = Pressure(kPa * 1_000.0)
    }
}

/**
 * Mole fraction — dimensionless, range [0.0, 1.0].
 *
 * Represents the fraction of a component in a phase or system.
 * Enforced at construction time — physically invalid values are rejected.
 *
 * ```kotlin
 * val x = MoleFraction(0.02)
 * val x = 0.02.moleFraction
 * ```
 */
@JvmInline
value class MoleFraction(val value: Double) : Comparable<MoleFraction> {

    init {
        require(value in 0.0..1.0) { "Mole fraction must be in [0, 1]: $value" }
    }

    operator fun plus(other: MoleFraction)  = MoleFraction(value + other.value)
    operator fun minus(other: MoleFraction) = MoleFraction(value - other.value)
    operator fun times(factor: Double)      = MoleFraction(value * factor)

    val complement get() = MoleFraction(1.0 - value)

    override fun compareTo(other: MoleFraction) = value.compareTo(other.value)

    override fun toString() = "$value"
}

/** Convenience extensions. */
val Double.pascal      get() = Pressure(this)
val Double.moleFraction get() = MoleFraction(this)
