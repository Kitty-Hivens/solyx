package dev.hivens.solyx.thermo

import dev.hivens.solyx.core.units.JoulePerMole
import dev.hivens.solyx.core.units.Kelvin
import dev.hivens.solyx.core.units.PhysicalConstants.R
import kotlin.math.ln
import kotlin.math.pow

/**
 * Magnetic ordering contribution to Gibbs energy — Inden-Hillert-Jarl model.
 *
 * Many transition metals (Fe, Ni, Co, Mn) exhibit magnetic ordering that
 * significantly affects their thermodynamic properties. Below the Curie
 * temperature [tc] the material is ferromagnetic; above it — paramagnetic.
 *
 * This contribution is added on top of the non-magnetic Gibbs energy:
 * ```
 * G_total = G_non_magnetic + G_magnetic
 * G_magnetic = R * T * ln(β + 1) * f(τ)
 * ```
 *
 * Where:
 * - τ = T / T_c — reduced temperature
 * - β — average magnetic moment per atom (Bohr magnetons)
 * - f(τ) — polynomial function, different below and above T_c
 *
 * Parameters for common phases (from CALPHAD databases):
 * - BCC Fe (ferrite): T_c = 1043 K, β = 2.22
 * - FCC Fe (austenite): T_c = 201 K,  β = -2.1  (antiferromagnetic)
 * - FCC Ni: T_c = 633 K,  β = 0.52
 *
 * Reference: Hillert & Jarl, CALPHAD 2 (1978) 227-238
 *
 * ```kotlin
 * val magneticBccFe = MagneticContribution(
 *     tc   = Kelvin(1043.0),
 *     beta = 2.22,
 *     p    = 0.28   // BCC structure constant
 * )
 * val gMag = magneticBccFe.compute(Kelvin(800.0))
 * ```
 *
 * @param tc   Curie temperature (ferromagnetic) or Néel temperature (antiferromagnetic)
 * @param beta average magnetic moment per atom in Bohr magnetons
 * @param p    structure-dependent constant — 0.28 for BCC, 0.25 for FCC/HCP
 */
class MagneticContribution(
    val tc: Kelvin,
    val beta: Double,
    val p: Double
) {
    init {
        require(tc.value > 0.0)   { "Curie temperature must be positive: ${tc.value} K" }
        require(p in 0.0..1.0)    { "Structure constant p must be in (0, 1): $p" }
    }

    /**
     * Compute magnetic Gibbs energy contribution at [temperature].
     *
     * Returns zero if [beta] is zero — no magnetic moment means no contribution.
     */
    fun compute(temperature: Kelvin): JoulePerMole {
        if (beta == 0.0) return JoulePerMole(0.0)

        val tau = temperature.value / tc.value
        val lnBeta1 = ln(beta.coerceAtLeast(0.0) + 1.0)
        val f = if (tau < 1.0) fBelow(tau) else fAbove(tau)

        return JoulePerMole(R * temperature.value * lnBeta1 * f)
    }

    /**
     * f(τ) below T_c — ferromagnetic region.
     *
     * ```
     * f(τ) = 1 - (1/A) * [79τ⁻¹/(140p) + (474/497)(1/p - 1)(τ³/6 + τ⁹/135 + τ¹⁵/600)]
     * ```
     */
    private fun fBelow(tau: Double): Double {
        val a = 518.0 / 1125.0 + 11692.0 / 15975.0 * (1.0 / p - 1.0)
        val term1 = 79.0 * tau.pow(-1) / (140.0 * p)
        val term2 = (474.0 / 497.0) * (1.0 / p - 1.0) *
                    (tau.pow(3) / 6.0 + tau.pow(9) / 135.0 + tau.pow(15) / 600.0)
        return 1.0 - (term1 + term2) / a
    }

    /**
     * f(τ) above T_c — paramagnetic region.
     *
     * ```
     * f(τ) = -(1/A) * [τ⁻⁵/10 + τ⁻¹⁵/315 + τ⁻²⁵/1500]
     * ```
     */
    private fun fAbove(tau: Double): Double {
        val a = 518.0 / 1125.0 + 11692.0 / 15975.0 * (1.0 / p - 1.0)
        return -(tau.pow(-5) / 10.0 + tau.pow(-15) / 315.0 + tau.pow(-25) / 1500.0) / a
    }

    companion object {
        /** BCC Fe (ferrite) — ferromagnetic below 1043 K. */
        val BCC_FE = MagneticContribution(tc = Kelvin(1043.0), beta = 2.22,  p = 0.28)

        /** FCC Fe (austenite) — antiferromagnetic, negative beta. */
        val FCC_FE = MagneticContribution(tc = Kelvin(201.0),  beta = -2.1,  p = 0.25)

        /** FCC Ni — ferromagnetic below 633 K. */
        val FCC_NI = MagneticContribution(tc = Kelvin(633.0),  beta = 0.52,  p = 0.25)

        /** FCC Co — ferromagnetic below 1388 K. */
        val FCC_CO = MagneticContribution(tc = Kelvin(1388.0), beta = 1.35,  p = 0.25)
    }
}
