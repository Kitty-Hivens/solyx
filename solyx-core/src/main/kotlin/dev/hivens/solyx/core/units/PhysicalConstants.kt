package dev.hivens.solyx.core.units

/**
 * Fundamental physical constants used in thermodynamic calculations.
 *
 * All values are in SI units (CODATA 2018 recommended values).
 */
object PhysicalConstants {

    /**
     * Universal gas constant: 8.314462618 J/(mol·K)
     *
     * Appears in every entropy of mixing calculation:
     * `S_mix = -R * Σ(x_i * ln(x_i))`
     */
    const val R: Double = 8.314462618

    /**
     * Avogadro constant: 6.02214076 × 10²³ mol⁻¹
     */
    const val N_A: Double = 6.02214076e23

    /**
     * Boltzmann constant: 1.380649 × 10⁻²³ J/K
     *
     * k_B = R / N_A
     */
    const val K_B: Double = 1.380649e-23

    /**
     * Faraday constant: 96485.33212 C/mol
     * Required for electrochemical calculations.
     */
    const val F: Double = 96485.33212
}
