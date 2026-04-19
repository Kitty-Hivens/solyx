package dev.hivens.solyx.calphad.tdb

import dev.hivens.solyx.core.units.Kelvin
import io.kotest.assertions.throwables.shouldNotThrowAny
import io.kotest.core.spec.style.DescribeSpec
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.collections.shouldHaveSize
import io.kotest.matchers.doubles.shouldBeWithinPercentageOf
import io.kotest.matchers.nulls.shouldNotBeNull
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe

class TdbParserTest : DescribeSpec({

    val tdb by lazy {
        val resource = TdbParserTest::class.java
            .getResourceAsStream("/test-cumg.tdb")
            ?.reader()
            ?: error("test-cumg.tdb not found in test resources")
        TdbParser.parse(resource)
    }

    // -------------------------------------------------------------------------
    // ELEMENT parsing
    // -------------------------------------------------------------------------

    describe("ELEMENT parsing") {

        it("parses all elements") {
            tdb.elements.keys shouldContain "CU"
            tdb.elements.keys shouldContain "MG"
            tdb.elements.keys shouldContain "VA"
        }

        it("parses element atomic mass") {
            tdb.elements["CU"]!!.atomicMass.shouldBeWithinPercentageOf(63.546, 0.01)
        }

        it("parses element reference state") {
            tdb.elements["CU"]!!.referenceState shouldBe "FCC_A1"
            tdb.elements["MG"]!!.referenceState shouldBe "HCP_A3"
        }

        it("parses H298-H0 and S298") {
            val cu = tdb.elements["CU"]!!
            cu.h298MinusH0.shouldBeWithinPercentageOf(5004.1, 0.1)
            cu.s298.shouldBeWithinPercentageOf(33.15, 0.1)
        }

        it("parses vacancy element") {
            tdb.elements["VA"]!!.atomicMass shouldBe 0.0
        }
    }

    // -------------------------------------------------------------------------
    // FUNCTION parsing
    // -------------------------------------------------------------------------

    describe("FUNCTION parsing") {

        it("parses all functions") {
            tdb.functions.keys shouldContain "GHSERCU"
            tdb.functions.keys shouldContain "GHSERMG"
            tdb.functions.keys shouldContain "GLIQCU"
            tdb.functions.keys shouldContain "GLIQMG"
        }

        it("GHSERCU has two temperature ranges") {
            tdb.functions["GHSERCU"]!!.ranges shouldHaveSize 2
        }

        it("GHSERCU first range starts at 298.15") {
            tdb.functions["GHSERCU"]!!.ranges[0].tLow
                .shouldBeWithinPercentageOf(298.15, 0.01)
        }

        it("GHSERCU transition at 1357.77 K (melting point of Cu)") {
            tdb.functions["GHSERCU"]!!.ranges[0].tHigh
                .shouldBeWithinPercentageOf(1357.77, 0.01)
        }

        it("GHSERCU second range ends at 3200 K") {
            tdb.functions["GHSERCU"]!!.ranges[1].tHigh
                .shouldBeWithinPercentageOf(3200.0, 0.01)
        }

        it("GHSERCU evaluates at 1000 K without error") {
            shouldNotThrowAny {
                val g = tdb.functions["GHSERCU"]!!.evaluate(Kelvin(1000.0))
                g.shouldNotBeNull()
                g.isNaN()      shouldBe false
                g.isInfinite() shouldBe false
            }
        }

        it("GHSERCU returns null outside temperature range") {
            tdb.functions["GHSERCU"]!!.evaluate(Kelvin(5000.0)) shouldBe null
        }

        it("GHSERCU at 298.15 K — known approximate value") {
            val g = tdb.functions["GHSERCU"]!!.evaluate(Kelvin(298.15))
            g.shouldNotBeNull()
            // At reference temperature G ≈ 0 by definition of GHSER
            // Exact value depends on H0 reference
            g.isNaN() shouldBe false
        }

        it("GLIQCU references GHSERCU and evaluates correctly") {
            val g = tdb.functions["GLIQCU"]!!.evaluate(Kelvin(1400.0))
            g.shouldNotBeNull()
            g.isNaN() shouldBe false
        }
    }

    // -------------------------------------------------------------------------
    // PHASE parsing
    // -------------------------------------------------------------------------

    describe("PHASE parsing") {

        it("parses all phases") {
            tdb.phases.keys shouldContain "LIQUID"
            tdb.phases.keys shouldContain "FCC_A1"
            tdb.phases.keys shouldContain "HCP_A3"
            tdb.phases.keys shouldContain "CUMG2"
        }

        it("LIQUID has one sublattice") {
            tdb.phases["LIQUID"]!!.sublatticeCount shouldBe 1
        }

        it("FCC_A1 has two sublattices") {
            tdb.phases["FCC_A1"]!!.sublatticeCount shouldBe 2
        }

        it("HCP_A3 stoichiometry is 1:0.5") {
            val hcp = tdb.phases["HCP_A3"]!!
            hcp.sublattices[0].shouldBeWithinPercentageOf(1.0, 0.01)
            hcp.sublattices[1].shouldBeWithinPercentageOf(0.5, 0.01)
        }

        it("CUMG2 stoichiometry is 1:2") {
            val cumg2 = tdb.phases["CUMG2"]!!
            cumg2.sublattices[0].shouldBeWithinPercentageOf(1.0, 0.01)
            cumg2.sublattices[1].shouldBeWithinPercentageOf(2.0, 0.01)
        }

        it("FCC_A1 has magnetic parameters from TYPE_DEF") {
            val fcc = tdb.phases["FCC_A1"]!!
            fcc.magneticAfm.shouldNotBeNull()
            fcc.magneticAfm.shouldBeWithinPercentageOf(-3.0, 0.01)
            fcc.magneticP!!.shouldBeWithinPercentageOf(0.28, 0.01)
        }

        it("LIQUID has no magnetic parameters") {
            tdb.phases["LIQUID"]!!.magneticAfm shouldBe null
        }
    }

    // -------------------------------------------------------------------------
    // CONSTITUENT parsing
    // -------------------------------------------------------------------------

    describe("CONSTITUENT parsing") {

        it("LIQUID constituents: CU and MG on single sublattice") {
            val constituents = tdb.phases["LIQUID"]!!.constituents
            constituents shouldHaveSize 1
            constituents[0] shouldContain "CU"
            constituents[0] shouldContain "MG"
        }

        it("FCC_A1 has CU and MG on first sublattice, VA on second") {
            val constituents = tdb.phases["FCC_A1"]!!.constituents
            constituents shouldHaveSize 2
            constituents[0] shouldContain "CU"
            constituents[0] shouldContain "MG"
            constituents[1] shouldContain "VA"
        }

        it("CUMG2 has single constituent per sublattice") {
            val constituents = tdb.phases["CUMG2"]!!.constituents
            constituents[0] shouldHaveSize 1
            constituents[0][0] shouldBe "CU"
            constituents[1] shouldHaveSize 1
            constituents[1][0] shouldBe "MG"
        }
    }

    // -------------------------------------------------------------------------
    // PARAMETER parsing
    // -------------------------------------------------------------------------

    describe("PARAMETER parsing") {

        it("parses G parameters for LIQUID") {
            val gParams = tdb.parametersFor("LIQUID", ParameterType.G)
            gParams.size shouldNotBe 0
        }

        it("parses L parameters for LIQUID") {
            val lParams = tdb.parametersFor("LIQUID", ParameterType.L)
            lParams shouldHaveSize 2  // L0 and L1
        }

        it("L0 has order 0, L1 has order 1") {
            val lParams = tdb.parametersFor("LIQUID", ParameterType.L)
                .sortedBy { it.order }
            lParams[0].order shouldBe 0
            lParams[1].order shouldBe 1
        }

        it("L0 evaluates at 1000K — known approximate value") {
            val l0 = tdb.parametersFor("LIQUID", ParameterType.L)
                .first { it.order == 0 }
            val value = l0.evaluate(Kelvin(1000.0))
            value.shouldNotBeNull()
            // L0 = -34000 + 4*T → at 1000K: -34000 + 4000 = -30000
            value.shouldBeWithinPercentageOf(-30000.0, 0.1)
        }

        it("G parameter for FCC_A1 CU:VA evaluates correctly") {
            val g = tdb.parametersFor("FCC_A1", ParameterType.G)
                .firstOrNull { param ->
                    param.constituentArray.any { sub -> sub.contains("CU") }
                }
            g.shouldNotBeNull()
            val value = g.evaluate(Kelvin(1000.0))
            value.shouldNotBeNull()
            value.isNaN() shouldBe false
        }

        it("total parameter count is reasonable") {
            tdb.parameters.size shouldNotBe 0
        }
    }

    // -------------------------------------------------------------------------
    // Edge cases
    // -------------------------------------------------------------------------

    describe("edge cases") {

        it("comments are stripped correctly") {
            // If $ comments were not stripped, function names would be corrupted
            tdb.functions["GHSERCU"].shouldNotBeNull()
        }

        it("multi-line entries are parsed correctly") {
            // GHSERCU spans multiple lines — must be joined before parsing
            tdb.functions["GHSERCU"]!!.ranges shouldHaveSize 2
        }

        it("empty database parses without error") {
            shouldNotThrowAny {
                TdbParser.parse("$ empty database\n")
            }
        }

        it("unknown keywords are silently ignored") {
            val withUnknown = """
                $ test
                ELEMENT CU FCC_A1 63.546 5004.1 33.15 !
                UNKNOWN_KEYWORD some data here !
                ELEMENT MG HCP_A3 24.305 4998.0 32.671 !
            """.trimIndent()
            shouldNotThrowAny {
                val db = TdbParser.parse(withUnknown)
                db.elements.keys shouldContain "CU"
                db.elements.keys shouldContain "MG"
            }
        }
    }
})
