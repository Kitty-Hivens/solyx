plugins {
    alias(libs.plugins.kotlin.jvm)
}

description = "CALPHAD engine for Solyx"

dependencies {
    api(project(":solyx-thermo"))
}
