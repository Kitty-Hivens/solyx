plugins {
    alias(libs.plugins.kotlin.jvm)
}

description = "Thermodynamic models for Solyx"

dependencies {
    api(project(":solyx-core"))
}
