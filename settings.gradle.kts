rootProject.name = "solyx"

dependencyResolutionManagement {
    repositories {
        mavenCentral()
    }
}

include(
    ":solyx-core",
    ":solyx-thermo",
    ":solyx-calphad"
)
