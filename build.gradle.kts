import org.jetbrains.kotlin.gradle.dsl.JvmTarget

plugins {
    alias(libs.plugins.kotlin.jvm) apply false
    alias(libs.plugins.vanniktech)  apply false
}

subprojects {
    apply(plugin = "org.jetbrains.kotlin.jvm")
    apply(plugin = "com.vanniktech.maven.publish")

    group   = "dev.hivens.solyx"
    version = "0.3.0-SNAPSHOT"

    extensions.configure<org.jetbrains.kotlin.gradle.dsl.KotlinJvmProjectExtension> {
        jvmToolchain(25)
        compilerOptions {
            jvmTarget.set(JvmTarget.JVM_25)
            allWarningsAsErrors.set(true)
            freeCompilerArgs.add("-Xcontext-parameters")
        }
    }

    dependencies {
        val libs = rootProject.extensions
            .getByType<VersionCatalogsExtension>()
            .named("libs")

        "testImplementation"(libs.findLibrary("kotest-runner").get())
        "testImplementation"(libs.findLibrary("kotest-assertions").get())
        "testImplementation"(libs.findLibrary("kotest-property").get())
        "testImplementation"(libs.findLibrary("junit-api").get())
        "testRuntimeOnly"(libs.findLibrary("junit-engine").get())
    }

    tasks.withType<Test> {
        useJUnitPlatform()
    }

    extensions.configure<com.vanniktech.maven.publish.MavenPublishBaseExtension> {
        publishToMavenCentral(com.vanniktech.maven.publish.SonatypeHost.CENTRAL_PORTAL)
        signAllPublications()

        pom {
            name.set(project.name)
            description.set("Computational thermodynamics library for the JVM")
            url.set("https://solyx.hivens.dev")
            inceptionYear.set("2026")

            licenses {
                license {
                    name.set("MIT License")
                    url.set("https://opensource.org/licenses/MIT")
                }
            }

            developers {
                developer {
                    id.set("hivens")
                    name.set("Haru")
                    url.set("https://hivens.dev")
                }
            }

            scm {
                url.set("https://github.com/Kitty-Hivens/solyx")
                connection.set("scm:git:git://github.com/Kitty-Hivens/solyx.git")
                developerConnection.set("scm:git:ssh://git@github.com/Kitty-Hivens/solyx.git")
            }
        }
    }
}
