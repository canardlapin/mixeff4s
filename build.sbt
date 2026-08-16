import org.scalajs.linker.interface.ModuleKind
import sbtcrossproject.CrossPlugin.autoImport.{crossProject, CrossType}
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

ThisBuild / organization := "io.github.canardlapin"
ThisBuild / scalaVersion := "3.7.4"
ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / versionScheme := Some("early-semver")
ThisBuild / licenses := Seq("Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt"))
ThisBuild / homepage := Some(url("https://github.com/canardlapin/mixeff4s"))
ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/canardlapin/mixeff4s"),
    "scm:git:https://github.com/canardlapin/mixeff4s.git",
    Some("scm:git:git@github.com:canardlapin/mixeff4s.git")
  )
)
ThisBuild / developers := List(
  Developer(
    "canardlapin",
    "canardlapin",
    "307091466+canardlapin@users.noreply.github.com",
    url("https://github.com/canardlapin")
  )
)

// Gale is source-only. Pin an immutable revision; a local checkout is admitted
// only through -Dmixeff4s.gale.build=/path/to/gale.
lazy val galeRevision = "f869613cec0a89e57b6c995b0a02cf471ac7127c"
lazy val galeBuild: java.net.URI =
  sys.props
    .get("mixeff4s.gale.build")
    .map(path => file(path).getCanonicalFile.toURI)
    .getOrElse(uri(s"https://github.com/canardlapin/gale.git#$galeRevision"))
lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")
lazy val galeCoreJS = ProjectRef(galeBuild, "coreJS")

lazy val commonSettings = Seq(
  name := "mixeff4s",
  description := "Linear and generalized linear mixed-effects models for Scala 3.",
  scalacOptions ++= Seq(
    "-deprecation",
    "-feature",
    "-unchecked",
    "-Wunused:all",
    "-Wvalue-discard"
  ),
  libraryDependencies += "org.scalameta" %%% "munit" % "1.2.1" % Test,
  Test / fork := false
)

lazy val mixeff4s =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Pure)
    .in(file("."))
    .settings(commonSettings)
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(
      scalaJSLinkerConfig ~= (_.withModuleKind(ModuleKind.CommonJSModule))
    )

lazy val mixeff4sJVM = mixeff4s.jvm
lazy val mixeff4sJS = mixeff4s.js

lazy val root = project
  .in(file("."))
  .aggregate(mixeff4sJVM, mixeff4sJS)
  .settings(
    name := "mixeff4s-root",
    publish / skip := true,
    Compile / sources := Nil,
    Test / sources := Nil
  )
