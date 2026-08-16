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

lazy val mixeff4s = project
  .in(file("."))
  .settings(
    name := "mixeff4s",
    description := "Linear and generalized linear mixed-effects models for Scala 3.",
    scalacOptions ++= Seq(
      "-deprecation",
      "-feature",
      "-unchecked",
      "-Wunused:all",
      "-Wvalue-discard"
    ),
    libraryDependencies += "org.scalameta" %% "munit" % "1.2.1" % Test,
    Test / fork := false
  )
