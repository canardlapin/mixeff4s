package mixeff4s

import java.nio.file.{Files, Path}

class ArchitectureSuite extends munit.FunSuite:
  test("layer imports stay acyclic"):
    val root = Path.of("src/main/scala/mixeff4s")
    assume(Files.exists(root), "source tree is present")
    val violations = Vector.newBuilder[String]
    Files
      .walk(root)
      .forEach: path =>
        if path.toString.endsWith(".scala") then
          val layer = layerOf(root.relativize(path).toString)
          val text = Files.readString(path)
          forbidden(layer).foreach: banned =>
            if text.contains(s"import $banned") || text.contains(s"import $banned.") then
              violations += s"${path.getFileName}: $layer must not import $banned"
    val found = violations.result()
    assert(found.isEmpty, clues(found))

  private def layerOf(relative: String): String =
    relative.split("/").headOption.getOrElse("unknown")

  private def forbidden(layer: String): Set[String] =
    layer match
      case "error"   => Set("mixeff4s.data", "mixeff4s.lmm", "mixeff4s.model", "mixeff4s.design", "mixeff4s.glmm")
      case "formula" => Set("mixeff4s.lmm", "mixeff4s.model", "mixeff4s.design", "mixeff4s.glmm")
      case "data"    => Set("mixeff4s.formula", "mixeff4s.lmm", "mixeff4s.model", "mixeff4s.design", "mixeff4s.glmm")
      case "model"   => Set("mixeff4s.lmm", "mixeff4s.formula", "mixeff4s.data", "mixeff4s.design", "mixeff4s.glmm")
      case "design"  => Set("mixeff4s.lmm", "mixeff4s.glmm", "mixeff4s.compiler", "mixeff4s.stats")
      case "linalg"  =>
        Set(
          "mixeff4s.lmm",
          "mixeff4s.formula",
          "mixeff4s.data",
          "mixeff4s.model",
          "mixeff4s.design",
          "mixeff4s.optimizer",
          "mixeff4s.glmm"
        )
      case "optimizer" =>
        Set("mixeff4s.lmm", "mixeff4s.design", "mixeff4s.formula", "mixeff4s.data", "mixeff4s.model", "mixeff4s.glmm")
      case "lmm"  => Set("mixeff4s.compiler", "mixeff4s.stats", "mixeff4s.glmm")
      case "glmm" => Set("mixeff4s.compiler", "mixeff4s.stats")
      case _      => Set.empty
