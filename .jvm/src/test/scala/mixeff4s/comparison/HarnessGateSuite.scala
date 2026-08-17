package mixeff4s.comparison

import java.nio.file.{Files, Path}

/** JVM-only: the on-disk comparison tree must match the embedded catalog. */
class HarnessGateSuite extends munit.FunSuite:
  test("scorecard, fixtures, and frozen files match the embedded catalogs"):
    val root = comparisonRoot
    val scorecardDisk = Scorecard.parse(Files.readString(root.resolve("parity_scorecard.toml")))
    val scorecard = Scorecard.loadEmbedded
    assertEquals(scorecardDisk.rows.map(_.key), scorecard.rows.map(_.key))
    assertEquals(scorecardDisk.rows.map(_.classification), scorecard.rows.map(_.classification))
    assertEquals(scorecardDisk.rows.map(_.reference), scorecard.rows.map(_.reference))

    val fixturesDisk = FixtureCatalog.parse(Files.readString(root.resolve("fixtures.toml")))
    val fixtures = FixtureCatalog.loadEmbedded
    assertEquals(fixturesDisk, fixtures)

    val frozenDisk = FrozenCatalog.parse(Files.readString(root.resolve("frozen/references.json")))
    val frozen = FrozenCatalog.loadEmbedded
    assertEquals(frozenDisk.results.map(r => (r.key, r.engine)), frozen.results.map(r => (r.key, r.engine)))
    assertEquals(frozen.mixeffRsRevision, Some(fixtures.mixeffRsRevision))

  test("scorecard datasets exist and vendored directories match the catalog"):
    val root = comparisonRoot
    val fixtures = FixtureCatalog.loadEmbedded
    val scorecard = Scorecard.loadEmbedded
    val missing = scorecard.rows.map(_.key.dataset).distinct.filterNot(fixtures.available.contains)
    assertEquals(missing, Vector.empty)

    val onDisk =
      val stream = Files.list(root.resolve("datasets"))
      try
        stream.toArray.toVector
          .map(_.asInstanceOf[Path])
          .filter(Files.isDirectory(_))
          .map(_.getFileName.toString)
          .sorted
      finally stream.close()
    assertEquals(onDisk, fixtures.vendored.sorted)

    fixtures.vendored.foreach: name =>
      val csv = Files.readString(root.resolve(s"datasets/$name/data.csv"))
      assertEquals(csv, EmbeddedFrames.csv(name), clues(name))
      val meta = Files.readString(root.resolve(s"datasets/$name/meta.json"))
      assertEquals(Datasets.parseMeta(meta), Datasets.parseMeta(EmbeddedFrames.meta(name)), clues(name))

  test("every scorecard row names a loadable fixture and a frozen claim"):
    val scorecard = Scorecard.loadEmbedded
    val frozen = FrozenCatalog.loadEmbedded
    scorecard.rows.foreach: row =>
      val dataset = Datasets.load(row.key.dataset)
      assert(dataset.frame.nRows > 0, clues(row.key.id))
      assert(frozen.claimed(row).isDefined, clues(row.key.id, row.reference))
      if row.key.estimator == "fast_pirls" then assertEquals(row.classification, ScorecardClass.DocumentedDivergence)

  private def comparisonRoot: Path =
    val here = Path.of("").toAbsolutePath.normalize()
    val candidates = Vector(here.resolve("comparison"), here.getParent.resolve("comparison"))
    candidates
      .find(path => Files.exists(path.resolve("parity_scorecard.toml")))
      .getOrElse(fail(s"comparison/ not found from $here"))
