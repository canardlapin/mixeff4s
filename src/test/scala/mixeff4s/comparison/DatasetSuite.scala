package mixeff4s.comparison

class DatasetSuite extends munit.FunSuite:
  test("catalog names match the fixture manifest"):
    assertEquals(Datasets.names, FixtureCatalog.loadEmbedded.available)

  test("every catalog dataset loads with the declared row count"):
    Datasets.loadAll.foreach: dataset =>
      assertEquals(dataset.frame.nRows, dataset.meta.nRows, clues(dataset.name))

  test("rail keeps the nlme Rail level order"):
    val rail = Datasets.load("rail")
    assertEquals(rail.frame.factor("Rail").map(_.levels), Some(Vector("2", "5", "1", "6", "3", "4")))

  test("ergostool keeps the nlme Subject level order"):
    val stool = Datasets.load("ergostool")
    assertEquals(
      stool.frame.factor("Subject").map(_.levels),
      Some(Vector("8", "5", "4", "9", "6", "3", "7", "1", "2"))
    )

  test("oats Variety includes a spaced level"):
    val oats = Datasets.load("oats")
    assertEquals(
      oats.frame.factor("Variety").map(_.levels),
      Some(Vector("Golden Rain", "Marvellous", "Victory"))
    )
    assertEquals(oats.frame.nRows, 72)
