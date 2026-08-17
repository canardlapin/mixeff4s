package mixeff4s.comparison

object Csv:
  def parse(text: String): (Vector[String], Vector[Vector[String]]) =
    val rows = text
      .replace("\r\n", "\n")
      .replace("\r", "\n")
      .split("\n", -1)
      .toVector
      .map(_.trim)
      .filter(_.nonEmpty)
      .map(parseRow)
    if rows.isEmpty then throw IllegalArgumentException("CSV is empty")
    (rows.head, rows.tail)

  private def parseRow(line: String): Vector[String] =
    val fields = Vector.newBuilder[String]
    val buf = new StringBuilder
    var inQuotes = false
    var i = 0
    while i < line.length do
      val c = line.charAt(i)
      if inQuotes then
        if c == '"' then
          if i + 1 < line.length && line.charAt(i + 1) == '"' then
            buf += '"'
            i += 1
          else inQuotes = false
        else buf += c
      else if c == '"' then inQuotes = true
      else if c == ',' then
        fields += buf.toString
        buf.clear()
      else buf += c
      i += 1
    if inQuotes then throw IllegalArgumentException(s"unterminated CSV quote in `$line`")
    fields += buf.toString
    fields.result()
