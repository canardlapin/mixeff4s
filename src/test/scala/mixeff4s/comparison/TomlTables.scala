package mixeff4s.comparison

/** Minimal TOML reader for scorecard / fixture catalogs. Not a general parser. */
object TomlTables:
  final case class Document(root: Map[String, Value], rows: Vector[Map[String, Value]])

  enum Value:
    case Str(value: String)
    case Num(value: Double)
    case Arr(values: Vector[Value])

    def asString: String =
      this match
        case Str(value) => value
        case other      => throw IllegalArgumentException(s"expected string, got $other")

    def asStrings: Vector[String] =
      this match
        case Arr(values) => values.map(_.asString)
        case other       => throw IllegalArgumentException(s"expected string array, got $other")

    def asDouble: Double =
      this match
        case Num(value) => value
        case other      => throw IllegalArgumentException(s"expected number, got $other")

  def parse(text: String): Document =
    val root = scala.collection.mutable.LinkedHashMap.empty[String, Value]
    val rows = Vector.newBuilder[scala.collection.mutable.LinkedHashMap[String, Value]]
    var current: scala.collection.mutable.LinkedHashMap[String, Value] = root
    val logical = coalesce(text)
    logical.foreach: (line, idx) =>
      if line == "[[row]]" then
        val next = scala.collection.mutable.LinkedHashMap.empty[String, Value]
        rows += next
        current = next
      else if line.startsWith("[") && line.endsWith("]") && !line.startsWith("[[") then current = root
      else
        line.split("=", 2) match
          case Array(k, v) => current.update(k.trim, parseValue(v.trim, idx))
          case _           =>
            throw IllegalArgumentException(s"line $idx: expected key = value")
    Document(root.toMap, rows.result().map(_.toMap))

  private def coalesce(text: String): Vector[(String, Int)] =
    val out = Vector.newBuilder[(String, Int)]
    val buf = new StringBuilder
    var start = 0
    var depth = 0
    text.linesIterator.zipWithIndex.foreach: (raw, idx) =>
      val line = stripComment(raw).trim
      if line.nonEmpty then
        if buf.isEmpty then start = idx + 1
        if buf.nonEmpty then buf += ' '
        buf ++= line
        depth += line.count(_ == '[') - line.count(_ == ']')
        if depth < 0 then throw IllegalArgumentException(s"line ${idx + 1}: unmatched `]`")
        if depth == 0 then
          out += ((buf.toString, start))
          buf.clear()
    if buf.nonEmpty then throw IllegalArgumentException("unterminated TOML array")
    out.result()

  private def stripComment(line: String): String =
    val out = new StringBuilder
    var inString = false
    var i = 0
    while i < line.length do
      val c = line.charAt(i)
      if c == '"' && (i == 0 || line.charAt(i - 1) != '\\') then inString = !inString
      if c == '#' && !inString then i = line.length
      else
        out += c
        i += 1
    out.toString

  private def parseValue(raw: String, line: Int): Value =
    if raw.startsWith("[") && raw.endsWith("]") then
      val inner = raw.substring(1, raw.length - 1).trim
      if inner.isEmpty then Value.Arr(Vector.empty)
      else Value.Arr(splitTop(inner).map(part => parseValue(part.trim, line)))
    else if raw.startsWith("\"") && raw.endsWith("\"") && raw.length >= 2 then
      Value.Str(unescape(raw.substring(1, raw.length - 1)))
    else
      raw.toDoubleOption match
        case Some(n) => Value.Num(n)
        case None    => throw IllegalArgumentException(s"line $line: unsupported value `$raw`")

  private def splitTop(inner: String): Vector[String] =
    val parts = Vector.newBuilder[String]
    val buf = new StringBuilder
    var inString = false
    inner.foreach: c =>
      if c == '"' then inString = !inString
      if c == ',' && !inString then
        parts += buf.toString
        buf.clear()
      else buf += c
    if buf.nonEmpty then parts += buf.toString
    parts.result()

  private def unescape(value: String): String =
    value.replace("\\\"", "\"").replace("\\\\", "\\")
