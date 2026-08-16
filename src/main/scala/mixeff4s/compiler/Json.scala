package mixeff4s.compiler

/** Deterministic JSON for compiler snapshots. No third-party encoder. */
private[compiler] object Json:
  def str(value: String): String =
    val escaped = value.flatMap:
      case '"'  => "\\\""
      case '\\' => "\\\\"
      case '\n' => "\\n"
      case '\r' => "\\r"
      case '\t' => "\\t"
      case c    => c.toString
    s"\"$escaped\""

  def num(value: Int): String = value.toString

  def arr(values: Iterable[String]): String =
    values.mkString("[", ", ", "]")

  def obj(fields: (String, String)*): String =
    if fields.isEmpty then "{}"
    else
      fields
        .map((key, value) => s"${str(key)}: $value")
        .mkString("{", ", ", "}")

  def pretty(fields: Vector[(String, String)], indent: Int = 0): String =
    val pad = "  " * indent
    val inner = "  " * (indent + 1)
    if fields.isEmpty then "{}"
    else
      fields
        .map((key, value) => s"$inner${str(key)}: $value")
        .mkString("{\n", ",\n", s"\n$pad}")
