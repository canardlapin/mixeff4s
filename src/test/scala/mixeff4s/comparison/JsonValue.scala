package mixeff4s.comparison

/** Small JSON reader for frozen references and dataset meta. Not a general encoder. */
enum JsonValue:
  case Obj(fields: Vector[(String, JsonValue)])
  case Arr(values: Vector[JsonValue])
  case Str(value: String)
  case Num(value: Double)
  case Bool(value: Boolean)
  case Null

  def field(name: String): Option[JsonValue] =
    this match
      case Obj(fields) => fields.collectFirst { case (k, v) if k == name => v }
      case _           => None

  def req(name: String): JsonValue =
    field(name).getOrElse(throw IllegalArgumentException(s"missing JSON field `$name`"))

  def asString: String =
    this match
      case Str(value) => value
      case other      => throw IllegalArgumentException(s"expected string, got $other")

  def asDouble: Double =
    this match
      case Num(value) => value
      case other      => throw IllegalArgumentException(s"expected number, got $other")

  def asInt: Int =
    val n = asDouble
    if n != n.floor then throw IllegalArgumentException(s"expected integer, got $n")
    n.toInt

  def asBoolean: Boolean =
    this match
      case Bool(value) => value
      case other       => throw IllegalArgumentException(s"expected boolean, got $other")

  def asArray: Vector[JsonValue] =
    this match
      case Arr(values) => values
      case other       => throw IllegalArgumentException(s"expected array, got $other")

  def asObject: Vector[(String, JsonValue)] =
    this match
      case Obj(fields) => fields
      case other       => throw IllegalArgumentException(s"expected object, got $other")

object JsonValue:
  def parse(text: String): JsonValue =
    val parser = Parser(text)
    val value = parser.parseValue()
    parser.skipWs()
    if !parser.done then throw IllegalArgumentException("trailing JSON content")
    value

  private final class Parser(text: String):
    private var i = 0

    def done: Boolean = i >= text.length

    def skipWs(): Unit =
      while i < text.length && text.charAt(i).isWhitespace do i += 1

    def parseValue(): JsonValue =
      skipWs()
      if done then throw IllegalArgumentException("unexpected end of JSON")
      text.charAt(i) match
        case '{'                        => parseObject()
        case '['                        => parseArray()
        case '"'                        => JsonValue.Str(parseString())
        case 't'                        => keyword("true", JsonValue.Bool(true))
        case 'f'                        => keyword("false", JsonValue.Bool(false))
        case 'n'                        => keyword("null", JsonValue.Null)
        case c if c == '-' || c.isDigit => parseNumber()
        case c                          => throw IllegalArgumentException(s"unexpected JSON at $i: `$c`")

    private def parseObject(): JsonValue =
      i += 1
      val fields = Vector.newBuilder[(String, JsonValue)]
      skipWs()
      if peek('}') then
        i += 1
        JsonValue.Obj(fields.result())
      else
        var continue = true
        while continue do
          skipWs()
          val key = parseString()
          skipWs()
          expect(':')
          fields += (key -> parseValue())
          skipWs()
          if peek('}') then
            i += 1
            continue = false
          else expect(',')
        JsonValue.Obj(fields.result())

    private def parseArray(): JsonValue =
      i += 1
      val values = Vector.newBuilder[JsonValue]
      skipWs()
      if peek(']') then
        i += 1
        JsonValue.Arr(values.result())
      else
        var continue = true
        while continue do
          values += parseValue()
          skipWs()
          if peek(']') then
            i += 1
            continue = false
          else expect(',')
        JsonValue.Arr(values.result())

    private def parseString(): String =
      expect('"')
      val out = new StringBuilder
      var closed = false
      while !closed do
        if done then throw IllegalArgumentException("unterminated JSON string")
        val c = text.charAt(i)
        i += 1
        c match
          case '"'  => closed = true
          case '\\' =>
            if done then throw IllegalArgumentException("unterminated JSON escape")
            val e = text.charAt(i)
            i += 1
            e match
              case '"' | '\\' | '/' => out += e
              case 'b'              => out += '\b'
              case 'f'              => out += '\f'
              case 'n'              => out += '\n'
              case 'r'              => out += '\r'
              case 't'              => out += '\t'
              case 'u'              =>
                if i + 4 > text.length then throw IllegalArgumentException("bad unicode escape")
                out += Integer.parseInt(text.substring(i, i + 4), 16).toChar
                i += 4
              case other => throw IllegalArgumentException(s"bad escape `$other`")
          case other => out += other
      out.toString

    private def parseNumber(): JsonValue =
      val start = i
      if peek('-') then i += 1
      while i < text.length && text.charAt(i).isDigit do i += 1
      if peek('.') then
        i += 1
        while i < text.length && text.charAt(i).isDigit do i += 1
      if peek('e') || peek('E') then
        i += 1
        if peek('+') || peek('-') then i += 1
        while i < text.length && text.charAt(i).isDigit do i += 1
      JsonValue.Num(text.substring(start, i).toDouble)

    private def keyword(expected: String, value: JsonValue): JsonValue =
      if i + expected.length > text.length || text.substring(i, i + expected.length) != expected then
        throw IllegalArgumentException(s"expected `$expected` at $i")
      i += expected.length
      value

    private def expect(c: Char): Unit =
      skipWs()
      if done || text.charAt(i) != c then throw IllegalArgumentException(s"expected `$c` at $i")
      i += 1

    private def peek(c: Char): Boolean =
      i < text.length && text.charAt(i) == c
