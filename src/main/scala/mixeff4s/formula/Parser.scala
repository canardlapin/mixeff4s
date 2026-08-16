package mixeff4s.formula

/** Recursive-descent parser for lme4-style mixed-model formulas. */
object Parser:
  def parse(input: String): Either[FormulaError, Formula] =
    val trimmed = input.trim
    if trimmed.isEmpty then Left(FormulaError.Empty)
    else
      val derived = Vector.newBuilder[DerivedColumn]
      tokenize(trimmed, derived).flatMap: tokens =>
        if tokens.isEmpty then Left(FormulaError.Empty)
        else
          val parser = new FormulaParser(tokens, trimmed)
          for
            response <- parser.parseResponse()
            _ <- parser.expect(Token.Tilde)
            parsed <- parser.parseRhs()
            (rawFixed, random) = parsed
            _ <-
              if rawFixed.isEmpty && random.isEmpty then Left(FormulaError.EmptyRhs)
              else Right(())
          yield
            val fixed = canonicalizeFixed(rawFixed)
            val referenced = referencedNames(response, fixed, random)
            Formula(
              response,
              fixed,
              random,
              derived.result().filter(d => referenced.contains(d.label))
            )

  private def canonicalizeFixed(fixed: Vector[FixedTerm]): Vector[FixedTerm] =
    val hasExplicit =
      fixed.exists(t => t == FixedTerm.Intercept || t == FixedTerm.NoIntercept)
    val withIntercept =
      if hasExplicit then fixed else FixedTerm.Intercept +: fixed
    if withIntercept.exists(_ == FixedTerm.NoIntercept) then
      withIntercept.filterNot(t => t == FixedTerm.Intercept || t == FixedTerm.NoIntercept)
    else withIntercept

  private def referencedNames(
      response: String,
      fixed: Vector[FixedTerm],
      random: Vector[RandomTerm]
  ): Set[String] =
    val names = scala.collection.mutable.Set(response)
    fixed.foreach(t => names ++= t.columnNames)
    random.foreach(rt => rt.terms.foreach(t => names ++= t.columnNames))
    names.toSet

  private enum Token:
    case Tilde, Plus, Minus, Star, Colon, Slash, Pipe, DoublePipe, Ampersand, LParen, RParen, Zero, One
    case Ident(name: String)

  private final case class Spanned(token: Token, pos: Int)

  private def matchingParen(chars: Vector[Char], openIdx: Int): Either[FormulaError, Int] =
    var depth = 0
    var j = openIdx
    var inTick = false
    var found = -1
    while j < chars.length && found < 0 do
      chars(j) match
        case '`' =>
          inTick = !inTick
        case '(' if !inTick =>
          depth += 1
        case ')' if !inTick =>
          depth -= 1
          if depth == 0 then found = j + 1
        case _ =>
          ()
      j += 1
    if found >= 0 then Right(found)
    else if inTick then
      Left(
        FormulaError.Other(
          "unterminated backtick-quoted identifier inside an in-formula transform"
        )
      )
    else Left(FormulaError.UnmatchedParen)

  private def tokenize(
      input: String,
      derived: scala.collection.mutable.Builder[DerivedColumn, Vector[DerivedColumn]]
  ): Either[FormulaError, Vector[Spanned]] =
    val chars = input.toVector
    val tokens = Vector.newBuilder[Spanned]
    var i = 0
    var error: Option[FormulaError] = None
    val seen = scala.collection.mutable.Set.empty[String]
    while i < chars.length && error.isEmpty do
      val c = chars(i)
      if c.isWhitespace then i += 1
      else
        val pos = i
        c match
          case '~' =>
            tokens += Spanned(Token.Tilde, pos)
            i += 1
          case '+' =>
            tokens += Spanned(Token.Plus, pos)
            i += 1
          case '-' =>
            tokens += Spanned(Token.Minus, pos)
            i += 1
          case '*' =>
            tokens += Spanned(Token.Star, pos)
            i += 1
          case ':' =>
            tokens += Spanned(Token.Colon, pos)
            i += 1
          case '/' =>
            tokens += Spanned(Token.Slash, pos)
            i += 1
          case '&' =>
            tokens += Spanned(Token.Ampersand, pos)
            i += 1
          case '`' =>
            val start = i + 1
            var j = start
            while j < chars.length && chars(j) != '`' do j += 1
            if j >= chars.length then
              error = Some(
                FormulaError.Other(s"unterminated backtick-quoted identifier starting at position $pos")
              )
            else
              val name = chars.slice(start, j).mkString
              if name.isEmpty then
                error = Some(FormulaError.Other(s"empty backtick-quoted identifier at position $pos"))
              else
                tokens += Spanned(Token.Ident(name), pos)
                i = j + 1
          case '(' =>
            tokens += Spanned(Token.LParen, pos)
            i += 1
          case ')' =>
            tokens += Spanned(Token.RParen, pos)
            i += 1
          case '|' =>
            if i + 1 < chars.length && chars(i + 1) == '|' then
              tokens += Spanned(Token.DoublePipe, pos)
              i += 2
            else
              tokens += Spanned(Token.Pipe, pos)
              i += 1
          case '0' | '1' =>
            if i + 1 < chars.length && identContinue(chars(i + 1)) then
              val start = i
              var k = i
              while k < chars.length && identContinue(chars(k)) do k += 1
              tokens += Spanned(Token.Ident(chars.slice(start, k).mkString), pos)
              i = k
            else
              tokens += Spanned(if c == '0' then Token.Zero else Token.One, pos)
              i += 1
          case ch if ch.isLetter || ch == '_' || ch == '.' =>
            val start = i
            var k = i
            while k < chars.length && identContinue(chars(k)) do k += 1
            val word = chars.slice(start, k).mkString
            if k < chars.length && chars(k) == '(' && covarianceWrapper(word).isDefined then
              tokens += Spanned(Token.Ident(word), pos)
              i = k
            else if k < chars.length && chars(k) == '(' then
              matchingParen(chars, k) match
                case Left(err) =>
                  error = Some(err)
                case Right(end) =>
                  val inner = chars.slice(k + 1, end - 1).mkString
                  val parsed =
                    if word == "I" then Transform.parseArith(inner)
                    else if TransformFn.fromName(word).isDefined then Transform.parseBareCall(word, inner)
                    else
                      Left(
                        FormulaError.Other(
                          s"in-formula construct `$word(...)` at position $pos is not in the engine's stateless transform subset (allowed: `I(<+ - * / ^, unary -, parens, literals, columns>)` and pointwise `log`/`log2`/`log10`/`exp`/`sqrt`/`abs`). Stateful transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`, `center`, …) carry fitting-time state and must be precomputed as data columns or handled by the host wrapper."
                        )
                      )
                  parsed match
                    case Left(err) =>
                      error = Some(err)
                    case Right(expr) =>
                      val dc = DerivedColumn.from(expr)
                      if !seen.contains(dc.label) then
                        derived += dc
                        seen += dc.label
                      tokens += Spanned(Token.Ident(dc.label), pos)
                      i = end
            else
              tokens += Spanned(Token.Ident(word), pos)
              i = k
          case ch if ch.isDigit =>
            val start = i
            var k = i
            while k < chars.length && identContinue(chars(k)) do k += 1
            tokens += Spanned(Token.Ident(chars.slice(start, k).mkString), pos)
            i = k
          case ch if "^%=!<>".contains(ch) =>
            error = Some(
              FormulaError.Other(
                s"unexpected '$ch' at position $pos: bare formula-level arithmetic is not supported — wrap a stateless expression in `I(...)` (e.g. `I(x^2)`, `I(a*b)`, `I(1/x)`) or use a pointwise transform (`log`/`log2`/`log10`/`exp`/`sqrt`/`abs`, e.g. `sqrt(I(x + 1))`). Stateful transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`, `center`, …) carry fitting-time state and must be precomputed as data columns or handled by the host wrapper. If the column name itself contains unusual characters, quote it with backticks (e.g. `` `log x` ``)."
              )
            )
          case other =>
            error = Some(FormulaError.UnexpectedToken(other.toString, pos))
    error.toLeft(tokens.result())

  private def identContinue(c: Char): Boolean =
    c.isLetterOrDigit || c == '_' || c == '.'

  private def covarianceWrapper(name: String): Option[RandomCovariance] =
    name match
      case "us"   => Some(RandomCovariance.Full)
      case "diag" => Some(RandomCovariance.Diagonal)
      case "cs"   => Some(RandomCovariance.CompoundSymmetry)
      case "ar1"  => Some(RandomCovariance.Ar1)
      case _      => None

  private final case class ParsedGrouping(
      groupings: Vector[GroupingFactor],
      expansion: Option[RandomTermExpansion]
  )

  private final class FormulaParser(tokens: Vector[Spanned], input: String):
    private var cursor = 0

    def parseResponse(): Either[FormulaError, String] =
      peek match
        case Some(Token.Ident(name)) =>
          advance()
          Right(name)
        case _ =>
          Left(FormulaError.MissingResponse)

    def expect(expected: Token): Either[FormulaError, Unit] =
      peek match
        case Some(t) if t == expected =>
          advance()
          Right(())
        case Some(t) =>
          Left(FormulaError.UnexpectedToken(t.toString, pos))
        case None =>
          Left(FormulaError.UnexpectedToken("end of input", pos))

    def parseRhs(): Either[FormulaError, (Vector[FixedTerm], Vector[RandomTerm])] =
      val fixed = Vector.newBuilder[FixedTerm]
      val random = Vector.newBuilder[RandomTerm]
      var negate = false
      var expectTerm = true
      var pendingOperator = false
      var error: Option[FormulaError] = None
      while error.isEmpty && !atEnd do
        peek match
          case Some(Token.Plus) =>
            advance()
            negate = false
            expectTerm = true
            pendingOperator = true
          case Some(Token.Minus) =>
            advance()
            negate = true
            expectTerm = true
            pendingOperator = true
          case Some(Token.LParen) =>
            if !expectTerm then error = Some(FormulaError.MissingTermSeparator(pos))
            else if negate then error = Some(FormulaError.NegatedRandomEffect(pos))
            else
              parseRandomTerm() match
                case Left(err) => error = Some(err)
                case Right(rts) =>
                  random ++= rts
                  expectTerm = false
                  pendingOperator = false
          case Some(Token.One) =>
            if !expectTerm then error = Some(FormulaError.MissingTermSeparator(pos))
            else
              advance()
              if negate then fixed += FixedTerm.NoIntercept else fixed += FixedTerm.Intercept
              negate = false
              expectTerm = false
              pendingOperator = false
          case Some(Token.Zero) =>
            if !expectTerm then error = Some(FormulaError.MissingTermSeparator(pos))
            else
              advance()
              fixed += FixedTerm.NoIntercept
              negate = false
              expectTerm = false
              pendingOperator = false
          case Some(Token.Ident(name))
              if peekNext.contains(Token.LParen) && covarianceWrapper(name).isDefined =>
            if !expectTerm then error = Some(FormulaError.MissingTermSeparator(pos))
            else if negate then error = Some(FormulaError.NegatedRandomEffect(pos))
            else
              val sourceStart = pos
              advance()
              covarianceWrapper(name) match
                case None =>
                  error = Some(FormulaError.Other(s"unknown random-effect covariance wrapper `$name`"))
                case Some(cov) =>
                  parseRandomTermWithCovariance(sourceStart, cov) match
                    case Left(err) => error = Some(err)
                    case Right(rts) =>
                      random ++= rts
                      negate = false
                      expectTerm = false
                      pendingOperator = false
          case Some(Token.Ident(_)) =>
            if !expectTerm then error = Some(FormulaError.MissingTermSeparator(pos))
            else
              parseTermExpr() match
                case Left(err) => error = Some(err)
                case Right(terms) =>
                  if negate then
                    val current = fixed.result()
                    fixed.clear()
                    fixed ++= current.filterNot(terms.contains)
                  else fixed ++= terms
                  negate = false
                  expectTerm = false
                  pendingOperator = false
          case Some(other) =>
            error = Some(FormulaError.UnexpectedToken(other.toString, pos))
          case None =>
            ()
      if error.isEmpty && pendingOperator then Left(FormulaError.TrailingOperator(pos))
      else error.toLeft((fixed.result(), random.result()))

    private def parseTermExpr(): Either[FormulaError, Vector[FixedTerm]] =
      parseAtom().flatMap: first =>
        peek match
          case Some(Token.Colon) =>
            val names = Vector.newBuilder[String]
            names += first
            var err: Option[FormulaError] = None
            while err.isEmpty && peek.contains(Token.Colon) do
              advance()
              parseAtom() match
                case Left(e)     => err = Some(e)
                case Right(name) => names += name
            err.toLeft(Vector(FixedTerm.Interaction(names.result())))
          case Some(Token.Star) =>
            val names = Vector.newBuilder[String]
            names += first
            var err: Option[FormulaError] = None
            while err.isEmpty && peek.contains(Token.Star) do
              advance()
              parseAtom() match
                case Left(e)     => err = Some(e)
                case Right(name) => names += name
            err.toLeft(expandStarTerms(names.result()))
          case Some(Token.Slash) =>
            val names = Vector.newBuilder[String]
            names += first
            val terms = Vector.newBuilder[FixedTerm]
            terms += FixedTerm.Column(first)
            var err: Option[FormulaError] = None
            while err.isEmpty && peek.contains(Token.Slash) do
              advance()
              parseAtom() match
                case Left(e) => err = Some(e)
                case Right(name) =>
                  names += name
                  terms += FixedTerm.Interaction(names.result())
            err.toLeft(terms.result())
          case _ =>
            Right(Vector(FixedTerm.Column(first)))

    private def parseAtom(): Either[FormulaError, String] =
      peek match
        case Some(Token.Ident(name)) =>
          val at = pos
          advance()
          name.toDoubleOption match
            case Some(_) => Left(FormulaError.NumericLiteralTerm(name, at))
            case None    => Right(name)
        case Some(other) =>
          Left(FormulaError.UnexpectedToken(other.toString, pos))
        case None =>
          Left(FormulaError.UnexpectedToken("end of input", pos))

    private def parseRandomTerm(): Either[FormulaError, Vector[RandomTerm]] =
      parseRandomTermWithCovariance(pos, RandomCovariance.Full)

    private def parseRandomTermWithCovariance(
        sourceStart: Int,
        requested: RandomCovariance
    ): Either[FormulaError, Vector[RandomTerm]] =
      expect(Token.LParen).flatMap: _ =>
        val terms = Vector.newBuilder[FixedTerm]
        var negate = false
        var bar: Option[Boolean] = None
        var error: Option[FormulaError] = None
        while error.isEmpty && bar.isEmpty do
          peek match
            case Some(Token.Pipe) =>
              advance()
              bar = Some(false)
            case Some(Token.DoublePipe) =>
              advance()
              bar = Some(true)
            case Some(Token.RParen) =>
              error = Some(FormulaError.MissingBar)
            case None =>
              error = Some(FormulaError.UnmatchedParen)
            case Some(Token.Plus) =>
              advance()
              negate = false
            case Some(Token.Minus) =>
              advance()
              negate = true
            case Some(Token.One) =>
              advance()
              if negate then terms += FixedTerm.NoIntercept else terms += FixedTerm.Intercept
              negate = false
            case Some(Token.Zero) =>
              advance()
              terms += FixedTerm.NoIntercept
              negate = false
            case Some(Token.Ident(_)) =>
              parseTermExpr() match
                case Left(err) => error = Some(err)
                case Right(exprTerms) =>
                  terms ++= exprTerms
                  negate = false
            case Some(other) =>
              error = Some(FormulaError.UnexpectedToken(other.toString, pos))
        error match
          case Some(err) => Left(err)
          case None =>
            val collected = terms.result()
            if collected.isEmpty then Left(FormulaError.EmptyRandomTerms)
            else
              val withIntercept =
                if collected.exists(t => t == FixedTerm.Intercept || t == FixedTerm.NoIntercept) then collected
                else FixedTerm.Intercept +: collected
              val zerocorr = bar.contains(true)
              val covariance =
                if requested == RandomCovariance.Full && zerocorr then RandomCovariance.Diagonal
                else requested
              parseGrouping().flatMap: parsed =>
                val sourceEnd = pos
                expect(Token.RParen).map: _ =>
                  val written = sourceSpan(sourceStart, sourceEnd)
                  parsed.groupings.map: grouping =>
                    RandomTerm(
                      withIntercept,
                      grouping,
                      zerocorr,
                      covariance,
                      Some(RandomTermSource(written, parsed.expansion))
                    )

    private def parseGrouping(): Either[FormulaError, ParsedGrouping] =
      peek match
        case Some(Token.Ident(_)) =>
          parseAtom().flatMap: first =>
            peek match
              case Some(Token.Ampersand) =>
                val names = Vector.newBuilder[String]
                names += first
                var err: Option[FormulaError] = None
                while err.isEmpty && peek.contains(Token.Ampersand) do
                  advance()
                  parseAtom() match
                    case Left(e)     => err = Some(e)
                    case Right(name) => names += name
                err.toLeft(ParsedGrouping(Vector(GroupingFactor.Interaction(names.result())), None))
              case Some(Token.Colon) =>
                val names = Vector.newBuilder[String]
                names += first
                var err: Option[FormulaError] = None
                while err.isEmpty && peek.contains(Token.Colon) do
                  advance()
                  parseAtom() match
                    case Left(e)     => err = Some(e)
                    case Right(name) => names += name
                err.toLeft(ParsedGrouping(Vector(GroupingFactor.Cell(names.result())), None))
              case Some(Token.Slash) =>
                val names = Vector.newBuilder[String]
                names += first
                var err: Option[FormulaError] = None
                while err.isEmpty && peek.contains(Token.Slash) do
                  advance()
                  parseAtom() match
                    case Left(e)     => err = Some(e)
                    case Right(name) => names += name
                err.toLeft(
                  ParsedGrouping(expandNestedGrouping(names.result()), Some(RandomTermExpansion.NestedGrouping))
                )
              case Some(Token.Star) =>
                val names = Vector.newBuilder[String]
                names += first
                var err: Option[FormulaError] = None
                while err.isEmpty && peek.contains(Token.Star) do
                  advance()
                  parseAtom() match
                    case Left(e)     => err = Some(e)
                    case Right(name) => names += name
                err.toLeft(
                  ParsedGrouping(expandCrossedGrouping(names.result()), Some(RandomTermExpansion.CrossedGrouping))
                )
              case _ =>
                Right(ParsedGrouping(Vector(GroupingFactor.Single(first)), None))
        case _ =>
          Left(FormulaError.EmptyGrouping)

    private def peek: Option[Token] =
      tokens.lift(cursor).map(_.token)

    private def peekNext: Option[Token] =
      tokens.lift(cursor + 1).map(_.token)

    private def pos: Int =
      tokens.lift(cursor).map(_.pos).getOrElse(tokens.lastOption.map(_.pos + 1).getOrElse(0))

    private def advance(): Unit =
      if cursor < tokens.length then cursor += 1

    private def atEnd: Boolean = cursor >= tokens.length

    private def sourceSpan(start: Int, endInclusive: Int): String =
      input.drop(start).take(endInclusive - start + 1).trim

  private def expandStarTerms(names: Vector[String]): Vector[FixedTerm] =
    val terms = Vector.newBuilder[FixedTerm]
    (1 to names.length).foreach: size =>
      combinations(names, size).foreach: combo =>
        if combo.length == 1 then terms += FixedTerm.Column(combo.head)
        else terms += FixedTerm.Interaction(combo)
    terms.result()

  private def expandNestedGrouping(names: Vector[String]): Vector[GroupingFactor] =
    names.indices.map: end =>
      val slice = names.take(end + 1)
      if slice.length == 1 then GroupingFactor.Single(slice.head)
      else GroupingFactor.Cell(slice)
    .toVector

  private def expandCrossedGrouping(names: Vector[String]): Vector[GroupingFactor] =
    (1 to names.length).toVector.flatMap: size =>
      combinations(names, size).map: combo =>
        if combo.length == 1 then GroupingFactor.Single(combo.head)
        else GroupingFactor.Cell(combo)

  private def combinations(names: Vector[String], size: Int): Vector[Vector[String]] =
    def loop(start: Int, remaining: Int, current: Vector[String]): Vector[Vector[String]] =
      if remaining == 0 then Vector(current)
      else
        (start to names.length - remaining).toVector.flatMap: idx =>
          loop(idx + 1, remaining - 1, current :+ names(idx))
    loop(0, size, Vector.empty)
