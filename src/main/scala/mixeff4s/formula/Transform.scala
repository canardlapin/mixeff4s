package mixeff4s.formula

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.MixedModelError

/** Whitelisted single-argument pointwise function. */
enum TransformFn:
  case Ln, Log2, Log10, Exp, Sqrt, Abs

  def label: String =
    this match
      case TransformFn.Ln    => "log"
      case TransformFn.Log2  => "log2"
      case TransformFn.Log10 => "log10"
      case TransformFn.Exp   => "exp"
      case TransformFn.Sqrt  => "sqrt"
      case TransformFn.Abs   => "abs"

  def apply(x: Double): Double =
    this match
      case TransformFn.Ln    => math.log(x)
      case TransformFn.Log2  => math.log(x) / math.log(2.0)
      case TransformFn.Log10 => math.log10(x)
      case TransformFn.Exp   => math.exp(x)
      case TransformFn.Sqrt  => math.sqrt(x)
      case TransformFn.Abs   => math.abs(x)

object TransformFn:
  def fromName(name: String): Option[TransformFn] =
    name match
      case "log"   => Some(TransformFn.Ln)
      case "log2"  => Some(TransformFn.Log2)
      case "log10" => Some(TransformFn.Log10)
      case "exp"   => Some(TransformFn.Exp)
      case "sqrt"  => Some(TransformFn.Sqrt)
      case "abs"   => Some(TransformFn.Abs)
      case _       => None

enum BinOp:
  case Add, Sub, Mul, Div, Pow

  def symbol: Char =
    this match
      case BinOp.Add => '+'
      case BinOp.Sub => '-'
      case BinOp.Mul => '*'
      case BinOp.Div => '/'
      case BinOp.Pow => '^'

  def precedence: Int =
    this match
      case BinOp.Add | BinOp.Sub => 1
      case BinOp.Mul | BinOp.Div => 2
      case BinOp.Pow             => 3

  def apply(a: Double, b: Double): Double =
    this match
      case BinOp.Add => a + b
      case BinOp.Sub => a - b
      case BinOp.Mul => a * b
      case BinOp.Div => a / b
      case BinOp.Pow => math.pow(a, b)

enum Expr:
  case Lit(value: Double)
  case Col(name: String)
  case Neg(inner: Expr)
  case Bin(op: BinOp, left: Expr, right: Expr)
  case Call(fn: TransformFn, arg: Expr)

final case class DerivedColumn(label: String, expr: Expr)

object DerivedColumn:
  def from(expr: Expr): DerivedColumn =
    new DerivedColumn(Transform.canonicalLabel(expr), expr)

object Transform:
  private val MaxDepth = 64

  def canonicalLabel(expr: Expr): String =
    expr match
      case Expr.Col(name) => name
      case Expr.Call(fn, arg) =>
        s"${fn.label}(${innerArgLabel(arg)})"
      case _ =>
        val out = StringBuilder("I(")
        writeExpr(out, expr, 0)
        out.append(')')
        out.toString

  def parseArith(src: String): Either[FormulaError, Expr] =
    for
      toks <- lex(src)
      expr <- TParser(toks).parseExpr(0)
    yield expr

  def parseBareCall(name: String, argSrc: String): Either[FormulaError, Expr] =
    TransformFn.fromName(name) match
      case None => Left(refuse(s"$name(…)"))
      case Some(fn) =>
        parseCallArgument(argSrc).map(arg => Expr.Call(fn, arg))

  def eval(expr: Expr, frame: ModelFrame): Either[MixedModelError, Vector[Double]] =
    val n = frame.nRows
    val buf = Vector.newBuilder[Double]
    var row = 0
    var error: Option[MixedModelError] = None
    while row < n && error.isEmpty do
      evalRow(expr, frame, row) match
        case Left(err) => error = Some(err)
        case Right(value) =>
          buf += value
          row += 1
    error.toLeft(buf.result())

  def materializeColumn(derived: DerivedColumn, frame: ModelFrame): Either[MixedModelError, Vector[Double]] =
    eval(derived.expr, frame).flatMap: values =>
      values.indexWhere(v => !v.isFinite) match
        case -1 => Right(values)
        case pos =>
          Left(
            MixedModelError.InvalidArgument(
              s"in-formula transform `${derived.label}` produced a non-finite value (${values(pos)}) at row $pos; the transform is undefined there (e.g. log/sqrt of a non-positive value, or division by zero) — clean or restrict the data before fitting"
            )
          )

  private def innerArgLabel(expr: Expr): String =
    expr match
      case Expr.Col(name)   => name
      case Expr.Call(_, _)  => canonicalLabel(expr)
      case _ =>
        val out = StringBuilder("I(")
        writeExpr(out, expr, 0)
        out.append(')')
        out.toString

  private def fmtLit(v: Double): String =
    if v.isFinite && v == math.rint(v) && math.abs(v) < 1e15 then v.toLong.toString
    else v.toString

  private def writeExpr(out: StringBuilder, expr: Expr, parentPrec: Int): Unit =
    expr match
      case Expr.Lit(v) =>
        out.append(fmtLit(v))
        ()
      case Expr.Col(name) =>
        out.append(name)
        ()
      case Expr.Neg(inner) =>
        val need = parentPrec > BinOp.Mul.precedence
        if need then out.append('(')
        out.append('-')
        writeExpr(out, inner, BinOp.Pow.precedence)
        if need then out.append(')')
        ()
      case Expr.Bin(op, a, b) =>
        val p = op.precedence
        val need = p < parentPrec
        if need then out.append('(')
        writeExpr(out, a, p)
        out.append(op.symbol)
        writeExpr(out, b, p + 1)
        if need then out.append(')')
        ()
      case Expr.Call(fn, arg) =>
        out.append(fn.label)
        out.append('(')
        out.append(innerArgLabel(arg))
        out.append(')')
        ()

  private def evalRow(expr: Expr, frame: ModelFrame, row: Int): Either[MixedModelError, Double] =
    expr match
      case Expr.Lit(v) => Right(v)
      case Expr.Col(name) =>
        frame.column(name) match
          case Some(Column.Numeric(values)) =>
            Right(values(row))
          case Some(_) =>
            Left(
              MixedModelError.InvalidArgument(
                s"in-formula transform references categorical column `$name`; stateless transforms operate on numeric columns only — precompute a numeric encoding or use the host wrapper for factor handling"
              )
            )
          case None =>
            Left(
              MixedModelError.InvalidArgument(
                s"in-formula transform references column `$name`, which is not present in the data"
              )
            )
      case Expr.Neg(inner) =>
        evalRow(inner, frame, row).map(v => -v)
      case Expr.Bin(op, a, b) =>
        for
          lhs <- evalRow(a, frame, row)
          rhs <- evalRow(b, frame, row)
        yield op(lhs, rhs)
      case Expr.Call(fn, arg) =>
        evalRow(arg, frame, row).map(fn.apply)

  private def refuse(construct: String): FormulaError =
    FormulaError.Other(
      s"in-formula construct `$construct` is not in the engine's stateless transform subset (allowed: `I(<+ - * / ^, unary -, parens, literals, columns>)` and pointwise `log`/`log2`/`log10`/`exp`/`sqrt`/`abs`). Stateful transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`, `center`, …) carry fitting-time state and must be precomputed as data columns or handled by the host wrapper."
    )

  private def parseCallArgument(src: String): Either[FormulaError, Expr] =
    val trimmed = src.trim
    if trimmed.contains(',') then
      Left(
        FormulaError.Other(
          s"multi-argument call `…($trimmed)` is not a stateless pointwise transform and is out of scope for the engine — base changes like `log(x, base)` are stateful; precompute the column or handle it in the host wrapper"
        )
      )
    else
      for
        toks <- lex(trimmed)
        parser = TParser(toks)
        expr <- parser.parsePrimary
        _ <-
          if parser.atEnd then Right(())
          else
            Left(
              FormulaError.Other(
                s"argument `$trimmed` to a pointwise transform must be a column, a nested whitelisted call, or an `I(...)` arithmetic expression — wrap arithmetic in `I(...)` (e.g. `sqrt(I(x + 1))`)"
              )
            )
      yield expr

  private enum Tok:
    case Num(value: Double)
    case Ident(name: String)
    case Op(op: BinOp)
    case LParen
    case RParen
    case Call(fn: TransformFn)
    case IOpen

  private def lex(src: String): Either[FormulaError, Vector[Tok]] =
    val chars = src.toVector
    val out = Vector.newBuilder[Tok]
    var i = 0
    var error: Option[FormulaError] = None
    while i < chars.length && error.isEmpty do
      val c = chars(i)
      if c.isWhitespace then i += 1
      else
        c match
          case '+' =>
            out += Tok.Op(BinOp.Add)
            i += 1
          case '-' =>
            out += Tok.Op(BinOp.Sub)
            i += 1
          case '*' =>
            out += Tok.Op(BinOp.Mul)
            i += 1
          case '/' =>
            out += Tok.Op(BinOp.Div)
            i += 1
          case '^' =>
            out += Tok.Op(BinOp.Pow)
            i += 1
          case '(' =>
            out += Tok.LParen
            i += 1
          case ')' =>
            out += Tok.RParen
            i += 1
          case ',' =>
            error = Some(
              FormulaError.Other(
                "multi-argument calls are not stateless pointwise transforms and are out of scope for the engine — precompute the column or handle it in the host wrapper"
              )
            )
          case '`' =>
            val start = i + 1
            var j = start
            while j < chars.length && chars(j) != '`' do j += 1
            if j >= chars.length then
              error = Some(FormulaError.Other("unterminated backtick-quoted identifier inside transform"))
            else
              out += Tok.Ident(chars.slice(start, j).mkString)
              i = j + 1
          case ch if ch.isDigit || ch == '.' =>
            val start = i
            var k = i
            while k < chars.length && (
                chars(k).isDigit || chars(k) == '.' || chars(k) == 'e' || chars(k) == 'E'
                  || ((chars(k) == '+' || chars(k) == '-') && k > start && (chars(k - 1) == 'e' || chars(k - 1) == 'E'))
              )
            do k += 1
            val lit = chars.slice(start, k).mkString
            lit.toDoubleOption match
              case Some(v) =>
                out += Tok.Num(v)
                i = k
              case None =>
                error = Some(FormulaError.Other(s"invalid numeric literal `$lit` in transform"))
          case ch if ch.isLetter || ch == '_' || ch == '.' =>
            val start = i
            var k = i
            while k < chars.length && (chars(k).isLetterOrDigit || chars(k) == '_' || chars(k) == '.') do k += 1
            val word = chars.slice(start, k).mkString
            var look = k
            while look < chars.length && chars(look).isWhitespace do look += 1
            if look < chars.length && chars(look) == '(' then
              if word == "I" then
                out += Tok.IOpen
                out += Tok.LParen
                i = look + 1
              else
                TransformFn.fromName(word) match
                  case Some(fn) =>
                    out += Tok.Call(fn)
                    out += Tok.LParen
                    i = look + 1
                  case None =>
                    error = Some(refuse(s"$word(…)"))
            else
              out += Tok.Ident(word)
              i = k
          case other =>
            error = Some(
              FormulaError.Other(
                s"unexpected character `$other` in in-formula transform; allowed inside `I(...)`: `+ - * / ^`, unary `-`, parentheses, numeric literals, and column references"
              )
            )
    error.toLeft(out.result())

  private final class TParser(toks: Vector[Tok]):
    var pos: Int = 0
    private var depth: Int = 0

    def atEnd: Boolean = pos >= toks.length

    def parseExpr(minPrec: Int): Either[FormulaError, Expr] =
      enter().flatMap: _ =>
        val result = parseExprInner(minPrec)
        leave()
        result

    def parsePrimary: Either[FormulaError, Expr] =
      enter().flatMap: _ =>
        val result = parsePrimaryInner
        leave()
        result

    private def peek: Option[Tok] = toks.lift(pos)

    private def bump(): Option[Tok] =
      val tok = toks.lift(pos)
      if tok.isDefined then pos += 1
      tok

    private def enter(): Either[FormulaError, Unit] =
      depth += 1
      if depth > MaxDepth then
        Left(
          FormulaError.Other(
            s"in-formula transform nesting exceeds the maximum supported depth ($MaxDepth); deeply nested parentheses/calls are rejected to bound parser recursion"
          )
        )
      else Right(())

    private def leave(): Unit =
      depth -= 1

    private def parseExprInner(minPrec: Int): Either[FormulaError, Expr] =
      parseUnary().flatMap: start =>
        var lhs = start
        var error: Option[FormulaError] = None
        var done = false
        while !done && error.isEmpty do
          peek match
            case Some(Tok.Op(op)) if op.precedence >= minPrec =>
              bump()
              val nextMin = if op == BinOp.Pow then op.precedence else op.precedence + 1
              parseExpr(nextMin) match
                case Left(err) =>
                  error = Some(err)
                case Right(rhs) =>
                  lhs = Expr.Bin(op, lhs, rhs)
            case _ =>
              done = true
        error.toLeft(lhs)

    private def parseUnary(): Either[FormulaError, Expr] =
      peek match
        case Some(Tok.Op(BinOp.Sub)) =>
          bump()
          parseExpr(BinOp.Pow.precedence).map(Expr.Neg.apply)
        case _ =>
          while peek.contains(Tok.Op(BinOp.Add)) do bump(): Unit
          parsePrimary

    private def parsePrimaryInner: Either[FormulaError, Expr] =
      bump() match
        case Some(Tok.Num(v)) =>
          Right(Expr.Lit(v))
        case Some(Tok.Ident(name)) =>
          Right(Expr.Col(name))
        case Some(Tok.LParen) =>
          for
            e <- parseExpr(0)
            _ <- bump() match
              case Some(Tok.RParen) => Right(())
              case _ =>
                Left(FormulaError.Other("unbalanced parentheses in in-formula transform"))
          yield e
        case Some(Tok.IOpen) =>
          bump() match
            case Some(Tok.LParen) =>
              peek match
                case Some(Tok.RParen) =>
                  Left(
                    FormulaError.Other(
                      "empty `I(...)` — expected an arithmetic expression inside the parentheses (e.g. `I(x^2)`, `I(a*b)`, `I(1/x)`)"
                    )
                  )
                case _ =>
                  for
                    e <- parseExpr(0)
                    _ <- bump() match
                      case Some(Tok.RParen) => Right(())
                      case _ =>
                        Left(FormulaError.Other("unbalanced parentheses in `I(...)`"))
                  yield e
            case _ =>
              Left(FormulaError.Other("malformed `I(...)` in in-formula transform"))
        case Some(Tok.Call(fn)) =>
          bump() match
            case Some(Tok.LParen) =>
              for
                arg <- parsePrimary
                _ <- bump() match
                  case Some(Tok.RParen) => Right(())
                  case _ =>
                    Left(
                      FormulaError.Other(
                        s"pointwise `${fn.label}` takes exactly one argument (a column, a nested whitelisted call, or `I(...)`); a second argument is stateful and out of scope"
                      )
                    )
              yield Expr.Call(fn, arg)
            case _ =>
              Left(FormulaError.Other("malformed pointwise call in in-formula transform"))
        case other =>
          Left(FormulaError.Other(s"unexpected token $other in in-formula transform"))
