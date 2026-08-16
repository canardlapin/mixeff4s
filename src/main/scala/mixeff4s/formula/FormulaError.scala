package mixeff4s.formula

/** Errors produced while parsing a formula string. */
enum FormulaError:
  case Empty
  case MissingTilde
  case MissingResponse
  case UnexpectedToken(token: String, position: Int)
  case UnmatchedParen
  case MissingBar
  case EmptyGrouping
  case EmptyRandomTerms
  case EmptyRhs
  case TrailingOperator(position: Int)
  case MissingTermSeparator(position: Int)
  case NumericLiteralTerm(literal: String, position: Int)
  case NegatedRandomEffect(position: Int)
  case Other(details: String)

  def code: String =
    this match
      case FormulaError.Empty                    => "empty"
      case FormulaError.MissingTilde             => "missing_tilde"
      case FormulaError.MissingResponse          => "missing_response"
      case FormulaError.UnexpectedToken(_, _)    => "unexpected_token"
      case FormulaError.UnmatchedParen           => "unmatched_paren"
      case FormulaError.MissingBar               => "missing_bar"
      case FormulaError.EmptyGrouping            => "empty_grouping"
      case FormulaError.EmptyRandomTerms         => "empty_random_terms"
      case FormulaError.EmptyRhs                 => "empty_rhs"
      case FormulaError.TrailingOperator(_)      => "trailing_operator"
      case FormulaError.MissingTermSeparator(_)  => "missing_term_separator"
      case FormulaError.NumericLiteralTerm(_, _) => "numeric_literal_term"
      case FormulaError.NegatedRandomEffect(_)   => "negated_random_effect"
      case FormulaError.Other(_)                 => "other"

  def message: String =
    this match
      case FormulaError.Empty =>
        "empty formula string"
      case FormulaError.MissingTilde =>
        "formula must contain '~' separating response from predictors"
      case FormulaError.MissingResponse =>
        "missing response variable on the left-hand side of '~'"
      case FormulaError.UnexpectedToken(token, position) =>
        s"unexpected token '$token' at position $position"
      case FormulaError.UnmatchedParen =>
        "unmatched opening parenthesis — expected ')'"
      case FormulaError.MissingBar =>
        "random-effect term is missing '|' or '||' separator"
      case FormulaError.EmptyGrouping =>
        "random-effect term has an empty grouping factor"
      case FormulaError.EmptyRandomTerms =>
        "random-effect term has no model terms before '|'"
      case FormulaError.EmptyRhs =>
        "formula has no terms on the right-hand side of '~'"
      case FormulaError.TrailingOperator(position) =>
        s"formula ends with a dangling '+'/'-' operator at position $position"
      case FormulaError.MissingTermSeparator(position) =>
        s"expected '+' or '-' separating model terms at position $position"
      case FormulaError.NumericLiteralTerm(literal, position) =>
        s"numeric literal '$literal' at position $position is not a valid model term (only the 0/1 intercept literals are allowed)"
      case FormulaError.NegatedRandomEffect(position) =>
        s"'-' cannot remove a random-effect term at position $position"
      case FormulaError.Other(details) =>
        details

  override def toString: String = message
