package mixeff4s.error

import mixeff4s.formula.FormulaError

/** Recoverable mixed-model failure. Display text may improve; [[code]] is the binding contract. */
enum MixedModelError:
  case Formula(error: FormulaError)
  case LinAlg(error: LinAlgError)
  case Optimization(details: String)
  case DimensionMismatch(details: String)
  case NotFitted
  case AlreadyFitted
  case Interrupted(details: String)
  case ConstantResponse
  case NoRandomEffects
  case InvalidArgument(details: String)
  case Unsupported(details: String)
  case UnsupportedFamilyLink(family: String, link: String)
  case ProblemTooLarge(details: String)
  case Singular(details: String)
  case RankSaturatedFixedEffects(rank: Int, nobs: Int)
  case PosDefException
  case InferenceUnavailable(reasonCode: String, details: String)

  def code: String =
    this match
      case MixedModelError.Formula(_)                      => "formula"
      case MixedModelError.LinAlg(_)                       => "linalg"
      case MixedModelError.Optimization(_)                 => "optimization"
      case MixedModelError.DimensionMismatch(_)            => "dimension_mismatch"
      case MixedModelError.NotFitted                       => "not_fitted"
      case MixedModelError.AlreadyFitted                   => "already_fitted"
      case MixedModelError.Interrupted(_)                  => "interrupted"
      case MixedModelError.ConstantResponse                => "constant_response"
      case MixedModelError.NoRandomEffects                 => "no_random_effects"
      case MixedModelError.InvalidArgument(_)              => "invalid_argument"
      case MixedModelError.Unsupported(_)                  => "unsupported"
      case MixedModelError.UnsupportedFamilyLink(_, _)     => "unsupported_family_link"
      case MixedModelError.ProblemTooLarge(_)              => "problem_too_large"
      case MixedModelError.Singular(_)                     => "singular_model"
      case MixedModelError.RankSaturatedFixedEffects(_, _) =>
        "rank_saturated_fixed_effects"
      case MixedModelError.PosDefException               => "positive_definite_exception"
      case MixedModelError.InferenceUnavailable(code, _) => code

  def message: String =
    this match
      case MixedModelError.Formula(error) =>
        s"Formula error: ${error.message}"
      case MixedModelError.LinAlg(error) =>
        s"Linear algebra error: ${error.message}"
      case MixedModelError.Optimization(details) =>
        s"Optimization error: $details"
      case MixedModelError.DimensionMismatch(details) =>
        s"Dimension mismatch: $details"
      case MixedModelError.NotFitted =>
        "Model not fitted: call fit() first"
      case MixedModelError.AlreadyFitted =>
        "Model already fitted: use refit() instead"
      case MixedModelError.Interrupted(details) =>
        s"Operation interrupted: $details"
      case MixedModelError.ConstantResponse =>
        "Constant response: model fitting failed"
      case MixedModelError.NoRandomEffects =>
        "No random effects in formula: this is not a mixed model"
      case MixedModelError.InvalidArgument(details) =>
        s"Invalid argument: $details"
      case MixedModelError.Unsupported(details) =>
        s"Unsupported model: $details"
      case MixedModelError.UnsupportedFamilyLink(family, link) =>
        s"Unsupported family/link combination: $family/$link"
      case MixedModelError.ProblemTooLarge(details) =>
        s"Problem too large: $details"
      case MixedModelError.Singular(details) =>
        s"Singular model: $details"
      case MixedModelError.RankSaturatedFixedEffects(rank, nobs) =>
        s"Fixed-effect design is rank-saturated: rank(X) = $rank and n = $nobs, leaving zero residual degrees of freedom."
      case MixedModelError.PosDefException =>
        "Positive definite exception during Cholesky"
      case MixedModelError.InferenceUnavailable(_, details) =>
        details

  override def toString: String = message

enum LinAlgError:
  case NotPositiveDefinite
  case DimensionMismatch(details: String)
  case Singular
  case RankDeficient(rank: Int, expected: Int)

  def code: String =
    this match
      case LinAlgError.NotPositiveDefinite  => "matrix_not_positive_definite"
      case LinAlgError.DimensionMismatch(_) => "dimension_mismatch"
      case LinAlgError.Singular             => "singular_matrix"
      case LinAlgError.RankDeficient(_, _)  => "rank_deficient"

  def message: String =
    this match
      case LinAlgError.NotPositiveDefinite =>
        "Matrix is not positive definite"
      case LinAlgError.DimensionMismatch(details) =>
        s"Dimension mismatch: $details"
      case LinAlgError.Singular =>
        "Singular matrix"
      case LinAlgError.RankDeficient(rank, expected) =>
        s"Rank deficient matrix (rank $rank, expected $expected)"

type FitResult[+A] = Either[MixedModelError, A]
