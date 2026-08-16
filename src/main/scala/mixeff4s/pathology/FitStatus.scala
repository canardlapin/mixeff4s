package mixeff4s.pathology

/** Contract fit status. Design-time certificates use `NotAssessed` until a fit is run. */
enum FitStatus:
  case ConvergedInterior, ConvergedBoundary, ConvergedReducedRank, ConvergedPenalised,
    NotIdentifiable, NotOptimized, NotAssessed

  def code: String =
    this match
      case FitStatus.ConvergedInterior    => "converged_interior"
      case FitStatus.ConvergedBoundary    => "converged_boundary"
      case FitStatus.ConvergedReducedRank => "converged_reduced_rank"
      case FitStatus.ConvergedPenalised   => "converged_penalised"
      case FitStatus.NotIdentifiable      => "not_identifiable"
      case FitStatus.NotOptimized         => "not_optimized"
      case FitStatus.NotAssessed          => "not_assessed"

enum Stratum:
  case Easy, Boundary, ReducedRank, Refusal, NotAssessed

  def code: String =
    this match
      case Stratum.Easy        => "easy"
      case Stratum.Boundary    => "boundary"
      case Stratum.ReducedRank => "reduced_rank"
      case Stratum.Refusal     => "refusal"
      case Stratum.NotAssessed => "not_assessed"

enum BoundaryKind:
  case ZeroVariance(index: Int)
  case UnitCorrelation(i: Int, j: Int)

  def code: String =
    this match
      case BoundaryKind.ZeroVariance(index) => s"zero_variance:$index"
      case BoundaryKind.UnitCorrelation(i, j) =>
        s"unit_correlation:$i,$j"

enum StructuralIssue:
  case RankSaturated(rank: Int, n: Int)
  case FewLevels(grouping: String, nLevels: Int)
  case SingletonsWithSlope(grouping: String, minGroupSize: Int)
  case InformationSaturated(nParams: Int, n: Int)
  case MalformedSpec(reason: String)
  case Separation(kind: SeparationKind)
  case CollinearFixedEffects(rank: Int, requested: Int)

  def code: String =
    this match
      case StructuralIssue.RankSaturated(_, _)        => "rank_saturated_fixed_effects"
      case StructuralIssue.FewLevels(_, _)            => "few_random_effect_levels"
      case StructuralIssue.SingletonsWithSlope(_, _)  => "singletons_with_slope"
      case StructuralIssue.InformationSaturated(_, _) => "information_saturated"
      case StructuralIssue.MalformedSpec(_)           => "malformed_spec"
      case StructuralIssue.Separation(_)              => "separation"
      case StructuralIssue.CollinearFixedEffects(_, _) =>
        "collinear_fixed_effects"

  def details: String =
    this match
      case StructuralIssue.RankSaturated(rank, n) =>
        s"fixed-effect rank $rank meets or exceeds n = $n"
      case StructuralIssue.FewLevels(grouping, nLevels) =>
        s"grouping `$grouping` has $nLevels level(s)"
      case StructuralIssue.SingletonsWithSlope(grouping, minGroupSize) =>
        s"grouping `$grouping` has a random slope and minimum group size $minGroupSize"
      case StructuralIssue.InformationSaturated(nParams, n) =>
        s"$nParams free parameters for n = $n observations"
      case StructuralIssue.MalformedSpec(reason) =>
        reason
      case StructuralIssue.Separation(kind) =>
        kind.details
      case StructuralIssue.CollinearFixedEffects(rank, requested) =>
        s"fixed-effect predictor rank $rank < requested $requested"

enum FeSeparationKind:
  case Complete, QuasiComplete

  def code: String =
    this match
      case FeSeparationKind.Complete      => "complete"
      case FeSeparationKind.QuasiComplete => "quasi_complete"

enum SeparationKind:
  case FixedEffect(kind: FeSeparationKind)
  case Conditional(nGroups: Int)
  case Both(feKind: FeSeparationKind, nGroups: Int)

  def details: String =
    this match
      case SeparationKind.FixedEffect(kind) =>
        s"fixed-effect ${kind.code} separation"
      case SeparationKind.Conditional(nGroups) =>
        s"conditional separation in $nGroups grouping level(s)"
      case SeparationKind.Both(feKind, nGroups) =>
        s"fixed-effect ${feKind.code} separation and conditional separation in $nGroups grouping level(s)"
