package mixeff4s.pathology

/** Contract fit status. Design-time certificates use `NotAssessed` until a fit is run. */
enum FitStatus:
  case ConvergedInterior, ConvergedBoundary, ConvergedReducedRank, NotIdentifiable, NotOptimized,
    NotAssessed

  def code: String =
    this match
      case FitStatus.ConvergedInterior    => "converged_interior"
      case FitStatus.ConvergedBoundary    => "converged_boundary"
      case FitStatus.ConvergedReducedRank => "converged_reduced_rank"
      case FitStatus.NotIdentifiable      => "not_identifiable"
      case FitStatus.NotOptimized         => "not_optimized"
      case FitStatus.NotAssessed          => "not_assessed"

enum Stratum:
  case Easy, Refusal, NotAssessed

  def code: String =
    this match
      case Stratum.Easy        => "easy"
      case Stratum.Refusal     => "refusal"
      case Stratum.NotAssessed => "not_assessed"

enum StructuralIssue:
  case RankSaturated(rank: Int, n: Int)
  case FewLevels(grouping: String, nLevels: Int)
  case SingletonsWithSlope(grouping: String, minGroupSize: Int)
  case InformationSaturated(nParams: Int, n: Int)

  def code: String =
    this match
      case StructuralIssue.RankSaturated(_, _)         => "rank_saturated_fixed_effects"
      case StructuralIssue.FewLevels(_, _)             => "few_random_effect_levels"
      case StructuralIssue.SingletonsWithSlope(_, _)   => "singletons_with_slope"
      case StructuralIssue.InformationSaturated(_, _)  => "information_saturated"

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
