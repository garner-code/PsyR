#' Plant spacing data
#'
#' A small dataset for testing PsyR functionality. Data has plant yield
#' as the DV, and comes from a 4 (group) x 3 (spacing) mixed design.
#'
#' @format ## `spacing`
#' A data frame with 48 rows and 4 columns:
#' \describe{
#'   \item{subj}{subject}
#'   \item{group}{between subjects factor - 1:4}
#'   \item{spacing}{within subjects factor - TWENTY, FORTY, or SIXTY}
#'   \item{yield}{number of plants yielded}
#'   ...
#' }
#' @source <https://www.unsw.edu.au/science/our-schools/psychology/our-research/research-tools/psy-statistical-program>
"spacing"

#' Patient experience data
#'
#' A study examines the effect of therapist experience and therapy type on
#' patient outcome, where 132 patients are allocated to one of 12 groups in a
#' 3 x 4 two-way factorial design (J=3, K=4, n=11, N=132).
#'
#' @format ## `experience`
#' A data frame with 132 rows and 5 columns:
#' \describe{
#'   \item{Experience}{between subjects factor - Years experience as a therapist, 10 yrs, 3 yrs or <1yr}
#'   \item{Therapy}{between subjects factor - behavioural therapy, cognitive, psychotherapy, group}
#'   \item{Group}{a factor denoting group membership based on the factorial crossover of experience x therapy}
#'   \item{GrpNum}{same as Group, but numerical}
#'   \item{Outcome}{outcome measure. higher score is better}
#'   ...
#' }
"experience"
