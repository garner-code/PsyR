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
#' experience - patient experience data
#'
#' A study examining the effect of experience of therapist and type of therapy
#' on patient outcome. 132 patients are randomly allocated to one of 12 groups
#' in a 3 x 4 two-way factorial design (between-groups).
#'
#' @format ## `experience`
#' A data frame with 132 rows and 5 columns:
#' \describe{
#'   \item{Experience}{Yrs experience of therapist: 10, 3, or <1}
#'   \item{Therapy}{Therapy type: Behavioural, Cognitive, Psychotherapy or Group}
#'   \item{GrpNum}{Group number: 1-12}
#'   \item{Outcome}{Patient wellbeing score: higher is better}
#'   \item{SubjID}{Subject ID}
#'   ...
#' }
#' @source <UNSW: Dr. Melanie Gleitzman>
"experience"
#' depression - patient happiness scores from pre to post and to follow up
#'
#' Data showing happiness scores from 15 patients. Scores taken at Pre, Post and
#' Follow Up (FU) to a treatment programme. Patients are allocated to one of three
#' Treatment Groups (Tmt1, Tmt2, Ctrl).
#'
#' @format ## `depression`
#' A data frame with 45 rows and 5 columns:
#' \describe{
#'   \item{Group}{Treatment Group: Tmt1 (treatment 1), Tmt2, Ctrl}
#'   \item{GrpNum}{Numerical vrsn of Group (1, 2, or 3)}
#'   \item{Time}{Timepoint of happiness measurement - Pre, Post, FU}
#'   \item{Happiness}{Happiness score: higher is better}
#'   \item{Subject}{Subject number}
#'   ...
#' }
#' @source <UNSW: Melanie Gleitzman>
"depression"
