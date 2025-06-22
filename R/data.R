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

#' Priming lexical decision task data
#'
#' A synthetic dataset modeling a lexical decision task used to investigate automatic activation of stereotypes through priming. Participants are subliminally primed with a social category label (e.g., ‘unemployed’ or ‘elderly’) before responding to target words. The target words are either stereotypical attributes, atypical attributes, or neutral words. The dependent variable is reaction time (in tenths of a second).
#'
#' Priming is hypothesized to enhance recognition for stereotypical attributes (i.e., faster reaction times), and inhibit recognition for atypical attributes, relative to neutral words.
#'
#' Contrasts of interest include:
#' \itemize{
#'   \item{"Other-Stereo" = c(-2, 1, 1)/2} — comparing stereotypical words against the average of atypical and neutral words
#'   \item{"Atyp-Neut" = c(0, 1, -1)} — directly comparing atypical and neutral words
#' }
#'
#' @format ## `priming`
#' A data frame with 15 rows and 4 columns:
#' \describe{
#'   \item{subj}{subject ID}
#'   \item{Condition}{within-subjects factor with 3 levels: STEREO, ATYPICAL, NEUTRAL}
#'   \item{WordType}{same as Condition, possibly labeled differently}
#'   \item{RT}{reaction time in tenths of a second}
#' }
#' @source <https://www.unsw.edu.au/science/our-schools/psychology/our-research/research-tools/psy-statistical-program>
"priming"

