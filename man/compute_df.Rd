% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_df.R
\name{compute_df}
\alias{compute_df}
\title{Compute degrees of freedom for the within or between family of contrasts}
\usage{
compute_df(model, data = NA, fctrs)
}
\arguments{
\item{model}{an anova model object of the afex_aov class}

\item{data}{instead of the anova model you can use a long form data frame}

\item{fctrs}{a list of characters or strings - list of relevant factor names}
}
\value{
a single value which is the degrees of freedom for that family
}
\description{
Given a list of within or between factors, and an object of the afex_aov
model class, this function will compute the degrees of freedom for that
family For further details see p. 216-217 from Bird (2002),
https://doi.org/10.1177/0013164402062002001.
}
\examples{
data("spacing")
spacing$group <- as.factor(spacing$group)
spacing$spacing <- as.factor(spacing$spacing)
 # get the degrees of freedom for the between family
 data("spacing")
 spacing$group <- as.factor(spacing$group)
 spacing$spacing <- as.factor(spacing$spacing)
 compute_df(model=NA, data=spacing, fctrs=list("group"))
# here is an example for multiple factor inputs, which is illegal
 # due to the factors being from both between and within families
# use the example with caution
compute_df(model=NA, data=spacing, fctrs=list("group", "spacing"))
}
