\name{log_diameter_peak_boxplot}
\alias{log_diameter_peak_boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot for Log Diameter, split by Peak
%%  ~~function to do ... ~~
}
\description{Boxplot for log diameter split into variable of choice for a specific peak, colour coded based on variable of interest.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
log_diameter_peak_boxplot(data, peak, xvar, variable, target, original, indent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data to be visualised.
  \item{peak}{Peak to be visualised.
  \item{xvar}{Variable to split the boxplot along the x-axis (grouping variable).
  \item{variable}{Variable of interest to colour-code and split the data.
  \item{target}{Adds a line equivalent to 10nm (target diameter of SMA nanodiscs).
  \item{original}{Adds a line equivalent to 200nm (starting size of DMPC vesicles).
  \item{indent}{Function to shift the labels of the original and the target lines horizontally to not impede the data.
}
}
}
}
}
}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
