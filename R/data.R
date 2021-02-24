#' Baseball hitters career and 1986 season total statistics
#'
#' 16 statistics of 263 major league hitters some of the
#' overall career and others relative to the 1986 season. Available at StatLib.
#' The matrix has a block structure, defined by season offensive play, career
#' offensive play and season defensive play.
#'
#'
#' @name hitters
#' @docType data
#' @format A \emph{263} by \emph{16} matrix. The variables are centered
#' and standardized to unit norm.
#' \describe{
#' \item{TAB_86}{times at bat in 1986}
#' \item{HIT_86}{hits in 1986}
#' \item{HR_86}{home runs in 1986}
#' \item{RUN_86}{runs in 1986}
#' \item{RB_86}{runs batted-in in 1986}
#' \item{WAL_86}{walks in 1986}
#' \item{WAL}{walks during his career}
#' \item{PO_86}{put outs in 1986}
#' \item{ASS_86}{assists in 1986}
#' \item{ERR_86}{errors in 1986}
#' \item{YC}{years in the major leagues}
#' \item{TAB}{times at bat during his career}
#' \item{HIT}{hits during his career}
#' \item{HR}{home runs during his career}
#' \item{RUN}{runs during his career}
#' \item{RUNB}{runs batted-in during his career}
#' }
#' @source \url{http://lib.stat.cmu.edu/datasets/baseball.data}
#' @keywords datasets
NULL


#' Baseball hitters statistics labels reference table
#'
#' This data frame provides descriptive labels for the variables in the
#' hitters datasets matching the short labels.
#'
#' @name hitters_labels
#' @docType data
#' @format A dataframe with columns
#' \describe{
#' \item{short.namethe}{labels ised in plots and tables}
#' \item{label}{the explanatory name}
#' \item{type}{whether offensive year 1986, defensive year 86 or offensive over the whole career}
#' }
#' @source \url{http://lib.stat.cmu.edu/datasets/baseball.data}
#' @keywords datasets
NULL

#' Holzinger student's ability data
#'
#'Holzinger and Swineford (1939) data is widely cited. This dataset
#' includes only the 12 tests on the Grant-White School students used in Ferrara et al .

#' @name holzinger
#' @docType data
#' @format A \emph{263} by \emph{17} matrix. The first 16 variables are centered
#' and standardized to unit norm.
#'
#' \strong{spatial} (SPL)
#' \describe{
#' \item{visual}{Visual perception test, a nonlanguage multiple-choice test of spatial relations }
#' \item{cubes}{Cubes test, spatial relations}
#' \item{flags}{Lozenges test, a visual imagery test in two or three dimensions}
#' }
#' \strong{verbal} (VBL)
#' \describe{
#' \item{paragraph}{Paragraph comprehension test, comprehension as measured by completion and multiple-choice questions}
#' \item{sentence}{Sentence completion test, a multiple-choice test in which “correct” answers reflect good judgment on the part of the subject}
#' \item{wordm}{Word meaning test, a multiple-choice vocabulary test}
#' }
#' \strong{speed} (SPD)
#' \describe{
#' \item{addition}{Addition test, speed of adding pairs of one-digit numbers}
#' \item{counting}{Counting groups of dots test, 4–7 dots, arranged in random patterns to be counted by subject}
#' \item{straight}{Straight and curved capitals test, a series of capital letters to be distinguished between those composed of straight lines only and those containing curved lines}
#' }
#' \strong{mathematical} (MTH)
#' \describe{
#' \item{deduct}{Deduction test, logical deduction test using the symbols and the letters}
#' \item{numeric}{Numerical puzzles test, a numerical deduction test, the object being to supply four numbers which will produce four given answers employing the operations of addition, multiplication, or division}
#' \item{series}{Series completion test, from a series of five numbers, the subject is supposed to deduce the rule of procedure from one number to the next and thus supply the sixth number in the series}
#' }
#' @source The complete set of data is included, for example in the R package MBESS.

#' @references
#' Holzinger, K. J. and Swineford, F. A. (1939). A study in factor analysis: The stability of a bi-factor solution. Supplementary Education Monographs, 48. University of Chicago.
#'
#' Ferrara, C., Martella, F., and Vichi, M. (2019). Probabilistic disjoint principal component analysis. Multivariate Behavioral Research, 54(1):47–61.
#' @keywords datasets
NULL


#' Holzinger student's ability data labels reference table
#'
#' This data frame provides descriptive labels for the variables in the
#' Holzinger student's ability data matching the short labels used.
#'
#' @name holzinger_labels
#' @docType data
#' @format A dataframe with columns
#' \describe{
#' \item{Test}{Name of each test}
#' \item{ Ability}{Short label for the abilities}
#' \item{Long_Ability}{Long name of th abilities}
#' }
#' @keywords datasets
NULL

