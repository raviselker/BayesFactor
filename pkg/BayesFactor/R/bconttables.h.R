
# This file is automatically generated, you probably don't want to edit this

bcontTablesOptions <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bcontTablesOptions",
    inherit = jmvcore::Options,
    public = list(
        initialize = function(
            rows = NULL,
            cols = NULL,
            counts = NULL,
            layers = NULL,
            sampling = "indepMultiRowsFixed",
            hypothesis = "different",
            bfType = "BF10",
            logOdds = FALSE,
            pp = FALSE,
            priorWidth = 1,
            ciWidth = 95,
            pcRow = FALSE,
            pcCol = FALSE,
            pcTot = FALSE, ...) {

            super$initialize(
                package='BayesFactor',
                name='bcontTables',
                requiresData=TRUE,
                ...)

            private$..rows <- jmvcore::OptionVariable$new(
                "rows",
                rows,
                suggested=list(
                    "nominal",
                    "ordinal"))
            private$..cols <- jmvcore::OptionVariable$new(
                "cols",
                cols,
                suggested=list(
                    "nominal",
                    "ordinal"))
            private$..counts <- jmvcore::OptionVariable$new(
                "counts",
                counts,
                suggested=list(
                    "continuous"),
                permitted=list(
                    "numeric"),
                default=NULL)
            private$..layers <- jmvcore::OptionVariables$new(
                "layers",
                layers,
                default=NULL)
            private$..sampling <- jmvcore::OptionList$new(
                "sampling",
                sampling,
                options=list(
                    "poisson",
                    "jointMulti",
                    "indepMultiRowsFixed",
                    "indepMultiColsFixed",
                    "hypergeom"),
                default="indepMultiRowsFixed")
            private$..hypothesis <- jmvcore::OptionList$new(
                "hypothesis",
                hypothesis,
                options=list(
                    "different",
                    "oneGreater",
                    "twoGreater"),
                default="different")
            private$..bfType <- jmvcore::OptionList$new(
                "bfType",
                bfType,
                options=list(
                    "BF10",
                    "BF01"),
                default="BF10")
            private$..logOdds <- jmvcore::OptionBool$new(
                "logOdds",
                logOdds,
                default=FALSE)
            private$..pp <- jmvcore::OptionBool$new(
                "pp",
                pp,
                default=FALSE)
            private$..priorWidth <- jmvcore::OptionNumber$new(
                "priorWidth",
                priorWidth,
                default=1,
                min=0.5,
                max=2)
            private$..ciWidth <- jmvcore::OptionNumber$new(
                "ciWidth",
                ciWidth,
                min=50,
                max=99.9,
                default=95)
            private$..pcRow <- jmvcore::OptionBool$new(
                "pcRow",
                pcRow,
                default=FALSE)
            private$..pcCol <- jmvcore::OptionBool$new(
                "pcCol",
                pcCol,
                default=FALSE)
            private$..pcTot <- jmvcore::OptionBool$new(
                "pcTot",
                pcTot,
                default=FALSE)

            self$.addOption(private$..rows)
            self$.addOption(private$..cols)
            self$.addOption(private$..counts)
            self$.addOption(private$..layers)
            self$.addOption(private$..sampling)
            self$.addOption(private$..hypothesis)
            self$.addOption(private$..bfType)
            self$.addOption(private$..logOdds)
            self$.addOption(private$..pp)
            self$.addOption(private$..priorWidth)
            self$.addOption(private$..ciWidth)
            self$.addOption(private$..pcRow)
            self$.addOption(private$..pcCol)
            self$.addOption(private$..pcTot)
        }),
    active = list(
        rows = function() private$..rows$value,
        cols = function() private$..cols$value,
        counts = function() private$..counts$value,
        layers = function() private$..layers$value,
        sampling = function() private$..sampling$value,
        hypothesis = function() private$..hypothesis$value,
        bfType = function() private$..bfType$value,
        logOdds = function() private$..logOdds$value,
        pp = function() private$..pp$value,
        priorWidth = function() private$..priorWidth$value,
        ciWidth = function() private$..ciWidth$value,
        pcRow = function() private$..pcRow$value,
        pcCol = function() private$..pcCol$value,
        pcTot = function() private$..pcTot$value),
    private = list(
        ..rows = NA,
        ..cols = NA,
        ..counts = NA,
        ..layers = NA,
        ..sampling = NA,
        ..hypothesis = NA,
        ..bfType = NA,
        ..logOdds = NA,
        ..pp = NA,
        ..priorWidth = NA,
        ..ciWidth = NA,
        ..pcRow = NA,
        ..pcCol = NA,
        ..pcTot = NA)
)

bcontTablesResults <- if (requireNamespace('jmvcore')) R6::R6Class(
    inherit = jmvcore::Group,
    active = list(
        freqs = function() private$.items[["freqs"]],
        tests = function() private$.items[["tests"]],
        odds = function() private$.items[["odds"]]),
    private = list(),
    public=list(
        initialize=function(options) {
            super$initialize(
                options=options,
                name="",
                title="Bayesian Contingency Tables")
            self$add(jmvcore::Table$new(
                options=options,
                name="freqs",
                title="Contingency Tables",
                columns=list(),
                clearWith=list(
                    "rows",
                    "cols",
                    "counts",
                    "layers")))
            self$add(jmvcore::Table$new(
                options=options,
                name="tests",
                title="Bayesian Tests",
                clearWith=list(
                    "rows",
                    "cols",
                    "counts",
                    "layers"),
                columns=list(
                    list(
                        `name`="test[bf10]", 
                        `title`="", 
                        `type`="text", 
                        `content`="BF\u2081\u2080", 
                        `visible`="(bfType:BF10)"),
                    list(
                        `name`="value[bf10]", 
                        `title`="Value", 
                        `visible`="(bfType:BF10)"),
                    list(
                        `name`="test[bf01]", 
                        `title`="", 
                        `type`="text", 
                        `content`="BF\u2080\u2081", 
                        `visible`="(bfType:BF01)"),
                    list(
                        `name`="value[bf01]", 
                        `title`="Value", 
                        `visible`="(bfType:BF01)"),
                    list(
                        `name`="test[N]", 
                        `title`="", 
                        `type`="text", 
                        `content`="N"),
                    list(
                        `name`="value[N]", 
                        `title`="Value", 
                        `type`="integer"))))
            self$add(jmvcore::Table$new(
                options=options,
                name="odds",
                title="Log Odds Ratio",
                visible="(logOdds)",
                clearWith=list(
                    "rows",
                    "cols",
                    "counts",
                    "layers",
                    "ciWidth"),
                columns=list(
                    list(
                        `name`="t[lo]", 
                        `title`="", 
                        `type`="text", 
                        `content`="Log odds ratio"),
                    list(
                        `name`="v[lo]", 
                        `title`="Value"),
                    list(
                        `name`="cil[lo]", 
                        `title`="Lower", 
                        `superTitle`="Confidence Intervals"),
                    list(
                        `name`="ciu[lo]", 
                        `title`="Upper", 
                        `superTitle`="Confidence Intervals"),
                    list(
                        `name`="t[f]", 
                        `title`="", 
                        `type`="text", 
                        `content`="Fisher's exact test"),
                    list(
                        `name`="v[f]", 
                        `title`="Value"),
                    list(
                        `name`="cil[f]", 
                        `title`="Lower", 
                        `superTitle`="Confidence Intervals"),
                    list(
                        `name`="ciu[f]", 
                        `title`="Upper", 
                        `superTitle`="Confidence Intervals"))))}))

bcontTablesBase <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bcontTablesBase",
    inherit = jmvcore::Analysis,
    public = list(
        initialize = function(options, data=NULL, datasetId="", analysisId="", revision=0) {
            super$initialize(
                package = 'BayesFactor',
                name = 'bcontTables',
                version = c(1,0,0),
                options = options,
                results = bcontTablesResults$new(options=options),
                data = data,
                datasetId = datasetId,
                analysisId = analysisId,
                revision = revision,
                pause = NULL,
                completeWhenFilled = TRUE)
        }))

#' Bayesian Contingency Tables
#'
#' 
#'
#' @examples
#' data('HairEyeColor')
#' dat <- as.data.frame(HairEyeColor)
#'
#' bcontTables(dat, rows = 'Hair', cols = 'Eye', counts = 'Freq')
#'
#' @param data the data as a data frame
#' @param rows a string naming the variable to use as the rows in the
#'   contingency table
#' @param cols a string naming the variable to use as the columns in the
#'   contingency table
#' @param counts a string naming the variable to use as counts, or NULL if
#'   each row represents a single observation
#' @param layers a character vector naming variables to split the contingency
#'   table across
#' @param sampling .
#' @param hypothesis \code{'different'} (default), \code{'oneGreater'} or
#'   \code{'twoGreater'}, the alternative hypothesis; group 1 different to group
#'   2, group 1 greater than group 2, and group 2 greater than group 1
#'   respectively
#' @param bfType \code{'BF10'} (default) or \code{'BF01'}
#' @param logOdds \code{TRUE} or \code{FALSE} (default), provide the log odds
#'   ratio (only available for 2x2 tables)
#' @param pp .
#' @param priorWidth .
#' @param ciWidth a number between 50 and 99.9 (default: 95), width of the
#'   confidence intervals to provide
#' @param pcRow \code{TRUE} or \code{FALSE} (default), provide row percentages
#' @param pcCol \code{TRUE} or \code{FALSE} (default), provide column
#'   percentages
#' @param pcTot \code{TRUE} or \code{FALSE} (default), provide total
#'   percentages
#' @return A results object containing:
#' \tabular{llllll}{
#'   \code{results$freqs} \tab \tab \tab \tab \tab a table of proportions \cr
#'   \code{results$tests} \tab \tab \tab \tab \tab a table \cr
#'   \code{results$odds} \tab \tab \tab \tab \tab a table of odds ratio results \cr
#' }
#'
#' Tables can be converted to data frames with \code{asDF} or \code{\link{as.data.frame}}. For example:
#'
#' \code{results$freqs$asDF}
#'
#' \code{as.data.frame(results$freqs)}
#'
#' @export
bcontTables <- function(
    data,
    rows,
    cols,
    counts = NULL,
    layers = NULL,
    sampling = "indepMultiRowsFixed",
    hypothesis = "different",
    bfType = "BF10",
    logOdds = FALSE,
    pp = FALSE,
    priorWidth = 1,
    ciWidth = 95,
    pcRow = FALSE,
    pcCol = FALSE,
    pcTot = FALSE) {

    if ( ! requireNamespace('jmvcore'))
        stop('bcontTables requires jmvcore to be installed (restart may be required)')

    if (missing(data))
        data <- jmvcore:::marshalData(
            parent.frame(),
            `if`( ! missing(rows), rows, NULL),
            `if`( ! missing(cols), cols, NULL),
            `if`( ! missing(counts), counts, NULL),
            `if`( ! missing(layers), layers, NULL))

    options <- bcontTablesOptions$new(
        rows = rows,
        cols = cols,
        counts = counts,
        layers = layers,
        sampling = sampling,
        hypothesis = hypothesis,
        bfType = bfType,
        logOdds = logOdds,
        pp = pp,
        priorWidth = priorWidth,
        ciWidth = ciWidth,
        pcRow = pcRow,
        pcCol = pcCol,
        pcTot = pcTot)

    results <- bcontTablesResults$new(
        options = options)

    analysis <- bcontTablesClass$new(
        options = options,
        data = data)

    analysis$run()

    analysis$results
}
