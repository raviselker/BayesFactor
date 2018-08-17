
# This file is a generated template, your changes will not be overwritten

blinRegClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "blinRegClass",
    inherit = blinRegBase,
    private = list(
      .init=function() {

        dep <- self$options$dep
        fixed <- self$options$covs
        modelTerms <- self$options$modelTerms

        anovaTable <- self$results$main

        anovaTable$addRow(rowKey='', list(name='Null model'))

        if (length(modelTerms) > 0) {
          data <- private$.cleanData()
          fmla <- jmvcore::constructFormula('.', modelTerms)
          fmla <- as.formula(fmla)
          models <- enumerateGeneralModels(fmla=fmla, whichModels='withmain', data=data)

          for (model in models) {
            terms <- attr(terms.formula(model), 'term.labels')
            terms <- jmvcore::decomposeTerms(terms)
            niceTerms <- vapply(terms, jmvcore::stringifyTerm, '', USE.NAMES=FALSE)
            name <- paste(niceTerms, collapse=' + ')
            anovaTable$addRow(rowKey=terms, list(name=name))
          }

        } else {
          anovaTable$addRow(rowKey='.', list(name='.'))
        }

      },
      .run = function() {

        if (is.null(self$options$dep) || length(self$options$covs) == 0)
          return()


      },
      .cleanData=function() {

        dep <- self$options$dep
        covs <- self$options$covs

        data <- self$data

        if ( ! is.null(dep))
          data[[dep]] <- jmvcore::toNumeric(data[[dep]])

        for (covariate in covs)
          data[[covariate]] <- jmvcore::toNumeric(data[[covariate]])

        data <- na.omit(data)

        data
      })
)
