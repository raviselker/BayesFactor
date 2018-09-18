
bttestOneSClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bttestOneSClass",
    inherit = bttestOneSBase,
    private = list(
        #### Member variables ----
        postPointEst = list(),
        hyp = NULL,
        nullInterval = NULL,

        #### Init + run functions ----
        .init = function() {

            private$.initTTestTable()

        },
        .run = function() {

            if (is.null(self$options$vars))
                return()

            data <- private$.cleanData()

            results <- private$.compute(data)

            private$.populateDescTable(results)
            private$.populateTTestTable(results)

            if (self$options$robust)
                private$.prepareRobust(results, data)

            if (self$options$pp)
                private$.preparePP(results)

        },

        #### Compute results ----
        .compute = function(data) {

            vars <- self$options$vars
            testValue <- self$options$testValue

            if (self$options$hypothesis == 'greater') {
                alternative <- 'greater'
                nullInterval <- c(0, Inf)
            }
            else if (self$options$hypothesis == 'less') {
                alternative <- 'less'
                nullInterval <- c(-Inf, 0)
            }
            else {
                alternative <- 'two.sided'
                nullInterval <- NULL
            }

            private$hyp <- alternative
            private$nullInterval <- nullInterval

            r <- list()
            for (var in vars) {

                col <- jmvcore::toNumeric(data[[var]])

                if (self$options$miss != 'listwise')
                    col <- na.omit(col)

                colTest <- col - testValue

                tres <- t.test(x = colTest, var.equal = TRUE, alternative = alternative)
                bfres <- ttestBF(x = colTest, nullInterval = nullInterval, rscale = self$options$bfPrior)
                bfo <- extractBF(bfres)

                m <- jmvcore::tryNaN(mean(col))
                v <- jmvcore::tryNaN(var(col))
                n <- length(col)
                sd <- sqrt(v)

                n[is.na(n)] <- 0
                m[is.na(m)] <- NaN
                sd[is.na(sd)] <- NaN

                desc <- list(m=m, v=v, n=n, sd=sd)

                r[[var]] <- list(tres=tres, bfres=bfres, bfo=bfo, desc=desc)
            }

            return(r)
        },

        #### Init tables/plots functions ----
        .initTTestTable = function() {

            table <- self$results$get('ttest')

            ci <- paste0(self$options$ciWidth, '% Credible Interval')
            table$getColumn('cil')$setSuperTitle(ci)
            table$getColumn('ciu')$setSuperTitle(ci)

            hypothesis <- self$options$hypothesis
            testValue <- self$options$testValue

            if (hypothesis == 'greater')
                table$setNote("hyp", jmvcore::format("H\u2090 Measure > {}", testValue))
            else if (hypothesis == 'less')
                table$setNote("hyp", jmvcore::format("H\u2090 Measure < {}", testValue))
            else
                table$setNote("hyp", NULL)

        },

        #### Populate tables functions ----
        .populateTTestTable = function(results) {

            table <- self$results$ttest
            lls <- private$postPointEst

            for (var in self$options$vars) {

                r <- results[[var]]

                row <- list(
                    bf10 = r$bfo$bf[1],
                    bf01 = 1 / r$bfo$bf[1],
                    err = r$bfo$error[1],
                    md = r$tres$estimate
                )

                table$setRow(rowKey=var, row)

                if (self$options$effectSize || self$options$ci) {

                    private$.checkpoint()

                    if (is.null(lls[[var]])) {
                        ll <- getPostPointEst(r, type = 'OS', self$options$bfPrior)
                        lls[[var]] <- ll
                    } else {
                        ll <- lls[[var]]
                    }

                    ciWidth <- self$options$ciWidth

                    es <- ll[[1]]
                    cil <- qtScaled(1 - (ciWidth / 200 + 0.5), mean=ll[1], sd=ll[2], df=ll[3], hyp=private$hyp)
                    ciu <- qtScaled(ciWidth / 200 + 0.5, mean=ll[1], sd=ll[2], df=ll[3], hyp=private$hyp)
                    table$setRow(rowKey=var, list(es=es, cil=cil, ciu=ciu))

                }
            }
        },
        .populateDescTable = function(results) {

            table <- self$results$desc

            for (var in self$options$vars) {

                r <- results[[var]]

                row <- list(
                    num = r$desc$n,
                    mean = r$desc$m,
                    sd = r$desc$sd
                )

                table$setRow(rowKey=var, row)
            }
        },

        #### Plot functions ----
        .prepareRobust = function(results, data) {

            plots <- self$results$plots

            vars <- self$options$vars
            nSteps <- 400
            steps <- c(seq(0.001, 2, length.out = nSteps), self$options$bfPrior)

            bfType <- self$options$bfType

            for (var in vars) {

                image <- plots$get(key=var)$robust

                r <- results[[var]]

                col <- jmvcore::toNumeric(data[[var]])

                if (self$options$miss != 'listwise')
                    col <- na.omit(col)

                colTest <- col - self$options$testValue

                BF <- numeric(length(steps))
                for (i in seq_along(steps))
                    BF[i] <- BayesFactor::ttestBF(x = colTest, nullInterval = private$nullInterval, rscale = steps[i])@bayesFactor$bf

                if (bfType == 'BF01')
                    BF <- -BF

                df <- data.frame('rscale'=c(0, steps), 'BF'=c(0, BF))
                userBF <- data.frame('rscale'=self$options$bfPrior,
                                     'BF'=ifelse(bfType == 'BF10', log(r$bfo$bf[1]), -log(r$bfo$bf[1])))

                image$setState(list(df=df, userBF=userBF))

            }
        },
        .robust = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            p <- robustPlot(image$state, self$options$bfType, self$options$bfPrior, ggtheme)

            print(p)

            return(TRUE)
        },
        .preparePP = function(results) {

            plots <- self$results$plots

            vars <- self$options$vars

            lls <- private$postPointEst

            for (var in vars) {

                image <- plots$get(key=var)$pp

                r <- results[[var]]

                if (is.null(lls[[var]])) {
                    ll <- getPostPointEst(r, type = 'OS', self$options$bfPrior)
                    lls[[var]] <- ll
                } else {
                    ll <- lls[[var]]
                }

                plotData <- ppPlotPrepare(ll, r$bfo$bf[1], self$options$bfType, self$options$bfPrior, private$hyp)

                image$setState(plotData)
            }
        },
        .pp = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            p <- ppPlot(image$state, self$options$bfPrior, ggtheme)

            print(p)

            return(TRUE)
        },

        #### Helper functions ----
        .cleanData=function() {

            data <- self$data
            if (self$options$miss == 'listwise')
                data <- na.omit(data)

            data
        })
)
