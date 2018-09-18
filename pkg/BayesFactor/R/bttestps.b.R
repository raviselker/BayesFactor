
bttestPSClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bttestPSClass",
    inherit = bttestPSBase,
    private = list(
        #### Member variables ----
        postPointEst = list(),
        hyp = NULL,
        nullInterval = NULL,

        #### Init + run functions ----
        .init = function() {

            private$.initTTestTable()
            private$.initDescTable()
            private$.initPlots()

        },
        .run = function() {

            if (length(self$options$pairs) == 0 || is.null(self$options$pairs[[1]]$i2))
                return()

            data <- private$.cleanData()

            results <- private$.compute(data)

            private$.populateTTestTable(results)
            private$.populateDescTable(results)

            if (self$options$robust)
                private$.prepareRobust(results, data)

            if (self$options$pp)
                private$.preparePP(results)

        },

        #### Compute results ----
        .compute = function(data) {

            pairs <- self$options$pairs

            if (self$options$hypothesis == 'oneGreater') {
                alternative <- 'greater'
                nullInterval <- c(0, Inf)
            }
            else if (self$options$hypothesis == 'twoGreater') {
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
            for (i in seq_along(pairs)) {

                pair <- pairs[[i]]

                if (is.null(pair$i2)) {
                    r[[i]] <- list()
                    next()
                }

                dataA <- data.frame(
                    x1 = jmvcore::toNumeric(data[[pair$i1]]),
                    x2 = jmvcore::toNumeric(data[[pair$i2]]))

                if (self$options$miss != 'listwise')
                    dataA <- na.omit(dataA)

                tres <- t.test(x = dataA$x1, y = dataA$x2, var.equal = TRUE, alternative = alternative, paired = TRUE)
                bfres <- ttestBF(x = dataA$x1, y = dataA$x2, nullInterval = nullInterval, rscale = self$options$bfPrior, paired = TRUE)
                bfo <- extractBF(bfres)

                m <- apply(dataA, 2, jmvcore::tryNaN(mean))
                v <- apply(dataA, 2, jmvcore::tryNaN(var))
                n <- apply(dataA, 2, length)
                sd <- sqrt(v)

                n[is.na(n)] <- 0
                m[is.na(m)] <- NaN
                sd[is.na(sd)] <- NaN

                desc <- list(m=m, n=n, sd=sd)

                r[[i]] <- list(tres=tres, bfres=bfres, bfo=bfo, desc=desc)
            }

            return(r)
        },

        #### Init tables/plots functions ----
        .initTTestTable = function() {

            table <- self$results$ttest

            ci <- paste0(self$options$ciWidth, '% Credible Interval')
            table$getColumn('cil')$setSuperTitle(ci)
            table$getColumn('ciu')$setSuperTitle(ci)

            hypothesis <- self$options$hypothesis

            if (hypothesis == 'oneGreater')
                table$setNote("hyp", "H\u2090 Measure 1 > Measure 2")
            else if (hypothesis == 'twoGreater')
                table$setNote("hyp", "H\u2090 Measure 1 < Measure 2")
            else
                table$setNote("hyp", NULL)

            pairs <- self$options$pairs

            for (i in seq_along(pairs)) {

                pair <- pairs[[i]]

                table$setRow(rowKey=pair, list(
                    `var1`=pair$i1,
                    `var2`=pair$i2))
            }
        },
        .initDescTable = function() {

            table <- self$results$desc

            pairs <- self$options$pairs

            for (i in seq_along(pairs)) {

                pair <- pairs[[i]]

                row1Key <- paste0(pair$i1, i)
                row2Key <- paste0(pair$i2, i)

                table$addRow(row1Key, list(name=pair$i1))
                table$addFormat(rowKey=row1Key, col=1, jmvcore::Cell.BEGIN_GROUP)
                table$addRow(row2Key, list(name=pair$i2))
                table$addFormat(rowKey=row2Key, col=2, jmvcore::Cell.END_GROUP)

            }
        },
        .initPlots = function() {

            pairs <- self$options$pairs
            plots <- self$results$plots

            for (i in seq_along(pairs)) {

                pair <- pairs[[i]]
                plots$get(pair)$setTitle(paste0(pair, collapse=' - '))

            }
        },

        #### Populate tables functions ----
        .populateTTestTable = function(results) {

            table <- self$results$ttest
            lls <- private$postPointEst

            pairs <- self$options$pairs

            for (i in seq_along(pairs)) {

                r <- results[[i]]

                if (length(r) == 0)
                    next()

                row <- list(
                    bf10 = r$bfo$bf[1],
                    bf01 = 1 / r$bfo$bf[1],
                    err = r$bfo$error[1],
                    md = r$tres$estimate
                )

                table$setRow(rowNo=i, row)

                if (self$options$effectSize || self$options$ci) {

                    private$.checkpoint()

                    if (is.null(lls[[ paste0(pairs[[i]], collapse = '') ]])) {
                        ll <- getPostPointEst(r, type = 'PS', self$options$bfPrior)
                        lls[[ paste0(pairs[[i]], collapse = '') ]] <- ll
                    } else {
                        ll <- lls[[ paste0(pairs[[i]], collapse = '') ]]
                    }

                    ciWidth <- self$options$ciWidth

                    es <- ll[[1]]
                    cil <- qtScaled(1 - (ciWidth / 200 + 0.5), mean=ll[1], sd=ll[2], df=ll[3], hyp=private$hyp)
                    ciu <- qtScaled(ciWidth / 200 + 0.5, mean=ll[1], sd=ll[2], df=ll[3], hyp=private$hyp)
                    table$setRow(rowNo=i, list(es=es, cil=cil, ciu=ciu))

                }
            }
        },
        .populateDescTable = function(results) {

            table <- self$results$desc

            pairs <- self$options$pairs

            iter <- 1

            for (i in seq_along(pairs)) {

                r <- results[[i]]

                if (length(r) == 0)
                    next()

                for (j in 1:2) {

                    row <- list(
                        num = r$desc$n[j],
                        mean = r$desc$m[j],
                        sd = r$desc$sd[j]
                    )

                    table$setRow(rowNo=iter, row)

                    iter <- iter + 1
                }
            }
        },

        #### Plot functions ----
        .prepareRobust = function(results, data) {

            plots <- self$results$plots

            pairs <- self$options$pairs
            nSteps <- 400
            steps <- c(seq(0.001, 2, length.out = nSteps), self$options$bfPrior)

            bfType <- self$options$bfType

            for (i in seq_along(pairs)) {

                image <- plots$get(key=pairs[[i]])$robust

                r <- results[[i]]

                pair <- pairs[[i]]

                if (is.null(pair$i2)) {
                    r[[i]] <- list()
                    next()
                }

                dataA <- data.frame(
                    x1 = jmvcore::toNumeric(data[[pair$i1]]),
                    x2 = jmvcore::toNumeric(data[[pair$i2]]))

                if (self$options$miss != 'listwise')
                    dataA <- na.omit(dataA)

                BF <- numeric(length(steps))
                for (i in seq_along(steps))
                    BF[i] <- BayesFactor::ttestBF(x=dataA$x1, y=dataA$x2, nullInterval = private$nullInterval, rscale = steps[i], paired = TRUE)@bayesFactor$bf

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
            pairs <- self$options$pairs

            lls <- private$postPointEst

            for (i in seq_along(pairs)) {

                image <- plots$get(key=pairs[[i]])$pp

                r <- results[[i]]

                if (is.null(lls[[ paste0(pairs[[i]], collapse = '') ]])) {
                    ll <- getPostPointEst(r, type = 'PS', self$options$bfPrior)
                    lls[[ paste0(pairs[[i]], collapse = '') ]] <- ll
                } else {
                    ll <- lls[[ paste0(pairs[[i]], collapse = '') ]]
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
