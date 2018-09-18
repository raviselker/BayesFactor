
bttestISClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bttestISClass",
    inherit = bttestISBase,
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

            if (is.null(self$options$group) || length(self$options$vars) == 0)
                return()

            data <- private$.cleanData()
            private$.errorCheck(data)

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

            group <- self$options$group
            groupLevels <- base::levels(data[[group]])

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
            for (var in self$options$vars) {

                dataA <- data.frame(
                    dep = jmvcore::toNumeric(data[[var]]),
                    group = data[[group]])

                if (self$options$miss != 'listwise')
                    dataA <- na.omit(dataA)

                tres <- t.test(formula = dep ~ group, data = dataA, var.equal = TRUE, alternative = alternative)
                bfres <- ttestBF(formula = dep ~ group, data = dataA, nullInterval = nullInterval, rscale = self$options$bfPrior)
                bfo <- extractBF(bfres)

                m <- tapply(dataA$dep, dataA$group, function(x) jmvcore::tryNaN(mean(x)))
                v <- tapply(dataA$dep, dataA$group, function(x) jmvcore::tryNaN(var(x)))
                n <- tapply(dataA$dep, dataA$group, length)
                sd <- sqrt(v)

                n[is.na(n)] <- 0
                m[is.na(m)] <- NaN
                sd[is.na(sd)] <- NaN

                desc <- list(m=m, n=n, sd=sd, groupLevels=groupLevels)

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
            groupName <- self$options$group

            groups <- NULL
            if ( ! is.null(groupName))
                groups <- base::levels(self$data[[groupName]])
            if (length(groups) != 2)
                groups <- c('Group 1', 'Group 2')

            if (hypothesis == 'oneGreater')
                table$setNote("hyp", jmvcore::format("H\u2090 {} > {}", groups[1], groups[2]))
            else if (hypothesis == 'twoGreater')
                table$setNote("hyp", jmvcore::format("H\u2090 {} < {}", groups[1], groups[2]))
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
                    md = r$tres$estimate[1] - r$tres$estimate[2]
                )

                table$setRow(rowKey=var, row)

                if (self$options$effectSize || self$options$ci) {

                    private$.checkpoint()

                    if (is.null(lls[[var]])) {
                        ll <- getPostPointEst(r, type = 'IS', self$options$bfPrior)
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

            private$.checkpoint()
            table$setStatus('complete')
        },
        .populateDescTable = function(results) {

            table <- self$results$desc

            for (var in self$options$vars) {

                r <- results[[var]]

                row <- list(
                    "dep" = var,
                    "group[1]" = r$desc$groupLevels[1],
                    "group[2]" = r$desc$groupLevels[2],
                    "num[1]" = r$desc$n[1],
                    "num[2]" = r$desc$n[2],
                    "mean[1]" = r$desc$m[1],
                    "mean[2]" = r$desc$m[2],
                    "sd[1]" = r$desc$sd[1],
                    "sd[2]" = r$desc$sd[2]
                )

                table$setRow(rowKey=var, row)
            }

            private$.checkpoint()
            table$setStatus('complete')
        },

        #### Plot functions ----
        .prepareRobust = function(results, data) {

            plots <- self$results$plots

            vars <- self$options$vars
            group <- self$options$group
            nSteps <- 400
            steps <- c(seq(0.001, 2, length.out = nSteps), self$options$bfPrior)

            bfType <- self$options$bfType

            for (var in vars) {

                image <- plots$get(key=var)$robust

                r <- results[[var]]

                dataA <- data.frame(
                    dep = jmvcore::toNumeric(data[[var]]),
                    group = data[[group]])

                if (self$options$miss != 'listwise')
                    dataA <- na.omit(dataA)

                BF <- numeric(length(steps))
                for (i in seq_along(steps))
                    BF[i] <- BayesFactor::ttestBF(formula = dep ~ group, data = dataA, nullInterval = private$nullInterval, rscale = steps[i])@bayesFactor$bf

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
                    ll <- getPostPointEst(r, type = 'IS', self$options$bfPrior)
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
        },
        .errorCheck = function(data) {

            group <- self$options$group
            groupLevels <- base::levels(data[[group]])
            if (length(groupLevels) != 2)
                jmvcore::reject("Grouping variable '{a}' must have exactly 2 levels", code="grouping_var_must_have_2_levels", a=group)


        })
)
