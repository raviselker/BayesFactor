
bttestISClass <- if (requireNamespace('jmvcore')) R6::R6Class(
  "bttestISClass",
  inherit = bttestISBase,
  private = list(
    .init = function() {

      ci <- paste0(self$options$ciWidth, '% Credible Interval')
      self$results$ttest$getColumn('cil')$setSuperTitle(ci)
      self$results$ttest$getColumn('ciu')$setSuperTitle(ci)
      
      hypothesis <- self$options$hypothesis
      groupName <- self$options$group
      
      groups <- NULL
      if ( ! is.null(groupName))
        groups <- base::levels(self$data[[groupName]])
      if (length(groups) != 2)
        groups <- c('Group 1', 'Group 2')
      
      if (hypothesis == 'oneGreater')
        self$results$ttest$setNote("hyp", jmvcore::format("H\u2090 {} > {}", groups[1], groups[2]))
      else if (hypothesis == 'twoGreater')
        self$results$ttest$setNote("hyp", jmvcore::format("H\u2090 {} < {}", groups[1], groups[2]))
      else
        self$results$ttest$setNote("hyp", NULL)

    },
    .run = function() {
      
      if (is.null(self$options$group))
        return()
      if (length(self$options$vars) == 0)
        return()
      
      if (self$options$hypothesis == 'oneGreater') {
        alternative <- 'less'
        nullInterval <- c(0, Inf)
      }
      else if (self$options$hypothesis == 'twoGreater') {
        alternative <- 'greater'
        nullInterval <- c(-Inf, 0)
      }
      else {
        alternative <- 'two.sided'
        nullInterval <- NULL
      }
      
      data <- self$data
      if (self$options$miss == 'listwise')
        data <- na.omit(data)
      
      group <- self$options$group
      groupLevels <- base::levels(data[[group]])
      if (length(groupLevels) != 2)
        jmvcore::reject("Grouping variable '{a}' must have exactly 2 levels", code="grouping_var_must_have_2_levels", a=groupVarName)

              
      for (var in self$options$vars) {
        
        dataA <- data.frame(
          dep=jmvcore::toNumeric(data[[var]]),
          group=data[[group]])
        
        if (self$options$miss != 'listwise')
          dataA <- na.omit(dataA)
        
        #if (is.factor(dataTTest$dep))
        #  res <- createError('Variable is not numeric')
        #else if (any(is.infinite(dataTTest$dep)))
        #  res <- createError('Variable contains infinite values')
        #else
        #  res <- try(t.test(dep ~ group, data=dataTTest, var.equal=TRUE, paired=FALSE, alternative=Ha, conf.level=confInt), silent=TRUE)
        
        tres <- t.test(
          formula = dep ~ group,
          data = dataA,
          var.equal = TRUE,
          alternative = alternative)
        
        bfres <- ttestBF(
          formula = dep ~ group,
          data = dataA,
          nullInterval = nullInterval,
          rscale = self$options$bfPrior)
        
        bfo <- extractBF(bfres)
        
        m <- tapply(dataA$dep, dataA$group, function(x) jmvcore::tryNaN(mean(x)))
        v <- tapply(dataA$dep, dataA$group, function(x) jmvcore::tryNaN(var(x)))
        n <- tapply(dataA$dep, dataA$group, length)
        se <- sqrt(v/n)
        sd <- sqrt(v)
        pooledVar <- jmvcore::tryNaN(((n[1]-1)*v[1]+(n[2]-1)*v[2])/(n[1]+n[2]-2))
        sediff <- jmvcore::tryNaN(sqrt((pooledVar/n[1])+(pooledVar/n[2])))
        
        n[is.na(n)] <- 0
        m[is.na(m)] <- NaN
        se[is.na(se)] <- NaN
        sd[is.na(sd)] <- NaN
        sediff[is.na(sediff)] <- NaN
        pooledVar[is.na(pooledVar)] <- NaN
        
        row <- list(
          t = tres$statistic,
          df = tres$parameter,
          bf10=bfo$bf[1],
          bf01=1/bfo$bf[1],
          err=bfo$error[1],
          md=tres$estimate[1]-tres$estimate[2],
          sed=sediff
        )
        
        self$results$ttest$setRow(rowKey=var, row)
        
        #result <- extractBF(result)
        
        self$results$desc$setRow(rowKey=var, list(
          "dep"=var,
          "group[1]"=groupLevels[1],
          "group[2]"=groupLevels[2],
          "num[1]"=n[1],
          "num[2]"=n[2],
          "mean[1]"=m[1],
          "mean[2]"=m[2],
          "sd[1]"=sd[1],
          "sd[2]"=sd[2],
          "se[1]"=se[1],
          "se[2]"=se[2]))
        
        
      }
    })
)
