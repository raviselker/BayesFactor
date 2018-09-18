
tScaledLL = function(x, par) {

    -sum(log( dtScaled(x, df=par[3], mean=par[1], sd=par[2]) ))

}

dtScaled = function(x, df, mean = 0, sd = 1, hyp = 'twoSided') {

    # Inspired by the 'dt.scaled' function in the 'metRology' package
    if (hyp == 'less') {

        ifelse(x <= 0,
               (dt((x - mean) / sd, df) / sd) / pt(-mean / sd, df, lower.tail=TRUE),
               0)

    } else if (hyp == 'greater') {

        ifelse(x >= 0,
               (dt((x - mean) / sd, df) / sd) / pt(-mean / sd, df, lower.tail=FALSE),
               0)

    } else {

        dt((x - mean) / sd, df) / sd

    }
}

qtScaled = function(p, df, mean = 0, sd = 1, lower.tail = TRUE, hyp = 'twoSided') {

    # Inspired by the 'qt.scaled' function in the 'metRology' package
    if (hyp == 'less') {

        below0 <- pt((0 - mean) / sd, df, lower.tail=TRUE)
        mean + sd * qt(p * below0, df)

    } else if (hyp == 'greater') {

        below0 <- pt((0 - mean) / sd, df, lower.tail=TRUE)
        mean + sd * qt(below0 + p * (1 - below0), df)

    } else {

        mean + sd * qt(p, df, lower.tail = lower.tail)

    }
}

dCauchy = function(x, scale = 0.707, hyp = 'twoSided') {

    # Inspired by the 'dhalfcauchy' function in the 'LaplacesDemon' package
    if (hyp == 'less') {

        ifelse(x <= 0, 2 * scale / (pi * (x*x + scale*scale)), 0)

    } else if (hyp == 'greater') {

        ifelse(x >= 0, 2 * scale / (pi * (x*x + scale*scale)), 0)

    } else {

        dcauchy(x, scale=scale)

    }
}

getPostPointEst = function(results, type, bfPrior = 0.707, iterations = 10000) {

    t <- results$tres$statistic
    df <- results$tres$parameter
    n1 <- results$desc$n[1]

    if (length(results$desc$n) == 1)
        n2 <- NULL
    else
        n2 <- results$desc$n[2]

    bfres <- BayesFactor::meta.ttestBF(t = t, n1 = n1, n2 = n2, rscale = bfPrior)
    samples <- BayesFactor::posterior(model = bfres, index = 1, iterations = iterations)
    samplesDelta <- as.numeric(samples[,'delta'])

    if (type == 'OS' || type == 'PS') {
        delta <- t * sqrt(1 / n1)
        sigma <- 1 / n1
    } else {
        delta <- t * sqrt((n1 + n2) / (n1 * n2))
        sigma <- sqrt((n1 * n1) / (n1 + n2))
    }

    if (sigma < 0.01)
        sigma <- 0.01

    ll <- optim(par=c(delta, sigma, df), fn=tScaledLL, x=samplesDelta, method='BFGS')$par

    return(ll)
}

ppPlotPrepare = function(ll, bf, bfType = 'BF01', bfPrior = 0.707, hyp = 'twoSided') {

    pointPrior <- dCauchy(0, scale = bfPrior, hyp=hyp)
    pointPosterior <- dtScaled(0, mean=ll[1], sd=ll[2], df=ll[3], hyp=hyp)


    points <- data.frame(x=c(0,0),
                         y=c(pointPrior, pointPosterior),
                         color=factor(c('Prior', 'Posterior'), levels=c('Prior', 'Posterior')))

    bf <- ifelse(bfType == 'BF10', bf, 1 / bf)

    xlim <- c(min(qtScaled(0.01, mean=ll[1], sd=ll[2], df=ll[3], hyp=hyp) - 0.2, -2),
              max(qtScaled(0.99, mean=ll[1], sd=ll[2], df=ll[3], hyp=hyp) + 0.2, 2))

    if (bf > 1000000 || bf < 0.000001)
        bfTitle <- formatC(bf, 3, format = "e")
    else
        bfTitle <- formatC(bf , 3, format = "f")

    if (bfType == 'BF10')
        title <- bquote(BF['10'] == .(bfTitle))
    else
        title <- bquote(BF['01'] == .(bfTitle))

    return(list(ll=ll, points=points, hyp=hyp, title=title, xlim=xlim))
}

ppPlot = function(data, rscale, ggtheme = NULL) {

    ll <- data$ll
    points <- data$points
    hyp <- data$hyp
    xlim <- data$xlim
    title <- data$title

    p <- ggplot2::ggplot(data.frame(x = c(-2, 2)), ggplot2::aes(x)) +
        ggplot2::stat_function(fun = dCauchy, args = list(scale = rscale, hyp=hyp), ggplot2::aes(linetype='dashed'), n=10000) +
        ggplot2::stat_function(fun = dtScaled, args = list(mean=ll[1], sd=ll[2], df=ll[3], hyp=hyp), ggplot2::aes(linetype='solid'), n=10000) +
        ggplot2::geom_point(size=3, data=points, shape=21, ggplot2::aes(x=x, y=y, fill=color)) +
        ggplot2::labs(list(x='Effect size', y='Density', title=title)) +
        ggplot2::scale_linetype_manual(name = 'color', values=c('dashed', 'solid'), labels = c('Prior','Posterior')) +
        ggplot2::xlim(xlim[1], xlim[2]) + ggtheme +
        ggplot2::theme(legend.justification = c(ifelse(ll[1] < 0, 1, 0), 1), legend.key.width=grid::unit(2.5,"line"),
                       legend.position = c(ifelse(ll[1] < 0, 1, 0.05), 1), legend.title=ggplot2::element_blank())

    return(p)
}

robustPlotAxes = function(bf, bfType = 'BF01') {

    labelsY <- c('1/100', '1/30', '1/10', '1/3', '1', '3', '10', '30', '100')

    bs <- log(c(100, 30, 10, 3))
    breaksY <- c(-bs, 0, rev(bs))

    minBoundX <- as.numeric(cut(min(bf) * 1.02, c(-Inf, breaksY, Inf))) - 1
    maxBoundX <- as.numeric(cut(max(bf) * 1.02, c(-Inf, breaksY, Inf)))

    if (minBoundX == 0 | minBoundX > length(breaksY)) {
        minBound <- min(bf) * 1.02
    } else {
        minBound <- breaksY[minBoundX]
    }

    if (maxBoundX > length(breaksY)) {
        maxBound <- max(bf) * 1.02
    } else {
        maxBound <- breaksY[maxBoundX]
    }

    if (minBound > log(1/3)) minBound <- log(1/3)
    if (minBound > min(bf)) minBound <- min(bf) * 1.02
    if (maxBound < log(3)) maxBound <- log(3)
    if (maxBound < max(bf)) maxBound <- max(bf) * 1.02

    ylim <- c(minBound, maxBound)

    if (bfType == 'BF01')
        H1 <- FALSE
    else
        H1 <- TRUE

    labelsCrit <- c(
        paste0("~~Very~strong~H[", ifelse(H1, 0, 1), "]"),
        paste0("~~Strong~H[", ifelse(H1, 0, 1), "]"),
        paste0("~~Moderate~H[", ifelse(H1, 0, 1), "]"),
        paste0("~~Anectodal~H[", ifelse(H1, 0, 1), "]"),
        paste0("~~Anectodal~H[", ifelse(H1, 1, 0), "]"),
        paste0("~~Moderate~H[", ifelse(H1, 1, 0), "]"),
        paste0("~~Strong~H[", ifelse(H1, 1, 0), "]"),
        paste0("~~Very~strong~H[", ifelse(H1, 1, 0), "]")
    )

    breaksTemp <- breaksY

    if (ylim[1] < min(breaksTemp)) {
        breaksTemp <- c(ylim[1], breaksTemp)
        labelsCrit <- c(paste0("~~Extreme~H[", ifelse(H1, 0, 1), "]"), labelsCrit)
    }


    if (ylim[2] > max(breaksTemp)) {
        breaksTemp <- c(breaksTemp, ylim[2])
        labelsCrit <- c(labelsCrit, paste0("~~Extreme~H[", ifelse(H1, 1, 0), "]"))
    }

    breaksCrit <- numeric()
    for (i in 1:(length(breaksTemp) - 1))
        breaksCrit[i] <- (breaksTemp[i] + breaksTemp[i + 1]) / 2

    labelsCrit <- labelsCrit[breaksCrit > ylim[1] & breaksCrit < ylim[2]]
    breaksCrit <- breaksCrit[breaksCrit > ylim[1] & breaksCrit < ylim[2]]

    crit <- data.frame(labels=labelsCrit, breaks=breaksCrit)

    return(list(ylim=ylim, breaksY=breaksY, labelsY=labelsY, crit=crit))
}

robustPlot = function(data, bfType = 'BF01', rscale = 0.707, ggtheme = NULL) {

    axis <- robustPlotAxes(data$df$BF, bfType)
    ylim <- axis$ylim
    breaks <- axis$breaksY
    labels <- axis$labelsY
    crit <- axis$crit

    labels <- labels[breaks >= ylim[1] & breaks <= ylim[2]]
    breaks <- breaks[breaks >= ylim[1] & breaks <= ylim[2]]

    p <- ggplot2::ggplot(data=data$df, ggplot2::aes(x=rscale, y=BF)) +
        ggplot2::geom_hline(yintercept=breaks, linetype='dotted') +
        ggplot2::geom_hline(yintercept=0, linetype='solid') +
        ggplot2::geom_line() +
        ggplot2::geom_point(data=data$userBF, size=3, shape=21, ggplot2::aes(fill=paste0('User prior width: ', round(rscale, 3)))) +
        ggplot2::geom_text(data=crit, ggplot2::aes(y=breaks, label=labels), x=Inf,
                           hjust=0, vjust=.5, size=4.5, parse = TRUE) +
        ggplot2::labs(list(x='Cauchy prior width', y=ifelse(bfType == 'BF10', expression(BF['10']), expression(BF['01'])))) +
        ggplot2::coord_cartesian(ylim=ylim, clip = 'off') +
        ggplot2::scale_y_continuous(breaks=breaks, labels=labels) +
        ggtheme +
        ggplot2::theme(legend.justification = 0.5, legend.position = 'top', legend.title=ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = 14), legend.margin = ggplot2::margin(0, 0, 0, 0),
                       axis.text.y = ggplot2::element_text(hjust = 1), plot.margin = ggplot2::margin(15, 100, 15, 15))

    return(p)
}
