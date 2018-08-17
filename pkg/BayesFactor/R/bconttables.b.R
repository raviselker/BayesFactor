
bcontTablesClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "bcontTablesClass",
    inherit = bcontTablesBase,
    private = list(
      #### Init + run functions ----
      .init=function() {

        rowVarName <- self$options$rows
        colVarName <- self$options$cols
        layerNames <- self$options$layers
        countsName <- self$options$counts

        freqs <- self$results$freqs
        tests <- self$results$tests
        odds  <- self$results$odds

        data <- private$.cleanData()

        reversed <- rev(layerNames)
        for (i in seq_along(reversed)) {
          layer <- reversed[[i]]
          freqs$addColumn(name=layer, type='text', combineBelow=TRUE)
          tests$addColumn(index=i, name=layer, type='text', combineBelow=TRUE)
          odds$addColumn(index=i, name=layer, type='text', combineBelow=TRUE)
        }

        # add the row column, containing the row variable
        # fill in dots, if no row variable specified

        if ( ! is.null(rowVarName))
          title <- rowVarName
        else
          title <- '.'

        freqs$addColumn(
          name=title,
          title=title,
          type='text')

        # add the column columns (from the column variable)
        # fill in dots, if no column variable specified

        if ( ! is.null(colVarName)) {
          superTitle <- colVarName
          levels <- base::levels(data[[colVarName]])
        }
        else {
          superTitle <- '.'
          levels <- c('.', '.')
        }

        subNames  <- c('[count]', '[pcRow]', '[pcCol]', '[pcTot]')
        subTitles <- c('Observed', '% within row', '% within column', '% of total')
        visible   <- c('TRUE', '(pcRow)', '(pcCol)', '(pcTot)')
        types     <- c('integer', 'number', 'number', 'number')
        formats   <- c('', 'pc', 'pc', 'pc')

        # iterate over the sub rows

        for (j in seq_along(subNames)) {
          subName <- subNames[[j]]
          if (j == 1)
            v <- '(pcRow || pcCol || pcTot)'
          else
            v <- visible[j]

          freqs$addColumn(
            name=paste0('type', subName),
            title='',
            type='text',
            visible=v)
        }

        for (i in seq_along(levels)) {
          level <- levels[[i]]

          for (j in seq_along(subNames)) {
            subName <- subNames[[j]]
            freqs$addColumn(
              name=paste0(i, subName),
              title=level,
              superTitle=superTitle,
              type=types[j],
              format=formats[j],
              visible=visible[j])
          }
        }

        # add the Total column

        freqs$addColumn(
          name='.total[count]',
          title='Total',
          type='integer')

        # populate the first column with levels of the row variable

        values <- list()
        for (i in seq_along(subNames))
          values[[paste0('type', subNames[i])]] <- subTitles[i]

        rows <- private$.grid(data=data, incRows=TRUE)

        for (i in seq_len(nrow(rows))) {
          for (name in dimnames(rows)[[2]]) {
            value <- as.character(rows[i, name])
            if (value == '.total')
              value <- 'Total'
            values[[name]] <- value
          }
          key <- paste0(rows[i,], collapse='`')
          freqs$addRow(rowKey=key, values=values)

          if (i == 1)
            freqs$addFormat(rowNo=i, 1, jmvcore::Cell.BEGIN_GROUP)
          else if (i == nrow(rows) - 1)
            freqs$addFormat(rowNo=i, 1, jmvcore::Cell.END_GROUP)
          else if (i == nrow(rows))
            freqs$addFormat(rowNo=i, 1, jmvcore::Cell.BEGIN_END_GROUP)
        }

        rows <- private$.grid(data=data, incRows=FALSE)
        values <- list()

        if (length(rows) == 0) {

          tests$addRow(rowKey=1, values=list())
          odds$addRow(rowKey=1, values=list())

        } else {

          for (i in seq_len(nrow(rows))) {

            for (name in dimnames(rows)[[2]]) {
              value <- as.character(rows[i, name])
              if (value == '.total')
                value <- 'Total'
              values[[name]] <- value
            }

            tests$addRow(rowKey=i, values=values)
            odds$addRow(rowKey=i, values=values)
          }
        }

        ciText <- paste0(self$options$ciWidth, '% Confidence Intervals')
        odds$getColumn('cil[lo]')$setSuperTitle(ciText)
        odds$getColumn('ciu[lo]')$setSuperTitle(ciText)
        odds$getColumn('cil[f]')$setSuperTitle(ciText)
        odds$getColumn('ciu[f]')$setSuperTitle(ciText)

      },
      .run=function() {

        rowVarName <- self$options$rows
        colVarName <- self$options$cols
        countsName <- self$options$counts

        if (is.null(rowVarName) || is.null(colVarName))
          return()

        data <- private$.cleanData()

        if (nlevels(data[[rowVarName]]) < 2)
          jmvcore::reject("Row variable '{}' contains less than 2 levels", code='', rowVarName)
        if (nlevels(data[[colVarName]]) < 2)
          jmvcore::reject("Column variable '{}' contains less than 2 levels", code='', colVarName)

        if ( ! is.null(countsName)) {
          countCol <- data[[countsName]]
          if (any(countCol < 0, na.rm=TRUE))
            jmvcore::reject('Counts may not be negative')
          if (any(is.infinite(countCol)))
            jmvcore::reject('Counts may not be infinite')
        }

        freqs <- self$results$freqs
        tests <- self$results$tests
        odds  <- self$results$odds

        freqRowNo <- 1
        othRowNo <- 1

        mats <- private$.matrices(data)

        nRows  <- base::nlevels(data[[rowVarName]])
        nCols  <- base::nlevels(data[[colVarName]])
        nCells <- nRows * nCols

        ciWidth <- self$options$ciWidth / 100

        for (mat in mats) {

          suppressWarnings({

            test <- try(extractBF(contingencyTableBF(mat, sampleType='poisson')))
            n <- sum(mat)

          }) # suppressWarnings

          total <- sum(mat)
          colTotals <- apply(mat, 2, sum)

          for (rowNo in seq_len(nRows)) {

            values <- mat[rowNo,]
            rowTotal <- sum(values)

            pcRow <- values / rowTotal

            values <- as.list(values)
            names(values) <- paste0(1:nCols, '[count]')
            values[['.total[count]']] <- rowTotal

            pcRow <- as.list(pcRow)
            names(pcRow) <- paste0(1:nCols, '[pcRow]')

            pcCol <- as.list(mat[rowNo,] / colTotals)
            names(pcCol) <- paste0(1:nCols, '[pcCol]')

            pcTot <- as.list(mat[rowNo,] / total)
            names(pcTot) <- paste0(1:nCols, '[pcTot]')

            values <- c(values, pcRow, pcCol, pcTot)

            freqs$setRow(rowNo=freqRowNo, values=values)
            freqRowNo <- freqRowNo + 1
          }

          values <- apply(mat, 2, sum)
          rowTotal <- sum(values)
          values <- as.list(values)
          names(values) <- paste0(1:nCols, '[count]')
          values[['.total[count]']] <- rowTotal

          pcRow <- apply(mat, 2, sum) / rowTotal
          pcRow <- as.list(pcRow)
          names(pcRow) <- paste0(1:nCols, '[pcRow]')

          pcCol <- rep(1, nCols)
          pcCol <- as.list(pcCol)
          names(pcCol) <- paste0(1:nCols, '[pcCol]')

          pcTot <- apply(mat, 2, sum) / total
          pcTot <- as.list(pcTot)
          names(pcTot) <- paste0(1:nCols, '[pcTot]')

          values <- c(values, pcRow, pcCol, pcTot)

          freqs$setRow(rowNo=freqRowNo, values=values)
          freqRowNo <- freqRowNo + 1

          # populate chi squared table

          if (base::inherits(test, 'try-error'))
            values <- list(
              `value[bf]`=NaN,
              `df[bf]`='',
              `p[bf]`='',
              `value[N]`=n)
          else
            values <- list(
              `value[bf]`=unname(test$statistic),
              `df[bf]`=unname(test$parameter),
              `p[bf]`=unname(test$p.value),
              `value[N]`=n)
          tests$setRow(rowNo=othRowNo, values=values)

          if (TRUE) {  # should test 2x2


          } else {
            odds$setRow(rowNo=othRowNo, list(
              `v[lo]`=NaN, `cil[lo]`='', `ciu[lo]`='',
              `v[f]`=NaN, `cil[f]`='', `ciu[f]`=''))
            odds$addFootnote(rowNo=othRowNo, 'v[lo]', 'Available for 2x2 tables only')
            odds$addFootnote(rowNo=othRowNo, 'v[f]', 'Available for 2x2 tables only')
          }

          othRowNo <- othRowNo + 1
        }

      },

      #### Helper functions ----
      .cleanData = function() {

        data <- self$data

        rowVarName <- self$options$rows
        colVarName <- self$options$cols
        layerNames <- self$options$layers
        countsName <- self$options$counts

        if ( ! is.null(rowVarName))
          data[[rowVarName]] <- as.factor(data[[rowVarName]])
        if ( ! is.null(colVarName))
          data[[colVarName]] <- as.factor(data[[colVarName]])
        for (layerName in layerNames)
          data[[layerName]] <- as.factor(data[[layerName]])
        if ( ! is.null(countsName))
          data[[countsName]] <- jmvcore::toNumeric(data[[countsName]])

        data
      },
      .matrices=function(data) {

        matrices <- list()

        rowVarName <- self$options$rows
        colVarName <- self$options$cols
        layerNames <- self$options$layers
        countsName <- self$options$counts

        if (length(layerNames) == 0) {

          subData <- jmvcore::select(data, c(rowVarName, colVarName))

          if (is.null(countsName))
            .COUNTS <- rep(1, nrow(subData))
          else
            .COUNTS <- jmvcore::toNumeric(data[[countsName]])

          matrices <- list(ftable(xtabs(.COUNTS ~ ., data=subData)))

        } else {

          layerData <- jmvcore::select(data, layerNames)
          dataList <- do.call(split, list(data, layerData))

          tables <- lapply(dataList, function(x) {

            xTemp <- jmvcore::select(x, c(rowVarName, colVarName))

            if (is.null(countsName))
              .COUNTS <- rep(1, nrow(xTemp))
            else
              .COUNTS <- jmvcore::toNumeric(x[[countsName]])

            ftable(xtabs(.COUNTS ~ ., data=xTemp))
          })

          rows <- private$.grid(data=data, incRows=FALSE)

          expand <- list()

          for (layerName in layerNames)
            expand[[layerName]] <- c(base::levels(data[[layerName]]))

          tableNames <- rev(expand.grid(expand))

          matrices <- list()
          for (i in seq_along(rows[,1])) {

            indices <- c()
            for (j in seq_along(tableNames[,1])) {

              row <- as.character(unlist((rows[i,])))
              tableName <- as.character(unlist(tableNames[j,]))

              if (all(row == tableName | row == '.total'))
                indices <- c(indices, j)
            }

            matrices[[i]] <- Reduce("+", tables[indices])
          }

        }

        matrices
      },
      .grid=function(data, incRows=FALSE) {

        rowVarName <- self$options$rows
        layerNames <- self$options$layers

        expand <- list()

        if (incRows) {
          if (is.null(rowVarName))
            expand[['.']] <- c('.', '. ', 'Total')
          else
            expand[[rowVarName]] <- c(base::levels(data[[rowVarName]]), '.total')
        }

        for (layerName in layerNames)
          expand[[layerName]] <- c(base::levels(data[[layerName]]), '.total')

        rows <- rev(expand.grid(expand))

        rows
      })
)
