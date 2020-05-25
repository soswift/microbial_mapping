
# function that reshapes a data frame into the correct format for plotly sunburst
# https://stackoverflow.com/questions/57395424/how-to-format-data-for-plotly-sunburst-diagram
as.sunburstDF <- function(DF, valueCol = NULL) {
  require(data.table)
  # TEST CASE
  DF =troph_data[c(1,3,2,4)]
  valueCol = "tot"
  ###############
  
  
  DT <- data.table(DF, stringsAsFactors = FALSE)
  DT[, root := "total"]
  setcolorder(DT, c("root", names(DF)))
  
  hierarchyList <- list()
  if (!is.null(valueCol)) {
    setnames(DT, valueCol, "values", skip_absent = TRUE)
  }
  hierarchyCols <- setdiff(names(DT), "values")
  
  for (i in seq_along(hierarchyCols)) {
      currentCols <- names(DT)[1:i]
    if (is.null(valueCol)) {
      currentDT <-
        unique(DT[, ..currentCols][, values := .N, by = currentCols], by = currentCols)
    } else {
      currentDT <-
        DT[, lapply(.SD, sum, na.rm = TRUE), by = currentCols, .SDcols = "values"]
    }
    setnames(currentDT, length(currentCols), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <-
    rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parentCols <-
    setdiff(names(hierarchyDT), c("labels", "values", valueCol))
  hierarchyDT[, parents := apply(.SD, 1, function(x) {
    fifelse(all(is.na(x)),
            yes = NA_character_,
            no = paste(x[!is.na(x)], sep = ":", collapse = " - "))
  }), .SDcols = parentCols]
  hierarchyDT[, ids := apply(.SD, 1, function(x) {
    paste(x[!is.na(x)], collapse = " - ")
  }), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parentCols) := NULL]
  return(hierarchyDT)
}