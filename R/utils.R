#' Title
#'
#' @param blob
#' @param type
#'
#' @return
#' @export
#'
#' @examples
blob_extract_1d <- function(blob, type = "double") {
  numbers <- packBits(rawToBits(unlist(blob)), type = type)

  return(numbers)
}

#' Title
#'
#' @param blob
#' @param type
#'
#' @return
#' @export
#'
#' @examples
blob_extract_2d <- function(blob, type = "double") {
  numbers <- packBits(rawToBits(unlist(blob)), type = type)

  df <- data.frame(
    # get every odd entry
    x = numbers[seq(1, length(numbers), by = 2)],
    # get every even entry
    y = numbers[seq(2, length(numbers), by = 2)]
  )

  return(df)
}


#' Extract Table From .moc File
#'
#' @param file file path for the `.moc` file to extract from.
#' @param table
#'
#' @return a [tibble][tibble::tibble-package]
#' @export
#'
#' @examples
extract_table <- function(file, table) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), file)

  df <- RSQLite::dbReadTable(con, name = table) |>
    tibble::as_tibble()

  RSQLite::dbDisconnect(con)

  df
}
