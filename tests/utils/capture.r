# adapted from answers on https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
capture <- function(expr) {
  wrn <- err <- NULL
  value <- withCallingHandlers(tryCatch(expr,
                                        error=function(e) {
                                          err <<- e
                                          NULL
                                        }),
                                        warning=function(w) {
                                          wrn <<- append(wrn, conditionMessage(w))
                                          invokeRestart("muffleWarning")
                                        })
  list(value=value, warning=wrn, error=err$message)
}
