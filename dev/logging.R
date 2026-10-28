# Minimal lifecycle logging for runnable scripts.

log_script_start <- function(script_name) {
  logger::log_info(paste("🛫 Starting:", script_name))
}

log_script_end <- function(script_name) {
  logger::log_info(paste("🛬 Completed:", script_name))
}
