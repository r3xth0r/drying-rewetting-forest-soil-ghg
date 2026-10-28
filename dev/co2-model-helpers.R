group_key <- function(site, depth_cm, treatment) {
  stringr::str_c(site, depth_cm, treatment, sep = "_")
}

select_group <- function(data, site, depth_cm, treatment) {
  dplyr::filter(
    data,
    Site == .env$site,
    Depth_cm == .env$depth_cm,
    Treatment == .env$treatment
  )
}
