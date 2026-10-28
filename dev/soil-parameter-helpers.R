# Shared soil-parameter analysis helpers

treatment_codes <- c("C", "C+N", "S", "S+N")
alpha <- 0.05

diagnostic_plot <- list(
  mfrow = c(2, 2),
  outer_margin = c(0, 0, 3, 0),
  which = c(1, 2, 3, 5),
  title_line = 1,
  title_font = 2
)

diagnostic_dir <- file.path(dir_interim, "diagnostics")
diagnostic_pdf_width_mm <- 250
diagnostic_pdf_height_mm <- 200
diagnostic_png_width_mm <- 250
diagnostic_png_height_mm <- 200
diagnostic_png_resolution <- 180

clean_text <- function(x) trimws(as.character(x))

check_treatments <- function(x) {
  if (anyNA(x) || any(!x %in% treatment_codes)) {
    stop("Unexpected or missing treatment codes.")
  }
}

as_treatment <- function(x) {
  factor(x, levels = treatment_codes, labels = treatment_levels)
}

strict_numeric <- function(x, column) {
  original <- clean_text(x)
  original[!is.na(original) & original %in% c("", "NA", "NaN")] <- NA_character_

  converted <- suppressWarnings(as.numeric(original))
  bad <- !is.na(original) & is.na(converted)

  if (any(bad)) {
    stop(
      "Non-numeric entries in ", column, ": ",
      paste(unique(original[bad]), collapse = ", ")
    )
  }

  converted
}

make_letters <- function(p_matrix, means, alpha = 0.05) {
  groups <- names(means)[order(means, decreasing = TRUE)]
  p_matrix <- p_matrix[groups, groups, drop = FALSE]

  if (!any(p_matrix[upper.tri(p_matrix)] < alpha)) {
    return(setNames(rep("", length(groups)), groups))
  }

  multcompLetters(
    p_matrix,
    threshold = alpha
  )$Letters
}

collect_results <- function(results, component) {
  imap_dfr(results, function(result, key) {
    result[[component]] |>
      mutate(Model = key, .before = 1)
  })
}

print_section <- function(title, object) {
  cat("\n", title, "\n", sep = "")
  print(object)
}

print_model_heading <- function(key, n_excluded) {
  cat(
    "\n\n", strrep("=", 72), "\n", key, "\n",
    strrep("=", 72), "\n",
    "Missing/non-finite observations excluded: ", n_excluded, "\n",
    sep = ""
  )
}

draw_diagnostics <- function(model, title, hide_annotations = FALSE) {
  par(mfrow = diagnostic_plot$mfrow, oma = diagnostic_plot$outer_margin)

  if (hide_annotations) {
    plot(model, which = diagnostic_plot$which, ask = FALSE, ann = FALSE)
  } else {
    plot(model, which = diagnostic_plot$which, ask = FALSE, sub.caption = "")
  }

  mtext(
    title,
    outer = TRUE,
    line = diagnostic_plot$title_line,
    font = diagnostic_plot$title_font
  )
}

save_diagnostic_pdf <- function(results, filename, hide_annotations = FALSE) {
  if (!save_pdf_figures) {
    return(invisible(NULL))
  }

  ensure_dir(diagnostic_dir)
  pdf(
    file.path(diagnostic_dir, filename),
    width = mm_to_inches(diagnostic_pdf_width_mm),
    height = mm_to_inches(diagnostic_pdf_height_mm)
  )
  on.exit(dev.off(), add = TRUE)

  for (key in names(results)) {
    draw_diagnostics(results[[key]]$model, key, hide_annotations)
  }
}

save_diagnostic_png <- function(result, key) {
  ensure_dir(diagnostic_dir)
  png(
    file.path(
      diagnostic_dir,
      paste0(
        "figure-s09-diagnostic-",
        str_to_lower(str_replace_all(key, "_", "-")), ".png"
      )
    ),
    width = mm_to_inches(diagnostic_png_width_mm),
    height = mm_to_inches(diagnostic_png_height_mm),
    units = "in",
    res = diagnostic_png_resolution,
    bg = "white"
  )
  on.exit(dev.off(), add = TRUE)

  draw_diagnostics(result$model, key)
}
