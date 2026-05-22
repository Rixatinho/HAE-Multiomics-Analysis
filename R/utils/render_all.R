#!/usr/bin/env Rscript
# =============================================================================
# render_all.R -- Orchestrator for Modular Panel Rendering
# =============================================================================
# Usage:
#   conda run -n multiomics Rscript render_all.R --all
#   conda run -n multiomics Rscript render_all.R --figure 2
#   conda run -n multiomics Rscript render_all.R --panel Fig2a_volcano_TC
#   conda run -n multiomics Rscript render_all.R --list
# =============================================================================

# Robust script-dir detection (works with Rscript and source())
.get_script_dir <- function() {
  # Method 1: commandArgs (Rscript invocation)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }
  # Method 2: sys.frame ofile (source() invocation)
  for (i in seq_len(sys.nframe())) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(normalizePath(dirname(ofile)))
  }
  # Fallback: working directory
  return(getwd())
}
SCRIPT_DIR <- .get_script_dir()

# --- Source shared modules (order matters) ---
source(file.path(SCRIPT_DIR, "nc_theme.R"))
source(file.path(SCRIPT_DIR, "panel_factories.R"))
source(file.path(SCRIPT_DIR, "figure_layouts.R"))

PANELS_DIR <- file.path(SCRIPT_DIR, "panels")

# =============================================================================
# CLI argument parsing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || args[1] == "--help") {
  cat("render_all.R -- Panel Orchestrator\n\n")
  cat("Usage:\n")
  cat("  --all            Render all panel scripts in panels/\n")
  cat("  --figure N       Render panels for Figure N (e.g., --figure 2)\n")
  cat("  --suppfig N      Render panels for SuppFig N (e.g., --suppfig 03)\n")
  cat("  --panel NAME     Render a single panel (e.g., --panel Fig2a_volcano_TC)\n")
  cat("  --list           List all registered panels and their slots\n")
  cat("  --review NAME    Print review card for a panel\n")
  quit(status = 0)
}

# =============================================================================
# Panel discovery: find .R scripts in panels/ directory
# =============================================================================
discover_panels <- function(pattern = NULL) {
  scripts <- list.files(PANELS_DIR, pattern = "\\.R$", full.names = TRUE)
  if (!is.null(pattern)) {
    scripts <- scripts[grepl(pattern, basename(scripts), ignore.case = TRUE)]
  }
  scripts
}

# =============================================================================
# Render a single panel script (source it in current session)
# =============================================================================
render_panel <- function(script_path) {
  panel_name <- gsub("\\.R$", "", basename(script_path))
  slot <- get_panel_slot(paste0(panel_name, ".pdf"))

  cat(sprintf("\n--- %s ---\n", panel_name))
  if (!is.null(slot)) {
    cat(sprintf("  Slot: %s  |  %d x %d mm\n", slot$figure, slot$w, slot$h))
  }

  t0 <- proc.time()
  tryCatch(
    source(script_path, local = new.env(parent = globalenv())),
    error = function(e) cat("  ERROR:", e$message, "\n")
  )
  dt <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  Done in %.1fs\n", dt))
}

# =============================================================================
# Dispatch
# =============================================================================
mode <- args[1]

if (mode == "--list") {
  all_panels <- list_all_panels()
  cat(sprintf("Registered panels: %d\n\n", length(all_panels)))
  for (p in all_panels) {
    s <- get_panel_slot(p)
    script_exists <- file.exists(file.path(PANELS_DIR, gsub("\\.pdf$", ".R", p)))
    status <- if (script_exists) "OK" else "--"
    cat(sprintf("  [%s] %-45s %3dx%3d mm  (%s)\n", status, p, s$w, s$h, s$figure))
  }

} else if (mode == "--panel") {
  if (length(args) < 2) stop("Usage: --panel NAME (without .R extension)")
  target <- args[2]
  script <- file.path(PANELS_DIR, paste0(target, ".R"))
  if (!file.exists(script)) stop("Panel script not found: ", script)
  render_panel(script)

} else if (mode == "--figure") {
  if (length(args) < 2) stop("Usage: --figure N")
  fig_key <- paste0("Figure_", args[2])
  if (!fig_key %in% names(LAYOUTS)) stop("Unknown figure: ", fig_key)
  panels_in_fig <- unlist(lapply(LAYOUTS[[fig_key]]$rows, function(row)
    sapply(row, function(e) gsub("\\.pdf$", "", e$name))))
  scripts <- file.path(PANELS_DIR, paste0(panels_in_fig, ".R"))
  scripts <- scripts[file.exists(scripts)]
  if (length(scripts) == 0) {
    cat("No panel scripts found for", fig_key, "\n")
  } else {
    cat(sprintf("Rendering %d panels for %s\n", length(scripts), fig_key))
    for (s in scripts) render_panel(s)
  }

} else if (mode == "--suppfig") {
  if (length(args) < 2) stop("Usage: --suppfig N (e.g., 03)")
  fig_key <- paste0("SuppFig_", sprintf("%02d", as.integer(args[2])))
  if (!fig_key %in% names(LAYOUTS)) stop("Unknown supp figure: ", fig_key)
  panels_in_fig <- unlist(lapply(LAYOUTS[[fig_key]]$rows, function(row)
    sapply(row, function(e) gsub("\\.pdf$", "", e$name))))
  scripts <- file.path(PANELS_DIR, paste0(panels_in_fig, ".R"))
  scripts <- scripts[file.exists(scripts)]
  if (length(scripts) == 0) {
    cat("No panel scripts found for", fig_key, "\n")
  } else {
    cat(sprintf("Rendering %d panels for %s\n", length(scripts), fig_key))
    for (s in scripts) render_panel(s)
  }

} else if (mode == "--all") {
  scripts <- discover_panels()
  if (length(scripts) == 0) {
    cat("No panel scripts found in", PANELS_DIR, "\n")
  } else {
    cat(sprintf("Rendering ALL %d panel scripts\n", length(scripts)))
    t0 <- proc.time()
    for (s in scripts) render_panel(s)
    dt <- (proc.time() - t0)["elapsed"]
    cat(sprintf("\n=== All done. Total: %.1fs ===\n", dt))
  }

} else if (mode == "--review") {
  if (length(args) < 2) stop("Usage: --review NAME")
  target <- args[2]
  slot <- get_panel_slot(paste0(target, ".pdf"))
  if (is.null(slot)) {
    cat("Panel not found in layouts:", target, "\n")
  } else {
    cat("## Review Card:", target, "\n")
    cat("- Figure:", slot$figure, "\n")
    cat("- Slot:", slot$w, "x", slot$h, "mm\n")
    cat("- Font: Arial (unified)\n")
    cat("- Theme: theme_nc (theme_bw base, rectangular border)\n")
    script <- file.path(PANELS_DIR, paste0(target, ".R"))
    if (file.exists(script)) {
      n_lines <- length(readLines(script))
      cat("- Script:", n_lines, "lines\n")
    } else {
      cat("- Script: NOT YET CREATED\n")
    }
  }

} else {
  cat("Unknown mode:", mode, "\n")
  cat("Run with --help for usage.\n")
}
