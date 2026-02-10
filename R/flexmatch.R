library(data.table)
library(MatchIt)
library(ggplot2)

# ==============================================================================
# 1. SCAN: Label Candidates
# ==============================================================================
#' Label Treatment and Control Stack Candidates (Corrected First-Treatment Logic)
#'
#' @description
#' Identifies valid stack candidates.
#' Fix: Strictly enforces that 'first' treated stacks must match the ID's
#' first-ever event time, discarding subsequent events even if they are valid.
#'
#' @export
label_candidates <- function(dt,
                             id_col,
                             time_col,
                             event_col,
                             k_pre,
                             k_post,
                             k_freeze,
                             match_by = -1,
                             treated_type = c("first", "all"),
                             control_type = c("all", "nyt", "nt", "nevert")) {

  treated_type <- match.arg(treated_type)
  control_type <- match.arg(control_type)

  if (!inherits(dt, "data.table")) setDT(dt)

  # 1. Setup & Sorting
  calc_dt <- dt[, .SD, .SDcols = c(id_col, time_col, event_col)]
  setnames(calc_dt, c("id", "period", "event"))

  # Ensure strict ordering by Time
  setorderv(calc_dt, c("id", "period"))

  # --- CRITICAL: Create Time Index (1..N) ---
  # This abstracts "Quarters/Days" into simple integers 1, 2, 3...
  calc_dt[, time_idx := seq_len(.N), by = id]

  # 2. Identify "Ever Treated" & "First Event Time"
  global_stats <- calc_dt[, .(
    is_ever_treated = max(event, na.rm=TRUE) == 1,
    first_treat_idx = if(max(event, na.rm=TRUE)==1) min(time_idx[event==1], na.rm=TRUE) else NA_integer_
  ), by = id]

  calc_dt <- merge(calc_dt, global_stats, by = "id", all.x = TRUE)

  # 3. Vectorized Distance Logic (Using Index)
  calc_dt[, evt_loc := ifelse(event == 1, time_idx, NA_integer_)]
  calc_dt[, last_evt := nafill(evt_loc, type = "locf"), by = id]
  calc_dt[, next_evt := nafill(evt_loc, type = "nocb"), by = id]

  calc_dt[, dist_last := time_idx - last_evt]
  calc_dt[, dist_next := next_evt - time_idx]

  # Dirty Statuses
  calc_dt[, is_dirty_past := !is.na(dist_last) & dist_last >= 0 & dist_last <= (k_post + k_freeze)]
  calc_dt[, is_dirty_future := !is.na(dist_next) & dist_next <= k_post]

  # Continuity
  calc_dt[, has_history := !is.na(shift(time_idx, n = k_pre, type="lag")), by = id]
  calc_dt[, has_future  := !is.na(shift(time_idx, n = k_post, type="lead")), by = id]

  # Valid Treated
  calc_dt[, is_dirty_past_int := as.integer(is_dirty_past)]
  calc_dt[, sum_dirty_pre := frollsum(is_dirty_past_int, n = k_pre, align = "right", fill = NA), by = id]

  calc_dt[, is_valid_treated := (event == 1) &
            (has_history == TRUE) &
            (has_future == TRUE) &
            (shift(sum_dirty_pre, 1) == 0)]

  # Valid Control
  calc_dt[, is_dirty_any := is_dirty_past | is_dirty_future]
  calc_dt[, is_dirty_any_int := as.integer(is_dirty_any)]
  win_len <- k_pre + 1L + k_post
  calc_dt[, sum_dirty_full := frollsum(is_dirty_any_int, n = win_len, align = "center", fill = NA), by = id]

  calc_dt[, is_valid_control := (has_history == TRUE) &
            (has_future == TRUE) &
            (sum_dirty_full == 0)]

  # 4. Filter Candidates
  treats <- calc_dt[is_valid_treated == TRUE, .(id, cohort_year = period, time_idx, treat = 1, first_treat_idx)]

  if (treated_type == "first") {
    treats <- treats[time_idx == first_treat_idx]
    treats <- unique(treats, by = "id")
  }
  treats[, c("first_treat_idx") := NULL]

  ctrls_raw <- calc_dt[is_valid_control == TRUE, .(id, cohort_year = period, time_idx, treat = 0, is_ever_treated, first_treat_idx)]

  if (control_type == "nevert") {
    ctrls <- ctrls_raw[is_ever_treated == FALSE]
  } else if (control_type == "nt") {
    ctrls <- ctrls_raw[is_ever_treated == TRUE]
  } else if (control_type == "nyt") {
    ctrls <- ctrls_raw[(is_ever_treated == TRUE) & (time_idx < first_treat_idx)]
  } else {
    ctrls <- ctrls_raw
  }
  ctrls[, c("is_ever_treated", "first_treat_idx") := NULL]

  candidates <- rbind(treats, ctrls)

  # --- FIX: EXACT MATCH PERIOD LOOKUP ---
  # We calculate the Target Index (e.g. Index 10 + (-1) = Index 9)
  candidates[, match_time_idx := time_idx + match_by]

  # We map Index 9 back to the REAL Date (e.g., 2015-04-01) using the original data
  time_map <- calc_dt[, .(id, time_idx, real_match_period = period)]

  # Inner Join serves two purposes:
  # 1. Retrieves correct Date
  # 2. Filters out cases where lag is out of bounds
  candidates <- merge(candidates, time_map,
                      by.x = c("id", "match_time_idx"),
                      by.y = c("id", "time_idx"),
                      all.x = FALSE)

  candidates[, match_period := real_match_period]
  # Cleanup
  candidates[, c("match_time_idx", "real_match_period", "time_idx") := NULL]

  # Restore user's ID column name
  setnames(candidates, "id", id_col)

  return(candidates)
}

# ==============================================================================
# 2. MATCH: Cross-Sectional Matching
# ==============================================================================
#' Cross-Sectional Matching for Stacked DiD
#'
#' @description
#' Performs exact cohort matching + 1-to-K matching using MatchIt.
#'
#' @export
match_candidates <- function(dt, stack_candidates, id_col, time_col, match_vars, ...) {

  if (!inherits(dt, "data.table")) setDT(dt)

  # 1. Prepare Covariate Data
  cols_to_fetch <- c(id_col, time_col, match_vars)
  cov_data <- dt[, .SD, .SDcols = cols_to_fetch]

  # Rename for internal safety
  setnames(cov_data, c(id_col, time_col), c("id", "match_period"))

  # 2. Prepare Candidate Data
  match_input_base <- copy(stack_candidates)
  setnames(match_input_base, id_col, "id")

  # 3. Merge
  match_input <- merge(match_input_base, cov_data,
                       by = c("id", "match_period"),
                       all.x = FALSE)

  match_input <- na.omit(match_input, cols = match_vars)

  # 4. MatchIt
  fm <- as.formula(paste("treat ~", paste(match_vars, collapse = " + ")))
  args <- list(...)
  args$formula <- fm
  args$data    <- match_input

  if ("exact" %in% names(args)) {
    cur <- args$exact
    if (inherits(cur, "formula")) {
      args$exact <- as.formula(paste("~", as.character(cur)[2], "+ cohort_year"))
    } else {
      args$exact <- as.formula(paste("~", paste(c(cur, "cohort_year"), collapse=" + ")))
    }
  } else {
    args$exact <- ~ cohort_year
  }

  m_out <- do.call(MatchIt::matchit, args)

  # 5. Output
  matched_data <- as.data.table(MatchIt::get_matches(m_out, id = "matchit_row_id"))
  matched_data[, stack_id := paste(id, cohort_year, sep = "_")]

  treated_map <- matched_data[treat == 1, .(subclass, cohort_id = stack_id)]
  matched_data <- merge(matched_data, treated_map, by = "subclass", all.x = TRUE)

  # Restore User ID Name
  setnames(matched_data, "id", id_col)

  cols_to_keep <- c("stack_id", "cohort_id", id_col, "cohort_year",
                    "treat", "match_period", "weights", match_vars)
  cols_to_keep <- intersect(cols_to_keep, names(matched_data))

  final_output <- matched_data[, .SD, .SDcols = cols_to_keep]

  return(list(matched_cross_section = final_output, match_object = m_out))
}

# ==============================================================================
# 3. BUILD: Expand Stacks
# ==============================================================================
#' Expand Matched Cross-Section to Full Stacks
#'
#' @export
expand_matched_stacks <- function(dt, matched_data, id_col, time_col, k_pre, k_post, match_by = -1) {

  if (!inherits(dt, "data.table")) setDT(dt)
  if (!inherits(matched_data, "data.table")) setDT(matched_data)

  # 1. Build Time Index Map from ORIGINAL data
  # This tells us: For ID "A", Index 5 = "2015-01-01", Index 6 = "2015-04-01"
  # This abstracts away the "Days vs Quarters" problem.
  time_map <- dt[, .SD, .SDcols = c(id_col, time_col)]
  setnames(time_map, c(id_col, time_col), c("id", "period"))
  setorderv(time_map, c("id", "period"))
  time_map[, time_idx := seq_len(.N), by = id]

  # 2. Get Cohort Index
  # We need to know: What Index corresponds to the 'cohort_year' (Event Time)?
  matched_base <- copy(matched_data)
  setnames(matched_base, id_col, "id")

  # Join map to get the Index of the Cohort Date
  matched_base <- merge(matched_base, time_map[, .(id, period, cohort_idx = time_idx)],
                        by.x = c("id", "cohort_year"),
                        by.y = c("id", "period"),
                        all.x = TRUE, sort = FALSE)

  # 3. Expand Rel Times
  rel_times <- seq(-k_pre, k_post)

  meta_cols <- c("stack_id", "cohort_id", "id", "cohort_idx", "treat", "weights", "cohort_year")
  meta_cols <- intersect(meta_cols, names(matched_base))

  stack_defs <- matched_base[, ..meta_cols]

  n_stacks <- nrow(stack_defs)
  expanded_stacks <- stack_defs[rep(1:n_stacks, each = length(rel_times))]
  expanded_stacks[, rel_time := rep(rel_times, times = n_stacks)]

  # 4. Calculate Target Index
  # Example: Cohort is Index 10. Rel_time -2. Target = Index 8.
  expanded_stacks[, target_idx := cohort_idx + rel_time]

  # 5. Retrieve Actual Date from Map using Target Index
  # This effectively does: "What is the date 2 rows before the event?"
  # It works perfectly for Quarters, Days, Years, etc.
  expanded_stacks <- merge(expanded_stacks, time_map[, .(id, target_idx = time_idx, period)],
                           by = c("id", "target_idx"),
                           all.x = TRUE, sort = FALSE)

  # 6. Merge Variables from DT using the retrieved Date
  final_df <- merge(expanded_stacks, dt,
                    by.x = c("id", "period"),
                    by.y = c(id_col, time_col),
                    all.x = TRUE, sort = FALSE)

  final_df[, post := as.integer(rel_time >= 0)]

  # Cleanup & Rename
  final_df[, c("cohort_idx", "target_idx") := NULL]
  setnames(final_df, "period", time_col)
  setnames(final_df, "id", id_col)

  first_cols <- c("stack_id", "cohort_id", id_col, time_col, "rel_time", "treat", "post", "weights")
  rest_cols <- setdiff(names(final_df), first_cols)
  setcolorder(final_df, c(first_cols, rest_cols))

  return(final_df)
}


# ==============================================================================
# 4. AUDIT: Quality Assessment
# ==============================================================================
#' Assess Matching Quality (Academic Standard)
#'
#' @export
assess_match_quality <- function(match_res, stack_candidates, dt, match_vars, id_col, time_col) {

  if (!inherits(dt, "data.table")) setDT(dt)

  m_out <- match_res$match_object
  matched_df <- match_res$matched_cross_section

  covariates <- intersect(match_vars, names(matched_df))
  if (length(covariates) == 0) stop("None of the 'match_vars' found in matched results.")

  # Unmatched Pool
  cols_fetch <- c(id_col, time_col, covariates)
  dt_subset <- dt[, ..cols_fetch]

  sc_copy <- copy(stack_candidates)
  setnames(sc_copy, id_col, "id")

  dt_sub_copy <- copy(dt_subset)
  setnames(dt_sub_copy, c(id_col, time_col), c("id", "match_period"))

  unmatched_pool <- merge(sc_copy, dt_sub_copy,
                          by = c("id", "match_period"),
                          all.x = TRUE)
  unmatched_pool <- na.omit(unmatched_pool, cols = covariates)

  # Labels
  get_lab_safe <- function(v, df) {
    if (!v %in% names(df)) return(v)
    lbl <- attr(df[[v]], "label")
    if (is.null(lbl) || length(lbl) == 0) return(as.character(v))
    return(as.character(lbl))
  }
  var_labs <- vapply(covariates, get_lab_safe, df = matched_df, FUN.VALUE = character(1))
  names(var_labs) <- covariates

  # Plot Data
  dt_matched <- melt(matched_df[, c("treat", "weights", covariates), with=FALSE],
                     measure.vars = covariates, variable.name = "term")
  dt_matched[, Stage := "Matched"]

  dt_unmatched <- melt(unmatched_pool[, c("treat", covariates), with=FALSE],
                       measure.vars = covariates, variable.name = "term")
  dt_unmatched[, weights := 1]
  dt_unmatched[, Stage := "Unmatched"]

  dt_plot <- rbind(dt_matched, dt_unmatched, fill = TRUE)
  dt_plot[, Variable_Label := var_labs[as.character(term)]]
  dt_plot[, Group := factor(treat, levels = c(0, 1), labels = c("Control", "Treated"))]

  # Density Plot
  p_density <- ggplot() +
    geom_density(data = dt_plot[Stage == "Unmatched"],
                 aes(x = value, color = Group, linetype = "Unmatched"),
                 size = 0.6, alpha = 0.6) +
    geom_density(data = dt_plot[Stage == "Matched"],
                 aes(x = value, color = Group, weight = weights, linetype = "Matched"),
                 size = 1) +
    facet_wrap(~ Variable_Label, scales = "free", ncol = 3) +
    scale_linetype_manual(name = "Stage", values = c("Unmatched" = "dashed", "Matched" = "solid")) +
    scale_color_manual(values = c("Control" = "#377EB8", "Treated" = "#E41A1C")) +
    labs(title = "Covariate Balance: Full Pool (Dashed) vs Matched Set (Solid)",
         subtitle = "Solid lines should overlap. Dashed lines show original bias.",
         x = "Covariate Value", y = "Density") +
    theme_minimal() +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))

  # Table
  s_out <- summary(m_out, standardize = TRUE)
  stats <- as.data.table(s_out$sum.matched, keep.rownames = "term")
  stats <- stats[term %in% covariates]

  clean_stats <- stats[, .(Variable = var_labs[term],
                           Mean_T = `Means Treated`, Mean_C = `Means Control`,
                           SMD = `Std. Mean Diff.`, VR = `Var. Ratio`)]
  clean_stats[, Quality_SMD := ifelse(abs(SMD) < 0.1, "Excellent", ifelse(abs(SMD) < 0.25, "Acceptable", "Poor"))]
  clean_stats[, Quality_VR  := ifelse(VR >= 0.8 & VR <= 1.25, "Balanced", "Distorted")]
  clean_stats[is.na(VR), Quality_VR := "N/A"]

  return(list(balance_table = clean_stats, density_panel = p_density))
}

# ==============================================================================
# 5. VISUALIZE: Event Dynamics
# ==============================================================================
#' Plot Event Dynamics (Lasagna Plot with Overlay)
#'
#' @export
plot_event_dynamics <- function(dt, dt_matched = NULL, id_col, time_col, event_col, n_bins = 4L, sort_by_intensity = FALSE) {

  if (!inherits(dt, "data.table")) setDT(dt)

  d <- copy(dt[, .SD, .SDcols = c(id_col, time_col, event_col)])
  setnames(d, c(id_col, time_col, event_col), c("id", "time", "event"))
  setkey(d, NULL)

  t_levels <- sort(unique(d$time))

  id_stats <- d[, {
    tr_times <- time[event == 1L]
    tr_idx <- match(tr_times, t_levels)
    list(n_event = length(tr_idx), tr_idx = list(tr_idx))
  }, by = id]

  id_stats[, is_treated := n_event > 0]

  if (isTRUE(sort_by_intensity)) {
    treated_mask <- id_stats$is_treated == TRUE
    if (sum(treated_mask) > 0) {
      qs <- unique(quantile(id_stats$n_event[treated_mask], probs = seq(0, 1, length.out = n_bins + 1L), type = 1))
      if (length(qs) < 2) {
        id_stats[treated_mask, bin := 1L]
      } else {
        id_stats[treated_mask, bin := cut(n_event, breaks = qs, include.lowest = TRUE, labels = FALSE)]
      }
    }
    id_stats[is.na(bin), bin := 0L]
  } else {
    id_stats[, bin := ifelse(is_treated, 1L, 0L)]
  }

  K <- max(id_stats$n_event, na.rm = TRUE)
  if (!is.finite(K) || K < 1L) K <- 1L

  pad_func <- function(v) {
    out <- rep(-1L, K)
    if (length(v) > 0) out[1:length(v)] <- v
    as.integer(out)
  }

  padded_list <- lapply(id_stats$tr_idx, pad_func)
  keys <- as.data.table(do.call(rbind, padded_list))
  time_key_cols <- paste0("t", seq_len(K))
  setnames(keys, time_key_cols)
  id_stats <- cbind(id_stats, keys)

  sort_order <- c("bin", time_key_cols, "id")
  order_directions <- c(1, rep(-1, length(time_key_cols)), 1)
  setorderv(id_stats, cols = sort_order, order = order_directions)

  id_levels_sorted <- id_stats$id

  p <- ggplot() +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(title = "Event Dynamics & Matched Stacks",
         subtitle = if(!is.null(dt_matched)) "Green: Events | Red: Matched Treated (Solid) | Blue: Matched Control (Freq)" else "Green: Events",
         x = time_col, y = id_col) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  plot_data_base <- d[event == 1L]
  if (nrow(plot_data_base) > 0) {
    plot_data_base[, id_factor := factor(id, levels = id_levels_sorted)]
    plot_data_base[, time_factor := factor(time, levels = t_levels)]
    p <- p + geom_raster(data = plot_data_base, aes(x = time_factor, y = id_factor), fill = "green4")
  }

  if (!is.null(dt_matched)) {
    if (!inherits(dt_matched, "data.table")) setDT(dt_matched)
    use_id <- if (id_col %in% names(dt_matched)) id_col else "id"
    use_time <- if (time_col %in% names(dt_matched)) time_col else "period"

    mat_trt <- dt_matched[treat == 1, .N, by = c(use_id, use_time)]
    setnames(mat_trt, c("id", "time", "N"))
    mat_trt <- mat_trt[id %in% id_levels_sorted & time %in% t_levels]
    if (nrow(mat_trt) > 0) {
      mat_trt[, id_factor := factor(id, levels = id_levels_sorted)]
      mat_trt[, time_factor := factor(time, levels = t_levels)]
      p <- p + geom_raster(data = mat_trt, aes(x = time_factor, y = id_factor), fill = "firebrick", alpha = 1)
    }

    mat_ctrl <- dt_matched[treat == 0, .N, by = c(use_id, use_time)]
    setnames(mat_ctrl, c("id", "time", "N"))
    mat_ctrl <- mat_ctrl[id %in% id_levels_sorted & time %in% t_levels]
    if (nrow(mat_ctrl) > 0) {
      mat_ctrl[, id_factor := factor(id, levels = id_levels_sorted)]
      mat_ctrl[, time_factor := factor(time, levels = t_levels)]
      p <- p + geom_raster(data = mat_ctrl, aes(x = time_factor, y = id_factor, alpha = N), fill = "steelblue")
      p <- p + scale_alpha_continuous(range = c(0.3, 1), name = "Ctrl Freq.")
    }
  }

  return(list(plot = p, order_table = id_stats))
}

# ==============================================================================
# 6. MASTER WRAPPER: flexmatch()
# ==============================================================================
#' Flexmatch: The All-in-One Stacked DiD Wrapper
#'
#' @description
#' Orchestrates the entire Stacked DiD pipeline:
#' Scans candidates, Matches them, Expands to stacks, Assesses Quality, and Visualizes Events.
#'
#' @param dt The original panel data (data.table).
#' @param id_col Character. Name of ID column.
#' @param time_col Character. Name of time column (integer/numeric).
#' @param event_col Character. Name of event dummy (1 = start of treatment).
#' @param k_pre Integer. Required pre-periods.
#' @param k_post Integer. Required post-periods.
#' @param k_freeze Integer. Freeze/Cool-down periods after k_post.
#' @param match_vars Character vector. Covariates to match on.
#' @param match_by Integer. Relative period for matching (default -1).
#' @param assess_quality Logical. If TRUE, runs balance checks. Default TRUE.
#' @param plot_events Logical. If TRUE, runs the event dynamics plot. Default TRUE.
#' @param verbose Logical. If TRUE, prints progress. Default TRUE.
#' @param ... Additional arguments passed directly to `MatchIt` (e.g., `ratio`, `method`, `replace`, `caliper`).
#'
#' @return A list containing:
#'   - `final_expanded_data`: The analysis-ready stacked dataset.
#'   - `matched_cross_section`: The matched pairs at time t.
#'   - `stack_candidates`: The full pool of eligible units.
#'   - `quality_stats`: List with Balance Table and Density Plots.
#'   - `event_plot`: The Lasagna plot object.
#'   - `match_object`: The raw MatchIt model.
#' @export
flexmatch <- function(dt,
                      id_col,
                      time_col,
                      event_col,
                      k_pre,
                      k_post,
                      k_freeze,
                      match_vars,
                      match_by = -1,
                      treated_type = "first",
                      control_type = "all",
                      assess_quality = TRUE,
                      plot_events = TRUE,
                      verbose = TRUE,
                      ...) {

  if (verbose) message("\n=== Starting Flexmatch Pipeline ===")
  if (!inherits(dt, "data.table")) setDT(dt)

  if (verbose) message(sprintf("1. Scanning Candidates (Type: Treated='%s', Control='%s')...", treated_type, control_type))

  cands <- label_candidates(
    dt = dt, id_col = id_col, time_col = time_col, event_col = event_col,
    k_pre = k_pre, k_post = k_post, k_freeze = k_freeze, match_by = match_by,
    treated_type = treated_type, control_type = control_type
  )

  n_treat <- nrow(cands[treat == 1])
  if (n_treat == 0) stop("No valid treated candidates found.")
  if (verbose) message(sprintf("   Candidates found: %d Treated, %d Control.", n_treat, nrow(cands[treat == 0])))

  if (verbose) message(sprintf("2. Matching on [%s]...", paste(match_vars, collapse = ", ")))
  match_res <- match_candidates(dt, cands, id_col, time_col, match_vars, ...)

  # FIX: Count matches using get(id_col)
  n_matched <- uniqueN(match_res$matched_cross_section[treat==1, get(id_col)])
  if (verbose) message(sprintf("   Matched %d / %d Treated units (%.1f%%).", n_matched, n_treat, 100 * n_matched / n_treat))

  if (verbose) message("3. Expanding Stacks...")
  final_data <- expand_matched_stacks(dt, match_res$matched_cross_section, id_col, time_col, k_pre, k_post, match_by)

  quality_stats <- NULL
  if (assess_quality) {
    if (verbose) message("4. Assessing Quality...")
    tryCatch({
      quality_stats <- assess_match_quality(match_res, cands, dt, match_vars, id_col, time_col)
    }, error = function(e) warning(paste("Quality Check Failed:", e$message)))
  }

  event_plot <- NULL
  if (plot_events) {
    if (verbose) message("5. Generating Event Dynamics Plot...")
    tryCatch({
      ep <- plot_event_dynamics(
        dt = dt,
        dt_matched = final_data,
        id_col = id_col,
        time_col = time_col,
        event_col = event_col,
        sort_by_intensity = FALSE
      )
      event_plot <- ep$plot
    }, error = function(e) warning(paste("Event Plot Failed:", e$message)))
  }

  if (verbose) message("=== Pipeline Complete ===\n")

  return(list(
    final_expanded_data   = final_data,
    matched_cross_section = match_res$matched_cross_section,
    stack_candidates      = cands,
    quality_stats         = quality_stats,
    event_plot            = event_plot,
    match_object          = match_res$match_object
  ))
}
