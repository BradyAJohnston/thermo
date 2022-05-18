library(tidyverse)

# fl <-
#   "~/Dropbox/BondLab/Data/MST data/Brady/13rna_com_exp/20200201/20200201_CompBind_assay.moc"
#
fl <- "~/Dropbox/BondLab/Data/MST data/Brady/Paper Bindg/PPRatph.moc"

tables_to_extract <- c("mCapScan",
                       "mMst",
                       "Annotation",
                       "tCapillary",
                       "tContainer")

test <- function(x, fl) {
  thermo::extract_table(fl, x) |>
    janitor::clean_names()
}

tables <- tables_to_extract |>
  purrr::set_names() |>
  purrr::map2(.y = fl, purrr::possibly(test, NA))

mst_values <- tables$mMst |>
  dplyr::group_by(id) |>
  dplyr::mutate(
    trace_mst = purrr::map(mst_trace, purrr::possibly(blob_extract_1d, NA)),
    step_times = list(
      c(
        nominal_duration_of_phase1,
        nominal_duration_of_phase2,
        nominal_duration_of_phase3
      )
    ),
    trace_mst = purrr::map2(trace_mst, step_times, purrr::possibly(function(x, y) {
      df <- tibble(fnorm = x / x[1],
                   time = seq_along(x) - 1)
      df |>
        dplyr::mutate(time = (time) / max(time) * sum(y) - 0.5)
    }, NA))
  ) |>
  dplyr::select(id, container, trace_mst, step_times)

capillary_values <- tables$tCapillary |>
  dplyr::select(id,
                annotations,
                index_on_parent_container,
                container_type,
                parent_container) |>
  dplyr::rename(id_capillary = id,
                cap = index_on_parent_container,
                cap_type = container_type) |>
  dplyr::mutate(cap = cap + 1) |>
  tidyr::separate_rows(annotations, sep = ";")

annotation_values <- tables$Annotation |>
  dplyr::select(id,
                caption,
                annotation_role,
                annotation_type,
                numeric_value,
                text_value)

tables$Annotation |>
  )
pull(tag) |>
  unique()

ext_values <- function(x) {
  start <- x |>
    dplyr::filter(time < 5 & time > 3) |>
    dplyr::pull(fnorm) |>
    mean()

  end <- x |>
    dplyr::filter(time > 20 & time < 22) |>
    dplyr::pull(fnorm) |>
    mean()

  diff <- end / start

  diff
}


df <- mst_values |>
  dplyr::left_join(capillary_values, by = c('container' = 'id_capillary')) |>
  dplyr::left_join(annotation_values, by = c('annotations' = 'id'))

df |>
  filter(
    annotation_type != "text"
  ) |>
  select(id, cap, caption, annotation_role, numeric_value) |>

  pivot_wider(
    names_from = annotation_role,
    values_from = c(caption, numeric_value)
  ) |>
  rename(
    ligand = caption_ligand,
    target = caption_target,
    conc_ligand = numeric_value_ligand,
    conc_target = numeric_value_target
  )
pivot_wider(
  names_from = annotation_role,
  values_from = caption,
) |>
  mutate(
    type = if_else(
      is.na(ligand),
      "target",
      "ligand"
    )
  ) |>
  pivot_wider(
    names_from = type,
    values_from = numeric_value,
    names_prefix = "conc_",
    names_repair = "unique"
  ) |>
  fill(ligand, target, conc_ligand, conc_target, .direction = c("downup")) |>
  unique()

pivot_wider(id_cols = c(id, container, trace_mst, cap, caption),
            names_from = c(annotation_role),
            values_from = numeric_value)
tidyr::pivot

values <- df |>
  dplyr::filter(!is.na(map(trace_mst, `[[`, 1))) |>
  dplyr::filter(annotation_type != "text") |>
  pivot_wider(
    names_from = annotation_role,
    values_from = c(caption, numeric_value, text_value, annotations)
  ) |>
  rename(
    ligand = caption_ligand,
    target = caption_target,
    conc_ligand = numeric_value_ligand,
    conc_target = numeric_value_target
  ) |>
  dplyr::select(
    id,
    container,
    trace_mst,
    annotation_type,

    parent_container,
    target,
    ligand,
    conc_ligand,
    conc_target
  ) |>
  # dplyr::filter(annotation_role == "ligand") |>
  group_by(id) |>
  dplyr::mutate(
    diff = purrr::map_dbl(trace_mst, purrr::possibly(ext_values, NA))
  )

values |>
  mutate(group = interaction(target)) |>
  group_by(group) |>
  biochemr::bio_binding(conc_ligand, diff) |>
  biochemr::bio_plot() +
  ggplot2::scale_x_log10() +
  labs(x = "[RNA] nm",
       y = "Fnorm")
