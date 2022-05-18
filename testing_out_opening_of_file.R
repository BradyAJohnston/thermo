library(RSQLite)
library(tidyverse)

fl <- "~/Dropbox/BondLab/Data/MST data/Brady/13rna_com_exp/20200201/20200201_CompBind_assay.moc"

con <- dbConnect(RSQLite::SQLite(), fl)


df <- dbReadTable(con, name = "mCapScan") |>
  as_tibble() |>
  # select(ID, CapScanTrace) |>
  filter(!is.na(CapScanTrace))




blob_extract_1d <- function(blob, type = "double") {
  numbers <- packBits(rawToBits(unlist(blob)), type = type)
  df <- data.frame(
    x = numbers
  )

  return(df)
}

blob_extract_2d <- function(blob, type = "double") {
  numbers <- packBits(rawToBits(unlist(blob)), type = type)

  df <- data.frame(
    x = numbers[seq(1, length(numbers), by = 2)],
    y = numbers[seq(2, length(numbers), by = 2)]
  )

  return(df)
}





df |>
  mutate(
    trace = map(CapScanTrace, blob_extract_2d)
  ) |>
  unnest(trace) |>
  ggplot(aes(x, y)) +
  geom_line(aes(group = CenterPosition))

dbReadTable(con, name = "mMST") |>
  as_tibble() |>
  select(ID, MstTrace) |>
  filter(!is.na(MstTrace)) |>
  mutate(
    trace = map(MstTrace, blob_extract_1d)
  ) |>
  select(ID, trace) |>
  unnest(trace) |>
  group_by(ID) |>
  mutate(int = row_number()) |>
  mutate(trace = x / x[1]) |>
  ggplot(aes(int, trace)) +
  geom_line(aes(group = ID))



fl <- "~/Dropbox/BondLab/Data/MST data/Brady/13rna_com_exp/20200201/20200201_CompBind_assay.moc"

con <- dbConnect(RSQLite::SQLite(), fl)

df_mst <- dbReadTable(con, name = "mMST") |>
  as_tibble() |>
  # select(ID, MstTrace) |>
  filter(!is.na(MstTrace)) |>
  mutate(
    trace = map(MstTrace, ~packBits(rawToBits(unlist(.x)), type = "double"))
  )

df_mst <- df_mst |>
  group_by(ID) |>
  mutate(
    time_sum = list(c(NominalDurationOfPhase1, NominalDurationOfPhase2, NominalDurationOfPhase3))
  ) |>
  # select(ID, time_sum) |
  ungroup() |>

  mutate(trace = map2(trace, time_sum, function(x, y) {
    df <- tibble(
      fnorm = x / x[1],
      time = seq_along(x) - 1
    )
    df |>
      mutate(time = (time) / max(time) * sum(y) - 0.5)
  }))

df_mst |>
  select(ID, trace, NominalDurationOfPhase1) |>
  unnest(trace) |>
  ggplot(aes(time, fnorm)) +
  geom_vline(aes(xintercept = NominalDurationOfPhase1)) +
  geom_line(aes(group = ID), alpha = 0.05) +
  # coord_cartesian(xlim = c(4, 6)) +
  theme_light()


df_analysis <- dbReadTable(con, name = "AffinityAnalysisMstTrace") |>
  as_tibble() |>
  select(ID, RawMst, ParentRun)
df_analysis |>
  left_join(df_mst |> select(ID, trace), by = c("RawMst" = "ID")) |>
  unnest(trace) |>
  ggplot(aes(time, fnorm, group = RawMst)) +
  geom_line() +
  facet_wrap(~ParentRun) +
  theme_light()


df_anno <- dbReadTable(con, name = "Annotation") |>
  as_tibble()
# select(ID, AnnotationRole, AnnotationType, Caption, NumericValue, TextValue) |>
# left_join(df_analysis, by = c("ID" = "ID"))



df_values <- dbReadTable(con, name = "tCapillary") |>
  as_tibble() |>
  select(ID, Annotations, Caption, ContainerType, ParentContainer) |>
  separate_rows(
    Annotations, sep = ";"
  ) |>
  left_join(df_anno, by = c("Annotations" = "ID"))


df_values |>
  select(
    ID,
    Annotations,
    ParentContainer,
    Caption.x,
    AnnotationRole,
    AnnotationType,
    Caption.y,
    NumericValue,
    Tag,
    TextValue
  ) |>
  left_join(df_mst, by = c("ParentContainer" = "Container"))

dbReadTable(con, name = "tContainer") |>
  as_tibble() |>

  df_mst |>
  select(ID, Container, ExcitationPower, MstPower, trace, time_sum) |>
  left_join(df_values, by = )
df_mst |> select(matches("Table"))
dbReadTable(con, "ExpertModeCapillarySettings") |>
  as_tibble()

dbReadTable(con, name = "tCapillary") |>
  as_tibble() |>
  select(ID, Annotations, ParentContainer) |>
  left_join(dbReadTable(con, name = "tContainer") |>
              as_tibble(),
            by = c("ParentContainer" = "ID")) |>
  left_join(df_values, by = c("ParentContainer" = "ParentContainer")) |>
  select(ID.x, ID.y, ParentContainer, NumericValue) |>
  left_join(df_mst, by = c("ParentContainer" = "Container")) |>
  drop_na()

df_mst |>
  select(ID, Container, trace, time_sum) |>
  left_join(df_values, by = c("Container" = "ID")) |>
  unnest(trace) |>
  filter(AnnotationRole == "ligand") |>
  ggplot(aes(time, fnorm, gorup = interaction(ParentContainer, Caption.x))) +
  geom_line(alpha = 0.6, mapping = aes(colour = log10(NumericValue))) +
  facet_wrap(~Caption.y~ParentContainer) +
  scale_colour_viridis_c() +
  theme_light()
