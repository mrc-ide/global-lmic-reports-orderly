date_0 <- as.Date(date, "%Y-%m-%d")

#get open data from next strain
download.file("https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz",
              "next_strain.gz")
#read in chunkwise and simplify
next_strain <- read_tsv_chunked(gzfile("next_strain.gz"), chunk_size = 5000, callback = DataFrameCallback$new(function(x, pos){
  group_by(x, date, country, Nextstrain_clade) %>%
    summarise(
      count = n(),
      .groups = "keep"
    )
})) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  )
unlink("next_strain.gz")

#convert clades to variants
clade_to_variant <- function(clades){
  variants <- case_when(
    clades == "recombinant" | is.na(clades) ~ as.character(NA),
    clades == "22A (Omicron)" | clades == "22B (Omicron)" ~ "Omicron Sub-Variant",
    TRUE ~ stringr::str_sub(clades, 5, -1) %>%
      stringr::str_remove_all("[\\(\\)V\\d, ]")
  )
  #restrict to the who variants
  who_variants <- c("Alpha", "Beta", "Gamma", "Delta", "Omicron", "Omicron Sub-Variant")

  variants <- case_when(
    is.na(variants) ~ as.character(NA),
    !variants %in% who_variants ~ "Wild",
    TRUE ~ variants
  )
  ordered(variants, levels = c("Wild", who_variants))
}
variants_df <- next_strain %>%
  mutate(variant = clade_to_variant(Nextstrain_clade)) %>%
  filter(!is.na(variant)) %>%
  group_by(date, country, variant) %>%
  summarise(count = sum(count),
            .groups = "drop")

#format data
iso3cs <- unique(squire::population$iso3c)
variants_df <- variants_df %>%
  filter(nchar(date) == 10 & date != "XXXX-XX-XX") %>%
  transmute(iso3c = countrycode::countrycode(country, "country.name", "iso3c",
                                          custom_match = c(Kosovo = "XKX")),
         date = as_date(date), variant = variant, count = count) %>%
  #ensure we have an entry for all countries dates and variants
  complete(
    iso3c = iso3cs, variant = unique(variant), date = seq(min(date), max(date), by = 1), fill = list(count = 0)
  ) %>%
  group_by(iso3c) %>%
  arrange(date, desc(count)) %>%
  filter(cumsum(count) > 0) %>%
  arrange(date, variant)

#function to find the start/end times of a variant
linearly_interpolate <- function(vector){
  if_else(
    is.na({vector}),
    approx(1:length({vector}), {vector}, 1:length({vector}), yright = NA)$y,
    {vector}
  )
}
find_timings <- function(df){
  variants <- setdiff(unique(df$variant), "Wild")
  #convert counts into weekly percentages
  df <- df %>%
    mutate(week = floor_date(date, "week")) %>%
    group_by(week, variant) %>%
    summarise(count = sum(count), .groups = "drop_last") %>%
    filter(sum(count) > 20) %>%
    mutate(p_var = count/sum(count)) %>%
    ungroup() %>%
    complete(week = seq(min(week), max(week), 1), variant = unique(variant),
             fill = list(p_var = NA)) %>%
    select(!count) %>%
    #linearlly interpolate the percentages
    group_by(variant) %>%
    arrange(week) %>%
    mutate(p_var = linearly_interpolate(p_var))
  map_dfr(variants, function(this_variant){
    temp <- filter(df, variant == this_variant)
    #find the mid point, (i.e first date where p_var > 0.5)
    mid_date <- suppressWarnings(min(temp$week[temp$p_var > 0.5]))
    if(is.infinite(mid_date)){
      mid_date <- temp$week[which.max(temp$p_var)]
    }
    #start_date is last then where p_var < 10 before the mid_date
    start_date <- suppressWarnings(temp %>%
      filter(week < mid_date & p_var < 0.1) %>%
      pull(week) %>%
      max)
    if(is.infinite(start_date)){
      start_date <- mid_date - 30
    }
    #end_date is the first time p_var >80% after mid_date
    end_date <- temp %>%
      filter(week > mid_date & p_var > 0.8)
    if(nrow(end_date) > 0){
      end_date <- end_date %>%
        pull(week) %>%
        max
    } else {
      end_date <- mid_date + as.numeric(mid_date - start_date)
    }
    tibble(start_date = start_date, mid_date = mid_date, end_date = end_date)
  }) %>%
    mutate(variant = variants)
}
can_fit_to <- function(df){
  df %>%
    mutate(week = floor_date(date, "week")) %>%
    group_by(week, .add = TRUE) %>%
    summarise(count = sum(count), .groups = "drop_last") %>%
    filter(count > 20) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 20) %>%
    select(!count)
}

#Determine which countries have enough data to fit to (more than 20 sequences
#on more than 20 weeks)
variants_df <- variants_df %>%
  mutate(
    un_region = countrycode::countrycode(iso3c, "iso3c", "un.region.name",
                                         custom_match = c(TWN = "Asia", XKX = "Europe")),
    un_sub_region = countrycode::countrycode(iso3c, "iso3c", "un.regionsub.name",
                                             custom_match = c(TWN = "Eastern Asia", XKX = "Eastern Europe"))
  )
fitting_info <- tibble(
  iso3c = iso3cs,
  fitting_type = "World"
) %>%
  left_join(
    variants_df %>%
      right_join(
        variants_df %>%
          group_by(un_region) %>%
          can_fit_to(),
        by = "un_region"
      ) %>%
      select(iso3c) %>%
      unique() %>%
      mutate(region_level = TRUE),
    by = "iso3c"
  ) %>%
  left_join(
    variants_df %>%
      right_join(
        variants_df %>%
          group_by(un_sub_region) %>%
          can_fit_to(),
        by = "un_sub_region"
      ) %>%
      select(iso3c) %>%
      unique() %>%
      mutate(sub_region_level = TRUE),
    by = "iso3c"
  ) %>%
  left_join(
    variants_df %>%
      group_by(iso3c) %>%
      can_fit_to() %>%
      mutate(country_level = TRUE),
    by = "iso3c"
  ) %>%
  transmute(
    iso3c = iso3c,
    fitting_type = case_when(
      country_level ~ "Country",
      sub_region_level ~ "Sub Region",
      region_level ~ "Region",
      TRUE ~ "World"
    )
  ) %>%
  left_join(
    variants_df %>%
      select(iso3c, un_region, un_sub_region) %>%
      unique(),
    by = "iso3c"
  )

#now figure out which ones we need to fit

if("World" %in% fitting_info$fitting_type){
  world_fit <- variants_df %>%
    group_by(variant, date) %>%
    summarise(
      count = sum(count),
      .groups = "drop"
    ) %>%
    find_timings
}

if("Region" %in% fitting_info$fitting_type){
  regions_to_fit <- fitting_info %>%
    filter(fitting_type == "Region") %>%
    pull(un_region) %>%
    unique()
  names(regions_to_fit) <- regions_to_fit
  region_fits <- map(
    regions_to_fit,
    ~variants_df %>%
      filter(un_region == .x) %>%
      group_by(variant, date) %>%
      summarise(
        count = sum(count),
        .groups = "drop"
      ) %>%
      find_timings
  )
  rm(regions_to_fit)
}

if("Sub Region" %in% fitting_info$fitting_type){
  regions_to_fit <- fitting_info %>%
    filter(fitting_type == "Sub Region") %>%
    pull(un_sub_region) %>%
    unique()
  names(regions_to_fit) <- regions_to_fit
  sub_region_fits <- map(
    regions_to_fit,
    ~variants_df %>%
      filter(un_sub_region == .x) %>%
      group_by(variant, date) %>%
      summarise(
        count = sum(count),
        .groups = "drop"
      ) %>%
      find_timings
  )
  rm(regions_to_fit)
}

if("Country" %in% fitting_info$fitting_type){
  regions_to_fit <- fitting_info %>%
    filter(fitting_type == "Country") %>%
    pull(iso3c) %>%
    unique()
  names(regions_to_fit) <- regions_to_fit
  country_fits <- map(
    regions_to_fit,
    ~variants_df %>%
      filter(iso3c == .x) %>%
      group_by(variant, date) %>%
      summarise(
        count = sum(count),
        .groups = "drop"
      ) %>%
      find_timings
  )
  rm(regions_to_fit)
}

#now compile these into a list
variant_timings <- map(purrr::transpose(fitting_info), function(x){
  if(x$fitting_type == "World"){
    world_fit
  } else if(x$fitting_type == "Country"){
    country_fits[[x$iso3c]]
  } else if(x$fitting_type == "Sub Region"){
    sub_region_fits[[x$un_sub_region]]
  } else if(x$fitting_type == "Region"){
    region_fits[[x$un_region]]
  }
})
names(variant_timings) <- iso3cs

saveRDS(variant_timings, "variant_timings.Rds")
