
income_Rt <- function(date_0) {

  # parse for latest date
  url <- "https://github.com/mrc-ide/global-lmic-reports/tree/master/data"
  html <- xml2::read_html(url)
  links <- rvest::html_nodes(html, ".js-navigation-open.Link--primary")
  text <- rvest::html_text(links, "title")
  latest <- which.max(as.Date(head(substr(text, 1, 10), -1)))
  d_link <- paste0("https://github.com/mrc-ide/global-lmic-reports/raw/master/data/",text[latest])

  # download and read in data file
  tf <- tempfile()
  suppressMessages(download.file(d_link, tf, quiet = TRUE))
  tf2 <- tempfile()
  extract <- unzip(tf, exdir = tf2)
  dat <- as.data.frame(data.table::fread(extract))

  # subset to Rt values
  dat <- dat[dat$compartment == "Rt" & dat$scenario == "Maintain Status Quo",]
  dat$income <- squire.page::get_income_group(dat$iso3c)
  filt_date <- min(date_0, max(dat$date))
  income_rt <- group_by(dat %>% filter(date == filt_date), income) %>% summarise(Rt = mean(y_mean))
  return(income_rt)
}


income_R0 <- function() {
  wb_metadata <- read.csv("gdp_income_group.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)

  R0s <- c(AFG = 1.94681197654148, AGO = 1.61007163196684, ALB = 1.96763103318809,
           ARE = 2.50253413533096, ARG = 2.75871385432299, ARM = 3.22817027616877,
           ATG = 2.4613728783852, AUS = 1.844035890742, AUT = 2.63323525599435,
           AZE = 2.55504023020025, BDI = 1.87623494135328, BEL = 3.20242750866164,
           BEN = 1.69995516225645, BFA = 2.03584781926971, BGD = 2.26914584142685,
           BGR = 1.96372896156003, BHR = 2.66280471201419, BHS = 2.92948926838559,
           BIH = 2.20660321948102, BLR = 2.12087726238117, BLZ = 3.98731571431526,
           BOL = 2.65462386446251, BRA = 3.35790560869017, BRB = 3.22416560168968,
           BRN = 4.29093542863153, BTN = 2.62477398120918, BWA = 3.07113496710246,
           CAF = 2.35794847845685, CAN = 2.96655337671002, CHE = 2.92733794711689,
           CHL = 3.58909774924995, CHN = 4.36822374260214, CIV = 1.98787174032542,
           CMR = 2.14549213886904, COD = 2.21116650066373, COG = 4.10238375616632,
           COL = 2.79077829877009, COM = 1.99125818369441, CPV = 1.76940551504754,
           CRI = 3.99413824030293, CUB = 2.54351237003809, CYP = 2.33473835346989,
           CZE = 3.30720884685794, DEU = 4.53855743696443, DJI = 1.70595258323997,
           DNK = 2.76264281150239, DOM = 2.7648322484037, DZA = 2.33881363433516,
           ECU = 3.03023855624653, EGY = 2.35744999249924, ERI = 3.39014842143898,
           ESP = 3.10464499150531, EST = 2.42305701041133, ETH = 2.01609626269111,
           FIN = 4.01051380376075, FJI = 2.70213252496874, FRA = 2.82518271889552,
           GAB = 2.54335510763381, GBR = 3.10616184489041, GEO = 2.1877677124125,
           GHA = 1.78844534217843, GIN = 2.58945541041892, GMB = 2.09904868844098,
           GNB = 3.8016943777235, GNQ = 2.1911241113226, GRC = 2.2470272203058,
           GRD = 2.90248979859742, GTM = 1.97984313982912, GUF = 2.6821201467294,
           GUY = 2.26564488915524, HKG = 2.54670715898017, HND = 2.71982595692518,
           HRV = 2.0074705489574, HTI = 3.09777895923623, HUN = 2.51314386302673,
           IDN = 2.3163696484572, IND = 3.00255927889085, IRL = 2.80712551223447,
           IRN = 3.51489818616122, IRQ = 2.28583764086555, ISL = 3.19597949565307,
           ISR = 3.05805319252858, ITA = 3.31298432164246, JAM = 3.02271712854222,
           JOR = 2.36607781951645, JPN = 2.12976088501996, KAZ = 2.70664459437929,
           KEN = 2.25569491579805, KGZ = 3.02620712531101, KHM = 3.77732364079532,
           KOR = 3.1257066910769, KWT = 2.96381415836947, LAO = 3.26126173664687,
           LBN = 1.94867905965103, LBR = 2.98352105956304, LBY = 1.69941521798719,
           LCA = 2.62260767845799, LKA = 2.40389495706905, LSO = 3.08772003655257,
           LTU = 2.50541574983208, LUX = 2.44278340646689, LVA = 1.94147700505066,
           MAR = 2.82596368511724, MDA = 2.43155455077294, MDG = 2.28192164307764,
           MDV = 2.63304801841058, MEX = 2.92179761664674, MKD = 5.12090391036425,
           MLI = 2.24740888084917, MLT = 2.96556528387917, MMR = 3.5553485736226,
           MNE = 3.49198279700695, MNG = 3.6002021677823, MOZ = 2.89698831663373,
           MRT = 2.79875462868845, MUS = 4.16921814519984, MWI = 2.60945091543225,
           MYS = 2.51215741924466, NAM = 3.84374660756451, NER = 1.85994735802755,
           NGA = 2.39687673804528, NIC = 2.45447448075758, NLD = 2.92658176336903,
           NOR = 2.59317735344656, NPL = 2.08650179431455, NZL = 2.12854821640393,
           OMN = 2.4807487638753, PAK = 2.40877929934856, PAN = 2.55544128138134,
           PER = 3.9299165943542, PHL = 1.83730758332291, PNG = 2.26177757599129,
           POL = 2.51008357140915, PRT = 2.97631636852076, PRY = 3.09753893131212,
           PSE = 2.09458916546166, QAT = 2.22672750732284, ROU = 2.96374255744107,
           RUS = 2.37746795980324, RWA = 2.06489410624638, SAU = 2.74011647848861,
           SDN = 2.1175149397001, SEN = 2.14953784249528, SGP = 2.59016673257639,
           SLE = 3.54783639015446, SLV = 4.33788960528107, SOM = 1.66303921637251,
           SRB = 2.51119258901574, SSD = 2.83342659159218, STP = 2.80832893411712,
           SUR = 2.75450148412546, SVK = 2.0614957459155, SVN = 2.41043614357078,
           SWE = 3.58787970888471, SWZ = 3.67394060609172, SYC = 3.01937945039727,
           SYR = 1.75617667335523, TCD = 2.13429521935772, TGO = 2.04553313964911,
           THA = 2.44888779825846, TJK = 3.08353013739756, TLS = 3.06982806631476,
           TTO = 3.34905776289758, TUN = 3.14937614855191, TUR = 2.96807878049415,
           TWN = 2.34753106748836, TZA = 1.73654113918716, UGA = 2.32953149631786,
           UKR = 2.24977562051465, URY = 3.00467163271078, USA = 3.03720462901783,
           UZB = 2.20365838608623, VCT = 2.65146649800355, VEN = 2.534477035193,
           VNM = 2.72443455557906, VUT = 1.79761052578649, YEM = 2.53331036695775,
           ZAF = 2.28563969296165, ZMB = 1.93282569548775, ZWE = 2.15991499009875,
           SLB = 3.73118501420733)
  inc_r0s <- data.frame("R0" = R0s, "iso3c" = names(R0s), "income" = as.character(squire.page::get_income_group(names(R0s))))
  R0s_income <- group_by(inc_r0s, income) %>% summarise(R0 = median(R0))

  return(R0s_income)

}

init_state <- function(deaths_removed, iso3c, seeding_cases = 5) {

  # get an initial
  pop <- squire::get_population(iso3c = iso3c, simple_SEIR = FALSE)
  init <- squire:::init_check_explicit(NULL, pop$n, seeding_cases = seeding_cases)

  if(deaths_removed > 0) {

    # work out how many deaths and where
    probs <- (squire:::probs$prob_hosp * squire:::probs$prob_severe * squire:::probs$prob_severe_death_treatment) +
      (squire:::probs$prob_hosp * (1-squire:::probs$prob_severe * squire:::probs$prob_non_severe_death_treatment))
    probs <- probs*pop$n
    probs <- probs/sum(probs)
    deaths <- as.numeric(t(rmultinom(1, deaths_removed, probs)))

    # approximate IFR for income group
    wb_metadata <- read.csv("gdp_income_group.csv", fileEncoding="UTF-8-BOM", stringsAsFactors = TRUE)
    income <- wb_metadata$income_group[match(iso3c, wb_metadata$country_code)]
    ifrs <- data.frame("income" = c("Low income", "Lower middle income", "Upper middle income", "High income"),
                       "ifr" = c(0.17, 0.31, 0.51, 1.02))
    ifr <- ifrs$ifr[ifrs$income == income]
    R <- rpois(1, deaths_removed*1/ifr/0.01)
    R <- as.numeric(t(rmultinom(1, R, rep(1/length(probs), length(probs)))))
    R <- R - deaths

    # and update the inital to reflect
    init$D <- deaths
    init$S <- init$S - R - deaths
    init$R1 <- R

  }

  return(init)
}

