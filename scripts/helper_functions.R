######### Regression Functions ########

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

listing <- function(...) {
  rlist::list.flatten(rlist::list.append(...))
}

create_col <- function(.data, lol) {
  lol <- enquo(lol)
  
  if (any(.data %>% 
          select(!!lol) %>% 
          str_detect(., "engage10"))) {
    mutate(.data, type = "engage10")
  } else if (any(.data %>% 
                 select(!!lol) %>% 
                 str_detect(., "delib10"))) {
    mutate(.data, type = "delib10")
  } else if (any(.data %>% 
                 select(!!lol) %>% 
                 str_detect(., "consult10"))) {
    mutate(.data, type = "consult10")
  } else if (any(.data %>% 
                 select(!!lol) %>% 
                 str_detect(., "reason10"))) {
    mutate(.data, type = "reason10")
  } else if (any(.data %>% 
                 select(!!lol) %>% 
                 str_detect(., "common10"))) {
    mutate(.data, type = "common10")
  } else if (any(.data %>% 
                 select(!!lol) %>% 
                 str_detect(., "countr10"))) {
    mutate(.data, type = "countr10")
  } else{
    
    .data <- 2
    
  }
}

get_p_values <- function (pval) {
  dplyr::case_when(
    is.na(pval) ~ "", 
    pval < 0.001 ~ "***", 
    pval < 0.01 ~ "**", 
    pval < 0.05 ~ "*",     
    pval < 0.10 ~ "^", 
#    pval < 0.10 ~ expression(phantom(0)^"\u2020"), 
    TRUE ~ "")
}

prep_dat <- function(model, bias) {
  model %>% 
    purrr::map(broom::tidy) %>%
    purrr::map(~create_col(., term)) %>% 
    bind_rows() %>% 
    mutate(bias = bias) %>% 
    na.omit() %>% 
    left_join(model %>% 
                purrr::map(p_value) %>% 
                purrr::map(~create_col(., term)) %>% 
                bind_rows()) %>% 
    left_join(model %>% 
                purrr::map(icc) %>% 
                purrr::map("cntry") %>% 
                unlist() %>% 
                data.frame(icc = ., type = c("delib10", "consult10", "reason10", "common10", "countr10", "engage10")))
}

prep_reml <- function(model, bias) {
  model %>% 
    purrr::map(broom::glance) %>% 
    bind_rows() %>% 
    mutate(type = c("delib10", "consult10", "reason10", 
                    "common10", "countr10", "engage10")) %>% 
    mutate(bias = bias)
}

clean_icc_with_controls <- function(with_controls) {
  with_controls %>% 
    select(term, estimate, type, bias, icc) %>% 
    group_by(type, bias) %>% 
    #  slice(1) %>% 
    mutate(icc = round(icc, 2)) %>% 
    filter(term %nin% c("sex","income","educ", 
                        "work", "age", "sex",
                        "afro", "latino", "americas", "asian", "ESS",
                        "ethnic10", "pop10", "gdp10", "lifeexp10", 
                        "corrupt10", "urbanratio10",
                        "(Intercept)")) %>% 
    mutate(term = case_when(
      term == "delib10" ~ "Deliberation Component Index",
      term == "reason10" ~ "Reasoned Justification",
      term == "common10" ~ "Common Good",
      term == "countr10" ~ "Counter-Arguments",
      term == "consult10" ~ "Range of Consultation",
      term == "engage10" ~ "Engaged Society"
    )) %>% 
    mutate(term = factor(term, 
                         levels = rev(c("Deliberation Component Index",
                                        "Reasoned Justification",
                                        "Common Good",
                                        "Counter-Arguments",
                                        "Range of Consultation",
                                        "Engaged Society"))))
}

clean_icc_with_pols <- function(with_pols) {
  with_pols %>% 
    select(term, estimate, type, bias, icc) %>% 
    group_by(type, bias) %>% 
    #  slice(1) %>% 
    mutate(icc = round(icc, 2)) %>% 
    filter(term %nin% c("sex","income","educ", 
                        "work", "age", "sex",
                        "afro", "latino", "americas", "asian", "ESS",
                        "ethnic10", "pop10", "gdp10", "lifeexp10", 
                        "corrupt10", "urbanratio10",
                        "(Intercept)", "polity_demdummy", "polity_autodummy")) %>% 
    mutate(term = case_when(
      term == "delib10" ~ "Deliberation Component Index",
      term == "reason10" ~ "Reasoned Justification",
      term == "common10" ~ "Common Good",
      term == "countr10" ~ "Counter-Arguments",
      term == "consult10" ~ "Range of Consultation",
      term == "engage10" ~ "Engaged Society",
    )) %>% 
    mutate(term = factor(term, 
                         levels = rev(c("Deliberation Component Index",
                                        "Reasoned Justification",
                                        "Common Good",
                                        "Counter-Arguments",
                                        "Range of Consultation",
                                        "Engaged Society")))) 
}

get_anovies <- function(models1, models2, model_number) {
  anova(models1[[model_number]], models2[[model_number]]) -> anovies
  
  anovies %>% 
    data.frame() %>% 
    rownames_to_column("comparison") %>% 
    janitor::clean_names() %>% 
    mutate(stars = get_p_values(pr_chisq)) %>% 
    mutate(report = paste0("(", chi_df, ")==", sprintf('%.2f', chisq, 2))) %>% 
    select(comparison, report, stars) %>% 
    .[2,]
}


get_anovies_pol <- function(models1, models2, model_number) {
  anova(models1[[model_number]], models2) -> anovies
  
  anovies %>% 
    data.frame() %>% 
    rownames_to_column("comparison") %>% 
    janitor::clean_names() %>% 
    mutate(stars = get_p_values(pr_chisq)) %>% 
    mutate(report = paste0("(", chi_df, ")==", sprintf('%.2f', chisq, 2))) %>% 
    select(comparison, report, stars) %>% 
    .[2,]
}

######### plot_models2 ########

plot_models2 <- function(..., exponentiate, std.est = NULL, rm.terms = NULL, 
                         title = NULL, m.labels = NULL, legend.title = "Dependent Variables", 
                         legend.pval.title = "p-level", axis.labels = NULL, axis.title = NULL, 
                         axis.lim = NULL, wrap.title = 50, wrap.labels = 25, wrap.legend.title = 20, 
                         grid.breaks = NULL, geom.size = 3, geom.spacing = 0.4, geom.colors = "Set1", 
                         show.values = FALSE, show.legend = TRUE, show.intercept = FALSE, 
                         show.p = TRUE, p.shape = FALSE, vline.type = 2, vline.color = "grey70", 
                         digits = 2, facet.grid = FALSE, thatorder = NULL, plot_colors = c("green","green","green","green","green","green",
                                                                                           "blue","blue","blue","blue","blue","blue",
                                                                                           "red","red","red","red","red","red"), col_mods = NULL, labelv = -0.5, labelh = -1.1, labelsize = 2) 
{
  
  is_merMod <- function(fit) {
    inherits(fit, c("lmerMod", "glmerMod", "nlmerMod", "merModLmerTest"))
  }
  
  #' @importFrom sjmisc str_contains
  get_glm_family <- function(fit) {
    c.f <- class(fit)
    
    # do we have glm? if so, get link family. make exceptions
    # for specific models that don't have family function
    if (any(c.f %in% c("lme", "plm"))) {
      fitfam <- ""
      logit_link <- FALSE
    } else {
      fitfam <- stats::family(fit)$family
      logit_link <- stats::family(fit)$link == "logit"
    }
    
    # create logical for family
    binom_fam <- fitfam %in% c("binomial", "quasibinomial")
    poisson_fam <- fitfam %in% c("poisson", "quasipoisson") ||
      sjmisc::str_contains(fitfam, "negative binomial", ignore.case = T)
    
    list(is_bin = binom_fam, is_pois = poisson_fam, is_logit = logit_link)
  }
  
  # return names of objects passed as ellipses argument
  dot_names <- function(dots) unname(unlist(lapply(dots, as.character)))
  
  
  #' @importFrom tidyr nest
  #' @importFrom dplyr select filter
  #' @importFrom stats complete.cases
  #' @importFrom rlang .data
  get_grouped_data <- function(x) {
    # nest data frame
    grps <- tidyr::nest(x)
    
    # remove NA category
    cc <- grps %>%
      dplyr::select(-.data$data) %>%
      stats::complete.cases()
    
    # select only complete cases
    grps <- grps %>% dplyr::filter(!! cc)
    
    # arrange data
    if (length(attr(x, "vars", exact = T)) == 1)
      reihe <- order(grps[[1]])
    else
      reihe <- order(grps[[1]], grps[[2]])
    grps <- grps[reihe, ]
    
    grps
  }
  
  
  str_ends_with <- function(x, .match) {
    l <- nchar(x)
    n <- nchar(.match)
    m <- substr(x, pmax(1, l - n + 1), l)
    which(m == .match)
  }
  
  input_list <- tibble::lst(...)
  if (length(input_list) == 1 && class(input_list[[1]]) == 
      "list") 
    input_list <- purrr::map(input_list[[1]], ~.x)
  if (missing(exponentiate)) 
    exponentiate <- inherits(input_list[[1]], c("glm", "glmerMod", 
                                                "glmmTMB"))
  if (!any(inherits(input_list[[1]], c("lm", "lmerMod", "lme", "merModLmerTest"), 
                    which = TRUE) == 1)) 
    std.est <- NULL
  
  #purrr::map(input_list, ~summary(.x))
  #summary(input_list)
  #s$coefficients[,5]
  #if (!show.intercept) 
  #  fl <- purrr::map(fl, ~dplyr::slice(.x, -1))
  #if (!is.null(std.est)) {
  #  
  get_p_stars <- function(pval) {
    dplyr::case_when(
      is.na(pval) ~ "",
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      pval < 0.10 ~ "^",
      TRUE ~ ""
    )
  }
  
  #    fl <- fl %>% purrr::map(~sjstats::std_beta(.x, 
  #                                 type = std.est)) %>% 
  #                    purrr::map(~sjmisc::var_rename(.x, 
  #                   std.estimate = "estimate")) %>% 
  #      purrr::map2(input_list,  ~tieele::add_column(.x, 
  #                    p.value = sjstats::get_model_pval(.y)[["p.value"]][-1]))
  #  }
  input_list2 <- input_list %>% 
    purrr::map(~sjstats::std_beta(.x, type = std.est)) %>% 
    purrr::map(~sjmisc::var_rename(.x,  std.estimate = "estimate")) 
  
  input_list2  <- purrr::map(input_list2 , ~dplyr::slice(.x, -1))
  
  nineties <- function(x) {
    x$conf.high2 <- x$estimate + 1.645 * x$std.error
    x$conf.low2 <- x$estimate - 1.645 * x$std.error
    s <- tibble(x$conf.high2,x$conf.low2)
    names(s) <- c("conf.high2","conf.low2")
    s
  }
  
  input_list2 <- purrr::map(input_list2 , ~nineties(.x)) %>%  purrr::map2(input_list2, cbind) 
  
  
  fl <- input_list2 %>% 
    purrr::map2(input_list, ~tibble::add_column(.x, p.value = sjstats::p_value(.y)[["p.value"]][-1]))
  
  get_axis_limits_and_ticks <- function(axis.lim, min.val, max.val, grid.breaks, exponentiate) {
    # axis limits
    if (is.null(axis.lim)) {
      lower_lim <- min.val * .95
      upper_lim <- max.val * 1.05
    } else {
      lower_lim <- axis.lim[1]
      upper_lim <- axis.lim[2]
    }
    
    # determine gridbreaks
    if (is.null(grid.breaks)) {
      if (exponentiate) {
        # use pretty distances for log-scale
        ticks <- grDevices::axisTicks(log(c(lower_lim, upper_lim)), log = TRUE)
        # make sure that scale is not too wide
        ticks <- ticks[1:which(ticks > max.val)[1]]
      } else {
        ticks <- pretty(c(lower_lim, upper_lim))
      }
    } else {
      ticks <- seq(lower_lim, upper_lim, by = grid.breaks)
    }
    
    # save proper axis limits
    list(axis.lim = c(min(ticks), max(ticks)), ticks = ticks)
  }
  #s$coefficients[,5]
  if (!show.intercept) 
    fl <- purrr::map(fl, ~dplyr::slice(.x, -1))
  
  if (exponentiate) 
    fl <- purrr::map(fl, function(x) {
      x[["estimate"]] <- exp(x[["estimate"]])
      x[["conf.low"]] <- exp(x[["conf.low"]])
      x[["conf.high"]] <- exp(x[["conf.high"]])
      x
    })
  for (i in 1:length(fl)) fl[[i]] <- tibble::add_column(fl[[i]], 
                                                        group = as.character(i))
  ff <- dplyr::bind_rows(fl)
  if (!is.null(std.est) && std.est == "std2") 
    ff$term <- substring(ff$term, first = 3)
  if (!is.null(rm.terms)) 
    ff <- subset(ff,!(ff$term %in% rm.terms))
  if (is.null(m.labels)) 
    m.labels <- sjlabelled::get_dv_labels(input_list)
  m.labels <- sjmisc::word_wrap(m.labels, wrap = wrap.labels)
  if (anyDuplicated(m.labels) > 0) 
    m.labels <- suppressMessages(tibble::tidy_names(m.labels))
  ff$group <- as.factor(ff$group)
  levels(ff$group) <- m.labels
  #  ff$group <- forcats::fct_rev(ff$group)
  
  if (!is.null(thatorder)) ff$term <- forcats::fct_relevel(ff$term, thatorder)
  
  ff$term <- forcats::fct_rev(ff$term)
  
  ff$p.stars <- get_p_stars(ff$p.value)
  ff$p.label <- sprintf("%.*f", digits, ff$estimate)
  if (show.p) 
    ff$p.label <- sprintf("%s %s", ff$p.label, ff$p.stars)
  if (p.shape) 
    p <- ggplot(ff, aes_string(x = "term", y = "estimate", 
                               colour = "group", shape = "p.stars"))
  else p <- ggplot(ff, aes_string(x = "term", y = "estimate", 
                                  colour = "group"))
  p <- p + geom_hline(yintercept = 0, 
                      linetype = vline.type, 
                      color = vline.color) + 
    geom_point(position = position_dodge(geom.spacing), size = geom.size) + 
    geom_errorbar(aes_string(ymin = "conf.low2", 
                             ymax = "conf.high2"), 
                  position = position_dodge(geom.spacing), width = 0) + 
    coord_flip() + 
    guides(colour = guide_legend(reverse = TRUE))
  
  get_estimate_axis_title <- function(fit, axis.title) {
    # check if we have a linear model
    is.linear.model <- any(inherits(fit, c("lm", "lmerMod", "lme"), which = TRUE) == 1)
    
    if (!is.linear.model) {
      # get information of glm
      fitfam <- get_glm_family(fit)
      
      # create logical for family
      poisson_fam <- fitfam$is_pois
      binom_fam <- fitfam$is_bin
      logit_link <- fitfam$is_logit
    }
    
    # check default label and fit family
    if (is.null(axis.title)) {
      if (is.linear.model)
        axis.title <- "Estimates"
      else if (poisson_fam)
        axis.title <- "Incident Rate Ratios"
      else if (binom_fam && !logit_link)
        axis.title <- "Risk Ratios"
      else
        axis.title <- "Odds Ratios"
    }
    
    axis.title
  }
  
  col_check2 <- function(geom.colors, collen) {
    # --------------------------------------------
    # check color argument
    # --------------------------------------------
    # check for corrct color argument
    if (!is.null(geom.colors)) {
      # check for color brewer palette
      #if (is.brewer.pal(geom.colors[1])) {
      if (TRUE) {
        geom.colors <- scales::brewer_pal(palette = geom.colors[1])(collen)
      } 
      else if (geom.colors[1] == "gs") {
        geom.colors <- scales::grey_pal()(collen)
        # do we have correct amount of colours?
      } else if (geom.colors[1] == "bw") {
        geom.colors <- rep("black", times = collen)
        # do we have correct amount of colours?
      } else if (length(geom.colors) > collen) {
        # shorten palette
        geom.colors <- geom.colors[1:collen]
      } else if (length(geom.colors) < collen) {
        # warn user abount wrong color palette
        warning(
          sprintf(
            "Insufficient length of color palette provided. %i color values needed.",
            collen
          ),
          call. = F
        )
        # set default palette
        geom.colors <- scales::brewer_pal(palette = "Set1")(collen)
      }
    } else {
      geom.colors <- scales::brewer_pal(palette = "Set1")(collen)
    }
    
    geom.colors
  }
  
  if (p.shape) 
    p <- p + scale_shape_manual(values = c(1, 16, 17, 15), 
                                labels = c("n.s.", "*", "**", "***"))
  if (show.values) 
    p <- p + geom_text(aes_string(label = "p.label"), position = position_dodge(geom.spacing), 
                       vjust = geom.spacing * labelv, hjust = labelh, show.legend = FALSE, size = labelsize)
  if (is.null(axis.labels)) 
    axis.labels <- sjlabelled::get_term_labels(input_list)
  p <- p + scale_x_discrete(labels = sjmisc::word_wrap(axis.labels, 
                                                       wrap = wrap.labels))
  if (!show.legend) 
    p <- p + guides(colour = "none", shape = "none")
  if (facet.grid) 
    p <- p + facet_grid(~group)
  axis.scaling <- get_axis_limits_and_ticks(axis.lim = axis.lim, 
                                            min.val = min(ff$conf.low2), max.val = max(ff$conf.high2), 
                                            grid.breaks = grid.breaks, exponentiate = exponentiate)
  #if (exponentiate) {
  #  p <- p + scale_y_continuous(trans = "log10", limits = axis.scaling$axis.lim, 
  #                              breaks = axis.scaling$ticks, labels = prettyNum)
  #}
  #else {
  #  p <- p + scale_y_continuous(limits = axis.scaling$axis.lim, 
  #                              breaks = axis.scaling$ticks, 
  #                              labels = axis.scaling$ticks)
  p <- p + scale_y_continuous(limits = axis.scaling$axis.lim)
  #}
  # p <- p + scale_colour_manual(values = col_check2(geom.colors, 
  #                                                   length(m.labels)))
  # p <- p + scale_colour_manual(values = palette(rainbow(18)))
  if (!is.null(col_mods)) {
    p <- p + scale_colour_manual(values = c(rep(col_mods[1], 4), 
                                            col_mods[2], 
                                            rep(col_mods[1], 2), 
                                            rep(col_mods[2], 5), 
                                            rep(col_mods[3], 6)))  
  } else {
    p <- p + scale_colour_manual(values = plot_colors)
  }
  
  p <- p + labs(x = NULL, 
                y = sjmisc::word_wrap(
                  get_estimate_axis_title(
                    input_list[[1]], 
                    axis.title), 
                  wrap = wrap.title), 
                title = sjmisc::word_wrap(title, 
                                          wrap = wrap.title), 
                colour = sjmisc::word_wrap(legend.title, 
                                           wrap = wrap.legend.title), 
                shape = sjmisc::word_wrap(legend.pval.title, 
                                          wrap = wrap.legend.title))
  p
}



######### W/O Controls ########

multilevels_no <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

multilevels_low <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_low ~ income + educ + 
                              work + age + sex + i +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}


multilevels_high <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_high ~ income + educ + 
                              work + age + sex + i +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

######### With Controls ########

multilevels_no_c <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}



multilevels_low_c <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_low ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

multilevels_high_c <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_high ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

multilevels_no_c_dem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_dem, weights = weight, ...)
  return(result)
}

multilevels_no_c_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}

multilevels_low_c_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_low ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}


multilevels_high_c_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_high ~ income + educ + 
                              work + age + sex + i +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}



######### With Polity Controls ########

multilevels_no_pol <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i + polity_demdummy + polity_autodummy +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

multilevels_low_pol <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_low ~ income + educ + 
                              work + age + sex + i + polity_demdummy + polity_autodummy +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

multilevels_high_pol <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_high ~ income + educ + 
                              work + age + sex + i + polity_demdummy + polity_autodummy +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}



multilevels_no_pol_dem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i + polity10 +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_dem, weights = weight, ...)
  return(result)
}

multilevels_no_pol_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i + polity10 +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}

multilevels_low_pol_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_low ~ income + educ + 
                              work + age + sex + i + polity10 +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}


multilevels_high_pol_nondem <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust_high ~ income + educ + 
                              work + age + sex + i + polity10 +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined_nondem, weights = weight, ...)
  return(result)
}


######### Polity Models ########

multilevels_no_pol <- function(delibs, ...) {
  result <- lmer(substitute(gov_trust ~ income + educ + 
                              work + age + sex + i + polity_demdummy + polity_autodummy +
                              pop10 + gdp10 + lifeexp10 + corrupt10 + urbanratio10 +
                              afro + latino + americas + asian + ESS +
                              (1|cntry), list(i = as.name(delibs))) , data = combined, weights = weight, ...)
  return(result)
}

