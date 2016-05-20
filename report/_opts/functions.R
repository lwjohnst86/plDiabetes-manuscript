
##' Functions
##' =========
##'
##' Custom functions I use for my analysis (those that don't warrant
##' their own package developed).
##'

med <- function(x) {
    round(median(x, na.rm = TRUE), 1)
}

iqr <- function(x) {
    paste0(round(quantile(x, 0.25, na.rm = TRUE), 1),
           '-', round(quantile(x, 0.75, na.rm = TRUE), 1))
}

correlate <- function(x, y) {
    round(cor(x, y, use = "complete.obs",
        method = 'spearman'), 2)
}

medianIQR <- function(x) {
    paste0(med(x), ' (', iqr(x), ')')
}

average <- function(x) {
    round(mean(x, na.rm = TRUE), 1) %>% format(., nsmall = 1)
}

stddev <- function(x) {
    round(sd(x, na.rm = TRUE), 1) %>% format(., nsmall = 1)
}

meanSD <- function(x) {
    paste0(average(x), ' (', stddev(x), ')')
}

anovaPval <- function(continuous, factor) {
    ## Uses dplyr
    fit <- aov(continuous ~ factor) %>% summary() %>% unlist()
    pval <- ifelse(fit['Pr(>F)1'] <= 0.01, '<0.01',
                   paste0(fit['Pr(>F)1'] %>% round(2) %>% format(., nsmall = 2)))
    return(pval)
}

renamePool <- function(x) {
    x %>%
      gsub('CE|ce', 'Cholesteryl esters', .) %>%
      gsub('PL|pl', 'Phospholipids', .)
}

renameSpecies <- function(x) {
    x %>%
      gsub('^pct_', '', .) %>%
      gsub('_Q', 'q', .) %>%
      gsub('_', '.', .) %>%
      gsub('^(\\D\\D)(\\d\\d)(\\d)', '\\1.\\2:\\3', .) %>%
      gsub('n(\\d)$', 'n-\\1', .) %>%
      gsub('D(\\d\\d)$', 'D-\\1', .) %>%
      gsub('indexDNL', 'DNL Index', .) %>%
      gsub('Total', 'Total (mole)', .)
}

splitSpeciesPools <- function(data, fa.variable) {
    data %>%
      mutate_(fa = lazyeval::interp(~renameSpecies(x),
                                    x = as.name(fa.variable))) %>%
      separate(fa, c('pool', 'species'), sep = '\\.') %>%
      mutate(species = species %>% factor(., levels = unique(.)),
             pool = renamePool(pool) %>% factor(., levels = unique(.)))
}

geePval <- function(continuous, time, data = ds, round.values = TRUE) {
    data <- data %>%
      select_(continuous = continuous, time = time, 'SID')

    pvalue <- geeglm(continuous ~ as.numeric(time), data = data,
                     id = SID, corstr = 'ar1') %>%
      tidy()

    if (round.values) {
        ifelse(pvalue[2,5] < 0.01, '<0.01', format(round(pvalue[2,5], 2), nsmall = 2))
    } else {
        pvalue[2,5]
    }
}

tertile <- function(x) {
    ## Uses dplyr
    x %>%
      quantile(., c(0, .333, .666, 1), na.rm = TRUE) %>%
      cut(x, ., include.lowest = TRUE)
}

chiSqPval <- function(factor1, factor2) {
    ## Uses dplyr and broom packages.
    table(factor1, factor2, useNA = 'no') %>%
      chisq.test() %>%
      tidy() %>%
      mutate(p = ifelse(p.value <= 0.05, '<0.01',
                        paste0(format(round(p.value, 2), nsmall = 2)))) %>%
      head(1)['p']
}


geeDF <- function(data, exposure, outcomes, covariates, renaming, pool.unit = '', p.adj = TRUE) {
    data %>%
      loopGEE(., outcomes, exposure, 'SID', covariates, corstr = 'ar1',
              filter.indep = TRUE, adjust.p.value = p.adj) %>%
      mutate(estimate = (exp(estimate) - 1) * 100,
             conf.low = (exp(conf.low) - 1) * 100,
             conf.high = (exp(conf.high) - 1) * 100) %>%
      mutate(dep = factor(dep, levels = unique(dep),
                          labels = renaming)) %>%
      splitSpeciesPools(., 'indep') %>%
      select(-indep, -term) %>%
      rename(indep = species) %>%
      mutate(pool = paste(pool, pool.unit) %>% factor(., levels = unique (.)))
}

plotForestGEE <- function(data, ylab, xlab.var.types = 'fatty acids.') {
    data %>%
      plotForest(., pvalue.factor = 'f.pvalue', groups = 'pool ~ dep',
                 x.axis.label = paste0('Percent difference with 95% CI in the outcomes\nfor each SD increase in ', xlab.var.types),
                 y.axis.label = ylab) +
      theme_tufte(base_family = 'Arial') +
      theme(axis.text.y = element_text(face=ifelse(grepl('Total', data$indep), "bold", "plain")),
            legend.key.width = unit(0.7, "line"),
            legend.key.height = unit(0.7, "line"),
            strip.background = element_blank(), #element_rect(fill = 'grey90', colour = 'grey90'),
            plot.margin = unit(c(0.5, 0, 0, 0), "cm")
            ) +
      scale_x_continuous(breaks = c(-10, -5, 0, 5, 10))
}

theme_simple <- function(base_size = 12, base_family = "Arial")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
            strip.background = element_rect(fill = 'grey90', colour = 'grey90'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = 'grey95'),
            panel.border = element_rect(fill = NA, colour = 'grey80', linetype = 'solid')
           )
  }

theme_nothing <- function(base_size = 12, base_family = "Arial")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
            rect             = element_blank(),
            line             = element_blank(),
            text             = element_blank(),
            axis.ticks.margin = unit(0, "lines")
           )
  }

geeInteractDF <- function(data, outcomes, exposures,
                             covariates, interaction,
                             renaming, pool.unit = '') {
    data %>%
      loopGEE(., outcomes, exposures, 'SID', covariates, interaction, corstr = 'ar1',
              filter.interact = TRUE, adjust.p.value = TRUE) %>%
      mutate(dep = factor(dep, levels = unique(dep),
                          labels = renaming)) %>%
      splitSpeciesPools(., 'indep') %>%
      select(-indep) %>%
      rename(indep = species) %>%
      mutate(pool = paste(pool, pool.unit) %>% factor(., levels = unique (.)))
}

plotManhattanGEE <- function(data, ylab) {
    data %>%
        mutate(fraction = factor(fraction, unique(fraction), ordered = TRUE)) %>%
      plotManhattanStyle(., 'Xterms', 'p.value', groups = 'fraction ~ Yterms',
                     y.axis.label = ylab) +
      theme_tufte(10, base_family = 'Arial') +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),
            strip.background = element_rect(fill = 'grey95', colour = 'grey95'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.y = element_text(face=ifelse(grepl('Total', levels(data$Xterms)), "bold", "plain"))) +
      theme(legend.key.width = unit(0.85, "line"),
            legend.key.height = unit(0.85, "line"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            )
}

table_geeDF <- function(data) {
    data %>%
      mutate(estCI = paste0(format(round(estimate, 1), nsmall = 1), ' (',
                            trim(format(round(conf.low, 1), nsmall = 1)), ', ',
                            trim(format(round(conf.high, 1), nsmall = 1)), ')',
                            ifelse(p.value < 0.05, '*', ' '))) %>%
      select(orderno, Yterms, Xterms, estCI) %>%
      spread(Yterms, estCI) %>%
      arrange(orderno, desc(Xterms)) %>%
      select(-orderno) %>%
      rename(Species = Xterms) %>%
      pander(justify = c('left', rep('right', 4)))
}

tableRowDescriptiveStat <- function(data.subset, func = meanSD) {
    data.subset %>%
      select(group = which(sapply(., is.factor)), which(sapply(., is.numeric))) %>%
      group_by(group) %>%
      summarise_each(funs(func)) %>%
      gather(Measure, val, -group) %>%
      spread(group, val)
}

tableRowSampleSize <- function(data, row.variable, column.variable) {
    data %>%
      select_(group = column.variable, row.variable) %>%
      na.omit() %>%
      group_by(group) %>%
      summarise(n = n(),
                Measure = 'n') %>%
      spread(group, n)
}

tableRowNPercent <- function(data.subset, column.variable) {
    data.subset %>%
      rename_(group = column.variable) %>%
      gather(Demographic, val, -group) %>%
      group_by(group, Demographic, val) %>%
      tally(sort = TRUE) %>%
      group_by(group, Demographic) %>%
      mutate(npct = paste0(n, ' (', ((100*n) / sum(n)) %>% round(0), '%)')) %>%
      select(-n) %>%
      spread(group, npct) %>%
      rename(Measure = val)
}

tableColumnPvalueStat <- function(data.subset, group, func) {
    data.subset %>%
      rename_(group = group) %>%
      gather(Measure, val, -group) %>%
      group_by(Measure) %>%
      summarise(p = func(val, group))
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

tableModelSelection <- function(selection.data) {
    selection.data %>%
      add_rownames() %>%
      select(Model = rowname, QIC, Delta = delta) %>%
      pander(digits = 2)
}

figRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("figcap.prefix"),
        sep = options("figcap.sep"), prefix.highlight = options("figcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight,
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})

supFigRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("supfigcap.prefix"),
        sep = options("figcap.sep"), prefix.highlight = options("figcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight,
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})

supTabRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("suptabcap.prefix"),
        sep = options("tabcap.sep"), prefix.highlight = options("tabcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight,
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})

tabRef <- local({
    tag <- numeric()
    created <- logical()
    used <- logical()
    function(label, caption, prefix = options("tabcap.prefix"),
        sep = options("tabcap.sep"), prefix.highlight = options("tabcap.prefix.highlight")) {
        i <- which(names(tag) == label)
        if (length(i) == 0) {
            i <- length(tag) + 1
            tag <<- c(tag, i)
            names(tag)[length(tag)] <<- label
            used <<- c(used, FALSE)
            names(used)[length(used)] <<- label
            created <<- c(created, FALSE)
            names(created)[length(created)] <<- label
        }
        if (!missing(caption)) {
            created[label] <<- TRUE
            paste0(prefix.highlight, prefix, " ", i, sep, prefix.highlight,
                " ", caption)
        } else {
            used[label] <<- TRUE
            paste(prefix, tag[label])
        }
    }
})

checkFigTabNumbering <- function(ref.function) {
    if (!all(environment(ref.function)$used)) {
        missingRef <- which(!environment(ref.function)$used)
        warning("Figure(s) ", paste(missingRef, sep = ", "), " with label(s) '",
                paste(names(environment(ref.function)$used)[missingRef], sep = "', '"),
                "' are present in the document but are never referred to in the text.")
    }
    if (!all(environment(ref.function)$used)) {
        missingRef <- which(!environment(ref.function)$used)
        warning("Figure(s) ", paste(missingRef, sep = ", "), " with label(s) '",
                paste(names(environment(ref.function)$used)[missingRef], sep = "', '"),
                "' are present in the document but are never referred to in the text.")
    }
}

calcCumPctFA <- function(data, variable.pattern, filter.value) {
    data %>%
      filter(VN == 0) %>%
      select(matches(variable.pattern)) %>%
      gather(fat, value) %>%
      group_by(fat) %>%
      summarise(pct = mean(value, na.rm = TRUE)) %>%
      arrange(desc(pct)) %>%
      mutate(cumpct = cumsum(pct) %>% round(1) %>% format(nsmall = 1),
             pct = round(pct, 1) %>% format(nsmall = 1)) %>%
      filter(cumpct <= filter.value)
}

plotInteract <- function(data, y, x,
                         fatty.acid.levels = c(2, 3, 4),
                         covars = covariates,
                         title = 'A',
                         interaction.term = 'VN',
                         is.fatty.acid = TRUE,
                         y.axis.label = NULL) {
    fa.lvl <- fatty.acid.levels
    if (!fa.lvl %in% c(2, 3, 4)) {
        print('Please select either 2, 3, or 4 for the fatty.acid.levels')
    } else if (fa.lvl == 2) {
        lvls <- 0:2/2
        lvl.labs <- c('<50%', '>50%')
    } else if (fa.lvl == 3) {
        lvls <- 0:3/3
        lvl.labs <- c('<33%', '33-66%', '>66%')
    } else if (fa.lvl == 4) {
        lvls <- 0:4/4
        lvl.labs <- c('Lowest', 'Low', 'High', 'Highest')
    }

    if (is.null(y.axis.label)) {
        y.axis.label <- paste0('log(', substring(y, 2), ')')
    }

    fatty.acid.name <- renameSpecies(x) %>%
      gsub('(\\D\\D)\\.', '\\U\\1 ', ., perl = TRUE)

    if (is.fatty.acid) {
        fatty.acid.name <- paste(fatty.acid.name,
                                 ifelse(grepl('pct_', x), '(mol%)', '\n(nmol/mL)'))
    }

    geeFormula <- paste(y, 'f.fat', sep = ' ~ ') %>%
      paste(., paste(covars, collapse = ' + '),
            paste('f.fat', interaction.term, sep = ':'), sep = ' + ') %>%
      as.formula()

    prep.fit <- data %>%
      select_(.dots = c(y, x, covars, interaction.term, 'SID')) %>%
      filter(complete.cases(.)) %>%
      rename_('fa' = x, 'VN' = interaction.term) %>%
      ungroup() %>%
      mutate(f.fat = cut(fa, quantile(fa, probs = lvls),
                         include.lowest=TRUE, labels = lvl.labs),
             VN = factor(VN, ordered = TRUE))

    fit <-
        geeglm(geeFormula, id = SID, family = gaussian,
               corstr = 'ar1', data = prep.fit)

    dodge <- position_dodge(width = 0.1)

    lsmeans(fit, ~ f.fat:VN) %>%
      summary() %>%
      as.data.frame() %>%
      ggplot(aes(VN, lsmean, group = f.fat, shape = f.fat, ymax = max(asymp.UCL))) +
      geom_line(aes(linetype = f.fat), position = dodge) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, linetype = f.fat), width = 0.1,
                    position = dodge) +
      geom_point(position = dodge) +
      scale_shape_discrete('Percentile') +
      scale_linetype_discrete('Percentile') +
      scale_x_discrete('Visit number', labels = c('0-yr', '3-yrs', '6-yrs')) +
      theme_tufte(10, base_family = 'Arial') +
      theme(legend.position = 'bottom',
            strip.background = element_rect(fill = 'grey95', colour = 'grey95'),
            plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
            legend.key.width = unit(0.75, "line"),
            legend.key.height = unit(0.75, "line")
            ) +
      ## theme(legend.justification = c(1, 1),
      ##       legend.position = c(1, 1)) +
      ylab(y.axis.label) +
      ggtitle(paste0(title, ': ', fatty.acid.name)) +
      theme(plot.title = element_text(hjust = 0, vjust = 1))
}

# Analyze -----------------------------------------------------------------

analyze_gee <- function(data, y, x, covar, unit,
                        adj.p = TRUE, adj.p.method = 'BH', int = FALSE) {
    if (int) {
        extract_term <- ':'
    } else {
        extract_term <- 'Xterm$'
    }
    gee_results <- data %>%
        mason::design('gee', family = gaussian, corstr = 'ar1') %>%
        {
            if (int) {
                mason::lay_base(., 'SID', y, x, covar, intvar = 'VN')

            } else {
                mason::lay_base(., 'SID', y, x, covar)
            }
        } %>%
        mason::build() %>%
        mason::polish(
            extract_term, adjust.p = FALSE,
            transform.beta.funs = function(x)
                (exp(x) - 1) * 100,
            rename.vars.funs = renaming_list
        ) %>%
        dplyr::mutate(
            unit = unit,
            order1 = substr(Xterms, nchar(Xterms), nchar(Xterms)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == ')', 20, order1),
            order2 = substr(Xterms, 1, 2),
            order1 = as.integer(order1),
            Yterms = factor(Yterms,
                            levels = c('log(1/HOMA-IR)', 'log(ISI)',
                                       'log(IGI/IR)', 'log(ISSI-2)'),
                            ordered = TRUE)
        ) %>%
        dplyr::arrange(desc(order1), Yterms) %>%
        tidyr::separate(Xterms, into = c('fraction', 'Xterms'), sep = '_') %>%
        dplyr::mutate(fraction = renaming_fraction(fraction),
                      Xterms = factor(Xterms, unique(Xterms)))

    if (adj.p) {
        gee_results <- gee_results %>%
            dplyr::mutate(p.value = p.adjust(p.value, method = adj.p.method))
    }

    return(gee_results)
}

plot_gee_results <- function(results.gee) {
    data <- results.gee %>%
        mutate(Xterms = Xterms %>% factor(., unique(.)),
               fraction = factor(fraction, unique(fraction), ordered = TRUE))
    data %>%
        seer::trance('main_effect') %>%
        seer::visualize(groups = 'fraction~Yterms',
                        xlab = 'Percent difference with 95% CI in the outcomes\nfor each SD increase in fatty acid',
                        ylab = 'Fatty acids') %>%
        seer::vision_simple(10, legend.position = 'right') +
        ggplot2::theme(
            axis.text.y = element_text(face = ifelse(
                levels(data$Xterms) == 'Total (mole)', "bold", "plain"
            )),
            legend.key.width = grid::unit(0.75, "line"),
            legend.key.height = grid::unit(0.75, "line"),
            panel.margin = grid::unit(0.75, "lines")
        ) +
        ggplot2::scale_alpha_discrete(name = 'FDR-adjusted\np-value',
                                      range = c(1.0, 1.0)) +
        ggplot2::scale_size_discrete(name = 'FDR-adjusted\np-value',
                                     range = c(0.75, 4)) +
        ggplot2::facet_grid(fraction~Yterms, scale = 'free_y') +
        scale_color_grey(start = 0.75, end = 0, name = 'FDR-adjusted\np-value') +
        scale_x_continuous(breaks = c(-10, -5, 0, 5, 10))
}

renaming_outcomes <- function(x) {
    x %>%
        gsub('linvHOMA', 'log(1/HOMA-IR)', .) %>%
        gsub('lISI', 'log(ISI)', .) %>%
        gsub('lIGIIR', 'log(IGI/IR)', .) %>%
        gsub('lISSI2', 'log(ISSI-2)', .)
}

renaming_fa <- function(x) {
    x %>%
        gsub('(.*)(\\d\\d)(\\d)', '\\1_\\2:\\3', .) %>%
        gsub('n(\\d)$', 'n-\\1', .) %>%
        gsub('D(\\d\\d)$', 'D-\\1', .) %>%
        gsub('^pct_', '', .) %>%
        gsub('Total', 'Total (mole)', .)
}

renaming_fraction <- function(x) {
    x %>%
        gsub('CE|ce', 'Cholesteryl esters', .) %>%
        gsub('PL|pl', 'Phospholipids', .)
}

renaming_list <- function(x) {
    x %>%
        renaming_fa() %>%
        renaming_outcomes()
}
