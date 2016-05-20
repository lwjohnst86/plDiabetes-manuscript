## ----libraryData, message = FALSE, echo = FALSE--------------------------

source('../_opts/options.R')
source('../_opts/functions.R')
load('../../data/dysgly_data.RData')
ds.dysgly <- ds 
# %>% 
#     filter(ifelse(VN == 0 & DM == 1, 'drop', 'keep') == 'keep')
load('../../data/faDiabetes_data.RData')
ds <- ds %>% 
    rename(TotalNE = totalNE)

# ds <- ds %>% 

ds0yr <- filter(ds, VN == 0, !is.na(TotalNE))


## ----geeData, cache = TRUE-----------------------------------------------

gee.df <- ## Merge the generated datasets into one dataset.
    full_join(ds %>%
                select(-matches('pct_(pl|ce)|(pl|ce)'), -TotalNE,
                       -pl_Total, -ce_Total),
              ds %>%
                filter(VN == 0) %>%
                select(SID, TotalNE, pl_Total, ce_Total,
                       matches('pct_(pl|ce)|(pl|ce)')) %>%
                mutate_each(funs(as.numeric(scale(.))), -SID),
              by = 'SID') %>%
      group_by(VN) %>%
      mutate(Waist = as.numeric(scale(Waist)),
             ALT = as.numeric(scale(ALT)),
             BaseAge = as.numeric(scale(BaseAge)),
             FamHistDiab = as.factor(FamHistDiab)) %>%
      ungroup() %>%
      arrange(SID, VN)

outcomes <- c('linvHOMA', 'lISI', 'lIGIIR', 'lISSI2')
covariates <- c('Waist', 'VN', 'Sex', 'Ethnicity', 'TotalNE', 'ALT',
                'BaseAge', 'FamHistDiab')
plPool.pct <- c('pl_Total', grep('^pct_pl', names(ds), value = TRUE))
cePool.pct <- c('ce_Total', grep('^pct_ce', names(ds), value = TRUE))
plPool.conc <- grep('^pl\\d\\d', names(ds), value = TRUE)
cePool.conc <- grep('^ce\\d\\d', names(ds), value = TRUE)

fitPL.pct <- analyze_gee(gee.df, outcomes, plPool.pct, covariates, 'mol%')
fitCE.pct <- analyze_gee(gee.df, outcomes, cePool.pct, covariates, 'mol%')
fitPL.conc <- analyze_gee(gee.df, outcomes, plPool.conc, covariates, 'nmol/mL')
fitCE.conc <- analyze_gee(gee.df, outcomes, cePool.conc, covariates, 'nmol/mL')

results.gee <- rbind(fitPL.pct, fitCE.pct, fitPL.conc, fitCE.conc)

intPL.pct <- analyze_gee(gee.df, outcomes, plPool.pct, covariates, 'mol%', int = TRUE)
intCE.pct <- analyze_gee(gee.df, outcomes, cePool.pct, covariates, 'mol%', int = TRUE)
intPL.conc <- analyze_gee(gee.df, outcomes, plPool.conc, covariates, 'nmol/mL', int = TRUE)
intCE.conc <- analyze_gee(gee.df, outcomes, cePool.conc, covariates, 'nmol/mL', int = TRUE)

results.int.gee <- rbind(intPL.pct, intCE.pct, intPL.conc, intCE.conc)


## ----dysglyData----------------------------------------------------------

# ds <- ds %>% 
#     filter(ifelse(VN == 0 & DM == 1, 'drop', 'keep') == 'keep')

ds.logreg <-
    left_join(
    ds.dysgly %>%
        filter(VN == 0),
    ds.dysgly %>%
        filter(!is.na(TotalNE)) %>%
        select(SID, VN, IFG, IGT, DM) %>%
        mutate_each(funs(ifelse(is.na(.), 0, .)), IFG, IGT) %>%
        mutate(PreDM = as.numeric(rowSums(.[c('IFG', 'IGT')], na.rm = TRUE))) %>%
        select(SID, VN, DM, PreDM) %>%
        mutate(FactorDysgly = ifelse(PreDM == 1, 'PreDM',
                                     ifelse(DM == 1, 'DM',
                                            'NGT'))) %>%
        select(SID, VN, FactorDysgly) %>%
        spread(VN, FactorDysgly) %>%
        mutate(
            DetailedConvert = as.factor(paste0(`0`, '-', `1`, '-', `2`)),
            ConvertDysgly = as.numeric(!grepl('NGT$|NGT-NA$|NGT-NGT-NA$|NGT-NA-NA$',
                                              DetailedConvert)), 
            ConvertDM = as.numeric(grepl('-DM$|-DM-DM$|-DM-NA$|-DM-NGT$',
                                         DetailedConvert)),
            ConvertPreDM = as.numeric(grepl('-PreDM$|-PreDM-PreDM$|-PreDM-NA$',
                                            DetailedConvert))
            #Drop = grepl('-PreDM$|-PreDM-PreDM$|-PreDM-NA$', DetailedConvert)
            )
    ) %>%
    filter(!is.na(pl_Total)) %>%
    mutate_each(funs(as.numeric(scale(.))),
                matches('pct_(pl|ce)|(pl|ce)'),
                Waist, ALT, BaseAge, TotalNE)

## ----codeForInline, dependson = 'geeData'--------------------------------

pctChg.o <- ds %>%
  select(f.VN, HOMA, ISI, IGIIR, ISSI2) %>%
  gather(Measure, Value, -f.VN) %>%
  na.omit() %>%
  group_by(Measure, f.VN) %>%
  summarise(med = median(Value)) %>%
  spread(f.VN, med) %>%
  mutate(pctChg = abs(((`6-yrs` - `0-yr`) / `0-yr`) * 100) %>% 
             round(1)) %>%
  ungroup() %>% 
  select(pctChg) %>% 
  {paste0(min(.), '% to ', max(.), '%')}

cumPctPL <- calcCumPctFA(ds, 'pct_pl', 90)
cumPctCE <- calcCumPctFA(ds, 'pct_ce', 90)

corDS <- ds %>%
  filter(VN == 0) %>%
  select(pl_Total, matches('^pct_pl\\d'),
         ce_Total, matches('^pct_ce\\d'),
         Waist, Age, ALT, NEFA = TotalNE, `1/HOMA-IR` = invHOMA, ISI,
         `IGI/IR` = IGIIR, `ISSI-2` = ISSI2) %>%
  do(data.frame(cor(.[, grep('pl|ce', names(.), invert = TRUE)],
                    .[, grep('pl|ce', names(.))], use = "complete.obs", 
                    method = "spearman")) %>% add_rownames()) %>%
  melt() %>%
  splitSpeciesPools(., 'variable') %>%
  mutate(rowname = factor(rowname, levels = unique(rowname)),
         value = round(value, 2)) %>%
  select(-variable) %>%
  tbl_df()

extCorr <- function(data, fa.pool, fa.pool.abbrev){
    data %>%
      filter(value >= 0.3 | value <= -0.3, pool == fa.pool) %>%
      mutate(species = gsub('Total \\(mole\\)', paste('total', fa.pool.abbrev), species),
             basic.char = gsub('Chol', 'cholesterol', rowname) %>%
               gsub('Waist', 'waist', .)) %>%
      mutate(cor.fa = paste0(species, ' (r~s~=', round(value, 2), ')'),
             cor.basic = paste0(basic.char, ' (r~s~=', round(value, 2), ')')) %>%
      arrange(desc(species)) %>%
      select(cor.fa, cor.basic)
}

corPL <- extCorr(corDS, 'Phospholipids', 'PL')
corCE <- extCorr(corDS, 'Cholesteryl esters', 'CE')

extractBeta <- function(data) {
    data %>%
      mutate(estCI = paste0('(beta: ', trim(format(round(estimate, 1), nsmall = 1)),
                            ' CI: ',
                            trim(format(round(conf.low, 1), nsmall = 1)),
                            ', ',
                            trim(format(round(conf.high, 1), nsmall = 1)),
                            ')'),
             #indep.estCI = paste(indep, estCI),
             #dep = paste(dep %>% gsub('log\\((.*)\\)', '\\1', .), estCI),
             dep = paste(Yterms %>% gsub('log\\((.*)\\)', '\\1', .))
             ) %>%
      filter(p.value <= 0.05) %>%
      arrange(Yterms, desc(estimate)) %>%
      select(Yterms, Xterms, estCI)
}

plgee <- extractBeta(fitPL.pct)
 ##extractBeta(fitPL.conc)
cegee <- extractBeta(fitCE.pct)
 ##extractBeta(fitCE.conc)

# Conversion to dysglycemia
dm.convert <- table(ds.logreg$ConvertDM)[2] %>% 
    {paste0(., ' (', round((. / 477)*100, 0), '%)')}
predm.convert <- table(ds.logreg$ConvertPreDM)[2] %>% 
    {paste0(., ' (', round((. / 477)*100, 0), '%)')}
#table(ds.logreg$ConvertDysgly)[2]

## ----inlineCodeDisc------------------------------------------------------

betagee <- rbind(fitPL.pct %>% mutate(unit = 'pct'),
                 fitPL.conc %>% mutate(unit = 'conc'),
                 fitCE.pct %>% mutate(unit = 'conc'),
                 fitCE.conc %>% mutate(unit = 'conc')) %>%
  filter(p.value < 0.05) %>%
  mutate(est = abs(round(estimate, 2)),
         dep = ifelse(as.numeric(Yterms) %in% 1:2, 'IR', 'BCD') %>% factor,
         frac = ifelse(grepl('Phospholipids', fraction), 'PL', 'CE') %>% factor) %>%
  select(indep = Xterms, dep, unit, frac, est)

exB <- function(stat, fa, y = c('IR', 'BCD'), ## exB = extract Beta
                 pool = c('PL', 'CE'), u = c('pct'),
                 data=betagee) {
    data %>%
      filter(indep %in% fa, dep %in% y, frac %in% pool,
             unit %in% u) %>%
      select(est) %>%
      .[[1]] %>%
      as.numeric() %>%
      stat
}


## ----tabBasicChar, results = 'asis'--------------------------------------

tab1SampleRow <- tableRowSampleSize(ds0yr, 'pl_Total', 'f.VN')

tab1MeanSDRows <- ds0yr %>%
  select(f.VN, Age, BMI, Waist, TAG, Chol, LDL, HDL, ALT, TotalNE,
         pl_Total, ce_Total) %>%
  tableRowDescriptiveStat(., meanSD)

tab1MedianIQRRows <- ds0yr %>%
  select(f.VN, HOMA, ISI, IGIIR, ISSI2) %>%
  tableRowDescriptiveStat(., medianIQR)

tab1NPctRows <- ds0yr %>%
  select(f.VN, Sex, Ethnicity, IFG, IGT, FamHistDiab) %>%
  mutate(FamHistDiab =
           plyr::mapvalues(FamHistDiab, c('No', 'Yes'),
                           c('No family history', 'Family history of diabetes'))) %>%
  tableRowNPercent(., 'f.VN') %>%
  ungroup() %>% 
  filter(Measure != '?', !is.na(Measure), Measure != 0, Measure != 'M',
         Measure != 'No family history') %>%
  mutate(Measure = ifelse(Measure == 1, paste0(Demographic), Measure)) %>%
  select(-Demographic)

set.caption(tabRef('tabBasicChar', 'Baseline characteristics of PROMISE participants at the baseline
examination (2004-2006).'))
rbind(tab1SampleRow, tab1NPctRows, tab1MeanSDRows, tab1MedianIQRRows) %>%
  mutate(Measure = Measure %>%
           factor(., levels = unique(.), labels =
                    c('n',
                      'Female', 'European', 'Latino/a',
                      'Other ethnicity', 'South Asian', 'IFG', 'IGT',
                      'Family members with diabetes',
                      'Age (yrs)', 'BMI (kg/m^2^)', 'WC (cm)',
                      'TAG (mmol/L)', 'Chol (mmol/L)', 'LDL (mmol/L)',
                      'HDL (mmol/L)', 'ALT (U/L)', 'NEFA (nmol/mL)', 'FA in PL (nmol/mL)',
                      'FA in CE (nmol/mL)',
                      'HOMA-IR', 'ISI', 'IGI/IR', 'ISSI-2'
                      )
                  )) %>%
  rename(`Baseline visit` = `0-yr`) %>%
  pander()


## ----figDistFA, cache = TRUE, fig.cap = figRef('figDistFA', 'Distribution of fatty acids in the phospholipid and cholesteryl ester fractions as a percentage of the total pool.  Each dot represents a data point and each box shows the median and interquartile range.')----

ds0yr %>%
  select(matches('^pct_')) %>%
  gather(FA, conc) %>%
  splitSpeciesPools(., 'FA') %>%
    # spread(species, conc) %>% 
    # seer::trance('boxes_dots') %>% 
    # seer::visualize(dots = FALSE)
    ggplot(aes(species, conc)) +
    geom_boxplot(outlier.shape = NA) + coord_flip() +
  #plotBoxWithJitter(., 'species', 'conc') +
  facet_wrap( ~ pool, ncol = 2) +
  ylab('Proportion (mol%)') +
  xlab('Fatty acids') +
  theme_tufte(10, base_family = 'Arial') +
  theme(strip.background = element_rect(fill = 'grey95', colour = 'grey95'))


## ----figHeatmap, cache = TRUE, fig.cap = figRef('figHeatmap', 'Spearman correlation coefficient heatmap of the phospholipid and cholesteryl ester fatty acids with covariates and outcomes of the PROMISE participants at the baseline visit (2004-2006).')----

corDS %>%
  plotHeatmap(., x = 'rowname', y = 'species',
              x.axis.label = 'Covariates and outcomes',
              y.axis.label = 'Fatty acids (mol%)',
              corr.text.size = 3.5) +
  facet_grid(pool ~ .) +
  theme_tufte(10, base_family = 'Arial') +
  theme(legend.key.width = unit(0.7, "line"),
        strip.background = element_rect(fill = 'grey95', colour = 'grey95'),
        legend.key.height = unit(0.7, "line"),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm")
        )


## ----figForest-pct, cache = TRUE, dependson = 'geeData', fig.cap = figRef('geeFP_pct', 'Longitudinal associations of individual phospholipid and cholesteryl ester fatty acids (mol%) with insulin sensitivity and beta-cell function using generalized estimating equations over the 6 years in the PROMISE cohort.  Adjusted for time, sex, ethnicity, baseline age, waist circumference, family history of diabetes, total non-esterified fatty acids, and alanine aminotransferase.  Outcome variables were log-transformed, predictor variables were scaled, and x-axis values were exponentiated to represent percent difference per SD increase in the fatty acid.  P-values were adjusted for the Benjamini-Hochberg False Discovery Rate.')----

#plotForestGEE(rbind(fitPL.pct, fitCE.pct), ylab = 'Fatty acids (mol%)')
results.gee %>% 
    filter(unit == 'mol%') %>% 
    plot_gee_results() +
    ylab('Fatty acids (mol%)')


## ----figManhattan-pct, cache = TRUE, dependson = c('geeData', 'interactData'), fig.cap = figRef('geeInt_pct', 'P-values from a test for a time interaction by phospholipid and cholesteryl ester fatty acids (mol%) in the fully adjusted GEE model, including the interaction term.  P-values were adjusted for the Benjamini-Hochberg False Discovery Rate.')----
results.int.gee %>% 
    filter(unit == 'mol%') %>% 
    plotManhattanGEE(ylab = 'Fatty acids (mol%)')

## ----figInt-pl, fig.height = 5, fig.width = 8, dependson = 'geeData', cache = TRUE, fig.cap = figRef('intPlots_pl_pct', 'Calculated least-squares means (with 95% CI) of the outcomes from fully adjusted GEE models based on phospholipid fatty acids (mol%), split by median for the groups, that had a significant interaction by time.')----

multiPlot(plotInteract(gee.df, 'lISI', 'pct_pl202n6', 2, title = 'A'),
          plotInteract(gee.df, 'lISSI2', 'pct_pl202n6', 2, title = 'B',
                       y.axis.label = 'log(ISSI-2)'),
          plotInteract(gee.df, 'lISSI2', 'pct_pl160', 2, title = 'C',
                       y.axis.label = 'log(ISSI-2)'),
          cols = 3)


## ----tabQIC, results = 'asis'--------------------------------------------

subds <- ds %>%
  select(SID, VN, lISSI2, pl_Total, Sex, Ethnicity, MET, ALT,
         BaseAge, Waist, BMI, AlcoholPerWk, TobaccoUse, SelfEdu,
         Occupation, TotalNE, FamHistDiab) %>%
  mutate(SelfEdu = as.factor(SelfEdu),
         Occupation = as.factor(Occupation),
         TobaccoUse = as.factor(TobaccoUse),
         AlcoholPerWk = as.factor(AlcoholPerWk)) %>%
  filter(complete.cases(.))

levels(subds$SelfEdu)[1:2] <- c('2', '2')
levels(subds$SelfEdu)[5] <- '2'
subds$SelfEdu <- droplevels(subds$SelfEdu)
#table(subds$SelfEdu, useNA = 'ifany')

M0 <- geeglm(lISSI2 ~ pl_Total + VN, data = subds,
                id = SID, family = gaussian,
                corstr = 'exchangeable')
M1 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge)
M2 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + Waist)
#M3 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + BMI)
M4 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + Waist + MET +
               AlcoholPerWk + TobaccoUse + SelfEdu + Occupation +
               FamHistDiab)
M5 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + MET + Waist +
               AlcoholPerWk + TobaccoUse + SelfEdu + Occupation +
               FamHistDiab + TotalNE + ALT)
M6 <- update(M0, . ~ . + Sex + Ethnicity + Waist + TotalNE)
M7 <- update(M0, . ~ . + Sex + Ethnicity + Waist + TotalNE + ALT)
M8 <- update(M0, . ~ . + Sex + Ethnicity + Waist + TotalNE + ALT + BaseAge)
M9 <- update(M0, . ~ . + Sex + Ethnicity + Waist + TotalNE + ALT + BaseAge + FamHistDiab)

set.caption(supTabRef('tabQIC', 'Comparing generalized estimating equation
models adjusting for different covariates using
Quasi-Likelihood Information Criterion.'))
model.sel(M0, M1, M2, #M3, 
          M4, M5, M6, M7, M8, M9,
          rank = QIC) %>%
  tableModelSelection()


## ----figForest-conc, cache = TRUE, dependson = 'geeData', fig.cap = supFigRef('geeFP_conc', 'Longitudinal associations of individual phospholipid and cholesteryl ester fatty acids as a concentration (nmol/mL) with insulin sensitivity and beta-cell function using generalized estimating equations for the 6 years in the PROMISE cohort.  Adjusted for visit, sex, ethnicity, baseline age, waist circumference, family history of diabetes, total free fatty acids, and alanine aminotransferase.  Outcome variables were log-transformed and x-axis values were exponentiated to represent percent difference per SD increase in the fatty acids.  P-values were adjusted for the False Discovery Rate.')----
#plotForestGEE(rbind(fitPL.conc, fitCE.conc), ylab = 'Fatty acid\n(nmol/mL)')
results.gee %>% 
    filter(unit == 'nmol/mL') %>% 
    plot_gee_results() +
    ylab('Fatty acid (nmol/mL)')

## ----figManhattan-conc, cache = TRUE, dependson = c('geeData', 'interactData'), fig.cap = supFigRef('geeInt_conc', 'P-values from a test for a time interaction by phospholipid and cholesteryl ester fatty acids (as a concentration, nmol/mL) in the fully adjusted GEE model, including the interaction term.  P-values were adjusted for the False Discovery Rate.')----

results.int.gee %>% 
    filter(unit == 'nmol/mL') %>% 
    plotManhattanGEE(ylab = 'Fatty acids (nmol/mL)')
#plotManhattanGEE(rbind(intPL.conc, intCE.conc), ylab = 'Fatty acids (nmol/mL)')


## ----tabDistriFA, results = 'asis'---------------------------------------

set.caption(supTabRef('supTabDistFA', "Mean concentration and
percent of each fatty acid with SD
in the phospholipid (PL) and cholesteryl ester (CE) fraction in PROMISE
participants without diabetes at the baseline visit only (2004-2006)."))
ds %>% filter(VN == 0) %>%
  select(matches('pl'), matches('ce')) %>%
  gather(FA, value) %>%
  group_by(FA) %>%
  summarize(meanSD = meanSD(value)) %>%
  mutate(pct = ifelse(grepl('^pct_', FA), '(mol%)', '(nmol/mL)')) %>%
  splitSpeciesPools('FA') %>%
  mutate(pool = ifelse(pool == 'Phospholipids', 'PL', 'CE') %>%
           paste(., pct) %>%
         factor(., levels = unique(.))) %>%
  select(-pct, -FA) %>%
  spread(pool, meanSD) %>%
  arrange(desc(species)) %>%
  rename('Fatty acid' = species) %>%
  pander(justify = c('left', rep('center', 4)),
         emphasize.strong.cells = which(.[, 1] == 'Total', arr.ind = TRUE),
         missing = '')


## ----tabEstCI-GEE, results = 'asis', dependson = 'geeData'---------------

set.caption(supTabRef('supTabGEE_pct',
'Exponentiated GEE beta coefficients and 95% confidence intervals for the
associations of phospholipid and cholesteryl ester fatty acids **(mol%)** with
the outcome variables. Calculated 95% CI were *not* adjusted for multiple
testing using the False Discovery Rate (FDR).  Astericks indicates p<0.05
adjusted by the FDR.'))
table_geeDF(rbind(fitPL.pct %>%
                    mutate(Xterms = paste('PL', Xterms) %>%
                             factor(., levels = unique(.)),
                           orderno = 1),
                fitCE.pct %>%
                  mutate(Xterms = paste('CE', Xterms) %>%
                           factor(., levels = unique(.)),
                         orderno = 2)))

set.caption(supTabRef('supTabGEE_conc',
'Exponentiated GEE beta coefficients and 95% confidence intervals for the
associations of phospholipid and cholesteryl ester fatty acids **(nmol/mL)** with
the outcome variables.  Calculated 95% CI were *not* adjusted for multiple
testing using the False Discovery Rate (FDR).  Astericks indicates p<0.05
adjusted by the FDR.'))
table_geeDF(rbind(fitPL.conc %>%
                    mutate(Xterms = paste('PL', Xterms) %>%
                             factor(., levels = unique(.)),
                           orderno = 1),
                fitCE.conc %>%
                  mutate(Xterms = paste('CE', Xterms) %>%
                           factor(., levels = unique(.)),
                         orderno = 2)))


## ----figInt-ce, fig.height = 9, fig.width = 8, dependson = 'geeData', cache = TRUE, fig.cap = supFigRef('intPlots_ce', 'Calculated least-squares means (with 95% CI) of the outcomes from fully adjusted GEE models based on cholesteryl ester fatty acids (mol% for A-F and nmol/mL for G-H), split by tertile for the groups, that had a significant interaction by time.')----

multiPlot(## percent
          ## 1
          plotInteract(gee.df, 'linvHOMA', 'pct_ce203n6', 3, title = 'A',
                       y.axis.label = 'log(1/HOMA-IR)'), 
          ## 4
          plotInteract(gee.df, 'linvHOMA', 'pct_ce183n6', 3, title = 'D',
                       y.axis.label = 'log(1/HOMA-IR)'),
          ## 7
          plotInteract(gee.df, 'lISI', 'ce183n6', 3, title = 'G'),
          ## 2
          plotInteract(gee.df, 'lIGIIR', 'pct_ce203n6', 3, title = 'B',
                       y.axis.label = 'log(IGI/IR)'),
          ## 5
          plotInteract(gee.df, 'lISI', 'pct_ce183n6', 3, title = 'E'),
          ## 8
          plotInteract(gee.df, 'lIGIIR', 'ce203n6', 3, title = 'H'),
          ## 3
          plotInteract(gee.df, 'lISSI2', 'pct_ce203n6', 3, title = 'C',
                       y.axis.label = 'log(ISSI-2)'), 
          ## 6
          plotInteract(gee.df, 'lISSI2', 'pct_ce160', 3, title = 'F',
                       y.axis.label = 'log(ISSI-2)'),
          cols = 3)


## ----figInt-pl-conc, fig.height = 5, fig.width = 8, dependson = 'geeData', cache = TRUE, fig.cap = supFigRef('intPlots_pl_conc', 'Calculated least-squares means (with 95% CI) of the outcomes from fully adjusted GEE models based on phospholipid fatty acids (nmol/mL), split by median for the groups, that had a significant interaction by time.')----

multiPlot(plotInteract(gee.df, 'lISI', 'pl202n6', 2, title = 'A'),
          plotInteract(gee.df, 'lISSI2', 'pl202n6', 2, title = 'B',
                       y.axis.label = 'log(ISSI-2)'),
          plotInteract(gee.df, 'lISI', 'pl183n6', 2, title = 'C'),
          cols = 3)


## ----figLogReg, fig.cap=supFigRef('dysglyLogReg', 'Odds ratios for the association of conversion to dysglycemia over the 6 years and fatty acids within the phospholipid and cholesteryl ester serum fractions (mol%). Models were adjusted using the covariates from the primary GEE model, not including time. P-values were adjusted for the Benjamini-Hochberg False Discovery Rate.'), cache=TRUE----

covariates <- c('Waist', 'Sex', 'Ethnicity', 'TotalNE', 'ALT',
                'BaseAge', 'FamHistDiab')
ds.names <- names(ds.logreg)
plPool.pct <- c('pl_Total', grep('^pct_pl', ds.names, value = TRUE))
cePool.pct <- c('ce_Total', grep('^pct_ce', ds.names, value = TRUE))
plPool.conc <- grep('^pl\\d\\d', ds.names, value = TRUE)
cePool.conc <- grep('^ce\\d\\d', ds.names, value = TRUE)

logistic_regression_loop <- function(data, y, x, covars) {
    data <- data %>%
        dplyr::select_(.dots = c(y, x, covars)) %>%
        tidyr::gather_('Yterms', 'Yterm', y) %>%
        tidyr::gather_('Xterms', 'Xterm', x)

    logreg_formula <- reformulate(c('Xterm', covars),
                                  response = 'Yterm')
    ds.plot <- data %>%
        dplyr::group_by(Yterms, Xterms) %>%
        dplyr::do(glm(logreg_formula, data = .,
                      family = binomial()) %>%
                      broom::tidy(., conf.int = TRUE,
                                  exponentiate = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(grepl('Xterm$', term)) %>%
        dplyr::mutate(p.value = p.adjust(p.value, 'BH')) %>%
        dplyr::mutate(Xterms = renameSpecies(Xterms) %>%
                          renamePool()) %>%
        tidyr::separate(Xterms, into = c('Fraction', 'Xterms'),
                        sep = '\\.') %>%
        dplyr::mutate(
            order1 = substr(Xterms, nchar(Xterms), nchar(Xterms)),
            order1 = ifelse(order1 == 0, 10, order1),
            order1 = ifelse(order1 == ')', 20, order1),
            order1 = as.integer(order1)
        ) %>%
        dplyr::arrange(desc(order1), Yterms) %>%
        dplyr::select(-order1) %>%
        dplyr::mutate(Xterms = Xterms %>% factor(., unique(.)),
               Fraction = factor(Fraction, rev(unique(Fraction)), ordered = TRUE)) 
    
    ds.plot %>%
        seer::trance('main_effect') %>%
        seer::visualize(groups = 'Fraction~.', center.line = '1',
                        ylab = 'Fatty acids (mol%)',
                        xlab = 'Odds ratio with 95% CI per\nSD increase in fatty acid') %>%
        seer::vision_simple(base.size = 10, legend.position = 'right') +
        ggplot2::theme(
            axis.text.y = element_text(face = ifelse(
                grepl('Total', levels(ds.plot$Xterms)), "bold", "plain"
            )),
            legend.key.width = grid::unit(0.75, "line"),
            legend.key.height = grid::unit(0.75, "line"),
            panel.margin = grid::unit(0.75, "lines")
        ) +
        ggplot2::scale_alpha_discrete(name = 'FDR-adjusted\np-value',
                                      range = c(1.0, 1.0)) +
        ggplot2::scale_size_discrete(name = 'FDR-adjusted\np-value',
                                     range = c(0.75, 4)) +
        scale_color_grey(start = 0.75, end = 0, name = 'FDR-adjusted\np-value')
        
}

logistic_regression_loop(ds.logreg, 'ConvertDysgly',
                         c(plPool.pct, cePool.pct),
                         covariates)

## ----testFigRefs---------------------------------------------------------

checkFigTabNumbering(figRef)
checkFigTabNumbering(supFigRef)
checkFigTabNumbering(tabRef)
checkFigTabNumbering(supTabRef)


## ----sessioninfo---------------------------------------------------------

 #sessionInfo()


