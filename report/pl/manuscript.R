
source('../_opts/options.R')
source('../_opts/functions.R')
load('../../data/faDiabetes_data.RData')

ds <- ds %>%
  rename(pl_Total = totalPL, ce_Total = totalCE) %>%
  mutate(f.VN = factor(VN, c(0, 3, 6), c('0-yr', '3-yrs', '6-yrs')))

ds0yr <- filter(ds, VN == 0, !is.na(totalNE)) %>%
  mutate(FHDiab = as.numeric(FamHistDiab >= 1))


gee.df <- ## Merge the generated datasets into one dataset.
    full_join(ds %>%
                select(-matches('pct_(pl|ce)|(pl|ce)'), -totalNE,
                       -pl_Total, -ce_Total),
              ds %>%
                filter(VN == 0) %>%
                select(SID, totalNE, pl_Total, ce_Total,
                       matches('pct_(pl|ce)|(pl|ce)')) %>%
                mutate_each(funs(as.numeric(scale(.))), -SID),
              by = 'SID') %>%
      mutate(FamHistDiab =
               plyr::mapvalues(FamHistDiab, c(0, 1:12),
                               c('No', rep('Yes', 12)))) %>%
      group_by(VN) %>%
      mutate(Waist = as.numeric(scale(Waist)),
             ALT = as.numeric(scale(ALT)),
             BaseAge = as.numeric(scale(BaseAge))) %>%
      ungroup() %>%
      arrange(SID, VN)

outcomes <- c('linvHOMA', 'lISI', 'lIGIIR', 'lISSI2')
outcomes.rename <- c('log(1/HOMA-IR)', 'log(ISI)', 'log(IGI/IR)', 'log(ISSI-2)')
covariates <- c('Waist', 'VN', 'Sex', 'Ethnicity', 'totalNE', 'ALT',
                'BaseAge', 'FamHistDiab')
plPool.pct <- c('pl_Total', grep('^pct_pl', names(ds), value = TRUE))
cePool.pct <- c('ce_Total', grep('^pct_ce', names(ds), value = TRUE))
plPool.conc <- grep('^pl\\d\\d', names(ds), value = TRUE)
cePool.conc <- grep('^ce\\d\\d', names(ds), value = TRUE)

fitPL.pct <- geeDF(gee.df, plPool.pct, outcomes, covariates,
               renaming = outcomes.rename)
fitCE.pct <- geeDF(gee.df, cePool.pct, outcomes, covariates,
               renaming = outcomes.rename)
fitPL.conc <- geeDF(gee.df, plPool.conc, outcomes, covariates,
               renaming = outcomes.rename)
fitCE.conc <- geeDF(gee.df, cePool.conc, outcomes, covariates,
               renaming = outcomes.rename)

intPL.pct <- geeInteractDF(gee.df, outcomes, plPool.pct, covariates, interaction = 'VN',
                           renaming = outcomes.rename)
intPL.conc <- geeInteractDF(gee.df, outcomes, plPool.conc, covariates, interaction = 'VN',
                            renaming = outcomes.rename)
intCE.pct <- geeInteractDF(gee.df, outcomes, cePool.pct, covariates, interaction = 'VN',
                           renaming = outcomes.rename)
intCE.conc <- geeInteractDF(gee.df, outcomes, cePool.conc, covariates, interaction = 'VN',
                            renaming = outcomes.rename)


pctChg.o <- ds %>%
  select(f.VN, HOMA, ISI, IGIIR, ISSI2) %>%
  gather(Measure, Value, -f.VN) %>%
  na.omit() %>%
  group_by(Measure, f.VN) %>%
  summarise(med = median(Value)) %>%
  spread(f.VN, med) %>%
  mutate(pctChg = ((`6-yrs` - `0-yr`) / `0-yr`) * 100) %>%
  select(pctChg) %>%
  abs() %>%
  round(1) %>%
  {paste0(min(.), '% to ', max(.), '%')}

cumPctPL <- calcCumPctFA(ds, 'pct_pl', 90)
cumPctCE <- calcCumPctFA(ds, 'pct_ce', 90)

corDS <- ds %>%
  filter(VN == 0) %>%
  select(pl_Total, matches('^pct_pl\\d'),
         ce_Total, matches('^pct_ce\\d'),
         Waist, Age, ALT, NEFA = totalNE, `1/HOMA-IR` = invHOMA, ISI,
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
             dep = paste(dep %>% gsub('log\\((.*)\\)', '\\1', .))
             ) %>%
      filter(p.value <= 0.05) %>%
      arrange(dep, desc(estimate)) %>%
      select(dep, indep, estCI)
}

plgee <- extractBeta(fitPL.pct)
 ##extractBeta(fitPL.conc)
cegee <- extractBeta(fitCE.pct)
 ##extractBeta(fitCE.conc)


betagee <- rbind(fitPL.pct %>% mutate(unit = 'pct'),
                 fitPL.conc %>% mutate(unit = 'conc'),
                 fitCE.pct %>% mutate(unit = 'conc'),
                 fitCE.conc %>% mutate(unit = 'conc')) %>%
  filter(p.value < 0.05) %>%
  mutate(est = abs(round(estimate, 2)),
         dep = ifelse(as.numeric(dep) %in% 1:2, 'IR', 'BCD') %>% factor,
         frac = ifelse(grepl('Phospholipids', pool), 'PL', 'CE') %>% factor) %>%
  select(indep, dep, unit, frac, est)

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


tab1SampleRow <- tableRowSampleSize(ds0yr, 'pl_Total', 'f.VN')

tab1MeanSDRows <- ds0yr %>%
  select(f.VN, Age, BMI, Waist, TAG, Chol, LDL, HDL, ALT, totalNE,
         pl_Total, ce_Total) %>%
  tableRowDescriptiveStat(., meanSD)

tab1MedianIQRRows <- ds0yr %>%
  select(f.VN, HOMA, ISI, IGIIR, ISSI2) %>%
  tableRowDescriptiveStat(., medianIQR)

tab1NPctRows <- ds0yr %>%
  select(f.VN, Sex, Ethnicity, IFG, IGT, FamHistDiab) %>%
  mutate(FamHistDiab =
           plyr::mapvalues(FamHistDiab, c(0, 1:6),
                           c('No family history', rep('Family history of diabetes', 6)))) %>%
  tableRowNPercent(., 'f.VN') %>%
  filter(Measure != '?', !is.na(Measure), Measure != 0, Measure != 'M',
         Measure != 'No family history') %>%
  mutate(Measure = gsub('^F$', 'Female', Measure),
         Measure = ifelse(Measure == 1, paste0(Demographic), Measure)) %>%
  select(-Demographic)

set.caption(tabRef('tabBasicChar', 'Basic characteristics of PROMISE participants at the baseline
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


ds0yr %>%
  select(matches('^pct_')) %>%
  gather(FA, conc) %>%
  splitSpeciesPools(., 'FA') %>%
  plotBoxWithJitter(., 'species', 'conc') +
  facet_wrap( ~ pool, ncol = 2) +
  ylab('Proportion (mol%)') +
  xlab('Fatty acids') +
  theme_tufte(10) +
  theme(strip.background = element_rect(fill = 'grey90', colour = 'grey90'))


corDS %>%
  plotHeatmap(., x = 'rowname', y = 'species',
              x.axis.label = 'Covariates and outcomes',
              y.axis.label = 'Fatty acids (mol%)',
              corr.text.size = 3.5) +
  facet_grid(pool ~ .) +
  theme_tufte(10) +
  theme(legend.key.width = unit(0.7, "line"),
        strip.background = element_rect(fill = 'grey90', colour = 'grey90'),
        legend.key.height = unit(0.7, "line"),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm")
        )


plotForestGEE(rbind(fitPL.pct, fitCE.pct), ylab = 'Fatty acids (mol%)')


plotManhattanGEE(rbind(intPL.pct, intCE.pct), ylab = 'Fatty acids (mol%)')


multiPlot(plotInteract(gee.df, 'lISI', 'pct_pl202n6', 2, title = 'A'),
          plotInteract(gee.df, 'lISSI2', 'pct_pl202n6', 2, title = 'B',
                       y.axis.label = 'log(ISSI-2)'),
          plotInteract(gee.df, 'lISSI2', 'pct_pl160', 2, title = 'C',
                       y.axis.label = 'log(ISSI-2)'),
          cols = 3)


subds <- ds %>%
  select(SID, VN, lISSI2, pl_Total, Sex, Ethnicity, MET, ALT,
         BaseAge, Waist, BMI, AlcoholPerWk, TobaccoUse, SelfEdu,
         Occupation, totalNE, FamHistDiab) %>%
  mutate(SelfEdu = as.factor(SelfEdu),
         Occupation = as.factor(Occupation),
         TobaccoUse = as.factor(TobaccoUse),
         AlcoholPerWk = as.factor(AlcoholPerWk)) %>%
  filter(complete.cases(.))

levels(subds$SelfEdu)[1:2] <- c('2', '2')
subds$SelfEdu <- droplevels(subds$SelfEdu)

M0 <- geeglm(lISSI2 ~ pl_Total + VN, data = subds,
                id = SID, family = gaussian,
                corstr = 'exchangeable')
M1 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge)
M2 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + Waist)
M3 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + BMI)
M4 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + Waist + MET +
               AlcoholPerWk + TobaccoUse + SelfEdu + Occupation +
               FamHistDiab)
M5 <- update(M0, . ~ . + Sex + Ethnicity + BaseAge + MET + Waist +
               AlcoholPerWk + TobaccoUse + SelfEdu + Occupation +
               FamHistDiab + totalNE + ALT)
M6 <- update(M0, . ~ . + Sex + Ethnicity + Waist + totalNE)
M7 <- update(M0, . ~ . + Sex + Ethnicity + Waist + totalNE + ALT)
M8 <- update(M0, . ~ . + Sex + Ethnicity + Waist + totalNE + ALT + BaseAge)
M9 <- update(M0, . ~ . + Sex + Ethnicity + Waist + totalNE + ALT + BaseAge + FamHistDiab)

set.caption(supTabRef('tabQIC', 'Comparing generalized estimating equation
models adjusting for different covariates using
Quasi-Likelihood Information Criterion.'))
model.sel(M0, M1, M2, M3, M4, M5, M6, M7, M8, M9,
          rank = QIC) %>%
  tableModelSelection()


plotForestGEE(rbind(fitPL.conc, fitCE.conc), ylab = 'Fatty acid\n(nmol/mL)')


plotManhattanGEE(rbind(intPL.conc, intCE.conc), ylab = 'Fatty acids (nmol/mL)')


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


set.caption(supTabRef('supTabGEE_pct',
'Exponentiated GEE beta coefficients and 95% confidence intervals for the
associations of phospholipid and cholesteryl ester fatty acids **(mol%)** with
the outcome variables. Calculated 95% CI were *not* adjusted for multiple
testing using the False Discovery Rate (FDR).  Astericks indicates p<0.05
adjusted by the FDR.'))
table_geeDF(rbind(fitPL.pct %>%
                    mutate(indep = paste('PL', indep) %>%
                             factor(., levels = unique(.)),
                           orderno = 1),
                fitCE.pct %>%
                  mutate(indep = paste('CE', indep) %>%
                           factor(., levels = unique(.)),
                         orderno = 2)))

set.caption(supTabRef('supTabGEE_conc',
'Exponentiated GEE beta coefficients and 95% confidence intervals for the
associations of phospholipid and cholesteryl ester fatty acids **(nmol/mL)** with
the outcome variables.  Calculated 95% CI were *not* adjusted for multiple
testing using the False Discovery Rate (FDR).  Astericks indicates p<0.05
adjusted by the FDR.'))
table_geeDF(rbind(fitPL.conc %>%
                    mutate(indep = paste('PL', indep) %>%
                             factor(., levels = unique(.)),
                           orderno = 1),
                fitCE.conc %>%
                  mutate(indep = paste('CE', indep) %>%
                           factor(., levels = unique(.)),
                         orderno = 2)))


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


multiPlot(plotInteract(gee.df, 'lISI', 'pl202n6', 2, title = 'A'),
          plotInteract(gee.df, 'lISSI2', 'pl202n6', 2, title = 'B',
                       y.axis.label = 'log(ISSI-2)'),
          plotInteract(gee.df, 'lISI', 'pl183n6', 2, title = 'C'),
          cols = 3)


checkFigTabNumbering(figRef)
checkFigTabNumbering(supFigRef)
checkFigTabNumbering(tabRef)
checkFigTabNumbering(supTabRef)


 #sessionInfo()

