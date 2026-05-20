########### Packages Needed   ###############
library(data.table)
##########

##Note: This code will only provide the models and results needed to later make the plots
######### Data Set Up #############
shandiv_df_all <- read.csv("figure3.csv") %>% as.data.table

dt <- shandiv_df_all; rm(shandiv_df_all)

## re-encode bfdur_rec
dt.bf.dur <- rbind(
  dt[
    i = month.alias == "M1",
    j = .(month.alias = month.alias, bfdur_rec = bfdur_rec_m1_pct,participant_id)
  ],
  dt[
    i = month.alias == "M6",
    j = .(month.alias = month.alias, bfdur_rec = bfdur_rec_m6_pct,participant_id)
  ]
)
dt <- merge(dt,dt.bf.dur,all.x = T,sort = F)
## re-factor highest_education_mother
dt[,highest_education_mother := factor(
  highest_education_mother,
  levels = c("Post graduate degree", 
             "<HS or HS", 
             "Some college", 
             "Bachelors")
)]
## id variables for melting
vars <- c('month.alias','shannon','ga_best','rec_pca_weeks','feeding_method')
## factors for melting
measures.factor <- c(
  'infant_sex',
  'mother_ethrace_2',
  'medicaid_in_preg',
  'exposed_bm',
  'highest_education_mother',
  'feeding_method'
)
## numeric for melting
measures.numeric <- c(
  'pre.bmi',
  'bfdur_rec'
)

# Models 
dt.melt.list <- lapply(list(factors = measures.factor,numerics = measures.numeric),function(i){
  data.table::melt(
    dt[,j = .SD,.SDcols = c(vars,i)],
    id.vars = vars,
    value.factor = TRUE
  )[!is.na(value)]
})
## for linear models; results (list) stored in data.table
models.lm <- data.table::rbindlist(lapply(dt.melt.list,function(dt.melt){
  dt.melt[
    ,
    j = .(model.lm = {
      f <- paste("shannon","~","value")
      f <- ifelse(
        .BY$month.alias=="Birth",
        paste(f,"+ ga_best"),
        paste(f,"+ rec_pca_weeks + feeding_method")
      )
      if(.BY$variable %in% c("feeding_method", "bfdur_rec")){
        f <- "shannon ~ value + rec_pca_weeks"
      }
      list(lm(formula = as.formula(f), data = .SD))
    }),
    keyby = .(variable,month.alias)
  ]
}))

## coefficients; results (as.data.table) stored in data.table
models.lm <- merge(
  models.lm,
  models.lm[
    ,
    j = .(model.coefficients = {
      res <- summary(model.lm[[1]])$coefficients
      dt <- as.data.table(res)
      dt[,measure.type := rownames(res)]
      dt[,p.val.alias := {
        p_label_text <- as.numeric(sprintf("%.4f",round(`Pr(>|t|)`,4)))
        p_label_text <- gsub("\\.", "\u00B7", p_label_text)
        ifelse(p_label_text == "0","p < 0\u00B70001",paste("p =",p_label_text))
      }]
      list(dt)
    }),
    by=.(variable,month.alias)
  ]
)

## lsmeans; results (list) stored in data.table
models.lm <- merge(
  models.lm,
  models.lm[
    i = variable %in% c("highest_education_mother","feeding_method"),
    j = .(lsmeans = list(lsmeans::lsmeans(model.lm[[1]],~value))),
    by=.(variable,month.alias)
  ],
  all.x = TRUE
)
data.table::set(
  x = models.lm,
  i = models.lm[,.I[sapply(lsmeans,is.null)]],
  j = 'lsmeans',
  value = NA
)

## contrasts; results (as.data.table) stored in data.table
models.lm <- merge(
  models.lm,
  models.lm[
    i = !is.na(lsmeans),
    j = .(contrasts = {
      res <- data.table::as.data.table(
        lsmeans::contrast(lsmeans[[1]],"trt.vs.ctrl", adjust = "BH")
      )
      res[,p.val.alias := {
        p_label_text <- as.numeric(sprintf("%.4f",round(p.value,4)))
        p_label_text <- gsub("\\.", "\u00B7", p_label_text)
        ifelse(p_label_text == "0","p < 0\u00B70001",paste("p =",p_label_text))
      }]
      list(res)
      # list(data.table::as.data.table(
      # lsmeans::contrast(lsmeans[[1]],"trt.vs.ctrl", adjust = "BH")))
    }),
    by=.(variable,month.alias)
  ],
  all.x = TRUE
)
data.table::set(
  x = models.lm,
  i = models.lm[,.I[sapply(contrasts,is.null)]],
  j = 'contrasts',
  value = NA
)
