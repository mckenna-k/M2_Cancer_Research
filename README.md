# M2_Cancer_Research

## Les vignettes test 

Lancer en bash les commandes suivantes
avant de faire un push sur le git.

Il ne doit pas y avoir d'erreur.

```
# test
echo 'rmarkdown::render("stress_test.Rmd")' | Rscript -

# launch msm and mstate models
echo 'rmarkdown::render("data_msm_updated.Rmd")' | Rscript -

# screen censor_param
echo 'rmarkdown::render("censor_effect.Rmd")' | Rscript -

# adding significant covariate
echo 'rmarkdown::render("censor_effect_covar.Rmd")' | Rscript -

# all models on BRCA cancer
echo 'rmarkdown::render("full_kc_stage_study.Rmd")' | Rscript -


```
