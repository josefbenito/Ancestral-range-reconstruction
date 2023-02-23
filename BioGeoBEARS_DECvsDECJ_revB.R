library(optimx)
library(GenSA)
library(FD) 
library(snow)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
setwd("~/UAH_Docs/Fall_2021/Cave_beetles_proposal/BioGeoBEARS")
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

trfn = "Beast_DivTime_NoOutgroup_M0.0012SD0.059_S500_25M_combined_summary_B10.newick"
moref(trfn)
pdffn = "Cave_beetle_tree.pdf"
pdf(file=pdffn, width=12, height=9)
tr = read.tree(trfn)
tr
plot(tr)
title("Cave beetle Abyss Iden70 Cov75 75p")
axisPhylo()
dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

geogfn = "Cave_beetle_geog_no_outgroup_karst_subregion.data"
moref(geogfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 4
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)
numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)
numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50 
BioGeoBEARS_run_object$speedup = TRUE 
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
runslow = TRUE
resfn = "Cave_beetle_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDEC = res
} else {
  load(resfn)
  resDEC = res
}

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE 
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE 
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resfn = "Cave_beetle_DEC+J_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECj = res
} else {
  load(resfn)
  resDECj = res
}

pdffn = "Cave_beetle_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=12, height=9)

analysis_titletxt ="BioGeoBEARS DEC on Cave beetle M0_unconstrained"
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.01, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.01, tipcex=0.8, statecex=0.4, splitcex=0.4, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="BioGeoBEARS DEC+J on Cave beetle M0_unconstrained"
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.01, tipcex=0.8, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.01, tipcex=0.8, statecex=0.4, splitcex=0.4, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

restable = NULL
teststable = NULL
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)
teststable$alt = c("DEC+J")
teststable$null = c("DEC")
row.names(restable) = c("DEC", "DEC+J")
restable = put_jcol_after_ecol(restable)
teststable
write.table(conditional_format_table(teststable), file="DEC vs DEC+J teststable.txt", quote=FALSE, sep="\t")
restable2 = restable
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike
write.table(conditional_format_table(restable_AIC_rellike), file="DEC vs DEC+J restable.txt", quote=FALSE, sep="\t")

model_name = "DEC"
res = resDEC
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  load(BSM_inputs_fn)
}
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))
runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  load(file="RES_clado_events_tables.Rdata")
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} 
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)
include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 4
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

pdffn = paste0(model_name, "_Cave_beetle_stochastic_map.pdf")
pdf(file=pdffn, width=12, height=9)
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC Cave beetle stochastic map", addl_params=list("j"), label.offset=0.05, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=2, lty=par("lty"), root.edge=TRUE, stratified=stratified)
plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC Cave beetle stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

model_name = "DEC+J"
res = resDECj
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  load(BSM_inputs_fn)
}
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))
runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  load(file="RES_clado_events_tables.Rdata")
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
}
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)
include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 4
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

pdffn = paste0(model_name, "_Cave_beetle_stochastic_map.pdf")
pdf(file=pdffn, width=12, height=9)
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC+J Cave beetle stochastic map", addl_params=list("j"), label.offset=0.01, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=2, lty=par("lty"), root.edge=TRUE, stratified=stratified)
plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="DEC+J cave beetle stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)