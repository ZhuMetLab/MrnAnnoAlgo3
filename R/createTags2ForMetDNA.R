
# createTags2ForMetDNA ------------------------------------------------------------------------

#' createTags2ForMetDNA
#'
#' @param ms1_data A tibble containing the MS1 data in MRN3 format with columns 'name', 'mz', 'rt', and 'ccs'.
#' @param result_df A table of MRN3 annotation result ("result_final.csv")
#'
#' @return Tags2 of MetDNA
#' @export
#'
#' @examples createTags2ForMetDNA(ms1_data = ms1_data, result_df = result_df)
#'
createTags2ForMetDNA <- function(ms1_data, result_df) {

    cat('MetDNA Tags2 creation from MRN3.', '\n')
    cat('\n')

    tags2_mrn3 = lapply(seq_len(nrow(ms1_data)), function(idx) {

        # get feature annotation info
        f_name <- ms1_data$name[idx]
        table_anno <- result_df[result_df$feature_name == f_name, ]
        # setup (tags2) list for MetDNA
        list_anno <- list()
        if (nrow(table_anno) > 0) {
            for (i in seq_len(nrow(table_anno))) {
                list_anno[[i]] <- list(
                    # tag(MRN3) == type(MetDNA): seed, adductAnnotation or metAnnotation
                    # type = ifelse(table_anno$tag[i] == 'seed', yes = 'seed', no = 'metAnnotation'),
                    type = table_anno$tag[i],
                    # feature from info
                    From = table_anno$id_from[i],
                    From.peak = table_anno$feature_from[i],
                    # metabolite annotation result
                    to = table_anno$id[i],
                    step = ifelse(table_anno$tag[i] == 'seed', yes = NA, no = 1), # step in MRN3 only 1
                    level = as.numeric(stringr::str_extract(table_anno$confidence[i], '\\d')),
                    as.seed = TRUE, # no redundancy in MRN3
                    as.seed.round = table_anno$seed_round[i],
                    isotope = "[M]", # no isotope annotation in MRN3
                    adduct = table_anno$adduct[i],
                    charge = 1, # single charge in MRN3 & MetDNA
                    formula = table_anno$formula[i],
                    # match error and score
                    mz.error = as.numeric(table_anno$mz_error[i]),
                    rt.error = as.numeric(table_anno$rt_error_rel[i]),
                    ccs.error = as.numeric(table_anno$ccs_error[i]),
                    int.error = NA,
                    ms2.sim = as.numeric(table_anno$ms2_score[i]),
                    nfrag = as.numeric(table_anno$matched_frag[i]),
                    score = unname(table_anno$total_score[i])
                )
            }
        }

        # setup MetDNA PeakInfo object
        data_PeakInfo = new(
            Class = "PeakInfo",
            name = ms1_data$name[idx] %>% as.character(),
            mz = ms1_data$mz[idx] %>% as.numeric() %>% round(x = ., digits = 4),
            rt = ms1_data$rt[idx] %>% as.numeric() %>% round(x = ., digits = 2),
            ccs = -1,
            ms2 = data.frame(),
            annotation = list_anno
        )
    })

    return(tags2_mrn3)

}

# library(tidyverse)
# library(MRN3)
# wd <- 'I:/Data/projects/MetDNA3_R/test_sample/NIST_plasma_hilic/POS/'
# ms1_data = getMs1Data4Mrn3(wd = wd)
# ms2_data = getMs2Data4Mrn3(wd = wd)
# # result_df = readr::read_csv(file.path(wd, '02_result_MRN3', 'result_mrn_anno.csv'))
# result_df = readr::read_csv(file.path(wd, '02_result_MRN3', 'result_final.csv'))
# # test
# f_name <- 'M72T316'
# idx <- which(ms1_data$name == f_name)
# tags2_after_annotation[[idx]]
# tags2_mrn3[[idx]]
# # load tags2
# load('S:/memberdata/zhanghaosong/04_test/metdna_mrn3/20240110_MetDNA3_param_test/NIST_plasma_hilic/POS/02_result_MRN_annotation/00_intermediate_data/tags2_after_redundancy_remove')
# tags2 <- lapply(tags2_after_annotation, function(x) {
#     temp_anno <- x@annotation
#     temp_anno <- data.frame(do.call(cbind, temp_anno)) %>% t() %>% as.data.frame()
# })
# tags2_len <- sapply(tags2, length)
# idx_tags2_keep <- which(tags2_len > 0)
# names(tags2) <- ms1_data$name
# tags2 <- tags2[idx_tags2_keep]
# table_tags2 <- lapply(seq_along(tags2), function(x) {
#     f_name <- names(tags2)[x]
#     temp_anno <- tags2[[x]]
#     temp_anno$feature_name <- f_name
#     return(temp_anno)
# }) %>% dplyr::bind_rows()
# table_tags2 <- tibble::as_tibble(table_tags2)
# colnames(table_tags2)
# table_tags2$step %>% unlist() %>% unique()

# MetDNA method ----
# in function: changeTags
#
# for(i in 1:length(kegg.id)){
#     score1 <- coefficients(lm.reg.mz)[1] +
#         coefficients(lm.reg.mz)[2] * mz.error[i]
#     score2 <- coefficients(lm.reg.rt)[1] +
#         coefficients(lm.reg.rt)[2] * rt.error[i]
#     score3 <- coefficients(lm.reg.dp)[1] +
#         coefficients(lm.reg.dp)[2] * ms2.sim[i]
#     score4 <- coefficients(lm.reg.ccs)[1] +
#         coefficients(lm.reg.ccs)[2] * ccs.error[i]
#
#     score <-
#         score1*weight_mz + score2*weight_rt  + score4*weight_ccs + score3*weight_dp
#
#     temp@annotation[[i]] <-
#         list(type = "seed",
#              From = NA,
#              From.peak = NA,
#              to = stringr::str_trim(kegg.id[i]),
#              step = NA,
#              level = 1,
#              as.seed = FALSE,
#              as.seed.round = NA,
#              isotope = "[M]",
#              adduct = stringr::str_trim(adduct[i]),
#              charge = 1,
#              # as.numeric(adduct_table$charge[match(adduct[i], adduct_table$name)]),
#              formula = stringr::str_trim(formula[i]),
#              mz.error = as.numeric(mz.error[i]),
#              rt.error = as.numeric(rt.error[i]),
#              ccs.error = as.numeric(ccs.error[i]),
#              int.error = NA,
#              ms2.sim = as.numeric(ms2.sim[i]),
#              nfrag = as.numeric(nfrag[i]),
#              score = unname(score)
#         )
# }
