
# rmRedunInMetdna3 ----------------------------------------------------------------------------

#' rmRedunInMetdna3
#'
#' @param result_df result_mrn_anno table (feature-metabolite pairs)
#' @param num_candi number of candidates (topN) for redundancy removal
#' @param is_known_prioritization priority: level3.1 (knowns) > level3.2 (unknowns) for redundancy removal
#' @param wd_output file.path(wd, '02_result_MRN3')
#'
#' @return result_mrn_anno table with redundancy removal
#' @export
#'
rmRedunInMetdna3 <- function(
        result_df,
        num_candi = 3, # number of candidates
        is_known_prioritization = T, # priority: level3.1 (knowns) > level3.2 (unknowns)
        wd_output
        # info_mrn,
        # md_mrn
    ) {

    cat('Redundancy removal by MetDNA3 (MrnAnnoAlgo3):', '\n')
    cat('\n')

    # remove ID metAnnotations which have seed annotations -----------------------------------------
    ids_seed <- result_df$id[which(result_df$tag == 'seed')] %>% unique() %>% sort()
    idx_metanno <- which(result_df$tag == 'metAnnotation')
    idx_metanno_idx_rm <- which(result_df$id[idx_metanno] %in% ids_seed)
    if (length(idx_metanno[idx_metanno_idx_rm]) > 0) {
        result_df_rm_seedMetAnno <- result_df[-idx_metanno[idx_metanno_idx_rm], ]
        cat('Remove ID metAnnotations which have seed annotations:', length(idx_metanno[idx_metanno_idx_rm]), '\n')
        cat('\n')
    } else {
        result_df_rm_seedMetAnno <- result_df
    }

    # remove redundancy: 1 to N --------------------------------------------------------------------
    result_df_rm_seedMetAnno <- result_df_rm_seedMetAnno %>% arrange(feature_mz, feature_rt)
    fs <- result_df_rm_seedMetAnno$feature_name %>% unique()
    check_seed_filter <- c()
    check_known_filter <- c()
    check_candi_filter <- c()
    list_res_anno <- lapply(fs, function(temp_f) {
        temp_res <- result_df_rm_seedMetAnno %>% filter(feature_name == temp_f) %>% arrange(desc(total_score), id)
        # keep seed first (same as grade1)
        if ('seed' %in% temp_res$tag) {
            temp_res <- temp_res %>% filter(tag == 'seed')
            check_seed_filter <<- append(check_seed_filter, T)
        } else {
            check_seed_filter <<- append(check_seed_filter, F)
        }
        # keep level 3.1 first (grade 3.1)
        if (is_known_prioritization & any(startsWith(temp_res$id, 'MRN'))) {
            temp_res <- temp_res[startsWith(temp_res$id, 'MRN'), ]
            check_known_filter <<- append(check_known_filter, T)
        } else {
            check_known_filter <<- append(check_known_filter, F)
        }
        # keep top N candidates (manually)
        if (nrow(temp_res) > num_candi) {
            temp_res <- temp_res[1:num_candi, ]
            check_candi_filter <<- append(check_candi_filter, T)
        } else {
            check_candi_filter <<- append(check_candi_filter, F)
        }
        return(temp_res)
    })
    names(list_res_anno) <- fs
    result_df_rm_1toN <- dplyr::bind_rows(list_res_anno)
    result_df_rm_1toN <- result_df_rm_1toN %>% arrange(feature_mz, feature_rt)

    stat_df <- tibble::tibble(
        feature_name = fs,
        filtered_by_seed = check_seed_filter,
        filtered_by_known = check_known_filter,
        filtered_by_topN = check_candi_filter
    )
    readr::write_csv(stat_df, file.path(wd_output, 'stat_rm_redun.csv'), na = '')

    # only confidence assignment (same with MetDNA)
    result_df_confi_assign <- MrnAnnoAlgo3::rmRedunInMrn3(
        result_df = result_df_rm_1toN,
        is_only_confidence_assignment = T
    )

    ### calculate redundancy (same with MetDNA)
    ids_anno <- sort(unique(result_df_confi_assign$id))
    peaks_anno <- sort(unique(result_df_confi_assign$feature_name))
    redun <- calculateRedundancyInMrn3(result_df_confi_assign)
    cat('Annotated peaks:', length(peaks_anno), '\n')
    cat('Annotated metabolites:', length(ids_anno), '\n')
    cat("Peak redundancy:", round(redun[1], 2), '\n')
    cat("Metabolite redundancy:", round(redun[2], 2), '\n')
    cat('\n')

    return(result_df_confi_assign)


    # # No application !
    # # remove redundancy: N to 1 --------------------------------------------------------------------
    # # reconstruct RT prediction model
    # table_mrn_rt <- result_df_rm_1toN %>% select(id, feature_rt) %>% arrange(id)
    # table_mrn_rt <- table_mrn_rt %>%
    #     dplyr::group_by(id) %>%
    #     dplyr::summarize(feature_rt = median(feature_rt, na.rm = T))
    # table_mrn_rt$lib_pred_rt <- match(table_mrn_rt$id, info_mrn$id) %>% info_mrn$rt[.]
    # idx_na <- which(is.na(table_mrn_rt$lib_pred_rt))
    # if (length(idx_na) > 0) table_mrn_rt <- table_mrn_rt[-idx_na, ]
    # plot(table_mrn_rt$lib_pred_rt, table_mrn_rt$feature_rt)
    # cat("There are ", nrow(table_mrn_rt), " metabolites are used for RT prediction.\n", sep = "")
    #
    # ### filter ID MRN
    # md <- md_mrn[match(table_mrn_rt$id, rownames(md_mrn)), ]
    # ### remove NA which appear in more than 50% metabolites
    # md1 <- md
    # remove.idx1 <- which(apply(md, 2, function(x) {sum(is.na(x)/nrow(md))}) > 0.5)
    # if (length(remove.idx1) > 0) md1 <- md[, -remove.idx1]
    # ### impute NA
    # md2 <- t(impute::impute.knn(data = t(md1))[[1]])
    # ### remove MD which are same in all metabolites
    # md3 <- md2
    # remove.idx2 <- which(apply(md2, 2, sd) == 0)
    # if (length(remove.idx2) > 0) md3 <- md2[, -remove.idx2]
    #
    # # construct RF model
    # train.y <- table_mrn_rt$feature_rt
    # train.x <- md3
    #
    # rm(list = c('md', 'md1', 'md2', 'md3')); gc()
    #
    # cat('Tuning RF function using caret approach.\n')
    # train.data <- data.frame(RT = train.y, train.x, stringsAsFactors = FALSE)
    # para_control <- caret::trainControl(method = 'cv', number = 10, search = 'random', verbose = F)
    # set.seed(100)
    # rf.reg <- caret::train(
    #     RT ~ ., data = train.data,
    #     method = 'rf',
    #     metric = 'Rsquared',
    #     tuneLength = 10,
    #     trControl = para_control,
    #     importance = T,
    #     allowParallel = T
    # )
    #
    # # predict RT
    # test.x <- md_mrn[match(table_mrn_rt$id, rownames(md_mrn)), ]
    # temp.rt <- predict(object = rf.reg, newdata = test.x)
    # table_mrn_rt$lib_pred_rt_2 <- temp.rt
    # plot(table_mrn_rt$lib_pred_rt_2, table_mrn_rt$feature_rt)

}
