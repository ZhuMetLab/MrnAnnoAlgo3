
# rmRedunInMetdna3 ----------------------------------------------------------------------------

#' rmRedunInMetdna3
#'
#' @param result_df result_mrn_anno table (feature-metabolite pairs)
#' @param wd_output file.path(wd, '02_result_MRN_annotation')
#' @param num_candi number of candidates (topN) for redundancy removal
#' @param is_known_prioritization priority: level3.1 (knowns) > level3.2 (unknowns) for redundancy removal
#' @param is_clean_Nto1 re-regression by seeds & level3 results
#' @param md_mrn md_mrn for RT re-regression
#' @param column column for RT re-regression
#'
#' @return result_mrn_anno table with redundancy removal
#' @export
#'
rmRedunInMetdna3 <- function(
        result_df,
        wd_output,
        num_candi = 3, # number of candidates
        is_known_prioritization = T, # priority: level3.1 (knowns) > level3.2 (unknowns)
        is_clean_Nto1 = T, # re-regression by seeds & level3 results
        md_mrn = NULL, # for RT re-regression
        column = NULL # for RT re-regression
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

    # refine metabolite name with consensus name ---------------------------------------------------
    tbl_cons_name <- readr::read_csv(system.file('mrn_consensus_name.csv', package = 'MrnAnnoAlgo3'))
    idx_in_tbl <- which(result_df_rm_seedMetAnno$id %in% tbl_cons_name$id_mrn)
    if (length(idx_in_tbl) > 0) {
        cons_name_in_tbl <- match(result_df_rm_seedMetAnno$id[idx_in_tbl], tbl_cons_name$id_mrn) %>%
            tbl_cons_name$consensus_name[.]
        result_df_rm_seedMetAnno$name[idx_in_tbl] <- cons_name_in_tbl
    }

    # remove redundancy: 1 to N --------------------------------------------------------------------
    result_df_rm_seedMetAnno <- result_df_rm_seedMetAnno %>% arrange(feature_mz, feature_rt)
    fs <- result_df_rm_seedMetAnno$feature_name %>% unique()
    check_seed_filter <- c()
    check_known_filter <- c()
    check_candi_filter <- c()
    list_res_anno <- pbapply::pblapply(fs, function(temp_f) {
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
        # keep name unique by consensus name
        temp_res <- temp_res %>% arrange(id) %>% distinct(name, .keep_all = T) %>%
            arrange(desc(total_score))
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

    # remove redundancy: N to 1 --------------------------------------------------------------------
    if (is_clean_Nto1) {
        # reconstruct RT prediction model
        table_mrn_rt <- result_df_rm_1toN %>% select(id, feature_rt) %>% arrange(id)
        # remove duplicated ids
        dup_ids <- table_mrn_rt$id[which(duplicated(table_mrn_rt$id))] %>% unique()
        table_mrn_rt <- table_mrn_rt %>%
            dplyr::filter(!(id %in% dup_ids))
        table_mrn_rt <- table_mrn_rt[table_mrn_rt$id %in% rownames(md_mrn), ]
        cat("There are ", nrow(table_mrn_rt), " metabolites are used for RT re-prediction.\n", sep = "")
        ### filter ID MRN
        md <- md_mrn[match(table_mrn_rt$id, rownames(md_mrn)), ]
        ### remove NA which appear in more than 50% metabolites
        md1 <- md
        remove.idx1 <- which(apply(md, 2, function(x) {sum(is.na(x)/nrow(md))}) > 0.5)
        if (length(remove.idx1) > 0) md1 <- md[, -remove.idx1]
        ### impute NA
        md2 <- t(impute::impute.knn(data = t(md1))[[1]])
        ### remove MD which are same in all metabolites
        md3 <- md2
        remove.idx2 <- which(apply(md2, 2, sd) == 0)
        if (length(remove.idx2) > 0) md3 <- md2[, -remove.idx2]
        # construct RF model
        train.y <- table_mrn_rt$feature_rt
        train.x <- md3
        rm(list = c('md', 'md1', 'md2', 'md3')); gc()
        switch(
            column,
            "hilic" = {
                marker.name <- c('XLogP', "tpsaEfficiency", "WTPT.5", "khs.dsCH", "MLogP", "nAcid", "nBase", "BCUTp.1l")
            },
            "rp" = {
                marker.name <- c('XLogP', "WTPT.4", "WTPT.5", "ALogp2", "BCUTp.1l")
            }
        )
        idx <- match(marker.name, colnames(train.x))
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0) stop("Your markers are not in MD data.\n")
        train.x <- train.x[, idx]
        cat('Tuning RF function using caret approach.\n')
        train.data <- data.frame(RT = train.y, train.x, stringsAsFactors = FALSE)
        para_control <- caret::trainControl(method = 'cv', number = 10, search = 'random', verbose = F)
        set.seed(100)
        rf.reg <- caret::train(
            RT ~ ., data = train.data,
            method = 'rf',
            metric = 'Rsquared',
            tuneLength = 10,
            trControl = para_control,
            importance = T,
            allowParallel = T
        )
        # test RT prediction results
        test.x <- md_mrn[match(table_mrn_rt$id, rownames(md_mrn)), idx]
        temp.rt <- predict(object = rf.reg, newdata = test.x)
        table_mrn_rt$lib_pred_rt_2 <- temp.rt
        rt_exp <- table_mrn_rt$feature_rt
        rt_pred <- table_mrn_rt$lib_pred_rt_2
        r_squared <- 1 - sum((rt_exp - rt_pred)^2) / sum((rt_exp - mean(rt_pred))^2)
        mse <- mean((rt_exp - rt_pred)^2)
        rmse <- sqrt(mse)
        mae <- mean(abs(rt_exp - rt_pred))
        medae <- median(abs(rt_exp - rt_pred))
        mre <- mean(abs((rt_exp - rt_pred) / rt_exp)) * 100
        medre <- median(abs((rt_exp - rt_pred) / rt_exp)) * 100
        dir.create(path = file.path(wd_output, '00_intermediate_data'), showWarnings = F, recursive = T)
        pdf(file.path(wd_output, '00_intermediate_data/rt_prediction_plot_2.pdf'), width = 6, height = 6)
        max_val <- max(c(rt_exp, rt_pred))
        plot(
            x = rt_exp,
            y = rt_pred,
            xlim = c(0, max_val),
            ylim = c(0, max_val),
            main = "RT prediction of seeds & MRN annotations",
            xlab = "Experimental RT (s)",
            ylab = "Predicted RT (s)",
            pch = 19, col = "blue"
        )
        abline(a = 0, b = 1, lty = 2, col = "red")
        text(x = max_val * 0.05, y = max_val * 0.9, labels = paste0("R-squared: ", round(r_squared, 4)), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.85, labels = paste0("MeanAE: ", round(mae, 2), " s"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.8, labels = paste0("MedianAE: ", round(medae, 2), " s"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.75, labels = paste0("MeanRE: ", round(mre, 2), "%"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.7, labels = paste0("MedianRE: ", round(medre, 2), "%"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.65, labels = paste0("RMSE: ", round(rmse, 2), " s"), pos = 4)
        dev.off()

        # predict RT for MRN results
        test.x <- md_mrn[sort(unique(result_df_rm_1toN$id)), idx]
        mrn_rt <- rep(NA, nrow(test.x))
        names(mrn_rt) <- rownames(test.x)
        ### remove unpredictable metabolite and impute NA in text.x
        idx1 <- which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
        test.x1 <- test.x[idx1, ]
        test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
        ### predict RT
        temp.rt <- predict(object = rf.reg, newdata = test.x1)
        mrn_rt[idx1] <- temp.rt
        ### remove NA
        mrn_rt <- mrn_rt[!is.na(mrn_rt)]

        # check delta RT before-after re-prediction
        # exclude seeds (feature RT = predicted RT)
        idx_level3 <- which(result_df_rm_1toN$tag != 'seed')
        rt_after <- mrn_rt[match(result_df_rm_1toN$id[idx_level3], names(mrn_rt))]
        rt_before <- result_df_rm_1toN$rt[idx_level3]
        # remove NA
        idx_na <- unname(which(is.na(rt_after)))
        if (length(idx_na) > 0) {
            rt_after[idx_na] <- 0
            rt_before[idx_na] <- 0
        }
        mae <- mean(abs(rt_before - rt_after))
        medae <- median(abs(rt_before - rt_after))
        mre <- mean(abs((rt_before - rt_after) / rt_before)) * 100
        medre <- median(abs((rt_before - rt_after) / rt_before)) * 100
        pdf(file.path(wd_output, '00_intermediate_data/rt_prediction_comparsion.pdf'), width = 6, height = 6)
        max_val <- max(c(rt_after, rt_before))
        plot(
            x = rt_after,
            y = rt_before,
            xlim = c(0, max_val),
            ylim = c(0, max_val),
            main = "RT prediction comparsion",
            xlab = "Re-predicted RT (s)",
            ylab = "Predicted RT (s)",
            pch = 19, col = "blue"
        )
        abline(a = 0, b = 1, lty = 2, col = "red")
        text(x = max_val * 0.05, y = max_val * 0.85, labels = paste0("MeanAE: ", round(mae, 2), " s"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.8, labels = paste0("MedianAE: ", round(medae, 2), " s"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.75, labels = paste0("MeanRE: ", round(mre, 2), "%"), pos = 4)
        text(x = max_val * 0.05, y = max_val * 0.7, labels = paste0("MedianRE: ", round(medre, 2), "%"), pos = 4)
        dev.off()

        # keep only 1 smallest RT error for N TO 1 results
        ids <- result_df_rm_1toN$id %>% unique() %>% sort()
        keep_features <- sapply(ids, function(temp_id) {
            temp_res <- result_df_rm_1toN %>% filter(id == temp_id)
            if (nrow(temp_res) <= 1) return(paste0(temp_res$feature_name, '--', temp_res$id))
            if (temp_id %in% names(rt_after)) {
                temp_feature_rt <- temp_res$feature_rt
                temp_check_rt <- rt_after[temp_id]
                idx_closed <- which.min(abs(temp_feature_rt - temp_check_rt))[1]
                return(paste0(temp_res$feature_name[idx_closed], '--', temp_res$id[idx_closed]))
            } else return(paste0(temp_res$feature_name, '--', temp_res$id))
        }) %>% unlist()
        all_f_m_pairs <- paste0(result_df_rm_1toN$feature_name, '--', result_df_rm_1toN$id)
        # keep seed results
        keep_features <- c(keep_features, all_f_m_pairs[result_df_rm_1toN$tag == 'seed']) %>% unique()
        result_df_rm_1toN_Nto1 <- result_df_rm_1toN[all_f_m_pairs %in% keep_features, ]
        removed_Nto1_results <- result_df_rm_1toN[!(all_f_m_pairs %in% keep_features), ]
        readr::write_csv(removed_Nto1_results, file.path(wd_output, 'removed_redun_results_Nto1_part.csv'), na = '')

    } else {
        result_df_rm_1toN_Nto1 <- result_df_rm_1toN
    }

    # only confidence assignment (same with MetDNA)
    result_df_confi_assign <- MrnAnnoAlgo3::rmRedunInMrn3(
        result_df = result_df_rm_1toN_Nto1,
        is_only_confidence_assignment = T
    )
    result_df_confi_assign <- result_df_confi_assign %>% arrange(feature_mz, feature_rt)
    all_f_m_pairs_raw <- paste0(result_df$feature_name, '--', result_df$id)
    all_f_m_pairs_final <- paste0(result_df_confi_assign$feature_name, '--', result_df_confi_assign$id)
    removed_results <- result_df[!(all_f_m_pairs_raw %in% all_f_m_pairs_final), ]
    readr::write_csv(removed_results, file.path(wd_output, 'removed_redun_results.csv'), na = '')

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

}
