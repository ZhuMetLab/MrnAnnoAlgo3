
# rmRedunInMetdna3 ----------------------------------------------------------------------------

#' rmRedunInMetdna3
#'
#' @param result_df result_mrn_anno table (feature-metabolite pairs)
#' @param wd_output file.path(wd, '02_result_MRN_annotation')
#' @param num_candi number of candidates (topN) for redundancy removal
#' @param is_known_prioritization priority: level3.1 (knowns) > level3.2 (unknowns) for redundancy removal
#' @param is_clean_Nto1 re-regression by seeds & level3 results
#' @param md_mrn 'md_mrn' for RT re-regression
#' @param column 'column' for RT re-regression
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
        md_mrn = NULL, # for RT re-regression!
        column = NULL # for RT re-regression!
    ) {

    cat('Redundancy removal by MetDNA3 (MrnAnnoAlgo3):', '\n')
    cat('\n')

    # refine metabolite name with consensus name ---------------------------------------------------
    tbl_cons_name <- readr::read_csv(system.file('mrn_consensus_name.csv', package = 'MrnAnnoAlgo3'))
    idx_in_tbl <- which(result_df$id %in% tbl_cons_name$id_mrn)
    if (length(idx_in_tbl) > 0) {
        cons_name_in_tbl <- match(result_df$id[idx_in_tbl], tbl_cons_name$id_mrn) %>%
            tbl_cons_name$consensus_name[.]
        result_df$name[idx_in_tbl] <- cons_name_in_tbl
    }

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

    # remove redundancy by feature -----------------------------------------------------------------
    fs <- result_df_rm_seedMetAnno$feature_name %>% unique()
    list_res_anno <- lapply(fs, function(temp_f) {
        temp_res <- result_df_rm_seedMetAnno %>% filter(feature_name == temp_f) %>% arrange(desc(total_score), id)
        # keep seed first (same as grade1)
        if ('seed' %in% temp_res$tag) temp_res <- temp_res %>% filter(tag == 'seed')
        # keep level 3.1 first (grade 3.1)
        if (is_known_prioritization & any(startsWith(temp_res$id, 'MRN'))) temp_res <- temp_res[startsWith(temp_res$id, 'MRN'), ]
        # keep name unique by consensus name
        temp_res <- temp_res %>% arrange(id) %>% distinct(name, .keep_all = T) %>% arrange(desc(total_score))
        return(temp_res)
    })
    uni_fs <- fs[which(sapply(list_res_anno, nrow) == 1)]
    result_df_rm_by_f <- bind_rows(list_res_anno)
    cat('Remove annotation redundancy by feature-level:', nrow(result_df_rm_seedMetAnno) - nrow(result_df_rm_by_f), '\n')
    cat('\n')
    f_m_pairs_1 <- paste0(result_df_rm_by_f$feature_name, '--', result_df_rm_by_f$id)
    f_m_pairs_all <- paste0(result_df$feature_name, '--', result_df$id)
    removed_redun_results_feature_part <- result_df[!f_m_pairs_all %in% f_m_pairs_1, ]
    readr::write_csv(removed_redun_results_feature_part, file.path(wd_output, 'removed_redun_results_feature_part.csv'), na = '')

    # reconstruct RT prediction model (for 1 to 1 results) -----------------------------------------
    mrn_rts <- reconstructRtPredModel(
        result_df_rm_by_f = result_df_rm_by_f,
        md_mrn = md_mrn,
        column = column,
        wd_output = wd_output
    )
    res_df_seed <- result_df_rm_by_f %>% filter(tag == 'seed') # seedAnno no change !!!
    res_df_mrnAnno <- result_df_rm_by_f %>% filter(tag != 'seed') # check MRN annotation results !!!

    # remove redundancy: N to 1 --------------------------------------------------------------------
    # check before
    to_check_res_df <- result_df_rm_by_f %>% filter(tag != 'seed') %>%
        select(feature_name, feature_rt, id, rt)
    to_check_res_df$pred_rt <- mrn_rts[to_check_res_df$id]
    dup_ids <- to_check_res_df$id[duplicated(to_check_res_df$id)]
    to_check_res_df$group <- 'others'
    to_check_res_df$group[to_check_res_df$id %in% dup_ids] <- 'N to 1'
    max_val <- max(c(to_check_res_df$pred_rt, to_check_res_df$feature_rt)) * 1.05
    library(ggplot2)
    p_raw <- ggplot2::ggplot(to_check_res_df, aes(x = feature_rt, y = pred_rt, group = group, color = group)) +
        geom_point(size = 1.5) +
        geom_abline(intercept = 0, slope = 1.3, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1.1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
        geom_abline(intercept = 0, slope = 0.9, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 0.7, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(x = "Feature RT (s)", y = "Re-predicted RT (s)", title = "RT plot of metAnnotations") +
        xlim(c(0, max_val)) +
        ylim(c(0, max_val)) +
        coord_fixed() +
        theme_bw()
    ggsave(filename = file.path(wd_output, '00_intermediate_data/rt_plot_rm_Nto1_before.pdf'), plot = p_raw, height = 6, width = 6)
    # execution
    ids <- res_df_mrnAnno$id %>% unique() %>% sort()
    list_metAnno <- lapply(ids, function(temp_id) {
        temp_res <- res_df_mrnAnno %>% filter(id == temp_id)
        if (nrow(temp_res) == 1) return(temp_res)
        temp_rt_repred <- mrn_rts[temp_id]
        idx_keep <- which.min(abs(temp_res$feature_rt - temp_rt_repred))[1]
        return(temp_res[idx_keep, ])
    })
    res_df_mrnAnno_rm_Nto1 <- bind_rows(list_metAnno)
    res_df_mrnAnno_rm_Nto1 <- res_df_mrnAnno_rm_Nto1 %>% arrange(feature_mz, feature_rt)
    cat('Remove annotation redundancy of N to 1:', nrow(res_df_mrnAnno) - nrow(res_df_mrnAnno_rm_Nto1), '\n')
    cat('\n')
    f_m_pairs_mrnAnno_rm_Nto1 <- paste0(res_df_mrnAnno_rm_Nto1$feature_name, '--', res_df_mrnAnno_rm_Nto1$id)
    f_m_pairs_mrnAnno <- paste0(res_df_mrnAnno$feature_name, '--', res_df_mrnAnno$id)
    removed_redun_results_Nto1_part <- res_df_mrnAnno[!f_m_pairs_mrnAnno %in% f_m_pairs_mrnAnno_rm_Nto1, ]
    readr::write_csv(removed_redun_results_Nto1_part, file.path(wd_output, 'removed_redun_results_Nto1_part.csv'), na = '')
    # check after rm. N to 1
    to_check_res_df <- res_df_mrnAnno_rm_Nto1 %>% filter(tag != 'seed') %>%
        select(feature_name, feature_rt, id, rt)
    to_check_res_df$pred_rt <- mrn_rts[to_check_res_df$id]
    to_check_res_df$group <- 'others'
    to_check_res_df$group[to_check_res_df$id %in% dup_ids] <- 'N to 1' # old 'dup_ids'
    p_rm_Nto1 <- ggplot2::ggplot(to_check_res_df, aes(x = feature_rt, y = pred_rt, group = group, color = group)) +
        geom_point(size = 1.5) +
        geom_abline(intercept = 0, slope = 1.3, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1.1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
        geom_abline(intercept = 0, slope = 0.9, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 0.7, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(x = "Feature RT (s)", y = "Re-predicted RT (s)", title = "RT plot of metAnnotations (rm. Nto1)") +
        xlim(c(0, max_val)) +
        ylim(c(0, max_val)) +
        coord_fixed() +
        theme_bw()
    ggsave(filename = file.path(wd_output, '00_intermediate_data/rt_plot_rm_Nto1_after.pdf'), plot = p_rm_Nto1, height = 6, width = 6)

    # remove redundancy: 1 to N --------------------------------------------------------------------
    # check before
    to_check_res_df <- to_check_res_df %>% arrange(feature_name) # same as before rm. N to 1
    to_check_res_df$group <- 'others'
    dup_fs <- to_check_res_df$feature_name[duplicated(to_check_res_df$feature_name)] %>% unique() %>% sort()
    to_check_res_df$group[to_check_res_df$feature_name %in% dup_fs] <- '1 to N'
    p_raw <- ggplot2::ggplot(to_check_res_df, aes(x = feature_rt, y = pred_rt, group = group, color = group)) +
        geom_point(size = 1.5) +
        geom_abline(intercept = 0, slope = 1.3, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1.1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
        geom_abline(intercept = 0, slope = 0.9, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 0.7, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(x = "Feature RT (s)", y = "Re-predicted RT (s)", title = "RT plot of metAnnotations") +
        xlim(c(0, max_val)) +
        ylim(c(0, max_val)) +
        coord_fixed() +
        theme_bw()
    ggsave(filename = file.path(wd_output, '00_intermediate_data/rt_plot_rm_1toN_before.pdf'), plot = p_raw, height = 6, width = 6)
    # execution
    fs <- res_df_mrnAnno_rm_Nto1$feature_name %>% unique() %>% sort()
    list_metAnno <- lapply(fs, function(temp_f) {
        temp_res <- res_df_mrnAnno_rm_Nto1 %>% filter(feature_name == temp_f) %>% arrange(desc(total_score), id)
        if (nrow(temp_res) <= num_candi) return(temp_res)
        return(temp_res[1:num_candi,])
    })
    res_df_mrnAnno_rm_1toN <- bind_rows(list_metAnno)
    cat('Remove annotation redundancy of 1 to N:', nrow(res_df_mrnAnno_rm_Nto1) - nrow(res_df_mrnAnno_rm_1toN), '\n')
    cat('\n')
    f_m_pairs_mrnAnno_rm_Nto1 <- paste0(res_df_mrnAnno_rm_Nto1$feature_name, '--', res_df_mrnAnno_rm_Nto1$id)
    f_m_pairs_mrnAnno_rm_1toN <- paste0(res_df_mrnAnno_rm_1toN$feature_name, '--', res_df_mrnAnno_rm_1toN$id)
    removed_redun_results_1toN_part <- res_df_mrnAnno_rm_Nto1[!f_m_pairs_mrnAnno_rm_Nto1 %in% f_m_pairs_mrnAnno_rm_1toN, ]
    readr::write_csv(removed_redun_results_1toN_part, file.path(wd_output, 'removed_redun_results_1toN_part.csv'), na = '')
    # check after rm. 1 to N
    to_check_res_df <- res_df_mrnAnno_rm_1toN %>% filter(tag != 'seed') %>%
        select(feature_name, feature_rt, id, rt)
    to_check_res_df$pred_rt <- mrn_rts[to_check_res_df$id]
    to_check_res_df$group <- 'others'
    to_check_res_df$group[to_check_res_df$feature_name %in% dup_fs] <- '1 to N' # old 'dup_fs'
    p_rm_1toN <- ggplot2::ggplot(to_check_res_df, aes(x = feature_rt, y = pred_rt, group = group, color = group)) +
        geom_point(size = 1.5) +
        geom_abline(intercept = 0, slope = 1.3, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1.1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
        geom_abline(intercept = 0, slope = 0.9, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 0.7, linetype = "dashed", color = "black", linewidth = 0.5) +
        labs(x = "Feature RT (s)", y = "Re-predicted RT (s)", title = "RT plot of metAnnotations (rm. 1toN)") +
        xlim(c(0, max_val)) +
        ylim(c(0, max_val)) +
        coord_fixed() +
        theme_bw()
    ggsave(filename = file.path(wd_output, '00_intermediate_data/rt_plot_rm_1toN_after.pdf'), plot = p_rm_1toN, height = 6, width = 6)

    # combine seedAnno and metAnno -----------------------------------------------------------------
    result_df_rm_redun_ok <- bind_rows(res_df_seed, res_df_mrnAnno_rm_1toN) %>% arrange(feature_mz, feature_rt)

    # only confidence assignment (same with MetDNA) ------------------------------------------------
    result_df_confi_assign <- MrnAnnoAlgo3::rmRedunInMrn3(
        result_df = result_df_rm_redun_ok,
        is_only_confidence_assignment = T
    )
    result_df_confi_assign <- result_df_confi_assign %>% arrange(feature_mz, feature_rt)

    ### calculate redundancy (same with MetDNA) ----------------------------------------------------
    ids_anno <- sort(unique(result_df_confi_assign$id))
    peaks_anno <- sort(unique(result_df_confi_assign$feature_name))
    redun <- calculateRedundancyInMrn3(result_df_confi_assign)
    cat('Annotated peaks:', length(peaks_anno), '\n')
    cat('Annotated metabolites:', length(ids_anno), '\n')
    cat('Annotations:', nrow(result_df_confi_assign), '\n')
    cat("Peak redundancy:", round(redun[1], 2), '\n')
    cat("Metabolite redundancy:", round(redun[2], 2), '\n')
    cat('\n')

    return(result_df_confi_assign)

}


# reconstructRtPredModel ----------------------------------------------------------------------

#' reconstructRtPredModel
#'
#' @param result_df_rm_by_f 'result_df_rm_by_f' of rmRedunInMetdna3 function
#' @param md_mrn 'md_mrn' of rmRedunInMetdna3 function
#' @param column 'column' of rmRedunInMetdna3 function
#' @param wd_output 'wd_output' of rmRedunInMetdna3 function
#'
#' @return MRN RTs
#' @export
#'
#' @examples reconstructRtPredModel(result_df_rm_by_f)
reconstructRtPredModel <- function(result_df_rm_by_f, md_mrn, column, wd_output) {
    table_mrn_rt <- result_df_rm_by_f %>%
        # filter(feature_name %in% uni_fs) %>%
        select(id, feature_rt)
    # remove duplicated ids
    dup_ids <- table_mrn_rt$id[which(duplicated(table_mrn_rt$id))] %>% unique()
    table_mrn_rt <- table_mrn_rt %>% dplyr::filter(!(id %in% dup_ids))
    # remove duplicated rts (n > 5)
    rt_summary <- table_mrn_rt %>% group_by(feature_rt) %>% summarise(n = n())
    rm_rts <- rt_summary$feature_rt[which(rt_summary$n > 5)]
    table_mrn_rt <- table_mrn_rt %>% dplyr::filter(!(feature_rt %in% rm_rts))
    # keep valid id in md_mrn
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
    test.x <- md_mrn[sort(unique(result_df_rm_by_f$id)), idx]
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
    return(mrn_rt)

}
