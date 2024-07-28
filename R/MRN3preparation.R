
# getMs1Data4Mrn3 -----------------------------------------------------------------------------

#' getMs1Data4Mrn3
#'
#' @param wd Working directory (default is the current directory).
#' @param ms1_filename Filename of the MS1 data (default is 'data.csv').
#' @param is_output output ms1_data?
#'
#' @return A tibble containing the MS1 data in MRN3 format with columns 'name', 'mz', 'rt', and 'ccs'.
#' @export
#'
#' @examples
#' getMs1Data4Mrn3(wd = 'TO/YOUR/PATH/', ms1_filename = 'data.csv')
#'
getMs1Data4Mrn3 <- function(
        wd = '.',
        ms1_filename = 'data.csv',
        is_output = TRUE
        ) {

    library(dplyr)

    wd_output <- file.path(wd, '02_result_MRN_annotation')
    dir.create(wd_output, showWarnings = F, recursive = T)

    # convert to MRN3 ms1_data format [, c('name', 'mz', 'rt', 'ccs')]
    ms1_data <- read.csv(file.path(wd, ms1_filename))
    if (!('ccs' %in% colnames(ms1_data))) ms1_data[, 'ccs'] <- -1
    ms1_data <- ms1_data[, c('name', 'mz', 'rt', 'ccs')] %>% tibble::as_tibble()
    if (is_output) readr::write_csv(ms1_data, file.path(wd_output, 'ms1_data.csv'))

    cat('Total features:', nrow(ms1_data), '\n')

    return(ms1_data)

}

# getMs2Data4Mrn3 -----------------------------------------------------------------------------

#' getMs2Data4Mrn3
#'
#' Converts MS2 spectra data to MRN3 format.
#'
#' @param wd Working directory (default is the current directory).
#' @param ms2_filename Filename of the MS2 spectra data (default is 'spectra.msp').
#' @param is_output output ms2_data?
#' @param is_use_01_ms2 use MS2 in 01_result_initial_seed_annotation/00_intermediate_data/ms2
#'
#' @return A list containing MS2 spectra data in MRN3 format.
#' @export
#'
#' @examples
#' getMs2Data4Mrn3(wd = 'TO/YOUR/PATH/', ms2_filename = 'spectra.msp')
#'
getMs2Data4Mrn3 <- function(
        wd = '.',
        ms2_filename = 'spectra.msp',
        is_output = TRUE,
        is_use_01_ms2 = TRUE
        ) {

    library(dplyr)

    wd_output <- file.path(wd, '02_result_MRN_annotation')
    dir.create(wd_output, showWarnings = F, recursive = T)

    # convert to MRN3 ms2_data format (feature name matched)
    if (is_use_01_ms2) {
        load(file.path(wd, '01_result_initial_seed_annotation/00_intermediate_data/ms2'))
        ms2_data <- lapply(ms2, function(temp_ms2) {
            temp_ms2$info <- t(temp_ms2$info)
            return(temp_ms2)
        })
    } else {
        ms2_data <- SpectraTools::ParseMSP(file.path(wd, ms2_filename))
    }
    if (is_output) save(ms2_data, file = file.path(wd_output, 'ms2_data.rda'))

    cat('Total MS2 spectra:', length(ms2_data), '\n')

    return(ms2_data)

}

# getTableSeedLong4Mrn3 -----------------------------------------------------------------------

#' getTableSeedLong4Mrn3
#'
#' Reads the seed annotation file and converts it to a long-format table for MRN3 annotation.
#'
#' @param wd Working directory (default is the current directory).
#' @param seedAnno_filename Filename of the seed annotation file (default is '01_result_initial_seed_annotation/ms2_match_annotation_result.csv').
#' @param id_type ID conversion to "id_mrn" or "id_kegg"
#' @param is_output output seed table (long-format)?
#'
#' @return A tibble containing the long-format seed table for MRN3 annotation with columns 'feature_name' and 'id_mrn'.
#' @export
#'
#' @examples
#' getTableSeedLong4Mrn3(wd = 'TO/YOUR/PATH/', seedAnno_filename = '01_result_initial_seed_annotation/ms2_match_annotation_result.csv')
#'
getTableSeedLong4Mrn3 <- function(
        wd = '.',
        seedAnno_filename = '01_result_initial_seed_annotation/ms2_match_annotation_result.csv',
        id_type = c('id_mrn', 'id_kegg'),
        is_output = TRUE
        ) {

    id_type <- match.arg(id_type)

    library(dplyr)

    wd_output <- file.path(wd, '02_result_MRN_annotation')
    dir.create(wd_output, showWarnings = F, recursive = T)

    # get seed table (long) for MRN3 annotation
    seed_anno <- read.csv(file.path(wd, seedAnno_filename)) %>% tibble::as_tibble()
    idx_feature_seed <- which(!is.na(seed_anno$id_reverse_summary))
    id_raw_seed <- lapply(idx_feature_seed, function(x) {
        temp_result <- seed_anno$id_reverse_summary[x]
        labids <- stringr::str_extract_all(temp_result, "labid\\{(.*?)\\}")[[1]] %>%
            sapply(., function(x) stringr::str_replace_all(x, "labid\\{|\\}", "")) %>% unname()
        # cat(x, labids, '\n')
        return(labids)
    })
    metinfo <- MetLib::metinfo #### MetLib: id_zhulab to id_mrn or id_kegg!!!
    id_mrn_seed <- lapply(id_raw_seed, function(x) unique(sort(unlist(metinfo[, id_type])[match(x, metinfo$id)])))
    names(id_mrn_seed) <- seed_anno$name[idx_feature_seed]
    id_mrn_seed_len <- unname(sapply(id_mrn_seed, length))
    id_mrn_seed <- id_mrn_seed[which(id_mrn_seed_len > 0)]
    id_mrn_seed_pair <- sort(unlist(id_mrn_seed))
    id_mrn_seed_all <- unique(id_mrn_seed_pair)
    cat('Raw annotated features:', length(id_raw_seed), '\n')
    cat('Annotated features:', length(id_mrn_seed), '\n')
    cat('Annotated feature-metabolite pairs:', length(id_mrn_seed_pair), '\n')
    cat('Annotated unique metabolties:', length(id_mrn_seed_all), '\n')
    table_seed_wide <- tibble::tibble(
        feature_name = names(id_mrn_seed),
        id_mrn = sapply(id_mrn_seed, function(x) stringr::str_c(x, collapse = ';'))
    )
    table_seed_long <- table_seed_wide %>%
        tidyr::separate_rows(id_mrn, sep = ";") %>%
        tidyr::separate_rows(id_mrn, sep = ", ") %>% # for MetLib/MetDNA2 id_kegg
        tidyr::separate_rows(id_mrn, sep = ",") # for MetLib/MetDNA2 id_kegg

    if (is_output) readr::write_csv(table_seed_long, file.path(wd_output, 'table_seed_long.csv'))

    return(table_seed_long)

}

# predictRtInMrn3 -----------------------------------------------------------------------------

#' predictRtInMrn3
#'
#' @param wd Working directory
#' @param ms1_data ms1_data
#' @param table_seed_long table_seed_long
#' @param md_mrn md_mrn
#' @param column column
#'
#' @return "mrn_rt" (Named predicted RT vector of rownames(md_mrn))
#' @export
#'
#' @examples predictRtInMrn3(wd, ms1_data, table_seed_long, md_mrn, column = 'hilic')
#'
predictRtInMrn3 <- function(
        wd,
        ms1_data,
        table_seed_long,
        md_mrn,
        column = c('hilic', 'rp')
        ) {

    # MetDNA2 method!

    library(dplyr)
    library(randomForest)

    column <- match.arg(column)

    # prepare RT training data set
    cat("RT prediction.", "\n")
    cat("\n")
    table_exp_rt <- table_seed_long
    table_exp_rt$rt_exp <- ms1_data$rt[match(table_exp_rt$feature_name, ms1_data$name)]
    table_exp_rt <- table_exp_rt[table_exp_rt$id_mrn %in% rownames(md_mrn), ]
    table_exp_rt <- table_exp_rt %>% dplyr::arrange(id_mrn, rt_exp)
    table_exp_rt <- table_exp_rt %>%
        dplyr::group_by(id_mrn) %>%
        dplyr::summarize(rt_exp = median(rt_exp, na.rm = T))
    if (nrow(table_exp_rt) == 0) stop("ERROR: No metabolites are identified.\n")
    if (nrow(table_exp_rt) < 20) stop("ERROR: < 20 metabolites are identified.\n")
    cat("There are ", nrow(table_exp_rt), " metabolites are used for RT prediction.\n", sep = "")

    # MD cleaning
    ### filter ID MRN
    md <- md_mrn[match(table_exp_rt$id_mrn, rownames(md_mrn)), ]
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
    train.y <- table_exp_rt$rt_exp
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
    # rf.reg <- MetDNA2:::trainRandomForestWithCaret(x = train.data)
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

    # predict RT for MRN database
    test.x <- md_mrn[, idx]
    mrn_rt <- rep(NA, nrow(test.x))
    names(mrn_rt) <- rownames(test.x)
    ### remove unpredictable metabolite and impute NA in text.x
    idx1 <- which(apply(test.x, 1, function(x) {sum(is.na(x))/ncol(test.x) < 0.5}))
    test.x1 <- test.x[idx1, ]
    test.x1 <- t(impute::impute.knn(data = t(test.x1))[[1]])
    ### predict RT
    temp.rt <- predict(object = rf.reg, newdata = test.x1)
    mrn_rt[idx1] <- temp.rt
    ### NA give the median RT of all peaks
    mrn_rt[is.na(mrn_rt)] <- median(ms1_data$rt)
    table_exp_rt$rt_pred <- temp.rt[match(table_exp_rt$id_mrn, names(temp.rt))]

    ### summary
    rt_exp <- table_exp_rt$rt_exp
    rt_pred <- table_exp_rt$rt_pred
    r_squared <- 1 - sum((rt_exp - rt_pred)^2) / sum((rt_exp - mean(rt_pred))^2)
    mse <- mean((rt_exp - rt_pred)^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(rt_exp - rt_pred))
    medae <- median(abs(rt_exp - rt_pred))
    mre <- mean(abs((rt_exp - rt_pred) / rt_exp)) * 100
    medre <- median(abs((rt_exp - rt_pred) / rt_exp)) * 100
    cat('RT prediction result:', '\n')
    cat('R-squared:', round(r_squared, 4), '\n')
    cat('Mean Absolute Error (MAE):', round(mae, 4), '\n')
    cat('Median Absolute Error (MedAE):', round(medae, 4), '\n')
    cat('Mean Relative Error (MRE):', round(mre, 4), '\n')
    cat('Median Relative Error (MedRE):', round(medre, 4), '\n')
    cat('Root Mean Squared Error (RMSE):', round(rmse, 4), '\n')
    cat('\n')

    cat('RT prediction finish.\n')
    cat('\n')

    # save "rt_result" for MetDNA result generation
    # TODO: modify MetDNA functions (structure issue)
    rt_result <- list(
        'KEGG.rt' = data.frame('RT' = mrn_rt, row.names = names(mrn_rt))
    )
    wd_output <- file.path(wd, "02_result_MRN_annotation", "00_intermediate_data")
    dir.create(wd_output, showWarnings = FALSE, recursive = TRUE)
    save(
        rt_result,
        file = file.path(wd_output, "rt_result"),
        compress = "xz", version = 2
    )
    pdf(file.path(wd_output, 'rt_prediction_plot.pdf'), width = 6, height = 6)
    max_val <- max(c(rt_exp, rt_pred))
    plot(
        x = rt_exp,
        y = rt_pred,
        xlim = c(0, max_val),
        ylim = c(0, max_val),
        main = "RT prediction of seeds",
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

    return(mrn_rt)

}
