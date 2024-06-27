
# getAdductAnno -------------------------------------------------------------------------------

#' getAdductAnno
#'
#' @param table_anno table_mrn_anno (total)
#' @param idx_seeds idx_seeds
#' @param round round
#' @param rt_tol_adductAnno rt_tol_adductAnno = 3 s
#'
#' @return table_mrn_anno (after adduct annotation)
#' @export
#'
#' @examples table_mrn_anno <- getAdductAnno(table_anno = table_mrn_anno, idx_seeds = c(1, 2), round = 1, rt_tol_adductAnno = 3)
#'
getAdductAnno <- function(
        table_anno,
        idx_seeds,
        round,
        rt_tol_adductAnno = 3
    ) {

    # f-m pair info
    comb_f_m_all <- paste(table_anno$feature_name, table_anno$id, sep = '--')
    comb_f_m_annotated <- comb_f_m_all[!is.na(table_anno$tag)]

    list_adductAnno <- lapply(idx_seeds, function(idx) {

        # get info
        temp_f_name <- table_anno$feature_name[idx]
        temp_f_rt <- table_anno$feature_rt[idx]
        temp_id <- table_anno$id[idx]
        temp_ms2_score <- table_anno$ms2_score[idx]

        # filter:
        # 1. same id
        # 2. RT diff. <= rt_tol_adductAnno
        logic_adduct <- (table_anno$id == temp_id) &
            (abs(table_anno$feature_rt - temp_f_rt) <= rt_tol_adductAnno)
        idx_adduct <- which(logic_adduct)
        # idx_adduct <- idx_adduct[-which(idx_adduct == idx)] (check after result adding)
        if (length(idx_adduct) == 0) return(NULL)

        # save results
        temp_result <- table_anno[idx_adduct, ]
        # temp_result$tag <- paste0('annoRound', sprintf("%02d", round))
        temp_result$tag <- 'adductAnnotation'
        temp_result$level <- round # (level = round)
        temp_result$feature_from <- temp_f_name
        temp_result$id_from <- temp_id
        temp_result$ms2_score <- temp_ms2_score # adductAnno use from MS2 score
        temp_result$matched_frag <- NA # adductAnno
        temp_result$matched_nl <- NA # adductAnno

        return(temp_result)
    })

    table_adductAnno <- list_adductAnno %>% dplyr::bind_rows()

    # remove annotated f.-m. pairs
    comb_f_m_adductAnno <- paste(table_adductAnno$feature_name, table_adductAnno$id, sep = '--')
    idx_f_m_unannotated <- which(!(comb_f_m_adductAnno %in% comb_f_m_annotated))
    if (length(idx_f_m_unannotated) == 0) return(table_anno) # no new adduct annotations
    table_adductAnno <- table_adductAnno[idx_f_m_unannotated, ]

    # remove duplicated f.-m. pairs
    comb_f_m_adductAnno <- paste(table_adductAnno$feature_name, table_adductAnno$id, sep = '--')
    idx_keep_adductAnno <- !duplicated(comb_f_m_adductAnno)
    table_adductAnno_uni <- table_adductAnno[idx_keep_adductAnno, ]
    comb_f_m_adductAnno_uni <- comb_f_m_adductAnno[idx_keep_adductAnno]

    # mark to final result: table_anno
    idx_to_mark <- match(comb_f_m_adductAnno_uni, comb_f_m_all)
    table_anno[idx_to_mark, ] <- table_adductAnno_uni

    return(table_anno)

}


# getMetAnno ----------------------------------------------------------------------------------

#' getMetAnno
#'
#' @param table_anno table_mrn_anno (total)
#' @param idx_seeds idx_seeds
#' @param round round
#' @param obj_mrn MRN object
#' @param table_ms2_network_filter table_ms2_network_filter
#' @param step_mrn step_mrn = 3 (MetDNA2/MRN2), = 1 (MetDNA3/MRN3)
#' @param topN_neighbor topN_neighbor = Inf (100 in MetDNA2)
#' @param lib_type_mrn "MRN3" or "MRN2"
#' @param thread thread = 4
#'
#' @return table_mrn_anno (after metabolite annotation)
#' @export
#'
#' @examples table_mrn_anno <- getMetAnno(...)
#'
getMetAnno <- function(
        table_anno,
        idx_seeds,
        round,
        obj_mrn,
        table_ms2_network_filter,
        step_mrn = 1,
        topN_neighbor = Inf,
        lib_type_mrn = 'MRN3',
        thread = 4
    ) {

    # f-m pair info
    comb_f_m_all <- paste(table_anno$feature_name, table_anno$id, sep = '--')
    comb_f_m_annotated <- comb_f_m_all[!is.na(table_anno$tag)]

    if (length(idx_seeds) < 64) thread <- min(4, thread)
    if (length(idx_seeds) < 16) thread <- min(2, thread)
    if (length(idx_seeds) < 4) thread <- 1

    library(parallel)
    cl <- parallel::makeCluster(thread)
    clusterExport(cl, c(
        # 'table_anno',
        # 'obj_mrn',
        # 'table_ms2_network_filter',
        # 'step_mrn',
        '%>%'
    ), envir = environment())

    # list_metAnno <- pbapply::pblapply(idx_seeds, function(idx) {
    tempMetAnnoFun <- function(idx) {

        cat(idx, '')

        # get basic info
        temp_f_name <- table_anno$feature_name[idx]
        temp_id <- table_anno$id[idx]

        # if id is not in MRN: next
        if (!(temp_id %in% igraph::V(obj_mrn)$name)) return(NULL)
        # get neighbor features (by MS2 similarity network) !!!
        temp_anno <- table_ms2_network_filter[(
            table_ms2_network_filter$from == temp_f_name |
                table_ms2_network_filter$to == temp_f_name
        ), ]
        if (nrow(temp_anno) == 0) return(NULL)
        # rerank 'from' and 'to' features
        idx_rerank <- which(temp_anno$from != temp_f_name)
        if (length(idx_rerank) > 0) {
            temp_anno$to[idx_rerank] <- temp_anno$from[idx_rerank]
            temp_anno$from[idx_rerank] <- temp_f_name
        }
        # add feature itself (MetDNA method)
        temp_anno <- temp_anno %>% dplyr::add_row(
            'from' = temp_f_name,
            'to' = temp_f_name,
            'ms2_score' = 1,
            'matchFrag' = 99,
            'matchNl' = 99
        )
        # get & separate neighbor features' ids
        temp_anno$from_id <- temp_id
        temp_anno$to_id <- sapply(temp_anno$to, function(temp_f) {
            table_anno$id[which(table_anno$feature_name == temp_f)] %>%
                paste(., collapse = ';')
        })
        temp_anno <- temp_anno %>% tidyr::separate_rows(to_id, sep = ";")

        # match ids by from id (MRN guided) !!!
        # step_mrn = 3 in MetDNA2(MRN2)
        # step_mrn = 1 in MrnAnnoAlgo3
        round_mrn <- 1
        round_id <- temp_id
        # cat('round: ')
        while (round_mrn <= step_mrn) {

            # cat(round_mrn, '')
            temp_to_id <- MrnAnnoAlgo3::getMrnNeighbor(id = temp_id, step = round_mrn, graph = obj_mrn)
            # use 'step == round_mrn' neighbors (MetDNA method) (by round_mrn)
            temp_to_id <- temp_to_id$id[temp_to_id$step == round_mrn]
            # 100 in MetDNA2
            if (length(temp_to_id) > topN_neighbor) temp_to_id <- sort(temp_to_id)[1:topN_neighbor]

            # # version 0.2.14.2.4: FCMN (Fully-connected metabolite network) # neighbors = all nodes temp!
            # temp_to_id <- sprintf("MRN%06d", 1:53583)

            # check valid result:
            temp_anno_to_check <- temp_anno %>% dplyr::filter(to_id %in% temp_to_id)
            if (nrow(temp_anno_to_check) == 0) {
                round_id <- temp_to_id
                round_mrn <- round_mrn + 1 # next round
            } else {
                # ONLY CHECK VALID ADDUCT & RT ERROR!!! (SAME WITH METDNA)
                comb_f_m_to_check <- paste(temp_anno_to_check$to, temp_anno_to_check$to_id, sep = '--')
                idx_check <- match(comb_f_m_to_check, comb_f_m_all)
                temp_res_to_check <-  table_anno[idx_check, ]
                valid_adducts <- c(
                    "[M+H]+",      "[M+NH4]+",    "[M+Na]+",   "[M-H+2Na]+",  "[M+K]+",
                    "[M-H+2K]+",   "[2M+H]+",     "[2M+NH4]+", "[2M+Na]+",    "[2M+K]+",
                    "[M-H2O+H]+",
                    "[M-H]-",      "[M+Na-2H]-",  "[M+K-2H]-", "[M+NH4-2H]-", "[2M-H]-",
                    "[M+CH3COO]-", "[2M+Na-2H]-", "[M-H2O-H]-"
                ) # 11 pos., 8 neg.
                valid_rt_error_rel <- 30 # (%)
                idx_adduct_ok <- temp_res_to_check$adduct %in% valid_adducts
                idx_rt_error_ok <- temp_res_to_check$rt_error_rel < valid_rt_error_rel
                idx_ok <- idx_adduct_ok & idx_rt_error_ok
                temp_anno_to_check <- temp_anno_to_check[idx_ok, ]
                if (nrow(temp_anno_to_check) == 0) {
                    round_id <- temp_to_id
                    round_mrn <- round_mrn + 1 # next round
                } else {
                    round_mrn <- step_mrn + 1 # end
                }
                # NO CHECK IN METDNA! (SAME WITH METDNA)
                # comb_f_m_to_check <- paste(temp_anno_to_check$to, temp_anno_to_check$to_id, sep = '--')
                # idx_f_m_annotated <- which(comb_f_m_to_check %in% comb_f_m_annotated)
                # if (length(idx_f_m_annotated) > 0) {
                #     comb_f_m_to_check <- comb_f_m_to_check[-idx_f_m_annotated]
                #     temp_anno_to_check <- temp_anno_to_check[-idx_f_m_annotated, ]
                # }
                # if (length(comb_f_m_to_check) == 0) {
                #     round_id <- temp_to_id
                #     round_mrn <- round_mrn + 1 # next round
                # } else {
                #     round_mrn <- step_mrn + 1 # end
                # }
            }
        }
        temp_anno <- temp_anno_to_check
        if (nrow(temp_anno) == 0) return(NULL)

        # mark feature-id combination to [table_anno]
        temp_comb_f_m <- paste(temp_anno$to, temp_anno$to_id, sep = '--')
        temp_idx_to_mark <- match(temp_comb_f_m, comb_f_m_all)

        # save results
        temp_result <- table_anno[temp_idx_to_mark, ]
        # temp_result$tag <- paste0('annoRound', sprintf("%02d", round))
        temp_result$tag <- 'metAnnotation'
        temp_result$level <- round # level = round
        temp_result$feature_from <- temp_f_name
        temp_result$id_from <- temp_id
        temp_result$ms2_score <- temp_anno$ms2_score
        temp_result$matched_frag <- temp_anno$matchFrag
        temp_result$matched_nl <- temp_anno$matchNl

        return(temp_result)

    }
    # })

    system.time(list_metAnno <- parallel::parLapply(cl, idx_seeds, tempMetAnnoFun))
    stopCluster(cl)

    table_metAnno <- list_metAnno %>% dplyr::bind_rows()

    if (nrow(table_metAnno) == 0) return(table_anno)

    # remove annotated f.-m. pairs (after MRN search for MetDNA method)
    comb_f_m_metAnno <- paste(table_metAnno$feature_name, table_metAnno$id, sep = '--')
    idx_f_m_unannotated <- which(!(comb_f_m_metAnno %in% comb_f_m_annotated))
    if (length(idx_f_m_unannotated) == 0) return(table_anno) # no new met. annotations
    table_metAnno <- table_metAnno[idx_f_m_unannotated, ]

    # remove 'peak -> peak' same peak problem (MetDNA method)
    idx_rm <- which(table_metAnno$feature_name == table_metAnno$feature_from)
    if (length(idx_rm) > 0) table_metAnno <- table_metAnno[-idx_rm, ]

    # remove duplicated f.-m. pairs
    comb_f_m_metAnno <- paste(table_metAnno$feature_name, table_metAnno$id, sep = '--')
    idx_keep_metAnno <- !duplicated(comb_f_m_metAnno)
    table_metAnno_uni <- table_metAnno[idx_keep_metAnno, ]
    comb_f_m_metAnno_uni <- comb_f_m_metAnno[idx_keep_metAnno]

    # mark to final result: table_anno
    idx_to_mark <- match(comb_f_m_metAnno_uni, comb_f_m_all)
    table_anno[idx_to_mark, ] <- table_metAnno_uni

    return(table_anno)

}
