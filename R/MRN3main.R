
# MrnAnnoAlgo ---------------------------------------------------------------------------------

#' MrnAnnoAlgo
#'
#' @param wd Working directory (default is the current directory).
#' @param ms1_data MS1 data in MRN3 format with columns 'name', 'mz', 'rt', and 'ccs'
#' @param ms2_data A list containing MS2 spectra data in MRN3 format (feature name matched with MS1).
#' @param table_seed_long A tibble containing the long-format seed table for MRN3 annotation with columns 'feature_name' and 'id_mrn'.
#' @param instrument_type instrument_type: "Orbitrap" or "TOF"
#' @param column LC type: "hilic" or "rp"
#' @param method_lc LC mwethod: 'Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min', 'RP12min_v2'
#' @param polarity polarity: "positive" or "negative"
#' @param thread thread (default is 4).
#' @param lib_type_mrn 'MRN3', 'MRN2'
#' @param ex_step extension step: '0', '1', '2'
#' @param is_rt_pred_by_metdna2 TRUE!
#' @param adduct_for_MRN_anno NULL or limited to some types
#' @param ms2_sim_method dp or hybrid
#' @param num_candi number of candidates (topN) for redundancy removal
#' @param is_known_prioritization priority: level3.1 (knowns) > level3.2 (unknowns) for redundancy removal
#'
#' @return null
#' @export
#'
#' @examples
#' MrnAnnoAlgo(
#' wd = 'TO/YOUR/PATH/',
#' ms1_data = ms1_data,
#' ms2_data = ms2_data,
#' table_seed_long = table_seed_long,
#' thread = 8)
#'
MrnAnnoAlgo <- function(
        wd = '.',
        ms1_data,
        ms2_data,
        table_seed_long,
        # experimental settings
        instrument_type = c('Orbitrap', 'TOF'),
        column = c('hilic', 'rp'),
        method_lc = c('Amide12min', 'Amide23min', 'Other', 'MetlinRP', 'RP12min', 'RP12min_v2'),
        polarity = c('positive', 'negative'),
        # MRN annotation settings
        lib_type_mrn = c('MRN3', 'MRN2'), # MRN2 = KEGG_MRN_ex2 in MetDNA2_function
        ex_step = c('0', '1', '2'), # BioTransformer extension step for unknown
        is_rt_pred_by_metdna2 = TRUE, # RT prediction by MetDNA2_function
        adduct_for_MRN_anno = NULL, # adduct for MS1 annotation
        ms2_sim_method = c('dp', 'hybrid'), # DP or HSS
        # redundancy removal settings (for MetDNA3, no old MetDNA!)
        num_candi = 3, # number of candidates
        is_known_prioritization = T, # priority: level3.1 (knowns) > level3.2 (unknowns)
        # other settings
        thread = 4
        ) {

    wd_output <- file.path(wd, '02_result_MRN3')
    dir.create(wd_output, showWarnings = F, recursive = T)

    library(dplyr)
    library(igraph)

    cat('MRN3 annotation start.', '\n')
    cat('\n')

    # param. preparation ---------------------------------------------------------------------------
    cat('Prepare parameters.', '\n')
    cat('\n')

    # fixed parameters
    is_rt_valid <- TRUE
    is_ccs_valid <- FALSE
    mz_tol <- 15
    mz_tol_ms2 <- 25
    rt_tol_rel <- 30
    rt_tol_adductAnno <- 3
    ccs_tol_rel <- 4
    ms2_cutoff <- 0.5
    w_mz <- 0.25
    w_rt <- 0.25
    w_ccs <- 0
    w_ms2 <- 0.5

    # instrument setting
    switch(
        instrument_type,
        'Orbitrap' = {
            mz_ppm_thr <- 0
        },
        'TOF' = {
            mz_ppm_thr <- 400
        }
    )

    # LC and method setting
    column <- match.arg(column)
    method_lc <- match.arg(method_lc)
    ex_step <- match.arg(ex_step)

    # polarity & adduct setting
    switch(
        polarity,
        'positive' = {
            adduct_table_all <- MetDNA3::lib_adduct_nl$positive
            if (lib_type_mrn == 'MRN2') {
                adduct_table <- adduct_table_all %>% dplyr::filter(annotation == 'Yes')
            } else {
                adduct_table <- adduct_table_all %>% dplyr::filter(annotation == 'Yes')
                # adduct_table <- adduct_table_all %>%
                #     dplyr::filter(adduct %in% c("[M+H]+", "[M-H2O+H]+", "[M+Na]+", "[M+NH4]+"))
            }
        },
        'negative' = {
            adduct_table_all <- MetDNA3::lib_adduct_nl$negative
            if (lib_type_mrn == 'MRN2') {
                adduct_table <- adduct_table_all %>% dplyr::filter(annotation == 'Yes')
            } else {
                adduct_table <- adduct_table_all %>% dplyr::filter(annotation == 'Yes')
                # adduct_table <- adduct_table_all %>%
                #     dplyr::filter(adduct %in% c("[M-H]-", "[M-H2O-H]-", "[M+Na-2H]-", "[M+CH3COO]-"))
            }
        },
        stop("Error: Invalid variable value of 'polarity'!")
    )
    # adduct_for_MRN_anno setting
    if (is.null(adduct_for_MRN_anno)) {
        adduct_for_MRN_anno <- adduct_table$adduct
    } else {
        adduct_table <- adduct_table %>%
            dplyr::filter(adduct %in% adduct_for_MRN_anno)
    }

    # MS2 match setting
    ms2_sim_method <- match.arg(ms2_sim_method)
    matchParam <- MrnAnnoAlgo3::getMatchParam(
        mz_tol_ms2 = mz_tol_ms2,
        mz_ppm_thr = mz_ppm_thr,
        methodScore = ms2_sim_method
    )

    st <- Sys.time()

    # load MRN objects and MRN info ----------------------------------------------------------------
    cat('Load MRN objects and MRN info.', '\n')
    cat('\n')
    switch(
        lib_type_mrn,
        'MRN3' = {
            switch(
                ex_step,
                '0' = {
                    load(system.file("obj_mrn.rda", package = "MrnAnnoAlgo3")) # obj_mrn
                    load(system.file("info_mrn.rda", package = "MrnAnnoAlgo3")) # info_mrn
                    load(system.file("md_mrn.rda", package = "MrnAnnoAlgo3")) # md_mrn
                },
                '1' = {
                    load(system.file("obj_mrn_3x1.rda", package = "MrnAnnoAlgo3")) # obj_mrn
                    load(system.file("info_mrn_3x.rda", package = "MrnAnnoAlgo3")) # info_mrn
                    load(system.file("md_mrn_3x.rda", package = "MrnAnnoAlgo3")) # md_mrn
                    info_mrn <- info_mrn_3x
                    info_mrn <- info_mrn %>% filter(min_reaction_step <= 1)
                    md_mrn <- md_mrn_3x
                    md_mrn <- md_mrn[rownames(md_mrn) %in% info_mrn$id, ]
                    obj_mrn <- obj_mrn_3x1
                    rm(list = c('info_mrn_3x', 'md_mrn_3x', 'obj_mrn_3x1')); gc()
                    gc()
                },
                '2' = {
                    load(system.file("obj_mrn_3x2.rda", package = "MrnAnnoAlgo3")) # obj_mrn
                    load(system.file("info_mrn_3x.rda", package = "MrnAnnoAlgo3")) # info_mrn
                    load(system.file("md_mrn_3x.rda", package = "MrnAnnoAlgo3")) # md_mrn
                    info_mrn <- info_mrn_3x
                    md_mrn <- md_mrn_3x
                    obj_mrn <- obj_mrn_3x2
                    rm(list = c('info_mrn_3x', 'md_mrn_3x', 'obj_mrn_3x2')); gc()
                    gc()
                }
            )
            step_mrn <- 1
            id_type <- 'id_mrn'
        },
        'MRN2' = {
            info_mrn <- MetDNA3::cpd_emrn
            md_mrn <- MetDNA3::md_mrn_emrn$version2
            obj_mrn <- MetDNA3::reaction_pair_network$version2$step2 # default ex_step = 2
            step_mrn <- 3
            id_type <- 'id_kegg'
        }
    )


    # add predicted RT & CCS info
    cat('Add predicted RT & CCS info.', '\n')
    cat('\n')
    if (is_rt_pred_by_metdna2) {
        cat('Predict RT by MetDNA2.', '\n')
        cat('\n')
        rt_result_path <- file.path(wd, "02_result_MRN_annotation", "00_intermediate_data")
        dir.create(rt_result_path, showWarnings = F, recursive = T)
        # MetDNA2 RT prediction method
        if (all(c('rt_result') %in% list.files(rt_result_path))) {
            load(file.path(rt_result_path, 'rt_result'))
        } else {
            inHouse_compound = as.data.frame(MetLib::metinfo) %>%
                dplyr::filter(lib_type %in% c('BioLib', 'zhulib', 'NIST_2017'))
            inHouse_compound$id_kegg <- inHouse_compound$id_mrn
            ms1_data_metdna2 <- MetDNA3::readInitialAnnotation(
                data = "ms2_match_annotation_result.csv",
                direction = 'reverse',
                rt_filter = FALSE,
                inHouse_compound = inHouse_compound,
                instrument = 'ThermoExploris',
                path = file.path(wd, "01_result_initial_seed_annotation")
            )
            rt_result <- MetDNA3::predictRT(
                data = ms1_data_metdna2,
                prefer_adduct = ifelse(polarity == "positive", "M+H", "M-H"),
                threads = thread,
                md_inHouse_cpd = MetDNA3::md_zhumetlib,
                md_kegg_cpd = md_mrn,
                use_default_md = TRUE,
                column = column,
                method_lc = method_lc,
                metdna_version = 'version2',
                use_pretrained_model = FALSE
            )
            save(
                rt_result,
                file = file.path(rt_result_path, 'rt_result'),
                compress = "xz", version = 2
            )
        }
        info_mrn$rt <- rt_result[[2]]$RT[match(info_mrn$id, rownames(rt_result[[2]]))]
    } else {
        rt_result_path <- file.path(wd, "02_result_MRN_annotation", "00_intermediate_data")
        dir.create(rt_result_path, showWarnings = F, recursive = T)
        if (all(c('rt_result') %in% list.files(rt_result_path))) {
            load(file.path(rt_result_path, 'rt_result'))
            info_mrn$rt <- rt_result$KEGG.rt$RT
        } else {
            # MRN3 prediction (simplified)
            info_mrn$rt <- MrnAnnoAlgo3::predictRtInMrn3(
                wd = wd,
                ms1_data = ms1_data,
                table_seed_long = table_seed_long,
                md_mrn = md_mrn,
                column = column
            )
        }
    }
    # info_mrn$rt <- -1
    info_mrn$ccs <- -1
    # reassign RT with seed exp RT
    table_seed_long$rt <- ms1_data$rt[match(table_seed_long$feature_name, ms1_data$name)]
    table_seed_exp_rt <- table_seed_long %>%
        dplyr::group_by(id_mrn) %>%
        summarise(rt = median(rt))
    table_seed_exp_rt <- table_seed_exp_rt[table_seed_exp_rt$id_mrn %in% info_mrn$id, ]
    plot(
        x = info_mrn$rt[match(table_seed_exp_rt$id_mrn, info_mrn$id)],
        y = table_seed_exp_rt$rt
    )
    info_mrn$rt[match(table_seed_exp_rt$id_mrn, info_mrn$id)] <- table_seed_exp_rt$rt

    # get seed annotation with MetLib --------------------------------------------------------------
    lib_levels = c('zhulib', 'BioLib', 'AllCCS_set2', 'NIST_2017', 'fiehnHilicLib')
    metinfo <- MetLib::metinfo %>%
        mutate(lib_type = factor(lib_type, lib_levels)) %>%
        arrange(lib_type)
    idx_match_f <- match(table_seed_long$feature_name, ms1_data$name)
    idx_match_m <- sapply(table_seed_long$id_mrn, function(x) {
        temp_idx <- grep(x, unlist(metinfo[, id_type]))
        if (length(temp_idx) == 0) return(NA)
        return(temp_idx[1])
    })
    seed_anno <- dplyr::add_row(
        MrnAnnoAlgo3::getEmptyAnno(),
        'feature_name' = table_seed_long$feature_name,
        'feature_mz' = ms1_data$mz[idx_match_f],
        'feature_rt' = ms1_data$rt[idx_match_f],
        'feature_ccs' = ms1_data$ccs[idx_match_f],
        'id' = table_seed_long$id_mrn,
        'name' = metinfo$name[idx_match_m],
        'formula' = metinfo$formula[idx_match_m],
        'em' = metinfo$monoisotopic_mass[idx_match_m]
    )
    if (length(which(is.na(seed_anno$name))) > 0) seed_anno <- seed_anno[-which(is.na(seed_anno$name)), ]
    # assign seed anno.: adduct, mz, rt, ccs, tag, level, seed_round
    temp_delta_mz <- (seed_anno$feature_mz - seed_anno$em)
    closest_index <- sapply(temp_delta_mz, function(x) which.min(abs(adduct_table_all$delta_mz - x)))
    seed_anno$adduct <- adduct_table_all$adduct[closest_index]
    seed_anno$mz <- seed_anno$em + adduct_table_all$delta_mz[closest_index]
    idx_m <- match(seed_anno$id, info_mrn$id)
    seed_anno$rt <- info_mrn$rt[idx_m]
    seed_anno$ccs <- info_mrn$ccs[idx_m]
    if (sum(is.na(idx_m)) > 0) {
        seed_anno$rt[which(is.na(idx_m))] <- 9999
        seed_anno$ccs[which(is.na(idx_m))] <- -1
    }
    seed_anno$tag = 'seed'
    seed_anno$level <- 0
    seed_anno$seed_round = 1
    # calculate m/z, rt, ccs error and score
    seed_anno <- MrnAnnoAlgo3::getErrorAndScore(
        table_anno = seed_anno,
        is_rt_valid = is_rt_valid,
        is_ccs_valid = is_ccs_valid,
        mz_ppm_thr = mz_ppm_thr,
        mz_tol = mz_tol,
        rt_tol_rel = rt_tol_rel,
        ccs_tol_rel = ccs_tol_rel
    )
    idx_fix <- which(seed_anno$rt_score < 0)
    if (length(idx_fix) > 0) seed_anno$rt_score[idx_fix] <- 0
    seed_anno$ms2_score <- 1
    readr::write_csv(seed_anno, file.path(wd_output, 'table_seed_anno.csv'), na = '')

    # simplify MRN object & MS1 data ---------------------------------------------------------------
    cat('Simplify MRN object & MS1 data.', '\n')
    cat('\n')

    # simply info_mrn !
    info_mrn <- info_mrn[info_mrn$id %in% igraph::V(obj_mrn)$name, ]
    # create global potential annotation results
    if (all(c('table_ms1_anno.csv') %in% list.files(wd_output))) {
        table_ms1_anno <- readr::read_csv(file.path(wd_output, 'table_ms1_anno.csv'))
    } else {
        # search by adducts
        library(parallel)
        cl <- parallel::makeCluster(min(thread, length(adduct_table$adduct)))
        clusterExport(cl, c(
            # 'adduct_table',
            # 'info_mrn',
            # 'mz_ppm_thr',
            # 'mz_tol',
            # 'ms1_data',
            '%>%'
        ), envir = environment())
        tempFunMs1Anno <- function(i) {
            adduct <- adduct_table$adduct[i]
            cat('Match adduct:', adduct, '\n')
            delta_mz <- adduct_table$delta_mz[i]
            temp_M <- regmatches(adduct, regexec("\\[(\\d+)", adduct))[[1]]
            if (length(temp_M) > 1) {temp_M <- as.numeric(temp_M[2])} else {temp_M <- 1}
            mrn_adduct_mz <- info_mrn$monoisotopic_mass * temp_M + delta_mz
            mrn_adduct_mz_ranges <- MetDNA3:::getMzRange(mz = mrn_adduct_mz, ppm = mz_tol, mz_ppm_thr = mz_ppm_thr)
            # match ms1 data and MRN: remove unmatched features & metabolites
            # [sub_table_ms1_anno]
            sub_table_ms1_anno <- lapply(seq_along(ms1_data$mz), function(i) {
                idx_mrn <- which(ms1_data$mz[i] >= mrn_adduct_mz_ranges[, 1] &
                                     ms1_data$mz[i] <= mrn_adduct_mz_ranges[, 2])
                if (length(idx_mrn) > 0) {
                    idx_features <- rep(i, length(idx_mrn))
                    new_anno <- dplyr::add_row(
                        MrnAnnoAlgo3::getEmptyAnno(),
                        'feature_name' = ms1_data$name[idx_features],
                        'feature_mz' = ms1_data$mz[idx_features],
                        'feature_rt' = ms1_data$rt[idx_features],
                        'feature_ccs' = ms1_data$ccs[idx_features],
                        'id' = info_mrn$id[idx_mrn],
                        'name' = info_mrn$name[idx_mrn],
                        'formula' = info_mrn$formula[idx_mrn],
                        'em' = info_mrn$monoisotopic_mass[idx_mrn],
                        'adduct' = adduct,
                        'mz' = info_mrn$monoisotopic_mass[idx_mrn] * temp_M + delta_mz,
                        'rt' = info_mrn$rt[idx_mrn],
                        'ccs' = info_mrn$ccs[idx_mrn]
                    )
                    return(new_anno)
                } else return(MrnAnnoAlgo3::getEmptyAnno())
            }) %>% dplyr::bind_rows()
            return(sub_table_ms1_anno)
        }

        system.time(list_ms1_anno_adducts <- parallel::parLapply(cl, seq_along(adduct_table$adduct), tempFunMs1Anno))
        stopCluster(cl)
        table_ms1_anno_adducts <- dplyr::bind_rows(list_ms1_anno_adducts)

        # add potential annotations to table_ms1_anno
        table_ms1_anno <- table_ms1_anno_adducts

        # calculate m/z, rt, ccs error and score
        table_ms1_anno <- MrnAnnoAlgo3::getErrorAndScore(
            table_anno = table_ms1_anno,
            is_rt_valid = is_rt_valid,
            is_ccs_valid = is_ccs_valid,
            mz_ppm_thr = mz_ppm_thr,
            mz_tol = mz_tol,
            rt_tol_rel = rt_tol_rel,
            ccs_tol_rel = ccs_tol_rel
        )

        readr::write_csv(table_ms1_anno, file.path(wd_output, 'table_ms1_anno.csv'), na = '')
    }

    # get simplified annotation table of features (with MS2)
    # MS2 index preparation
    ms2_name <- unname(sapply(ms2_data, function(x) x$info[,1]))

    # RT, CCS filter ?!
    # mark seed anno score = 1!
    comb_seed <- paste(table_seed_long$feature_name, table_seed_long$id_mrn, sep = '--')
    comb_f_m_all <- paste(table_ms1_anno$feature_name, table_ms1_anno$id, sep = '--')
    idx_seed_mark <- which(comb_f_m_all %in% comb_seed)
    idx_fix <- which(table_ms1_anno$rt_score[idx_seed_mark] < 0)
    if (length(idx_fix) > 0) table_ms1_anno$rt_score[idx_seed_mark[idx_fix]] <- 0 # level 2 seed no RT
    # # ignore seed "ID" metabolite RT error (< version 0.2)
    # ids_seed <- sort(unique(table_seed_long$id_mrn))
    # idx_rt_remark <- which(table_ms1_anno$id %in% ids_seed)
    # idx_fix <- which(table_ms1_anno$rt_score[idx_rt_remark] < 0)
    # table_ms1_anno$rt_score[idx_rt_remark[idx_fix]] <- 0
    # filter
    if (is_rt_valid) {
        table_ms1_anno <- table_ms1_anno[table_ms1_anno$rt_score >= 0, ]
    }
    if (is_ccs_valid) {
        table_ms1_anno <- table_ms1_anno[table_ms1_anno$ccs_score >= 0, ]
    }

    # features for MRN annotation:
    # 1. with MS2
    # 2. in adducts for MRN annotation
    idx_for_anno <- which(
        table_ms1_anno$feature_name %in% ms2_name &
            table_ms1_anno$adduct %in% adduct_for_MRN_anno
    )
    # table_ms1_anno_with_ms2 <- table_ms1_anno[table_ms1_anno$feature_name %in% ms2_name, ]
    table_ms1_anno_with_ms2 <- table_ms1_anno[idx_for_anno, ]
    # table_ms1_anno_only <- table_ms1_anno[!(table_ms1_anno$feature_name %in% ms2_name), ]
    table_ms1_anno_only <- table_ms1_anno[-idx_for_anno, ]
    readr::write_csv(table_ms1_anno_with_ms2, file.path(wd_output, 'table_ms1_anno_with_ms2.csv'), na = '')
    readr::write_csv(table_ms1_anno_only, file.path(wd_output, 'table_ms1_anno_only.csv'), na = '')

    cat('Putative annotations:', nrow(table_ms1_anno), '\n')
    cat('Putative annotations (features with MS2, for MRN annotation):', nrow(table_ms1_anno_with_ms2), '\n')

    # get simplified MRN
    # After MRN2 is simplified, it will reduce the propagation efficiency and will not be implemented.
    if (lib_type_mrn == 'MRN3') {
        simplified_nodes <- table_ms1_anno$id %>% unique() %>% sort()
        simplified_nodes <- simplified_nodes[simplified_nodes %in% V(obj_mrn)$name]
        obj_mrn <- induced_subgraph(obj_mrn, simplified_nodes)
    }
    cat('Simplified MRN metabolites:', length(obj_mrn), '\n')
    cat('Simplified MRN pairs:', nrow(igraph::get.edgelist(obj_mrn)), '\n')
    # cat('MRN metabolites:', length(obj_mrn), '\n')
    # cat('MRN pairs:', nrow(igraph::get.edgelist(obj_mrn)), '\n')
    cat('\n')

    # get links of features with MS2 (MS2 network edges) -------------------------------------------
    if ('table_ms2_edges.rda' %in% list.files(wd_output)) {
        cat('Load MS2 network edges (feature MS2 pairs).', '\n')
        cat('\n')
        load(file.path(wd_output, 'table_ms2_edges.rda'))
    } else {
        cat('Get MS2 network edges (feature MS2 pairs).', '\n')

        library(parallel)
        cl <- parallel::makeCluster(thread)
        clusterExport(cl, c(
            # 'table_ms1_anno_with_ms2',
            # 'obj_mrn',
            # 'step_mrn',
            '%>%'
        ), envir = environment())

        getMs2NameLink <- function(feature_from) {

            table_from <- table_ms1_anno_with_ms2[table_ms1_anno_with_ms2$feature_name == feature_from, ]
            if (nrow(table_from) == 0) return(NULL)
            id_from <- table_from$id
            id_from <- id_from[id_from %in% igraph::V(obj_mrn)$name]
            if (length(id_from) == 0) return(NULL)

            node_to <- MrnAnnoAlgo3::getMrnNeighbor(id = id_from, step = step_mrn, graph = obj_mrn)
            # use all-step ids
            node_to <- node_to$id

            # # version 0.2.14.2.4: FCMN (Fully-connected metabolite network) # neighbors = all nodes temp!
            # node_to <- sprintf("MRN%06d", 1:53583)

            if (length(node_to) == 0) return(NULL)
            table_to <- table_ms1_anno_with_ms2[table_ms1_anno_with_ms2$id %in% node_to, ]
            feature_to <- unique(table_to$feature_name)
            feature_to <- feature_to[feature_to != feature_from]
            if (length(feature_to) == 0) return(NULL)
            return(data.frame(from = feature_from, to = feature_to))

        }
        system.time(list_ms2_name_link <- parallel::parLapply(cl, ms2_name, getMs2NameLink))
        stopCluster(cl)
        table_ms2_edges <- dplyr::bind_rows(list_ms2_name_link) %>% tibble::as_tibble()

        comb_from_to <- apply(table_ms2_edges, 1, function(row)
            stringr::str_c(sort(c(row[1], row[2])), collapse = '--'))
        table_ms2_edges <- table_ms2_edges[which(!duplicated(comb_from_to)),]
        cat('MS2 network edges:', nrow(table_ms2_edges), '\n')
        cat('\n')
        rm(list = c('list_ms2_name_link'))
        save(table_ms2_edges, file = file.path(wd_output, 'table_ms2_edges.rda'))
    }

    # calculate similarity of MS2 network edges (MS2 pairs) ----------------------------------------
    if ('table_ms2_network.rda' %in% list.files(wd_output)) {
        cat('Load similarity of MS2 network edges  (feature MS2 pairs).', '\n')
        cat('\n')
        load(file.path(wd_output, 'table_ms2_network.rda'))
    } else {
        cat('Calculate similarity of MS2 network edges (feature MS2 pairs).', '\n')
        cat('\n')
        library(parallel)
        cl <- parallel::makeCluster(thread)
        clusterExport(cl, c(
            'table_ms2_edges',
            'ms2_name',
            'ms2_data',
            'matchParam',
            '%>%'
        ), envir = environment())
        tempFunMs2 <- function(x) {
            feature_from <- table_ms2_edges$from[x]
            feature_to <- table_ms2_edges$to[x]
            idx_from <- match(feature_from, ms2_name)
            idx_to <- match(feature_to, ms2_name)
            # the lib spec always use spec from the peak with smaller mz (MetDNA2)
            mz_from <- as.data.frame(ms2_data[[idx_from]]$info)$PRECURSORMZ %>% as.numeric()
            mz_to <- as.data.frame(ms2_data[[idx_to]]$info)$PRECURSORMZ %>% as.numeric()
            if (mz_to <= mz_from) {
                match_result <- SpectraTools::MatchSpectra(
                    dataExp = MrnAnnoAlgo3::convertSpectraData(ms2_data = ms2_data[[idx_from]]),
                    dataRef = MrnAnnoAlgo3::convertSpectraData(ms2_data = ms2_data[[idx_to]]),
                    matchParam)
            } else {
                match_result <- SpectraTools::MatchSpectra(
                    dataExp = MrnAnnoAlgo3::convertSpectraData(ms2_data = ms2_data[[idx_to]]),
                    dataRef = MrnAnnoAlgo3::convertSpectraData(ms2_data = ms2_data[[idx_from]]),
                    matchParam)
            }
            if (!is.null(match_result)) {
                match_info <- as.data.frame(match_result@info)
                if (length(match_result@matchedFragments) > 0) {
                    temp_matchedFragments <- match_result@matchedFragments[[1]]
                    if (is.null(temp_matchedFragments)) {
                        n_frag_match <- 0
                    } else {
                        n_frag_match <- temp_matchedFragments %>%
                            dplyr::filter(intensity > 0 & intensityExp > 0) %>%
                            dplyr::count() %>% dplyr::pull()
                    }
                } else n_frag_match <- 0
                if (length(match_result@nlFragments) > 0) {
                    n_nl_match <- match_result@nlFragments[[1]] %>%
                        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
                        dplyr::count() %>% dplyr::pull()
                } else n_nl_match <- 0
                match_info$matchFrag <- n_frag_match
                match_info$matchNl <- n_nl_match
            } else {
                match_info <- data.frame(name = NA)
            }
            #### msentropy
            # ms2_entropy_score <- msentropy::calculate_entropy_similarity(
            #     peaks_a = ms2_data[[idx_from]]$spec,
            #     peaks_b = ms2_data[[idx_to]]$spec,
            #     ms2_tolerance_in_da = ms2_tolerance_in_da,
            #     ms2_tolerance_in_ppm = mz_tol_ms2,
            #     clean_spectra = T,
            #     min_mz = -1,
            #     max_mz = -1,
            #     noise_threshold = 0.01,
            #     max_peak_num = -1
            # )
            # match_info <- data.frame(
            #     ms2_score = ms2_entropy_score,
            #     matchFrag = NA,
            #     matchNl = NA
            # )
            return(match_info)
        }
        system.time(ms2_scores <- parallel::parLapply(cl, 1:nrow(table_ms2_edges), tempFunMs2))
        stopCluster(cl)
        ms2_scores <- dplyr::bind_rows(ms2_scores)
        table_ms2_network <- dplyr::bind_cols(table_ms2_edges, ms2_scores)
        if ('scoreReverse' %in% colnames(table_ms2_network)) table_ms2_network <- table_ms2_network %>%
            dplyr::rename(ms2_score = scoreReverse) %>% dplyr::select(-name, -mz)
        if ('score' %in% colnames(table_ms2_network)) table_ms2_network <- table_ms2_network %>%
            dplyr::rename(ms2_score = score) %>% dplyr::select(-name, -mz)
        save(table_ms2_network, file = file.path(wd_output, 'table_ms2_network.rda'))
        rm(list = c('ms2_scores')); gc()
    }

    table_ms2_network_filter <- table_ms2_network %>%
        # filter(scoreReverse >= ms2_cutoff) # 20240131 before
        # filter(scoreReverse >= ms2_cutoff | matchFrag >= 4) # 20240401 before
        filter(ms2_score >= ms2_cutoff)
    save(table_ms2_network_filter, file = file.path(wd_output, 'table_ms2_network_filter.rda'))
    rm(list = c('table_ms2_edges')); gc()

    # recursive annotation with seeds --------------------------------------------------------------
    if ('table_mrn_anno.csv' %in% list.files(wd_output)) {
        cat('Load recursive annotation result.', '\n')
        cat('\n')
        table_mrn_anno <- readr::read_csv(file.path(wd_output, 'table_mrn_anno.csv'))
    } else {
        cat('Do recursive annotation.', '\n')
        cat('\n')
        # get MRN annotation table (seed anno + table_ms1_anno_with_ms2)
        table_mrn_anno <- table_ms1_anno_with_ms2
        seed_f_m_pair <- paste0(seed_anno$feature_name, '--', seed_anno$id)
        potential_f_m_pair <- paste0(table_mrn_anno$feature_name, '--', table_mrn_anno$id)
        idx_in_seed <- which(potential_f_m_pair %in% seed_f_m_pair)
        if (length(idx_in_seed) > 0) table_mrn_anno <- table_mrn_anno[-idx_in_seed, ]
        table_mrn_anno <- dplyr::bind_rows(seed_anno, table_mrn_anno) %>% dplyr::arrange(id)

        # start recursive annotation by seeds
        round <- 1
        idx_seeds <- which(table_mrn_anno$seed_round == round)
        # comb_f_m_all <- paste(table_mrn_anno$feature_name, table_mrn_anno$id, sep = '--')
        # comb_f_m_annotated <- comb_f_m_all[which(table_mrn_anno$tag == 'seed')]

        while (length(idx_seeds) > 0) {

            # summary before recursive annotation
            len_seed_feature <- table_mrn_anno$feature_name[idx_seeds] %>% unique() %>% length()
            cat('\n')
            cat('Round:', round, '\n')
            cat('Seed feature number:', len_seed_feature, '\n')
            cat('Seed metabolite number:', length(idx_seeds), '\n')
            cat('\n')

            cat('Isotope & Adduct Annotation', '\n')
            cat('\n')

            # skip isotopc annotation (meaningless)
            # input 'table_mrn_anno' and 'idx_seeds' for adduct search
            table_mrn_anno <- getAdductAnno(
                table_anno = table_mrn_anno,
                idx_seeds = idx_seeds,
                round = round,
                rt_tol_adductAnno = rt_tol_adductAnno
            )

            cat('Metabolite Annotation', '\n')
            cat('\n')

            table_mrn_anno <- getMetAnno(
                table_anno = table_mrn_anno,
                idx_seeds = idx_seeds,
                round = round,
                obj_mrn = obj_mrn,
                table_ms2_network_filter = table_ms2_network_filter,
                step_mrn = step_mrn,
                lib_type_mrn = lib_type_mrn,
                thread = thread
            )

            # mark new seeds
            # 1. tag: Annotated
            # 2. seed_round: NA
            round <- round + 1
            idx_seeds_new <- which(!is.na(table_mrn_anno$tag) & is.na(table_mrn_anno$seed_round))
            table_mrn_anno$seed_round[idx_seeds_new] <- round

            idx_seeds <- idx_seeds_new

        }

        readr::write_csv(table_mrn_anno, file.path(wd_output, 'table_mrn_anno.csv'), na = '')

    }

    # summary
    result_mrn_anno <- table_mrn_anno[!is.na(table_mrn_anno$tag), ]
    result_mrn_anno$total_score <-
        result_mrn_anno$mz_score * w_mz +
        result_mrn_anno$rt_score * w_rt +
        result_mrn_anno$ccs_score * w_ccs +
        result_mrn_anno$ms2_score * w_ms2
    # seed total score is 1
    result_mrn_anno$total_score[which(result_mrn_anno$tag == 'seed')] <- 1

    # add adduct annotation result in table_ms1_anno_only (same with MetDNA) -----------------------
    # feature without MS2 should be annotated again
    # method modified by 'getAdductAnno' function
    # filter
    # 1. same id
    # 2. RT diff. <= rt_tol_adductAnno
    table_ms1_anno_only <- readr::read_csv(file.path(wd_output, 'table_ms1_anno_only.csv'))
    ids_annotated <- sort(unique(result_mrn_anno$id))
    max_round <- max(result_mrn_anno$seed_round)
    list_adduct_anno_only <- lapply(ids_annotated, function(temp_id) {
        temp_table_ms1_anno_only <- table_ms1_anno_only[which(table_ms1_anno_only$id == temp_id), ]
        # RT check
        temp_mrn_anno <- result_mrn_anno[which(result_mrn_anno$id == temp_id), ]
        rt_ranges <- data.frame(
            minus = temp_mrn_anno$feature_rt - rt_tol_adductAnno,
            plus = temp_mrn_anno$feature_rt + rt_tol_adductAnno
        )
        idx_rt_check <- sapply(
            temp_table_ms1_anno_only$feature_rt,
            function(x) which(x >= rt_ranges[, 1] & x <= rt_ranges[, 2])[1]
        )
        idx_na_check <- which(!is.na(idx_rt_check))
        if (length(idx_na_check) == 0) return(NULL)
        temp_table_ms1_anno_only <- temp_table_ms1_anno_only[idx_na_check, ]
        idx_rt_check <- idx_rt_check[idx_na_check]
        temp_table_ms1_anno_only$tag <- 'adductAnnotation'
        temp_table_ms1_anno_only$level <- max_round
        temp_table_ms1_anno_only$feature_from <- temp_mrn_anno$feature_name[idx_rt_check]
        temp_table_ms1_anno_only$id_from <- temp_id
        temp_table_ms1_anno_only$total_score <- temp_mrn_anno$total_score[idx_rt_check]
        return(temp_table_ms1_anno_only)
    })
    table_adduct_anno_only <- list_adduct_anno_only %>% dplyr::bind_rows()
    # readr::write_csv(table_mrn_anno, file.path(wd_output, 'table_mrn_anno_before_ms1_adduct_anno.csv'))
    result_mrn_anno <- dplyr::bind_rows(result_mrn_anno, table_adduct_anno_only)
    result_mrn_anno <- result_mrn_anno %>% dplyr::arrange(id)
    readr::write_csv(table_adduct_anno_only, file.path(wd_output, 'table_adduct_anno_only.csv'), na = '')
    readr::write_csv(result_mrn_anno, file.path(wd_output, 'result_mrn_anno_raw.csv'), na = '')

    # redundancy removal (MetDNA method) -----------------------------------------------------------
    # result_mrn_anno_with_confidence <- MrnAnnoAlgo3::rmRedunInMrn3(result_df = result_mrn_anno, is_only_confidence_assignment = T)
    result_mrn_anno_rm_redun_by_metdna <- MrnAnnoAlgo3::rmRedunInMrn3(result_df = result_mrn_anno)
    readr::write_csv(result_mrn_anno_rm_redun_by_metdna, file.path(wd_output, 'result_mrn_anno_rm_redun_by_old_metdna.csv'), na = '')

    # redundancy removal ---------------------------------------------------------------------------
    result_mrn_anno_rm_redun_by_metdna3 <- MrnAnnoAlgo3::rmRedunInMetdna3(
        result_df = result_mrn_anno,
        num_candi = num_candi, # number of candidates
        is_known_prioritization = is_known_prioritization, # priority: level3.1 (knowns) > level3.2 (unknowns)
        wd_output = wd_output
    )

    readr::write_csv(result_mrn_anno_rm_redun_by_metdna3, file.path(wd_output, 'result_mrn_anno_with_confidence.csv'), na = '')

    et <- Sys.time()
    et - st

    cat('MRN3 annotation finished.', '\n')
    cat('\n')

}
