
# assignConfidenceInMrn3 ----------------------------------------------------------------------

#' assignConfidenceInMrn3
#'
#' @param result_df annotation result data frame with column: "tag" & "adduct" ! (same peak group!)
#'
#' @return "confidence" of this peak group (a value vector)
#' @export
#'
#' @examples assignConfidenceInMrn3(result_df)
#'
assignConfidenceInMrn3 <- function(result_df) {
    temp_tag <- result_df$tag
    temp_adduct <- result_df$adduct
    ### seed -> [grade 1]
    if (any(temp_tag == 'seed')) return('grade1')
    ### no isotope (no grade 2)
    # temp.idx1 <- which(x$isotope %in% c("[M+1]","[M+2]","[M+3]","[M+4]"))
    # if (length(temp.idx1) > 0) return("grade2")
    ### prefer adduct -> [grade 3]
    prefer_adduct <- c(
        '[M+H]+', '[M+Na]+', '[M+NH4]+',
        '[M-H]-', '[M+CH3COO]-', '[M+Cl]-'
    ) # prefer_adduct same with MetDNA
    temp.idx2 <- which(temp_adduct %in% prefer_adduct)
    if (length(temp.idx2) > 0) return('grade3')
    ### other -> [grade 4]
    return('grade4')
}

# calculateRedundancyInMrn3 -------------------------------------------------------------------

#' calculateRedundancyInMrn3
#'
#' @param result_df annotation result data frame (all annotation results!)
#'
#' @return "Peak Redundancy" & "Metabolite Redundancy" (two value vector)
#' @export
#'
#' @examples calculateRedundancyInMrn3(result_df)
#'
calculateRedundancyInMrn3 <- function(result_df) {
    ### calculate redundancy
    # Peak redundancy: 1 feature -> x metabolites (peak number / annotated metabolite number)
    # Metabolite redundancy: 1 metabolite -> x peak group (metabolite peak group number / unique metabolite number)
    uni_f_name <- unique(result_df$feature_name)
    uni_id <- unique(result_df$id)
    redun_peak <- nrow(result_df) / length(uni_f_name) # 1 feature -> x metabolites
    redun_id <- length(unique(result_df$group)) / length(uni_id) # MetDNA method (modified)
    # group_num_by_id <- sapply(uni_id, function(x) max(result_df$group[result_df$id == x]))
    # redun_id <- mean(group_num_by_id) # 1 metabolite -> x peak group
    return(c(redun_peak, redun_id))
}

# rmRedunInMrn3 -------------------------------------------------------------------------------

#' rmRedunInMrn3
#'
#' @param result_df annotation result data frame (all annotation results!)
#' @param is_only_confidence_assignment FALSE
#'
#' @return annotation result data frame after redundancy removal (with column: "group" & "confidence" !)
#' @export
#'
#' @examples rmRedunInMrn3(result_df)
#'
rmRedunInMrn3 <- function(result_df, is_only_confidence_assignment = FALSE) {

    # RT tolerance of peak grouping (s)
    rt_tol_group <- 3

    # Peak grouping & confidence assignment & redundancy removal (MetDNA method) ----
    # 1. RT group
    # 2. annotation selection
    # 3. confidence calculation
    #
    # (NOTE IN METDNA)
    # # for each unique id:
    #    1. group all annotation into peak group: RT = 3 (with same annotation)
    #    2. if seed exist, only keep the 1st annotation records
    #    3. if not seed annotation, only keep the highest score records
    # (Note: the metAnnotation records with highest score is prior to reserve after version 0.6.62)
    #
    # f-m pair is unique in MRN3 design, so anno. result selection is not needed!

    # result_df_bak <- result_df
    # result_df <- result_df_bak

    # Peak grouping & confidence assignment ----
    cat('Peak grouping & confidence assignment.', '\n')
    cat('\n')
    ids_anno <- sort(unique(result_df$id))
    peaks_anno <- sort(unique(result_df$feature_name))
    result_df <- lapply(ids_anno, function(temp_id) {
        temp_res <- result_df[result_df$id == temp_id, ]
        rt_class <- MetDNA3::groupRT2(rt_vec = temp_res$feature_rt, rt_tol = rt_tol_group)
        # assign confidence: column 'tag' & 'adduct'
        temp_confidence <- unlist(lapply(rt_class, function(x) assignConfidenceInMrn3(temp_res[x, ])))
        names(rt_class) <- temp_confidence
        temp_confidence <- stringr::str_extract(names(sort(unlist(rt_class))), 'grade\\d')
        temp_group <- c()
        for (i in seq_along(rt_class)) {
            temp_group[rt_class[[i]]] <- i
        }
        temp_res$group <- stringr::str_c(temp_id, "_pg", temp_group)
        temp_res$confidence <- temp_confidence
        return(temp_res)
    }) %>% dplyr::bind_rows()

    # CHECK POINT! is_only_confidence_assignment ----
    if (is_only_confidence_assignment) return(result_df)

    cat('Redundancy removal.', '\n')
    cat('\n')

    redun_record <- list()
    round_redun <- 1

    while (round_redun <= 5) {

        ### calculate redundancy
        ids_anno <- sort(unique(result_df$id))
        peaks_anno <- sort(unique(result_df$feature_name))
        redun <- calculateRedundancyInMrn3(result_df)
        cat('Round:', round_redun, '\n')
        cat('Annotated peaks:', length(peaks_anno), '\n')
        cat('Annotated metabolites:', length(ids_anno), '\n')
        cat("Peak redundancy:", round(redun[1], 2), '\n')
        cat("Metabolite redundancy:", round(redun[2], 2), '\n')
        cat('\n')

        # CHECK POINT! delta_redun ----
        redun_record[[round_redun]] <- redun
        if (round_redun > 1) {
            delta_redun <- abs(mean(redun_record[[round_redun]]) - mean(redun_record[[round_redun - 1]]))
            if (delta_redun == 0) break
        }

        # remove redundancy (MetDNA method) ----

        ### remove [grade 4] annotations (by metabolite!)
        result_df_rm_redun_met <- lapply(ids_anno, function(temp_id) {
            temp_res <- result_df[result_df$id == temp_id, ]
            temp_confidence <- temp_res$confidence
            # (NOTE IN METDNA)
            # if the metabolite has more than [one] groups, and one group confidence is grade 4,
            # then the group with the grade 4 is removed.
            if (length(unique(temp_confidence)) > 1 & any(temp_confidence == "grade4")) {
                temp_res <- temp_res[-which(temp_confidence == "grade4"),]
                temp_confidence <- temp_confidence[-which(temp_confidence == "grade4")]
            }
            return(temp_res)
        }) %>% dplyr::bind_rows()

        ### remove peak [1 to n] annotations (by feature!)
        peaks_anno <- unique(result_df_rm_redun_met$feature_name)
        remain_idx <- lapply(peaks_anno, function(temp_peak){
            temp.idx <- which(temp_peak == result_df_rm_redun_met$feature_name)
            if (length(temp.idx) == 1) return(temp.idx)
            temp_res <- result_df_rm_redun_met[temp.idx, ]
            temp_confidence <- temp_res$confidence
            if (length(unique(temp_confidence)) > 1) {
                temp_grade <- as.numeric(substr(temp_confidence, 6, 6))
                remain_temp_idx <- temp.idx[which(temp_grade == temp_grade[which.min(temp_grade)])]
                return(remain_temp_idx)
            }
            return(temp.idx)
        }) %>% unlist()
        result_df_rm_redun_peak <- result_df_rm_redun_met[sort(remain_idx), ]

        ### assign confidence to result again
        groups_anno <- unique(result_df_rm_redun_peak$group)
        result_df_reassign <- lapply(groups_anno, function(temp_group) {
            temp_res <- result_df_rm_redun_peak[result_df_rm_redun_peak$group == temp_group, ]
            # assign confidence: column 'tag' & 'adduct'
            temp_confidence <- assignConfidenceInMrn3(temp_res)
            temp_res$confidence <- temp_confidence
            return(temp_res)
        }) %>% dplyr::bind_rows()

        ### remark
        result_df <- result_df_reassign
        round_redun <- round_redun + 1

    }

    return(result_df)
}
