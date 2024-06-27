# getEmptyAnno --------------------------------------------------------------------------------

#' getEmptyAnno
#'
#' @return A tibble of 'empty_anno' (for MRN3)
#' @export
#'
#' @examples getEmptyAnno()
getEmptyAnno <- function() {
    empty_anno <- tibble::tibble(
        'feature_name' = character(),
        'feature_mz' = numeric(),
        'feature_rt' = numeric(),
        'feature_ccs' = numeric(),
        'id' = character(),
        'name' = character(),
        'formula' = character(),
        'em' = numeric(),
        'adduct' = character(),
        'mz' = numeric(),
        'rt' = numeric(),
        'ccs' = numeric(),
        'tag' = character(),
        'level' = numeric(),
        'seed_round' = numeric(),
        'feature_from' = character(),
        'id_from' = character(),
        'matched_frag' = numeric(),
        'matched_nl' = numeric(),
        'mz_error' = numeric(),
        'rt_error_abs' = numeric(),
        'rt_error_rel' = numeric(),
        'ccs_error' = numeric(),
        'mz_score' = numeric(),
        'rt_score' = numeric(),
        'ccs_score' = numeric(),
        'ms2_score' = numeric(),
        'total_score' = numeric()
    )
    return(empty_anno)
}

# getMatchParam -------------------------------------------------------------------------------

#' getMatchParam
#'
#' @param mz_tol_ms2 mz_tol_ms2 25 ppm
#' @param mz_ppm_thr mz_ppm_thr: 400 Da in TOF, 0 in Orbitrap
#' @param includePrecursor FALSE in MRN annotation
#' @param ppmPrecursorFilter ppmPrecursorFilter 30 ppm
#' @param thrIntensityAbs thrIntensityAbs 50
#' @param thrIntensityRel thrIntensityRel 0.01 = 1 pctg.
#' @param methodScore methodScore: 'dp' or 'hybrid'
#'
#' @return SpectraTools::MatchParam
#' @export
#'
#' @examples getMatchParam(mz_tol_ms2 = 25, mz_ppm_thr = 0)
getMatchParam <- function(
        mz_tol_ms2 = 25,
        mz_ppm_thr = 400,
        includePrecursor = FALSE,
        ppmPrecursorFilter = 30,
        thrIntensityAbs = 50,
        thrIntensityRel = 0.01,
        methodScore = c("dp", "hybrid")
    ) {
    methodScore <- match.arg(methodScore)
    matchParam <- SpectraTools::MatchParam(
        ppm = mz_tol_ms2,
        resDefineAt = mz_ppm_thr,
        cutoff = 0,
        weightIntensity = 1,
        weightMZ = 0,
        normIntensity = TRUE,
        tuneLibSpectra = TRUE,
        intensityExpNormed = TRUE,
        intensityLibNormed = TRUE,
        includePrecursor = includePrecursor,
        ppmPrecursorFilter = ppmPrecursorFilter,
        thrIntensityAbs = thrIntensityAbs,
        thrIntensityRel = thrIntensityRel,
        intensityNormedMethod = 'maximum',
        methodMatch = 'direct',
        methodScore = methodScore
    ) %>% new(Class = 'MatchParam')
    return(matchParam)
}

# convertSpectraData --------------------------------------------------------------------------

#' convertSpectraData
#'
#' @param ms2_data ms2_data
#'
#' @return SpectraData
#' @export
#'
#' @examples convertSpectraData(ms2_data)
convertSpectraData <- function(ms2_data) {
    options(readr.num_columns = 0)
    temp_info <- ms2_data$info %>%
        as.data.frame() %>%
        dplyr::rename(
            name = NAME,
            mz = PRECURSORMZ) %>%
        dplyr::select(name:mz) %>%
        readr::type_convert()
    temp_ms2_data <- ms2_data$spec
    result <- new(
        'SpectraData',
        info = temp_info,
        spectra = list(temp_ms2_data))
    return(result)
}


# getErrorAndScore ----------------------------------------------------------------------------

#' getErrorAndScore
#'
#' @param table_anno table_anno
#' @param is_rt_valid is_rt_valid (TRUE)
#' @param is_ccs_valid is_ccs_valid (FALSE)
#' @param mz_ppm_thr mz_ppm_thr (0 Da for Orbitrap; 400 Da for TOF)
#' @param mz_tol mz_tol (MS1, 15 ppm)
#' @param rt_tol_rel rt_tol_rel (30 pctg.)
#' @param ccs_tol_rel ccs_tol_rel (4 pctg.)
#'
#' @return table_anno
#' @export
#'
#' @examples getErrorAndScore(table_anno)
getErrorAndScore <- function(
        table_anno,
        is_rt_valid = TRUE,
        is_ccs_valid = FALSE,
        mz_ppm_thr = 0,
        mz_tol = 15,
        rt_tol_rel = 30,
        ccs_tol_rel = 4
    ) {
    table_anno$mz_error <- round(
        abs(table_anno$feature_mz - table_anno$mz) / max(table_anno$mz, mz_ppm_thr) * 1e6, 2)
    if (is_rt_valid) {
        table_anno$rt_error_abs <- round(abs(table_anno$feature_rt - table_anno$rt), 2)
        table_anno$rt_error_rel <- (abs(table_anno$feature_rt - table_anno$rt) / table_anno$rt * 100) %>% round(., 2)
    } else {
        table_anno$rt_error_abs <- -1
        table_anno$rt_error_rel <- -1
    }
    if (is_ccs_valid) {
        table_anno$ccs_error <- (abs(table_anno$feature_ccs - table_anno$ccs) / table_anno$ccs * 100) %>% round(., 2)
    } else table_anno$ccs_error <- -1
    # m/z, rt, ccs score
    table_anno$mz_score <- round(1 - table_anno$mz_error / mz_tol, 4)
    table_anno$rt_score <- round(1 - table_anno$rt_error_rel / rt_tol_rel, 4)
    if (is_ccs_valid) {
        table_anno$ccs_score <- round(1 - table_anno$ccs_error / ccs_tol_rel, 4)
    } else table_anno$ccs_score <- 0

    return(table_anno)
}

