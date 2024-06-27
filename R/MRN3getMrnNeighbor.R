
# getMrnNeighbor ------------------------------------------------------------------------------

#' getMrnNeighbor
#'
#' @param id A character string specifying the metabolite identifier.
#' @param step An integer indicating the number of steps to expand the neighborhood in the network.
#' @param graph The metabolic reaction network represented as an igraph object.
#'
#' @return A data frame containing the metabolite IDs and the corresponding steps in the network.
#'
#' @import igraph
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' # Assuming 'graph' is an igraph object representing the metabolic reaction network
#' # Get neighbors of metabolite with ID 'C03793' within 3 steps
#' getMrnNeighbor(id = 'C03793', step = 3, graph = MetDNA3::reaction_pair_network$version2$step2)
#'
#' @seealso
#' \code{\link{igraph::neighbors}}
#'
getMrnNeighbor <- function(id, step, graph) {

    # check valid ids
    valid_ids <- id[id %in% igraph::V(graph)$name]
    if (length(valid_ids) == 0) return(NULL)

    # get neighbors by step
    temp_step <- 1
    temp_ids <- valid_ids
    sum_ids <- valid_ids # step 0
    names(sum_ids) <- rep(0, length(sum_ids))
    while (temp_step <= step) {
        temp_neighbors <- lapply(temp_ids, function(x) names(igraph::neighbors(graph, x))) %>%
            unlist() %>% unique() %>% sort()
        names(temp_neighbors) <- rep(temp_step, length(temp_neighbors))
        sum_ids <- c(sum_ids, temp_neighbors)
        temp_step <- temp_step + 1
        temp_ids <- temp_neighbors
    }
    sum_ids <- sum_ids[!duplicated(sum_ids)]
    # if (length(which(sum_ids %in% valid_ids)) > 0) sum_ids <- sum_ids[-which(sum_ids %in% valid_ids)]
    if (length(sum_ids) == 0) return(NULL)

    # return data frame
    res <- data.frame(
        id = sum_ids,
        step = names(sum_ids)
    )
    return(res)
}
