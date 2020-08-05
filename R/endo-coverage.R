#' Compute coverage over a Ranges object
#'
#' @param x a `Ranges` object
#' @param shift shift how much should each range in x be shifted by? (default = 0L)
#' @param width width how long should the returned coverage score be?
#' This must be either a  positive integer or NULL (default = NULL)
#' @param weight weight how much weight should be assigned to each range? Either
#' an integer or numeric vector or a column in x. (default = 1L)
#' @param ... other optional parameters to pass to coverage
#'
#' @return An expanded Ranges object with a score column corresponding to
#' the coverage value over that interval. Note that compute_coverage
#' drops metadata associated with the orginal ranges.
#' @seealso \code{IRanges::\link[IRanges:coverage-methods]{coverage()}}, 
#'          \code{GenomicRanges::\link[GenomicRanges:coverage-methods]{coverage()}}
#' @examples
#' rng <- as_iranges(data.frame(start = 1:10, width = 5))
#' compute_coverage(rng)
#' compute_coverage(rng, shift = 14L)
#' compute_coverage(rng, width = 10L)
#' @importFrom IRanges coverage ranges
#' @importFrom S4Vectors runValue
#' @export
compute_coverage <- function(x, shift, width, weight, ...) {
  UseMethod("compute_coverage")
}

#' @export
compute_coverage.default <- function(x, shift = 0L, width = NULL, weight = 1L, ...) {
  as_ranges(coverage(x, shift, width, weight, ...))
}



#' @export
compute_coverage.GroupedGenomicRanges <- function(.data, shift = 0L, width = NULL, weight = 1L, ...) {
    #split out the dataframe by group
    splitted_list <- S4Vectors::split(.data@delegate, .data@group_indices@values)
    seqinfo = .data@delegate@seqinfo #store the sequence information for us to reset

    #compute the coverage for each individual group
    coverage_list <- lapply(splitted_list, function(x) as_ranges(coverage(x, shift, width, weight, ...)))

    #convenience function for setting the sequence information 
    set_seqinfo <- function(x, info) {
    seqinfo(x) <- info
    x
    }
    #push the lists back together making sure they keep the same seqinfo
    merged_list <- as_granges(GRangesList(lapply(coverage_list, function(x) set_sequences(x, seqinfo))))

    #remove the group_name column and return the grouped Genomic Ranges
    merged_list$group_name <- NULL
    group_by(merged_list, group)
}


setMethod("coverage", "DelegatingGenomicRanges",
          function(x, shift = 0L, width = NULL, weight = 1L, ...) {
            coverage(load_delegate(x), shift, width, weight, ...)
          })

