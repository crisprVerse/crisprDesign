#' Add representative coordinates from amino acid position annotations
#'
#' Convert amino acid position annotations into representative coordinates for
#' downstream one-dimensional smoothing or clustering analyses.
#'
#' The column annotation of the input \code{GuideSet} must contain a 
#' column named \code{aaPos}. This column
#' typically stores one or more amino acid positions per row, encoded as
#' semi-colon-separated integers (for example, \code{"120;121"}).
#'
#' When \code{positionMethod = "average"}, each row is assigned a single
#' coordinate equal to the floor of the mean of the positions listed in
#' \code{aaPos}.
#'
#' When \code{positionMethod = "max"}, the function uses the \code{aaChanges}
#' column, if available, to identify the amino acid position associated with
#' the largest editing score. The returned integer vector is named by the
#' corresponding amino acid substitution (for example, \code{"R314K"}). If
#' \code{aaChanges} is missing, the function falls back to \code{aaPos} and
#' returns an unnamed vector.
#'
#' When \code{positionMethod = "all"}, the function returns all amino acid
#' positions with editing weights greater than or equal to
#' \code{minEditingWeight}, based on \code{aaChanges}. If no position passes
#' the threshold, the first listed substitution is returned. Each element of
#' the returned list is a named integer vector, where names are amino acid
#' substitutions (for example, \code{"R314K"}). If \code{aaChanges} is
#' missing, the function falls back to \code{aaPos}.
#'
#' The \code{aaChanges} column is expected to contain strings such as
#' \code{"R314K(0.41);V315I(0.29);A311T(0.10)"}, where the amino acid position
#' is extracted from each substitution and the value in parentheses is treated
#' as its editing weight.
#'
#' @param guideSet A \code{GuideSet} object containing at least a column named
#'   \code{aaPos}. Optionally, it may also contain a column named
#'   \code{aaChanges}.
#' @param positionMethod Character string specifying how plotting coordinates
#'   should be assigned. Must be one of \code{"average"}, \code{"max"}, or
#'   \code{"all"}.
#' @param minEditingWeight Numeric scalar giving the minimum editing weight
#'   required for an edited amino acid position to be retained when
#'   \code{positionMethod = "all"}. Default is \code{0.3}.
#'
#' @return
#' If \code{positionMethod = "average"}, an integer vector of length
#' \code{nrow(df)} is added to the mcols of the GuideSet object. 
#'
#' If \code{positionMethod = "max"}, a named integer vector of length
#' \code{nrow(df)} when \code{aaChanges} is available; otherwise an unnamed
#' integer vector.
#'
#' If \code{positionMethod = "all"}, a list of length \code{nrow(df)}, where
#' each element is an integer vector containing one or more plotting
#' coordinates for that row. When \code{aaChanges} is available, each integer
#' vector is named by the corresponding amino acid substitution.
#'
#'
#' @importFrom stringr str_extract str_extract_all
#' @export
addRepCoordinates <- function(
    guideSet,
    positionMethod = c("average", "max", "all"),
    minEditingWeight = 0.3
){
    cols <- c("aaPos", "aaChanges")
    gsCols <- colnames(mcols(guideSet))
    cols <- intersect(cols, gsCols)
    guides <- mcols(guideSet)[, cols, drop=FALSE]
    guides <- as.data.frame(guides)
    coords <- .getRepCoordinates(guides, positionMethod=positionMethod)
    mcols(guideSet)$repCoordinate <- unname(coords)
    mcols(guideSet)$repSub <- names(coords)
    return(guideSet)
}








.getRepCoordinates <- function(
    df,
    positionMethod = c("average", "max", "all"),
    minEditingWeight = 0.3
) {

    positionMethod <- match.arg(positionMethod)

    if (!is.data.frame(df)) {
        stop("`df` must be a data.frame.")
    }
    if (!"aaPos" %in% colnames(df)) {
        stop("`df` must contain a column named `aaPos`.")
    }
    if (!is.numeric(minEditingWeight) || length(minEditingWeight) != 1L ||
        is.na(minEditingWeight)) {
        stop("`minEditingWeight` must be a single non-missing numeric value.")
    }

    aaPos <- as.character(df$aaPos)

    .parse_aa_pos_string <- function(x) {
        if (is.na(x) || !nzchar(x)) {
            return(NA_integer_)
        }
        vals <- suppressWarnings(as.integer(strsplit(x, split = ";", fixed = TRUE)[[1]]))
        vals <- as.integer(vals[!is.na(vals)])
        if (length(vals) == 0L) {
            return(NA_integer_)
        }
        vals
    }

    .average_position <- function(x) {
        pos <- .parse_aa_pos_string(x)
        if (length(pos) == 1L && is.na(pos)) {
            return(NA_integer_)
        }
        as.integer(floor(mean(pos)))
    }

    .extract_positions_scores_subs <- function(x) {
        if (is.na(x) || !nzchar(x)) {
            return(list(
                pos = integer(0),
                scores = numeric(0),
                subs = character(0)
            ))
        }

        subs <- stringr::str_extract_all(x, "[^;(]+(?=\\()")[[1]]
        scores <- as.numeric(
            stringr::str_extract_all(x, "(?<=\\()[^)]*(?=\\))")[[1]]
        )
        pos <- as.integer(stringr::str_extract(subs, "[0-9]+"))

        keep <- !is.na(pos)
        list(
            pos = pos[keep],
            scores = scores[keep],
            subs = subs[keep]
        )
    }

    .max_position_from_changes <- function(x, fallback) {
        parsed <- .extract_positions_scores_subs(x)
        pos <- parsed$pos
        scores <- parsed$scores
        subs <- parsed$subs

        if (length(pos) == 0L) {
            out <- .average_position(fallback)
            return(stats::setNames(out, NA_character_))
        }

        if (length(scores) == 0L || all(is.na(scores))) {
            return(stats::setNames(pos[1], subs[1]))
        }

        idx <- which.max(replace(scores, is.na(scores), -Inf))
        stats::setNames(pos[idx], subs[idx])
    }

    .all_positions_from_changes <- function(x, fallback, minEditingWeight = 0.3) {
        parsed <- .extract_positions_scores_subs(x)
        pos <- parsed$pos
        scores <- parsed$scores
        subs <- parsed$subs

        if (length(pos) == 0L) {
            fallback_pos <- .parse_aa_pos_string(fallback)
            fallback_pos <- fallback_pos[!is.na(fallback_pos)]
            if (length(fallback_pos) == 0L) {
                return(stats::setNames(NA_integer_, NA_character_))
            }
            return(stats::setNames(fallback_pos[1], NA_character_))
        }

        good <- scores >= minEditingWeight
        good[is.na(good)] <- FALSE

        if (any(good)) {
            return(stats::setNames(pos[good], subs[good]))
        }

        stats::setNames(pos[1], subs[1])
    }

    if (positionMethod == "average") {
        plottingPos <- vapply(aaPos, .average_position, FUN.VALUE = integer(1))
        return(plottingPos)
    }

    hasAaChanges <- "aaChanges" %in% colnames(df)

    # For indel screens / no aaChanges available
    if (!hasAaChanges) {
        fallback <- lapply(aaPos, .parse_aa_pos_string)

        if (positionMethod == "max") {
            return(vapply(fallback, function(x) x[1], FUN.VALUE = integer(1)))
        }

        if (positionMethod == "all") {
            return(lapply(fallback, function(x) x[!is.na(x)]))
        }
    }

    aaChanges <- as.character(df$aaChanges)

    if (positionMethod == "max") {
        plottingPos <- lapply(
            seq_len(nrow(df)),
            function(i) {
                if (is.na(aaChanges[i]) || !nzchar(aaChanges[i])) {
                    stats::setNames(.average_position(aaPos[i]), NA_character_)
                } else {
                    .max_position_from_changes(aaChanges[i], aaPos[i])
                }
            }
        )
        plottingPos <- unlist(plottingPos, use.names = TRUE)
        return(plottingPos)
    }

    if (positionMethod == "all") {
        plottingPos <- lapply(
            seq_len(nrow(df)),
            function(i) {
                if (is.na(aaChanges[i]) || !nzchar(aaChanges[i])) {
                    pos <- .parse_aa_pos_string(aaPos[i])
                    pos <- pos[!is.na(pos)]
                    return(pos)
                } else {
                    .all_positions_from_changes(
                        aaChanges[i],
                        aaPos[i],
                        minEditingWeight = minEditingWeight
                    )
                }
            }
        )
        return(plottingPos)
    }
}





