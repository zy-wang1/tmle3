#' Defines a tmle (minus the data)
#'
#' Current limitations:
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_vimPAR <- R6Class(
  classname = "tmle3_Spec_risk",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    make_params = function(tmle_task, likelihood) {
      # todo: export and use sl3:::get_levels
      A_vals <- tmle_task$get_tmle_node("A")
      if (is.factor(A_vals)) {
        A_levels <- sort(unique(A_vals))
        A_levels <- factor(A_levels, levels(A_vals))
      } else {
        A_levels <- sort(unique(A_vals))
      }
      tsm_params <- lapply(A_levels, function(A_level) {
        intervention <- define_lf(LF_static, "A", value = A_level)
        tsm <- Param_TSM$new(likelihood, intervention)
        return(tsm)
      })
        mean_param <- Param_mean$new(likelihood)

        ate_params <- lapply(tsm_params, function(comparison_param){
          Param_delta$new(likelihood, delta_param_ATE, list(mean_param,
            comparison_param))
        })

        tmle_params <- c(tsm_params, mean_param, ate_params)

      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' Estimates PAR using TSM and observed mean
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Variables to Consider for VIM
#' Y=Outcome
#' @importFrom sl3 make_learner Lrnr_mean
#' @export
tmle_vimPAR <- function() {
  tmle3_Spec_vimPAR$new()
}
