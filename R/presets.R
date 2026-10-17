#' Disease-specific LOT rule presets
#'
#' Starting points for common BPC cohorts. Each preset is an ordinary
#' [lot_rules()] object, so any setting can be changed afterwards. Drug names
#' match by prefix and ignore case ("Goserelin" matches "Goserelin Acetate").
#'
#' * `"NSCLC"`: carboplatin/cisplatin swaps stay in the same line. Dropping
#'   chemo and continuing pemetrexed or immunotherapy (maintenance) stays in
#'   the same line.
#' * `"CRC"`: fluorouracil/capecitabine swaps stay in the same line; adding
#'   bevacizumab or leucovorin later does not start a new line.
#' * `"BrCa"`: switching between aromatase inhibitors stays in the same line;
#'   adding ovarian suppression or HER2 antibodies later does not start a new
#'   line; longer 180-day gap, since hormone therapy runs for years.
#' * `"Custom"`: package defaults.
#'
#' These are reasonable defaults, not clinical guidance. Review them with a
#' clinician for your analysis.
#'
#' @param disease One of "NSCLC", "CRC", "BrCa", "Custom".
#' @return A [lot_rules()] object.
#' @examples
#' lot_preset("CRC")
#' @export
lot_preset <- function(disease = c("NSCLC", "CRC", "BrCa", "Custom")) {
  disease <- match.arg(disease)
  switch(disease,
    NSCLC = lot_rules(window = 28, gap = 120,
                      swaps = list(c("Carboplatin", "Cisplatin"))),
    CRC = lot_rules(window = 28, gap = 120,
                    swaps = list(c("Fluorouracil", "Capecitabine")),
                    add_without_new_line = c("Bevacizumab", "Leucovorin")),
    BrCa = lot_rules(window = 28, gap = 180,
                     swaps = list(c("Anastrozole", "Letrozole", "Exemestane")),
                     add_without_new_line = c("Goserelin", "Leuprolide",
                                              "Trastuzumab", "Pertuzumab")),
    Custom = lot_rules()
  )
}
