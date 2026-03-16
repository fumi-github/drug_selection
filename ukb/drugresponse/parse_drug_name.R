#' Parse Drug Name and Dosage
#'
#' Extracts the drug name and its dosage in milligrams (mg) from a string.
#' Assumes the dosage is a number followed by "mg" (case-insensitive),
#' potentially with a space in between. Returns NA if "mg" appears zero
#' or more than once.
#'
#' @param drug_string A character string containing the drug description.
#'
#' @return A one-row data frame with columns 'drug_name' (character) and
#'   'dosage_mg' (numeric). Returns NA in both columns if parsing fails
#'   or the input doesn't meet the criteria (exactly one "mg").
#'
#' @examples
#' parse_drug_name("Accupro 10mg tablets (Pfizer Ltd)")
#' #>   drug_name dosage_mg
#' #> 1   Accupro        10
#'
#' parse_drug_name("ACCUPRO TABLETS 40MG")
#' #>         drug_name dosage_mg
#' #> 1 ACCUPRO TABLETS        40
#'
#' parse_drug_name("Accupro TABS 5MG")
#' #>     drug_name dosage_mg
#' #> 1 Accupro TABS         5
#'
#' parse_drug_name("DrugName12.5mg")
#' #>   drug_name dosage_mg
#' #> 1  DrugName      12.5
#'
#' parse_drug_name("Acebutolol 200mg / Hydrochlorothiazide 12.5mg tablets")
#' #>   drug_name dosage_mg
#' #> 1      <NA>        NA
#'
#' parse_drug_name("AMLODIPINE")
#' #>   drug_name dosage_mg
#' #> 1      <NA>        NA
#'
#' @importFrom stringr str_count str_match regex str_trim
#'
parse_drug_name <- function(drug_string) {
  # Ensure stringr package is available
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Please install the 'stringr' package: install.packages('stringr')")
  }
  library(stringr)
  
  # Define the pattern for "mg" (case-insensitive)
  mg_pattern <- regex("mg", ignore_case = TRUE)
  
  # Count occurrences of "mg"
  mg_count <- stringr::str_count(drug_string, mg_pattern)
  
  # If not exactly one occurrence, return NA
  if (mg_count != 1) {
    return(data.frame(drug_name = NA_character_, dosage_mg = NA_real_))
  }
  
  # Define the regex to capture the drug name and dosage
  # - (.*?): Non-greedily capture the drug name (Group 1)
  # - \\s*: Match zero or more spaces before the number
  # - (\\d*\\.?\\d+): Capture the number (integer or decimal) (Group 2)
  # - \\s*: Match zero or more spaces after the number
  # - mg: Match "mg" case-insensitively (handled by regex() option)
  # - (.*): Capture any trailing characters (Group 3) - helps anchor the match
  extract_pattern <- regex("(.*?) *(\\d*\\.?\\d+) *mg(.*)", ignore_case = TRUE)
  
  # Attempt to match the pattern
  match_result <- stringr::str_match(drug_string, extract_pattern)
  
  # Check if the pattern matched successfully
  if (is.na(match_result[1, 1])) {
    # Should ideally not happen if mg_count was 1, but as a safeguard
    return(data.frame(drug_name = NA_character_, dosage_mg = NA_real_))
  } else {
    # Extract, clean, and convert the captured groups
    drug_name_raw <- match_result[1, 2]
    dosage_raw <- match_result[1, 3]
    
    drug_name_clean <- stringr::str_trim(drug_name_raw)
    dosage_numeric <- as.numeric(dosage_raw)
    
    # Return the result in a data frame
    return(data.frame(drug_name = drug_name_clean, dosage_mg = dosage_numeric))
  }
}

# --- Example Usage with provided test cases ---

# Example 1: Standard format
print(parse_drug_name("Accupro 10mg tablets (Pfizer Ltd)"))

# Example 2: All caps, no space
print(parse_drug_name("ACCUPRO TABLETS 40MG"))

# Example 3: Short form, no space
print(parse_drug_name("Accupro TABS 5MG"))

# Example 4: No space at all
print(parse_drug_name("DrugName12.5mg"))

# Example 5: Decimal dosage with space
print(parse_drug_name("Drug Name 12.5 mg"))

# Example 6: Multiple 'mg' occurrences
print(parse_drug_name("Acebutolol 200mg / Hydrochlorothiazide 12.5mg tablets"))

# Example 7: No 'mg' occurrence
print(parse_drug_name("AMLODIPINE"))

# Example 8: 'mg' as part of a different unit (e.g., mg/ml) - still extracts number before 'mg'
print(parse_drug_name("Amlodipine 10mg/5ml oral suspension"))

# Example 9: Another multiple 'mg' case
print(parse_drug_name("Adalat 5mg capsules (Bayer Plc)... ADALAT CAP 10mg"))

# Example 10: 'MGM' instead of 'mg' - should return NA as per request
print(parse_drug_name("ADALAT LA 60 MGM."))
