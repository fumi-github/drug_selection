# Install and load the stringr package if you haven't already
# install.packages("stringr")
library(stringr)

# 2025.06.27 added sach sachet SACHET for bile acid sequestrants

#' Calculate Total Tablets/Capsules/Days from Prescription String
#'
#' Parses various string formats to extract the total number of units (tablets/capsules/days).
#' Handles simple numbers, numbers with units, multipliers (e.g., N*M, NxM or N packs of M),
#' removes surrounding parentheses or square brackets, and removes a leading 'x'.
#'
#' @param input_string A character string describing the quantity.
#'
#' @return The total number of units as a numeric value, or NA if parsing fails.
#'
#' @examples
#' calculate_total_units("28 tab")
#' calculate_total_units("2*28 tablet(s)")
#' calculate_total_units("28x2")
#' calculate_total_units("3 packs of 28 capsule(s)")
#' calculate_total_units("56.000")
#' calculate_total_units("(28)")
#' calculate_total_units("[56]")
#' calculate_total_units("84 DAYS")
#' calculate_total_units("x28")
#' calculate_total_units("X56 capsules") # Test leading X with unit

calculate_total_units <- function(input_string) {
  # Return NA for empty or invalid input
  if (is.null(input_string) || is.na(input_string) || nchar(trimws(input_string)) == 0) {
    return(NA_real_)
  }
  
  # --- Preprocessing ---
  # 1. Remove surrounding brackets/parentheses
  processed_string <- str_replace(input_string, "^\\((.*)\\)$", "\\1")
  processed_string <- str_replace(processed_string, "^\\[(.*)\\]$", "\\1")
  
  # 2. Remove a single 'x' or 'X' if it's the very first character
  #    Use ^[xX] to match 'x' or 'X' only at the beginning (^)
  processed_string <- str_replace(processed_string, "^[xX]", "")
  
  # 3. Convert to lowercase and trim whitespace
  text <- str_to_lower(trimws(processed_string))
  
  # --- Pattern 1: Multiplication ---
  # Case A: "N * M unit" or "N x M unit" (e.g., "2*28 tablet(s)", "28x2")
  match_mult_op <- str_match(text, "^(\\d+)\\s*[\\*x]\\s*(\\d+)")
  if (!is.na(match_mult_op[1, 1])) {
    num1 <- as.numeric(match_mult_op[1, 2])
    num2 <- as.numeric(match_mult_op[1, 3])
    if (!is.na(num1) && !is.na(num2)) {
      return(num1 * num2)
    }
  }
  
  # Case B: "N packs of M unit" (e.g., "2 packs of 28 tablet(s)")
  match_mult_pack <- str_match(text, "^(\\d+)\\s+packs?\\s+of\\s+(\\d+)")
  if (!is.na(match_mult_pack[1, 1])) {
    num1 <- as.numeric(match_mult_pack[1, 2])
    num2 <- as.numeric(match_mult_pack[1, 3])
    if (!is.na(num1) && !is.na(num2)) {
      return(num1 * num2)
    }
  }
  
  # --- Pattern 2: Simple Number before Unit ---
  # Case C: "N unit" (e.g., "28 tab", "30 tablets - 2.5 mg", "84 days")
  match_simple <- str_match(text, "^(\\d+(\\.\\d+)?)\\s*[-–]?\\s*(tab|cap|tablet|capsule|day|days|sach|sachet)")
  if (!is.na(match_simple[1, 1])) {
    num <- as.numeric(match_simple[1, 2])
    if (!is.na(num)) {
      return(num)
    }
  }
  
  # --- Pattern 3: Simple Number Only ---
  # Case D: "N" or "N.N" (e.g., "56.000", handles output of "x28" -> "28")
  match_numeric_only <- str_match(text, "^(\\d+(\\.\\d+)?)$")
  if (!is.na(match_numeric_only[1, 1])) {
    num <- as.numeric(match_numeric_only[1, 2])
    if (!is.na(num)) {
      return(num)
    }
  }
  
  # --- No Match Found ---
  return(NA_real_)
}

# --- Example Usage ---

# Define the input strings based on your examples, including the new ones
input_examples <- c(
  "56.000",
  "28 tab",
  "30 tablet",
  "28 Tabs",
  "28 TABLET(S)",
  "28 - tablet(s)",
  "30 tablets - 2.5 mg",
  "28 tablet(s) - 2.5 mg",
  "2*28 tablet(s) - 2.5 mg",
  "28x2",
  " 5 x 28 ",
  "1 pack of 28 tablet(s)",
  "2 packs of 28 tablet(s)",
  "28 capsules",
  "28 – capsule",
  "3 packs of 28 capsule(s)",
  "84 DAYS",
  "28 day supply",
  "(28)",
  "[56]",
  " ( 28 tabs ) ",
  "x28",               # New example with leading 'x'
  "X56 capsules",      # New example with leading 'X' and unit
  " x 14 days",      # New example with leading 'x' and space (should NOT remove x)
  "Text starting with x", # Should not remove 'x' here
  "Some unrelated text",
  NA,
  ""
)

# Apply the function to each example
results <- sapply(input_examples, calculate_total_units, USE.NAMES = FALSE)

# Display the results alongside the inputs
output_df <- data.frame(Input = input_examples, Calculated_Units = results)
print(output_df)
