subset_diagnoses <- function(icd10_diagnoses, 
                             icd10_class, 
                             icd10_subclass=NULL, 
                             diagnosis_level=NULL,
                             keep_level=NULL) {
  
  # Check if AID column exists in the dataframe
  if (!"AID" %in% names(icd10_diagnoses)) {
    stop("Error: 'AID' column not found in the provided dataframe.\n\n")
  }
  
  # Warn and convert to uppercase if icd10_class is not in uppercase
  if(icd10_class != toupper(icd10_class)) {
    warning(paste("Converted", icd10_class, "to uppercase. Please provide ICD-10 codes in uppercase in the future.\n\n"))
    icd10_class <- toupper(icd10_class)
  }
  
  # If subclasses are specified, create a full list of ICD10 codes to match
  if (!is.null(icd10_subclass)) {
    subclasses <- sort(icd10_subclass)  # Sort subclasses in ascending order
    full_classes <- paste0(icd10_class, subclasses)  # Combine main class with each sub-class
  } else {
    full_classes <- icd10_class
  }
  
  # Use grep to find matching column names
  matching_cols <- unique(grep(paste(full_classes, collapse="|"), names(icd10_diagnoses), value = TRUE))
  
  # Check if there are columns that match the specified ICD10 class
  if (length(matching_cols) == 0) {
    stop(paste("Error: No columns found matching the specified ICD10 class(es):\n\n", paste(full_classes, collapse=", ")))
  }
  
  # Ensure AID column is included
  if (!("AID" %in% matching_cols)) {
    matching_cols <- c("AID", matching_cols)
  }
  
  # Subset the dataframe by columns
  subset_df <- icd10_diagnoses[, matching_cols, drop = FALSE]
  
  # If diagnosis_level is specified, subset rows based on this diagnosis level
  # diagnosis_level is the number 4,2,1 Xioran codes the 1st, 2nd, 3rd diagnoses
  if (!is.null(diagnosis_level)) {
    matching_rows <- apply(subset_df[-1], 1, function(row) any(diagnosis_level %in% row))
    subset_df <- subset_df[matching_rows, ]
  }
  
  # If keep_level is specified, filter duplicates by retaining only the specified level
  if (!is.null(keep_level)) {
    # Identify duplicated AIDs
    duplicate_AIDs <- subset_df$AID[duplicated(subset_df$AID, fromLast = TRUE) | duplicated(subset_df$AID)]
    
    # Filter out rows not having the keep_level for these duplicated AIDs
    filter_cond <- !(subset_df$AID %in% duplicate_AIDs & !apply(subset_df[-1], 1, function(row) any(keep_level %in% row)))
    subset_df <- subset_df[filter_cond, ]
  }
  
  # Adding secondary label based on icd10_class and icd10_subclass
  #if (!is.null(icd10_subclass)) {
  #  subset_df$secondary_label <- paste0(icd10_class, icd10_subclass)
  #} else {
  #  subset_df$secondary_label <- icd10_class
  #}
  
  return(subset_df)
}



###############################################################################
# An individual may have multiple assessments indicated by their ID followed by 
# different dates--therefore, an individual may appear on multiple rows. 
# Function lets the user choose the assessment number. 
# Using the assessment number, the variable cans_scores in the extract_cans 
# function will be subset to include individuals nth assessment
###############################################################################

filter_by_nth_assessment <- function(df, n) {
  # Extract the ID and date components of the AID
  df$ID_component <- sub(":.*", "", df$AID)
  df$Date_component <- as.Date(sub(".*:", "", df$AID))
  
  # Order by ID and Date
  df <- df[order(df$ID_component, df$Date_component), ]
  
  # Use ave() to generate a sequence number for each ID group
  df$seq_num <- ave(seq_along(df$AID), df$ID_component, FUN = seq_along)
  
  # Check which IDs have all the assessments in the range
  valid_ids <- unique(df[df$seq_num %in% max(n), "ID_component"])
  
  # Subset the data based on valid IDs and the nth assessments
  df <- df[df$ID_component %in% valid_ids & df$seq_num %in% n, ]
  
  # Drop temporary columns
  df$ID_component <- NULL
  df$Date_component <- NULL
  df$seq_num <- NULL
  
  return(df)
}



match_rows_by_AID <- function(icd10_diagnoses, cans_scores, n=NULL) {
  
  # Check if AID column exists in both dataframes
  if (!"AID" %in% names(icd10_diagnoses) || !"AID" %in% names(cans_scores)) {
    stop("Error: 'AID' column not found in one of the provided dataframes.")
  }
  
  # If n is provided, filter cans_scores by nth assessment
  if (!is.null(n)) {
    cans_scores <- filter_by_nth_assessment(cans_scores, n)
  }
  
  # Subset cans_scores based on AID values in icd10_diagnoses
  matched_rows <- cans_scores[cans_scores$AID %in% icd10_diagnoses$AID, ]
  matched_rows$group_label <- 1
  
  # Return the matched rows
  return(matched_rows)
}


# ORIGINAL
generate_control_group <- function(icd10_diagnoses, 
                                   cans_scores, 
                                   icd10_controls, 
                                   control_sample_size) {
  
  # Identify the rows to exclude based on the provided ICD10 classes
  rows_to_exclude <- apply(icd10_diagnoses, 1, function(row) {
    any(sapply(icd10_controls, function(code) grepl(paste0("^", code), row)))
  })
  
  # Subset icd10_diagnoses to remove rows with the specified ICD10 classes
  control_diagnoses_subset <- icd10_diagnoses[!rows_to_exclude, ]
  
  # If the desired sample size is greater than the number of available rows, adjust the sample size
  if (control_sample_size > nrow(control_diagnoses_subset)) {
    control_sample_size <- nrow(control_diagnoses_subset)
  }
  
  # Sample rows to get the control group
  control_sample_indices <- sample(1:nrow(control_diagnoses_subset), control_sample_size)
  control_diagnoses_sample <- control_diagnoses_subset[control_sample_indices, ]
  
  # Extract matching rows from cans_scores using the AID column from control_diagnoses_sample
  control_cans_scores_subset <- match_rows_by_AID(control_diagnoses_sample, 
                                                  cans_scores)
  
  return(control_cans_scores_subset)
}


filter_cans_by_cases <- function(cans_scores, cases) {
  
  # Check if ID column exists in cans_scores; if not, create it
  if (!"ID" %in% colnames(cans_scores)) {
    cans_scores$ID <- sapply(strsplit(as.character(cans_scores$AID), ":"), `[`, 1)
  }
  
  # Check if ID column exists in cases; if not, create it
  if (!"ID" %in% colnames(cases)) {
    cases$ID <- sapply(strsplit(as.character(cases$AID), ":"), `[`, 1)
  }
  
  # Get common columns between cans_scores and cases
  common_columns <- intersect(colnames(cans_scores), colnames(cases))
  
  # Subset cans_scores to retain only the common columns
  cans_scores <- cans_scores[, common_columns, drop=FALSE]
  
  # Filter out rows from cans_scores where ID is found in cases dataframe
  cans_minus_cases <- cans_scores[!(cans_scores$ID %in% cases$ID), ]
  
  # Take a sample from the filtered cans_scores
  control_sample_size <- nrow(cases)
  if (control_sample_size > nrow(cans_minus_cases)) {
    control_sample_size <- nrow(cans_minus_cases)
  }
  control_sample_indices <- sample(1:nrow(cans_minus_cases), control_sample_size)
  control_diagnoses_sample <- cans_minus_cases[control_sample_indices, ]
  
  return(control_diagnoses_sample)
}



# Main function
# icd10_diagnoses: df of diagnoses for all assessments 
# cans_scores: df of cans assessment scores
# icd10_class=NULL: choose ICD10 class ("F32")
# icd10_subclass=NULL: choose ICD10 subclass (.1) 
# return_diagnoses=FALSE: Return a dataframe with filtered diagnoses
# n=NULL: choose to filter by 1st, 2nd, nth assessment number, or 1 to n assessments
# diagnosis_level=NULL: Filter based on whether icd10 class is 1,2,4
# keep_level=NULL: choose which levels to keep when there are duplicates
# generate_control_group=FALSE: for case/control analyses
# icd10_controls=c()): choose specific icd10 classes for control, or leave blank
#                      to randomly choose.
extract_cans_by_icd10 <- function(icd10_diagnoses, 
                                  cans_scores, 
                                  icd10_class=NULL, 
                                  icd10_subclass=NULL, 
                                  return_diagnoses=FALSE, 
                                  n=NULL, 
                                  diagnosis_level=NULL,
                                  keep_level=NULL,
                                  generate_control_group=FALSE,
                                  icd10_controls=c()) {
  
  # If diagnosis_class is specified, subset icd10_diagnoses
  if (!is.null(icd10_class)) {
    icd10_diagnoses_subset <- subset_diagnoses(icd10_diagnoses, 
                                               icd10_class, 
                                               icd10_subclass, 
                                               diagnosis_level,
                                               keep_level)
  } else {
    icd10_diagnoses_subset <- icd10_diagnoses
  }
  
  # Extract matching rows from cans_scores using the AID column from icd10_diagnoses_subset
  cans_scores_subset <- match_rows_by_AID(icd10_diagnoses_subset, 
                                          cans_scores, 
                                          n)
  
  # Match the secondary label to the cans_scores_subset
  secondary_labels <- icd10_diagnoses_subset[icd10_diagnoses_subset$AID %in% cans_scores_subset$AID, "secondary_label"]
  cans_scores_subset$secondary_label <- secondary_labels[match(cans_scores_subset$AID, icd10_diagnoses_subset$AID)]
  
  # If generate_control_group is TRUE, generate the control group and bind it to cans_scores_subset
  if (generate_control_group) {
    control_data <- generate_control_group(icd10_diagnoses, 
                                           cans_scores, 
                                           icd10_controls, 
                                           nrow(cans_scores_subset))
    control_data$group_label <- 0
    
    cans_scores_subset <- rbind(cans_scores_subset, 
                                control_data)
  }
  
  # Depending on the user's request, return one or both dataframes
  if (return_diagnoses) {
    return(list(icd10_diagnoses=icd10_diagnoses_subset, cans_scores=cans_scores_subset))
  } else {
    return(as.data.frame(cans_scores_subset))
  }
}


generate_diagnosis_df <- function(cans_extracted, diagnoses) {
  
  # Matching rows based on AID
  matching_rows <- diagnoses[diagnoses$AID %in% cans_extracted$AID, ]
  
  # Function to generate the column based on value
  generate_column <- function(df, value) {
    value_matrix <- df == value
    
    # Replace TRUE values in the matrix with the column names
    value_names_matrix <- matrix(ncol=ncol(value_matrix), nrow=nrow(value_matrix))
    
    for (i in 1:ncol(value_matrix)) {
      value_names_matrix[,i] <- ifelse(value_matrix[,i], colnames(df)[i], NA)
    }
    
    # Extract the first non-NA value for each row
    value_col <- apply(value_names_matrix, 1, function(x) {
      val <- x[!is.na(x)]
      if(length(val) > 0) {
        return(val[1])
      } else {
        return(NA)
      }
    })
    
    return(value_col)
  }
  
  # Generate the three columns
  primary_col <- generate_column(matching_rows, 4)
  secondary_col <- generate_column(matching_rows, 2)
  tertiary_col <- generate_column(matching_rows, 1)
  
  # Create and return the new dataframe
  result_df <- data.frame(
    AID = cans_extracted$AID,
    Primary = primary_col,
    Secondary = secondary_col,
    Tertiary = tertiary_col
  )
  
  return(result_df)
}


# Use for setting up dataframe for analysis at the individual level
# Input raw cans extraction which will be imputed and filtered
# Each individual has 'threshold' number of assessments
# AID is split into ID and Date

filter_and_impute <- function(cans_extracted, threshold, group_label_col) {
  # Extract the IDs as before
  split_ids <- strsplit(cans_extracted$AID, ":")
  ids_only <- sapply(split_ids, function(x) x[1])
  
  # Create a table of counts
  id_counts <- table(ids_only)
  
  # Filter the IDs based on the threshold
  filtered_ids <- names(id_counts[id_counts >= threshold])
  
  # Filter the original dataframe using these IDs
  filtered_data <- cans_extracted[ids_only %in% filtered_ids,]
  
  # Filter cans_extracted using AIDs from unique_filtered_data
  final_data <- cans_extracted[cans_extracted$AID %in% filtered_data$AID,]
  
  # Impute the dataframe without the AID and group_label columns
  imputed_data <- impute_with_mode(final_data[, !names(final_data) %in% c("AID", group_label_col)], threshold = 0.2)
  
  # Bind the AID and group_label columns from the original dataframe to the imputed data
  result_df <- cbind(final_data[, c("AID", group_label_col)], imputed_data)
  
  # Split the AID column into two new columns: ID and Date
  result_df$ID <- sapply(strsplit(result_df$AID, ":"), `[`, 1)
  result_df$Date <- sapply(strsplit(result_df$AID, ":"), `[`, 2)
  
  return(result_df)
}

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

impute_with_mode <- function(df, threshold = 0.2) {
  
  # Calculate the proportion of NA values for each column
  na_proportions <- sapply(df, function(col) sum(is.na(col)) / length(col))
  
  # Identify columns to drop based on the threshold
  columns_to_drop <- names(na_proportions[na_proportions > threshold])
  
  # Drop columns that exceed the threshold of NA values
  df <- df[, !names(df) %in% columns_to_drop, drop = FALSE]
  
  # Impute NAs with mode
  df[] <- lapply(df, function(col) {
    if(any(is.na(col))) {
      mode_value <- get_mode(col[!is.na(col)])
      col[is.na(col)] <- mode_value
    }
    return(col)
  })
  
  # Test for ordinal values between 0 and 3
  if (any(df < 0 | df > 3, na.rm = TRUE)) {
    stop("Data contains values outside the ordinal range [0,3].")
  }
  
  # Test if all values are integers
  if (any(df != round(df), na.rm = TRUE)) {
    stop("Data contains non-integer values.")
  }
  
  return(df)
}