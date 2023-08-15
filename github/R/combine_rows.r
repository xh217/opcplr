#' Function to combine rows based on shared elements
#'
#' This function takes a data frame and combines rows that have shared elements. It iterates through each row and checks for shared elements with subsequent rows. If shared elements are found, the rows are combined.
#'
#' @param data A data frame where rows need to be combined.
#' @param sep The separator used to split and combine elements (default is ";").
#' @return A list of combined rows.
#' @examples
#' data <- data.frame(V1 = c("A;B", "C;D", "B;E", "F"), stringsAsFactors = FALSE)
#' combined <- combine_rows(data, sep = ";")
#'
#' @importFrom base nrow
#' @importFrom base paste
#' @importFrom base unique
#' @importFrom base any
#' @importFrom base strsplit
#' @export

combine_rows <- function(data,sep = ";") {
  # Initialize a list to store the combined rows
  tab <- sep
  combined_rows <- list()
  
  # Initialize a vector to keep track of processed rows
  processed <- logical(nrow(data))
  
  # Iterate through each row
  for (i in 1:(nrow(data)-1)) {
  	
    if (!processed[i]) {
      # Get the elements from the current row
      elements <- unlist(strsplit(data[i,], tab))
      
      # Initialize a vector to store the combined elements
      combined_elements <- elements
      
      # Check if other rows have shared elements
      for (j in (i + 1):nrow(data)) {
        if (!processed[j]) {
          # Get the elements from the current row being compared
          compare_elements <- unlist(strsplit(data[j, ], tab))
          
          # Check if there are shared elements
          if (any(compare_elements %in% combined_elements)) {
            # Combine the elements
            combined_elements <- unique(c(combined_elements, compare_elements))
            
            # Mark the row as processed
            processed[j] <- TRUE
          }
        }
      }
      
      # Store the combined row
      combined_rows[[length(combined_rows) + 1]] <- paste(combined_elements, collapse = tab)
    }
  }
  
  # Return the combined rows
  return(combined_rows)
}
