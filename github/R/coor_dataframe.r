#' Convert junction bin to a dataframe including the first and second junction positions
#'
#' This function takes a character vector of junction coordinates and converts it into a dataframe. For each junction, it calculates the positions of the first and second junction blocks based on a specified window size.
#'
#' @param coordinates A character vector of junction coordinates in the format "chrN:start-chrM:end".
#' @param window The window size used for creating genomic regions around junction blocks (default is 10000).
#' @return A list of dataframes, where each dataframe contains information about the first and second junction positions.
#' @examples
#' coords <- c("chr1:1000-chr2:2000", "chr3:3000-chr4:4000")
#' coor_dataframe(coords, window = 5000)
#'
#' @importFrom base regmatches
#' @importFrom base gregexpr
#' @importFrom base as.numeric
#' @export
coor_dataframe<-function(coordinates,window=10000){
matches <- regmatches(coordinates, gregexpr("\\d+", coordinates))[[1]]

# Extract chromosome number and positions for first junction block
first_chromosome <- matches[1]
first_start_position <- as.numeric(matches[2])*1000
first_end_position <- as.numeric(matches[2])*1000 + window

# Extract chromosome number and positions for second junction block
second_chromosome <- matches[3]
second_start_position <- as.numeric(matches[4])*1000
second_end_position<- as.numeric(matches[4])*1000 + window

# Construct dataframe
first_junction <- data.frame(chr = paste0("chr",first_chromosome),
                    start = first_start_position,
                    end = first_end_position
                    )
second_junction <- data.frame(chr = paste0("chr",second_chromosome),
                    start = second_start_position,
                    end = second_end_position
                    )  
res<-list(first_junction, second_junction)
                     
return(res)                                    
}                    
