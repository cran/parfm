################################################################################
#  Print method for class 'select.parfm'                                       #
################################################################################
#                                                                              #
#  This function prints the objects of class 'parfm'                           #
#                                                                              #
#  Its parameters are                                                          #
#   - x         : the object of class 'select.parfm'                           #
#   - digits    : number of significant digits                                 #
#   - na.prints : character string indicating NA values in printed output      #
#                                                                              #
#                                                                              #
#   Date:                 January, 10, 2012                                    #
#   Last modification on: December  2, 2016                                    #
################################################################################

print.select.parfm <- function(x,
                               digits=3,
                               na.print="----",
                               ...) {
    if (missing(digits))
        digits <- 3
    
    if (!is.null(x)) {
        x <- lapply(x, function(c) {
            C <- cbind(paste(' ', rownames(c)), 
                         format(c, digits = digits, justify = "right", 
                                width = max(nchar(colnames(c)))))
            rownames(C) <- rep("", nrow(C))
            C[grepl('NA', C)] <- na.print
            C[, -1] <- format(C[, -1], nsmall = max(nchar(C[, -1])),
                                justify = "right")
            colnames(C) <- format(colnames(C), justify = 'right')
            return(C)
        })

        # AIC
        colnames(x[[1]])[1] <- 'AIC:'
        print(x[[1]], quote = FALSE)
        
        cat('\n')
        
        # BIC
        colnames(x[[2]])[1] <- 'BIC:'
        print(x[[2]], quote = FALSE)
    } 
}
