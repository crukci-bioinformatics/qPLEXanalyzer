# Mean variance plot
plotMeanVar <- function(MSnSetObj, title="") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    
    intensities <- log(exprs(MSnSetObj) + 0.001)
    mvDat <- data.frame(
        Mean = rowMeans(intensities),
        Variance = apply(intensities, 1, var)
    )
    
    ssDat <- smooth.spline(x = mvDat$Mean, y = mvDat$Variance, spar = 1) %$%
        data.frame(x = x, y = y)
    
    ggplot(mvDat, aes(x = Mean, y = Variance)) +
        geom_point(size = 0.5, alpha = 0.6, colour = "darkblue") +
        geom_line(data = ssDat, aes(x = x, y = y), colour = "red", size = 0.5) +
        theme_bw() +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5))
}

