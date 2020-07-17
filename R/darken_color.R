  .darken_color <- function(color, factor=1.4){
    col <- grDevices::col2rgb(color)
    col <- col/factor
    col <- grDevices::rgb(t(col), maxColorValue=255)
  }
