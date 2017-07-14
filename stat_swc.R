#Joe Snider
#7/17
#
#A stat for ggplot2 to plot a 2d swc file.
# Decodes the parent points and creates segments.

StatSWC <- 
  ggproto("StatSWC", Stat,
          
          compute_group = function(data, scales) {
            kp.parent <- data$parent > 0
            data.frame(x=data[kp.parent > 0, "x"],
                       y=data[kp.parent > 0, "y"],
                       xend=data[data[kp.parent,]$parent, "x"],
                       yend=data[data[kp.parent,]$parent, "y"])
          },
          
          required_aes = c("x", "y", "parent")
)

stat_swc <- function(mapping = NULL, data = NULL, geom = "segment",
                     position = "identity", na.rm = FALSE, show.legend = NA, 
                     inherit.aes = TRUE, ...) {
  layer(
    stat = StatSWC, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
