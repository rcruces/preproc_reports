library(farver)

# Function - Color Palette Visualization
color.test <- function(color.vec) {
  n=length(color.vec)
  barplot(rep(1,n),space = 0,border = NA, axes = F, col = color.vec, main=deparse(substitute(color.vec)))
}

# Function that plot the circular interpolation of a given color vector
plot.palette <- function(color.vec,N,rad=2) {
  # col.vector - vector of colors as input
  # N - the number of out colors desired
  color.pal <- colorRampPalette(c(color.vec))
  pie(rep(1,N),col = color.pal(N), labels = NA,border = NA,radius = rad, main=deparse(substitute(color.vec)))
}

# -----------------------------------------------
# Vector with the colors of the BuGrRd palette
BuGrRd.vec <- c("navy", "royalblue3", "#5B83FF", "#73B5E3", "skyblue", "slategray1", 
               "gray60", "#F4E287", "#F1CB00", "#FF8F07", "#FB4200", "#CC0000", "#530000")

# Number of points to interpolate
N <- 500

# Colormap BuGrRd
BuGrRd <- colorRampPalette(BuGrRd.vec, space = "Lab")(N)

# CIELAB representation of the BuGrRd colormap
cielab <- decode_colour(BuGrRd, to = "lab")

# -----------------------------------------------
# Plots of the colormap
par(mfrow=c(2,1))
plot.palette(rev(BuGrRd.vec), N, rad=4)
color.test(BuGrRd.vec)

# cieL* lightness profile of the BuGrRd colormap
layout(matrix(c(1,2), nrow=2), heights=c(5,1))
# L*
plot(
  cielabF[,1],
  type="p",
  xlab="Palette index",
  ylab="Lightness L*",
  main="Lightness Profile",
  pch=20,
  cex=9, frame.plot=FALSE, col=BuGrRd, ylim = c(0,100)
)
abline(v=N/2, col="black", lty=2, lwd=2)
grid()
color.test(BuGrRd)

colorRampPalette(BuGrRd.vec, space = "Lab")(13)
