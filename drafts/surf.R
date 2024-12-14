# -------------------------------------------------------------------------------------
#### Cortical surface plotting ####
require('fsbrain')           # version 0.4.2
require('rgl')               # version 0.1.54
library('freesurferformats')  # read freesurfer data as arrays #requires giftit package
library("RColorBrewer")

# Helper function
plot_surface <-function(brainMesh, legend='', view_angles=c('sd_lateral_lh', 'sd_medial_lh', 'sd_medial_rh', 'sd_lateral_rh'), img_only=FALSE, Out) {
  try(img <- vis.export.from.coloredmeshes(brainMesh, colorbar_legend = legend, grid_like = FALSE, view_angles = view_angles, img_only = img_only, horizontal=TRUE, output_img = Out, silent = TRUE), silent = TRUE)
  while (rgl.cur() > 0) { close3d() }; file.remove(list.files(path = getwd(), pattern = 'fsbrain'))
  return(img)
}

# Colormap
BlGrRd_g90 <- colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "gray90", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")))
BlGyRd_centered <- colorRampPalette(c(rev(brewer.pal(9,"Blues")[3:9]), rep("gray90",3), brewer.pal(9,"Reds")[3:9]))
BlGyRd <- colorRampPalette(c(rev(brewer.pal(9,"Blues")[2:7]), rep("gray90",3), brewer.pal(9,"Reds")[3:9]))


# Set the path to the surface
fsLR32k.lh <- read.fs.surface(filepath = "/Users/rcruces/git_here/micapipe/surfaces/fsLR-32k.L.inflated.surf.gii")
fsLR32k.rh <- read.fs.surface(filepath = "/Users/rcruces/git_here/micapipe/surfaces/fsLR-32k.R.inflated.surf.gii")

# Load labels
aparc <- read.csv("/Users/rcruces/git_here/micapipe/parcellations/aparc_conte69.csv", header = FALSE)

map_to_labels <- function(source_val, target_lab, mask = NULL, fill = 0) {
 
  # Find all the unique values
  uniq <- unique(target_lab)[[1]]
  
  # Create an empty array to store the mapped values
  mapped <- rep(fill, length(target_lab))
  
  for (i in uniq) { 
    # Find the indices of the target labels that match the current unique value
    idx <- which(target_lab == i)
    
    # If a mask is provided, only map the values that are in the mask
    if (!is.null(mask)) {
      idx <- idx[mask[idx]]
    }
    
    # Map the values
    mapped[idx] <- source_val[i+1]
  }
  
  return(mapped)
}

mk6240.diff.H.surf <- map_to_labels(source_val = as.vector(mk6240.diff.H), target_lab = aparc)
mk6240.diff.P.surf <- map_to_labels(source_val = as.vector(mk6240.diff.P), target_lab = aparc)

mk6240.P <- map_to_labels(source_val = as.vector(colMeans(mk.all[mk.all$type=="P",82:152])), target_lab = aparc)

# Create the coloredmeshes
lf= limit_fun(-0.05, 0.05)
cml = coloredmesh.from.preloaded.data(fsLR32k.lh, morph_data = lf(mk6240.diff.P.surf[1:32492]), hemi = 'lh', makecmap_options = list('range'=c(-0.05, 0.05),'colFn'=BlGrRd_g90))
cmr = coloredmesh.from.preloaded.data(fsLR32k.rh, morph_data = lf(mk6240.diff.P.surf[32493:64984]), hemi = 'rh', makecmap_options = list('range'=c(-0.05, 0.05), 'colFn'=BlGrRd_g90))
sph.nat <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), rglactions = list('no_vis'=T))
# Plot the surface
plot_surface(sph.nat, 'mk6240.diff-PX', Out = "~/Desktop/mk6240.diff.P.surf-brainviews.png")

# Call high-level API for live plot.
library(viridis)
inferno_r <- colorRampPalette(inferno(256))
surfaces = hemilist(fsLR32k.lh, fsLR32k.rh);
pvd  = hemilist(mk6240.P[1:32492], mk6240.P[32493:64984]);
cm = vis.subject.pre(surfaces, pvd , draw_colorbar = T, rglactions = list('trans_fun'=limit_fun(0.5, 1.4)), makecmap_options = list('range'=c(0.5, 1.4), 'colFn'=inferno_r));
cm = vis.subject.pre(surfaces, pvd , draw_colorbar = T, rglactions = list('trans_fun'=limit_fun(0.5, 1.4)), makecmap_options = list('range'=c(0.5, 1.4), 'colFn'=inferno_r));

# Save plot rho values of correlations
export(cm, output_img = "~/Desktop/mk6240.P.png", grid_like = FALSE, colorbar_legend="rho", horizontal = TRUE, 
       view_angles = c('sd_lateral_lh', 'sd_medial_lh', 'sd_medial_rh', 'sd_lateral_rh'));
while (rgl.cur() > 0) { close3d() }; file.remove(list.files(path = getwd(), pattern = 'fsbrain'))



