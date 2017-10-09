####### Make flowchart Figure 1

## Install and source necessary packages
rm(list = ls())

## Load and attach necessary add-on packages
library(diagram)

# Load in Rdata file created by ElNinoMetaAnalysis_code.R 
load(file="data/ElNino.RData")

## You need to set your working directory to be within the GitHub repo "ElNino_MetaAnalysis" to proceed with the following filepaths.

tiff(filename = "figures/flowchart.tiff", # Open tiff file
     width = 2250, height = 2651, units = "px", # Set tiff size
     compression = "lzw", res = 300) # Set compression and resolution
par(mar = c(0.1, 0.1, 0.1, 0.1)) # Set graphical parameters
openplotmat() # Open plot matrix
elpos <- coordinates(c(1, 1, 1, 2, 2, 3, 3, 3, 3)) # Define coordinates
fromto <- matrix(ncol = 2, byrow = TRUE, data = c(1, 2, 2, 3, 3, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9, 7, 10, 8, 11, 9, 18, 10, 13, 11, 14, 13, 16, 14, 17, 16, 19)) # Define lines connecting boxes
nr <- nrow(fromto) # Count how many connections there are
arrpos <- matrix(ncol = 2, nrow = nr) # Create matrix to place the arrows
for (i in 1:nr) arrpos[i, ] <- straightarrow (to = elpos[fromto[i, 2], ], from = elpos[fromto[i, 1], ], lwd = 3, arr.pos = 1, arr.length = 0.6) # Place arrows
textrect(elpos[1,], 0.1, 0.04, lab = paste("All 
n =",NTotal_CanAccess), cex = 1.25, font=2, family="Times", shadow.size = 0) # Add text
textrect (elpos[2,], 0.1, 0.04,lab = paste("Articles that are 
primary source of data
n =",NTotal_NoPseudo), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect(elpos[3,], 0.1, 0.04, lab = paste("El Niño/La Niña Warming
& Coral
n =", NTotal_ElNino), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[4,], 0.1, 0.04,lab = paste("Cover
n =", NTotal_Cover), cex = 1.25, font=2, family="Times", shadow.size = 0) # Add text
textrect (elpos[5,], 0.1, 0.04,lab = paste("Bleaching
n =", NTotal_Bleaching), cex = 1.25, family="Times", font=2, shadow.size = 0) # Add text
textrect (elpos[6,], 0.1, 0.04,lab = paste("Quantitative Studies
n =", NTotal_Cover_quant), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[7,], 0.1, 0.04,lab = paste("Quantitative Studies
n =", NTotal_Bleaching_quant), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[8,], 0.1, 0.04,lab = paste("Includes Data 
Before El Niño/La Niña 
n =", NTotal_CovWBef), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[9,], 0.1, 0.04,lab = paste("Includes Data Before 
El Niño/La Niña 
n =", NTotal_BlWBef), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[10,], 0.1, 0.04,lab = paste("After El Niño/La Niña 
Only (No Before Data)
n =", NTotal_BlAfterOnly), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[11,], 0.1, 0.04,lab = paste("Includes sample 
size and error
n =", NTotal_HasNandSD_cov), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[13,], 0.1, 0.04,lab = paste("Includes sample 
size and error
n =", NTotal_HasNandSD_BlAfterOnly), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[14,], 0.1, 0.04,lab = paste("Includes Percent
Cover Only
n =",NTotal_ElNino_Excel_cov_pcovonly), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[16,], 0.1, 0.04,lab = paste("Includes Percent
Bleaching Only
n =",NTotal_ElNino_Excel_BlAfterOnly_pblonly), cex = 0.75, family="Times", shadow.size = 0) # Add text
textrect (elpos[17,], 0.1, 0.04,lab = paste("n =",NTotal_ElNino_Excel_cov_LATE), cex = 1.5, font=2, family="Times", shadow.size = 0) # Add text
textrect (elpos[18,], 0.1, 0.04,lab = paste("n =",NTotal_ElNino_Excel_BlWBef_LATE), cex = 1.5, font=2, family="Times", shadow.size = 0) # Add text
textrect (elpos[19,], 0.1, 0.04,lab = paste("n =",NTotal_ElNino_Excel_BlAfterOnly_LATE), cex = 1.5, font=2, family="Times", shadow.size = 0) # Add text

dd <- c(0.0, 0.025)
text(arrpos[1, 1] + 0.165, arrpos[1, 2] + 0.055, cex = 0.75, family="Times", "Remove   Review/Meta-analysis, Secondary Lit, Repeated Data") # Add text
text(arrpos[2, 1] + 0.059, arrpos[2, 2] + 0.055, cex = 0.75, family="Times", " Remove Not   Relevant/Not El Niño/La Niña") # Add text
text(arrpos[5, 1] + 0.03, arrpos[5, 2] + 0.055, cex = 0.75, family="Times", "Remove   Qualitative Data") # Add text
text(arrpos[6, 1] + 0.03, arrpos[6, 2] + 0.055, cex = 0.75, family="Times", "Remove   Qualitative Data") # Add text
text(arrpos[7, 1] + 0.012, arrpos[7, 2] - 0.055, cex = 0.75, family="Times", "Remove   No Error/N") # Add text
text(arrpos[12, 1] + 0.012, arrpos[12, 2] + 0.055, cex = 0.75, family="Times", "Remove   No Error/N") # Add text
text(arrpos[14, 1] + 0.028, arrpos[14, 2] + 0.055, cex = 0.75, family="Times", "Remove   Density Studies") # Add text
text(arrpos[15, 1] + 0.04, arrpos[15, 2] + 0.052, cex = 0.75, family="Times", "Remove   Long-term Studies") # Add text
text(arrpos[11, 1] + 0.04, arrpos[11, 2] + 0.052, cex = 0.75, family="Times", "Remove   Long-term Studies") # Add text
text(arrpos[16, 1] + 0.04, arrpos[16, 2] + 0.052, cex = 0.75, family="Times", "Remove   Long-term Studies") # Add text
text(arrpos[13, 1] + 0.028, arrpos[13, 2] + 0.055, cex = 0.75, family="Times", "Remove   Density Studies") # Add text

dev.off() # Close tiff file to complete, and reset graphical parameters
 