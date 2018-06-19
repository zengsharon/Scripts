setwd("C:/Users/sharo/Google Drive/Savage Lab - Drive/R Files/CM/Data")
library("epade", lib.loc="~/R/win-library/3.4")
library("plyr", lib.loc="~/R/win-library/3.4")
library("questionr", lib.loc="~/R/win-library/3.4")
Naive_TCRs_weighted <- read.csv('Naive_TCRS_weighted.csv')
Naive_TCRs <- c(as.character.factor(Naive_TCRs_weighted[, 1]))
MP_TCRs_weighted <- read.csv('VM_TCRs_weighted.csv')
MP_TCRs <- as.character(MP_TCRs_weighted[, 1])

MP_residues <- data.frame(substr(MP_TCRs, 1, 1))
for (i in 2:16) {
  blah <- data.frame(substr(MP_TCRs, i, i))
  MP_residues <- cbind(MP_residues, blah)
}

Naive_residues <- data.frame(substr(Naive_TCRs, 1, 1))
for (i in 2:16) {
  blah <- data.frame(substr(Naive_TCRs, i, i))
  Naive_residues <- cbind(Naive_residues, blah)
}

Naive_weights = (Naive_TCRs_weighted$Mean_Freq/sum(Naive_TCRs_weighted$Mean_Freq))
MP_weights = (MP_TCRs_weighted$Mean_Freq/sum(MP_TCRs_weighted$Mean_Freq))

amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

barplot(wtd.table(factor(Naive_residues[, 1], levels = amino_acids), weights = Naive_weights), ylim = c(0, 1), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 1")

for (i in 2:11) {
  barplot(wtd.table(factor(Naive_residues[, i], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.5), ylab = "Frequency", col = "gray", main = paste("CD8-Naive: Residue", i))
}  

barplot(wtd.table(factor(Naive_residues[, 12], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.25), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 12")

barplot(wtd.table(factor(Naive_residues[, 13], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.25), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 13")

barplot(wtd.table(factor(Naive_residues[, 14], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.1), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 14")

barplot(wtd.table(factor(Naive_residues[, 15], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.01), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 15")

barplot(wtd.table(factor(Naive_residues[, 16], levels = amino_acids), weights = Naive_weights), ylim = c(0, 0.001), ylab = "Frequency", col = "gray", main = "CD8-Naive: Residue 16")

barplot(wtd.table(factor(MP_residues[, 1], levels = amino_acids), weights = MP_weights), ylim = c(0, 1), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 1")

for (i in 2:11) {
  barplot(wtd.table(factor(MP_residues[, i], levels = amino_acids), weights = MP_weights), ylab = "Frequency", ylim = c(0, 0.5), col = "gray", main = paste("CD8-MP: Residue", i))
}  

barplot(wtd.table(factor(MP_residues[, 12], levels = amino_acids), weights = MP_weights), ylim = c(0, 0.25), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 12")

barplot(wtd.table(factor(MP_residues[, 13], levels = amino_acids), weights = MP_weights), ylim = c(0, 0.25), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 13")

barplot(wtd.table(factor(MP_residues[, 14], levels = amino_acids), weights = MP_weights), ylim = c(0, 0.1), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 14")

barplot(wtd.table(factor(MP_residues[, 15], levels = amino_acids), weights = MP_weights), ylim = c(0, 0.01), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 15")

barplot(wtd.table(factor(MP_residues[, 16], levels = amino_acids), weights = MP_weights), ylim = c(0, 0.001), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 16")

# barplot(wtd.table(Naive_residues[, 7], weights = Naive_TCRs_weighted$Mean_Freq)[-1], ylab = "Frequency (%)", col = "gray", main="CD8-Naive: Residue 7", ylim = c(0, 40), axes = TRUE)

# bar.plot.wtd(substr(Naive_TCRs, 7, 7), w = Naive_TCRs_weighted[, 2], ylab = "Frequency (%)", col = "gray", main=paste("CD8-Naive: Residue" , 7), vnames.x = c("<7", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

# bar.plot.wtd(substr(Naive_TCRs, 11, 11), w = Naive_TCRs_weighted[, 2], ylab = "Frequency (%)", col = "gray", main = "CD8-Naive: Residue 11", vnames.x = c("<11", "A", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))

# for (i in 1:16) {
#   png(file = paste("CD8-Naive Residue", i, sep = " ", ".png"), bg = "transparent")
#   bar.plot.wtd(substr(Naive_TCRs, i, i), w = Naive_TCRs_weighted[, 2], ylab = "Frequency (%)", col = "gray", main=paste("CD8-Naive: Residue" , i))
#   dev.off()
# }
# for (i in 1:16) {
#   png(file = paste("CD8-MP Residue", i, sep = " ", ".png"), bg = "transparent")
#   bar.plot.wtd(substr(MP_TCRs, i, i), w = MP_TCRs_weighted[, 2], ylab = "Frequency (%)", col = "gray", main=paste("CD8-MP: Residue" , i))
#   dev.off()
# }

# barplot(wtd.table(factor(MP_residues[, 16], levels = amino_acids), weights = MP_weights)[-1], ylim = c(0, 0.001), ylab = "Frequency", col = "gray", main = "CD8-MP: Residue 16")


