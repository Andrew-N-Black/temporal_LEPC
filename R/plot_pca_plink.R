library(tidyverse)

eigenvec <- read.table("joint_germline.pca.eigenvec", header = FALSE)
eigenval <- scan("joint_germline.pca.eigenval")

colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigenvec) - 2)))
pve <- round(eigenval / sum(eigenval) * 100, 1)

ggplot(eigenvec, aes(PC1, PC2)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pve[1], "%)"),
    y = paste0("PC2 (", pve[2], "%)")
  ) +
  theme_minimal()
