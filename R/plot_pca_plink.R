library(tidyverse)
library(dplyr)

eigenvec <- read.table("joint_germline.pca.eigenvec", header = FALSE)
eigenval <- scan("joint_germline.pca.eigenval")

colnames(eigenvec) <- c("IID", paste0("PC", 1:(ncol(eigenvec) - 1)))
pve <- round(eigenval / sum(eigenval) * 100, 1)


eigenvec <- eigenvec %>%
    mutate(ID = sub("^normal_", "", IID)) %>%
    left_join(metadata %>% select(ID, GRP), by = "ID")


ggplot(eigenvec, aes(PC1, PC2)) +
    geom_point(size = 6, shape = 21, color = "black", aes(fill = GRP)) +
    labs(
        x = paste0("PC1 (", pve[1], "%)"),
        y = paste0("PC2 (", pve[2], "%)")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual("Group", values = c("cadetblue", "black")) +
    theme(legend.position = "right") +
    guides(
        fill = guide_legend(override.aes = list(shape = 21, size = 5, stroke = 0.5))
    )+theme_classic()
