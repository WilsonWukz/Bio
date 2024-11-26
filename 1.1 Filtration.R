#Subgroup
Preneo$subgroup <- sub("_.*", "", colnames(Preneo))

#Count cells' amount before filtration
Preneo_before <- table(Preneo$subgroup) %>% as.data.frame()
colnames(Preneo_before) <- c("Subgroup", "CellCount_Before")

#Filtration
Preneo_filtered <- subset(Preneo, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
Preneo_filtered$subgroup <- sub("_.*", "", colnames(Preneo_filtered))

#Count cells' amount after filtration
Preneo_after <- table(Preneo_filtered$subgroup) %>% as.data.frame()
colnames(Preneo_after) <- c("Subgroup", "CellCount_After")

#Combine the counting data
Preneo_stats <- merge(Preneo_before, Preneo_after, by = "Subgroup")
Preneo_stats <- Preneo_stats %>%
  mutate(Group = "Preneo",
         Percentage = (CellCount_After / CellCount_Before) * 100)

#Subgroup
Tumor$subgroup <- sub("_.*", "", colnames(Tumor))

#Count cells' amount before filtration
Tumor_before <- table(Tumor$subgroup) %>% as.data.frame()
colnames(Tumor_before) <- c("Subgroup", "CellCount_Before")

#Filtration
Tumor_filtered <- subset(Tumor, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
Tumor_filtered$subgroup <- sub("_.*", "", colnames(Tumor_filtered))

#Count cells' amount after filtration
Tumor_after <- table(Tumor_filtered$subgroup) %>% as.data.frame()
colnames(Tumor_after) <- c("Subgroup", "CellCount_After")

#Combine the counting data
Tumor_stats <- merge(Tumor_before, Tumor_after, by = "Subgroup")
Tumor_stats <- Tumor_stats %>%
  mutate(Group = "Tumor",
         Percentage = (CellCount_After / CellCount_Before) * 100)

#Combine Preneo and Tumor
combined_stats <- rbind(Preneo_stats, Tumor_stats)

# Draw the plot
ggplot(combined_stats, aes(x = Subgroup, y = CellCount_Before, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", aes(y = CellCount_Before, fill = "Before")) +
  geom_bar(stat = "identity", position = "dodge", aes(y = CellCount_After, fill = "After")) +
  geom_text(aes(y = CellCount_After + 100, label = paste0(round(Percentage, 1), "%")), 
            position = position_dodge(0.9), size = 4) +
  facet_wrap(~Group, scales = "free_x") +
  labs(title = "Cell Count Before and After Filtration in Preneo and Tumor Subgroups",
       x = "Subgroup",
       y = "Cell Count",
       fill = "Stage") +
  theme_minimal() +
  theme(text = element_text(size = 14))
