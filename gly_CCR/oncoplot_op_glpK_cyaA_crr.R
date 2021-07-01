library(ComplexHeatmap)

mat = read.csv(file = './glpK_cyaA_crr_op_oncoprint_mat_gly_genes.csv', header=TRUE, row.names="X")
mat = as.matrix(mat)
# Using "Table 10" colors
alter_fun_list = list(
  # The order in which these are defined is the order in which their corresponding rects are drawn within the plot.
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  MOB = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E45756", col = NA))
  },
  # INS = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.70, gp = gpar(fill = "#76b7b2", col = NA))
  # },
  CNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.70, "mm"), gp = gpar(fill = "#B279A2", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "#F58518", col = NA))
  },
  SNP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#4C78A8", col = NA))
  }
)
get_type_fun = function(x) strsplit(x, ";")[[1]]


cond_mat = read.csv(file = './glyK_cyaA_crr_op_oncoprint_mat_conditions.csv', header=TRUE, row.names="X")
cond_mat = as.matrix(cond_mat)
# cond_mat <- cond_mat.frame(lapply(cond_mat, function(x){
#   gsub("pre-evolved", "evolved", x)
# }))

study = c()
temperature = c()
carbon_source = c()
strain = c()
supplement = c()
sulfur_source = c()
nitrogen_source = c()
phosphorous_source = c()
calcium_source = c()

for (sample_name in colnames(mat)) {
  study = append(study, cond_mat["exp id", sample_name])
  temperature = append(temperature, cond_mat["temperature", sample_name])
  carbon_source = append(carbon_source, cond_mat["carbon-source", sample_name])
  strain = append(strain, cond_mat["strain-description", sample_name])
  supplement = append(supplement, cond_mat["supplement", sample_name])
  sulfur_source = append(sulfur_source, cond_mat["sulfur-source", sample_name])
  nitrogen_source = append(nitrogen_source, cond_mat["nitrogen-source", sample_name])
  phosphorous_source = append(phosphorous_source, cond_mat["phosphorous-source", sample_name])
  calcium_source = append(calcium_source, cond_mat["calcium-source", sample_name])
}

cond_color_mat = read.csv(file = '../data/condition_clustermap_colors.csv', header=TRUE, row.names="X")

temperature_colors = c()
uniq_conds = unique(temperature)
cond_type = "temperature"
for (x in uniq_conds) {
  color = cond_color_mat[with(cond_color_mat, condition.category == cond_type & condition == x), ]$"color"[1]
  color = toString(color)
  temperature_colors = append(temperature_colors, color)
}
names(temperature_colors) = uniq_conds

library(stringr)
# carbon_source_colors = c('#ededed', '#d1d1d1', '#adadad', '#828282', '#5c5c5c', '#2b2b2b')  # Don't want to use original clustermap colors
carbon_source_colors = c('#ededed', '#d1d1d1', '#adadad', '#828282', '#5c5c5c')  # Don't want to use original clustermap colors
uniq_conds = unique(carbon_source)
uniq_conds= str_sort(uniq_conds)  # get the strains in ascending color order
names(carbon_source_colors) = uniq_conds

library(stringr)
strain_colors = c('#e7e1ef', '#d4b9da', '#c993c7', '#df64af', '#e72989', '#cd1256', '#970042')  # Don't want to use original clustermap colors
uniq_conds = unique(strain)
uniq_conds = str_sort(uniq_conds)  # get the strains in ascending color order
names(strain_colors) = uniq_conds

library(stringr)
# supplement_colors = c('#c6dbef', '#6aaed6', '#2070b4')
supplement_colors = c('#c6dbef')
uniq_conds = unique(supplement)
uniq_conds = str_sort(uniq_conds)  # get the strains in ascending color order
names(supplement_colors) = uniq_conds

sulfur_source_colors = c()
uniq_conds = unique(sulfur_source)
cond_type = "sulfur source"
for (x in uniq_conds) {
  color = cond_color_mat[with(cond_color_mat, condition.category == cond_type & condition == x), ]$"color"[1]
  color = toString(color)
  sulfur_source_colors = append(sulfur_source_colors, color)
}
names(sulfur_source_colors) = uniq_conds

nitrogen_source_colors = c()
uniq_conds = unique(nitrogen_source)
cond_type = "nitrogen source"
for (x in uniq_conds) {
  color = cond_color_mat[with(cond_color_mat, condition.category == cond_type & condition == x), ]$"color"[1]
  color = toString(color)
  nitrogen_source_colors = append(nitrogen_source_colors, color)
}
names(nitrogen_source_colors) = uniq_conds

phosphorous_source_colors = c()
uniq_conds = unique(phosphorous_source)
cond_type = "phosphorous source"
for (x in uniq_conds) {
  color = cond_color_mat[with(cond_color_mat, condition.category == cond_type & condition == x), ]$"color"[1]
  color = toString(color)
  phosphorous_source_colors = append(phosphorous_source_colors, color)
}
names(phosphorous_source_colors) = uniq_conds

calcium_source_colors = c()
uniq_conds = unique(calcium_source)
cond_type = "calcium source"
for (x in uniq_conds) {
  color = cond_color_mat[with(cond_color_mat, condition.category == cond_type & condition == x), ]$"color"[1]
  color = toString(color)
  calcium_source_colors = append(calcium_source_colors, color)
}
names(calcium_source_colors) = uniq_conds


# Matplotlib Tab20 colormap used for experiment colors.
# Getting colors of conditions from original condition clustermap for this
ha = HeatmapAnnotation(
  annotation_label = c(
    "sulfur source",
    "nitrogen source",
    "phosphorous source",
    "calcium source",
    "temperature",
    "carbon source",
    "strain",
    "supplement",
    "study"),
  sulfur_source = sulfur_source,
  nitrogen_source = nitrogen_source,
  phosphorous_source = phosphorous_source,
  calcium_source = calcium_source,
  temperature = temperature,
  carbon_source = carbon_source,
  strain = strain,
  supplement = supplement,
  study = study,
  col = list(
    sulfur_source = sulfur_source_colors,
    nitrogen_source = nitrogen_source_colors,
    phosphorous_source = phosphorous_source_colors,
    calcium_source = calcium_source_colors,
    temperature = temperature_colors,
    carbon_source = carbon_source_colors,
    strain = strain_colors,
    supplement = supplement_colors,
    study = c(
      "GYD"="#AEC7E8",
      "SSW_GLY"="#2CA02C",
      "SSW_GLU_GLY"="#97DE94",
      "SSW_GLU_XYL"="#B34142",
      "GLU"="#FCB887",
      "m-tartrate2"="#C5B0D5",
      "pgiKO_1"="#E377C2",
      "pgiKO_2"="#9467BD",
      "pgiBME"="#F7B6D2",
      "pgiHSA"="#2AB0C1",
      "pgiPAE"="#B2B35A",
      "gentamycin"="#D7D7A8",
      "TOL_hexamethylenediamine"="#C49C94",
      "ptsHI-crrKO"="#C7C7C7")
    ),
  annotation_name_side = "left"
)

oncoplot_col = c("SNP" = "#4C78A8",
                 "DEL" = "#F58518",
                 "MOB"="#E45756",
                 "CNV"="#B279A2"
                 )
options(repr.plot.width = 7, repr.plot.height = 2, repr.plot.res = 200)
p = oncoPrint(
  mat,
  alter_fun = alter_fun_list,
  get_type=get_type_fun,
  col = oncoplot_col,
  remove_empty_columns = TRUE,
  pct_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  pct_side = "right",
  row_names_side = "left",
  top_annotation = ha,
  # bottom_annotation = ha,
  right_annotation = NULL,
  heatmap_legend_param = list(title = "mutation type",
                              at = c("SNP", "DEL", "MOB", "CNV"),
                              labels = c("single nucleotide polymorphism (SNP)",
                                         "deletion (DEL)",
                                         "mobile insertion element (MOB)",
                                         "copy number variant (CNV)"))
)

draw(p,
     merge_legend = TRUE,
     # heatmap_legend_side = "bottom"
)
