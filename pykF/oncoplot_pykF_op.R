library(ComplexHeatmap)

mat = read.csv(file = './pykF_op_oncoprint_mat_genes.csv', header=TRUE, row.names="X")
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
  INS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.70, gp = gpar(fill = "#76b7b2", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "#F58518", col = NA))
  },
  SNP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#4C78A8", col = NA))
  }
)
get_type_fun = function(x) strsplit(x, ";")[[1]]

cond_mat = read.csv(file = './pykF_op_oncoprint_mat_conditions.csv', header=TRUE, row.names="X")
cond_mat = as.matrix(cond_mat)

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


ha = HeatmapAnnotation(
  temperature = temperature,
  carbon_source = carbon_source,
  strain = strain,
  supplement = supplement,
  sulfur_source = sulfur_source,
  nitrogen_source = nitrogen_source,
  phosphorous_source = phosphorous_source,
  calcium_source = calcium_source,
  study = study,
  col = list(
    study = c(
      # "10.1093/molbev/msu209"="#4C78A8",
      # "10.1128/AEM.02246-14"="#F58518",
      # "10.1016/j.ymben.2016.11.008"="#E45756",
      # "10.1128/AEM.00410-17"="#72B7B2",
      # "10.1101/634105"="#54A24B",
      # "10.1371/journal.pgen.1005715"="#EECA3B",
      # "10.1016/j.ymben.2018.05.012"="#B279A2",
      # "10.1038/s41559-020-1271-x"="#FF9DA6",
      # "10.3389/fmicb.2018.00427"="#9D755D"
      "42C"="#AEC7E8",
      "GLU"="#FCB887",
      "SER"="#2CA02C",
      "SSW_GLU_AC"="#97DE94",
      "SSW_GLU_XYL"="#B34142",
      "TOL_hexanoic_acid"="#FF9896",
      "TOL_isobutyric_acid"="#9467BD",
      "TOL_putrescine"="#C5B0D5",
      "pyocyanin"="#E377C2",
      "tpiAKO"="#F7B6D2",
      "tpiHSA"="#2AB0C1",
      "tpiPAE"="#B2B35A",
      "gentamycin"="#D7D7A8"
    ),
    temperature = c("42 Celsius" = "#E74C3C",
                    "37 Celsius" = "#E67E23"),
    carbon_source = c(
      "glucose(2)" = "#e9e9e9",
      "glucose(3)" = "#c6c6c6",
      "glucose(4)" = "#959595",
      "glucose(4) or acetate(4)" = "#686868",
      "glucose(4) or xylose(4)" = "#333333"
                      ),
    strain = c(
      "H. sapien tpiA" = "#e1d4e8",
      "P. aerophilum tpiA" = "#cda0cd",
      "WT" = "#df64af",
      "ΔsdaA ΔsdaB ΔtdcG ΔglyA" = "#df2179",
      "ΔtpiA glucose M9 pre-evolved" = "#a90649"
               ),
    supplement = c(
      "C13H10N2O(0.168184)" = "#deebf7",
      "gentamycin(0.03)" = "#c6dbef",
      "glycine(2mM) & L-serine" = "#9dcae1",
      "hexanoic acid" = "#6aaed6",
      "isobutyric acid" = "#4191c6",
      "NaCl(0.5) trace elements" = "#2070b4",
      "putrescine" = "#08509b"
                   ),
    sulfur_source = c(
      "MgSO4(0.24)" = "#F69B58",
      "(NH4)2SO4(2) MgSO4(0.12)" = "#AC501A"
      ),
    nitrogen_source = c(
      "NH4Cl(1)" = "#F2735D",
      "(NH4)2SO4(2)" = "#C03029"
    ),
    phosphorous_source = c(
      "KH2PO4(3) Na2HPO4(6.8)" = "#7EC686",
      "KH2PO4(13.6)" = "#368A4E"
    ),
    calcium_source = c(
      "CaCl2(0.1)" = "#7BC4C3",
      "none" = "#3B6C7F"
      )
    ),
  annotation_label = c(
    "temperature",
    "carbon source",
    "strain",
    "supplement",
    "sulfur source",
    "nitrogen source source",
    "phosphorous source",
    "calcium source",
    "study"),
  annotation_name_side = "left"
)

oncoplot_col = c("SNP" = "#4C78A8", "DEL" = "#F58518", "INS"="#76b7b2", "MOB"="#E45756")
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
                              at = c("SNP", "DEL", "MOB", "INS"),
                              labels = c("single nucleotide polymorphism (SNP)", "deletion (DEL)", "mobile insertion element (MOB)", "insertion (INS)"))
)

draw(p,
     merge_legend = TRUE,
     # heatmap_legend_side = "bottom"
     )
