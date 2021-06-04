library(trackViewer)

UNKNOWN_EFFECT = "#ABABAB"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

feats = read.csv(file = './glpK_pub_aa_feats.csv')

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]

muts = read.csv(file = './glpK_pub_aa_muts.csv')

muts$clr <- "black"
muts$clr <- muts$color

mutations <- GRanges("r", IRanges(muts$AA.pos, width=1, names=muts$name))
muts$border_color <- muts$clr
muts$border_color[muts$border_color==UNKNOWN_EFFECT] <- "black"
mutations$score <- muts$mutation.count
mutations$score <- muts$mutation.count
mutations$border <- muts$border_color
muts$text_colors = muts$clr
muts$text_colors[muts$text_colors==UNKNOWN_EFFECT] <- "black"
mutations$label.parameter.gp <- gpar(col=muts$text_colors)
mutations$dashline.col <- muts$clr
mutations$color <- muts$study.color

features <- GRanges("r", IRanges(feats$start, end=feats$end))
feats$height <- rep(0.03, length(features))
features$height = feats$height
# 
feats$layer <- 0
feats$layer[feats$feature=="Carbohydrate kinase N-terminal domain"] <- 0
feats$layer[feats$feature=="Carbohydrate kinase C-terminal domain"] <- 0
feats$layer[feats$feature=="Helix"] <- 1
feats$layer[feats$feature=="Beta strand"] <- 1
feats$layer[feats$feature=="Turn"] <- 1
feats$layer[feats$feature=="GlpK subunit binding site"] <- 2
feats$layer[feats$feature=="GlpK subunit interface"] <- 3
feats$layer[feats$feature=="Crr interface"] <- 4
feats$layer[feats$feature=="Crr and Zinc binding site"] <- 5
feats$layer[feats$feature=="Allosteric FBP inhibitor binding site"] <- 5
feats$layer[feats$feature=="Substrate binding site"] <- 5
feats$layer[feats$feature=="ATP binding site"] <- 6
feats$layer[feats$feature=="N6-malonyllysine site"] <- 6
features$featureLayerID <- paste(feats$layer)
names(features) <- paste(feats$feature)
features$fill <- feats$color
features$color <- feats$color

legend <- list(
  labels=c(
    "GYD",
    "SSW_GLY",
    "SSW_GLU_GLY",
    "unknown effect",
    "deleterious (SIFT < 0.05)",
    "structurally destabilizing (ΔΔG > 2)",
    "deleterious and structurally destabilizing"
  ),
  fill=c(
    "#AEC7E8",
    "#2CA02C",
    "#97DE94",
    "white",
    "white",
    "white",
    "white"
  ),
  col=c(
    "white",
    "white",
    "white",
    "black",
    FUNC_DISRUPT_COLOR,
    STRUCT_DISRUPT_COLOR,
    FUNC_AND_STRUCT_DISRUPT_COLOR
    )
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$start, aa_chain$end)),
          yaxis = FALSE
)
# grid.text("Mutations across GlpK's amino acid chain", x=.5, y=0.99, just="top", gp=gpar(cex=1.5, fontface="bold"))
