library(trackViewer)

UNKNOWN_EFFECT = "#ABABAB"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

feats = read.csv(file = './cyaA_aa_feats.csv', stringsAsFactors=FALSE)

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]

muts = read.csv(file = './cyaA_aa_muts.csv', stringsAsFactors=FALSE)

muts$clr <- "black"
muts$clr <- muts$color
mutations <- GRanges("r", IRanges(muts$AA.pos, width=1, names=muts$name))
muts$border_color <- muts$clr
muts$border_color[muts$border_color==UNKNOWN_EFFECT] <- "#6e6e6e"
mutations$score <- muts$mutation.count
mutations$score <- muts$mutation.count
mutations$border <- muts$border_color
muts$text_colors = muts$clr
muts$text_colors[muts$text_colors==UNKNOWN_EFFECT] <- "#6e6e6e"
mutations$label.parameter.gp <- gpar(col=muts$text_colors)
mutations$dashline.col <- muts$clr
mutations$color <- muts$study.color

features <- GRanges("r", IRanges(feats$start, end=feats$end))
feats$height <- rep(0.03, length(features))
features$height = feats$height
feats$layer <- 0
feats$layer[feats$feature=="Catalytic region"] <- 0
feats$layer[feats$feature=="Regulatory region"] <- 0
feats$layer[feats$feature=="Nucleotidyltransferase PFAM domain"] <- 1
feats$layer[feats$feature=="Adenylate cyclase class-I PFAM domain"] <- 1
feats$layer[feats$feature=="G3P associated inhibition"] <- 2
feats$layer[feats$feature=="Phosphohistidine by Crr"] <- 2
features$featureLayerID <- feats$layer
names(features) <- feats$feature
features$fill <- feats$color
features$color <- feats$color

legend <- list(
  labels=
    c(
      "m-tartrate2",
      "GYD",
      "pgiBME",
      "pgiKO_1",
      "pgiKO_2",
      "ptsHI-crrKO",
      "pgiHSA",
      "unknown effect",
      "truncation",
      "deleterious (SIFT < 0.05)"
    ),
  fill=c(
    "#C5B0D5",
    "#AEC7E8",
    "#F7B6D2",
    "#E377C2",
    "#9467BD",
    "#C7C7C7",
    "#2AB0C1",
    "white",
    "white",
    "white"),
  col=c("white",
        "white",
        "white",
        "white",
        "white",
        "white",
        "white",
        "#6e6e6e",
        CODING_SEQ_DISRUPT_COLOR,
        FUNC_DISRUPT_COLOR)
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$start, aa_chain$end)),
          yaxis = FALSE
)
# grid.text("Mutations across CyaA's amino acid chain", x=.5, y=0.7, just="top", gp=gpar(cex=1.5, fontface="bold"))
