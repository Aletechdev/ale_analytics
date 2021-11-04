library(trackViewer)

UNKNOWN_EFFECT = "black"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

feats = read.csv(file = './Crr_feats.csv', stringsAsFactors=FALSE)

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]

muts = read.csv(file = './crr_pub_aa_muts.csv', stringsAsFactors=FALSE)

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

features <- GRanges("r", IRanges(feats$AA.start, end=feats$AA.end))
feats$height <- rep(0.03, length(features))
features$height = feats$height

feats$layer <- 0
feats$layer[feats$feature=="Crr 1 PFAM domain"] <- 0
feats$layer[feats$feature=="PTS Crr type-1 domain"] <- 1
feats$layer[feats$feature=="Membrane interface"] <- 1
feats$layer[feats$feature=="Helix"] <- 2
feats$layer[feats$feature=="Beta strand"] <- 2
feats$layer[feats$feature=="Turn"] <- 2
feats$layer[feats$feature=="GlpK Zinc binding site"] <- 4
feats$layer[feats$feature=="Crr activity tele-phosphohistidine intermediate"] <- 5
feats$layer[feats$feature=="Phospho-donor active site"] <- 6
feats$layer[feats$feature=="Phosphohistidine by HPr"] <- 6
feats$layer[feats$feature=="N6-acetyllysine site"] <- 6
feats$layer[feats$feature=="MalK Crr interface"] <- 20
feats$layer[feats$feature=="FrsA Crr interface"] <- 21
feats$layer[feats$feature=="PtsI Crr interface"] <- 22
feats$layer[feats$feature=="PtsG Crr interface"] <- 23
feats$layer[feats$feature=="PtsH Crr interface"] <- 24
feats$layer[feats$feature=="GlpK Crr interface"] <- 25
feats$layer[feats$feature=="Crr Crr interface"] <- 26
features$featureLayerID <- paste(feats$layer)
names(features) <- paste(feats$feature)
features$fill <- feats$color
features$color <- feats$color

legend <- list(
  labels=c(
    "GYD",
    "pgiKO_1",
    "pgiKO_2",
    "unknown effect",
    "truncation",
    "deleterious (SIFT < 0.05)",
    "structurally destabilizing (ΔΔG > 2)",
    "deleterious and structurally destabilizing"
  ),
  fill=c(
    "#AEC7E8",
    "#E377C2",
    "#9467BD",
    "white",
    "white",
    "white",
    "white",
    "white"
  ),
  col=c(
    "white",
    "white",
    "white",
    "#6e6e6e",
    CODING_SEQ_DISRUPT_COLOR,
    FUNC_DISRUPT_COLOR,
    STRUCT_DISRUPT_COLOR,
    FUNC_AND_STRUCT_DISRUPT_COLOR
  )
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$AA.start, aa_chain$AA.end)),
          yaxis = FALSE
)
