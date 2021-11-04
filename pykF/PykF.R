library(trackViewer)

UNKNOWN_EFFECT = "black"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

# Read in feats features and remove those not interested in
feats = read.csv(file = './PykF_feats.csv', stringsAsFactors=FALSE)

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]

muts = read.csv(file = './pykF_pub_aa_muts.csv', stringsAsFactors=FALSE)

unique(muts$study)


muts$clr <- "black"
muts$clr <- muts$color

mutations <- GRanges("r", IRanges(muts$AA.pos, width=1, names=muts$name))
muts$border_color <- muts$clr
muts$border_color[muts$border_color==UNKNOWN_EFFECT] <- "#6e6e6e"
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
feats$layer[feats$feature=="Helix"] <- 0
feats$layer[feats$feature=="Beta strand"] <- 0
feats$layer[feats$feature=="Turn"] <- 0
feats$layer[feats$feature=="Barrel domain"] <- 1
feats$layer[feats$feature=="Alpha/beta domain"] <- 1
feats$layer[feats$feature=="Transition state stabilizer"] <- 2
feats$layer[feats$feature=="N6-acetyllysine site"] <- 2
feats$layer[feats$feature=="Potassium binding site"] <- 3
feats$layer[feats$feature=="ATP binding site"] <- 4
feats$layer[feats$feature=="Magnesium binding site"] <- 4
feats$layer[feats$feature=="Substrate binding site"] <- 5
feats$layer[feats$feature=="PykF subunit interface"] <- 6

features$featureLayerID <- paste(feats$layer)
names(features) <- paste(feats$feature)
features$fill <- feats$color
features$color <- feats$color


legend <- list(labels=
                 c(
                   "TOL_isobutyric_acid",
                   "tpiAKO",
                   "TOL_putrescine",
                   "tpiHSA",
                   "tpiPAE",
                   "SER",
                   "SSW_GLU_AC",
                   "42C",
                   "TOL_hexanoic_acid",
                   "pyocyanin",
                   "SSW_GLU_XYL",
                   "unknown effect",
                   "truncation",
                   "deleterious (SIFT < 0.05)",
                   "structurally destabilizing (ΔΔG > 2)",
                   "deleterious and structurally destabilizing"
                   ),
               fill=c(
                 "#9467BD",
                 "#F7B6D2",
                 "#C5B0D5",
                 "#2AB0C1",
                 "#B2B35A",
                 "#2CA02C",
                 "#97DE94",
                 "#AEC7E8",
                 "#FF9896",
                 "#E377C2",
                 "#B34142",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white"),
               col=c(
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "white",
                 "#6e6e6e",
                 CODING_SEQ_DISRUPT_COLOR,
                 FUNC_DISRUPT_COLOR,
                 STRUCT_DISRUPT_COLOR,
                 FUNC_AND_STRUCT_DISRUPT_COLOR)
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$start, aa_chain$end)),
          yaxis = FALSE
)
# grid.text("Mutations across PykF's amino acid chain", x=.5, y=0.98, just="top", gp=gpar(cex=1.5, fontface="bold"))
