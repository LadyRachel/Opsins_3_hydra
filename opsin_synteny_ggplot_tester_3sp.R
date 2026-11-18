#!/usr/bin/env Rscript
# opsin_synteny_ggplot_tester_3sp.R   
# TESTER: three-species opsin synteny, AEP -> HV105 -> HVIRI. HOLI removed.
#
# Background thinning is by UNIFORM SUBSAMPLING (BG_SUBSAMPLE), which lowers
# ribbon density without biasing the diagonal:off-diagonal (movement) ratio.
# Do NOT thin the background by tightening DIAMOND/alignment stringency:
# that preferentially keeps conserved syntenic genes and makes the genome look
# more conserved than it is, inflating the opsin contrast.


library(ggplot2)

# ── INPUT ────────────────────────────────────────────────────────────────────
CHR_SIZES  <- "chr_sizes.tsv"
MARKER_DIR <- "opsin_markers"
LINK_DIR   <- "opsin_links"
BG_LINK_DIR <- "synteny_links"      # whole-genome 1:1 RBH ortholog links

# ── Gray background controls ─────────────────────────────────────────────────
DRAW_BACKGROUND <- TRUE
BG_SUBSAMPLE    <- 20      # keep every Nth RBH ribbon. HONEST thinning: uniform,
                          #   preserves the movement:conserved ratio. Raise for a
                          #   sparser wash, lower (down to 1) for the full field.
BG_MOVERS_ONLY  <- FALSE   # FALSE = full background (conserved diagonal + movers).
                          #   TRUE  = draw ONLY chromosome-jumping orthologs. Much
                          #   cleaner, but drops the conserved "denominator" — the
                          #   reader can no longer see what FRACTION of the genome
                          #   stayed put. Decide this framing with Aide (see notes).
# <<< BACKGROUND DARKNESS KNOB >>> lower grey number = darker wash, higher =
# lighter. grey72 here; try grey65 (darker) or grey80 (lighter) to taste.
BG_COLOR <- "grey72"
BG_ALPHA <- 0.5
BG_LW    <- 0.12

SPECIES_ORDER  <- c("AEP", "HV105", "HVIRI")
SPECIES_LABELS <- c(AEP = "H. vulgaris AEP", HV105 = "H. vulgaris HV105",
                    HVIRI = "H. viridissima")
# Muted/pastel versions of the tree strain colors (each hex blended ~40% toward
# white). Same hue family as the phylogeny, but light enough that the
# clade-colored opsin markers read against them. Full-saturation tree hues, for
# reference: AEP #5AB4E5, HV105 #0073B2, HVIRI #E5A024.
SPECIES_FILL   <- c(AEP = "#9CD2EF", HV105 = "#66ABD1", HVIRI = "#EFC67C")
MARK_COLOR <- "#1F4E79"
LINK_COLOR <- "#1F4E79"
CLADE_COLORS <- c(CI = "#E69F00", CII = "#009E73", CIII = "#0072B2",
                  CIV = "#D62728", CV = "#7570B3")

PAIRS <- list(c("AEP","HV105"), c("HV105","HVIRI"))

# Column i is the same homologous chromosome across all rows, so dropping the
# HOLI row leaves the 15-column correspondence intact.
SYNTENIC_ORDER <- list(
  AEP = c("HVAEP1","HVAEP2","HVAEP3","HVAEP4","HVAEP5","HVAEP6","HVAEP7",
          "HVAEP8","HVAEP9","HVAEP10","HVAEP11","HVAEP12","HVAEP13","HVAEP14",
          "HVAEP15"),
  HV105 = c("NC_061159.1","NC_061161.1","NC_061158.1","NC_061167.1",
            "NC_061170.1","NC_061163.1","NC_061169.1","NC_061166.1",
            "NC_061160.1","NC_061165.1","NC_061162.1","NC_061157.1",
            "NC_061164.1","NC_061168.1","NC_061156.1"),
  HVIRI = c("NC_136456.1","NC_136448.1","NC_136447.1","NC_136454.1",
            "NC_136452.1","NC_136451.1","NC_136457.1","NC_136450.1",
            "NC_136449.1","NC_136459.1","NC_136455.1","NC_136446.1",
            "NC_136453.1","NC_136460.1","NC_136458.1")
)

# ── Layout ───────────────────────────────────────────────────────────────────
COL_W       <- 10
COL_GAP     <- 2
BAR_H       <- 0.30
ROW_SPACING <- 1.5
# ─────────────────────────────────────────────────────────────────────────────

chr_sizes <- read.delim(CHR_SIZES, stringsAsFactors = FALSE)
n_sp <- length(SPECIES_ORDER)
ycenters <- setNames((n_sp - seq_len(n_sp)) * ROW_SPACING, SPECIES_ORDER)

chrtab <- do.call(rbind, lapply(SPECIES_ORDER, function(sp) {
  ord   <- SYNTENIC_ORDER[[sp]]
  sizes <- chr_sizes[chr_sizes$species == sp, ]
  len   <- sizes$length[match(ord, sizes$accession)]
  data.frame(species = sp, accession = ord, col_i = seq_along(ord),
             length = len, stringsAsFactors = FALSE)
}))
chrtab$col_left <- (chrtab$col_i - 1) * (COL_W + COL_GAP)
chrtab$y        <- ycenters[chrtab$species]
rownames(chrtab) <- paste(chrtab$species, chrtab$accession)

xmap <- function(species, chr, bp) {
  k <- match(paste(species, chr), rownames(chrtab))
  frac <- bp / chrtab$length[k]
  over <- !is.na(frac) & (frac > 1 | frac < 0)
  if (any(over)) {
    cat(sprintf("xmap: clamped %d marker(s) past chromosome bounds:\n", sum(over)))
    bad <- data.frame(species = species[over], chr = chr[over], bp = bp[over],
                       length = chrtab$length[k][over], frac = round(frac[over], 3))
    print(bad)
    frac[over] <- pmin(pmax(frac[over], 0), 1)
  }
  chrtab$col_left[k] + frac * COL_W
}

# ── Bars ─────────────────────────────────────────────────────────────────────
bars <- data.frame(
  xmin = chrtab$col_left, xmax = chrtab$col_left + COL_W,
  ymin = chrtab$y - BAR_H / 2, ymax = chrtab$y + BAR_H / 2,
  species = SPECIES_LABELS[chrtab$species]
)

# ── Markers ──────────────────────────────────────────────────────────────────
# Marker files are per-strain and already on disk for all four strains; we just
# read the three we need. Opsins on unplaced scaffolds get x = NA (their chr is
# not in SYNTENIC_ORDER) and are dropped below.
markers <- do.call(rbind, lapply(SPECIES_ORDER, function(sp) {
  d <- read.delim(file.path(MARKER_DIR, sprintf("%s_opsin_markers.tsv", sp)),
                  stringsAsFactors = FALSE)
  d$species <- sp
  d
}))
markers$mid  <- (markers$start + markers$end) / 2
markers$x    <- xmap(markers$species, markers$chr, markers$mid)
markers$ymin <- ycenters[markers$species] - BAR_H / 2
markers$ymax <- ycenters[markers$species] + BAR_H / 2
markers <- markers[!is.na(markers$x), ]
markers$clade[is.na(markers$clade) | markers$clade == ""] <- NA

# ── Links ────────────────────────────────────────────────────────────────────
links <- do.call(rbind, lapply(PAIRS, function(p) {
  path <- file.path(LINK_DIR, sprintf("%s_vs_%s_opsin_links.tsv", p[1], p[2]))
  if (!file.exists(path)) {
    cat(sprintf("  [opsin] missing %s -- skipping this pair\n", path)); return(NULL)
  }
  d <- read.delim(path, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(NULL)
  data.frame(
    x1 = xmap(p[1], d$chr1, (d$start1 + d$end1) / 2),
    y1 = ycenters[p[1]] - BAR_H / 2,
    x2 = xmap(p[2], d$chr2, (d$start2 + d$end2) / 2),
    y2 = ycenters[p[2]] + BAR_H / 2,
    clade = d$clade)
}))
links <- links[!is.na(links$x1) & !is.na(links$x2), ]
links$clade[is.na(links$clade) | links$clade == ""] <- NA

# ── Gray background: whole-genome RBH orthologs ──────────────────────────────
# Same xmap() and column system as the opsins, so conserved orthologs draw as
# near-vertical gray ribbons and chromosome-jumps draw as gray diagonals — the
# baseline the opsin ribbons are meant to be read against.
bg_links <- NULL
if (DRAW_BACKGROUND) {
  bg_links <- do.call(rbind, lapply(PAIRS, function(p) {
    path <- file.path(BG_LINK_DIR, sprintf("%s_vs_%s_links_bg.tsv", p[1], p[2]))
    if (!file.exists(path)) {
      cat(sprintf("  [bg] missing %s -- skipping this pair\n", path)); return(NULL)
    }
    d <- read.delim(path, stringsAsFactors = FALSE)
    if (nrow(d) == 0) return(NULL)
    # Keep only links whose chromosomes are in the syntenic set (this is what
    # drops the unplaced NW_ scaffold links), then flag movers by column mismatch.
    col1 <- match(d$chr1, SYNTENIC_ORDER[[p[1]]])
    col2 <- match(d$chr2, SYNTENIC_ORDER[[p[2]]])
    keep <- !is.na(col1) & !is.na(col2)
    cat(sprintf("  [bg] %s vs %s: %d links, %d on syntenic chromosomes\n",
                p[1], p[2], nrow(d), sum(keep)))
    d <- d[keep, ]; col1 <- col1[keep]; col2 <- col2[keep]
    if (BG_MOVERS_ONLY) {
      mv <- col1 != col2
      d <- d[mv, ]
    }
    if (nrow(d) == 0) return(NULL)
    d <- d[seq(1, nrow(d), by = BG_SUBSAMPLE), ]   # uniform, unbiased thinning
    data.frame(
      x1 = xmap(p[1], d$chr1, (d$start1 + d$end1) / 2),
      y1 = ycenters[p[1]] - BAR_H / 2,
      x2 = xmap(p[2], d$chr2, (d$start2 + d$end2) / 2),
      y2 = ycenters[p[2]] + BAR_H / 2)
  }))
  if (!is.null(bg_links)) bg_links <- bg_links[!is.na(bg_links$x1) & !is.na(bg_links$x2), ]
  cat(sprintf("background ribbons drawn: %d  (subsample 1/%d, movers_only = %s)\n",
              if (is.null(bg_links)) 0L else nrow(bg_links),
              BG_SUBSAMPLE, BG_MOVERS_ONLY))
}

# ── Annotation ───────────────────────────────────────────────────────────────
sp_labels <- data.frame(x = -2.5, y = ycenters[SPECIES_ORDER],
                        label = SPECIES_LABELS[SPECIES_ORDER])
col_labels <- data.frame(
  x = (seq_len(15) - 1) * (COL_W + COL_GAP) + COL_W / 2,
  y = max(ycenters) + BAR_H / 2 + 0.45, label = as.character(seq_len(15)))

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot()

# Gray background FIRST so everything else sits on top of it.
if (DRAW_BACKGROUND && !is.null(bg_links) && nrow(bg_links) > 0) {
  p <- p + geom_curve(data = bg_links,
                      aes(x = x1, y = y1, xend = x2, yend = y2),
                      curvature = 0.03, angle = 90, ncp = 8,
                      color = BG_COLOR, alpha = BG_ALPHA, linewidth = BG_LW)
}

p <- p +
  geom_curve(data = links, aes(x = x1, y = y1, xend = x2, yend = y2,
                               color = clade),
             curvature = 0.03, angle = 90, ncp = 10,
             alpha = 0.55, linewidth = 0.35) +
  geom_rect(data = bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = species), color = "grey35", linewidth = 0.2) +
  geom_segment(data = markers, aes(x = x, xend = x, y = ymin, yend = ymax,
                                   color = clade), linewidth = 0.6) +
  geom_text(data = sp_labels, aes(x = x, y = y, label = label),
            hjust = 1, fontface = "italic", size = 3.6) +
  geom_text(data = col_labels, aes(x = x, y = y, label = label),
            size = 2.8, color = "grey45") +
  scale_fill_manual(values = setNames(SPECIES_FILL[SPECIES_ORDER],
                                      SPECIES_LABELS[SPECIES_ORDER]),
                    guide = "none") +
  scale_color_manual(values = CLADE_COLORS, na.value = "grey75",
                     name = "Clade", breaks = names(CLADE_COLORS),
                     na.translate = FALSE) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "right", plot.margin = margin(12, 24, 24, 70))
OUT_PDF <- "opsin_synteny_3species_tester.pdf"
OUT_PNG <- "opsin_synteny_3species_tester.png"

ggsave(OUT_PDF, p, width = 16, height = 4.2, device = cairo_pdf)
ggsave(OUT_PNG, p, width = 16, height = 4.2, dpi = 200, bg = "white")

cat(sprintf("markers drawn: %d   links drawn: %d\n", nrow(markers), nrow(links)))
cat(sprintf("wrote %s and %s\n", OUT_PDF, OUT_PNG))

# ── off-diagonal table ─────────────────────────────────────────────────────────────────────
opsin_cols <- do.call(rbind, lapply(PAIRS, function(p) {
  path <- file.path(LINK_DIR, sprintf("%s_vs_%s_opsin_links.tsv", p[1], p[2]))
  d <- read.delim(path, stringsAsFactors = FALSE)
  col1 <- match(d$chr1, SYNTENIC_ORDER[[p[1]]])
  col2 <- match(d$chr2, SYNTENIC_ORDER[[p[2]]])
  data.frame(pair = paste(p[1], p[2], sep = "-"),
             chr1 = d$chr1, start1 = d$start1, end1 = d$end1,
             chr2 = d$chr2, start2 = d$start2, end2 = d$end2,
             clade = d$clade, col1 = col1, col2 = col2,
             diagonal = !is.na(col1) & !is.na(col2) & col1 == col2,
             row.names = NULL, stringsAsFactors = FALSE)
}))

table(opsin_cols$pair, opsin_cols$diagonal, useNA = "ifany")
mk2 <- do.call(rbind, lapply(SPECIES_ORDER, function(sp) {
  d <- read.delim(file.path(MARKER_DIR, sprintf("%s_opsin_markers.tsv", sp)),
                  stringsAsFactors = FALSE)
  data.frame(species = sp, chr = d$chr, start = d$start, end = d$end,
             gene_id = d$gene_id, row.names = NULL, stringsAsFactors = FALSE)
}))

hit <- function(sp, chr, s, e) {
  vapply(seq_along(chr), function(i) {
    j <- which(mk2$species == sp[i] & mk2$chr == chr[i] &
                 mk2$start <= e[i] & mk2$end >= s[i])
    if (length(j) == 1) mk2$gene_id[j] else if (length(j) == 0) NA_character_ else "MULTI"
  }, character(1))
}

sp1 <- sub("-.*", "", opsin_cols$pair); sp2 <- sub(".*-", "", opsin_cols$pair)
opsin_cols$gene1 <- hit(sp1, opsin_cols$chr1, opsin_cols$start1, opsin_cols$end1)
opsin_cols$gene2 <- hit(sp2, opsin_cols$chr2, opsin_cols$start2, opsin_cols$end2)

table(gene1 = ifelse(is.na(opsin_cols$gene1), "none", opsin_cols$gene1 == "MULTI"),
      gene2 = ifelse(is.na(opsin_cols$gene2), "none", opsin_cols$gene2 == "MULTI"))

amb <- which(opsin_cols$gene1 == "MULTI" | opsin_cols$gene2 == "MULTI")

cands <- function(sp, chr, s, e) {
  j <- which(mk2$species == sp & mk2$chr == chr & mk2$start <= e & mk2$end >= s)
  ov <- pmin(e, mk2$end[j]) - pmax(s, mk2$start[j])
  paste(sprintf("%s(%dbp)", mk2$gene_id[j], ov), collapse = " | ")
}

for (i in amb) {
  r <- opsin_cols[i, ]
  side <- if (r$gene1 == "MULTI") 1 else 2
  sp <- if (side == 1) sub("-.*", "", r$pair) else sub(".*-", "", r$pair)
  chr <- if (side == 1) r$chr1 else r$chr2
  s   <- if (side == 1) r$start1 else r$start2
  e   <- if (side == 1) r$end1 else r$end2
  cat(sprintf("[%s] %s side%d %s:%d-%d clade=%s diag=%s -> %s\n",
              i, r$pair, side, chr, s, e, r$clade, r$diagonal, cands(sp, chr, s, e)))
}
cw <- read.delim("hv105_loc_crosswalk.tsv", stringsAsFactors = FALSE)
head(cw); nrow(cw); colnames(cw)

m105 <- subset(m, Strain == "Hv105")
acc2gene <- do.call(rbind, lapply(seq_len(nrow(m105)), function(i) {
  a <- trimws(strsplit(m105$Legacy_ID[i], ";")[[1]])
  data.frame(acc = a[a != ""], gene = m105$Gene_ID_phylo[i],
             row.names = NULL, stringsAsFactors = FALSE)
}))

nrow(acc2gene); length(unique(acc2gene$gene))

nine <- c("XP_012562617.2","XP_047139602.1","XP_047139603.1","XP_047139604.1",
          "NP_001274706.1","XP_012562830.2","XP_047140615.1",
          "XP_002163358.2","XP_047126525.1")
table(acc2gene$gene[match(nine, acc2gene$acc)], useNA = "ifany")
m105$key <- paste(m105$chromosome_scaffold, m105$start, m105$end, sep = ":")
mk2$key   <- paste(mk2$chr, mk2$start, mk2$end, sep = ":")
mk105 <- mk2[mk2$species == "HV105", ]
mk105$gene_phylo <- m105$Gene_ID_phylo[match(mk105$key, m105$key)]

sum(is.na(mk105$gene_phylo)); nrow(mk105); length(unique(mk105$gene_phylo))
subset(mk105, key %in% c("NC_061164.1:42824078:42824992",
                         "NC_061157.1:57638486:57639415",
                         "NC_061157.1:57822440:57823417"))

mkall <- mk2
mkall$gene_phylo <- NA_character_

for (i in seq_len(nrow(mkall))) {
  sp <- c(AEP = "HvAEP", HV105 = "Hv105", HVIRI = "Hviri")[mkall$species[i]]
  cand <- m[m$Strain == sp, ]
  j <- which(cand$chromosome_scaffold == mkall$chr[i] &
               cand$start <= mkall$end[i] & cand$end >= mkall$start[i])
  if (length(j) == 1) mkall$gene_phylo[i] <- cand$Gene_ID_phylo[j]
  else if (length(j) > 1) mkall$gene_phylo[i] <- "MULTI"
}

table(mkall$species, ifelse(is.na(mkall$gene_phylo), "none",
                            ifelse(mkall$gene_phylo == "MULTI", "multi", "ok")))

aep <- m[m$Strain == "HvAEP", ]
aep$chr_mapped <- sub("^chr", "HVAEP", aep$chromosome_scaffold)

idx <- which(mkall$species == "AEP")
for (i in idx) {
  j <- which(aep$chr_mapped == mkall$chr[i] &
               aep$start <= mkall$end[i] & aep$end >= mkall$start[i])
  mkall$gene_phylo[i] <- if (length(j) == 1) aep$Gene_ID_phylo[j]
  else if (length(j) > 1) "MULTI" else NA_character_
}

table(mkall$species, ifelse(is.na(mkall$gene_phylo), "none",
                            ifelse(mkall$gene_phylo == "MULTI", "multi", "ok")))

m$lo <- pmin(m$start, m$end); m$hi <- pmax(m$start, m$end)
strain_of <- c(AEP = "HvAEP", HV105 = "Hv105", HVIRI = "Hviri")

mkall$gene_phylo <- NA_character_
for (i in seq_len(nrow(mkall))) {
  cand <- m[m$Strain == strain_of[mkall$species[i]] & m$in_final_tree != "FALSE", ]
  cchr <- if (mkall$species[i] == "AEP")
    sub("^chr", "HVAEP", cand$chromosome_scaffold) else cand$chromosome_scaffold
  j <- which(cchr == mkall$chr[i] & cand$lo <= mkall$end[i] & cand$hi >= mkall$start[i])
  mkall$gene_phylo[i] <- if (length(j) == 1) cand$Gene_ID_phylo[j]
  else if (length(j) > 1) "MULTI" else NA_character_
}

table(mkall$species, ifelse(is.na(mkall$gene_phylo), "none",
                            ifelse(mkall$gene_phylo == "MULTI", "multi", "ok")))


lookup <- function(sp, chr, s, e) {
  vapply(seq_along(chr), function(i) {
    j <- which(mkall$species == sp[i] & mkall$chr == chr[i] &
                 mkall$start <= e[i] & mkall$end >= s[i])
    g <- unique(mkall$gene_phylo[j])
    if (length(g) == 1) g else if (!length(g)) NA_character_ else "MULTI"
  }, character(1))
}

opsin_cols$g1 <- lookup(sp1, opsin_cols$chr1, opsin_cols$start1, opsin_cols$end1)
opsin_cols$g2 <- lookup(sp2, opsin_cols$chr2, opsin_cols$start2, opsin_cols$end2)
sum(is.na(opsin_cols$g1) | opsin_cols$g1 == "MULTI",
    is.na(opsin_cols$g2) | opsin_cols$g2 == "MULTI")

off <- subset(opsin_cols, !diagonal)[, c("pair","clade","g1","col1","g2","col2")]
off[order(off$pair, off$clade, off$g1), ]
write.csv(opsin_cols, "opsin_links_labelled.csv", row.names = FALSE)