library(readr)
pairs <- read_table2("F5selectednodesfilteredsam.txt", col_names = FALSE)
colnames(pairs) <- c("count", "from", "to", "flags")

skip <- sum(substr(readLines("graph.gfa"),1,1)=="S")
graph <- read_table2("graph.gfa", col_names = FALSE, skip = skip, comment = "#")


names <- unique(c(pairs$from, pairs$to))
ltr <- c(LETTERS, letters)[1:length(names)]
names(ltr) <- names


p <- data.frame(fr=ltr[pairs$from], se=ltr[pairs$to], count=pairs$count, flags=bitAnd(48, pairs$flags))
for(i in 1:nrow(p)) {
  if(as.character(p[i,"fr"]) > as.character(p[i,"se"])) {
    x <- p[i,"se"]
    p[i,"se"] <- p[i,"fr"]
    p[i,"fr"] <-  x
    p[i,"flags"] <- 48 - p[i,"flags"]
  }
}

up <- unique(p[,-3]) # no duplicated, ignring count
up$count <- sapply(1:nrow(up),
              function(i) {sum(subset(p, fr==up$fr[i] & se==up$se[i] & flags==up$flags[i])$count)})
up$pair <- paste(up$fr, up$se)
# subset(up, count > 1)



names <- stringi::stri_replace(names, "", fixed="NODE_")
names <- stringi::stri_replace(names, "", regex="[+]_length.*")
ltr2 <- ltr
names(ltr2) <- names

g <- data.frame(fr=ltr2[as.character(graph$X2)], se=ltr2[as.character(graph$X4)])
g <- subset(g, rowSums(is.na(g))==0)

for(i in 1:nrow(g)) {
  if(as.character(g[i,"fr"]) > as.character(g[i,"se"])) {
    x <- g[i,"se"]
    g[i,"se"] <- g[i,"fr"]
    g[i,"fr"] <-  x
  }
}

g$pair <- paste(g$fr, g$se)

new_pair <- subset(up, !(pair %in% g$pair) & count> 2)
names(names) <- ltr2

new_pair$fr <- names[as.character(new_pair$fr)]
new_pair$se <- names[as.character(new_pair$se)]
