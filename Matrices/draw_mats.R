calc_score <- function(a, gap=-2, mismatch=-1, match=1) {
  sc <- matrix(0, nrow(a), ncol(a))
  a <- ifelse(a==1, match, mismatch)
  for(i in 2:nrow(a)) {
    for(j in 2:ncol(a)) {
      sc[i,j] <- max(0,sc[i-1,j-1]+a[i,j], sc[i,j-1]+gap, sc[i-1,j]+gap)
    }
  }
  return(max(sc))
}

draw_one <- function(fname) {
  m <- read.csv(fname, header=FALSE)
  n <- as.matrix(m[-1,-1])
  a <- matrix(as.numeric(n), nrow(n))
  image(a, breaks=(0:2)/2, col=c("white","black"))
  title(fname)
}

draw_pair <-function(r,c) {
  fname <- paste0(r,"vs",c,".csv")
  transpose <- FALSE
  if(!file.exists(fname)) {
    fname <- paste0(c,"vs",r,".csv")
    transpose <- TRUE
  }
  m <- read.csv(fname, header=FALSE)
  n <- as.matrix(m[-1,-1])
  a <- matrix(as.numeric(n), nrow(n))
  if(transpose) {
    a <- t(a)
  }
  image(a, breaks=(0:2)/2, col=c("white","black"),
        xlab=r, ylab=c)
  title(calc_score(a))
}


fnames <- dir(pattern = "*csv")
fnames <- grep("KJ170699", fnames, invert=TRUE, value=TRUE)
ref_names <- unique(c(substr(fnames,1,8),substr(fnames,11,18)))


pdf(file="3x3.pdf", width=8, height=8)
par(mfcol=c(3,3), mar=c(4, 4, 3, 1))
for(i in 1:3) {
  for(j in 1:3) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 4:6) {
  for(j in 4:6) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 7:9) {
  for(j in 7:9) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 1:3) {
  for(j in 4:6) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 1:3) {
  for(j in 7:9) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 4:6) {
  for(j in 7:9) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
dev.off()


ref_names <- ref_names[ref_names!="LR025097"]

pdf(file="4x4.pdf", width=8, height=8)
par(mfcol=c(4,4), mar=c(4, 4, 3, 1))
for(i in 1:4) {
  for(j in 1:4) {
      draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 5:8) {
  for(j in 5:8) {
      draw_pair(ref_names[i], ref_names[j])
  }
}
for(i in 1:4) {
  for(j in 5:8) {
    draw_pair(ref_names[i], ref_names[j])
  }
}
dev.off()