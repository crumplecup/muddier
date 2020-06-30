# assign workspace
setwd('/home/crumplecup/work')

# load required packages
library(roxygen2)
library(data.table)
library(muddier)
library(magrittr)
library(parallel)

# set cores for parallel processing
options(mc.cores=detectCores())

# load sediment residence time data
char <- fread('charcoal.csv')
setwd('/home/crumplecup/work/muddier')
facies <- char[1,-c(1:2,74,129,339,376:377)]

# strip out stephens summary rows, keeping pmf values
pmfs <- char[-c(1,4432:nrow(char)),-c(1:2,74,129,339,376:377)]
# replace character entries with zero
pmfs[pmfs=='' | pmfs=='out of range' | pmfs=='(note: there is no DFK_84a)'] <- 0

# convert pmfs to numeric array
npmfs <- apply(pmfs,2,function(x) as.numeric(x))
npmfs <- npmfs*5		#values are annual values, but for 5 year increments

# convert years (row index of pmfs array) to numeric vector
years <- char[-c(1,4432:nrow(char)),2] %>% unlist %>% as.numeric
# mean expected age of samples
mns <- apply(npmfs,2,function(x) weighted.mean(years,x))
# record site ids
site_id <- names(pmfs)

rank <- order(mns)


# subset inherited age samples-

# convert sample id to character from factor
all_letters <- c(letters,letters %>% toupper)

# strip the last letter off of ids with ending letters
ids <- 0
for (i in 1:length(site_id)) {
  lab <- site_id[i] %>% strsplit('')   #split name into vector of single letters
  lab <- lab[[1]]     #extract character vector from list
  #if there is a letter on the end ignore it
  if (lab[length(lab)] %in% all_letters) ids[i] <- paste(lab[1:(length(lab)-1)], collapse='')
  #else copy the full sample id
  else ids[i] <- paste(lab, collapse='')
}


# build site by var matrix for selecting sites based on characteristics

charcoal <- data.table(site_id = site_id, family = ids,
      facies = unlist(facies), mn = mns)
charcoal <- charcoal[order(mns)]
charcoal$rank <- 1:nrow(charcoal)
rownames(charcoal) <- site_id[order(mns)]

setwd('/home/crumplecup/work/muddier')
usethis::use_data(charcoal, overwrite = T)

#reorder pmfs to match charcoal
char_pmfs <- array(0,dim(npmfs))
for (i in 1:ncol(npmfs))  {
  char_pmfs[,i] <- npmfs[,names(pmfs) == charcoal$site_id[i]]
}
char_pmfs <- as.data.table(char_pmfs)
colnames(char_pmfs) <- charcoal$site_id
rownames(char_pmfs) <- as.character(years)
usethis::use_data(char_pmfs, overwrite = T)


# define test var = 'more than 3 samples per site'
inher_age <- charcoal[,.(test = .N > 3, site_id, rank), keyby = .(facies, family)]
# subset sites that pass test
inher_age <- inher_age[test == 'TRUE', .(site_id, rank, count = .N), keyby = .(facies, family)]
usethis::use_data(inher_age)


# change rownames of char_pmfs to change from year 2000 instead of 1950
library(magrittr)
new_pmfs <- char_pmfs
rownames(char_pmfs) <- seq(0, (nrow(new_pmfs)*5)-1, 5) %>% as.character %>% rev
setwd('/home/crumplecup/work/muddier')
usethis::use_data(char_pmfs, overwrite = T)
