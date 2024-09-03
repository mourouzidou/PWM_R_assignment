#install.packages("stringdist")
library(stringdist)
library(ggseqlogo)
library(gplots)
library(ggplot2)
library(seqinr)




#read the downloaded file containig the IDs of all human genes
geneIDs <- read.table("allgeneids.txt", sep = "\t")
IDs <- unique(geneIDs$V1)

#pick 1000random genes and save it to a file to get the entries from Ensembl
random_gene_ids <- unique(sample(IDs, size = 1000, replace = FALSE))

write.table(random_gene_ids, "rand1000IDs.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)


glycolysis <- "glycolysis.txt"
gluconeogenesis <- "gluconeogeneis.txt"
rand1000genes <- "1000randgenes.txt"

#readLines for every file-sample
glyc_lines <- readLines(glycolysis)
glucon_lines <- readLines(gluconeogenesis)
all_lines <- readLines(rand1000genes)

# find the gene IDs fir each pathway
glycolysisIDs <- grep("^>", glyc_lines, value = TRUE)
gluconeogenesisIDs <- grep("^>", glucon_lines, value = TRUE)

commonGenes <- intersect(gluconeogenesisIDs, glycolysisIDs)

ngly <- length(glycolysisIDs)
nglu <- length(gluconeogenesisIDs)
sample_size <- ngly + nglu

#rnd67gly <- sample(IDs, size = ngly)
#rnd35glu <- sample(IDs, size = nglu)
#rndCommons<- intersect(rnd35, rnd67gly)


#calculate the p-value to check the significance of common occurancies
common_prob <- phyper(length(commonGenes)-1, nglu, length(IDs)-nglu, ngly, lower.tail = FALSE)

format(common_prob, scientific = FALSE)

getData = function(filepath) {
  con = file(filepath, "r")
  data = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(length(grep(">", line)) > 0){
      name = gsub("^>(\\w+)", replacement="\\1", x=line)
      data[[name]] = ""
    }
    else{
      data[[name]] = paste(data[[name]], line, sep="")
    }
  }
  close(con)
  return(data)
}

glycolysis_data <- getData(glycolysis)
gluconeogenesis_data <- getData(gluconeogenesis)
background <- getData(rand1000genes)




#lapply(gluconeogenesis_data, nchar)
#lapply(glycolysis_data, nchar)
#lapply(background, nchar)


processSequences = function(seqList, len=5){
  x = sapply(seqList, function(i){
    v = strsplit(i,"")[[1]]
    sapply(1:(length(v)-len+1), function(j){paste(v[j:(j+len-1)], collapse="")})
  })
  table(x)
}

motifLength=8

fgGluconeogenesisDataset <- processSequences(gluconeogenesis_data, motifLength)
fgGlycolysisDataset <- processSequences(glycolysis_data, motifLength)
bgDataset <- processSequences(background, motifLength)
bgDataset

fgGluconeogenesisDataset
fgGlycolysisDataset

sort(fgGluconeogenesisDataset, decreasing=TRUE)
sort(fgGlycolysisDataset, decreasing=TRUE)


fgGlycolCounts <- sum(fgGlycolysisDataset)
fgGluconCounts <- sum(fgGluconeogenesisDataset)
bgCounts <- sum(bgDataset)

getProb = function(foreground, background){
  probs = vector("numeric", length=length(foreground))
  sumforground = sum(foreground)
  sumbackground = sum(background)
  for(i in 1:length(foreground)){
    bcounts = 0
    if( names(foreground)[i] %in% names(background)){
      bcounts = background[[names(foreground)[i]]]
      prob = phyper(q=foreground[i]-1,  m = bcounts, n = sumbackground - bcounts, k = sumforground, lower.tail = FALSE)
      probs[i] = prob
    }else{
      bcounts = 1
      prob = phyper(q=foreground[i]-1,  m = bcounts, n = sumbackground - bcounts, k = sumforground, lower.tail = FALSE)
      probs[i] = prob
    }
  }
  names(probs) = names(foreground)
  return(sort(probs, decreasing=FALSE))
}

hyperGluconprobs = getProb(fgGluconeogenesisDataset, bgDataset)
# count the number of 8-mer substrings of gluconeogenesis pathway with p-value<0.001
overGluc <- names(hyperGluconprobs)[which(hyperGluconprobs<0.001)]
length(overGluc)
#res = 175
hyperGlycolprobs = getProb(fgGlycolysisDataset, bgDataset)
# count the number of 8-mer substrings of glycolysis pathway with p-value<0.001
overGly <- names(hyperGlycolprobs)[which(hyperGlycolprobs<0.001)]
length(overGly)
 #res = 487

#find the common overrepresenting substring in both pathways
length(intersect(overGluc, overGly))


# get the 8mer with the smallest p-value for each pathway
mostfreqGlu <- names(hyperGluconprobs[1])
#                 _____AATCTCGC____ 

mostfreqGly <- names(hyperGlycolprobs[1])
#                 _____AAAGACGC____ 



# find sequences similar to the substring with the smallest p-value for each pathway
getAllInstances = function(candidate, foreground, threshold){
  allnames = names(foreground)
  motifstrings = c()
  for(i in 1:length(allnames)){
    if( stringdist(candidate, allnames[i], method = "hamming") < threshold){
      motifstrings = c(motifstrings, rep(allnames[i], foreground[i]))
    }
  }
  return(motifstrings)
}

stringGlycolMotifs = getAllInstances(mostfreqGly, fgGlycolysisDataset, 3)
stringGluconMotifs = getAllInstances(mostfreqGlu, fgGluconeogenesisDataset, 3)




getACGT = function(dataset, alphabet=c("A", "C", "G", "T")){
  counts = vector("numeric", length=length(alphabet))
  names(counts) = alphabet
  for(i in 1:length(names(dataset))){
    nam = strsplit(names(dataset)[i], "")[[1]]
    ntall = rep(nam, each=dataset[i])
    ntcounts = table(factor(ntall, levels=alphabet))
    counts[names(ntcounts)] = counts[names(ntcounts)] + ntcounts 
  }
  counts/sum(counts)
}


basefreqs = getACGT(bgDataset)
basefreqs


getPWM = function(stringMotifs, length=6, alphabet =c("A", "C", "G", "T"),  freqs = rep(0.25, 4)){
  pfm = matrix(0, nrow=4, ncol=length)
  row.names(pfm) = alphabet
  for(i in 1:length(stringMotifs)){
    v = unlist(strsplit(stringMotifs[i], "")) ## or strsplit(stringMotifs[i], "")[[1]]
    for(j in 1:length(v)){
      pfm[v[j], j] = pfm[v[j], j] + 1
    }
  }
  ppm = pfm/colSums(pfm)
  pwm = pwm = log2((ppm+1e-4)/freqs)
  return(list(pwm=pwm, ppm=ppm))
}

psGlucon <- getPWM(stringMotifs = stringGluconMotifs, length=motifLength, freqs = basefreqs)
psGlycol <- getPWM(stringMotifs = stringGlycolMotifs, length=motifLength, freqs = basefreqs)
psGlucon
psGlycol

ggseqlogo(psGlucon$ppm) + ggtitle("Gluconeogenesis")
ggseqlogo(psGlycol$ppm) + ggtitle("Glycolysis")

getScoreSimple  = function(vstring, pwm){
  score = 0
  v = strsplit(vstring, "")[[1]]
  scores = vector("numeric", length=length(v)-ncol(pwm)+1)
  for(i in 1:(length(v)-ncol(pwm)+1)){
    score = 0
    for(j in 1:ncol(pwm)){
      letter = v[i+j-1]
      score = score + pwm[letter, j] 
    }
    scores[i] = score
  }
  return(scores)
}

#Gluconeogenesis
resGluc = t(sapply(gluconeogenesis_data, getScoreSimple, psGlucon$pwm))
matrix(as.numeric(resGluc > 3), nrow=nrow(resGluc))


#glycolysis
resGly = t(sapply(glycolysis_data, getScoreSimple, psGlycol$pwm))

#   find the maximum score for each sequence in gluconeogenesis and glycolysis
maxScoresGluc <- apply(resGluc, 1, max)
maxScoresGly <- apply(resGly, 1, max)

# glycolysis matrix to scan gluconeogenesis
glyScanGlu = t(sapply(gluconeogenesis_data, getScoreSimple, psGlycol$pwm))
maxGluScannedByGly <- apply(glyScanGlu, 1, max)
gluScanGly = t(sapply(glycolysis_data, getScoreSimple, psGlucon$pwm))
maxGlyScannedByGlu <- apply(gluScanGly, 1, max)



logicalresGluc = matrix(as.numeric(resGluc > -4), nrow=nrow(resGluc))
heatmap.2(logicalresGluc, dendrogram='none', Rowv=F, Colv=F, trace = 'none', main = "Gluconeogenesis")






logicalresGly = matrix(as.numeric(resGly > -4), nrow=nrow(resGly))
logicalresGlu = matrix(as.numeric(resGluc > -4), nrow=nrow(resGluc))
heatmap.2(logicalresGly, dendrogram='none', Rowv=F, Colv=F, trace = 'none', main = "Glycolysis")

heatmap.2(logicalresGlu, dendrogram='none', Rowv=F, Colv=F, trace = 'none', main = "Gluconeogenesis")









