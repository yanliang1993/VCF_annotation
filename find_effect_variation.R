codon_to_effect <- function(txid, codon_start, codon_end, alt, alt_idx){  ## helper function: define the effect using altered codon
  effect_temp <- c()
  s <- extractTranscriptSeqs(genome, cdsBy(txdb)[txid])
  codon <- DNAString(substring(s, codon_start, codon_end))
  codon_alt <- codon
  for(i in 1:length(alt)){
    codon_alt[alt_idx + i - 1] <- alt[i]
  }
  ## find reference and altered amino acid (DO NOT consider alternative initiation codons)
  ref_aa <- translate(codon, genetic.code = GENETIC_CODE)
  alt_aa <- translate(codon_alt, genetic.code = GENETIC_CODE)
  ## check if encoding amino acid is changed
  if (ref_aa == alt_aa){
    effect_temp <- c(effect_temp, "Silent mutation") 
  } else if (as.character(ref_aa) == 'ATG') {
    effect_temp <- c(effect_temp, "Start codon lost")
  } else if (as.character(ref_aa) == '*') {
    effect_temp <- c(effect_temp, "Stop codon lost")
  } else if (as.character(alt_aa) == '*') {
    effect_temp <- c(effect_temp, "Nonsense mutation")
  } else { 
    effect_temp <- c(effect_temp, "Missense mutation")
  }
  return(effect_temp)
}


find_effect_snp <- function(chr, pos, ref, alt){  ## define the effect of query snp
  if (nchar(alt) != 1){
    ref_split <- strsplit(ref, '')[[1]]
    alt_split <- strsplit(alt, '')[[1]]
    ref <- ref_split[ref_split != alt_split]
    alt <- alt_split[ref_split != alt_split]
  }
  gr <- GRanges(seqnames = chr, strand = "+", ranges = IRanges(start = pos, end = pos))
  ## check if query position is in cds or not
  if_cds <- suppressWarnings(subsetByOverlaps(cds(txdb), gr))
  if (length(if_cds) == 0) return("Intergenic/Intragenic variant")
  ## find codon
  map <- mapToTranscripts(gr, cdsBy(txdb), ignore.strand=F)
  effect_temp <- c()
  for (i in 1:length(map)){  ##map to multiple transcripts
    codon_start <- (start(ranges(map))[i] - 1) %/% 3 * 3 + 1
    codon_end <- codon_start + 2
    idx <- start(ranges(map))[i] - codon_start + 1   #record the position within codon
    effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start, codon_end, alt, idx))
  }
  ## find the most deleterious effect
  return(effect_temp[which.max(effect_severity[effect_temp])])
}


find_effect_insdel <- function(chr, pos, ref, alt){ ## define the effect of query insert or deletion
  gr <- GRanges(seqnames = chr, strand = "+", ranges = IRanges(start = pos, end = pos+max(nchar(ref), nchar(alt))))
  ## check if query position is in cds or not
  if_cds <- suppressWarnings(subsetByOverlaps(cds(txdb), gr))
  if (length(if_cds) == 0) return("Intergenic/Intragenic variant")
  ## check the change of length
  length_diff <- nchar(alt) - nchar(ref)
  if (length_diff %% 3 != 0) return("Frameshift")
  else if (nchar(alt) > nchar(ref)) return("Inframe insertion")
  else return("Inframe deletion")
}


find_effect_mnp <- function(chr, pos, alt){ ## define the effect of query mnp (n=2)
  alt <- strsplit(alt, '')[[1]]
  gr <- GRanges(seqnames = chr, strand = "+", ranges = IRanges(start = pos, end = pos+1))
  ## check if query position is in cds or not
  if_cds <- suppressWarnings(subsetByOverlaps(cds(txdb), gr))
  if (length(if_cds) == 0) return("Intergenic/Intragenic variant")
  ## find codon
  map <- mapToTranscripts(gr, cdsBy(txdb), ignore.strand=F)
  effect_temp <- c()
  for (i in 1:length(map)){
    codon_start <- (start(ranges(map))[i] - 1) %/% 3 * 3 + 1
    codon_end <- codon_start + 2
    idx <- start(ranges(map))[i] - codon_start + 1   ## record the position within codon of first base
    if (idx == 3){  ## two bases are NOT in the same codon
      ## for the first base:
      effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start, codon_end, alt[1], 3))
      ## for the second base:
      effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start+3, codon_end+3, alt[2], 1))
    } else{
      effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start, codon_end, alt, idx))
    }
  }
  ## find the most deleterious effect
  return(effect_temp[which.max(effect_severity[effect_temp])])
}


find_effect_complex <- function(chr, pos, ref, alt){
  gr <- GRanges(seqnames = chr, strand = "+", ranges = IRanges(start = pos, end = pos+max(nchar(ref), nchar(alt))))
  ## check if query position is in cds or not
  if_cds <- suppressWarnings(subsetByOverlaps(cds(txdb), gr))
  if (length(if_cds) == 0) return("Intergenic/Intragenic variant")
  ## check if the length was changed
  length_diff <- nchar(alt) - nchar(ref)
  if (length_diff != 0){
    if (length_diff %% 3 != 0) return("Frameshift")
    else if (nchar(alt) > nchar(ref)) return("Inframe insertion")
    else return("Inframe deletion")
  }
  ## find codon (when REF and ALT have the same length)
  map <- mapToTranscripts(gr, cdsBy(txdb), ignore.strand=F)
  effect_temp <- c()
  alt <- strsplit(alt, '')[[1]]
  for (i in 1:length(map)){
    codon_start <- (start(ranges(map))[i] - 1) %/% 3 * 3 + 1
    codon_end <- codon_start + 2
    idx <- start(ranges(map))[i] - codon_start + 1   ## record the position within codon of first base
    effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start, codon_end, alt[1:(4-idx)], idx))
    if (length(alt) > 3 | idx != 1){
      n <- 1 + (length(alt) - 5 + idx)%/%3  ## number of affected codons left
      for (j in 1:n){
        effect_temp <- c(effect_temp, codon_to_effect(map@seqnames[i], codon_start+3*j, codon_end+3*j, 
                                                      alt[(5-idx+3*(n-1)):min(length(alt), (5-idx+3*(n-1)+2))], 1))
      }
    }
  }
  ## find the most deleterious effect
  return(effect_temp[which.max(effect_severity[effect_temp])])
}