write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  fasta_lines <- character()
  for (i in seq_along(sequences)) {
    chr <- bed_file[i, "seqnames"]
    start <- bed_file[i, "start"]
    end <- bed_file[i, "end"]
    seq_name <- paste0(chr, ":", start, "-", end)
    fasta_lines <- c(fasta_lines, paste0(">", seq_name), as.character(sequences[[i]]))
  }
  writeLines(fasta_lines, con = filepath)
}
