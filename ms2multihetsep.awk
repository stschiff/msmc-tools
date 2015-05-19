BEGIN {
  OFS = "\t"
}

/^\[[0-9]+\]/ {
  split($0, f, "]")
  printf("%s %s\n", substr(f[1], 2), f[2]) > treeFile
}

/^positions/ {
  for(i = 2; i <= NF; i += 1) {
    positions[i - 1] = $i
    alleles[i - 1] = ""
  }
  while(getline) {
    sites = $0
    for(i = 1; i <= length(positions); i += 1) {
      alleles[i] = alleles[i] substr(sites, i, 1)
    }
  }
  lastPos = 0
  for(i = 1; i <= length(positions); i += 1) {
    realPos = int(positions[i] * L)
    if(realPos > lastPos) {
      print chr, realPos, realPos - lastPos, alleles[i] > sitesFile
    }
    lastPos = realPos
  }
}