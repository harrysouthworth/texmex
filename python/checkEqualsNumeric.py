def checkEqualsNumeric(fname):
  inpath = "../R/"
  outpath = "Rout/"

  infile = file(inpath + fname, "r")
  outfile = file(outpath + fname, "w")
  
  for line in infile:
    sline = line.strip()
    if sline[0:18] == "checkEqualsNumeric":
      outfile.write(" expect_that(")
      sline = sline[19:]
      i = 0
      for char in sline:
        if char != ",":
          outfile.write(char)
          i += 1
        else:
          break
      outfile.write(",  equals(")
      for char in sline[(i+2):]:
        if char != ",":
          outfile.write(char)
          i += 1
        else:
          break
      outfile.write("), label")
      for char in sline[(i+1):]:
        if char != "\"":
          i += 1
          pass
        else:
          break
      for char in sline[i:]:
        outfile.write(char)
      outfile.write("\n")
    else:
      outfile.write(line)
      
  
  
  infile.close()
  outfile.close()

checkEqualsNumeric("test.bootmex.R")
