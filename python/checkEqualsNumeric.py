def checkEqualsNumeric(fname):
  inpath = "../R/"
  outpath = "Rout/"

  infile = file(inpath + fname, "r")
  outfile = file(outpath + fname, "w")
  
  for line in infile:
    sline = line.strip()
    
    if sline[0:18] == "checkEqualsNumeric":
      index = 19
      translateArgs(sline, index, outfile)
    elif sline[0:11] == "checkEquals":
      index = 12
      translateArgs(sline, index, outfile)
    else:
      outfile.write(line)

  infile.close()
  outfile.close()

def translateArgs(txt, i, outfile):
    txt = txt[i:]
    outfile.write(" expect_that(")
    i = 0
    for char in txt:
      if char != ",":
        outfile.write(char.strip())
        i += 1
      else:
        break
    outfile.write(", equals(")
    for char in txt[(i+1):]:
      if char != ",":
        outfile.write(char.strip())
        i += 1
      else:
        break
    outfile.write("), label")
    for char in txt[(i+1):]:
      if char != "\"":
        i += 1
        pass
      else:
        break
    for char in txt[i:]:
      outfile.write(char.strip())
    outfile.write("\n")


checkEqualsNumeric("test.bootmex.R")
