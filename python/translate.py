def getTestFiles():
  from os import listdir
  inpath = "../R/"
  f = listdir(inpath)
  
  testFiles = []
  for fi in f:
    if fi[0:5] == "test.":
      testFiles.append(fi[5:(-2)])
  
  return(testFiles)

def topMatter(testname, outfile):
  outfile.write("context(\"" + testname + "\")\n\n")
  outfile.write("test_that(\"" + testname + " behaves as it should\", {\n")

def checkEquals(fname):
  inpath = "../R/"
  outpath = "Rout/"
  
  files = getTestFiles()
  
  for f in files:
    infile = file(inpath + "test." + f + ".R", "r")
    outfile = file(outpath + "test." + f + ".R", "w")

    topMatter(f, outfile)

    for line in infile:
      line = line.replace("msg=", "label=")
      line = line.replace("texmexPst", "texmex:::texmexPst")
      sline = line.strip()

      if "test." + f + " <-" in sline:
        pass
      elif "function(){" in sline:
        sline = sline.split("(){")[1]
        outfile.write("  " + sline)
      elif sline[0:18] == "checkEqualsNumeric":
        index = 19
        translateArgs(sline, index, outfile)
      elif sline[0:11] == "checkEquals":
        index = 12
        translateArgs(sline, index, outfile)
      elif sline[0:14] == "checkException":
        index = 15
        translateArgs(sline, index, outfile)
      elif sline[0:9] == "checkTrue":
        index=10
        translateArgs(sline, index, outfile)
      else:
        outfile.write(line)
  
    outfile.write(")\n")
  
    infile.close()
    outfile.close()

  

def translateArgs(txt, i, outfile):
    checkTrue = False
    if i == 10:
      checkTrue = True
    txt = txt[i:]
    outfile.write("  expect_that(")
    i = 0

    for char in txt:
      if char != ",":
        outfile.write(char.strip())
        i += 1
      else:
        break

    if not checkTrue:
      outfile.write(", equals(")
      for char in txt[(i+1):]:
        if char != ",":
          outfile.write(char.strip())
          i += 1
        else:
          break
      outfile.write("), ")

    else: # if checkTrue
      outfile.write(", is_true(), ")
      for char in txt[(i+2):]:
        outfile.write(char.strip())
      outfile.write("\n")


checkEquals("wh")
#getTestFiles()
