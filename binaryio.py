from struct import unpack

# A library for reading binary files written with fortran.

# Rewind NUM entries
def rewind(f,num):
  for nrec in range(num):
    f.seek(-4,1)
    foot = f.read(4)
    f1 = unpack('i',foot)[0]
    f.seek(-4,1)
    f.seek(-f1,1)
    f.seek(-4,1)
    head = f.read(4)
    h1 = unpack('i',head)[0]
    f.seek(-4,1)
    if h1 != f1 :
      print "Error during skip ", nrec, ": h1 = ", h1, " and f1 = ", f1
      return -1


# Reset pointer to first bit in file.
def restart(f):
  f.seek(0,0)


# Set pointer to last bit in file.
def end(f):
  f.seek(0,2)


# Move forward NUM entries.
def skip(f,num):
  for nrec in range(num):
    head = f.read(4)
    h1 = unpack('i',head)[0]
    f.seek(h1,1)
    foot = f.read(4)
    f1 = unpack('i',foot)[0]
    if h1 != f1 :
      print "Error during skip ", nrec, ": h1 = ", h1, " and f1 = ", f1
      return -1


# Read NUM integer variables in a single write entry
def bread_int(f,num,rd_head=True):
  if rd_head:
    head = f.read(4)
    h1 = unpack('i',head)[0]
    vals = f.read(h1)
  
    foot = f.read(4)
    f1 = unpack('i',foot)[0]

    if h1 != f1 :
      print "Error: h1 = ", h1, " and f1 = ", f1
      return -1

    nint = h1/4
    if nint*4 != h1 :
      print "Error: h1 isn't multiple of 4. Weird output: h1=",h1
      return -1

    if num == 1:
      return unpack('i'*nint,vals)[0]
    else :
      return unpack('i'*nint,vals)[0:num]
  else :
    vals = f.read(4*num)
    if num == 1:
      return unpack('i'*num,vals)[0]
    else :
      return unpack('i'*num,vals)[0:num]


# Read NUM real*4 variables in a single write entry
def bread_float(f,num,rd_head=True):
  if rd_head:
    head = f.read(4)
    h1 = unpack('i',head)[0]
    vals = f.read(h1)

    foot = f.read(4)
    f1 = unpack('i',foot)[0]

    if h1 != f1 :
      print "Error: h1 = ", h1, " and f1 = ", f1
      return -1

    nfloat = h1/4
    if nfloat*4 != h1 :
      print "Error: h1 isn't multiple of 4. Weird output. h1=",h1
      return -1

    if num == 1:
      return unpack('f'*nfloat,vals)[0]
    else :
      return unpack('f'*nfloat,vals)[0:num]
  else :
    vals = f.read(4*num)
    if num == 1:
      return unpack('f'*num,vals)[0]
    else :
      return unpack('f'*num,vals)[0:num]


# Read NUM real*8 variables in a single write entry
def bread_dbl(f,num,rd_head=True):
  if rd_head:
    head = f.read(4)
    h1 = unpack('i',head)[0]
    vals = f.read(h1)

    foot = f.read(4)
    f1 = unpack('i',foot)[0]

    if h1 != f1 :
      print "Error: h1 = ", h1, " and f1 = ", f1
      return -1

    ndbl = h1/8
    if ndbl*8 != h1 :
      print "Error: h1 isn't multiple of 8. Weird output. h1=",h1
#    return -1

    if num == 1:
      return unpack('d'*ndbl,vals)[0]
    else :
      return unpack('d'*ndbl,vals)[0:num]
  else :
    vals = f.read(8*num)
    if num == 1:
      return unpack('d'*num,vals)[0]
    else :
      return unpack('d'*num,vals)[0:num]


# Read NUM char variables in a single write entry
def bread_char(f,num):
  head = f.read(4)
  h1 = unpack('i',head)[0]
  vals = f.read(h1)

  foot = f.read(4)
  f1 = unpack('i',foot)[0]

  if h1 != f1 :
    print "Error: h1 = ", h1, " and f1 = ", f1
    return -1

  nchar = h1

  if num == 1:
    return unpack('c'*nchar,vals)[0]
  else :
    return unpack('c'*nchar,vals)[0:num]


# Read NUM variables of mixed types in a single write entry
def bread_mixed(f,types):
  head = f.read(4)
  h1 = unpack('i',head)[0]
  vals = f.read(h1)

  foot = f.read(4)
  f1 = unpack('i',foot)[0]

  if h1 != f1 :
    print "Error: h1 = ", h1, " and f1 = ", f1
    return -1

  length_req = 0
  returned_vals = []
  count = 0
  for t in types:
    count = count + 1
    if t == 'd':
      length_req = length_req + 8
    if t == 'i':
      length_req = length_req + 4
    if t == 'f':
      length_req = length_req + 4
    if t == 'c':
      length_req = length_req + 1

  if length_req > h1:
    print "Error: h1 is smaller than desired record length: h1=",h1
    return -1

  extra = (h1 - length_req)/4

  if extra*4+length_req != h1 :
    print "Warning: Leftover space isn't divisible by 4."

  extra = h1 - length_req
  return unpack('='+types+'c'*extra,vals)[0:count]

