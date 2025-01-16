import math,os,sys,numpy,random

ang_to_bohr=1.889725989
auTOcm=219474.63137
atTOaumass=1822.8886

def get_natom(filename):
 natom=0
 with open(filename, 'r+') as ffile:
  for line in ffile:
   if 'NAtom' in line:
    natom=int((line.split())[1])
    break
 nfreq=3*natom-6
 return (natom,nfreq)

def readfile(filename,natom,nfreq):
 with open(filename, 'r+') as ffile:
  atmass=[]
  wavenum=[]
  eigenvec=numpy.zeros((nfreq,natom*3))
  for line in ffile:
   if '                         Standard orientation:' in line: # Getting last geometry
    geom=[]
    for i in range(4): line=ffile.next()
    for i in range(0,natom):
     line=ffile.next()
     for j in range(3):
      geom.append(float((line.split())[j+3]))
   if 'Frequencies --' in line:
    for i in range(2,len(line.split())):
     wavenum.append(float((line.split())[i]))
   if '  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z' in line:
    for i in range(natom):
     line=ffile.next()
     for j in range(3):
      eigenvec[len(wavenum)-3][i*3+j]=float((line.split())[j+2])
      eigenvec[len(wavenum)-2][i*3+j]=float((line.split())[j+5])
      eigenvec[len(wavenum)-1][i*3+j]=float((line.split())[j+8])
   if 'Temperature' in line and 'Pressure' in line:
    for i in range(natom):
     line=ffile.next()
     atmass.append(float((line.split())[8]))
 return (geom,atmass,wavenum,eigenvec)

# MAIN CODE STARTS HERE

if (len(sys.argv)<4):
  print ('Please run with the following arguments: 1) name of the freq output; 2,3) two atoms defining the dissociation vector')
  sys.exit(0)
filename = str(sys.argv[1])
atom1 = int(sys.argv[2])-1
atom2 = int(sys.argv[3])-1

debug=open("debug", 'w+')

if os.path.exists(filename):
 (natom,nfreq)=get_natom(filename)
# print (>> debug, get_natom(filename))
 (geom,atmass,wavenum,eigenvec)=readfile(filename,natom,nfreq)
 # print (>> debug, "Geometry:", geom )
# print (>> debug, "Atomic mass:", atmass )
# print (>> debug, "Wavenumbers:", wavenum )
# print (>> debug, "Eigenvectors:", eigenvec)
 vector1=[]
 product=[]
 for i in range(3):
  vector1.append(geom[atom2*3+i]-geom[atom1*3+i])  # change this if you want to change the dissociation coordinate
 for ifreq in range(len(wavenum)):
  vector2=[]
  for i in range(3):
   vector2.append(eigenvec[ifreq][atom2*3+i]-eigenvec[ifreq][atom1*3+i])
  product.append(abs(numpy.dot(vector1,vector2)))
 for ifreq in range(len(wavenum)):
  if(product[ifreq]>max(product)/10.):
   print (ifreq+1,wavenum[ifreq],product[ifreq]/max(product))
 print
# vector=numpy.zeros(natom*3)
# for i in range(3):
#  vector[atom1*3+i]=-geom[atom1*3+i]*atmass[atom1]
#  vector[atom2*3+i]=geom[atom2*3+i]*atmass[atom2]
#  vector[atom2*3+i]=geom[atom2*3+i]-geom[atom1*3+i]
# for i in range(len(wavenum)):
#  for j in range(3):
#   tot+=(eigenvec[i][3*atom2+j]-eigenvec[i][3*atom1+j])**2
#  print i+1,wavenum[i],numpy.dot(vector,eigenvec[i])
# print eigenvec[i]
else:
 print (filename, "does not exist")

debug.close()
