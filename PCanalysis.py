import numpy as np
from matplotlib.mlab import PCA
import numpy.linalg as LA
import math

def getList(files, version):
    f = open("/home/ichproject/Expression_Profiles2/"+files,'rU')
    myList = []

#    totalPerc = []
#    counter = 0
    if "Perc" in version:
        totalPerc = []
        counter = 0
        for line in f:
            line = line.split()
            #fline = [float(x) for x in line[1:]] #for not transposed data
            fline = [float(x) for x in line] # for transposed data
            for h in range(len(fline)):
                if counter == 0:
                    if "Log" in version:
                        totalPerc.append(math.log(fline[h]))
                    else:
                        totalPerc.append(fline[h])
                else:
                    if "Log" in version:
                        totalPerc[h] = totalPerc[h] + math.log(fline[h])
                    else:
                        totalPerc[h] = totalPerc[h] + fline[h]
            counter = counter + 1
        #print totalPerc, "percentages"
        f.seek(0,0)

    for line in f:
        hourList = []
        line = line.split()
#        fline = [float(x) for x in line[1:]]  #for not transposed data
        fline = [float(x) for x in line] #this is when data is transposed
        if "raw" in version:
            for hour in range(len(fline)):
                if version == "raw":
                    hourList.append(fline[hour])
                if version == "rawPerc":
                    hourList.append(fline[hour]/totalPerc[hour]*100)
                if version == "rawPercLog":
                    hourList.append(math.log(fline[hour])/totalPerc[hour]*100)
                if version == "rawLog":
                    hourList.append(math.log(fline[hour]))

        if "diff" in version:
            for hour in range(len(fline)-1):
                if version == "diff":
                    hourList.append(fline[hour+1] - fline[hour])
                if version == "diffPerc":
                    hourList.append((fline[hour+1]/totalPerc[hour+1]*100)-(fline[hour]/totalPerc[hour]*100))
                if version == "diffPercLog":
                    hourList.append((math.log(fline[hour+1])/totalPerc[hour+1]*100)-(math.log(fline[hour])/totalPerc[hour]*100))
                if version == "diffLog":
                    hourList.append(math.log(fline[hour+1]) - math.log(fline[hour]))
        myList.append(hourList)


    return myList
    
def getPCA(version):    
    #fileList = ["ConjugationA1.txt","ConjugationA2.txt", "StarvationA1.txt","StarvationA2.txt","StarvationA3.txt", "GrowthA1.txt", "GrowthA2.txt","Gro
wthA3.txt"]
    #fileList = ["GrowthA1.txt","GrowthA2.txt","GrowthA3.txt"]
    fileList = ["CombinedStarvationT.txt","CombinedConjugationT.txt"]
    conjList =[]
    starvList = []
    growthList = []
    comboStarvList = []
    comboConjList = []
    for files in fileList:
        print files
        eigenvector = []
        myList = getList(files,version)
        data = np.array(myList)
        #autoResults = PCA(data)
        #print autoResults.Wt, "results"
        mu = data.mean(axis=0) #other option
        data = data - mu #other option
# data = (data - mu)/data.std(axis=0)  # Uncomment this reproduces mlab.PCA results
        eigvec, eigval, V = np.linalg.svd(data.T, full_matrices=False) #other option
#        print eigvec, "eigvec"
#        data = data.transpose()  #trying THIS OUT
        #print data, "before"
        #data -= np.mean(data, axis=0)  #PUT BACK
        #print data
        #c = np.cov(data.transpose()) #put BACK
        #print c
        #print "covariance matrix"
#        c = np.cov(data,None,0)
#        print data.transpose(), "transposed"
#        print data, "not transposed"

#        print c, "covariance"
#        np.savetxt("covariance"+files,c,'%.14f')
        #eigval, eigvec = np.linalg.eig(c) #PUT BACK!!
        print len(eigvec), "len of eigenvector"
        print eigvec, "eigvec"

        eig_pairs = [(np.abs(eigval[i]), eigvec[:,i]) for i in range(len(eigval))] #PUT BACK

#        print eig_pairs
        eig_pairs.sort() #PUT BACK
        eig_pairs.reverse() #PUT BACK
        
#        print files
        #for i in range(len(eigvec)):
            #print eigvec[i], "eigvec", i
#            print eigval[i], "eigval", i
            #eigenvector = eigvec[i].tolist()
        for i in eig_pairs: #PUT BACK
            eigenvector = i[1].tolist() #PUT BACK
            print(i[1]), "eig vec" #PUT BACK
            #print eigenvector, "list form"
            if files[:-6] == "Conjugation":
                conjList.append(eigenvector)
            elif files[:-6] == "Starvation":
                starvList.append(eigenvector)
            elif files[:-6] == "Growth":
                growthList.append(eigenvector)
            elif files[:-5] == "CombinedStarvation": #-5 FOR transpose, -4 otherwise
                comboStarvList.append(eigenvector)
            elif files[:-5] == "CombinedConjugation": #-5 FOR transpose, -4 otherwise
                comboConjList.append(eigenvector)

    return conjList, starvList, growthList, comboStarvList, comboConjList

def main():
#    fileList = ["ConjugationA1.txt","ConjugationA2.txt", "StarvationA1.txt","StarvationA2.txt","StarvationA3.txt"]

#    for files in fileList:
#        print files
#        myList = getList(files,"raw")
    versions = ["raw","rawPerc","rawPercLog","diff","diffPerc","diffPercLog", "rawLog","diffLog"]
    for version in versions:
        conjList, starvList, growthList, comboStarvList, comboConjList = getPCA(version)
#        print starvList, "starvation"

        for x in range(len(conjList)/2):
            f = open("conjugationPCA"+version+str(x+1)+".txt",'w')
            f.write("rep1,rep2")
            for j in range(len(conjList[x])):
                if "diff" in version:
                    f.write("\n"+str(conjList[x][j])+","+str(conjList[x+9][j]))
                    print conjList[x][j], "conj1", conjList[x+9][j], "conj2"
                else:
                    f.write("\n"+str(conjList[x][j])+","+str(conjList[x+10][j]))
                    print conjList[x][j], "conj1", conjList[x+10][j], "conj2"
            print "NEXT PCA"

        for z in range(len(starvList)/3):
            re = open("starvationPCA"+version+str(z+1)+".txt",'w')
            re.write("rep1,rep2,rep3")
            for y in range(len(starvList[z])):
                if "diff" in version:
                    re.write("\n"+str(starvList[z][y])+","+str(starvList[z+6][y])+","+str(starvList[z+12][y]))
                    print starvList[z][y], "starv1", starvList[z+6][y], "starv2", starvList[z+12][y], "starv3"
                else:
                    re.write("\n"+str(starvList[z][y])+","+str(starvList[z+7][y])+","+str(starvList[z+14][y]))
                    print starvList[z][y], "starv1", starvList[z+7][y], "starv2", starvList[z+14][y], "starv3"
            print "NEXT PCA STARVATION"

        for w in range(len(growthList)/3):
            re = open("growthPCA"+version+str(w+1)+".txt",'w')
            re.write("rep1,rep2,rep3")
            for v in range(len(growthList[w])):
                if "diff" in version:
                    re.write("\n"+str(growthList[w][v]) +"," +str(growthList[w+2][v]) + "," +str(growthList[w+4][v]))
                    print growthList[w][v], "growthv1",growthList[w+2][v],"growthv2",growthList[w+4][v],"growthv3"
                else:
                    re.write("\n"+str(growthList[w][v]) +"," +str(growthList[w+3][v]) + "," +str(growthList[w+6][v]))
                    print growthList[w][v], "growthv1",growthList[w+3][v],"growthv2",growthList[w+6][v],"growthv3"
            print "NEXT PCA GROWTH"

#        for a in range(len(comboStarvList)):
        for a in range(0,2): #I only want 2 PCs
            re = open("comboStarvPCA"+version+str(a+1)+".txt",'w')
            re.write("rep1")
            for b in range(len(comboStarvList[a])):
                re.write("\n"+str(comboStarvList[a][b]))
#                print comboStarvList[a][b], "comboStarv"
#            print "NEXT PCA comboStarv"
        
#        for c in range(len(comboConjList)):
        for c in range(0,2): #I only want 2 PCs
            re = open("comboConjPCA" + version+str(c+1)+".txt",'w')
            re.write("rep1")
            for d in range(len(comboConjList[c])):
                re.write("\n"+str(comboConjList[c][d]))
#                print comboConjList[c][d], "comboConj"
#            print "NEXT PCA comboConj"
if __name__ == "__main__":
    main()
