import numpy as np
from matplotlib.mlab import PCA
import numpy.linalg as LA
import math

#put the data in the correct format for PCA analysis. 
def getList(files, version):
    f = open("/home/ichproject/Expression_Profiles2/"+files,'rU')
    myList = []

    if "Perc" in version: #if we look at the expression data in terms of percentages
        totalPerc = []
        counter = 0
        for line in f:
            line = line.split()
            fline = [float(x) for x in line[1:]] #for not transposed data
            for h in range(len(fline)):
                if counter == 0:
                    if "Log" in version:
                        totalPerc.append(math.log(fline[h]))
                    else:
                        totalPerc.append(fline[h])
                else:
                    if "Log" in version: #if logging all expression data
                        totalPerc[h] = totalPerc[h] + math.log(fline[h])
                    else:
                        totalPerc[h] = totalPerc[h] + fline[h]
            counter = counter + 1
        f.seek(0,0)

    for line in f:
        hourList = []
        line = line.split()
        fline = [float(x) for x in line[1:]]  #for not transposed data
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

        if "diff" in version: #if taking the difference between time intervals
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
    fileList = ["ConjugationA1.txt","ConjugationA2.txt", "StarvationA1.txt","StarvationA2.txt","StarvationA3.txt", "GrowthA1.txt", "GrowthA2.txt","GrowthA3.txt"]
    conjList =[] #list of eigenvectors for conjugation
    starvList = []# for starvation data
    growthList = [] #for growth data
    for files in fileList: # get eigenvector for each data file
        eigenvector = []
        myList = getList(files,version) #get data for PCA analysis
        data = np.array(myList)
        data -= np.mean(data, axis=0)
        c = np.cov(data.transpose())
        eigval, eigvec = np.linalg.eig(c) # create eigenvectors
        print len(eigvec), "len of eigenvector"
        print eigvec, "eigvec"

        eig_pairs = [(np.abs(eigval[i]), eigvec[:,i]) for i in range(len(eigval))] #pair eigenvectors with there values

        eig_pairs.sort() #order the eigenvectors from most variablility to least
        eig_pairs.reverse()
        
        for i in eig_pairs: #put eigenvectors in lists to return
            eigenvector = i[1].tolist()
            print(i[1]), "eig vec"
            if files[:-6] == "Conjugation":
                conjList.append(eigenvector)
            elif files[:-6] == "Starvation":
                starvList.append(eigenvector)
            elif files[:-6] == "Growth":
                growthList.append(eigenvector)

    return conjList, starvList, growthList

def main():
    #versions are the different types of PCA analysis done
    versions = ["raw","rawPerc","rawPercLog","diff","diffPerc","diffPercLog", "rawLog","diffLog"]
    for version in versions:
        conjList, starvList, growthList = getPCA(version) #do PCA of each file for each version

        for x in range(len(conjList)/2): #put conjugation data into a table
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

        for z in range(len(starvList)/3): #put starvation data into a table
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

        for w in range(len(growthList)/3): #put growth data into a table
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


if __name__ == "__main__":
    main()
