from math import sqrt

file1 = open("freq_RMS_5.000000_250.000000_600.000000_3.000000<-.dat","r")
file2 = open("freq_RMS_5.000000_250.000000_600.000000_4.000000<-.dat","r")
file3 = open("freq_RMS_5.000000_250.000000_600.000000_5.000000<-.dat","r")

data1 = []
data2 = []
data3 = []
data4 = []
data5 = []

for line in file1:
    temp = line.split()
    data1.append([float(temp[0]),float(temp[1])])

for line in file2:
    temp = line.split()
    data2.append([float(temp[0]),float(temp[1])])

for line in file3:
    temp = line.split()
    data3.append([float(temp[0]),float(temp[1])])

##for line in file4:
##    temp = line.split()
##    data4.append([float(temp[0]),float(temp[1])])
##
##for line in file5:
##    temp = line.split()
##    data5.append([float(temp[0]),float(temp[1])])

def diference(d1,d2,shift):
    total = 0
    for i in range(min(len(d1),len(d2))-20):
        if((d1[i][0]>300)&(d2[i+shift][0]>300)&(d1[i][0]<590)&(d2[i+shift][0]<590)):
            total += (d1[i][1]-d2[i+shift][1])**2
    return(total)

s12 = -10
s13 = -2
s14 = 0
s15 = 0

print(diference(data1,data2,s12))
print(diference(data1,data3,s13))
##print(diference(data1,data4,s14))
##print(diference(data1,data5,s15))

out_file = open("final.dat","w")
out_average = open("average_5.0<-.dat","w")

aver_data = []
error_data = []

for i in range(min(len(data1),len(data2),len(data3))-20):
    if((data1[i][0]>300)&(data2[i + s12][0]>300)&(data3[i + s13][0]>300)&(data1[i][0]<590)&(data2[i + s12][0]<590)&(data3[i + s13][0]<590)):
        out_file.write(str(data1[i][0])+"\t"+str(data1[i][1])+"\t"+str(data2[i + s12][1])+"\t"+str(data3[i + s13][1])+"\n")
        f_av = (data1[i][0] + data2[i + s12][0] + data3[i + s13][0])/3
        v_av = (data1[i][1] + data2[i + s12][1] + data3[i + s13][1])/3
        sigma_f = sqrt(((data1[i][0] - f_av)**2 + (data2[i + s12][0] - f_av)**2 + (data3[i + s13][0] - f_av)**2)/2)
        sigma_v = sqrt(((data1[i][1] - v_av)**2 + (data2[i + s12][1] - v_av)**2 + (data3[i + s13][1] - v_av)**2)/2)
        out_average.write(str(f_av)+"\t"+str(v_av)+"\t"+str(3.0*sigma_f/2)+"\t"+str(3.0*sigma_v/2)+"\n")
        
out_file.close()
out_average.close()

            
