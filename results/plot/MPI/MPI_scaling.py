import matplotlib.pyplot as plt
import pandas as pd
import csv

#--------------------------- Use csv.reader ------------------------------#
# file = open("compute_time.csv")
# # print(f"file type: {type(file)}")
# csvreader = csv.reader(file)

# x = []
# y = {}
# for i,row in enumerate(csvreader):
#     # first row is the x axis(# cores)
#     if (i==0):
#         x = [int(core) for core in row[1:]]
#     # 5 - 7th rows are the speedup for 128x128, 256x256, 1024x1024 respectively
#     # if (i==5 or i==6 or i==7):
#     if (row[0]=="128x128 speedup" or row[0]=="256x256 speedup" or row[0]=="1024x1024 speedup"):
#         y[row[0]] = [float(speedup) for speedup in row[1:]]
# file.close()

#--------------------------- Alternatively use pandas ------------------------------#
data = pd.read_csv("compute_time.csv")
x = [eval(x) for x in data.columns.tolist()[1:]]
# x = list(map(int,data.columns.tolist()[1:]))
y = {}

for row in data.values:
    k = row[0]
    v = row[1:].tolist()
    if (k=="128x128 speedup" or k=="256x256 speedup" or k=="1024x1024 speedup"):
        y[k] = v

# for row in range(data.shape[0]):
#     if (list(data.loc[row])[0]=="128x128 speedup" or list(data.loc[row])[0]=="256x256 speedup" or list(data.loc[row])[0]=="1024x1024 speedup"):
#         y[list(data.loc[row])[0]] = list(data.loc[row])

#--------------------------- Start plotting ------------------------------#
color=''
for k,v in y.items():
    if k=="128x128 speedup":
        color='c'
    elif k=="256x256 speedup":
        color='g'
    else:
        color='m'
    plt.plot(x,v,'o-',markersize=4,linewidth=1,color=color,label=k)
plt.legend(loc="upper left")
# plt.title('MPI speedup')
plt.xlabel('#cores',fontsize=12)
plt.ylabel('speed-up',fontsize=12)
plt.tight_layout()
# plt.show() #this cannot be used together with plt.savefig() to produce the correct svg files
plt.savefig("svg/MPI_speedup.svg")
plt.close()