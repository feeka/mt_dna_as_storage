from workflow_manager import WorkflowManager
from copy import *
import numpy as np
import matplotlib.pyplot as plt
from bioinformatics.work_without_ecc import *
from helper import *
e = [0,0.05,0.1,0.15,0.2,0.25,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.8,0.9,1.0]
x_values= []
y_values_for_cw= []
y_values_for_dist_r= []
#######RAW DATA FOR TESTING########
plt.figure(figsize=(12, 25))
for t in range(3000):
    print(t)

    if t==2999:
        break
    try:
        proto_messages = generate_messages(7,5)
        k = len(proto_messages[0])
        polynom = [1,1,1,2,3]
        degree = 3
        s_part=[1,1,1,1,1,0]
        n=15

        flag = True
        states = []
        wm = None
        wm = WorkflowManager(proto_messages,polynom,s_part,degree,n)
        wm.generate_messages()
        wm.encode()
        wm.map_and_synthesize()

        y_for_cw = []
        y_for_dist_reads = []
        x_data = []

        for cnt1 in e:
            x_data.append(cnt1)
            wm.create_reads_sequence(cnt1)#reads are created with 2 distorted positions
            wm.detect_codewords()
            Temp = wm.assemble()
            lal = wm.compare_data()

        #cnt1=cnt1+0.05
            y_for_cw.append(lal[0]/len(wm.original))
        #print(y_for_cw)
    #x_values.append(x_data)

        count=0
        y_values_for_cw.append(y_for_cw)
#print(from_dist_reads)
        for i in x_data:
            wwe=WorkWithoutECC(proto_messages,i)
            wwe.concatenate_to_one_string()
            wwe.distort()
    #wwe.distort_reads_weak()
            wwe.make_reads()
    #from_dist_reads = greedy_scs(list(reads_from_original),len(reads_from_original[0])-2)
    #from_dist_reads = greedy_scs(list(reads_from_original),len(reads_from_original)-2)
            y_for_dist_reads.append((wwe.compare_data())/len(wwe.one_str))

        #print(y_for_dist_reads)

        y_values_for_dist_r.append(y_for_dist_reads)
        #plt.plot(x_data, y_for_dist_reads,'r',label='y_for_dist_reads')
        #plt.plot(x_data, y_for_cw,'g',label='y_for_cw')
        #plt.legend()
        #plt.ylabel("Errors")
        #plt.xlabel("P_err")
        #plt.ylim(top=1)
        #plt.savefig("results_two/reads_"+str(t)+".png")
        #print("graph drawn")
    except RecursionError:
        continue

#plt.legend()

# And plot it

"""
    plt.figure(figsize=(12, 9))

    plt.plot(x_data, y_for_dist_reads,'r',label='y_for_dist_reads')
    plt.plot(x_data, y_for_cw,'g',label='y_for_cw')
    plt.legend()if self.probability_a>0.5:
    plt.ylabel("Errors")
    plt.xlabel("P_err")
    plt.ylim(top=1)
    plt.savefig("results_one/reads_"+str(t)+".png")
"""
average_y_for_cw=[sum(x)/3000 for x in zip(*y_values_for_cw)]
average_y_for_dist_reads=[sum(x)/3000 for x in zip(*y_values_for_dist_r)]
print("cw",average_y_for_cw)
print("dist_reads",average_y_for_dist_reads)
plt.figure(figsize=(12, 9))

plt.plot(e, average_y_for_dist_reads,'r',label='y_for_dist_reads')
plt.plot(e, average_y_for_cw,'g',label='y_for_cw')
plt.legend()
plt.ylabel("Errors")
plt.xlabel("P_err")
plt.ylim(top=1)
plt.savefig("results_one/reads_big_pic3000"+".png")
# And plot it
