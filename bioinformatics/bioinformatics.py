from sequencer import *
from termcolor import colored

""""Input data"""
repeats_experiment = "TTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAAACGTTGACACGTACGTACGTTGACACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGACACGTACGTACGTACGTACGTACGTACGTACGTACTGACGTACGTACGTACGTACGTACGTACGTACGTACGTTGACTTCTCGAGTTCCGTCCGGTTAAGCGTGACAGTCCCAGTGAACCCACAAACCGTGATGGCTGTCCTTGGAGTCATACGCAA"

"""Count the repetitions"""
count_of_acgt={"ACGTACGTACGTACGT":0 }
for i in range(len(repeats_experiment)/16):
	if repeats_experiment[i:i+16] == "ACGTACGTACGTACGT":
		count_of_acgt["ACGTACGTACGTACGT"]+=1

print(count_of_acgt)
with open("count_of_genes.txt","w") as f:
    for item in count_of_acgt:
        f.write("%s" % item+"->"+str(count_of_acgt[item])+"\n")
generated_reads=generate_reads(repeats_experiment,16)

with open('randomized_reads.txt', 'w') as f:
    for item in generated_reads:
        f.write("%s\n" % item)
reconstructed_genome = genome_reconstruction(generated_reads)

distorted_posels = {}
for i in range(0,len(reconstructed_genome)):
    if reconstructed_genome[i]!=repeats_experiment[i]:
        distorted_posels[i]=repeats_experiment[i]+"->"+reconstructed_genome[i]

with open("distorted_posels.txt",'w') as f:
    for i in range(0,len(reconstructed_genome)):
        if reconstructed_genome[i]!=repeats_experiment[i]:
            distorted_posels[i]=repeats_experiment[i]+"->"+reconstructed_genome[i]
            f.write("%s" % str(i)+": "+distorted_posels[i]+"\n")

k=37
generated_reads=generate_reads(repeats_experiment,k)
reconstructed_genome = genome_reconstruction(generated_reads)
print "------------------"+ str(k)+"-mers--------------------"
for i in range(len(reconstructed_genome)):
	if reconstructed_genome[i]==repeats_experiment[i]:
		print colored(reconstructed_genome[i],"green"),
	else:
		print colored(reconstructed_genome[i],"red") +"("+colored(repeats_experiment[i],"white")+")",
