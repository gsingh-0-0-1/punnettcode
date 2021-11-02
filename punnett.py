#Gurmehar Singh, 1 Nov 2021

import numpy as np

#Let's say we have two "parents" that we are crossing
#in order to find the genotypes of their children.
#Let's also say that we don't know how many genes we
#are dealing with. So, the relevant genotypes for our
#parents might look like Aa and AA, or they might look
#like SsYy and SSYy - or even something as complex as
#GgJJBbMm and GGJjBBmm.

#The goal of this script will be to take two parent
#genotypes with n genes each, and produce a comprehensive
#output of the possible genotypes of their children.

#We can view the pairing of two alleles as a 2-dimensional
#vector - for example, "Aa" can be represented as [A a].
#Thus, crossing two genotypes works as the Kronecker product - 
#or more specifically, the outer product - of those two vectors.

#In other words - crossing [A a] and [A a], treating A and a as
#mathematical quantities yields a matrix that looks like this:
#[AA Aa]
#[Aa aa]
#which should not be surprising.

#Now imagine that we have a parent with a genes like SsYy, who
#will be crossed with another parent with the same set of genes.
#To create a Punnett square to represent the offspring of such a
#pairing, we first need to obtain the possible combinations of
#genes that will be placed along the row and the column of the 
#Punnett square. This can be done by splitting up the genes of the
#two parents gene-wise into vectors: [S s] and [Y y], for both the
#parents in this case. We can then take the outer product of these
#two vectors for each parent, and then flatten the resulting matrix
#to get the following output:
#Parent 1: [SY Sy sY sy]
#Parent 2: [SY Sy sY sy]
#which is to be expected.
#Conveniently, given that we have defined a generalized pseudo-outer-product
#function, we can simply plug those two vectors back into the function to
#obtain the full Punnett Square. Try it below! This script has these
#two example gene sets assigned to two hypothetical parents - and
#you can change these as you like, with any number of two-allele genes.

#definition for the pseudo-outer-product that I will be using here.
def pseudo_outer_prod(l):
	output_dimension_num = len(l)
	shape = []
	indshape = []
	for i in range(output_dimension_num):
		shape.append(len(l[i]))
		indshape.append(1)
	shape = tuple(shape)
	indshape = tuple(indshape)

	fillstr = ''

	out_array = np.full(shape, fillstr, dtype=object)

	for i in range(output_dimension_num):
		for j in range(len(l[i])):
			inds = np.full(indshape, j)
			temparr = np.full(shape, fillstr, dtype=object)
			np.put_along_axis(temparr, inds, l[i][j], i)
			out_array = out_array + temparr


	return out_array
#Change these as you like
p1 = 'SsYy'
p2 = 'SsYy'

#split into lists
l1 = [p1[ind] + p1[ind + 1] for ind in range(0, len(p1), 2)]
l2 = [p2[ind] + p2[ind + 1] for ind in range(0, len(p2), 2)]

#make sure that the number of genes in each of the lists is the same
assert len(l1) == len(l2)

#some other checks
assert len(l1) >= 1
n = len(l1)

#obtain the possible genotype pairings for the Punnett Square we want to set up
#we can skip this for the case of n = 1 - since we'll just cross the contents of
#l1 and l2.
if n == 1:
	punnett1 = [l1[0][0], l1[0][1]]
	punnett2 = [l2[0][0], l2[0][1]]
if n > 1:
	#create a list of lists - or a list of vectors - which separate the genotypes
	#for each gene, for each parent
	cross1 = [[l1[ind][0], l1[ind][1]] for ind in range(len(l1))]
	cross2 = [[l2[ind][0], l2[ind][1]] for ind in range(len(l2))]

	#get the genotype combinations for each parent
	punnett1 = pseudo_outer_prod(cross1).flatten()
	punnett2 = pseudo_outer_prod(cross2).flatten()

#get the final output
#we can conveniently do this by plugging in our
#Punnett row and column into the outer product function:
output = pseudo_outer_prod([punnett1, punnett2])
shape = output.shape

#sort the strings alphabetically...
output = output.flatten()
for ind in range(len(output)):
	output[ind] = "".join(sorted(output[ind], key = str.casefold))
	templ = [output[ind][i] + output[ind][i + 1] for i in range(0, len(output[ind]), 2)]
	templ = ["".join(sorted(el)) for el in templ]
	sorted_output = "".join(templ)
	output[ind] = sorted_output

output = np.reshape(output, shape)

print(output)


