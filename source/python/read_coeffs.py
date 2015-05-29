import numpy as np

filename = "morten.txt"
with open(filename, "r") as f: # infile
	with open ("coeff_out.cpp", "w") as out: # outfile
		# reading from coefficient file
		coeffs = np.zeros([9,5]) # coefficient matrix
		i = 0
		for line in f: 
			elements = line.split()
			if len(elements) > 0: # if line not empty
				for j, element in enumerate(elements):
					coeffs[i][j] = float(element)
				i += 1
		
		#print coeffs

		# writing to file
		out.write("if (j == 0) {\n")
		for j in range(5):
			if j != 0 : out.write("}\nelse if (j == %d) {\n" %(j))
			for i in range(9):
				#print j, i
				out.write("    coeffs(%d) = %e;\n" %(i, coeffs[i][j]))


				

				