import sys
import difflib

def main(output1, output2):
	with open(output1, 'r') as rfile1, open(output2, 'r') as rfile2:
		lines1 = rfile1.readlines()
		lines2 = rfile2.readlines()
		nlines1 = len(lines1)
		nlines2 = len(lines2)
		minlines = min(nlines1, nlines2)
		no_differences = True
		for i in range(minlines):
			if i % 10000 == 0:
				print (str(i)+" lines checked")
			l1split = lines1[i].split()
			l2split = lines2[i].split()
			if (l1split[0].strip() != l2split[0].strip()) or (l1split[1].strip() != l2split[1].strip()):
				#if int(l1split[1].strip()) >= 255:
				#	continue
				print("Difference at line " + str(i))
				print(lines1[i])
				print(lines2[i])
				print("=============================================\n")
				no_differences = False
		if no_differences:
			print("Outputs were identical")
		else:
			print("Outputs were different")
		


if __name__=="__main__":
	main(sys.argv[1], sys.argv[2])
