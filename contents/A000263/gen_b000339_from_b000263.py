import os

content = []
with open("b000263.txt") as f :
	for line in f.readlines() :
		line = line.strip()
		if line == "": continue
		x, y = map(int, line.split(' '))
		content.append((x, y))

if os.path.exists("b000339.txt") :
	print("The output file b000339.txt exists. Please rename or remove it before running this script.")
	exit(1)

with open("b000339.txt", "w") as f :
	for x, y in content :
		f.write("{} {}\n".format(x, y + (x * x // 4)))
