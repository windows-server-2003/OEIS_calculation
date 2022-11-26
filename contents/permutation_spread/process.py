

spread0 = ""
spread1 = ""
spread2 = ""
spread3 = ""
table_all = ""
n = 1
table_all_n = 1
with open("stdout.txt") as f :
	for line in f.readlines() :
		values = list(map(int, line.replace('[', '').replace(']', '').split(' ')))
		if n >= 1 : spread0 += str(n) + " " + str(values[0]) + "\n"
		if n >= 2 : spread1 += str(n) + " " + str(values[1]) + "\n"
		if n >= 3 : spread2 += str(n) + " " + str(values[2]) + "\n"
		if n >= 4 : spread3 += str(n) + " " + str(values[3]) + "\n"
		for j in values :
			table_all += str(table_all_n) + " " + str(j) + "\n"
			table_all_n += 1
		n += 1

with open("b004204.txt", "w") as f : f.write(spread0)
with open("b004205.txt", "w") as f : f.write(spread1)
with open("b004206.txt", "w") as f : f.write(spread2)
with open("b004246.txt", "w") as f : f.write(spread3)
with open("b147679.txt", "w") as f : f.write(table_all)

