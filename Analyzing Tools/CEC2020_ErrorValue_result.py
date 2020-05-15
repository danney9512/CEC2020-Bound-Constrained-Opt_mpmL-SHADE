import os
import csv


filename = ""
D = ["05", "10", "15", "20"]
fun_num = 10

alg_name = input("請輸入演算法名稱: ")
alg_exp = input("請輸入演算法實驗:")

for d in D:

    output_filename = alg_name + "+" + alg_exp + "_result_D" + d + ".csv"
    output_file = open(output_filename, 'a', newline='')
    writer = csv.writer(output_file)
    writer.writerow(["Func.", "Best", "Worst", "Median", "Mean", "Std"])

    for i in range(1,fun_num+1, 1):
        if d=="05" and (i==6 or i==7):
            continue

        if i <10:
            filename = alg_name + "_F0" + str(i) + "_D" + d + ".csv"
        else:
            filename = alg_name + "_F" + str(i) + "_D" + d + ".csv"
        
        try:
            file = open(filename, 'r')
            data = csv.reader(file)
            for row in data:
                writer.writerow(row)
        except IOError:
            print(filename + "讀取失敗")
        finally:
            file.close()


        filename = ""
    
    output_file.close()


