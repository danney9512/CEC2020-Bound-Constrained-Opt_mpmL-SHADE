import os
import pandas as pd


def filename_proccess(o_str):
    
    text = o_str.split("_")
    
    F = text[1]
    D = text[2]
    D = D.replace("D", "")
    return F, D


#
target_filename = "ErrorValue.csv"

 # get all files' and folders' names in the current directory
filenames= os.listdir (".")

folders = []

for folder in filenames:
    if os.path.isdir(os.path.join(os.path.abspath("."), folder)): # check whether the current object is a folder or not
        folders.append(folder)


# Create a new directory folder for output
dirName = '30Run results'
 
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")


# Loop all folders to merge the data
for folder in folders:
    
    runs_col = []
    for numofRun in range(0,30):
        runs_col.append("Run" + str(numofRun + 1))
    
    data_frame = pd.DataFrame(columns=['NFE'].append(runs_col))
    
    for numofRun in range(0,30):
        path = folder + "/" + str(numofRun) + "/"
        
        """        
        try:
            file = open(path + target_filename, "r")
            data = csv.reader(file)
            for row in data:
                print(row)
        except:
            print(path + target_filename + " 讀取失敗")
        finally:
            file.close()
        """
        
        data_frame_tmp = pd.read_csv(path + target_filename, header=None)
        
        if numofRun == 0:
            data_frame["NFE"] = data_frame_tmp[0]
            
        data_frame["Run" + str(numofRun+1)] = data_frame_tmp[1]
        
    data_frame.set_index("NFE",inplace=True)
    print(data_frame)
    
    tab = '\t'

    F, D = filename_proccess(folder)
    data_frame.to_csv("./" + dirName + "/" + "mpmL-SHADE_" + F + "_" + D + ".csv")
    data_frame.to_csv("./" + dirName + "/" + "mpmL-SHADE_" + F + "_" + D + ".txt", header=False, index=False, sep=tab, float_format='%19.8f')
        
        