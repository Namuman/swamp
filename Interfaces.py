import os

def get_comppropDB_names(directory: str):
    all_filenames = os.listdir(directory)
    comppropDB_names_list = []
    for filename in all_filenames:
        if filename[-5:] == '.xlsx' and 'CompProp' in filename:
            comppropDB_names_list.append(filename)
    return comppropDB_names_list
## update to use DB files outside from main.py directory

def get_binarycoefDB_names(directory: str):
    all_filenames = os.listdir(directory)
    binarycoefDB__names_list = []
    for filename in all_filenames:
        if filename[-5:] == '.xlsx' and 'BinaryCoef' in filename:
            binarycoefDB__names_list.append(filename)
    return binarycoefDB__names_list
## update to use DB files outside from main.py directory

def get_streamcomp_names(directory: str):
    all_filenames = os.listdir(directory)
    streamcomp_names_list = []
    for filename in all_filenames:
        if filename[-5:] == '.xlsx' and 'StreamComposition' in filename:
            streamcomp_names_list.append(filename)
    return streamcomp_names_list

def select_file(filenames_list: list):
   print('\tSelect file from listed below:')
   for i in range(len(filenames_list)):
       print('\t{} -\t{}'.format(i + 1, filenames_list[i]))
   print('\tInput file name number:')
   selected_number = int(input())
   selected_file = filenames_list[selected_number - 1]
   print('\tSELECTED FILE:\t{}\n'.format(selected_file))
   return selected_file
