import csv
import numpy as np

def readFunc(path, layout = "colMajor", delimiter = ","):
    '''
    read a file and store it in an np.array. default input file is in col major

    input:
    path: the file path

    output:
    funcX(n * m): n * m matrix with the first column label, second row represents x and the rests are y(x)
    '''
    func = np.array([])
    dataTable = []
    with open(path, 'rb') as csvfile:
        dataRaw = csv.reader(csvfile, delimiter = delimiter) # read from file
        # Read the whole file
        for row in dataRaw:
            dataTable.append(row)
        # If rowMajor, transpose
        if(layout == "rowMajor"):
            data = list(map(list, zip(*data)))
        else: data = dataTable
        # read
        label = data[0][:]
        func = np.array(data[1:][:], dtype = np.float64)
    return (label, func)

def readFuncMul(path, delimiter = ","):
    xs = []
    ys = []
    labels = []
    with open(path, 'rb') as csvfile:
        data = csv.reader(csvfile, delimiter = delimiter) # read from file
        i = 0
        for row in data:
            if(i % 3 == 0): labels.append(row)
            elif(i % 3 == 1): xs.append(np.array(row, dtype = np.float64))
            else: ys.append(np.array(row, dtype = np.float64))
    return (labels, xs, ys)
