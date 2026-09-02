import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

user_input = input("Enter a time value:")
time = float(user_input)
matrix_size = 3
data = []
dt = 0
filename = "test3"

def check_empty():
    with open(filename, "r") as f:
        line = f.readline()
        isNotEmpty1 = False
        for entries in line.split():
            if entries == "results":
                isNotEmpty1 = True
        if not isNotEmpty1:
            print("There isn't a header or data. Run the program")
            sys.exit()
        while "RDM results" not in line:
            line = f.readline()
        line = next(f)
        line = next(f)
        line = next(f)
        test = line.split()
        isNotEmpty = False
        for entries in test:
            if (entries == "time"):
                isNotEmpty = True
        if not isNotEmpty:
            print("The data is empty but there is a header. Did you run the program?")
            sys.exit()

check_empty()

def get_dt() :
     with open(filename, "r") as f:
        line = f.readline()
        while "dt" not in line: #finds Dt and writes dt
                line = f.readline()
        temp2 = line.split()
        dt1 = float(temp2[2])
        return dt1

def get_size():
    with open(filename, "r") as f:
        line = f.readline()
        while "levels" not in line: #finds levels and writes levels as Matrix_size
                line = f.readline()
        temp = line.split()
        matrix_size = int(temp[2])
        return matrix_size
    
def get_data(dt, time):
    if not os.path.exists(filename):
        print("Error: File not found")
        return None
    else:
        with open(filename, "r") as f:
            line = f.readline()
            while "RDM results" not in line: #looks for RDM results andn skips the header of the file
                line = f.readline()
            line = next(f)
            line = next(f) 
            line = next(f)
            line = next(f)
            line = next(f)
            times = line.split()
            error = abs(float(times[0]) - time)
            lastrow = times
            for line in f:
                if not times:
                    break
                times = line.split()
                if error > abs(float(times[0]) -  time):
                    error = abs(float(times[0]) - time)
                    lastrow = times
                else:
                    return lastrow
            return lastrow


def get_ReIm(matrix_size, time):
    if not os.path.exists(filename):
        print("Error: File not found")
        return None
    else:
        with open(filename, "r") as f:
            line = f.readline()
            while "RDM results" not in line: #looks for RDM results andn skips the header of the file
                line = f.readline()
            line = next(f)
            line = next(f)  # might need to remove one of these
            line = next(f)
            parts = line.split()
            reVsIm = []
            for part in parts:   # make a map that tells me if i x mod n^2 + 2 is imag or real or junk
                if "imag" in part:   #imag is 1
                    reVsIm.append(1)        
                elif "norm" in part:
                    reVsIm.append(2)
                elif "time" in part: #junk is 2
                    reVsIm.append(2)
                elif "real" in part:  #real is 0
                    reVsIm.append(0)
    return reVsIm

matrix_size = get_size()
dt = get_dt()
ReIm = get_ReIm(matrix_size, time)
data =  get_data(dt, time)

def sort_real(data, indexer):
    real = []
    for x in range(len(data)):
        if (indexer[x] == 0):
            real.append(float(data[x]))
    return real
def sort_im(data, indexer):
    im = []
    for x in range(len(data)):
        if (indexer[x] == 1):
            im.append(float(data[x]))
    return im

real = sort_real(data, ReIm)
im = sort_im(data, ReIm)

def sorterReal(real, matrix_size):
    real_Matrix = np.zeros((matrix_size,matrix_size))
    y = matrix_size
    i = 0
    for x in range(y):
        for z in range(x + 1):
            real_Matrix[x,z] = real[i]
            real_Matrix[z,x] = real_Matrix[x,z]
            i = i + 1
    return real_Matrix
    

def sorterIm(im, matrix_size):
    im_Matrix = np.zeros((matrix_size,matrix_size))
    y = matrix_size
    i = 0
    for x in range(y):
        for z in range(x + 1):
            if (x == z) :
                im_Matrix[x,z] = 0 
            else:
                im_Matrix[x ,z] = im[i]
                im_Matrix[z,x] = -1 * im_Matrix[x,z]
                i = i + 1
    return im_Matrix

def find_largestRe(real):
    largest = real[0]
    for x in range(len(real)):
        if (largest < real[x]):
            largest = real[x]
    return largest
def find_LargestIm(im):
    largest = im[0]
    for x in range(len(im)):
        if (largest < im[x]):
            largest = im[x]
    return largest 

t = str(data[0])

im_Matrix = sorterIm(im, matrix_size)
real_Matrix = sorterReal(real, matrix_size)

indices = np.arange(matrix_size)
x, y = np.meshgrid(indices, indices)
a = 1*10**-7
LRE = -find_largestRe(real) - a
UPRE = find_largestRe(real) + a
LIM = -find_LargestIm(im) - a
UPIM = find_LargestIm(im) + a

size = (1/matrix_size) * 2500   #this define dots size larger = bigger can be changed if wanted to 

fig, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6))  #The Real part of RDM is symmetric wherea the Im part is antisymmetric
ax1.set_title("Real Time:" + t)
ax2.set_title("Imaginary Time:" + t)

scat1 = ax1.scatter(x,y, c = real_Matrix, s = size, cmap = 'bwr_r', vmin = -find_largestRe(real), vmax = find_largestRe(real), alpha =1)#This is for real. here if you want to make the bounds flaconstant change vmin and vmax
plt.tight_layout()
scat2 = ax2.scatter(x,y, c = im_Matrix, s = size, cmap = 'PRGn_r', vmin = -find_LargestIm(im), vmax = find_LargestIm(im), alpha =1) #This is for real. here if you want to make the bounds flaconstant change vmin and vmax
plt.tight_layout()
plt.colorbar(scat2, ax = ax2)
plt.colorbar(scat1, ax = ax1)
plt.show()
