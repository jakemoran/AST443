import numpy as np
import matplotlib
import requests
import math
import csv

from matplotlib import pyplot as plt

def get_header_index(array, header_name):

    count = 0
    for element in array:
        if element == header_name:
            return count
        count += 1


'''
We use this function to determine if an entry has a valid number. We cannot use
isnumeric() or isdigit() as those only work for integers. If the last character
in the string is a number, consider the whole string a number.
'''
def is_number(string):
    if string[-1:].isdigit():
        return True
    return False


# Open the exoplanet data
with open('exoplanet_catalog.csv') as exoplanets:
    csv_reader = csv.reader(exoplanets, delimiter=',')
    exoplanets_array = []
    for line in csv_reader:
        exoplanets_array.append(line)

# Number of planets: 4837

'''
First we filter by detection type, filtering out all planets that
were not first detected by transit. By doing this we only include planets that
are capable of being observed transiting.
'''
detection_index = get_header_index(exoplanets_array[0], 'detection_type')
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    if exoplanets_array[count - num_removed + 1][detection_index] != 'Primary Transit':
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 3451

'''
Next we have to ensure that the host star is bright enough. It must have an apparent
magnitude in the V-band of at most 12.
'''

mag_index = get_header_index(exoplanets_array[0], 'mag_v')
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    if not is_number(exoplanets_array[count - num_removed + 1][mag_index]) or \
        float(exoplanets_array[count - num_removed + 1][mag_index]) > 12.0:
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 432

'''
Only consider planets that have a max altitude of 40 degrees as seen from Stony Brook (corresponds
to a declination between -9 and 91).
'''

dec_index = get_header_index(exoplanets_array[0], 'dec')
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    if not is_number(exoplanets_array[count - num_removed + 1][dec_index]) or \
        -9.0 > float(exoplanets_array[count - num_removed + 1][dec_index]) < 91.0:
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 263

'''
Only consider stars that culminate between 10:30 pm and 3:30 am local time.
'''

ra_index = get_header_index(exoplanets_array[0], 'ra')
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    if not is_number(exoplanets_array[count - num_removed + 1][ra_index]) or \
        300 > float(exoplanets_array[count - num_removed + 1][ra_index]) > 10:
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 47

for planet in exoplanets_array:
    print(planet[0], planet[detection_index], planet[mag_index], planet[dec_index], planet[ra_index])

print(len(exoplanets_array)-1)
