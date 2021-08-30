import numpy as np
import math
import csv
import datetime as dt

from datetime import datetime, timezone

now = datetime.now(timezone.utc)
this_year = now.year
this_month = now.month
this_day = now.day + (now.hour * 3600 + now.minute * 60 + now.second) / 86400

# Function that converts a time in UTC to Julian time
def getJD(year=this_year, month=this_month, day=this_day):

    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month

    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
            (year == 1582 and month < 10) or
            (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)

    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)

    D = math.trunc(30.6001 * (monthp + 1))

    jd = B + C + D + day + 1720994.5

    return jd


# Function that converts Julian time to a date in UTC
def jd_to_UTC(jd):

    jd = jd + 0.5

    F, I = math.modf(jd)
    I = int(I)

    A = math.trunc((I - 1867216.25) / 36524.25)

    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I

    C = B + 1524

    D = math.trunc((C - 122.1) / 365.25)

    E = math.trunc(365.25 * D)

    G = math.trunc((C - E) / 30.6001)

    day = C - E + F - math.trunc(30.6001 * G)

    if G < 13.5:
        month = G - 1
    else:
        month = G - 13

    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    return year, month, day

def calc_transit_time(sma, star_radius, planet_radius, period):

    arc_length = 2*sma * math.asin(star_radius/sma)
    circumference = 2*math.pi*sma
    transit_proportion = (arc_length + 2*planet_radius)/circumference

    return transit_proportion * period

# Get the index of a certain header in the exoplanet data
def get_header_index(array, header_name):

    count = 0
    for element in array:
        if element == header_name:
            return count
        count += 1


'''
We use this function to determine if an entry has a valid number (some are blank). We cannot use
isnumeric() or isdigit() as those only work for integers. If the last character in the string is 
a number, consider the whole string a number.
'''
def is_number(string):
    if string[-1:].isdigit():
        return True
    return False


# Open the exoplanet data and store it all in exoplanets_array
with open('exoplanet_catalog.csv') as exoplanets:
    csv_reader = csv.reader(exoplanets, delimiter=',')
    exoplanets_array = []
    for line in csv_reader:
        exoplanets_array.append(line)

# Number of planets: 4837

# Save important header indecies for later
detection_index = get_header_index(exoplanets_array[0], 'detection_type')
mag_index = get_header_index(exoplanets_array[0], 'mag_v')
dec_index = get_header_index(exoplanets_array[0], 'dec')
ra_index = get_header_index(exoplanets_array[0], 'ra')
radius_index = get_header_index(exoplanets_array[0], 'radius')
star_radius_index = get_header_index(exoplanets_array[0], 'star_radius')
tconj_index = get_header_index(exoplanets_array[0], 'tconj')
tzero_index = get_header_index(exoplanets_array[0], 'tzero_tr')
period_index = get_header_index(exoplanets_array[0], 'orbital_period')
sma_index = get_header_index(exoplanets_array[0], 'semi_major_axis')

for planet in exoplanets_array:
    if planet[0] == 'WASP-52 b':
        print(planet[period_index])

'''
First we filter by detection type, filtering out all planets that
were not first detected by transit. By doing this we only include planets that
are capable of being observed transiting.
'''
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    # If the planet was not detected by transit, remove it from the array
    if exoplanets_array[count - num_removed + 1][detection_index] != 'Primary Transit':
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 3451

'''
Next we have to ensure that the host star is bright enough. It must have an apparent
magnitude in the V-band of at most 12.
'''
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    # Remove the planet if the v-band magnitude is not in the database or is outside the range
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
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    # Remove the planet if the declination is not in the database or is outside the range
    if not is_number(exoplanets_array[count - num_removed + 1][dec_index]) or \
        -9.0 > float(exoplanets_array[count - num_removed + 1][dec_index]) < 91.0:
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 263

'''
Only consider stars that culminate in the middle of the night (approximately 11:30 pm to
3:30 am), corresponding to a right ascension between 270 and 330 degrees.
'''
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    # Remove the planet if the RA is not in the database or is outside the range
    if not is_number(exoplanets_array[count - num_removed + 1][ra_index]) or \
        7.5 < float(exoplanets_array[count - num_removed + 1][ra_index]) < 295:
        del exoplanets_array[count - num_removed + 1]
        num_removed += 1
    count += 1

# Number of planets remaining: 100

'''
Next we filter by the relative sizes of the planet and star to ensure that enough
of the star's light will be blocked. The cross section ratio must be at least 0.008.
'''
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

while count < num_planets:
    # Remove planet if the radius of the planet or star is not in the database
    if not is_number(exoplanets_array[count - num_removed + 1][radius_index]) or \
        not is_number(exoplanets_array[count - num_removed + 1][star_radius_index]):

        del exoplanets_array[count - num_removed + 1]
        num_removed += 1

    else:
        # Convert both radii to meters
        radius = float(exoplanets_array[count - num_removed + 1][radius_index]) * 6.991e7
        star_radius = float(exoplanets_array[count - num_removed + 1][star_radius_index]) * 6.96e8

        # Remove planet if it does not block enough of the star
        if (radius / star_radius)**2 < 0.008:
            del exoplanets_array[count - num_removed + 1]
            num_removed += 1

    count += 1

# Number of planets remaining: 23

'''
Finally we predict the transit times based on the Julian time of the first transit and the
period of the orbit. We compare these times to a list of observing times throughout the
month of September.
'''
num_planets = len(exoplanets_array) - 1
count = 0
num_removed = 0

observing_times = []
list_of_transits = []

# Consider the first observing start time to be the night of September 1
first_start_jd = getJD(2021, 9, 2.067)

# Define the length of an observation session (4 hours)
time_range = 0.16712
start_jd = first_start_jd

for i in range(1, 30):
    observing_times.append([start_jd, start_jd+time_range])
    start_jd += 1.00272

while count < num_planets:

    # Discard planets that don't have available transit time, period, or semi major axis
    if not is_number(exoplanets_array[count - num_removed + 1][tconj_index]) and \
        not is_number(exoplanets_array[count - num_removed + 1][tzero_index]) or \
        not is_number(exoplanets_array[count - num_removed + 1][period_index]) or \
        not is_number(exoplanets_array[count - num_removed + 1][sma_index]):

        del exoplanets_array[count - num_removed + 1]
        valid = False
        num_removed += 1

    elif is_number(exoplanets_array[count - num_removed + 1][tconj_index]):
        first_transit = float(exoplanets_array[count - num_removed + 1][tconj_index])
        valid = True
    else:
        first_transit = float(exoplanets_array[count - num_removed + 1][tzero_index])
        valid = True

    # If there is an available transit time, calculate when the last transit was and
    # append the next 30 transits to a list
    if valid:

        period = float(exoplanets_array[count - num_removed + 1][period_index])
        name = exoplanets_array[count - num_removed + 1][0]
        ra = float(exoplanets_array[count - num_removed + 1][ra_index])
        dec = float(exoplanets_array[count - num_removed + 1][dec_index])
        radius = float(exoplanets_array[count - num_removed + 1][radius_index]) * 6.991e7
        star_radius = float(exoplanets_array[count - num_removed + 1][star_radius_index]) * 6.96e8
        sma = float(exoplanets_array[count - num_removed + 1][sma_index]) * 1.5e11
        transit_signature = (radius / star_radius)**2

        last_transit = (math.floor((getJD() - first_transit) / period) * period) + first_transit

        for i in range(1, 50):
            list_of_transits.append([name, ra, dec, last_transit + i*period, transit_signature,
                                     calc_transit_time(sma, star_radius, radius, period*24)])

    count += 1

usable_transits = []

for transit in list_of_transits:
    for time in observing_times:
        if time[0] < transit[3] < time[1]:
            usable_transits.append(transit)

for transit in usable_transits:
    print(f'Name: {transit[0]}, RA: {transit[1]}, DEC: {transit[2]}, '
          f'Transit Time (JD): {transit[3]}, Transit Time (UTC): {jd_to_UTC(transit[3])}, '
          f'Transit Signature: {transit[4]}, Transit Duration (hr): {transit[5]}')

print(len(usable_transits))

print(getJD())
