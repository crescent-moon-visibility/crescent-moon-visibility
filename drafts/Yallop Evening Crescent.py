#!/usr/bin/env python3
# MIT, @ebraminio and @hidp123

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.image as mpimg
import math
import numpy as np
import astronomy # pip install astronomy-engine
#import pandas
#from tqdm import tqdm
from mpire import WorkerPool # pip install mpire



KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit.

def calculate(base_time, latitude, longitude, STEPS):
    
    if not -60 <= latitude <= 60: return {}
    
    observer = astronomy.Observer(latitude, longitude)
    time = base_time.AddDays(-observer.longitude / 360) # this corrects the base time based on timezone
    sunset   = astronomy.SearchRiseSet(astronomy.Body.Sun,  observer, astronomy.Direction.Set, time, 1)
    moonset  = astronomy.SearchRiseSet(astronomy.Body.Moon, observer, astronomy.Direction.Set, time, 1)
    
    if sunset is None or moonset is None: return {}
    #if sunset.ut > moonset.ut: q_code = 'G'

    # https://astro.ukho.gov.uk/moonwatch/background.html
    # lag time: The time interval between sunset and moonset. The lag time is usually
    # given in minutes. It can be negative, indicating that the Moon sets before the Sun.
    if lag_time < 0:
        best_time = sunset
    else: best_time = astronomy.Time(sunset.ut + lag_time * 4/9)

    # best time: an empirical prediction of the time which gives the observer the best opportunity
    # to see the new crescent Moon (Sunset time + (4/9)*Lag time).
    best_time = astronomy.Time(sunset.ut + lag_time * 4/9)

    new_moon_prev = astronomy.SearchMoonPhase(0, sunset, -35)
    new_moon_next = astronomy.SearchMoonPhase(0, sunset, 35)
    
    if (sunset.ut - new_moon_prev.ut) <= (new_moon_next.ut - sunset.ut):
       new_moon_nearest = new_moon_prev
    else: new_moon_nearest = new_moon_next

    moon_age_to_next_moon = best_time.ut - new_moon_next.ut # moon age at best time.
    moon_age_to_prev_moon = best_time.ut - new_moon_prev.ut # moon age at best time.
    moon_age_to_nearest_new_moon = best_time.ut - new_moon_nearest.ut # moon age in days.

    # if (moon_age_to_nearest_moon * 24) % 1 < STEPS / 15: return {"q_code": 'I'} # to plot moon age lines at every hour.
    if (moon_age_to_nearest_moon * 24) % 1 < STEPS / 30 or (moon_age_to_nearest_moon * 24) % 1 > (1 - STEPS / 30): return {"q_code": 'I'} # to plot moon age lines at every hour.

    if lag_time < 0 and sunset.ut < new_moon_nearest.ut: return {"q_code": 'J'}
    if lag_time < 0: return {"q_code": 'G'} # i.e. moonset is before sunset
    if sunset.ut < new_moon_nearest.ut: return {"q_code": 'H'} # i.e. conjunction after sunset

    sun_equator = astronomy.Equator(astronomy.Body.Sun, best_time, observer, True, True)
    #sun_distance = KM_PER_AU * sun_equator.dist #topocentric
    sun_horizon = astronomy.Horizon(best_time, observer, sun_equator.ra, sun_equator.dec, astronomy.Refraction.Airless)
    sun_alt = sun_horizon.altitude
    sun_az = sun_horizon.azimuth

    moon_elongation_event = astronomy.Elongation(astronomy.Body.Moon, best_time) #geocentric elongation
    moon_equator = astronomy.Equator(astronomy.Body.Moon, best_time, observer, True, True) #RA is in h.dd (hours.degrees)
    #moon_distance = KM_PER_AU * moon_equator.dist #topocentric
    #moon_distance2 = astronomy.Libration(best_time).dist_km #AU_IN_M * moon_equator.vec.Length()
    moon_horizon = astronomy.Horizon(best_time, observer, moon_equator.ra, moon_equator.dec, astronomy.Refraction.Airless)
    moon_alt = moon_horizon.altitude
    moon_az = moon_horizon.azimuth

    # https://github.com/rob-blackbourn/PyFinance/blob/2bbad39b/py_calendrical/location.py#L217
    #lunar_parallax = 6378140 / moon_distance * math.cos(math.radians(moon_alt)) #lunar_parallax in radians.
    #lunar_parallax = math.degrees(lunar_parallax_RAD)

    # https://github.com/abdullah-alhashim/prayer_calculator/blob/8abe558/moon_sighting.py#L54-L62
    #HP = lunar_parallax / math.cos(math.radians(moon_alt))
    SD = astronomy.Libration(best_time).diam_deg * 60 / 2 #semi-diameter of the Moon in arcminutes, geocentric
    lunar_parallax = SD/0.27245 #in arcminutes
    #SD = 0.27245 * HP * (180 * 60 / math.pi)        # semi-diameter of the Moon
    SD_topo = SD * (1 + (math.sin(math.radians(moon_alt)) * math.sin(math.radians(lunar_parallax/60)))) #in arcminutes. Here SD is in arcminutes, moon_alt in degrees, lunar_parallax in degrees (that's why it has been divided by 60).

    # https://github.com/abdullah-alhashim/prayer_calculator/blob/8abe558/moon_sighting.py#L71-L77
    ARCL = moon_elongation_event.elongation #in degrees, geocentric
    DAZ = sun_az - moon_az

    #COSARCV = math.cos(math.radians(ARCL))/math.cos(math.radians(DAZ))

    #if -1 <= COSARCV <= 1: ARCV = math.degrees(math.acos(COSARCV)) #math.degrees(math.acos(math.cos(math.radians(ARCL))/math.cos(math.radians(DAZ)))) #moon_alt - sun_alt
    #elif COSARCV < -1: ARCV = math.degrees(math.acos(-1))
    #elif COSARCV > 1: ARCV = math.degrees(math.acos(1)) 

    #ARCV
    geomoon = astronomy.GeoVector(astronomy.Body.Moon, best_time, True)
    geosun = astronomy.GeoVector(astronomy.Body.Sun, best_time, True)
    rot = astronomy.Rotation_EQJ_EQD(best_time)
    rotmoon = astronomy.RotateVector(rot, geomoon)
    rotsun  = astronomy.RotateVector(rot, geosun)
    meq = astronomy.EquatorFromVector(rotmoon)
    seq = astronomy.EquatorFromVector(rotsun)
    mhor = astronomy.Horizon(best_time, observer, meq.ra, meq.dec, astronomy.Refraction.Airless)
    shor = astronomy.Horizon(best_time, observer, seq.ra, seq.dec, astronomy.Refraction.Airless)
    ARCV = mhor.altitude - shor.altitude

    W_topo = SD_topo * (1 - (math.cos(math.radians(ARCL)))) #in arcminutes
    
    q = (ARCV - (11.8371 - 6.3226*W_topo + 0.7319*W_topo**2 - 0.1018*W_topo**3)) / 10

    if q > +0.216: q_code = 'A' # Crescent easily visible
    elif +0.216 >= q > -0.014: q_code = 'B' # Crescent visible under perfect conditions
    elif -0.014 >= q > -0.160: q_code = 'C' # May need optical aid to find crescent
    elif -0.160 >= q > -0.232: q_code = 'D' # Will need optical aid to find crescent
    elif -0.232 >= q > -0.293: q_code = 'E' # Crescent not visible with telescope
    elif -0.293 >= q: q_code = 'F'

    return {
        "lat": observer.latitude,
        "long": observer.longitude,
        "sunset": str(sunset.Utc()),
        "moonset": str(moonset.Utc()),
        "lag time": lag_time,
        "best time": str(best_time.Utc()),
        #"sun dist": sun_distance,
        "sun alt": sun_alt,
        "sun az": sun_az,
        #"moon dist": moon_distance,
        "moon alt": moon_alt,
        "moon az": moon_az,
        "lunar parllax": lunar_parallax,
        #"HP": HP,
        "SD": SD,
        "SD_topo": SD_topo,
        "ARCV": ARCV,
        "DAZ": DAZ,
        "ARCL": ARCL,
        "W_topo": W_topo,
        "q_code": q_code,
        "q": q
    }

def run(base_time):
    result = []
    STEPS = 5
    latitudes = np.arange(-90, 90, STEPS)
    longitudes = np.arange(-180, 180, STEPS)
    args_list = [(base_time, - lat, lng, STEPS) for lat in latitudes for lng in longitudes]

    with WorkerPool(n_jobs=4) as pool:
        result = list(pool.map(calculate, args_list, progress_bar=True))

    # Define the color mapping
    color_mapping = {
        'A': "#83c702",
        'B': "#709a08",
        'C': "#416100",
        'D': "gold",
        'E': "orange",
        'F': (0, 0, 0, 0),
        'G': "red",
        'H': "purple",
        'I': "gray",
        'J': (1, 0, 0, 0.5),  # semi-transparent red        
    }

    # Create a colormap
    cmap = ListedColormap(list(color_mapping.values()))

    H = np.zeros((len(latitudes), len(longitudes), 4))  # Use 4 instead of 3
    for i, r in enumerate(result):
        if 'q_code' in r:
            color = matplotlib.colors.to_rgba(color_mapping[r['q_code']])  # Use to_rgba instead of to_rgb
            H[i // len(longitudes), i % len(longitudes)] = color

    plt.figure(figsize=(20, 10))  # Size is in inches

    # Load the image
    img = mpimg.imread('..map.png') # enter your file path for map image.

    # Display the image
    plt.imshow(img, extent=[-180, 180, -90, 90])

    # Plot the data on the map with some transparency
    plt.imshow(H, cmap=cmap, extent=[-180, 180, -90, 90], alpha=0.6)

    # Add ticks and labels on both sides
    plt.tick_params(axis='both', direction='inout', left=True, right=True, top=True, bottom=True, labelleft=True, labelright=True, labeltop=True, labelbottom=True)

    # Add horizontal lines at y=60 and y=-60
    plt.axhline(60, color='gray', linewidth=1.0)
    plt.axhline(-60, color='gray', linewidth=1.0)
    
    # plt.xticks([])  # Remove x-axis ticks
    # plt.yticks([])  # Remove y-axis ticks

    plt.show()

if __name__ == '__main__':
    run(astronomy.Time.Make(2022, 8, 27, 0, 0, 0))

# def run(base_time):
#     result = []
#     STEPS = 10
#     H = numpy.ndarray(shape=(int(180 / STEPS), int(360 / STEPS)))
#     for lng in tqdm(numpy.arange(0, 360, STEPS)):
#         for lat in numpy.arange(0, 180, STEPS):
#             r = calculate(base_time, 90 - lat, lng - 180)
#             if 'q_code' in r:
#                 result.append(r)
#                 H[int(lat / STEPS), int(lng / STEPS)] = ord('F') - ord(r['q_code'])
#             else:
#                 H[int(lat / STEPS), int(lng / STEPS)] = 0

#     #result = pandas.DataFrame.from_records(result)
#     #result.to_excel("output_" + str(base_time.Utc()) + ".xlsx")
#     plt.imshow(H)
#     #plt.savefig(str(base_time.Utc()) + '.png')
#     plt.show() # uncomment this when running inside Jupyter Notebook / Google Colab

# run(astronomy.Time.Make(2022, 8, 27, 0, 0, 0))
# run(astronomy.Time.Make(2022, 8, 28, 0, 0, 0))
# run(astronomy.Time.Make(2022, 8, 29, 0, 0, 0))

# run(astronomy.Time.Make(2022, 6, 29, 0, 0, 0))
#run(astronomy.Time.Make(2022, 6, 30, 0, 0, 0))
#run(astronomy.Time.Make(2022, 8, 29, 0, 0, 0))
