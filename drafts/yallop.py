#!/usr/bin/env python3
# MIT, @ebraminio and @hidp123

import matplotlib.pyplot as plt
import math
import numpy
import astronomy # pip install astronomy-engine
import pandas
from tqdm import tqdm

#KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit.

def calculate(base_time, latitude, longitude):
    observer = astronomy.Observer(latitude, longitude)
    time = base_time.AddDays(-longitude / 360) # this corrects the base time based on timezone
    sunset   = astronomy.SearchRiseSet(astronomy.Body.Sun,  observer, astronomy.Direction.Set, time, True)
    moonset  = astronomy.SearchRiseSet(astronomy.Body.Moon, observer, astronomy.Direction.Set, time, True)
    if sunset is None or moonset is None: return

    # https://astro.ukho.gov.uk/moonwatch/background.html
    # lag time: The time interval between sunset and moonset. The lag time is usually
    # given in minutes. It can be negative, indicating that the Moon sets before the Sun.
    lag_time = moonset.ut - sunset.ut
    if lag_time < 0: return {"q_code": "G"}

    # best time: a semi empirical prediction of the time which gives the observer the best opportunity
    # to see the new crescent Moon according to Yallop (Sunset time + (4/9)*Lag time).
    best_time = astronomy.Time(sunset.ut + lag_time * 4/9)

    sun_equator = astronomy.Equator(astronomy.Body.Sun, best_time, observer, True, True)
    #sun_distance = KM_PER_AU * sun_equator.dist #topocentric
    sun_horizon = astronomy.Horizon(best_time, observer, sun_equator.ra, sun_equator.dec, astronomy.Refraction.Airless)
    sun_alt = sun_horizon.altitude
    sun_az = sun_horizon.azimuth

    moon_elongation_event = astronomy.Elongation(astronomy.Body.Moon, best_time) #geocentric elongation event. Gives three outputs.
    moon_equator = astronomy.Equator(astronomy.Body.Moon, best_time, observer, True, True) #RA is in h.dd (hours.degrees)
    #moon_distance = KM_PER_AU * moon_equator.dist #topocentric
    #moon_distance2 = astronomy.Libration(best_time).dist_km #AU_IN_M * moon_equator.vec.Length()
    moon_horizon = astronomy.Horizon(best_time, observer, moon_equator.ra, moon_equator.dec, astronomy.Refraction.Airless)
    moon_alt = moon_horizon.altitude
    moon_az = moon_horizon.azimuth

    # https://github.com/rob-blackbourn/PyFinance/blob/2bbad39b/py_calendrical/location.py#L217
    #lunar_parallax = 6378140 / moon_distance * math.cos(math.radians(moon_alt)) #lunar_parallax in radians.
    #lunar_parallax = math.degrees(lunar_parallax_RAD)


    #HP = lunar_parallax / math.cos(math.radians(moon_alt))
    SD = astronomy.Libration(best_time).diam_deg * 60 / 2 #semi-diameter of the Moon in arcminutes, geocentric
    lunar_parallax = SD/0.27245 #in arcminutes
    #SD = 0.27245 * HP * (180 * 60 / math.pi)        # semi-diameter of the Moon
    SD_topo = SD * (1 + (math.sin(math.radians(moon_alt)) * math.sin(math.radians(lunar_parallax/60)))) #in arcminutes. Here SD is in arcminutes, moon_alt in degrees, lunar_parallax in degrees (that's why it has been divided by 60).


    ARCL = moon_elongation_event.elongation #in degrees
    DAZ = sun_az - moon_az
    ARCV = math.degrees(math.acos(math.cos(math.radians(ARCL))/math.cos(math.radians(DAZ))))

    W_topo = SD_topo * (1 - (math.cos(math.radians(ARCL)))) #in arcminutes
    q = (ARCV - (11.8371 - 6.3226*W_topo + 0.7319*W_topo**2 - 0.1018*W_topo**3)) / 10

    if q > +0.216: q_code = 'A' # Crescent easily visible
    elif +0.216 >= q > -0.014: q_code = 'B' # Crescent visible under perfect conditions
    elif -0.014 >= q > -0.160: q_code = 'C' # May need optical aid to find crescent
    elif -0.160 >= q > -0.232: q_code = 'D' # Will need optical aid to find crescent
    elif -0.232 >= q > -0.293: q_code = 'E' # Crescent not visible with telescope
    elif -0.293 >= q: q_code = 'F'

    return q_code
    #return {
    #    "lat": observer.latitude,
    #    "long": observer.longitude,
    #    "sunset": str(sunset.Utc()),
    #    "moonset": str(moonset.Utc()),
    #    "lag time": lag_time,
    #    "best time": str(best_time.Utc()),
    #    "sun dist": sun_distance,
    #    "sun alt": sun_alt,
    #    "sun az": sun_az,
    #    "moon dist": moon_distance,
    #    "moon alt": moon_alt,
    #    "monn az": moon_az,
    #    "lunar parllax": lunar_parallax,
    #    "HP": HP,
    #    "SD": SD,
    #    "SD_topo": SD_topo,
    #    "ARCV": ARCV,
    #    "DAZ": DAZ,
    #    #"ARCL": ARCL,
    #    "W_topo": W_topo,
    #    "q_code": q_code,
    #    "q": q
    #}

def run(base_time):
    #result = []
    STEPS = 5
    H = numpy.ndarray(shape=(180 // STEPS, 360 // STEPS))
    for lng in tqdm(range(0, 360, STEPS)):
        for lat in range(0, 180, STEPS):
            q_code = calculate(base_time, 90 - lat, lng - 180) or 'F'
            #result.append(r)
            H[lat // STEPS, lng // STEPS] = ord('F') - ord(q_code)

    #result = pandas.DataFrame.from_records(result)
    #result.to_excel("output_" + str(base_time.Utc()) + ".xlsx")

    plt.imshow(H)
    #plt.savefig(str(base_time.Utc()) + '.png')
    plt.show() # uncomment this when running inside Jupyter Notebook / Google Colab

run(astronomy.Time.Make(2022, 8, 27, 0, 0, 0))
run(astronomy.Time.Make(2022, 8, 28, 0, 0, 0))
run(astronomy.Time.Make(2022, 8, 29, 0, 0, 0))

run(astronomy.Time.Make(2022, 6, 29, 0, 0, 0))
#run(astronomy.Time.Make(2022, 6, 30, 0, 0, 0))
#run(astronomy.Time.Make(2022, 8, 29, 0, 0, 0))
