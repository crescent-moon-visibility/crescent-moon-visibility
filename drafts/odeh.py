import matplotlib.pyplot as plt
import math
import numpy
import astronomy
#import pandas
from tqdm import tqdm


KM_PER_AU = 1.4959787069098932e+8   #<const> The number of kilometers per astronomical unit.

def calculate(base_time, latitude, longitude):
    observer = astronomy.Observer(latitude, longitude)
    time = base_time.AddDays(-observer.longitude / 360) # this corrects the base time based on timezone
    sunset   = astronomy.SearchRiseSet(astronomy.Body.Sun,  observer, astronomy.Direction.Set, time, 1)
    moonset  = astronomy.SearchRiseSet(astronomy.Body.Moon, observer, astronomy.Direction.Set, time, 1)
    if sunset is None or moonset is None: return {}
    #print(latitude, longitude)

    # https://astro.ukho.gov.uk/moonwatch/background.html
    # lag time: The time interval between sunset and moonset. The lag time is usually
    # given in minutes. It can be negative, indicating that the Moon sets before the Sun.
    lag_time = moonset.ut - sunset.ut

    # best time: an empirical prediction of the time which gives the observer the best opportunity
    # to see the new crescent Moon (Sunset time + (4/9)*Lag time).
    best_time = astronomy.Time(sunset.ut + lag_time * 4/9)

    sun_equator = astronomy.Equator(astronomy.Body.Sun, best_time, observer, True, True)
    #sun_distance = KM_PER_AU * sun_equator.vec.Length()
    sun_horizon = astronomy.Horizon(best_time, observer, sun_equator.ra, sun_equator.dec, astronomy.Refraction.Airless)
    sun_alt = sun_horizon.altitude
    sun_az = sun_horizon.azimuth

    moon_equator = astronomy.Equator(astronomy.Body.Moon, best_time, observer, True, True) #RA is in h.dd (hours.degrees)
    moon_elongation_geo = astronomy.Elongation(astronomy.Body.Moon, best_time) #geocentric elongation
    moon_elongation_topo = astronomy.AngleBetween(sun_equator.vec, moon_equator.vec) #topocentric elongation
    #moon_distance = KM_PER_AU * moon_equator.vec.Length()
    moon_horizon = astronomy.Horizon(best_time, observer, moon_equator.ra, moon_equator.dec, astronomy.Refraction.Airless)
    moon_alt = moon_horizon.altitude
    moon_az = moon_horizon.azimuth
    
    #print(moon_elongation_topo)

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
    ARCL = moon_elongation_topo #in degrees
    DAZ = sun_az - moon_az
    DALT = moon_alt - sun_alt
    COSARCV = math.cos(math.radians(ARCL))/math.cos(math.radians(DAZ))
    if -1 <= COSARCV <= 1: ARCV = math.degrees(math.acos(COSARCV)) #math.degrees(math.acos(math.cos(math.radians(ARCL))/math.cos(math.radians(DAZ)))) #moon_alt - sun_alt
    elif COSARCV < -1: ARCV = math.degrees(math.acos(-1))
    elif COSARCV > 1: ARCV = math.degrees(math.acos(1)) 
    #print(ARCV)

    W_topo = SD_topo * (1 - (math.cos(math.radians(ARCL)))) #in arcminutes
    #q = (ARCV - (11.8371 - 6.3226*W_topo + 0.7319*W_topo**2 - 0.1018*W_topo**3)) / 10
    V = ARCV - (7.1651 - 6.3226 * W_topo + 0.7319 * math.pow(W_topo, 2) - 0.1018 * math.pow(W_topo, 3))

    if V >= 5.65: q_code = 'A' # Crescent is visible by naked eye
    elif +5.65 > V >= 2: q_code = 'B' # Crescent is visible by optical aid
    elif +2 > V >= -0.96: q_code = 'C' # Crescent is visible only by optical aid
    elif -0.96 > V: q_code = 'D'

    #if q > +0.216: q_code = 'A' # Crescent easily visible
    #elif +0.216 >= q > -0.014: q_code = 'B' # Crescent visible under perfect conditions
    #elif -0.014 >= q > -0.160: q_code = 'C' # May need optical aid to find crescent
    #elif -0.160 >= q > -0.232: q_code = 'D' # Will need optical aid to find crescent
    #elif -0.232 >= q > -0.293: q_code = 'E' # Crescent not visible with telescope
    #elif -0.293 >= q: q_code = 'F'

    return {
        "lat": observer.latitude,
        "long": observer.longitude,
        "sunset": sunset.Utc(),
        "moonset": moonset.Utc(),
        "lag time": lag_time,
        "best time": best_time.Utc(),
        #"sun dist": sun_distance,
        "sun alt": sun_alt,
        "sun az": sun_az,
        "moon elong geo": moon_elongation_geo.elongation,
        "moon elong topo": moon_elongation_topo,
        #"moon dist": moon_distance,
        "moon alt": moon_alt,
        "monn az": moon_az,
        "lunar parllax": lunar_parallax,
        #"HP": HP,
        "SD": SD,
        "SD_topo": SD_topo,
        "COSARCV": COSARCV,
        "ARCV": ARCV,
        "DALT": DALT,
        "DAZ": DAZ,
        "ARCL": ARCL,
        "W_topo": W_topo,
        "q_code": q_code,
        "V": V
    }

def run(base_time):
    result = []
    STEPS = 3
    H = numpy.ndarray(shape=(180 // STEPS, 360 // STEPS))
    for lng in tqdm(range(0, 360, STEPS)):
        for lat in range(0, 180, STEPS):
            r = calculate(base_time, 90 - lat, lng - 180)
            if 'q_code' in r:
                result.append(r)
                H[lat // STEPS, lng // STEPS] = ord('D') - ord(r['q_code'])
            else:
                H[lat // STEPS, lng // STEPS] = 0

    #result = pandas.DataFrame.from_records(result)
    #result.to_excel("output_" + str(base_time.Utc()).replace(":", "-") + " (" + str(STEPS) + " steps)" + ".xlsx") # + str(base_time.Utc()) + 

    plt.imshow(H)
    plt.show()

run(astronomy.Time.Make(2022, 6, 29, 0, 0, 0))
#run(astronomy.Time.Make(2022, 6, 30, 0, 0, 0))
#run(astronomy.Time.Make(2022, 8, 29, 0, 0, 0))

# result = [calculate(d, c[1], c[0]) 
#           for d in [
#               astronomy.Time.Make(2022, 6, 29, 0, 0, 0),
#               astronomy.Time.Make(2022, 8, 28, 0, 0, 0),
#               astronomy.Time.Make(2022, 8, 29, 0, 0, 0),
#           ] 
#           for c in [
#               [52, -1.5],
#               [32, -8],
#               [-29, 26],
#           ]
#          ]
# result = pandas.DataFrame.from_records(result)
# result.to_excel("output_selection.xlsx")
