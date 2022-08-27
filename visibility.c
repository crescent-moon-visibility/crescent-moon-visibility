// NOTE: It is not yet accurate as we expect
// MIT, @ebraminio and @hidp123

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "thirdparty/astro_demo_common.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "thirdparty/stb_image_write.h"

const int PixelsPerDegree = 2;
const int MinLatitude = -90;
const int MaxLatitude = +90;
const int MinLongitude = -180;
const int MaxLongitude = +180;
const int PixelsWide = (MaxLongitude - MinLongitude) * PixelsPerDegree;
const int PixelsHigh = (MaxLatitude - MinLatitude) * PixelsPerDegree;

void render(uint32_t *image, unsigned width, unsigned height, astro_time_t base_time);

int main(int argc, const char **argv) {
    const char *outFileName;
    astro_time_t time;

    if (argc == 1) {
        printf("Run this like, ./visibility.out out.png [2022-08-27T00:00:00Z]");
        return 1;
    } else if (argc == 2) {
        outFileName = argv[1];
        time = Astronomy_CurrentTime();
    } else {
        outFileName = argv[1];
        if (ParseTime(argv[2], &time)) return 1;
    }

    // Create a world map image in memory for the given time.
    uint32_t image[PixelsWide * PixelsHigh] = {0};
    render(image, PixelsWide, PixelsHigh, time);

    return !stbi_write_png(outFileName, PixelsWide, PixelsHigh, 4, (uint8_t *) image, PixelsWide * 4);
}

void render(uint32_t *image, unsigned width, unsigned height, astro_time_t base_time) {
    astro_observer_t observer;
    observer.height = 0.0;
    const double AU_IN_M = 149597871000.0;

    for (unsigned i = 0; i < PixelsWide; ++i) {
        observer.longitude = (i / (double)PixelsPerDegree) + MinLongitude;
        for (unsigned j = 0; j < PixelsHigh; ++j) {
            observer.latitude = ((PixelsHigh - (j + 1)) / (double)PixelsPerDegree) + MinLatitude;

            astro_time_t time = Astronomy_AddDays(base_time, -observer.longitude / 360.0);
            astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
            astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
            if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) continue;
            // https://astro.ukho.gov.uk/moonwatch/background.html
            // lag time: The time interval between sunset and moonset. The lag time is usually
            // given in minutes. It can be negative, indicating that the Moon sets before the Sun.
            double lag_time = moonset.time.ut - sunset.time.ut;
            // best time: an empirical prediction of the time which gives the observer the best opportunity
            // to see the new crescent Moon (Sunset time + (4/9)*Lag time).
            astro_time_t best_time = Astronomy_TimeFromDays(sunset.time.ut + lag_time * 4.0/9);
            //
            astro_equatorial_t sun_equator = Astronomy_Equator(BODY_SUN, &best_time, observer, EQUATOR_J2000, ABERRATION);
            // double sun_distance = AU_IN_M * Astronomy_VectorLength(sun_equator.vec);
            astro_horizon_t sun_horizon = Astronomy_Horizon(&best_time, observer, sun_equator.ra, sun_equator.dec, REFRACTION_NORMAL);
            double sun_alt = sun_horizon.altitude;
            double sun_az = sun_horizon.azimuth;
            //
            astro_equatorial_t moon_equator = Astronomy_Equator(BODY_MOON, &best_time, observer, EQUATOR_J2000, ABERRATION);
            double moon_distance = AU_IN_M * Astronomy_VectorLength(moon_equator.vec);
            astro_horizon_t moon_horizon = Astronomy_Horizon(&best_time, observer, moon_equator.ra, moon_equator.dec, REFRACTION_NORMAL);
            double moon_alt = moon_horizon.altitude;
            double moon_az = moon_horizon.azimuth;

            // https://github.com/rob-blackbourn/PyFinance/blob/2bbad39b/py_calendrical/location.py#L217
            double lunar_parallax = 6378140.0 / moon_distance * cos(moon_alt * DEG2RAD);

            // https://github.com/abdullah-alhashim/prayer_calculator/blob/8abe558/moon_sighting.py#L54-L62
            double HP = lunar_parallax / cos(moon_alt * DEG2RAD);
            double SD = 0.27245 * HP * (180.0 * 60 / M_PI); // semi-diameter of the Moon
            double SD_topo = SD * (1 + (sin(moon_alt * DEG2RAD) * sin(HP)));

            // https://github.com/abdullah-alhashim/prayer_calculator/blob/8abe558/moon_sighting.py#L71-L77
            double ARCV = moon_alt - sun_alt;
            double DAZ = sun_az - moon_az;
            // double ARCL = acos(cos(ARCV * DEG2RAD) * cos(DAZ * DEG2RAD)) * RAD2DEG;
            double W_topo = SD_topo * (1 - (cos(ARCV * DEG2RAD) * cos(DAZ * DEG2RAD)));
            double q = (ARCV - (11.8371 - 6.3226*W_topo + 0.7319*pow(W_topo, 2) - 0.1018*pow(W_topo, 3))) / 10;

            unsigned char q_code = 'F';
            if (q > +0.216) q_code = 'A'; // Crescent easily visible
            else if (+0.216 >= q && q > -0.014) q_code = 'B'; // Crescent visible under perfect conditions
            else if (-0.014 >= q && q > -0.160) q_code = 'C'; // May need optical aid to find crescent
            else if (-0.160 >= q && q > -0.232) q_code = 'D'; // Will need optical aid to find crescent
            else if (-0.232 >= q && q > -0.293) q_code = 'E'; // Crescent not visible with telescope
            else if (-0.293 >= q) q_code = 'F';

            uint8_t value = ('F' - q_code) * 255 / 5;
            uint32_t *pixel = &image[i + j*width];
            assert(i >= 0 || i < width);
            assert(i >= 0 || i < height);
            *pixel = value + (value << 8) + (value << 16) + 0xFF000000;
        }
    }
}
