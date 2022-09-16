// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cmath>
#include <cstdint>
#include "thirdparty/astro_demo_common.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "thirdparty/stb_image_write.h"

const int pixelsPerDegree = 4;
const int minLatitude = -90;
const int maxLatitude = +90;
const int minLongitude = -180;
const int maxLongitude = +180;
const int width = (maxLongitude - minLongitude) * pixelsPerDegree;
const int height = (maxLatitude - minLatitude) * pixelsPerDegree;

void render(uint32_t *image, astro_time_t base_time);

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
    uint32_t *image = (uint32_t *) calloc(width * height, 4);
    render(image, time);

    return !stbi_write_png(outFileName, width, height, 4, image, width * 4);
}

void render(uint32_t *image, astro_time_t base_time) {
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned i = 0; i < width; ++i) {
        for (unsigned j = 0; j < height; ++j) {
            double longitude = (i / (double)pixelsPerDegree) + minLongitude;
            astro_observer_t observer = {
                .latitude = ((height - (j + 1)) / (double)pixelsPerDegree) + minLatitude,
                .longitude = longitude,
                .height = .0,
            };

            astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
#if 0
            astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
            astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
            if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) continue;
            double lag_time = moonset.time.ut - sunset.time.ut;
            if (lag_time < 0) {
                image[i + j * width] = 0xFF0000FF;
                continue;
            }
            astro_time_t best_time = Astronomy_AddDays(sunset.time, lag_time * 4.0/9);
#else
            astro_search_result_t sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 1);
            astro_search_result_t moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 1);
            if (sunrise.status != ASTRO_SUCCESS || moonrise.status != ASTRO_SUCCESS) continue;
            double lag_time = sunrise.time.ut - moonrise.time.ut;
            if (lag_time < 0) {
                image[i + j * width] = 0xFF0000FF;
                continue;
            }
            astro_time_t best_time = Astronomy_AddDays(sunrise.time, -lag_time * 4.0/9);
#endif
            //
            astro_equatorial_t sun_equator = Astronomy_Equator(BODY_SUN, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
            // double sun_distance = AU_IN_M * Astronomy_VectorLength(sun_equator.vec);
            astro_horizon_t sun_horizon = Astronomy_Horizon(&best_time, observer, sun_equator.ra, sun_equator.dec, REFRACTION_NONE);
            // double sun_alt = sun_horizon.altitude;
            double sun_az = sun_horizon.azimuth;
            //
            astro_equatorial_t moon_equator = Astronomy_Equator(BODY_MOON, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
            astro_libration_t liberation = Astronomy_Libration(best_time);
            // moon_elongation_geo = Astronomy_Elongation(BODY_MOON, best_time); // geocentric elongation

            // double moon_distance = liberation.dist_km;
            astro_horizon_t moon_horizon = Astronomy_Horizon(&best_time, observer, moon_equator.ra, moon_equator.dec, REFRACTION_NONE);
            double moon_alt = moon_horizon.altitude;
            double moon_az = moon_horizon.azimuth;

            // https://github.com/rob-blackbourn/PyFinance/blob/2bbad39b/py_calendrical/location.py#L217

            double SD = liberation.diam_deg * 60 / 2; // semi-diameter of the Moon in arcminutes, geocentric
            double lunar_parallax = SD/0.27245; // in arcminutes
            double SD_topo = SD * (1 + (sin(moon_alt * DEG2RAD) * sin(lunar_parallax/60 * DEG2RAD))); // in arcminutes. Here SD is in arcminutes, moon_alt in degrees, lunar_parallax in degrees (that's why it has been divided by 60).

#define YALLOP 0
#if YALLOP
            astro_elongation_t moon_elongation = Astronomy_Elongation(BODY_MOON, best_time);
            double ARCL = moon_elongation.elongation;
#else // Odeh
            astro_angle_result_t moon_elongation_topo = Astronomy_AngleBetween(sun_equator.vec, moon_equator.vec); // topocentric elongation
            double ARCL = moon_elongation_topo.angle; // in degrees
            //assert(moon_elongation_topo.status == ASTRO_SUCCESS);
#endif
            double DAZ = sun_az - moon_az;
            // double DALT = moon_alt - sun_alt;
            double COSARCV = cos(ARCL * DEG2RAD) / cos(DAZ * DEG2RAD);
            if (COSARCV < -1) COSARCV = -1;
            else if (COSARCV > 1) COSARCV = 1;
            double ARCV = acos(COSARCV) * RAD2DEG;
            double W_topo = SD_topo * (1 - (cos(ARCL * DEG2RAD))); // in arcminutes
#if YALLOP
            double q = (ARCV - (11.8371 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3))) / 10;

            unsigned char q_code;
            if (q > +.216) q_code = 'A'; // Crescent easily visible
            else if (q > -.014) q_code = 'B'; // Crescent visible under perfect conditions
            else if (q > -.160) q_code = 'C'; // May need optical aid to find crescent
            else if (q > -.232) q_code = 'D'; // Will need optical aid to find crescent
            else if (q > -.293) q_code = 'E'; // Crescent not visible with telescope
            else q_code = 'F';
#else // Odeh
            unsigned char q_code;
            double V = ARCV - (7.1651 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3));
            if (V >= 5.65) q_code = 'A'; // Crescent is visible by naked eye
            else if (V >= 2) q_code = 'C'; // Crescent is visible by optical aid
            else if (V >= -.96) q_code = 'E'; // Crescent is visible only by optical aid
            else q_code = 'F';
#endif

            uint32_t color = 0x00000000;
            if (q_code == 'A') color = 0xFF3EFF00;
            else if (q_code == 'B') color = 0xFF3EFF6D;
            else if (q_code == 'C') color = 0xFF00FF9E;
            else if (q_code == 'D') color = 0xFF00FFFA;
            else if (q_code == 'E') color = 0xFF3C78FF;
            else if (q_code == 'F') color = 0x00000000;
            // assert(i >= 0 || i < width);
            // assert(i >= 0 || i < height);
            image[i + j * width] = color;
        }
    }
}
