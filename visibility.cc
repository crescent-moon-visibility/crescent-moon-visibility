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

void render(uint32_t *image, bool waning, bool odeh, astro_time_t base_time);
char calculate(double latitude, double longitude, bool waning, bool odeh, astro_time_t base_time);

int main(int argc, const char **argv) {
    if (argc == 1) {
        printf("Run this like,\n./visibility.out 2022-08-27T00:00:00Z waning yallop map out.png\n./visibility.out 2022-08-27T00:00:00Z waxing odeh calculate 34.23 23.3");
        return 1;
    }

    astro_time_t time;
    if (ParseTime(argv[1], &time)) return 1;

    bool waning;
    if (strcmp(argv[2], "waning") == 0) waning = true;
    else if (strcmp(argv[2], "waxing") == 0) waning = false;
    else return 1;

    bool odeh;
    if (strcmp(argv[3], "odeh") == 0) odeh = true;
    else if (strcmp(argv[3], "yallop") == 0) odeh = false;
    else return 1;

    if (strcmp(argv[4], "map") == 0) {
        uint32_t *image = (uint32_t *) calloc(width * height, 4);
        render(image, waning, odeh, time);
        return !stbi_write_png(argv[5], width, height, 4, image, width * 4);
    } else if (strcmp(argv[4], "calculate") == 0) {
        printf("%c", calculate(atof(argv[5]), atof(argv[6]), waning, odeh, time));
        return 0;
    } else return 1;
}

char calculate(double latitude, double longitude, bool waning, bool odeh, astro_time_t base_time) {
    astro_observer_t observer = {
        .latitude = latitude,
        .longitude = longitude,
        .height = .0,
    };

    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);

    astro_time_t best_time;
    if (!waning) {
        astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
        astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
        if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) return 'F';
        double lag_time = moonset.time.ut - sunset.time.ut;
        if (lag_time < 0) return 'G';
        best_time = Astronomy_AddDays(sunset.time, lag_time * 4.0/9);
    } else {
        astro_search_result_t sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 1);
        astro_search_result_t moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 1);
        if (sunrise.status != ASTRO_SUCCESS || moonrise.status != ASTRO_SUCCESS) return 'F';
        double lag_time = sunrise.time.ut - moonrise.time.ut;
        if (lag_time < 0) return 'G';
        best_time = Astronomy_AddDays(sunrise.time, -lag_time * 4.0/9);
    }

    astro_equatorial_t sun_equator = Astronomy_Equator(BODY_SUN, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
    astro_horizon_t sun_horizon = Astronomy_Horizon(&best_time, observer, sun_equator.ra, sun_equator.dec, REFRACTION_NONE);
    double sun_az = sun_horizon.azimuth;
    astro_equatorial_t moon_equator = Astronomy_Equator(BODY_MOON, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
    astro_libration_t liberation = Astronomy_Libration(best_time);

    astro_horizon_t moon_horizon = Astronomy_Horizon(&best_time, observer, moon_equator.ra, moon_equator.dec, REFRACTION_NONE);
    double moon_alt = moon_horizon.altitude;
    double moon_az = moon_horizon.azimuth;

    // https://github.com/rob-blackbourn/PyFinance/blob/2bbad39b/py_calendrical/location.py#L217

    double SD = liberation.diam_deg * 60 / 2; // semi-diameter of the Moon in arcminutes, geocentric
    double lunar_parallax = SD/0.27245; // in arcminutes
    double SD_topo = SD * (1 + (sin(moon_alt * DEG2RAD) * sin(lunar_parallax/60 * DEG2RAD))); // in arcminutes. Here SD is in arcminutes, moon_alt in degrees, lunar_parallax in degrees (that's why it has been divided by 60).

    double ARCL;
    if (!odeh) { // Yallop
        astro_elongation_t moon_elongation = Astronomy_Elongation(BODY_MOON, best_time);
        ARCL = moon_elongation.elongation;
    } else {
        astro_angle_result_t moon_elongation_topo = Astronomy_AngleBetween(sun_equator.vec, moon_equator.vec); // topocentric elongation
        ARCL = moon_elongation_topo.angle; // in degrees
        //assert(moon_elongation_topo.status == ASTRO_SUCCESS);
    }

    double DAZ = sun_az - moon_az;
    // double DALT = moon_alt - sun_alt;
    double COSARCV = cos(ARCL * DEG2RAD) / cos(DAZ * DEG2RAD);
    if (COSARCV < -1) COSARCV = -1;
    else if (COSARCV > 1) COSARCV = 1;
    double ARCV = acos(COSARCV) * RAD2DEG;
    double W_topo = SD_topo * (1 - (cos(ARCL * DEG2RAD))); // in arcminutes
    
    if (!odeh) { // Yallop
        double q = (ARCV - (11.8371 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3))) / 10;
        if (q > +.216) return 'A'; // Crescent easily visible
        else if (q > -.014) return 'B'; // Crescent visible under perfect conditions
        else if (q > -.160) return 'C'; // May need optical aid to find crescent
        else if (q > -.232) return 'D'; // Will need optical aid to find crescent
        else if (q > -.293) return 'E'; // Crescent not visible with telescope
        else return 'F';
    } else { // Odeh
        double V = ARCV - (7.1651 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3));
        if (V >= 5.65) return 'A'; // Crescent is visible by naked eye
        else if (V >= 2) return 'C'; // Crescent is visible by optical aid
        else if (V >= -.96) return 'E'; // Crescent is visible only by optical aid
        else return 'F';
    }
}

void render(uint32_t *image, bool waning, bool odeh, astro_time_t base_time) {
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned i = 0; i < width; ++i) {
        for (unsigned j = 0; j < height; ++j) {
            double latitude = ((height - (j + 1)) / (double)pixelsPerDegree) + minLatitude;
            double longitude = (i / (double)pixelsPerDegree) + minLongitude;
            char q_code = calculate(latitude, longitude, waning, odeh, base_time);
            uint32_t color = 0x00000000;
            if      (q_code == 'A') color = 0xFF3EFF00;
            else if (q_code == 'B') color = 0xFF3EFF6D;
            else if (q_code == 'C') color = 0xFF00FF9E;
            else if (q_code == 'D') color = 0xFF00FFFA;
            else if (q_code == 'E') color = 0xFF3C78FF;
            else if (q_code == 'F') color = 0x00000000;
            else if (q_code == 'G') color = 0xFF0000FF;
            image[i + j * width] = color;
        }
    }
}
