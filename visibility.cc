// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cstdlib>
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

template<bool evening, bool yallop>
static void render(uint32_t *image, astro_time_t base_time);
template<bool evening, bool yallop>
static char calculate(double latitude, double longitude, astro_time_t base_time);

int main(int argc, const char **argv) {
    if (argc == 1) {
        printf("Run this like,\n"
               "./visibility.out 2022-08-27T00:00:00Z evening yallop map out.png\n"
               "./visibility.out 2022-08-27T00:00:00Z morning odeh calculate 34.23 23.3");
        return 1;
    }

    astro_time_t time;
    if (ParseTime(argv[1], &time)) return 1;

    bool evening;
    if      (strcmp(argv[2], "evening") == 0) evening = true;
    else if (strcmp(argv[2], "morning") == 0) evening = false;
    else return 1;

    bool yallop;
    if      (strcmp(argv[3], "yallop") == 0) yallop = true;
    else if (strcmp(argv[3], "odeh")   == 0) yallop = false;
    else return 1;

    if (strcmp(argv[4], "map") == 0) {
        uint32_t *image = (uint32_t *) calloc(width * height, 4);
        evening
            ? (yallop ? render<true,  true>(image, time) : render<true,  false>(image, time))
            : (yallop ? render<false, true>(image, time) : render<false, false>(image, time));
        return !stbi_write_png(argv[5], width, height, 4, image, width * 4);
    } else if (strcmp(argv[4], "calculate") == 0) {
        char c = evening
            ? (yallop ? calculate<true,  true>(atof(argv[5]), atof(argv[6]), time) : calculate<true,  false>(atof(argv[5]), atof(argv[6]), time))
            : (yallop ? calculate<false, true>(atof(argv[5]), atof(argv[6]), time) : calculate<false, false>(atof(argv[5]), atof(argv[6]), time));
        printf("%c\n", c);
        return 0;
    } else return 1;
}

template<bool evening, bool yallop>
static char calculate(double latitude, double longitude, astro_time_t base_time) {
    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
    astro_observer_t observer = { .latitude = latitude, .longitude = longitude, .height = .0 };
    astro_time_t best_time;
    if (evening) {
        astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
        astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
        if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) return 'F';
        double lag_time = moonset.time.ut - sunset.time.ut;
        if (lag_time < 0) return 'G'; // Moonset before sunset
        best_time = Astronomy_AddDays(sunset.time, lag_time * 4 / 9);
    } else {
        astro_search_result_t sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 1);
        astro_search_result_t moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 1);
        if (sunrise.status != ASTRO_SUCCESS || moonrise.status != ASTRO_SUCCESS) return 'F';
        double lag_time = sunrise.time.ut - moonrise.time.ut;
        if (lag_time < 0) return 'G'; // Moonrise after sunrise
        best_time = Astronomy_AddDays(sunrise.time, -lag_time * 4 / 9);
    }

    astro_equatorial_t sun_equator = Astronomy_Equator(BODY_SUN, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
    astro_horizon_t sun_horizon = Astronomy_Horizon(&best_time, observer, sun_equator.ra, sun_equator.dec, REFRACTION_NONE);
    astro_equatorial_t moon_equator = Astronomy_Equator(BODY_MOON, &best_time, observer, EQUATOR_OF_DATE, ABERRATION);
    astro_horizon_t moon_horizon = Astronomy_Horizon(&best_time, observer, moon_equator.ra, moon_equator.dec, REFRACTION_NONE);
    astro_libration_t liberation = Astronomy_Libration(best_time);

    double SD = liberation.diam_deg * 60 / 2; // Semi-diameter of the Moon in arcminutes, geocentric
    double lunar_parallax = SD / 0.27245; // In arcminutes
    // As SD_topo should be in arcminutes as SD, but moon_alt and lunar_parallax are in degrees, it is divided by 60.
    double SD_topo = SD * (1 + sin(moon_horizon.altitude * DEG2RAD) * sin(lunar_parallax / 60 * DEG2RAD));

    double ARCL = yallop
        ? Astronomy_Elongation(BODY_MOON, best_time).elongation // Geocentric elongation in Yallop
        : Astronomy_AngleBetween(sun_equator.vec, moon_equator.vec).angle; // Topocentric elongation in Odeh

    double DAZ = sun_horizon.azimuth - moon_horizon.azimuth;
    double COSARCV = cos(ARCL * DEG2RAD) / cos(DAZ * DEG2RAD);
    if      (COSARCV < -1) COSARCV = -1;
    else if (COSARCV > +1) COSARCV = +1;
    double ARCV = acos(COSARCV) * RAD2DEG;
    double W_topo = SD_topo * (1 - cos(ARCL * DEG2RAD)); // in arcminutes

    if (yallop) {
        double q = (ARCV - (11.8371 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3))) / 10;
        if      (q > +.216) return 'A'; // Crescent easily visible
        else if (q > -.014) return 'B'; // Crescent visible under perfect conditions
        else if (q > -.160) return 'C'; // May need optical aid to find crescent
        else if (q > -.232) return 'D'; // Will need optical aid to find crescent
        else if (q > -.293) return 'E'; // Crescent not visible with telescope
        else return 'F';
    } else { // Odeh
        double V = ARCV - (7.1651 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3));
        if      (V >= 5.65) return 'A'; // Crescent is visible by naked eye
        else if (V >= 2.00) return 'C'; // Crescent is visible by optical aid
        else if (V >= -.96) return 'E'; // Crescent is visible only by optical aid
        else return 'F';
    }
}

template<bool evening, bool yallop>
static void render(uint32_t *image, astro_time_t base_time) {
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned i = 0; i < width; ++i) {
        for (unsigned j = 0; j < height; ++j) {
            double latitude = ((height - (j + 1)) / (double) pixelsPerDegree) + minLatitude;
            double longitude = (i / (double) pixelsPerDegree) + minLongitude;
            char q_code = calculate<evening, yallop>(latitude, longitude, base_time);
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
