// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#include "thirdparty/astronomy.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wno-deprecated-declarations"
#include "thirdparty/stb_image_write.h"
#pragma GCC diagnostic pop

// To be passed compiled time
#ifndef PIXEL_PER_DEGREE
#define PIXEL_PER_DEGREE 4
#endif

const unsigned pixelsPerDegree = PIXEL_PER_DEGREE;
const int minLatitude = -90;
const int maxLatitude = +90;
const int minLongitude = -180;
const int maxLongitude = +180;
const unsigned width = (maxLongitude - minLongitude) * pixelsPerDegree;
const unsigned height = (maxLatitude - minLatitude) * pixelsPerDegree;

struct details_t {
    astro_time_t sunset_sunrise, moonset_moonrise, best_time;
    double lag_time, moon_age_prev, moon_age_next;
    double sd, lunar_parallax, arcl, arcv, daz, w_topo, sd_topo, value;
    double moon_azimuth, moon_altitude, moon_ra, moon_dec;
    double sun_azimuth, sun_altitude, sun_ra, sun_dec;
};

template<bool evening, bool yallop>
static char calculate(
    double latitude, double longitude, double altitude, astro_time_t base_time,
    /* optional, used for table extra results */ details_t *details = nullptr,
    /* optional, used for moon ages lines in map */ bool *draw_moon_line = nullptr,
    /* optional, used for first visibility points in map */ double *result_time = nullptr,
    /* optional, used for red line in map */ double *q_value = nullptr
) {
    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
    astro_observer_t observer = { .latitude = latitude, .longitude = longitude, .height = altitude };

    astro_direction_t direction = evening ? DIRECTION_SET : DIRECTION_RISE;
    astro_search_result_t sunset_sunrise   = Astronomy_SearchRiseSet(BODY_SUN,  observer, direction, time, 1);
    astro_search_result_t moonset_moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, direction, time, 1);
    if (sunset_sunrise.status != ASTRO_SUCCESS || moonset_moonrise.status != ASTRO_SUCCESS) return 'H'; // No sun{set,rise} or moon{set,rise}
    double lag_time = (moonset_moonrise.time.ut - sunset_sunrise.time.ut) * (evening ? 1 : -1);
    if (details) { details->lag_time = lag_time; details->moonset_moonrise = moonset_moonrise.time; details->sunset_sunrise = sunset_sunrise.time; }
    astro_time_t best_time = lag_time < 0 ? sunset_sunrise.time : Astronomy_AddDays(sunset_sunrise.time, lag_time * 4 / 9 * (evening ? 1 : -1));
    if (result_time) *result_time = best_time.ut;
    astro_time_t new_moon_prev = Astronomy_SearchMoonPhase(0, sunset_sunrise.time, -35).time;
    astro_time_t new_moon_next = Astronomy_SearchMoonPhase(0, sunset_sunrise.time, +35).time;
    astro_time_t new_moon_nearest = (sunset_sunrise.time.ut - new_moon_prev.ut) <= (new_moon_next.ut - sunset_sunrise.time.ut)
        ? new_moon_prev : new_moon_next;
    if (draw_moon_line) *draw_moon_line = ((int) round((best_time.ut - new_moon_nearest.ut) * 24 * 20) % 20) == 0;
    if (details) {
        details->moon_age_prev = best_time.ut - new_moon_prev.ut;
        details->moon_age_next = best_time.ut - new_moon_next.ut;
    }
    bool before_new_moon = (sunset_sunrise.time.ut - new_moon_nearest.ut) * (evening ? 1 : -1) < 0;
    if (lag_time < 0 && before_new_moon) return 'J'; // Checks both of the conditions on the two below lines, shows a mixed color
    if (lag_time < 0) return 'I'; // Moonset before sunset / Moonrise after sunrise
    if (before_new_moon) return 'G'; // Sunset is before new moon / Sunrise is after new moon

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
    double W_topo = SD_topo * (1 - cos(ARCL * DEG2RAD)); // In arcminutes

    char result = ' ';
    double value;
    if (yallop) {
        value = (ARCV - (11.8371 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3))) / 10;
        if      (value > +.216) result = 'A'; // Crescent easily visible
        else if (value > -.014) result = 'B'; // Crescent visible under perfect conditions
        else if (value > -.160) result = 'C'; // May need optical aid to find crescent
        else if (value > -.232) result = 'D'; // Will need optical aid to find crescent
        else if (value > -.293) result = 'E'; // Crescent not visible with telescope
        else result = 'F';
    } else { // Odeh
        value = ARCV - (7.1651 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3));
        if      (value >= 5.65) result = 'A'; // Crescent is visible by naked eye
        else if (value >= 2.00) result = 'C'; // Crescent is visible by optical aid
        else if (value >= -.96) result = 'E'; // Crescent is visible only by optical aid
        else result = 'F';
    }
    if (q_value) *q_value = value;

    if (details) {
        details->best_time = best_time;
        details->sd = SD; details->lunar_parallax = lunar_parallax; details->arcl = ARCL; details->arcv = ARCV;
        details->daz = DAZ; details->w_topo = W_topo; details->sd_topo = SD_topo; details->value = value;
        details->moon_azimuth = moon_horizon.azimuth, details->moon_altitude = moon_horizon.altitude;
        details->moon_ra = moon_horizon.ra; details->moon_dec = moon_horizon.dec;
        details->sun_azimuth = sun_horizon.azimuth; details->sun_altitude = sun_horizon.altitude;
        details->sun_ra = sun_horizon.ra; details->sun_dec = sun_horizon.dec;
    }

    return result;
}

template<bool evening, bool yallop>
static void render(uint32_t *image, astro_time_t base_time) {
    double min_naked_eye_time = INFINITY; unsigned min_naked_eye_x = 0, min_naked_eye_y = 0;
    double min_telescope_time = INFINITY; unsigned min_telescope_x = 0, min_telescope_y = 0;

#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned i = 0; i < width; ++i) {
        // double max_q_value = -INFINITY; unsigned max_q_value_x = 0, max_q_value_y = 0;
        for (unsigned j = 0; j < height; ++j) {
            double latitude = ((height - (j + 1)) / (double) pixelsPerDegree) + minLatitude;
            double longitude = (i / (double) pixelsPerDegree) + minLongitude;
            bool draw_moon_line = false;
            double result_time = 0;
            double q_value = -INFINITY;
            char q_code = calculate<evening, yallop>(latitude, longitude, 0, base_time, nullptr, &draw_moon_line, &result_time, &q_value);
            uint32_t color = 0x00000000;
            if      (q_code == 'A') color = 0xFF3EFF00; // These color codes are in 0xAAGGBBRR format
            else if (q_code == 'B') color = 0xFF3EFF6D;
            else if (q_code == 'C') color = 0xFF00FF9E;
            else if (q_code == 'D') color = 0xFF00FFFA;
            else if (q_code == 'E') color = 0xFF3C78FF;
            else if (q_code == 'F') color = 0x00000000;
            else if (q_code == 'G') color = 0xFFAD0D6A;
            else if (q_code == 'H') color = 0x00000000;
            else if (q_code == 'I') color = 0xFF0000FF;
            else if (q_code == 'J') color = 0xFF5707B5;
            if (draw_moon_line) color = 0xFFB0B0B0;
            image[i + j * width] = color;

            if ((q_code == 'A' || q_code == 'B') && result_time < min_naked_eye_time)
#if defined(_OPENMP)
                #pragma omp critical
#endif
            { min_naked_eye_x = i; min_naked_eye_y = j; min_naked_eye_time = result_time; }
            if ((q_code == 'C' || q_code == 'D') && result_time < min_telescope_time)
#if defined(_OPENMP)
                #pragma omp critical
#endif
            { min_telescope_x = i; min_telescope_y = j; min_telescope_time = result_time; }

//             if (q_value > max_q_value)
// #if defined(_OPENMP)
//                 #pragma omp critical
// #endif
//             { max_q_value_x = i; max_q_value_y = j; max_q_value = q_value; }
        }

        // if (max_q_value_x != 0 && max_q_value_y != 0)
        //     image[max_q_value_x + max_q_value_y * width] = 0xFF0000FF;
    }

    #define DIAMOND_SIZE 7
    if (min_naked_eye_x != 0 && min_naked_eye_y != 0) {
        for (int i = -DIAMOND_SIZE; i <= DIAMOND_SIZE; ++i) {
            for (int j = -DIAMOND_SIZE; j <= DIAMOND_SIZE; ++j) {
                if (abs(i) + abs(j) > DIAMOND_SIZE) continue;
                unsigned naked_eye = min_naked_eye_x + i + (min_naked_eye_y + j) * width;
                if (naked_eye < width * height) image[naked_eye] = 0xFF0000FF;
            }
        }
    }
    if (min_telescope_x != 0 && min_telescope_y != 0) {
        for (int i = -DIAMOND_SIZE; i <= DIAMOND_SIZE; ++i) {
            for (int j = -DIAMOND_SIZE; j <= DIAMOND_SIZE; ++j) {
                if (abs(i) + abs(j) > DIAMOND_SIZE) continue;
                unsigned telescope = min_telescope_x + i + (min_telescope_y + j) * width;
                if (telescope < width * height) image[telescope] = 0xFF0000FF;
            }
        }
    }
    #undef DIAMOND_SIZE
}

int main(int argc, const char **argv) {
    if (argc == 1) {
        printf("Run like this:\n"
               "./visibility 2022-08-27 map evening yallop out.png\n"
               "./visibility 2022-08-27 table 34.23,23.3,0 100 > results.tsv");
        return 1;
    }

    int year = atoi(strtok((char *) argv[1], "-"));
    int month = atoi(strtok(nullptr, "-"));
    int day = atoi(strtok(nullptr, "-"));
    astro_time_t time = Astronomy_MakeTime(year, month, day, 0, 0, 0);

    if (!strcmp(argv[2], "map")) {
        bool evening;
        if      (!strcmp(argv[3], "evening")) evening = true;
        else if (!strcmp(argv[3], "morning")) evening = false;
        else return 1;

        bool yallop;
        if      (!strcmp(argv[4], "yallop")) yallop = true;
        else if (!strcmp(argv[4], "odeh"  )) yallop = false;
        else return 1;

        uint32_t *image = (uint32_t *) calloc(width * height, 4);
        evening
            ? (yallop ? render<true,  true>(image, time) : render<true,  false>(image, time))
            : (yallop ? render<false, true>(image, time) : render<false, false>(image, time));
        return !stbi_write_png(argv[5], width, height, 4, image, width * 4);
    } else if (!strcmp(argv[2], "table")) {
        details_t details;
        double latitude = atof(strtok((char *) argv[3], ","));
        double longitude = atof(strtok(nullptr, ","));
        double altitude = atof(strtok(nullptr, ","));
        unsigned days = atoi(argv[4]);
        printf("UTC Date\tLatitude\tLongitude\tAltitude\t");

        printf("Sunset\tMoonset\tBest time\tMoon age from prev\tMoon age from next\tlag time\t");
        printf("Evening/Yallop\tq value\t");
        printf("Moon azimuth\tMoon altitude\tMoon ra\tMoon dec\t");
        printf("Sun azimuth\tSun altitude\tSun ra\tSun dec\t");
        printf("Moon sd\tlunar parallax\tarcl geo\tarcv yallop\tdaz\tw topo\tsd topo\t");

        printf("Evening/Odeh\tV value\t");
        printf("Moon sd\tlunar parallax\tarcl topo\tarcv odeh\tdaz\tw topo\tsd topo\t");

        printf("Sunrise\tMoonrise\tBest time\tMoon age from prev\tMoon age from next\tlag time\t");
        printf("Morning/Yallop\tq value\t");
        printf("Moon azimuth\tMoon altitude\tMoon ra\tMoon dec\t");
        printf("Sun azimuth\tSun altitude\tSun ra\tSun dec\t");
        printf("Moon sd\tlunar parallax\tarcl geo\tarcv yallop\tdaz\tw topo\tsd topo\t");

        printf("Morning/Odeh\tV value\t");
        printf("Moon sd\tlunar parallax\tarcl topo\tarcv odeh\tdaz\tw topo\tsd topo\t");

        printf("\n");
        for (unsigned i = 0; i < days; ++i) {
            char result;
            astro_utc_t utc = Astronomy_UtcFromTime(time);
            printf("%d-%d-%d\t%f\t%f\t%f\t", utc.year, utc.month, utc.day, latitude, longitude, altitude);
#define LOG(v) printf("%f\t", details.v)
#define TIME(t) utc = Astronomy_UtcFromTime(details.t); printf("%d-%02d-%02d %02d:%02d:%02.2f\t", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
            memset(&details, 0, sizeof (details_t));
            result = calculate<true,  true >(latitude, longitude, altitude, time, &details);
            TIME(sunset_sunrise); TIME(moonset_moonrise); TIME(best_time); LOG(moon_age_prev); LOG(moon_age_next); LOG(lag_time);
            printf("%c\t", result); LOG(value);
            LOG(moon_azimuth); LOG(moon_altitude); LOG(moon_ra); LOG(moon_dec);
            LOG(sun_azimuth); LOG(sun_altitude); LOG(sun_ra); LOG(sun_dec);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            printf("%c\t", calculate<true,  false>(latitude, longitude, altitude, time, &details)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            result = calculate<false,  true >(latitude, longitude, altitude, time, &details);
            TIME(sunset_sunrise); TIME(moonset_moonrise); TIME(best_time); LOG(moon_age_prev); LOG(moon_age_next); LOG(lag_time);
            printf("%c\t", result); LOG(value);
            LOG(moon_azimuth); LOG(moon_altitude); LOG(moon_ra); LOG(moon_dec);
            LOG(sun_azimuth); LOG(sun_altitude); LOG(sun_ra); LOG(sun_dec);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            printf("%c\t", calculate<false, false>(latitude, longitude, altitude, time, &details)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);
#undef TIME
#undef LOG
            printf("\n");
            time = Astronomy_AddDays(time, 1);
        }
        return 0;
    } else return 1;
}
