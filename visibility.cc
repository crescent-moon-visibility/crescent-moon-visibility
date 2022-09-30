// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#include "thirdparty/astronomy.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "thirdparty/stb_image_write.h"

const unsigned pixelsPerDegree = PIXEL_PER_DEGREE;
const int minLatitude = -90;
const int maxLatitude = +90;
const int minLongitude = -180;
const int maxLongitude = +180;
const unsigned width = (maxLongitude - minLongitude) * pixelsPerDegree;
const unsigned height = (maxLatitude - minLatitude) * pixelsPerDegree;

struct details_t {
    astro_time_t sun_rise, moon_rise, best_time;
    double lag_time, moon_age_prev, moon_age_next;
    double sd, lunar_parallax, arcl, arcv, daz, w_topo, sd_topo, value;
    double moon_azimuth, moon_altitude, moon_ra, moon_dec;
    double sun_azimuth, sun_altitude, sun_ra, sun_dec;
};

template<bool evening, bool yallop>
static char calculate(double latitude, double longitude, double altitude, astro_time_t base_time, details_t *details) {
    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
    astro_observer_t observer = { .latitude = latitude, .longitude = longitude, .height = altitude };
    astro_time_t best_time;
    if (evening) {
        astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
        astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
        if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) return 'H'; // No sunset or moonset
        double lag_time = moonset.time.ut - sunset.time.ut;
        if (details) { details->lag_time = lag_time; details->moon_rise = moonset.time; details->sun_rise = sunset.time; }
        if (lag_time < 0) return 'I'; // Moonset before sunset
        best_time = Astronomy_AddDays(sunset.time, lag_time * 4 / 9);

        astro_time_t new_moon_prev = Astronomy_SearchMoonPhase(0, time, -35).time;
        astro_time_t new_moon_next = Astronomy_SearchMoonPhase(0, time, +35).time;
        astro_time_t new_moon_nearest = (sunset.time.ut - new_moon_prev.ut) <= (new_moon_next.ut - sunset.time.ut)
            ? new_moon_prev : new_moon_next;
        if (details) {
            details->moon_age_prev = best_time.ut - new_moon_prev.ut;
            details->moon_age_next = best_time.ut - new_moon_next.ut;
        }
        if (sunset.time.ut < new_moon_nearest.ut) return 'G'; // sunset is before new moon
    } else {
        astro_search_result_t sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 1);
        astro_search_result_t moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 1);
        if (sunrise.status != ASTRO_SUCCESS || moonrise.status != ASTRO_SUCCESS) return 'H'; // No sunrise or moonrise
        double lag_time = sunrise.time.ut - moonrise.time.ut;
        if (details) { details->lag_time = lag_time; details->moon_rise = moonrise.time; details->sun_rise = sunrise.time; }
        if (lag_time < 0) return 'I'; // Moonrise after sunrise
        best_time = Astronomy_AddDays(sunrise.time, -lag_time * 4 / 9);

        astro_time_t new_moon_prev = Astronomy_SearchMoonPhase(0, time, -35).time;
        astro_time_t new_moon_next = Astronomy_SearchMoonPhase(0, time, +35).time;
        astro_time_t new_moon_nearest = (sunrise.time.ut - new_moon_prev.ut) <= (new_moon_next.ut - sunrise.time.ut)
            ? new_moon_prev : new_moon_next;
        if (details) {
            details->moon_age_prev = best_time.ut - new_moon_prev.ut;
            details->moon_age_next = best_time.ut - new_moon_next.ut;
        }
        if (sunrise.time.ut > new_moon_nearest.ut) return 'G'; // sunrise is after new moon
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
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned i = 0; i < width; ++i) {
        for (unsigned j = 0; j < height; ++j) {
            double latitude = ((height - (j + 1)) / (double) pixelsPerDegree) + minLatitude;
            double longitude = (i / (double) pixelsPerDegree) + minLongitude;
            char q_code = calculate<evening, yallop>(latitude, longitude, 0, base_time, nullptr);
            uint32_t color = 0x00000000;
            if      (q_code == 'A') color = 0xFF3EFF00; // These color codes are in 0xAAGGBBRR format
            else if (q_code == 'B') color = 0xFF3EFF6D;
            else if (q_code == 'C') color = 0xFF00FF9E;
            else if (q_code == 'D') color = 0xFF00FFFA;
            else if (q_code == 'E') color = 0xFF3C78FF;
            else if (q_code == 'F') color = 0x00000000;
            else if (q_code == 'G') color = 0xFF9988CC;
            else if (q_code == 'H') color = 0x00000000;
            else if (q_code == 'I') color = 0xFF0000FF;
            image[i + j * width] = color;
        }
    }
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

        printf("Sunset\tMoonset\tBest time\tMoon age prev\tMoon age next\tlag time\t");
        printf("Evening/Yallop\tq value\t");
        printf("Moon azimuth\tMoon altitude\tMoon ra\tMoon dec\t");
        printf("Sun azimuth\tSun altitude\tSun ra\tSun dec\t");
        printf("Moon sd\tlunar parallax\tarcl geo\tarcv yallop\tdaz\tw topo\tsd topo\t");

        printf("Evening/Odeh\tV value\t");
        printf("Moon sd\tlunar parallax\tarcl topo\tarcv odeh\tdaz\tw topo\tsd topo\t");

        printf("Sunrise\tMoonrise\tBest time\tMoon age prev\tMoon age next\tlag time\t");
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
            TIME(sun_rise); TIME(moon_rise); TIME(best_time); LOG(moon_age_prev); LOG(moon_age_next); LOG(lag_time);
            printf("%c\t", result); LOG(value);
            LOG(moon_azimuth); LOG(moon_altitude); LOG(moon_ra); LOG(moon_dec);
            LOG(sun_azimuth); LOG(sun_altitude); LOG(sun_ra); LOG(sun_dec);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            printf("%c\t", calculate<true,  false>(latitude, longitude, altitude, time, &details)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            result = calculate<false,  true >(latitude, longitude, altitude, time, &details);
            TIME(sun_rise); TIME(moon_rise); TIME(best_time); LOG(moon_age_prev); LOG(moon_age_next); LOG(lag_time);
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
