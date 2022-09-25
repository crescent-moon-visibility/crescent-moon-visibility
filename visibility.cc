// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#include "thirdparty/astronomy.h"
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

struct details_t {
    astro_time_t sun_rise, moon_rise;
    double lag_time;
    double sd, lunar_parallax, arcl, arcv, daz, w_topo, sd_topo, value, cosarcv;
    double moon_horizon_azimuth, moon_horizon_altitude, moon_horizon_ra, moon_horizon_dec;
    double sun_horizon_azimuth, sun_horizon_altitude, sun_horizon_ra, sun_horizon_dec;
};
template<bool evening, bool yallop>
static char calculate(double latitude, double longitude, double altitude, astro_time_t base_time, details_t *details);

int main(int argc, const char **argv) {
    if (argc == 1) {
        printf("Run this like,\n"
               "./visibility 2022-08-27 map evening yallop out.png\n"
               "./visibility 2022-08-27 csv 34.23,23.3,0 100");
        return 1;
    }

    int year = atoi(strtok((char *) argv[1], "-"));
    int month = atoi(strtok(nullptr, "-"));
    int day = atoi(strtok(nullptr, "-"));
    astro_time_t time = Astronomy_MakeTime(year, month, day, 0, 0, 0);

    if (strcmp(argv[2], "map") == 0) {
        bool evening;
        if      (strcmp(argv[3], "evening") == 0) evening = true;
        else if (strcmp(argv[3], "morning") == 0) evening = false;
        else return 1;

        bool yallop;
        if      (strcmp(argv[4], "yallop") == 0) yallop = true;
        else if (strcmp(argv[4], "odeh")   == 0) yallop = false;
        else return 1;

        uint32_t *image = (uint32_t *) calloc(width * height, 4);
        evening
            ? (yallop ? render<true,  true>(image, time) : render<true,  false>(image, time))
            : (yallop ? render<false, true>(image, time) : render<false, false>(image, time));
        return !stbi_write_png(argv[5], width, height, 4, image, width * 4);
    } else if (strcmp(argv[2], "csv") == 0) {
        details_t details;
        double latitude = atof(strtok((char *) argv[3], ","));
        double longitude = atof(strtok(nullptr, ","));
        double altitude = atof(strtok(nullptr, ","));
        unsigned days = atoi(argv[4]);
        printf("UTC Date,Latitude,Longitude,Altitude,");
        printf("Evening/Yallop,q,");
        printf("Sunset,Moonset,lag time,");
        printf("sd,lunar parallax,arcl,arcv,daz,w topo,sd topo,cosarcv,");
        printf("moon horizon azimuth,moon horizon altitude,moon horizon ra,moon horizon dec,");
        printf("sun horizon azimuth,sun horizon altitude,sun horizon ra,sun horizon dec,");

        printf("Evening/Odeh,V,");
        printf("sd,lunar parallax,arcl,arcv,daz,w topo,sd topo,cosarcv,");

        printf("Morning/Yallop,q,");
        printf("Sunrise,Moonrise,lag time,");
        printf("sd,lunar parallax,arcl,arcv,daz,w topo,sd topo,cosarcv,");
        printf("moon horizon azimuth,moon horizon altitude,moon horizon ra,moon horizon dec,");
        printf("sun horizon azimuth,sun horizon altitude,sun horizon ra,sun horizon dec,");

        printf("Morning/Odeh,V,");
        printf("sd,lunar parallax,arcl,arcv,daz,w topo,sd topo,cosarcv");

        printf("\n");
        for (unsigned i = 0; i < days; ++i) {
            astro_utc_t utc = Astronomy_UtcFromTime(time);
            printf("%d-%d-%d,%f,%f,%f,", utc.year, utc.month, utc.day, latitude, longitude, altitude);
#define LOG(v) printf("%f,", details.v)
#define TIME(t) utc = Astronomy_UtcFromTime(details.t); printf("%d-%d-%d-%d-%d-%f,", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
            memset(&details, 0, sizeof (details_t));
            printf("%c,", calculate<true,  true >(latitude, longitude, altitude, time, &details)); LOG(value);
            TIME(sun_rise); TIME(moon_rise); LOG(lag_time);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo); LOG(cosarcv);
            LOG(moon_horizon_azimuth); LOG(moon_horizon_altitude); LOG(moon_horizon_ra); LOG(moon_horizon_dec);
            LOG(sun_horizon_azimuth); LOG(sun_horizon_altitude); LOG(sun_horizon_ra); LOG(sun_horizon_dec);

            memset(&details, 0, sizeof (details_t));
            printf("%c,", calculate<true,  false>(latitude, longitude, altitude, time, &details)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo); LOG(cosarcv);

            memset(&details, 0, sizeof (details_t));
            printf("%c,", calculate<false, true >(latitude, longitude, altitude, time, &details)); LOG(value);
            TIME(sun_rise); TIME(moon_rise); LOG(lag_time);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo); LOG(cosarcv);
            LOG(moon_horizon_azimuth); LOG(moon_horizon_altitude); LOG(moon_horizon_ra); LOG(moon_horizon_dec);
            LOG(sun_horizon_azimuth); LOG(sun_horizon_altitude); LOG(sun_horizon_ra); LOG(sun_horizon_dec);

            memset(&details, 0, sizeof (details_t));
            printf("%c,", calculate<false, false>(latitude, longitude, altitude, time, &details)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo); LOG(cosarcv);
#undef TIME
#undef LOG
            printf("\n");
            time = Astronomy_AddDays(time, 1);
        }
        return 0;
    } else return 1;
}

template<bool evening, bool yallop>
static char calculate(double latitude, double longitude, double altitude, astro_time_t base_time, details_t *details) {
    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
    astro_observer_t observer = { .latitude = latitude, .longitude = longitude, .height = altitude };
    astro_time_t best_time;
    if (evening) {
        astro_search_result_t sunset  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_SET, time, 1);
        astro_search_result_t moonset = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_SET, time, 1);
        if (sunset.status != ASTRO_SUCCESS || moonset.status != ASTRO_SUCCESS) return 'G'; // No sunset or moonset
        double lag_time = moonset.time.ut - sunset.time.ut;
        if (details) { details->lag_time = lag_time; details->moon_rise = time; details->sun_rise = time; }
        if (lag_time < 0) return 'H'; // Moonset before sunset
        best_time = Astronomy_AddDays(sunset.time, lag_time * 4 / 9);
    } else {
        astro_search_result_t sunrise  = Astronomy_SearchRiseSet(BODY_SUN,  observer, DIRECTION_RISE, time, 1);
        astro_search_result_t moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, DIRECTION_RISE, time, 1);
        if (sunrise.status != ASTRO_SUCCESS || moonrise.status != ASTRO_SUCCESS) return 'G'; // No sunrise or moonrise
        double lag_time = sunrise.time.ut - moonrise.time.ut;
        if (details) { details->lag_time = lag_time; details->moon_rise = time; details->sun_rise = time; }
        if (lag_time < 0) return 'H'; // Moonrise after sunrise
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
        details->sd = SD; details->lunar_parallax = lunar_parallax; details->arcl = ARCL; details->arcv = ARCV;
        details->daz = DAZ; details->w_topo = W_topo; details->sd_topo = SD_topo; details->value = value;
        details->cosarcv = COSARCV;
        details->moon_horizon_azimuth = moon_horizon.azimuth, details->moon_horizon_altitude = moon_horizon.altitude;
        details->moon_horizon_ra = moon_horizon.ra; details->moon_horizon_dec = moon_horizon.dec;
        details->sun_horizon_azimuth = sun_horizon.azimuth; details->sun_horizon_altitude = sun_horizon.altitude;
        details->sun_horizon_ra = sun_horizon.ra; details->sun_horizon_dec = sun_horizon.dec;
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
            if      (q_code == 'A') color = 0xFF3EFF00;
            else if (q_code == 'B') color = 0xFF3EFF6D;
            else if (q_code == 'C') color = 0xFF00FF9E;
            else if (q_code == 'D') color = 0xFF00FFFA;
            else if (q_code == 'E') color = 0xFF3C78FF;
            else if (q_code == 'F') color = 0x00000000;
            else if (q_code == 'G') color = 0x00000000;
            else if (q_code == 'H') color = 0xFF0000FF;
            image[i + j * width] = color;
        }
    }
}
