// MIT, @ebraminio and @hidp123

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>

#include "thirdparty/astronomy.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
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
    astro_time_t sunset_sunrise, moonset_moonrise, best_time, new_moon_prev, new_moon_next;
    double lag_time, moon_age_prev, moon_age_next;
    double sd, lunar_parallax, arcl, arcv, daz, w_topo, sd_topo, value;
    double moon_azimuth, moon_altitude, moon_ra, moon_dec;
    double sun_azimuth, sun_altitude, sun_ra, sun_dec;
};

template<bool evening, bool yallop>
static char calculate(
    double latitude, double longitude, double altitude, astro_time_t base_time,
    /* optional, used for table extra results */ details_t *details = nullptr,
    /* optional, used as an option in table results */ bool ignore_besttime = false,
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
    astro_time_t best_time = (lag_time < 0 || ignore_besttime)
                           ? sunset_sunrise.time
                           : Astronomy_AddDays(sunset_sunrise.time, lag_time * 4 / 9 * (evening ? 1 : -1));
    if (result_time) *result_time = best_time.ut;
    astro_time_t new_moon_prev = Astronomy_SearchMoonPhase(0, sunset_sunrise.time, -35).time;
    astro_time_t new_moon_next = Astronomy_SearchMoonPhase(0, sunset_sunrise.time, +35).time;
    astro_time_t new_moon_nearest = (sunset_sunrise.time.ut - new_moon_prev.ut) <= (new_moon_next.ut - sunset_sunrise.time.ut)
        ? new_moon_prev : new_moon_next;
    if (details) { details->new_moon_prev = new_moon_prev; details->new_moon_next = new_moon_next; }
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
    double ARCV;
    if (yallop) {
        astro_vector_t geomoon = Astronomy_GeoVector(BODY_MOON, best_time, ABERRATION);
        astro_vector_t geosun = Astronomy_GeoVector(BODY_SUN, best_time, ABERRATION);
        astro_rotation_t rot = Astronomy_Rotation_EQJ_EQD(&best_time);
        astro_vector_t rotmoon = Astronomy_RotateVector(rot, geomoon);
        astro_vector_t rotsun  = Astronomy_RotateVector(rot, geosun);
        astro_equatorial_t meq = Astronomy_EquatorFromVector(rotmoon);
        astro_equatorial_t seq = Astronomy_EquatorFromVector(rotsun);
        astro_horizon_t mhor = Astronomy_Horizon(&best_time, observer, meq.ra, meq.dec, REFRACTION_NONE);
        astro_horizon_t shor = Astronomy_Horizon(&best_time, observer, seq.ra, seq.dec, REFRACTION_NONE);
        ARCV = mhor.altitude - shor.altitude;
    } else { // Odeh
        double COSARCV = cos(ARCL * DEG2RAD) / cos(DAZ * DEG2RAD);
        if      (COSARCV < -1) COSARCV = -1;
        else if (COSARCV > +1) COSARCV = +1;
        ARCV = acos(COSARCV) * RAD2DEG;
    }
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

// Lightweight result of the expensive per-location astronomy, used by the map
// renderer. Everything except `value` is expressed in UTC days.
struct geo_point_t {
    bool ok = false;         // Both the sun and moon rise/set events were found.
    double sunset = 0;       // Sun set (evening) or rise (morning) time, UTC days.
    double moonset = 0;      // Moon set (evening) or rise (morning) time, UTC days.
    double value = 0;        // Yallop q value or Odeh V value.
};

// Runs the expensive part of `calculate`: the two rise/set searches and the
// astrometric value computation. It stops short of classifying the result,
// which is cheap and deferred to `classify` so it can be done after bilinear
// interpolation on the coarse grid.
template<bool evening, bool yallop>
static void compute_geo(
    double latitude, double longitude, double altitude, astro_time_t base_time,
    bool ignore_besttime, geo_point_t *out
) {
    astro_time_t time = Astronomy_AddDays(base_time, -longitude / 360);
    astro_observer_t observer = { .latitude = latitude, .longitude = longitude, .height = altitude };

    astro_direction_t direction = evening ? DIRECTION_SET : DIRECTION_RISE;
    astro_search_result_t sunset_sunrise   = Astronomy_SearchRiseSet(BODY_SUN,  observer, direction, time, 1);
    astro_search_result_t moonset_moonrise = Astronomy_SearchRiseSet(BODY_MOON, observer, direction, time, 1);
    out->ok = (sunset_sunrise.status == ASTRO_SUCCESS) && (moonset_moonrise.status == ASTRO_SUCCESS);
    if (!out->ok) return;

    out->sunset  = sunset_sunrise.time.ut;
    out->moonset = moonset_moonrise.time.ut;
    double lag_time = (out->moonset - out->sunset) * (evening ? 1 : -1);
    astro_time_t best_time = (lag_time < 0 || ignore_besttime)
                           ? sunset_sunrise.time
                           : Astronomy_AddDays(sunset_sunrise.time, lag_time * 4 / 9 * (evening ? 1 : -1));

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
    double ARCV;
    if (yallop) {
        astro_vector_t geomoon = Astronomy_GeoVector(BODY_MOON, best_time, ABERRATION);
        astro_vector_t geosun = Astronomy_GeoVector(BODY_SUN, best_time, ABERRATION);
        astro_rotation_t rot = Astronomy_Rotation_EQJ_EQD(&best_time);
        astro_vector_t rotmoon = Astronomy_RotateVector(rot, geomoon);
        astro_vector_t rotsun  = Astronomy_RotateVector(rot, geosun);
        astro_equatorial_t meq = Astronomy_EquatorFromVector(rotmoon);
        astro_equatorial_t seq = Astronomy_EquatorFromVector(rotsun);
        astro_horizon_t mhor = Astronomy_Horizon(&best_time, observer, meq.ra, meq.dec, REFRACTION_NONE);
        astro_horizon_t shor = Astronomy_Horizon(&best_time, observer, seq.ra, seq.dec, REFRACTION_NONE);
        ARCV = mhor.altitude - shor.altitude;
    } else { // Odeh
        double COSARCV = cos(ARCL * DEG2RAD) / cos(DAZ * DEG2RAD);
        if      (COSARCV < -1) COSARCV = -1;
        else if (COSARCV > +1) COSARCV = +1;
        ARCV = acos(COSARCV) * RAD2DEG;
    }
    double W_topo = SD_topo * (1 - cos(ARCL * DEG2RAD)); // In arcminutes

    if (yallop) {
        out->value = (ARCV - (11.8371 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3))) / 10;
    } else { // Odeh
        out->value = ARCV - (7.1651 - 6.3226 * W_topo + .7319 * pow(W_topo, 2) - .1018 * pow(W_topo, 3));
    }
}

// Cheap classification of an (interpolated) point into a visibility code.
template<bool evening, bool yallop>
static char classify(
    const geo_point_t &p, double new_moon_prev, double new_moon_next,
    bool ignore_besttime, bool *draw_moon_line, double *result_time
) {
    if (!p.ok) return 'H'; // No sun{set,rise} or moon{set,rise}

    double lag_time = (p.moonset - p.sunset) * (evening ? 1 : -1);
    double best_time = (lag_time < 0 || ignore_besttime)
                     ? p.sunset
                     : p.sunset + (p.moonset - p.sunset) * 4 / 9;
    if (result_time) *result_time = best_time;

    double new_moon_nearest = (p.sunset - new_moon_prev) <= (new_moon_next - p.sunset)
        ? new_moon_prev : new_moon_next;
    if (draw_moon_line) *draw_moon_line = ((int) round((best_time - new_moon_nearest) * 24 * 20) % 20) == 0;
    bool before_new_moon = (p.sunset - new_moon_nearest) * (evening ? 1 : -1) < 0;

    if (lag_time < 0 && before_new_moon) return 'J'; // Checks both of the conditions on the two below lines, shows a mixed color
    if (lag_time < 0) return 'I'; // Moonset before sunset / Moonrise after sunrise
    if (before_new_moon) return 'G'; // Sunset is before new moon / Sunrise is after new moon

    double value = p.value;
    if (yallop) {
        if      (value >  .216) return 'A'; // Crescent easily visible
        else if (value > -.014) return 'B'; // Crescent visible under perfect conditions
        else if (value > -.160) return 'C'; // May need optical aid to find crescent
        else if (value > -.232) return 'D'; // Will need optical aid to find crescent
        else if (value > -.293) return 'E'; // Crescent not visible with telescope
        else                    return 'F';
    } else { // Odeh
        if      (value >= 5.65) return 'A'; // Crescent is visible by naked eye
        else if (value >= 2.00) return 'C'; // Crescent is visible by optical aid
        else if (value >= -.96) return 'E'; // Crescent is visible only by optical aid
        else                    return 'F';
    }
}

template<bool evening, bool yallop>
static void render(uint32_t *image, astro_time_t base_time) {
    double min_naked_eye_time = INFINITY; unsigned min_naked_eye_x = 0, min_naked_eye_y = 0;
    double min_telescope_time = INFINITY; unsigned min_telescope_x = 0, min_telescope_y = 0;

    const unsigned step = 4; // Coarse grid spacing in pixels: heavy astronomy runs at 1/step of the full resolution.
    const unsigned GW = width / step + 1;
    const unsigned GH = height / step + 1;

    // The new moons are global events, identical for every pixel; find the four
    // that bracket the map's date window once instead of searching per pixel.
    astro_time_t t_center = Astronomy_AddDays(base_time, 0.75);
    astro_time_t nm_prev = Astronomy_SearchMoonPhase(0, t_center, -35).time;
    astro_time_t nm_next = Astronomy_SearchMoonPhase(0, t_center, +35).time;
    astro_time_t nm_prev2 = Astronomy_SearchMoonPhase(0, Astronomy_AddDays(nm_prev, -1), -35).time;
    astro_time_t nm_next2 = Astronomy_SearchMoonPhase(0, Astronomy_AddDays(nm_next, +1), +35).time;
    const double nm_prev_ut = nm_prev.ut, nm_next_ut = nm_next.ut;
    const double nm_prev2_ut = nm_prev2.ut, nm_next2_ut = nm_next2.ut;

    geo_point_t *grid = (geo_point_t *) calloc(GW * GH, sizeof(geo_point_t));

    // Heavy astronomy is only evaluated on the coarse grid.
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned g = 0; g < GW * GH; ++g) {
        unsigned gx = g % GW, gy = g / GW;
        unsigned i = gx * step; if (i > width) i = width;
        unsigned j = gy * step; if (j > height) j = height;
        double latitude  = ((height - (j + 1)) / (double) pixelsPerDegree) + minLatitude;
        double longitude = (i / (double) pixelsPerDegree) + minLongitude;
        compute_geo<evening, yallop>(latitude, longitude, 0, base_time, false, &grid[g]);
    }

    // Cheap per-pixel work: bilinear interpolation over the coarse grid, then classification.
#if defined(_OPENMP)
    #pragma omp parallel for
#endif
    for (unsigned j = 0; j < height; ++j) {
        unsigned gy0 = j / step;
        double fy = (j - gy0 * step) / (double) step;
        const geo_point_t *row0 = &grid[gy0 * GW];
        const geo_point_t *row1 = &grid[(gy0 + 1) * GW];
        for (unsigned i = 0; i < width; ++i) {
            unsigned gx0 = i / step;
            double fx = (i - gx0 * step) / (double) step;
            const geo_point_t &g00 = row0[gx0];
            const geo_point_t &g10 = row0[gx0 + 1];
            const geo_point_t &g01 = row1[gx0];
            const geo_point_t &g11 = row1[gx0 + 1];
            double w00 = (1 - fx) * (1 - fy), w10 = fx * (1 - fy), w01 = (1 - fx) * fy, w11 = fx * fy;

            geo_point_t p;
            double oksum = w00 * g00.ok + w10 * g10.ok + w01 * g01.ok + w11 * g11.ok;
            p.ok = oksum >= 0.5;
            if (oksum > 0) {
                p.sunset  = (w00 * g00.sunset  + w10 * g10.sunset  + w01 * g01.sunset  + w11 * g11.sunset ) / oksum;
                p.moonset = (w00 * g00.moonset + w10 * g10.moonset + w01 * g01.moonset + w11 * g11.moonset) / oksum;
                p.value   = (w00 * g00.value   + w10 * g10.value   + w01 * g01.value   + w11 * g11.value  ) / oksum;
            } else {
                p.sunset = p.moonset = p.value = 0;
            }

            double new_moon_prev, new_moon_next;
            if (p.sunset < nm_prev_ut)      { new_moon_prev = nm_prev2_ut; new_moon_next = nm_prev_ut; }
            else if (p.sunset < nm_next_ut) { new_moon_prev = nm_prev_ut;  new_moon_next = nm_next_ut; }
            else                            { new_moon_prev = nm_next_ut;  new_moon_next = nm_next2_ut; }

            bool draw_moon_line = false;
            double result_time = 0;
            char q_code = classify<evening, yallop>(p, new_moon_prev, new_moon_next, false, &draw_moon_line, &result_time);
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
        }
    }

    free(grid);

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
               "./visibility 2022-08-27 table 34.23,23.3,0 100 > results.tsv\n"
               "./visibility 2022-08-27 table-ignore-besttime 34.23,23.3,0 100 > results.tsv");
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
        
    } else if (!strcmp(argv[2], "table") || !strcmp(argv[2], "table-ignore-besttime")) {
        details_t details;
        bool ignore_besttime = !strcmp(argv[2], "table-ignore-besttime");
        double latitude = atof(strtok((char *) argv[3], ","));
        double longitude = atof(strtok(nullptr, ","));
        double altitude = atof(strtok(nullptr, ","));
        unsigned days = atoi(argv[4]);
        printf("UTC Date\tLatitude\tLongitude\tAltitude\t");

        printf("Sunset\tMoonset%s\tPrev New Moon\tNext New Moon\tMoon age from prev\tMoon age to next\tLag time\t",
               ignore_besttime ? "" : "\tBest time");
        printf("Evening (Yallop)\tq value\t");
        printf("Moon azimuth\tMoon altitude\tMoon ra\tMoon dec\t");
        printf("Sun azimuth\tSun altitude\tSun ra\tSun dec\t");
        printf("Moon sd\tlunar parallax\tarcl geo\tarcv yallop\tdaz\tw topo\tsd topo\t");

        printf("Evening (Odeh)\tV value\t");
        printf("Moon sd\tlunar parallax\tarcl topo\tarcv odeh\tdaz\tw topo\tsd topo\t");

        printf("Sunrise\tMoonrise%s\tPrev New Moon\tNext New Moon\tMoon age from prev\tMoon age to next\tlag time\t",
               ignore_besttime ? "" : "\tBest time");
        printf("Morning (Yallop)\tq value\t");
        printf("Moon azimuth\tMoon altitude\tMoon ra\tMoon dec\t");
        printf("Sun azimuth\tSun altitude\tSun ra\tSun dec\t");
        printf("Moon sd\tlunar parallax\tarcl geo\tarcv yallop\tdaz\tw topo\tsd topo\t");

        printf("Morning (Odeh)\tV value\t");
        printf("Moon sd\tlunar parallax\tarcl topo\tarcv odeh\tdaz\tw topo\tsd topo\t");

        printf("\n");
        for (unsigned i = 0; i < days; ++i) {
            char result;
            astro_utc_t utc = Astronomy_UtcFromTime(time);
            printf("%d-%d-%d\t%f\t%f\t%f\t", utc.year, utc.month, utc.day, latitude, longitude, altitude);
#define LOG(v) printf("%f\t", details.v)
#define TIME(t) utc = Astronomy_UtcFromTime(details.t); printf("%d-%02d-%02d %02d:%02d:%02.2f\t", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
#define TIME_DIFF(t) printf("%s%d:%02d:%02d\t", details.t < 0 ? "-" : "", (int) floor(abs(details.t) * 24), (int) floor(abs(details.t) * 24 * 60 - floor(abs(details.t) * 24) * 60), (int) floor(abs(details.t) * 24 * 60 * 60 - floor(abs(details.t) * 24 * 60) * 60))
            memset(&details, 0, sizeof (details_t));
            result = calculate<true,  true >(latitude, longitude, altitude, time, &details, ignore_besttime);
            TIME(sunset_sunrise); TIME(moonset_moonrise);
            if (!ignore_besttime) { TIME(best_time); }
            TIME(new_moon_prev); TIME(new_moon_next); TIME_DIFF(moon_age_prev); TIME_DIFF(moon_age_next); TIME_DIFF(lag_time);
            printf("%c\t", result); LOG(value);
            LOG(moon_azimuth); LOG(moon_altitude); LOG(moon_ra); LOG(moon_dec);
            LOG(sun_azimuth); LOG(sun_altitude); LOG(sun_ra); LOG(sun_dec);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            printf("%c\t", calculate<true,  false>(latitude, longitude, altitude, time, &details, ignore_besttime)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            result = calculate<false,  true >(latitude, longitude, altitude, time, &details, ignore_besttime);
            TIME(sunset_sunrise); TIME(moonset_moonrise);
            if (!ignore_besttime) { TIME(best_time); }
            TIME(new_moon_prev); TIME(new_moon_next); TIME_DIFF(moon_age_prev); TIME_DIFF(moon_age_next); TIME_DIFF(lag_time);
            printf("%c\t", result); LOG(value);
            LOG(moon_azimuth); LOG(moon_altitude); LOG(moon_ra); LOG(moon_dec);
            LOG(sun_azimuth); LOG(sun_altitude); LOG(sun_ra); LOG(sun_dec);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);

            memset(&details, 0, sizeof (details_t));
            printf("%c\t", calculate<false, false>(latitude, longitude, altitude, time, &details, ignore_besttime)); LOG(value);
            LOG(sd); LOG(lunar_parallax); LOG(arcl); LOG(arcv); LOG(daz); LOG(w_topo); LOG(sd_topo);
#undef TIME
#undef LOG
            printf("\n");
            time = Astronomy_AddDays(time, 1);
        }
        return 0;
    } else {
        printf("Invalid command\n");
        return 1;
    }
}
