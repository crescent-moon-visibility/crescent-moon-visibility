#ifndef __ASTRONOMY_DEMO_COMMON_H
#define __ASTRONOMY_DEMO_COMMON_H
#include "astronomy.h"

#ifdef __cplusplus
extern "C" {
#endif
int ParseArgs(int argc, const char *argv[], astro_observer_t *observer, astro_time_t *time);
int ParseTime(const char *text, astro_time_t *time);
void PrintTime(astro_time_t time);
#ifdef __cplusplus
}
#endif

#endif /* __ASTRONOMY_DEMO_COMMON_H */
