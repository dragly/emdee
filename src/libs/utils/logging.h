#ifndef LOGGING_H
#define LOGGING_H

#ifdef MD_USE_GLOG
#include <glog/logging.h>
#else
#include <glogfallback.h>
#endif

#endif // LOGGING_H
