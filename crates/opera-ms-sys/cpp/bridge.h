// C++ side of the opera-ms-sys bridge: runs the main() of sigma, OPERA-LG
// and short_read_analysis in the calling process (an opera-ms worker
// process).
#pragma once

#include <cstdint>
#include "rust/cxx.h"

namespace opera_ms {

int32_t run_sigma(const rust::Vec<rust::String>& args);
int32_t run_opera_lg(const rust::Vec<rust::String>& args, bool fast);
int32_t run_short_read_analysis(const rust::Vec<rust::String>& args);

}  // namespace opera_ms
