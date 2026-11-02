#include "opera-ms-sys/cpp/bridge.h"

#include <string>
#include <vector>

#include "Configure.h"

// The programs' main() functions, renamed by build.rs (-Dmain=...).
int sigma_main(int argc, char** argv);
int opera_lg_main(int argc, char* argv[]);
int short_read_analysis_main(int argc, char* argv[]);

namespace opera_ms {

namespace {

// Calls a main() with argv = {program, args...}.
int32_t call_main(int (*main_function)(int, char**), const char* program,
                  const rust::Vec<rust::String>& args) {
    std::vector<std::string> owned;
    owned.emplace_back(program);
    for (const auto& arg : args) {
        owned.emplace_back(std::string(arg));
    }
    std::vector<char*> argv;
    for (auto& arg : owned) {
        argv.push_back(&arg[0]);
    }
    argv.push_back(nullptr);
    return main_function(static_cast<int>(owned.size()), argv.data());
}

}  // namespace

int32_t run_sigma(const rust::Vec<rust::String>& args) {
    return call_main(sigma_main, "sigma", args);
}

int32_t run_opera_lg(const rust::Vec<rust::String>& args, bool fast) {
    if (fast) {
        // OPERA-LG-fast: the values the Makefile builds in with -D
        // (OPERA-LG/src/Makefile, FAST_DEFINES); no config-file key sets them.
        Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 100000;
        Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 5;
    }
    return call_main(opera_lg_main, fast ? "OPERA-LG-fast" : "OPERA-LG", args);
}

int32_t run_short_read_analysis(const rust::Vec<rust::String>& args) {
    return call_main(short_read_analysis_main, "short_read_analysis", args);
}

}  // namespace opera_ms
