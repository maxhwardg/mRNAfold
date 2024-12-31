#include "logging.hpp"

#include <iostream>

namespace mwmrna {
bool allow_logging = true;
void log(const std::string &msg, const std::experimental::source_location &loc) {
    if (!allow_logging) return;
    std::clog << loc.file_name() << ":" << loc.line() << ":" << loc.function_name() << " " << msg
              << std::endl;
}
}  // namespace mwmrna