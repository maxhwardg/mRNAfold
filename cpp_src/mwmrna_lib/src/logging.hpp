#ifndef MWMRNA_LIB_LOGGING_HPP
#define MWMRNA_LIB_LOGGING_HPP

#include <experimental/source_location>
#include <string>

namespace mwmrna {
extern bool allow_logging;
void log(const std::string &msg, const std::experimental::source_location &loc =
                                     std::experimental::source_location::current());
}  // namespace mwmrna

#endif  // MWMRNA_LIB_LOGGING_HPP