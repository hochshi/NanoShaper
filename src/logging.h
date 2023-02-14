#ifndef logging_h
#define logging_h

#include <string>
#include <utility>

#ifdef SPDLOG_ENABLED
#include <spdlog/spdlog.h>

// template <typename... Args>
// using format_string_t = spdlog::format_string_t<Args...>;

// #else
// template <typename... Args> using format_string_t = std::string;
// template <typename... Args> using format_string_t = std::string_view;
#endif  // SPDLOG_ENABLED

namespace nanoshaper {
namespace logging {

#ifdef SPDLOG_ENABLED
using level = spdlog::level::level_enum;
#else
enum level : int {
  trace = 0,
  debug = 1,
  info = 2,
  warn = 3,
  err = 4,
  critical = 5,
  off = 6,
  n_levels
};
#endif

template <level lvl, typename... Args>
inline void log(std::string fmt, Args&&... args) {
#ifdef SPDLOG_ENABLED
  spdlog::log(lvl, fmt, std::forward<Args>(args)...);
#endif
}

}  // namespace logging
}  // namespace nanoshaper
#endif  // logging_h
