// ConfigFile.cpp

#include <ConfigFile.h>
#include <main_functions.h>
#include <istream>
#include <stdexcept>
#include <string>

#ifdef JSON_ENABLED
#include <nlohmann/json.hpp>
#endif

#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW

using std::string;

namespace nanoshaper {

ConfigFile::ConfigFile(string filename, string delimiter, string comment,
                       string sentry, string format)
    : myDelimiter(delimiter), myComment(comment), mySentry(sentry) {
  // Construct a ConfigFile, getting keys and values from given file

  std::ifstream in(filename.c_str());

  if (!in) {
    std::string err_msg = "File not found ";
    err_msg += filename;
    throw std::invalid_argument(filename.c_str());
  }

  if ("plain" == format) {
    fileParser(in, *this);
  } else if ("json" == format) {
    jsonParser(in, *this);
  } else {
    std::string err_msg = "Invalid format ";
    err_msg += format;
    throw std::invalid_argument(err_msg.c_str());
  }
}

ConfigFile::ConfigFile()
    : myDelimiter(string(1, '=')), myComment(string(1, '#')) {
  // Construct a ConfigFile without a file; empty
}

void ConfigFile::remove(const string& key) {
  // Remove key and its value
  myContents.erase(myContents.find(key));
  return;
}

bool ConfigFile::keyExists(const string& key) const {
  // Indicate whether key is found
  mapci p = myContents.find(key);
  return (p != myContents.end());
}

/* static */
void ConfigFile::trim(string& s) {
  // Remove leading and trailing whitespace
  static const char whitespace[] = " \n\t\v\r\f";
  s.erase(0, s.find_first_not_of(whitespace));
  s.erase(s.find_last_not_of(whitespace) + 1U);
}

std::ostream& operator<<(std::ostream& os, const ConfigFile& cf) {
  // Save a ConfigFile to os
  for (ConfigFile::mapci p = cf.myContents.begin(); p != cf.myContents.end();
       ++p) {
    os << p->first << " " << cf.myDelimiter << " ";
    os << p->second << std::endl;
  }
  return os;
}

std::istream& jsonParser(std::istream& is, ConfigFile& cf) {
#ifdef JSON_ENABLED
  nlohmann::json data = nlohmann::json::parse(is);
  auto contents = data.get<std::map<std::string, std::string>>();
  cf.setContents(contents);
  return is;
#else
  throw std::logic_error(
      "JSON parsing is not available. To enable recompile with "
      "ENABLE_JSON=ON.");
#endif
}

// std::istream& operator>>( std::istream& is, ConfigFile& cf )
std::istream& fileParser(std::istream& is, ConfigFile& cf) {
  // Load a ConfigFile from is
  // Read in keys and values, keeping internal whitespace
  typedef string::size_type pos;
  const string& delim = cf.myDelimiter;  // separator
  const string& comm = cf.myComment;     // comment
  const string& sentry = cf.mySentry;    // end of file sentry
  const pos skip = delim.length();       // length of separator

  string nextline = "";  // might need to read ahead to see where value ends
  std::map<std::string, std::string> contents;

  while (is || nextline.length() > 0) {
    // Read an entire line at a time
    string line;
    if (nextline.length() > 0) {
      line = nextline;  // we read ahead; use it now
      nextline = "";
    } else {
      std::getline(is, line);
    }

    // Ignore comments
    line = line.substr(0, line.find(comm));

    // Check for end of file sentry
    if (sentry != "" && line.find(sentry) != string::npos)
      return is;

    // Parse the line if it contains a delimiter
    pos delimPos = line.find(delim);
    if (delimPos < string::npos) {
      // Extract the key
      string key = line.substr(0, delimPos);
      line.replace(0, delimPos + skip, "");

      // See if value continues on the next line
      // Stop at blank line, next line with a key, end of stream,
      // or end of file sentry
      bool terminate = false;
      while (!terminate && is) {
        std::getline(is, nextline);
        terminate = true;

        string nlcopy = nextline;
        ConfigFile::trim(nlcopy);
        if (nlcopy == "")
          continue;

        nextline = nextline.substr(0, nextline.find(comm));
        if (nextline.find(delim) != string::npos)
          continue;
        if (sentry != "" && nextline.find(sentry) != string::npos)
          continue;

        nlcopy = nextline;
        ConfigFile::trim(nlcopy);
        if (nlcopy != "")
          line += "\n";
        line += nextline;
        terminate = false;
      }

      // Store key and value
      ConfigFile::trim(key);
      ConfigFile::trim(line);
      contents[key] = line;  // overwrites if key is repeated
    }
  }

  cf.setContents(contents);
  return is;
}
}  // namespace nanoshaper
