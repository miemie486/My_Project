#include "msg.h"

int verbose_level_fdv {_VBSLVL_HIGH_};

void show_message_fdv(std::string message, int vbs_level) {
    if (! (verbose_level_fdv < vbs_level)) {
    std::cout << message << std::endl << std::flush;
  }
}

// void show_message_timestamp_fdv(std::string &message, int vbs_level) {
void show_message_timestamp_fdv(std::string message, int vbs_level) {

  if (! (verbose_level_fdv < vbs_level)) {
    show_message_fdv(message, vbs_level);
    time_t now = time(0);
    std::cout << _TIMER_COL_ << ctime(&now) << std::flush;
  }
}

void set_verbose_level_fdv(int vbs_level) {
  verbose_level_fdv = vbs_level;
}
