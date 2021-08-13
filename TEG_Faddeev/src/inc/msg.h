#ifndef MSGINC
#define MSGINC
#include <iostream>
#include <string>

#define _VBSLVL_DEBUG_ 20
#define _VBSLVL_HIGH_ 10
#define _VBSLVL_NORMAL_ 5
#define _VBSLVL_LOW_ 1
#define _VBSLVL_URGENT_ -100
#define _TIMER_COL_ "                                    ["

extern int verbose_level_fdv;

// show message if verbose_level_fdv >= vbs_level
// void show_message_fdv(std::string &message, int vbs_level);
void show_message_fdv(std::string message, int vbs_level);
// show message and time stamp if verbose_level_fdv >= vbs_level
// void show_message_timestamp_fdv(std::string &message, int vbs_level);
void show_message_timestamp_fdv(std::string message, int vbs_level);
void set_verbose_level_fdv(int vbs_level);

#endif
