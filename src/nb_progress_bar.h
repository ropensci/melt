// A custom ProgressBar class to display a progress bar with time estimation.
// Modified to display elapsed time.
//
// Original version available from:
// https://github.com/kforner/rcpp_progress/tree/master/inst/examples/RcppProgressETA/
#ifndef NB_PROGRESS_BAR_HPP_
#define NB_PROGRESS_BAR_HPP_

#include <R_ext/Print.h>
#include <ctime>
#include <stdio.h>
#include <sstream>
#include <string.h>

#include "progress_bar.hpp"

// for unices only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

class NB_ProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====
    /**
     * Main constructor
     */
    NB_ProgressBar()  {
      _max_ticks = 40;
      _finalized = false;
      _timer_flag = true;
    }
    ~NB_ProgressBar() {
    }
    public: // ===== main methods =====
      void display() {
        REprintf("\n");
        flush_console();
      }
      // update display
      void update(float progress) {
        // stop if already finalized
        if (_finalized) return;
        // start time measurement when update() is called the first time
        if (_timer_flag) {
          _timer_flag = false;
          // measure start time
          time(&start);
        } else {
          // measure current time
          time(&end);
          // calculate passed time and remaining time (in seconds)
          double pas_time = std::difftime(end, start);
          double rem_time = (pas_time / progress) * (1 - progress);
          // convert seconds to time string
          std::string rem_time_string = _time_to_string(rem_time);
          std::string pas_time_string = _time_to_string(pas_time);
          // create progress bar string
          std::string progress_bar_string = _current_ticks_display(progress);
          // ensure overwriting of old time info
          int rem_empty_length = rem_time_string.length();
          int pas_empty_length = pas_time_string.length();
          std::string rem_empty_space = std::string(rem_empty_length, ' ');
          std::string pas_empty_space = std::string(pas_empty_length, ' ');
          // merge progress bar and time string
          std::stringstream strs;
          strs << "[" << progress_bar_string << "] " <<
            std::round(100 * progress) <<
              "% | eta: "  << rem_time_string <<
                "| elapsed: " << pas_time_string << pas_empty_space;
          std::string temp_str = strs.str();
          char const* char_type = temp_str.c_str();
          // print: remove old and replace with new
          REprintf("\r");
          REprintf("%s", char_type);
          // finalize display when ready
          if(progress == 1) {
            _finalize_display();
          }
        }
      }
      void end_display() {
        update(1);
      }

      protected: // ==== other instance methods =====

        // convert double with seconds to time string
        std::string _time_to_string(double seconds) {

          int time = (int) seconds;

          int hour = 0;
          int min = 0;
          int sec = 0;

          hour = time / 3600;
          time = time % 3600;
          min = time / 60;
          time = time % 60;
          sec = time;

          std::stringstream time_strs;
          if (hour != 0) time_strs << hour << "h ";
          if (min != 0) time_strs << min << "m ";
          if (sec != 0) time_strs << sec << "s ";
          std::string time_str = time_strs.str();
          return time_str;
        }
        // update the ticks display corresponding to progress
        std::string _current_ticks_display(float progress) {
          int nb_ticks = _compute_nb_ticks(progress);
          std::string cur_display = _construct_ticks_display_string(nb_ticks);
          return cur_display;
        }
        // construct progress bar display
        std::string _construct_ticks_display_string(int nb) {
          std::stringstream ticks_strs;
          for (int i = 0; i < (_max_ticks - 1); ++i) {
            if (i < nb) {
              ticks_strs << "#";
            } else {
              ticks_strs << "-";
            }
          }
          std::string tick_space_string = ticks_strs.str();
          return tick_space_string;
        }
        // finalize
        void _finalize_display() {
          if (_finalized) return;
          // REprintf("\n");
          flush_console();
          _finalized = true;
        }
        // compute number of ticks according to progress
        int _compute_nb_ticks(float progress) {
          return int(progress * _max_ticks);
        }
        // N.B: does nothing on windows
        void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
#endif
        }
        private: // ===== INSTANCE VARIABLES ====
          int _max_ticks;   		// the total number of ticks to print
          bool _finalized;
          bool _timer_flag;
          time_t start,end;

};
#endif
