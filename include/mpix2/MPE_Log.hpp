/***
 *  $Id$
 **
 *  File: MPE_Log.hpp
 *  Created: Oct 15, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_MPE_LOG_HPP
#define MPIX2_MPE_LOG_HPP

#include <mpe.h>

namespace mpix {

  class MPE_Log {
  public:
      MPE_Log() : active_(false), running_(false) { }

      MPE_Log(const std::string& name, const std::string& color)
	  : active_(false), running_(false) {
	  init(name, color);
      } // MPE_Log

      ~MPE_Log() { stop(); }

      bool init(const std::string& name, const std::string& color) {
	  if (running_) return false;
	  start_ = MPE_Log_get_event_number();
	  stop_ = MPE_Log_get_event_number();
	  MPE_Describe_state(start_, stop_, name.c_str(), color.c_str());
	  active_ = true;
      } // init

      bool start() {
	  if ((!active_) || (running_)) return false;
	  MPE_Log_event(start_, 0, 0);
	  running_ = true;
	  return true;
      } // start

      bool stop() {
	  if (!running_) return false;
	  MPE_Log_event(stop_, 0, 0);
	  running_ = false;
	  return true;
      } // stop

  private:
      MPE_Log(const MPE_Log&);
      void operator=(const MPE_Log&);

      bool active_;
      bool running_;

      int start_;
      int stop_;

  }; // class MPE_Log

} // namespace mpix

#endif // MPIX2_MPE_LOG_HPP
