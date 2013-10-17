/***
 *  $Id$
 **
 *  File: MPE_Log.hpp
 *  Created: Oct 15, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_MPE_LOG_HPP
#define MPIX2_MPE_LOG_HPP

#include <mpe.h>

namespace mpix {

  class MPE_Log {
  public:
      explicit MPE_Log(const std::string& name, std::string color) : active_(false) {
	  start(name, color);
      } // MPE_Log

      ~MPE_Log() { stop(); }

      bool start(const std::string& name, std::string color) {
	  if (active_) return false;
	  name_ = name;
	  start_ = MPE_Log_get_event_number();
	  stop_ = MPE_Log_get_event_number();
	  MPE_Describe_state(start_, stop_, name_.c_str(), color.c_str());
	  MPE_Log_event(start_, 0, (name_ + " start").c_str());
	  active_ = true;
	  return true;
      } // start

      void stop() {
	  if (active_) MPE_Log_event(stop_, 0, (name_ + " stop").c_str());
	  active_ = false;
      } // stop

  private:
      MPE_Log(const MPE_Log&);
      void operator=(const MPE_Log&);

      std::string name_;
      bool active_;
      int start_;
      int stop_;

  }; // class MPE_Log

} // namespace mpix

#endif // MPIX2_MPE_LOG_HPP
