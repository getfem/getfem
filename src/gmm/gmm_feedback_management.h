/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2017 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/** @file gmm_feedback_management.h 
    @date July 03, 2017.
    @brief Support for run time management of trace, warning and assert 
           feedback.
*/

#ifndef GMM_FEEDBACK_MANAGEMENT_H__
#define GMM_FEEDBACK_MANAGEMENT_H__

namespace gmm {

/* *********************************************************************** */
/*	GetFEM++ feedback management                  			   */
/* *********************************************************************** */

enum class FeedbackType {
  TRACE = 0,
  WARNING,
  ASSERT
};

// Abstract class providing feedback management interface.
// The interface consist of three functions:
//   * for sending feedback message
//   * for getting traces level
//   * for getting warning level
//   * for action to be done after feedback is handled
struct base_feedback_handler {
  virtual ~base_feedback_handler() = default;
  virtual void send(const std::string &message, FeedbackType messageType, size_t level) = 0;
  virtual size_t traces_level() { return get_traces_level(); }
  virtual size_t warning_level() { return get_warning_level(); }
  virtual void terminating_action() = 0;
};


// Provides the default implementation of feedback handling.
struct default_feedback_handler final : public base_feedback_handler {
  void send(const std::string &message, FeedbackType, size_t) override {
    std::cerr << message << std::endl;
  }
  void terminating_action() override {
    std::exit(1);
  }
};

// This class acts as a run-time dispatcher for sending feedback
// messages and getting trace and warning levels.
class feedback_manager {
public:
    // Steals the pointer to a messenger object that provides
    // feedback handling implementation.
    //
    // Example:
    //   feedback_manager::manage(new default_feedback_handler);
    //
    static base_feedback_handler* manage(base_feedback_handler *pHandler=nullptr);
    static void send(const std::string &message, FeedbackType type, size_t level);
    static size_t traces_level();
    static size_t warning_level();
    // Action to be taken when feedback handling is done
    static void terminating_action();
};

inline base_feedback_handler* feedback_manager::manage(base_feedback_handler *pHandler) {
  static std::unique_ptr<base_feedback_handler> pHandler_ =
    std::move(std::unique_ptr<base_feedback_handler>(new default_feedback_handler));
  if (pHandler != nullptr) {
    pHandler_.reset(pHandler);
  }
  return pHandler_.get();
}

inline void feedback_manager::send(const std::string &message, FeedbackType type, size_t level) {
  feedback_manager::manage()->send(message, type, level);
}

inline void feedback_manager::terminating_action() {
  feedback_manager::manage()->terminating_action();
}

inline size_t feedback_manager::traces_level() {
  return feedback_manager::manage()->traces_level();
}

inline size_t feedback_manager::warning_level() {
  return feedback_manager::manage()->warning_level();
}

} // namespace gmm
#endif /* GMM_FEEDBACK_MANAGEMENT_H__ */
