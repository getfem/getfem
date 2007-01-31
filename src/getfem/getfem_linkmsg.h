// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**@file linkmsg.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date March 21, 1999.
   @brief structure for managing Messages between objects.
*/

#ifndef GETFEM_LINKMSG_H__
#define GETFEM_LINKMSG_H__

#include "dal_basic.h"
#include "dal_tas.h"
#include "getfem_config.h"
#include <queue>

namespace getfem {

  /* ******************************************************************** */
  /*	MASK TYPE.                                		          */
  /* ******************************************************************** */

  class mask {

    protected :

      gmm::uint32_type tab; /* uint 64 does'nt exist every where.         */

    public :

      bool operator[] (int i) const { return ((tab & (1 << i)) != 0); }
      const mask &operator &= (const mask &m) { tab &= m.tab; return *this; }
      const mask &operator |= (const mask &m) { tab |= m.tab; return *this; }
      mask operator &(const mask &m) const { mask r(*this); r &= m; return r; }
      mask operator |(const mask &m) const { mask r(*this); r |= m; return r; }
      mask(void) { tab = 0xFFFFFFFF; }
      mask(int i) { tab = (1 << i); }
  };


  /* ******************************************************************** */
  /*	SENDER CLASS.                                		          */
  /* ******************************************************************** */

  class virtual_linkmsg_receiver;

  class virtual_linkmsg_sender {
  public : 
    virtual void sup_receiver(virtual_linkmsg_receiver *) = 0;
    virtual ~virtual_linkmsg_sender() {}
  };


  template<class RECEIVER>
  class linkmsg_sender : public virtual_linkmsg_sender {
  public :
    
    typedef int msg_id_type;
    typedef size_t size_type;
    
  protected :
    
    typedef dal::dynamic_tas<RECEIVER *, 3> RECEIVERTAB;
    RECEIVERTAB receivers;
    dal::dynamic_array<mask, 3> masks;
    
  public :
    
    void add_receiver(RECEIVER &re, mask m = mask())
    { masks[receivers.add(&re)] = m; }
    void sup_receiver(virtual_linkmsg_receiver *);
    template<class T> void send(const T &) const;
    linkmsg_sender(void) {};
    virtual ~linkmsg_sender();
    const linkmsg_sender& operator =(const linkmsg_sender &) {
      DAL_THROW(internal_error, "The copy of this object doesn't work");
      return *this;
    }
    linkmsg_sender(const linkmsg_sender &v) { *this = v; }
  };

  template<class RECEIVER>
    void linkmsg_sender<RECEIVER>::sup_receiver(virtual_linkmsg_receiver *p) {
    typename RECEIVERTAB::tas_iterator it = receivers.tas_begin(),
      ite = receivers.tas_end();
    // cerr << "sup_receiver " << p << " " << receivers.card() << endl;
    for (; it != ite; ++it)
      if (*it == p) {
	receivers.sup(it.index());
	// cerr << " sup " << it.index() << endl;
      }
  }

  template<class RECEIVER> template<class T> 
    void linkmsg_sender<RECEIVER>::send(const T &msg) const {
    typename RECEIVERTAB::const_tas_iterator it = receivers.tas_begin(),
      ite = receivers.tas_end();
    for (; it != ite; ++it)
      if ((masks[it.index()])[T::ID]) (*it)->receipt(msg);
  }

  template<class RECEIVER>
     linkmsg_sender<RECEIVER>::~linkmsg_sender() {
    typename RECEIVERTAB::tas_iterator it = receivers.tas_begin(),
      ite = receivers.tas_end();
    for (; it != ite; ++it)
      (*it)->out_sender(*this);
  }

  /* ******************************************************************** */
  /*	VIRTUAL RECEIVER CLASS, declaration                               */
  /* ******************************************************************** */

  class virtual_linkmsg_receiver {
  public :
    
    typedef int msg_id_type;
    
  protected :
    
    dal::dynamic_tas<virtual_linkmsg_sender *> senders;
    
    template<class ITER>
    void sup_sender_(virtual_linkmsg_sender *s, ITER b, const ITER &e,
		     bool t = true);

  public :
    
    template <class SENDER, class RECEIVER>
    void add_sender(SENDER &s, RECEIVER &r, mask m = mask())
      { senders.add(&s); s.add_receiver(r, m); }
    template <class SENDER>
    void sup_sender(SENDER &s)
      { sup_sender_(&s, senders.tas_begin(), senders.tas_end()); }
    template <class SENDER>
    void out_sender(SENDER &s)
      { sup_sender_(&s, senders.tas_begin(), senders.tas_end(), false); }
    virtual_linkmsg_receiver(void) {}
    template <class SENDER>
    virtual_linkmsg_receiver(SENDER &s, mask m = mask())
      { senders.add(&s); s.add_receiver(*this, m); }
    ~virtual_linkmsg_receiver()
      { sup_sender_(NULL, senders.tas_begin(), senders.tas_end()); }
    const virtual_linkmsg_receiver& operator =(const
					       virtual_linkmsg_receiver &) {
      DAL_THROW(internal_error, "The copy of this object doesn't work");
      return *this;
    }
    virtual_linkmsg_receiver(const virtual_linkmsg_receiver &v) { *this = v; }
  };


  /* ******************************************************************** */
  /*	VIRTUAL RECEIVER CLASS. Members functions.          	          */
  /* ******************************************************************** */

  template <class ITER>
    void virtual_linkmsg_receiver::sup_sender_(virtual_linkmsg_sender *s,
					       ITER b, const ITER &e, bool t) {
    for ( ; b != e; ++b)
      if (((*b) == s) || (s == NULL)) {
	if (t) senders[b.index()]->sup_receiver(this);
	senders.sup(b.index());
      }
  }

}

#endif /* GETFEM_LINKMSG_H__ */

