/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  LINK MeSsaGE library (linkmsg)                               */
/* File    :  linkmsg.h : structure for managing Messages between objects. */
/*     									   */
/*                                                                         */
/* Date : March 21, 1999.                                                  */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __LINKMSG_H
#define __LINKMSG_H

#include <dal_basic.h>
#include <dal_tas.h>
#include <queue>

namespace lmsg
{
  /* ******************************************************************** */
  /*	MASK TYPE.                                		          */
  /* ******************************************************************** */

  class mask
  {
    protected :

      dal::uint32_type tab; /* uint 64 does'nt exist every where.         */

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

  class virtual_linkmsg_sender {
  public : 
    virtual void sup_receiver(void *) = 0;
    virtual ~virtual_linkmsg_sender() {}
  };


  template<class RECEIVER> class linkmsg_sender : public virtual_linkmsg_sender
  {
    public :
      
      typedef int msg_id_type;
      typedef size_t size_type;

    protected :

      dal::dynamic_tas<RECEIVER *, 3> receivers;
      dal::dynamic_array<mask, 3> masks;

    public :

      void add_receiver(RECEIVER &re, mask m = mask())
      { masks[receivers.add(&re)] = m; }
      void sup_receiver(void *);
      template<class T> void send(const T &) const;
      linkmsg_sender(void) {};
      virtual ~linkmsg_sender();
  };

  template<class RECEIVER>
    void linkmsg_sender<RECEIVER>::sup_receiver(void *p)
  {
    for (int i = receivers.card() - 1; i >= 0; --i)
      if (receivers[i] == p) receivers.sup(i);
    receivers.compact();
  }

  template<class RECEIVER> template<class T> 
    void linkmsg_sender<RECEIVER>::send(const T &msg) const
  {
    for (int i = receivers.card() - 1; i >= 0; --i)
      if ((masks[i])[int(msg)]) receivers[i]->receipt(msg);
  }

  template<class RECEIVER>
     linkmsg_sender<RECEIVER>::~linkmsg_sender()
  {
    for (int i = receivers.card() - 1; i >= 0; --i)
      receivers[i]->out_sender(*this);
  }

  /* ******************************************************************** */
  /*	VIRTUAL RECEIVER CLASS, declaration                               */
  /* ******************************************************************** */

  class virtual_linkmsg_receiver
  {
    public :
      
      typedef int msg_id_type;

    protected :
      
      dal::dynamic_tas<virtual_linkmsg_sender *> senders;

      template<class ITER>
	void _sup_sender(virtual_linkmsg_sender *s, ITER b, const ITER &e,
			 bool t = true);

    public :

      template <class SENDER, class RECEIVER>
        void add_sender(SENDER &s, RECEIVER &r, mask m = mask())
      { senders.add(&s); s.add_receiver(r, m); }
      template <class SENDER>
        void sup_sender(SENDER &s)
      { _sup_sender(&s, senders.tas_begin(), senders.tas_end()); }
      template <class SENDER>
        void out_sender(SENDER &s)
      { _sup_sender(&s, senders.tas_begin(), senders.tas_end(), false); }
      virtual_linkmsg_receiver(void) {}
      template <class SENDER>
	virtual_linkmsg_receiver(SENDER &s, mask m = mask())
      { senders.add(&s); s.add_receiver(*this, m); }
      ~virtual_linkmsg_receiver()
      { _sup_sender(NULL, senders.tas_begin(), senders.tas_end()); }
  };


  /* ******************************************************************** */
  /*	VIRTUAL RECEIVER CLASS. Members functions.          	          */
  /* ******************************************************************** */

  template <class ITER>
    void virtual_linkmsg_receiver::_sup_sender(virtual_linkmsg_sender *s,
					       ITER b, const ITER &e, bool t)
  {
    for ( ; b != e; ++b)
      if (((*b) == s) || (s == NULL))
      {
	if (t) senders[b.index()]->sup_receiver(this);
	senders.sup(b.index());
      }
  }

}

#endif /* __LINKMSG_H */
