/*===========================================================================

 Copyright (C) 2005-2017 Yves Renard

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

===========================================================================*/

#include "getfem/getfem_mesh_region.h"
#include "getfem/getfem_mesh.h"
#include "getfem/getfem_omp.h"

namespace getfem {
  typedef mesh_region::face_bitset face_bitset;

  mesh_region::mesh_region(const mesh_region &other)
    : p(std::make_shared<impl>()), id_(size_type(-2)), parent_mesh(0),
      index_updated(false)
  {
    this->operator=(other);
  }


  mesh_region::mesh_region() : p(std::make_shared<impl>()), id_(size_type(-2)),
                               type_(size_type(-1)),
    partitioning_allowed(true), parent_mesh(0), index_updated(false)
  {
    if (me_is_multithreaded_now()) prohibit_partitioning();
  }

  mesh_region::mesh_region(size_type id__) : id_(id__), type_(size_type(-1)),
    partitioning_allowed(true), parent_mesh(0), index_updated(false)
  { }

  mesh_region::mesh_region(mesh& m, size_type id__, size_type type) :
    p(std::make_shared<impl>()), id_(id__), type_(type), partitioning_allowed(true), parent_mesh(&m),
    index_updated(false)
  {
    if (me_is_multithreaded_now()) prohibit_partitioning();
  }

  mesh_region::mesh_region(const dal::bit_vector &bv) :
    p(std::make_shared<impl>()), id_(size_type(-2)), type_(size_type(-1)),
    partitioning_allowed(true), parent_mesh(0), index_updated(false)
  {
    if (me_is_multithreaded_now()) prohibit_partitioning();
    add(bv);
  }

  void mesh_region::touch_parent_mesh()
  {
    if (parent_mesh) parent_mesh->touch_from_region(id_);
  }

  const mesh_region& mesh_region::from_mesh(const mesh &m) const
  {
    if (!p.get())
    {
      mesh_region *r = const_cast<mesh_region*>(this);
      if (id_ == size_type(-1))
      {
        r->p = std::make_shared<impl>();
        r->add(m.convex_index());
      }
      else if (id_ != size_type(-2))
      {
        *r = m.region(id_);
      }
    }
    index_updated = false;
    return *this;
  }

  mesh_region& mesh_region::operator=(const mesh_region &from)
  {
    if (!this->parent_mesh && !from.parent_mesh)
    {
      this->id_ = from.id_;
      this->type_ = from.type_;
      this->partitioning_allowed = from.partitioning_allowed;
      if (from.p.get()) {
        if (!this->p.get()) this->p = std::make_shared<impl>();
        this->wp() = from.rp();
      }
      else
        this->p.reset();
    }
    else if (!this->parent_mesh) {
      this->p = from.p;
      this->id_ = from.id_;
      this->type_ = from.type_;
      this->parent_mesh = from.parent_mesh;
      this->partitioning_allowed = from.partitioning_allowed;
    }
    else {
      if (from.p.get())
      {
        this->wp() = from.rp();
        this->type_= from.get_type();
        this->partitioning_allowed = from.partitioning_allowed;
      }
      else if (from.id_ == size_type(-1)) {
        this->clear();
        this->add(this->parent_mesh->convex_index());
        this->type_ = size_type(-1);
        this->partitioning_allowed = true;
      }
      touch_parent_mesh();
    }
    index_updated = false;
    return *this;
  }

  bool mesh_region::compare(const mesh &m1, const mesh_region &mr,
                            const mesh &m2) const {
    if (&m1 != &m2) return false;
    if (!p.get() && !mr.p.get()) return (id_ == mr.id_);
    this->from_mesh(m1);
    mr.from_mesh(m2);
    if (this->p.get() && !(mr.p.get())) return false;
    if (!(this->p.get()) && (mr.p.get())) return false;
    if (this->p.get())
      if (this->p.get()->m != mr.p.get()->m) return false;
    return true;
  }

  face_bitset mesh_region::operator[](size_t cv) const
  {
    map_t::const_iterator it = rp().m.find(cv);
    if (it != rp().m.end()) return (*it).second;
    else return face_bitset();
  }

  void mesh_region::update_partition_iterators() const
  {
    if (index_updated) return;
    itbegin = partition_begin();
    itend   = partition_end  ();
    index_updated = true;
  }
  mesh_region::const_iterator
    mesh_region::partition_begin( ) const
  {
    size_type region_size = rp().m.size();
    if (region_size < num_threads())
    { //for small regions: put the whole region into zero thread
      if (this_thread() == 0) return rp().m.begin();
      else return rp().m.end();
    }
    size_type partition_size = static_cast<size_type>
      (std::ceil(static_cast<scalar_type>(region_size)/
       static_cast<scalar_type >(num_threads())));
    size_type index_begin = partition_size * this_thread();
    if (index_begin >= region_size ) return  rp().m.end();

    const_iterator it = rp().m.begin();
    for (size_type i=0;i<index_begin;++i) ++it;
    return it;
  }

  mesh_region::const_iterator
    mesh_region::partition_end( ) const
  {
    size_type region_size = rp().m.size();
    if (region_size < num_threads()) return rp().m.end();

    size_type partition_size = static_cast<size_type>
      (std::ceil(static_cast<scalar_type>(region_size)/
       static_cast<scalar_type >(num_threads())));
    size_type index_end = partition_size * (this_thread() + 1);
    if (index_end >= region_size ) return  rp().m.end();

    const_iterator it = rp().m.begin();
    for (size_type i=0;i<index_end && it!=rp().m.end();++i) ++it;
    return it;
  }

  mesh_region::const_iterator mesh_region::begin( ) const
  {
    GMM_ASSERT1(p != 0, "Internal error");
    if (me_is_multithreaded_now() && partitioning_allowed)
    {
      update_partition_iterators();
      return itbegin;
    }
    else { return rp().m.begin(); }
  }

  mesh_region::const_iterator mesh_region::end  ( ) const
  {
    if (me_is_multithreaded_now() && partitioning_allowed)
    {
      update_partition_iterators();
      return itend;
    }
    else return rp().m.end();
  }

  void  mesh_region::allow_partitioning()
  {
    if (me_is_multithreaded_now()) partitioning_allowed = true;
    else partitioning_allowed.all_threads() = true;
  }

  void mesh_region::bounding_box(base_node& Pmin, base_node& Pmax) const {
    auto &mesh = *this->get_parent_mesh();
    for (auto cv : dal::bv_iterable_c(index())) {
      for (const auto &pt : mesh.points_of_convex(cv)) {
        for (auto j = 0; j < Pmin.size(); j++){
          Pmin[j] = std::min(Pmin[j], pt[j]);
          Pmax[j] = std::max(Pmax[j], pt[j]);
        }
      }
    }
  }

  void  mesh_region::prohibit_partitioning()
  {
    if (me_is_multithreaded_now()) partitioning_allowed = false;
    else partitioning_allowed.all_threads() = false;
  }

  bool mesh_region::is_partitioning_allowed() const
  {
    return partitioning_allowed;
  }

  /* may be optimized .. */
  const dal::bit_vector&  mesh_region::index() const
  {
    GMM_ASSERT1(p.get(), "Use from_mesh on that region before");
    dal::bit_vector& convex_index = rp().index_.thrd_cast();
    convex_index.clear();
    for (const_iterator it = begin(); it != end(); ++it)
    {
      if (it->second.any()) convex_index.add(it->first);
    }
    return convex_index;
  }

  void mesh_region::add(const dal::bit_vector &bv)
  {
    for (dal::bv_visitor i(bv); !i.finished(); ++i)
    {
      wp().m[i].set(0,1);
    }
    touch_parent_mesh();
    index_updated = false;
  }

  void mesh_region::add(size_type cv, short_type f)
  {
    wp().m[cv].set(short_type(f+1),1);
    touch_parent_mesh();
    index_updated = false;
  }

  void mesh_region::sup_all(size_type cv)
  {
    map_t::iterator it = wp().m.find(cv);
    if (it != wp().m.end()) {
      wp().m.erase(it);
      touch_parent_mesh();
    }
    index_updated = false;
  }

  void mesh_region::sup(size_type cv, short_type f)
  {
    map_t::iterator it = wp().m.find(cv);
    if (it != wp().m.end()) {
      (*it).second.set(short_type(f+1),0);
      if ((*it).second.none()) wp().m.erase(it);
      touch_parent_mesh();
    }
    index_updated = false;
  }

  void mesh_region::clear()
  {
    wp().m.clear(); touch_parent_mesh();
    index_updated = false;
  }

  void mesh_region::clean()
  {
    for (map_t::iterator it = wp().m.begin(), itn;
      it != wp().m.end(); it = itn)
    {
      itn = it;
      ++itn;
      if ( !(*it).second.any() ) wp().m.erase(it);
    }
    touch_parent_mesh();
    index_updated = false;
  }


  void mesh_region::swap_convex(size_type cv1, size_type cv2)
  {
    map_t::iterator it1 = wp().m.find(cv1), it2 = wp().m.find(cv2),
      ite = wp().m.end();
    face_bitset f1, f2;

    if (it1 != ite) f1 = it1->second;
    if (it2 != ite) f2 = it2->second;
    if (!f1.none()) wp().m[cv2] = f1;
    else if (it2 != ite) wp().m.erase(it2);
    if (!f2.none()) wp().m[cv1] = f2;
    else if (it1 != ite) wp().m.erase(it1);
    touch_parent_mesh();
    index_updated = false;
  }

  bool mesh_region::is_in(size_type cv, short_type f) const
  {
    GMM_ASSERT1(p.get(), "Use from mesh on that region before");
    map_t::const_iterator it = rp().m.find(cv);
    if (it == rp().m.end() || short_type(f+1) >= MAX_FACES_PER_CV) return false;
    return ((*it).second)[short_type(f+1)];
  }

  bool mesh_region::is_in(size_type cv, short_type f, const mesh &m) const
  {
    if (p.get()) {
      map_t::const_iterator it = rp().m.find(cv);
      if (it == rp().m.end() || short_type(f+1) >= MAX_FACES_PER_CV)
        return false;
      return ((*it).second)[short_type(f+1)];
    }
    else
    {
      if (id() == size_type(-1)) return true;
      else return m.region(id()).is_in(cv, f);
    }
  }



  bool mesh_region::is_empty() const
  {
    return rp().m.empty();
  }

  bool mesh_region::is_only_convexes() const
  {
    return is_empty() ||
      (or_mask()[0] == true && or_mask().count() == 1);
  }

  bool mesh_region::is_only_faces() const
  {
    return is_empty() || (and_mask()[0] == false);
  }

  face_bitset mesh_region::faces_of_convex(size_type cv) const
  {
    map_t::const_iterator it = rp().m.find(cv);
    if (it != rp().m.end()) return ((*it).second) >> 1;
    else return face_bitset();
  }

  face_bitset mesh_region::and_mask() const
  {
    face_bitset bs;
    if (rp().m.empty()) return bs;
    bs.set();
    for (map_t::const_iterator it = rp().m.begin(); it != rp().m.end(); ++it)
      if ( (*it).second.any() )  bs &= (*it).second;
    return bs;
  }

  face_bitset mesh_region::or_mask() const
  {
    face_bitset bs;
    if (rp().m.empty()) return bs;
    for (map_t::const_iterator it = rp().m.begin(); it != rp().m.end(); ++it)
      if ( (*it).second.any() )  bs |= (*it).second;
    return bs;
  }

  size_type mesh_region::size() const
  {
    size_type sz=0;
    for (map_t::const_iterator it = begin(); it != end(); ++it)
      sz += (*it).second.count();
    return sz;
  }

  size_type mesh_region::unpartitioned_size() const
  {
    size_type sz=0;
    for (map_t::const_iterator it = rp().m.begin(); it != rp().m.end(); ++it)
      sz += (*it).second.count();
    return sz;
  }


  mesh_region mesh_region::intersection(const mesh_region &a,
                                        const mesh_region &b)
  {
    GMM_TRACE4("intersection of "<<a.id()<<" and "<<b.id());
    mesh_region r;
    /* we do not allow the "all_convexes" kind of regions
    for these operations as there are not intended to be manipulated
    (they only exist to provide a default argument to the mesh_region
    parameters of assembly procedures etc. */
    GMM_ASSERT1(a.id() !=  size_type(-1)||
                b.id() != size_type(-1), "the 'all_convexes' regions "
                "are not supported for set operations");
    if (a.id() == size_type(-1))
    {
      for (const_iterator it = b.begin();it != b.end(); ++it) r.wp().m.insert(*it);
      return r;
    }
    else if (b.id() == size_type(-1))
    {
      for (const_iterator it = a.begin();it != a.end(); ++it) r.wp().m.insert(*it);
      return r;
    }

    map_t::const_iterator
      ita = a.begin(), enda = a.end(),
      itb = b.begin(), endb = b.end();

    while (ita != enda && itb != endb) {
      if (ita->first < itb->first) ++ita;
      else if (ita->first > itb->first) ++itb;
      else {
        face_bitset maska = ita->second, maskb = itb->second, bs;
        if (maska[0] && !maskb[0]) bs = maskb;
        else if (maskb[0] && !maska[0]) bs = maska;
        else bs = maska & maskb;
        if (bs.any()) r.wp().m.insert(r.wp().m.end(), std::make_pair(ita->first,bs));
        ++ita; ++itb;
      }
    }
    return r;
  }

  mesh_region mesh_region::merge(const mesh_region &a,
                                 const mesh_region &b)
  {
    GMM_TRACE4("Merger of " << a.id() << " and " << b.id());
    mesh_region r;
    GMM_ASSERT1(a.id() != size_type(-1) &&
      b.id() != size_type(-1), "the 'all_convexes' regions "
      "are not supported for set operations");
    for (const_iterator it = a.begin();it != a.end(); ++it)
    {
      r.wp().m.insert(*it);
    }
    for (const_iterator it = b.begin();it != b.end(); ++it)
    {
      r.wp().m[it->first] |= it->second;
    }
    return r;
  }


  mesh_region mesh_region::subtract(const mesh_region &a,
                                    const mesh_region &b)
  {
    GMM_TRACE4("subtraction of "<<a.id()<<" and "<<b.id());
    mesh_region r;
    GMM_ASSERT1(a.id() != size_type(-1) &&
      b.id() != size_type(-1), "the 'all_convexes' regions "
      "are not supported for set operations");
    for (const_iterator ita = a.begin();ita != a.end(); ++ita)
      r.wp().m.insert(*ita);

    for (const_iterator itb = b.begin();itb != b.end(); ++itb)
    {
      size_type cv = itb->first;
      map_t::iterator it = r.wp().m.find(cv);
      if (it != r.wp().m.end()){
                it->second &= ~(itb->second);
                if (it->second.none()) r.wp().m.erase(it);
          }
    }
    return r;
  }

  int mesh_region::region_is_faces_of(const getfem::mesh& m1,
                                      const mesh_region &rg2,
                                      const getfem::mesh& m2) const {
    if (&m1 != &m2) return 0;
    int r = 1, partially = 0;
    for (mr_visitor cv(*this, m1); !cv.finished(); cv.next())
      if (cv.is_face() && rg2.is_in(cv.cv(),short_type(-1), m2))
        partially = -1;
      else
        r = 0;
    if (r == 1) return 1; else return partially;
  }

  size_type mesh_region::free_region_id(const getfem::mesh& m)
  {
    return m.regions_index().last_true()+1;
  }


  void mesh_region::error_if_not_faces() const
  {
    GMM_ASSERT1(is_only_faces(), "Expecting a set of faces, not convexes");
  }

  void mesh_region::error_if_not_convexes() const
  {
    GMM_ASSERT1(is_only_convexes(), "Expecting a set of convexes, not faces");
  }

  void mesh_region::error_if_not_homogeneous() const
  {
    GMM_ASSERT1(is_only_faces() || is_only_convexes(), "Expecting a set "
                "of convexes or a set of faces, but not a mixed set");
  }




#if GETFEM_PARA_LEVEL > 1

  mesh_region::visitor::visitor(const mesh_region &s, const mesh &m,
                                bool intersect_with_mpi) :
    cv_(size_type(-1)), f_(short_type(-1)), finished_(false)
  {
    if ((me_is_multithreaded_now() && s.partitioning_allowed)) {
      s.from_mesh(m);
      init(s);
    } else {
      if (s.id() == size_type(-1)) {
        if (intersect_with_mpi)
          init(m.get_mpi_region());
        else
          init(m.convex_index());
      } else if (s.id() == size_type(-2) && s.p.get()) {
        if (intersect_with_mpi) {
          mpi_rg = std::make_unique<mesh_region>(s);
          mpi_rg->from_mesh(m);
          m.intersect_with_mpi_region(*mpi_rg);
          init(*mpi_rg);
        } else
          init(s);
      } else {
        GMM_ASSERT1(s.id() != size_type(-2), "Internal error");
        if (intersect_with_mpi)
          init(m.get_mpi_sub_region(s.id()));
        else
          init(m.region(s.id()));
      }
    }
  }

#else

  mesh_region::visitor::visitor(const mesh_region &s, const mesh &m, bool)
    :cv_(size_type(-1)), f_(short_type(-1)), finished_(false)
  {
    if ((me_is_multithreaded_now() && s.partitioning_allowed)) {
      s.from_mesh(m);
      init(s);
    } else {
      if (s.id() == size_type(-1)) {
        init(m.convex_index());
      } else if (s.id() == size_type(-2) && s.p.get()) {
        init(s);
      } else {
        GMM_ASSERT1(s.id() != size_type(-2), "Internal error");
        init(m.region(s.id()));
      }
    }
  }

#endif


  bool mesh_region::visitor::next()
  {
    if (whole_mesh) {
      if (itb == iteb) { finished_ = true; return false; }
      cv_ = itb.index();
      c = 0;
      f_ = 0;
      ++itb; while (itb != iteb && !(*itb)) ++itb;
      return true;
    }
    while (c.none())
      {
        if (it == ite) { finished_=true; return false; }
        cv_ = it->first;
        c   = it->second;
        f_ = short_type(-1);
        ++it;
        if (c.none()) continue;
      }
    next_face();
    return true;
  }

  mesh_region::visitor::visitor(const mesh_region &s) :
    cv_(size_type(-1)), f_(short_type(-1)), finished_(false)
  {
    init(s);
  }

  void mesh_region::visitor::init(const dal::bit_vector &bv)
  {
    whole_mesh = true;
    itb = bv.begin(); iteb = bv.end();
    while (itb != iteb && !(*itb)) ++itb;
    next();
  }

  void mesh_region::visitor::init(const mesh_region &s)
  {
    whole_mesh = false;
    it  = s.begin();
    ite = s.end();
    next();
  }

  std::ostream & operator <<(std::ostream &os, const mesh_region &w)
  {
    if (w.id() == size_type(-1))
      os << " ALL_CONVEXES";
    else if (w.p.get())
    {
      for (mr_visitor cv(w); !cv.finished(); cv.next())
        {
          os << cv.cv();
          if (cv.is_face()) os << "/" << cv.f();
          os << " ";
        }
    }
    else
    {
      os << " region " << w.id();
    }
    return os;
  }

  struct dummy_mesh_region_ {
    mesh_region mr;
    dummy_mesh_region_() : mr() {}
  };

  const mesh_region &dummy_mesh_region()
  { return dal::singleton<dummy_mesh_region_>::instance().mr; }
}
