/*===========================================================================

Copyright (C) 2002-2012 Yves Renard

This file is a part of GETFEM++

Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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


#include "getfem/dal_static_stored_objects.h"
#include "getfem/dal_singleton.h"
#include <map>
#include <list>
#include <set>
#include <algorithm>
#include <deque>

namespace dal {

	void collect_static_stored_objects_garbage(){
		if (getfem::me_is_multithreaded_now()) return;
		std::set<const static_stored_object*> total_garbage;
		for(int thread=0; thread < int(getfem::num_threads()); thread++) {
                  static_stored_objects_garbage& garbage = 
                    dal::singleton<static_stored_objects_garbage>::instance(thread);
                  total_garbage.insert(garbage.garbage_pointers.begin(),garbage.garbage_pointers.end());
                  garbage.garbage_pointers.clear();
		}
	        for(std::set<const static_stored_object*>::iterator it = total_garbage.begin();
		    it!=total_garbage.end();it++)
		    if ((**it).ref_sum()==0) delete *it;
		}

	// Pointer to an object with the dependencies
	struct enr_static_stored_object {
		pstatic_stored_object p;
		bool valid;
		permanence perm;
		std::set<pstatic_stored_object> dependent_object;
		std::set<pstatic_stored_object> dependencies;
		enr_static_stored_object(pstatic_stored_object o, permanence perma)
			: p(o), valid(true), perm(perma) {}
		enr_static_stored_object(void)
			: p(0), valid(true), perm(STANDARD_STATIC_OBJECT) {}
	};

	// Pointer to a key with a coherent order
	struct enr_static_stored_object_key {
		pstatic_stored_object_key p;
		bool operator < (const enr_static_stored_object_key &o) const
		{ return (*p) < (*(o.p)); }
		enr_static_stored_object_key(pstatic_stored_object_key o) : p(o) {}
	};

	// Storing array types
	typedef std::map<enr_static_stored_object_key, enr_static_stored_object>
		stored_object_tab;
	struct stored_key_tab : public std::map<pstatic_stored_object,
		pstatic_stored_object_key> {
			~stored_key_tab() {
				for (iterator it = begin(); it != end(); ++it) delete it->second;
			}
	};

    // Gives a pointer to a key of an object from its pointer
	pstatic_stored_object_key key_of_stored_object(pstatic_stored_object o, size_t thread) {
          stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance(int(thread));
		stored_key_tab::iterator it = stored_keys.find(o);
		if (it != stored_keys.end()) return it->second;
		return 0;
	}

	/** Gives a pointer to a key of an object from its pointer 
    (searches in the storage of all threads) */
	pstatic_stored_object_key key_of_stored_object(pstatic_stored_object o) {
		for(size_t thread = 0; thread<getfem::num_threads();thread++){
			pstatic_stored_object_key key = key_of_stored_object(o,thread);
			if (key) return key;
		}
		return 0;
	}



	// Test if an object is stored (in current thread storage).
	bool exists_stored_object(pstatic_stored_object o) {
		stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance();
		return (stored_keys.find(o) != stored_keys.end());
	}

	// Test if an object is stored (in any of the thread's storage).
	bool exists_stored_object_all_threads(pstatic_stored_object o) {
          for(int thread = 0; thread < int(getfem::num_threads()); thread++){
			stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance(thread);
			if (stored_keys.find(o) != stored_keys.end()) return true;
		}
		return false;
	}


	/* Gives a pointer to an object from a key pointer (by looking in the
	current thread storage)*/
	pstatic_stored_object search_stored_object(pstatic_stored_object_key k) {
		stored_object_tab& stored_objects
			= dal::singleton<stored_object_tab>::instance();
		stored_object_tab::iterator it
			= stored_objects.find(enr_static_stored_object_key(k));
		if (it != stored_objects.end()) return it->second.p;
		return 0;
	}

	/* Search for an object in the storage of all threads*/
	pstatic_stored_object
        search_stored_object_all_threads(pstatic_stored_object_key k) {
          for(int thread = 0; thread<int(getfem::num_threads());thread++){
            stored_object_tab& stored_objects
            = dal::singleton<stored_object_tab>::instance(thread);
            stored_object_tab::iterator it
            = stored_objects.find(enr_static_stored_object_key(k));
            if (it != stored_objects.end()) return it->second.p;
          }
          return 0;
	}

	/** Gives an iterator on stored object from a pointer object 
	also indicates in which thread storage the object is found*/
	static inline stored_object_tab::iterator 
        iterator_of_object(pstatic_stored_object o, size_t& thread_found) {
          thread_found=0;
          pstatic_stored_object_key k = key_of_stored_object(o);
          if (k) {
            for(int thread = 0; thread<int(getfem::num_threads());thread++){
              stored_object_tab& stored_objects
                = dal::singleton<stored_object_tab>::instance(thread);
              stored_object_tab::iterator it
                = stored_objects.find(enr_static_stored_object_key(k));
              if (it != stored_objects.end()) {thread_found=thread; return it;}
            }
            GMM_ASSERT1(false,"Object has key but cannot be found");
          }
          return dal::singleton<stored_object_tab>::instance().end();
	}

	/* Gives an iterator on stored object from a pointer object*/
	static inline stored_object_tab::iterator 
        iterator_of_object(pstatic_stored_object o) {
          size_t thread;
          return iterator_of_object(o,thread);
	}

	// Test the validity of arrays
	void test_stored_objects(void) {
		stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance();
		for (stored_key_tab::iterator it = stored_keys.begin();
			it != stored_keys.end(); ++it)
			iterator_of_object(it->first);
		stored_object_tab& stored_objects
			= dal::singleton<stored_object_tab>::instance();
		for (stored_object_tab::iterator it = stored_objects.begin();
			it != stored_objects.end(); ++it)
			GMM_ASSERT1(iterator_of_object(it->second.p) != stored_objects.end(),
			"Object has key but cannot be found");
	}

	/* Add a dependency, object o1 will depend on object o2 */
	void add_dependency(pstatic_stored_object o1,
                            pstatic_stored_object o2) {
          stored_object_tab& stored_objects
          = dal::singleton<stored_object_tab>::instance();    
          std::vector<stored_object_tab*> all_stored_objects;
          for(size_t i= 0; i<getfem::num_threads();i++) 
            all_stored_objects.push_back
              (&(dal::singleton<stored_object_tab>::instance(int(i))));
          size_t thread1, thread2;
          stored_object_tab::iterator it1 = iterator_of_object(o1,thread1);
          stored_object_tab::iterator it2 = iterator_of_object(o2,thread2);
          if (it1 != (*all_stored_objects[thread1]).end() && 
              it2 != (*all_stored_objects[thread2]).end()) {
            getfem::omp_guard local_lock;
            it2->second.dependent_object.insert(o1);
            it1->second.dependencies.insert(o2);
          }
          else {
            cerr << "Problem adding dependency between " << o1 << " of type "
            << typeid(*o1).name() << " and " << o2 << " of type "
            << typeid(*o2).name() << ". ";
            if (it1 == (*all_stored_objects[thread1]).end()) 
              cerr << "First object does not exist.";
            if (it2 == (*all_stored_objects[thread2]).end()) 
              cerr << "Second object does not exist.";
            cerr<<" thread N = "<<getfem::this_thread();
            cerr << endl;
            assert(false);
            GMM_ASSERT1(false, "Add_dependency : Inexistent object");
          }
	}
        
	/*remove a dependency (from storages of all threads). 
	Return true if o2 has no more dependent object. */
	bool del_dependency(pstatic_stored_object o1, pstatic_stored_object o2) {
		stored_object_tab& stored_objects
			= dal::singleton<stored_object_tab>::instance();
		stored_object_tab::iterator it1 = iterator_of_object(o1);
		stored_object_tab::iterator it2 = iterator_of_object(o2);
		if (it1 != stored_objects.end() && it2 != stored_objects.end()) {
			getfem::omp_guard local_lock; 
			it2->second.dependent_object.erase(o1);
			it1->second.dependencies.erase(o2);
			return it2->second.dependent_object.empty();
		}
		return true;
	}

	// Add an object (local thread storage) with two optional dependencies
	void add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
		permanence perm) {

			stored_object_tab& stored_objects
				= dal::singleton<stored_object_tab>::instance();
			stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance();
			GMM_ASSERT1(stored_keys.find(o) == stored_keys.end(),
				"This object has already been stored, "
				"possibly with another key");
			stored_keys[o] = k;
			stored_objects[enr_static_stored_object_key(k)]
			= enr_static_stored_object(o, perm);
	}



	/* Only delete the list of objects but not the dependencies */
	static void basic_delete(std::list<pstatic_stored_object> &to_delete){
          for(int thread=0; thread < int(getfem::num_threads()); thread++) {
            stored_object_tab& stored_objects
            = dal::singleton<stored_object_tab>::instance(thread);
            stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance(thread);
            std::list<pstatic_stored_object>::iterator it;
            for (it = to_delete.begin(); it != to_delete.end(); ++it) {
              // cout << "delete object " << (*it).get() << " of type "
              //      << typeid(*(*it)).name() << endl;
              pstatic_stored_object_key k = key_of_stored_object(*it,thread);
              stored_object_tab::iterator ito = stored_objects.end();
              if (k) ito = stored_objects.find(k);
              if (k) {getfem::omp_guard local_lock; stored_keys.erase(*it);}
              if (ito != stored_objects.end()) {
                getfem::omp_guard local_lock;
                delete ito->first.p;
                stored_objects.erase(ito);
              }
            }
          }
	}

	// Delete a list of objects and their dependencies
	void del_stored_objects_immediate(std::list<pstatic_stored_object> &to_delete,
		bool ignore_unstored) {
			stored_object_tab& stored_objects
				= dal::singleton<stored_object_tab>::instance();
			std::list<pstatic_stored_object>::iterator it, itnext;
			for (it = to_delete.begin(); it != to_delete.end(); it = itnext) {
				itnext = it; itnext++;
				stored_object_tab::iterator ito = iterator_of_object(*it);
				if (ito == stored_objects.end()) {
					if (ignore_unstored)
						to_delete.erase(it);
					else
						GMM_ASSERT1(false, "This object is not stored : " << it->get()
						<< " typename: " << typeid(*it->get()).name());
				}
				else
					iterator_of_object(*it)->second.valid = false;
			}
			std::set<pstatic_stored_object>::iterator itd;
			for (it = to_delete.begin(); it != to_delete.end(); ++it) {
				if (*it) {
					stored_object_tab::iterator ito = iterator_of_object(*it);
					GMM_ASSERT1(ito != stored_objects.end(), "An object disapeared !");
					ito->second.valid = false;
					std::set<pstatic_stored_object> dep = ito->second.dependencies;
					for (itd = dep.begin(); itd != dep.end(); ++itd) {
						if (del_dependency(*it, *itd)) {
							stored_object_tab::iterator itod=iterator_of_object(*itd);
							if (itod->second.perm == AUTODELETE_STATIC_OBJECT
								&& itod->second.valid) {
									itod->second.valid = false;
									to_delete.push_back(*itd);
							}
						}
					}
					for (itd = ito->second.dependent_object.begin();
						itd != ito->second.dependent_object.end(); ++itd) {
							stored_object_tab::iterator itod=iterator_of_object(*itd);
							if (itod != stored_objects.end()) {
								GMM_ASSERT1(itod->second.perm != PERMANENT_STATIC_OBJECT,
									"Trying to delete a permanent object " << *itd);
								if (itod->second.valid) {
									itod->second.valid = false;
									to_delete.push_back(itod->second.p);
								}
							}
					}
				}
			}
			basic_delete(to_delete);
	}


	class object_terminator{
		struct deletion_unit{
			bool ignore_unstored;
			std::list<pstatic_stored_object> to_delete;
			deletion_unit(bool ignore, const std::list<pstatic_stored_object>& list) : 
				ignore_unstored(ignore), to_delete(list){}
		};

	public:

		void add_to_deletion(std::list<pstatic_stored_object> &to_delete,
			bool ignore_unstored){
				deletion_unit u(ignore_unstored,to_delete);
				deletion_list.push_back(u);
		}

		void clear(){deletion_list.clear();}
		void delete_stored_content(){
			for(std::list<deletion_unit>::iterator it=deletion_list.begin(); it!=deletion_list.end();it++)
			    del_stored_objects_immediate(it->to_delete,it->ignore_unstored);

			deletion_list.clear();
		}
	private:
		std::list<deletion_unit> deletion_list;
	};


	void del_stored_objects(std::list<pstatic_stored_object> &to_delete,
		bool ignore_unstored){
			//if (getfem::me_is_multithreaded_now()){
			//	object_terminator& terminator = dal::singleton<object_terminator>::instance();
			//	terminator.add_to_deletion(to_delete,ignore_unstored);
			//} else 
				del_stored_objects_immediate(to_delete,ignore_unstored);
	}


	void flush_deleted_objects(){
          GMM_ASSERT1(!getfem::me_is_multithreaded_now(), 
                      "Actual object deletion should be done outside " 
                      "the parallel region (preferably right after)");
          for(int thread=0; thread < int(getfem::num_threads()); thread++){
            object_terminator& terminator = 
              dal::singleton<object_terminator>::instance(thread);
            terminator.delete_stored_content();
          }
	}


	// Delete an object and its dependencies
	void del_stored_object(pstatic_stored_object o, bool ignore_unstored) {
		std::list<pstatic_stored_object> to_delete;
		to_delete.push_back(o);
		del_stored_objects(to_delete, ignore_unstored);
	}

	// Delete all the object whose perm is greater or equal to perm
	void del_stored_objects(permanence perm) {
		stored_object_tab& stored_objects
			= dal::singleton<stored_object_tab>::instance();
		if (perm == PERMANENT_STATIC_OBJECT) perm = STRONG_STATIC_OBJECT;
		std::list<pstatic_stored_object> to_delete;
		stored_object_tab::iterator it;
		for (it = stored_objects.begin(); it != stored_objects.end(); ++it)
			if (it->second.perm >= perm)
				to_delete.push_back(it->second.p);
		del_stored_objects(to_delete, false);
	}

	// List the stored objects for debugging purpose
	void list_stored_objects(std::ostream &ost) {
          for(int thread=0; thread<int(getfem::num_threads());thread++){
            stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance(thread);
            if (stored_keys.begin() == stored_keys.end())
              ost << "No static stored objects" << endl;
            else
              ost << "Static stored objects" << endl;
            for (stored_key_tab::iterator it = stored_keys.begin();
                 it != stored_keys.end(); ++it) {
              ost << "Object: " << it->first << " typename: "
                  << typeid(*it->first).name() << endl;
            }
          }
	}

	// Number of stored objects
	size_t nb_stored_objects(void) {
		long num_objects=0;
		for(size_t thread=0;thread<getfem::num_threads();thread++){
			stored_key_tab& stored_keys = dal::singleton<stored_key_tab>::instance();
			num_objects+=stored_keys.size();
		}
		return num_objects;
	}

}
