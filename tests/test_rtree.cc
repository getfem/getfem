#include <bgeot_rtree.h>
#include <ftool.h>
#include <dal_bit_vector.h>
using bgeot::base_node;
using bgeot::size_type;
using bgeot::dim_type;
using bgeot::rtree;

static bool quick = false;


static bool r1_ge_r2(const base_node& min1, const base_node& max1,
		     const base_node& min2, const base_node& max2) {
  for (size_type i=0; i < min1.size(); ++i) { 
    if (!(min1[i] <= min2[i] && max1[i] >= max2[i])) return false; 
  }
  return true;
}
static bool r1_inter_r2(const base_node& min1, const base_node& max1,
			const base_node& min2, const base_node& max2) {
  for (size_type i=0; i < min1.size(); ++i) 
    if (max1[i] < min2[i] || min1[i] > max2[i]) return false; 
  return true;
}

struct intersection_p {
  const base_node min,max;
  void print(std::ostream &o) { o << "intersects(" << min << ".." << max << ")"; }
  intersection_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
  bool operator()(const base_node& min2, const base_node& max2)
  { return r1_inter_r2(min,max,min2,max2); }
  bool accept(const base_node& min2, const base_node& max2) 
  { return operator()(min2,max2); }
};

struct contains_p {
  const base_node min,max;
  contains_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
  void print(std::ostream &o) { o << "contains(" << min << ".." << max << ")"; }
  bool operator()(const base_node& min2, const base_node& max2)
  { return r1_ge_r2(min2,max2,min,max); }
  bool accept(const base_node& min2, const base_node& max2) 
  { return operator()(min2,max2); }
};

struct contained_p {
  const base_node min,max;
  void print(std::ostream &o) { o << "contained(" << min << ".." << max << ")"; }
  contained_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
  bool accept(const base_node& min2, const base_node& max2)
  { return r1_inter_r2(min,max,min2,max2); }
  bool operator()(const base_node& min2, const base_node& max2) 
  { return r1_ge_r2(min,max,min2,max2); }
};

struct has_point_p {
  const base_node P;
  void print(std::ostream &o) { o << "has_point(" << P << ")"; }
  has_point_p(const base_node& P_) : P(P_) {}
  bool operator()(const base_node& min2, const base_node& max2) {
    for (size_type i=0; i < P.size(); ++i) 
      if (P[i] < min2[i] || P[i] > max2[i]) return false;
    return true;
  }
  bool accept(const base_node& min2, const base_node& max2) 
  { return operator()(min2,max2); }
};
  
template <typename Predicate>
static void brute_force_check(const std::vector<base_node>& rmin,
		       const std::vector<base_node>& rmax,			 
		       std::vector<size_type>& pbset, Predicate p) {    
  //cout << "brute_force_check("; p.print(cout); cout << ")\n";
  dal::bit_vector iset; 
  for (size_type i=0; i < rmin.size(); ++i) {
    if (p(rmin[i],rmax[i])) { 
      //cout << "brute_force_check: found match " << i << ":" << rmin[i] << "," << rmax[i] << "\n";
      iset.add(i);
    }
  }
  //cout << "brute force: iset = " << iset << "\n";
  dal::bit_vector bv; bv.merge_from(pbset);
  //cout << "      rtree: iset = " << bv << "\n";
  for (size_type i=0; i < pbset.size(); ++i) {
    assert(iset[pbset[i]]);
    iset[pbset[i]] = false;
  }
  if (iset.card()) { 
    cout << "the rtree found " << bv << " but failed to find " << iset << "\n";
    assert(iset.card() == 0);
  }
}

static void verify(const std::vector<base_node>& rmin, const std::vector<base_node>& rmax, bgeot::rtree& tree) {
  size_type N=rmin.front().size();
  std::vector<size_type> pbset;
  //tree.dump();
  base_node extent(rmin[0].size()); 
  for (size_type i=0; i < rmin.size(); ++i)
    for (size_type k=0; k < N; ++k) { 
      assert(rmax[i][k]-rmin[i][k] >= 0.);
      extent[k] = std::max(extent[k], rmax[i][k]-rmin[i][k]);
    }

  for (size_type i=0; i < 100; ++i) {
    base_node min(N), max(N);
    for (size_type k=0; k < N; ++k) { min[k] = gmm::random(double()*1.3); max[k] = min[k]+gmm::random()*0.1; }
    tree.find_containing_boxes(min,max,pbset);
    brute_force_check(rmin,rmax,pbset,contains_p(min,max));

    tree.find_intersecting_boxes(min,max,pbset);
    brute_force_check(rmin,rmax,pbset,intersection_p(min,max));

    tree.find_contained_boxes(min,max,pbset);
    brute_force_check(rmin,rmax,pbset,contained_p(min,max));

    tree.find_boxes_at_point(min,pbset);
    brute_force_check(rmin,rmax,pbset,has_point_p(min));

    tree.find_boxes_at_point(max,pbset);
    brute_force_check(rmin,rmax,pbset,has_point_p(max));
  }
  for (size_type i=0; i < rmin.size(); ++i) {
    base_node min2(rmin[i]); for (size_type k=0; k < N; ++k) { min2[k] -= extent[k]*gmm::random()*0.1; }
    base_node max2(rmax[i]); for (size_type k=0; k < N; ++k) { max2[k] += extent[k]*gmm::random()*0.1; }
    base_node min3(rmin[i]); for (size_type k=0; k < N; ++k) { min3[k] += extent[k]*gmm::random()*0.001; }
    base_node max3(rmax[i]); for (size_type k=0; k < N; ++k) { max3[k] -= extent[k]*gmm::random()*0.00001; }
    tree.find_boxes_at_point(rmin[i],pbset);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
    tree.find_boxes_at_point(rmax[i],pbset);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
    tree.find_contained_boxes(min2,max2,pbset);
    assert(pbset.size() >= 1);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
    //cout << "i=" << i << " find_contained_boxes " << rmin[i] << rmax[i] << "\n";
    tree.find_contained_boxes(rmin[i],rmax[i],pbset);
    assert(pbset.size() >= 1);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
    tree.find_containing_boxes(rmin[i],rmax[i],pbset);
    assert(pbset.size() >= 1);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
    { 
      contains_p p(min3,max3);
      tree.find_containing_boxes(min3,max3,pbset);
      //cout << "i=" << i << ", " << rmin[i] << rmax[i] << ", min3=" << min3 << ", max3=" << max3 << "\n";
      if (p(rmin[i],rmax[i])) {
	assert(pbset.size() >= 1);
	assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
      } else {
	assert(std::find(pbset.begin(), pbset.end(), i) == pbset.end());
      }
    }
    //cout << "i=" << i << " find_intersecting_boxes " << rmin[i] << rmax[i] << ", " << min2 << max3 << "\n";
    tree.find_intersecting_boxes(min2,max3,pbset);
    assert(pbset.size() >= 1);
    assert(std::find(pbset.begin(), pbset.end(), i) != pbset.end());
  }
}

static void check_tree() {
  bgeot::rtree tree;
  tree.add_box(base_node(1.0,0.),base_node(1.5,0.));
  tree.add_box(base_node(2.0,0.),base_node(3.0,0.),2);
  tree.add_box(base_node(1.5,0.),base_node(2.2,0.),1);
  std::vector<size_type> pi; tree.find_boxes_at_point(base_node(2.8,0.),pi);
  tree.dump();
  assert(pi.size()==1 && pi[0] == 2);  
  tree.clear();
  std::vector<base_node> rmin, rmax;
  cout << "1D checks\n";
  for (int C=5; C < 15; ++C) {
    static double dec[] = {.1, .5, .7423, .99, 1.,1.0001,1.25,1.5,2.,4.,4.1};
    for (size_type d=0; d < sizeof(dec)/sizeof(dec[0]); ++d) {
      //cout << "C=" << C << ",  DEC=" << dec[d] << "\n";
      for (int i=-C; i < C-1; ++i) {
	base_node a(1),b(1); a[0] = i/double(C); b[0]=(i+dec[d])/double(C);
	rmin.push_back(a); rmax.push_back(b);
	tree.add_box(rmin.back(),rmax.back());
      }
      verify(rmin, rmax, tree);
      tree.clear(); rmin.clear(); rmax.clear();
    }
  }

  cout << "2D random check\n";
  for (size_type i=0; i < 600; ++i) {
    rmin.push_back(base_node(gmm::random(double()), gmm::random(double())));
    rmax.push_back(rmin.back() + base_node(1.+gmm::random(), 1.+gmm::random())/10.);
    tree.add_box(rmin.back(),rmax.back());
  }
  verify(rmin, rmax, tree);

  cout << "2D/1D random check\n";
  tree.clear(); rmin.clear(); rmax.clear();
  for (size_type i=0; i < 600; ++i) {
    rmin.push_back(base_node(i/10000.,0.));
    rmax.push_back(base_node((i+1)/10000.,0.));
    tree.add_box(rmin.back(),rmax.back());
  }
  verify(rmin, rmax, tree);

  cout << "3D/2D random check\n";
  tree.clear(); rmin.clear(); rmax.clear();
  for (size_type i=0; i < 600; ++i) {
    rmin.push_back(base_node(gmm::random(double()), 0, gmm::random(double())));
    rmax.push_back(rmin.back() + base_node(.1+gmm::random(), 0.1, .1+gmm::random())/10.);
    tree.add_box(rmin.back(),rmax.back());
  }
  verify(rmin, rmax, tree);
  cout << "\nthe rtree is ok!\n";
}

int main(int argc, char **argv) {
  if (argc == 2 && strcmp(argv[1],"-quick")==0) quick = true;
  try {
    check_tree();
    /*if (!quick)
      speed_test(3,300000,20000);
      else speed_test(2,10000,100);*/
  } DAL_STANDARD_CATCH_ERROR;
  return 0;  
}
