#include <bgeot_geometric_trans.h>

using bgeot::size_type;

int main(void)
{
  try {
    

    bgeot::pgeometric_trans p1 =
      bgeot::geometric_trans_descriptor("GT_QK(2,2)");
    bgeot::pgeometric_trans p2 =
      bgeot::geometric_trans_descriptor("GT_PRODUCT(GT_PK(1,2),  GT_PK(1,2))");

    cout << "p1 = " << p1 << " p2 = " << p2 << endl;
    if (p1 != p2) DAL_THROW(dal::internal_error, "internal_error");

    cout << "Name of geometric trans : "
	 << bgeot::name_of_geometric_trans(p1) << endl;

    cout << "Name of the product of the two : " 
	 << bgeot::name_of_geometric_trans(product_geotrans(p1, p2)) << endl;


    bgeot::pgeometric_trans pai;

    char meth[500];

    for (int i = 1; i < 4; ++i)
      {
	pai = bgeot::simplex_geotrans(i,i);
	
	cout << "Simplexe de dimension " << i << " et de degre " << i << endl;
	
	for (size_type k = 0; k < pai->nb_points(); ++k)
	  {
	    cout << "Poly " << k << " : " << pai->poly_vector()[k];
	    cout << " point : " << pai->geometric_nodes()[k] << endl;
	    
	  }
	
	cout << endl << endl;
	
      }
    
    {
      sprintf(meth, "GT_PRODUCT(GT_PK(1,1), GT_PK(1,1))");
      pai = bgeot::geometric_trans_descriptor(meth);
      
      cout << "Produit de simplexes " << endl;
      
      for (size_type k = 0; k < pai->nb_points(); ++k)
	{
	  cout << "Poly " << k << " : " << pai->poly_vector()[k];
	  cout << " point : " << pai->geometric_nodes()[k] << endl;
	  
	}
      
      cout << endl << endl;
      
    }
    
    {
      sprintf(meth, "GT_LINEAR_PRODUCT(GT_PK(1,1), GT_PK(1,1))");
      pai = bgeot::geometric_trans_descriptor(meth);
      
      cout << "produit lineaire de simplexes " << endl;
      
      for (size_type k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

      cout << endl << endl;
      
    }
    
    {
      sprintf(meth, "GT_LINEAR_PRODUCT(GT_PK(2,1), GT_PK(1,1))");
      pai = bgeot::geometric_trans_descriptor(meth);
      
      cout << "produit lineaire de simplexes " << endl;
      
      for (size_type k = 0; k < pai->nb_points(); ++k)
	{
	  cout << "Poly " << k << " : " << pai->poly_vector()[k];
	  cout << " point : " << pai->geometric_nodes()[k] << endl;
	  
	}
      
      cout << endl << endl;
      
    }
    
    {
      sprintf(meth, 
   "GT_LINEAR_PRODUCT(GT_LINEAR_PRODUCT(GT_PK(1,1), GT_PK(1,1)), GT_PK(1,1))");
      pai = bgeot::geometric_trans_descriptor(meth);
      
      cout << "produit lineaire de simplexes " << endl;
      
      for (size_type k = 0; k < pai->nb_points(); ++k)
	{
	  cout << "Poly " << k << " : " << pai->poly_vector()[k];
	  cout << " point : " << pai->geometric_nodes()[k] << endl;
	  
	}
      
      cout << endl << endl;
      
    }
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
