#include <bgeot_geometric_trans.h>


int main(void)
{
  bgeot::pgeometric_trans pai;

  for (int i = 1; i < 4; ++i)
  {
    pai = bgeot::simplex_trans(i,i);

    cout << "Simplexe de dimension " << i << " et de degre " << i << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

    getchar(); cout << endl << endl;

  }

  {
    pai = bgeot::convex_product_trans(bgeot::simplex_trans(1,1), bgeot::simplex_trans(1,1));

    cout << "Produit de simplexes " << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

    getchar(); cout << endl << endl;

  }

  {
    pai = bgeot::linear_product_trans(bgeot::simplex_trans(1,1), bgeot::simplex_trans(1,1));

    cout << "produit lineaire de simplexes " << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

    getchar(); cout << endl << endl;

  }

  {
    pai = bgeot::linear_product_trans(bgeot::simplex_trans(2,1), bgeot::simplex_trans(1,1));

    cout << "produit lineaire de simplexes " << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

    getchar(); cout << endl << endl;

  }

{
    pai = bgeot::linear_product_trans(bgeot::linear_product_trans(bgeot::simplex_trans(1,1), bgeot::simplex_trans(1,1)), bgeot::simplex_trans(1,1));

    cout << "produit lineaire de simplexes " << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Poly " << k << " : " << pai->poly_vector()[k];
      cout << " point : " << pai->geometric_nodes()[k] << endl;

    }

    getchar(); cout << endl << endl;

  }

}
