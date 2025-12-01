#include"M_function.h"

surface s_init(std::vector<double> &all_currents, std::vector<double> &all_bks, std::vector<double> &all_vs){
  surface m;
  // Insérer les points
  for (size_t i = 0; i < all_currents.size(); i++) {
     Point p(all_currents[i], all_bks[i]);
     m.dt.insert(p);
     m.function_values[p] = all_vs[i];
  }

  return m;
}

static inline void barycentric_coords(
    const Point& A, const Point& B, const Point& C,
    const Point& P,
    double &u, double &v, double &w)
{
    double x = P.x(), y = P.y();
    double x1 = A.x(), y1 = A.y();
    double x2 = B.x(), y2 = B.y();
    double x3 = C.x(), y3 = C.y();

    double detT = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);

    u = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / detT;
    v = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / detT;
    w = 1.0 - u - v;
}

double interpolate(double x, double y, const surface& m) {
    Point query(x, y);
    auto fh = m.dt.locate(query);

    if (m.dt.is_infinite(fh)) return 0.0;

    // Les 3 sommets
    const Point& p0 = fh->vertex(0)->point();
    const Point& p1 = fh->vertex(1)->point();
    const Point& p2 = fh->vertex(2)->point();

    // Valeurs
    double f0 = m.function_values.at(p0);
    double f1 = m.function_values.at(p1);
    double f2 = m.function_values.at(p2);

    // Coordonnées barycentriques
    double l0, l1, l2;
    barycentric_coords(p0, p1, p2, query, l0, l1, l2);

    return l0*f0 + l1*f1 + l2*f2;
}


// Fonction pour exporter la surface vers un fichier CSV
void exportSurfaceToCSV(const std::string& filename, 
  double x_min, double x_max, int x_steps,
  double y_min, double y_max, int y_steps,const surface& M) {
  std::cout << "Début de l'export..." << std::endl;
    
  // Vérification des paramètres
  if (x_steps <= 0 || y_steps <= 0) {
    std::cerr << "Erreur: nombre de pas invalide" << std::endl;
    return;
  }
    
  if (x_steps == 1) x_steps = 2;  // Éviter division par zéro
  if (y_steps == 1) y_steps = 2;
    
  std::ofstream file(filename);
  
  if (!file.is_open()) {
    std::cerr << "Erreur: impossible d'ouvrir le fichier " << filename << std::endl;
    return;
  }
    
  // Écriture de l'en-tête
  file << "x,y,z" << std::endl;
  file.flush();  // Force l'écriture
  
  std::cout << "En-tête écrit, calcul des points..." << std::endl;
  
  // Calcul des pas
  double x_step = (x_max - x_min) / (x_steps - 1);
  double y_step = (y_max - y_min) / (y_steps - 1);
  
  int count = 0;
  
  // Itération sur la grille
  for (int i = 0; i < x_steps; ++i) {
    double x = x_min + i * x_step;
    
    for (int j = 0; j < y_steps; ++j) {
      double y = y_min + j * y_step;
            
      try {
        //std::cout << "calcul de z en " << x << " " << y << std::endl;
        double z = interpolate(x, y,M);
                
        // Vérifier que z est un nombre valide
        if (std::isnan(z) || std::isinf(z)) {
          std::cerr << "Valeur invalide à (" << x << "," << y << ")" << std::endl;
          z = 0.0;
        }
                
        // Écriture de la ligne: x,y,z
        file << std::fixed << std::setprecision(6) 
             << x << "," << y << "," << z << std::endl;
                
        count++;
                
        // Afficher la progression tous les 1000 points
        if (count % 1000 == 0) {
          std::cout << "Points écrits: " << count << std::endl;
        }
                
      } catch (const std::exception& e) {
        std::cerr << "Exception lors du calcul de M(" << x << "," << y << "): " << e.what() << std::endl;
      }
    }
  }
    
  file.close();
  std::cout << "Surface exportée dans " << filename << " avec succès!" << std::endl;
  std::cout << "Points générés: " << count << std::endl;
}
