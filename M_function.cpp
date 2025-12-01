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

double interpolate(double x, double y,const surface& m){
  // Vérifier que la triangulation n'est pas vide
  if (m.dt.number_of_vertices() < 3) {
    std::cerr << "Erreur: triangulation insuffisante" << std::endl;
    return 0.0;
  }

  Point query(x, y);

  // Vérifier si le point est dans la triangulation
  // Chercher la face contenant le point
  typename Delaunay::Face_handle face = m.dt.locate(query);

  // Si locate retourne nullptr ou si on est en dehors
  if (m.dt.is_infinite(face)) {
    //std::cerr << "Point (" << x << ", " << y << ") hors de la triangulation" << std::endl;
    return 0.0;
  }

  std::vector<std::pair<Point, double>> coords;
  double norm = CGAL::natural_neighbor_coordinates_2(
    m.dt, query, std::back_inserter(coords)).second;

  // Vérifier la validité de la normalisation
  if (norm <= 0.0 || coords.empty()) {
    std::cerr << "Coordonnées naturelles invalides pour (" << x << ", " << y << ")" << std::endl;
    std::cerr << "Norm: " << norm << ", Nb coords: " << coords.size() << std::endl;
    return 0.0;
  }

  //std::cout << "Interpolation en (" << x << ", " << y << ") avec " << coords.size() << " voisins, norm=" << norm << std::endl;

  return CGAL::linear_interpolation(
    coords.begin(), coords.end(), norm,
    CGAL::Data_access<std::map<Point, double>>(m.function_values));

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


