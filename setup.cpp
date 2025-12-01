#include "setup.h"

std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            row.push_back(std::stod(cell));
        }

        data.push_back(row);
    }

    file.close();
    return data;
}

void load_discharging_data(const std::string discharging_filename, std::vector<double>& all_currents, std::vector<double>& all_bks, std::vector<double>& all_vs, const float nominal_capacity, const float R_id, std::vector<double>& a_1, std::vector<double>& a_1_I, std::vector<double>& a_2, std::vector<double>& a_2_I, int S, int P){
    std::vector<std::vector<double>> discharging_data = readCSV(discharging_filename);

    std::vector<double> C_rates = {discharging_data[0][0]};
    std::vector<double> C_rates_start_index = {0};
    std::vector<double> C_rates_end_index;
    for(int i = 1; i < discharging_data.size(); i++){
        if ( discharging_data[i-1][0] != discharging_data[i][0]){
            C_rates.push_back(discharging_data[i][0]);
            C_rates_start_index.push_back(i);
            C_rates_end_index.push_back(i-1);
        }
    }
    C_rates_end_index.push_back(discharging_data.size() - 1);

    std::vector<double> I;
    std::vector<double> V;
    std::vector<double> E;

    for(int i = 0; i < C_rates.size(); i++){

        I.push_back(C_rates[i] * nominal_capacity * P * -1);
        V.push_back(discharging_data[C_rates_start_index[i]][2] * S);
        E.push_back(0);

        for(int j = C_rates_start_index[i]; j < C_rates_end_index[i]; j++){
            I.push_back( C_rates[i] * nominal_capacity * P * -1);
            V.push_back( discharging_data[j+1][2] * S);
            E.push_back( E.back() - ( discharging_data[j+1][1] - discharging_data[j][1] ) * discharging_data[j][2] * S * P * ( 1 + C_rates[i] * nominal_capacity * P * R_id / ( discharging_data[j][2] * S ) ) );
        }
    }

    double max_content = *std::max_element(E.begin(), E.end());

    for( int i = 0; i < C_rates.size(); i++){
        a_1.push_back(max_content - *std::max_element(&E[C_rates_start_index[i]], &E[C_rates_end_index[i]]));
        a_1_I.push_back( C_rates[i] * nominal_capacity * P * -1);
    }
/*
    for(int i = 0; i < E.size(); i++ ){
        E[i] = max_content - E[i];
    }
*/
    for(int i = 0; i < C_rates.size(); i++){
        a_2.push_back(*std::max_element(&E[C_rates_start_index[i]], &E[C_rates_end_index[i]]));
        a_2_I.push_back( C_rates[i] * nominal_capacity * P * -1);
    }
/*
    for(int i = 0; i < (a_1.size() / 2); i++){
        double temp = a_1[i];
        a_1[i] = a_1[a_1.size() - i - 1];
        a_1[a_1.size() - i - 1] = temp;
        temp = a_1_I[i];
        a_1_I[i] = a_1_I[a_1_I.size() - i - 1];
        a_1_I[a_1_I.size() - i - 1] = temp;
        temp = a_2[i];
        a_2[i] = a_2[a_2.size() - i - 1];
        a_2[a_2.size() - i - 1] = temp;
        temp = a_2_I[i];
        a_2_I[i] = a_2_I[a_2_I.size() - i - 1];
        a_2_I[a_2_I.size() - i - 1] = temp;
    }
*/

    all_currents.insert(all_currents.end(), I.begin(), I.end());
    all_vs.insert(all_vs.end(), V.begin(), V.end());
    all_bks.insert(all_bks.end(), E.begin(), E.end());
}

surface M_init(const std::string discharging_filename, const double nominal_capacity, const double R_id, std::vector<double>& a_1, std::vector<double>& a_1_I, std::vector<double>& a_2, std::vector<double>& a_2_I, int S, int P){
    
    // The data must be grouped by C_rate in the file and they must be sorted by voltage inside a group.
    
    std::vector<double> all_currents;
    std::vector<double> all_bks;
    std::vector<double> all_vs;

    load_discharging_data(discharging_filename, all_currents, all_bks, all_vs, nominal_capacity, R_id, a_1, a_1_I, a_2, a_2_I, S, P);

    return s_init(all_currents, all_bks, all_vs);
}

// Fonction pour exporter la surface vers un fichier CSV
void exportcurveToCSV(const std::string& filename, std::vector<double> x, std::vector<double> y){
  std::cout << "Début de l'export..." << std::endl;
    
  // Vérification des paramètres
  if (x.size() != y.size()) {
    std::cerr << "Erreur: |x| =/= |y|" << std::endl;
    return;
  }
    
  std::ofstream file(filename);
  
  if (!file.is_open()) {
    std::cerr << "Erreur: impossible d'ouvrir le fichier " << filename << std::endl;
    return;
  }
    
  // Écriture de l'en-tête
  file << "x,y" << std::endl;
  file.flush();  // Force l'écriture
  
  std::cout << "En-tête écrit, calcul des points..." << std::endl;
  
  // Itération sur la grille
  for (int i = 0; i < x.size(); ++i) {
    // Écriture de la ligne: x,y,z
    file << std::fixed << std::setprecision(6) 
         << x[i] << "," << y[i] << std::endl;
                
  }
    
  file.close();
  std::cout << "courbe exportée dans " << filename << " avec succès!" << std::endl;
  std::cout << "Points générés: " << x.size() << std::endl;
}


