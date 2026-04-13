#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]

List cpp_model(IntegerVector population, NumericVector beta,
               IntegerVector init_inf, int time_max, 
               double period_infectious, double period_latent,
               double connection_strength) {

  ////////////////////////////////////////////////
  // Define objects
  int n_patch = population.size();
  
  IntegerMatrix sus(n_patch, time_max); // Susceptible
  IntegerMatrix exp(n_patch, time_max); // Exposed
  IntegerMatrix inf(n_patch, time_max); // Infectious
  IntegerMatrix rec(n_patch, time_max); // Recovered
  IntegerMatrix new_recoveries(n_patch, time_max);  // New recoveries
  for (int i = 0; i < n_patch; ++i) {
    for (int j = 0; j < time_max; ++j) {
      new_recoveries(i, j) = 0;  
    }
  }
  IntegerMatrix exp_s(n_patch, time_max); // Start of exposure (newly exposed)
  IntegerMatrix inf_s(n_patch, time_max); // Start of infectiousness (newly infectious)
  IntegerVector num(n_patch); // Number of infectious individuals
  NumericMatrix movement(n_patch, n_patch); // Air movement
  NumericMatrix beta_matrix(n_patch, n_patch); // Beta matrix
  unordered_map<string, vector<int>> data; // Dataset of each unique ID
  vector<vector<vector<int>>> inf_ids; // [patch][time][unique ID]
  inf_ids.resize(n_patch);
  for (int patch = 0; patch < n_patch; ++patch) {
    inf_ids[patch].resize(time_max);
  }
  
  ////////////////////////////////////////////////
  // Set initial conditions
  for (int i = 0; i < n_patch; ++i) {
    sus(i, 0) = population[i] - init_inf[i];
    std::cout << "sus[" << i <<  "] = " << sus(i, 0) << std::endl;
  }
  inf(_, 0) = init_inf;
  rec(_, 0) = rep(0, n_patch);

  int curr_time = 0;
  
  for (int patch = 0; patch < n_patch; ++patch) {
    num[patch] = init_inf[patch];
    
    for (int id = 0; id < init_inf[patch]; ++id) {
      
      //int inf_period = sample(inf_period_dist.size(), 1, false, inf_period_dist, true)[0];
      int inf_period = max(int(round(rgamma(1, period_infectious*1.5, 1/1.5), 0)[0]), 1);
      int stop = curr_time + inf_period;
      int hard_stop = min(stop, time_max);
      
      for (int step = curr_time; step < hard_stop; ++step) {
        inf_ids[patch][step].push_back(id);
      }
      
      string unique_id = to_string(patch) + "_" + to_string(id);
      data[unique_id] = {-1, -1, -1, -1, curr_time,
                         inf_period, -1, stop, int(hard_stop != stop)};
    }
  }
  
  ////////////////////////////////////////////////
  // Update loop
  NumericMatrix lambda(n_patch, n_patch);
  NumericVector proba_patch(n_patch);
  IntegerVector pop(n_patch);
  NumericMatrix foi(n_patch, n_patch);

  int duration = 0;
  
  for (int curr_time = 1; curr_time < time_max; ++curr_time) { // 

    // Load movement matrix directly in the loop
    std::string filename = "./matrices1/matrix_" + std::to_string(curr_time + 0) + ".txt"; // 2 Jan 2009
    std::ifstream file(filename);
    if (!file.is_open()) {
      stop("Could not open file: " + filename);
    }
    
    std::vector<std::vector<double>> matrixData;
    std::string line;
    
    while (std::getline(file, line)) {
      std::stringstream lineStream(line);
      std::string cell;
      std::vector<double> rowData;
      
      while (lineStream >> cell) {
        rowData.push_back(std::stod(cell));
      }
      matrixData.push_back(rowData);
    }
    
    // Convert the NumericMatrix
    int rows = matrixData.size();
    int cols = rows > 0 ? matrixData[0].size() : 0;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        movement(i, j) = matrixData[i][j];
      }
    }

    if (curr_time == 1) {
      std::cout << "movement matrix [1, 1]: " << movement(1-1, 1-1) <<  std::endl;
      std::cout << "movement matrix [1, 2]: " << movement(1-1, 2-1) <<  std::endl;
      std::cout << "movement matrix [1, 3]: " << movement(1-1, 3-1) <<  std::endl;
      std::cout << "movement matrix [2, 1]: " << movement(2-1, 1-1) <<  std::endl;
      std::cout << "movement matrix [3, 1]: " << movement(3-1, 1-1) <<  std::endl;
    }


    if (curr_time >= 114 && curr_time <= 127 ) { // from 25 Apr 2009 to 8 May 2009 (116)
      for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (j != i) {
          beta_matrix(i,j) = beta[j] * connection_strength * movement(i,j) * 0.7; 
        } else if (j == i){
          beta_matrix(i,j) = beta[i] * 0.7;
        }      
      }
    }
    } else if (curr_time >= 128) { // after 9 May 2009 (117)
      for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (j != i) {
          beta_matrix(i,j) = beta[j] * connection_strength * movement(i,j) * 1; 
        } else if (j == i){
          beta_matrix(i,j) = beta[i] * 1;
        }      
      }
    }
    } else { // before 24 Apr 2009 
      for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (j != i) {
          beta_matrix(i,j) = beta[j] * connection_strength * movement(i,j); 
        } else if (j == i){
          beta_matrix(i,j) = beta[i];
        }      
      }
    }
    }

/*
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (j != i) {
          beta_matrix(i,j) = beta[j] * connection_strength * movement(i,j); 
        } else if (j == i){
          beta_matrix(i,j) = beta[i];
        }      
      }
    }
  */

    // Total population vector
    pop = sus(_, curr_time - 1) + exp(_, curr_time - 1) + inf(_, curr_time - 1) + rec(_, curr_time - 1);
    // int total_population = sum(sus(_, curr_time-1)) + sum(exp(_, curr_time-1)) + sum(inf(_, curr_time-1)) + sum(rec(_, curr_time-1));
    std::cout << "Time: " << curr_time << " || total_population: " << pop << endl; // DEBUG

    // Compute lambda
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        // if (pop[j] == 0) {
        //   lambda(i,j) = 0;
        // } else {
          lambda(i, j) = beta_matrix(i, j) * sus(j, curr_time - 1) * inf(i, curr_time - 1); // / pop[j];
        // Rcpp::Rcout << "lambda(" << i << "," << j << ") = " << lambda(i, j) << std::endl;
        
        // }
      }
    }
    
    if (curr_time == 1) {
      foi = clone(lambda);
    }
    
    // Update
    for (int patch = 0; patch < n_patch; ++patch) {
      
      double sum_lambda = sum(lambda(_, patch));
      int s_to_e = min(int(round(rpois(1, sum_lambda), 0)[0]), sus(patch, curr_time - 1));
      // Rcout << "Time: " << curr_time << " || Patch: " << patch << " || s_to_e: " << s_to_e << endl; // DEBUG
      
      proba_patch = lambda(_, patch);
      if (sum_lambda > 0.0) {
        proba_patch = proba_patch / sum_lambda;
      }
      
      for (int l = 0; l < n_patch; ++l) {
        if (proba_patch[l] < 0) {
          IntegerVector S = sus(_, curr_time - 1);
          IntegerVector E = exp(_, curr_time - 1);
          IntegerVector I = inf(_, curr_time - 1);
          IntegerVector R = rec(_, curr_time - 1);
          Rcout << "Time = " << curr_time << " || " << S << " || " << I << " || "<< E << " || " << R << endl;
          stop("Negative lambda!");
        } 
      }

      // Loop on new exposed in patch patch and at time curr_time
      for (int id = num[patch]; id < num[patch] + s_to_e; ++id) {
        // Rcout << "id: " << id;
        
        // ID
        string unique_id = to_string(patch) + "_" + to_string(id);
        
        // Times of events
        int lat_period = max(int(round(rgamma(1, period_latent*2.5, 1/2.5), 0)[0]), 1);  // Latent period
        int inf_period = max(int(round(rgamma(1, period_infectious*1.5, 1/1.5), 0)[0]), 1); // int inf_period = sample(inf_period_dist.size(), 1, false, inf_period_dist, true)[0];

        // Source 
        int patch_source = sample(n_patch, 1, false, proba_patch, false)[0];
        vector<int> ids = inf_ids[patch_source][curr_time - 1];
        IntegerVector vec(ids.begin(), ids.end());
        int id_source = sample(vec, 1)[0];
        
        //  Generation time
        string unique_id_source = to_string(patch_source) + "_" + to_string(id_source);
        int ser_interval_time = curr_time + lat_period - data[unique_id_source][4];

        // Parameters to updata data
        int stop = curr_time + lat_period + inf_period;
        if (curr_time + lat_period < time_max) inf_s(patch, curr_time + lat_period) += 1;

        int hard_stop = min(stop, time_max);
        if (curr_time + lat_period + inf_period < time_max) new_recoveries(patch, curr_time + lat_period + inf_period) +=  1;
        for (int step = curr_time + lat_period; step < hard_stop; ++step) {
          inf_ids[patch][step].push_back(id);
        }
          
        data[unique_id] = {patch_source, id_source, curr_time, lat_period,
                             curr_time + lat_period, inf_period,
                             ser_interval_time, stop,
                             int(hard_stop != stop)};
      }
        
      // Counter of exposed
      num[patch] += s_to_e;
      
      //  Update counts
      int e_to_i = inf_s(patch, curr_time);
      sus(patch, curr_time) = sus(patch, curr_time - 1) - s_to_e;
      exp(patch, curr_time) = exp(patch, curr_time - 1) + s_to_e - e_to_i;
      exp_s(patch, curr_time) = s_to_e;

      inf(patch, curr_time) = inf_ids[patch][curr_time].size();
      rec(patch, curr_time) = rec(patch, curr_time - 1) + new_recoveries(patch, curr_time);

      // Update duration
      if (sum(inf(_, curr_time)) != 0 && sum(exp(_, curr_time)) != 0) duration = curr_time + 1;
      
    }
  }

  //////////////////////////////////////////////
  // Outputs
  List dynamics = List::create(Named("sus") = sus,
                               Named("exp") = exp,
                               Named("inf") = inf,
                               Named("rec") = rec,
                               Named("inf_s") = inf_s,
                               Named("exp_s") = exp_s);

  List sourcing = List::create(Named("data") = data,
                               Named("inf_ids") = inf_ids);
  int epi_size = sum(inf_s);

  return List::create(Named("dynamics") = dynamics,
                      Named("sourcing") = sourcing,
                      Named("size") = epi_size,
                      Named("foi") = foi,
                      Named("beta_mat") = beta_matrix,
                      Named("duration") = duration);

}
