#include <algorithm>
#include <cmath>
#include <cstdlib> // For std::getenv, std::stod
#include <filesystem>   // For file system operations (C++17 and later)
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers> // For std::numbers::pi
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h> // OpenMP for parallelization

#include <boost/math/quadrature/gauss_kronrod.hpp>

constexpr double PI_CONST = std::numbers::pi;
constexpr double SQRT_2PI = std::sqrt(2.0 * PI_CONST);

int get_int_from_env(const char* env_var_name, int default_value) {
    const char* env_val_str = std::getenv(env_var_name);
    if (env_val_str == nullptr) {
        std::cout << "INFO: Environment variable '" << env_var_name
                  << "' not set. Using default value: " << default_value << std::endl;
        return default_value;
    }
    try {
        // std::string s(env_val_str); // Create string for stoi
        // size_t processed_chars;
        // int value = std::stoi(s, &processed_chars);
        // if (processed_chars < s.length()) { // Check if the whole string was processed
        //     std::cerr << "WARNING: Trailing characters for environment variable '" << env_var_name
        //               << "'. Value: '" << env_val_str << "'. Using default value: "
        //               << default_value << std::endl;
        //     return default_value;
        // }
        // More robust: use strtol and check for errors and unprocessed characters
        char* end_ptr;
        long val_long = std::strtol(env_val_str, &end_ptr, 10);

        if (env_val_str == end_ptr || *end_ptr != '\0') { // No conversion or trailing chars
            std::cerr << "WARNING: Invalid integer format for environment variable '" << env_var_name
                      << "'. Value: '" << env_val_str << "'. Using default value: "
                      << default_value << std::endl;
            return default_value;
        }
        if (val_long > std::numeric_limits<int>::max() || val_long < std::numeric_limits<int>::min()) {
             std::cerr << "WARNING: Value out of range for int for environment variable '" << env_var_name
                      << "'. Value: '" << env_val_str << "'. Using default value: "
                      << default_value << std::endl;
            return default_value;
        }
        int value = static_cast<int>(val_long);
        std::cout << "INFO: Using value from environment variable '" << env_var_name
                  << "': " << value << std::endl;
        return value;
    } catch (const std::invalid_argument& ia) { // Should be caught by strtol checks
        std::cerr << "WARNING: Invalid argument for environment variable '" << env_var_name
                  << "'. Value: '" << env_val_str << "'. Not a valid integer. Using default value: "
                  << default_value << std::endl;
        return default_value;
    } catch (const std::out_of_range& oor) { // Should be caught by strtol checks
        std::cerr << "WARNING: Value out of range for environment variable '" << env_var_name
                  << "'. Value: '" << env_val_str << "'. Using default value: "
                  << default_value << std::endl;
        return default_value;
    }
}
// Retrieves a double value from an environment variable.
// If the variable is not set or invalid, returns a default value.
double get_double_from_env(const char* env_var_name, double default_value) {
    const char* env_val_str = std::getenv(env_var_name);
    if (env_val_str == nullptr) {
        std::cout << "INFO: Environment variable '" << env_var_name
                  << "' not set. Using default value: " << default_value << std::endl;
        return default_value;
    }
    try {
        double value = std::stod(std::string(env_val_str)); // Convert to string then to double
        std::cout << "INFO: Using value from environment variable '" << env_var_name
                  << "': " << value << std::endl;
        return value;
    } catch (const std::invalid_argument& ia) {
        std::cerr << "WARNING: Invalid argument for environment variable '" << env_var_name
                  << "'. Value: '" << env_val_str << "'. Not a valid double. Using default value: "
                  << default_value << std::endl;
        return default_value;
    } catch (const std::out_of_range& oor) {
        std::cerr << "WARNING: Value out of range for environment variable '" << env_var_name
                  << "'. Value: '" << env_val_str << "'. Using default value: "
                  << default_value << std::endl;
        return default_value;
    }
}

// Calculates the logarithm of the standard normal probability density function (PDF).
double log_standard_normal_pdf(double t) {
    double t_squared_half = 0.5 * t * t;
    if (t_squared_half > 700.0) { // exp(-700) is already very small, and exp(700) would overflow
        return -std::numeric_limits<double>::infinity();
    }
    return -0.5 * std::log(2.0 * PI_CONST) - t_squared_half;
}

// Calculates the standard normal probability density function (PDF) value.
double standard_normal_pdf_val(double t) {
    double log_pdf = log_standard_normal_pdf(t);
    if (std::isinf(log_pdf) && log_pdf < 0) {
        return 0.0;
    }
    return std::exp(log_pdf);
}

// Calculates the first term of function h in a numerically stable way.
double calculate_term1_stable(double x, double y, double q_val,
                              const std::function<double(double)>& sigma_func,
                              const std::function<double(double)>& mu_func) {
    if ((1.0 - q_val) <= std::numeric_limits<double>::epsilon()) {
        return 0.0;
    }
    double current_sigma = sigma_func(x);
    double current_mu = mu_func(x);

    if (current_sigma <= std::numeric_limits<double>::epsilon()) {
        // If sigma is very small, this is a Dirac-delta like peak at mu(x).
        // (1-q) * (1/sigma) * f_std((y-mu)/sigma)
        // For numerical stability, we work on a logarithmic scale as long as possible.
        // If y is close to mu, the value can be large, otherwise 0.
        // For more accurate handling, an analytical limit would be needed, or sigma could be shifted by a small epsilon.
        // We are now handling a simplified case: if y is very close to mu, we return a large value, otherwise 0.
        // But log_standard_normal_pdf already handles the t->inf case.
        double t_val = (y - current_mu) / (current_sigma + std::numeric_limits<double>::epsilon()); // Avoid division by 0
        double log_val = std::log(1.0 - q_val) - std::log(current_sigma + std::numeric_limits<double>::epsilon()) + log_standard_normal_pdf(t_val);
        
        if (log_val > std::log(std::numeric_limits<double>::max() / 2.0)) {
            return std::numeric_limits<double>::infinity();
        }
        if (log_val < std::log(std::numeric_limits<double>::min() * 2.0)) {
            return 0.0;
        }
        return std::exp(log_val);
    }

    double t_val = (y - current_mu) / current_sigma;
    double log_scaled_f_val = -std::log(current_sigma) + log_standard_normal_pdf(t_val);
    double log_term1_val = std::log(1.0 - q_val) + log_scaled_f_val;

    if (log_term1_val < std::log(std::numeric_limits<double>::min() * 2.0)) {
        return 0.0;
    } else if (log_term1_val > std::log(std::numeric_limits<double>::max() / 2.0)) {
        return std::numeric_limits<double>::infinity();
    } else {
        return std::exp(log_term1_val);
    }
}

// Calculates the second term of function h (where u is standard normal) in a numerically stable way.
double calculate_term2_stable_u_is_std_normal(double x, double y, double q_val,
                                            const std::function<double(double)>& sigma_func,
                                            const std::function<double(double)>& mu_func,
                                            int quadrature_points_param) { // quadrature_points_param is currently not used due to fixed QUADRATURE_POINTS_CONST
    if (q_val <= std::numeric_limits<double>::epsilon()) {
        return 0.0;
    }

    double current_sigma = sigma_func(x);
    double current_mu = mu_func(x);

    if (current_sigma <= std::numeric_limits<double>::epsilon()) {
        double log_val_term2 = std::log(q_val) + log_standard_normal_pdf(y - current_mu);
        if (log_val_term2 < std::log(std::numeric_limits<double>::min() * 2.0)) {
            return 0.0;
        }
        if (log_val_term2 > std::log(std::numeric_limits<double>::max() / 2.0)) {
            return std::numeric_limits<double>::infinity();
        }
        return std::exp(log_val_term2);
    }

    auto integrand_transformed_lambda = [&](double w_val) {
        double arg1 = y - current_mu - current_sigma * w_val;
        double log_f1 = log_standard_normal_pdf(arg1);
        double log_f2 = log_standard_normal_pdf(w_val);

        if ((std::isinf(log_f1) && log_f1 < 0) || (std::isinf(log_f2) && log_f2 < 0)) {
            return 0.0;
        }
        double log_product_integrand = log_f1 + log_f2;
        if (log_product_integrand < std::log(std::numeric_limits<double>::min() * 2.0)) {
            return 0.0;
        } else if (log_product_integrand > std::log(std::numeric_limits<double>::max() / 2.0)) {
            return std::numeric_limits<double>::infinity();
        } else {
            return std::exp(log_product_integrand);
        }
    };

    double transformed_lower_limit = -10.0; // Integration limits for w
    double transformed_upper_limit = 10.0;
    constexpr unsigned int QUADRATURE_POINTS_CONST = 15; // Number of points for Gauss-Kronrod
    
    double error; // Newer versions of Boost K-G integrator might require the error variable
    double integral_val = boost::math::quadrature::gauss_kronrod<double, QUADRATURE_POINTS_CONST>::integrate(
        integrand_transformed_lambda, transformed_lower_limit, transformed_upper_limit,
        5,    // max_depth_recursion - according to Boost quadrature docs, this can be a good value
        1e-9, &error); // tolerance and error output

    return q_val * integral_val;
}

// func_h: Core function combining term1 and term2.
double func_h_final_u_is_std_normal(double x, double y, double q_val,
                                  const std::function<double(double)>& sigma_func,
                                  const std::function<double(double)>& mu_func,
                                  int quadrature_points_param) {
    double term1 = calculate_term1_stable(x, y, q_val, sigma_func, mu_func);
    double term2 = calculate_term2_stable_u_is_std_normal(x, y, q_val, sigma_func, mu_func, quadrature_points_param);
    return term1 + term2;
}

// Calculates xtilde based on x and gamma boundaries.
double calculate_xtilde(double x, double gamma_minus_val, double gamma_plus_val) {
    if (x >= -gamma_minus_val && x <= gamma_plus_val) {
        return x;
    } else if (x > gamma_plus_val) {
        return gamma_plus_val;
    } else { // x is less than -gamma_minus_val
        return -gamma_minus_val;
    }
}

// Calculates the q function used in the iteration.
double calculate_q_func(double x0, double x_curr, double p_prob,
                        double gamma_minus_val, double gamma_plus_val,
                        double q_for_h_calc,
                        const std::function<double(double)>& sigma_func,
                        const std::function<double(double)>& mu_func,
                        int quadrature_points_h) {
    double xtilde_val = calculate_xtilde(x0, gamma_minus_val, gamma_plus_val);
    double h1 = func_h_final_u_is_std_normal(x0, x_curr - x0, q_for_h_calc, sigma_func, mu_func, quadrature_points_h);
    double h2 = func_h_final_u_is_std_normal(xtilde_val, x_curr - xtilde_val, q_for_h_calc, sigma_func, mu_func, quadrature_points_h);
    return h1 * (1.0 - p_prob) + h2 * p_prob;
}

// Integrates a function on a grid using the trapezoidal rule.
double integrate_on_grid(const std::vector<double>& x_grid_vals,
                         const std::vector<double>& integrand_values_on_grid) {
    if (x_grid_vals.size() != integrand_values_on_grid.size() || x_grid_vals.size() < 2) {
        return 0.0;
    }
    double integral_sum = 0.0;
    for (size_t i = 0; i < x_grid_vals.size() - 1; ++i) {
        double dx = x_grid_vals[i + 1] - x_grid_vals[i];
        integral_sum += 0.5 * (integrand_values_on_grid[i] + integrand_values_on_grid[i + 1]) * dx;
    }
    return integral_sum;
}

// Helper function to normalize a density function on the grid.
void normalize_density_on_grid(const std::vector<double>& x_grid_vals,
                               std::vector<double>& f_values_to_normalize) {
    // return; // Uncomment to disable normalization for debugging
    if (x_grid_vals.size() != f_values_to_normalize.size() || x_grid_vals.empty()) {
        std::cerr << "WARNING: Cannot normalize density. Grid and values mismatch or empty." << std::endl;
        return;
    }

    double integral_val = integrate_on_grid(x_grid_vals, f_values_to_normalize);

    if (std::abs(integral_val) < std::numeric_limits<double>::epsilon() * 100.0) { // *100 due to floating point errors
        std::cerr << "WARNING: Integral of the function is close to zero (" << integral_val 
                  << "). Cannot normalize. Values might all be zero or cancel out." << std::endl;
        // Optionally, all values could be set to 0 here, or left as is.
        // If all values are 0, the integral is 0, and normalization would be 0/0.
        // If the function is e.g. symmetric around 0 and the grid is too, but values are not all 0,
        // then the trapezoidal rule might yield 0. This case should be investigated more thoroughly.
        // For now, we don't change the values if the integral is 0.
        return;
    }

    for (size_t i = 0; i < f_values_to_normalize.size(); ++i) {
        f_values_to_normalize[i] /= integral_val;
    }
}

// Initializes f0 as a uniform distribution (window function) between -gamma_minus and +gamma_plus.
void initialize_f0_as_window(const std::vector<double>& x_grid_vals,
                            std::vector<double>& f0_values,
                            double neg_gamma_boundary, // This will be -gamma_minus_val
                            double pos_gamma_boundary) { // This will be gamma_plus_val
    f0_values.assign(x_grid_vals.size(), 0.0); // Initially zero everywhere

    if (pos_gamma_boundary <= neg_gamma_boundary) {
        std::cerr << "WARNING: Invalid gamma boundaries for initial window function (pos_gamma <= neg_gamma). "
                  << "Initializing f0 to all zeros." << std::endl;
        // Since f0_values is already initialized with zeros, there's nothing more to do here.
        // Normalization will later indicate an error if the integral is 0.
        return;
    }

    double interval_length = pos_gamma_boundary - neg_gamma_boundary;
    if (interval_length <= std::numeric_limits<double>::epsilon()) {
        std::cerr << "WARNING: Gamma interval length is zero or negative for initial window function. "
                  << "Initializing f0 to all zeros." << std::endl;
        return;
    }

    double function_value_in_window = 1.0 / interval_length;
    std::cout << "Initializing f0 as a window function between " << neg_gamma_boundary
              << " and " << pos_gamma_boundary
              << " with value " << function_value_in_window << std::endl;

    for (size_t i = 0; i < x_grid_vals.size(); ++i) {
        if (x_grid_vals[i] >= neg_gamma_boundary && x_grid_vals[i] <= pos_gamma_boundary) {
            f0_values[i] = function_value_in_window;
        }
        // Elsewhere it's already 0 (due to assign)
    }

    // Normalization is not strictly necessary here if it was analytically designed to integrate to 1,
    // BUT numerical integration on the grid might yield a different value.
    // Therefore, it's better to keep numerical normalization so it actually integrates to 1 on the grid.
    std::cout << "Normalizing initial f0 (window function)..." << std::endl;
    normalize_density_on_grid(x_grid_vals, f0_values); // Use the normalize_density_on_grid function
    
    double check_integral = integrate_on_grid(x_grid_vals, f0_values);
    std::cout << "Integral of initial f0 (window) after numerical normalization: " << check_integral << std::endl;
}

// Runs the iterative process to calculate f_k functions.
void run_iteration(int num_iterations_param,
                   const std::vector<double>& x_grid_vals,
                   const std::vector<double>& f_initial, // This is f0
                   std::vector<std::vector<double>>& f_history_output,
                   double p_prob, double gamma_minus_val, double gamma_plus_val,
                   double q_for_h_calc,
                   const std::function<double(double)>& sigma_func,
                   const std::function<double(double)>& mu_func,
                   int quadrature_points_h) {
    
    std::vector<double> f_current_iter = f_initial;
    f_history_output.push_back(f_current_iter); // Store f0

    int n_grid_points = x_grid_vals.size(); // Number of grid points

    // Start from iter_idx = 1 because f0 is already in f_history_output
    for (int iter_idx = 1; iter_idx < num_iterations_param; ++iter_idx) {
        std::cout << "Calculating iteration f_" << iter_idx << " (based on f_"
                  << (iter_idx - 1) << ") "
                  << "(OMP threads: " << omp_get_max_threads() << ")" << std::endl;
        
        std::vector<double> f_k_plus_1_next(n_grid_points); // Result of the next iteration

        // This loop is parallelized. Calculation of f_k_plus_1_next[i] is independent for each i.
        #pragma omp parallel
        {
            // Each thread needs its own copy of the integrand vector,
            // to avoid data races when passing it to integrate_on_grid.
            std::vector<double> integrand_for_f_k_plus_1_thread_local(n_grid_points);

            #pragma omp for schedule(dynamic) // Dynamic scheduling might be good if work per i varies
            for (int i = 0; i < n_grid_points; ++i) {
                double current_x = x_grid_vals[i];
                // The integrand: q(x0, current_x) * f_current_iter(x0)
                // This needs to be evaluated at every x0 grid point
                for (int j = 0; j < n_grid_points; ++j) {
                    double x0_val = x_grid_vals[j];
                    double q_val_at_x0_current_x = calculate_q_func(x0_val, current_x, p_prob,
                                                                    gamma_minus_val, gamma_plus_val,
                                                                    q_for_h_calc, sigma_func, mu_func,
                                                                    quadrature_points_h);
                    integrand_for_f_k_plus_1_thread_local[j] = q_val_at_x0_current_x * f_current_iter[j];
                }
                // Integration over x0 on the grid
                f_k_plus_1_next[i] = integrate_on_grid(x_grid_vals, integrand_for_f_k_plus_1_thread_local);
            }
        } // end omp parallel

        // NORMALIZATION for the newly calculated f_k_plus_1_next
        // std::cout << "Normalizing f_" << iter_idx << "..." << std::endl; // Optional log
        normalize_density_on_grid(x_grid_vals, f_k_plus_1_next);
        // double check_integral = integrate_on_grid(x_grid_vals, f_k_plus_1_next); // Optional log
        // std::cout << "Integral of f_" << iter_idx << " after normalization: " << check_integral << std::endl; // Optional log

        f_current_iter = f_k_plus_1_next; // Update for the next iteration

        // Writing to f_history_output is done on the main thread, after the parallel section
        f_history_output.push_back(f_current_iter);
    }
}

// Global parameters (consider passing them or using a struct if they grow)
double q_for_h_param;
double p_probability_param;
double gamma_plus_param;
double gamma_minus_param;
double const_sigma_param;
double const_mu_param;

// Can also be defined as a global constant if needed in multiple places
const double SECONDS_IN_A_DAY = 86400.0;

// Sigma function (constant, scaled for time unit).
double my_sigma_function(double x) {
    // The x parameter is not used here, as a constant value is returned
    (void)x; // To prevent compiler warnings about unused parameter

    // Scaling the read const_sigma_param
    // Assume const_sigma_param is a daily standard deviation, and we want a "base unit" standard deviation,
    // where the base unit is what the rest of the model (e.g., dt in SDE) is calibrated to.
    // If this is the sqrt(dt) factor, and dt=1 day, then the std dev is const_sigma_param itself.
    // If dt in the model is 1 second, then the daily std dev must be divided by sqrt(seconds_per_day).
    if (const_sigma_param < 0) {
        std::cerr << "WARNING: const_sigma_param (" << const_sigma_param
                  << ") is negative. Using its absolute value for scaling." << std::endl;
        return std::abs(const_sigma_param) / std::sqrt(SECONDS_IN_A_DAY);
    }
    if (SECONDS_IN_A_DAY <= 0) { // Safety check, though it's a constant here
        std::cerr << "ERROR: SECONDS_IN_A_DAY is not positive in my_sigma_function. Returning unscaled sigma." << std::endl;
        return const_sigma_param;
    }

    return const_sigma_param / std::sqrt(SECONDS_IN_A_DAY);
}

// Mu function (constant, scaled for time unit).
double my_mu_function(double x) {
    (void)x; // x is not used

    // Assume const_mu_param is a "daily" drift, and we need the drift per base unit of time.
    // If the base unit is 1 second, we divide by the number of seconds in a day.
    if (SECONDS_IN_A_DAY <= 0) { // Safety check
        std::cerr << "ERROR: SECONDS_IN_A_DAY is not positive in my_mu_function. Returning unscaled mu." << std::endl;
        return const_mu_param;
    }
    return const_mu_param / SECONDS_IN_A_DAY; // Linear scaling with time
}


int main() {
    // 1. Set parameters from environment variables or default values
    std::cout << "--- Reading parameters ---" << std::endl;
    q_for_h_param       = get_double_from_env("Q_FOR_H", 0.5);
    p_probability_param = get_double_from_env("P_PROBABILITY", 0.3);
    gamma_plus_param    = get_double_from_env("GAMMA_PLUS", 0.002); // Adjusted default
    gamma_minus_param   = get_double_from_env("GAMMA_MINUS", 0.002); // Adjusted default

    // NEW: Reading constant values for Sigma and Mu
    const_sigma_param   = get_double_from_env("CONST_SIGMA", 0.1); // Example: 10% daily volatility
    const_mu_param      = get_double_from_env("CONST_MU", 0.0);    // Example: 0 daily drift
    std::cout << "--------------------------" << std::endl;

    std::function<double(double)> sigma_f = my_sigma_function;
    std::function<double(double)> mu_f = my_mu_function;

    // The number of iterations determines how many f_k functions are generated.
    // f_history[0] will be f0, f_history[num_iterations-1] will be f_(num_iterations-1).

    // 2. Define grid (x_grid)
    std::vector<double> x_grid_values;
    double x_min   = get_double_from_env("X_MIN", -0.002);
    double x_max   = get_double_from_env("X_MAX", 0.002);
    int n_points_grid = get_int_from_env("N_POINTS_GRID", 3201);
    int num_iterations = get_int_from_env("NUM_ITERATIONS", 100);
    double dx_grid = (x_max - x_min) / (n_points_grid - 1);
    for (int i = 0; i < n_points_grid; ++i) {
        x_grid_values.push_back(x_min + i * dx_grid);
    }

    // 3. Initialize f_0(x)
    std::vector<double> f0_values;
    // Initialize as a window function between the negative of gamma_minus_param and gamma_plus_param
    // Note: gamma_minus_param is expected to be positive. The lower boundary is -gamma_minus_param.
    initialize_f0_as_window(x_grid_values, f0_values, -gamma_minus_param, gamma_plus_param);


    // 4. Container for the f_k series
    std::vector<std::vector<double>> f_history;

    // 5. Run iteration
    // f0 is already handled by initialize_f0_as_window and will be added to f_history at the start of run_iteration
    run_iteration(num_iterations, x_grid_values, f0_values, f_history,
                  p_probability_param, gamma_minus_param, gamma_plus_param,
                  q_for_h_param, sigma_f, mu_f, 15); // 15 for quadrature_points_h

    // 6. Process results: write to CSV file
    std::cout << "Iterations finished. " << f_history.size() << " f functions generated (from f0 to f"
              << (f_history.size() - 1) << ")." << std::endl;

    std::string csv_filename = "f_functions_output.csv";
    std::ofstream outfile(csv_filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Failed to open " << csv_filename << " for writing." << std::endl;
        return 1; // Indicate error
    }

    std::cout << "Writing CSV file: " << csv_filename << std::endl;

    // Write header
    outfile << "x_value";
    for (size_t k = 0; k < f_history.size(); ++k) {
        outfile << ",f" << k;
    }
    outfile << std::endl;

    // Write data rows
    if (!x_grid_values.empty() && !f_history.empty()) {
        if (f_history[0].size() != x_grid_values.size()) {
            std::cerr << "Error: Dimensions of x_grid and f_history[0] do not match!" << std::endl;
            outfile.close();
            return 1;
        }
    } else if (x_grid_values.empty() || f_history.empty()) {
        std::cerr << "Error: x_grid or f_history is empty, nothing to write." << std::endl;
        outfile.close();
        return 1;
    }
    
    outfile << std::fixed << std::setprecision(10); // Increased precision for CSV
    for (size_t i = 0; i < x_grid_values.size(); ++i) { // Iterate over x values (rows)
        outfile << x_grid_values[i];
        for (size_t k = 0; k < f_history.size(); ++k) { // Iterate over f_k functions (columns)
            outfile << "," << f_history[k][i];
        }
        outfile << std::endl;
    }
    outfile.close();
    std::cout << "CSV file successfully created: " << csv_filename << std::endl;
    
    const auto& x_grid = x_grid_values; // Alias for clarity
    const auto& last_fk_values = f_history.back(); // Values of the last function

    // Optional: Print some values of the last function to console for checking
    if (!f_history.empty()) {
        std::cout << "\nSome values of the last generated function (f" << f_history.size() - 1 << ") (for checking):" << std::endl;
        std::cout << std::fixed << std::setprecision(5);
        for (size_t i = 0; i < std::min(x_grid.size(), (size_t)5); ++i) { // Only the first 5 x values
            std::cout << "x = " << x_grid[i]
                      << ", f" << (f_history.size() - 1) << "(x) = " << last_fk_values[i] << std::endl;
        }
        std::cout.unsetf(std::ios_base::floatfield); // Reset precision printing for cout
    }

    // --- Integration analysis on the last function ---
    if (!f_history.empty()) {
        std::cout << "\n--- Integration analysis on the last function (f"
                  << f_history.size() - 1 << ") ---" << std::endl;
        std::cout << "Used gamma_minus: " << gamma_minus_param
                  << " (so the lower boundary for xtilde is: " << -gamma_minus_param << ")" << std::endl;
        std::cout << "Used gamma_plus: " << gamma_plus_param 
                  << " (so the upper boundary for xtilde is: " << gamma_plus_param << ")" << std::endl;

        // 1. Total integral on the grid
        double total_integral = integrate_on_grid(x_grid, last_fk_values);
        std::cout << "Total integral on the grid (" << x_grid.front() << " to " << x_grid.back()
                  << "): " << std::fixed << std::setprecision(10) << total_integral << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);

        if (std::abs(total_integral) < std::numeric_limits<double>::epsilon() * 100.0) {
            std::cout << "The total integral is close to zero, ratios cannot be meaningfully interpreted." << std::endl;
        } else {
            // 2. Integral on the segment (-inf, -gamma_minus_param] (within the grid limits)
            // Note: gamma_minus_param is stored as positive in the C++ code.
            double lower_tail_integral = 0.0;
            std::vector<double> x_segment_lower;
            std::vector<double> f_segment_lower;
            for (size_t i = 0; i < x_grid.size(); ++i) {
                if (x_grid[i] <= -gamma_minus_param) {
                    x_segment_lower.push_back(x_grid[i]);
                    f_segment_lower.push_back(last_fk_values[i]);
                } else {
                    // To integrate exactly up to -gamma_minus_param if it's not a grid point,
                    // the last trapezoid would need to be split or interpolated.
                    // Current simplification: only include grid points strictly <= -gamma_minus_param.
                    // If the first point is already > -gamma_minus_param, this segment will be empty or have one point.
                    if (!x_segment_lower.empty() && x_grid[i-1] < -gamma_minus_param && x_grid[i] > -gamma_minus_param) {
                        // Add the boundary point and interpolated f value for better accuracy
                        x_segment_lower.push_back(-gamma_minus_param);
                        double f_at_boundary = f_segment_lower.back() + (last_fk_values[i] - f_segment_lower.back()) * 
                                               (-gamma_minus_param - x_grid[i-1]) / (x_grid[i] - x_grid[i-1]);
                        f_segment_lower.push_back(f_at_boundary);
                    }
                    break; // The grid is sorted, no more relevant points
                }
            }
            if (x_segment_lower.size() >= 2) {
                lower_tail_integral = integrate_on_grid(x_segment_lower, f_segment_lower);
            } else if (x_segment_lower.size() == 1 && x_segment_lower.front() == x_grid.front()) {
                // If only the FIRST grid point is in the range, trapezoidal integral is 0.
                lower_tail_integral = 0.0; 
            }
            std::cout << "Integral from grid left (" << (x_segment_lower.empty() ? "N/A" : std::to_string(x_segment_lower.front()))
                      << ") to -gamma_minus (" << -gamma_minus_param << "): "
                      << std::fixed << std::setprecision(10) << lower_tail_integral << std::endl;
            std::cout.unsetf(std::ios_base::floatfield);

            // 3. Integral on the segment [gamma_plus_param, +inf) (within the grid limits)
            double upper_tail_integral = 0.0;
            std::vector<double> x_segment_upper;
            std::vector<double> f_segment_upper;
            bool first_upper_point_added = false;
            for (size_t i = 0; i < x_grid.size(); ++i) {
                if (x_grid[i] >= gamma_plus_param) {
                    if (!first_upper_point_added && i > 0 && x_grid[i-1] < gamma_plus_param) {
                        // Add the boundary point and interpolated f value for better accuracy
                        x_segment_upper.push_back(gamma_plus_param);
                        double f_at_boundary = last_fk_values[i-1] + (last_fk_values[i] - last_fk_values[i-1]) *
                                               (gamma_plus_param - x_grid[i-1]) / (x_grid[i] - x_grid[i-1]);
                        f_segment_upper.push_back(f_at_boundary);
                    }
                    x_segment_upper.push_back(x_grid[i]);
                    f_segment_upper.push_back(last_fk_values[i]);
                    first_upper_point_added = true;
                } else {
                    if (first_upper_point_added) { // Should not happen with sorted grid if logic is correct
                        break;
                    }
                }
            }
            if (x_segment_upper.size() >= 2) {
                upper_tail_integral = integrate_on_grid(x_segment_upper, f_segment_upper);
            } else if (x_segment_upper.size() == 1 && x_segment_upper.back() == x_grid.back()) {
                // If only the LAST grid point is in the range.
                upper_tail_integral = 0.0;
            }
            std::cout << "Integral from gamma_plus (" << gamma_plus_param << ") to grid right ("
                      << (x_segment_upper.empty() ? "N/A" : std::to_string(x_segment_upper.back()))
                      << "): " << std::fixed << std::setprecision(10) << upper_tail_integral << std::endl;
            std::cout.unsetf(std::ios_base::floatfield);

            // 4. Aggregated "tail" integral and ratio
            double tails_integral_sum = lower_tail_integral + upper_tail_integral;
            double tails_ratio = (total_integral == 0) ? 0 : (tails_integral_sum / total_integral);

            std::cout << "Aggregated 'tail' integral (left + right): " << std::fixed << std::setprecision(10) << tails_integral_sum << std::endl;
              std::cout << "Ratio of 'tail' areas to the total integral: "
                        << tails_ratio << " (i.e., " << tails_ratio * 100.0 << "%)" << std::endl;
              std::cout.unsetf(std::ios_base::floatfield);
  
              double tails_ratio_percentage = tails_ratio * 100.0;
              std::ostringstream filename_body_sstr;
              filename_body_sstr << std::fixed << std::setprecision(8); // Precision for filename parts
              filename_body_sstr << "tails_ratio_Q" << q_for_h_param
                                 << "_P" << p_probability_param
                                 << "_GP" << gamma_plus_param
                                 << "_GM" << gamma_minus_param
                                 << "_CS" << const_sigma_param
                                 << "_CM" << const_mu_param;
              std::string filename_body = filename_body_sstr.str();
              std::replace(filename_body.begin(), filename_body.end(), '.', 'p');
              std::string ratio_filename = filename_body + ".txt";

              // --- Step 1: Write and close the file ---
              std::ofstream outfile_ratio(ratio_filename);
              if (outfile_ratio.is_open()) {
                  outfile_ratio << std::fixed << std::setprecision(10);
                  outfile_ratio << tails_ratio_percentage << std::endl;
                  outfile_ratio.close(); // Close the file
                  std::cout << "Tails ratio successfully written to ./" << ratio_filename << std::endl; // Clarify current dir

                  // --- Step 2: Now attempt to move the closed file ---
                  try {
                      std::filesystem::path source_file_path = ratio_filename; // Path relative to CWD
                      std::filesystem::path target_directory_path = "tails_ratios";
      
                      if (!std::filesystem::exists(target_directory_path)) {
                          if (std::filesystem::create_directories(target_directory_path)) {
                              std::cout << "INFO: Created directory: " << target_directory_path.string() << std::endl;
                          } else {
                              std::cerr << "ERROR: Could not create directory: " << target_directory_path.string() << std::endl;
                              // Optionally, you might want to skip the move if directory creation fails
                          }
                      }
      
                      std::filesystem::path destination_file_path = target_directory_path / source_file_path.filename();
      
                      std::filesystem::rename(source_file_path, destination_file_path);
                      std::cout << "Successfully moved " << source_file_path.string()
                                << " to " << destination_file_path.string() << std::endl;
      
                  } catch (const std::filesystem::filesystem_error& e) {
                      std::cerr << "ERROR: Filesystem operation failed. Could not move file '" << ratio_filename
                                << "' to 'tails_ratios/' directory. Reason: " << e.what() << std::endl;
                      std::cerr << "Details: error code " << e.code() << ", path1: '" << e.path1().string() << "', path2: '" << e.path2().string() << "'" << std::endl;
                  }
              } else {
                  std::cerr << "Error: Failed to open " << ratio_filename << " for writing. File cannot be moved." << std::endl;
              }

            if (outfile_ratio.is_open()) {
                outfile_ratio << std::fixed << std::setprecision(10); // Precision for file content
                outfile_ratio << tails_ratio_percentage << std::endl;
                outfile_ratio.close();
                std::cout << "Tails ratio successfully written to " << ratio_filename << std::endl;
            } else {
                std::cerr << "Error: Failed to open " << ratio_filename << " for writing." << std::endl;
            }
        }
    }

    // --- NEW PART: Calculation of c1, c2, and the two extra terms ---
    if (!f_history.empty()) {
        double L_const = 1.0;
        double W_const = 1.0;
        double theta_const = 0.5;
        double p_val = p_probability_param; // Use the global p_probability

        double gamma_abs_boundary;
        if (std::abs(gamma_plus_param - gamma_minus_param) > 1e-9 * std::max(std::abs(gamma_plus_param), std::abs(gamma_minus_param))) {
            std::cout << "\nWARNING: gamma_plus_param (" << gamma_plus_param
                      << ") and gamma_minus_param (" << gamma_minus_param
                      << ") differ significantly. Using gamma_plus_param (" << gamma_plus_param 
                      << ") as the absolute gamma boundary for c1/c2 related term calculations."
                      << std::endl;
        }
        gamma_abs_boundary = gamma_plus_param; // Use gamma_plus_param as the representative |gamma|

        double c1 = std::pow(theta_const / (1.0 - theta_const), 1.0 - theta_const);
        double c2 = std::pow((1.0 - theta_const) / theta_const, theta_const);

        std::cout << "\n--- Calculation of new terms (c1, c2, etc.) ---" << std::endl;
        std::cout << "L = " << L_const << ", W = " << W_const << ", theta = " << theta_const
                  << ", p = " << p_val << ", gamma_abs_boundary = " << gamma_abs_boundary << std::endl;
        std::cout << "Calculated c1 = " << c1 << ", c2 = " << c2 << std::endl;

        auto integrand_upper_lambda = [&](double t_val) {
            if (std::abs(theta_const - 0.5) < 1e-9) { // Special case theta = 0.5
                return (std::exp(0.5 * t_val) + std::exp(gamma_abs_boundary - 0.5 * t_val) - 2.0 * std::exp(0.5 * gamma_abs_boundary));
            } else {
                double term_c1_part = c1 * (std::exp(t_val * (1.0 - theta_const)) - std::exp(gamma_abs_boundary * (1.0 - theta_const)));
                double term_c2_part = c2 * (std::exp(gamma_abs_boundary - t_val * theta_const) - std::exp(gamma_abs_boundary * (1.0 - theta_const))); // This seems to be gamma_abs * (theta_const) as per typical forms
                // Recheck formula for term_c2_part's second exp argument if issues arise. Assuming original was correct.
                // Original form from user had std::exp(gamma_abs * (1.0 - theta_const)) for the second part of c2 term as well. Let's stick to it.
                 return (term_c1_part + term_c2_part);
            }
        };

        auto integrand_lower_lambda = [&](double t_val) { // t_val here is x, which is negative in this domain
            if (std::abs(theta_const - 0.5) < 1e-9) { // Special case theta = 0.5
                 // For t_val < 0, often form is exp(-0.5*t) + exp(0.5*t - gamma_abs) - 2*exp(-0.5*gamma_abs)
                 // Original: (std::exp(0.5 * t_val) + std::exp(-0.5 * t_val - gamma_abs) - 2.0 * std::exp(-0.5 * gamma_abs));
                 // This assumes t_val is negative, and gamma_abs is positive boundary.
                return (std::exp(0.5 * t_val) + std::exp(-gamma_abs_boundary - 0.5 * t_val) - 2.0 * std::exp(-0.5 * gamma_abs_boundary));

            } else {
                // Original:
                // double exp_gamma_theta_minus_1 = std::exp(gamma_abs_boundary * (theta_const - 1.0));
                // double term_c1_part = c1 * (std::exp(t_val * (1.0 - theta_const)) - exp_gamma_theta_minus_1);
                // double term_c2_part = c2 * (std::exp(-t_val * theta_const - gamma_abs_boundary) - exp_gamma_theta_minus_1);
                // return (term_c1_part + term_c2_part);
                // This form implies boundary at -gamma_abs.
                // Let's assume the form: c1(e^{(1-th)t} - e^{-(1-th)gamma}) + c2(e^{-th*t - gamma} - e^{-(1-th)gamma})
                // This structure seems to be for t < -gamma.
                // The provided code has:
                 double exp_neg_gamma_one_minus_theta = std::exp(-gamma_abs_boundary * (1.0 - theta_const)); // Equivalent to exp(gamma * (theta-1))
                 double term_c1_part = c1 * (std::exp(t_val * (1.0 - theta_const)) - exp_neg_gamma_one_minus_theta);
                 double term_c2_part = c2 * (std::exp(-t_val * theta_const - gamma_abs_boundary * theta_const) - exp_neg_gamma_one_minus_theta); // Modified c2 term, using -gamma*theta for second part of exponent
                return (term_c1_part + term_c2_part);
            }
        };
        
        double upper_tail_term_value = 0.0;
        std::vector<double> x_segment_upper_calc;
        std::vector<double> integrand_times_f_upper;

        for (size_t i = 0; i < x_grid.size(); ++i) {
            if (x_grid[i] >= gamma_abs_boundary) {
                if (x_segment_upper_calc.empty() && i > 0 && x_grid[i-1] < gamma_abs_boundary && std::abs(x_grid[i] - x_grid[i-1]) > 1e-9) {
                    // Interpolate f at gamma_abs_boundary for better accuracy at the start of the segment
                    double f_star_at_gamma_abs = last_fk_values[i-1] + (last_fk_values[i] - last_fk_values[i-1]) * 
                                                 (gamma_abs_boundary - x_grid[i-1]) / (x_grid[i] - x_grid[i-1]);
                    x_segment_upper_calc.push_back(gamma_abs_boundary);
                    integrand_times_f_upper.push_back(integrand_upper_lambda(gamma_abs_boundary) * f_star_at_gamma_abs);
                }
                x_segment_upper_calc.push_back(x_grid[i]);
                integrand_times_f_upper.push_back(integrand_upper_lambda(x_grid[i]) * last_fk_values[i]);
            }
        }
        if (x_segment_upper_calc.size() >= 2) {
            upper_tail_term_value = integrate_on_grid(x_segment_upper_calc, integrand_times_f_upper);
        }
        upper_tail_term_value *= p_val; // The p_val multiplier is applied to the integral only at the end
        std::cout << "Upper tail term (p * integral from " << gamma_abs_boundary << " to " << x_grid.back() << "): "
                  << std::fixed << std::setprecision(10) << upper_tail_term_value << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);

        double lower_tail_term_value = 0.0;
        std::vector<double> x_segment_lower_calc;
        std::vector<double> integrand_times_f_lower;
        double neg_gamma_abs_boundary = -gamma_abs_boundary;

        for (size_t i = 0; i < x_grid.size(); ++i) {
            if (x_grid[i] <= neg_gamma_abs_boundary) {
                x_segment_lower_calc.push_back(x_grid[i]);
                integrand_times_f_lower.push_back(integrand_lower_lambda(x_grid[i]) * last_fk_values[i]);
            } else { // x_grid[i] > neg_gamma_abs_boundary
                if (!x_segment_lower_calc.empty() && i > 0 && x_grid[i-1] < neg_gamma_abs_boundary && std::abs(x_grid[i] - x_grid[i-1]) > 1e-9) {
                     // Interpolate f at neg_gamma_abs_boundary for better accuracy at the end of the segment
                     double f_star_at_neg_gamma_abs = last_fk_values[i-1] + (last_fk_values[i] - last_fk_values[i-1]) * 
                                                      (neg_gamma_abs_boundary - x_grid[i-1]) / (x_grid[i] - x_grid[i-1]);
                     x_segment_lower_calc.push_back(neg_gamma_abs_boundary);
                     integrand_times_f_lower.push_back(integrand_lower_lambda(neg_gamma_abs_boundary) * f_star_at_neg_gamma_abs);
                }
                break; // Past the relevant segment
            }
        }
        // integrate_on_grid assumes increasing x order. x_segment_lower_calc is already sorted.
        if (x_segment_lower_calc.size() >= 2) {
            lower_tail_term_value = integrate_on_grid(x_segment_lower_calc, integrand_times_f_lower);
        }
        lower_tail_term_value *= p_val; // The p_val multiplier is applied to the integral only at the end
        std::cout << "Lower tail term (p * integral from " << x_grid.front() << " to " << neg_gamma_abs_boundary << "): "
                  << std::fixed << std::setprecision(10) << lower_tail_term_value << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);
        
        // The factor L_const * std::pow(W_const, theta_const) simplifies to 1 since L=1, W=1.
        double total_extra_term = upper_tail_term_value + lower_tail_term_value;

        std::cout << "Total extra term value (L*W^theta * (upper_term + lower_term), with L*W^theta=1): "
                  << std::fixed << std::setprecision(10) << total_extra_term << std::endl;
        std::cout.unsetf(std::ios_base::floatfield);

        // Writing extra terms to file, using the same filename convention as tails_ratio
        std::ostringstream extra_filename_body_sstr;
        extra_filename_body_sstr << std::fixed << std::setprecision(8); // Precision for filename parts
        extra_filename_body_sstr << "extra_calc_Q" << q_for_h_param
                                 << "_P" << p_probability_param
                                 << "_GP" << gamma_plus_param
                                 << "_GM" << gamma_minus_param
                                 << "_CS" << const_sigma_param
                                 << "_CM" << const_mu_param
                                 << "_Th" << theta_const; // Include Theta as well
        std::string extra_filename_body = extra_filename_body_sstr.str();
        std::replace(extra_filename_body.begin(), extra_filename_body.end(), '.', 'p');
        std::string extra_term_filename = extra_filename_body + ".txt";

        std::ofstream outfile_extra_terms(extra_term_filename);
        if (outfile_extra_terms.is_open()) {
            outfile_extra_terms << std::fixed << std::setprecision(10); // Precision for file content
            outfile_extra_terms << "c1: " << c1 << std::endl;
            outfile_extra_terms << "c2: " << c2 << std::endl;
            outfile_extra_terms << "upper_tail_term_calc: " << upper_tail_term_value << std::endl;
            outfile_extra_terms << "lower_tail_term_calc: " << lower_tail_term_value << std::endl;
            outfile_extra_terms << "total_extra_term_calc: " << total_extra_term << std::endl;
            outfile_extra_terms.close();
            std::cout << "Extra calculated terms successfully written to " << extra_term_filename << std::endl;
        } else {
            std::cerr << "Error: Could not open " << extra_term_filename << " for writing." << std::endl;
        }
    }
    return 0;
}
