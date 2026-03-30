#include "mcfost_interface.hpp"
#include <iostream>
#include <random>

/**
 * @brief Example usage of the MCFOST C++ interface
 */

void exampleUsage() {
    try {
        // Create MCFOST interface
        MCFOST::Interface mcfost;

        // Configure MCFOST
        MCFOST::MCFOSTConfig config;
        config.parameter_file = "disk.para";
        config.use_sph_limits_file = false;
        config.fix_star = false;
        config.turn_on_lacc = true;
        config.use_ism_heating = false;
        config.compute_frad = true;
        config.write_t_files = false;
        config.ism_heating_mode = 0;

        // Number of dust types
        int ndusttypes = 2;

        // Initialize MCFOST
        std::cout << "Initializing MCFOST..." << std::endl;
        mcfost.initialize(config, ndusttypes);
        std::cout << "MCFOST initialized successfully!" << std::endl;

        // Create example particle data
        int np = 10000;  // number of particles
        int nptmass = 1; // one star
        int ntypes = 2;  // gas + dust

        auto particles = MCFOST::createParticleData(np, ndusttypes, nptmass);

        // Set up example particle data
        setupExampleParticles(particles, np, ndusttypes, nptmass, ntypes);

        // Run MCFOST calculation
        std::cout << "Running MCFOST calculation..." << std::endl;
        auto results = mcfost.run(particles, config);
        std::cout << "MCFOST calculation completed!" << std::endl;

        // Process results
        processResults(results, np);

        // Example of getting diffusion opacity
        double temp = 100.0; // K
        int icell = 1;
        double kappa = mcfost.getDiffusionOpacity(temp, icell);
        std::cout << "Diffusion opacity at " << temp << " K: " << kappa << " cm²/g" << std::endl;

        // Reset for next calculation
        mcfost.reset();

    } catch (const MCFOST::MCFOSTException& e) {
        std::cerr << "MCFOST Error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void setupExampleParticles(MCFOST::ParticleData& particles, int np, int ndusttypes, int nptmass, int ntypes) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(-10.0, 10.0); // AU
    std::uniform_real_distribution<> vel_dist(-1.0, 1.0);   // km/s
    std::uniform_real_distribution<> temp_dist(10.0, 1000.0); // K
    std::uniform_real_distribution<> dust_dist(0.0, 0.1);    // dust fraction

    // Set up particle positions, velocities, and properties
    for (int i = 0; i < np; ++i) {
        // Positions (x, y, z) and smoothing length h
        particles.xyzh[0][i] = pos_dist(gen); // x
        particles.xyzh[1][i] = pos_dist(gen); // y
        particles.xyzh[2][i] = pos_dist(gen); // z
        particles.xyzh[3][i] = 0.1;           // smoothing length in AU

        // Velocities (vx, vy, vz) and internal energy u
        particles.vxyzu[0][i] = vel_dist(gen); // vx
        particles.vxyzu[1][i] = vel_dist(gen); // vy
        particles.vxyzu[2][i] = vel_dist(gen); // vz
        particles.vxyzu[3][i] = 1e12;         // internal energy (erg/g)

        // Phase information (0 = gas, 1 = dust, etc.)
        particles.iphase[i] = (i < np/2) ? 0 : 1;

        // Gas temperature
        particles.t_gas[i] = temp_dist(gen);

        // Dust fractions for each dust type
        for (int j = 0; j < ndusttypes; ++j) {
            particles.dustfrac[j][i] = dust_dist(gen);
        }
    }

    // Set up dust grain properties
    if (ndusttypes > 0) {
        particles.grainsize[0] = 1.0e-5;  // 0.1 micron
        particles.grainsize[1] = 1.0e-4;  // 1 micron
        particles.graindens[0] = 3.0;     // g/cm³
        particles.graindens[1] = 3.0;     // g/cm³
    }

    // Set up central star
    if (nptmass > 0) {
        particles.xyzmh_ptmass[0][0] = 0.0; // x
        particles.xyzmh_ptmass[1][0] = 0.0; // y
        particles.xyzmh_ptmass[2][0] = 0.0; // z
        particles.xyzmh_ptmass[3][0] = 1.0; // mass in solar masses
        particles.xyzmh_ptmass[4][0] = 0.1; // stellar radius in AU

        particles.vxyz_ptmass[0][0] = 0.0; // vx
        particles.vxyz_ptmass[1][0] = 0.0; // vy
        particles.vxyz_ptmass[2][0] = 0.0; // vz
    }

    // Set up particle type information
    particles.npoftype.resize(ntypes);
    particles.massoftype.resize(ntypes);
    particles.npoftype[0] = np/2;    // gas particles
    particles.npoftype[1] = np/2;    // dust particles
    particles.massoftype[0] = 1e-6;  // mass per gas particle (solar masses)
    particles.massoftype[1] = 1e-8;  // mass per dust particle (solar masses)

    // Physical units (CGS-based)
    particles.hfact = 1.2;           // smoothing length factor
    particles.umass = 1.989e33;      // solar mass in g
    particles.utime = 3.154e7;       // year in s
    particles.udist = 1.496e13;      // AU in cm

    // Energy derivatives (set to zero for this example)
    particles.dudt.resize(np, 0.0);
}

void processResults(const MCFOST::MCFOSTResults& results, int np) {
    std::cout << "\n--- MCFOST Results ---" << std::endl;
    std::cout << "Mean molecular weight: " << results.mu_gas << std::endl;

    // Calculate temperature statistics
    double min_temp = 1e10, max_temp = 0.0, mean_temp = 0.0;
    int valid_particles = 0;

    for (int i = 0; i < np; ++i) {
        if (results.temperature[i] > 0) {
            min_temp = std::min(min_temp, static_cast<double>(results.temperature[i]));
            max_temp = std::max(max_temp, static_cast<double>(results.temperature[i]));
            mean_temp += results.temperature[i];
            valid_particles++;
        }
    }

    if (valid_particles > 0) {
        mean_temp /= valid_particles;
        std::cout << "Temperature statistics:" << std::endl;
        std::cout << "  Valid particles: " << valid_particles << "/" << np << std::endl;
        std::cout << "  Min temperature: " << min_temp << " K" << std::endl;
        std::cout << "  Max temperature: " << max_temp << " K" << std::endl;
        std::cout << "  Mean temperature: " << mean_temp << " K" << std::endl;
    }

    // Calculate packet statistics
    double total_packets = 0.0;
    for (int i = 0; i < np; ++i) {
        total_packets += results.n_packets[i];
    }
    std::cout << "Total photon packets: " << total_packets << std::endl;
    std::cout << "Mean packets per particle: " << total_packets / np << std::endl;
}

/**
 * @brief Example integration with a simple hydrodynamics code
 */
class SimpleHydroCode {
private:
    MCFOST::Interface mcfost_;
    MCFOST::ParticleData particles_;
    bool mcfost_initialized_ = false;

public:
    void initializeMCFOST(const std::string& parameter_file, int ndusttypes) {
        MCFOST::MCFOSTConfig config;
        config.parameter_file = parameter_file;
        config.compute_frad = true;
        config.turn_on_lacc = true;

        mcfost_.initialize(config, ndusttypes);
        mcfost_initialized_ = true;

        std::cout << "MCFOST integrated into hydro code" << std::endl;
    }

    void updateRadiativeHeating() {
        if (!mcfost_initialized_) {
            throw std::runtime_error("MCFOST not initialized");
        }

        MCFOST::MCFOSTConfig config;
        config.parameter_file = "disk.para"; // This should be stored

        // Run MCFOST calculation
        auto results = mcfost_.run(particles_, config);

        // Update particle temperatures and heating rates
        for (size_t i = 0; i < particles_.t_gas.size(); ++i) {
            if (results.temperature[i] > 0) {
                // Update gas temperature based on dust temperature
                // This is a simplified example - real implementation would be more sophisticated
                particles_.t_gas[i] = std::max(particles_.t_gas[i],
                                             static_cast<double>(results.temperature[i]) * 0.8);
            }
        }

        std::cout << "Updated radiative heating for " << particles_.t_gas.size() << " particles" << std::endl;
    }

    void setParticleData(const MCFOST::ParticleData& data) {
        particles_ = data;
    }

    ~SimpleHydroCode() {
        if (mcfost_initialized_) {
            mcfost_.reset();
        }
    }
};

int main() {
    std::cout << "MCFOST C++ Interface Example" << std::endl;
    std::cout << "=============================" << std::endl;

    // Basic usage example
    exampleUsage();

    std::cout << "\n--- Integration Example ---" << std::endl;

    // Example of integration with hydro code
    try {
        SimpleHydroCode hydro;

        // Set up particle data
        auto particles = MCFOST::createParticleData(1000, 1, 1);
        setupExampleParticles(particles, 1000, 1, 1, 2);
        hydro.setParticleData(particles);

        // Initialize MCFOST
        hydro.initializeMCFOST("disk.para", 1);

        // Simulate time evolution with radiative heating
        for (int timestep = 0; timestep < 3; ++timestep) {
            std::cout << "\nTimestep " << timestep + 1 << std::endl;
            hydro.updateRadiativeHeating();
        }

    } catch (const std::exception& e) {
        std::cerr << "Integration example failed: " << e.what() << std::endl;
    }

    return 0;
}
