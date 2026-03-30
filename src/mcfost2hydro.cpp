#ifndef MCFOST_INTERFACE_HPP
#define MCFOST_INTERFACE_HPP

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

// Forward declarations for Fortran interface
extern "C" {
    // Fortran subroutine declarations with C binding
    void init_mcfost_phantom_(
        const char* mcfost_para_filename,
        const int* ndusttypes,
        const bool* use_sph_limits_file,
        const char* sph_limits_file,
        double* sph_limits,
        int* ierr,
        const float* keep_particles,
        const bool* fix_star,
        const bool* turn_on_lacc,
        const bool* turn_on_dust_subl,
        const bool* use_ism_heating,
        int filename_len,
        int limits_file_len
    );

    void run_mcfost_phantom_(
        const int* np,
        const int* nptmass,
        const int* ntypes,
        const int* ndusttypes,
        const int* dustfluidtype,
        const int* npoftype,
        const double* xyzh,
        const double* vxyzu,
        const char* iphase,
        const double* grainsize,
        const double* graindens,
        const double* dustfrac,
        const double* massoftype,
        const double* xyzmh_ptmass,
        const double* vxyz_ptmass,
        const double* hfact,
        const double* umass,
        const double* utime,
        const double* udist,
        const int* ndudt,
        const double* dudt,
        const bool* compute_frad,
        const double* sph_limits,
        float* tphantom,
        float* n_packets,
        double* mu_gas,
        int* ierr,
        const bool* write_t_files,
        const int* ism,
        const double* t_gas
    );

    void reset_mcfost_phantom_();

    void diffusion_opacity_(
        const double* temp,
        const int* icell,
        double* kappa_diffusion
    );
}

namespace MCFOST {

/**
 * @brief Exception class for MCFOST errors
 */
class MCFOSTException : public std::runtime_error {
public:
    explicit MCFOSTException(const std::string& message)
        : std::runtime_error("MCFOST Error: " + message) {}
};

/**
 * @brief Structure to hold particle data for MCFOST calculations
 */
struct ParticleData {
    std::vector<std::vector<double>> xyzh;    // positions and smoothing lengths [4][np]
    std::vector<std::vector<double>> vxyzu;   // velocities and internal energy [4][np]
    std::vector<char> iphase;                 // particle phase information [np]
    std::vector<double> t_gas;                // gas temperature [np]

    // Dust properties
    std::vector<std::vector<double>> dustfrac; // dust fractions [ndusttypes][np]
    std::vector<double> grainsize;             // grain sizes [ndusttypes]
    std::vector<double> graindens;             // grain densities [ndusttypes]

    // Point masses (stars)
    std::vector<std::vector<double>> xyzmh_ptmass; // star positions, masses [5][nptmass]
    std::vector<std::vector<double>> vxyz_ptmass;  // star velocities [3][nptmass]

    // Type information
    std::vector<int> npoftype;      // number of particles of each type
    std::vector<double> massoftype; // mass of each particle type

    // Energy derivatives (for live calculations)
    std::vector<double> dudt;       // energy derivative [ndudt]

    // Physical units
    double hfact;  // smoothing length factor
    double umass;  // mass unit
    double utime;  // time unit
    double udist;  // distance unit
};

/**
 * @brief Structure to hold MCFOST results
 */
struct MCFOSTResults {
    std::vector<float> temperature;  // dust temperature for each particle
    std::vector<float> n_packets;    // number of photon packets per particle
    double mu_gas;                   // mean molecular weight
};

/**
 * @brief Configuration options for MCFOST
 */
struct MCFOSTConfig {
    std::string parameter_file;      // path to MCFOST parameter file
    std::string sph_limits_file;     // path to SPH limits file (optional)
    bool use_sph_limits_file = false;

    // Optional physics settings
    std::optional<float> keep_particles;
    bool fix_star = false;
    bool turn_on_lacc = false;
    bool turn_on_dust_subl = false;
    bool use_ism_heating = false;

    // Calculation options
    bool compute_frad = false;       // compute radiation pressure
    bool write_t_files = false;      // write temperature files
    int ism_heating_mode = 0;        // ISM heating mode (0=off, 1=ProDiMo, 2=Bate&Keto)
};

/**
 * @brief Main MCFOST interface class
 */
class Interface {
private:
    bool initialized_ = false;
    int ndusttypes_ = 0;
    std::vector<double> sph_limits_;

public:
    /**
     * @brief Initialize MCFOST with given configuration
     * @param config Configuration parameters
     * @param ndusttypes Number of dust types
     */
    void initialize(const MCFOSTConfig& config, int ndusttypes) {
        if (initialized_) {
            throw MCFOSTException("MCFOST already initialized");
        }

        ndusttypes_ = ndusttypes;
        sph_limits_.resize(6, 0.0);

        // Prepare optional parameters
        const float* keep_particles_ptr = config.keep_particles ? &(*config.keep_particles) : nullptr;

        int ierr = 0;

        // Call Fortran initialization routine
        init_mcfost_phantom_(
            config.parameter_file.c_str(),
            &ndusttypes_,
            &config.use_sph_limits_file,
            config.sph_limits_file.c_str(),
            sph_limits_.data(),
            &ierr,
            keep_particles_ptr,
            &config.fix_star,
            &config.turn_on_lacc,
            &config.turn_on_dust_subl,
            &config.use_ism_heating,
            static_cast<int>(config.parameter_file.length()),
            static_cast<int>(config.sph_limits_file.length())
        );

        if (ierr != 0) {
            throw MCFOSTException("Initialization failed with error code " + std::to_string(ierr));
        }

        initialized_ = true;
    }

    /**
     * @brief Run MCFOST radiative transfer calculation
     * @param particles Input particle data
     * @param config Configuration options
     * @return MCFOST calculation results
     */
    MCFOSTResults run(const ParticleData& particles, const MCFOSTConfig& config) {
        if (!initialized_) {
            throw MCFOSTException("MCFOST not initialized");
        }

        // Validate input data
        validateParticleData(particles);

        // Prepare input parameters
        int np = particles.xyzh[0].size();
        int nptmass = particles.xyzmh_ptmass.empty() ? 0 : particles.xyzmh_ptmass[0].size();
        int ntypes = particles.npoftype.size();
        int dustfluidtype = 1; // Assuming single dust fluid type
        int ndudt = particles.dudt.size();

        // Flatten 2D arrays for Fortran interface
        std::vector<double> xyzh_flat = flatten2D(particles.xyzh);
        std::vector<double> vxyzu_flat = flatten2D(particles.vxyzu);
        std::vector<double> dustfrac_flat = flatten2D(particles.dustfrac);
        std::vector<double> xyzmh_ptmass_flat = flatten2D(particles.xyzmh_ptmass);
        std::vector<double> vxyz_ptmass_flat = flatten2D(particles.vxyz_ptmass);

        // Prepare output arrays
        MCFOSTResults results;
        results.temperature.resize(np, -1.0f);
        results.n_packets.resize(np, 0.0f);

        int ierr = 0;

        // Call Fortran computation routine
        run_mcfost_phantom_(
            &np,
            &nptmass,
            &ntypes,
            &ndusttypes_,
            &dustfluidtype,
            particles.npoftype.data(),
            xyzh_flat.data(),
            vxyzu_flat.data(),
            particles.iphase.data(),
            particles.grainsize.data(),
            particles.graindens.data(),
            dustfrac_flat.data(),
            particles.massoftype.data(),
            xyzmh_ptmass_flat.data(),
            vxyz_ptmass_flat.data(),
            &particles.hfact,
            &particles.umass,
            &particles.utime,
            &particles.udist,
            &ndudt,
            particles.dudt.data(),
            &config.compute_frad,
            sph_limits_.data(),
            results.temperature.data(),
            results.n_packets.data(),
            &results.mu_gas,
            &ierr,
            &config.write_t_files,
            &config.ism_heating_mode,
            particles.t_gas.data()
        );

        if (ierr != 0) {
            throw MCFOSTException("Calculation failed with error code " + std::to_string(ierr));
        }

        return results;
    }

    /**
     * @brief Reset MCFOST memory state
     */
    void reset() {
        if (initialized_) {
            reset_mcfost_phantom_();
        }
    }

    /**
     * @brief Calculate diffusion opacity for a given temperature and cell
     * @param temperature Temperature in K
     * @param icell Cell index
     * @return Diffusion opacity in cm^2/g
     */
    double getDiffusionOpacity(double temperature, int icell) {
        if (!initialized_) {
            throw MCFOSTException("MCFOST not initialized");
        }

        double kappa_diffusion;
        diffusion_opacity_(&temperature, &icell, &kappa_diffusion);
        return kappa_diffusion;
    }

    /**
     * @brief Destructor - cleanup MCFOST state
     */
    ~Interface() {
        if (initialized_) {
            reset();
        }
    }

private:
    /**
     * @brief Validate particle data consistency
     */
    void validateParticleData(const ParticleData& particles) {
        if (particles.xyzh.size() != 4) {
            throw MCFOSTException("xyzh must have 4 components");
        }
        if (particles.vxyzu.size() != 4) {
            throw MCFOSTException("vxyzu must have 4 components");
        }

        size_t np = particles.xyzh[0].size();
        for (const auto& component : particles.xyzh) {
            if (component.size() != np) {
                throw MCFOSTException("Inconsistent particle array sizes in xyzh");
            }
        }
        for (const auto& component : particles.vxyzu) {
            if (component.size() != np) {
                throw MCFOSTException("Inconsistent particle array sizes in vxyzu");
            }
        }

        if (particles.iphase.size() != np) {
            throw MCFOSTException("iphase size mismatch");
        }
        if (particles.t_gas.size() != np) {
            throw MCFOSTException("t_gas size mismatch");
        }

        if (ndusttypes_ > 0) {
            if (particles.dustfrac.size() != ndusttypes_) {
                throw MCFOSTException("dustfrac first dimension must equal ndusttypes");
            }
            for (const auto& dust_component : particles.dustfrac) {
                if (dust_component.size() != np) {
                    throw MCFOSTException("dustfrac particle dimension mismatch");
                }
            }
        }
    }

    /**
     * @brief Flatten 2D vector for Fortran interface (column-major order)
     */
    std::vector<double> flatten2D(const std::vector<std::vector<double>>& vec2d) {
        if (vec2d.empty()) return {};

        size_t rows = vec2d.size();
        size_t cols = vec2d[0].size();
        std::vector<double> flat(rows * cols);

        // Fortran uses column-major order
        for (size_t j = 0; j < cols; ++j) {
            for (size_t i = 0; i < rows; ++i) {
                flat[j * rows + i] = vec2d[i][j];
            }
        }
        return flat;
    }
};

/**
 * @brief Convenience function to create particle data structure
 */
ParticleData createParticleData(int np, int ndusttypes = 0, int nptmass = 0) {
    ParticleData data;

    // Initialize particle arrays
    data.xyzh.resize(4);
    data.vxyzu.resize(4);
    for (int i = 0; i < 4; ++i) {
        data.xyzh[i].resize(np);
        data.vxyzu[i].resize(np);
    }

    data.iphase.resize(np);
    data.t_gas.resize(np);

    // Initialize dust arrays if needed
    if (ndusttypes > 0) {
        data.dustfrac.resize(ndusttypes);
        for (int i = 0; i < ndusttypes; ++i) {
            data.dustfrac[i].resize(np);
        }
        data.grainsize.resize(ndusttypes);
        data.graindens.resize(ndusttypes);
    }

    // Initialize star arrays if needed
    if (nptmass > 0) {
        data.xyzmh_ptmass.resize(5);
        data.vxyz_ptmass.resize(3);
        for (int i = 0; i < 5; ++i) {
            data.xyzmh_ptmass[i].resize(nptmass);
        }
        for (int i = 0; i < 3; ++i) {
            data.vxyz_ptmass[i].resize(nptmass);
        }
    }

    return data;
}

} // namespace MCFOST

#endif // MCFOST_INTERFACE_HPP
