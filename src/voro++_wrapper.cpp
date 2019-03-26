/*
 C Wrapper of voro++ C++ libraries for mcfost
*/
#include "voro++.hh"
#include <iostream>
#include <iomanip>

using namespace voro;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void progress_bar(float progress) {

  int barWidth = 50;

  std::cout  << " " << std::setfill(' ') << std::setw(3) << int(progress * 100.0) << "%";
  std::cout << " |";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << " "; //">";
    else std::cout << " ";
  }
  std::cout << "| \r";
  std::cout.flush();

  if (progress >= 1.0) std::cout << std::endl;
}
// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

extern "C" {
  void voro_C(int n, int max_neighbours, double limits[6], double x[], double y[], double z[], double h[],  double threshold, int n_vectors, double cutting_vectors[][3], double cutting_distance_o_h, int icell_start, int icell_end, int cpu_id, int n_cpu, int n_points_per_cpu,
              int &n_in, double volume[], double delta_edge[], double delta_centroid[], int first_neighbours[], int last_neighbours[], int n_neighbours[], int neighbours_list[], bool was_cell_cut[],int &ierr,
              int n_stars, double stars_radius[], int stars_cell[]) {

    /*
     voro_C is called by mcfost via call voro()

     input int n 			: number of cells inside the model limits
     input int max_neighbours		: maxium number neighbours per cell
     input double limits		: array of size 6 containing the upper and lower limits of the model in each direction
     input double x, y, z 		: array of size n, spatial coordinates of the grid
     input double h 			: array of size n, cells larger than threshold * h are cut up to 2 * cutting_distance_o_h * h
     input double threshold		: set the maximum elongation of cell beyond which it is cut
     input in n_vectors 		: number of faces of the polyhedron used to cut cells
     input double cutting_vectors 	: array of size (n_vectors, 3) defining each face of the polyhedron used to cut cells
     input double cutting_distance_o_h  : distance of a face to the center of the polyhedron divided by h
     input int icell_start, icell_end   : For each thread call of voro_C, the loop over particles starts and ends at icell_start, icell_end
     input int cpu_id			: id of the current thread. 0 if no OpenMP
     input integer n_cpu		: Number of cores for OpenMP (OPENMP_NUM_THREADS)
     output integer n_in		: Number of cells that were part of the tesselation
     output double volume		: array of size n containing the volume of each cell
     output double delta_edge		: array of size n containing the maximum distance to a vertex for each cell
     output double delta_centroid	: array of size n containing the distance of particule to the center of a voronoi cell
     output integer first_neighbours    : array of size n, for each cell stores its first neighbour
     output integer last_neighbours     : array of size n, for each cell stores its last neighbour
     output integer n_neighbours	: size(n), total number of neighbours for each thread
     output integer neighbours_list	: size(n*max_neihbours*n_cpu), contains the neighbours for each cell
     output bool was_cell_cut		: size(n), array of .true. value is the cell has been cut because too elongated
     input int n_points_per_cpu		: number of cells per cpu, n_cells if only one cpu
     output integer ierr		: a integer for passive exit error message. > 0 if any problem during the tesselation
     input integer n_stars		: Number of stars in the model
     input int stars_radius		: size(n_stars) array of stellar radius for each star
     input int stars_cell		: size(n_stars) array of stars index on the whole grid
    */

    ierr = 0 ;

    // intent(in)
    double ax,bx, ay,by, az, bz ;
    ax = limits[0] ; bx = limits[1] ;
    ay = limits[2] ; by = limits[3] ;
    az = limits[4] ; bz = limits[5] ;

    int i, k, nx,ny,nz,init_mem(8);
    int istar, ijk;
    double rx, ry, rz, rs,rsq,r;

    int row = n_points_per_cpu * max_neighbours ;

    //wall_list wl; // I am not adding extra walls yet
    bool xperiodic,yperiodic,zperiodic ;
    pre_container *pcon=NULL; // pre-container pointer

    xperiodic = false ;
    yperiodic = false ;
    zperiodic = false ;

    // Define a pre-container to determine the optimal block size
    pcon=new pre_container(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
    for(i=0;i<n;i++) {
      //std::cout << i << " " << x[i] << std::endl;
      pcon->put(i,x[i],y[i],z[i]);
    }
    pcon->guess_optimal(nx,ny,nz);

    // define the proper container and point the pre-containet toward it
    particle_order vo;
    container con(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);

    //con.add_wall(wl);
    pcon->setup(vo,con);
    delete pcon;

    c_loop_order vlo(con,vo);
    c_loop_subset cls(con);


    // Perform the Voronoi tesselation
    voronoicell_neighbor c(con);
    int pid ; //particle index
    std::vector<int> vi;

    int n_neighbours_cell, first_neighbour, last_neighbour ;
    int max_size_list = max_neighbours * n ;
    double cx, cy, cz, cutting_distance ;


    n_neighbours[cpu_id] = 0 ;
    last_neighbour = -1 ;
    n_in = 0 ;

    float progress = 0.0 ;
    float progress_bar_step = 0.01 ;
    float pb_threshold = progress_bar_step*(1.*n)/n_cpu ;

    if (!vlo.start()) {
      std::cout << "Error : voro++ did not manage to initialize" << std::endl;
      std::cout << "Exiting" << std::endl;
      ierr = 1 ;
      exit(1) ;
    }

    // First loop on stars
    ijk=0;
    for (istar=0; istar < n_stars; istar++) {
      cls.setup_sphere(0., 0., 0., 2.001*stars_radius[istar], false);
      rs = stars_radius[istar];
      /* Search cells which are neighbours of the star etoile(istar)
         And cut them at the stellar surface.
         Their volume will be computed later.
         Store a specific arrays if a cell has been cut at the stellar surface.
         The later is needed in cross_voronoi_cell() to compute the proper ray path
         in a cell. */

      //if (con.find_voronoi_cell(0.0, 0.0, 0.0, cx, cy, cz, pid))
      //printf("%d %f %f %f\n", pid, cx, cy, cz);

      if (cls.start()) do {
          if (cls.pid()==stars_cell[istar] && con.compute_cell(c, cls)) {
            puts("Cutting star cell");
            c.init(-rs,rs,-rs,rs,-rs,rs);
            //c.plane(0,Phi,1,rs);c.plane(0,-Phi,1,rs);c.plane(0,Phi,-1,rs);
            //c.plane(0,-Phi,-1,rs);c.plane(1,0,Phi,rs);c.plane(-1,0,Phi,rs);
            //c.plane(1,0,-Phi,rs);c.plane(-1,0,-Phi,rs);c.plane(Phi,1,0,rs);
            //c.plane(-Phi,1,0,rs);c.plane(Phi,-1,0,rs);c.plane(-Phi,-1,0,rs);
           for(i=0;i<250;i++) {
             rx=rs*(2*rnd()-1);
             ry=rs*(2*rnd()-1);
             rz=rs*(2*rnd()-1);
             rsq=rx*rx+ry*ry+rz*rz;
             if(rsq>0.&&rsq<rs*rs) {
               r=1/sqrt(rsq);rx*=r;ry*=r;rz*=r;
               c.plane(rx,ry,rz,rs);
             }
           }
          }
        } while(cls.inc());
    } //stars loop
    if (ijk) printf("Cut %d cells which overlap the stars\n", ijk);

    //float V_old, V_new ;
    do {
      pid = vlo.pid() ; // id of the current cell in the c_loop

      if ((pid >= icell_start) && (pid <= icell_end)) { // We only compute a fraction of the cells in a given thread
        if (con.compute_cell(c,vlo)) { // return false if the cell was removed
          n_in++ ;

          if (cpu_id == n_cpu-1) {
            if (n_in > pb_threshold) {
              progress  += progress_bar_step ;
              pb_threshold += progress_bar_step*(1.*n)/n_cpu ;
              progress_bar(progress) ;
            }
          }

          // Store the neighbours list
          n_neighbours_cell = c.number_of_faces() ;
          n_neighbours[cpu_id] = n_neighbours[cpu_id] + n_neighbours_cell ;

          first_neighbour = last_neighbour+1 ; first_neighbours[pid] = first_neighbour ;
          last_neighbour  = last_neighbour + n_neighbours_cell ; last_neighbours[pid]  = last_neighbour ;

          if (n_neighbours[cpu_id] > max_size_list) {ierr = 1 ; exit(1) ;}

          c.neighbors(vi) ;
          for (i=0 ; i<n_neighbours_cell ; i++) {
            //std::cout << "test0 " << row  << " " << n_points_per_cpu << " " << max_neighbours << std::endl;
            if (vi[i] >=0) {
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i] + 1 ;
            } else {
              // Wall
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i] ;
            }
          }

          // Compute the maximum distance to a vertex and the distance to the centroid from the particule position
          delta_edge[pid] = sqrt(c.max_radius_squared()) ; // does not cost anything

          c.centroid(cx,cy,cz) ;
          delta_centroid[pid] = sqrt(cx*cx + cy*cy + cz*cz) ;


          // If the Voronoi cell is elongated, we intersect it with a dodecahedron
          was_cell_cut[pid] = false ;
          if (delta_edge[pid] > threshold * h[pid]) {
            cutting_distance = cutting_distance_o_h * h[pid] ;

            //V_old = c.volume() ;

            // Adding the n-planes to cut the cell
            for (k=0 ; k<n_vectors ; k++) {
              // All the planes are a distance of 0.5 x vector length away from the center of the cell --> extra factor 2
              c.plane(cutting_vectors[k][0], cutting_vectors[k][1], cutting_vectors[k][2], 2*cutting_distance) ;
            }
            was_cell_cut[pid] = true ;

            //V_new = c.volume() ;
            // std::cout << "e/3h= " << delta_edge[pid] / ( threshold * h[pid])   << " V_ratio=" << V_new / V_old << std::endl ; // Ok, numbers make sense
          }
          // Volume of the cell (computed after eventual cut)
          volume[pid] = c.volume() ;
        }  // con.compute_cell
      } // pid test

    } while(vlo.inc()); //Finds the next particle to test

    if (cpu_id == n_cpu-1) progress_bar(1.0) ;
    con.draw_particles("Voro.gnu");
  }
}
