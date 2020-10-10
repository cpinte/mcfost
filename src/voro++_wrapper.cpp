#include <voro++/voro++.hh>
#include <iostream>
#include <iomanip>
using namespace voro;

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


extern "C" {
  void voro_C(int n, int max_neighbours, double limits[6], double x[], double y[], double z[], double h[],  double threshold, int n_vectors, double cutting_vectors[][3], double cutting_distance_o_h, int icell_start, int icell_end, int cpu_id, int n_cpu, int n_points_per_cpu,
              int &n_in, double volume[], double delta_edge[], double delta_centroid[], int first_neighbours[], int last_neighbours[], int n_neighbours[], int neighbours_list[], bool was_cell_cut[],int &ierr) {

    ierr = 0 ;

    double ax,bx, ay,by, az, bz;
    ax = limits[0]; bx = limits[1];
    ay = limits[2]; by = limits[3];
    az = limits[4]; bz = limits[5];

    int i, k, nx,ny,nz,init_mem(8);

    int row = n_points_per_cpu * max_neighbours;

    //wall_list wl; // I am not adding extra walls yet
    bool xperiodic,yperiodic,zperiodic;

    xperiodic = false;
    yperiodic = false;
    zperiodic = false;

    // Makes a guess at the optimal grid of blocks to use, as in voro++ pre_container
    double dx=bx-ax, dy=by-ay, dz=bz-az;
    double optimal_particles = 5.6;
    double ilscale=pow(n/(optimal_particles*dx*dy*dz),1/3.0);
    nx=int(dx*ilscale+1);
    ny=int(dy*ilscale+1);
    nz=int(dz*ilscale+1);

    // Define the container
    particle_order po;
    container con(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
    for(i=0;i<n;i++) {
      con.put(po,i,x[i],y[i],z[i]);
    }

    // Perform the Voronoi tesselation
    c_loop_order vlo(con,po);

    voronoicell_neighbor c(con);
    int pid;
    std::vector<int> vi;

    int n_neighbours_cell, first_neighbour, last_neighbour;
    int max_size_list = max_neighbours * n;
    double cx, cy, cz, cutting_distance;

    n_neighbours[cpu_id] = 0;
    last_neighbour = -1;
    n_in = 0;

    float progress = 0.0;
    float progress_bar_step = 0.01;
    float pb_threshold = progress_bar_step*(1.*n)/n_cpu;

    if (!vlo.start()) {
      std::cout << "Error : voro++ did not manage to initialize" << std::endl;
      std::cout << "Exiting" << std::endl;
      ierr = 1;
      exit(1);
    }

    do {
      pid = vlo.pid(); // id of the current cell in the c_loop

      if ((pid >= icell_start) && (pid <= icell_end)) { // We only compute a fraction of the cells in a given thread
        if (con.compute_cell(c,vlo)) { // return false if the cell was removed
          n_in++;

          if (cpu_id == n_cpu-1) {
            if (n_in > pb_threshold) {
              progress  += progress_bar_step;
              pb_threshold += progress_bar_step*(1.*n)/n_cpu;
              progress_bar(progress);
            }
          }

          // Store the neighbours list
          n_neighbours_cell = c.number_of_faces();
          n_neighbours[cpu_id] = n_neighbours[cpu_id] + n_neighbours_cell;

          first_neighbour = last_neighbour+1; first_neighbours[pid] = first_neighbour;
          last_neighbour  = last_neighbour + n_neighbours_cell; last_neighbours[pid]  = last_neighbour;

          if (n_neighbours[cpu_id] > max_size_list) {ierr = 1; exit(1);}

          c.neighbors(vi);
          for (i=0; i<n_neighbours_cell; i++) {
            if (vi[i] >=0) {
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i] + 1;
            } else { // Wall
              neighbours_list[row * cpu_id + first_neighbour+i] = vi[i];
            }
          }

          // Compute the maximum distance to a vertex and the distance to the centroid from the particule position
          delta_edge[pid] = sqrt(c.max_radius_squared()); // does not cost anything

          c.centroid(cx,cy,cz);
          delta_centroid[pid] = sqrt(cx*cx + cy*cy + cz*cz);

          // If the Voronoi cell is elongated, we intersect it with a dodecahedron
          was_cell_cut[pid] = false;
          if (delta_edge[pid] > threshold * h[pid]) {
            cutting_distance = cutting_distance_o_h * h[pid];

            // Adding the n-planes to cut the cell
            for (k=0; k<n_vectors; k++) {
              // All the planes are a distance of 0.5 x vector length away from the center of the cell --> extra factor 2
              c.plane(cutting_vectors[k][0], cutting_vectors[k][1], cutting_vectors[k][2], 2*cutting_distance);
            }
            was_cell_cut[pid] = true;
          }
          // Volume of the cell (computed after eventual cut)
          volume[pid] = c.volume();
        }  // con.compute_cell
      } // pid test
    } while(vlo.inc()); //Finds the next particle to test

    if (cpu_id == n_cpu-1) progress_bar(1.0);
  }
}
